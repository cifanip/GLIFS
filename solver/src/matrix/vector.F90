module vector_mod
  
  use matrix_mod
  
  implicit none

  type, extends(matrix), public :: vector
    private
    
    !number of proc row
    integer, public :: nprow

    !proc row (grid coordinates)
    integer, public :: prow
    
    !local number of rows
    integer, public :: nrow

    !block size
    integer, public :: nblk
    
    !vector descriptor
    integer, dimension(7), public :: desc
    
    complex(double_p), allocatable, dimension(:), public :: v
    
    contains
    
    procedure :: vector_ctor_opt_pgrid
    procedure :: vector_ctor_set_pgrid
    generic, public :: ctor => vector_ctor_opt_pgrid,&
                               vector_ctor_set_pgrid
    procedure, public :: delete
    procedure, public :: default_parameters
    procedure :: init_descriptor => init_vector_descriptor
    procedure :: set_local_size
    procedure :: set_block_size
    procedure, public :: write_to_disk

  end type
  
  private :: vector_ctor_opt_pgrid,&
             vector_ctor_set_pgrid,&
             delete,&
             default_parameters,&
             init_vector_descriptor,&
             set_local_size,&
             set_block_size,&
             write_to_disk

contains

!========================================================================================!
  subroutine delete(this,delete_mpi)
    class(vector), intent(inout) :: this
    logical, intent(in), optional :: delete_mpi
    
    call this%default_parameters()
    call deallocate_array(this%v)
    call this%delete_matrix(delete_mpi)
        
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine vector_ctor_opt_pgrid(this,mpic,fname,n)
    class(vector), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    character(len=*), intent(in) :: fname
    integer, intent(in) :: n
    integer :: tmp
    
    call this%init_matrix(mpic,fname,n)
    
    call set_1d_pgrid(this%n,mpic%nprocs,this%nprow)
    
    call set_blacs_pgrid(1,this%nprow,tmp,this%prow,this%ctxt)
    
    if (this%ctxt==-1) then
      this%is_allocated = .FALSE.
    else
      this%is_allocated = .TRUE.
    end if

    call init_aux_comm(this%is_allocated,this%comm,this%rank,this%nprocs)

    if (.not.this%is_allocated) then
      call this%default_parameters()
      this%desc(2)=-1
      call allocate_array(this%v,1,this%nrow)
      return
    end if
    
    call this%set_block_size()
    call this%set_local_size()
    
    call this%init_descriptor()
    
    call allocate_array(this%v,1,this%nrow)
    this%v=(0.d0,0.d0)
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine vector_ctor_set_pgrid(this,mpic,fname,n,nprow)
    class(vector), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    character(len=*), intent(in) :: fname
    integer, intent(in) :: n,nprow
    integer :: tmp
    
    call this%init_matrix(mpic,fname,n)
    
    this%nprow = nprow
    
    call check_pgrid(this%n,this%nprow,1,BLOCK_ROW)
    
    call set_blacs_pgrid(1,this%nprow,tmp,this%prow,this%ctxt)

    if (this%ctxt==-1) then
      this%is_allocated = .FALSE.
    else
      this%is_allocated = .TRUE.
    end if

    call init_aux_comm(this%is_allocated,this%comm,this%rank,this%nprocs)

    if (.not.this%is_allocated) then
      call this%default_parameters()
      this%desc(2)=-1
      call allocate_array(this%v,1,this%nrow)
      return
    end if
    
    call this%set_block_size()
    call this%set_local_size()
    
    call this%init_descriptor()
    
    call allocate_array(this%v,1,this%nrow)
    this%v=(0.d0,0.d0)
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine default_parameters(this)
    class(vector), intent(inout) :: this
    
    this%nrow  = -1
    this%prow  = -1

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_vector_descriptor(this)
    class(vector), intent(inout) :: this

    !set descriptor
    this%desc(1) = 502
    this%desc(2) = this%ctxt
    this%desc(3) = this%n
    this%desc(4) = this%nblk
    this%desc(5) = 0
    this%desc(6) = this%lld

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine set_block_size(this)
    class(vector), intent(inout) :: this
    
    this%nblk = ceiling(real(this%n)/real(this%nprow))

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine set_local_size(this)
    class(vector), intent(inout) :: this

    this%nrow = NUMROC(this%n,this%nblk,this%prow,0,this%nprow)
    
    !local leading dimension
    this%lld = this%nrow

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine write_to_disk(this,time_level,mat_dir)
    class(vector), intent(in) :: this
    integer, intent(in) :: time_level
    character(len=*), intent(in) :: mat_dir
    character(len=20) :: time_dir
    character(len=:), allocatable :: path
    integer, allocatable, dimension(:) :: buff_size
    integer, dimension(2) :: grid_info,idx
    integer :: status(MPI_STATUS_SIZE)
    integer(mpi_offset_kind), allocatable, dimension(:) :: disp
    integer :: f_handle,err,local_buff_size,prow,nprow,i,&
               mpi_value_type,data_size
    
    mpi_value_type = MPI_DOUBLE_COMPLEX
    data_size = complex_dp_size
    
    if (this%is_allocated) then
    
      write(time_dir,int_format_) time_level
      path=trim(adjustl(time_dir))//'/'//mat_dir//'/'//this%fname
    
      prow=this%prow
      nprow=this%nprow
    
      grid_info(1)=this%n
      grid_info(2)=nprow

      idx(1)=get_gidx(1,this%nblk,prow,nprow)
      idx(2)=get_gidx(this%nrow,this%nblk,prow,nprow)
    
      !gather local buff size
      call allocate_array(buff_size,0,nprow-1)
      local_buff_size=size(this%v)
      call MPI_ALLGATHER(local_buff_size, 1, MPI_INTEGER, buff_size, 1,&
                         MPI_INTEGER, this%comm, err)
    
      !compute displacement
      allocate(disp(0:nprow-1))
      disp(0)=0
      do i=1,nprow-1
        if (i==1) then
          disp(i)=data_size*buff_size(i-1)+&
                  integer_size*(size(idx)+size(grid_info))
        else
          disp(i)=disp(i-1)+data_size*buff_size(i-1)+integer_size*size(idx)
        end if
      end do

      call mpi_file_open(this%comm,path,MPI_MODE_CREATE+MPI_MODE_WRONLY,&
                         MPI_INFO_NULL,f_handle,err)
                       
      if (err.ne.MPI_SUCCESS) then
        call abort_run('Write file '//this%fname//': mpi file open failed',err)
      end if

      call mpi_file_seek(f_handle,disp(prow),MPI_SEEK_SET,err)
      if (err.ne.MPI_SUCCESS) then
        call abort_run('Write file '//this%fname//': mpi file seek failed',err)
      end if
    
      if (IS_MASTER) then
        call mpi_file_write(f_handle,grid_info,size(grid_info),MPI_INTEGER,status,err)
        if (err.ne.MPI_SUCCESS) then
          call abort_run('Write file '//this%fname//&
                         ': mpi file write info_grid failed',err)
        end if
      end if

      call mpi_file_write(f_handle,idx,size(idx),MPI_INTEGER,status,err)
      if (err.ne.MPI_SUCCESS) then
        call abort_run('Write file '//this%fname//': mpi file write idx failed',err)
      end if

      call mpi_file_write(f_handle,this%v,size(this%v),mpi_value_type,status,err)
      if (err.ne.MPI_SUCCESS) then
        call abort_run('Write file '//this%fname//': mpi file write field failed',err)
      end if

      call mpi_file_close(f_handle,err)
      if (err.ne.MPI_SUCCESS) then
        call abort_run('Write file '//this%fname//': mpi file close failed',err)
      end if
  
      call MPI_BARRIER(this%comm,err)
      
    end if

  end subroutine
!========================================================================================!

end module vector_mod