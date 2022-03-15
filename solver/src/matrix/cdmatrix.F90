module cdmatrix_mod
  
  use dmatrix_mod
  
  implicit none

  type, extends(dmatrix), public :: cdmatrix
    private

    !matrix (local)
    complex(double_p), allocatable, dimension(:,:), public :: m
    
    !global indexes
    integer, allocatable, dimension(:,:,:), public :: gidx
    
    !type pointer
    type(cdmatrix), pointer, public :: ptrm => NULL() 
    
    contains

    procedure, public :: cdmatrix_ctor_opt_pgrid
    procedure, public :: cdmatrix_ctor_set_pgrid
    generic, public :: ctor => cdmatrix_ctor_opt_pgrid,&
                               cdmatrix_ctor_set_pgrid
    
    procedure, public :: delete => delete_cdmatrix
    procedure, public :: cyclic_to_column
    procedure, public :: column_to_cyclic
    procedure, public :: write_to_disk
    procedure, public :: read_from_disk
    procedure, public :: make_skewh
    procedure, public :: hconj
    procedure, public :: copy_values
    procedure, public :: multiply
    procedure, public :: multiply_h
    procedure, public :: inf_norm
    procedure, public :: l2_norm
    procedure, public :: subtract
    procedure, public :: sum
    procedure, public :: set_to_zero
    procedure, public :: is_stored
    procedure :: set_global_indexes
    procedure, public :: allocate_ptrm
    procedure :: s_mul
    
    procedure :: equal
    generic, public :: assignment(=) => equal

  end type
  
  private :: cdmatrix_ctor_opt_pgrid,&
             cdmatrix_ctor_set_pgrid,&
             delete_cdmatrix,&
             cyclic_to_column,&
             column_to_cyclic,&
             write_to_disk,&
             read_from_disk,&
             equal,&
             make_skewh,&
             hconj,&
             copy_values,&
             multiply,&
             inf_norm,&
             l2_norm,&
             subtract,&
             sum,&
             set_to_zero,&
             set_global_indexes,&
             is_stored,&
             allocate_ptrm,&
             s_mul
             
  interface
    real(double_p) function pzlange(norm,m,n,a,ia,ja,desca,work)
      import
      character :: norm
      integer :: m,n,ia,ja
      integer :: desca(*)
      complex(double_p) :: a(*)
      real(double_p) :: work(*)
    end function pzlange
  end interface
  
  !out of type
  private :: check_grid_info
  private :: reset_uppert

contains

!========================================================================================!
  recursive subroutine delete_cdmatrix(this,delete_mpi)
    class(cdmatrix), intent(inout) :: this
    logical, intent(in), optional :: delete_mpi
 
    call deallocate_array(this%m)
    call deallocate_array(this%gidx)
    if (associated(this%ptrm)) then
      call this%ptrm%delete(delete_mpi=.FALSE.)
    end if
    call this%delete_dmatrix(delete_mpi)
        
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine equal(lhs,rhs)
    class(cdmatrix), intent(inout) :: lhs
    type(cdmatrix), intent(in) :: rhs
    
    if (lhs%is_allocated) then
      call abort_run('lhs already allocated found in assignment(=) cdmatrix')
    end if
    call allocate_array(lhs%m,1,rhs%nrow,1,rhs%ncol)
    call allocate_array(lhs%gidx,1,2,1,rhs%nrow,1,rhs%ncol)
    
    if (rhs%is_allocated) then
      call copy_array(rhs%m,lhs%m)
      call copy_array(rhs%gidx,lhs%gidx)
    end if
    
    call lhs%equal_dmatrix(rhs)
    
    call lhs%allocate_ptrm()
    
  end subroutine   
!========================================================================================!

!========================================================================================!
  subroutine cdmatrix_ctor_opt_pgrid(this,mpic,distr_type,fname,read_pgrid,n)
    class(cdmatrix), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    integer, intent(in) :: distr_type
    character(len=*), intent(in) :: fname
    logical, intent(in) :: read_pgrid
    integer, intent(in), optional :: n
    
    call this%init_dmatrix_opt_pgrid(mpic,distr_type,fname,read_pgrid,n)

    call allocate_array(this%m,1,this%nrow,1,this%ncol)
    call allocate_array(this%gidx,1,2,1,this%nrow,1,this%ncol)
    if (this%is_allocated) then
      this%m=(0.d0,0.d0)
      call this%set_global_indexes()
      call this%allocate_ptrm()
    end if
     
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine cdmatrix_ctor_set_pgrid(this,mpic,distr_type,fname,n,nprow,npcol)
    class(cdmatrix), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    integer, intent(in) :: distr_type
    character(len=*), intent(in) :: fname
    integer, intent(in) :: n,nprow,npcol
    
    call this%init_dmatrix_set_pgrid(mpic,distr_type,fname,n,nprow,npcol)

    call allocate_array(this%m,1,this%nrow,1,this%ncol)
    call allocate_array(this%gidx,1,2,1,this%nrow,1,this%ncol)
    if (this%is_allocated) then
      this%m=(0.d0,0.d0)
      call this%set_global_indexes()
      call this%allocate_ptrm()
    end if
     
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine set_global_indexes(this)
    class(cdmatrix), intent(inout) :: this
    integer :: i,j,nrow,ncol,nrblk,prow,nprow,ncblk,pcol,npcol
    
    if (this%is_allocated) then
    
      nrow  = this%nrow
      ncol  = this%ncol
      nrblk = this%nrblk
      ncblk = this%ncblk
      prow  = this%prow
      pcol  = this%pcol
      nprow = this%nprow
      npcol = this%npcol

      do j=1,ncol
        do i=1,nrow
          this%gidx(1,i,j) = get_gidx(i,nrblk,prow,nprow)
          this%gidx(2,i,j) = get_gidx(j,ncblk,pcol,npcol)
        end do
      end do

    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine cyclic_to_column(this,q)
    class(cdmatrix), intent(in) :: this
    type(cdmatrix), intent(inout) :: q
    integer :: n

    if ((.not.this%is_allocated).AND.(.not.q%is_allocated)) then
      return
    end if

    if (this%dtype.ne.BLOCK_CYCLIC) then
      call abort_run('Input matrix not BLOCK_CYLCIC in cyclic_to_column')
    end if

    if (q%dtype.ne.BLOCK_COLUMN) then
      call abort_run('Input matrix not BLOCK_COLUMN in cyclic_to_column')
    end if
    
    n=this%n

    if (this%nprocs >= q%nprocs) then
      call pzgemr2d(n,n,this%m,1,1,this%desc,q%m,1,1,q%desc,this%ctxt)
    else
      call pzgemr2d(n,n,this%m,1,1,this%desc,q%m,1,1,q%desc,q%ctxt)
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine column_to_cyclic(this,q)
    class(cdmatrix), intent(in) :: this
    type(cdmatrix), intent(inout) :: q
    integer :: n

    if ((.not.this%is_allocated).AND.(.not.q%is_allocated)) then
      return
    end if
    
    if (this%dtype.ne.BLOCK_COLUMN) then
      call abort_run('Input matrix not BLOCK_COLUMN in column_to_cyclic')
    end if

    if (q%dtype.ne.BLOCK_CYCLIC) then
      call abort_run('Output matrix not BLOCK_CYLCIC in column_to_cyclic')
    end if
    
    n=this%n

    if (this%nprocs >= q%nprocs) then
      call pzgemr2d(n,n,this%m,1,1,this%desc,q%m,1,1,q%desc,this%ctxt)
    else
      call pzgemr2d(n,n,this%m,1,1,this%desc,q%m,1,1,q%desc,q%ctxt)
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine make_skewh(this)
    class(cdmatrix), intent(inout) :: this
    complex(double_p) :: alpha,beta
    
    if (.not.this%is_allocated) then
      return
    end if

    call reset_uppert(this%gidx,this%m)

    alpha=(-1.d0,0.0)
    beta=(1.d0,0.d0)
    call copy_values(this%ptrm,this)
     
    call set_num_threads()

    call PZTRANC(this%n,this%n,alpha,this%ptrm%m,1,1,this%ptrm%desc,beta,&
                 this%m,1,1,this%desc)
     
    call set_sigle_thread()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine hconj(this,a)
    class(cdmatrix), intent(inout) :: this
    type(cdmatrix), intent(in) :: a
    complex(double_p) :: alpha,beta

    if (.not.this%is_allocated) then
      return
    end if
    
    alpha=(1.d0,0.0)
    beta=(0.d0,0.d0)
    
    call set_num_threads()
    
    call PZTRANC(this%n,this%n,alpha,a%m,1,1,a%desc,beta,this%m,1,1,this%desc) 
    
    call set_sigle_thread()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine multiply(this,a,b)
    class(cdmatrix), intent(inout) :: this
    type(cdmatrix), intent(in) :: a,b
    complex(double_p) :: alpha,beta
    integer :: m,n,k

    if (.not.this%is_allocated) then
      return
    end if
    
    m=this%n
    n=m
    k=m
    
    alpha=(1.d0,0.d0)
    beta=(0.d0,0.d0)
    
    call set_num_threads()

    call pzgemm('N','N',m,n,k,alpha,a%m,1,1,a%desc,b%m,1,1,b%desc,beta,&
                 this%m,1,1,this%desc)

    call set_sigle_thread()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine multiply_h(this,a,b)
    class(cdmatrix), intent(inout) :: this,a
    type(cdmatrix), intent(in) :: b
    complex(double_p) :: alpha,beta,s
    integer :: m,n,k

    if (.not.this%is_allocated) then
      return
    end if
    
    m=this%n
    n=m
    k=m
    
    alpha=(0.d0,-1.d0)
    beta=(0.d0,0.d0)
    
    a%ptrm%m=a%m
    s = (0.d0,1.d0)
    call a%s_mul(s)
    
    call set_num_threads()

    call pzhemm('R','L',m,n,alpha,a%m,1,1,a%desc,b%m,1,1,b%desc,beta,&
                this%m,1,1,this%desc)

    call set_sigle_thread()
    
    a%m=a%ptrm%m

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine s_mul(this,s)
    class(cdmatrix), intent(inout) :: this
    complex(double_p), intent(in) :: s
    integer :: i,j,nrow,ncol
    
    if (.not.this%is_allocated) then
      return
    end if
    
    nrow = this%nrow
    ncol = this%ncol

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nrow,ncol,this,s) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        this%m(i,j) = s*this%m(i,j)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine subtract(this,a,b)
    class(cdmatrix), intent(inout) :: this
    type(cdmatrix), intent(in) :: a,b
    integer :: i,j,nrow,ncol
    
    if (.not.this%is_allocated) then
      return
    end if
    
    nrow = this%nrow
    ncol = this%ncol

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nrow,ncol,this,a,b) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        this%m(i,j) = a%m(i,j) - b%m(i,j)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine sum(this,a,b)
    class(cdmatrix), intent(inout) :: this
    type(cdmatrix), intent(in) :: a,b
    integer :: i,j,nrow,ncol
    
    if (.not.this%is_allocated) then
      return
    end if
    
    nrow = this%nrow
    ncol = this%ncol

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nrow,ncol,this,a,b) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        this%m(i,j) = a%m(i,j) + b%m(i,j)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine set_to_zero(this)
    class(cdmatrix), intent(inout) :: this
    this%m=(0.d0,0.d0)
  end subroutine
!========================================================================================!

!========================================================================================!
  function inf_norm(this) result(r)
    class(cdmatrix), intent(in) :: this
    real(double_p), allocatable, dimension(:) :: work
    integer :: m,n,iarow,mp0,ierror
    real(double_p) :: r

    if (this%is_allocated) then
    
      m=this%n
      n=m

      iarow=indxg2p(1,this%nrblk,this%prow,0,this%nprow)
      mp0=NUMROC(m,this%nrblk,this%prow,iarow,this%nprow)
      call allocate_array(work,1,mp0)
      
      call set_num_threads()

      r=pzlange('I',m,n,this%m,1,1,this%desc,work)
      
      call set_sigle_thread()
      
    end if
    
    call mpi_bcast(r,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)

  end function
!========================================================================================!

!========================================================================================!
  function l2_norm(this) result(r)
    class(cdmatrix), intent(in) :: this
    integer :: n,nrow,ncol,i,j,ierror
    real(double_p) :: r

    if (this%is_allocated) then
      
      n=this%n
      nrow=this%nrow
      ncol=this%ncol
      
      r=0.d0
      
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP SHARED(nrow,ncol,this) &
      !$OMP PRIVATE(i,j) &
      !$OMP REDUCTION(+:r)
      do j=1,ncol
        do i=1,nrow
          r = r + this%m(i,j)*conjg(this%m(i,j))
        end do
      end do
      !$OMP END PARALLEL DO
      
      r = r/(n*n)
      
    end if
    
    call mpi_bcast(r,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
    r = sqrt(r)

  end function
!========================================================================================!

!========================================================================================!
  subroutine copy_values(lhs,rhs)
    class(cdmatrix), intent(inout) :: lhs
    type(cdmatrix), intent(in) :: rhs
    integer :: i,j,nrow,ncol
    
    if ((lhs%is_allocated).AND.(rhs%is_allocated)) then
    
      nrow = lhs%nrow
      ncol = lhs%ncol

      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP SHARED(nrow,ncol,lhs,rhs) &
      !$OMP PRIVATE(i,j)
      do j=1,ncol
        do i=1,nrow
          lhs%m(i,j) = rhs%m(i,j)
        end do
      end do
      !$OMP END PARALLEL DO

    end if

  end subroutine   
!========================================================================================!

!========================================================================================!
  subroutine reset_uppert(gidx,m)
    integer, allocatable, dimension(:,:,:), intent(in) :: gidx
    complex(double_p), allocatable, dimension(:,:), intent(inout) :: m
    integer :: i,j,ie,je,gi,gj
    
    ie = size(m,1)
    je = size(m,2)
    
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(ie,je,gidx,m) &
    !$OMP PRIVATE(i,j,gi,gj)
    do j=1,je
      do i=1,ie
        
        gi = gidx(1,i,j)
        gj = gidx(2,i,j)
        
        if (gj>gi) then
          m(i,j) = (0.d0,0.d0)
        else if (gj==gi) then
          m(i,j) = 0.5d0*m(i,j)
        end if

      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine write_to_disk(this,time_level,mat_dir)
    class(cdmatrix), intent(in) :: this
    integer, intent(in) :: time_level
    character(len=*), intent(in) :: mat_dir
    type(cdmatrix) :: q
    character(len=20) :: time_dir
    character(len=:), allocatable :: path
    integer, allocatable, dimension(:) :: buff_size
    integer, dimension(2) :: grid_info,idx
    integer :: status(MPI_STATUS_SIZE)
    integer(mpi_offset_kind), allocatable, dimension(:) :: disp
    integer :: f_handle,err,local_buff_size,pcol,npcol,i,&
               mpi_value_type,data_size
               
    if (this%dtype.ne.BLOCK_CYCLIC) then
      call abort_run('Input matrix not BLOCK_CYCLIC in write_to_disk')
    end if
    
    mpi_value_type = MPI_DOUBLE_COMPLEX
    data_size = complex_dp_size
    
    call q%ctor(this%mpic,BLOCK_COLUMN,'q',this%n,1,this%nprocs)
    
    call this%cyclic_to_column(q)
    
    if (q%is_allocated) then
    
      write(time_dir,int_format_) time_level
      path=trim(adjustl(time_dir))//'/'//mat_dir//'/'//this%fname
    
      pcol=q%pcol
      npcol=q%npcol
    
      grid_info(1)=this%n
      grid_info(2)=npcol

      idx(1)=get_gidx(1,q%ncblk,pcol,npcol)
      idx(2)=get_gidx(q%ncol,q%ncblk,pcol,npcol)
    
      !gather local buff size
      call allocate_array(buff_size,0,npcol-1)
      local_buff_size=size(q%m)
      call MPI_ALLGATHER(local_buff_size, 1, MPI_INTEGER, buff_size, 1,&
                         MPI_INTEGER, q%comm, err)
    
      !compute displacement
      allocate(disp(0:npcol-1))
      disp(0)=0
      do i=1,npcol-1
        if (i==1) then
          disp(i)=data_size*buff_size(i-1)+&
                  integer_size*(size(idx)+size(grid_info))
        else
          disp(i)=disp(i-1)+data_size*buff_size(i-1)+integer_size*size(idx)
        end if
      end do

      call mpi_file_open(q%comm,path,MPI_MODE_CREATE+MPI_MODE_WRONLY,&
                         MPI_INFO_NULL,f_handle,err)
                       
      if (err.ne.MPI_SUCCESS) then
        call abort_run('Write file '//this%fname//': mpi file open failed',err)
      end if

      call mpi_file_seek(f_handle,disp(pcol),MPI_SEEK_SET,err)
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

      call mpi_file_write(f_handle,q%m,size(q%m),mpi_value_type,status,err)
      if (err.ne.MPI_SUCCESS) then
        call abort_run('Write file '//this%fname//': mpi file write field failed',err)
      end if

      call mpi_file_close(f_handle,err)
      if (err.ne.MPI_SUCCESS) then
        call abort_run('Write file '//this%fname//': mpi file close failed',err)
      end if
  
      call MPI_BARRIER(q%comm,err)
      
    end if
    
    call q%delete()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine read_from_disk(this,time_level,mat_dir)
    class(cdmatrix), intent(inout) :: this
    integer, intent(in) :: time_level
    character(len=*), intent(in) :: mat_dir
    type(cdmatrix) :: q
    character(len=20) :: time_dir
    character(len=:), allocatable :: path
    integer, allocatable, dimension(:) :: buff_size
    integer, dimension(2) :: grid_info,idx
    integer :: status(MPI_STATUS_SIZE)
    integer(mpi_offset_kind), allocatable, dimension(:) :: disp
    integer :: f_handle,err,local_buff_size,pcol,npcol,i,&
               mpi_value_type,data_size
               
    if (this%dtype.ne.BLOCK_CYCLIC) then
      call abort_run('Input matrix not BLOCK_CYCLIC in read_to_disk')
    end if
    
    mpi_value_type = MPI_DOUBLE_COMPLEX
    data_size = complex_dp_size
    
    call q%ctor(this%mpic,BLOCK_COLUMN,'q',this%n,1,this%nprocs)

    if (q%is_allocated) then
    
      write(time_dir,int_format_) time_level
      path=trim(adjustl(time_dir))//'/'//mat_dir//'/'//this%fname
    
      pcol=q%pcol
      npcol=q%npcol
    
      !gather local buff size
      call allocate_array(buff_size,0,npcol-1)
      local_buff_size=size(q%m)
      call MPI_ALLGATHER(local_buff_size, 1, MPI_INTEGER, buff_size, 1,&
                         MPI_INTEGER, q%comm, err)
    
      !compute displacement
      allocate(disp(0:npcol-1))
      disp(0)=0
      do i=1,npcol-1
        if (i==1) then
          disp(i)=data_size*buff_size(i-1)+&
                  integer_size*(size(idx)+size(grid_info))
        else
          disp(i)=disp(i-1)+data_size*buff_size(i-1)+integer_size*size(idx)
        end if
      end do

      call mpi_file_open(q%comm,path,MPI_MODE_RDONLY,MPI_INFO_NULL,f_handle,err)

      if (err.ne.MPI_SUCCESS) then
        call abort_run('Read file '//this%fname//': mpi file open failed',err)
      end if

      call mpi_file_seek(f_handle,disp(pcol),MPI_SEEK_SET,err)
      if (err.ne.MPI_SUCCESS) then
        call abort_run('Read file '//this%fname//': mpi file seek failed',err)
      end if
    
      if (IS_MASTER) then
        call mpi_file_read(f_handle,grid_info,size(grid_info),MPI_INTEGER,status,err)
        if (err.ne.MPI_SUCCESS) then
          call abort_run('Read file '//this%fname//&
                         ': mpi file read info_grid failed',err)
        end if
        call check_grid_info(grid_info,this%n,q%npcol)
      end if

      call mpi_file_read(f_handle,idx,size(idx),MPI_INTEGER,status,err)
      if (err.ne.MPI_SUCCESS) then
        call abort_run('Read file '//this%fname//': mpi file read idx failed',err)
      end if

      call mpi_file_read(f_handle,q%m,size(q%m),mpi_value_type,status,err)
      if (err.ne.MPI_SUCCESS) then
        call abort_run('Write file '//this%fname//': mpi file read field failed',err)
      end if

      call mpi_file_close(f_handle,err)
      if (err.ne.MPI_SUCCESS) then
        call abort_run('Read file '//this%fname//': mpi file close failed',err)
      end if
  
      call MPI_BARRIER(q%comm,err)
    
    end if
    
    call q%column_to_cyclic(this)
    
    call q%delete()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine check_grid_info(grid_info,n,npcol)
    integer, dimension(2), intent(in) :: grid_info
    integer, intent(in) :: n,npcol
    
    if (grid_info(1).ne.n) then
      call abort_run('Read matrix size not equal to actual matrix size')
    end if

   if (grid_info(2).ne.npcol) then
      call abort_run('Read npcol not equal to actual npcol')
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  function is_stored(this,time_level,mat_dir)
    class(cdmatrix), intent(inout) :: this
    integer, intent(in) :: time_level
    character(len=*), intent(in) :: mat_dir
    character(len=20) :: time_dir
    character(len=:), allocatable :: path
    logical :: is_stored
    integer :: ierror
    
    if (IS_MASTER) then
      
      write(time_dir,int_format_) time_level
      path=trim(adjustl(time_dir))//'/'//mat_dir//'/'//this%fname
      
      inquire(file=path,exist=is_stored)
      
    end if
    
    call mpi_bcast(is_stored,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
      
  end function
!========================================================================================!

!========================================================================================!
  subroutine allocate_ptrm(this) 
    class(cdmatrix), intent(inout) :: this
    integer :: err

    if (.not.associated(this%ptrm)) then
      
      allocate(this%ptrm,STAT=err)
      if (err.ne.0) then
        call abort_run('Allocation of type pointer cdmatrix failed ') 
      end if

      call allocate_array(this%ptrm%m,1,this%nrow,1,this%ncol)
      call allocate_array(this%ptrm%gidx,1,2,1,this%nrow,1,this%ncol)
    
      if (this%is_allocated) then
        call copy_array(this%m,this%ptrm%m)
        call copy_array(this%gidx,this%ptrm%gidx)
      end if
    
      call this%ptrm%equal_dmatrix(this)

    else
      call abort_run('Attempt to allocate an associated type pointer cdmatrix ')
    end if
	
  end subroutine
!========================================================================================!

end module cdmatrix_mod