module dmatrix_mod
  
  use matrix_mod
  
  implicit none

  type, abstract, extends(matrix), public :: dmatrix
    private

    !number of proc row/columns
    integer, public :: nprow,npcol

    !proc row/column (grid coordinates)
    integer, public :: prow,pcol
    
    !local number of rows/columns
    integer, public :: nrow,ncol

    !block size
    integer, public :: nrblk,ncblk
    
    !matrix descriptor
    integer, dimension(9), public :: desc
    
    !distribution type
    integer, public :: dtype
    
    contains

    procedure, public :: init_dmatrix_opt_pgrid
    procedure, public :: init_dmatrix_set_pgrid
    procedure, public :: delete_dmatrix
    
    procedure :: set_local_size
    procedure :: set_block_size
    procedure :: default_parameters
    procedure :: init_descriptor

    procedure, public :: equal_dmatrix

  end type
  
  private :: init_dmatrix_opt_pgrid,&
             init_dmatrix_set_pgrid,&
             delete_dmatrix,&
             set_block_size,&
             default_parameters,&
             equal_dmatrix,&
             init_descriptor

contains

!========================================================================================!
  subroutine delete_dmatrix(this,delete_mpi)
    class(dmatrix), intent(inout) :: this
    logical, intent(in), optional :: delete_mpi

    call this%default_parameters
    this%dtype = -1
    this%desc = -1
    call this%delete_matrix(delete_mpi)
        
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine equal_dmatrix(lhs,rhs)
    class(dmatrix), intent(inout) :: lhs
    class(dmatrix), intent(in) :: rhs
    
    lhs%nprow = rhs%nprow
    lhs%npcol = rhs%npcol
    lhs%prow  = rhs%prow
    lhs%pcol  = rhs%pcol
    lhs%nrow  = rhs%nrow
    lhs%ncol  = rhs%ncol
    lhs%nrblk = rhs%nrblk
    lhs%ncblk = rhs%ncblk
    lhs%desc  = rhs%desc
    lhs%dtype = rhs%dtype
    
    call lhs%equal_matrix(rhs)
    
  end subroutine   
!========================================================================================!

!========================================================================================!
  subroutine default_parameters(this)
    class(dmatrix), intent(inout) :: this

    this%prow  = -1
    this%pcol  = -1
    this%nrow  = -1
    this%ncol  = -1

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_dmatrix_opt_pgrid(this,mpic,distr_type,fname,read_pgrid,n)
    class(dmatrix), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    integer, intent(in) :: distr_type
    character(len=*), intent(in) :: fname
    logical, intent(in) :: read_pgrid
    integer, intent(in), optional :: n
    integer :: nprow,npcol
    type(par_file) :: pfile
    
    call this%init_matrix(mpic,fname,n)
    
    this%dtype=distr_type
    
    if (read_pgrid) then
      
      call pfile%ctor('input_parameters','specs')
      call pfile%read_parameter(this%nprow,'nprow')
      call pfile%read_parameter(this%npcol,'npcol')
      
      call check_pgrid(this%n,this%nprow,this%npcol,this%dtype)   
      
    else
    
      select case(this%dtype)
        case(BLOCK_CYCLIC)
          call set_2d_pgrid(this%n,mpic%nprocs,this%nprow,this%npcol)
        case(BLOCK_COLUMN)
          call set_1d_pgrid(this%n,mpic%nprocs,this%npcol)
          this%nprow=1
        case(BLOCK_ROW)
          call set_1d_pgrid(this%n,mpic%nprocs,this%nprow)
          this%npcol=1
        case default
          call abort_run('wrong distribution type found in init_dmatrix_opt_pgrid')
      end select
      
    end if

    call set_blacs_pgrid(this%nprow,this%npcol,this%prow,this%pcol,this%ctxt)

    if (this%ctxt==-1) then
      this%is_allocated = .FALSE.
    else
      this%is_allocated = .TRUE.
    end if

    call init_aux_comm(this%is_allocated,this%comm,this%rank,this%nprocs)

    !exit if context is invalid
    if (.not.this%is_allocated) then
      call this%default_parameters()
      this%desc(2)=-1
      return
    end if
    
    call this%set_block_size()
    call this%set_local_size()

    call this%init_descriptor()
    
    if ((IS_MASTER).AND.(fname=='w')) then
      write(*,'(A,I4,I4)') 'Processor gird:   ', this%nprow,this%npcol
      write(*,'(A,I4,I4)') 'Blocking factors: ', this%nrblk,this%ncblk
    end if
     
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_dmatrix_set_pgrid(this,mpic,distr_type,fname,n,nprow,npcol)
    class(dmatrix), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    integer, intent(in) :: distr_type
    character(len=*), intent(in) :: fname
    integer, intent(in) :: n,nprow,npcol
    
    call this%init_matrix(mpic,fname,n)
    
    this%dtype=distr_type
    
    this%nprow = nprow
    this%npcol = npcol
    
    call check_pgrid(this%n,this%nprow,this%npcol,this%dtype)

    call set_blacs_pgrid(this%nprow,this%npcol,this%prow,this%pcol,this%ctxt)    

    if (this%ctxt==-1) then
      this%is_allocated = .FALSE.
    else
      this%is_allocated = .TRUE.
    end if

    call init_aux_comm(this%is_allocated,this%comm,this%rank,this%nprocs)

    !exit if context is invalid
    if (.not.this%is_allocated) then
      call this%default_parameters()
      this%desc(2)=-1
      return
    end if
    
    call this%set_block_size()
    call this%set_local_size()

    call this%init_descriptor()
     
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine set_block_size(this)
    class(dmatrix), intent(inout) :: this
    integer :: q,r,c
    
    r = ceiling(real(this%n)/real(this%nprow))
    c = ceiling(real(this%n)/real(this%npcol))
    
    select case(this%dtype)
      
      case(BLOCK_CYCLIC)
        
        if (ref_blk_size<=r) then
          this%nrblk = ref_blk_size
        else
          this%nrblk = r
        end if
        
        if (ref_blk_size<=c) then
          this%ncblk = ref_blk_size
        else
          this%ncblk = c
        end if
        
      case(BLOCK_COLUMN)
      
        this%ncblk = c
        this%nrblk = this%n
        
      case(BLOCK_ROW)
        
        this%ncblk = this%n
        this%nrblk = r
        
    end select

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine set_local_size(this)
    class(dmatrix), intent(inout) :: this
    integer :: lld,glld,ierror
    
    this%nrow = NUMROC(this%n,this%nrblk,this%prow,0,this%nprow)
    this%ncol = NUMROC(this%n,this%ncblk,this%pcol,0,this%npcol)
    
    !local leading dimension
    this%lld = this%nrow

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_descriptor(this)
    class(dmatrix), intent(inout) :: this
    integer :: info
    
    call DESCINIT(this%desc,this%n,this%n,this%nrblk,this%ncblk,0,0,&
                  this%ctxt,this%lld,info)
    
    if (info<0) then
      call abort_run('call DESCINIT in dmatrix failed')
    end if

  end subroutine
!========================================================================================!

end module dmatrix_mod