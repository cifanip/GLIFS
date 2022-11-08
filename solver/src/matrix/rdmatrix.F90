module rdmatrix_mod
  
  use dmatrix_mod
  
  implicit none

  type, extends(dmatrix), public :: rdmatrix
    private

    !matrix (local)
    real(double_p), allocatable, dimension(:,:), public :: m
    
    contains

    procedure, public :: rdmatrix_ctor_opt_pgrid
    procedure, public :: rdmatrix_ctor_set_pgrid
    generic, public :: ctor => rdmatrix_ctor_opt_pgrid,&
                               rdmatrix_ctor_set_pgrid
    procedure, public :: delete
    procedure, public :: cyclic_to_column
    procedure, public :: column_to_cyclic
    
    procedure :: equal
    generic, public :: assignment(=) => equal

  end type
  
  private :: rdmatrix_ctor_opt_pgrid,&
             rdmatrix_ctor_set_pgrid,&
             delete,&
             cyclic_to_column,&
             column_to_cyclic,&
             equal

contains

!========================================================================================!
  subroutine delete(this,delete_mpi)
    class(rdmatrix), intent(inout) :: this
    logical, intent(in), optional :: delete_mpi
    
    call deallocate_array(this%m)
    call this%delete_dmatrix(delete_mpi)
        
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine equal(lhs,rhs)
    class(rdmatrix), intent(inout) :: lhs
    type(rdmatrix), intent(in) :: rhs
    
    if (lhs%is_allocated) then
      call abort_run('lhs already allocated found in assignment(=) rdmatrix')
    end if
    
    if (rhs%is_allocated) then
      call copy_array(rhs%m,lhs%m)
    end if
    
    call lhs%equal_dmatrix(rhs)
    
  end subroutine   
!========================================================================================!

!========================================================================================!
  subroutine rdmatrix_ctor_opt_pgrid(this,mpic,distr_type,read_pgrid,n)
    class(rdmatrix), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    integer, intent(in) :: distr_type
    logical, intent(in) :: read_pgrid
    integer, intent(in), optional :: n
    
    call this%init_dmatrix_opt_pgrid(mpic,distr_type,'',read_pgrid,n)
    
    if (this%is_allocated) then
      call allocate_array(this%m,1,this%nrow,1,this%ncol)
      this%m=0.d0
    else
      call allocate_array(this%m,0,0,0,0)
    end if
     
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine rdmatrix_ctor_set_pgrid(this,mpic,distr_type,n,nprow,npcol)
    class(rdmatrix), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    integer, intent(in) :: distr_type
    integer, intent(in) :: n,nprow,npcol
    
    call this%init_dmatrix_set_pgrid(mpic,distr_type,'',n,nprow,npcol)
    
    if (this%is_allocated) then
      call allocate_array(this%m,1,this%nrow,1,this%ncol)
      this%m=0.d0
    else
      call allocate_array(this%m,0,0,0,0)
    end if
     
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine cyclic_to_column(this,q)
    class(rdmatrix), intent(in) :: this
    type(rdmatrix), intent(inout) :: q
    integer :: n
    
    if ((.not.this%is_allocated).AND.(.not.q%is_allocated)) then
      return
    end if
    
    if (this%dtype.ne.BLOCK_CYCLIC) then
      call abort_run('Input matrix not BLOCK_CYLCIC in cyclic_to_column')
    end if

    if (q%dtype.ne.BLOCK_COLUMN) then
      call abort_run('Output matrix not BLOCK_COLUMN in cyclic_to_column')
    end if
    
    n=this%n

    if (this%nprocs >= q%nprocs) then
      call pdgemr2d(n,n,this%m,1,1,this%desc,q%m,1,1,q%desc,this%ctxt)
    else
      call pdgemr2d(n,n,this%m,1,1,this%desc,q%m,1,1,q%desc,q%ctxt)
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine column_to_cyclic(this,q)
    class(rdmatrix), intent(in) :: this
    type(rdmatrix), intent(inout) :: q
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
      call pdgemr2d(n,n,this%m,1,1,this%desc,q%m,1,1,q%desc,this%ctxt)
    else
      call pdgemr2d(n,n,this%m,1,1,this%desc,q%m,1,1,q%desc,q%ctxt)
    end if

  end subroutine
!========================================================================================!

end module rdmatrix_mod