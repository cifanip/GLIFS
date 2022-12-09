module dbracket_mod

  use laplacian_mod
  
  implicit none

  type, extends(laplacian), public :: dbracket
    private
    
    type(cdmatrix) :: z1,z2,z3
    
    real(double_p) :: alpha
    
    contains
    
    procedure, public :: delete
    procedure, public :: ctor
    
    procedure, public :: apply_inv
    
  end type
  
  private :: delete,&
             ctor,&
             apply_inv

contains

!========================================================================================!
  subroutine delete(this) 
    class(dbracket), intent(inout) :: this

    call this%laplacian%delete()
    call this%z1%delete(delete_mpi=.FALSE.)
    call this%z2%delete(delete_mpi=.FALSE.)
    call this%z3%delete(delete_mpi=.FALSE.)

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ctor(this,mpic,w) 
    class(dbracket), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    type(cdmatrix), intent(in) :: w
    type(par_file) :: pfile
    
    call this%laplacian%ctor(mpic,w)
    
    this%z1 = w
    this%z2 = w
    this%z3 = w

    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(this%alpha,'alpha_db')

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine apply_inv(this,w)
    class(dbracket), intent(inout) :: this
    type(cdmatrix), intent(inout) :: w
    integer :: i,j,ncol,nrow
    
    this%z1%m = w%m
    
    call this%laplacian%apply_inv(w)
    
    call this%z2%multiply(w,this%z1)
    call this%z3%multiply(this%z1,w)
    
    ncol=w%ncol
    nrow=w%nrow
    
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(ncol,nrow,w,this) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        w%m(i,j) = w%m(i,j) + this%alpha*(this%z2%m(i,j) - this%z3%m(i,j))
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine
!========================================================================================!

end module dbracket_mod