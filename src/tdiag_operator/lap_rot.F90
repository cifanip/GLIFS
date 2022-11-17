module lap_rot_mod

  use dsph_procs_mod
  
  implicit none

  type, extends(laplacian), public :: lap_rot
    private
    
    !coriolis matrix
    type(cdmatrix), public :: f
    
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
    class(lap_rot), intent(inout) :: this
    integer :: i,ierror
    
    call this%f%delete(delete_mpi=.FALSE.)  
    call this%laplacian%delete()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ctor(this,mpic,w) 
    class(lap_rot), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    type(cdmatrix), intent(in) :: w
    
    call this%laplacian%ctor(mpic,w) 
    
    call init_coriolis_matrix(this%laplacian,w,this%f)

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine apply_inv(this,w)
    class(lap_rot), intent(inout) :: this
    type(cdmatrix), intent(inout) :: w
    
    !subtract coriolis matrix
    call w%subtract(w,this%f)
    
    call this%laplacian%apply_inv(w)

  end subroutine
!========================================================================================!

end module lap_rot_mod