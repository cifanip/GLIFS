module integrator_mod

  use tdiag_operator_mod
  
  implicit none

  type, abstract, public :: integrator
    private
    
    !keep a copy to time-step
    real(double_p), public :: dt
    
    contains

    procedure, public :: delete_integrator
    procedure, public :: init_integrator
    
    procedure(solve_interface_integrator), public, deferred :: solve

  end type

  interface 
    subroutine solve_interface_integrator(this,top,w,psi,dt_opt)
      import
      class(integrator), intent(inout) :: this
      class(tdiag_operator), intent(inout) :: top
      type(cdmatrix), intent(inout) :: w,psi
      real(double_p), intent(in), optional :: dt_opt
    end subroutine solve_interface_integrator
  end interface
  
  private :: delete_integrator,&
             init_integrator,&
             solve_interface_integrator

contains

!========================================================================================!
  subroutine delete_integrator(this) 
    class(integrator), intent(inout) :: this
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_integrator(this,dt) 
    class(integrator), intent(out) :: this
    real(double_p), optional :: dt
    type(par_file) :: pfile

    this%dt = dt

  end subroutine
!========================================================================================!

end module integrator_mod