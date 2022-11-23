module integrators_mod

  use iso_amp_mod
  use iso_gmp_mod
  use heun_mod
  
  implicit none

  integer, parameter :: ISO_AMP_INTEG = 0
  integer, parameter :: ISO_GMP_INTEG = 1
  integer, parameter :: HEUN_INTEG = 2

contains

!========================================================================================!
  subroutine allocate_integrator(integ,w,dt)
    class(integrator), allocatable, intent(inout) :: integ
    type(cdmatrix), intent(in) :: w
    real(double_p), intent(in) :: dt
    integer :: integ_type
    type(par_file) :: pfile
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(integ_type,'integrator')
    
    select case(integ_type)
      case(ISO_AMP_INTEG)
        allocate(iso_amp::integ)
      case(ISO_GMP_INTEG)
        allocate(iso_gmp::integ)
      case(HEUN_INTEG)
        allocate(heun::integ)
      case default
        call abort_run('Wrong integrator found')
    end select

    select type(integ)
      type is(iso_amp)
        call integ%ctor(w,dt)
      type is(iso_gmp)
        call integ%ctor(w,dt)
      type is(heun)
        call integ%ctor(w,dt)
    end select

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine deallocate_integrator(integ)
    class(integrator), intent(inout) :: integ
    
    select type(integ)
      type is(iso_amp)
        call integ%delete()
      type is(iso_gmp)
        call integ%delete()
      type is(heun)
        call integ%delete()
    end select
    
  end subroutine
!========================================================================================!

end module integrators_mod