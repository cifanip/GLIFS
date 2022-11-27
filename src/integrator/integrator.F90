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
    
    procedure(print_gmap_interface_integrator), public, deferred :: print_gmap

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

  interface 
    subroutine print_gmap_interface_integrator(this,output_dir,fields_dir)
      import
      class(integrator), intent(inout) :: this
      integer, intent(in) :: output_dir
      character(len=*), intent(in) :: fields_dir
    end subroutine print_gmap_interface_integrator
  end interface
  
  private :: delete_integrator,&
             init_integrator,&
             solve_interface_integrator,&
             print_gmap_interface_integrator
             
  !outside of type
   public :: apply_Ads,&
             update_gt_map,&
             compute_commutator

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

!========================================================================================!
  subroutine apply_Ads(g,gd,w_gd,w) 
    type(cdmatrix), intent(in) :: g
    type(cdmatrix), intent(inout) :: gd,w_gd,w
    
    call gd%hconj(g)
    
    call w_gd%multiply(w,gd)

    call w%multiply(g,w_gd)

    call w%make_skewh()
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine update_gt_map(g,aux,gt)
    type(cdmatrix), intent(in) :: g
    type(cdmatrix), intent(inout) :: aux,gt
    
    call aux%copy_values(gt)    
    call gt%multiply(g,aux)

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine compute_commutator(a,b,c)
    type(cdmatrix), intent(in) :: a,b
    type(cdmatrix), intent(inout) :: c
    integer :: i,j,ncol,nrow
    
    call c%multiply(a,b)
    call c%ptrm%multiply(b,a)
    
    ncol=c%ncol
    nrow=c%nrow
    
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(a,b,c,ncol,nrow) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        c%m(i,j) = c%m(i,j) - c%ptrm%m(i,j)
      end do
    end do
    !$OMP END PARALLEL DO 

  end subroutine
!========================================================================================!

end module integrator_mod