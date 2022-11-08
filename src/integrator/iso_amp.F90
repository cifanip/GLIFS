module iso_amp_mod

  use integrator_mod
  
  implicit none

  type, extends(integrator), public :: iso_amp

    !work matrices
    type(cdmatrix) :: wt0,wt,psi_wt,psi_wth,&
                      psi_wt_psi,delta_w
    
    !iteration implicit step
    real(double_p) :: tol
    integer :: max_iter
    
    contains
    
    procedure, public :: delete
    procedure, public :: ctor
    procedure, public :: solve
    procedure :: clean_up
    procedure :: expl
    procedure :: fpiter
    
  end type
  
  private :: delete,&
             ctor,&
             solve,&
             clean_up,&
             expl,&
             fpiter
             
  private :: fpiter_add,&
             expl_add

contains

!========================================================================================!
  subroutine delete(this) 
    class(iso_amp), intent(inout) :: this

    call this%wt0%delete(delete_mpi=.FALSE.)
    call this%wt%delete(delete_mpi=.FALSE.)
    call this%psi_wt%delete(delete_mpi=.FALSE.)
    call this%psi_wth%delete(delete_mpi=.FALSE.)
    call this%psi_wt_psi%delete(delete_mpi=.FALSE.)
    call this%delta_w%delete(delete_mpi=.FALSE.)
    
    call this%delete_integrator()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ctor(this,w,dt) 
    class(iso_amp), intent(out) :: this
    type(cdmatrix), intent(in) :: w
    real(double_p), intent(in) :: dt
    type(par_file) :: pfile
    
    call this%init_integrator(dt)

    this%wt0        = w
    this%wt         = w
    this%psi_wt     = w
    this%psi_wth    = w
    this%psi_wt_psi = w
    this%delta_w    = w
    
    !clean-up some memory
    call this%clean_up()
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(this%tol,'tol')
    call pfile%read_parameter(this%max_iter,'max_iter')

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine solve(this,top,w,psi,dt_opt)
    class(iso_amp), intent(inout) :: this
    class(tdiag_operator), intent(inout) :: top
    type(cdmatrix), intent(inout) :: w,psi
    real(double_p), intent(in), optional :: dt_opt
    real(double_p) :: ts,te
    
    ts = MPI_Wtime()
    
    call this%fpiter(top,w,psi,dt_opt)    
    call this%expl(top,w,psi,dt_opt)   
      
    te = MPI_Wtime() 

    if (IS_MASTER) then
      write(*,'(A,'//output_format_(2:9)//')') '    ISO_AMP CPU TIME: ', te-ts
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine expl(this,top,w,psi,dt_opt)
    class(iso_amp), intent(inout) :: this
    class(tdiag_operator), intent(inout) :: top
    type(cdmatrix), intent(inout) :: w,psi
    real(double_p), intent(in), optional :: dt_opt
    real(double_p) :: dt
    
    if (present(dt_opt)) then
      dt = dt_opt
    else
      dt = this%dt
    end if
    
    !scale commutator
    dt=dt*(w%n**(1.5d0))/sqrt(16.d0*pi)
    
    !apply inverse operator: laplacian/helmholtz
    call psi%copy_values(this%wt)
    call top%apply_inv(psi)
    
    !(\Delta^{-1} Wt)Wt
    call this%psi_wt%multiply(psi,this%wt)

    !Wt(\Delta^{-1} Wt)
    call this%psi_wth%hconj(this%psi_wt)
    
    !(\Delta^{-1} Wt)Wt(\Delta^{-1} Wt)
    call this%psi_wt_psi%multiply(this%psi_wt,psi)
    
    !add-up
    call expl_add(this%wt%m,         &
                  psi%m,             &
                  this%psi_wt%m,     &
                  this%psi_wth%m,    &
                  this%psi_wt_psi%m, &
                  dt,                &
                  this%wt%nrow,      &
                  this%wt%ncol,      &
                  w%m)
                            
    call w%make_skewh()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine fpiter(this,top,w,psi,dt_opt)
    class(iso_amp), intent(inout) :: this
    class(tdiag_operator), intent(inout) :: top
    type(cdmatrix), intent(inout) :: w,psi
    real(double_p), intent(in), optional :: dt_opt
    real(double_p) :: dt,err
    integer :: iter
    logical :: converged

    if (present(dt_opt)) then
      dt = dt_opt
    else
      dt = this%dt
    end if

    !scale commutator
    dt=dt*(w%n**(1.5d0))/sqrt(16.d0*pi)
    
    iter = 0
    converged = .FALSE.
    
    call this%wt0%copy_values(w)
    
    do
    
      !apply inverse operator: laplacian/helmholtz
      call psi%copy_values(this%wt0)
      call top%apply_inv(psi)

      !(\Delta^{-1} Wt)Wt
      call this%psi_wt%multiply(psi,this%wt0)
    
      !Wt(\Delta^{-1} Wt)
      call this%psi_wth%hconj(this%psi_wt)
    
      !(\Delta^{-1} Wt)Wt(\Delta^{-1} Wt)
      call this%psi_wt_psi%multiply(this%psi_wt,psi)
      
      !add-up
      call fpiter_add(w%m,               &
                      psi%m,             &
                      this%psi_wt%m,     &
                      this%psi_wth%m,    &
                      this%psi_wt_psi%m, &
                      dt,                &
                      this%wt%nrow,      &
                      this%wt%ncol,      &
                      this%wt%m)
      
      !compute infinity norm
      call this%delta_w%subtract(this%wt,this%wt0)
      err = this%delta_w%inf_norm()
    
      if (err <= this%tol) then
        converged = .TRUE.
        if (IS_MASTER) then
          write(*,'(A,I3,A,'//output_format_(2:9)//')') '    converged: iter ',iter,&
                                                   ';  err ',err
        end if
        return
      end if
      
      call this%wt0%copy_values(this%wt)
      
      iter = iter + 1
    
      if (iter == this%max_iter) then
        call abort_run('Maximum number of iterations reached in isomp_fpiter')
      end if
    
    end do

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine fpiter_add(wn,psi,psi_wt,psi_wth,psi_wt_psi,dt,nrow,ncol,wt)
    complex(double_p), allocatable, dimension(:,:), intent(in) ::  wn,psi,psi_wt,&
                                                                   psi_wth,psi_wt_psi
    real(double_p), intent(in) :: dt
    integer, intent(in) :: nrow,ncol
    complex(double_p), allocatable, dimension(:,:), intent(inout) :: wt
    real(double_p) x,y
    integer :: i,j
    
    x=0.5d0*dt
    y=0.25d0*dt*dt

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(ncol,nrow,wt,wn,x,y,psi_wt,psi_wth,psi_wt_psi) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        wt(i,j) = wn(i,j) + x*(psi_wt(i,j)-psi_wth(i,j)) + y*psi_wt_psi(i,j)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine expl_add(wt,psi,psi_wt,psi_wth,psi_wt_psi,dt,nrow,ncol,w)
    complex(double_p), allocatable, dimension(:,:), intent(in) ::  wt,psi,psi_wt,&
                                                                   psi_wth,psi_wt_psi
    real(double_p), intent(in) :: dt
    integer, intent(in) :: nrow,ncol
    complex(double_p), allocatable, dimension(:,:), intent(inout) :: w
    real(double_p) x,y
    integer :: i,j
    
    x=0.5d0*dt
    y=0.25d0*dt*dt

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(ncol,nrow,w,wt,x,y,psi_wt,psi_wth,psi_wt_psi) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        w(i,j) = wt(i,j) + x*(psi_wt(i,j)-psi_wth(i,j)) - y*psi_wt_psi(i,j)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine clean_up(this)
    class(iso_amp), intent(inout) :: this
    
    call this%wt0%ptrm%delete(delete_mpi=.FALSE.)
    call this%wt%ptrm%delete(delete_mpi=.FALSE.)
    call this%psi_wt%ptrm%delete(delete_mpi=.FALSE.)
    call this%psi_wth%ptrm%delete(delete_mpi=.FALSE.)
    call this%psi_wt_psi%ptrm%delete(delete_mpi=.FALSE.)
    call this%delta_w%ptrm%delete(delete_mpi=.FALSE.)  

  end subroutine
!========================================================================================!

end module iso_amp_mod