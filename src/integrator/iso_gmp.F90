module iso_gmp_mod

  use integrator_mod
  
  implicit none

  type, extends(integrator), public :: iso_gmp
    private

    !Lie group element
    type(cdmatrix) :: gt
    
    !identity matrix
    type(cdmatrix) :: id
    
    !work matrices
    type(cdmatrix) :: g,g0,z1,z2,z3
    
    !iteration implicit step
    real(double_p) :: tol
    integer :: max_iter
    
    contains
    
    procedure, public :: delete
    procedure, public :: ctor
    procedure, public :: solve
    procedure, public :: print_gmap
    procedure :: expl_mp
    procedure :: euler_fwd
    procedure :: clean_up
    
  end type
  
  private :: delete,&
             ctor,&
             solve,&
             print_gmap,&
             expl_mp,&
             euler_fwd,&
             clean_up
             
  private :: compute_gh

contains

!========================================================================================!
  subroutine delete(this) 
    class(iso_gmp), intent(inout) :: this

    call this%gt%delete(delete_mpi=.FALSE.)
    call this%g%delete(delete_mpi=.FALSE.)
    call this%g0%delete(delete_mpi=.FALSE.)
    call this%id%delete(delete_mpi=.FALSE.)
    call this%z1%delete(delete_mpi=.FALSE.)
    call this%z2%delete(delete_mpi=.FALSE.)
    call this%z3%delete(delete_mpi=.FALSE.)
    
    call this%delete_integrator()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ctor(this,w,dt) 
    class(iso_gmp), intent(out) :: this
    type(cdmatrix), intent(in) :: w
    real(double_p), intent(in) :: dt
    type(par_file) :: pfile
    
    call this%init_integrator(dt)

    this%gt = w
    call this%gt%rename('gt')
    
    this%g  = w
    this%id = w
    this%z1 = w
    this%z2 = w
    this%z3 = w
    this%g0 = w
    
    call this%g%set_identity()
    call this%gt%set_identity()
    call this%id%set_identity()
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(this%tol,'tol')
    call pfile%read_parameter(this%max_iter,'max_iter')
    
    call this%clean_up()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine solve(this,top,w,psi,dt_opt)
    class(iso_gmp), intent(inout) :: this
    class(tdiag_operator), intent(inout) :: top
    type(cdmatrix), intent(inout) :: w,psi
    real(double_p), intent(in), optional :: dt_opt
    real(double_p) :: dt,err,ts,te
    integer :: iter,ncol,nrow,i,j
    logical :: converged
    
    ts = MPI_Wtime()

    if (present(dt_opt)) then
      dt = dt_opt
    else
      dt = this%dt
    end if
    
    !initial guess of group element
    !-- call this%euler_fwd(top,w,psi,dt_opt)
    !-- call this%expl_mp(top,w,psi,dt_opt)
    this%g0%m=this%g%m

    !scale commutator
    dt=dt*(w%n**(1.5d0))/sqrt(16.d0*pi)
    
    iter = 0
    converged = .FALSE.
    
    nrow=w%nrow
    ncol=w%ncol

    do
    
      !z1 = gh = 0.5*(I+Ut)
      call compute_gh(this%id,this%g,this%z1)
    
      !z2 = gh dagger = z1 dagger
      call this%z2%hconj(this%z1)
      
      !z3 = Wn*gh' = w*z2
      call this%z3%multiply(w,this%z2)

      !z2 = Wh = gh*Wn*gh' = z1*z3
      call this%z2%multiply(this%z1,this%z3)
    
      !apply inverse operator: laplacian/helmholtz
      call psi%copy_values(this%z2)
      call top%apply_inv(psi)

      !z2 = Psi(Wh)*gh = psi*z1
      call this%z2%multiply(psi,this%z1)
      
      !update group element
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP SHARED(this,dt,ncol,nrow) &
      !$OMP PRIVATE(i,j)
      do j=1,ncol
        do i=1,nrow
          this%g%m(i,j) = this%id%m(i,j) + dt*this%z2%m(i,j)
        end do
      end do
      !$OMP END PARALLEL DO
      
      !compute infinity norm
      call this%z1%subtract(this%g,this%g0)
      err = this%z1%inf_norm()
    
      if (err <= this%tol) then
        converged = .TRUE.
         if (IS_MASTER) then
           write(*,'(A,I3,A,'//output_format_(2:9)//')') '    converged: iter ',iter,&
                                                    ';  err ',err
         end if
         exit
      end if
      
      call this%g0%copy_values(this%g)
      
      iter = iter + 1
    
      if (iter == this%max_iter) then
        call abort_run('Maximum number of iterations reached in iso_gmp%solve()')
      end if
    
    end do
    
    !Ad^* operator
    call apply_Ads(this%g,this%z1,this%z2,w)
    
    !update g_t map
    call update_gt_map(this%g,this%gt)

    te = MPI_Wtime() 

    if (IS_MASTER) then
      write(*,'(A,'//output_format_(2:9)//')') '    ISO_GMP CPU TIME: ', te-ts
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine expl_mp(this,top,w,psi,dt_opt)
    class(iso_gmp), intent(inout) :: this
    class(tdiag_operator), intent(inout) :: top
    type(cdmatrix), intent(in) :: w
    type(cdmatrix), intent(inout) :: psi
    real(double_p), intent(in), optional :: dt_opt
    real(double_p) :: dt,h
    integer :: i,j,ncol,nrow
    
    if (present(dt_opt)) then
      dt = dt_opt
    else
      dt = this%dt
    end if
    
    !scale commutator
    dt=dt*(w%n**(1.5d0))/sqrt(16.d0*pi)
    
    call psi%copy_values(w)
    call top%apply_inv(psi)
    
    h=0.5*dt
    
    ncol=w%ncol
    nrow=w%nrow

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(this,dt,psi,h,ncol,nrow) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        this%g%m(i,j) = this%id%m(i,j) + h*psi%m(i,j)
      end do
    end do
    !$OMP END PARALLEL DO
    
    call this%z3%copy_values(w)
    
    call apply_Ads(this%g,this%z1,this%z2,this%z3)

    call psi%copy_values(this%z3)
    call top%apply_inv(psi)
    
    call this%z1%multiply(psi,this%g)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(this,dt,psi,h,ncol,nrow) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        this%g%m(i,j) = this%id%m(i,j) + dt*this%z1%m(i,j)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine euler_fwd(this,top,w,psi,dt_opt)
    class(iso_gmp), intent(inout) :: this
    class(tdiag_operator), intent(inout) :: top
    type(cdmatrix), intent(in) :: w
    type(cdmatrix), intent(inout) :: psi
    real(double_p), intent(in), optional :: dt_opt
    real(double_p) :: dt
    integer :: i,j,ncol,nrow
    
    if (present(dt_opt)) then
      dt = dt_opt
    else
      dt = this%dt
    end if
    
    !scale commutator
    dt=dt*(w%n**(1.5d0))/sqrt(16.d0*pi)
    
    ncol=w%ncol
    nrow=w%nrow
    
    !apply inverse operator: laplacian/helmholtz
    call psi%copy_values(w)
    call top%apply_inv(psi)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(this,dt,psi,ncol,nrow) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        this%g%m(i,j) = this%id%m(i,j) + dt*psi%m(i,j)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine compute_gh(id,g,gh)
    type(cdmatrix), intent(in) :: id,g
    type(cdmatrix), intent(inout) :: gh
    integer :: i,j,ncol,nrow
    
    ncol=gh%ncol
    nrow=gh%nrow
    
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(gh,id,g,ncol,nrow) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        gh%m(i,j) = 0.5d0*(id%m(i,j)+g%m(i,j))
      end do
    end do
    !$OMP END PARALLEL DO    

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine clean_up(this)
    class(iso_gmp), intent(inout) :: this
    
    call this%gt%ptrm%delete(delete_mpi=.FALSE.)
    call this%g%ptrm%delete(delete_mpi=.FALSE.)
    call this%id%ptrm%delete(delete_mpi=.FALSE.)
    call this%z1%ptrm%delete(delete_mpi=.FALSE.)
    call this%z2%ptrm%delete(delete_mpi=.FALSE.)
    call this%g0%ptrm%delete(delete_mpi=.FALSE.)
    
    !Note: z3 calls make_skewh, thus ptrm has to be allocated

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine print_gmap(this,output_dir,fields_dir)
    class(iso_gmp), intent(inout) :: this
    integer, intent(in) :: output_dir
    character(len=*), intent(in) :: fields_dir
    
    call this%gt%write_to_disk(output_dir,fields_dir)

  end subroutine
!========================================================================================!

end module iso_gmp_mod
