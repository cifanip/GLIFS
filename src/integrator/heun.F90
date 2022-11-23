module heun_mod

  use integrator_mod
  
  implicit none

  type, extends(integrator), public :: heun
    private
    
    !work matrices
    type(cdmatrix) :: z1,z2,z3
    
    contains
    
    procedure, public :: delete
    procedure, public :: ctor
    procedure, public :: solve
    procedure, public :: print_gmap
    procedure :: clean_up
    
  end type
  
  private :: delete,&
             ctor,&
             solve,&
             print_gmap,&
             clean_up

contains

!========================================================================================!
  subroutine delete(this) 
    class(heun), intent(inout) :: this

    call this%z1%delete(delete_mpi=.FALSE.)
    call this%z2%delete(delete_mpi=.FALSE.)
    call this%z3%delete(delete_mpi=.FALSE.)
    
    call this%delete_integrator()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ctor(this,w,dt) 
    class(heun), intent(out) :: this
    type(cdmatrix), intent(in) :: w
    real(double_p), intent(in) :: dt
    type(par_file) :: pfile
    
    call this%init_integrator(dt)

    this%z1 = w
    this%z2 = w
    this%z3 = w
    
    call this%clean_up()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine solve(this,top,w,psi,dt_opt)
    class(heun), intent(inout) :: this
    class(tdiag_operator), intent(inout) :: top
    type(cdmatrix), intent(inout) :: w,psi
    real(double_p), intent(in), optional :: dt_opt
    real(double_p) :: dt,ts,te
    integer :: ncol,nrow,i,j
    logical :: converged
    
    ts = MPI_Wtime()

    if (present(dt_opt)) then
      dt = dt_opt
    else
      dt = this%dt
    end if

    !scale commutator
    dt=dt*(w%n**(1.5d0))/sqrt(16.d0*pi)
    
    nrow=w%nrow
    ncol=w%ncol

    call psi%copy_values(w)
    call top%apply_inv(psi)
    
    call compute_commutator(psi,w,this%z1)
    
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(w,this,dt,ncol,nrow) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        this%z2%m(i,j) = w%m(i,j) + dt*this%z1%m(i,j)
      end do
    end do
    !$OMP END PARALLEL DO

    call psi%copy_values(this%z2)
    call top%apply_inv(psi)
    
    call compute_commutator(psi,this%z2,this%z3)

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(this,w,dt,ncol,nrow) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        w%m(i,j) = w%m(i,j) + 0.5*dt*(this%z1%m(i,j)+this%z3%m(i,j))
      end do
    end do
    !$OMP END PARALLEL DO

    te = MPI_Wtime() 

    if (IS_MASTER) then
      write(*,'(A,'//output_format_(2:9)//')') '    HEUN CPU TIME: ', te-ts
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine clean_up(this)
    class(heun), intent(inout) :: this

    call this%z2%ptrm%delete(delete_mpi=.FALSE.)
    
    !Note: z1,z3 call compute_commutator, thus ptrm has to be allocated

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine print_gmap(this,output_dir,fields_dir)
    class(heun), intent(inout) :: this
    integer, intent(in) :: output_dir
    character(len=*), intent(in) :: fields_dir
    
    !nothing to be done
    return

  end subroutine
!========================================================================================!

end module heun_mod