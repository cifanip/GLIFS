module ns_operator_mod

  use laplacian_mod
  use forcing_mod
  
  implicit none
  
  type, private :: tdiags
    private
    
    !diagonal element Laplacian
    real(double_p), allocatable, dimension(:) :: d
    
    !off diagonals and rhs
    complex(double_p), allocatable, dimension(:) :: e,rhs
    
  end type

  type, extends(tdiag_operator), public :: ns_operator
    private
    
    type(tdiags), allocatable, dimension(:) :: tds
    
    !viscosity
    real(double_p) :: nu

    !damping factor
    real(double_p) :: alpha0

    !damping factor (accounts for 2\nu)
    real(double_p) :: alpha

    !copy of time-step (accounts for strang splitting)
    real(double_p) :: dt
    
    !theta method
    real(double_p) :: theta
    
    !keep a pointer to laplacian and forcing
    type(laplacian), pointer :: lap => NULL()
    type(forcing), pointer   :: forc => NULL()
    
    !matrix used to correct damping at mode (1,0)
    type(cdmatrix) :: d10
    
    contains
    
    procedure, public :: delete
    procedure, public :: ctor
    procedure, public :: apply_inv
    procedure :: apply_forced
    
    procedure :: init_blks
    
  end type
  
  private :: delete,&
             ctor,&
             init_blks,&
             apply_inv,&
             apply_forced

contains

!========================================================================================!
  subroutine delete(this) 
    class(ns_operator), intent(inout) :: this
    
    deallocate(this%tds)
    this%forc => NULL()
    this%lap => NULL()
    call this%delete_tdiag_operator()
    call this%d10%delete(delete_mpi=.FALSE.)

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ctor(this,mpic,w,lap,forc)
    class(ns_operator), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    type(cdmatrix), intent(in) :: w
    type(laplacian), intent(inout), target :: lap
    type(forcing), intent(in), target :: forc
    complex(double_p) :: t
    
    call this%init_tdiag_operator(mpic,w)
    
    call this%init_blks()
    
    this%lap  => lap
    
    this%forc => forc
    
    call init_T10_matrix(lap,w,this%d10)
    call this%d10%ptrm%multiply(w,this%d10)
    t=-this%d10%ptrm%trace()    
    this%d10%m = t*this%dt*this%alpha0*this%d10%m
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine apply_inv(this,w)
    class(ns_operator), intent(inout) :: this
    type(cdmatrix), intent(inout) :: w
    
    call this%apply_forced(this%lap,w,this%forc,this%dt)

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine apply_forced(this,lap,w,forc,dt)
    class(ns_operator), intent(inout) :: this
    type(laplacian), intent(inout) :: lap
    type(cdmatrix), intent(inout) :: w
    type(forcing), intent(inout) :: forc
    real(double_p), intent(in) :: dt
    real(double_p) :: r,s,ts,te
    integer :: i,j,k,n,n_diags,nrow,ncol
    
    ts = MPI_Wtime()
    
    n=this%n
    n_diags=size(this%d_midx)
    
    call forc%update(dt)
    
    call forc%f%ptrm%copy_values(w)
    call lap%apply(w)
    
    nrow=w%nrow
    ncol=w%ncol
    
    r=this%theta*this%dt*this%nu
    s=this%theta*this%dt*this%alpha

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(this,nrow,ncol,r,s,w,forc) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        w%m(i,j)=r*w%m(i,j) + (1.d0-s)*forc%f%ptrm%m(i,j) + forc%f%m(i,j)+&
                 this%d10%m(i,j)
      end do
    end do
    !$OMP END PARALLEL DO
    
    call w%cyclic_to_column(this%q)
    call this%diags_to_buffer()
    
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n_diags,n,this) &
    !$OMP PRIVATE(i,k)
    do i=1,n_diags
    
      k=n-this%d_midx(i)
      this%tds(i)%rhs=this%d_buff(i,1:k)
      call solve_inv_lap(this%tds(i)%d,this%tds(i)%e,this%tds(i)%rhs)
      this%d_buff(i,1:k)=this%tds(i)%rhs
    
    end do
    !$OMP END PARALLEL DO
    
    call this%buffer_to_diags()
    call this%q%column_to_cyclic(w)
    call w%make_skewh()

    te = MPI_Wtime() 

    if (IS_MASTER) then
      write(*,'(A,'//output_format_(2:9)//')') '    APPLY H_TURB CPU TIME: ', te-ts
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_blks(this)
    class(ns_operator), intent(inout) :: this
    real(double_p), allocatable, dimension(:) :: e_aux
    integer :: n,m,i,j,err
    type(par_file) :: pfile
    real(double_p) :: dt,nu,alpha,theta
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(dt,'dt')
    call pfile%read_parameter(nu,'nu')
    call pfile%read_parameter(alpha,'alpha')
    theta = 0.5d0
    !symmetric strang splitting
    dt=0.5d0*dt
    
    this%nu = nu
    this%alpha0 = alpha
    this%alpha = alpha - 2.d0*nu
    this%dt = dt
    this%theta = theta
    
    if (.not.this%q%is_allocated) then
      return
    end if
    
    n=this%n
    
    allocate(this%tds(1:size(this%d_midx)),stat=err)
    if (err /= 0) then
      call abort_run('Allocation of ns_operator block systems failed')
    end if    

    do i=1,size(this%d_midx)
    
      m=this%d_midx(i)

      call compute_lap_blk(n,m,this%tds(i)%d,e_aux)
      this%tds(i)%d=-this%tds(i)%d
      e_aux=-e_aux  
        
      call reallocate_array(this%tds(i)%e,1,n-m)
      call reallocate_array(this%tds(i)%rhs,1,n-m)
      do j=1,n-m
        this%tds(i)%e(j)%re=-(1.d0-theta)*dt*nu*e_aux(j)
        this%tds(i)%e(j)%im=0.d0
      end do
      this%tds(i)%d = -(1.d0-theta)*dt*nu*this%tds(i)%d +1.d0 +&
                       (1.d0-theta)*dt*this%alpha
      call store_factorization_lap(this%tds(i)%d,this%tds(i)%e)
      
    end do

  end subroutine
!========================================================================================!

end module ns_operator_mod
