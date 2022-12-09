module qg_operator_mod

  use dsph_procs_mod
  
  implicit none
  
  type, private :: tdiags
    private

    !diagonal element Helmholtz
    complex(double_p), allocatable, dimension(:) :: d

    !off diagonals and rhs
    complex(double_p), allocatable, dimension(:) :: e,rhs

    !additional Lapack storage needed for Helmholtz operator
    complex(double_p), allocatable, dimension(:) :: u,du2
    integer, allocatable, dimension(:) :: ipiv
    
  end type

  type, extends(tdiag_operator), public :: qg_operator
    private
    
    type(tdiags), allocatable, dimension(:) :: tds
    
    !coriolis matrix
    type(cdmatrix), public :: f
    
    contains
    
    procedure, public :: delete
    procedure, public :: ctor
    
    procedure, public :: apply_inv
    procedure :: init_blks
    
  end type
  
  private :: delete,&
             ctor,&
             apply_inv,&
             init_blks
             
  private :: store_factorization,&
             solve_inv_hel,&
             build_T20,&
             compute_sin_sq

contains

!========================================================================================!
  subroutine delete(this) 
    class(qg_operator), intent(inout) :: this
    
    deallocate(this%tds)
    call this%f%delete(delete_mpi=.FALSE.)  
    call this%delete_tdiag_operator()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ctor(this,mpic,lap,w) 
    class(qg_operator), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    type(laplacian), intent(inout) :: lap
    type(cdmatrix), intent(in) :: w
    
    call this%init_tdiag_operator(mpic,w)
    
    call this%init_blks()
    
    call init_coriolis_matrix(lap,w,this%f)

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine apply_inv(this,w)
    class(qg_operator), intent(inout) :: this
    type(cdmatrix), intent(inout) :: w
    complex(double_p), allocatable, dimension(:) :: rhs
    integer :: i,n,m,s,n_diags
    
    n=this%n
    n_diags=size(this%d_midx)
    
    !subtract coriolis matrix
    call w%subtract(w,this%f)

    call w%cyclic_to_column(this%q)    
    call this%diags_to_buffer()
    
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n_diags,n,this) &
    !$OMP PRIVATE(i,m,s)
    do i=1,n_diags
    
      m=this%d_midx(i)
      s=n-m
      this%tds(i)%rhs=this%d_buff(i,1:s)
      
      call solve_inv_hel(this%tds(i)%d,&
                         this%tds(i)%e,&
                         this%tds(i)%u,&
                         this%tds(i)%du2,&
                         this%tds(i)%ipiv,&
                         this%tds(i)%rhs)
      this%d_buff(i,1:s)=this%tds(i)%rhs
    
    end do
    !$OMP END PARALLEL DO
    
    call this%buffer_to_diags()
    call this%q%column_to_cyclic(w)
    call w%make_skewh() 

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_blks(this)
    class(qg_operator), intent(inout) :: this
    real(double_p), allocatable, dimension(:) :: d_aux,e_aux
    real(double_p) :: k_hel
    integer :: i,j,n,m,err
    complex(double_p), allocatable, dimension(:) :: t
    type(par_file) :: pfile
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(k_hel,'k_hel')
    
    if (.not.this%q%is_allocated) then
      return
    end if
    
    n=this%n
    
    call build_T20(n,t)
    call compute_sin_sq(n,t)
    
    allocate(this%tds(1:size(this%d_midx)),stat=err)
    if (err /= 0) then
      call abort_run('Allocation of qg_operator block systems failed')
    end if    
    
    do i=1,size(this%d_midx)
    
      m=this%d_midx(i)
      
      call compute_lap_blk(n,m,d_aux,e_aux)
      call reallocate_array(this%tds(i)%e,1,n-m)
      call reallocate_array(this%tds(i)%rhs,1,n-m)
      do j=1,n-m
        this%tds(i)%e(j)%re=e_aux(j)
        this%tds(i)%e(j)%im=0.d0
      end do
      
      call reallocate_array(this%tds(i)%d,1,n-m)
      call reallocate_array(this%tds(i)%u,1,n-m)
      do j=1,n-m
        this%tds(i)%d(j)%re = d_aux(j) - k_hel*0.5d0*(t(j)%re+t(j+m)%re)
        this%tds(i)%d(j)%im = 0.d0
      end do
      this%tds(i)%u=this%tds(i)%e
      call store_factorization(this%tds(i)%d,&
                               this%tds(i)%e,&
                               this%tds(i)%u,&
                               this%tds(i)%du2,&
                               this%tds(i)%ipiv)
      
    end do

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine store_factorization(d,e,u,du2,ipiv)
    complex(double_p), allocatable, dimension(:), intent(inout) :: d,e,u
    complex(double_p), allocatable, dimension(:), intent(inout) :: du2
    integer, allocatable, dimension(:), intent(inout) :: ipiv
    integer :: n,info
    
    n=size(d)

    call reallocate_array(du2,1,n-2)
    call reallocate_array(ipiv,1,n)
    
    call zgttrf(n,e,d,u,du2,ipiv,info)

    if (info.ne.0) then
      call abort_run('zgttrf in lap_blk failed',info)
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine solve_inv_hel(d,e,u,du2,ipiv,rhs)
    complex(double_p), allocatable, dimension(:), intent(in) :: d,e,u,du2
    integer, allocatable, dimension(:), intent(in) :: ipiv
    complex(double_p), allocatable, dimension(:), intent(inout) :: rhs
    integer :: n,ldb,nrhs,info
    
    n=size(d)
    ldb=n
    nrhs=1
    
    call zgttrs('N',n,nrhs,e,d,u,du2,ipiv,rhs,ldb,info)

    if (info.ne.0) then
      call abort_run('zgttrs in lap_blk failed',info)
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine compute_sin_sq(n,t)
    integer, intent(in) :: n
    complex(double_p), allocatable, dimension(:), intent(inout) :: t
    integer :: i
    complex(double_p) :: im
    real(double_p) :: x,y
    
    im=(0.d0,1.d0)
    
    do i=1,n
      x=(4.d0/3.d0)*sqrt(pi/n)
      y=-(4.d0/3.d0)*sqrt(pi/5.d0)
      t(i) = im*(x+y*t(i))
    end do
    
    t=-im*t*sqrt(n/(4.d0*pi))

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine build_T20(n,t)
    integer, intent(in) :: n
    complex(double_p), allocatable, dimension(:), intent(out) :: t
    integer :: i
    real(double_p) :: s,c,j
    
    call allocate_array(t,1,n)
    
    s = real(n-1)/real(2.d0)
    
    c = sqrt((2.d0*s+3.d0)*(2.d0*s+2.d0)*(2.d0*s+1.d0)*2.d0*s*(2.d0*s-1.d0))
    
    do i=1,n
      j = i-s-1.d0
      t(i)%re = sqrt(5.d0)*(2.d0*(3*j*j-s*(s+1.d0)))/c
      t(i)%im = 0.d0
    end do

  end subroutine
!========================================================================================!

end module qg_operator_mod
