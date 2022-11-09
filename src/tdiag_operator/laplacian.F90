module laplacian_mod

  use tdiag_operator_mod
  use rdmatrix_mod
  
  implicit none
  
  type, private :: tdiags
    
    !diagonal element Laplacian
    real(double_p), allocatable, dimension(:) :: d
    
    !off diagonals and rhs
    complex(double_p), allocatable, dimension(:) :: e,rhs
    
  end type

  type, extends(tdiag_operator), public :: laplacian
    private
    
    !used to apply inverse
    type(tdiags), allocatable, dimension(:) :: tdsi

    !used to apply
    type(tdiags), allocatable, dimension(:) :: tds
    
    !auxiliary eigenvector matrix
    type(rdmatrix), public :: v
    
    contains
    
    procedure, public :: delete
    procedure, public :: ctor
    
    procedure, public :: apply
    procedure, public :: apply_inv
    procedure :: init_blks
    
  end type
  
  private :: delete,&
             ctor,&
             apply,&
             apply_inv,&
             init_blks
             
  !outside of type
  public :: store_factorization_lap,&
            solve_inv_lap,&
            compute_lap_blk

  private :: solve_inv_m0

contains

!========================================================================================!
  subroutine delete(this) 
    class(laplacian), intent(inout) :: this
    integer :: i,ierror
    
    deallocate(this%tdsi)
    deallocate(this%tds)
    if (this%v%is_allocated) then
      call this%v%delete()
    end if
    call this%delete_tdiag_operator()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ctor(this,mpic,w) 
    class(laplacian), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    type(cdmatrix), intent(in) :: w
    
    call this%init_tdiag_operator(mpic,w)
    
    call this%init_blks()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine apply_inv(this,w)
    class(laplacian), intent(inout) :: this
    type(cdmatrix), intent(inout) :: w
    complex(double_p), allocatable, dimension(:) :: rhs
    integer :: i,n,m,s,n_diags
    
    n=this%n
    n_diags=size(this%d_midx)

    call w%cyclic_to_column(this%q)    
    call this%diags_to_buffer()
    
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n_diags,n,this) &
    !$OMP PRIVATE(i,m,s)
    do i=1,n_diags
    
      m=this%d_midx(i)
      s=n-m
      this%tdsi(i)%rhs=this%d_buff(i,1:s)
      if (m==0) then
        call solve_inv_m0(this%tdsi(i)%d,this%tdsi(i)%e,this%tdsi(i)%rhs)
      else
        call solve_inv_lap(this%tdsi(i)%d,this%tdsi(i)%e,this%tdsi(i)%rhs)
      end if
      this%d_buff(i,1:s)=this%tdsi(i)%rhs
    
    end do
    !$OMP END PARALLEL DO
    
    call this%buffer_to_diags()
    call this%q%column_to_cyclic(w)
    call w%make_skewh()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine apply(this,w)
    class(laplacian), intent(inout) :: this
    type(cdmatrix), intent(inout) :: w
    real(double_p), allocatable, dimension(:) :: aux_e
    complex(double_p), allocatable, dimension(:) :: aux_v
    integer :: i,j,n,n_diags,s
    
    n=this%n
    n_diags=size(this%d_midx)

    call w%cyclic_to_column(this%q)
    call this%diags_to_buffer()

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n_diags,n,this) &
    !$OMP PRIVATE(i,j,s,aux_v,aux_e)
    do i=1,n_diags
    
      s=n-this%d_midx(i)
      
      call allocate_array(aux_v,0,s+1)
      aux_v(0)=(0.d0,0.d0)
      aux_v(s+1)=(0.d0,0.d0)
      aux_v(1:s)=this%d_buff(i,1:s)

      call allocate_array(aux_e,0,s)
      aux_e(0)=0.d0
      aux_e(1:s)=dble(this%tds(i)%e)

      do j=1,s
        this%tds(i)%rhs(j)=aux_v(j-1)*aux_e(j-1) + aux_v(j)*this%tds(i)%d(j)+&
                           aux_v(j+1)*aux_e(j)
      end do
      
      this%d_buff(i,1:s)=this%tds(i)%rhs
      
      call deallocate_array(aux_v)
      call deallocate_array(aux_e)
    
    end do
    !$OMP END PARALLEL DO
    
    call this%buffer_to_diags()
    call this%q%column_to_cyclic(w)
    call w%make_skewh()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_blks(this)
    class(laplacian), intent(inout) :: this
    real(double_p), allocatable, dimension(:) :: e_aux
    integer :: n,m,i,j,err
    
    if (.not.this%q%is_allocated) then
      return
    end if
    
    n=this%n
    
    allocate(this%tdsi(1:size(this%d_midx)),stat=err)
    if (err /= 0) then
      call abort_run('Allocation of laplacian block systems failed')
    end if
    
    allocate(this%tds(1:size(this%d_midx)),stat=err)
    if (err /= 0) then
      call abort_run('Allocation of laplacian block systems failed')
    end if   
    
    do i=1,size(this%d_midx)
    
      m=this%d_midx(i)
      
      call compute_lap_blk(n,m,this%tdsi(i)%d,e_aux)
      call reallocate_array(this%tdsi(i)%e,1,n-m)
      call reallocate_array(this%tdsi(i)%rhs,1,n-m)
      do j=1,n-m
        this%tdsi(i)%e(j)%re=e_aux(j)
        this%tdsi(i)%e(j)%im=0.d0
      end do
      
      call reallocate_array(this%tds(i)%d,1,n-m)
      call reallocate_array(this%tds(i)%e,1,n-m)
      call reallocate_array(this%tds(i)%rhs,1,n-m)
      this%tds(i)%d = - this%tdsi(i)%d
      this%tds(i)%e = - this%tdsi(i)%e
      
      if (m>0) then
        call store_factorization_lap(this%tdsi(i)%d,this%tdsi(i)%e)
      end if
      
    end do

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine store_factorization_lap(d,e)
    real(double_p), allocatable, dimension(:), intent(inout) :: d
    complex(double_p), allocatable, dimension(:), intent(inout) :: e
    integer :: n,info
    
    n=size(d)
    
    call zpttrf(n,d,e,info)

    if (info.ne.0) then
      call abort_run('zpttrf in lap_blk failed',info)
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine solve_inv_lap(d,e,rhs)
    real(double_p), allocatable, dimension(:), intent(in) :: d
    complex(double_p), allocatable, dimension(:), intent(in) :: e
    complex(double_p), allocatable, dimension(:), intent(inout) :: rhs
    integer :: n,ldb,nrhs,info
    
    n=size(d)
    ldb=n
    nrhs=1
    
    call zpttrs('L',n,nrhs,d,e,rhs,ldb,info)

    if (info.ne.0) then
      call abort_run('zpttrs in lap_blk failed',info)
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine solve_inv_m0(d,e,rhs)
    real(double_p), allocatable, dimension(:), intent(in) :: d
    complex(double_p), allocatable, dimension(:), intent(in) :: e
    complex(double_p), allocatable, dimension(:), intent(inout) :: rhs
    complex(double_p), allocatable, dimension(:) :: aux
    real(double_p) :: r,r0,w
    integer :: i,n
    
    n=size(d)
    call allocate_array(aux,1,n)
    
    rhs=rhs-sum(rhs)/n

    r=1.d0/d(1)
    rhs(1)=r*rhs(1)
    aux(1)=1.d0
    do i=2,n
      r0=1.d0/d(i-1)
      r=1.d0/d(i)
      w=e(i-1)*r/aux(i-1)
      aux(i)=1.d0-w*r0*e(i-1)
      rhs(i)=r*rhs(i)-w*rhs(i-1)
    end do
    rhs(n)=0.d0    
    do i=n-1,1,-1
      r=1.d0/d(i)
      rhs(i)=(rhs(i)-r*e(i)*rhs(i+1))/aux(i)
    end do
    rhs=rhs-sum(rhs)/n

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine compute_lap_blk(n,m,d,e)
    integer, intent(in) :: n,m
    real(double_p), allocatable, dimension(:), intent(out) :: d,e
    integer, allocatable, dimension(:) :: j
    real(double_p) :: s
    integer :: i,p
    
    s=real(n-1)/real(2)
    
    call allocate_array(j,1,n-m)
    call allocate_array(d,1,size(j))
    call allocate_array(e,1,size(j))
    d=0.d0
    e=0.d0
      
    !fill up diagonal
    j=(/(i,i=0,n-m-1)/)
    do p=1,size(j)
      d(p)=2.d0*(s*(2.d0*j(p)+1.d0+m)-j(p)*(j(p)+m))
    end do
      
    !fill up off-diagonal
    j=(/(i,i=1,n-m)/)
    do p=1,size(j)-1
      e(p)=-sqrt(dble((j(p)+m)*(n-j(p)-m)))*sqrt(dble(j(p)*(n-j(p))))
    end do
 
  end subroutine
!========================================================================================!

end module laplacian_mod
