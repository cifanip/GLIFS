module matrix_mod
  
  use mpi_control_mod
  
  implicit none

  integer, parameter :: BLOCK_CYCLIC = 0
  integer, parameter :: BLOCK_COLUMN = 1
  integer, parameter :: BLOCK_ROW    = 2
  
  integer, parameter :: min_mat_size = 1000
  integer, parameter :: ref_blk_size = 128

  type, abstract, public :: matrix
    private
    
    !global matrix size (assumed square global matrix)
    integer, public :: n
    
    !keep a pointer to mpi_control
    type(mpi_control), pointer, public :: mpic => NULL()
    
    !BLACS context
    integer, public :: ctxt
    
    !local leading dimension
    integer, public :: lld
    
    logical, public :: is_allocated = .FALSE.

    !file name
    character(len=:), allocatable, public :: fname
    
    !matrix MPI communicator
    integer, public :: comm   = MPI_COMM_NULL
    integer, public :: rank   = -1
    integer, public :: nprocs
    
    contains
    
    procedure, public :: init_matrix
    procedure, public :: delete_matrix
    procedure, public :: equal_matrix

  end type
  
  private :: init_matrix,&
             delete_matrix
             
  public :: set_1d_pgrid,&
            set_2d_pgrid,&
            delete_ctxt,&
            init_aux_comm,&
            set_blacs_pgrid,&
            equal_matrix,&
            check_pgrid
             
  interface
    integer function numroc(n,nb,iproc,isrcproc,nprocs)
      integer :: n,nb,iproc,isrcproc,nprocs
    end function numroc
  end interface

  interface
    integer function indxg2l(ig,nblk,pcoord,isrcproc,nprocs)
      integer :: ig,nblk,pcoord,isrcproc,nprocs
    end function indxg2l
    
    integer function indxl2g(il,nblk,pcoord,isrcproc,nprocs)
      integer :: il,nblk,pcoord,isrcproc,nprocs
    end function indxl2g

    integer function indxg2p(ig,nblk,tmp,isrcproc,nprocs)
      integer :: il,nblk,pcoord,tmp,isrcproc,nprocs
    end function indxg2p
  end interface

contains

!========================================================================================!
  subroutine delete_matrix(this,delete_mpi)
    class(matrix), intent(inout) :: this
    logical, intent(in), optional :: delete_mpi
    logical :: delete_mpi_opt
    integer :: ierror
    
    if (present(delete_mpi)) then
      delete_mpi_opt = delete_mpi
    else
      delete_mpi_opt = .TRUE.
    end if
        
    this%mpic => NULL()
    this%n     = -1
    this%lld   = -1
    if ((this%comm.ne.MPI_COMM_NULL).AND.(delete_mpi_opt)) then
      call mpi_comm_free(this%comm,ierror)
    end if
    if ((this%is_allocated).AND.(delete_mpi_opt)) then
      call delete_ctxt(this%ctxt)
    end if
    this%is_allocated = .FALSE.
    this%rank=-1
    this%nprocs=-1
        
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine equal_matrix(lhs,rhs)
    class(matrix), intent(inout) :: lhs
    class(matrix), intent(in) :: rhs
    
    lhs%n            = rhs%n
    lhs%mpic         => rhs%mpic
    lhs%ctxt         = rhs%ctxt
    lhs%lld          = rhs%lld
    lhs%is_allocated = rhs%is_allocated
    lhs%comm         = rhs%comm
    lhs%rank         = rhs%rank
    lhs%nprocs       = rhs%nprocs
    lhs%fname        = rhs%fname
    
  end subroutine   
!========================================================================================!

!========================================================================================!
  subroutine init_matrix(this,mpic,fname,n)
    class(matrix), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    character(len=*), intent(in) :: fname
    integer, intent(in), optional :: n
    type(par_file) :: pfile
    integer, allocatable, dimension(:,:) :: distr
    integer :: ctxt,nl
    
    this%mpic => mpic
    this%fname = fname
    
    !read matrix size
    if (present(n)) then
      this%n=n
    else
      call pfile%ctor('input_parameters','specs')
      call pfile%read_parameter(this%n,'N')
    end if
     
  end subroutine
!========================================================================================!

!========================================================================================!
  function get_lidx(ig,nblk,pcoord,nprocs) result(il)
    integer, intent(in) :: ig,nblk,pcoord,nprocs
    integer :: il
    
    il = indxg2l(ig,nblk,pcoord,0,nprocs)

  end function
!========================================================================================!

!========================================================================================!
  function get_gidx(il,nblk,pcoord,nprocs) result(ig)
    integer, intent(in) :: il,nblk,pcoord,nprocs
    integer :: ig
    
    ig = indxl2g(il,nblk,pcoord,0,nprocs)

  end function
!========================================================================================!

!========================================================================================!
  function get_pcoord(ig,nblk,nprocs) result(pcoord)
    integer, intent(in) :: ig,nblk,nprocs
    integer :: pcoord,tmp
    
    pcoord = indxg2p(ig,nblk,tmp,0,nprocs)

  end function
!========================================================================================!

!========================================================================================!
  subroutine set_1d_pgrid(n,p,ngrid)
    integer, intent(in) :: n,p
    integer, intent(out) :: ngrid
    integer :: nl,pt,nc
    
    !target number of procs
    pt=n/min_mat_size
    
    !more procs needed than available
    if (pt>=p) then
      ngrid=p
    else
      ngrid=max(pt,1)
    end if

    !use only one core
    if (ngrid==1) then
      return
    end if
    
    !all but the last proc gets same amount of points
    ngrid = ceiling(dble(n)/ceiling(dble(n)/dble(ngrid)))

    if (ceiling(dble(n)/dble(ngrid))<1) then
      call abort_run('nroc < 1 found in set_1d_pgrid')
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine set_2d_pgrid(n,p,nprow,npcol)
    integer, intent(in) :: n,p
    integer, intent(out) :: nprow,npcol
    integer :: nl,pt,q,nc
    
    !target number of procs
    pt=n*n/(min_mat_size*min_mat_size)
    
    !more procs needed than available
    if (pt>=p) then
      q=p
    else
      q=max(pt,1)
    end if
    
    !use only one core
    if (q==1) then
      nprow=1
      npcol=1
      return
    end if
    
    !split grid equally
    nprow=floor(sqrt(real(q)))

    !all but the last proc gets same amount of points
    nprow = ceiling(dble(n)/ceiling(dble(n)/dble(nprow)))
    npcol = nprow

    if (ceiling(dble(n)/dble(nprow))<1) then
      call abort_run('nroc < 1 found in set_2d_pgrid')
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine check_pgrid(n,pr,pc,grid_type)
    integer, intent(in) :: n,pr,pc,grid_type
    integer :: qr,qc
    
    select case(grid_type)
      case(BLOCK_COLUMN)
        if (pr > 1) then
          call abort_run('pr > 1 found in BLOCK_COLUMN check_pgrid')
        end if
      case(BLOCK_ROW)
        if (pc > 1) then
          call abort_run('pc > 1 found in BLOCK_ROW check_pgrid')
        end if
    end select
    
    qr = (pr-1)*ceiling(dble(n)/dble(pr))
    if (qr>n) then
      call abort_run('Invalid number of processor rows found in check_pgrid')
    end if
    if (ceiling(dble(n)/dble(pr))<1) then
      call abort_run('nrow < 1 found in set_2d_pgrid')
    end if

    qc = (pc-1)*ceiling(dble(n)/dble(pc))
    if (qc>n) then
      call abort_run('Invalid number of processor columns found in check_pgrid')
    end if
    if (ceiling(dble(n)/dble(pc))<1) then
      call abort_run('ncol < 1 found in set_2d_pgrid')
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_aux_comm(is_allocated,comm,rank,size)
    logical, intent(in) :: is_allocated
    integer, intent(out) :: comm,rank,size
    integer, allocatable, dimension(:) :: g_color,ranks
    integer :: group,new_group,ierror,n,color,i,j,nprocs
    
    call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierror)

    call mpi_comm_group(MPI_COMM_WORLD,group,ierror)
    
    call allocate_array(g_color,0,nprocs-1)
    g_color=0
    
    if (is_allocated) then
      color=1
    else
      color=0
    end if
    
    call mpi_allgather(color,1,MPI_INTEGER,g_color,1,MPI_INTEGER,MPI_COMM_WORLD,ierror)
    
    n=sum(g_color)
    call allocate_array(ranks,1,n)
    
    j=0
    do i=0,nprocs-1
      if (g_color(i)==1) then
        j=j+1
        ranks(j)=i
      end if
    end do
    
    call mpi_group_incl(group, n, ranks, new_group, ierror)
    
    call mpi_comm_create(MPI_COMM_WORLD, new_group, comm, ierror)

    if (ierror.ne.0) then
      call abort_run('mpi_comm_create in init_aux_comm failed')
    end if
    
    if (is_allocated) then
      call mpi_comm_rank(comm,rank,ierror)
      call mpi_comm_size(comm,size,ierror)
    end if
    
    call mpi_bcast(size,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine set_blacs_pgrid(nprow,npcol,prow,pcol,ctxt)
    integer, intent(in) :: nprow,npcol
    integer, intent(out) :: prow,pcol,ctxt
    integer :: i,j
    
    call BLACS_GET(0,0,ctxt)
    call BLACS_GRIDINIT(ctxt,'R',nprow,npcol)
    call BLACS_GRIDINFO(ctxt,i,j,prow,pcol)

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine delete_ctxt(ctxt)
    integer, intent(inout) :: ctxt
    call blacs_gridexit(ctxt)
  end subroutine
!========================================================================================!

end module matrix_mod