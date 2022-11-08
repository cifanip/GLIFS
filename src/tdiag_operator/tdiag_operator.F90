module tdiag_operator_mod

  use cdmatrix_mod
  
  implicit none

  !auxiliary type (diagonal points info)
  type :: dinfo
    integer, allocatable, dimension(:,:) :: idx
    integer :: rank
  end type

  type, abstract, public :: tdiag_operator
    private
    
    integer, public :: n
    type(mpi_control), pointer, public :: mpic => NULL()
    
    type(dinfo), allocatable, dimension(:), public :: dw
    
    !MPI diagonals communication
    integer, allocatable, dimension(:) :: dmat_type,dbuff_type
    integer, allocatable, dimension(:), public :: d_midx
    logical, allocatable, dimension(:,:) :: sr_diag
    complex(double_p), allocatable, dimension(:,:), public :: d_buff
    
    type(cdmatrix), public :: q
    
    contains

    procedure, public :: delete_tdiag_operator
    procedure, public :: init_tdiag_operator
    
    procedure(init_blks_interface_tdiag_operator), public, deferred :: init_blks
    procedure(apply_inv_interface_tdiag_operator), public, deferred :: apply_inv
    
    procedure :: init_dinfo
    procedure :: init_diags
    procedure :: diags_to_buffer
    procedure :: buffer_to_diags

  end type

  interface 
    subroutine init_blks_interface_tdiag_operator(this)
      import
      class(tdiag_operator), intent(inout) :: this
    end subroutine init_blks_interface_tdiag_operator
  end interface

  interface 
    subroutine apply_inv_interface_tdiag_operator(this,w)
      import
      class(tdiag_operator), intent(inout) :: this
      type(cdmatrix), intent(inout) :: w
    end subroutine apply_inv_interface_tdiag_operator
  end interface
  
  private :: delete_tdiag_operator,&
             init_tdiag_operator,&
             init_dinfo,&
             init_diags,&
             diags_to_buffer,&
             buffer_to_diags,&
             init_blks_interface_tdiag_operator,&
             apply_inv_interface_tdiag_operator

contains

!========================================================================================!
  subroutine delete_tdiag_operator(this) 
    class(tdiag_operator), intent(inout) :: this
    integer :: i,ierror

    this%mpic => NULL()
    
    deallocate(this%dw)
    
    if (allocated(this%dmat_type)) then
      do i=lbound(this%dmat_type,1),ubound(this%dmat_type,1)
        if (this%dmat_type(i).ne.MPI_DATATYPE_NULL) then
          call mpi_type_free(this%dmat_type(i),ierror)
        end if
      end do
    end if

    if (allocated(this%dbuff_type)) then
      do i=lbound(this%dbuff_type,1),ubound(this%dbuff_type,1)
        if (this%dbuff_type(i).ne.MPI_DATATYPE_NULL) then
          call mpi_type_free(this%dbuff_type(i),ierror)
        end if
      end do
    end if
    
    call deallocate_array(this%d_midx)
    call deallocate_array(this%sr_diag)
    call deallocate_array(this%d_buff)

    call this%q%delete()
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_tdiag_operator(this,mpic,w) 
    class(tdiag_operator), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    type(cdmatrix), intent(in) :: w

    this%mpic => mpic
    
    this%n = w%n
    
    call this%q%ctor(this%mpic,BLOCK_COLUMN,'q',w%n,1,w%nprocs)
    
    call this%init_dinfo(this%q)

    call this%init_diags()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_diags(this)
    class(tdiag_operator), intent(inout) :: this
    integer :: n,np,ncol,n_diags,s,i,j,p,l,nblocks,gcs,count,rank,ierror,res,displ_0
    integer, allocatable, dimension(:) :: block_lenghts,displ,gncol
    
    if (.not.this%q%is_allocated) then
      return
    end if
    
    n=this%n
    np=this%q%nprocs
    rank=this%q%rank
    ncol=this%q%ncol
    call allocate_array(gncol,0,np-1)
    call mpi_allgather(ncol,1,MPI_INTEGER,gncol,1,MPI_INTEGER,this%q%comm,ierror)
    
    !init logical send/recv
    call allocate_array(this%sr_diag,1,2,0,np-1)
    gcs=get_gidx(1,this%q%ncblk,this%q%pcol,this%q%npcol)
    this%sr_diag=.FALSE.
    !send
    do i=0,np-1
      
      l=gcs+i
      if (l<=n) then
        this%sr_diag(1,i)=.TRUE.
      end if
      
    end do

    !recv
    do i=0,np-1
      
      gcs=get_gidx(1,this%q%ncblk,i,this%q%npcol)
    
      do j=0,np-1
        
        if (gcs+j<=n) then
        
          if (j==rank) then
            this%sr_diag(2,i)=.TRUE.
            exit
          end if
          
        end if
      
      end do
    
    end do
    
    !init mpi derived types
    !dmat_type
    gcs=get_gidx(1,this%q%ncblk,this%q%pcol,this%q%npcol)
    call allocate_array(this%dmat_type,0,np-1)
    this%dmat_type=MPI_DATATYPE_NULL
    do p=0,np-1
    
      if (this%sr_diag(1,p)) then
    
        l=gcs+p
      
        !set block size
        count=0
        do
          nblocks = min(ncol,n-l+1)
          if (nblocks>=1) then
            count=count+nblocks
            l=l+np
          else
            exit
          end if
        end do
      
        !set block lengths and displ
        call reallocate_array(block_lenghts,1,count)
        call reallocate_array(displ,1,count)
        block_lenghts=1
        displ(1)=0
        
        count=1
        do i=1,ncol
        
          l=gcs+p+(i-1)
          
          do
            l=l+np
            if (l<=n) then
              if (count<size(displ)) then
                count=count+1
                displ(count)=displ(count-1)+np
              end if
            else
              if (count<size(displ)) then
                count=count+1
                res=n-(l-np)
                displ(count)=displ(count-1)+res+gcs+p+i
              end if
              exit
            end if
          end do
          
        end do

        call mpi_type_indexed(size(block_lenghts),block_lenghts,displ,&
                              MPI_DOUBLE_COMPLEX,this%dmat_type(p),ierror)
        call mpi_type_commit(this%dmat_type(p),ierror)
        
      end if

    end do
    
    !dbuff_type
    call allocate_array(this%dbuff_type,0,np-1)
    this%dbuff_type=MPI_DATATYPE_NULL
    n_diags=1
    count=rank+1
    do
      count=count+np
      if (count<=n) then
        n_diags=n_diags+1
      else
        exit
      end if
    end do
    call allocate_array(this%d_midx,1,n_diags)
    do i=1,n_diags
      s=n-rank-(i-1)*np
      this%d_midx(i)=n-s
    end do
    call allocate_array(this%d_buff,1,n_diags,1,maxval(n-this%d_midx(:)))
    this%d_buff=(-1.d0,-1.d0)
    
    do p=0,np-1
    
      if (this%sr_diag(2,p)) then
        
        gcs=get_gidx(1,this%q%ncblk,p,this%q%npcol)
        
        count=0
        do i=1,gncol(p)
          
          l=gcs+rank+(i-1)
          if (l<=n) then
            count=count+1
          end if

        end do
        
        call reallocate_array(block_lenghts,1,count)
        call reallocate_array(displ,1,count)
        
        do i=1,size(block_lenghts)
          
          count=0
          l=gcs+rank+(i-1)
          
          do
            if (l<=n) then
              count=count+1
            else
              block_lenghts(i)=count
              if (p>0) then
                displ(i)=(i-1)*n_diags + n_diags*sum(gncol(0:p-1))
              else
                displ(i)=(i-1)*n_diags
              end if
              exit
            end if
            l=l+np
          end do

        end do
        
        call mpi_type_indexed(size(block_lenghts),block_lenghts,displ,&
                              MPI_DOUBLE_COMPLEX,this%dbuff_type(p),ierror)
        call mpi_type_commit(this%dbuff_type(p),ierror)
        
      end if
      
    end do

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine diags_to_buffer(this)
    class(tdiag_operator), intent(inout) :: this
    integer :: n,np,tag,gcs,p,ireq,nreq,ierror,is
    integer, allocatable, dimension(:,:) :: status
    integer, allocatable, dimension(:) :: requests

    if (.not.this%q%is_allocated) then
      return
    end if
     
    n=this%n
    np=this%q%nprocs
    tag=0
    
    nreq=0
    do p=0,np-1
      if (this%sr_diag(1,p)) then
        nreq=nreq+1
      end if
      if (this%sr_diag(2,p)) then
        nreq=nreq+1
      end if
    end do
    call allocate_array(status,1,MPI_STATUS_SIZE,1,nreq)
    call allocate_array(requests,1,nreq)
    
    ireq=0
     
    !recv diags from matrix
    do p=0,np-1
        
      if (this%sr_diag(2,p)) then
      
        ireq=ireq+1
        call MPI_IRECV(this%d_buff,1,this%dbuff_type(p),p,tag,&
                       this%q%comm,requests(ireq),ierror)
      end if
      
    end do
    
    !send diags do buffer
    gcs=get_gidx(1,this%q%ncblk,this%q%pcol,this%q%npcol)
    do p=0,np-1
    
      if (this%sr_diag(1,p)) then
        
        is=gcs+p
        
        ireq=ireq+1
        call MPI_ISSEND(this%q%m(is,1),1,this%dmat_type(p),p,tag,&
                        this%q%comm,requests(ireq),ierror)
      end if
    
    end do
    
    call MPI_WAITALL(nreq,requests,status,ierror)

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine buffer_to_diags(this)
    class(tdiag_operator), intent(inout) :: this
    integer :: n,np,tag,gcs,p,ireq,nreq,ierror,is
    integer, allocatable, dimension(:,:) :: status
    integer, allocatable, dimension(:) :: requests

    if (.not.this%q%is_allocated) then
      return
    end if
     
    n=this%n
    np=this%q%nprocs
    tag=0
    
    nreq=0
    do p=0,np-1
      if (this%sr_diag(1,p)) then
        nreq=nreq+1
      end if
      if (this%sr_diag(2,p)) then
        nreq=nreq+1
      end if
    end do
    call allocate_array(status,1,MPI_STATUS_SIZE,1,nreq)
    call allocate_array(requests,1,nreq)
    
    ireq=0
     
    !recv diags from buffer
    gcs=get_gidx(1,this%q%ncblk,this%q%pcol,this%q%npcol)
    do p=0,np-1
        
      if (this%sr_diag(1,p)) then
      
        is=gcs+p
      
        ireq=ireq+1
        call MPI_IRECV(this%q%m(is,1),1,this%dmat_type(p),p,tag,&
                       this%q%comm,requests(ireq),ierror)
      end if
      
    end do
    
    !send diags to matrix
    do p=0,np-1
    
      if (this%sr_diag(2,p)) then
        
        ireq=ireq+1
        call MPI_ISSEND(this%d_buff,1,this%dbuff_type(p),p,tag,&
                        this%q%comm,requests(ireq),ierror)
      end if
    
    end do
    
    call MPI_WAITALL(nreq,requests,status,ierror)

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_dinfo(this,w)
    class(tdiag_operator), intent(inout) :: this
    type(cdmatrix), intent(in) :: w
    integer :: n,npcol,pcol,ncblk,ncol,m,i,j,gc_s,gc_e,gis,gie,d_size,err
    
    if (.not.w%is_allocated) then
      return
    end if
    
    n=w%n
    npcol=w%npcol
    pcol=w%pcol
    ncblk=w%ncblk
    ncol=w%ncol
    
    allocate(this%dw(0:n-1),stat=err)
    if (err /= 0) then
      call abort_run('Allocation of dinfo in tdiag_operator failed ')
    end if 

    !store diagonal info
    gc_s=get_gidx(1,ncblk,pcol,npcol)
    gc_e=get_gidx(ncol,ncblk,pcol,npcol)
    
    do i=1,size(w%m,1)
      
      m = i-1
      gis = gc_s + m
      
      if (gis > n) then
        this%dw(m)%rank = -1
        cycle
      end if
      
      this%dw(m)%rank = w%rank
      
      gie = min(gis+ncol-1,n)
      d_size = gie - gis +1
      call reallocate_array(this%dw(m)%idx,1,4,1,d_size)
      
      do j=1,d_size
        this%dw(m)%idx(1,j)=gis+(j-1)
        this%dw(m)%idx(2,j)=j
        this%dw(m)%idx(3,j)=this%dw(m)%idx(1,j)-m
      end do
      
    end do

  end subroutine
!========================================================================================!

end module tdiag_operator_mod