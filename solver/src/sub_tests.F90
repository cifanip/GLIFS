!========================================================================================!
  subroutine set_blacs_pgrid(nprow,npcol,prow,pcol,ctxt)
    integer, intent(in) :: nprow,npcol
    integer, intent(out) :: prow,pcol,ctxt
    integer :: i,j,ierror,mctxt
    logical :: invalid,ginvalid
    
    call BLACS_GET(0,0,ctxt)
    call BLACS_GRIDINIT(ctxt,'R',nprow,npcol)
    
    do
    
      call mpi_allreduce(ctxt,mctxt,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierror)
    
      if ((ctxt > -1).AND.(ctxt < mctxt)) then
        invalid = .TRUE.
      else
        invalid = .FALSE.
      end if
    
      call mpi_allreduce(invalid,ginvalid,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierror)
    
      if (ginvalid) then
      
        call BLACS_GET(0,0,ctxt)
        call BLACS_GRIDINIT(ctxt,'R',nprow,npcol)
        if ((ctxt > -1).AND.(ctxt > mctxt)) then
          call delete_ctxt(ctxt)
          ctxt = ctxt - 1
        end if
        
      else
      
        call BLACS_GRIDINFO(ctxt,i,j,prow,pcol)
        return
      
      end if
      
    end do

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