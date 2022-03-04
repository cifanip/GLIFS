module prun_mod

  use mpi
  use omp_lib

  implicit none
  
  logical, protected :: IS_MASTER, IS_PAR 
  integer, protected :: N_THREADS
  
  public :: init_run
  public :: abort_run
  
contains

!========================================================================================!
  subroutine init_run() 
    integer :: rank, nproc, ierror
    
    call BLACS_PINFO(rank,nproc)
        
    !set isMaster
    if (rank == 0) then
      IS_MASTER = .TRUE.
    else
      IS_MASTER = .FALSE.
    end if
        
    !set isParallel
    if (nproc == 1) then
      IS_PAR = .FALSE.
    else
      IS_PAR = .TRUE.
    end if
        
    !set number of threads
    !$OMP PARALLEL DEFAULT(none) &
    !$OMP SHARED(N_THREADS)
     N_THREADS = omp_get_num_threads() 
    !$OMP END PARALLEL
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine abort_run(msg,opt) 
    character(len=*), intent(in) :: msg
    integer, optional :: opt
    integer :: error_code, ierror

    if (present(opt)) then
      error_code = opt
    else
      error_code = -1
    end if
    
    !if (IS_MASTER) then
      write(*,*) '/***************************************************************/'
      write(*,*) 'mpiABORT CALLED with ERROR CODE: ', error_code
      write(*,*) 'ERROR MESSAGE: ', msg
      write(*,*) '/***************************************************************/'
    !end if
    
    call MPI_ABORT(MPI_COMM_WORLD,error_code,ierror)
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine end_run()
    call BLACS_EXIT(0)
  end subroutine
!========================================================================================!
  
end module prun_mod