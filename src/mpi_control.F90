module mpi_control_mod

  use allocate_arrays_mod
  use par_file_mod
  use aux_mod
  
  implicit none

  type, public :: mpi_control
    private
    
    integer, public :: rank
    integer, public :: nprocs
    
    contains

    procedure, public :: delete
    procedure, public :: ctor

  end type
  
  private :: delete,&
             ctor


contains

!========================================================================================!
  subroutine delete(this) 
    class(mpi_control), intent(inout) :: this
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ctor(this) 
    class(mpi_control), intent(out) :: this
    type(par_file) :: pfile
    integer :: ierr,nprow,npcol,n

    call pfile%ctor('input_parameters','specs')

    !set proc rank
    call MPI_COMM_RANK(MPI_COMM_WORLD, this%rank, ierr)
    
    !set number of procs
    call MPI_COMM_SIZE(MPI_COMM_WORLD, this%nprocs, ierr)
    
    !check read-in number of procs
    call pfile%read_parameter(nprow,'nprow')
    call pfile%read_parameter(npcol,'npcol')
    if ((this%nprocs /= nprow*npcol).AND.(IS_PAR)) then
      call abort_run('Number of read in procs not equal to MPI_SIZE')
    end if
    
    call pfile%read_parameter(n,'N')
    if (this%nprocs>n) then
      call abort_run('nprocs > N found in mpi_control')
    end if
 
  end subroutine
!========================================================================================!

end module mpi_control_mod