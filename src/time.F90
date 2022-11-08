module time_mod

  use allocate_arrays_mod
  use par_file_mod
  
  implicit none
  
  type, public :: time
    private
  
    real(double_p), public :: t
    real(double_p) :: Tf
    real(double_p), public :: dt
    real(double_p) :: dtout,tout
    integer :: write_interval
    integer :: write_iter
    
    !counter time iterations
    integer, public :: iter
    
    !input-output folder
    integer, public :: input_dir
    integer, public :: output_dir

    contains
      
    procedure, public :: delete
    procedure, public :: ctor
    procedure, public :: loop
    procedure, public :: output
    procedure, public :: write_out
    procedure :: update
    procedure :: info
    procedure :: write_folder
    procedure :: init_time_level
    
  end type
  
  private :: delete,&
             ctor,&
             update,&
             loop,&
             info,&
             output,&
             write_out,&
             write_folder,&
             init_time_level

contains

!========================================================================================!
  subroutine delete(this)
    class(time), intent(inout) :: this
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ctor(this)
    class(time), intent(out) :: this
    type(par_file) :: pfile
    
    call pfile%ctor('input_parameters','specs')

    call pfile%read_parameter(this%Tf,'Tf')
    call pfile%read_parameter(this%dt,'dt')
    call pfile%read_parameter(this%input_dir,'input_folder')
    call pfile%read_parameter(this%write_interval,'write_interval')
    call pfile%read_parameter(this%dtout,'dtout')
    
    this%iter = 0
    this%write_iter = 0
    this%output_dir = this%input_dir

    call this%init_time_level() 

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine update(this)
    class(time), intent(inout) :: this
        
    this%iter = this%iter + 1
    this%write_iter = this%write_iter + 1
  
    this%t = this%t + this%dt
    
  end subroutine
!========================================================================================!

!========================================================================================!
  function loop(this) result(is_run)
    class(time), intent(inout) :: this
    logical :: is_run
    real(double_p) :: small
        
    small = 1.d-1*this%dt
        
    call this%update()
        
    if (this%t >= this%Tf + small) then
      is_run = .FALSE.
      this%t=this%t-this%dt
    else
      call this%info()
      is_run = .TRUE.
    end if
    
  end function
!========================================================================================!

!========================================================================================!
  function output(this) result(is_output)
    class(time), intent(inout) :: this
    logical :: is_output
    real(double_p) :: t,tout,dt,small
        
    t=this%t
    tout=this%tout
    dt=this%dt
    small=dt/4.d0

    if ((this%write_interval == this%write_Iter).OR. &
        ((t >= tout).OR.((tout >= t-small).AND.(tout <= t+small)))) then
        is_output = .TRUE.
        this%write_iter = 0
        this%tout = this%tout + this%dtout
    else
        is_output = .FALSE.
    end if

  end function
!========================================================================================!

!========================================================================================!
  subroutine write_out(this,fdir)
    class(time), intent(inout) :: this
    character(len=*), intent(in) :: fdir
      
    call this%write_folder(fdir)
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine write_folder(this,fdir)
    class(time), intent(inout) :: this
    character(len=*), intent(in) :: fdir
    character(len=10) :: char_tout
    character(len=20) :: dir
    integer :: ierror
        
    this%output_dir = this%output_dir + 1
        
    write(char_tout,int_format_) this%output_dir
    dir = trim(adjustl(char_tout))//'/'//fdir
        
    if (IS_MASTER) then   
          
      !mkdir
      call system('mkdir -p ./' // adjustl(trim(dir)))
          
      call system(adjustl('touch ./'//trim(char_tout)//'/info_restart'))

      !write time
      open(UNIT=IO_unit_,FILE=trim(adjustl(char_tout))//'/info_restart',&
           STATUS='REPLACE',ACTION='WRITE')
        write(IO_unit_,double_format_) this%t
      close(IO_unit_)

    end if
        
    call mpi_barrier(MPI_COMM_WORLD,ierror)
        
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine info(this)
    class(time), intent(inout) :: this

    if (IS_MASTER) then
      write(*,*) ''
      
      write(*,*) 'SOLVING FOR TIME:'
      
      write(*,'(A,'//double_format_(2:10)//')') '    t  =', this%t
      write(*,'(A,'//output_format_(2:9)//')')  '    dt =  ', this%dt

    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_time_level(this)
    class(time), intent(inout) :: this
    character(len=10) :: dir_name  
    integer :: ierror
      
    if (IS_MASTER) then
      
      write(dir_name,int_format_) this%input_dir
      
      if (this%input_dir == 0) then
        this%t = 0.d0
        this%tout = this%dtout
      else
        open(UNIT=IO_unit_,FILE=adjustl(trim(dir_name)//'/info_restart'),&
             STATUS='old',ACTION='read')
          read(IO_unit_,double_format_) this%t
        close(IO_unit_)  
        this%tout = this%dtout+this%t
      end if
    end if
      
    call MPI_BCAST(this%t, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(this%tout, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

  end subroutine
!========================================================================================!
  
end module time_mod