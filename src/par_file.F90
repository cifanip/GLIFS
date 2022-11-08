module par_file_mod

  use prun_mod
  use kinds_mod
  use formats_mod
  
  implicit none
  
  type, public :: par_file
    private
  
    character(len=:), allocatable :: file_name
    character(len=:), allocatable :: file_path
    
    contains
    
    procedure, public :: ctor   => par_file_CTOR

    procedure :: read_par_int,&
                 read_par_bool,&
                 read_par_real
    generic, public :: read_parameter => read_par_int,&
                                         read_par_bool,&
                                         read_par_real
    
  end type
  
  private :: par_file_CTOR,&
             read_par_int,&
             read_par_bool,&
             read_par_real
  
  !outside type        
  private :: frstnb

contains

!========================================================================================!
  subroutine par_file_CTOR(this,file_name,file_dir)
    class(par_file), intent(out) :: this
    character(len=*), intent(in) :: file_name
    character(len=*), intent(in) :: file_dir 
    
    this%file_name = file_name
    this%file_path = file_dir//'/'//file_name
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine read_par_int(this,x,name,bcast) 
    class(par_file), intent(in) :: this
    integer, intent(out) :: x
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: bcast
    logical :: opt_bcast
    character(len=100) :: whole
    integer :: i1, i2, ios, ierror 
    
    if (present(bcast)) then
      opt_bcast = bcast
    else
      opt_bcast = .TRUE.
    end if
    
    if (IS_MASTER) then
    
      open(UNIT=IO_unit_,FILE=this%file_path,STATUS='OLD',ACTION='READ')
      do 
        read(IO_unit_,char_format_,IOSTAT=ios) whole
      
          if (ios /= 0) then
            call abort_run('Parameter '//name//' not found in file '//this%file_name)
          else
            !parse whole
            i1 = index(whole,' ')
              if (whole(1:i1-1) == name) then
                i2 = frstnb(whole(i1+1:),name)
                read(whole(i1+i2:),int_format_) x
                exit
              end if
          end if
      end do
      close(IO_unit_)
    
    end if
    
    if (opt_bcast) then
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_BCAST(x, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    end if
  
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine read_par_bool(this,x,name,bcast) 
    class(par_file), intent(in) :: this
    logical, intent(out) :: x
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: bcast
    logical :: opt_bcast
    character(len=100) :: whole
    integer :: i1, i2, ios, ierror 
    
    if (present(bcast)) then
      opt_bcast = bcast
    else
      opt_bcast = .TRUE.
    end if
    
    if (IS_MASTER) then

      open(UNIT=IO_unit_,FILE=this%file_path,STATUS='OLD',ACTION='READ')
      do 
        read(IO_unit_,char_format_,IOSTAT=ios) whole
      
          if (ios /= 0) then
            call abort_run('Parameter '//name//' not found in file '//this%file_name)
          else
            !parse whole
            i1 = index(whole,' ')
              if (whole(1:i1-1) == name) then
                i2 = frstnb(whole(i1+1:),name)
                read(whole(i1+i2:),logical_format_) x
                exit
              end if
          end if
      end do
      close(IO_unit_)
    
    end if
    
    if (opt_bcast) then
      call MPI_BCAST(x, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine read_par_real(this,x,name,bcast) 
    class(par_file), intent(in) :: this
    real(double_p), intent(out) :: x
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: bcast
    logical :: opt_bcast
    character(len=100) :: whole
    integer :: i1, i2, ios, ierror 
    
    if (present(bcast)) then
      opt_bcast = bcast
    else
      opt_bcast = .TRUE.
    end if
    
    if (IS_MASTER) then
    
      open(UNIT=IO_unit_,FILE=this%file_path,STATUS='OLD',ACTION='READ')
      do 
        read(IO_unit_,char_format_,IOSTAT=ios) whole
      
          if (ios /= 0) then
            call abort_run('Parameter '//name//' not found in file '//this%file_name)
          else
            !parse whole
            i1 = index(whole,' ')
              if (whole(1:i1-1) == name) then
                i2 = frstnb(whole(i1+1:),name)
                read(whole(i1+i2:),double_format_) x
                exit
              end if
          end if
      end do
      close(IO_unit_)
    
    end if

    if (opt_bcast) then
      call MPI_BCAST(x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    end if
  
  end subroutine
!========================================================================================!

!========================================================================================!
  function frstnb(str,par) result(ind)
    character(len=*), intent(in) :: str, par
    integer :: ind, i
    
    do i=1,len(str)
      if (str(i:i) /= ' ') then
        ind = i
        return 
      end if
    end do
    
    call abort_run('Invalid entry while parsing parameter: ' // par)
  
  end function
!========================================================================================!

end module par_file_mod