module forcing_mod

  use dsph_procs_mod
  
  implicit none
  
  integer, parameter :: HTURB_DECAYING = 0
  integer, parameter :: HTURB_FORCING  = 1

  type, public :: forcing
    
    !turb. forcing matrices
    type(cdmatrix), allocatable, dimension(:) :: f0

    !resultant forcing matrix
    type(cdmatrix) :: f
    
    !auxiliary matrix
    type(cdmatrix) :: q
    
    !forcing sph index
    integer :: lf

    !forcing box size
    integer :: dlf

    !selector forcing
    integer :: forc_type
    
    contains
    
    procedure, public :: delete
    procedure, public :: ctor
    procedure, public :: update
    
  end type
  
  private :: delete,&
             ctor,&
             update

  private :: update_hturb_forcing_rphase,&
             update_hturb_forcing_dW         

contains

!========================================================================================!
  subroutine delete(this) 
    class(forcing), intent(inout) :: this

    call this%q%delete()
    call this%f%delete(delete_mpi=.FALSE.)
    if (allocated(this%f0)) then
      deallocate(this%f0)
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ctor(this,lap,w)
    class(forcing), intent(out) :: this
    type(laplacian), intent(inout) :: lap
    type(cdmatrix), intent(in) :: w
    type(par_file) :: pfile
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(this%forc_type,'forcing')
    
    call this%q%ctor(lap%mpic,BLOCK_COLUMN,'q',w%n,1,w%nprocs)
    
    select case(this%forc_type)
      case(HTURB_FORCING)
        call init_hturb_forcing(lap,w,this%lf,this%dlf,this%f,this%f0)
      case(HTURB_DECAYING)
        this%f = w
        call this%f%rename('f')
        this%f%m = (0.d0,0.d0)
      case default
        call abort_run('Wrong forcing type found')
    end select

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine update(this,dt)
    class(forcing), intent(inout) :: this
    real(double_p), intent(in) :: dt
    integer :: nmat,i,lfi,lf,dlf
    
    if (this%forc_type==HTURB_DECAYING) then
      return
    end if
    
    call this%f%set_to_zero()
    
    nmat = size(this%f0)
    
    lf = this%lf
    dlf = this%dlf
    
    do i=1,nmat
      
      lfi = lf-dlf+i-1
      
      !call update_hturb_forcing_rphase(this%f,this%q,this%f0(i),lfi,dlf,dt)
      call update_hturb_forcing_dW(this%f,this%q,this%f0(i),lfi,dlf,dt)
      
    end do
    
    call this%f%make_skewh()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine update_hturb_forcing_rphase(f,q,f0,lf,dlf,dt)
    type(cdmatrix), intent(inout) :: f,q
    type(cdmatrix), intent(in) :: f0
    integer, intent(in) :: lf,dlf
    real(double_p), intent(in) :: dt
    real(double_p), allocatable, dimension(:) :: r
    integer :: j,n,ierror,gcs,m,row,ncol,l_max
    complex(double_p) :: aux
    
    n=f%n
    call allocate_array(r,0,lf)
    
    if (IS_MASTER) then
      call random_number(r)
    end if
    call MPI_Bcast(r,size(r),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
    
    r=r*2.d0*pi
    
    call q%copy_values(f0)
    
    gcs=get_gidx(1,q%ncblk,q%pcol,q%npcol)
    ncol=q%ncol
    
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,q,r,gcs,ncol,lf,dt) &
    !$OMP PRIVATE(m,row,j,aux)
    do m=0,lf
      
      row = gcs + m
      
      if (row <= n) then
      
        j=1
        do
        
          if ((row>n).OR.(j>ncol)) then
            exit
          end if
          
          if (m==0) then
            if (r(m)>=0.5d0) then
              aux%re = sqrt(2.d0)
            else
              aux%re = -sqrt(2.d0)
            end if
            aux%im = 0.d0
          else
            aux%re = cos(r(m))
            aux%im = sin(r(m))
          end if
          
          q%m(row,j) = q%m(row,j)*aux*dt
          
          j = j + 1
          row = row + 1
        
        end do
        
      end if
      
    end do
    !$OMP END PARALLEL DO

    call q%column_to_cyclic(f%ptrm)
    call f%sum(f,f%ptrm)
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine update_hturb_forcing_dW(f,q,f0,lf,dlf,dt)
    type(cdmatrix), intent(inout) :: f,q
    type(cdmatrix), intent(in) :: f0
    integer, intent(in) :: lf,dlf
    real(double_p), intent(in) :: dt
    real(double_p), allocatable, dimension(:) :: r1,r2,u1,u2
    integer :: j,n,ierror,gcs,m,row,ncol
    complex(double_p) :: aux
    
    n=f%n
    call allocate_array(r1,0,lf)
    call allocate_array(r2,0,lf)
    call allocate_array(u1,0,lf)
    call allocate_array(u2,0,lf)
    
    if (IS_MASTER) then
      call random_number(u1)
      call random_number(u2)
    end if

    call MPI_Bcast(u1,size(u1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
    call MPI_Bcast(u2,size(u2),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
    
    r1=sqrt(-2.d0*log(u1))*cos(2.d0*pi*u2)*sqrt(dt)
    r2=sqrt(-2.d0*log(u1))*sin(2.d0*pi*u2)*sqrt(dt)
    
    call q%copy_values(f0)
    
    gcs=get_gidx(1,q%ncblk,q%pcol,q%npcol)
    ncol=q%ncol
    
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,q,u1,r1,r2,gcs,ncol,lf,dt) &
    !$OMP PRIVATE(m,row,j,aux)
    do m=0,lf
      
      row = gcs + m
      
      if (row <= n) then
      
        j=1
        do
        
          if ((row>n).OR.(j>ncol)) then
            exit
          end if
          
          !case m==0
          if (m==0) then
            aux%re = sqrt(2.d0)*r1(m)
            aux%im = 0.d0
          else
            aux%re = r1(m)
            aux%im = r2(m)
          end if
          
          q%m(row,j) = q%m(row,j)*aux
          
          j = j + 1
          row = row + 1
        
        end do
        
      end if
      
    end do
    !$OMP END PARALLEL DO

    call q%column_to_cyclic(f%ptrm)
    call f%sum(f,f%ptrm)
    
  end subroutine
!========================================================================================!

end module forcing_mod