module solver_mod

  use lap_blk_mod
  use time_mod
  
  implicit none
  
  type, public :: solver
    private
    
    type(mpi_control) :: mpic
    type(cdmatrix) :: w
    type(lap_blk) :: lap
    type(time) :: run_time
    
    !work matrices
    type(cdmatrix) :: wt0,wt,psi,psi_wt,psi_wth,&
                      psi_wt_psi,delta_w
    
    !iteration implicit step
    real(double_p) :: tol
    integer :: max_iter
    
    !fields output folder
    character(len=:), allocatable, public :: fields_dir
    
    !type of flow
    integer :: flow
    
    !viscosity
    real(double_p) :: nu
    
    !damping
    real(double_p) :: alpha
    
    !magnitude forcing
    real(double_p) :: fmag
    
    !turb. forcing matrix
    type(cdmatrix) :: f0,f
    
    !forcing sph index
    integer :: lf
    
    !theta-method hturb
    real(double_p) :: theta

    contains
      
    procedure, public :: delete
    procedure, public :: ctor
    procedure, public :: solve_euler
    procedure, public :: solve_hturb
    procedure :: isomp
    procedure :: isomp_fpiter
    procedure :: isomp_expl

  end type
  
  private :: delete,&
             ctor,&
             solve_euler,&
             solve_hturb,&
             isomp_fpiter,&
             isomp_expl,&
             isomp,&
             info_run_cpu_time
             
  !outside of type
  private :: isomp_fpiter_add,&
             isomp_expl_add

contains

!========================================================================================!
  subroutine delete(this)
    class(solver), intent(inout) :: this
    
    call this%mpic%delete()
    call this%w%delete()
    call this%lap%delete()
    call this%run_time%delete()
    
    call this%wt0%delete(delete_mpi=.FALSE.)
    call this%wt%delete(delete_mpi=.FALSE.)
    call this%psi%delete(delete_mpi=.FALSE.)
    call this%psi_wt%delete(delete_mpi=.FALSE.)
    call this%psi_wth%delete(delete_mpi=.FALSE.)
    call this%psi_wt_psi%delete(delete_mpi=.FALSE.)
    call this%delta_w%delete(delete_mpi=.FALSE.)
    if (this%flow==H_TURB) then
      call this%f%delete(delete_mpi=.FALSE.)
      call this%f0%delete(delete_mpi=.FALSE.)
    end if
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ctor(this)
    class(solver), intent(out) :: this
    type(par_file) :: pfile
    
    call set_sigle_thread()
    
    call this%mpic%ctor()
    call this%w%ctor(this%mpic,BLOCK_CYCLIC,'w',read_pgrid=.TRUE.)
    
    if (IS_MASTER) then
      write(*,'(A,I4)') 'Threads:     ', N_THREADS
      write(*,'(A,I4)') 'MPI procs:   ', this%w%nprocs
      write(*,'(A,I4)') 'Total cores: ', this%w%nprocs*N_THREADS
    end if

    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(this%tol,'tol')
    call pfile%read_parameter(this%max_iter,'max_iter') 
    call pfile%read_parameter(this%flow,'flow') 
    
    if (this%flow==H_TURB) then
      call pfile%read_parameter(this%nu,'nu')
      call pfile%read_parameter(this%alpha,'alpha')
      call pfile%read_parameter(this%fmag,'fmag')
      call pfile%read_parameter(this%lf,'lf')
      call pfile%read_parameter(this%theta,'theta')
      !account for term \nu \omega
      this%alpha = this%alpha - this%nu
      if (this%lf>this%w%n-1) then
        call abort_run('lf > N-1')
      end if
    end if
    
    this%fields_dir = 'fields'

    call this%lap%ctor(this%mpic,this%w,this%flow)
    call this%run_time%ctor()
    
    this%wt0        = this%w
    this%wt         = this%w
    this%psi        = this%w
    this%psi_wt     = this%w
    this%psi_wth    = this%w
    this%psi_wt_psi = this%w
    this%delta_w    = this%w

    if (this%flow==H_TURB) then
      this%f  = this%w
      this%f0 = this%w
      call this%lap%init_hturb_forcing(this%fmag,this%lf,this%f0)
    end if
    
    !init vorticity
    if (this%w%is_stored(this%run_time%output_dir,this%fields_dir)) then
    
      call this%w%read_from_disk(this%run_time%output_dir,this%fields_dir)
      if (IS_MASTER) then
        write(*,'(A)') 'Vorticity read from disk'
      end if
      
    else
      
      if (this%flow==EULER) then
        call this%lap%compute_ic(this%w)        
      end if

      if (this%flow==H_TURB) then
        this%w%m=(0.d0,0.d0)
      end if

      if (IS_MASTER) then
        write(*,'(A)') 'Vorticity initialized from given i.c.'
      end if

    end if
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine solve_hturb(this)
    class(solver), intent(inout) :: this
    real(double_p) :: ts,te,l2W

    do while (this%run_time%loop())

      ts = MPI_Wtime()
      
      call this%lap%apply_hturb(this%w,&
                                this%f,&
                                this%f0,&
                                this%lf,&
                                0.5d0*this%run_time%dt,&
                                this%nu,&
                                this%alpha,&
                                this%theta)

     call this%isomp()

     call this%lap%apply_hturb(this%w,&
                               this%f,&
                               this%f0,&
                               this%lf,&
                               0.5d0*this%run_time%dt,&
                               this%nu,&
                               this%alpha,&
                               this%theta)                
       
      if (this%run_time%output()) then
        call this%run_time%write_out(this%fields_dir)
        call this%w%write_to_disk(this%run_time%output_dir,this%fields_dir)
        call this%lap%compute_sph_coeff(this%w,this%run_time%output_dir,this%fields_dir)
      end if
      
      l2W = this%w%l2_norm()
      if (IS_MASTER) then
        write(*,'(A,'//double_format_(2:10)//')') '    l2(W): ', l2W
      end if 
      
      te = MPI_Wtime()
    
      call info_run_cpu_time(ts,te)

    end do

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine solve_euler(this)
    class(solver), intent(inout) :: this
    real(double_p) :: ts,te
    
    !cycle time loop
    do while (this%run_time%loop())

      ts = MPI_Wtime()

      call this%isomp()
       
      if (this%run_time%output()) then
        call this%run_time%write_out(this%fields_dir)
        call this%w%write_to_disk(this%run_time%output_dir,this%fields_dir)
        call this%lap%compute_sph_coeff(this%w,this%run_time%output_dir,this%fields_dir)
      end if
      
      te = MPI_Wtime()
    
      call info_run_cpu_time(ts,te)
    
    end do

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine isomp(this,dt_opt)
    class(solver), intent(inout) :: this
    real(double_p), intent(in), optional :: dt_opt
    real(double_p) :: ts,te
    
    ts = MPI_Wtime()
    
    call this%isomp_fpiter(dt_opt)
    call this%isomp_expl(dt_opt)   
      
    te = MPI_Wtime() 

    if (IS_MASTER) then
      write(*,'(A,'//output_format_(2:9)//')') '    ISOMP CPU TIME: ', te-ts
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine isomp_expl(this,dt_opt)
    class(solver), intent(inout) :: this
    real(double_p), intent(in), optional :: dt_opt
    real(double_p) :: dt
    
    if (present(dt_opt)) then
      dt = dt_opt
    else
      dt = this%run_time%dt
    end if

    !apply inverse laplacian
    call this%psi%copy_values(this%wt)
    call this%lap%apply_inv(this%psi)
    
    !(\Delta^{-1} Wt)Wt
    call this%psi_wt%multiply(this%psi,this%wt)
    !call this%psi_wt%multiply_h(this%wt,this%psi)

    !Wt(\Delta^{-1} Wt)
    call this%psi_wth%hconj(this%psi_wt)
    
    !(\Delta^{-1} Wt)Wt(\Delta^{-1} Wt)
    call this%psi_wt_psi%multiply(this%psi_wt,this%psi)
    !call this%psi_wt_psi%multiply_h(this%psi,this%psi_wt)
    
    !add-up
    call isomp_expl_add(this%wt%m,             &
                            this%psi%m,        &
                            this%psi_wt%m,     &
                            this%psi_wth%m,    &
                            this%psi_wt_psi%m, &
                            dt,                &
                            this%wt%nrow,      &
                            this%wt%ncol,      &
                            this%w%m)
                            
    call this%w%make_skewh()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine isomp_fpiter(this,dt_opt)
    class(solver), intent(inout) :: this
    real(double_p), intent(in), optional :: dt_opt
    real(double_p) :: dt,err
    integer :: iter
    logical :: converged

    if (present(dt_opt)) then
      dt = dt_opt
    else
      dt = this%run_time%dt
    end if
    
    iter = 0
    converged = .FALSE.
    
    call this%wt0%copy_values(this%w)
    
    do
    
      !apply inverse laplacian
      call this%psi%copy_values(this%wt0)
      call this%lap%apply_inv(this%psi)

      !(\Delta^{-1} Wt)Wt
      call this%psi_wt%multiply(this%psi,this%wt0)
      !call this%psi_wt%multiply_h(this%wt0,this%psi)
    
      !Wt(\Delta^{-1} Wt)
      call this%psi_wth%hconj(this%psi_wt)
    
      !(\Delta^{-1} Wt)Wt(\Delta^{-1} Wt)
      call this%psi_wt_psi%multiply(this%psi_wt,this%psi)
      !call this%psi_wt_psi%multiply_h(this%psi,this%psi_wt)
      
      !add-up
      call isomp_fpiter_add(this%w%m,          &
                            this%psi%m,        &
                            this%psi_wt%m,     &
                            this%psi_wth%m,    &
                            this%psi_wt_psi%m, &
                            dt,                &
                            this%wt%nrow,      &
                            this%wt%ncol,      &
                            this%wt%m)
      
      !compute infinity norm
      call this%delta_w%subtract(this%wt,this%wt0)
      err = this%delta_w%inf_norm()
    
      if (err <= this%tol) then
        converged = .TRUE.
        if (IS_MASTER) then
          write(*,'(A,I3,A,'//output_format_(2:9)//')') '    converged: iter ',iter,&
                                                   ';  err ',err
        end if
        return
      end if
      
      call this%wt0%copy_values(this%wt)
      
      iter = iter + 1
    
      if (iter == this%max_iter) then
        call abort_run('Maximum number of iterations reached in isomp_fpiter')
      end if
    
    end do

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine isomp_fpiter_add(wn,psi,psi_wt,psi_wth,psi_wt_psi,dt,nrow,ncol,wt)
    complex(double_p), allocatable, dimension(:,:), intent(in) ::  wn,psi,psi_wt,&
                                                                   psi_wth,psi_wt_psi
    real(double_p), intent(in) :: dt
    integer, intent(in) :: nrow,ncol
    complex(double_p), allocatable, dimension(:,:), intent(inout) :: wt
    real(double_p) x,y
    integer :: i,j
    
    x=0.5d0*dt
    y=0.25d0*dt*dt

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(ncol,nrow,wt,wn,x,y,psi_wt,psi_wth,psi_wt_psi) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        wt(i,j) = wn(i,j) + x*(psi_wt(i,j)-psi_wth(i,j)) + y*psi_wt_psi(i,j)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine isomp_expl_add(wt,psi,psi_wt,psi_wth,psi_wt_psi,dt,nrow,ncol,w)
    complex(double_p), allocatable, dimension(:,:), intent(in) ::  wt,psi,psi_wt,&
                                                                   psi_wth,psi_wt_psi
    real(double_p), intent(in) :: dt
    integer, intent(in) :: nrow,ncol
    complex(double_p), allocatable, dimension(:,:), intent(inout) :: w
    real(double_p) x,y
    integer :: i,j
    
    x=0.5d0*dt
    y=0.25d0*dt*dt

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(ncol,nrow,w,wt,x,y,psi_wt,psi_wth,psi_wt_psi) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        w(i,j) = wt(i,j) + x*(psi_wt(i,j)-psi_wth(i,j)) - y*psi_wt_psi(i,j)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine info_run_cpu_time(ts,te)
    real(double_p), intent(in) :: ts,te

    if (IS_MASTER) then
      write(*,'(A,'//output_format_(2:9)//')') '    TOTAL CPU TIME: ', te-ts
    end if
  
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine info_run_end(ts,te)
    real(double_p), intent(in) :: ts,te
  
    if (IS_MASTER) then
      write(*,*) ''
      write(*,'(A,'//output_format_(2:9)//')') &
          'EXIT RUN NORMAL. SIMULATION TIME: ', te-ts
    end if    
  
  end subroutine
!========================================================================================!
  
end module solver_mod