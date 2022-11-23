PROGRAM main

  use time_mod
  use qg_operator_mod
  use integrators_mod

  IMPLICIT NONE
  
  !laplacian operator
  type(laplacian) :: lap
  
  !QG operator
  type(qg_operator) :: qg_op
  
  !convection integrator
  class(integrator), allocatable :: integ
  
  type(time) :: run_time
  type(mpi_control) :: mpic
  
  !vorticity and stream matrix
  type(cdmatrix) :: w,psi

  !fields output folder
  character(len=:), allocatable :: fields_dir
  
  real(double_p) :: ts,te
  
  ! ------------ SOLVER MAIN ------------  !
  call init_run()
  
  ts = MPI_Wtime()
  
  call init_solver()

  call solve_main()
  
  call delete_solver()
  
  te = MPI_Wtime()
  
  call info_run_end(ts,te)
  
  call end_run()
  ! -------------------------------------  !
  
contains

!========================================================================================!
  subroutine solve_main()
    
    !cycle time loop
    do while (run_time%loop())

      ts = MPI_Wtime()

      call integ%solve(qg_op,w,psi)
       
      if (run_time%output()) then
        call write_out()
      end if
      
      te = MPI_Wtime()
    
      call info_run_cpu_time(ts,te)
    
    end do

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_solver()
    
    call set_blas_sigle_thread()

    !call random_seed()
    
    call mpic%ctor()
    call w%ctor(mpic,BLOCK_CYCLIC,'w',read_pgrid=.TRUE.)
  
    if (IS_MASTER) then
      write(*,'(A,I4)') 'Threads:     ', N_THREADS
      write(*,'(A,I4)') 'MPI procs:   ', w%nprocs
      write(*,'(A,I4)') 'Total cores: ', w%nprocs*N_THREADS
    end if
    
    fields_dir = 'fields'
    
    call lap%ctor(mpic,w)
    call qg_op%ctor(mpic,lap,w) 
    call run_time%ctor()
    
    psi = w
    call psi%rename('psi')
    
    call allocate_integrator(integ,w,run_time%dt)
    
    !init vorticity
    if (w%is_stored(run_time%output_dir,fields_dir)) then
    
      call w%read_from_disk(run_time%output_dir,fields_dir)
      if (IS_MASTER) then
        write(*,'(A)') 'Vorticity read from disk'
      end if
      
    else

      call compute_ic(lap,w)
      
      if (IS_MASTER) then
        write(*,'(A)') 'Vorticity initialized by computing i.c.'
      end if
    
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine delete_solver()

    call deallocate_integrator(integ)
    call psi%delete(delete_mpi=.FALSE.)
    call run_time%delete()
    call lap%delete()
    call qg_op%delete()
    call w%delete()
    call mpic%delete()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine write_out()

    call run_time%write_out(fields_dir)
        
    !write vorticity
    call w%write_to_disk(run_time%output_dir,fields_dir)
        
    !write stream-function
    call psi%copy_values(w)
    call qg_op%apply_inv(psi)
    call psi%write_to_disk(run_time%output_dir,fields_dir)
        
    !write sph coefficients of vorticity
    call compute_sph_coeff(lap,w,run_time%output_dir,fields_dir)
        
    !write sph coefficients of vorticity
    call compute_sph_coeff(lap,psi,run_time%output_dir,fields_dir)

    !write discrete diffeomorphism (if computed)
    call integ%print_gmap(run_time%output_dir,fields_dir)

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

END PROGRAM main