PROGRAM main

  USE solver_mod

  IMPLICIT NONE
  
  type(solver) :: qturb
  real(double_p) :: ts,te
  type(par_file) :: pfile
  integer :: flow
  
  call init_run()
  
  ts = MPI_Wtime()
  
  call pfile%ctor('input_parameters','specs')
  call pfile%read_parameter(flow,'flow')
  
  call qturb%ctor()
  
  select case(flow)
    case(EULER)
      call qturb%solve_euler()
    case(H_TURB)
      call qturb%solve_hturb()
    case(D_TURB)
      call qturb%solve_dturb()
    case(QG_FLOW)
      call qturb%solve_qg_flow()
    case default
      call abort_run('Wrong flow solver')
  end select
  
  te = MPI_Wtime()
  
  call info_run_end(ts,te)
  
  call end_run()

END PROGRAM main