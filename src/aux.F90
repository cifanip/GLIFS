module aux_mod

  use allocate_arrays_mod

  implicit none

  real(double_p), parameter :: pi=4.d0*datan(1.d0)
  
contains

!========================================================================================!
  subroutine check_odd(n)
    integer, intent(in) :: n
    integer :: r

      r = modulo(n,2)
      if (r==0) then
        call abort_run('N cannot be even')
      end if

  end subroutine
!========================================================================================!

!========================================================================================!
  function sph_size(n) result(r)
    integer, intent(in) :: n
    integer :: q,r
    
    q=n-1
    r = q*(q+1)/2+q

  end function
!========================================================================================!

!========================================================================================!
  function sph_idx(l,m) result(r)
    integer, intent(in) :: l,m
    integer :: q,r
    
    q=l-1
    r=q*(q+1)/2+q+m+1

  end function
!========================================================================================!

!========================================================================================!
  subroutine set_blas_multiple_threads()
#ifdef CMP_INTEL
    call mkl_set_num_threads(N_THREADS)
#endif
#ifdef CMP_GNU
    call openblas_set_num_threads(N_THREADS)
#endif
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine set_blas_sigle_thread()
#ifdef CMP_INTEL
    call mkl_set_num_threads(1)
#endif
#ifdef CMP_GNU
    call openblas_set_num_threads(1)
#endif
  end subroutine
!========================================================================================!
  
end module aux_mod