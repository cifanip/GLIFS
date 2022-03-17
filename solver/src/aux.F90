module aux_mod

  use allocate_arrays_mod

  implicit none

  real(double_p), parameter :: pi=4.d0*datan(1.d0)
  
contains

!========================================================================================!
  function ic_gen_test(m,l,l_cut) result(w)
    integer, intent(in) :: m,l,l_cut
    complex(double_p) :: w,w0,im
    real(double_p) :: k
    
    if (l>l_cut) then
      w=(0.d0,0.d0)
      return
    end if
    
    w0=(1.d0,1.d0)
    
    k=sqrt(dble(l*l+m*m))
    w=w0/(abs(w0)*k)
    
    if (m==0) then
      w%im = 0.d0
    end if

  end function
!========================================================================================!

!========================================================================================!
  function ic_gen_hturb(m,l,l_min,l_max) result(w)
    integer, intent(in) :: m,l,l_min,l_max
    real(double_p), dimension(2) :: r
    real(double_p) :: w_ref,w_mag
    complex(double_p) :: w
    
    if ((l>l_max).OR.(l<l_min)) then
      w=(0.d0,0.d0)
      return
    end if
    
    w_ref = 1.d-2
    w_mag = sqrt(2.d0*w_ref/(2.d0*l+1.d0))
    
    !randomize module
    call random_number(r)
    r(1) = r(1) - 0.5d0
    w_mag = (1.d0 + r(1)*0.4d0)*w_mag
    
    !randomize phase
    w%re = w_mag*cos(2.d0*pi*r(2))
    w%im = w_mag*sin(2.d0*pi*r(2))
    
    if (m==0) then
      w%re = w_mag
      w%im = 0.d0
    end if

  end function
!========================================================================================!

!========================================================================================!
  function hturbf_gen(m,l,dlf,lf) result(f)
    integer, intent(in) :: m,l,dlf,lf
    real(double_p) :: a,b
    complex(double_p) :: f
    
    a = (l-lf)*(l-lf)
    b = 0.5d0*dlf*dlf/log(10.d0)
    
    f%re=exp(-a/b)
    f%im=0.d0

  end function
!========================================================================================!

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