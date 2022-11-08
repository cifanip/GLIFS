module ic_mod

  use vector_mod

  implicit none
  
  type, public :: ic
    private
    
    !i.c. of spherical harmonic coefficients
    type(vector), public :: ic_vec
    
    !range of spherical h. index l
    integer, public :: l_min,l_max
    
    contains
    
    procedure, public :: ic_gen_test
    procedure, public :: ic_gen_rnd
    procedure, public :: ic_gen_hturbf
    procedure, public :: ic_gen_Tbase
    
    procedure, public :: delete
  
  end type
  
  private :: ic_gen_test,&
             ic_gen_rnd,&
             ic_gen_hturbf,&
             ic_gen_Tbase,&
             delete

  
contains

!========================================================================================!
  subroutine delete(this) 
    class(ic), intent(inout) :: this
    
    call this%ic_vec%delete()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ic_gen_test(this,mpic,n)
    class(ic), intent(inout) :: this
    type(mpi_control), intent(in) :: mpic
    integer, intent(in) :: n
    complex(double_p) :: w,w0
    real(double_p) :: k
    integer :: l,m,ic_size,gidx,lidx,r_min,r_max,l_cut
    type(par_file) :: pfile
    
    ic_size=sph_size(n)
    call this%ic_vec%ctor(mpic,'ic_vec',ic_size,mpic%nprocs)
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(l_cut,'l_cut')

    if (l_cut>n-1) then
      call abort_run('l_cut > N-1')
    end if
    
    this%l_min = 1
    this%l_max = l_cut
    
    w0=(1.d0,1.d0)
    w=(0.d0,0.d0)
    
    r_min = get_gidx(1,this%ic_vec%nblk,this%ic_vec%prow,this%ic_vec%nprocs)
    r_max = get_gidx(this%ic_vec%nrow,this%ic_vec%nblk,this%ic_vec%prow,&
                     this%ic_vec%nprocs)
    
    do l=1,n-1
    
      do m=0,l

        gidx = sph_idx(l,m)
        
        if ((gidx>=r_min).AND.(gidx<=r_max)) then
        
          ! --------
          if (l>l_cut) then
            w=(0.d0,0.d0)
            cycle
          end if
        
          k=sqrt(dble(l*l+m*m))
          w=w0/(abs(w0)*k)
    
          if (m==0) then
            w%im = 0.d0
          end if
          ! --------
          
          lidx = get_lidx(gidx,this%ic_vec%nblk,this%ic_vec%prow,this%ic_vec%nprocs)
        
          this%ic_vec%v(lidx)%re=w%re
          this%ic_vec%v(lidx)%im=w%im            
        
        end if
        
      end do
    
    end do
    
    call this%ic_vec%all_gather()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ic_gen_rnd(this,mpic,n)
    class(ic), intent(inout) :: this
    type(mpi_control), intent(in) :: mpic
    integer, intent(in) :: n
    complex(double_p) :: w
    real(double_p) :: w_ref,w_mag
    integer :: l,m,ic_size,gidx,lidx,r_min,r_max,l_min,l_max
    type(par_file) :: pfile
    real(double_p), allocatable, dimension(:,:) :: rnd_num
    
    ic_size=sph_size(n)
    call this%ic_vec%ctor(mpic,'ic_vec',ic_size,mpic%nprocs)
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(l_min,'l_min')
    call pfile%read_parameter(l_max,'l_max')

    if (l_max>n-1) then
      call abort_run('l_max > N-1')
    end if
    if (l_min<1) then
      call abort_run('l_min < 1')
    end if

    this%l_min = l_min
    this%l_max = l_max
    
    !array of random numbers
    call allocate_array(rnd_num,1,2,1,this%ic_vec%nrow)
    call random_number(rnd_num)
    
    r_min = get_gidx(1,this%ic_vec%nblk,this%ic_vec%prow,this%ic_vec%nprocs)
    r_max = get_gidx(this%ic_vec%nrow,this%ic_vec%nblk,this%ic_vec%prow,&
                     this%ic_vec%nprocs)
                     
    w_ref = 4.d0
    w = (0.d0,0.d0)
    
    do l=1,n-1
    
      do m=0,l

        gidx = sph_idx(l,m)
        lidx = get_lidx(gidx,this%ic_vec%nblk,this%ic_vec%prow,this%ic_vec%nprocs)
        
        if ((gidx>=r_min).AND.(gidx<=r_max)) then
        
          ! --------
          if ((l>l_max).OR.(l<l_min)) then
            w=(0.d0,0.d0)
            cycle
          end if
    
          w_mag = sqrt(2.d0*w_ref/((2.d0*l+1.d0)))
    
          !randomize module
          rnd_num(1,lidx) = rnd_num(1,lidx) - 0.5d0
          w_mag = (1.d0 + rnd_num(1,lidx)*0.4d0)*w_mag
    
          !randomize phase
          w%re = (w_mag/sqrt(2.d0))*cos(2.d0*pi*rnd_num(2,lidx))
          w%im = (w_mag/sqrt(2.d0))*sin(2.d0*pi*rnd_num(2,lidx))
    
          if (m==0) then
            if (rnd_num(1,lidx)>=0.d0) then
              w%re = w_mag
            else
              w%re = -w_mag
            end if
              w%im = 0.d0
          end if
          ! --------
        
          this%ic_vec%v(lidx)%re=w%re
          this%ic_vec%v(lidx)%im=w%im            
        
        end if
        
      end do
    
    end do
    
    call this%ic_vec%all_gather()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ic_gen_hturbf(this,mpic,n,lf,lref,dlf)
    class(ic), intent(inout) :: this
    type(mpi_control), intent(in) :: mpic
    integer, intent(in) :: n,lf,lref,dlf
    integer :: l,m,gidx,lidx,r_min,r_max,ic_size
    real(double_p) :: a,b,f_mag
    type(par_file) :: pfile
    
    ic_size=sph_size(n)
    call this%ic_vec%ctor(mpic,'ic_vec',ic_size,mpic%nprocs)
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(f_mag,'fmag')

    this%l_min = lref
    this%l_max = lref
    
    r_min = get_gidx(1,this%ic_vec%nblk,this%ic_vec%prow,this%ic_vec%nprocs)
    r_max = get_gidx(this%ic_vec%nrow,this%ic_vec%nblk,this%ic_vec%prow,&
                     this%ic_vec%nprocs)
    
    if (dlf>0) then
      b = 0.5d0*dlf*dlf/log(10.d0)
    else
      b = 1.d0
    end if
    
    do l=this%l_min,this%l_max
    
      do m=0,l

        gidx = sph_idx(l,m)
        lidx = get_lidx(gidx,this%ic_vec%nblk,this%ic_vec%prow,this%ic_vec%nprocs)
        
        if ((gidx>=r_min).AND.(gidx<=r_max)) then
        
          ! --------    
          a = (l-lf)*(l-lf)
    
          this%ic_vec%v(lidx)%re=f_mag*exp(-a/b)
          this%ic_vec%v(lidx)%im=0.d0
          ! --------        
        
        end if
        
      end do
    
    end do
    
    call this%ic_vec%all_gather()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine ic_gen_Tbase(this,mpic,n,l,m)
    class(ic), intent(inout) :: this
    type(mpi_control), intent(in) :: mpic
    integer, intent(in) :: n,l,m
    integer :: ic_size,r_min,r_max,gidx,lidx
    
    ic_size=sph_size(n)
    call this%ic_vec%ctor(mpic,'ic_vec',ic_size,mpic%nprocs)

    if (l>n-1) then
      call abort_run('l > N-1')
    end if
    if (l<1) then
      call abort_run('l < 1')
    end if
    if (m>l) then
      call abort_run('m > l')
    end if

    this%l_min = l
    this%l_max = l
    
    r_min = get_gidx(1,this%ic_vec%nblk,this%ic_vec%prow,this%ic_vec%nprocs)
    r_max = get_gidx(this%ic_vec%nrow,this%ic_vec%nblk,this%ic_vec%prow,&
                     this%ic_vec%nprocs)

    gidx = sph_idx(l,m)
    lidx = get_lidx(gidx,this%ic_vec%nblk,this%ic_vec%prow,this%ic_vec%nprocs)
        
    if ((gidx>=r_min).AND.(gidx<=r_max)) then
    
      this%ic_vec%v(lidx)%re=1.d0
      this%ic_vec%v(lidx)%im=0.d0
        
    end if
    
    call this%ic_vec%all_gather()

  end subroutine
!========================================================================================!
  
end module ic_mod