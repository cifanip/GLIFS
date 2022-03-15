module lap_blk_mod

  use rdmatrix_mod
  use cdmatrix_mod
  use vector_mod
  use flow_types_mod
  
  implicit none
  
  integer, parameter :: SPH_IC = 1
  integer, parameter :: SPH_F  = 2

  !auxiliary type (diagonal points info)
  type :: dinfo
    integer, allocatable, dimension(:,:) :: idx
    integer :: rank
  end type

  !auxiliary type (tridiagonal systems)
  type :: tdiags
    !used for apply inverse Laplacian
    real(double_p), allocatable, dimension(:) :: d
    complex(double_p), allocatable, dimension(:) :: e,rhs
    !used for apply Laplacian
    real(double_p), allocatable, dimension(:) :: dl,el
    !used for viscosity and damping
    real(double_p), allocatable, dimension(:) :: dn
    complex(double_p), allocatable, dimension(:) :: en
  end type

  type, public :: lap_blk
    private
    
    integer, public :: n
    type(rdmatrix) :: v
    type(mpi_control), pointer :: mpic => NULL()
    
    type(dinfo), allocatable, dimension(:) :: dw
    
    !MPI diagonals communication
    integer, allocatable, dimension(:) :: dmat_type,dbuff_type
    integer, allocatable, dimension(:) :: d_midx
    logical, allocatable, dimension(:,:) :: sr_diag
    complex(double_p), allocatable, dimension(:,:) :: d_buff
    
    type(tdiags), allocatable, dimension(:) :: tds
    
    type(cdmatrix) :: q
    
    contains

    procedure, public :: delete => delete_lap_blk
    procedure, public :: ctor   => lap_blk_CTOR
    
    procedure, public :: apply
    procedure, public :: apply_inv
    procedure, public :: apply_hturb
    procedure, public :: compute_ic
    procedure, public :: compute_sph_coeff
    procedure, public :: init_hturb_forcing
    procedure :: add_q_component
    procedure :: compute_vort_component
    procedure :: init_dinfo
    procedure :: init_diags
    procedure :: diags_to_buffer
    procedure :: buffer_to_diags
    procedure :: init_block_systems

  end type
  
  private :: delete_lap_blk,&
             compute_ic,&
             init_hturb_forcing,&
             add_q_component,&
             init_dinfo,&
             apply,&
             apply_inv,&
             apply_hturb,&
             compute_sph_coeff,&
             init_diags,&
             diags_to_buffer,&
             buffer_to_diags,&
             init_block_systems
  
  !outside of type
  private :: assemble_ic,&
             compute_tmat,&
             compute_eigv,&
             adjust_basis_sign,&
             compute_vort_component,&
             frobenius_prod,&
             store_factorization,&
             solve_inv_lap,&
             update_hturb_forcing_rphase,&
             update_hturb_forcing_dW

contains

!========================================================================================!
  subroutine delete_lap_blk(this) 
    class(lap_blk), intent(inout) :: this
    integer :: i,ierror

    if (this%v%is_allocated) then
      call this%v%delete()
    end if
    this%mpic => NULL()
    
    deallocate(this%dw)
    
    if (allocated(this%dmat_type)) then
      do i=1,size(this%dmat_type)
        call mpi_type_free(this%dmat_type(i),ierror)
      end do
    end if

    if (allocated(this%dbuff_type)) then
      do i=1,size(this%dbuff_type)
        call mpi_type_free(this%dbuff_type(i),ierror)
      end do
    end if
    
    call deallocate_array(this%d_midx)
    call deallocate_array(this%sr_diag)
    call deallocate_array(this%d_buff)
    
    deallocate(this%tds)

    call this%q%delete()
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine lap_blk_CTOR(this,mpic,w,flow) 
    class(lap_blk), intent(out) :: this
    type(mpi_control), intent(in), target :: mpic
    type(cdmatrix), intent(in) :: w
    integer, intent(in) :: flow
    type(par_file) :: pfile

    this%mpic => mpic
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(this%n,'N')
    
    call this%q%ctor(this%mpic,BLOCK_COLUMN,'q',w%n,1,w%nprocs)
    
    call this%init_dinfo(this%q)

    call this%init_diags()
    
    call this%init_block_systems(flow)

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine apply_hturb(this,w,f,f0,lf,dlf,dt,nu,alpha,theta)
    class(lap_blk), intent(inout) :: this
    type(cdmatrix), intent(inout) :: w,f
    type(cdmatrix), intent(in), allocatable, dimension(:) :: f0
    integer, intent(in) :: lf,dlf
    real(double_p) :: dt,nu,alpha,theta,r,s,ts,te
    integer :: i,j,k,n,n_diags,nrow,ncol
    
    ts = MPI_Wtime()
    
    n=this%n
    n_diags=size(this%d_midx)
    
    call update_hturb_forcing(f,this%q,f0,lf,dlf,dt)
    
    call f%ptrm%copy_values(w)
    call this%apply(w)
    
    nrow=w%nrow
    ncol=w%ncol
    
    r=theta*dt*nu
    s=theta*dt*alpha

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nrow,ncol,dt,r,s,w,f) &
    !$OMP PRIVATE(i,j)
    do j=1,ncol
      do i=1,nrow
        w%m(i,j)=r*w%m(i,j) + (1.d0-s)*f%ptrm%m(i,j) + f%m(i,j)
      end do
    end do
    !$OMP END PARALLEL DO
    
    call w%cyclic_to_column(this%q)
    call this%diags_to_buffer()
    
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n_diags,n,this) &
    !$OMP PRIVATE(i,k)
    do i=1,n_diags
    
      k=n-this%d_midx(i)
      this%tds(i)%rhs=this%d_buff(i,1:k)
      call solve_inv_lap(this%tds(i)%dn,this%tds(i)%en,this%tds(i)%rhs)
      this%d_buff(i,1:k)=this%tds(i)%rhs
    
    end do
    !$OMP END PARALLEL DO
    
    call this%buffer_to_diags()
    call this%q%column_to_cyclic(w)
    call w%make_skewh()

    te = MPI_Wtime() 

    if (IS_MASTER) then
      write(*,'(A,'//output_format_(2:9)//')') '    APPLY H_TURB CPU TIME: ', te-ts
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine apply_inv(this,w)
    class(lap_blk), intent(inout) :: this
    type(cdmatrix), intent(inout) :: w
    complex(double_p), allocatable, dimension(:) :: rhs
    integer :: i,n,s,n_diags
    
    n=this%n
    n_diags=size(this%d_midx)

    call w%cyclic_to_column(this%q)
    call this%diags_to_buffer()
    
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n_diags,n,this) &
    !$OMP PRIVATE(i,s)
    do i=1,n_diags
    
      s=n-this%d_midx(i)
      this%tds(i)%rhs=this%d_buff(i,1:s)
      call solve_inv_lap(this%tds(i)%d,this%tds(i)%e,this%tds(i)%rhs)
      this%d_buff(i,1:s)=this%tds(i)%rhs
    
    end do
    !$OMP END PARALLEL DO
    
    call this%buffer_to_diags()
    call this%q%column_to_cyclic(w)
    call w%make_skewh()  

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine apply(this,w)
    class(lap_blk), intent(inout) :: this
    type(cdmatrix), intent(inout) :: w
    real(double_p), allocatable, dimension(:) :: aux_e
    complex(double_p), allocatable, dimension(:) :: aux_v
    integer :: i,j,n,n_diags,s
    
    n=this%n
    n_diags=size(this%d_midx)

    call w%cyclic_to_column(this%q)
    call this%diags_to_buffer()

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n_diags,n,this) &
    !$OMP PRIVATE(i,j,s,aux_v,aux_e)
    do i=1,n_diags
    
      s=n-this%d_midx(i)
      
      call reallocate_array(aux_v,0,s+1)
      aux_v(0)=(0.d0,0.d0)
      aux_v(s+1)=(0.d0,0.d0)
      aux_v(1:s)=this%d_buff(i,1:s)

      call reallocate_array(aux_e,0,s)
      aux_e(0)=0.d0
      aux_e(1:s)=this%tds(i)%el

      do j=1,s
        this%tds(i)%rhs(j)=aux_v(j-1)*aux_e(j-1) + aux_v(j)*this%tds(i)%dl(j) +&
                           aux_v(j+1)*aux_e(j)
      end do
      
      this%d_buff(i,1:s)=this%tds(i)%rhs
    
    end do
    !$OMP END PARALLEL DO
    
    call this%buffer_to_diags()
    call this%q%column_to_cyclic(w)
    call w%make_skewh()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine update_hturb_forcing(f,q,f0,lf,dlf,dt)
    type(cdmatrix), intent(inout) :: f,q
    type(cdmatrix), intent(in), allocatable, dimension(:) :: f0
    integer, intent(in) :: lf,dlf
    real(double_p), intent(in) :: dt
    integer :: nmat,i,lfi
    
    call f%set_to_zero()
    
    nmat = size(f0)
    
    do i=1,nmat
      
      lfi = lf-dlf+i-1
      
      call update_hturb_forcing_rphase(f,q,f0(i),lfi,dlf,dt)
      
    end do
    
    call f%make_skewh()

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
    call allocate_array(r,1,lf)
    
    if (IS_MASTER) then
      call random_number(r)
      r=r*2.d0*pi
    end if

    call MPI_Bcast(r,lf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
    
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
            aux%re = cos(r(m))
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
    call allocate_array(r1,1,lf)
    call allocate_array(r2,1,lf)
    call allocate_array(u1,1,lf)
    call allocate_array(u2,1,lf)
    
    if (IS_MASTER) then
      call random_number(u1)
      call random_number(u2)
      r1=sqrt(-2.d0*log(u1))*cos(2.d0*pi*u2)*sqrt(dt/3.d0)
      r2=sqrt(-2.d0*log(u1))*sin(2.d0*pi*u2)*sqrt(dt/3.d0)
    end if

    call MPI_Bcast(r1,lf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
    call MPI_Bcast(r2,lf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
    
    call q%copy_values(f0)
    
    gcs=get_gidx(1,q%ncblk,q%pcol,q%npcol)
    ncol=q%ncol
    
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,q,r1,r2,gcs,ncol,lf,dt) &
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
            aux%re = r1(m)
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

!========================================================================================!
  subroutine compute_ic(this,w)
    class(lap_blk), intent(inout) :: this
    type(cdmatrix), intent(inout) :: w
    integer :: m,l_cut
    type(par_file) :: pfile
    complex(double_p) :: alpha,beta
    real(double_p) :: ts,te
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(l_cut,'l_cut')
    
    if (l_cut>this%n-1) then
      call abort_run('l_cut > N-1')
    end if
    
    if (IS_MASTER) then
      write(*,'(A)') ' '
      write(*,'(A)') 'Start computation of i.c.'
    end if
    
    ts = MPI_Wtime()
    
    call this%q%set_to_zero()

    do m=0,this%n-1
    
      if (IS_MASTER) then
        write(*,'(A,I5)') 'Solving eigv problem: ', m
      end if

      !allocate memory for eigenvectors
      !call this%v%ctor(this%mpic,BLOCK_CYCLIC,.FALSE.,this%n-m)
      call this%v%ctor(this%mpic,BLOCK_CYCLIC,this%n-m,1,1)

      !solve eigenvalue problem
      call compute_eigv(this%n,this%v%n,this%v%npcol,this%v%nrow,this%v%ncol,&
                        this%v%is_allocated,this%v%desc,this%v%m)
      
      call this%add_q_component(this%q,l_cut,SPH_IC)

      call this%v%delete()

    end do
    
    call this%q%column_to_cyclic(w)
    call w%make_skewh()
    te = MPI_Wtime()
    
    if (IS_MASTER) then
      write(*,'(A)') 'End computation of i.c.'
      write(*,'(A,'//output_format_(2:9)//')') '    cpu time: ', te-ts
    end if
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine compute_sph_coeff(this,w,out_dir,fields_dir)
    class(lap_blk), intent(inout) :: this
    type(cdmatrix), intent(in) :: w
    integer, intent(in) :: out_dir
    character(len=*), intent(in) :: fields_dir
    type(vector) :: om
    integer :: m,v_size
    
    v_size=sph_size(this%n)
    call om%ctor(this%mpic,'sphc',v_size,w%nprocs)
    
    call w%cyclic_to_column(this%q)
    
    do m=0,this%n-1
      
      !call this%v%ctor(this%mpic,BLOCK_CYCLIC,.FALSE.,this%n-m)
      call this%v%ctor(this%mpic,BLOCK_CYCLIC,this%n-m,1,1)

      call compute_eigv(this%n,this%v%n,this%v%npcol,this%v%nrow,this%v%ncol,&
                        this%v%is_allocated,this%v%desc,this%v%m)
      
      call this%compute_vort_component(this%q,om)
      
      call this%v%delete()  
      
    end do
    
    call om%write_to_disk(out_dir,fields_dir)
    
    call om%delete()
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_hturb_forcing(this,fmag,lf,dlf,f)
    class(lap_blk), intent(inout) :: this
    real(double_p), intent(in) :: fmag
    integer, intent(in) :: lf,dlf
    type(cdmatrix), allocatable, dimension(:), intent(inout) :: f
    integer :: m,lfi,i,nmat
    type(par_file) :: pfile
    complex(double_p) :: alpha,beta
    real(double_p) :: ts,te
    
    !init random generator
    ! -- call random_seed()
    
    nmat = size(f)
    
    do i=1,nmat
      
      call f(i)%ctor(this%mpic,BLOCK_COLUMN,'htf',this%q%n,1,this%q%nprocs)
      
      lfi = lf-dlf+i-1
      
      do m=0,lfi

        !call this%v%ctor(this%mpic,BLOCK_CYCLIC,.FALSE.,this%n-m)
        call this%v%ctor(this%mpic,BLOCK_CYCLIC,this%n-m,1,1)

        call compute_eigv(this%n,this%v%n,this%v%npcol,this%v%nrow,this%v%ncol,&
                          this%v%is_allocated,this%v%desc,this%v%m)
      
        call this%add_q_component(f(i),lfi,SPH_F,dlf)

        call this%v%delete()

      end do
      
      f(i)%m = fmag*f(i)%m
      
    end do
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine add_q_component(this,w,l_ref,qbasis_op,dl_ref)
    class(lap_blk), intent(in) :: this
    type(cdmatrix), intent(inout) :: w
    integer, intent(in) :: l_ref,qbasis_op
    integer, intent(in), optional :: dl_ref
    type(rdmatrix) :: q
    real(double_p), allocatable, dimension(:,:) :: buff
    integer :: gn,n,m,i,ierror
    integer, dimension(4) :: info
    logical :: is_allocated
    
    gn=this%n
    n=this%v%n
    m=gn-n
    
    call q%ctor(this%mpic,BLOCK_COLUMN,.FALSE.,n)

    call this%v%cyclic_to_column(q)
    
    if (q%is_allocated) then
      call adjust_basis_sign(gn,q%ncblk,q%pcol,q%npcol,q%m)      
    end if
    
    do i=0,this%mpic%nprocs-1

      info(1)=q%nrow
      info(2)=q%ncol
      info(3)=q%npcol
      info(4)=q%ncblk
      is_allocated=q%is_allocated
    
      call MPI_BCAST(is_allocated,1,MPI_LOGICAL,i,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(info,4,MPI_INTEGER,i,MPI_COMM_WORLD,ierror)
      call reallocate_array(buff,1,info(1),1,info(2))
      
      if (is_allocated) then
        if (this%mpic%rank==i) then
          buff=q%m
        end if
      else      
        cycle
      end if

      call MPI_BCAST(buff,size(buff),MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierror)
      
      select case(qbasis_op)
        case(SPH_IC)
          call assemble_ic(m,l_ref,info(4),i,info(3),buff,this%dw,&
                           w%is_allocated,w%m)
        case(SPH_F)
          call assemble_hturb_f(m,l_ref,dl_ref,info(4),i,info(3),buff,this%dw,&
                                w%is_allocated,w%m)
        case default
      end select
           
    end do
     
    call q%delete()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine compute_vort_component(this,w,om)
    class(lap_blk), intent(in) :: this
    type(cdmatrix), intent(in) :: w
    type(vector), intent(inout) :: om
    type(rdmatrix) :: q
    real(double_p), allocatable, dimension(:,:) :: buff
    complex(double_p) :: om_val,gom_val
    integer :: gn,n,m,i,j,l,i_loc,i_glo,ierror,gis,gie
    integer, dimension(4) :: info
    logical :: is_allocated
    
    gn=this%n
    n=this%v%n
    m=gn-n
    
    call q%ctor(this%mpic,BLOCK_COLUMN,.FALSE.,n)

    call this%v%cyclic_to_column(q)
    
    if (q%is_allocated) then
      call adjust_basis_sign(gn,q%ncblk,q%pcol,q%npcol,q%m)      
    end if
    
    do i=0,this%mpic%nprocs-1
    
      info(1)=q%nrow
      info(2)=q%ncol
      info(3)=q%npcol
      info(4)=q%ncblk
      is_allocated=q%is_allocated
    
      call MPI_BCAST(is_allocated,1,MPI_LOGICAL,i,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(info,4,MPI_INTEGER,i,MPI_COMM_WORLD,ierror)
      call reallocate_array(buff,1,info(1),1,info(2))

      if (is_allocated) then
        if (this%mpic%rank==i) then
          buff=q%m
        end if
      else      
        cycle
      end if

      call MPI_BCAST(buff,size(buff),MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierror)
      
      gis = get_gidx(1,om%nblk,om%prow,om%nprow)
      gie = get_gidx(om%nrow,om%nblk,om%prow,om%nprow) 
  
      do j=1,size(buff,2)
      
        l = m + (get_gidx(j,info(4),i,info(3))-1)
        om_val=frobenius_prod(l,m,buff(:,j),this%dw,w%is_allocated,w%m)    
        call mpi_allreduce(om_val,gom_val,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
                           MPI_COMM_WORLD,ierror)
                           
        i_glo = sph_idx(l,m)                 
        i_loc = get_lidx(i_glo,om%nblk,om%prow,om%nprow)   
        
        if ( (i_glo>=gis).AND.(i_glo<=gie) ) then
          om%v(i_loc) = gom_val
        end if
        
      end do
           
    end do
     
    call q%delete()    

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine assemble_ic(m,l_cut,ncblk,pcol,npcol,eigv,dw,is_allocated,w)
    integer, intent(in) :: m,l_cut,ncblk,pcol,npcol
    real(double_p), allocatable, dimension(:,:), intent(in) :: eigv
    type(dinfo), allocatable, dimension(:), intent(in) :: dw
    logical, intent(in) :: is_allocated
    complex(double_p), allocatable, dimension(:,:), intent(inout) :: w
    complex(double_p) :: im,wh
    integer :: p,q,i,j,jg,l
    
    if (.not.is_allocated) then
      return
    end if
    
    if (dw(m)%rank==-1) then
      return
    end if
    
    im=(0.d0,1.d0)

    !assign new elements to the basis
    do q=1,size(eigv,2)
    
      l = m + (get_gidx(q,ncblk,pcol,npcol)-1)
      
      if ((l==0).AND.(m==0)) then
        cycle
      end if

      if (l>l_cut) then
        return
      end if
    
      wh=ic_gen(m,l,l_cut)

      do p=1,size(dw(m)%idx,2)
        i=dw(m)%idx(1,p)
        j=dw(m)%idx(2,p)
        jg=dw(m)%idx(3,p)
        w(i,j)=w(i,j)+im*wh*eigv(jg,q)
      end do
    
    end do

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine assemble_hturb_f(m,l_ref,dl_ref,ncblk,pcol,npcol,eigv,dw,is_allocated,w)
    integer, intent(in) :: m,l_ref,dl_ref,ncblk,pcol,npcol
    real(double_p), allocatable, dimension(:,:), intent(in) :: eigv
    type(dinfo), allocatable, dimension(:), intent(in) :: dw
    logical, intent(in) :: is_allocated
    complex(double_p), allocatable, dimension(:,:), intent(inout) :: w
    complex(double_p) :: im,wh
    integer :: p,q,i,j,jg,l
    
    if (.not.is_allocated) then
      return
    end if
    
    if (dw(m)%rank==-1) then
      return
    end if
    
    im=(0.d0,1.d0)

    !assign new elements to the basis
    do q=1,size(eigv,2)
    
      l = m + (get_gidx(q,ncblk,pcol,npcol)-1)
      
      if ((l==0).AND.(m==0)) then
        cycle
      end if
      
      if (l>l_ref) then
        return
      end if

      if (l==l_ref) then
    
        wh=hturbf_gen(m,l,dl_ref,l_ref)

        do p=1,size(dw(m)%idx,2)
          i=dw(m)%idx(1,p)
          j=dw(m)%idx(2,p)
          jg=dw(m)%idx(3,p)
          w(i,j)=w(i,j)+im*wh*eigv(jg,q)
        end do
        
      end if
    
    end do

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_block_systems(this,flow)
    class(lap_blk), intent(inout) :: this
    integer, intent(in) :: flow
    real(double_p), allocatable, dimension(:) :: d_aux,e_aux
    integer :: n,m,i,j,err
    type(par_file) :: pfile
    real(double_p) :: dt,nu,alpha,theta
    
    if (flow==H_TURB) then
      call pfile%ctor('input_parameters','specs')
      call pfile%read_parameter(dt,'dt')
      call pfile%read_parameter(nu,'nu')
      call pfile%read_parameter(alpha,'alpha')
      call pfile%read_parameter(theta,'theta')
      !symmetric strang splitting
      dt=0.5d0*dt
      !account for term \nu \omega
      alpha = alpha - nu
    end if
    
    if (.not.this%q%is_allocated) then
      return
    end if
    
    n=this%n
    
    allocate(this%tds(1:size(this%d_midx)),stat=err)
    if (err /= 0) then
      call abort_run('Allocation of Laplacian block systems failed')
    end if    
    
    do i=1,size(this%d_midx)
    
      m=this%d_midx(i)
      
      !init inverse Laplacian
      call compute_tmat(n,m,this%tds(i)%d,e_aux)
      call reallocate_array(this%tds(i)%e,1,n-m)
      call reallocate_array(this%tds(i)%rhs,1,n-m)
      do j=1,n-m
        this%tds(i)%e(j)%re=e_aux(j)
        this%tds(i)%e(j)%im=0.d0
      end do
      if (m==0) then
        this%tds(i)%d(1)=this%tds(i)%d(1)+0.5d0
      end if
      call store_factorization(this%tds(i)%d,this%tds(i)%e)
      
      if (flow==H_TURB) then
        !init Laplacian
        call compute_tmat(n,m,this%tds(i)%dl,this%tds(i)%el)
        this%tds(i)%dl=-this%tds(i)%dl
        this%tds(i)%el=-this%tds(i)%el

        !init viscosity and damping
        call compute_tmat(n,m,this%tds(i)%dn,e_aux)
        this%tds(i)%dn=-this%tds(i)%dn
        e_aux=-e_aux  
        
        call reallocate_array(this%tds(i)%en,1,n-m)
        do j=1,n-m
          this%tds(i)%en(j)%re=-(1.d0-theta)*dt*nu*e_aux(j)
          this%tds(i)%en(j)%im=0.d0
        end do
        this%tds(i)%dn = -(1.d0-theta)*dt*nu*this%tds(i)%dn +1.d0 +&
                          (1.d0-theta)*dt*alpha
        call store_factorization(this%tds(i)%dn,this%tds(i)%en)
      end if
      
    end do

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine store_factorization(d,e)
    real(double_p), allocatable, dimension(:), intent(inout) :: d
    complex(double_p), allocatable, dimension(:), intent(inout) :: e
    integer :: n,info
    
    n=size(d)
    
    call zpttrf(n,d,e,info)

    if (info.ne.0) then
      call abort_run('zpttrf in lap_blk failed',info)
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine solve_inv_lap(d,e,rhs)
    real(double_p), allocatable, dimension(:), intent(in) :: d
    complex(double_p), allocatable, dimension(:), intent(in) :: e
    complex(double_p), allocatable, dimension(:), intent(inout) :: rhs
    integer :: n,ldb,nrhs,info
    
    n=size(d)
    ldb=n
    nrhs=1
    
    call zpttrs('L',n,nrhs,d,e,rhs,ldb,info)

    if (info.ne.0) then
      call abort_run('zpttrs in lap_blk failed',info)
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  function frobenius_prod(l,m,eigv,dw,is_allocated,w) result(val)
    integer, intent(in) :: l,m
    real(double_p), dimension(:), intent(in) :: eigv
    type(dinfo), allocatable, dimension(:), intent(in) :: dw
    logical, intent(in) :: is_allocated
    complex(double_p), allocatable, dimension(:,:), intent(in) :: w
    complex(double_p) :: im,val
    integer :: p,i,j,jg

    val=(0.d0,0.d0)

    if (.not.is_allocated) then
      return
    end if
    
    if (dw(m)%rank==-1) then
      return
    end if

    if ((l==0).AND.(m==0)) then
      return
    end if
    
    im=(0.d0,1.d0)

    do p=1,size(dw(m)%idx,2)
      i=dw(m)%idx(1,p)
      j=dw(m)%idx(2,p)
      jg=dw(m)%idx(3,p)
      val=val+w(i,j)*(-im)*eigv(jg)
    end do

  end function
!========================================================================================!

!========================================================================================!
  subroutine compute_tmat(n,m,d,e)
    integer, intent(in) :: n,m
    real(double_p), allocatable, dimension(:), intent(out) :: d,e
    integer, allocatable, dimension(:) :: j
    real(double_p) :: s
    integer :: i,p
    
    s=real(n-1)/real(2)
    
    call allocate_array(j,1,n-m)
    call allocate_array(d,1,size(j))
    call allocate_array(e,1,size(j))
    d=0.d0
    e=0.d0
      
    !fill up diagonal
    j=(/(i,i=0,n-m-1)/)
    do p=1,size(j)
      d(p)=2.d0*(s*(2.d0*j(p)+1.d0+m)-j(p)*(j(p)+m))
    end do
      
    !fill up off-diagonal
    j=(/(i,i=1,n-m)/)
    do p=1,size(j)-1
      e(p)=-sqrt(dble((j(p)+m)*(n-j(p)-m)))*sqrt(dble(j(p)*(n-j(p))))
    end do
 
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine compute_eigv(gn,n,npcol,nrow,ncol,is_allocated,desc,v)
    integer, intent(in) :: gn,n,npcol,nrow,ncol
    logical, intent(in) :: is_allocated
    integer, dimension(9), intent(in) :: desc
    real(double_p), allocatable, dimension(:,:), intent(inout) :: v
    real(double_p), allocatable, dimension(:) :: work,d,e
    integer, allocatable, dimension(:) :: iwork,recvcount,displs
    integer :: lwork,liwork,iq,jq,info
    
    if (.not.is_allocated) then
      return
    end if

    call compute_tmat(gn,gn-n,d,e)

    liwork = 2 + 7*gn + 8*npcol
    call allocate_array(iwork,1,liwork)
    iq=1
    jq=1
    
    call allocate_array(work,1,1)
    call PDSTEDC('I',n,d,e,v,iq,jq,desc,work,-1,iwork,liwork,info)

    if (info.ne.0) then
      call abort_run('subroutine PDSTEDC (LWORK=-1) in lap_blk failed',info)
    end if
    
    lwork=int(work(1))
    call reallocate_array(work,1,lwork)
    
    call set_num_threads()
    
    call PDSTEDC('I',n,d,e,v,iq,jq,desc,work,lwork,iwork,liwork,info)
    
    call set_sigle_thread()

    if (info.ne.0) then
      call abort_run('subroutine PDSTEDC in lap_blk failed',info)
    end if

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine adjust_basis_sign(n,ncblk,pcol,npcol,v)
    integer, intent(in) :: n,ncblk,pcol,npcol
    real(double_p), allocatable, dimension(:,:) :: v
    real(double_p) :: val
    integer :: j,jg,jh,nrow
    logical :: is_n_odd
    
    if (modulo(n,2)==1) then
      is_n_odd = .TRUE.
    else
      is_n_odd = .FALSE.
    end if
    
    nrow = size(v,1)
    
    do j=1,size(v,2)
    
      jg = get_gidx(j,ncblk,pcol,npcol)
      
      if (modulo(nrow,2)==0) then
      
        val = v(nrow/2+1,j)
        if (val==0.d0) then
          call abort_run('could not determine sign quantised basis')
        end if
        
        if (modulo(jg,2)==0) then
          jh = jg/2
        else
          jh = jg/2 + 1
        end if
        
        if (modulo(jh,2)==0) then
          if (is_n_odd) then
            if (val<0.d0) then
              v(:,j)=-v(:,j)
            end if
          else
            if (val>0.d0) then
              v(:,j)=-v(:,j)
            end if            
          end if
        else
          if (is_n_odd) then
            if (val>0.d0) then
              v(:,j)=-v(:,j)
            end if
          else
            if (val<0.d0) then
              v(:,j)=-v(:,j)
            end if
          end if       
        end if
        
      else
    
        if (modulo(jg,2)==0) then
        
          val = v(nrow/2+2,j)
          if (val==0.d0) then
            call abort_run('could not determine sign quantised basis')
          end if
          jh = jg/2
          if (modulo(jh,2)==0) then
            if (is_n_odd) then
              if (val>0.d0) then
                v(:,j)=-v(:,j)
              end if
            else
              if (val<0.d0) then
                v(:,j)=-v(:,j)
              end if
            end if          
          else
            if (is_n_odd) then
              if (val<0.d0) then
                v(:,j)=-v(:,j)
              end if
            else
              if (val>0.d0) then
                v(:,j)=-v(:,j)
              end if
            end if
          end if
        
        else
      
          val = v(nrow/2+1,j)
          if (val==0.d0) then
            call abort_run('could not determine sign quantised basis')
          end if
          jh = jg/2 + 1
          if (modulo(jh,2)==0) then
            if (is_n_odd) then
              if (val>0.d0) then
                v(:,j)=-v(:,j)
              end if
            else
              if (val<0.d0) then
                v(:,j)=-v(:,j)
              end if
            end if        
          else
            if (is_n_odd) then
              if (val<0.d0) then
                v(:,j)=-v(:,j)
              end if
            else
              if (val>0.d0) then
                v(:,j)=-v(:,j)
              end if
            end if
          end if        
        
        end if

      end if

    end do

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_diags(this)
    class(lap_blk), intent(inout) :: this
    integer :: n,np,ncol,n_diags,s,i,j,p,l,nblocks,gcs,count,rank,ierror,res,displ_0
    integer, allocatable, dimension(:) :: block_lenghts,displ,gncol
    
    if (.not.this%q%is_allocated) then
      return
    end if
    
    n=this%n
    np=this%q%nprocs
    rank=this%q%rank
    ncol=this%q%ncol
    call allocate_array(gncol,0,np-1)
    call mpi_allgather(ncol,1,MPI_INTEGER,gncol,1,MPI_INTEGER,this%q%comm,ierror)
    
    !init logical send/recv
    call allocate_array(this%sr_diag,1,2,0,np-1)
    gcs=get_gidx(1,this%q%ncblk,this%q%pcol,this%q%npcol)
    this%sr_diag=.FALSE.
    !send
    do i=0,np-1
      
      l=gcs+i
      if (l<=n) then
        this%sr_diag(1,i)=.TRUE.
      end if
      
    end do

    !recv
    do i=0,np-1
      
      gcs=get_gidx(1,this%q%ncblk,i,this%q%npcol)
    
      do j=0,np-1
        
        if (gcs+j<=n) then
        
          if (j==rank) then
            this%sr_diag(2,i)=.TRUE.
            exit
          end if
          
        end if
      
      end do
    
    end do
    
    !init mpi derived types
    !dmat_type
    gcs=get_gidx(1,this%q%ncblk,this%q%pcol,this%q%npcol)
    call allocate_array(this%dmat_type,0,np-1)
    do p=0,np-1
    
      if (this%sr_diag(1,p)) then
    
        l=gcs+p
      
        !set block size
        count=0
        do
          nblocks = min(ncol,n-l+1)
          if (nblocks>=1) then
            count=count+nblocks
            l=l+np
          else
            exit
          end if
        end do
      
        !set block lengths and displ
        call reallocate_array(block_lenghts,1,count)
        call reallocate_array(displ,1,count)
        block_lenghts=1
        displ(1)=0
        
        count=1
        do i=1,ncol
        
          l=gcs+p+(i-1)
          
          do
            l=l+np
            if (l<=n) then
              if (count<size(displ)) then
                count=count+1
                displ(count)=displ(count-1)+np
              end if
            else
              if (count<size(displ)) then
                count=count+1
                res=n-(l-np)
                displ(count)=displ(count-1)+res+gcs+p+i
              end if
              exit
            end if
          end do
          
        end do

        call mpi_type_indexed(size(block_lenghts),block_lenghts,displ,&
                              MPI_DOUBLE_COMPLEX,this%dmat_type(p),ierror)
        call mpi_type_commit(this%dmat_type(p),ierror)
        
      end if

    end do
    
    !dbuff_type
    call allocate_array(this%dbuff_type,0,np-1)
    n_diags=1
    count=rank+1
    do
      count=count+np
      if (count<=n) then
        n_diags=n_diags+1
      else
        exit
      end if
    end do
    call allocate_array(this%d_midx,1,n_diags)
    do i=1,n_diags
      s=n-rank-(i-1)*np
      this%d_midx(i)=n-s
    end do
    call allocate_array(this%d_buff,1,n_diags,1,maxval(n-this%d_midx(:)))
    this%d_buff=(-1.d0,-1.d0)
    
    do p=0,np-1
    
      if (this%sr_diag(2,p)) then
        
        gcs=get_gidx(1,this%q%ncblk,p,this%q%npcol)
        
        count=0
        do i=1,gncol(p)
          
          l=gcs+rank+(i-1)
          if (l<=n) then
            count=count+1
          end if

        end do
        
        call reallocate_array(block_lenghts,1,count)
        call reallocate_array(displ,1,count)
        
        do i=1,size(block_lenghts)
          
          count=0
          l=gcs+rank+(i-1)
          
          do
            if (l<=n) then
              count=count+1
            else
              block_lenghts(i)=count
              if (p>0) then
                displ(i)=(i-1)*n_diags + n_diags*sum(gncol(0:p-1))
              else
                displ(i)=(i-1)*n_diags
              end if
              exit
            end if
            l=l+np
          end do

        end do
        
        call mpi_type_indexed(size(block_lenghts),block_lenghts,displ,&
                              MPI_DOUBLE_COMPLEX,this%dbuff_type(p),ierror)
        call mpi_type_commit(this%dbuff_type(p),ierror)
        
      end if
      
    end do

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine diags_to_buffer(this)
    class(lap_blk), intent(inout) :: this
    integer :: n,np,tag,gcs,p,ireq,nreq,ierror,is
    integer, allocatable, dimension(:,:) :: status
    integer, allocatable, dimension(:) :: requests

    if (.not.this%q%is_allocated) then
      return
    end if
     
    n=this%n
    np=this%q%nprocs
    tag=0
    
    nreq=0
    do p=0,np-1
      if (this%sr_diag(1,p)) then
        nreq=nreq+1
      end if
      if (this%sr_diag(2,p)) then
        nreq=nreq+1
      end if
    end do
    call allocate_array(status,1,MPI_STATUS_SIZE,1,nreq)
    call allocate_array(requests,1,nreq)
    
    ireq=0
     
    !recv diags from matrix
    do p=0,np-1
        
      if (this%sr_diag(2,p)) then
      
        ireq=ireq+1
        call MPI_IRECV(this%d_buff,1,this%dbuff_type(p),p,tag,&
                       this%q%comm,requests(ireq),ierror)
      end if
      
    end do
    
    !send diags do buffer
    gcs=get_gidx(1,this%q%ncblk,this%q%pcol,this%q%npcol)
    do p=0,np-1
    
      if (this%sr_diag(1,p)) then
        
        is=gcs+p
        
        ireq=ireq+1
        call MPI_ISSEND(this%q%m(is,1),1,this%dmat_type(p),p,tag,&
                        this%q%comm,requests(ireq),ierror)
      end if
    
    end do
    
    call MPI_WAITALL(nreq,requests,status,ierror)

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine buffer_to_diags(this)
    class(lap_blk), intent(inout) :: this
    integer :: n,np,tag,gcs,p,ireq,nreq,ierror,is
    integer, allocatable, dimension(:,:) :: status
    integer, allocatable, dimension(:) :: requests

    if (.not.this%q%is_allocated) then
      return
    end if
     
    n=this%n
    np=this%q%nprocs
    tag=0
    
    nreq=0
    do p=0,np-1
      if (this%sr_diag(1,p)) then
        nreq=nreq+1
      end if
      if (this%sr_diag(2,p)) then
        nreq=nreq+1
      end if
    end do
    call allocate_array(status,1,MPI_STATUS_SIZE,1,nreq)
    call allocate_array(requests,1,nreq)
    
    ireq=0
     
    !recv diags from buffer
    gcs=get_gidx(1,this%q%ncblk,this%q%pcol,this%q%npcol)
    do p=0,np-1
        
      if (this%sr_diag(1,p)) then
      
        is=gcs+p
      
        ireq=ireq+1
        call MPI_IRECV(this%q%m(is,1),1,this%dmat_type(p),p,tag,&
                       this%q%comm,requests(ireq),ierror)
      end if
      
    end do
    
    !send diags to matrix
    do p=0,np-1
    
      if (this%sr_diag(2,p)) then
        
        ireq=ireq+1
        call MPI_ISSEND(this%d_buff,1,this%dbuff_type(p),p,tag,&
                        this%q%comm,requests(ireq),ierror)
      end if
    
    end do
    
    call MPI_WAITALL(nreq,requests,status,ierror)

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_dinfo(this,w)
    class(lap_blk), intent(inout) :: this
    type(cdmatrix), intent(in) :: w
    integer :: n,npcol,pcol,ncblk,ncol,m,i,j,gc_s,gc_e,gis,gie,d_size,err
    
    if (.not.w%is_allocated) then
      return
    end if
    
    n=w%n
    npcol=w%npcol
    pcol=w%pcol
    ncblk=w%ncblk
    ncol=w%ncol
    
    allocate(this%dw(0:n-1),stat=err)
    if (err /= 0) then
      call abort_run('Allocation of dinfo in lap_blk failed ')
    end if 

    !store diagonal info
    gc_s=get_gidx(1,ncblk,pcol,npcol)
    gc_e=get_gidx(ncol,ncblk,pcol,npcol)
    
    do i=1,size(w%m,1)
      
      m = i-1
      gis = gc_s + m
      
      if (gis > n) then
        this%dw(m)%rank = -1
        cycle
      end if
      
      this%dw(m)%rank = w%rank
      
      gie = min(gis+ncol-1,n)
      d_size = gie - gis +1
      call reallocate_array(this%dw(m)%idx,1,4,1,d_size)
      
      do j=1,d_size
        this%dw(m)%idx(1,j)=gis+(j-1)
        this%dw(m)%idx(2,j)=j
        this%dw(m)%idx(3,j)=this%dw(m)%idx(1,j)-m
      end do
      
    end do

  end subroutine
!========================================================================================!

end module lap_blk_mod