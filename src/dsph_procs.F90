module dsph_procs_mod

  use laplacian_mod
  use ic_mod
  
  implicit none

  integer, parameter :: IC_TEST   = 0
  integer, parameter :: IC_ZERO   = 1
  integer, parameter :: IC_RANDOM = 2
  
  public :: compute_ic,&
            init_hturb_forcing,&
            init_coriolis_matrix,&
            compute_sph_coeff
             
  private :: add_q_component,&
             assemble_ic,&
             frobenius_prod,&
             compute_vort_component,&
             compute_eigv,&
             adjust_basis_sign

contains

!========================================================================================!
  subroutine compute_ic(lap,w)
    type(laplacian), intent(inout) :: lap
    type(cdmatrix), intent(inout) :: w
    type(ic) :: ic_field
    integer :: m,l_max,ic_type
    real(double_p) :: ts,te
    type(par_file) :: pfile
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(ic_type,'initial_conditions')
    
    if (IS_MASTER) then
      write(*,'(A)') ' '
      write(*,'(A)') 'Start computation of i.c.'
    end if
    
    select case(ic_type)
      case(IC_TEST)
        call ic_field%ic_gen_test(lap%mpic,w%n)
      case(IC_ZERO)
        w%m = (0.d0,0.d0)
        if (IS_MASTER) then
          write(*,'(A)') 'End computation of i.c.'
        end if
        return
      case(IC_RANDOM)
        call ic_field%ic_gen_rnd(lap%mpic,w%n)
      case default
        call abort_run('Wrong i.c. type found')
    end select
    
    ts = MPI_Wtime()
    
    l_max = ic_field%l_max
    
    call lap%q%set_to_zero()

    do m=0,l_max
    
      if (IS_MASTER) then
        write(*,'(A,I5)') 'Solving eigv problem: ', m
      end if

      !allocate memory for eigenvectors
      !call lap%v%ctor(lap%mpic,BLOCK_CYCLIC,.FALSE.,lap%n-m)
      call lap%v%ctor(lap%mpic,BLOCK_CYCLIC,lap%n-m,1,1)

      !solve eigenvalue problem
      call compute_eigv(lap%n,lap%v%n,lap%v%npcol,lap%v%nrow,lap%v%ncol,&
                        lap%v%is_allocated,lap%v%desc,lap%v%m)
      
      call add_q_component(lap,lap%q,ic_field)

      call lap%v%delete()

    end do
    
    call lap%q%column_to_cyclic(w)
    call w%make_skewh()
    te = MPI_Wtime()
    
    call ic_field%delete()
    
    if (IS_MASTER) then
      write(*,'(A)') 'End computation of i.c.'
      write(*,'(A,'//output_format_(2:9)//')') '    cpu time: ', te-ts
    end if
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_hturb_forcing(lap,w,lf,dlf,f,f0)
    type(laplacian), intent(inout) :: lap
    type(cdmatrix), intent(in) :: w
    integer, intent(out) :: lf,dlf
    type(cdmatrix), intent(out) :: f
    type(cdmatrix), allocatable, dimension(:), intent(out) :: f0
    type(ic), allocatable, dimension(:) :: ic_f
    integer :: m,lfi,i,nmat,err
    real(double_p) :: ts,te
    type(par_file) :: pfile
    
    f = w
    call f%rename('f')
    
    !set spectral space box width to 2
    dlf = 2
    
    allocate(f0(2*dlf+1),stat=err)
    if (err /= 0) then
      call abort_run('Allocation of f0 in init_hturb_forcing failed ')
    end if 
    allocate(ic_f(2*dlf+1),stat=err)
    if (err /= 0) then
      call abort_run('Allocation of ic_f in init_hturb_forcing failed ')
    end if 
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(lf,'lf')
    
    if (lf>lap%n-1) then
      call abort_run('lf > N-1')
    end if
    if (lf<1) then
      call abort_run('lf < 1')
    end if    
    if (lf - dlf<1) then
      call abort_run('lf - dlf < 1')
    end if
    if (lf + dlf>lap%n-1) then
      call abort_run('lf + dlf > n-1')
    end if
    
    nmat = size(f0)
    
    do i=1,nmat
      
      call f0(i)%ctor(lap%mpic,BLOCK_COLUMN,'htf',lap%q%n,1,lap%q%nprocs)
      
      lfi = lf-dlf+i-1
      call ic_f(i)%ic_gen_hturbf(lap%mpic,lap%n,lf,lfi,dlf)

      do m=0,lfi

        !call lap%v%ctor(lap%mpic,BLOCK_CYCLIC,.FALSE.,lap%n-m)
        call lap%v%ctor(lap%mpic,BLOCK_CYCLIC,lap%n-m,1,1)

        call compute_eigv(lap%n,lap%v%n,lap%v%npcol,lap%v%nrow,lap%v%ncol,&
                          lap%v%is_allocated,lap%v%desc,lap%v%m)
      
        call add_q_component(lap,f0(i),ic_f(i))

        call lap%v%delete()

      end do
      
      call ic_f(i)%delete()
      
    end do
    
    deallocate(ic_f)
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine init_coriolis_matrix(lap,w,f)
    type(laplacian), intent(inout) :: lap
    type(cdmatrix), intent(in) :: w
    type(cdmatrix), intent(inout) :: f
    real(double_p) :: cor_par
    type(ic) :: ic_cor
    integer :: l,m
    type(par_file) :: pfile
    
    f = w
    
    call pfile%ctor('input_parameters','specs')
    call pfile%read_parameter(cor_par,'cor_par')
    
    l=1
    m=0
    
    call ic_cor%ic_gen_Tbase(lap%mpic,lap%n,l,m)
    
    call lap%q%set_to_zero()
    call lap%v%ctor(lap%mpic,BLOCK_CYCLIC,lap%n-m,1,1)
    call compute_eigv(lap%n,lap%v%n,lap%v%npcol,lap%v%nrow,lap%v%ncol,&
                      lap%v%is_allocated,lap%v%desc,lap%v%m)
    call add_q_component(lap,lap%q,ic_cor)
    call lap%v%delete()
    
    !iT_{1,0}
    call lap%q%column_to_cyclic(f)
    
    f%m = cor_par*f%m
    
    call ic_cor%delete()
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine add_q_component(lap,w,ic_field)
    type(laplacian), intent(in) :: lap
    type(cdmatrix), intent(inout) :: w
    type(ic), intent(in) :: ic_field
    type(rdmatrix) :: q
    real(double_p), allocatable, dimension(:,:) :: buff
    integer :: gn,n,m,i,l_ref_min,l_ref_max,ierror
    integer, dimension(4) :: info
    logical :: is_allocated
    
    gn=lap%n
    n=lap%v%n
    m=gn-n
    
    call q%ctor(lap%mpic,BLOCK_COLUMN,.FALSE.,n)

    call lap%v%cyclic_to_column(q)
    
    if (q%is_allocated) then
      call adjust_basis_sign(gn,q%ncblk,q%pcol,q%npcol,q%m)      
    end if
    
    do i=0,lap%mpic%nprocs-1

      info(1)=q%nrow
      info(2)=q%ncol
      info(3)=q%npcol
      info(4)=q%ncblk
      is_allocated=q%is_allocated
    
      call MPI_BCAST(is_allocated,1,MPI_LOGICAL,i,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(info,4,MPI_INTEGER,i,MPI_COMM_WORLD,ierror)
      call reallocate_array(buff,1,info(1),1,info(2))
      
      if (is_allocated) then
        if (lap%mpic%rank==i) then
          buff=q%m
        end if
      else      
        cycle
      end if

      call MPI_BCAST(buff,size(buff),MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierror)

      l_ref_min = ic_field%l_min
      l_ref_max = ic_field%l_max

      call assemble_ic(m,l_ref_min,l_ref_max,info(4),i,info(3),buff,&
                       lap%dw,w%is_allocated,ic_field%ic_vec,w%m)
           
    end do
     
    call q%delete()

  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine assemble_ic(m,l_min,l_max,ncblk,pcol,npcol,eigv,dw,is_allocated,ic_v,w)
    integer, intent(in) :: m,l_min,l_max,ncblk,pcol,npcol
    real(double_p), allocatable, dimension(:,:), intent(in) :: eigv
    type(dinfo), allocatable, dimension(:), intent(in) :: dw
    logical, intent(in) :: is_allocated
    type(vector), intent(in) :: ic_v
    complex(double_p), allocatable, dimension(:,:), intent(inout) :: w
    complex(double_p) :: im,wh
    integer :: p,q,i,j,jg,l,gr,lr
    
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

      if ((l>=l_min).AND.(l<=l_max)) then
      
        !global sph index
        gr = sph_idx(l,m)
        
        !i.c. value
        wh = ic_v%vg(gr)

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
  subroutine compute_sph_coeff(lap,w,out_dir,fields_dir)
    type(laplacian), intent(inout) :: lap
    type(cdmatrix), intent(in) :: w
    integer, intent(in) :: out_dir
    character(len=*), intent(in) :: fields_dir
    type(vector) :: om
    integer :: m,v_size
    
    v_size=sph_size(lap%n)
    call om%ctor(lap%mpic,trim(adjustl(w%fname))//'_sphc',v_size,w%nprocs)
    
    call w%cyclic_to_column(lap%q)
    
    do m=0,lap%n-1
      
      !call lap%v%ctor(lap%mpic,BLOCK_CYCLIC,.FALSE.,lap%n-m)
      call lap%v%ctor(lap%mpic,BLOCK_CYCLIC,lap%n-m,1,1)

      call compute_eigv(lap%n,lap%v%n,lap%v%npcol,lap%v%nrow,lap%v%ncol,&
                        lap%v%is_allocated,lap%v%desc,lap%v%m)
      
      call compute_vort_component(lap,lap%q,om)
      
      call lap%v%delete()  
      
    end do
    
    call om%write_to_disk(out_dir,fields_dir)
    
    call om%delete()
    
  end subroutine
!========================================================================================!

!========================================================================================!
  subroutine compute_vort_component(lap,w,om)
    type(laplacian), intent(in) :: lap
    type(cdmatrix), intent(in) :: w
    type(vector), intent(inout) :: om
    type(rdmatrix) :: q
    real(double_p), allocatable, dimension(:,:) :: buff
    complex(double_p) :: om_val,gom_val
    integer :: gn,n,m,i,j,l,i_loc,i_glo,ierror,gis,gie
    integer, dimension(4) :: info
    logical :: is_allocated
    
    gn=lap%n
    n=lap%v%n
    m=gn-n
    
    call q%ctor(lap%mpic,BLOCK_COLUMN,.FALSE.,n)

    call lap%v%cyclic_to_column(q)
    
    if (q%is_allocated) then
      call adjust_basis_sign(gn,q%ncblk,q%pcol,q%npcol,q%m)      
    end if
    
    do i=0,lap%mpic%nprocs-1
    
      info(1)=q%nrow
      info(2)=q%ncol
      info(3)=q%npcol
      info(4)=q%ncblk
      is_allocated=q%is_allocated
    
      call MPI_BCAST(is_allocated,1,MPI_LOGICAL,i,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(info,4,MPI_INTEGER,i,MPI_COMM_WORLD,ierror)
      call reallocate_array(buff,1,info(1),1,info(2))

      if (is_allocated) then
        if (lap%mpic%rank==i) then
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
        om_val=frobenius_prod(l,m,buff(:,j),lap%dw,w%is_allocated,w%m)    
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

    call compute_lap_blk(gn,gn-n,d,e)

    liwork = 2 + 7*gn + 8*npcol
    call allocate_array(iwork,1,liwork)
    iq=1
    jq=1
    
    call allocate_array(work,1,1)
    call PDSTEDC('I',n,d,e,v,iq,jq,desc,work,-1,iwork,liwork,info)

    if (info.ne.0) then
      call abort_run('subroutine PDSTEDC (LWORK=-1) in laplacian failed',info)
    end if
    
    lwork=int(work(1))
    call reallocate_array(work,1,lwork)
    
    call set_blas_multiple_threads()
    
    call PDSTEDC('I',n,d,e,v,iq,jq,desc,work,lwork,iwork,liwork,info)
    
    call set_blas_sigle_thread()

    if (info.ne.0) then
      call abort_run('subroutine PDSTEDC in laplacian failed',info)
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
  
end module dsph_procs_mod
