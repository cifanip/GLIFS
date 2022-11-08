#if vec_dim==1
if (.not. allocated(v)) then
  allocate(v(st:en),STAT=err)
  if (err /= 0) then
    call abort_run('Allocation of v in allocateArray1D failed ') 
  end if  
else
  call abort_run('Attempt to allocate allocated array ') 
end if
#endif

#if vec_dim==2
if (.not. allocated(v)) then
  allocate(v(st1:en1,st2:en2),STAT=err)
  if (err /= 0) then
    call abort_run('Allocation of v in allocateArray2D failed ')
  end if    
else
  call abort_run('Attempt to allocate allocated array ') 
end if
#endif

#if vec_dim==3
if (.not. allocated(v)) then
  allocate(v(st1:en1,st2:en2,st3:en3),STAT=err)
  if (err /= 0) then
    call abort_run('Allocation of v in allocateArray3D failed ')
  end if
else
  call abort_run('Attempt to allocate allocated array ') 
end if
#endif

#if vec_dim==4
if (.not. allocated(v)) then
  allocate(v(st1:en1,st2:en2,st3:en3,st4:en4),STAT=err)
  if (err /= 0) then
    call abort_run('Allocation of v in allocateArray4D failed ')
  end if
else
  call abort_run('Attempt to allocate allocated array ') 
end if
#endif