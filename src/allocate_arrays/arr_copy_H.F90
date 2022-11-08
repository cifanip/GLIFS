if (.not.allocated(v)) then
  call abort_run('Attempt to copy from deallocated array')
end if

#if vec_dim==1
if (allocated(cpy)) then
  if ((lbound(v,1).ne.lbound(cpy,1)).OR.&
     (ubound(v,1).ne.ubound(cpy,1))) then
     call abort_run('Attempt to copy arrays of different size')
  end if 
else
  call allocate_array(cpy,lbound(v,1),ubound(v,1))
end if
cpy=v
#endif

#if vec_dim==2
if (allocated(cpy)) then
  if ((lbound(v,1).ne.lbound(cpy,1)).OR.&
     (ubound(v,1).ne.ubound(cpy,1)).OR.&
     (lbound(v,2).ne.lbound(cpy,2)).OR.&
     (ubound(v,2).ne.ubound(cpy,2))) then
     call abort_run('Attempt to copy arrays of different size')
  end if 
else
  call allocate_array(cpy,lbound(v,1),ubound(v,1),&
                         lbound(v,2),ubound(v,2))
end if
cpy=v
#endif

#if vec_dim==3
if (allocated(cpy)) then
  if ((lbound(v,1).ne.lbound(cpy,1)).OR.&
     (ubound(v,1).ne.ubound(cpy,1)).OR.&
     (lbound(v,2).ne.lbound(cpy,2)).OR.&
     (ubound(v,2).ne.ubound(cpy,2)).OR.&
     (lbound(v,3).ne.lbound(cpy,3)).OR.&
     (ubound(v,3).ne.ubound(cpy,3))) then
     call abort_run('Attempt to copy arrays of different size')
  end if 
else
  call allocate_array(cpy,lbound(v,1),ubound(v,1),&
                         lbound(v,2),ubound(v,2),&
                         lbound(v,3),ubound(v,3))
end if
cpy=v
#endif