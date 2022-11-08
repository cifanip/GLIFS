#if vec_dim==1
if (allocated(v)) then
  call deallocate_array(v)
  call allocate_array(v,st,en) 
else
  call allocate_array(v,st,en)
end if
#endif

#if vec_dim==2
if (allocated(v)) then
  call deallocate_array(v)
  call allocate_array(v,st1,en1,st2,en2) 
else
  call allocate_array(v,st1,en1,st2,en2)
end if
#endif

#if vec_dim==3
if (allocated(v)) then
  call deallocate_array(v)
  call allocate_array(v,st1,en1,st2,en2,st3,en3) 
else
  call allocate_array(v,st1,en1,st2,en2,st3,en3)
end if
#endif