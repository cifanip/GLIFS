module kinds_mod
  implicit none
  
  !kinds
  integer, parameter :: single_p = KIND(0.e0)
  integer, parameter :: double_p = KIND(0.d0)
  
  !size
  integer, parameter :: real_dp_size = SIZEOF(0.d0)
  integer, parameter :: complex_dp_size = 2*real_dp_size
  integer, parameter :: integer_size = SIZEOF(0)
  integer, parameter :: logical_size = SIZEOF(.TRUE.)

end module kinds_mod