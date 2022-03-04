module formats_mod
  
  implicit none
  
  integer, parameter :: IO_unit_ = 50
  character(len=4), parameter  :: int_format_     = '(I6)'
  character(len=3), parameter  :: char_format_    = '(A)'
  character(len=4), parameter  :: logical_format_ = '(L7)'
  character(len=11), parameter :: double_format_  = '(ES25.15E3)'
  character(len=10), parameter :: output_format_  = '(ES11.4E2)'
  
end module formats_mod