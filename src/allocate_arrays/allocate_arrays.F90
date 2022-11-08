module allocate_arrays_mod

  use kinds_mod
  use prun_mod

  implicit none
  
  interface allocate_array 
    module PROCEDURE allocate_array_1d_complex
    module PROCEDURE allocate_array_2d_complex
    module PROCEDURE allocate_array_3d_complex

    module PROCEDURE allocate_array_1d_double
    module PROCEDURE allocate_array_2d_double
    module PROCEDURE allocate_array_3d_double
    module PROCEDURE allocate_array_4d_double
      
    module PROCEDURE allocate_array_1d_integer
    module PROCEDURE allocate_array_2d_integer
    module PROCEDURE allocate_array_3d_integer
      
    module PROCEDURE allocate_array_1d_logical
    module PROCEDURE allocate_array_2d_logical
    module PROCEDURE allocate_array_3d_logical
  end interface
    
  interface deallocate_array
    module PROCEDURE deallocate_array_1d_complex
    module PROCEDURE deallocate_array_2d_complex
    module PROCEDURE deallocate_array_3d_complex
    
    module PROCEDURE deallocate_array_1d_double
    module PROCEDURE deallocate_array_2d_double
    module PROCEDURE deallocate_array_3d_double
      
    module PROCEDURE deallocate_array_1d_integer
    module PROCEDURE deallocate_array_2d_integer
    module PROCEDURE deallocate_array_3d_integer
      
    module PROCEDURE deallocate_array_1d_logical
    module PROCEDURE deallocate_array_2d_logical
    module PROCEDURE deallocate_array_3d_logical
  end interface
    
  interface reallocate_array
    module PROCEDURE reallocate_array_1d_complex
    module PROCEDURE reallocate_array_2d_complex
    module PROCEDURE reallocate_array_3d_complex
    
    module PROCEDURE reallocate_array_1d_double
    module PROCEDURE reallocate_array_2d_double
    module PROCEDURE reallocate_array_3d_double
      
    module PROCEDURE reallocate_array_1d_integer
    module PROCEDURE reallocate_array_2d_integer
    module PROCEDURE reallocate_array_3d_integer
      
    module PROCEDURE reallocate_array_1d_logical
    module PROCEDURE reallocate_array_2d_logical
    module PROCEDURE reallocate_array_3d_logical
  end interface

  interface copy_array
    module PROCEDURE copyArray_1d_complex
    module PROCEDURE copyArray_2d_complex
    module PROCEDURE copyArray_3d_complex
    
    module PROCEDURE copyArray_1d_double
    module PROCEDURE copyArray_2d_double
    module PROCEDURE copyArray_3d_double
      
    module PROCEDURE copyArray_1d_integer
    module PROCEDURE copyArray_2d_integer
    module PROCEDURE copyArray_3d_integer
      
    module PROCEDURE copyArray_1d_logical
    module PROCEDURE copyArray_2d_logical
    module PROCEDURE copyArray_3d_logical
  end interface
  
contains

!                     allocate 
!========================================================================================!
  subroutine allocate_array_1d_complex(v,st,en) 
    integer, intent(IN) :: st, en
    complex(double_p), allocatable, dimension(:), intent(INOUT) :: v
    integer :: err
#   define vec_dim 1
#   include "arr_alloc_H.F90"
#   undef vec_dim 
  end subroutine

  subroutine allocate_array_2d_complex(v,st1,en1,st2,en2) 
    integer, intent(IN) :: st1, en1, st2, en2
    complex(double_p), allocatable, dimension(:,:), intent(INOUT) :: v
    integer :: err
#   define vec_dim 2
#   include "arr_alloc_H.F90"
#   undef vec_dim
  end subroutine
  
  subroutine allocate_array_3d_complex(v,st1,en1,st2,en2,st3,en3) 
    integer, intent(IN) :: st1, en1, st2, en2, st3, en3
    complex(double_p), allocatable, dimension(:,:,:), intent(INOUT) :: v
    integer :: err
#   define vec_dim 3
#   include "arr_alloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine allocate_array_1d_double(v,st,en) 
    integer, intent(IN) :: st, en
    real(double_p), allocatable, dimension(:), intent(INOUT) :: v
    integer :: err
#   define vec_dim 1
#   include "arr_alloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine allocate_array_2d_double(v,st1,en1,st2,en2) 
    integer, intent(IN) :: st1, en1, st2, en2
    real(double_p), allocatable, dimension(:,:), intent(INOUT) :: v
    integer :: err
#   define vec_dim 2
#   include "arr_alloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine allocate_array_3d_double(v,st1,en1,st2,en2,st3,en3) 
    integer, intent(IN) :: st1, en1, st2, en2, st3, en3
    real(double_p), allocatable, dimension(:,:,:), intent(INOUT) :: v
    integer :: err
#   define vec_dim 3
#   include "arr_alloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine allocate_array_4d_double(v,st1,en1,st2,en2,st3,en3,st4,en4) 
    integer, intent(IN) :: st1, en1, st2, en2, st3, en3, st4, en4
    real(double_p), allocatable, dimension(:,:,:,:), intent(INOUT) :: v
    integer :: err
#   define vec_dim 4
#   include "arr_alloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine allocate_array_1d_integer(v,st,en) 
    integer, intent(IN) :: st, en
    integer, allocatable, dimension(:), intent(INOUT) :: v
    integer :: err
#   define vec_dim 1
#   include "arr_alloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine allocate_array_2d_integer(v,st1,en1,st2,en2) 
    integer, intent(IN) :: st1, en1, st2, en2
    integer, allocatable, dimension(:,:), intent(INOUT) :: v
    integer :: err
#   define vec_dim 2
#   include "arr_alloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine allocate_array_3d_integer(v,st1,en1,st2,en2,st3,en3) 
    integer, intent(IN) :: st1, en1, st2, en2, st3, en3
    integer, allocatable, dimension(:,:,:), intent(INOUT) :: v
    integer :: err
#   define vec_dim 3
#   include "arr_alloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine allocate_array_1d_logical(v,st,en) 
    integer, intent(IN) :: st, en
    logical, allocatable, dimension(:), intent(INOUT) :: v
    integer :: err
#   define vec_dim 1
#   include "arr_alloc_H.F90"
#   undef vec_dim
  end subroutine
  
  subroutine allocate_array_2d_logical(v,st1,en1,st2,en2) 
    integer, intent(IN) :: st1, en1, st2, en2
    logical, allocatable, dimension(:,:), intent(INOUT) :: v
    integer :: err
#   define vec_dim 2
#   include "arr_alloc_H.F90"
#   undef vec_dim
  end subroutine
  
  subroutine allocate_array_3d_logical(v,st1,en1,st2,en2,st3,en3) 
    integer, intent(IN) :: st1, en1, st2, en2, st3, en3
    logical, allocatable, dimension(:,:,:), intent(INOUT) :: v
    integer :: err
#   define vec_dim 3
#   include "arr_alloc_H.F90"
#   undef vec_dim
  end subroutine
!========================================================================================!

!                     deallocate 
!========================================================================================!
  subroutine deallocate_array_1d_complex(v) 
    complex(double_p), allocatable, dimension(:), intent(INOUT) :: v
#   include "deallocate_arr_H.F90"
  end subroutine

  subroutine deallocate_array_2d_complex(v) 
    complex(double_p), allocatable, dimension(:,:), intent(INOUT) :: v
#   include "deallocate_arr_H.F90"
  end subroutine
  
  subroutine deallocate_array_3d_complex(v) 
    complex(double_p), allocatable, dimension(:,:,:), intent(INOUT) :: v
#   include "deallocate_arr_H.F90"
  end subroutine

  subroutine deallocate_array_1d_double(v) 
    real(double_p), allocatable, dimension(:), intent(INOUT) :: v
#   include "deallocate_arr_H.F90"
  end subroutine

  subroutine deallocate_array_2d_double(v) 
    real(double_p), allocatable, dimension(:,:), intent(INOUT) :: v
#   include "deallocate_arr_H.F90"
  end subroutine

  subroutine deallocate_array_3d_double(v) 
    real(double_p), allocatable, dimension(:,:,:), intent(INOUT) :: v
#   include "deallocate_arr_H.F90"
  end subroutine

  subroutine deallocate_array_1d_integer(v) 
    integer, allocatable, dimension(:), intent(INOUT) :: v
#   include "deallocate_arr_H.F90"
  end subroutine
  
  subroutine deallocate_array_2d_integer(v) 
    integer, allocatable, dimension(:,:), intent(INOUT) :: v
#   include "deallocate_arr_H.F90"
  end subroutine

  subroutine deallocate_array_3d_integer(v) 
    integer, allocatable, dimension(:,:,:), intent(INOUT) :: v
#   include "deallocate_arr_H.F90"
  end subroutine

  subroutine deallocate_array_1d_logical(v) 
    logical, allocatable, dimension(:), intent(INOUT) :: v
#   include "deallocate_arr_H.F90"
  end subroutine
  
  subroutine deallocate_array_2d_logical(v) 
    logical, allocatable, dimension(:,:), intent(INOUT) :: v
#   include "deallocate_arr_H.F90"
  end subroutine

  subroutine deallocate_array_3d_logical(v) 
    logical, allocatable, dimension(:,:,:), intent(INOUT) :: v
#   include "deallocate_arr_H.F90"
  end subroutine
!========================================================================================!

!                     reAllocate 
!========================================================================================!
  subroutine reallocate_array_1d_complex(v,st,en) 
    integer, intent(IN) :: st, en
    complex(double_p), allocatable, dimension(:), intent(INOUT) :: v
#   define vec_dim 1
#   include "arr_realloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine reallocate_array_2d_complex(v,st1,en1,st2,en2)
    integer, intent(IN) :: st1, en1, st2, en2
    complex(double_p), allocatable, dimension(:,:), intent(INOUT) :: v
#   define vec_dim 2
#   include "arr_realloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine reallocate_array_3d_complex(v,st1,en1,st2,en2,st3,en3)
    integer, intent(IN) :: st1, en1, st2, en2, st3, en3
    complex(double_p), allocatable, dimension(:,:,:), intent(INOUT) :: v
#   define vec_dim 3
#   include "arr_realloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine reallocate_array_1d_double(v,st,en) 
    integer, intent(IN) :: st, en
    real(double_p), allocatable, dimension(:), intent(INOUT) :: v
#   define vec_dim 1
#   include "arr_realloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine reallocate_array_2d_double(v,st1,en1,st2,en2)
    integer, intent(IN) :: st1, en1, st2, en2
    real(double_p), allocatable, dimension(:,:), intent(INOUT) :: v
#   define vec_dim 2
#   include "arr_realloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine reallocate_array_3d_double(v,st1,en1,st2,en2,st3,en3)
    integer, intent(IN) :: st1, en1, st2, en2, st3, en3
    real(double_p), allocatable, dimension(:,:,:), intent(INOUT) :: v
#   define vec_dim 3
#   include "arr_realloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine reallocate_array_1d_integer(v,st,en) 
    integer, intent(IN) :: st, en
    integer, allocatable, dimension(:), intent(INOUT) :: v
#   define vec_dim 1
#   include "arr_realloc_H.F90"
#   undef vec_dim
  end subroutine
  
  subroutine reallocate_array_2d_integer(v,st1,en1,st2,en2)
    integer, intent(IN) :: st1, en1, st2, en2
    integer, allocatable, dimension(:,:), intent(INOUT) :: v
#   define vec_dim 2
#   include "arr_realloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine reallocate_array_3d_integer(v,st1,en1,st2,en2,st3,en3)
    integer, intent(IN) :: st1, en1, st2, en2, st3, en3
    integer, allocatable, dimension(:,:,:), intent(INOUT) :: v
#   define vec_dim 3
#   include "arr_realloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine reallocate_array_1d_logical(v,st,en)
    integer, intent(IN) :: st, en
    logical, allocatable, dimension(:), intent(INOUT) :: v
#   define vec_dim 1
#   include "arr_realloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine reallocate_array_2d_logical(v,st1,en1,st2,en2)
    integer, intent(IN) :: st1, en1, st2, en2
    logical, allocatable, dimension(:,:), intent(INOUT) :: v
#   define vec_dim 2
#   include "arr_realloc_H.F90"
#   undef vec_dim
  end subroutine

  subroutine reallocate_array_3d_logical(v,st1,en1,st2,en2,st3,en3)
    integer, intent(IN) :: st1, en1, st2, en2, st3, en3
    logical, allocatable, dimension(:,:,:), intent(INOUT) :: v
#   define vec_dim 3
#   include "arr_realloc_H.F90"
#   undef vec_dim
  end subroutine
!========================================================================================!

!                              copy 
!========================================================================================!
  subroutine copyArray_1d_complex(v,cpy) 
    complex(double_p), allocatable, dimension(:), intent(IN) :: v
    complex(double_p), allocatable, dimension(:), intent(INOUT) :: cpy
#   define vec_dim 1
#   include "arr_copy_H.F90"
#   undef vec_dim
  end subroutine

  subroutine copyArray_2d_complex(v,cpy)
    complex(double_p), allocatable, dimension(:,:), intent(IN) :: v
    complex(double_p), allocatable, dimension(:,:), intent(INOUT) :: cpy
#   define vec_dim 2
#   include "arr_copy_H.F90"
#   undef vec_dim
  end subroutine

  subroutine copyArray_3d_complex(v,cpy)
    complex(double_p), allocatable, dimension(:,:,:), intent(IN) :: v
    complex(double_p), allocatable, dimension(:,:,:), intent(INOUT) :: cpy
#   define vec_dim 3
#   include "arr_copy_H.F90"
#   undef vec_dim
  end subroutine

  subroutine copyArray_1d_double(v,cpy) 
    real(double_p), allocatable, dimension(:), intent(IN) :: v
    real(double_p), allocatable, dimension(:), intent(INOUT) :: cpy
#   define vec_dim 1
#   include "arr_copy_H.F90"
#   undef vec_dim
  end subroutine

  subroutine copyArray_2d_double(v,cpy)
    real(double_p), allocatable, dimension(:,:), intent(IN) :: v
    real(double_p), allocatable, dimension(:,:), intent(INOUT) :: cpy
#   define vec_dim 2
#   include "arr_copy_H.F90"
#   undef vec_dim
  end subroutine

  subroutine copyArray_3d_double(v,cpy)
    real(double_p), allocatable, dimension(:,:,:), intent(IN) :: v
    real(double_p), allocatable, dimension(:,:,:), intent(INOUT) :: cpy
#   define vec_dim 3
#   include "arr_copy_H.F90"
#   undef vec_dim
  end subroutine

  subroutine copyArray_1d_integer(v,cpy) 
    integer, allocatable, dimension(:), intent(IN) :: v
    integer, allocatable, dimension(:), intent(INOUT) :: cpy
#   define vec_dim 1
#   include "arr_copy_H.F90"
#   undef vec_dim
  end subroutine
  
  subroutine copyArray_2d_integer(v,cpy)
    integer, allocatable, dimension(:,:), intent(IN) :: v
    integer, allocatable, dimension(:,:), intent(INOUT) :: cpy
#   define vec_dim 2
#   include "arr_copy_H.F90"
#   undef vec_dim
  end subroutine

  subroutine copyArray_3d_integer(v,cpy)
    integer, allocatable, dimension(:,:,:), intent(IN) :: v
    integer, allocatable, dimension(:,:,:), intent(INOUT) :: cpy
#   define vec_dim 3
#   include "arr_copy_H.F90"
#   undef vec_dim
  end subroutine

  subroutine copyArray_1d_logical(v,cpy)
    logical, allocatable, dimension(:), intent(IN) :: v
    logical, allocatable, dimension(:), intent(INOUT) :: cpy
#   define vec_dim 1
#   include "arr_copy_H.F90"
#   undef vec_dim
  end subroutine

  subroutine copyArray_2d_logical(v,cpy)
    logical, allocatable, dimension(:,:), intent(IN) :: v
    logical, allocatable, dimension(:,:), intent(INOUT) :: cpy
#   define vec_dim 2
#   include "arr_copy_H.F90"
#   undef vec_dim
  end subroutine

  subroutine copyArray_3d_logical(v,cpy)
    logical, allocatable, dimension(:,:,:), intent(IN) :: v
    logical, allocatable, dimension(:,:,:), intent(INOUT) :: cpy
#   define vec_dim 3
#   include "arr_copy_H.F90"
#   undef vec_dim
  end subroutine
!========================================================================================!

  
end module allocate_arrays_mod