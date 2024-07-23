!! This module sets up the finite difference coefficients and relevenet matrix parameters
!!
Module maths_constants_diff_coeff
   use type_kinds, only: dp
   implicit none

   integer :: DiffOrder  ! order of the equation, 2, 3 or 4
   integer :: sub_diag   ! number of sub_diagonals of the given matrices
   integer :: sup_diag   ! number of super_diagonals of the given matrices
   integer :: nband      ! resulting bandwith of the banded matrix
   real(dp), dimension(1:5, 1:6) :: D1, D2

   ! D1 and D2 are the coefficient finite diff matrices for first and second derivatives respecitely
   ! rows 1 and 2 and forward differences
   ! 3rd row is a central difference
   ! rows 4 and 5 are backward differenced
   ! note that in D2(3,4) is the central node

contains

!!
! @brief      This subroutine sets up the matrices D1 and D2 (the first and second derivatives finite difference coefficinets)
!!
   Subroutine diff_initialisation
      integer :: i

      D1 = 0.d0; D2 = 0.d0

      !We first set the first two forward differences and central and then use symmetry to obtain the backward

      Select Case (DiffOrder)
      Case (4)

         sub_diag = 5; sup_diag = 5; nband = sup_diag + sub_diag + 1

!!!! fourth order accurate everywhere

    !!!! First Derivative
    !! forward 1
         D1(1, 1) = -25.d0/12.d0; D1(1, 2) = 4.d0; D1(1, 3) = -3.d0; D1(1, 4) = 4.d0/3.d0; D1(1, 5) = -1.d0/4.d0
    !! forward 2
         D1(2, 1) = -1.d0/5.d0; D1(2, 2) = -13.d0/12.d0; D1(2, 3) = 2.d0; 
         D1(2, 4) = -1.d0; D1(2, 5) = 1.d0/3.d0; D1(2, 6) = -1.d0/20.d0
    !! Central
         D1(3, 2) = 1.d0/12.d0; D1(3, 3) = -2.d0/3.d0; D1(3, 4) = 0.d0; D1(3, 5) = 2.d0/3.d0; D1(3, 6) = -1.d0/12.d0

!!!! Second Derivative
    !! forward 1
         D2(1, 1) = 15.d0/4.d0; D2(1, 2) = -77.d0/6.d0; D2(1, 3) = 107.d0/6.d0; 
         D2(1, 4) = -13.d0; D2(1, 5) = 61.d0/12.d0; D2(1, 6) = -5.d0/6.d0
    !! forward 2
         D2(2, 1) = 5.d0/6.d0; D2(2, 2) = -5.d0/4.d0; D2(2, 3) = -1.d0/3.d0; 
         D2(2, 4) = 7.d0/6.d0; D2(2, 5) = -1.d0/2.d0; D2(2, 6) = 1.d0/12.d0
    !! Central
         D2(3, 1) = 0.d0; D2(3, 2) = -1.d0/12.d0; D2(3, 3) = 4.d0/3.d0; 
         D2(3, 4) = -5.d0/2.d0; D2(3, 5) = 4.d0/3.d0; D2(3, 6) = -1.d0/12.d0

      Case (3)

         sub_diag = 4; sup_diag = 4; nband = sup_diag + sub_diag + 1

!!!! fourth order accurate in the entirior
!!!! second order accurate in the boundaries
!!!! makes banded form equations slightly smaller

!!!! (3rd order accurate first derivative)
!!!! First Derivative
    !! forward 1
         D1(1, 1) = -11.d0/6.d0; D1(1, 2) = 3.d0; D1(1, 3) = -3.d0/2.d0; D1(1, 4) = 1.d0/3.d0
    !! forward 2
         D1(2, 1) = -1.d0/5.d0; D1(2, 2) = -13.d0/12.d0; D1(2, 3) = 2.d0; 
         D1(2, 4) = -1.d0; D1(2, 5) = 1.d0/3.d0; D1(2, 6) = -1.d0/20.d0
    !! Central
         D1(3, 2) = 1.d0/12.d0; D1(3, 3) = -2.d0/3.d0; D1(3, 4) = 0.d0; D1(3, 5) = 2.d0/3.d0; D1(3, 6) = -1.d0/12.d0

!!!! Second Derivatives
         D2(1, 1) = 35.d0/12.d0; D2(1, 2) = -26.d0/3.d0; D2(1, 3) = 19.d0/2.d0; 
         D2(1, 4) = -14.d0/3.d0; D2(1, 5) = 11.d0/12.d0
     !! forward 2
         D2(2, 1) = 5.d0/6.d0; D2(2, 2) = -5.d0/4.d0; D2(2, 3) = -1.d0/3.d0; 
         D2(2, 4) = 7.d0/6.d0; D2(2, 5) = -1.d0/2.d0; D2(2, 6) = 1.d0/12.d0
     !! Central
         D2(3, 1) = 0.d0; D2(3, 2) = -1.d0/12.d0; D2(3, 3) = 4.d0/3.d0; 
         D2(3, 4) = -5.d0/2.d0; D2(3, 5) = 4.d0/3.d0; D2(3, 6) = -1.d0/12.d0

      Case (2)

         sub_diag = 4; sup_diag = 4; nband = sup_diag + sub_diag + 1

!!!! second order accurate everywhere

!!!! First Derivative
    !! forward 1
         D1(1, 1) = -11.d0/6.d0; D1(1, 2) = 3.d0; D1(1, 3) = -3.d0/2.d0; D1(1, 4) = 1.d0/3.d0
    !! forward 2
         D1(2, 1) = -1.d0/3.d0; D1(2, 2) = -1.d0/2.d0; D1(2, 3) = 1.d0; D1(2, 4) = -1.d0/6.d0
    !! Central
         D1(3, 1) = 0.d0; D1(3, 3) = -1.d0/2.d0; D1(3, 4) = 0.d0; D1(3, 5) = 1.d0/2.d0

!!!! Second Derivatives
         D2(1, 1) = 35.d0/12.d0; D2(1, 2) = -26.d0/3.d0; D2(1, 3) = 19.d0/2.d0; 
         D2(1, 4) = -14.d0/3.d0; D2(1, 5) = 11.d0/12.d0
     !! forward 2
         D2(2, 1) = 1.d0; D2(2, 2) = -2.d0; D2(2, 3) = 1.d0; D2(2, 4) = 0.d0; D2(2, 5) = 0.d0; d2(2, 6) = 0.d0
     !! Central
         D2(3, 3) = 1.d0; D2(3, 4) = -2.d0; D2(3, 5) = 1.d0

      End Select

!! We set the backward from the forward using symmetry

      Do i = 1, 6
!backward 2
         D1(4, i) = -D1(2, 7 - i)
         D2(4, i) = D2(2, 7 - i)
!backward 1
         D1(5, i) = -D1(1, 7 - i)
         D2(5, i) = D2(1, 7 - i)
      End Do


   End Subroutine diff_initialisation

End Module maths_constants_diff_coeff
