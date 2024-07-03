!!{This module sets up the test equations to solve - equation equation has three subroutines.
!  Equation, BC_Bot and BC_Top - representing the equation and both boundaries respectively}
!!
Module equations_definition_test
   use type_kinds, only: dp
   use maths_constants, only: pi, ex

contains

!!
! @brief      Sets up the first equation in the form A u_xx + B u_x + C u = D
!
! @param       x     physical domain input - can have variable coefficients with x
!
! @return      A     second order derivative coefficient
! @return      B     first order derivative coefficient
! @return      C     zeroth order derivative coefficient
! @return      D     RHS/inhomogenous term
!
!!
   Subroutine equation1_test(x, A, B, C, D)
      real(dp), intent(in) :: x ! xpoisiton in the domain
      real(dp), intent(out) :: A, B, C, D

      A = 1.d0
      B = 0.d0
      C = -1.d0
      D = 0.d0

   End Subroutine equation1_test

!!
! @brief      Sets up the bottom/lefthand boundary in the first equation in the form A u_xx + B u_x + C u = D
!
! @param       x     physical domain input - can have variable coefficients with x
!
! @return      A     second order derivative coefficient
! @return      B     first order derivative coefficient
! @return      C     zeroth order derivative coefficient
! @return      D     RHS/inhomogenous term
!
!!
   Subroutine equation1_BC_Bot_test(x, A, B, C, D)
      real(dp), intent(in) :: x ! xpoisiton in the domain
      real(dp), intent(out) :: A, B, C, D

      A = 0.d0
      B = 0.d0
      C = 1.d0
      D = 1.d0

   End Subroutine equation1_BC_Bot_test

!!
! @brief      Sets up the top/righthand boundary in the first equation in the form A u_xx + B u_x + C u = D
!
! @param       x     physical domain input - can have variable coefficients with x
!
! @return      A     second order derivative coefficient
! @return      B     first order derivative coefficient
! @return      C     zeroth order derivative coefficient
! @return      D     RHS/inhomogenous term
!
!!
   Subroutine equation1_BC_Top_test(x, A, B, C, D)
      real(dp), intent(in) :: x ! xpoisiton in the domain
      real(dp), intent(out) :: A, B, C, D

      A = 0.d0
      B = 0.d0
      C = 1.d0
      D = ex

   End Subroutine equation1_BC_Top_test

End Module equations_definition_test

