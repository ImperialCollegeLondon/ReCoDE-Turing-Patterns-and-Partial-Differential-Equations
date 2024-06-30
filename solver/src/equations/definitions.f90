!!{This module sets up the equations to solve - equation equation has three subroutines.
!  Equation, BC_Bot and BC_Top - representing the equation and both boundaries respectively
!
!  Sets equations
!  A u_xx + B u_x + C u + F(u) = D or D u_t
!
!
!!
Module equations_definition
   use type_kinds, only: dp
   use maths_constants, only: pi, ex

contains

!!
! @brief      Sets up the first equation in the form A u_xx + B u_x + C u = D or D u_t
!
! @param       x     physical domain input - can have variable coefficients with x
!
! @return      A     second order derivative coefficient
! @return      B     first order derivative coefficient
! @return      C     zeroth order derivative coefficient
! @return      D     RHS/inhomogenous term
!
!!
   Subroutine equation1_linear(x, A, B, C, D)
      real(dp), intent(in) :: x ! xpoisiton in the domain
      real(dp), intent(out) :: A, B, C, D
      real(dp) :: epsi

      epsi = 0.001

      A = epsi
      B = 0.d0
      C = 0.d0
      D = 1.d0

   End Subroutine equation1_linear

!!
! @brief      Sets up the nonlinear terms F in the equation A u_xx + B u_x + C u + F(u) = D or D u_t
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       u     non-linear variable
! @param       v     non-linear variable
!
! @return      F     Non-linear output (F is a function of u and v)
! @return      Fu    u - Derivative of Non-linear output (F is a function of u and v)
! @return      Fv    v - Derivative of Non-linear output (F is a function of u and v)
!
!!
   Subroutine equation1_non_linear(x, u, v, F, Fu, Fv)
      real(dp), intent(in) :: x, u, v
      real(dp), intent(out) :: F, Fu, Fv
!!! F is a function of u - Fu is F'(u)
!!!

      F = u - u**2.d0 - 1.2d0*v*u/(0.2d0 + u)
      Fu = 1.d0 - 2.d0*u - 1.2d0*(v/(0.2d0 + u) - v*u/((0.2d0 + u)**2.d0))
      Fv = -1.2d0*u/(0.2d0 + u)
   End Subroutine equation1_non_linear

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
   Subroutine equation1_BC_Bot(x, A, B, C, D)
      real(dp), intent(in) :: x ! xpoisiton in the domain
      real(dp), intent(out) :: A, B, C, D

      A = 0.d0
      B = 1.d0
      C = 0.d0
      D = 0.d0

   End Subroutine equation1_BC_Bot

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
   Subroutine equation1_BC_Top(x, A, B, C, D)
      real(dp), intent(in) :: x ! xpoisiton in the domain
      real(dp), intent(out) :: A, B, C, D

      A = 0.d0
      B = 1.d0
      C = 0.d0
      D = 0.d0

   End Subroutine equation1_BC_Top

!!
! @brief      Sets up the the initial conidiotn in the form u = F
!
! @param       x     physical domain input - can have variable coefficients with x
!
! @return      IC     Initial condition
!
!!
   Subroutine equation1_initial_condition(x, IC)
      real(dp), intent(in) :: x ! xpoisiton in the domain
      real(dp), intent(out) :: IC

      !IC = sin(2.d0*pi*x)
      IC = (cos(pi*x*4.d0)**2.d0)

   End Subroutine equation1_initial_condition

!!
! @brief      Sets up the second equation in the form A v_xx + B v_x + C v = D or D v_t
!
! @param       x     physical domain input - can have variable coefficients with x
!
! @return      A     second order derivative coefficient
! @return      B     first order derivative coefficient
! @return      C     zeroth order derivative coefficient
! @return      D     RHS/inhomogenous term
!
!!
   Subroutine equation2_linear(x, A, B, C, D)
      real(dp), intent(in) :: x ! xpoisiton in the domain
      real(dp), intent(out) :: A, B, C, D

      A = 0.1d0
      B = 0.d0
      C = 0.d0
      D = 1.d0

   End Subroutine equation2_linear

!!
! @brief      Sets up the nonlinear terms F in the equation A v_xx + B v_x + C v + G(u,v) = D or D v_t
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       u     non-linear variable
! @param       v     non-linear variable
!
! @return      G     Non-linear output (G is a function of u and v)
! @return      Gu    u - Derivative of Non-linear output (Gu is a function of u and v)
! @return      Gv    v - Derivative of Non-linear output (Gv is a function of u and v)
!
!!
   Subroutine equation2_non_linear(x, u, v, G, Gu, Gv)
      real(dp), intent(in) :: x, u, v
      real(dp), intent(out) :: G, Gu, Gv

      G = 1.2d0*v*(u/(0.2d0 + u) + 0.d0)
      Gu = 1.2d0*v/(0.2d0 + u) - v*u/((0.2d0 + u)**2.d0)
      Gv = 1.2d0*(u/(0.2d0 + u) + 0.d0)

   End Subroutine equation2_non_linear

!!
! @brief      Sets up the bottom/lefthand boundary in the second equation in the form A v_xx + B v_x + C v = D
!
! @param       x     physical domain input - can have variable coefficients with x
!
! @return      A     second order derivative coefficient
! @return      B     first order derivative coefficient
! @return      C     zeroth order derivative coefficient
! @return      D     RHS/inhomogenous term
!
!!
   Subroutine equation2_BC_Bot(x, A, B, C, D)
      real(dp), intent(in) :: x ! xpoisiton in the domain
      real(dp), intent(out) :: A, B, C, D

      A = 0.d0
      B = 1.d0
      C = 0.d0
      D = 0.d0

   End Subroutine equation2_BC_Bot

!!
! @brief      Sets up the top/righthand boundary in the second equation in the form A v_xx + B v_x + C v = D
!
! @param       x     physical domain input - can have variable coefficients with x
!
! @return      A     second order derivative coefficient
! @return      B     first order derivative coefficient
! @return      C     zeroth order derivative coefficient
! @return      D     RHS/inhomogenous term
!
!!
   Subroutine equation2_BC_Top(x, A, B, C, D)
      real(dp), intent(in) :: x ! xpoisiton in the domain
      real(dp), intent(out) :: A, B, C, D

      A = 0.d0
      B = 1.d0
      C = 0.d0
      D = 0.d0

   End Subroutine equation2_BC_Top

!!
! @brief      Sets up the the initial conidiotn in the form u = F
!
! @param       x     physical domain input - can have variable coefficients with x
!
! @return      IC     Initial condition
!
!!
   Subroutine equation2_initial_condition(x, IC)
      real(dp), intent(in) :: x ! xpoisiton in the domain
      real(dp), intent(out) :: IC

      ! IC = sin(2.d0*pi*x)
      If (x .gt. 0.5d0) Then
         IC = 0.d0
      Else
         IC = (cos(pi*x)**2.d0)
      End If

   End Subroutine equation2_initial_condition

End Module equations_definition

