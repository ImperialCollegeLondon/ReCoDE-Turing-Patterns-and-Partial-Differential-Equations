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
!
!
!!
   Subroutine equation1_linear(x, y, Ax, Bx, Ay, By, C, D)
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: Ax, Bx, Ay, By, C, D
      real(dp) :: epsi

      epsi = 0.001

      Ax = 1.d0
      Bx = 0.d0

      Ay = 1.d0
      By = 0.d0

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
   Subroutine equation1_non_linear(x, y, u, v, F, Fu, Fv)
      real(dp), intent(in) :: x, y, u, v
      real(dp), intent(out) :: F, Fu, Fv
!!! F is a function of u - Fu is F'(u)
!!!

      F = u - u**2.d0 !- 1.2d0*v*u/(0.2d0 + u)
      Fu = 1.d0 - 2.d0*u !- 1.2d0*(v/(0.2d0 + u) - v*u/((0.2d0 + u)**2.d0))
      Fv = 0.d0! -1.2d0*u/(0.2d0 + u)
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
   Subroutine equation1_BC_X_Bot(x, y, Ax, Bx, Ay, By, C, D)
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: Ax, Bx, Ay, By, C, D

      Ax = 0.d0
      Bx = 1.d0

      Ay = 0.d0
      By = 0.d0

      C = 0.d0
      D = 0.d0

   End Subroutine equation1_BC_X_Bot

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
   Subroutine equation1_BC_X_Top(x, y, Ax, Bx, Ay, By, C, D)
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: Ax, Bx, Ay, By, C, D

      Ax = 0.d0
      Bx = 1.d0

      Ay = 0.d0
      By = 0.d0

      C = 0.d0
      D = 0.d0

   End Subroutine equation1_BC_X_Top

   Subroutine equation1_BC_Y_Bot(x, y, Ax, Bx, Ay, By, C, D)
      !!! Note that the boundary conditions at the corners are governed by BC_Y
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: Ax, Bx, Ay, By, C, D

      Ax = 0.d0
      Bx = 0.d0

      Ay = 0.d0
      By = 1.d0

      C = 0.d0
      D = 0.d0

   End Subroutine equation1_BC_Y_Bot

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
   Subroutine equation1_BC_Y_Top(x, y, Ax, Bx, Ay, By, C, D)
   !!! Note that the boundary conditions at the corners are governed by BC_Y
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: Ax, Bx, Ay, By, C, D

      Ax = 0.d0
      Bx = 0.d0

      Ay = 0.d0
      By = 1.d0

      C = 0.d0
      D = 0.d0

   End Subroutine equation1_BC_Y_Top

!!
! @brief      Sets up the the initial conidiotn in the form u = F
!
! @param       x     physical domain input - can have variable coefficients with x
!
! @return      IC     Initial condition
!
!!
   Subroutine equation1_initial_condition(x, y, IC)
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: IC

      !IC = sin(2.d0*pi*x)
      ! IC = sin(pi*x) + sin(pi*y)
      IC = cos(2.d0*pi*y*x)**2.d0
      !IC = 5.d0*x*y

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
   Subroutine equation2_linear(x, y, Ax, Bx, Ay, By, C, D)
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: Ax, Bx, Ay, By, C, D

      Ax = 1.d0
      Bx = 0.d0

      Ay = 1.d0
      By = 0.d0

      C = 0.d0
      D = 1.d0

   End Subroutine equation2_linear

!!
! @brief      Sets up the nonlinear terms F in the equation A v_xx + B v_x + C v + F(u,v) = D or D v_t
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       u     non-linear variable
! @param       v     non-linear variable
!
! @return      F     Non-linear output (F is a function of u and v)
! @return      Fu    u - Derivative of Non-linear output (Fu is a function of u and v)
! @return      Fv    v - Derivative of Non-linear output (Fv is a function of u and v)
!
!!
   Subroutine equation2_non_linear(x, y, u, v, F, Fu, Fv)
      real(dp), intent(in) :: x, y, u, v
      real(dp), intent(out) :: F, Fu, Fv

      F = u - u**2.d0 !- 1.2d0*v*u/(0.2d0 + u)
      Fu = 1.d0 - 2.d0*u !- 1.2d0*(v/(0.2d0 + u) - v*u/((0.2d0 + u)**2.d0))
      Fv = 0.d0! -1.2d0*u/(0.2d0 + u)

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
   Subroutine equation2_BC_X_Bot(x, y, Ax, Bx, Ay, By, C, D)
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: Ax, Bx, Ay, By, C, D

      Ax = 0.d0
      Bx = 1.d0

      Ay = 0.d0
      By = 0.d0

      C = 0.d0
      D = 0.d0

   End Subroutine equation2_BC_X_Bot

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
   Subroutine equation2_BC_X_Top(x, y, Ax, Bx, Ay, By, C, D)
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: Ax, Bx, Ay, By, C, D

      Ax = 0.d0
      Bx = 1.d0

      Ay = 0.d0
      By = 0.d0

      C = 0.d0
      D = 0.d0

   End Subroutine equation2_BC_X_Top

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
   Subroutine equation2_BC_Y_Bot(x, y, Ax, Bx, Ay, By, C, D)
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: Ax, Bx, Ay, By, C, D

      Ax = 0.d0
      Bx = 0.d0

      Ay = 0.d0
      By = 1.d0

      C = 0.d0
      D = 0.d0
   End Subroutine equation2_BC_Y_Bot

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
   Subroutine equation2_BC_Y_Top(x, y, Ax, Bx, Ay, By, C, D)
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: Ax, Bx, Ay, By, C, D

      Ax = 0.d0
      Bx = 0.d0

      Ay = 0.d0
      By = 1.d0

      C = 0.d0
      D = 0.d0

   End Subroutine equation2_BC_Y_Top

!!
! @brief      Sets up the the initial conidiotn in the form u = F
!
! @param       x     physical domain input - can have variable coefficients with x
!
! @return      IC     Initial condition
!
!!
   Subroutine equation2_initial_condition(x, y, IC)
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: IC

      !IC = sin(2.d0*pi*x)
      ! IC = sin(pi*x) + sin(pi*y)
      IC = cos(2.d0*pi*x*y)**2.d0
      !IC = 5.d0*x*y

   End Subroutine equation2_initial_condition

End Module equations_definition

