!!{This module sets up the equations to solve - equation equation has three subroutines.
!  Equation, BC_Bot and BC_Top - representing the equation and both boundaries respectively
!
!  Sets equations
!  Ax1 u_xx + Bx1 u_x + Ay1 u_yy + By1 u_y + C1 u + F1(u, v) = D1 or D1 u_t
!  Ax2 v_xx + Bx2 v_x + Ay2 v_yy + By2 v_y + C1 v + F2(u, v) = D2 or D2 v_t
!
!!
Module equations_definition
   use type_kinds, only: dp
   use maths_constants, only: pi, ex

contains

!!
! @brief      Sets up the first equation in the form 
!              Ax u_xx + Bx u_x + Ay u_yy + By u_y + C u = D or D u_t
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
!
! @return      Ax     second order x derivative coefficient
! @return      Bx     first order x derivative coefficient
! @return      Ay     second order y derivative coefficient
! @return      By     first order y derivative coefficient
! @return      C     zeroth order derivative coefficient
! @return      D     RHS/inhomogenous term
!
!!
   Subroutine equation1_linear(x, y, Ax, Bx, Ay, By, C, D)
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: Ax, Bx, Ay, By, C, D
      real(dp) :: epsi

      epsi = 0.05d0

      Ax = epsi
      Bx = 0.d0

      Ay = epsi
      By = 0.d0
      
      C = 0.d0
      D = 1.d0

   End Subroutine equation1_linear

!!
! @brief      Sets up the nonlinear terms F in the equation 
!             Ax u_xx + Bx u_x + Ay u_yy + By u_y + C u + F(u,v) = D or D u_t
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
! @param       u     non-linear variable
! @param       v     non-linear variable
!
! @return      F     Non-linear output (F is a function of u and v)
! @return      Fu    u - Partial Derivative of Non-linear output (F is a function of u and v)
! @return      Fv    v - Partial Derivative of Non-linear output (F is a function of u and v)
!
!!
   Subroutine equation1_non_linear(x, y, u, v, F, Fu, Fv)
      real(dp), intent(in) :: x, y, u, v
      real(dp), intent(out) :: F, Fu, Fv
      !!! F is a function of u - Fu is dF(u,v)/du and Fv is dF(u,v)/dv
      !!!   Non linear terms do not effect boundaries

      F = u*(1.d0-u) - v*u/(1.d0+u)
      Fu = -2.d0*u + 1.d0 - v/(1.d0+u) + v*u/((1.d0+u)**2.d0)
      Fv = -u/(1.d0+u)
   End Subroutine equation1_non_linear

!!
! @brief      Sets up the bottom/lefthand boundary in the first equation in the form 
!             Ax u_xx + Bx u_x + Ay u_yy + By u_y + C u = D or D u_t
!             
!             When x = xl
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
!
! @return      Ax     second order x derivative coefficient
! @return      Bx     first order x derivative coefficient
! @return      Ay     second order y derivative coefficient
! @return      By     first order y derivative coefficient
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
! @brief      Sets up the top/righthand boundary in the first equation in the form 
!             Ax u_xx + Bx u_x + Ay u_yy + By u_y + C u = D or D u_t
!             
!             when x = xr
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
!
! @return      Ax     second order x derivative coefficient
! @return      Bx     first order x derivative coefficient
! @return      Ay     second order y derivative coefficient
! @return      By     first order y derivative coefficient
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

!!
! @brief      Sets up the top/righthand boundary in the first equation in the form 
!             Ax u_xx + Bx u_x + Ay u_yy + By u_y + C u = D or D u_t
!             
!             when y = yl
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
!
! @return      Ax     second order x derivative coefficient
! @return      Bx     first order x derivative coefficient
! @return      Ay     second order y derivative coefficient
! @return      By     first order y derivative coefficient
! @return      C     zeroth order derivative coefficient
! @return      D     RHS/inhomogenous term
!
!!
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
! @brief      Sets up the top/righthand boundary in the first equation in the form 
!             Ax u_xx + Bx u_x + Ay u_yy + By u_y + C u = D or D u_t
!             
!             when y = yr
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
!
! @return      Ax     second order x derivative coefficient
! @return      Bx     first order x derivative coefficient
! @return      Ay     second order y derivative coefficient
! @return      By     first order y derivative coefficient
! @return      C     zeroth order derivative coefficient
! @return      D     RHS/inhomogenous term
!
!!
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
! @brief      Sets up the the initial conidiotn in the form u = F. Note in 1D y = 1
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
!
! @return      IC     Initial condition
!
!!
   Subroutine equation1_initial_condition(x, y, IC)
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: IC
      real(dp) :: r

      !IC = sin(2.d0*pi*x)
      ! IC = sin(pi*x) + sin(pi*y)
       !IC = cos(2.d0*pi*x*y)**2.d0
      !IC = 5.d0*x*y
      IC = (cos(4.d0*x*pi))**2.d0*(ex**(-(x**2.d0)))
      
     ! Call random_seed()
     ! Call random_number(r)
      !IC = 1.d0+(2.d0*r-1.d0)*0.01d0!*abs(sin(5*x*pi))*0.1d0
      !If (x.gt.0.5d0) then
       !  IC = 1.d0 + 0.001*cos(6.d0*x*pi)
      !End If

   End Subroutine equation1_initial_condition


!!
! @brief      Sets up the second equation in the form 
!              Ax v_xx + Bx v_x + Ay v_yy + By v_y + C v = D or D v_t
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
!
! @return      Ax     second order x derivative coefficient
! @return      Bx     first order x derivative coefficient
! @return      Ay     second order y derivative coefficient
! @return      By     first order y derivative coefficient
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
! @brief      Sets up the nonlinear terms F in the equation 
!             Ax v_xx + Bx v_x + Ay v_yy + By v_y + C v + F(u,v) = D or D v_t
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
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
      real(dp) :: mu

      F = v*u/(1.d0+u)
      Fu = v/(1.d0+u) - v*u/((1.d0+u)**2.d0)
      Fv = u/(1.d0+u)
      
      !mu = 1.8d0
      !F = mu*((1.d0 - u*u*v))
      !Fu = mu*(-2.d0*u*v)
      !Fv = mu*(-u*u)

   End Subroutine equation2_non_linear
!!
! @brief      Sets up the bottom/lefthand boundary in the second equation in the form 
!             Ax v_xx + Bx v_x + Ay v_yy + By v_y + C v = D or D v_t
!             
!             x = xl
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
!
! @return      Ax     second order x derivative coefficient
! @return      Bx     first order x derivative coefficient
! @return      Ay     second order y derivative coefficient
! @return      By     first order y derivative coefficient
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
! @brief      Sets up the top/righthand boundary in the second equation in the form 
!             Ax v_xx + Bx v_x + Ay v_yy + By v_y + C v = D or D v_t
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
!
! @return      Ax     second order x derivative coefficient
! @return      Bx     first order x derivative coefficient
! @return      Ay     second order y derivative coefficient
! @return      By     first order y derivative coefficient
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
! @brief      Sets up the bottom/lefthand boundary in the second equation in the form
!             Ax v_xx + Bx v_x + Ay v_yy + By v_y + C v = D or D v_t
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
!
! @return      Ax     second order x derivative coefficient
! @return      Bx     first order x derivative coefficient
! @return      Ay     second order y derivative coefficient
! @return      By     first order y derivative coefficient
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
! @brief      Sets up the top/righthand boundary in the second equation in the form 
!             Ax v_xx + Bx v_x + Ay v_yy + By v_y + C v = D or D v_t
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
!
! @return      Ax     second order x derivative coefficient
! @return      Bx     first order x derivative coefficient
! @return      Ay     second order y derivative coefficient
! @return      By     first order y derivative coefficient
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
! @brief      Sets up the the initial conidiotn in the form v = F. Note in 1D y = 1
!
! @param       x     physical domain input - can have variable coefficients with x
! @param       y     physical domain input - can have variable coefficients with y
!
! @return      IC     Initial condition
!
!!
   Subroutine equation2_initial_condition(x, y, IC)
      real(dp), intent(in) :: x, y ! xpoisiton in the domain
      real(dp), intent(out) :: IC
      real(dp) :: r

      !Call random_number(r)
      !IC = 1.d0 !+ 0.01d0*r
      IC = (0.2 + 0.1*sin(pi*x*10.d0))**2.d0*(ex**(((x+0.8)**2.d0)))

   End Subroutine equation2_initial_condition

End Module equations_definition

