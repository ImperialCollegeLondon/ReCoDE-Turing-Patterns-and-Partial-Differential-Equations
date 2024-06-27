!!Program that builds and solves PDEs
!!
Program main
   use type_kinds
   use reader
   use domain
   use maths_constants
   use equations
   use linear_algebra
   use solve_ode, only : solver
   use temporal_marching, only : march
   implicit none
   external :: DGBMV

   integer :: i,j
   real(dp), dimension(:), allocatable :: X,X1
   real(dp), dimension(:,:), allocatable :: Soln
   real(dp), dimension(:, :), allocatable :: L !! banded form matrix
   real(dp), dimension(:), allocatable :: RHS,temp !! right hand side of equation

   Write (6, *)
   Write (6, *) '!!!!!!!!!!!!!!!!!!!!!!!!!'
   Write (6, *) '!!!!!!!!!!!!!!!!!!!!!!!!!'
   Write (6, *) '!!!!!! PDE SOLVER !!!!!!!'
   Write (6, *) '!!!!!!!!!!!!!!!!!!!!!!!!!'
   Write (6, *) '!!!!!!!!!!!!!!!!!!!!!!!!!'
   Write (6, *)

  !!!! Call intitial subroutines

   Call read_me  ! opens settings.input and reads and sets settings
   Call diff_initialisation  ! sets the finite difference coefficients
   Call initial_domain_settings !builds the domain and computational domains


   Select Case (Time_switch)
   Case (0)
      Call solver
      Stop
   Case (1)
      Call march
      stop
   End Select


End Program main






