!!Program that builds and solves PDEs
!!
Program main
   use type_kinds
   use reader
   use domain
   use maths_constants
   use equations
   use linear_algebra
   use solve_bvp, only : solve_runner
   use temporal_marching, only : march_runner
   implicit none


   Write (6, *)
   Write (6, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   Write (6, *) '!!! Partial Differentiation Solver !!!'
   Write (6, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   Write (6, *)

  !!!! Call intitial subroutines

   Call read_me  ! opens settings.input and reads and sets settings
   Call diff_initialisation  ! sets the finite difference coefficients
   Call initial_domain_settings !builds the domain and computational domains
   Call reset_domain_paramters ! rescales the domain terms

   Select Case (Time_switch)
   Case (0)
      Write (6, *) 'Ellipitc Solver'
      Call solve_runner
      Stop
   Case (1)
      Write (6, *) 'Parabolic Solver'
      Call march_runner
      Stop
   End Select

End Program main






