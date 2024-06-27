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
      WRITE(6,*) 'Ellipitc Solver'
      Call solve_runner
      Stop
   Case (1)
      WRITE(6,*) 'Parabolic Solver'
      Call march_runner
      stop
   End Select

End Program main





