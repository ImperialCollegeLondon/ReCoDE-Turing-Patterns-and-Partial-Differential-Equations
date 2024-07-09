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
   use solve_ibvp, only : march_runner
   implicit none

   real :: cpu_start, cpu_end

   call cpu_time(cpu_start)

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
      Write (6, *) 'Boundary Value Problem Solver'
      Call solve_runner
      !Stop
   Case (1)
      Write (6, *) 'Temporal Initial Boundary Value Solver'
      Call march_runner
      !Stop
   End Select

   call cpu_time(cpu_end)

     Write (6, '(A, F8.3, A)') 'Took: ', cpu_end - cpu_start, 's'


End Program main






