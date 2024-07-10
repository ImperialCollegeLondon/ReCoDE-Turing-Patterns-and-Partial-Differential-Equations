!!Program that builds and solves PDEs
!!
Program main
   use type_kinds
   use reader
   use domain
   use maths_constants
   use equations
   use linear_algebra
   use solve_bvp, only : bvp_runner
   use solve_ibvp, only : ibvp_runner
   implicit none

   real :: cpu_start, cpu_end

   !! set computational time begining
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


   Write (6, *)
   Write (6, *) 'size of x Domain::: ', nx, dxc
  
   Select Case(Domain_number)
   Case(2)
      Write (6, *) 'size of y Domain::: ', ny, dyc
   End Select
 
   Write (6, *) 'number of equations:::', Eqn_number
   Write (6, *)
   Write (6, *) 'size of total Domain::: ', idim
   Write (6, *)
   Write (6, *) 'size of t Domain::: ', nt, dt
   Write (6, *)
   Write (6, *) 'order of the finite differences', DiffOrder
   Write (6, *)
  
   Select Case(Non_Linear_switch)
   Case(1)
      Write (6, *) 'Non-linear iteration at error',Newton_Error
      Write (6, *)
   End Select


   !!! Decide which solver to call - temporal marching for BVP
   Select Case (Time_switch)
   Case (0)
      Write (6, *) '!!!!!! Boundary Value Problem Solver !!!!!!'
      Call bvp_runner(cpu_start)
      !Stop
   Case (1)
      Write (6, *) '!!!!!! Initial Boundary Value Problem Solver !!!!!!'
      Call ibvp_runner(cpu_start)
      !Stop
   End Select

   !! set final computational time
   call cpu_time(cpu_end)

   Write (6, '(A, F8.3, A)') ' Total time taken: ', cpu_end - cpu_start, 's'
   Write (6, *)
End Program main






