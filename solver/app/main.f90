!!Program that builds and solves PDEs
!!
Program main
   !use type_kinds, only: sp,dp
   use reader, only: read_me
   use domain
   use maths_constants
   use equations
   use linear_algebra, only: solver_banded_double_precision
   implicit none

   real(dp), dimension(:), allocatable :: X
   integer :: i
   real(dp), dimension(:, :), allocatable :: L !! banded form matrix
   real(dp), dimension(:), allocatable :: RHS !! right hand side of equation

   Write (6, *)
   Write (6, *) '!!!!!!!!!!!!!!!!!!!!!!!!!'
   Write (6, *) '!!!!!!!!!!!!!!!!!!!!!!!!!'
   Write (6, *) '!!!!!! PDE SOLVER !!!!!!!'
   Write (6, *) '!!!!!!!!!!!!!!!!!!!!!!!!!'
   Write (6, *) '!!!!!!!!!!!!!!!!!!!!!!!!!'
   Write (6, *)

  !!!! call intitial subroutines

   call read_me  ! opens settings.input and reads and sets settings
   call diff_initialisation  ! sets the finite difference coefficients
   call initial_domain_settings !builds the domain and computational domains

   Write (6, *) 'size of x domain::: ', nx
   Write (6, *) 'order of the finite differences', DiffOrder
   Write (6, *) '... Building matrix'

   call build_the_matrix(nx, dxc, dxcsq, xcdom, xmetric1, xmetric1sq, xmetric2, 1, L, RHS)
   call solver_banded_double_precision(nx, nband, sub_diag, sup_diag, L, RHS, X)
   deallocate (L, RHS)

   Write (6, *)
   Write (6, '(5(A20,x))') 'Physical domain', 'Comp Domain', 'Numerical Soln', 'Exact Soln', 'error'
   Do i = 1, nx
      Write (6, '(4(f20.14,1x),e20.10)') xdom(i), xcdom(i), x(i), ex**xdom(i), abs(x(i) - ex**xdom(i))
   End Do
   Write (6, *)
End Program main
