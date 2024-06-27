!!{The module that set up and runs the BVP Solver}
!!
Module solve_bvp
   use type_kinds
   use reader
   use domain
   use maths_constants
   use equations
   use linear_algebra

contains

!!
! @brief      {Subroutine that sets up a BVP solver for the equation given in definitions.f90}
!
! @return     Solves the BVP in the form A u_xx + B u_x + C = D. Boundary conditions can be general
!!
   Subroutine solve_runner

      integer :: i, j
      real(dp), dimension(:), allocatable :: X
      real(dp), dimension(:, :), allocatable :: L !! banded form matrix
      real(dp), dimension(:), allocatable :: RHS !! right hand side of equation

      Write (6, *) 'size of x domain::: ', nx
      Write (6, *) 'order of the finite differences', DiffOrder
      Write (6, *) '... Building matrix'
      Write (6, *)

      !!! domain has been built in main
      !!! next step is to build the discretised operators
      !!! equation is in the form L * X = RHS

      Call build_the_matrix(nx, dxc, dxcsq, xcdom, xmetric1, xmetric1sq, xmetric2, 1, L, RHS)

      !! Solve the equation

      Call solver_banded_double_precision(nx, nband, sub_diag, sup_diag, L, RHS, X)
      deallocate (L, RHS)

      !! Print the results

      Open (10, file='BVP.dat')

      Write (6, '(5(A20,x))') 'Physical domain', 'Comp Domain', 'Numerical Soln'!, 'Exact Soln', 'error'
      Write (10, '(5(A20,x))') 'Physical domain', 'Comp Domain', 'Numerical Soln'!, 'Exact Soln', 'error'
      Do i = 1, nx
         Write (6, '(4(f20.14,1x),e20.10)') xdom(i), xcdom(i), X(i)!, ex**xdom(i), abs(x(i)-ex**xdom(i))
         Write (10, '(4(f20.14,1x),e20.10)') xdom(i), xcdom(i), X(i)!, ex**xdom(i), abs(x(i)-ex**xdom(i))
      End Do
      Close (10)
      Write (6, *)

      deallocate (X)

      Return

   End Subroutine solve_runner

End Module solve_bvp

