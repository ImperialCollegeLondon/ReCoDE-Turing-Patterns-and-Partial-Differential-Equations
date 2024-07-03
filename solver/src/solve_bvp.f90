!!{The module that set up and runs the BVP Solver for both non-linear BVP and Linear BVP
!  Non-linear uses Newton iteration}
!!
Module solve_bvp
   use type_kinds
   use reader, only: Time_switch, Non_Linear_switch, Eqn_number
   use domain
   use linear_algebra
   use equations, only: equation_runner
   use maths_constants, only: DiffOrder, nband, sub_diag, sup_diag
   use Newtons_method, only: non_linear_iteration

contains

!!
! @brief      {Subroutine that sets up a BVP solver for the equation given in definitions.f90}
!
! @return     Solves the BVP in the form A u_xx + B u_x + C = D. Boundary conditions can be general
!!
   Subroutine solve_runner

      integer :: i, j, k, iteration
      real(dp), dimension(:), allocatable :: U
      real(dp), dimension(:,:,:), allocatable :: U_2d
      real(dp), dimension(:, :), allocatable :: L !! banded form matrix
      real(dp), dimension(:), allocatable :: RHS !! right hand side of equation

      Write (6, *) 'size of x domain::: ', nx
      Write (6, *) 'order of the finite differences', DiffOrder
      Write (6, *) '... Building linear matrix'
      Write (6, *)

      !!! domain has been built in main
      !!! next step is to build the discretised operators
      !!! equation is in the form L * X = RHS

      !! Set the equation
      Call equation_runner(L, RHS)

      !! Solve the equation
      Call solver_banded_double_precision(idim, nband, sub_diag, sup_diag, L, RHS, U)

      !! If non-linear BVP then go to non-linear solver

      Select Case (Non_Linear_switch)
      Case (1)
         Write (6, *)
         Write (6, *) 'Non linear BVP'
         Call non_linear_iteration(L, RHS, U, iteration)
      Case Default
         Write (6, *)
         Write (6, *) 'Linear BVP'
      End Select

      deallocate (L, RHS)
      !! Print the results





      If (Domain_number==1) then

            Open (10, file='BVP.dat')

            Write (6, '(5(A20,x))') 'Physical domain', 'Comp Domain', 'Numerical Soln'!, 'Exact Soln', 'error'
            Write (10, '(5(A20,x))') 'Physical domain', 'Comp Domain', 'Numerical Soln'!, 'Exact Soln', 'error'

      Select Case (Eqn_number)

      Case (1)
         Do i = 1, nx
            Write (6, '(4(f20.14,1x),e20.10)') xdom(i), xcdom(i), U(i),xdom(i)**2.d0!, ex**xdom(i), abs(x(i)-ex**xdom(i))
            Write (10, '(4(f20.14,1x),e20.10)') xdom(i), xcdom(i), U(i)!, ex**xdom(i), abs(x(i)-ex**xdom(i))
         End Do
      Case (2)
         Do i = 1, nx
            Write (6, '(4(f20.14,1x),e20.10)') xdom(i), xcdom(i), U(2*i - 1), U(2*i)!,!!, ex**xdom(i), abs(x(i)-ex**xdom(i))
            Write (10, '(4(f20.14,1x),e20.10)') xdom(i), xcdom(i), U(2*i - 1), U(2*i)!!, ex**xdom(i), abs(x(i)-ex**xdom(i))
         End Do
      End Select

      Else if(Domain_number==2)then



      Select Case (Eqn_number)

      Case (1)


            Open (10, file='BVP1.dat')
            Open (11, file='BVPx.dat')
            Open (12, file='BVPy.dat')

      allocate(U_2d(1:nx,1:ny,1))
        ! Convert the 1D solution vector into a 2D solution matrix
        do j = 1, ny
          do i = 1, nx
            k = (j-1)*nx + i
            U_2d(i, j, 1) = U(k)
          end do
        end do

         Do i = 1, nx
            Write (11, '(10000(f20.14,1x),e20.10)') (xdom(i),j=1,ny)
            Write (10, '(10000(f20.14,1x),e20.10)') (U_2d(i,j,1),j=1,ny)
            Write (12, '(10000(f20.14,1x),e20.10)') (ydom(j),j=1,ny)
         End Do


         Case (2)
         Open (9, file='BVP1.dat')
         Open (10, file='BVP2.dat')
         Open (11, file='BVPx.dat')
         Open (12, file='BVPy.dat')

      allocate(U_2d(1:nx,1:ny,2))
        ! Convert the 1D solution vector into a 2D solution matrix
        do j = 1, ny
          do i = 1, nx
            k = (j-1)*nx + i
            U_2d(i, j, 1) = U(2*k-1)
            U_2d(i, j, 2) = U(2*k)
          end do
        end do

         Do i = 1, nx
            Write (11, '(10000(f20.14,1x),e20.10)') (xdom(i),j=1,ny)
            Write (9, '(10000(f20.14,1x),e20.10)') (U_2d(i,j,1),j=1,ny)
            Write (10, '(10000(f20.14,1x),e20.10)') (U_2d(i,j,2),j=1,ny)
            Write (12, '(10000(f20.14,1x),e20.10)') (ydom(j),j=1,ny)
         End Do
      End Select
      End if

      



      Close (10)
      Write (6, *)
      deallocate (U)

      Return

   End Subroutine solve_runner

End Module solve_bvp

