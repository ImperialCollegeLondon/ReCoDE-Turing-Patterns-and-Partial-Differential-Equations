!!{Organises the temproal marching subroutines}
!!
Module solve_ibvp
   use type_kinds
   use reader, only: Time_switch, Non_Linear_switch, Eqn_number, Newton_Error
   use domain
   use equations, only: equation_runner, initial_condition, initial_condition2D
   use maths_constants, only: DiffOrder, nband, sub_diag, sup_diag
   use linear_algebra
   use Newtons_method, only: non_linear_iteration

   real(dp), dimension(:, :), allocatable :: Soln ! solution

contains

!!
! @brief      {Runs the temporal march}
! 
! @input      L    - LHS of equation, banded matrix form. Dimension (nband * idim)
! @input      RHS  - RHS of equation, vector. Dimension (idim)
!
! @return     Solves the IBVP in the form A u_xx + B u_x + C + non_linear terms = D u_t
!             Outputs to IVBP.dat
!!
   Subroutine ibvp_runner(L,RHS)
      real(dp), dimension(:, :), allocatable,intent(inout) :: L !! left hand side
      real(dp), dimension(:), allocatable,intent(inout) :: RHS !! right hand side of equation

      integer :: i, j, jk, k
      real(dp), dimension(:, :, :), allocatable :: U_2d
     

      allocate (Soln(1:idim, 1:nt))

      !! Set the initial conditions

      Select Case (Domain_number)
      Case (1)
         Call initial_condition(nx, xdom, Soln)
      Case (2)
         Call initial_condition2D(nx, ny, xdom, ydom, idim, idim_xy, Eqn_number, Soln)
      End Select

      !! March in time

      Call implicit_march(L, RHS)


      !!!! Print the solution
      If (Domain_number == 1) then

      !! Printing solution
         Select Case (Eqn_number)

         Case (1)

            Open (10, file='IVBP_1eqn_1D.dat')
            Write (10, '(2000000(f20.14,1x))') 0.d0, (tDom(j), j=1, nt)

            Do i = 1, nx
               Write (10, '(2000000(f20.14,1x))') xDom(i), (soln(i, j), j=1, nt)
            End Do
            Close (10)
            Write (6, *) 'March completed...'
            Write (6, *)

         Case (2)

            Open (10, file='IVBP1_2eqn_1D.dat')
            Open (11, file='IVBP2_2eqn_1D.dat')
            Write (10, '(20000000(f20.14,1x))') 0.d0, (tDom(j), j=1, nt)
            Write (11, '(20000000(f20.14,1x))') 0.d0, (tDom(j), j=1, nt)

            Do i = 1, nx
               Write (10, '(20000000(f20.14,1x))') xDom(i), (soln(2*i - 1, j), j=1, nt)
               Write (11, '(20000000(f20.14,1x))') xDom(i), (soln(2*i, j), j=1, nt)
            End Do
            Close (10)
            Close (11)
            Write (6, *) 'March completed...'
            Write (6, *)

         End Select

      Else If (Domain_Number == 2) then

         Select Case (Eqn_number)

         Case (1)
            Open (10, file='IBVP_1eqn_2D.dat')
            Open (11, file='IBVPx_1eqn_2D.dat')
            Open (12, file='IBVPy_1eqn_2D.dat')

            allocate (U_2d(1:nx, 1:ny, 1))

            Do i = 1, nx
               Write (11, '(20000000(f20.14,1x),e20.10)') (xdom(i), j=1, ny)
               Write (12, '(20000000(f20.14,1x),e20.10)') (ydom(j), j=1, ny)
            End Do

            Do jk = 1, nt

               ! Convert the 1D solution vector into a 2D solution matrix
               Do j = 1, ny
                  Do i = 1, nx
                     k = (j - 1)*nx + i
                     U_2d(i, j, 1) = Soln(k, jk)
                  end Do
               end Do

               Do i = 1, nx
                  Write (10, '(20000000(f20.14,1x),e20.10)') (U_2d(i, j, 1), j=1, ny)
               End Do

            end Do

         Case (2)
            Open (9, file='IBVP1_2eqn_2D.dat')
            Open (10, file='IBVP2_2eqn_2D.dat')
            Open (11, file='IBVPx_2eqn_2D.dat')
            Open (12, file='IBVPy_2eqn_2D.dat')

            allocate (U_2d(1:nx, 1:ny, 2))

            Do i = 1, nx
               Write (11, '(20000000(f20.14,1x),e20.10)') (xdom(i), j=1, ny)
               Write (12, '(20000000(f20.14,1x),e20.10)') (ydom(j), j=1, ny)
            End Do

            Do jk = 1, nt

               ! Convert the 1D solution vector into a 2D solution matrix
               Do j = 1, ny
                  Do i = 1, nx
                     k = (j - 1)*nx + i
                     U_2d(i, j, 1) = Soln(2*k - 1, jk)
                     U_2d(i, j, 2) = Soln(2*k, jk)
                  end Do
               end Do

               Do i = 1, nx
                  Write (9, '(20000000(f20.14,1x),e20.10)') (U_2d(i, j, 1), j=1, ny)
                  Write (10, '(20000000(f20.14,1x),e20.10)') (U_2d(i, j, 2), j=1, ny)
               End Do

            end Do

         End Select
      End If
      Return
   End Subroutine ibvp_runner

!!
! @brief      {Sets up and solves the implicit system
!              We want solve the system
!              RHS(u^{i+1} - u^{i}) = Lu^{i+1}  + boundary conditions
!
!              Outside of the boundaries:::
!              Ltemp = (RHS - dt*L)
!              Ltemp u^{i+1} = RHS u^{i}
!
!              At the boundaries we have Lu = RHS (independent of spatial march)
!
!
!              Ltemp   Implicit marching operator
!              temp    temporary RHS
!              X       temporary solution
!
!              }
!
! @input      L    - LHS of equation, banded matrix form. Dimension (nband * idim)
! @input      RHS  - RHS of equation, vector. Dimension (idim)
! 
!
! @return     Soln - contains the soln to the PDE
!!
   Subroutine implicit_march(L,RHS)
      real(dp), dimension(:, :), allocatable,intent(inout) :: L !! left hand side
      real(dp), dimension(:), allocatable,intent(inout) :: RHS !! right hand side of equation
      
      real(dp), dimension(:), allocatable :: U, temp ! temporary vectors
      !! operator from equations
      real(dp), dimension(:, :), allocatable :: L_March, L_march_unbanded !! implicit marching operator

      integer :: i, j, k, differ, iteration



     
!!!! Build the opeator
      allocate (temp(1:idim), U(idim))
      allocate (L_March(nband, 1:idim), L_march_unbanded(1:idim, 1:idim))


      !! Builds the operator
      !!! Interior points - L and L_March are in banded form - row sub_diag+1 is where the diagonals live

      Do i = 1, idim
         L_march_unbanded(i, i) = RHS(i)
      End Do
      Call band_the_matrix(idim, L_march_unbanded, sub_diag, sup_diag, nband, L_March)
      Deallocate (L_march_unbanded)

      !! Set the marching operator - differs from L !!!!
      L_March = L_March - dt*L

      !!!! March the operator
      !!! Begin the march in time
      Do j = 2, nt

      !!!! Set the boundaries
         Select Case (Domain_number)
         Case (1)

      !! set the boundaries

            temp(1:Eqn_number) = RHS(1:Eqn_number)
            temp(idim - Eqn_number:idim) = RHS(idim - Eqn_number:idim)

      !! set the interior points
            temp(1 + Eqn_number:idim - Eqn_number) = &
                      &RHS(1 + Eqn_number:idim - Eqn_number)*Soln(1 + Eqn_number:idim - Eqn_number, j - 1)

         Case (2)

      !! set the interior points
            temp(:) = RHS(:)*Soln(:, j - 1)

         !!! Set the boundary points
            differ = Eqn_number*nx
            temp(1:differ) = RHS(1:differ)
            temp(((ny - 1)*nx)*Eqn_number:idim) = RHS(((ny - 1)*nx)*Eqn_number:idim)


         !!! Here it's quite complicated to set the correct boundary conditions to temp 
            If (Eqn_number == 1) then

               Do k = 1, ny - 1
                  temp(1 + k*nx) = RHS(1 + k*nx)
                  temp((nx) + k*nx) = RHS((nx) + k*nx)
               End Do

            Else if (Eqn_number == 2) then

            !! This can probably be simplifed however it's a little complicated. Easy implementations over speed!
               Do k = 1, ny - 1
                  temp(1 + differ*k) = RHS(1 + differ*k)
                  temp(2 + differ*k) = RHS(2 + differ*k)
                  temp(2*nx - 1 + k*differ) = RHS(2*nx - 1 + k*differ)
                  temp(2*nx + k*differ) = RHS(2*nx + k*differ)
               End Do

            End If

         End Select

         !!! solve the system
         Call solver_banded_Double_precision(idim, nband, sub_diag, sup_diag, L_March, temp, U)

         !!! Non-linear solve
         Select Case (Non_Linear_switch)
         Case (1)
            Call non_linear_iteration(L_March, temp, U, iteration)
         End Select

         ! Set the solution
         Soln(:, j) = U(:)
         
         If (mod(j,100)==0) Then
            Write (6, *) 'Position done::', j,'/',nt
         End if
      End Do

      deallocate (temp, U, L_March)

   End Subroutine implicit_march

End Module solve_ibvp

