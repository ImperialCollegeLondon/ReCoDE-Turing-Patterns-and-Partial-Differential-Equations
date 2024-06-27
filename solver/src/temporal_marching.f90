!!{Organises the temproal marching subroutines}
!!
Module temporal_marching
   use omp_lib
   use type_kinds
   use reader, only: Time_switch, Non_Linear_switch, Eqn_number
   use domain
   use equations, only: equation_runner, initial_condition
   use maths_constants, only: DiffOrder, nband, sub_diag, sup_diag
   use linear_algebra
   use Newtons_method, only: non_linear_iteration

   real(dp), dimension(:, :), allocatable :: Soln ! solution

contains

!!
! @brief      {Runs the temporal march}
!
! @return     Solves the IBVP in the form A u_xx + B u_x + C + non_linear terms = D u_t
!             Outputs to IVBP.dat
!!
   Subroutine march_runner
      integer :: i, j

      Write (6, *) 'size of x Domain::: ', nx, dxc
      Write (6, *) 'size of total Domain::: ', idim
      Write (6, *) 'size of t Domain::: ', nt, dt
      Write (6, *) 'order of the finite differences', DiffOrder
      Write (6, *) '... Building matrix'
      Write (6, *)
      Write (6, *) 'Starting march...'

      allocate (Soln(1:idim, 1:nt))

      Call initial_condition(nx, xcdom, Soln)

      Call implicit_march

      Select Case (Eqn_number)

      Case (1)

         Open (10, file='IVBP.dat')
         Write (10, '(20000(f20.14,1x))') 0.d0, (tDom(j), j=1, nt)

         Do i = 1, nx
            Write (10, '(20000(f20.14,1x))') xDom(i), (soln(i, j), j=1, nt)
         End Do
         Close (10)
         Write (6, *) 'March completed...'
         Write (6, *)

      Case (2)

         Open (10, file='IVBP_1.dat')
         Open (11, file='IVBP_2.dat')
         Write (10, '(20000(f20.14,1x))') 0.d0, (tDom(j), j=1, nt)
         Write (11, '(20000(f20.14,1x))') 0.d0, (tDom(j), j=1, nt)

         Do i = 1, nx
            Write (10, '(20000(f20.14,1x))') xDom(i), (soln(2*i - 1, j), j=1, nt)
            Write (11, '(20000(f20.14,1x))') xDom(i), (soln(2*i, j), j=1, nt)
         End Do
         Close (10)
         Close (11)
         Write (6, *) 'March completed...'
         Write (6, *)

      End Select
      Return
   End Subroutine march_runner

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
! @return     Soln - contains the soln to the PDE
!!
   Subroutine implicit_march

      real(dp), dimension(:), allocatable :: U, temp ! temporary vectors
      real(dp), dimension(:, :), allocatable :: L !! operator from equations
      real(dp), dimension(:, :), allocatable :: L_March !! implicit marching operator
      real(dp), dimension(:), allocatable :: RHS !! right hand side of equation

      integer :: i, j, iteration

!!!! Build the opeator
      allocate (temp(1:idim), U(idim))
      allocate (L_March(nband, 1:idim))

      !! Obtain the original operator L
      Call equation_runner(L, RHS)

      !! Builds the operator
      !!! Interior points - L and L_March are in banded form - row sub_diag+1 is where the diagonals live
      L_March(sub_diag + 1, :) = RHS(:)
      L_March = L_March - dt*L

      !! Set the boundaries - to the equation Lu = RHS
      !! Boundaries are rows 1 and N in non-banded form
      !! These loops find those rows in banded form
      !! This needs much more explaining/simplifying?
      !$omp Parallel Do
      Do i = 0, Eqn_number - 1
      Do j = 1, sub_diag, Eqn_number
         L_March(sub_diag + 2 - j, j + i) = L(sub_diag + 2 - j, j + i)
         L_March(sub_diag + j, idim + 1 - j - i) = L(sub_diag + j, idim + 1 - j - i)
      End Do
      End do
      !$omp End Parallel Do

      !!!! March the operator
      !!! Begin the march in time
      Do j = 2, nt

      !! set the boundaries
         temp(1:Eqn_number) = RHS(1:Eqn_number)
         temp(idim - Eqn_number:idim) = RHS(idim - Eqn_number:idim)

      !! set the interior points
         temp(1+Eqn_number:idim - Eqn_number) = RHS(1+Eqn_number:idim - Eqn_number)*Soln(1+Eqn_number:idim - Eqn_number, j - 1)

         !!! solve the system
         Call solver_banded_Double_precision(idim, nband, sub_diag, sup_diag, L_March, temp, U)

         Select Case (Non_Linear_switch)
         Case (1)
            Call non_linear_iteration(L_March, temp, U, iteration)
         Case Default
         End Select

         Soln(:, j) = U(:)
      End Do

      deallocate (temp, U, L_March)

   End Subroutine implicit_march

End Module temporal_marching

