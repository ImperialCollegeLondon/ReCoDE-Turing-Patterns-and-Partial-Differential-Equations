!!{Organises the temproal marching subroutines}
!!
Module temporal_marching
   use omp_lib
   use type_kinds
   use reader
   use Domain
   use equations, only : build_the_matrix
   use maths_constants, only : DiffOrder, nband, sub_diag, sup_diag
   use linear_algebra

   real(dp), dimension(:), allocatable :: X, temp ! temporary vectors
   real(dp), dimension(:, :), allocatable :: Soln ! solution
   


   real(dp), dimension(:, :), allocatable :: L !! operator from equations
   real(dp), dimension(:, :), allocatable :: L_March !! implicit marching operator
   real(dp), dimension(:), allocatable :: RHS !! right hand side of equation

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
      Write (6, *) 'size of t Domain::: ', nt, dt
      Write (6, *) 'order of the finite differences', DiffOrder
      Write (6, *) '... Building matrix'
      Write (6, *)
      Write (6, *) 'Starting march...'

      allocate (Soln(1:nx, 1:nt))

      Call linear_implicit_march

      open (10, file='IVBP.dat')
      Write (10, '(20000(f20.14,1x))') 0.d0, (tDom(j), j=1, nt)

      Do i = 1, nx
         Write (10, '(20000(f20.14,1x))') xDom(i), (soln(i, j), j=1, nt)
      End Do

      Write (6, *) 'March completed...'
      Write (6, *)
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
   Subroutine linear_implicit_march
      integer :: i, j

      allocate (temp(1:nx), X(nx))
      allocate (L_March(nband, 1:nx))

      !! Obtain the original operator L
      Call build_the_matrix(nx, dxc, dxcsq, xcdom, xmetric1, xmetric1sq, xmetric2, 1, L, RHS)

      !! Builds the operator
      !!! Interior points - L and L_March are in banded form - row sub_diag+1 is where the diagonals live
      L_March(sub_diag + 1, :) = RHS(:)
      L_March = L_March - dt*L

      !! Set the boundaries - to the equation Lu = RHS
      !! Boundaries are rows 1 and N in non-banded form
      !! These loops find those rows in banded form

      !$omp Parallel Do
         Do j = 1, sub_diag
            L_March(sub_diag + 2 - j, j) = L(sub_diag + 2 - j, j)
            L_March(sub_diag + j, nx + 1 - j) = L(sub_diag + j, nx + 1 - j)
         End Do
      !$omp End Parallel Do

      !!! Begin the march in time
      Do j = 2, nt

      !! set the boundaries

         temp(1) = RHS(1) 
         temp(nx) = RHS(nx)

      !! set the interior points
         temp(2:nx - 1) = RHS(2:nx - 1)*Soln(2:nx - 1, j - 1)

         !!! solve the system
         Call solver_banded_Double_precision(nx, nband, sub_diag, sup_diag, L_March, temp, X)

         Soln(:, j) = X(:)
      End Do

      deallocate (temp, X, L_March)

   End Subroutine linear_implicit_march

End Module temporal_marching
