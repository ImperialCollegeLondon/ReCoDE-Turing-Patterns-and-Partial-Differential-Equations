!!{The module that set up and runs the BVP Solver for both non-linear BVP and Linear BVP
!  Non-linear uses Newton iteration}
!!
Module solve_bvp
   use type_kinds
   use reader, only : Time_switch, Non_Linear_switch
   use domain
   use linear_algebra
   use equations, only : build_the_matrix
   use maths_constants, only : DiffOrder, nband, sub_diag, sup_diag

contains

!!
! @brief      {Subroutine that sets up a BVP solver for the equation given in definitions.f90}
!
! @return     Solves the BVP in the form A u_xx + B u_x + C = D. Boundary conditions can be general
!!
   Subroutine solve_runner

      integer :: i, j
      real(dp), dimension(:), allocatable :: U
      real(dp), dimension(:, :), allocatable :: L !! banded form matrix
      real(dp), dimension(:), allocatable :: RHS !! right hand side of equation

      Write (6, *) 'size of x domain::: ', nx
      Write (6, *) 'order of the finite differences', DiffOrder
      Write (6, *) '... Building linear matrix'
      Write (6, *)

      !!! domain has been built in main
      !!! next step is to build the discretised operators
      !!! equation is in the form L * X = RHS

      Call build_the_matrix(nx, dxc, dxcsq, xcdom, xmetric1, xmetric1sq, xmetric2, 1, L, RHS)

      !! Solve the equation

      Call solver_banded_double_precision(nx, nband, sub_diag, sup_diag, L, RHS, U)
      
      !! If non-linear BVP then go to non-linear solver

      Select Case (Non_Linear_switch)
      Case (1)
         Write (6, *)
         Write (6, *) 'Non linear BVP'
         Call non_linear_BVP(L,RHS,U)
      Case Default
         Write (6, *)
         Write (6, *) 'Non linear BVP'
      End Select


      deallocate (L, RHS)
      !! Print the results

      Open (10, file='BVP.dat')

      Write (6, '(5(A20,x))') 'Physical domain', 'Comp Domain', 'Numerical Soln'!, 'Exact Soln', 'error'
      Write (10, '(5(A20,x))') 'Physical domain', 'Comp Domain', 'Numerical Soln'!, 'Exact Soln', 'error'
      Do i = 1, nx
         Write (6, '(4(f20.14,1x),e20.10)') xdom(i), xcdom(i), U(i)!, ex**xdom(i), abs(x(i)-ex**xdom(i))
         Write (10, '(4(f20.14,1x),e20.10)') xdom(i), xcdom(i), U(i)!, ex**xdom(i), abs(x(i)-ex**xdom(i))
      End Do
      Close (10)
      Write (6, *)
      deallocate (U)

      Return

   End Subroutine solve_runner

!!
! @brief      {Non-linear extension to BVP - use Newton Iteration
!              PDE can be written as: Lu  + F(u) = RHS, F has non-linear terms
!              
!              Newton iteration, given previous guess u1 
!              Solve for v: Lv + F'(u1)v = Lu1  + F(u1) - RHS
!              Then u2 = u1 + v
!
! @param      L     Linear opeartor (nband * nx)
! @param      RHS   RHS of the equation
! @param      U     Linear solution to L*U = RHS
!
! @return     U     Non-lienar Solution to L*U + F(u) = RHS
!!
   Subroutine non_linear_BVP(L,RHS,U)
      use equations, only : equation1_non_linear
      use domain, only : first_difference
      external :: DGBMV ! external banded matrix multiplication from Blas

      real(dp), dimension(:), allocatable ,intent(inout) :: U
      real(dp), dimension(:, :), allocatable ,intent(in) :: L !! banded form matrix
      real(dp), dimension(:), allocatable,intent(in)  :: RHS !! right hand side of equation
      real(dp), dimension(:,:), allocatable :: N_LHS !!
      real(dp), dimension(:), allocatable :: N_RHS !! 
      real(dp), dimension(:), allocatable :: F, Fu ! non-linear term and derivative
      real(dp), dimension(:), allocatable :: X ! current iteration
      integer :: iteration ! iteration level
      integer :: i,j
      real(dp) :: error

      allocate(X(1:nx),F(1:nx),Fu(1:nx))
      allocate(N_LHS(1:nband,1:nx),N_RHS(1:nx))
      N_LHS = 0.d0

      iteration = 0

      Do 
      !! Set the linear solution (when there are no non-linear terms) to current iteratin
      X(:) = U(:)

     
      !! Set the non-linear terms at the previous value
      !$omp Parallel Do
      Do i = 2, nx-1 
         Call equation1_non_linear(xcdom(i), U(i), F(i), Fu(i))
      End do
      !$omp End Parallel Do

      !! Set up RHS of Newton Iteration 
      !! Banded matrix multplication
      !! First we find L * u (the solution at the previous iteration)
      Call DGBMV('N', nx, nx, sub_diag, sup_diag, 1.d0, L, nband, X, 1, 0.d0, N_RHS, 1)
      !! Then we want the RHS and non-linear terms
      N_RHS = N_RHS - RHS + F
      N_RHS = - N_RHS

      !! Set up LHS of Newton Iteration
      ! The diagonal of N_LHS needs to be derivative of non-linear terms (diagonal) + L (note banded form) 
      ! Diagonal of a banded matrix is along the sub_diag + 1 row 
      N_LHS(sub_diag + 1, :) = Fu(:) 
      N_LHS = N_LHS + L 

     ! do i = 1,nx
     ! WRITE(6,*) Fu(i),u(i)
     ! end do
     ! WRITE(6,*) 
     !   do i = 1,nband
     ! WRITE(6,'(12f16.8)') (n_LHS(i,j),j=1,nx)
     ! end do

      !!! Solve for the newtown iteration error
      Call solver_banded_double_precision(nx, nband, sub_diag, sup_diag, N_LHS, N_RHS, U)
   
      !! calculate the error
      error = norm2(U)
      iteration = iteration+1
      
      !Write(6,'(i10,13(f15.6))') iteration,(U(j),j=1,nx),error

      ! reset the solution
      U = U + X
      
      ! if error is sufficiently small
      If (error.gt.1.d-5)then
      Else
         Goto 20
      End if

      If (iteration<100)then
      Else
         Write (6, *) 'Newton iteratin for non-linear BVP passed 100, current error::',error
         Write (6, *) 'Stopping'
         Stop 
      End If
      End do

      20 Continue
      Write (6, *) 'non-linear BVP completed, iterations::', iteration, 'error::', error
      Write (6, *)

      deallocate(N_LHS, N_RHS, X, F, Fu)
      Return
   End Subroutine non_linear_BVP

End Module solve_bvp













