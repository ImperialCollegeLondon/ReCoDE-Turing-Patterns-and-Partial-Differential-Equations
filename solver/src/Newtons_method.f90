!!{Module that contains the non_linear_iteratio}
!!
Module Newtons_method
  use type_kinds
  use reader, only: Time_switch, Non_Linear_switch, Newton_Error
  use equations, only: equation1_non_linear
  use domain
  use linear_algebra
  use maths_constants, only: DiffOrder, nband, sub_diag, sup_diag

  contains

!!
! @brief      {Non-linear extension to BVP and IVP - use Newton Iteration
!              PDE can be written as: Lu  + F(u) = RHS, F has non-linear terms
!
!              Will discuss in greater detail how this works... TODO
!
! @param      L          Linear opeartor (nband * nx)
! @param      RHS        RHS of the equation
! @param      U          Linear solution to L*U = RHS
!
! @return     U          Non-lienar Solution to L*U + F(u) = RHS
! @return     iteration  Number of iterations taken for convergence
!!
 Subroutine non_linear_iteration(L, RHS, U, iteration)
   use equations, only: equation1_non_linear

      external :: DGBMV ! external banded matrix multiplication from Blas

      real(dp), dimension(:), allocatable, intent(inout) :: U
      real(dp), dimension(:, :), allocatable, intent(in) :: L !! banded form matrix
      real(dp), dimension(:), allocatable, intent(in)  :: RHS !! right hand side of equation
      real(dp), dimension(:, :), allocatable :: N_LHS !!
      real(dp), dimension(:), allocatable :: N_RHS !!
      real(dp), dimension(:), allocatable :: F, Fu ! non-linear term and derivative
      real(dp), dimension(:), allocatable :: X ! current iteration
      integer :: iteration ! iteration level
      integer :: i, j
      real(dp) :: error

      iteration = 0

      Do
         allocate (X(1:nx), F(1:nx), Fu(1:nx))
         allocate (N_LHS(1:nband, 1:nx), N_RHS(1:nx))
         N_LHS = 0.d0; N_RHS = 0.d0


      !! Set the linear solution (when there are no non-linear terms) to current iteratin
         X(:) = U(:)

      !! Set the non-linear terms at the previous value
         !$omp Parallel Do
         Do i = 2, nx - 1
            Call equation1_non_linear(xcdom(i), U(i), F(i), Fu(i))
         End do
         !$omp End Parallel Do

         Select Case (Time_switch)
         Case (1)
           F = - F * dt
           Fu = - Fu * dt
         Case (0)
         End Select


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
         iteration = iteration + 1

         !Write(6,'(i10,13(f15.6))') iteration,(U(j),j=1,5),error

         ! reset the solution
         U = U + X

         ! if error is sufficiently small
         If (error .gt. Newton_Error) then
         Else
            Goto 20
         End if

         If (iteration < 100) then
         Else
            Write (6, *) 'Newton iteratin for non-linear BVP passed 100, current error::', error
            Write (6, *) 'Stopping'
            Stop
         End If

         deallocate (N_LHS, N_RHS, X, F, Fu)
      End do

20    Continue
!      Write (6, *) 'non-linear BVP completed, iterations::', iteration, 'error::', error
 !     Write (6, *)
      deallocate (N_LHS, N_RHS, X, F, Fu)
      Return
   End Subroutine non_linear_iteration

End Module Newtons_method