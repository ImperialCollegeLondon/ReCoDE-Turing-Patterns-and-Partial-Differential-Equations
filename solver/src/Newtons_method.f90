!!{Module that contains the non_linear_iteration subroutine}
!!
Module Newtons_method
   use type_kinds
   use reader, only: Time_switch, Non_Linear_switch, Newton_Error, Eqn_number, Max_iter, Domain_number
   use domain
   use linear_algebra
   use maths_constants, only: DiffOrder, nband, sub_diag, sup_diag

contains

!!
! @brief      {Non-linear extension to BVP and IVP - use Newton Iteration
!              PDE can be written as: Lu  + F(u) = RHS, F has non-linear terms
!
!              Will discuss in greater detail how this works
!
! @param      L          Linear opeartor (nband * idim)
! @param      RHS        RHS of the equation
! @param      U          Linear solution to L*U = RHS
!
! @return     U          Non-lienar Solution to L*U + F(u) = RHS
! @return     iteration  Number of iterations taken for convergence
!!
   Subroutine non_linear_iteration(L, RHS, U, iteration)
      use equations, only: non_linear_setup, non_linear_setup2D

      external :: DGBMV ! external banded matrix multiplication from Blas

      real(dp), dimension(:), allocatable, intent(inout) :: U !! Soln
      real(dp), dimension(:, :), allocatable, intent(in) :: L !! banded form matrix
      real(dp), dimension(:), allocatable, intent(in)  :: RHS !! right hand side of equation
      integer, intent(out) :: iteration

      real(dp), dimension(:, :), allocatable :: N_LHS !! Left hand side of equation
      real(dp), dimension(:), allocatable :: N_RHS !! Right hand side of equation
      real(dp), dimension(:), allocatable :: F !! Non linear terms
      real(dp), dimension(:, :), allocatable :: Fu, Fv! non-linear derivative terms
      real(dp), dimension(:), allocatable :: X ! current iteration solution
      integer :: i, j
      real(dp) :: error, blank

      iteration = 0

      Do iteration = 0, max_iter

         allocate (X(1:idim))
         allocate (F(1:idim))!, Fu(1:idim),Fv(1:idim))
         allocate (N_LHS(1:nband, 1:idim), N_RHS(1:idim))
         N_LHS = 0.d0; N_RHS = 0.d0; F = 0.d0

      !! Set the linear solution (when there are no non-linear terms) to current iteratin
         X = U

      !! Set the non-linear terms at the previous value
      !! F is a vector, Fu and Fv are banded matrices
         Select Case (Domain_number)
         Case (1)
            Call non_linear_setup(nx, xcdom, X, F, Fu, Fv)
         Case (2)
            !Write(6,*) 'here'
            Call non_linear_setup2D(nx, ny, xcdom, ycdom, idim, idim_xy, Eqn_number, nband, sub_diag, U, F, Fu, Fv)
         End Select

      !! If temporal marching, the non-linear terms need to multipled by -dt
         Select Case (Time_switch)
         Case (1)
            F = -F*dt
            Fu = -Fu*dt
            Fv = -Fv*dt
         Case (0)
         End Select

      !! Set up RHS of Newton Iteration
         !
      !! Banded matrix multplication
      !! First we find L * u (the solution at the previous iteration)
         Call DGBMV('N', idim, idim, sub_diag, sup_diag, 1.d0, L, nband, X, 1, 0.d0, N_RHS, 1)

         ! Do i = 1,idim
         ! Write(6,*) 'N_RHS', n_rhs(i)
         ! End do
         ! stop

      !! Then we want the RHS and non-linear terms

         N_RHS = RHS - N_RHS - F
         N_LHS = Fu + Fv + L

      !!! Solve for the newtown iteration error
         Call solver_banded_double_precision(idim, nband, sub_diag, sup_diag, N_LHS, N_RHS, X)

      !! calculate the error
         error = norm2(X)
         !  Write (6,*) 'iteration', iteration, error
      !! Add on the correction to the vector
         U = U + X

         deallocate (N_LHS, N_RHS, X, F, Fu, Fv)

         !Write(6,*) error, iteration, max_iter, Newton_Error

         ! if error is sufficiently small
         If (error .LT. Newton_Error) Then
            Return
         Else
         End if

         ! If iteration fails
         If (iteration == Max_iter) then
            !Do i = 1, idim
            !   Write (6, *) N_RHS(i), X(i), RHS(i), RHS(i) + N_RHS(i) + F(i), F(i)
            !End do
            Write (6, *) 'Newton iteratin reached maximum::', Max_iter
            Write (6, *) 'Last error', error
            Stop
         End if
      End do
   End Subroutine non_linear_iteration

End Module Newtons_method
