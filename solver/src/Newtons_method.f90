!!{Module that contains the non_linear_iteration}
!!
Module Newtons_method
   use type_kinds
   use reader, only: Time_switch, Non_Linear_switch, Newton_Error, Eqn_number
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
! @param      L          Linear opeartor (nband * idim)
! @param      RHS        RHS of the equation
! @param      U          Linear solution to L*U = RHS
!
! @return     U          Non-lienar Solution to L*U + F(u) = RHS
! @return     iteration  Number of iterations taken for convergence
!!
   Subroutine non_linear_iteration(L, RHS, U, iteration)
      use equations, only: non_linear_setup

      external :: DGBMV ! external banded matrix multiplication from Blas

      real(dp), dimension(:), allocatable, intent(inout) :: U !! Soln
      real(dp), dimension(:, :), allocatable, intent(in) :: L !! banded form matrix
      real(dp), dimension(:), allocatable, intent(in)  :: RHS !! right hand side of equation
      real(dp), dimension(:, :), allocatable :: N_LHS !!
      real(dp), dimension(:), allocatable :: N_RHS !!
      real(dp), dimension(:), allocatable :: F,Fu,Fv! non-linear term and derivative
      real(dp), dimension(:), allocatable :: X ! current iteration
      integer :: iteration ! iteration level
      integer :: i, j
      real(dp) :: error, blank

      iteration = 0

      Do

         allocate (X(1:idim))
         allocate( F(1:idim), Fu(1:idim),Fv(1:idim))
         allocate (N_LHS(1:nband, 1:idim), N_RHS(1:idim))
         N_LHS = 0.d0; N_RHS = 0.d0
         F = 0.d0; Fu = 0.d0; Fv = 0.d0

      !! Set the linear solution (when there are no non-linear terms) to current iteratin
         X(:) = U(:)



      !! Set the non-linear terms at the previous value
         Call non_linear_setup(nx, xcdom, U, F, Fu, Fv)

         Select Case (Time_switch)
         Case (1)
            F = -F*dt
            Fu = -Fu*dt
            Fv = -Fv*dt
         Case (0)
         End Select

      !! Set up RHS of Newton Iteration
      !! Banded matrix multplication
      !! First we find L * u (the solution at the previous iteration)
         Call DGBMV('N', idim, idim, sub_diag, sup_diag, 1.d0, L, nband, X, 1, 0.d0, N_RHS, 1)
      !! Then we want the RHS and non-linear terms


      !! Set up LHS of Newton Iteration
         ! The diagonal of N_LHS needs to be derivative of non-linear terms (diagonal) + L (note banded form)
         ! Diagonal of a banded matrix is along the sub_diag + 1 row
         
         Select Case (Eqn_number)
         
         Case(1)
            N_RHS = N_RHS - RHS + F
            N_RHS = -N_RHS
            N_LHS(sub_diag + 1, :) = Fu(:)

         Case(2)

         Do i = 1,nx
            N_RHS(2*i-1) = N_RHS(2*i-1) - RHS(2*i-1) + F(2*i-1)
            N_RHS(2*i) = N_RHS(2*i) - RHS(2*i) + F(2*i)

            N_LHS(sub_diag + 1, 2*i-1) = Fu(2*i-1)
            N_LHS(sub_diag + 1, 2*i) = Fv(2*i)
         End Do   

         N_RHS = -N_RHS

         End Select
         
         N_LHS = N_LHS + L

         ! do i = 1,idim
         ! WRITE(6,*) Fu(i),u(i)
         ! end do
         ! WRITE(6,*)
         !   do i = 1,nband
         ! WRITE(6,'(12f16.8)') (n_LHS(i,j),j=1,idim)
         ! end do

      !!! Solve for the newtown iteration error
         Call solver_banded_double_precision(idim, nband, sub_diag, sup_diag, N_LHS, N_RHS, U)

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

         deallocate (N_LHS, N_RHS, X, F, Fu,Fv)
      End do

20    Continue
!      Write (6, *) 'non-linear BVP completed, iterations::', iteration, 'error::', error
      !     Write (6, *)
      deallocate (N_LHS, N_RHS, X, F, Fu)
      Return
   End Subroutine non_linear_iteration

End Module Newtons_method
