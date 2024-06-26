!!{ This module sets up the discretised equation (including operator in banded form) and RHS }
!!
Module equations
   use omp_lib
   use type_kinds, only: dp
   use equations_builder
   !! here we allowing equations to read the equation definitions that we have set up:
   use equations_definition, only: equation1, equation1_BC_Top, equation1_BC_Bot
   use maths_constants, only: sub_diag, sup_diag, nband

!!! L and RHS form the equation
contains

  !!
   ! @brief      The subroutine sets up the equation in the form L * u = RHS
   !
   ! @param      n          dimension of problem
   ! @param      ch         step size
   ! @param      chsq       step size squared
   ! @param      cdom       computational domain
   ! @param      metric1    The metric 1 - determined by domains
   ! @param      metric1sq  The metric 1 sq - determined by domains
   ! @param      metric2    The metric 2 - determined by domains
   ! @param      which_equation    (integer - asks which equation you want to discretise)
   !
   ! @return     L     left hand side of the equation - returns a matrix in banded form (nband * n)
   ! @return     RHS   RHS of the equation - returns a vector (n)
  !!
   Subroutine build_the_matrix(n, ch, chsq, cdom, metric1, metric1sq, metric2, which_equation, L, RHS)

   !!! sets up a matrix equation in the form L * u = RHS
   !!! D is the RHS, Lx is the matrix

      integer, intent(in) :: n !! size of domain
      real(dp), intent(in) :: ch, chsq !computational distance/squared
      real(dp), dimension(:), allocatable, intent(in) :: cdom ! computational domain

    !! grid streching metric terms
      real(dp), dimension(:), allocatable, intent(in) :: metric1, metric1sq, metric2

    !! which equation?
      integer, intent(in) :: which_equation

    !! outputs
      real(dp), dimension(:, :), allocatable, intent(out) :: L !! banded form matrix
      real(dp), dimension(:), allocatable, intent(out) :: RHS !! right hand side of equation

      !!!! relevent vectors/infomation
      real(dp), dimension(:), allocatable :: A, B, C, D
      real(dp), dimension(:), allocatable :: At, Bt, Ct, Dt
      integer :: i, j

    !!! allocate
      allocate (L(nband, n), RHS(n))
      allocate (A(n), B(n), C(n), D(n))
      allocate (At(n), Bt(n), Ct(n), Dt(n))

    !! At,Bt,Ct and Dt form the coefficients of the discretised equation in the form Atu'' + Btu' + Ctu = Dt
    !! A,B,C and D are then the corrected coefficients given grid streching transfomations

      
    !! here we use OpenMP to parallelise the Do Loops - no infomation is shared between loops

      Select Case (which_equation)

      !!! Equation 1
      Case (1)

         !!! Boundaries
         Call equation1_BC_Bot(cdom(1), At(1), Bt(1), Ct(1), Dt(1))
         Call equation1_BC_Top(cdom(n), At(n), Bt(n), Ct(n), Dt(n))
         !!! Interior

         !$omp Parallel Do
           Do i = 2, n - 1
              Call equation1(cdom(i), At(i), Bt(i), Ct(i), Dt(i))
           End do
         !$omp End Parallel Do
      Case Default

         Write (6, *) 'Equation Error: which_equation in equations.f90 should be 1'

      End Select

   
    !! Apply the correct scallings
      !$omp Parallel Do
        Do i = 1, n 
           Call scales(A(i), B(i), C(i), D(i), At(i), Bt(i), Ct(i), Dt(i), ch, chsq, metric1(i), metric1sq(i), metric2(i))
        End Do
      !$omp End Parallel Do

    !!! Derivative runner moves the coefficients into a banded matrix L
      Call derivative_runner(n, A, B, C, L)
    !! set the fina output for RHS
      RHS = D

      deallocate (Dt, Ct, Bt, At, D, C, B, A)

      Return
   End Subroutine build_the_matrix

End Module equations

