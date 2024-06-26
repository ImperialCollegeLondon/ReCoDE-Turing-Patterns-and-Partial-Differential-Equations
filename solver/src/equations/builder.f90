!!{This module sets up the subroutines so an equation can be discretised onto a streched grid
!  band_the_matrix - is a banded matrix calculator}
!!
Module equations_builder
   use type_kinds, only: dp
   use maths_constants, only: D1, D2, sub_diag, sup_diag, nband

contains

!!
! @brief      { Takes a square matrix A (n * n) and returns a banded matrix AB (nband * n)
!               This saves memory and computation time when a calculation is conducted with A}
!
!                We use the method from lapack:
!                On entry, the matrix A in band storage, in rows 1 to sub_diag+sup_diag+1.
!                The j-th column of A is stored in the j-th column of the array AB as follows:
!                AB(sup_diag+1+i-j,j) = A(i,j) for max(1,j-sup_diag)<=i<=min(n,j+sub_diag)
!
!
! @param      n      dimension of A
! @param      A      square matrix to band (n * n)
! @param      sub_diag     number of subdiagonals in A
! @param      sup_diag     number of superdiagonals in A
!
! @return     nband  the new leading dimension of AB
! @return     AB     the banded matrix, dimension (nband * n)
!!
   Subroutine band_the_matrix(n, A, AB)

      integer, intent(in) :: n
      real(dp), dimension(:, :), allocatable, intent(in) :: A
      real(dp), dimension(:, :), allocatable, intent(out) :: AB
      integer :: i, j

      allocate (AB(nband, n))

      Do j = 1, n
      Do i = max(1, j - sup_diag), min(n, j + sub_diag)
         AB(sup_diag + 1 + i - j, j) = A(i, j)
      End Do
      End Do

      Return
   End Subroutine band_the_matrix

!!
! @brief      {This Subroutine takes the coefficient of the differential equation
!              and scales them with the correct grid streching terms and step sizes.
!              The main scaling is the new B term: this has contribution from both ai and bi.
!              We additionally multiply out all the scalings
!              so the second order derivative coefficients remains unchanged.
!             }
!
! @param      ai         second order derivative coefficient
! @param      bi         first order derivative coefficient
! @param      ci         zeroth order derivative coefficient
! @param      di         right hand side coefficient
! @param      ch         computational stepsize
! @param      chsq       computational stepsize
! @param      metric1    The metric 1 (defined in domains)
! @param      metric1sq  The metric 1 squared (defined in domains)
! @param      metric2    The metric 2 (defined in domains)
!
! @return      A          second order derivative coefficient
! @return      B          first order derivative coefficient
! @return      C          zeroth order derivative coefficient
! @return      D          right hand side coefficient
!!
   Subroutine scales(A, B, C, D, ai, bi, ci, di, ch, chsq, metric1, metric1sq, metric2)
      real(dp), intent(in) :: ai, bi, ci, di ! entries from the equation
      real(dp), intent(in) :: ch, chsq ! computational distance and computational distance squared
      real(dp), intent(in) :: metric1, metric1sq, metric2 ! metric terms
      real(dp), intent(out) :: A, B, C, D ! outputs

!!!! this Subroutine evaluates the terms of the equations in the
!!!! computational Domain and corrects them with the relevent metrics

      A = ai
      B = (bi*metric1 - ai*metric2/metric1)*ch
      C = ci*metric1sq*chsq
      D = di*metric1sq*chsq

   End Subroutine scales

!!
! @brief      This Subroutine builds the discretised differential operator out of the sum
!             of the zeroth, first and second derivative operators
!
! @param      n     dimension of the problem
! @param      A     second order operator coefficient (vector)
! @param      B     first order operator coefficient (vector)
! @param      C     zeroth order operator coefficient (vector)
!
! @return      L    Full discretised operator in banded form
   Subroutine derivative_runner(n, A, B, C, L)
      integer, intent(in) :: n
      real(dp), dimension(:), allocatable, intent(in) :: A, B, C
      real(dp), dimension(:, :), allocatable, intent(inout) :: L
      real(dp), dimension(:, :), allocatable :: eD0, eD1, eD2

!!! this builds the entire derivative L in banded form
      allocate (eD0(1:nband, 1:n), eD1(1:nband, 1:n), eD2(1:nband, 1:n))

      ed0 = 0.d0; ed1 = 0.d0; ed2 = 0.d0

!! each subroutine takes the coefficient vector and multiplies it by the operator
!! L is then the sum of these matrices
      Call zero(n, ed0, C)
      L = ed0
      Call first(n, ed1, B)
      L = L + ed1
      Call second(n, ed2, A)
      L = L + ed2

      deallocate (ed0, ed1, ed2)

   End Subroutine derivative_runner

!!
! @brief      { This function takes the zeroth derivative coefficient and populates it into a diagonal banded matrix}
!
! @param      n      Dimension of the problem
! @param      C      The coefficient of zeroth derivative
!
! @return     Deriv  The output (nband * n) matrix
!!
   Subroutine zero(n, Deriv, C)
      real(dp), dimension(:, :), allocatable, intent(inout) :: Deriv
      real(dp), dimension(:), allocatable, intent(in) :: C
      integer, intent(in) :: n
      integer :: ii
      real(dp), dimension(:, :), allocatable :: Dtemp

!!! this sets up the zeroth order derivative

      Deriv = 0.d0

!! Dtemp is a temporary operator that will be turned into banded form
!! The matrix Deriv could be populated directly.
      allocate (Dtemp(1:n, 1:n)); Dtemp = 0.d0

      Do ii = 1, n
         Dtemp(ii, ii) = C(ii)
      End Do

      Call band_the_matrix(n, Dtemp, Deriv)
      deallocate (Dtemp)
   End Subroutine zero

!!
! @brief      { This function takes the first derivative coefficient and populates it into a diagonal banded matrix
!               This subroutine is very similar to first_difference in domains.f90. However, here we do not sum the rows.}
!
! @param      n      Dimension of the problem
! @param      B      The coefficient of first order derivative
!
! @return     Deriv  The output (nband * n) matrix
!!
   Subroutine first(n, Deriv, B)
      use omp_lib
      real(dp), dimension(:, :), allocatable, intent(inout) :: Deriv
      real(dp), dimension(:), allocatable, intent(in) :: B
      integer, intent(in) :: n
      integer :: i, j
      real(dp), dimension(:, :), allocatable :: Dtemp

!!! this sets up the first order derivative

      Deriv = 0.d0

      allocate (Dtemp(1:n, 1:n)); Dtemp = 0.d0

      Dtemp(1, 1:6) = D1(1, 1:6)*B(1)
      Dtemp(2, 1:6) = D1(2, 1:6)*B(2)

!!! Here we wish to parallelize the following: specifically the j loop. However it is difficult to do with array notation as
!   j is an index in both arrays of Dtemp. Therefore we make use of OpenMP library
!   $omp parallel do: Tells the compiler to parallelize the loop.
!   private(i): Declares i as a private variable, meaning each thread gets its own copy of i to avoid race conditions.
!   note we have included the omp_lib library

      ! Parallelize the outer loop
      !$omp parallel do private(i)
      do j = 3, n - 2
         do i = -2, 2
            Dtemp(j, j + i) = D1(3, 4 + i)*B(j)
         end do
      end do
      !$omp end parallel do

      Dtemp(n - 1, n - 5:n) = D1(4, 1:6)*B(n - 1)
      Dtemp(n, n - 5:n) = D1(5, 1:6)*B(n)

      Call band_the_matrix(n, Dtemp, Deriv)
      deallocate (Dtemp)

   End Subroutine first

!!
! @brief      { This function takes the first derivative coefficient and populates it into a diagonal banded matrix
!               This subroutine is very similar to second_difference in domains.f90. However, here we do not sum the rows.}
!
! @param      n      Dimension of the problem
! @param      A      The coefficient of second order derivative
!
! @return     Deriv  The output (nband * n) matrix
!!
   Subroutine second(n, Deriv, A)
      use omp_lib
      real(dp), dimension(:, :), allocatable, intent(inout) :: Deriv
      real(dp), dimension(:), allocatable, intent(in) :: A
      integer, intent(in) :: n
      integer :: i, j
      real(dp), dimension(:, :), allocatable :: Dtemp

!!! this sets up the second order derivative

      Deriv = 0.d0

      allocate (Dtemp(1:n, 1:n)); Dtemp = 0.d0

      Dtemp(1, 1:6) = D2(1, 1:6)*A(1)
      Dtemp(2, 1:6) = D2(2, 1:6)*A(2)

      ! Parallelize the outer loop
      !$omp parallel do private(i)
      do j = 3, n - 2
         do i = -2, 2
            Dtemp(j, j + i) = D2(3, 4 + i)*A(j)
         end do
      end do
      !$omp end parallel do

      Dtemp(n - 1, n - 5:n) = D2(4, 1:6)*A(n - 1)
      Dtemp(n, n - 5:n) = D2(5, 1:6)*A(n)

      Call band_the_matrix(n, Dtemp, Deriv)
      deallocate (Dtemp)

   End Subroutine second

End Module equations_builder
