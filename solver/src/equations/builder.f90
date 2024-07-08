!!{This module sets up the subroutines so an equation can be discretised onto a streched grid}
!!
Module equations_builder
   use type_kinds, only: dp
   use maths_constants, only: D1, D2

contains

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
! @param      metric1    The metric 1 (defined in Domains)
! @param      metric1sq  The metric 1 squared (defined in Domains)
! @param      metric2    The metric 2 (defined in Domains)
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

      A = ai/(metric1sq*chsq)
      B = ((bi*metric1 - ai*metric2/metric1)*ch)/(metric1sq*chsq)
      C = ci
      D = di

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
      allocate (eD0(1:n, 1:n), eD1(1:n, 1:n), eD2(1:n, 1:n))

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

!!! this sets up the zeroth order derivative

      Deriv = 0.d0

      Do ii = 1, n
         Deriv(ii, ii) = C(ii)
      End Do

   End Subroutine zero

!!
! @brief      { This function takes the first derivative coefficient and populates it into a diagonal banded matrix
!               This subroutine is very similar to first_difference in Domains.f90. However, here we Do not sum the rows.}
!
! @param      n      Dimension of the problem
! @param      B      The coefficient of first order derivative
!
! @return     Deriv  The output (nband * n) matrix
!!
   Subroutine first(n, Deriv, B)
      real(dp), dimension(:, :), allocatable, intent(inout) :: Deriv
      real(dp), dimension(:), allocatable, intent(in) :: B
      integer, intent(in) :: n
      integer :: i, j

!!! this sets up the first order derivative

      Deriv = 0.d0

      Deriv(1, 1:6) = D1(1, 1:6)*B(1)
      Deriv(2, 1:6) = D1(2, 1:6)*B(2)

      Do j = 3, n - 2
         Do i = -2, 2
            Deriv(j, j + i) = D1(3, 4 + i)*B(j)
         End Do
      End Do

      Deriv(n - 1, n - 5:n) = D1(4, 1:6)*B(n - 1)
      Deriv(n, n - 5:n) = D1(5, 1:6)*B(n)

   End Subroutine first

!!
! @brief      { This function takes the first derivative coefficient and populates it into a diagonal banded matrix
!               This subroutine is very similar to second_difference in Domains.f90. However, here we Do not sum the rows.}
!
! @param      n      Dimension of the problem
! @param      A      The coefficient of second order derivative
!
! @return     Deriv  The output (nband * n) matrix
!!
   Subroutine second(n, Deriv, A)
      real(dp), dimension(:, :), allocatable, intent(inout) :: Deriv
      real(dp), dimension(:), allocatable, intent(in) :: A
      integer, intent(in) :: n
      integer :: i, j

!!! this sets up the second order derivative

      Deriv = 0.d0

      Deriv(1, 1:6) = D2(1, 1:6)*A(1)
      Deriv(2, 1:6) = D2(2, 1:6)*A(2)

      Do j = 3, n - 2
         Do i = -2, 2
            Deriv(j, j + i) = D2(3, 4 + i)*A(j)
         End Do
      End Do

      Deriv(n - 1, n - 5:n) = D2(4, 1:6)*A(n - 1)
      Deriv(n, n - 5:n) = D2(5, 1:6)*A(n)

   End Subroutine second

End Module equations_builder
