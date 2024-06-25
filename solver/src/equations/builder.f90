module equations_builder
use type_kinds, only : dp
use maths_constants, only : D1,D2,sub_diag, sup_diag, nband

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
! @brief      {This subroutine takes the coefficient of the differential equation 
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
! @brief      This subroutine builds the discretised differential operator out of the sum 
!             of the zeroth, first and second derivative operators
!
! @param      n     { parameter_description }
! @param      A     { parameter_description }
! @param      B     { parameter_description }
! @param      C     { parameter_description }
! @param      L     { parameter_description }
!
! @return     { description_of_the_return_value }
subroutine derivative_runner(n,A,B,C,L)
integer,intent(in) :: n
	real(dp),dimension(:),allocatable,intent(in) :: A,B,C
	real(dp),dimension(:,:),allocatable,intent(inout) :: L
	real(dp),dimension(:,:),allocatable :: eD0,eD1,eD2

!!! this builds the entire derivative L

allocate(eD0(1:nband,1:n),eD1(1:nband,1:n),eD2(1:nband,1:n))

ed0 = 0.d0;ed1=0.d0;ed2=0.d0
	
	Call zero(n,ed0,C)
		L=ed0
	Call first(n,ed1,B)
		L=L+ed1
	Call second(n,ed2,A)
		L=L+ed2

deallocate(ed0,ed1,ed2)

end subroutine derivative_runner

subroutine zero(n,Deriv,C)
	real(dp), dimension(:,:), allocatable,intent(inout) :: Deriv
	real(dp),dimension(:),allocatable,intent(in) :: C
	integer,intent(in) :: n
	integer :: ii
	real(dp),dimension(:,:),allocatable :: Dtemp

!!! this sets up the zeroth order derivative

Deriv=0.d0

allocate(Dtemp(1:n,1:n)); Dtemp = 0.d0

Do ii=1,n
	Dtemp(ii,ii) = C(ii)
End do

Call band_the_matrix(n, Dtemp, Deriv)
deallocate(Dtemp)
end subroutine zero


subroutine first(n,Deriv,B)
	real(dp), dimension(:,:), allocatable,intent(inout) :: Deriv
	real(dp),dimension(:),allocatable,intent(in) :: B
	Integer,intent(in) :: n
	integer :: i,jj	
	real(dp),dimension(:,:),allocatable :: Dtemp

!!! this sets up the first order derivative

Deriv=0.d0

allocate(Dtemp(1:n,1:n)); Dtemp = 0.d0

	Dtemp(1,1:6) = D1(1,1:6)*B(1)
	Dtemp(2,1:6) = D1(2,1:6)*B(2)

Do jj = 3,n-2
Do i=-2,2
	Dtemp(jj,jj+i) = D1(3,4+i)*B(jj)
End do
END DO

!Do i = 0,-5,-1
	Dtemp(n-1,n-5:n) = D1(4,1:6)*B(n-1)
	Dtemp(n,n-5:n) = D1(5,1:6)*B(n)
!End do

Call band_the_matrix(n, Dtemp, Deriv)
deallocate(Dtemp)

end subroutine first

subroutine second(n,Deriv,A)
	real(dp), dimension(:,:), allocatable,intent(inout) :: Deriv
	real(dp),dimension(:),allocatable,intent(in) :: A
	Integer,intent(in) :: n
	integer :: i,jj
	real(dp),dimension(:,:),allocatable :: Dtemp

!!! this sets up the second order derivative

!!! this sets up the first order derivative

Deriv=0.d0

allocate(Dtemp(1:n,1:n)); Dtemp = 0.d0

	Dtemp(1,1:6) = D2(1,1:6)*A(1)
	Dtemp(2,1:6) = D2(2,1:6)*A(2)

Do jj = 3,n-2
Do i=-2,2
	Dtemp(jj,jj+i) = D2(3,4+i)*A(jj)
End do
END DO

!Do i = 0,-5,-1
	Dtemp(n-1,n-5:n) = D2(4,1:6)*A(n-1)
	Dtemp(n,n-5:n) = D2(5,1:6)*A(n)
!End do

Call band_the_matrix(n, Dtemp, Deriv)
deallocate(Dtemp)


end subroutine second


end module equations_builder