module equations_builder
use type_kinds, only : dp
use maths_constants, only : D1,D2

contains 

subroutine band_the_matrix(n,A,kl,ku,LDAB,AB)
!!! puts the matrix A in banded form (AB) with dimensions (LDAB,N)

	integer,intent(in) :: n ! original size of A (nxn)
	real(dp),dimension(:,:),allocatable,intent(in) :: A
	integer,intent(in) :: kl,ku ! number of sub/super diagonals
	integer,intent(out) :: ldab ! dimension of AB
	real(dp),dimension(:,:),allocatable,intent(out) :: AB ! output matrix
	integer :: i,j,k

	LDAB = KL + KU + 1
	
	allocate(AB(LDAB,N))

	do j=1,n
	do i=max(1,j-KU),min(N,j+KL)
	  AB(KU+1+i-j,j) = A(i,j)
	end do
	end do
	
	return
end subroutine band_the_matrix


subroutine scales(A,B,C,D,ai,bi,ci,di,ch,chsq,metric1,metric1sq,metric2)
real(dp),intent(in) :: ai,bi,ci,di ! entries from the equation
real(dp),intent(in) :: ch,chsq ! computational distance and computational distance squared
real(dp),intent(in) :: metric1,metric1sq,metric2 ! metric terms
real(dp), intent(out) :: A,B,C,D ! outputs

!!!! this subroutine evaluates the terms of the equations in the
!!!! computational domain and corrects them with the relevent metrics

	A=  ai
	B = (bi*metric1-ai*metric2/metric1)*ch
	C = ci*metric1sq*chsq
	D = di*metric1sq*chsq

end subroutine scales

subroutine derivative_runner(n,A,B,C,L)
integer,intent(in) :: n
	real(dp),dimension(:),allocatable,intent(in) :: A,B,C
	real(dp),dimension(:,:),allocatable,intent(inout) :: L
	real(dp),dimension(:,:),allocatable :: eD0,eD1,eD2

!!! this builds the entire derivative L

allocate(eD0(1:n,1:n),eD1(1:n,1:n),eD2(1:n,1:n))

ed0 = 0.d0;ed1=0.d0;ed2=0.d0
	
	Call zero(n,ed0,C)
		L=ed0
	Call first(n,ed1,B)
		L=L+ed1
	Call second(n,ed2,A)
		L=L+ed2

deallocate(ed0,ed1,ed2)

end subroutine derivative_runner


!above are external
!below are internal

subroutine zero(n,Deriv,C)
	real(dp), dimension(:,:), allocatable,intent(inout) :: Deriv
	real(dp),dimension(:),allocatable,intent(in) :: C
	integer,intent(in) :: n
	integer :: ii

!!! this sets up the zeroth order derivative

Deriv=0.d0

Do ii=1,n
	Deriv(ii,ii) = C(ii)
End do

end subroutine zero


subroutine first(n,Deriv,B)
	real(dp), dimension(:,:), allocatable,intent(inout) :: Deriv
	real(dp),dimension(:),allocatable,intent(in) :: B
	Integer,intent(in) :: n
	integer :: i,jj	

!!! this sets up the first order derivative

Deriv=0.d0

Do i = 1,6
	Deriv(1,i) = D1(1,i)*B(1)
	Deriv(2,i) = D1(2,i)*B(2)
End do

Do jj = 3,n-2
Do i=-2,2
	Deriv(jj,jj+i) = D1(3,4+i)*B(jj)
End do
END DO

Do i = 0,-5,-1
	Deriv(n-1,n+i) = D1(4,6+i)*B(n-1)
	Deriv(n,n+i) = D1(5,6+i)*B(n)
End do


end subroutine first

subroutine second(n,Deriv,A)
	real(dp), dimension(:,:), allocatable,intent(inout) :: Deriv
	real(dp),dimension(:),allocatable,intent(in) :: A
	Integer,intent(in) :: n
	integer :: i,jj

!!! this sets up the second order derivative

Deriv=0.d0

Do i = 1,6
	Deriv(1,i) = D2(1,i)*A(1)
	Deriv(2,i) = D2(2,i)*A(2)
End do

Do jj = 3,n-2
Do i=-2,2
	Deriv(jj,jj+i) = D2(3,4+i)*a(jj)
End do
END DO

Do i = 0,-5,-1
	Deriv(n-1,n+i) = D2(4,6+i)*A(n-1)
	Deriv(n,n+i) = D2(5,6+i)*A(n)
End do


end subroutine second


end module equations_builder