Module Diff_Coeff
use type_kinds, only: dp
implicit none

integer :: DiffOrder

REAL(dp), dimension(1:5,1:6) :: D1,D2
!!!!! D1 and D2 are the coefficient finite diff matrices
!!!!! column 1 and 2 and forward, 3 is central, 4 and 5 are backward
contains
!!! note that in D2(3,4) is the central node


Subroutine First_Deriv(N, func, Deriv)
integer, intent(in) :: N ! dimension of vector
real(dp), dimension(:),allocatable,intent(in) ::    func
real(dp), dimension(:),allocatable,intent(inout) :: deriv
integer :: i,j

Deriv = 0.d0 

Do i = 1,6 
    deriv(1) = deriv(1) + D1(1,i)*func(i)
    deriv(2) = deriv(2) + D1(2,i)*func(i)


!!! because the y value is increasing dramatically towards the end
!!! use sided differences

    deriv(n-2) = deriv(n-2) + D1(4,i)*func((n-2)-5+i)
    deriv(n-1) = deriv(n-1) + D1(4,i)*func(n-6+i)
    deriv(n) = deriv(n) + D1(5,i)*func(n-6+i)
End do

Do j = 3,n-3
Do i = 2,6
    deriv(j) = deriv(j) + D1(3,i)*func(j-4+i)
END DO
END DO

!deriv(n-2) = (func(n-1)-func(n-3))/2.d0

return
End Subroutine First_Deriv

Subroutine Second_Deriv(N, func, Deriv)
integer, intent(in) :: N ! dimension of vector
real(dp), dimension(:),allocatable,intent(in) ::    func
real(dp), dimension(:),allocatable,intent(inout) :: deriv
integer :: i,j

!!! iteratively builds the derivatives

Deriv = 0.d0 

Do i = 1,6 
    deriv(1) = deriv(1) + D2(1,i)*func(i)
    deriv(2) = deriv(2) + D2(2,i)*func(i)

    deriv(n-2) = deriv(n-2) + D2(4,i)*func((n-2)-5+i)
    deriv(n-1) = deriv(n-1) + D2(4,i)*func(n-6+i)
    deriv(n) = deriv(n) + D2(5,i)*func(n-6+i)
End do

Do j = 3,n-3
Do i = 2,6
    deriv(j) = deriv(j) + D2(3,i)*func(j-4+i)
END DO
END DO

return
End Subroutine Second_Deriv


Subroutine diff_initialisation
integer :: i,j 

D1 = 0.d0; D2 = 0.d0
IF (DiffOrder==4) THEN
    !!!! First Derivative
    !! f1 
    D1(1,1) = -25.d0/12.d0; D1(1,2) = 4.d0; D1(1,3) = -3.d0; 
    D1(1,4) = 4.d0/3.d0; D1(1,5) = -1.d0/4.d0
    !! f2 
    D1(2,1) = -1.d0/5.d0; D1(2,2) = -13.d0/12.d0; D1(2,3) = 2.d0; 
    D1(2,4) = -1.d0; D1(2,5) = 1.d0/3.d0; D1(2,6) = -1.d0/20.d0
    !! C
    D1(3,2) = 1.d0/12.d0; D1(3,3) = -2.d0/3.d0; D1(3,4) = 0.d0; 
    D1(3,5) = 2.d0/3.d0; D1(3,6)=-1.d0/12.d0

!!!! Second Derivative
    !! f1 
    D2(1,1) = 15.d0/4.d0; D2(1,2) = -77.d0/6.d0; D2(1,3) = 107.d0/6.d0; 
    D2(1,4) = -13.d0; D2(1,5) = 61.d0/12.d0; D2(1,6) = -5.d0/6.d0
    !! f2 
    D2(2,1) = 5.d0/6.d0; D2(2,2) = -5.d0/4.d0; D2(2,3) = -1.d0/3.d0; 
    D2(2,4) = 7.d0/6.d0; D2(2,5) = -1.d0/2.d0; D2(2,6) = 1.d0/12.d0
    !! C
    D2(3,1) = 0.d0; D2(3,2) = -1.d0/12.d0; D2(3,3) = 4.d0/3.d0; 
    D2(3,4) = -5.d0/2.d0; D2(3,5)=4.d0/3.d0; D2(3,6) = -1.d0/12.d0

ELSE !diffOrder ==2 

!!!! First Derivative
    !! f1 
    D1(1,1) = -11.d0/6.d0; D1(1,2) = 3.d0; D1(1,3) = -3.d0/2.d0; D1(1,4) = 1.d0/3.d0
!   D1(1,1) = -3.d0/2.d0; D1(1,2) = 2.d0; D1(1,3) = -1.d0/2.d0
    !D1(1,1) = -25.d0/12.d0; D1(1,2) = 4.d0; D1(1,3) = -3.d0; D1(1,4) = 4.d0/3.d0; D1(1,5) = -1.d0/4.d0
    !! f2 
    D1(2,1) = -1.d0/3.d0; D1(2,2) = -1.d0/2.d0; D1(2,3) = 1.d0; D1(2,4) = -1.d0/6.d0
    !! C
    D1(3,1) = 0.d0; D1(3,3) = -1.d0/2.d0; D1(3,4) = 0.d0; D1(3,5) = 1.d0/2.d0

!!!! Second Derivatives
    D2(1,1) = 35.d0/12.d0; D2(1,2) = -26.d0/3.d0; 
    D2(1,3) = 19.d0/2.d0; D2(1,4) = -14.d0/3.d0; D2(1,5) = 11.d0/12.d0
        !! f2 
    D2(2,1) = 1.d0; D2(2,2) = -2.d0; D2(2,3) = 1.d0; D2(2,4) = 0.d0; D2(2,5) = 0.d0; d2(2,6) = 0.d0
        !! C
    D2(3,3) = 1.d0; D2(3,4) = -2.d0; D2(3,5) = 1.d0

END IF

Do i = 1,6
!b2
    D1(4,i) = -D1(2,6+1-i)
    D2(4,i) = D2(2,6+1-i)
!b1 
    D1(5,i) = -D1(1,6+1-i)
    D2(5,i) = D2(1,6+1-i)
End Do

!!!!! x domain
End Subroutine diff_initialisation
end Module Diff_Coeff