module equations_definition_test
use type_kinds, only : dp
use maths_constants,only : pi,ex
contains

subroutine equation1(x,A,B,C,D)
real(dp),intent(in) :: x ! xpoisiton in the domain
real(dp),intent(out) :: A,B,C,D

!!!! this equation is in the form
!!! A u_xx + B u_x + C u = D

A = 1.d0
B = 0.d0
C = -1.d0
D = 0.d0

end subroutine equation1

subroutine equation1_BC_Bot(x,A,B,C,D)
real(dp),intent(in) :: x ! xpoisiton in the domain
real(dp),intent(out) :: A,B,C,D

!!!! bottom/lefthand boundary

!!!! this equation is in the form
!!! A u_xx + B u_x + C u = D

	A = 0.d0
	B = 0.d0
	C = 1.d0
	D = 1.d0

end subroutine equation1_BC_Bot

subroutine equation1_BC_Top(x,A,B,C,D)
real(dp),intent(in) :: x ! xpoisiton in the domain
real(dp),intent(out) :: A,B,C,D

!!!! top/righthand boundary

!!!! this equation is in the form
!!! A u_xx + B u_x + C u = D

	A = 0.d0
	B = 0.d0
	C = 1.d0
	D = ex

end subroutine equation1_BC_Top

end module equations_definition_test


