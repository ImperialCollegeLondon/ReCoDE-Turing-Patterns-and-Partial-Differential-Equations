module equations
	use type_kinds, only : dp
	use equations_definition, only : equation1,equation1_BC_Top,equation1_BC_Bot
	use equations_builder, only : scales,derivative_runner,band_the_matrix
	use maths_constants, only : sub_diag,sup_diag

	real(dp),dimension(:,:),allocatable :: L !! banded form matrix
	integer :: nband 												 !! new dimension of L (nband x n)
	real(dp),dimension(:),allocatable :: RHS !! right hand side of equation

real(dp),dimension(:,:),allocatable :: Ltemp
contains

subroutine build_the_matrix(n,ch,chsq,cdom,metric1,metric1sq,metric2)

!!! sets up a matrix equation in the form L * u = RHS
!!! D is the RHS, Lx is the matrix

	integer,intent(in) :: n !! size of domain
	real(dp),intent(in) :: ch,chsq !computational distance/squared
	real(dp),dimension(:),allocatable,intent(in) :: cdom ! computational domain

!! grid streching metric terms
	real(dp),dimension(:),allocatable,intent(in) :: metric1,metric1sq,metric2


!!!! relevent vectors/infomation
	real(dp),dimension(:),allocatable :: A,B,C,D
	real(dp) :: Ain,Bin,Cin,Din
	integer :: ii,i,j


!!! allocate 
	allocate(A(n),B(n),C(n),D(n))
	allocate(Ltemp(n,n))

!!! Boundaries
!!! first call the equation
	call equation1_BC_Bot(cdom(1),ain,bin,cin,din)
!!! then apply the correct scalings - ch or chsq etc
	call scales(A(1),B(1),C(1),D(1),ain,bin,cin,din,ch,chsq,metric1(1),metric1sq(1),metric2(1))

	call equation1_BC_Top(cdom(n),ain,bin,cin,din)
	call scales(A(n),B(n),C(n),D(n),ain,bin,cin,din,ch,chsq,metric1(n),metric1sq(n),metric2(n))


!!! Interior
	do ii = 2,n-1
		call equation1(cdom(ii),ain,bin,cin,din)
		call scales(A(ii),B(ii),C(ii),D(ii),ain,bin,cin,din,ch,chsq,metric1(ii),metric1sq(ii),metric2(ii))
	end do

	call derivative_runner(n,A,B,C,Ltemp)
	RHS = D

	deallocate(c,b,d,a)

!!! put the matrix in banded form
	call band_the_matrix(n,Ltemp,sub_diag,sup_diag,nband,L)
	

	return
end subroutine build_the_matrix

end module equations