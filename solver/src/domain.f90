module domain
use type_kinds, only : dp

!!!!! domains
!!! xdomain 

integer :: nx ! size of domain in x
real(dp), dimension(:),allocatable :: xdom !physical xdomain
real(dp), dimension(:),allocatable :: xcdom !x computational xdomain
real(dp) :: xl,xr ! start and end points of x domain
real(dp), dimension(:),allocatable :: dx_vec,dxc_vec ! dx/dxc distances
real(dp) :: dxc,dxcsq ! computational distance
real(dp), dimension(:),allocatable :: xmetric1,xmetric1sq,xmetric2 ! xc derivative and xcc derivative
!!!!!

!!!!!! grid streching 
!!!!!! x grid streching
logical :: x_grid_strech_on
Character*1 :: x_singular_pertubation_location ! location of singular perturbation l or r
real(dp) :: xhalf ! puts halve the points intbetween left or right boundary and xhalf


contains 

subroutine initial_domain_settings
integer :: i

	!!! to be call first - sets xdom array
	call set_up_domain(nx,xl,xr,dxc,dx_vec,dxc_vec,xdom,xcdom,x_grid_strech_on,&
			&x_singular_pertubation_location,xhalf,xmetric1,xmetric1sq,xmetric2)
	dxcsq = dxc*dxc

end subroutine initial_domain_settings


subroutine set_up_domain(n,left,right,d,dx,dc,x,c,grid_strech_on,&
						&pertubation_location,grid_strech_parameter,metric1,metric1sq,metric2)

	integer,intent(in) :: n !size_of_domain
	real(dp),intent(in) :: left,right ! left and right ends of physical domain 
	real(dp),intent(out) :: d 			!! computational distance 
	real(dp),intent(out), dimension(:),allocatable :: dx,dc !physical distance/computational distance (vectors)
	real(dp),intent(out), dimension(:),allocatable :: x,c ! physical coordinate, computational coordiante

!!!! grid streching settings
	logical,intent(in) :: grid_strech_on
	character*1,intent(in) :: pertubation_location
	real(dp),intent(in) :: grid_strech_parameter
	real(dp),intent(out),dimension(:),allocatable :: metric1,metric2,metric1sq


real(dp),parameter :: error = 10.d-8
integer :: i,k
real(dp) :: sizer

!!! allocate the vectors
allocate(x(1:n),c(1:n),dx(1:n),dc(1:n),metric1(1:n),metric1sq(1:n),metric2(1:n))

!! set size of domain
	sizer = right - left

!! set computation distance
	d = 1.d0/(n-1.d0)
	dc = d

!!! set the left boundary
	x(1) = left
	c(1) = 0.d0
	dx(1) = 0.d0

! build other terms
	do i = 2,n
		k =  i - 1
		c(i) = (i-1)*d
		!x(i) = left + c(i)*sizer
		x(i) = mapping(grid_strech_on,c(i),left,right,sizer,pertubation_location,grid_strech_parameter)
		dx(i) = x(i) - x(k)
	end do
 

 !!!! find the derivatives of metrics
	call metrics(n,x,d,metric1,metric1sq,metric2)

	!do i = 1,n
	!	WRITE(6,*) i,c(i),x(i),metric1(i),metric1sq(i),metric2(i)
	!end do

return
end subroutine set_up_domain

subroutine metrics(n,x,ch,metric1,metric1sq,metric2)
integer,intent(in) :: n ! size of domain
real(dp),dimension(:),allocatable,intent(in) :: x
real(dp),intent(in) :: ch ! computatioanl distance
real(dp),dimension(:),allocatable,intent(inout) :: metric1,metric2,metric1sq

!! dx and dxx are dx/dc and d^2/dc^2 respectively
!! calculate using finite differences subroutine 

call first_deriv(n,x,metric1)
call second_deriv(n,x,metric2)

metric1 = metric1 /ch
metric2 = metric2 /(ch*ch)
metric1sq = metric1*metric1

end subroutine metrics

function mapping(on_off,c,left,right,size_of_domain,pertubation_location,half)
	logical,intent(in)    :: on_off
	real(dp),intent(in) 	:: c ! computational domain
	real(dp),intent(in)   :: left,right,half,size_of_domain
	character*1,intent(in) :: pertubation_location
	!! left boundary, right boundary, is the BL in the left or right boundary,half the points are below
	real(dp) ::  mapping
	!!! builds the mapping from the computational plane to the physical plane
	!!! dependendent on grid streching 

	real(dp) :: a,b

	select case(on_off)
	case(.TRUE.)

		select case(pertubation_location)
		case('L','l')
			!!!! boundary layer is on the left boundary 

			a = right*half/(right - 2.d0*half)
			b = 1.d0 + a/right
			mapping = a*c/(b-c)

		case('R','r')
			WRITE(6,*) 'Have not implemented yet pertubation at right boundary yet'
			stop
		end select

	case(.FALSE.)
	!!!! simple domain mapping

		mapping = left + size_of_domain*c
	end select

return
end function mapping



Subroutine First_Deriv(N, func, Deriv)
use maths_constants, only : D1,D2
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
use maths_constants, only : D1,D2
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



end module domain

















