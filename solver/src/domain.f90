module domain
use type_kinds, only : dp

!!!!! domains
integer :: nx ! size of domain in x
real(dp), dimension(:),allocatable :: xdom !physical xdomain
real(dp), dimension(:),allocatable :: xcdom !x computational xdomain
real(dp) :: xl,xr ! start and end points of x domain
real(dp), dimension(:),allocatable :: dx_vec,dxc_vec ! dx/dxc distances
real(dp) :: dxc,dxcsq ! computational distance
!!!!!

contains 

subroutine initial_domain_settings
integer :: i

!!! to be call first - sets xdom array
call set_up_domain(nx,xl,xr,dxc,dx_vec,dxc_vec,xdom,xcdom)
dxcsq = dxc*dxc

end subroutine initial_domain_settings


subroutine set_up_domain(n,left,right,d,dx,dc,x,c)
integer,intent(in) :: n !size_of_domain
real(dp),intent(in) :: left,right ! left and right ends of physical domain 
real(dp),intent(out) :: d 			!! computational distance 
real(dp),intent(out), dimension(:),allocatable :: dx,dc !physical distance/computational distance (vectors)
real(dp),intent(out), dimension(:),allocatable :: x,c ! physical coordinate, computational coordiante
real(dp),parameter :: error = 10.d-8
integer :: i,k
real(dp) :: sizer

!!! allocate the vectors
allocate(x(1:n),c(1:n),dx(1:n),dc(1:n))

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
		x(i) = left + c(i)*sizer
		dx(i) = x(i) - x(k)
	end do

return
end subroutine set_up_domain

end module domain

















