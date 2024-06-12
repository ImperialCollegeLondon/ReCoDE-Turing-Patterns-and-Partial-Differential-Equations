Module domain
use type_kinds, only : dp

!{ Builds/computes the domains of the problem - 
!        each domain has a computational domain and physical domain and corresponding metric terms
!         Computatioanl grid goes from 0 to 1
!         Physical grid goes from left to right (read in from reader)
!         
!        }
!!

!!!!! domains parameters
!!! xdomain parameters
!____________________________________________________
! nx                number of points in x domain
! xdom              physical x domain (vector)    
! xcdom             computational x domain
! xl                x left boundary           
! xr                x right boundary
! dx_vec            vector of x distances
! dxc_vec           vector of computational distance (always) contant ! might remove 
! dxc               scalar computational distance
! dxcsq             scalar computational distance square
! xmetric1          dx/dc vector
! xmetric1sq        dx/dc vector squared
! xmetric2          dx^2/dc^2 vector
! x_grid_strech_on                    turns on (TRUE) or off (FALSE) grid stretch function
! x_singular_pertubation_location     location of singular perturbation left boundary (L) or right (r)
! xhalf                               puts halve the points intbetween left or right boundary and xhalf


integer :: nx 
real(dp), dimension(:),allocatable :: xdom, xcdom 
real(dp) :: xl,xr 
real(dp), dimension(:),allocatable :: dx_vec,dxc_vec
real(dp) :: dxc,dxcsq 
real(dp), dimension(:),allocatable :: xmetric1,xmetric1sq,xmetric2 
logical :: x_grid_strech_on
character*1 :: x_singular_pertubation_location 
real(dp) :: xhalf                     
                    
contains 



Subroutine initial_domain_settings

!!
! @brief      { This builds the physical and computational domains given 
!                 above in the Module}
!!

  !!! sets up x and xc domains and metrics
  call set_up_domain(nx,xl,xr,dxc,dx_vec,dxc_vec,xdom,xcdom,x_grid_strech_on,&
          &x_singular_pertubation_location,xhalf,xmetric1,xmetric1sq,xmetric2)

  dxcsq = dxc*dxc


  !!! will add in y-domain in future
End Subroutine initial_domain_settings




Subroutine set_up_domain(n,left,right,d,dx,dc,x,c,gs_on,&
                    &pb_loc,gs,metric1,metric1sq,metric2)

!!
! @brief      { Set up domain builds a physical domain whihch is mapped to a 
! computational domain via grid streching. The computational domain will have 
! evenly space grid points }
!
! @param      n             How many points in domain
! @param      left          The left boundary 
! @param      right         The right bounday 
! @param      gs_on         Grid strech on - is the grid streching on? (logical)
! @param      pb_loc        The pertubation location (left or right)
!                            - if gs_on true then is the equation singular at left or right boundary (character)
! @param      gs            The grid strech parameter 
!                            - if gs_on true then half the grid points will be clustered between pb_loc and gs
! 
! @Return      d            computational distance (scalar + constant)
! @Return      dx           physical distance (vector)
! @Return      dc           computational distance (vector + constant)
! @Return      x            The physical domain (vector)
! @Return      c            Computational domain (vector)
! @Return      metric1      The metric dx/dc (vector)
! @Return      metric1sq    The metric dx/dc square vector
! @Return      metric2      The metric dx^2/dc^2 (vector)
!
!!

integer, intent(in)  :: n 
real(dp),intent(in)  :: left,right 
real(dp),intent(out) :: d              
real(dp),intent(out), dimension(:),allocatable :: dx,dc 
real(dp),intent(out), dimension(:),allocatable :: x,c 

!!!! grid streching settings
logical,intent(in) :: gs_on
character*1,intent(in) :: pb_loc
real(dp),intent(in) :: gs
real(dp),intent(out),dimension(:),allocatable :: metric1,metric2,metric1sq

integer :: i,k
real(dp) :: sizer 

  ! allocate the vectors
  allocate(x(1:n),c(1:n),dx(1:n),dc(1:n),metric1(1:n),metric1sq(1:n),metric2(1:n))

  ! set size of domain
  sizer = right - left

  !! set computation distance
  d = 1.d0/(n-1.d0)
  dc = d

  ! set the left boundary
  c(1) = 0.d0
  ! first distance is zero - only have one point
  dx(1) = 0.d0

  ! build the computational domain (evenly spaced by d)
  Do i = 2,n
    c(i) = (i-1)*d
  End Do

  ! build the physical domain 

  Select Case (gs_on)
  Case(.FALSE.)
    call mapping_no_strech(c,left,sizer,x)
  Case(.TRUE.)
    call mapping_strech(n,c,left,right,sizer,pb_loc,gs,x)
  End Select
   
  ! find the derivatives xdom.. i.e. the metrics
  call metrics(n,x,d,metric1,metric1sq,metric2)

Return
End Subroutine set_up_domain

Subroutine mapping_no_strech(c,left,size,mapping)
real(dp),dimension(:),allocatable,intent(in) :: c
real(dp),intent(in) :: left,size
real(dp),dimension(:),allocatable,intent(inout) ::  mapping

  !!!! simple linear domain mapping
  mapping(:) = left + size*c(:)

Return
End Subroutine

Subroutine mapping_strech(n,c,left,right,size,pb_loc,half,mapping)

!!
! @brief  {builds the mapping from the computational plane to the physical plane 
!                  dependendent on grid streching. 
!                  If grid streching is on then the function clusters 
!                  half the points between one boundary and pb_loc. 
!                  If grid streching is off then the function is a linear map
!                  
!                  The function is defined by:: 

!                  eta = a\lambda/(b-\lambda)
!                  b = 1 + a
!                  a = h/(1-2*h)
!                  
!                  where 
!                  \lambda is computational distance between 0 and 1
!                  \eta is the physical distance between 0 and 1
!                  h is the clustering point: half the grid points are placed between 0 and h
!                  (h must be between 0 and 1 and not 1/2)
!
! @param      n               Size of the domain
! @param      on_off          Turns on the grid streching if true
! @param      c               computational domain
! @param      left            The left boundary 
! @param      right           The right boundary
! @param      size            The size of the domain
! @param      pb_loc          The pertubation location (left or right)
! @param      half            Half the grid points will be clustered between a boundary and variable half
!
! @Return     physical domain coordinate corresponding to computational domain gird
!!

integer,intent(in) :: n ! size of domain
real(dp),dimension(:),allocatable,intent(in) :: c
real(dp),intent(in) :: left,right,half,size
character*1,intent(in) :: pb_loc
real(dp),dimension(:),allocatable,intent(inout)  ::  mapping
real(dp)  :: a,b,half_temporary

  ! we build the mapping for the domain size [0, 1] and then rescale back
  ! we therefore scale half to be half within [0,1]  
   
  half_temporary = (half - left)/size

  Select Case(pb_loc)
  Case('L','l')
   ! boundary layer is on the left boundary

    a = 1.d0*half_temporary/(1.d0 - 2.d0*half_temporary)
    b = 1.d0 + a/1.d0
    mapping(:) = a*c(:)/(b-c(:))

    !! scale back up

    mapping(:) = mapping(:)*size + left

  Case('R','r')
    WRITE(6,*) 'Have not implemented yet pertubation at right boundary yet'
    STOP
  End Select

Return
End Subroutine mapping_strech



Subroutine metrics(n,x,ch,metric1,metric1sq,metric2)

!!
! @brief      { Computes the metrics of the coordinate mapping}
!
! @param      n          Size of the domain
! @param      x          Physical coodinate (vector)
! @param      ch         Computational distance (scaler)
! 
! @Return      metric1      The metric dx/dc (vector)
! @Return      metric1sq    The metric dx/dc square vector
! @Return      metric2      The metric dx^2/dc^2 (vector)
!!

integer,intent(in) :: n ! size of domain
real(dp),dimension(:),allocatable,intent(in) :: x
real(dp),intent(in) :: ch ! computatioanl distance
real(dp),dimension(:),allocatable,intent(inout) :: metric1,metric2,metric1sq

  !! dx and dxx are dx/dc and d^2/dc^2 respectively
  !! calculate using finite differences Subroutine 

  !!!! find the first difference
  call first_difference(n,x,metric1)

  !!!! find the second difference
  call second_difference(n,x,metric2)

  !!! divide by computational distances to turn differences 

  metric1   = metric1/ch
  metric2   = metric2/(ch*ch)
  metric1sq = metric1*metric1

End Subroutine metrics



Subroutine first_difference(n, func, deriv)
use maths_constants, only : D1,D2
integer, intent(in) :: n ! dimension of vector
real(dp), dimension(:),allocatable,intent(in) ::    func
real(dp), dimension(:),allocatable,intent(inout) :: deriv
integer :: i

!!
! @brief      {This function computes the first difference of the function func
!               It uses the constants D1,D2 which contain the relevent finite difference
!               coefficients }
!
! @param      n      size of domain
! @param      func   function to difference 
!
! @return     deriv  the difference function 


  deriv = 0.d0 

  ! the first/second point in the difference is the 
  ! sum of the function * first/second row coefficient matrix for the first six points 
  ! note we are using sided difference here and fourth order accurate sided difference require 6 points

  deriv(1)   = Sum(D1(1,1:6)*func(1:6))
  deriv(2)   = Sum(D1(2,1:6)*func(1:6))

  ! the last/penultimate point in the difference is the 
  ! sum of the function * last/penultimate row coefficient matrix for the first six points
  deriv(n-1) = Sum(D1(4,1:6)*func(n-5:n))
  deriv(n)   = Sum(D1(5,1:6)*func(n-5:n))


  ! the middle points use central differences and only require 5 points (entires 2 to 6 in the D1(3,:) row)
  ! in the above cases we summed across 6 points and whereas now we sum across 6 points and n-4 times. 
  ! As n-4>6 we construct the vectorisation in n. This should be more efficient.

  Do i = 2,6
    deriv(3:n-2) = deriv(3:n-2) + D1(3,i)*func((3-4+i):(n-2-4+i))
  End Do 

Return
End Subroutine first_difference

Subroutine second_difference(n, func, deriv)
use maths_constants, only : D1,D2
integer, intent(in) :: n ! dimension of vector
real(dp), dimension(:),allocatable,intent(in) ::    func
real(dp), dimension(:),allocatable,intent(inout) :: deriv
integer :: i

!!
! @brief      {This function computes the first difference of the function func
!               It uses the constants D1,D2 which contain the relevent finite difference
!               coefficients }
!
! @param      n      size of domain
! @param      func   function to difference
!
! @return     deriv  the second difference function 


  deriv = 0.d0 

  ! the first/second point in the difference is the 
  ! sum of the function * first/second row coefficient matrix for the first six points 
  ! note we are using sided difference here and fourth order accurate sided difference require 6 points

  deriv(1)   = Sum(D2(1,1:6)*func(1:6))
  deriv(2)   = Sum(D2(2,1:6)*func(1:6))

  ! the last/penultimate point in the difference is the 
  ! sum of the function * last/penultimate row coefficient matrix for the first six points
  deriv(n-1) = Sum(D2(4,1:6)*func(n-5:n))
  deriv(n)   = Sum(D2(5,1:6)*func(n-5:n))


  ! the middle points use central differences and only require 5 points (entires 2 to 6 in the D1(3,:) row)
  ! in the above cases we summed across 6 points and whereas now we sum across 6 points and n-4 times. 
  ! As n-4>6 we construct the vectorisation in n. This should be more efficient.

  Do i = 2,6
    deriv(3:n-2) = deriv(3:n-2) + D2(3,i)*func((3-4+i):(n-2-4+i))
  End Do 

Return
End Subroutine second_difference

End Module domain

















