!!
!{ Builds/computes the domains of the problem -
!  each domain has a computational domain and physical domain and corresponding metric terms
!  Computatioanl grid goes from 0 to 1
!  Physical grid goes from left to right (read in from reader)
!  }
!!
!
! xdomain parameters
!____________________________________________________
! nx                number of points in x domain
! xdom              physical x domain (vector)
! xcdom             computational x domain
! xl                x left boundary
! xr                x right boundary
! dxc               scalar computational distance
! dxcsq             scalar computational distance square
! xmetric1          dx/dc vector
! xmetric1sq        dx/dc vector squared
! xmetric2          dx^2/dc^2 vector
! x_grid_strech_on  turns on (TRUE) or off (FALSE) grid stretch function
! xhalf             puts half the points intbetween left boundary and xhalf
! 
! 
! time domain parameters
!____________________________________________________
! nt                number of points in time
! tdom              time domain (vector)
! tl                initial time
! tr                final time
!!
Module domain
   use type_kinds, only: dp
   use omp_lib

!!! x domain
   integer :: nx
   real(dp), dimension(:), allocatable :: xdom, xcdom
   real(dp) :: xl, xr
   real(dp) :: dxc, dxcsq
   real(dp), dimension(:), allocatable :: xmetric1, xmetric1sq, xmetric2
   logical :: x_grid_strech_on
   real(dp) :: xhalf

!!! time domain
   integer :: nt
   real(dp), dimension(:), allocatable :: tdom
   real(dp) :: tl, tr, dt

contains

!!
! @brief      { This builds the physical and computational domains given
!               above in the module
!               }
!!
   Subroutine initial_domain_settings

      !Set up time domain
      Call time_domain

      ! sets up x and xc domains and metrics
      Call set_up_domain(nx, xl, xr, dxc, xdom, xcdom, x_grid_strech_on,&
              &xhalf, xmetric1, xmetric1sq, xmetric2)

      dxcsq = dxc*dxc

      ! will add in y-domain in future
   End Subroutine initial_domain_settings

   !!
   ! @brief      {Set up time domain}
   !!
   Subroutine time_domain
   integer :: i

      !!! set up time domain
      allocate(tdom(1:nt))
      tdom(1) = tl

      dt = (tr - tl)/(nt - 1.d0)

      !$omp Parallel Do
         Do i = 2, nt
            tdom(i) = tl + (i - 1)*dt
         End Do
      !$omp End Parallel Do

      If (tdom(nt)==tr) then
      Else
         Write (6,*) 'Error 1 in time_domain'
      End if

   End Subroutine time_domain

!!
! @brief      { Set up domain builds a physical domain whihch is mapped to a
!               computational domain via grid streching.
!               The computational domain will have evenly space grid points
!               }
!
! @param      n             How many points in domain
! @param      left          The left boundary
! @param      right         The right bounday
! @param      gs_on         Grid strech on - is the grid streching on? (logical)
! @param      gs            The grid strech parameter
!                            - If gs_on true then half the grid points will be clustered
!                            - between the left boundary and gs
!
! @Return      d            computational distance (scalar + constant)
! @Return      x            The physical domain (vector)
! @Return      c            Computational domain (vector)
! @Return      metric1      The metric dx/dc (vector)
! @Return      metric1sq    The metric dx/dc square vector
! @Return      metric2      The metric dx^2/dc^2 (vector)
!
!!
   Subroutine set_up_domain(n, left, right, d, x, c, gs_on, gs, metric1, metric1sq, metric2)

      integer, intent(in)  :: n
      real(dp), intent(in)  :: left, right
      real(dp), intent(out) :: d
      real(dp), intent(out), dimension(:), allocatable :: x, c

!!!! grid streching settings
      logical, intent(in)  :: gs_on
      real(dp), intent(in) :: gs
      real(dp), intent(out), dimension(:), allocatable :: metric1, metric2, metric1sq

      integer :: i, k
      real(dp) :: sizer

      ! allocate the vectors
      allocate (x(1:n), c(1:n), metric1(1:n), metric1sq(1:n), metric2(1:n))

      ! set size of domain
      sizer = right - left

      ! set computation distance
      d = 1.d0/(n - 1.d0)

      ! set the left boundary
      c(1) = 0

      ! build the computational domain (evenly spaced by d)
      !$omp Parallel Do
         Do i = 2, n
            c(i) = (i - 1)*d
         End Do
      !$omp End Parallel Do

      ! build the physical domain

      Select Case (gs_on)
      Case (.FALSE.)
         Call mapping_no_strech(c, left, sizer, x)
      Case (.TRUE.)
         Call mapping_strech(n, c, left, sizer, gs, x)
      End Select

      ! find the derivatives xdom.. i.e. the metrics
      Call metrics(n, x, d, metric1, metric1sq, metric2)

      !Do i = 1,n
      !  Write(6,*)i, metric1(i),metric2(i)
      !End Do

      Return
   End Subroutine set_up_domain

!!
! @brief      { builds the linear mapping function which maps the physical domain
!               (from left to right) to the computational domain [0,1] with evenly spaced
!               points
!               }
!
! @param      c        computatioanal domain
! @param      left     The left boundary
! @param      size     The size of the boundary
! @param      mapping  The mapping
!
! @Return     mapping  Physical domain coordinate corresponding to computational domain gird
!!
   Subroutine mapping_no_strech(c, left, size, mapping)
      real(dp), dimension(:), allocatable, intent(in) :: c
      real(dp), intent(in) :: left, size
      real(dp), dimension(:), allocatable, intent(inout) ::  mapping

      ! simple linear domain mapping
      mapping(:) = left + size*c(:)

      Return
   End Subroutine

!!
! @brief  {builds the mapping from the computational plane to the physical plane
!          depEndEndent on grid streching.
!          If grid streching is on then the function clusters
!          half the points between left boundary and the variable gs
!          If grid streching is off then the function is a linear map
!
!          The function is defined by::

!          eta = a\lambda/(b-\lambda)
!          b = 1 + a
!          a = h/(1-2*h)
!
!          where
!          \lambda is computational distance between 0 and 1
!          \eta is the physical distance between 0 and 1
!          h is the clustering point: half the grid points are placed between 0 and h
!          (h must be between 0 and 1 and not 1/2)
!          }
!
! @param      n          Size of the domain
! @param      on_off     Turns on the grid streching If true
! @param      c          Computational domain
! @param      left       The left boundary
! @param      size       The size of the domain
! @param      gs         Half the grid points will be clustered between a left boundary and variable gs
!
! @Return     mapping    Physical domain coordinate corresponding to computational domain gird
!!
   Subroutine mapping_strech(n, c, left, size, gs, mapping)

      integer, intent(in) :: n
      real(dp), dimension(:), allocatable, intent(in) :: c
      real(dp), intent(in) :: left, gs, size
      real(dp), dimension(:), allocatable, intent(inout)  ::  mapping
      real(dp)  :: a, b, half_temporary

      ! we build the mapping for the domain size [0, 1] and then rescale back
      ! we therefore scale gs to be half within [0,1]

      half_temporary = (gs - left)/size

      !half_temporary cannot equal 0.5 otherwise divide by 0
      If ((half_temporary - 0.5d0) == 0.d0) then
         half_temporary = 0.5000001d0
      End If

      ! build the mapping from inbetween 0 and 1

      a = 1.d0*half_temporary/(1.d0 - 2.d0*half_temporary)
      b = 1.d0 + a/1.d0
      mapping(:) = a*c(:)/(b - c(:))

  !! scale back up

      mapping(:) = mapping(:)*size + left

      Return
   End Subroutine mapping_strech

!!
! @brief      {Computes the metrics of the coordinate mapping
!              }
!
! @param      n             Size of the domain
! @param      x             Physical coodinate (vector)
! @param      ch            Computational distance (scaler)
!
! @Return      metric1      The metric dx/dc (vector)
! @Return      metric1sq    The metric dx/dc square vector
! @Return      metric2      The metric dx^2/dc^2 (vector)
!!
   Subroutine metrics(n, x, ch, metric1, metric1sq, metric2)

      integer, intent(in) :: n ! size of domain
      real(dp), dimension(:), allocatable, intent(in) :: x
      real(dp), intent(in) :: ch ! computatioanl distance
      real(dp), dimension(:), allocatable, intent(inout) :: metric1, metric2, metric1sq

      ! dx and dxx are dx/dc and d^2/dc^2 respectively
      ! calculate using finite differences Subroutine

  !!!! find the first difference
      Call first_difference(n, x, metric1)

  !!!! find the second difference
      Call second_difference(n, x, metric2)

      ! divide by computational distances to turn differences

      metric1 = metric1/ch
      metric2 = metric2/(ch*ch)
      metric1sq = metric1*metric1

   End Subroutine metrics

!!
! @brief      {This function computes the first difference of the function func
!              It uses the constants D1,D2 which contain the relevent finite difference coefficients
!              }
!
! @param      n      size of domain
! @param      func   function to difference
!
! @return     deriv  the difference function
   Subroutine first_difference(n, func, deriv)
      use maths_constants, only: D1, D2
      integer, intent(in) :: n ! dimension of vector
      real(dp), dimension(:), allocatable, intent(in) ::    func
      real(dp), dimension(:), allocatable, intent(inout) :: deriv
      integer :: i

      deriv = 0.d0

      ! the first/second point in the difference is the
      ! sum of the function * first/second row coefficient matrix for the first six points
      ! note we are using sided difference here and fourth order accurate sided difference require 6 points

      deriv(1) = Sum(D1(1, 1:6)*func(1:6))
      deriv(2) = Sum(D1(2, 1:6)*func(1:6))

      ! the last/penultimate point in the difference is the
      ! sum of the function * last/penultimate row coefficient matrix for the first six points
      deriv(n - 1) = Sum(D1(4, 1:6)*func(n - 5:n))
      deriv(n) = Sum(D1(5, 1:6)*func(n - 5:n))

      ! the middle points use central differences and only require 5 points (entires 2 to 6 in the D1(3,:) row)
      ! in the above cases we summed across 6 points and whereas now we sum across 6 points and n-4 times.
      ! As n-4>6 we construct the vectorisation in n. This should be more efficient.
      ! Constructing the vectorisation in n is more efficient

      Do i = 2, 6
         deriv(3:n - 2) = deriv(3:n - 2) + D1(3, i)*func((3 - 4 + i):(n - 2 - 4 + i))
      End Do

      Return
   End Subroutine first_difference

!!
! @brief      {This function computes the first difference of the function func
!              It uses the constants D1,D2 which contain the relevent finite difference coefficients
!              }
!
! @param      n      size of domain
! @param      func   function to difference
!
! @return     deriv  the second difference function
   Subroutine second_difference(n, func, deriv)
      use maths_constants, only: D1, D2
      integer, intent(in) :: n ! dimension of vector
      real(dp), dimension(:), allocatable, intent(in) ::    func
      real(dp), dimension(:), allocatable, intent(inout) :: deriv
      integer :: i

      deriv = 0.d0

      ! the first/second point in the difference is the
      ! sum of the function * first/second row coefficient matrix for the first six points
      ! note we are using sided difference here and fourth order accurate sided difference require 6 points

      deriv(1) = Sum(D2(1, 1:6)*func(1:6))
      deriv(2) = Sum(D2(2, 1:6)*func(1:6))

      ! the last/penultimate point in the difference is the
      ! sum of the function * last/penultimate row coefficient matrix for the first six points
      deriv(n - 1) = Sum(D2(4, 1:6)*func(n - 5:n))
      deriv(n) = Sum(D2(5, 1:6)*func(n - 5:n))

      ! the middle points use central differences and only require 5 points (entires 2 to 6 in the D1(3,:) row)
      ! in the above cases we summed across 6 points and whereas now we sum across 6 points and n-4 times.
      ! As n-4>6 we construct the vectorisation in n. This should be more efficient.
      ! Constructing the vectorisation in n is more efficient

      Do i = 2, 6
         deriv(3:n - 2) = deriv(3:n - 2) + D2(3, i)*func((3 - 4 + i):(n - 2 - 4 + i))
      End Do

      Return
   End Subroutine second_difference

End Module domain

