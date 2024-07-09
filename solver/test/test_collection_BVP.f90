!!{Collection of tests - banded matrix, domains and BVP solver}
!!
Module test_collection_BVP
   use testdrive, only: error_type, unittest_type, new_unittest, check
   implicit none

   public :: collect_tests1

contains

   !!
   ! @brief      { Collect all exported unit tests }
   !
   ! @return      testsuite  The testsuite value
   !!
   Subroutine collect_tests1(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
        new_unittest("test Domain builder (set_up_Domain) ", check_domain_builder), &
        new_unittest("test banded matrix calculations  ", check_banded_matrix), &
        new_unittest("test banded matrix multiplication ", check_banded_matrix_multiplication), &
        new_unittest("test bvp solver with no streching ", check_bvp_solver_4th_order_no_streching), &
        new_unittest("test bvp solver with streching ", check_bvp_solver_4th_order_streching) &
        &]

   End Subroutine collect_tests1

   !!
   ! @brief      {Tests the domain builder subroutine}
   !
   ! @return      error
   !!
   Subroutine check_domain_builder(error)
      use type_kinds, only: dp
      use Domain, only: set_up_Domain
      use maths_constants, only: pi

      type(error_type), allocatable, intent(out) :: error
      real(dp), dimension(:), allocatable :: x, c, metric1, metric1sq, metric2
      real(dp) :: check_value, error_value
      real(dp)  :: l, r, dc, check1, check2, check3, check4, check5
      integer  :: n, i

      error_value = 1.d-8

      n = 57
      l = 34.7d0
      r = pi/4.d0

!!! with this function the End points should match
!!! cDom should be 1 and zero
!!! dc should = dcv
      Call set_up_Domain(n, l, r, dc, x, c, .FALSE., 0.d0, metric1, metric1sq, metric2)

      check1 = abs(l - x(1))
      Call check(error, check1 < error_value)

      check2 = abs(r - x(n))
      Call check(error, check2 < error_value)

      check3 = abs(0.d0 - c(1))
      Call check(error, check3 < error_value)

      check4 = abs(1.d0 - c(n))
      Call check(error, check4 < error_value)

      Do i = 2, n
         check5 = abs(dc - (c(i) - c(i - 1)))
         Call check(error, check5 < error_value)
      End Do

      deallocate (x, c, metric1, metric1sq, metric2)
   End Subroutine check_domain_builder

   !!
   ! @brief      checks the banded matrix module
   !
   ! @return      error
   !!
   Subroutine check_banded_matrix(error)
      use type_kinds, only: dp
      use matrix_control, only: band_the_matrix
      use maths_constants, only: pi

      type(error_type), allocatable, intent(out) :: error
      real(dp), dimension(:, :), allocatable :: A, AB, AB_test
      real(dp) :: error_value
      integer  :: n, KL, KU, LDAB
      integer  :: i, j

      error_value = 1.d-8

      n = 5
      KL = 2
      KU = 1
      LDAB = KL + KU + 1
      allocate (A(n, n), AB(ldab, n))

      !! set up unbanded matrix

      A = reshape([2.d0, 3.d0, 2.d0, 0.d0, 0.d0, &
                  & 1.d0, 3.d0, 2.d0, 4.d0, 0.d0, &
                  & 0.d0, 2.d0, 2.d0, 12.d0, 8.d0, &
                  & 0.d0, 0.d0, 4.d0, 13.d0, 17.d0, &
                  & 0.d0, 0.d0, 0.d0, 7.d0, 1.d0], &
                  &[5, 5])

      Call band_the_matrix(n, A, KL, KU, LDAB, AB_test)

      AB = transpose(reshape([0.d0, 1.d0, 2.d0, 4.d0, 7.d0,&
                  &2.d0, 3.d0, 2.d0, 13.d0, 1.d0, &
                  &3.d0, 2.d0, 12.d0, 17.d0, 0.d0,&
                  &2.d0, 4.d0, 8.d0, 0.d0, 0.d0], &
                  &[5, 4]))

    !  Do i = 1, n
    !     Write(6, *) (A(i,j),j=1,n)
    !  End Do
    !  Write(6,*)
    !  Do i = 1, ldab
    !     Write(6, *) (AB(i,j),j=1,n)
    !  End Do


      Do i = 1, ldab
      Do j = 1, n
         Call check(error, ab(i, j), ab_test(i, j))
      End Do
      End Do

      deallocate (A, AB, AB_test)

      Return
   End Subroutine check_banded_matrix

      !!
   ! @brief      checks the banded matrix multiplication subroutine
   !
   ! @return      error
   !!
   Subroutine check_banded_matrix_multiplication(error)
      use type_kinds, only: dp
      use linear_algerba, only: band_the_matrix
      use maths_constants, only: pi

      external :: DGBMV
      type(error_type), allocatable, intent(out) :: error
      real(dp), dimension(:, :), allocatable :: A, AB
      real(dp), dimension(:), allocatable :: Vec,Sol1,Sol2
      real(dp) :: error_value
      integer  :: n, KL, KU, LDAB
      integer  :: i, j

      error_value = 1.d-8

      n = 5
      KL = 2
      KU = 1
      LDAB = KL + KU + 1
      allocate (A(n, n), Vec(n), Sol1(n),Sol2(n))

      !! set up unbanded matrix

      A = reshape([2.d0, 3.d0, 2.d0, 0.d0, 0.d0, &
                  & 1.d0, 3.d0, 2.d0, 4.d0, 0.d0, &
                  & 0.d0, 2.d0, 2.d0, 12.d0, 8.d0, &
                  & 0.d0, 0.d0, 4.d0, 13.d0, 17.d0, &
                  & 0.d0, 0.d0, 0.d0, 7.d0, 1.d0], &
                  &[5, 5])

      Vec = ([1.d0, 2.d0, 0.d0, 3.d0, -4.d0])

      Call band_the_matrix(n, A, KL, KU, LDAB, AB)

      Sol1 = matmul(A,vec)
      
      Call DGBMV('N', n, n, KL, KU, 1.d0, AB, LDAB, Vec, 1, 0.d0, Sol2, 1)

      Do j = 1, n
         Call check(error, Sol1(j), Sol2(j))
      End Do
      deallocate (A, AB, Vec, Sol1, Sol2)

      Return
   End Subroutine check_banded_matrix_multiplication

   !!
   ! @brief      {Checks the BVP solver is correct when there is no grid streching}
   !
   ! @return     error
   !!
   Subroutine check_bvp_solver_4th_order_no_streching(error)
      use type_kinds, only: dp
      use reader, only: read_me
      use domain
      use maths_constants
      use equations
      use linear_algebra, only: solver_banded_double_precision

      type(error_type), allocatable, intent(out) :: error

      real(dp) :: error_value, test_value
      real(dp), dimension(:), allocatable :: X
      integer :: j
      real(dp), dimension(:, :), allocatable :: L !! banded form matrix
      real(dp), dimension(:), allocatable :: RHS !! right hand side of equation

      !call read_me  ! opens settings.input and reads

      nx = 21; nt = 20
      xl = 0.d0; xr = 1.d0
      x_grid_strech_on = .FALSE.
      DiffOrder = 4

      call diff_initialisation  ! sets the finite difference coefficients
      call initial_domain_settings !builds the domain and computational domains
      call build_the_matrix(nx, dxc, dxcsq, xcdom, xmetric1, xmetric1sq, xmetric2, 1, L, RHS)
      call solver_banded_double_precision(nx, nband, sub_diag, sup_diag, L, RHS, X)
      deallocate (L, RHS)

      error_value = 1.d-7

      Do j = 1, nx
         test_value = abs(x(j) - ex**xdom(j))
      End Do

      deallocate (xdom, xcdom, X, tdom)
      Write (6, '(a,1x,e12.4)') '     Error Value::', error_value

      Return
   End Subroutine check_bvp_solver_4th_order_no_streching

   !!
   ! @brief      {Checks the BVP solver is correct when there is grid streching}
   !
   ! @return     error
   !!
   Subroutine check_bvp_solver_4th_order_streching(error)
      use type_kinds, only: dp
      use reader, only: read_me
      use domain
      use maths_constants
      use equations
      use linear_algebra, only: solver_banded_double_precision

      type(error_type), allocatable, intent(out) :: error

      real(dp) :: error_value, test_value
      real(dp), dimension(:), allocatable :: X
      integer :: j
      real(dp), dimension(:, :), allocatable :: L !! banded form matrix
      real(dp), dimension(:), allocatable :: RHS !! right hand side of equation

      nx = 21; nt = 20
      xl = 0.d0; xr = 1.d0
      x_grid_strech_on = .TRUE.
      xhalf = 0.3d0

      call diff_initialisation  ! sets the finite difference coefficients
      call initial_domain_settings !builds the domain and computational domains
      call build_the_matrix(nx, dxc, dxcsq, xcdom, xmetric1, xmetric1sq, xmetric2, 1, L, RHS)
      call solver_banded_double_precision(nx, nband, sub_diag, sup_diag, L, RHS, X)
      deallocate (L, RHS)

      error_value = 1.d-7

      Do j = 1, nx
         test_value = abs(x(j) - ex**xdom(j))
      End Do

      deallocate (xdom, xcdom, X, tdom)

      Write (6, '(a,1x,e12.4)') '     Error Value::', error_value

      Return
   End Subroutine check_bvp_solver_4th_order_streching

End Module test_collection_BVP
