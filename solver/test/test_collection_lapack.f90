!!{Module test to test LAPACK and the lacpack Caller Subroutine solver_banded_double_precision}
!!
Module test_collection_lapack
   use linear_algebra, only: band_the_matrix, solver_banded_double_precision
   use testdrive, only: error_type, unittest_type, new_unittest, check
   implicit none

   public :: collect_tests_lapack

contains

   ! @brief      Collect all exported unit tests
   !
   ! @return     testsuite - testvalue
   !!
   Subroutine collect_tests_lapack(testsuite)

      !! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("test lapack banded solver single percision - one vector ", test_solver_dp_1) &
                  ]

   End Subroutine collect_tests_lapack

  !!
   ! @brief      Check lapack solver_banded_double_precision is correct
   !
   ! @return     error
  !!
   Subroutine test_solver_dp_1(error)
      use type_kinds, only: dp

      type(error_type), allocatable, intent(out) :: error
      real(dp), dimension(:, :), allocatable :: amat, amat2
      real(dp), dimension(:), allocatable :: bmat, soln
      real(dp), dimension(:), allocatable  :: output
      real(dp) :: check_value, error_value
      integer  :: i, j
      integer  :: n, lda, ldb, nrhs, sub_diag, sup_diag, nband

      ! set up the equation

      n = 2
      lda = n
      ldb = n
      nrhs = 1

      sub_diag = 1
      sup_diag = 1
      nband = sub_diag + sup_diag + 1

      allocate (amat(LDA, N), bmat(n), soln(n))

      amat = reshape([2., 3., 2., 2.], [2, 2])
      bmat(:) = [3., 6.]
      soln(:) = [3., -1.5]

      ! band the equation

      Call band_the_matrix(n, amat, sub_diag, sup_diag, nband, amat2)

      Call solver_banded_double_precision(n, nband, sub_diag, sup_diag, amat2, bmat, output)

!!! is the error within float limit
      error_value = 1./(10.)**10.

      Do i = 1, n
         check_value = abs(soln(i) - output(i))
         Call check(error, check_value < error_value)
      End Do

      deallocate (soln, bmat, amat, output, amat2)

!!  Amat:
!    2.00000000       3.00000000
!    2.00000000       2.00000000
!  bmat:    3.00000000      -1.50000012
! solution: x vect  3.000000, -1.500000

   End Subroutine test_solver_dp_1

End Module test_collection_lapack
