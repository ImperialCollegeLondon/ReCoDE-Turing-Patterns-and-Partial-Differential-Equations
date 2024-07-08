module test_demo
    use fibonacci, only: fib
    use testdrive, only: error_type, unittest_type, new_unittest, check
    implicit none
    private

    public :: collect_tests

contains

    !> Collect all exported unit tests
    subroutine collect_tests(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
        new_unittest("test_negative", test_negative),  &
        new_unittest("test_0_or_1", test_0_or_1), &
        new_unittest("input_greater_than_1", test_greater_than_1) &
        ]
    end subroutine collect_tests

    !> Check input is negative
    subroutine test_negative(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call check(error, fib(-6), -1)

    end subroutine test_negative

    !> Check input is zero or one
    subroutine test_0_or_1(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        call check(error, fib(0), 0)
        call check(error, fib(1), 1)

    end subroutine test_0_or_1

    !> Check input is grerater than one
    subroutine test_greater_than_1(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      integer :: n

      do n = 2, 10
          call check(error, fib(n), fib(n-1) + fib(n-2))
      end do

  end subroutine test_greater_than_1

end module test_demo

program tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite
  use test_demo, only : collect_tests
  implicit none
  integer :: stat

  stat = 0
  call run_testsuite(collect_tests, error_unit, stat)

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if
end program tester
