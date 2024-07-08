module fibonacci
  implicit none
  private

  public :: fib
contains

  recursive function fib(n) result(f)
    use, intrinsic :: iso_fortran_env, only : error_unit
    integer, intent(in) :: n
    integer :: f

    if (n < 0) then
      write(error_unit, "(A)") "Error: n must be non-negative"
      ! This should end the program with `error stop 1` but that's untestable
      f = -1
    else if (n <= 1) then
      f = n
    else
      f = fib(n-1) + fib(n-2)
    end if
  end function fib

end module fibonacci
