program main
  use fibonacci, only: fib
  implicit none

  character(len=12) :: arg
  integer :: num

  !stop

  call get_command_argument(1, arg)
  read(arg, *) num
  print *, fib(num)
end program main
