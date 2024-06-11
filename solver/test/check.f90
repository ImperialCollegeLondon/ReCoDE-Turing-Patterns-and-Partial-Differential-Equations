program tester
  use type_kinds,only : error_unit
  use testdrive, only : run_testsuite
  use test_collection_lapack, only : collect_tests_lapack
  use test_collection1, only : collect_tests1
  implicit none
  integer :: stat

  stat = 0
  call run_testsuite(collect_tests_lapack, error_unit, stat)

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed! - LAPACK"
    error stop
  end if
  WRITE(6,*)
  
  stat = 0
  call run_testsuite(collect_tests1, error_unit, stat)
  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed! - LAPACK"
    error stop
  end if
  WRITE(6,*)
  
end program tester




