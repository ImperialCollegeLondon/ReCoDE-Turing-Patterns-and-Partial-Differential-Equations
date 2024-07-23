!!{Program that organises and runs the tests}
!!
Program tester
   use type_kinds, only: error_unit
   use testdrive, only: run_testsuite
   use test_collection_lapack, only: collect_tests_lapack
   use test_collection_BVP, only: collect_tests1
   implicit none
   integer :: stat

   stat = 0
   Call run_testsuite(collect_tests_lapack, error_unit, stat)

   If (stat > 0) then
      Write (error_unit, '(i0, 1x, a)') stat, "test(s) failed! - LAPACK"
      Error Stop
   End If
   Write (6, *)

   stat = 0
   Call run_testsuite(collect_tests1, error_unit, stat)
   If (stat > 0) then
      Write (error_unit, '(i0, 1x, a)') stat, "test(s) failed! - BVP"
      Error Stop
   End if
   Write (6, *)

End Program tester

