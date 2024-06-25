Module maths_constants_Diff_Coeff
   use type_kinds, only: dp
   implicit none

   integer :: DiffOrder
   integer :: sub_diag, sup_diag, nband ! number of sub/super_diagonals of the given matrices

   REAL(dp), dimension(1:5, 1:6) :: D1, D2
!!!!! D1 and D2 are the coefficient finite diff matrices
!!!!! column 1 and 2 and forward, 3 is central, 4 and 5 are backward
contains
!!! note that in D2(3,4) is the central node

   Subroutine diff_initialisation
      integer :: i, j

      D1 = 0.d0; D2 = 0.d0
      select case (difforder)
      case (4)

         sub_diag = 5; sup_diag = 5; nband = sup_diag + sub_diag + 1

!!!! fourth order accurate

    !!!! First Derivative
    !! f1
         D1(1, 1) = -25.d0/12.d0; D1(1, 2) = 4.d0; D1(1, 3) = -3.d0; D1(1, 4) = 4.d0/3.d0; D1(1, 5) = -1.d0/4.d0
    !! f2
         D1(2, 1) = -1.d0/5.d0; D1(2, 2) = -13.d0/12.d0; D1(2, 3) = 2.d0; 
         D1(2, 4) = -1.d0; D1(2, 5) = 1.d0/3.d0; D1(2, 6) = -1.d0/20.d0
    !! C
         D1(3, 2) = 1.d0/12.d0; D1(3, 3) = -2.d0/3.d0; D1(3, 4) = 0.d0; D1(3, 5) = 2.d0/3.d0; D1(3, 6) = -1.d0/12.d0

!!!! Second Derivative
    !! f1
         D2(1, 1) = 15.d0/4.d0; D2(1, 2) = -77.d0/6.d0; D2(1, 3) = 107.d0/6.d0; 
         D2(1, 4) = -13.d0; D2(1, 5) = 61.d0/12.d0; D2(1, 6) = -5.d0/6.d0
    !! f2
         D2(2, 1) = 5.d0/6.d0; D2(2, 2) = -5.d0/4.d0; D2(2, 3) = -1.d0/3.d0; 
         D2(2, 4) = 7.d0/6.d0; D2(2, 5) = -1.d0/2.d0; D2(2, 6) = 1.d0/12.d0
    !! C
         D2(3, 1) = 0.d0; D2(3, 2) = -1.d0/12.d0; D2(3, 3) = 4.d0/3.d0; 
         D2(3, 4) = -5.d0/2.d0; D2(3, 5) = 4.d0/3.d0; D2(3, 6) = -1.d0/12.d0

      case (3)

         sub_diag = 4; sup_diag = 4; nband = sup_diag + sub_diag + 1

!!!! fourth order accurate in the entirior
!!!! second order accurate in the boundaries

!!!! First Derivative
    !! f1
         D1(1, 1) = -11.d0/6.d0; D1(1, 2) = 3.d0; D1(1, 3) = -3.d0/2.d0; D1(1, 4) = 1.d0/3.d0
    !! f2
         D1(2, 1) = -1.d0/5.d0; D1(2, 2) = -13.d0/12.d0; D1(2, 3) = 2.d0; 
         D1(2, 4) = -1.d0; D1(2, 5) = 1.d0/3.d0; D1(2, 6) = -1.d0/20.d0
    !! C
         D1(3, 2) = 1.d0/12.d0; D1(3, 3) = -2.d0/3.d0; D1(3, 4) = 0.d0; D1(3, 5) = 2.d0/3.d0; D1(3, 6) = -1.d0/12.d0

!!!! Second Derivatives
         D2(1, 1) = 35.d0/12.d0; D2(1, 2) = -26.d0/3.d0; D2(1, 3) = 19.d0/2.d0; 
         D2(1, 4) = -14.d0/3.d0; D2(1, 5) = 11.d0/12.d0
        !! f2
         D2(2, 1) = 5.d0/6.d0; D2(2, 2) = -5.d0/4.d0; D2(2, 3) = -1.d0/3.d0; 
         D2(2, 4) = 7.d0/6.d0; D2(2, 5) = -1.d0/2.d0; D2(2, 6) = 1.d0/12.d0
        !! C
         D2(3, 1) = 0.d0; D2(3, 2) = -1.d0/12.d0; D2(3, 3) = 4.d0/3.d0; 
         D2(3, 4) = -5.d0/2.d0; D2(3, 5) = 4.d0/3.d0; D2(3, 6) = -1.d0/12.d0

      case (2)

         sub_diag = 4; sup_diag = 4; nband = sup_diag + sub_diag + 1

!!!! second order accurate

!!!! First Derivative
    !! f1
         D1(1, 1) = -11.d0/6.d0; D1(1, 2) = 3.d0; D1(1, 3) = -3.d0/2.d0; D1(1, 4) = 1.d0/3.d0
    !! f2
         D1(2, 1) = -1.d0/3.d0; D1(2, 2) = -1.d0/2.d0; D1(2, 3) = 1.d0; D1(2, 4) = -1.d0/6.d0
    !! C
         D1(3, 1) = 0.d0; D1(3, 3) = -1.d0/2.d0; D1(3, 4) = 0.d0; D1(3, 5) = 1.d0/2.d0

!!!! Second Derivatives
         D2(1, 1) = 35.d0/12.d0; D2(1, 2) = -26.d0/3.d0; D2(1, 3) = 19.d0/2.d0; 
         D2(1, 4) = -14.d0/3.d0; D2(1, 5) = 11.d0/12.d0
        !! f2
         D2(2, 1) = 1.d0; D2(2, 2) = -2.d0; D2(2, 3) = 1.d0; D2(2, 4) = 0.d0; D2(2, 5) = 0.d0; d2(2, 6) = 0.d0
        !! C
         D2(3, 3) = 1.d0; D2(3, 4) = -2.d0; D2(3, 5) = 1.d0

      end select

!!!! terms are symmetric

      Do i = 1, 6
!b2
         D1(4, i) = -D1(2, 6 + 1 - i)
         D2(4, i) = D2(2, 6 + 1 - i)
!b1
         D1(5, i) = -D1(1, 6 + 1 - i)
         D2(5, i) = D2(1, 6 + 1 - i)
      End Do

!!!!! x domain
   End Subroutine diff_initialisation

end Module maths_constants_Diff_Coeff
