!!Module that sets up the relevent maths constants for the problem
!!
Module maths_constants
   use type_kinds, only: dp
   use maths_constants_fundamental, only: pi, ex
   use maths_constants_diff_coeff, only: DiffOrder, D1, D2, diff_initialisation,&
                   &sub_diag, sup_diag, nband

!! note that::: 
   ! DiffOrder   order of the equation, 2, 3 or 4
   ! sub_diag    number of sub_diagonals of the given matrices
   ! sup_diag    number of super_diagonals of the given matrices
   ! nband       resulting bandwith of the banded matrix
   ! D1 and D2 are the coefficient finite diff matrices for first and second derivatives respecitely

End Module maths_constants
