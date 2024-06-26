!!{ Subroutines for manipluating matrices. Currently:: band_the_matrix which puts a matrix in banded form}
!!
Module matrix_control
  use type_kinds, only: dp

contains

!!
! @brief      { Takes a square matrix A (n * n) and returns a banded matrix AB (nband * n)
!               This saves memory and computation time when a calculation is conducted with A}
!
!                We use the method from lapack:
!                On entry, the matrix A in band storage, in rows 1 to sub_diag+sup_diag+1.
!                The j-th column of A is stored in the j-th column of the array AB as follows:
!                AB(sup_diag+1+i-j,j) = A(i,j) for max(1,j-sup_diag)<=i<=min(n,j+sub_diag)
!
!
! @param      n      dimension of A
! @param      A      square matrix to band (n * n)
! @param      sub_diag     number of subdiagonals in A
! @param      sup_diag     number of superdiagonals in A
! @param      nband  the new leading dimension of AB
! 
! @return     AB     the banded matrix, dimension (nband * n)
!!
   Subroutine band_the_matrix(n, A, sub_diag, sup_diag, nband, AB)

      integer, intent(in) :: n
      real(dp), dimension(:, :), allocatable, intent(in) :: A
      integer, intent(in) :: sub_diag, sup_diag,nband
      real(dp), dimension(:, :), allocatable, intent(out) :: AB
      integer :: i, j, btest

      btest = sub_diag + sup_diag + 1

      If (nband==btest) then
      Else
        Write (6,*) 'Error 1 in band_the_matrix...'
      End If

      allocate (AB(nband, n))

      Do j = 1, n
      Do i = max(1, j - sup_diag), min(n, j + sub_diag)
         AB(sup_diag + 1 + i - j, j) = A(i, j)
      End Do
      End Do

      Return
   End Subroutine band_the_matrix


End Module matrix_control