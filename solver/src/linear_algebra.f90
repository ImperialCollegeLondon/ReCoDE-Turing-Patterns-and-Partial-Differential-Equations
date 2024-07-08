!!
!{ Sets up the linear algebra solving Subroutines with Lapack integration
! }
!!
Module linear_algebra
   use type_kinds, only: dp
   implicit none
   external :: dgbsvx ! lapack double precision linear solve

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
      integer, intent(in) :: sub_diag, sup_diag, nband
      real(dp), dimension(:, :), allocatable, intent(out) :: AB
      integer :: i, j, btest

      btest = sub_diag + sup_diag + 1

      If (nband == btest) then
      Else
         Write (6, *) 'Error 1 in band_the_matrix...'
      End If

      allocate (AB(nband, n))

      Do j = 1, n
      Do i = max(1, j - sup_diag), min(n, j + sub_diag)
         AB(sup_diag + 1 + i - j, j) = A(i, j)
      End Do
      End Do

      Return
   End Subroutine band_the_matrix

!!
! @brief      { This subroutine sets up and solves linear systems with double precision
!               Solves the system l * soln = rhs
!               soln and rhs are vectors size nband
!               l is a banded or non-banded matrix (dimension nband * n_input)
!
!               All lapack subroutine terms have been taken from the website and copied directly
!               netlib.org/lapack/explore-html/d1/da6/group__gbsvx_ga38273d98ae4d598529fc9647ca847ce2.html#ga38273d98ae4d598529fc9647ca847ce2
!               }
!
! @param      n_input   Leading dimension of the system
! @param      nband     Second dimension of matrix l
! @param      sub_diag  Number of subdiagonal of l
! @param      sup_diag  Number of superdiagonals of l
! @param      l         Matrix dimension nband * n (matrix)
! @param      rhs       The known rhs (vector)
!
! @return     soln      The solution of the system
!!
   Subroutine solver_banded_double_precision(n_input, nband, sub_diag, sup_diag, l, rhs, soln)
      integer, intent(in) :: n_input, nband, sub_diag, sup_diag
      real(dp), dimension(:, :), allocatable, intent(in)  :: l
      real(dp), dimension(:), allocatable, intent(in) :: rhs
      real(dp), dimension(:), allocatable, intent(out)   :: soln
      integer :: i, j

    !! Lapack specific terms - lapack gives the required dimension of each term. Copy from website.
    !! X is soln
    !! B is the rhs
    !! AB is the matrix
    !! Everything else sets the correct settings
    !! Note that in dgbsvx B and X can be a matrices. We however Do not need this

      character*1 ::   FACT
      character*1 ::  TRANS
      integer   ::  N
      integer   ::  KL
      integer   ::  KU
      integer   ::  NRHS
      real(dp), dimension(:, :), allocatable :: AB
      integer   ::  LDAB
      real(dp), dimension(:, :), allocatable  ::   AFB
      integer   ::  LDAFB
      integer, dimension(:), allocatable  ::   IPIV
      character*1 ::   EQUED
      real(dp), dimension(:), allocatable  ::  R
      real(dp), dimension(:), allocatable  ::  C
      real(dp), dimension(:, :), allocatable :: B
      integer   ::  LDB
      real(dp), dimension(:, :), allocatable :: X
      integer   ::  LDX
      real(dp)  ::  RCOND
      real(dp), dimension(:), allocatable  :: FERR
      real(dp), dimension(:), allocatable  ::  BERR
      real(dp), dimension(:), allocatable  ::  WORK
      integer, dimension(:), allocatable  ::  IWORK
      integer   ::  INFO

      !! Lapack specific settings
      !! LU decomposition used to solve the system
      FACT = 'N'
      !! We don't want to use the transpose of l
      TRANS = 'N'

      !! Setting correct diagonals and dimensions

      KL = sub_diag
      KU = sup_diag
      N = n_input
      NRHS = 1

      LDAB = nband
      allocate (AB(LDAB, n_input))
      AB = l

      !! More lapack specific settings (copied from website)

      LDAFB = 2*KL + KU + 1
      allocate (AFB(LDAFB, N))
      allocate (IPIV(1:N))

      EQUED = 'N'

      allocate (R(1:N), C(1:N))

      LDB = N ! Dimension of B
      allocate (B(LDB, NRHS))
      B(:, 1) = rhs

      LDX = N
      allocate (X(LDX, NRHS))
      allocate (BERR(NRHS), FERR(NRHS))
      allocate (WORK(max(1, 3*N)), IWORK(N))

      !! X is the output
      Call dgbsvx(FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV,&
      &EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO)

      ! If INFO is not equal to zero, then dgbsvx has failed
      ! Writing messages allows us to determine what has gone wrong
      ! These conditions are from the lapack website

      If (INFO == 0) THEN
      !!! dgbsvx has succeded
      Else If (INFO .LT. 0) THEN
         Write (6, *) 'INFO < 0, dgbsvx, has an illegal input value'
         Write (6, *) 'INFO = ', INFO
         Stop
      Else If ((INFO .GT. 0) .AND. (INFO .LE. N)) THEN
         Write (6, *) 'Matrix l is not invertible'
         Write (6, *) 'INFO = ', INFO
         Stop
      Else If (INFO == (N + 1)) THEN
         Write (6, *) 'Matrix l is alomost singular - be careful with solution'
         Write (6, *) 'INFO = ', INFO
         Stop
      Else
         Write (6, *) 'Something is wrong with dgbsvx'
         Write (6, *) 'INFO = ', INFO
         Stop
      End If

      !! Setting up output and moving it into the soln
      allocate (soln(1:N))
      soln(:) = X(:, 1)

      deallocate (AB, AFB, IPIV, R, C, FERR, BERR, WORK, IWORK, X)
      Return
   End Subroutine solver_banded_double_precision

End Module linear_algebra

! Code taken from:::
!  https://github.com/ECoulter/Kronecker/blob/master/README.md

! Build module in here, with Block Matrix operations.

Module Kronecker
   use type_kinds, only: dp

contains
   
   !!
   ! @brief      Takes in Matrices A(i,j),B(k,l), assumed 2D, returns Kronecker Product C(i*k,j*l)
   !             Note a non-commutative solution
   !
   ! @param      A     Matrix (n*m)
   ! @param      B     Matrix (p*q)
   ! 
   ! @return     C     Matrix (np * mq)
   !!
   Function KronProd(A, B) result(C)
      IMPLICIT NONE
      real(dp), dimension(:, :), intent(in)  :: A, B
      real(dp), dimension(:, :), allocatable :: C
      integer :: i = 0, j = 0, k = 0, l = 0
      integer :: m = 0, n = 0, p = 0, q = 0

      allocate (C(size(A, 1)*size(B, 1), size(A, 2)*size(B, 2)))
      C = 0

      !!! note sure if this can be parallesied with omp?
      Do i = 1, size(A, 1)
         Do j = 1, size(A, 2)
            n = (i - 1)*size(B, 1) + 1
            m = n + size(B, 1) - 1
            p = (j - 1)*size(B, 2) + 1
            q = p + size(B, 2) - 1
            C(n:m, p:q) = A(i, j)*B
         End Do
      End Do

   End Function KronProd

End module Kronecker
