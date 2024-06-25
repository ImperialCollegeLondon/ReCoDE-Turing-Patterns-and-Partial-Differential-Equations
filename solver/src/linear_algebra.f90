!!
!{ Sets up the linear algebra solving Subroutines with Lapack integration
! }
!!
Module linear_algebra
   use type_kinds, only: sp, dp
   implicit none
   external :: dgbsvx ! lapack double precision linear solve

contains
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
    !! Note that in dgbsvx B and X can be a matrices. We however do not need this

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
      Do i = 1, N
         soln(i) = X(i, 1)
      End Do

      deallocate (AB, AFB, IPIV, R, C, FERR, BERR, WORK, IWORK, X)
      Return
   End Subroutine solver_banded_double_precision

End Module linear_algebra
