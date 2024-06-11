module linear_algebra
    use equations, only : band_the_matrix
	use type_kinds, only: sp,dp
	implicit none
  	external :: sgesv,dgbsvx ! lapack linear solve

   contains

subroutine solver_single_precision(n,lda,A,NRHS,ldb,B)
!https://www.netlib.org/lapack/single/sgesv.f
    real(sp),dimension(:,:),allocatable,intent(inout) :: A
    real(sp),dimension(:,:),allocatable,intent(inout) :: B
    integer,intent(in) :: n,lda,NRHS,ldb
    integer,dimension(:),allocatable :: ipiv
    integer :: INFO

!!! A  (LDA,N)
!!! B N-by-NRHS
!!! n The number of linear equations
!!! lda  The leading dimension of the array A.
!!! nrhs The number of right hand sides,
!!! ldb The leading dimension of the array B.  LDB >= max(1,N).

	 allocate(ipiv(n))

	 call sgesv(n, NRHS, a, lda, ipiv, b, ldb, info)

    if (info /= 0) then
        print '(a, i0)', 'Error: ', info
        stop
    end if

    deallocate(ipiv)

end subroutine solver_single_precision

subroutine solver_banded_double_precision(n_input,nband,sub_diag,sup_diag,L,RHS,soln,banded)
integer,intent(in) :: n_input, nband,sub_diag,sup_diag ! dimension, second dimension, sub and super diags
real(dp),dimension(:,:),allocatable,intent(in)  :: L !(nband*n)
real(dp),dimension(:),allocatable,intent(in) :: RHS !(n)
real(dp),dimension(:),allocatable,intent(out)   :: soln !(n)
logical,intent(in) :: banded
integer :: i,j

!! copied directly from: 
!https://netlib.org/lapack/explore-html/d1/da6/group__gbsvx_ga38273d98ae4d598529fc9647ca847ce2.html#ga38273d98ae4d598529fc9647ca847ce2

    character*1 ::   fact
    character*1 ::  trans
    integer   ::  n
    integer   ::  kl
    integer   ::  ku
    integer   ::  nrhs
    real(dp), dimension(:,:), allocatable :: ab
    integer   ::  ldab
    real(dp), dimension(:,:), allocatable  ::   afb
    integer   ::  ldafb
    integer, dimension(:), allocatable  ::   ipiv
    character*1 ::   equed
    real(dp), dimension(:), allocatable  ::  r
    real(dp), dimension(:), allocatable  ::  c
    real(dp), dimension(:,:), allocatable :: b
    integer   ::  ldb
    real(dp), dimension(:,:), allocatable :: x
    integer   ::  ldx
    real(dp)  ::  rcond
    real(dp), dimension(:), allocatable  :: ferr
    real(dp), dimension(:), allocatable  ::  berr
    real(dp), dimension(:), allocatable  ::  work
    integer, dimension(:),allocatable  ::  iwork
    integer   ::  info        


    fact = 'N'
    trans = 'N'

    kl = sub_diag
    ku = sup_diag
    n = n_input
    nrhs = 1

select case(banded)
case(.TRUE.) 
    ldab = nband
    allocate(ab(ldab,n_input))
    ab = L

    !!! FOR TESTING
case(.FALSE.)

    call band_the_matrix(n,L,kl,ku,LDAB,AB)

end select

LDAFB = 2*KL+KU+1
allocate(AFB(LDAFB,N)) !output 

allocate(ipiv(1:N)) !output

EQUED = 'N'

allocate(R(1:N),C(1:N))

ldb = N ! Dimension of B
allocate(B(LDB,NRHS))
B(:,1) = RHS 

LDX = N
allocate(X(LDX,NRHS))

allocate(BERR(NRHS),FERR(NRHS))

allocate(work(max(1,3*N)),iwork(N))

    
  call dgbsvx(FACT,TRANS,N,KL,KU,NRHS,AB,LDAB,AFB,LDAFB,IPIV,&
  &EQUED,R,C,B,LDB,X,LDX,RCOND,FERR,BERR,WORK,IWORK,INFO)

IF (INFO==0) THEN 
  ELSE IF(INFO.LT.0) THEN
    WRITE(6,*) 'Info < 0, dgbsvx, has an illegal input value'
  ELSE IF((info.GT.0).AND.(info.LE.N)) THEN
    WRITE(6,*) 'Matrix is not invertible'
  ELSE IF (info==(N+1)) THEN
    !Write(6,*) 'matrix is alomost singular - be careful with solution'
  ELSE
  WRITE(6,*) 'Other things are wrong in finding solution. Info = ',info
END IF

allocate(soln(1:n))
do i = 1,n
    soln(i) = x(i,1)
end do


deallocate(AB,AFB,IPIV,R,C,FERR,BERR,WORK,IWORK,X)


end subroutine solver_banded_double_precision


end module linear_algebra