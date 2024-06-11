module test_collection_lapack
    use linear_algebra, only: solver_single_precision,solver_banded_double_precision
    use testdrive, only: error_type, unittest_type, new_unittest, check
    implicit none
    private

    public :: collect_tests_lapack

contains

    !> Collect all exported unit tests
    subroutine collect_tests_lapack(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
        new_unittest("test lapack solver single percision - one vector ", test_solver_sp_1),  &
        new_unittest("test lapack solver single percision - two vector ", test_solver_sp_2),  &
        new_unittest("test lapack banded solver single percision - one vector ", test_solver_dp_1)  &
        ]

    end subroutine collect_tests_lapack

    !> Check lapack works correct for linear solve 1
    subroutine test_solver_dp_1(error)
      use type_kinds, only: dp

      type(error_type), allocatable, intent(out) :: error
       real(dp),dimension(:,:),allocatable :: amat
       real(dp),dimension(:),allocatable :: bmat,soln  
       real(dp),dimension(:),allocatable  :: output    
       real(dp) :: check_value,error_value
      integer  :: i, j
      integer  :: n,lda,ldb,nrhs,sub_diag,sup_diag

    n = 2
    lda=n
    ldb=n
    nrhs=1

    sub_diag = 1
    sup_diag = 1

    allocate(amat(LDA,N),bmat(n),soln(n))

    amat = reshape([ 2., 3., 2., 2. ], [ 2, 2 ])
    bmat(:) = [ 3., 6. ]
    soln(:) = [ 3., -1.5 ]

    call solver_banded_double_precision(n,lda,sub_diag,sup_diag,amat,bmat,output,.FALSE.)


!!! is the error within float limit
		error_value = 1./(10.)**10.
   
    do i = 1,n
    	check_value = abs(soln(i) - output(i) )
    	!call check(error, check_value<error_value )
    	call check(error, check_value<error_value )
    end do

    deallocate(soln,bmat,amat,output)

!!  Amat: 
!    2.00000000       3.00000000    
!    2.00000000       2.00000000    
!  bmat:    3.00000000      -1.50000012    
! solution: x vect  3.000000, -1.500000

    end subroutine test_solver_dp_1

    subroutine test_solver_sp_1(error)
      use type_kinds, only: sp

      type(error_type), allocatable, intent(out) :: error
       real(sp),dimension(:,:),allocatable :: amat
       real(sp),dimension(:,:),allocatable :: bmat,soln     
       real(sp) :: check_value,error_value
      integer  :: i, j
      integer  :: n,lda,ldb,nrhs

    n = 2
    lda=n
    ldb=n
    nrhs=1

    allocate(amat(LDA,N),bmat(ldb,nrhs),soln(ldb,nrhs))

    amat = reshape([ 2., 3., 2., 2. ], [ 2, 2 ])
    bmat(:,1) = [ 3., 6. ]
    soln(:,1) = [ 3., -1.5 ]

    call solver_single_precision(n,lda,amat,nrhs,ldb,bmat)

!!! is the error within float limit
        error_value = 1.d-6
    do i = 1,n
        check_value = abs(soln(i,1) - bmat(i,1) )
        !call check(error, check_value<error_value )
        call check(error, check_value<error_value )
    end do

    deallocate(soln,bmat,amat)

!!  Amat: 
!    2.00000000       3.00000000    
!    2.00000000       2.00000000    
!  bmat:    3.00000000      -1.50000012    
! solution: x vect  3.000000, -1.500000

    end subroutine test_solver_sp_1


        subroutine test_solver_sp_2(error)
      use type_kinds, only: sp

      type(error_type), allocatable, intent(out) :: error
       real(sp),dimension(:,:),allocatable :: amat
       real(sp),dimension(:,:),allocatable :: bmat,soln     
       real(sp) :: check_value,error_value
      integer  :: i, j
      integer  :: n,lda,ldb,nrhs

    n = 2
    lda=n
    ldb=n
    nrhs=2

    allocate(amat(LDA,N),bmat(ldb,nrhs),soln(ldb,nrhs))

    amat = reshape([ 2., 3., 2., 2. ], [ 2, 2 ])
    bmat =  reshape([ 3., 6.,3., 6. ], [ 2, 2 ]) 
    soln =  reshape([ 3., -1.5 ,3., -1.5  ], [ 2, 2 ]) 

    call solver_single_precision(n,lda,amat,nrhs,ldb,bmat)

!!! is the error within float limit
		error_value = 1.d-6
  
  	do j = 1,n
    do i = 1,n
    	check_value = abs(soln(i,j) - bmat(i,j) )
    	call check(error, check_value<error_value )
    end do
    end do

    deallocate(soln,bmat,amat)

!!  Amat: 
!    2.00000000       3.00000000    
!    2.00000000       2.00000000    
!  bmat:    3.00000000      -1.50000012    
! solution: x vect  3.000000, -1.500000

    end subroutine test_solver_sp_2

end module test_collection_lapack