module test_collection1
    use testdrive, only: error_type, unittest_type, new_unittest, check
    implicit none
    private

    public :: collect_tests1

contains

    !> Collect all exported unit tests
    subroutine collect_tests1(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
        new_unittest("test domain builder (set_up_domain) ", check_domain_builder),  &
        new_unittest("test banded matrix calculations  ", check_banded_matrix), &
        new_unittest("Test the second order BVP Solver 4th order accurate ", check_BVP_solver_4th_order), &
        new_unittest("Test the second order BVP Solver 3rd order accurate ", check_BVP_solver_3rd_order), &
        new_unittest("Test the second order BVP Solver 2nd order accurate ", check_BVP_solver_2nd_order) &
        &]

    end subroutine collect_tests1


    subroutine check_domain_builder(error)
      use type_kinds, only: dp
      use domain, only : set_up_domain
      use maths_constants, only : pi

      type(error_type), allocatable, intent(out) :: error
       real(dp),dimension(:),allocatable :: dv,dcv,x,c,metric1,metric1sq,metric2
       real(dp) :: check_value,error_value
       real(dp)  :: l,r,dc,check1,check2,check3,check4,check5
       integer  :: n,i

       error_value = 1.d-8

       n = 57
       l = 34.7d0
       r = pi/4.d0

!!! with this function the end points should match
!!! cdom should be 1 and zero
!!! dc should = dcv
       call set_up_domain(n,l,r,dc,dv,dcv,x,c,.FALSE.,'l',0.d0,metric1,metric1sq,metric2)

       check1 = abs(l-x(1))
       call check(error, check1<error_value )

       check2 = abs(r-x(n))
       call check(error, check2<error_value )

       check3 = abs(0.d0-c(1))
       call check(error, check3<error_value )

       check4 = abs(1.d0-c(n))
       call check(error, check4<error_value )     

       do i = 1,n
        check5 = abs(dc-dcv(n)) 
        call check(error, check5<error_value ) 
       end do

!!!! need to put tests in for the metrics - with streching on and off


       deallocate(dv,dcv,x,c,metric1,metric1sq,metric2)
    end subroutine check_domain_builder


    subroutine check_banded_matrix(error)
      use type_kinds, only: dp
      use equations, only : band_the_matrix
      use maths_constants, only : pi

      type(error_type), allocatable, intent(out) :: error
       real(dp),dimension(:,:),allocatable :: A,AB,AB_test
       real(dp) :: error_value
       integer  :: n,kl,ku,ldab,ldab_test
       integer  :: i,j

       error_value = 1.d-8

       n = 5
       kl = 2
       ku = 1
       LDAB = KL + KU + 1
       allocate(A(n,n),AB(ldab,n))

       A = reshape([ 2.d0, 3.d0, 2.d0, 0.d0, 0.d0, &
                   & 1.d0, 3.d0, 2.d0, 4.d0, 0.d0, &
                   & 0.d0, 2.d0, 2.d0, 12.d0, 8.d0, &
                   & 0.d0, 0.d0, 4.d0, 13.d0, 17.d0, &
                   & 0.d0, 0.d0, 0.d0, 7.d0, 1.d0], &
                   &[ 5, 5])

       call band_the_matrix(n,A,kl,ku,LDAB_test,AB_test)

    AB = transpose(reshape([0.d0,1.d0,2.d0,4.d0,7.d0,&
                &2.d0,3.d0,2.d0,13.d0,1.d0, &
                &3.d0,2.d0,12.d0,17.d0,0.d0,&
                &2.d0,4.d0,8.d0,0.d0,0.d0], &
                &[ 5, 4]))

    call check(error, ldab, ldab_test)

    do i = 1,ldab
    do j = 1,n
      call check(error, ab(i,j), ab_test(i,j))
    end do 
    end do

    deallocate(A,AB,Ab_test)

    return
    end subroutine check_banded_matrix



    subroutine check_BVP_solver_4th_order(error)
      use type_kinds, only: dp
      use equations_TEST, only : build_the_matrix_test,L,nband,RHS
      use maths_constants, only : diff_initialisation,pi,ex,sub_diag,sup_diag,difforder
      use reader, only : read_me
      use domain, only : initial_domain_settings,nx,xdom,xcdom,dxc,dxcsq,xmetric1,xmetric1sq,xmetric2
      use linear_algebra, only : solver_banded_double_precision

      type(error_type), allocatable, intent(out) :: error

       real(dp) :: error_value,test_value
       real(dp),dimension(:),allocatable :: X
       integer :: j


      call read_me  ! opens settings.input and reads
      difforder = 4   
      call diff_initialisation  ! sets the finite difference coefficients
      call initial_domain_settings !builds the domain and computational domains
      call build_the_matrix_test(nx,dxc,dxcsq,xcdom,xmetric1,xmetric1sq,xmetric2)
      call solver_banded_double_precision(nx,nband,sub_diag,sup_diag,L,RHS,X,.true.)

       error_value = 1.d-6

    do j = 1,nx
      !WRITE(6,*) X(j),ex**xdom(j)
      test_value = abs(x(j)-ex**xdom(j))
      !call check(error, test_value < error_value)
    end do 

    deallocate(X,L,xdom,xcdom)

        WRITE(6,'(a,1x,e12.4)') '     Error Value::',error_value  

      return
    end subroutine check_bvp_solver_4th_order

    subroutine check_BVP_solver_3rd_order(error)
      use type_kinds, only: dp
      use equations_TEST, only : build_the_matrix_test,L,nband,RHS
      use maths_constants, only : diff_initialisation,pi,ex,sub_diag,sup_diag,difforder
      use reader, only : read_me
      use domain, only : initial_domain_settings,nx,xdom,xcdom,dxc,dxcsq,xmetric1,xmetric1sq,xmetric2
      use linear_algebra, only : solver_banded_double_precision

      type(error_type), allocatable, intent(out) :: error

       real(dp) :: error_value,test_value
       real(dp),dimension(:),allocatable :: X
       integer :: j


      call read_me  ! opens settings.input and reads
      difforder = 3  
      call diff_initialisation  ! sets the finite difference coefficients
      call initial_domain_settings !builds the domain and computational domains
      call build_the_matrix_test(nx,dxc,dxcsq,xcdom,xmetric1,xmetric1sq,xmetric2)
      call solver_banded_double_precision(nx,nband,sub_diag,sup_diag,L,RHS,X,.true.)

       error_value = 1.d-4

    do j = 1,nx
      !WRITE(6,*) X(j),ex**xdom(j)
      test_value = abs(x(j)-ex**xdom(j))
      call check(error, test_value < error_value)
    end do 

    deallocate(X,L,xdom,xcdom)


    WRITE(6,'(a,1x,e12.4)') '     Error Value::',error_value  


    return
    end subroutine check_BVP_solver_3rd_order


        subroutine check_BVP_solver_2nd_order(error)
      use type_kinds, only: dp
      use equations_TEST, only : build_the_matrix_test,L,nband,RHS
      use maths_constants, only : diff_initialisation,pi,ex,sub_diag,sup_diag,difforder
      use reader, only : read_me
      use domain, only : initial_domain_settings,nx,xdom,xcdom,dxc,dxcsq,xmetric1,xmetric1sq,xmetric2
      use linear_algebra, only : solver_banded_double_precision

      type(error_type), allocatable, intent(out) :: error

       real(dp) :: error_value,test_value
       real(dp),dimension(:),allocatable :: X
       integer :: j


      call read_me  ! opens settings.input and reads
      difforder = 2   
      call diff_initialisation  ! sets the finite difference coefficients
      call initial_domain_settings !builds the domain and computational domains
      call build_the_matrix_test(nx,dxc,dxcsq,xcdom,xmetric1,xmetric1sq,xmetric2)
      call solver_banded_double_precision(nx,nband,sub_diag,sup_diag,L,RHS,X,.true.)

       error_value = 1.d-3

    do j = 1,nx
      !WRITE(6,*) X(j),ex**xdom(j)
      test_value = abs(x(j)-ex**xdom(j))
      call check(error, test_value < error_value)
    end do 

    deallocate(X,L,xdom,xcdom)

        WRITE(6,'(a,1x,e12.4)') '     Error Value::',error_value  

    return
    end subroutine check_BVP_solver_2nd_order


end module test_collection1