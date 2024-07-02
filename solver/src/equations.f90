!!{ This module sets up handles all the discretisation of the equations}
!!
Module equations
   use omp_lib
   use type_kinds, only: dp
   use equations_builder
   !! here we allowing equations to read the equation definitions that we have set up:
   use equations_definition
   use equations_definition_test
   use Equation_2D
   use maths_constants, only: sub_diag, sup_diag, nband
   use reader, only: Time_switch, Eqn_number
   use linear_algebra, only: band_the_matrix

contains

   Subroutine equation_runner(L, RHS)
      use domain
      use Kronecker, only : KronProd, KronProdMod
      real(dp), dimension(:, :), allocatable, intent(out) :: L !! banded form matrix
      real(dp), dimension(:), allocatable, intent(out) :: RHS !! right hand side of equation
      real(dp), dimension(:,:,:), allocatable :: L1x, L1y 
      real(dp), dimension(:,:), allocatable :: RHS_temp, L2x, L2y, Ix,Iy,Rx1,Ry1,rx2,ry2
      real(dp), dimension(:, :), allocatable :: Ltemp1,Ltemp2,Ltemp,rtemp1,rtemp2
      real(dp), dimension(:), allocatable :: RHSx1, RHSy1, RHSx2, RHSy2
      integer :: i,j,i1,i2,j1,j2


      Select Case (Domain_number)
      Case(1)

      Call equation_setup(Ltemp, RHS, nx, idim, nband, sub_diag, dxc, dxcsq, xdom, xcdom,&
                   &xmetric1, xmetric1sq, xmetric2)

      Case(2)

      Call equation_setup2D(Ltemp, RHS, nx, ny, idim_xy, idim, Eqn_number,&
                & dxc, dxcsq, xdom, xcdom,xmetric1, xmetric1sq, xmetric2,&
                & dyc, dycsq, ydom, ycdom,ymetric1, ymetric1sq, ymetric2)
      End Select

      open(10,file='f.dat')
      do i = 1,idim
        Write(10,'(100000(f9.3,1x))') Ltemp(i,:)
      end do
      close(10)


      Call band_the_matrix(idim, Ltemp, sub_diag, sup_diag, nband, L)      
      deallocate(Ltemp)




      Return
   End Subroutine equation_runner

  !!
   ! @brief      This subroutine sets up the linear equation in the form LHS * u = RHS
   !             Depending on how many equations, u may contain two different variables
   !
   ! @return     L    The linear spatial operator of the system. Built out of Eqn_number equations
   ! @return     RHS  Right hand side of the equations
   !
   !!
   Subroutine equation_setup(L, RHS, n, dim, band, diag, dc, dcsq, dom, cdom,&
                &metric1, metric1sq, metric2)

      real(dp), dimension(:, :), allocatable, intent(out) :: L !! banded form matrix
      real(dp), dimension(:), allocatable, intent(out) :: RHS !! right hand side of equation
      integer, intent(in) :: n, dim, band, diag
      real(dp), intent(in) :: dc,dcsq
      real(dp), dimension(:),allocatable,intent(in) :: dom, cdom, metric1, metric1sq, metric2

      real(dp), dimension(:, :), allocatable :: L1, L2, Ltemp !! banded form matrix
      real(dp), dimension(:), allocatable :: R1, R2 !! right hand side of equation
      integer :: i, j, i1, i2, j1, j2



      !!! L is the banded matrix, RHS is a vector
      allocate (L(1:band, 1:dim), RHS(dim))

      !!! Ltemp is the unbanded matrix
      allocate (Ltemp(1:dim, 1:dim))

      Select Case (Eqn_number)

      Case (1)
      !!! Build matrix sets up each equation
         Call build_the_matrix(n, dc, dcsq, cdom, metric1, metric1sq, metric2, 1, L1, R1)

         Ltemp = L1
         RHS = R1

         deallocate (L1, R1)
      Case (2)
         Call build_the_matrix(n, dc, dcsq, cdom, metric1, metric1sq, metric2, 1, L1, R1)
         Call build_the_matrix(n, dc, dcsq, cdom, metric1, metric1sq, metric2, 2, L2, R2)

         ! Here we fill Ltemp with the two linear equations L1 and L2
        !! Wonder how to parallize this?

         Do i = 1, n
           Do j = 1, n
              i1 = 2*i - 1
              i2 = 2*i - 0
              j1 = 2*j - 1
              j2 = 2*j - 0
              Ltemp(i1, j1) = L1(i, j)
              Ltemp(i2, j2) = L2(i, j)
           End Do
           RHS(i1) = R1(i)
           RHS(i2) = R2(i)
         End Do

         deallocate (L1, R1, L2, R2)
      End Select

      L = Ltemp

      !! Band the matrix
!      Call band_the_matrix(dim, Ltemp, diag, diag, band, L)
      deallocate (ltemp)

   End Subroutine equation_setup

  !!
   ! @brief      The subroutine sets up the equation in the form L * u = RHS
   !
   ! @param      n          dimension of problem
   ! @param      ch         step size
   ! @param      chsq       step size squared
   ! @param      cdom       computational domain
   ! @param      metric1    The metric 1 - determined by domains
   ! @param      metric1sq  The metric 1 sq - determined by domains
   ! @param      metric2    The metric 2 - determined by domains
   ! @param      which_equation    (integer - asks which equation you want to discretise)
   !
   ! @return     L     left hand side of the equation - returns a matrix in banded form (nband * n)
   ! @return     RHS   RHS of the equation - returns a vector (n)
  !!
   Subroutine build_the_matrix(n, ch, chsq, cdom, metric1, metric1sq, metric2, which_equation, L, RHS)

   !!! sets up a matrix equation in the form L * u = RHS
   !!! D is the RHS, Lx is the matrix

      integer, intent(in) :: n !! size of domain
      real(dp), intent(in) :: ch, chsq !computational distance/squared
      real(dp), dimension(:), allocatable, intent(in) :: cdom ! computational domain

    !! grid streching metric terms
      real(dp), dimension(:), allocatable, intent(in) :: metric1, metric1sq, metric2

    !! which equation?
      integer, intent(in) :: which_equation

    !! outputs
      real(dp), dimension(:, :), allocatable, intent(out) :: L !! banded form matrix
      real(dp), dimension(:), allocatable, intent(out) :: RHS !! right hand side of equation

      !!!! relevent vectors/infomation
      real(dp), dimension(:), allocatable :: A, B, C, D
      real(dp), dimension(:), allocatable :: At, Bt, Ct, Dt
      integer :: i, j
      real(dp) :: blank

    !!! allocate
      allocate (L(n, n), RHS(n))
      allocate (A(n), B(n), C(n), D(n))
      allocate (At(n), Bt(n), Ct(n), Dt(n))

    !! At,Bt,Ct and Dt form the coefficients of the discretised equation in the form Atu'' + Btu' + Ctu = Dt
    !! A,B,C and D are then the corrected coefficients given grid streching transfomations

    !! here we use OpenMP to parallelise the Do Loops - no infomation is shared between loops

      Select Case (which_equation)

      !!! Equation 1
      Case (1)

         !!! Boundaries
         Call equation1_BC_X_Bot(cdom(1), 0.d0, At(1), Bt(1), blank, blank, Ct(1), Dt(1))
         Call equation1_BC_X_Top(cdom(n), 0.d0, At(n), Bt(n), blank, blank, Ct(n), Dt(n))
         !!! Interior

         !$omp Parallel Do
         Do i = 2, n - 1
            Call equation1_linear(cdom(i), 0.d0, At(i), Bt(i), blank, blank, Ct(i), Dt(i))
         End do
         !$omp End Parallel Do

      !!! Equation 2
      Case (2)

         !!! Boundaries
         Call equation2_BC_X_Bot(cdom(1), 0.d0, At(1), Bt(1), blank, blank, Ct(1), Dt(1))
         Call equation2_BC_X_Top(cdom(n), 0.d0, At(n), Bt(n), blank, blank, Ct(n), Dt(n))
         !!! Interior

         !$omp Parallel Do
         Do i = 2, n - 1
            Call equation2_linear(cdom(i), 0.d0, At(i), Bt(i), blank, blank, Ct(i), Dt(i))
         End do
         !$omp End Parallel Do

      !!! Equation TEST
      Case (0)

         !!! Boundaries
         Call equation1_BC_Bot_test(cdom(1), At(1), Bt(1), Ct(1), Dt(1))
         Call equation1_BC_Top_test(cdom(n), At(n), Bt(n), Ct(n), Dt(n))
         !!! Interior

         !$omp Parallel Do
         Do i = 2, n - 1
            Call equation1_test(cdom(i), At(i), Bt(i), Ct(i), Dt(i))
         End do
         !$omp End Parallel Do

      Case Default

         Write (6, *) 'Equation Error: which_equation in equations.f90 should be 1 or 2'

      End Select

    !! Apply the correct scallings
      !$omp Parallel Do
      Do i = 1, n
         Call scales(A(i), B(i), C(i), D(i), At(i), Bt(i), Ct(i), Dt(i), ch, chsq, metric1(i), metric1sq(i), metric2(i))
      End Do
      !$omp End Parallel Do

    !!! Derivative runner moves the coefficients into a banded matrix L
      Call derivative_runner(n, A, B, C, L)
    !! set the final output for RHS
      RHS = D

      deallocate (Dt, Ct, Bt, At, D, C, B, A)

      Return
   End Subroutine build_the_matrix


 


   !!
   ! @brief      {Sets the initial condition for the temporal march}
   !
   ! @param      n     dimension of computational spatial domain
   ! @param      dom  computational spatial domain
   !
   ! @return      soln  solution
   !!
   Subroutine initial_condition(n, dom, soln)
      integer, intent(in) :: n
      real(dp), dimension(:), allocatable, intent(in) :: dom
      real(dp), dimension(:, :), allocatable, intent(inout) :: soln
      integer :: i



      Select Case (Eqn_number)

      Case (1)
         !$omp Parallel Do
         Do i = 1, n
            Call equation1_initial_condition(dom(i), 0.d0, soln(i, 1))
         End Do
         !$omp End Parallel Do

      Case (2)
         !$omp Parallel Do
         Do i = 1, n
            Call equation1_initial_condition(dom(i), 0.d0, soln(2*i - 1, 1))
            Call equation2_initial_condition(dom(i), 0.d0, soln(2*i, 1))
         End Do
         !$omp End Parallel Do
      End Select

   End Subroutine initial_condition

   !!
   ! @brief      {Sets up the non-linear components of the differential equation}
   !
   ! @param      n     dimension of the domain (not idim - nx)
   ! @param      cdom  computational domain
   !
   ! @return     U     Solution at previous state (split into two in cases)
   ! @return     F     Non linear term (Vector)
   ! @return     Fu    u derivative of non-linear term (matrix)
   ! @return     Fv    v derivative of non-linear term (matrix)
   !
   !!
   Subroutine non_linear_setup(n, cdom, U, F, Fu, Fv)
      use domain, only: idim
      integer, intent(in) :: n
      real(dp), dimension(:), allocatable, intent(in) :: cdom, U
      real(dp), dimension(:), allocatable, intent(inout) ::  F
      real(dp), dimension(:, :), allocatable, intent(inout) ::  Fu, Fv
      real(dp), dimension(:, :), allocatable ::  Fu_temp, Fv_temp
      integer :: i
      real(dp) :: blank

      allocate (Fu_temp(idim, idim), Fv_temp(idim, idim))
      Fu_temp = 0.d0; Fv_temp = 0.d0

      Select Case (Eqn_number)

      Case (1)
         !$omp Parallel Do
         Do i = 2, n - 1
            Call equation1_non_linear(cdom(i), 0.d0, U(i), 0.d0, F(i), Fu_temp(i, i), blank)
         End do
         !$omp End Parallel Do

      Case (2)
         !$omp Parallel Do
         Do i = 2, n - 1
            Call equation1_non_linear(cdom(i), 0.d0, U(2*i - 1), U(2*i), F(2*i - 1), Fu_temp(2*i - 1, 2*i - 1), &
                                                                                          &Fv_temp(2*i - 1, 2*i))
            Call equation2_non_linear(cdom(i), 0.d0, U(2*i - 1), U(2*i), F(2*i), Fu_temp(2*i, 2*i - 1), &
                                                                                          &Fv_temp(2*i, 2*i))
         End do
         !$omp End Parallel Do
      End Select

      Call band_the_matrix(idim, Fu_temp, sub_diag, sup_diag, nband, Fu)
      Call band_the_matrix(idim, Fv_temp, sub_diag, sup_diag, nband, Fv)

      deallocate (Fu_temp, Fv_temp)
   End Subroutine non_linear_setup



End Module equations

