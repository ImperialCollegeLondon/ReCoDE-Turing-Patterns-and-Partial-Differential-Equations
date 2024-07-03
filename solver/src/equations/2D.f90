
!!This module builds the equation in two Domains if Domain_number = 2
!!
Module Equation_2D
   use type_kinds, only: dp
   use Kronecker, only: KronProd
   use equations_builder
   use equations_definition
   use omp_lib
contains

   !!
   ! @brief      Sets up the 2D discretised operator size idim*idim unbanded and RHS size idim
   !

   ! @param      nx          Size of x Domain
   ! @param      ny          Size of y Domain
   ! @param      idim_xy     nx*ny
   ! @param      idim        nx*ny*Eqn_number
   ! @param      Eqn_number  Number of equations
   ! @param      dxc         Computational x distance
   ! @param      dxcsq       Computational xsq distance
   ! @param      xDom        xDomain
   ! @param      xcDom       x computational Domain
   ! @param      xmetric1    The xmetric 1
   ! @param      xmetric1sq  The xmetric 1 sq
   ! @param      xmetric2    The xmetric 2
   ! @param      dyc         Computational y distance
   ! @param      dycsq       Computational ysq distance
   ! @param      yDom        The yDom
   ! @param      ycDom       The ycDom
   ! @param      ymetric1    The ymetric 1
   ! @param      ymetric1sq  The ymetric 1 sq
   ! @param      ymetric2    The ymetric 2

   ! @return      L           Discretised operator 
   ! @return      RHS         RHS of equation
   !!
   Subroutine equation_setup2D(L, RHS, nx, ny, idim_xy, idim, Eqn_number,&
                   & dxc, dxcsq, xDom, xcDom, xmetric1, xmetric1sq, xmetric2,&
                   & dyc, dycsq, yDom, ycDom, ymetric1, ymetric1sq, ymetric2)

      !This subroutine builds the operator using the Kroncker Product method:
      !https://en.wikipedia.org/wiki/Kronecker_product

      real(dp), dimension(:, :), allocatable, intent(out) :: L !! unbanded form matrix
      real(dp), dimension(:), allocatable, intent(out) :: RHS !! right hand side of equation
      integer, intent(in) :: nx, ny, idim_xy, idim, eqn_number
      real(dp), intent(in) :: dxc, dxcsq, dyc, dycsq
      real(dp), dimension(:), allocatable, intent(in) :: xDom, xcDom, xmetric1, xmetric1sq, xmetric2
      real(dp), dimension(:), allocatable, intent(in) :: yDom, ycDom, ymetric1, ymetric1sq, ymetric2

      real(dp), dimension(:, :, :), allocatable :: Ix, Iy,L1x, L1y, L2x, L2y, ltemp1, ltemp2
      real(dp), dimension(:, :), allocatable :: R1, L1, R2, L2
      real(dp), dimension(:), allocatable :: Rtemp1, Rtemp2

      integer :: i, j, k, i1, i2, j1, j2

      allocate (RHS(idim), L(idim, idim))
      allocate (Ix(nx, nx, nx), Iy(ny, ny, ny))

      Ix = 0.d0; Iy = 0.d0

   
      ! Build identitiy tensor where the diagonal is taken in three dimensions
      !$omp Parallel Do
      Do i = 1, nx
         Ix(i, i, i) = 1.d0
      End Do
      !$omp End Parallel Do

      !$omp Parallel Do
      Do i = 1, ny
         Iy(i, i, i) = 1.d0
      End Do
      !$omp End Parallel Do


      ! How many equations?
      Select Case (Eqn_number)
      Case (1)


      !!! note in this case idim_xy = idim

         allocate (Ltemp1(idim_xy, idim_xy, ny + nx))
         Ltemp1 = 0.d0

         !! Build the equation - we have nx L1y's and ny L1x's - this explains the tensor form of Ix and Iy
         Call build_the_matrix2D(nx, dxc, dxcsq, xcDom, xmetric1, xmetric1sq, xmetric2, &
                                 &ny, dyc, dycsq, ycDom, ymetric1, ymetric1sq, ymetric2, &
                                 &1, L1x, L1y, R1)

        !! Use the Kroncker product to build operator
        !! first index 2 to ny - 1 to account for y BCs
         !$omp Parallel Do
         Do j = 2, ny - 1
            Ltemp1(:, :, j) = KronProd(Iy(:, :, j), L1x(:, :, j))
         End Do
         !$omp End Parallel Do

         Do j = 1, nx
            Ltemp1(:, :, j + ny) = KronProd(L1y(:, :, j), Ix(:, :, j))
         End Do
         !$omp End Parallel Do

         !! Turn off some y terms to account for x BCs
         !$omp Parallel Do
         Do i = 1, ny - 2
            Ltemp1(nx*i + 1, :, ny + 1) = 0.d0
            Ltemp1(nx*i + nx, :, ny + nx) = 0.d0
         End Do
         !$omp End Parallel Do


         ! Add up each of the contributions
         !$omp Parallel Do
         Do i = 1, idim_xy
         Do j = 1, idim_xy
            L1(i, j) = Sum(Ltemp1(i, j, :))
         End Do
         End Do
         !$omp End Parallel Do


         !! Reshape the RHS into a vector size idim_xy
         RHS = reshape(R1, (/idim_xy/))

         deallocate (ltemp1, L1x, L1y, R1)





      Case (2)

         allocate (Ltemp1(idim_xy, idim_xy, ny + nx), Ltemp2(idim_xy, idim_xy, ny + nx))
         Ltemp1 = 0.d0; Ltemp2 = 0.d0

         Call build_the_matrix2D(nx, dxc, dxcsq, xcDom, xmetric1, xmetric1sq, xmetric2, &
                                 &ny, dyc, dycsq, ycDom, ymetric1, ymetric1sq, ymetric2, &
                                 &1, L1x, L1y, R1)

         Call build_the_matrix2D(nx, dxc, dxcsq, xcDom, xmetric1, xmetric1sq, xmetric2, &
                                 &ny, dyc, dycsq, ycDom, ymetric1, ymetric1sq, ymetric2, &
                                 &2, L2x, L2y, R2)

        !! Build Ltemp1 operator and Rtemp1 with the Kronecker product
         !$omp Parallel Do
         Do j = 2, ny - 1
            Ltemp1(:, :, j) = KronProd(Iy(:, :, j), L1x(:, :, j))
            Ltemp2(:, :, j) = KronProd(Iy(:, :, j), L2x(:, :, j))
         End Do
         !$omp End Parallel Do

         !! Build Ltemp2 operator and Rtemp1 with the Kronecker product
         !$omp Parallel Do
         Do j = 1, nx
            Ltemp1(:, :, j + ny) = KronProd(L1y(:, :, j), Ix(:, :, j))
            Ltemp2(:, :, j + ny) = KronProd(L2y(:, :, j), Ix(:, :, j))
         End Do
         !$omp End Parallel Do

         !! Boundary conditions
         !$omp Parallel Do
         Do i = 1, ny - 2
            Ltemp1(nx*i + 1, :, ny + 1) = 0.d0; Ltemp1(nx*i + nx, :, ny + nx) = 0.d0
            Ltemp2(nx*i + 1, :, ny + 1) = 0.d0; Ltemp2(nx*i + nx, :, ny + nx) = 0.d0
         End Do
         !$omp End Parallel Do

         !Reshape RHS
         Rtemp1 = reshape(R1, (/idim_xy/))
         Rtemp2 = reshape(R2, (/idim_xy/))

         allocate (L1(idim_xy, idim_xy), L2(idim_xy, idim_xy))
         
         !Sum the operators
         !$omp Parallel Do
         Do i = 1, idim_xy
         Do j = 1, idim_xy
            L1(i, j) = Sum(Ltemp1(i, j, :))
            L2(i, j) = Sum(Ltemp2(i, j, :))
         End Do
         End Do
         !$omp End Parallel Do

         !  Set the Ltemp and RHS - two equations
         !$omp Parallel Do
         Do i = 1, idim_xy
            RHS(2*i - 1) = Rtemp1(i)
            RHS(2*i) = Rtemp2(i)
            Do j = 1, idim_xy
               i1 = 2*i - 1
               i2 = 2*i - 0
               j1 = 2*j - 1
               j2 = 2*j - 0

               L(i1, j1) = L1(i, j)
               L(i2, j2) = L2(i, j)

            End Do
         End Do
         !$omp End Parallel Do

         deallocate (Rtemp1, ltemp1, L1x, L1y, R1, L1)
         deallocate (Rtemp2, ltemp2, L2x, L2y, R2, L2)

      End Select

      deallocate (Ix, Iy)

      Return
   End Subroutine equation_setup2D

 !!

   !!
   ! @brief        The subroutine sets up equations in the form L * u = RHS where we discretised in x and y
   !
   ! @param      n_x             The size of x
   ! @param      ch_x            x computational distance
   ! @param      chsq_x          ch_x^2
   ! @param      cDom_x          computational xdomain
   ! @param      metric1_x       The metric 1 x
   ! @param      metric1sq_x     The metric 1 sq x
   ! @param      metric2_x       The metric 2 x
   ! @param      n_y             The size of y
   ! @param      ch_y            y computatioanl distance
   ! @param      chsq_y          ch_y^2
   ! @param      cDom_y          computational ydomain
   ! @param      metric1_y       The metric 1 y
   ! @param      metric1sq_y     The metric 1 sq y
   ! @param      metric2_y       The metric 2 y
   ! @param      which_equation  which equation 1 or 2
   ! 
   ! @return      Lx              { parameter_description }
   ! @return      Ly              { parameter_description }
   ! @return      RHS             The rhs
   !
   !!
   Subroutine build_the_matrix2D(n_x, ch_x, chsq_x, cDom_x, metric1_x, metric1sq_x, metric2_x,&
                               &n_y, ch_y, chsq_y, cDom_y, metric1_y, metric1sq_y, metric2_y,&
                               &which_equation, Lx, Ly, RHS)

      integer, intent(in) :: n_x !! size of Domain
      real(dp), intent(in) :: ch_x, chsq_x !computational distance/squared
      real(dp), dimension(:), allocatable, intent(in) :: cDom_x ! computational Domain

    !! grid streching metric terms
      real(dp), dimension(:), allocatable, intent(in) :: metric1_x, metric1sq_x, metric2_x

      integer, intent(in) :: n_y !! size of Domain
      real(dp), intent(in) :: ch_y, chsq_y !computational distance/squared
      real(dp), dimension(:), allocatable, intent(in) :: cDom_y ! computational Domain

    !! grid streching metric terms
      real(dp), dimension(:), allocatable, intent(in) :: metric1_y, metric1sq_y, metric2_y

    !! which equation?
      integer, intent(in) :: which_equation

    !! outputs
      real(dp), dimension(:, :, :), allocatable, intent(out) :: Lx, Ly !! banded form matrix
      real(dp), dimension(:, :), allocatable, intent(out) :: RHS !! right hand side of equation

      !!!! relevent vectors/infomation
      real(dp), dimension(:), allocatable :: A, B, C, D
      real(dp), dimension(:), allocatable :: At, Bt, Ct, Dt
      real(dp), dimension(:, :), allocatable :: Lxtemp, Lytemp
      integer :: i, j, jj
      real(dp) :: blank1 ,blank2

    !! Here we run through every point - building x operators and y operators
    !! We have ny x operators size Lx and nx y operators size ny
      allocate (Lx(n_x, n_x, n_y), Ly(n_y, n_y, n_x), RHS(n_x, n_y))

      !!! Set up the x operators at every point in y
      Do jj = 1, n_y

         allocate (A(n_x), B(n_x), C(n_x), D(n_x))
         allocate (At(n_x), Bt(n_x), Ct(n_x), Dt(n_x))
         allocate (Lxtemp(n_x, n_x))

         Select Case (which_equation)

         !! First fill the vectors of each equation - very similar to correpsonding equations in equations.f90
         Case (1)
         
         !!! Boundaries
            Call equation1_BC_X_Bot(cDom_x(1), cDom_y(jj), At(1), Bt(1), blank1, blank2, Ct(1), Dt(1))
            Call equation1_BC_X_Top(cDom_x(n_x), cDom_y(jj), At(n_x), Bt(n_x), blank1, blank2, Ct(n_x), Dt(n_x))
        
         !!! Interior
            !$omp Parallel Do
            Do i = 2, n_x - 1
               Call equation1_linear(cDom_x(i), cDom_y(jj), At(i), Bt(i), blank1, blank2, Ct(i), Dt(i))
            End Do
            !$omp End Parallel Do
         Case (2)
         
         !!! Boundaries
            Call equation2_BC_X_Bot(cDom_x(1), cDom_y(jj), At(1), Bt(1), blank1, blank2, Ct(1), Dt(1))
            Call equation2_BC_X_Top(cDom_x(n_x), cDom_y(jj), At(n_x), Bt(n_x), blank1, blank2, Ct(n_x), Dt(n_x))
         
         !!! Interior
            !$omp Parallel Do
            Do i = 2, n_x - 1
               Call equation2_linear(cDom_x(i), cDom_y(jj), At(i), Bt(i), blank1, blank2, Ct(i), Dt(i))
            End Do
            !$omp End Parallel Do
         End Select

    !! Apply the correct scallings
         !$omp Parallel Do
         Do i = 1, n_x
            Call scales(A(i), B(i), C(i), D(i), At(i), Bt(i), Ct(i), Dt(i), &
                      &ch_x, chsq_x, metric1_x(i), metric1sq_x(i), metric2_x(i))
         End Do
         !$omp End Parallel Do

    !!! Derivative runner moves the coefficients into a unbanded matrix 
         Call derivative_runner(n_x, A, B, C, Lxtemp)

         Lx(:, :, jj) = Lxtemp
         
         !! We set the RHS here
         RHS(:, jj) = D
         deallocate (Dt, Ct, Bt, At, D, C, B, A, Lxtemp)

      End Do

      !! Set up y operators
      
      Do jj = 1, n_x

         allocate (A(n_y), B(n_y), C(n_y), D(n_y))
         allocate (At(n_y), Bt(n_y), Ct(n_y), Dt(n_y))
         allocate (Lytemp(n_y, n_y))

         Select Case (which_equation)
         Case (1)

         !!! Boundaries
            Call equation1_BC_Y_Bot(cDom_x(jj), cDom_y(1), blank1, blank2, At(1), Bt(1), Ct(1), Dt(1))
            Call equation1_BC_Y_Top(cDom_x(jj), cDom_y(n_y), blank1, blank2, At(n_y), Bt(n_y), Ct(n_y), Dt(n_y))
         
         !!! Interior
            !$omp Parallel Do
            Do i = 2, n_y - 1
               Call equation1_linear(cDom_x(jj), cDom_y(i), blank1, blank2, At(i), Bt(i), Ct(i), Dt(i))
            End Do
            !$omp End Parallel Do

         Case (2)

         !!! Boundaries
            Call equation2_BC_Y_Bot(cDom_x(jj), cDom_y(1), blank1, blank2, At(1), Bt(1), Ct(1), Dt(1))
            Call equation2_BC_Y_Top(cDom_x(jj), cDom_y(n_y), blank1, blank2, At(n_y), Bt(n_y), Ct(n_y), Dt(n_y))
         
         !!! Interior
            !$omp Parallel Do
            Do i = 2, n_y - 1
               Call equation2_linear(cDom_x(jj), cDom_y(i), blank1, blank2, At(i), Bt(i), Ct(i), Dt(i))
            End Do
            !$omp End Parallel Do
         End Select

    !! Apply the correct scallings
         !$omp Parallel Do
         Do i = 1, n_y
            Call scales(A(i), B(i), C(i), D(i), At(i), Bt(i), Ct(i), Dt(i), &
               &ch_y, chsq_y, metric1_y(i), metric1sq_y(i), metric2_y(i))
         End Do
         !$omp End Parallel Do

    !!! Derivative runner moves the coefficients into a banded matrix
         Call derivative_runner(n_y, A, B, C, Lytemp)

         Ly(:, :, jj) = Lytemp

      !!! Note that the boundary conditions at the corners are governed by BC_Y
         RHS(jj, 1) = D(1)
         RHS(jj, n_y) = D(n_y)

         deallocate (Dt, Ct, Bt, At, D, C, B, A, Lytemp)

      End Do

      Return
   End Subroutine build_the_matrix2D

   !!
   ! @brief      Sets the 2D initial condition for the temporal march
   !
   ! @param      nx          size of x
   ! @param      ny          size of y
   ! @param      xDom        The x dom  !! NOT COMPUTATIONAL DOMAIN
   ! @param      yDom        The y dom  !! NOT COMPUTATIONAL DOMAIN
   ! @param      idim        nx*ny*Eqn_number
   ! @param      idim_xy     nx*ny
   ! @param      Eqn_number  number of equations
   !
   ! @return     soln        The solution
   !!
   Subroutine initial_condition2D(nx, ny, xDom, yDom, idim, idim_xy, Eqn_number, soln)
      integer, intent(in) :: nx, ny, idim, idim_xy, Eqn_number
      real(dp), dimension(:), allocatable, intent(in) :: xDom, yDom
      real(dp), dimension(:, :), allocatable, intent(inout) :: soln
      real(dp), dimension(:, :), allocatable :: temp1, temp2
      real(dp), dimension(:), allocatable :: Rtemp1, Rtemp2
      integer :: i, j

      Select Case (Eqn_number)

      Case (1)

         allocate (temp1(nx, ny))

         !!! temp1 is a nx * ny matrix - we fill it and then flatten it with reshape
         !$omp Parallel Do
         Do i = 1, nx
         Do j = 1, ny
            Call equation1_initial_condition(xDom(i), yDom(j), temp1(i, j))
         End Do
         End Do
         !$omp End Parallel Do

         !! set the solution as the flattened value
         soln(:, 1) = Reshape(temp1, (/idim_xy/))

         deallocate (temp1)

      Case (2)
         allocate (temp1(nx, ny), Rtemp1(idim_xy))
         allocate (temp2(nx, ny), Rtemp2(idim_xy))

         !! Call both initial conditions
         !$omp Parallel Do
         Do i = 1, nx
         Do j = 1, ny
            Call equation1_initial_condition(xDom(i), yDom(j), temp1(i, j))
            Call equation2_initial_condition(xDom(i), yDom(j), temp2(i, j))
         End Do
         End Do
         !$omp End Parallel Do

         !! Re shape and set
         Rtemp1 = Reshape(temp1, (/idim_xy/))
         Rtemp2 = Reshape(temp2, (/idim_xy/))

         !$omp Parallel Do
         Do i = 1, idim_xy
            soln(2*i - 1, 1) = Rtemp1(i)
            soln(2*i, 1) = Rtemp2(i)
         End Do
         !$omp End Parallel Do

         deallocate (Rtemp1, temp1)
         deallocate (Rtemp2, temp2)
      End Select

   End Subroutine initial_condition2D


   
   !!
   ! @brief      Sets up the non-linear components of the differential equation
   !
   ! @param      nx          size of x
   ! @param      ny          size of y
   ! @param      xcDom       The x computational dom
   ! @param      ycDom       The y computatioanl dom
   ! @param      idim        nx*ny*Eqn_number
   ! @param      idim_xy     nx*ny
   ! @param      Eqn_number  Number of equations
   ! @param      nband       Band of matrix
   ! @param      diag        Diagonals of matrix
   ! @param      U           Current State/Soln 
   ! 
   ! @return      F          Non-liner function F(U)
   ! @return      Fu         u derivative of F - matrix dimension (nband idim)
   ! @return      Fv         v derivative of F - matrix dimension (nband idim)
   !
   !!
   Subroutine non_linear_setup2D(nx, ny, xcDom, ycDom, idim, idim_xy, Eqn_number, nband, diag, U, F, Fu, Fv)
      use linear_algebra, only: band_the_matrix
      integer, intent(in) :: nx, ny, idim, idim_xy, Eqn_number, nband, diag
      real(dp), dimension(:), allocatable, intent(in) :: xcDom, ycDom
      real(dp), dimension(:), allocatable, intent(in) ::  U
      real(dp), dimension(:), allocatable, intent(inout) ::  F
      real(dp), dimension(:, :), allocatable, intent(inout) ::  Fu, Fv

      real(dp), dimension(:, :), allocatable ::  Fu_temp, Fv_temp
      integer :: i, j, k
      real(dp) :: blank

      allocate (Fu_temp(idim, idim), Fv_temp(idim, idim))
      Fu_temp = 0.d0; Fv_temp = 0.d0

      Select Case (Eqn_number)

      Case (1)
         !$omp Parallel Do

         Do i = 2, nx - 1
         Do j = 2, ny - 1
            !! term k redistributes the nx * ny grid into a nx*ny vector - similar to reshape
            k = (j - 1)*nx + i
            Call equation1_non_linear(xcDom(i), ycDom(j), U(k), 0.d0, F(k), Fu_temp(k, k), blank)
         End Do
         End Do

         !$omp End Parallel Do

      Case (2)
         ! When there are two equations need to split up each term
         !$omp Parallel Do
         Do i = 2, nx - 1
         Do j = 2, ny - 1
            k = (j - 1)*nx + i
            Call equation1_non_linear(xcDom(i), ycDom(j), U(2*k - 1), U(2*k), F(2*k - 1),  Fu_temp(2*k - 1, 2*k - 1), &
                                                                                          &Fv_temp(2*k - 1, 2*k))
            Call equation2_non_linear(xcDom(i), ycDom(j), U(2*k - 1), U(2*k), F(2*k), Fu_temp(2*k, 2*k - 1), &
                                                                                          &Fv_temp(2*k, 2*k))
         End Do
         End Do
         !$omp End Parallel Do
      End Select

      !! band the matrices
      Call band_the_matrix(idim, Fu_temp, diag, diag, nband, Fu)
      Call band_the_matrix(idim, Fv_temp, diag, diag, nband, Fv)

      deallocate (Fu_temp, Fv_temp)
   End Subroutine non_linear_setup2D

End Module Equation_2D
