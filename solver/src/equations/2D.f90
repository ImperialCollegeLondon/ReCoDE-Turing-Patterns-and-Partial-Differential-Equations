Module Equation_2D
use type_kinds, only: dp
use Kronecker, only : KronProdMod, KronProd
use equations_builder
use equations_definition
Contains

Subroutine equation_setup2D(L, RHS, nx, ny, idim_xy, idim, Eqn_number,&
                & dxc, dxcsq, xdom, xcdom,xmetric1, xmetric1sq, xmetric2,&
                & dyc, dycsq, ydom, ycdom,ymetric1, ymetric1sq, ymetric2)

      real(dp), dimension(:, :), allocatable, intent(out) :: L !! unbanded form matrix
      real(dp), dimension(:), allocatable, intent(out) :: RHS !! right hand side of equation
      integer, intent(in) :: nx, ny, idim_xy, idim, eqn_number
      real(dp), intent(in) :: dxc, dxcsq,dyc, dycsq
      real(dp), dimension(:), allocatable, intent(in) :: xdom, xcdom,xmetric1, xmetric1sq, xmetric2
      real(dp), dimension(:), allocatable, intent(in) :: ydom, ycdom,ymetric1, ymetric1sq, ymetric2

      real(dp), dimension(:,:,:), allocatable :: Ix,Iy
      real(dp), dimension(:,:,:), allocatable :: l1x, l1y
      real(dp), dimension(:,:), allocatable :: R1
      real(dp), dimension(:,:,:), allocatable :: l2x, l2y
      real(dp), dimension(:,:), allocatable ::  R2
      real(dp),dimension(:),allocatable :: Rtemp1, Rtemp2
      real(dp),dimension(:,:,:),allocatable :: ltemp1, ltemp2

      real(dp),dimension(3,3) :: hhh
      real(dp),dimension(9) :: lll

      integer :: i, j, k, i1,i2, j1, j2, differ

      allocate(RHS(idim),L(idim,idim))
      allocate(Ix(nx,nx,nx),Iy (ny,ny,ny))

      Ix = 0.d0; Iy = 0.d0


        Do i = 1,nx
          Ix(i,i,i) = 1.d0
        End do

        Do i = 1,ny
          Iy(i,i,i) = 1.d0
        End do


      Select Case (Eqn_number)
      Case(1)

      !!! note in this case idim_xy = idim

      allocate(Ltemp1(idim_xy,idim_xy,ny+nx))
      Ltemp1 = 0.d0

      Call build_the_matrix2D(nx, dxc, dxcsq, xcdom, xmetric1, xmetric1sq, xmetric2, &
                              &ny, dyc, dycsq, ycdom, ymetric1, ymetric1sq, ymetric2, &
                              &1, L1x, L1y, R1)


        !! Build Ltemp1 operator and Rtemp1 with the Kronecker product
        Do j = 2,ny-1
          Ltemp1(:,:,j) = KronProd(Iy(:,:,j),l1x(:,:,j)) 
        End do

        Do j = 1,nx
           Ltemp1(:,:,j+ny) =  KronProd(l1y(:,:,j),Ix(:,:,j))
        End do

        Do i = 1,ny-2
          Ltemp1(nx*i+1,:,ny+1) = 0.d0
          Ltemp1(nx*i+nx,:,ny+nx) = 0.d0
        End do

        !Rtemp1 = KronProd(Iy,R1) 
        RHS = reshape(R1, (/idim_xy/))


        Do i = 1,idim
        Do j = 1,idim
          L(i,j) = Sum(Ltemp1(i,j,:))
        End do
        End Do

        deallocate(ltemp1, l1x,l1y, R1)

      Case(2)

      allocate(Ltemp1(idim_xy,idim_xy,ny+nx),Ltemp2(idim_xy,idim_xy,ny+nx))
      Ltemp1 = 0.d0; Ltemp2 = 0.d0

      Call build_the_matrix2D(nx, dxc, dxcsq, xcdom, xmetric1, xmetric1sq, xmetric2, &
                              &ny, dyc, dycsq, ycdom, ymetric1, ymetric1sq, ymetric2, &
                              &1, L1x, L1y, R1)

      Call build_the_matrix2D(nx, dxc, dxcsq, xcdom, xmetric1, xmetric1sq, xmetric2, &
                              &ny, dyc, dycsq, ycdom, ymetric1, ymetric1sq, ymetric2, &
                              &2, L2x, L2y, R2)


        !! Build Ltemp1 operator and Rtemp1 with the Kronecker product
        Do j = 2,ny-1
          Ltemp1(:,:,j) = KronProd(Iy(:,:,j),l1x(:,:,j)) 
          Ltemp2(:,:,j) = KronProd(Iy(:,:,j),l2x(:,:,j)) 
        End do

        Do j = 1,nx
           Ltemp1(:,:,j+ny) =  KronProd(l1y(:,:,j),Ix(:,:,j))
           Ltemp2(:,:,j+ny) =  KronProd(l2y(:,:,j),Ix(:,:,j))
        End do

        Do i = 1,ny-2
          Ltemp1(nx*i+1,:,ny+1) = 0.d0; Ltemp1(nx*i+nx,:,ny+nx) = 0.d0
          Ltemp2(nx*i+1,:,ny+1) = 0.d0; Ltemp2(nx*i+nx,:,ny+nx) = 0.d0
        End do

        !Rtemp1 = KronProd(Iy,R1) 
        RHS = reshape(R1, (/idim_xy/))


        Do i = 1,idim
        Do j = 1,idim
          L(i,j) = Sum(Ltemp1(i,j,:))
        End do
        End Do

        deallocate(ltemp1, l1x,l1y, R1)
 
      End Select

      deallocate(Ix, Iy)

      Return
      End Subroutine equation_setup2D


 !!
   ! @brief      The subroutine sets up the equation in the form L * u = RHS
   !
   ! @param      n_x          dimension of problem
   ! @param      ch_x         step size
   ! @param      chsq_x       step size squared
   ! @param      cdom_x       computational domain
   ! @param      metric1_x    The metric 1 - determined by domains
   ! @param      metric1sq_x  The metric 1 sq - determined by domains
   ! @param      metric2_x    The metric 2 - determined by domains
   ! @param      which_equation    (integer - asks which equation you want to discretise)
   !
   ! @return     L     left hand side of the equation - returns a matrix in banded form (nband * n_x)
   ! @return     RHS   RHS of the equation - returns a vector (n_x)
  !!
   Subroutine build_the_matrix2D(n_x, ch_x, chsq_x, cdom_x, metric1_x, metric1sq_x, metric2_x,&
                               &n_y, ch_y, chsq_y, cdom_y, metric1_y, metric1sq_y, metric2_y,& 
                               &which_equation, Lx, Ly, RHS)

   !!! sets up a matrix equation in the form L * u = RHS
   !!! D is the RHS, Lx is the matrix

      integer, intent(in) :: n_x !! size of domain
      real(dp), intent(in) :: ch_x, chsq_x !computational distance/squared
      real(dp), dimension(:), allocatable, intent(in) :: cdom_x ! computational domain

    !! grid streching metric terms
      real(dp), dimension(:), allocatable, intent(in) :: metric1_x, metric1sq_x, metric2_x


      integer, intent(in) :: n_y !! size of domain
      real(dp), intent(in) :: ch_y, chsq_y !computational distance/squared
      real(dp), dimension(:), allocatable, intent(in) :: cdom_y ! computational domain

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
      real(dp) :: blank

    !!! allocate
      allocate (Lx(n_x, n_x, n_y), Ly(n_y, n_y, n_x), RHS(n_x, n_y))




      Do jj = 1, n_y

      allocate (A(n_x), B(n_x), C(n_x), D(n_x))
      allocate (At(n_x), Bt(n_x), Ct(n_x), Dt(n_x))
      allocate (Lxtemp(n_x, n_x))


      Select Case (which_equation)
      Case(1)
         !!! Boundaries
         Call equation1_BC_X_Bot(cdom_x(1), cdom_y(jj), At(1), Bt(1), blank, blank, Ct(1), Dt(1))
         Call equation1_BC_X_Top(cdom_x(n_x), cdom_y(jj), At(n_x), Bt(n_x),  blank, blank,  Ct(n_x), Dt(n_x))
         !!! Interior

         !$omp Parallel Do
         Do i = 2, n_x - 1
            Call equation1_linear(cdom_x(i), cdom_y(jj), At(i), Bt(i),  blank, blank,Ct(i), Dt(i))
         End do
         !$omp End Parallel Do
      Case(2)
         !!! Boundaries
         Call equation2_BC_X_Bot(cdom_x(1), cdom_y(jj), At(1), Bt(1),  blank, blank,Ct(1), Dt(1))
         Call equation2_BC_X_Top(cdom_x(n_x), cdom_y(jj), At(n_x), Bt(n_x), blank, blank,Ct(n_x), Dt(n_x))
         !!! Interior

         !$omp Parallel Do
         Do i = 2, n_x - 1
            Call equation2_linear(cdom_x(i), cdom_y(jj), At(i), Bt(i), blank, blank,Ct(i), Dt(i))
         End do
         !$omp End Parallel Do
      End Select


    !! Apply the correct scallings
      !$omp Parallel Do
      Do i = 1, n_x
         Call scales(A(i), B(i), C(i), D(i), At(i), Bt(i), Ct(i), Dt(i), ch_x, chsq_x, metric1_x(i), metric1sq_x(i), metric2_x(i))
      End Do
      !$omp End Parallel Do

    !!! Derivative runner moves the coefficients into a banded matrix L
      Call derivative_runner(n_x, A, B, C, Lxtemp)

      Lx(:,:,jj) = Lxtemp
      RHS(:, jj) = D
      deallocate (Dt, Ct, Bt, At, D, C, B, A, Lxtemp)


      End Do

      Do jj = 1, n_x

      allocate (A(n_y), B(n_y), C(n_y), D(n_y))
      allocate (At(n_y), Bt(n_y), Ct(n_y), Dt(n_y))
      allocate (Lytemp(n_y, n_y))

      Select Case (which_equation)
      Case(1)

         !!! Boundaries
         Call equation1_BC_Y_Bot(cdom_x(jj), cdom_y(1), blank, blank, At(1), Bt(1), Ct(1), Dt(1))
         Call equation1_BC_Y_Top(cdom_x(jj), cdom_y(n_y), blank, blank, At(n_y), Bt(n_y), Ct(n_y), Dt(n_y))
         !!! Interior

         !$omp Parallel Do
         Do i = 2, n_y - 1
            Call equation1_linear(cdom_x(jj), cdom_y(i), blank, blank,  At(i), Bt(i), Ct(i), Dt(i))
         End do
         !$omp End Parallel Do
        
      Case(2)
      
         !!! Boundaries
         Call equation2_BC_Y_Bot(cdom_x(jj), cdom_y(1), blank, blank, At(1), Bt(1), Ct(1), Dt(1))
         Call equation2_BC_Y_Top(cdom_x(jj), cdom_y(n_y), blank, blank, At(n_y), Bt(n_y), Ct(n_y), Dt(n_y))
         !!! Interior

         !$omp Parallel Do
         Do i = 2, n_y - 1
            Call equation2_linear(cdom_x(jj), cdom_y(i), blank, blank,  At(i), Bt(i), Ct(i), Dt(i))
         End do
         !$omp End Parallel Do
      End Select   


    !! Apply the correct scallings
      !$omp Parallel Do
      Do i = 1, n_y
         Call scales(A(i), B(i), C(i), D(i), At(i), Bt(i), Ct(i), Dt(i), ch_y, chsq_y, metric1_y(i), metric1sq_y(i), metric2_y(i))
      End Do
      !$omp End Parallel Do

    !!! Derivative runner moves the coefficients into a banded matrix L
      Call derivative_runner(n_y, A, B, C, Lytemp)

      Ly(:,:,jj) = Lytemp
    
    !! set the final output for RHS
      
    !  Do i = 1,n_y
    !  Write(6,*) jj, i, Ct(i), C(i),cdom_x(jj)
    !  end do

      !!! Note that the boundary conditions at the corners are governed by BC_Y 
      RHS(jj, 1) = D(1)
      RHS(jj, n_y) = D(n_y)

      deallocate (Dt, Ct, Bt, At, D, C, B, A, Lytemp)

      End Do

      Return
   End Subroutine build_the_matrix2D



   !!
   ! @brief      {Sets the initial condition for the temporal march}
   !
   ! @param      n     dimension of computational spatial domain
   ! @param      dom  computational spatial domain
   !
   ! @return      soln  solution
   !!
   Subroutine initial_condition2D(nx, ny, xdom,ydom,idim,idim_xy,Eqn_number, soln)
      integer, intent(in) :: nx,ny, idim, idim_xy,Eqn_number
      real(dp), dimension(:), allocatable, intent(in) :: xdom,ydom
      real(dp), dimension(:, :), allocatable, intent(inout) :: soln
      real(dp), dimension(:, :), allocatable :: temp1,temp2
      real(dp), dimension(:), allocatable :: Rtemp1, Rtemp2
      integer :: i, j


      Select Case (Eqn_number)

      Case (1)

      allocate(temp1(nx,ny))

         !$omp Parallel Do
         Do i = 1, nx
         Do j = 1, ny
            Call equation1_initial_condition(xdom(i), ydom(j), temp1(i, j))
         End Do
         End Do
         !$omp End Parallel Do


        soln(:,1) = Reshape(temp1, (/idim_xy/))

        open (10,file='t.dat')
        Do i = 1, nx
          WRite(10,*) (temp1(i,j),j=1,ny)
        End do
        close(10)

        deallocate(temp1)

      Case (2)
       allocate(temp1(nx,ny),Rtemp1(idim_xy))
       allocate(temp2(nx,ny),Rtemp2(idim_xy))

         !$omp Parallel Do
         Do i = 1, nx
         Do j = 1, ny
            Call equation1_initial_condition(xdom(i), ydom(j), temp1(i, j))
            Call equation2_initial_condition(xdom(i), ydom(j), temp2(i, j))
         End Do
         End Do
         !$omp End Parallel Do
        

        Rtemp1 = Reshape(temp1, (/idim_xy/))
        Rtemp2 = Reshape(temp2, (/idim_xy/))

        Do i = 1,idim_xy
          soln(2*i-1,1) = Rtemp1(i)
          soln(2*i,1) = Rtemp2(i)
        End do
        
        deallocate(Rtemp1,temp1)
        deallocate(Rtemp2,temp2)
      End Select

   End Subroutine initial_condition2D



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
   Subroutine non_linear_setup2D (nx, ny, xcdom,ycdom,idim,idim_xy,Eqn_number,nband,diag, U, F, Fu, Fv)
      use linear_algebra, only : band_the_matrix
      integer, intent(in) :: nx,ny, idim, idim_xy,Eqn_number,nband,diag
      real(dp), dimension(:), allocatable, intent(in) :: xcdom,ycdom
      real(dp), dimension(:), allocatable, intent(in) ::  U
      real(dp), dimension(:), allocatable, intent(inout) ::  F
      real(dp), dimension(:, :), allocatable, intent(inout) ::  Fu, Fv

      real(dp), dimension(:, :), allocatable ::  Fu_temp, Fv_temp
      integer :: i,j,k
      real(dp) :: blank



      allocate (Fu_temp(idim, idim), Fv_temp(idim, idim))
      Fu_temp = 0.d0; Fv_temp = 0.d0

      Select Case (Eqn_number)

      Case (1)
         !$omp Parallel Do
         
         Do i = 2, nx - 1
         Do j = 2, ny - 1
            k = (j-1)*nx + i
            Call equation1_non_linear(xcdom(i), ycdom(j), U(k), 0.d0, F(k), Fu_temp(k, k), blank)
         End Do
         End do


         !$omp End Parallel Do

      Case (2)
         !$omp Parallel Do
         Do i = 2, nx - 1
         Do j = 2, ny - 1
            k = (j-1)*nx + i
            Call equation1_non_linear(xcdom(i), ycdom(i), U(2*i - 1), U(2*i), F(2*i - 1), Fu_temp(2*i - 1, 2*i - 1), &
                                                                                          &Fv_temp(2*i - 1, 2*i))
            Call equation2_non_linear(xcdom(i), ycdom(i), U(2*i - 1), U(2*i), F(2*i), Fu_temp(2*i, 2*i - 1), &
                                                                                          &Fv_temp(2*i, 2*i))
         End do
         End do
         !$omp End Parallel Do
      End Select

      Call band_the_matrix(idim, Fu_temp, diag, diag, nband, Fu)
      Call band_the_matrix(idim, Fv_temp, diag, diag, nband, Fv)

      deallocate (Fu_temp, Fv_temp)
   End Subroutine non_linear_setup2D

End Module Equation_2D