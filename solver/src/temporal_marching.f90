!!{Organises the temproal marching subroutines}
!!
Module temporal_marching
   use omp_lib
   use type_kinds
   use reader, only: Time_switch, Non_Linear_switch, Eqn_number
   use domain
   use equations, only: equation_runner, initial_condition,initial_condition2D
   use maths_constants, only: DiffOrder, nband, sub_diag, sup_diag
   use linear_algebra
   use Newtons_method, only: non_linear_iteration

   real(dp), dimension(:, :), allocatable :: Soln ! solution

contains

!!
! @brief      {Runs the temporal march}
!
! @return     Solves the IBVP in the form A u_xx + B u_x + C + non_linear terms = D u_t
!             Outputs to IVBP.dat
!!
   Subroutine march_runner
      integer :: i, j, jk,k
      real(dp), dimension(:,:,:), allocatable :: U_2d

      Write (6, *) 'size of x Domain::: ', nx, dxc
      Write (6, *) 'size of total Domain::: ', idim
      Write (6, *)
      Write (6, *) 'size of t Domain::: ', nt, dt
      Write (6, *)
      Write (6, *) 'order of the finite differences', DiffOrder
      Write (6, *) '... Building matrix'
      Write (6, *)
      Write (6, *) 'Starting march...'

      allocate (Soln(1:idim, 1:nt))


      !! Set the initial conditions



      Select Case(Domain_number)
      Case(1)
         Call initial_condition(nx, xdom, Soln)
      Case(2)
         Call initial_condition2D(nx,ny,xdom,ydom,idim,idim_xy,Eqn_number,Soln)
      End Select

      !! March


      Call implicit_march

      

      If (Domain_number==1)then
      !! Printing Case
      Select Case (Eqn_number)

      Case (1)

         Open (10, file='IVBP1.dat')
         Write (10, '(20000(f20.14,1x))') 0.d0, (tDom(j), j=1, nt)

         Do i = 1, nx
            Write (10, '(20000(f20.14,1x))') xDom(i), (soln(i, j), j=1, nt)
         End Do
         Close (10)
         Write (6, *) 'March completed...'
         Write (6, *)

      Case (2)

         Open (10, file='IVBP1.dat')
         Open (11, file='IVBP2.dat')
         Write (10, '(200000(f20.14,1x))') 0.d0, (tDom(j), j=1, nt)
         Write (11, '(200000(f20.14,1x))') 0.d0, (tDom(j), j=1, nt)

         Do i = 1, nx
            Write (10, '(200000(f20.14,1x))') xDom(i), (soln(2*i - 1, j), j=1, nt)
            Write (11, '(200000(f20.14,1x))') xDom(i), (soln(2*i, j), j=1, nt)
         End Do
         Close (10)
         Close (11)
         Write (6, *) 'March completed...'
         Write (6, *)

      End Select


      Else If (Domain_Number==2) then

      Select Case (Eqn_number)

      Case (1)
            Open (10, file='IBVP1.dat')
            Open (11, file='IBVPx.dat')
            Open (12, file='IBVPy.dat')

      allocate(U_2d(1:nx,1:ny,1))
      
      Do i = 1, nx
            Write (11, '(200000(f20.14,1x),e20.10)') (xdom(i),j=1,ny)
            Write (12, '(200000(f20.14,1x),e20.10)') (ydom(j),j=1,ny)
      End Do



      Do jk = 1, nt

        ! Convert the 1D solution vector into a 2D solution matrix
        do j = 1, ny
          do i = 1, nx
            k = (j-1)*nx + i
            U_2d(i, j,1) = Soln(k,jk)
          end do
        end do

         Do i = 1, nx
            Write (10, '(200000(f20.14,1x),e20.10)') (U_2d(i,j,1),j=1,ny)
         End Do

      end do

         Case (2)
            Open (9, file='IBVP1.dat')
            Open (10, file='IBVP2.dat')
            Open (11, file='IBVPx.dat')
            Open (12, file='IBVPy.dat')

      allocate(U_2d(1:nx,1:ny,2))
      
      Do i = 1, nx
            Write (11, '(200000(f20.14,1x),e20.10)') (xdom(i),j=1,ny)
            Write (12, '(200000(f20.14,1x),e20.10)') (ydom(j),j=1,ny)
      End Do


      Do jk = 1, nt

        ! Convert the 1D solution vector into a 2D solution matrix
        do j = 1, ny
          do i = 1, nx
            k = (j-1)*nx + i
            U_2d(i, j,1) = Soln(2*k-1,jk)
            U_2d(i, j,2) = Soln(2*k,jk)
          end do
        end do

         Do i = 1, nx
            Write (9, '(200000(f20.14,1x),e20.10)') (U_2d(i,j,1),j=1,ny)
            Write (10, '(200000(f20.14,1x),e20.10)') (U_2d(i,j,2),j=1,ny)
         End Do

      end do

      End Select
      End If
      Return
   End Subroutine march_runner

!!
! @brief      {Sets up and solves the implicit system
!              We want solve the system
!              RHS(u^{i+1} - u^{i}) = Lu^{i+1}  + boundary conditions
!
!              Outside of the boundaries:::
!              Ltemp = (RHS - dt*L)
!              Ltemp u^{i+1} = RHS u^{i}
!
!              At the boundaries we have Lu = RHS (independent of spatial march)
!
!
!              Ltemp   Implicit marching operator
!              temp    temporary RHS
!              X       temporary solution
!
!              }
!
! @return     Soln - contains the soln to the PDE
!!
   Subroutine implicit_march

      real(dp), dimension(:), allocatable :: U, temp ! temporary vectors
      real(dp), dimension(:, :), allocatable :: L !! operator from equations
      real(dp), dimension(:, :), allocatable :: L_March, L_march_unbanded !! implicit marching operator
      real(dp), dimension(:), allocatable :: RHS !! right hand side of equation

      integer :: i, j, k, differ, iteration

!!!! Build the opeator
      allocate (temp(1:idim), U(idim))
      allocate (L_March(nband, 1:idim),L_march_unbanded(1;idim,1:idim))

      !! Obtain the original operator L
      Call equation_runner(L, RHS)

      !! Builds the operator
      !!! Interior points - L and L_March are in banded form - row sub_diag+1 is where the diagonals live
      
      Do i = 1,idim
         L_march_unbanded(i, i) = RHS(i)
      End Do
      Call band_the_matrix(idim, L_march_unbanded, sub_diag, sup_diag, nband, L_March)   
      Deallocate(L_march_unbanded)      


      L_March = L_March - dt*L



      !!!! March the operator
      !!! Begin the march in time
      Do j = 2, nt



      !!!! Set the boundaries
      Select Case (Domain_number)
      Case(1)

      !! set the boundaries

         temp(1:Eqn_number) = RHS(1:Eqn_number)
         temp(idim - Eqn_number:idim) = RHS(idim - Eqn_number:idim)

      !! set the interior points
         temp(1 + Eqn_number:idim - Eqn_number) = &
                   &RHS(1 + Eqn_number:idim - Eqn_number)*Soln(1 + Eqn_number:idim - Eqn_number, j - 1)

      Case(2)

      !! set the interior points
         temp(:) = RHS(:)*Soln(:, j - 1)


         !!! Set the boundary points
         differ = Eqn_number*nx
         temp(1:differ) = RHS(1:differ)
         temp(((ny-1)*nx)*Eqn_number:idim) = RHS(((ny-1)*nx)*Eqn_number:idim)
         
         If (Eqn_number == 1) then

         Do k = 1, ny-1
            temp(1 + k*nx) = RHS(1+k*nx) 
            temp((nx) + k*nx) = RHS((nx) + k*nx)
         End do

         Else if (Eqn_number == 2) then
         Do k = 1, ny-1
            temp(1 + differ*k) = RHS(1 + differ*k)
            temp(2 + differ*k) = RHS(2 + differ*k)
            temp(2*nx-1 + k*differ) = RHS(2*nx-1 + k*differ)
            temp(2*nx + k*differ) = RHS(2*nx + k*differ) 
         End do
         End if

      End Select



         !!! solve the system
         Call solver_banded_Double_precision(idim, nband, sub_diag, sup_diag, L_March, temp, U)



         Select Case (Non_Linear_switch)
         Case (1)
            Call non_linear_iteration(L_March, temp, U, iteration)
         Case Default
         End Select

         Soln(:, j) = U(:)
         Write (6, *) 'Position done::', j, U(20)
      End Do

      deallocate (temp, U, L_March)

   End Subroutine implicit_march

End Module temporal_marching

