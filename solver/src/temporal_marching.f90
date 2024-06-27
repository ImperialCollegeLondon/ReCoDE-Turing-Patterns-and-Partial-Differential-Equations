Module temporal_marching
   use type_kinds
   use reader
   use domain
   use maths_constants
   use equations
   use linear_algebra


   real(dp), dimension(:), allocatable :: X
   real(dp), dimension(:,:), allocatable :: Soln
   real(dp), dimension(:, :), allocatable :: L,Lmarch !! banded form matrix
   real(dp), dimension(:), allocatable :: RHS,temp !! right hand side of equation


contains
Subroutine linear_implicit_march
   integer :: i,j
      external :: DGBMV

      Write (6, *) 'size of x domain::: ', nx,dxc
      Write (6, *) 'size of t domain::: ', nt,dt
      Write (6, *) 'order of the finite differences', DiffOrder
      Write (6, *) '... Building matrix'
      Write (6, *)

      allocate(Soln(1:nx,1:nt),temp(1:nx),X(nx))
      allocate(Lmarch(1:nband,1:nx))

      Call initial_condition(nx,xcdom,Soln)
      Call build_the_matrix(nx, dxc, dxcsq, xcdom, xmetric1, xmetric1sq, xmetric2, 1, L, RHS)

      
      WRITE(6,*) RHS(1),RHS(nx)


      !! Builds the operator
      Lmarch(sub_diag+1,:) = RHS(:)
      Lmarch = Lmarch - dt*L
      do j = 1, sub_diag
        Lmarch(sub_diag + 2 - j, j) = L(sub_diag + 2 - j, j)
        Lmarch(sub_diag + j, nx + 1- j) = L(sub_diag + j, nx + 1- j)
      end do

      do j = 2,nt
         temp(1) = RHS(1)
         temp(2:nx-1) = RHS(2:nx-1)*Soln(2:nx-1,j-1)
         temp(nx) = RHS(nx)
         Call solver_banded_double_precision(nx, nband, sub_diag, sup_diag, Lmarch, temp, X)
         Soln(:,j) = X(:)
         !Call DGBMV('N', nx, nx, sub_diag, sub_diag, 1.d0, L, nband, temp, 1, 0.d0, X, 1)
         !soln(:,j) = soln(:,j-1) + dt*(X(:) - RHS(:))
      end do


open(10,file='1.dat')
WRITE(10,*) 0.d0, (tdom(j),j=1,nt)
do i = 1,nx
   WRITE(10,*) xdom(i), (soln(i,j),j=1,nt)
end do

End Subroutine linear_implicit_march

End Module temporal_marching