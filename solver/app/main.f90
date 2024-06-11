program main
  !use type_kinds, only: sp,dp
  use reader, only : read_me
  use domain
  use maths_constants
  use equations
  use linear_algebra, only : solver_banded_double_precision
  implicit none

  real(dp),dimension(:),allocatable :: X
  integer :: i

  WRITE(6,*)
  WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(6,*) '!!!!!! PDE SOLVER !!!!!!!'
  WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(6,*)
  
  !!!! call intitial subroutines

  call read_me  ! opens settings.input and reads 
  call diff_initialisation  ! sets the finite difference coefficients
  call initial_domain_settings !builds the domain and computational domains

  WRITE(6,*) 'size of x domain::: ', nx
  WRITE(6,*) 'order of the finite differences', DiffOrder

  WRITE(6,*) '... Building matrix'
  call build_the_matrix(nx,dxc,dxcsq,xcdom,xmetric1,xmetric1sq,xmetric2)

  call solver_banded_double_precision(nx,nband,sub_diag,sup_diag,L,RHS,X,.true.)
  deallocate(L)
  
  WRITE(6,*)
  do i = 1,nx
    WRITE(6,*) x(i),ex**xdom(i)
  end do

end program main