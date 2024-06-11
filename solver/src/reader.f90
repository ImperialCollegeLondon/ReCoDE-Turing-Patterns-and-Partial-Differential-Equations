module reader
implicit none

Logical :: file_exists

contains 

subroutine read_me
use domain, only : nx,xl,xr,x_singular_pertubation_location,xhalf,x_grid_strech_on
use maths_constants, only : DiffOrder
integer :: x_grid_strech_on_integer

!!! opens settings.input file and reads data

inquire (file='settings.input', exist=file_exists)
	select case(file_exists)
	case(.TRUE.)
	case(.FALSE.)
	WRITE(6,*) 'No settings.input file - stopping program'
	stop
	END Select

	open(1,file='settings.input',status='old',action='read')
	
	read(1,*) !first line is title
	read(1,*) ! x settings
	read(1,*) nx
	read(1,*) xl,xr
	read(1,*) x_grid_strech_on_integer,x_singular_pertubation_location,xhalf
	read(1,*) ! general settings
	read(1,*) DiffOrder

	close(1)

	select case(x_grid_strech_on_integer)
	case(1)
		x_grid_strech_on = .TRUE.
		WRITE(6,*) 'Grid streching ON'
	case(0)
		x_grid_strech_on = .FALSE.
		WRITE(6,*) 'Grid streching OFF'
	case default
		WRITE(6,*) 'Error in settings.input.... stopping program'
		WRITE(6,*) 'x_grid_strech_on_integer should be 1 or 0'
	end select


	select case(nx)
	case (7 :)
	case default
		WRITE(6,*)
		WRITE(6,*) 'Error in settings.input.... stopping program'
		WRITE(6,*)
		WRITE(6,*) 'nx must be greater than 6' 
		WRITE(6,*)
		stop
	end select


	select case(difforder)
	case (2,3,4)
	case default
		WRITE(6,*)
		WRITE(6,*) 'Error in settings.input.... stopping program'
		WRITE(6,*)
		WRITE(6,*) 'Diff order can be:::' 
		WRITE(6,*) '   2 - second order'
		WRITE(6,*) 'or 3 - second order boundaries, fourth order interior'
		WRITE(6,*) 'or 4 - fourth order'
		WRITE(6,*)
		stop
	end select


end subroutine read_me

end module reader