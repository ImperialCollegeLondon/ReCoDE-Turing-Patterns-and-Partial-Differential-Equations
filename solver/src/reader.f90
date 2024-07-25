Module reader
   use type_kinds
   implicit none

   logical :: file_exists
   integer :: Time_switch
   integer :: Non_Linear_switch
   integer :: Eqn_Number
   integer :: Max_iter
   real(dp) :: Newton_Error

   ! Domain Terms
   integer :: nx, ny, nt, Domain_number
   real(dp) :: xl, xr, xhalf, yl, yr, yhalf, tl, tr
   logical :: x_grid_stretch_on, y_grid_stretch_on

contains

   Subroutine read_me
      use maths_constants, only: DiffOrder
      integer :: x_grid_stretch_on_integer, y_grid_stretch_on_integer

!!! opens settings.input file and reads data

      Inquire (file='settings.input', exist=file_exists)
      Select Case (file_exists)
      Case (.TRUE.)
      Case (.FALSE.)
         Write (6, *) 'No settings.input file - stopping program'
         Stop
      End Select

      Open (1, file='settings.input', status='old', action='read')

      Read (1, *) !first line is title
      Read (1, *) Time_switch
      Read (1, *) Non_Linear_switch
      Read (1, *) Eqn_Number
      Read (1, *) Domain_number
      Read (1, *)
      Read (1, *) ! x settings
      Read (1, *) nx
      Read (1, *) xl, xr
      Read (1, *) x_grid_stretch_on_integer, xhalf
      Read (1, *)
      Read (1, *) ! y settings
      Read (1, *) ny
      Read (1, *) yl, yr
      Read (1, *) y_grid_stretch_on_integer, yhalf
      Read (1, *)
      Read (1, *) ! t settings
      Read (1, *) nt
      Read (1, *) tl, tr
      Read (1, *)
      Read (1, *) ! general settings
      Read (1, *) DiffOrder
      Read (1, *) Newton_Error
      Read (1, *) Max_iter

      Close (1)

      Select Case (Eqn_Number)
      Case (1)
      Case (2)
      Case default
         Write (6, *) 'Eqn_Number can be 1 or 2!!! Stopping Program'
         Stop
      End Select

      Select Case (Time_switch)
      Case (1)
      Case (0)
      Case default
         Write (6, *) 'Error in 1 settings.input.... stopping program'
      End Select

      Write (6, *)

      Select Case (x_grid_stretch_on_integer)
      Case (1)
         x_grid_stretch_on = .TRUE.
         Write (6, *) 'X Grid stretching ON'
      Case (0)
         x_grid_stretch_on = .FALSE.
         Write (6, *) 'X Grid stretching OFF'
      Case default
         Write (6, *) 'Error in settings.input.... stopping program'
         Write (6, *) 'x_grid_stretch_on_integer should be 1 or 0'
      End Select

      Select Case (nx)
      Case (1:5)
         Write (6, *)
         Write (6, *) 'Error in settings.input.... stopping program'
         Write (6, *)
         Write (6, *) 'nx must be greater than 5'
         Write (6, *)
         Stop
      End Select

      Select Case (y_grid_stretch_on_integer)
      Case (1)
         y_grid_stretch_on = .TRUE.
         Write (6, *) 'Y Grid stretching ON'
      Case (0)
         y_grid_stretch_on = .FALSE.
         Write (6, *) 'Y Grid stretching OFF'
      Case default
         Write (6, *) 'Error in settings.input.... stopping program'
         Write (6, *) 'y_grid_stretch_on_integer should be 1 or 0'
      End Select

      Select Case (ny)
      Case (1:5)
         Write (6, *)
         Write (6, *) 'Error in settings.input.... stopping program'
         Write (6, *)
         Write (6, *) 'ny must be greater than 5'
         Write (6, *)
         Stop
      End Select

      Select Case (difforder)
      Case (2, 3, 4)
      Case default
         Write (6, *)
         Write (6, *) 'Error in settings.input.... stopping program'
         Write (6, *)
         Write (6, *) 'Diff order can be:::'
         Write (6, *) '   2 - second order'
         Write (6, *) 'or 3 - second order boundaries, fourth order interior'
         Write (6, *) 'or 4 - fourth order'
         Write (6, *)
         Stop
      End Select

   End Subroutine read_me

End Module reader
