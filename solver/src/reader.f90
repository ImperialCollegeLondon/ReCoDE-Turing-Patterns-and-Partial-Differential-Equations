Module reader
   implicit none

   logical :: file_exists
   integer :: Time_switch

contains

   Subroutine read_me
      use domain, only: nx, xl, xr, xhalf, x_grid_strech_on, nt,tl,tr
      use maths_constants, only: DiffOrder
      integer :: x_grid_strech_on_integer

!!! opens settings.input file and reads data

      Inquire (file='settings.input', exist=file_exists)
      Select Case (file_exists)
      Case (.TRUE.)
      Case (.FALSE.)
         Write (6, *) 'No settings.input file - stopping program'
         Stop
      End Select

      open (1, file='settings.input', status='old', action='read')

      read (1, *) !first line is title
      read (1, *) Time_switch
      read (1, *) ! x settings
      read (1, *) nx
      read (1, *) xl, xr
      read (1, *) x_grid_strech_on_integer, xhalf
      read (1, *) ! t settings
      read (1, *) nt
      read (1, *) tl, tr
      read (1, *) ! general settings
      read (1, *) DiffOrder

      close (1)

      Select Case (Time_switch)
      Case (1)
         Write (6, *) 'Solving parabolic PDE in form  A u_xx + B u_x + C u = D u_t'
      Case (0)
         Write (6, *) 'Solving ODE in form  A u_xx + B u_x + C u = D'
      Case default
         Write (6, *) 'Error in settings.input.... stopping program'
      End Select

      Write (6, *)

      Select Case (x_grid_strech_on_integer)
      Case (1)
         x_grid_strech_on = .TRUE.
         Write (6, *) 'Grid streching ON'
      Case (0)
         x_grid_strech_on = .FALSE.
         Write (6, *) 'Grid streching OFF'
      Case default
         Write (6, *) 'Error in settings.input.... stopping program'
         Write (6, *) 'x_grid_strech_on_integer should be 1 or 0'
      End Select

      Select Case (nx)
      Case (7:)
      Case default
         Write (6, *)
         Write (6, *) 'Error in settings.input.... stopping program'
         Write (6, *)
         Write (6, *) 'nx must be greater than 6'
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
