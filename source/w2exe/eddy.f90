!*==eddy.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module EDDY
     use PREC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(r8), allocatable, dimension(:) :: arodi, e, erough, fric, strick,    &
           & tkelatprdconst, ustarbtke, wshy
     character(8), allocatable, dimension(:) :: azc, imptke
     real(r8), allocatable, dimension(:, :) :: azt, decay, dzt, fricbr
     integer, allocatable, dimension(:) :: firsti, lasti, tkebc, wallpnt
     real :: PREC
     real :: r8
     logical, allocatable, dimension(:) :: strickon, tkelatprd
     real(r8), allocatable, dimension(:, :, :) :: tke
!
!*** End of declarations rewritten by SPAG
!
     end module EDDY
