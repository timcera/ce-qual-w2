!*==initialvelocity.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module INITIALVELOCITY
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND), allocatable, dimension(:, :) :: bsave
     real(R8KIND), allocatable, dimension(:) :: elwss, qssi, uavg
     logical :: init_vel, once_through
     logical, allocatable, dimension(:) :: loop_branch
!
!*** End of declarations rewritten by SPAG
!
     end module INITIALVELOCITY
