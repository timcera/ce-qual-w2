!*==trans.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module TRANS
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND), allocatable, dimension(:, :) :: adx, adz, at, ct, dt, dx, dz,   &
           & dzq, vt
     real(R8KIND), pointer, dimension(:, :) :: cnew, cold, ssb, ssk
     real(R8KIND), allocatable, dimension(:) :: theta
!
!*** End of declarations rewritten by SPAG
!
     end module TRANS
