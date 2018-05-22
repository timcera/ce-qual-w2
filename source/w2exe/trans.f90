!*==trans.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module TRANS
     use PREC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(r8), allocatable, dimension(:, :) :: adx, adz, at, ct, dt, dx, dz,   &
           & dzq, vt
     real(r8), pointer, dimension(:, :) :: cnew, cold, ssb, ssk
     real :: PREC
     real :: r8
     real(r8), allocatable, dimension(:) :: theta
!
!*** End of declarations rewritten by SPAG
!
     end module TRANS
