!*==surfhe.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module SURFHE
     use PREC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(r8), allocatable, dimension(:) :: afw, bfw, cfw, cshe, et, lat,      &
           & longit, phi0, rb, rc, re, shade, wind, windh, wsc
     integer, allocatable, dimension(:) :: iwind
     real(r8) :: phiset, rhowcp
     real :: PREC
     real :: r8
     logical, allocatable, dimension(:) :: rh_evap
!
!*** End of declarations rewritten by SPAG
!
     end module SURFHE                                     !MLM 08/12/05
