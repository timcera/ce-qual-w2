!*==surfhe.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module SURFHE
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND), allocatable, dimension(:) :: afw, bfw, cfw, cshe, et, lat,      &
           & longit, phi0, rb, rc, re, shade, wind, windh, wsc
     integer, allocatable, dimension(:) :: iwind
     real(R8KIND) :: phiset, rhowcp
     logical, allocatable, dimension(:) :: rh_evap
!
!*** End of declarations rewritten by SPAG
!
     end module SURFHE                                     !MLM 08/12/05
