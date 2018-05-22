!*==selective1.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
 
     module SELECTIVE1
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     character(8), allocatable, dimension(:) :: dynsel, elcontspl, monctr,     &
               & tcelevcon, tcntr, tcyearly, tspltcntr, tsyearly
     integer, allocatable, dimension(:) :: jbmon, jsmon, kstrsplt, ncountcw,   &
          & nouts, seld, tciseg, tcjb, tcjs, tcnelev, tspltjb
     integer, allocatable, dimension(:, :) :: jstsplt, jstspltt, ncountc
     integer :: numtempc, numtsplt, tempn
     real, allocatable, dimension(:) :: nxsel, qwdfrac, tcklay, tctemp, tctend,&
                                      & tctsrt, temp2, tspltt, tstend, tstsrt, &
                                      & volm
     real :: nxtsplit, nxtstr, nxttcd, tcdfreq, tfrqtmp
     real, allocatable, dimension(:, :) :: qstrfrac, tcelev, tempcrit, volmc
     character(8) :: tempc, tspltc
!
!*** End of declarations rewritten by SPAG
!
     end module SELECTIVE1
