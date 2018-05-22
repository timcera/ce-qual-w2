!*==selective1usgs.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
!***********************************************************************************************************************************
!**                                                S E L E C T I V E   I N I T                                                    **
!***********************************************************************************************************************************
 
     module SELECTIVE1USGS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     character(8), allocatable, dimension(:) :: dynsel, elcontspl, tcelevcon,  &
               & tcntr, tcyearly, tsdynsel, tspltcntr, tsyearly
     real, allocatable, dimension(:, :) :: estrsav, qstrfrac, tcelev, tempcrit,&
       & tsdepth, tsmaxflow, tsmaxhead, tsminfrac, tsminhead, volmc
     real, allocatable, dimension(:) :: ewdsav, maxfrac1, maxfrac2, minfrac1,  &
                                      & minfrac2, nxsel, nxtssel, qwdfrac,     &
                                      & splt2e, splt2t, tcklay, tctemp, tctend,&
                                      & tctsrt, temp2, tspltt, tstemp2, tstend,&
                                      & tstsrt, volm
     integer, allocatable, dimension(:, :) :: jstsplt, ncountc, tsprior
     integer, allocatable, dimension(:) :: kstrsplt, ncountcw, nout0, nout1,   &
          & nout2, nouts, seld, tciseg, tcjb, tcjs, tcnelev, tspltjb, tsseld
     integer :: ng1, ng2, numtempc, numtsplt, tempn
     logical, allocatable, dimension(:, :) :: no_flow, str_active
     real :: nxtsplit, nxtstr, nxttcd, qfrac1, qfrac2, sum_minfrac1,           &
           & sum_minfrac2, tcdfreq, tfrqtmp, tsconv, tspltfreq
     logical, allocatable, dimension(:) :: share_flow, wd_active
     character(8) :: tempc, tspltc
     character(5), allocatable, dimension(:, :) :: tstype
!
!*** End of declarations rewritten by SPAG
!
     end module SELECTIVE1USGS
