!*==macrophytec.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module MACROPHYTEC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real, allocatable, dimension(:) :: anorm, armac, cddrag, dwsa, dwv, exm,  &
                                      & lrpmac, mbmp, mc, mg, mhsc, mhsn, mhsp,&
                                      & mk1, mk2, mk3, mk4, mm, mmax, mn, mp,  &
                                      & mpom, mr, msat, mt1, mt2, mt3, mt4,    &
                                      & nsed, o2mg, o2mr, psed
     real, allocatable, dimension(:, :) :: bic, cw, macmbrs, macmbrt, macwbci, &
       & ssmacmb
     character(10), allocatable, dimension(:, :) :: conv2
     real, pointer, dimension(:, :) :: domp, domr, ldommac, lpommac, nh4mg,    &
                                     & nh4mr, po4mg, po4mr, rpommac, ticmc
     real, allocatable, dimension(:, :, :) :: gammaj, mac, mact, mactrm,       &
       & mactrmf, mactrmr, mclim, mmr, mnlim, mplim, mrr, smac, smact
     logical, allocatable, dimension(:) :: kticol
     real, allocatable, dimension(:, :, :, :) :: macrc, macrm, macss, mgr,     &
       & mllim, smacrc, smacrm
     logical, allocatable, dimension(:, :) :: macrophyte_calc, print_macrophyte
     logical :: macrophyte_on
     character(8), allocatable, dimension(:, :) :: macwbc, mprwbc
     character(10), allocatable, dimension(:, :, :, :) :: mlfpr
!
!*** End of declarations rewritten by SPAG
!
     end module MACROPHYTEC                                             ! cb 8/24/15
