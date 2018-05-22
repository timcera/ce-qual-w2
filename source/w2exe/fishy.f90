!*==fishy.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***                                                                ****
!***        N U M E R I C A L   PARTICLE/FISH   S U R R O G A T E   ****
!***                                                                ****
!***********************************************************************
!***********************************************************************
!***********************************************************************
! Andy Goodwin F77 Version 1/2001
! Scott Wells F90 Version + Enhancements for Random Water Movement + Bug Fixes/Link to W2V3.1  1/2001
! Scott Wells Updates to latest version V3.7 5/1/2015
! Scott Wells Updates and particle tracking V4.1 8/1/2017
! Module Definitions
 
     module FISHY
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
     integer, parameter :: FPARA = 16, FFP = 4, WQP = 2
!
! Local variables
!
     real :: alphax, alphaz, aspratio, bbdysearch, botdo, bothvel, botk,       &
           & bottemp, botvvel, bxsensph, countfortruc, currtfishx, currtfishz, &
           & dbdysearch, depthint, diffromopt, distreactxdo, distreactxhv,     &
           & distreactxtp, distreactxvv, distreactzdo, distreactzhv,           &
           & distreactztp, distreactzvv, dostep1mult, dothres, dothres2,       &
           & doxweigt, dozweigt, dzsensph, epsilonrd, fage, fbdysearch, fnearx,&
           & fnearz, frdx, frdy, frdz, fromktbot, fsize, fxdospan, fxgraddo,   &
           & fxgradhv, fxgradtp, fxgradvv, fxhvspan, fxloc, fxrdspan,          &
           & fxreactdo, fxreactrnd, fxreacttmp, fxreactvel, fxsensph, fxtpspan,&
           & fxvvspan, fyloc, fyloctmp, fyvel, fzdospan, fzgraddo, fzgradhv,   &
           & fzgradtp, fzgradvv, fzhvspan, fzloc, fzrdspan, fzreactdo,         &
           & fzreactrnd, fzreacttmp, fzreactvel, fztpspan, fzvvspan,           &
           & hadepthint, hr, hvxweigt, hvzweigt, jdaydiff, lftdo, lfthvel,     &
           & lfttemp
     character(2) :: ampm
     integer :: barchrtxfn = 8003, barchrtzfn = 8004, datadebugfn = 8002,      &
              & diagfn = 8001, finalfn = 8005, initialfn = 8006
     integer :: botvvvel, dir, dnbp, fatpt, fcount, fi, fimp, fjr, fk, fkmp,   &
              & fkmptemp, fn, fnbp, fschl, fsnag, ilok, iswitch, jbp, jj,      &
              & klast, klok, ktwbf, locate, ncollector, net, nfish, nfishpcel, &
              & nfishseg, ngrpfh, nl, nr, numacoustics, numgillnets, nzones,   &
              & othfi, othfk, rd, seed, snd, surfcalc, tag, unbp, variable,    &
              & varydo, varyhvel, varytemp, varyvvel, wbrun
     real, allocatable, dimension(:, :, :) :: brchfish, flowfield, lastfield,  &
       & nodes, sndcatch, wqfield
     character(3) :: collector, dxtheory, partmod, parton
     integer, allocatable, dimension(:, :, :) :: corners, netcatch
     logical :: debug, hitstickbottom, hitstickside, line, linear, particle,   &
              & showsky, wbskip
     real, allocatable, dimension(:) :: delaydate, fxloci, fzloci, sedvel
     real(R8KIND), dimension(5) :: fishdo, fishtemp, fxvel, fzvel
     real, allocatable, dimension(:, :) :: fishes
     real(R8KIND) :: fx00count, fxdocount, fxhvcount, fxrdcount, fxtpcount,    &
                   & fxvvcount, fz00count, fzdocount, fzhvcount, fzrdcount,    &
                   & fztpcount, fzvvcount, nfsfreq, rundiff, temporarygrad,    &
                   & totfxwgt, totfzwgt, value, xgraddo, xgradhv, xgradtp,     &
                   & xgradvv, zgraddo, zgradhv, zgradtp, zgradvv
     real(R8KIND), dimension(4) :: graddo, gradtemp, gradxvel, gradzvel
     integer, allocatable, dimension(:) :: group, icoll, icollb, icollt, ifish,&
          & ifishb, ifisht, nbrf, snagcount, sndcount
     real :: lftvvel, ljday, lrunday, maxreactxdo, maxreactxhv, maxreactxtp,   &
           & maxreactxvv, maxreactzdo, maxreactzhv, maxreactztp, maxreactzvv,  &
           & maxufish, maxwfish, midid, midido, midih, midihvel, midit,        &
           & miditemp, midiv, midivvel, midkd, midkdo, midkh, midkhvel, midkt, &
           & midktemp, midkv, midkvvel, miltime, minreactxdo, minreactxhv,     &
           & minreactxtp, minreactxvv, minreactzdo, minreactzhv, minreactztp,  &
           & minreactzvv, mxxspdl, mxzspdl, oldfxloc, oldfyloc, oldfzloc,      &
           & otherfishx, otherfishz, outfreq, outfreqapr, outfreqaug,          &
           & outfreqdec, outfreqfeb, outfreqjan, outfreqjul, outfreqjun,       &
           & outfreqmar, outfreqmay, outfreqnov, outfreqoct, outfreqp,         &
           & outfreqsep, rdx, rdxweigt, rdy, rdyweigt, rdz, rdzweigt, rgtdo,   &
           & rgthvel, rgttemp, rgtvvel, rrr, skydawn, skyday, skydusk,         &
           & skynight, tempopt, tempoptd, tempoptn, tempthres, topdo, tophvel
     integer, allocatable, dimension(:, :) :: limpbr, ndinfo, rimpbr
     real, dimension(4) :: rmat
     real :: toptemp, topvvel, totfxspan, totfzspan, tpxweigt, tpzweigt,       &
           & tsetp2mult, tstep1, tstep1mult, tstep2, tstep2mult, ubdysearch,   &
           & ufish, uzsensph, vfish, vvelcap, vvxweigt, vvzweigt, wfish, wmax, &
           & xdispcoeff, xdotrucsum, xdotrunc, xhvtrucsum, xhvtrunc, xinzone,  &
           & xlok, xmdzone, xotzone, xrefl, xschzoneinn, xschzonemid,          &
           & xschzoneout, xsepcomfort, xtptrucsum, xtptrunc, xurgency,         &
           & xurgtrak, xvvtrucsum, xvvtrunc, yrefl, yurgency, zbotrefl,        &
           & zdispcoeff, zdotrucsum, zdotrunc, zhvtrucsum, zhvtrunc, zinzone,  &
           & zlok, zmax, zmdzone, zmin, zotzone, zschzoneinn, zschzonemid,     &
           & zschzoneout, zsepcomfort, zsurrefl, ztptrucsum, ztptrunc,         &
           & zurgency, zurgtrak, zvvtrucsum, zvvtrunc
     character(4) :: txtdorule, txtmultirule, txtnullff, txtnullwq,            &
                   & txtparticle, txtpasstran, txtrandrule, txtschoolg,        &
                   & txtstimrule, txttemprule, txtvelorule
!
!*** End of declarations rewritten by SPAG
!
                                     !   the number of Water Quality Parameters evaluated at nodes for use in the NFS
                                                    ! The number of fish parameters. The number of Flow Field Parameters evaluated &
                                ! SW 2/16/01
                                                                             ! SW 7/1/2017
 
                                                                                               !,DELAYDATE
 
 
                                                       ! SW 1/16/01
 
                                                                ! SW 2/16/01 7/25/2017
                                                                         ! SW 1/16/01
 
 
     end module FISHY                                              ! SW 2/16/01
