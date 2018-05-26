!*==main.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module MAIN
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
     character(72), parameter :: CONFN = 'w2_con.npt'
!
! Local variables
!
     real(R8KIND), allocatable, dimension(:) :: a, albedo, ax, betai, bhrho, bta,  &
           & c, cbhe, d, dlvol, dlvr, dlxrho, dxi, elbot, esr, etr, etsr, ev,  &
           & evbr, f, fetch, fi, fsed, fsedc1, fsedc2, fsod, gammai, gma, hwi, &
           & icemin, icesw, icet2, icethi, q, qdt, qinsum, qpr, qprbr, qssum,  &
           & qwdo, ranlw, rn, rs, sedci, sedci1, sedci2, srosh, t2i, tinsum,   &
           & tpb, tsed, v, volb, volg, wshx, xbr
     real(R8KIND) :: ab, bhrsum, bhsum, colb, coldep, del, depkti, dfc, dlmr,      &
               & dltcal, dlxmin, dtv, ea, effric, elt, eltms, es, gc2, heatex, &
               & hia, hice, hmin, hrad, iceth1, iceth2, icethu, qinfr, rhoicp, &
               & rhoin, rhoirl1, sroin, sronet, sroout, srosed, sstot, t2r4,   &
               & tairv, tau1, tau2, tflux, tice, udl, udr, v1, vqin, vqini,    &
               & wt1, wt2, wwt, zb
     logical :: add_layer, branch_found, derived_calc, downstream_outflow,     &
              & dsi_calc, end_run, error_open, fish_particle_exist,            &
              & ice_computation, new_page, n_calc, phbuff_exist, po4_calc,     &
              & pumps, restart_in, restart_out, retlog, sedcomp_exist,         &
              & spillway, sub_layer, surface_warning, tdgon, time_series,      &
              & update_kinetics, update_rates, volume_warning, warning_open,   &
              & water_age_active, weir_calc
     character(8) :: aeratec, ccc, closec, dltinter, envirpc, ext, habtatc,    &
                   & inituwl, limc, rsic, rsoc, selectc, tsrc, wdoc
     integer, save :: aeratefn, con, fishhabfn, flowbfn, massbfn, wlfn
     character(10), allocatable, dimension(:) :: alch, anch, apch, c2ch, cdch, &
                & cdwdoc, cwdoc, epch, kfch, macch
     logical, allocatable, dimension(:) :: alg_calc, allow_ice, bod_calc,      &
          & bod_calcn, bod_calcp, contour, detailed_ice, dn_head,              &
          & energy_balance, evaporation, fetch_calc, flux, head_boundary,      &
          & implicit_visc, iso_sediment, iso_sediment1, iso_sediment2,         &
          & iso_temp, long_profile, long_sediment, long_sediment1,             &
          & long_sediment2, long_temp, mass_balance, place_qin, place_qtr,     &
          & profile, pumpon, screen_output, sediment_calc, sediment_calc1,     &
          & sediment_calc2, snapshot, specify_qtr, spreadsheet, vector,        &
          & vert_profile, vert_sediment, vert_sediment1, vert_sediment2,       &
          & vert_temp, volume_balance, zero_slope
     character(8), allocatable, dimension(:) :: azslc, cdac, cdpltc, celc,     &
               & cplc, cpltc, dltadd, dtcac, dtrc, dtric, dynsedk, ebc, evc,   &
               & exc, exic, fetchc, flxc, fricc, gasgtc, gasspc, gridc, hdic,  &
               & heatc, hpltc, icec, incac, kfac, mbc, metic, pdgtc, pdpic,    &
               & pdspc, ppuc, pqc, prc, prcac, prfc, pugtc, pupic, puspc, qinc,&
               & qinic, qoutc, rhevc, scrc, sedcc, sedcc1, sedcc2, sedprc,     &
               & sedprc1, sedprc2, sedrc, seg, slhtc, slicec, sltrc, snpc,     &
               & sprc, sroc, tecplot, trc, trcac, tric, vbc, visc, vplc, wdic, &
               & windc, wtypec
     integer, allocatable, dimension(:, :) :: bl, cdn, iprf, isnp, ispr, kbswt,&
          & ktswt, wdo, wdo2
     character(10), save :: blank, blank1
     integer, allocatable, dimension(:) :: bth, ibpr, icpl, ilat, ilayer, itr, &
          & itsr, iwd, iwdo, jbdam, jbtr, jbwd, jss, kbdgt, kbdpi, kbdsp, kbi, &
          & kbmax, kbqin, kbr, kbtr, kbugt, kbupi, kbusp, ktdgt, ktdpi, ktdsp, &
          & ktqin, kttr, ktugt, ktupi, ktusp, kwd, lpr, nbl, ncpl, nflx, niprf,&
          & nisnp, nispr, nl, npoint, nprf, nscr, nsnp, nspr, nvpl, tsr, vpr
     character(72), allocatable, dimension(:) :: bthfn, cplfn, flxfn, flxfn2,  &
                & lprfn, prffn, snpfn, sprfn, vplfn, vprfn
     real(R8KIND), allocatable, dimension(:, :) :: c2i, epici, hseg, p, qtot, qtrf,&
           & tssdh1, tssuh1
     real, allocatable, dimension(:) :: cdsum, cdtot, csum, dltd, dltf, dltmax,&
                                      & ebdgt, ebdpi, ebdsp, ebpu, ebugt,      &
                                      & ebupi, ebusp, eltrb, eltrt, etdgt,     &
                                      & etdpi, etdsp, etpu, etugt, etupi,      &
                                      & etusp, qimxr, qoavr, qomxr, qtavb,     &
                                      & qtmxb, rsod, rsof, sedcic, sedcin,     &
                                      & sedcip, sedcis, tsedf, tsrd, tsrf,     &
                                      & wdod, wdof, x1
     character(8), allocatable, dimension(:, :) :: cdtbrc, cdwbc, cinbrc,      &
               & cprbrc, cprwbc, ctrtrc, epic, epiprc, hprwbc, kfwbc, sinkc,   &
               & sinkct, stric
     real, allocatable, dimension(:, :) :: cdwdo, cinsum, cmbrs, cout, cpb,    &
       & cpld, cplf, cwdo, estrt, flxd, flxf, hab, kfjw, prfd, prff, qinf,     &
       & scrd, scrf, sedvp, sedvp1, sedvp2, snpd, snpf, sprd, sprf, tvp, vpld, &
       & vplf, wstrt
     real :: celrty, dlxmax, hmax, jday1, jdayts, nxtvd, tmend, tmstrt, ttime
     character(10), allocatable, dimension(:, :) :: conv1
     real(R8KIND), save :: cp, ice_tol, rhoa, rhoi, rimt, rk1, rl1, thrkti, vtol
     real(R8KIND), allocatable, dimension(:, :, :) :: cssdh1, cssuh1
     character(4), allocatable, dimension(:) :: cunit1
     real, allocatable, dimension(:, :, :) :: cvp, epivp, macrclp, macrcvp
     character(2) :: deg
     logical, allocatable, dimension(:, :) :: epiphyton_calc, iso_conc,        &
          & iso_epiphyton, iso_macrophyte, long_conc, long_epiphyton,          &
          & long_macrophyte, tdg_gate, tdg_spillway, vert_conc, vert_epiphyton,&
          & vert_macrophyte
     character(1) :: esc
     character(3) :: gdch
     integer :: grf, ibod, idt, ie, ii, incr, incris, is, iut, j, ja, jac, jcb,&
              & jd, jdaynx, je, jf, jfile, jg, jg_age, jh, jjjb, jjz, jm, jp,  &
              & js, jsg, jt, jwd, jwr, k, kbp, kf_do_sed, kf_do_sod,           &
              & kf_nh4_sed, kf_nh4_sod, kf_po4_sed, kf_po4_sod, kf_sed_nburial,&
              & kf_sed_pburial, ktip, l, l1, l2, l3, m, nae, nalk, nas, nbodce,&
              & nbodcs, nbode, nbodne, nbodns, nbodpe, nbodps, nbods, nccs,    &
              & ndlt, ndo, ndsi, ndsp, ndt, nfe, ngc, ngctdg, ngn2, nit1, niw, &
              & niwdo, nldom, nlpom, nnh4, nno3, nnsg, npo4, npsi, nrdom,      &
              & nrpom, nrs
     real(8), allocatable, dimension(:) :: iceqss
     character(45), allocatable, dimension(:) :: kfname
     character(14), allocatable, dimension(:) :: kfname2
     character(72) :: line, rsofn, segnum, segnum2, text, tsrfn, tsrfn1, wdofn
     integer :: ndg = 16
     integer :: nrso, nsse, nsss, ntac, ntacmn, ntacmx, ntds, ntic, ntrt, ntsr,&
              & nwdt, rsi, sif, vsf
     character(10) :: sedcch, sedch, sednch, sedpch
!
!*** End of declarations rewritten by SPAG
!
!    Variable declaration
  !INTEGER       :: J,NIW,NGC,NGCS,NTDS,NCCS,NGCE,NSSS,NSSE,NPO4,NNH4
                                                                         ! CEMA -placed NGCS and NGCE in global module  SW 10/16/2015
                                                                                                                        ! SW 3/9/16
                                                                                   ! cb 6/6/10
                                                     ! SW 4/19/10
                                                                       ! SW 4/19/10
                                        ! SR 7/27/2017
                                                        ! 1/16/13
                                                          ! cb 10/12/11
                                                         ! cb 1/16/13
                                        ! SW 4/30/15
                                                                                                                                         ! SW 7/31/09; 8/24/09
 
!    Allocatable array declarations
 
  !REAL(R8),      ALLOCATABLE, DIMENSION(:)     :: SEDCI1,SEDCI2   ! Amaila
                                                                         ! SW 5/26/15  SR 7/27/2017
                                                                              ! cb 6/7/17, Amaila
                                                           ! QINT,   QOUTT,
                                                                     ! Amaila
                                                                                  ! cb 8/21/15
                                                                                           ! SW 1/23/06
                                                                                            ! cb 1/26/09
                                                                                                  ! ICE_IN,     RC/SW 4/28/11
                                                                                  ! Amaila
                                                                                                                                                     ! Amaila
                                                                                                                       ! cb 5/19/2011
                                                                                                      ! cb 8/21/15
                                                                                       ! Amaila
                                                                                                                       !SW 07/16/04
                                                                                                  !, KFNAME2
 
!    Data declarations
     data rk1/2.12D0/, rl1/333507.0D0/, rimt/0.0D0/, rhoa/1.25D0/,             &
        & rhoi/916.0D0/, vtol/1.0D3/, cp/4186.0D0/                                                       ! SW 10/20/15
     data ice_tol/0.0050D0/
                           !0.005
     data blank/'          '/, blank1/'    -99.00'/
     data con/10/                      !,  RSI /11/
     data thrkti/0.10D0/
     data flowbfn/9500/, wlfn/9510/, aeratefn/9520/, fishhabfn/9530/,          &
        & massbfn/9501/                                                                      ! SW 5/25/15 NOTE THAT FISHHABFN INCREMENTS 3 TIMES SO 9531,9532,9533 ARE RESERVED ALSO JWFILE1=9549+NWB  JBFILE1=9749+NBR
 
     end module MAIN
