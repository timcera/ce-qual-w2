!*==cemavars.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
     module CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     logical :: applybubbturb, bodtestout, bubbles_calculation, cao_method,    &
              & cemarelatedcode, cemasedimentprocessesinc,                     &
              & cema_pom_resuspension_processes, fftactive, firsttimeinbubbles,&
              & firsttimeincemamftseddiag, firsttimeinfftcode,                 &
              & includealkalinity, includebedconsolidation,                    &
              & includecemagenbodconstituents, includecemaseddiagenesis,       &
              & includedynamicph, includefftlayer, includeiron,                &
              & includemanganese, limbubbsize, movefftlayerdown, sd_global,    &
              & sediment_diagenesis, usereleasefraction, writebesnp,           &
              & writecemamftsedflx, writepwsnp
     logical, allocatable, dimension(:) :: applycemapwrelease, cemalayeradded, &
          & cemassapplied, crackopen, endbedconsolidation
     real(8), allocatable, dimension(:) :: bedconsolidrate, bedelevation,      &
          & bedelevationlayer, bedporosity, bottomturbulence, bubbleradiussed, &
          & c0sed, cemacumpwrelease, cemacumpwreleased, cemacumpwreleaserate,  &
          & cemacumpwtorelease, cgsed, ch4dis, ch4gas, co2dis, co2gas,         &
          & constconsolidrate, ctsed, fftlayconc, h2sdis, h2sgas, nh4dis,      &
          & nh4gas, porewaterrelrate, presbubbsed, prescritsed, sdnh4flux,     &
          & sdno3flux, sdpflux, sdregnae_ch4_co2, sdregnae_h2s_so4,            &
          & sdregnae_hs_nh4_nit, sdregnae_hs_o2_nit, sdregnae_nh3_no3_h,       &
          & sdregnae_nh3_no3_l, sdregnae_no3_n2_h, sdregnae_no3_n2_l,          &
          & sdregnalk_t, sdregnan_no3_n2, sdregnch4_t, sdregnfe2_t,            &
          & sdregnfeooh_t, sdregnh2s_t, sdregnmn2_t, sdregnmno2_t, sdregnnh3_t,&
          & sdregnox_threshold, sdregnpo4_t, sdregnpoc_l_fr, sdregnpoc_r_fr,   &
          & sdregnpoc_t, sdregnpon_l_fr, sdregnpon_r_fr, sdregnpon_t,          &
          & sdregnpop_l_fr, sdregnpop_r_fr, sdregnpop_t, sdregnpw_diffcoeff,   &
          & sdregnsul_t, sdregntic_t, sdregnt_t, sdregn_minrate_poc_ine,       &
          & sdregn_minrate_poc_lab, sdregn_minrate_poc_ref,                    &
          & sdregn_minrate_pon_ine, sdregn_minrate_pon_lab,                    &
          & sdregn_minrate_pon_ref, sdregn_minrate_pop_ine,                    &
          & sdregn_minrate_pop_lab, sdregn_minrate_pop_ref,                    &
          & sdregn_normconst_h2s_so4, sdregn_sulfate_ch4_h2s,                  &
          & sdregn_theta_ch4_co2, sdregn_theta_h2s_so4, sdregn_theta_nh3_no3,  &
          & sdregn_theta_no3_n2, sdregn_theta_poc_ine, sdregn_theta_poc_lab,   &
          & sdregn_theta_poc_ref, sdregn_theta_pon_ine
     real(8) :: bedelevationinit, bedporosityinit, bubbaccfraction,            &
              & bubbrelfraction, bubbrelfractionatm, bubbrelscale,             &
              & bubbwatgasexchrate, burialvel, calibparam_r1, cemaparticlesize,&
              & cemapwph, cemasedimentdensity, cemasedimentsvelocity,          &
              & cematurbulencescaling, coeffa_turb, coeffb_turb,               &
              & crackclosefraction, critstressif, crshields, dia_pom, fesetvel,&
              & fftlayersettvel, gasdiff_sed, gasreleasech4, henryconst_ch4,   &
              & henryconst_co2, henryconst_h2s, henryconst_nh3,                &
              & hs_h2s_eqb_const, icethicknesschange, initfftlayerconc, kdfe1, &
              & kdfe2, kdmn1, kdmn2, kdp1, kdp2, kfeooh_halfsat, kfe_oxid,     &
              & kfe_red, kmno2_halfsat, kmn_oxid, kmn_red, layeraddthkfrac,    &
              & maxbubbrad, mnsetvel, nh4_nh3_eqb_const, o2ch4, o2fe2, o2h2s,  &
              & o2mn2, partmixvel, spgrav_pom, taucrpom, totalporewatremoved,  &
              & totalporewatvolume, totalsedimentsinbed,                       &
              & volumeincreasedconsolid, youngmodulus
     real(8), allocatable, dimension(:, :, :) :: brrateagas, brvoluagas,       &
          & bubblesgasconc, sedgenbodconc
     real(8), allocatable, dimension(:, :) :: brrateagasnet, bubblerelwb,      &
          & bubblescarried, bubblesradius, bubblesreleaseallvalue,             &
          & bubblesrisev, cemasedconc, cematsscopy, cema_sd_vars,              &
          & dissolvedgassediments, mftsedflxvars, sconc, sedgenbodconsrate,    &
          & sedgenbodconstcoeff, sedgenboddecayrate, sedgenbodinit,            &
          & sedgenbodregnrate, sedgenbodregntcoeff, tconc, tconcp
     logical, allocatable, dimension(:, :) :: bubblesatsurface,                &
          & firstbubblesrelease
     integer(4), allocatable, dimension(:, :) :: bubbleslnumber, bubblesstatus
     integer(4), save :: cemabtmlayfiln, cemalogfiln, cemaoutfiln1,            &
                       & cemaoutfiln2, cemaoutfiln3, cemaoutfiln4,             &
                       & cemaoutfiln5, cemaoutfiln6, cemasedflxfiln1,          &
                       & cemasedflxfiln2, cemasedflxfiln3, cemasedflxfiln4,    &
                       & cemasedflxfiln5, cemasedflxfiln6, cemasedflxfiln7,    &
                       & cemasedflxfiln8, cemasedflxfiln9, cemasnpoutfiln,     &
                       & cematsr1outfiln, numgas
     integer(4) :: cemafiln, cemasedimenttype, fftactprd, nch4, nfe2, nfeooh,  &
                 & ngch4, ngfe2, ngfeooh, ngh2s, ngmft, ngmn2, ngmno2, ngso4,  &
                 & ngturb, nh2s, nmft, nmn2, nmno2, nso4, nturb, numbubrelarr, &
                 & numconsolidregns, numfftactiveprds,                         &
                 & numregnssedimentbedcomposition, numregnssedimentdiagenesis, &
                 & segnumi, tempcntr1
     integer(4), allocatable, dimension(:) :: cemamft_incond_regn,             &
             & cemamft_randc_regn, consolidationtype, consolidregnnum,         &
             & consregsegen, consregsegst, constporewtrrate, fftactprden,      &
             & fftactprdst, numcemapwinst, sedbeddiarcregsegen,                &
             & sedbeddiarcregsegst, sedbedinitregsegen, sedbedinitregsegst
     real, save :: cemaoutfilbub
     character(256), allocatable, dimension(:) :: consolidrateregnfil
     real(8), save :: gasconst_r
     integer(8), allocatable, dimension(:) :: lastdiffvolume, mftbubbreleased
     integer :: numgenbodconstituents, numgenbodconsumptionregions,            &
              & numgenbodinitregns
     real :: nxtsediag, sediagfreq
     integer, allocatable, dimension(:) :: sdregn_ch4compmethod,               &
          & sdregn_pomresuspmethod, sedgenbodconsregsegen,                     &
          & sedgenbodconsregsegst, sedgenbodregsegen, sedgenbodregsegst
     real(8), allocatable, dimension(:) :: sdregn_theta_pon_lab,               &
          & sdregn_theta_pon_ref, sdregn_theta_pop_ine, sdregn_theta_pop_lab,  &
          & sdregn_theta_pop_ref, sdregn_theta_pw, sd_aerlayerthick, sd_alk,   &
          & sd_ch4, sd_ch4p2, sd_denit, sd_epoc, sd_epon, sd_epop, sd_fe2,     &
          & sd_fe2t, sd_feooh, sd_fpoc, sd_fpon, sd_fpop, sd_hs, sd_hsp2,      &
          & sd_hst, sd_hstp, sd_hstp2, sd_jdenit, sd_jo2no3, sd_jpoc, sd_jpon, &
          & sd_jpop, sd_kdiapoc, sd_kdiapon, sd_kdiapop, sd_mn2, sd_mn2t,      &
          & sd_mno2, sd_nh3, sd_nh3p2, sd_nh3t, sd_nh3tp, sd_nh3tp2, sd_no3,   &
          & sd_no3p, sd_no3p2, sd_ph, sd_phvalue, sd_po4, sd_po4p2, sd_po4t,   &
          & sd_po4tp, sd_po4tp2, sd_poc2, sd_poc22, sd_pon2, sd_pon22, sd_pop2,&
          & sd_pop22, sd_so4conc, sd_t, sd_tds, sd_thtapoc, sd_thtapon,        &
          & sd_thtapop, sd_tic, sedgenbodss, volcema
     character(6), allocatable, dimension(:) :: sedgenbodname
!
!*** End of declarations rewritten by SPAG
!
 
   ! Real(8), Allocatable, Dimension(:) :: sedcellwidth
 
    !Real(8), Allocatable, Dimension(:) :: SDRegn_Theta_POC_Ine, SDRegn_CH4CompMethod
    ! SW 7/1/2017
 
 
 
 
    !Real(8), Allocatable, Dimension(:) :: IceQSS
 
 
    ! cb 10/8/13
    ! cb 2/5/13
    !,GasReleaseCO2   ! SW 10/10/2017   ! SW 10/19/2017
 
    ! SW 5/25/2017
    ! cb 5/22/15
 
    !Generic BOD constituent variables
    ! SW 10/20/2017
    !End generic BOD constituent variables
 
    ! CEMA testing variables start
!    real(8)
!    jnh4zz,odzz,jno3zz,jpo4zz,jsizz,o2zz,nh4zz,no3zz,po4zz,sizz,tempzz  !
!    CEMA testing variable only real(8) m1,m2,Dp,w2,Dd ,thtaDp,thtaDd         
!    !(A) real(8) kappnh4 ,pienh4 ,thtanh4 ,kmnh4,thtakmnh4 ,kmnh4o2 real(8)
!    kapp1no3,k2no3,thtano3 real(8) kappd1,kappp1,pie1s,pie2s,thtapd1,kmhso2
!    real(8) ksi,csisat,dpie1si,pie2si
!    real(8) h2ss,thtasi,kmpsi,o2critsi
!    real(8) dpie1po4,pie2po4,o2crit,kmo2Dp
!    real(8) frpon1,kpon1,thtapon1
!    real(8) frpon2,kpon2,thtapon2
!    real(8) frpon3,kpon3,thtapon3
!    real(8) frpoc1,kpoc1,thtapoc1
!    real(8) frpoc2,kpoc2,thtapoc2
!    real(8) frpoc3,kpoc3,thtapoc3
!    real(8) frpop1,kpop1,thtapop1
!    real(8) frpop2,kpop2,thtapop2
!    real(8) frpop3,kpop3,thtapop3
!    real(8) ratiocn,ratiocp,ratiocsi
!    real(8) xjnh4,jcinzz
!    real(8) SD_Jctest
    ! CEMA testing variables end
 
  !  Data CEMAFilN  /11/    ! SW 5/26/15
     data cemasnpoutfiln/2411/                                                        ! SW 8/31/2017
     data cematsr1outfiln/2412/
     data cemabtmlayfiln/2414/
     data cemasedflxfiln1/2415/
     data cemasedflxfiln2/2416/
     data cemasedflxfiln3/2417/
     data cemalogfiln/2418/
     data cemasedflxfiln4/2419/
     data cemasedflxfiln5/2420/
     data cemasedflxfiln6/2421/
     data cemasedflxfiln7/2422/
     data cemasedflxfiln8/2423/
     data cemasedflxfiln9/2424/
     data cemaoutfiln1/2426/
     data cemaoutfiln2/2429/
     data cemaoutfiln3/2435/
     data cemaoutfiln4/2437/
     data cemaoutfiln5/2438/
     data cemaoutfiln6/2439/
     data cemaoutfilbub/2440/
     data gasconst_r/0.0821/    !L.atm/mol/K
     data numgas/4/
 
     end module CEMAVARS
