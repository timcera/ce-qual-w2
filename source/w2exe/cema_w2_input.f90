!*==cema_w2_input.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     subroutine CEMA_W2_INPUT
 
    ! Type declarations
     use MAIN
     use GLOBAL
     use GEOMC
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: dayzz, monzz, ninp, yearzz
     logical :: file_exists, skiploop
     character(256) :: messagetemp
!
!*** End of declarations rewritten by SPAG
!
 
 
    ! start CEMA testing variables
 
     sd_global = .FALSE.
     includeiron = .FALSE.
     includemanganese = .FALSE.
     includedynamicph = .FALSE.
     includealkalinity = .FALSE.
     ngh2s = 0
     ngch4 = 0
     ngso4 = 0
     ngturb = 0
     ngfe2 = 0
     ngfeooh = 0
     ngmn2 = 0
     ngmno2 = 0
     ngmft = 0
 
    ! end CEMA testing variables
 
     inquire(file = "W2_CEMA_Input.npt", exist = file_exists)  ! file_exists will be TRUE if the file
 
     if(.NOT.file_exists)then
         cemarelatedcode = .FALSE.
         includebedconsolidation = .FALSE.
         return
     endif
     cemarelatedcode = .TRUE.
!
     cemafiln = nunit
!
     nunit = nunit + 1               ! SW 5/26/15
     open(cemafiln, file = "W2_CEMA_Input.npt")
     open(cemalogfiln, file = "CEMALogFile.opt")
!
!!Read Header
     skiploop = .FALSE.
     do while (.NOT.skiploop)
         read(cemafiln, '(a)')messagetemp
         if(INDEX(messagetemp, "$")==0)skiploop = .TRUE.
     enddo
 
     backspace(cemafiln)
     read(cemafiln, *)messagetemp, sd_global
     if(.NOT.sd_global)then
         cemarelatedcode = .FALSE.
         includebedconsolidation = .FALSE.
         return
     endif
     read(cemafiln, *)messagetemp, nh2s       ! cb 2/18/13  reading constituent #'s for h2s, ch4, so4, co2, turbidity, and mft
     read(cemafiln, *)messagetemp, nch4
     read(cemafiln, *)messagetemp, nso4
     read(cemafiln, *)messagetemp, nturb
     read(cemafiln, *)messagetemp, nfe2
     read(cemafiln, *)messagetemp, nfeooh
     read(cemafiln, *)messagetemp, nmn2
     read(cemafiln, *)messagetemp, nmno2
     read(cemafiln, *)messagetemp, nmft
     ngh2s = nh2s - ngcs + 1
     ngch4 = nch4 - ngcs + 1
     ngso4 = nso4 - ngcs + 1
     ngturb = nturb - ngcs + 1
     ngfe2 = nfe2 - ngcs + 1
     ngfeooh = nfeooh - ngcs + 1
     ngmn2 = nmn2 - ngcs + 1
     ngmno2 = nmno2 - ngcs + 1
     ngmft = nmft - ngcs + 1
 
     read(cemafiln, *)messagetemp, includebedconsolidation
     read(cemafiln, *)messagetemp, layeraddthkfrac
     read(cemafiln, *)messagetemp, numconsolidregns
!
     allocate(consolidationtype(numconsolidregns),                             &
            & constconsolidrate(numconsolidregns))
     allocate(constporewtrrate(numconsolidregns),                              &
            & consolidrateregnfil(numconsolidregns))
     allocate(consregsegst(numconsolidregns), consregsegen(numconsolidregns))
!
     read(cemafiln, *)messagetemp, (consregsegst(i), i = 1, numconsolidregns)
     read(cemafiln, *)messagetemp, (consregsegen(i), i = 1, numconsolidregns)
     read(cemafiln, *)messagetemp,                                             &
                    & (consolidationtype(i), i = 1, numconsolidregns)
     read(cemafiln, *)messagetemp,                                             &
                    & (constconsolidrate(i), i = 1, numconsolidregns)
     read(cemafiln, *)messagetemp,                                             &
                    & (consolidrateregnfil(i), i = 1, numconsolidregns)
!
     read(cemafiln, *)messagetemp, bedelevationinit
     read(cemafiln, *)messagetemp, bedporosityinit
     read(cemafiln, *)messagetemp, cemapwph
     read(cemafiln, *)messagetemp, includedynamicph
     read(cemafiln, *)messagetemp, includealkalinity
     read(cemafiln, *)messagetemp, cemaparticlesize
     cemaparticlesize = 1.D-6*cemaparticlesize      !Microns to m
     read(cemafiln, *)messagetemp, cemasedimenttype
     read(cemafiln, *)messagetemp, cemasedimentdensity
     read(cemafiln, *)messagetemp, cemasedimentsvelocity
     cemasedimentsvelocity = cemasedimentsvelocity/86400.D0     !m/d to m/s
     read(cemafiln, *)messagetemp, cemasedimentprocessesinc
!
     allocate(bedelevation(imx), bedelevationlayer(imx), bedporosity(imx))
    !allocate (sedcellwidth(imx))
     allocate(consolidregnnum(imx), bedconsolidrate(imx), porewaterrelrate(imx)&
            & )
     allocate(cemasedconc(imx, kmx))
     allocate(cemacumpwrelease(imx), cemalayeradded(imx), cemassapplied(imx))
     allocate(cemacumpwtorelease(imx), cemacumpwreleased(imx))
     allocate(numcemapwinst(imx))
     allocate(applycemapwrelease(imx))
     allocate(cemacumpwreleaserate(imx))
     allocate(endbedconsolidation(imx))
     allocate(cematsscopy(kmx, imx))
     allocate(volcema(nbr))
 
     bedelevation = bedelevationinit
     bedelevationlayer = 0.D00
    !sedcellwidth=0.d00
     bedporosity = bedporosityinit
     bedconsolidrate = 0.D00
     porewaterrelrate = 0.D00
     cemasedconc = 0.D00
     cemacumpwrelease = 0.D00
     cemacumpwreleaserate = 0.D00
     cemacumpwtorelease = 0.D00
     cemacumpwreleased = 0.D00
     endbedconsolidation = .FALSE.
     volcema = 0.D00
     numcemapwinst = 0
     applycemapwrelease = .FALSE.
 
     read(cemafiln, *)messagetemp, writebesnp
     read(cemafiln, *)messagetemp, writepwsnp
 
     read(cemafiln, *)messagetemp, includefftlayer
     read(cemafiln, *)messagetemp, numfftactiveprds
     if(includefftlayer)firsttimeinfftcode = .TRUE.
     allocate(fftactprdst(numfftactiveprds), fftactprden(numfftactiveprds))
     allocate(fftlayconc(imx))
     read(cemafiln, *)messagetemp, (fftactprdst(i), i = 1, numfftactiveprds)
     read(cemafiln, *)messagetemp, (fftactprden(i), i = 1, numfftactiveprds)
     read(cemafiln, *)messagetemp, initfftlayerconc
     read(cemafiln, *)messagetemp, fftlayersettvel
     fftlayersettvel = fftlayersettvel/86400.D00
                                                !m/d --> m/s
     fftlayconc = 0.D00
     fftactprd = 1
     movefftlayerdown = .FALSE.
     read(cemafiln, *)messagetemp, movefftlayerdown
 
     read(cemafiln, *)messagetemp, includecemaseddiagenesis
     if(includecemaseddiagenesis)sediment_diagenesis = .TRUE.
     if(.NOT.includecemaseddiagenesis)return
     if(includecemaseddiagenesis)firsttimeincemamftseddiag = .TRUE.
     includecemagenbodconstituents = .FALSE.
     read(cemafiln, *)messagetemp, numregnssedimentbedcomposition
 
     allocate(sdregnpoc_t(numregnssedimentbedcomposition),                     &
            & sdregnpon_t(numregnssedimentbedcomposition),                     &
            & sdregnsul_t(numregnssedimentbedcomposition))
     allocate(sdregnpop_t(numregnssedimentbedcomposition))
     allocate(sdregnh2s_t(numregnssedimentbedcomposition),                     &
            & sdregnnh3_t(numregnssedimentbedcomposition),                     &
            & sdregnch4_t(numregnssedimentbedcomposition))
     allocate(sdregntic_t(numregnssedimentbedcomposition),                     &
            & sdregnalk_t(numregnssedimentbedcomposition),                     &
            & sdregnpo4_t(numregnssedimentbedcomposition))
     allocate(sdregnfe2_t(numregnssedimentbedcomposition),                     &
            & sdregnfeooh_t(numregnssedimentbedcomposition))
     allocate(sdregnmn2_t(numregnssedimentbedcomposition),                     &
            & sdregnmno2_t(numregnssedimentbedcomposition))
     allocate(sdregnt_t(numregnssedimentbedcomposition))
     allocate(sedbedinitregsegst(numregnssedimentbedcomposition),              &
            & sedbedinitregsegen(numregnssedimentbedcomposition))
!
     read(cemafiln, *)messagetemp, (sedbedinitregsegst(i), i = 1,              &
                    & numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp, (sedbedinitregsegen(i), i = 1,              &
                    & numregnssedimentbedcomposition)
 
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnt_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnpoc_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnpon_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnpop_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnsul_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnnh3_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnpo4_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnh2s_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnch4_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregntic_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnalk_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnfe2_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnfeooh_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnmn2_t(i), i = 1, numregnssedimentbedcomposition)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnmno2_t(i), i = 1, numregnssedimentbedcomposition)
 
     read(cemafiln, *)messagetemp, includecemagenbodconstituents
    !If(IncludeCEMAGenBODConstituents)Then
     read(cemafiln, *)messagetemp, numgenbodconstituents
     allocate(sedgenbodname(numgenbodconstituents))
     read(cemafiln, *)messagetemp,                                             &
                    & (sedgenbodname(i), i = 1, numgenbodconstituents)
 
     read(cemafiln, *)messagetemp, numgenbodinitregns
     allocate(sedgenbodinit(numgenbodconstituents, numgenbodinitregns))
     allocate(sedgenbodregsegst(numgenbodinitregns),                           &
            & sedgenbodregsegen(numgenbodinitregns))
     read(cemafiln, *)messagetemp,                                             &
                    & (sedgenbodregsegst(i), i = 1, numgenbodinitregns)
     read(cemafiln, *)messagetemp,                                             &
                    & (sedgenbodregsegen(i), i = 1, numgenbodinitregns)
 
     do i = 1, numgenbodconstituents
         read(cemafiln, *)messagetemp,                                         &
                        & (sedgenbodinit(i, j), j = 1, numgenbodinitregns)
     enddo
 
     read(cemafiln, *)messagetemp, numgenbodconsumptionregions
     allocate(sedgenbodconsregsegst(numgenbodconsumptionregions),              &
            & sedgenbodconsregsegen(numgenbodconsumptionregions))
     allocate(sedgenbodregnrate(numgenbodconstituents,                         &
            & numgenbodconsumptionregions),                                    &
            & sedgenbodregntcoeff(numgenbodconstituents,                       &
            & numgenbodconsumptionregions))
     read(cemafiln, *)messagetemp, (sedgenbodconsregsegst(i), i = 1,           &
                    & numgenbodconsumptionregions)
     read(cemafiln, *)messagetemp, (sedgenbodconsregsegen(i), i = 1,           &
                    & numgenbodconsumptionregions)
 
     do i = 1, numgenbodconstituents
         read(cemafiln, *)messagetemp, (sedgenbodregnrate(i, j), j = 1,        &
                        & numgenbodconsumptionregions)
     enddo
     do i = 1, numgenbodconstituents
         read(cemafiln, *)messagetemp, (sedgenbodregntcoeff(i, j), j = 1,      &
                        & numgenbodconsumptionregions)
     enddo
    !End If
!
     read(cemafiln, *)messagetemp, numregnssedimentdiagenesis
 
     allocate(sdregnpoc_l_fr(numregnssedimentdiagenesis),                      &
            & sdregnpoc_r_fr(numregnssedimentdiagenesis),                      &
            & sdregnpon_l_fr(numregnssedimentdiagenesis))
     allocate(sdregnpon_r_fr(numregnssedimentdiagenesis),                      &
            & sdregnpw_diffcoeff(numregnssedimentdiagenesis),                  &
            & sdregnox_threshold(numregnssedimentdiagenesis))
     allocate(sdregnpop_l_fr(numregnssedimentdiagenesis),                      &
            & sdregnpop_r_fr(numregnssedimentdiagenesis))
     allocate(sdregnae_nh3_no3_l(numregnssedimentdiagenesis),                  &
            & sdregnae_nh3_no3_h(numregnssedimentdiagenesis),                  &
            & sdregnae_no3_n2_l(numregnssedimentdiagenesis))
     allocate(sdregnae_no3_n2_h(numregnssedimentdiagenesis),                   &
            & sdregnan_no3_n2(numregnssedimentdiagenesis),                     &
            & sdregnae_ch4_co2(numregnssedimentdiagenesis))
     allocate(sdregnae_hs_nh4_nit(numregnssedimentdiagenesis),                 &
            & sdregnae_hs_o2_nit(numregnssedimentdiagenesis),                  &
            & sdregn_theta_pw(numregnssedimentdiagenesis))
     allocate(sdregn_theta_nh3_no3(numregnssedimentdiagenesis),                &
            & sdregn_theta_no3_n2(numregnssedimentdiagenesis),                 &
            & sdregn_theta_ch4_co2(numregnssedimentdiagenesis))
     allocate(sdregn_sulfate_ch4_h2s(numregnssedimentdiagenesis),              &
            & sdregnae_h2s_so4(numregnssedimentdiagenesis),                    &
            & sdregn_theta_h2s_so4(numregnssedimentdiagenesis))
     allocate(sdregn_normconst_h2s_so4(numregnssedimentdiagenesis),            &
            & sdregn_minrate_pon_lab(numregnssedimentdiagenesis),              &
            & sdregn_minrate_pon_ref(numregnssedimentdiagenesis))
     allocate(sdregn_minrate_pon_ine(numregnssedimentdiagenesis),              &
            & sdregn_minrate_poc_lab(numregnssedimentdiagenesis),              &
            & sdregn_minrate_poc_ref(numregnssedimentdiagenesis))
     allocate(sdregn_minrate_poc_ine(numregnssedimentdiagenesis),              &
            & sdregn_theta_pon_lab(numregnssedimentdiagenesis),                &
            & sdregn_theta_pon_ref(numregnssedimentdiagenesis))
     allocate(sdregn_theta_pon_ine(numregnssedimentdiagenesis),                &
            & sdregn_theta_poc_lab(numregnssedimentdiagenesis),                &
            & sdregn_theta_poc_ref(numregnssedimentdiagenesis))
     allocate(sdregn_theta_poc_ine(numregnssedimentdiagenesis),                &
            & sdregn_ch4compmethod(numregnssedimentdiagenesis),                &
            & sdregn_pomresuspmethod(numregnssedimentdiagenesis))
     allocate(sdregn_theta_pop_lab(numregnssedimentdiagenesis),                &
            & sdregn_theta_pop_ref(numregnssedimentdiagenesis),                &
            & sdregn_theta_pop_ine(numregnssedimentdiagenesis))
     allocate(sdregn_minrate_pop_lab(numregnssedimentdiagenesis),              &
            & sdregn_minrate_pop_ref(numregnssedimentdiagenesis),              &
            & sdregn_minrate_pop_ine(numregnssedimentdiagenesis))
     allocate(sedbeddiarcregsegst(numregnssedimentdiagenesis),                 &
            & sedbeddiarcregsegen(numregnssedimentdiagenesis))
!
     read(cemafiln, *)messagetemp, (sedbeddiarcregsegst(i), i = 1,             &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sedbeddiarcregsegen(i), i = 1,             &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnpoc_l_fr(i), i = 1, numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnpoc_r_fr(i), i = 1, numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnpon_l_fr(i), i = 1, numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnpon_r_fr(i), i = 1, numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnpop_l_fr(i), i = 1, numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnpop_r_fr(i), i = 1, numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnpw_diffcoeff(i), i = 1, numregnssedimentdiagenesis&
                    & )
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnox_threshold(i), i = 1, numregnssedimentdiagenesis&
                    & )
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnae_nh3_no3_l(i), i = 1, numregnssedimentdiagenesis&
                    & )
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnae_nh3_no3_h(i), i = 1, numregnssedimentdiagenesis&
                    & )
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnae_no3_n2_l(i), i = 1, numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnae_no3_n2_h(i), i = 1, numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnan_no3_n2(i), i = 1, numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnae_ch4_co2(i), i = 1, numregnssedimentdiagenesis)              !Eq. 10.35
     read(cemafiln, *)messagetemp, (sdregnae_hs_nh4_nit(i), i = 1,             &
                    & numregnssedimentdiagenesis)                                              !Eq. 3.3
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnae_hs_o2_nit(i), i = 1, numregnssedimentdiagenesis&
                    & )                                                                       !Eq. 3.3
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregn_theta_pw(i), i = 1, numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_theta_nh3_no3(i), i = 1,            &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_theta_no3_n2(i), i = 1,             &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_theta_ch4_co2(i), i = 1,            &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_sulfate_ch4_h2s(i), i = 1,          &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp,                                             &
                    & (sdregnae_h2s_so4(i), i = 1, numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_theta_h2s_so4(i), i = 1,            &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_normconst_h2s_so4(i), i = 1,        &
                    & numregnssedimentdiagenesis)                                                   !Eq. 9.6
     read(cemafiln, *)messagetemp, (sdregn_minrate_pon_lab(i), i = 1,          &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_minrate_pon_ref(i), i = 1,          &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_minrate_pon_ine(i), i = 1,          &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_minrate_poc_lab(i), i = 1,          &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_minrate_poc_ref(i), i = 1,          &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_minrate_poc_ine(i), i = 1,          &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_minrate_pop_lab(i), i = 1,          &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_minrate_pop_ref(i), i = 1,          &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_minrate_pop_ine(i), i = 1,          &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_theta_pon_lab(i), i = 1,            &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_theta_pon_ref(i), i = 1,            &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_theta_pon_ine(i), i = 1,            &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_theta_poc_lab(i), i = 1,            &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_theta_poc_ref(i), i = 1,            &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_theta_poc_ine(i), i = 1,            &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_theta_pop_lab(i), i = 1,            &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_theta_pop_ref(i), i = 1,            &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, (sdregn_theta_pop_ine(i), i = 1,            &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, o2h2s        ! cb 5/22/15   moved from below
     read(cemafiln, *)messagetemp, o2ch4        ! cb 5/22/15
     read(cemafiln, *)messagetemp, kdp1      ! cb 5/22/15
     read(cemafiln, *)messagetemp, kdp2      ! cb 5/22/15
     read(cemafiln, *)messagetemp, (sdregn_ch4compmethod(i), i = 1,            &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, nh4_nh3_eqb_const
     read(cemafiln, *)messagetemp, hs_h2s_eqb_const
     read(cemafiln, *)messagetemp, henryconst_nh3
     read(cemafiln, *)messagetemp, henryconst_ch4
     read(cemafiln, *)messagetemp, henryconst_h2s
     read(cemafiln, *)messagetemp, henryconst_co2
     read(cemafiln, *)messagetemp, gasdiff_sed  ! in m^2/s
     read(cemafiln, *)messagetemp, calibparam_r1
     read(cemafiln, *)messagetemp, youngmodulus
     read(cemafiln, *)messagetemp, critstressif
     read(cemafiln, *)messagetemp, bubbrelscale
     read(cemafiln, *)messagetemp, crackclosefraction
     read(cemafiln, *)messagetemp, limbubbsize
     read(cemafiln, *)messagetemp, maxbubbrad
     read(cemafiln, *)messagetemp, usereleasefraction
     read(cemafiln, *)messagetemp, bubbrelfraction
     read(cemafiln, *)messagetemp, bubbaccfraction
     read(cemafiln, *)messagetemp, numbubrelarr
     read(cemafiln, *)messagetemp, bubbrelfractionatm
     read(cemafiln, *)messagetemp, bubbwatgasexchrate
     read(cemafiln, *)messagetemp, applybubbturb
     read(cemafiln, *)messagetemp, cematurbulencescaling
     read(cemafiln, *)messagetemp, coeffa_turb
     read(cemafiln, *)messagetemp, coeffb_turb
     read(cemafiln, *)messagetemp, writecemamftsedflx
     read(cemafiln, *)messagetemp, partmixvel
     read(cemafiln, *)messagetemp, burialvel
     read(cemafiln, *)messagetemp, includeiron      ! cb 5/22/15  moved from above
     read(cemafiln, *)messagetemp, includemanganese ! cb 5/22/15
     read(cemafiln, *)messagetemp, kfeooh_halfsat
     read(cemafiln, *)messagetemp, kfe_red
     read(cemafiln, *)messagetemp, kfe_oxid
     read(cemafiln, *)messagetemp, fesetvel
     read(cemafiln, *)messagetemp, kdfe1
     read(cemafiln, *)messagetemp, kdfe2
     read(cemafiln, *)messagetemp, o2fe2
     read(cemafiln, *)messagetemp, kmno2_halfsat
     read(cemafiln, *)messagetemp, kmn_red
     read(cemafiln, *)messagetemp, kmn_oxid
     read(cemafiln, *)messagetemp, mnsetvel
     read(cemafiln, *)messagetemp, kdmn1
     read(cemafiln, *)messagetemp, kdmn2
     read(cemafiln, *)messagetemp, o2mn2
     read(cemafiln, *)messagetemp, cema_pom_resuspension_processes
     read(cemafiln, *)messagetemp, (sdregn_pomresuspmethod(i), i = 1,          &
                    & numregnssedimentdiagenesis)
     read(cemafiln, *)messagetemp, taucrpom
     read(cemafiln, *)messagetemp, crshields
     read(cemafiln, *)messagetemp, cao_method
     read(cemafiln, *)messagetemp, spgrav_pom
     read(cemafiln, *)messagetemp, dia_pom
     read(cemafiln, *, end = 100)messagetemp, sediagfreq
                                                      ! FREQUENCY OF OUTPUT SW 5/25/2017
     read(cemafiln, *, end = 100)messagetemp, bubbles_calculation
                                                               ! BUBBLES_CALCULATION SW 5/25/2017
     goto 200
    ! IF INPUT FILE DOES NOT HAVE NEW LAST LINE - ERROR TRAPPING FOR NEW CODE
100   sediagfreq = 7.0
     bubbles_calculation = .FALSE.
200   close(cemafiln)
 
    !Allocate other variables
 
     allocate(cemamft_randc_regn(imx), cemamft_incond_regn(imx),               &
            & mftsedflxvars(imx, 47), cema_sd_vars(imx, 21))
     allocate(sd_no3p2(2), sd_nh3p2(2), sd_nh3tp2(2), sd_ch4p2(2), sd_po4p2(2),&
            & sd_po4tp2(2))
     allocate(sd_hsp2(2), sd_hstp2(2), sd_poc22(3), sd_pon22(3), sd_pop22(3))
     allocate(sd_poc2(3), sd_pon2(3), sd_pop2(3), sd_nh3tp(2), sd_no3p(2),     &
            & sd_po4tp(2), sd_hstp(2))
     allocate(sd_fpon(3), sd_fpoc(3), sd_kdiapon(3), sd_thtapon(3),            &
            & sd_kdiapoc(3), sd_thtapoc(3))
     allocate(sd_jpoc(3), sd_jpon(3), sd_jpop(3), sd_nh3(2), sd_tic(2),        &
            & sd_alk(2), sd_ph(2), sd_tds(2))
     allocate(sd_epoc(3), sd_epon(3), sd_epop(3))
     allocate(sd_t(2))
     allocate(sd_denit(2), sd_jdenit(2), sd_jo2no3(2), sd_ch4(2), sd_hs(2),    &
            & sd_po4t(2))
     allocate(sd_fe2t(2), sd_fe2(2), sd_feooh(2), sd_mn2t(2), sd_mn2(2),       &
            & sd_mno2(2))
     allocate(sd_kdiapop(3), sd_thtapop(3), sd_nh3t(2), sd_no3(2), sd_hst(2),  &
            & sd_po4(2), sd_fpop(3))
     allocate(sd_so4conc(imx), sd_phvalue(imx))
     allocate(sd_aerlayerthick(imx))
     allocate(h2sdis(imx), h2sgas(imx), ch4dis(imx), ch4gas(imx))
     allocate(nh4dis(imx), nh4gas(imx), co2dis(imx), co2gas(imx))
     allocate(bubbleradiussed(imx), presbubbsed(imx), prescritsed(imx))
     allocate(cgsed(imx), c0sed(imx), ctsed(imx))
     allocate(tconc(numgas, imx), tconcp(numgas, imx), sconc(numgas, imx))
     allocate(dissolvedgassediments(numgas, imx))
     allocate(crackopen(imx), mftbubbreleased(imx), lastdiffvolume(imx))
     allocate(bubblescarried(imx, numbubrelarr),                               &
            & bubblesradius(imx, numbubrelarr))
     allocate(bubbleslnumber(imx, numbubrelarr),                               &
            & bubblesstatus(imx, numbubrelarr))
     allocate(bubblesrisev(imx, numbubrelarr))
     allocate(bubblesgasconc(imx, numbubrelarr, numgas))
     allocate(brvoluagas(imx, numbubrelarr, numgas),                           &
            & brrateagas(imx, numbubrelarr, numgas))
     allocate(firstbubblesrelease(imx, numbubrelarr),                          &
            & bubblesreleaseallvalue(imx, numbubrelarr),                       &
            & bubblerelwb(nwb, numgas))                                                                                      ! SW 7/1/2017
     allocate(brrateagasnet(imx, numgas))
     allocate(bubblesatsurface(imx, numbubrelarr))
     allocate(bottomturbulence(imx))
!!Allocate(IceQSS(IMX)) Put in main code sw 6/17/15
!
!!Sediment Generic BOD
     allocate(sedgenbodconc(numgenbodconstituents, imx, 2),                    &
            & sedgenboddecayrate(numgenbodconstituents, imx))
     allocate(sedgenbodconsrate(numgenbodconstituents, imx),                   &
            & sedgenbodconstcoeff(numgenbodconstituents, imx))
     allocate(sedgenbodss(numgenbodconstituents), sdpflux(nwb), sdnh4flux(nwb),&
            & sdno3flux(nwb))
!
     sd_no3p2 = 0.D00
!
     sd_nh3p2 = 0.D00
!
     sd_nh3tp2 = 0.D00
!
     sd_ch4p2 = 0.D00
!
     bubblerelwb = 0.0                                                                         ! SW 7/1/2017
     sd_po4p2 = 0.D00
     sd_po4tp2 = 0.D00
     sd_hsp2 = 0.D00
     sd_hstp2 = 0.D00
     sd_poc22 = 0.D00
     sd_pon22 = 0.D00
     sd_pop22 = 0.D00
     sd_poc2 = 0.D00
     sd_pon2 = 0.D00
     sd_pop2 = 0.D00
     sd_nh3tp = 0.D00
     sd_no3p = 0.D00
     sd_po4tp = 0.D00
     sd_hstp = 0.D00
     sd_fpon = 0.D00
     sd_fpoc = 0.D00
     sd_kdiapon = 0.D00
     sd_thtapon = 0.D00
     sd_kdiapoc = 0.D00
     sd_thtapoc = 0.D00
     sd_jpoc = 0.D00
     sd_jpon = 0.D00
     sd_jpop = 0.D00
     sd_nh3 = 0.D00
     sd_epoc = 0.D00
     sd_epon = 0.D00
     sd_epop = 0.D00
     sd_denit = 0.D00
     sd_jdenit = 0.D00
     sd_jo2no3 = 0.D00
     sd_ch4 = 0.D00
     sd_hst = 0.D00
     sd_po4 = 0.D00
     sd_fpop = 0.D00
     sd_hs = 0.D00
     sd_po4t = 0.D00
     sd_fe2 = 0.D00
     sd_fe2t = 0.D00
     sd_feooh = 0.D00
     sd_mn2 = 0.D00
     sd_mn2t = 0.D00
     sd_mno2 = 0.D00
     sd_t = 0.D00
     sd_kdiapop = 0.D00
     sd_thtapop = 0.D00
     sd_nh3t = 0.D00
     sd_no3 = 0.D00
     sd_aerlayerthick = 0.D00
     h2sdis = 0.D00
     h2sgas = 0.D00
     ch4dis = 0.D00
     ch4gas = 0.D00
     nh4dis = 0.D00
     nh4gas = 0.D00
     co2dis = 0.D00
     co2gas = 0.D00
     bubbleradiussed = 0.D00
     presbubbsed = 0.D00
     prescritsed = 0.D00
     cgsed = 0.D00
     c0sed = 0.D00
     ctsed = 0.D00
     tconc = 0.D00
     tconcp = 0.D00
     sconc = 0.D00
     crackopen = .FALSE.
     mftbubbreleased = 0
     lastdiffvolume = 0.D00
     bubblescarried = 0
     bubbleslnumber = 0
     bubblesstatus = 0
     bubblesradius = 0.D00
     bubblesrisev = 0.D00
     bubblesgasconc = 0.D00
     cemamft_randc_regn = 0
     mftsedflxvars = 0.D00
     cema_sd_vars = 0.D00
     bubblesreleaseallvalue = 0.D00
     sd_so4conc = 0.D0
     brvoluagas = 0.D00
     brrateagas = 0.D00
     brrateagasnet = 0.D00
     bottomturbulence = 0.D00
    !IceQSS = 0.d00  Put in main code SW 6/17/15
     sedgenbodconc = 0.D00
    !IceQSS = 0.d00  Put in main code SW 6/17/15
     sedgenbodconsrate = 0.D00
     sedgenboddecayrate = 0.D00
     sedgenbodconstcoeff = 0.D00
     sedgenbodss = 0.D00
     dissolvedgassediments = 0.D00
     fesetvel = fesetvel/86400.0
     kfe_red = kfe_red/86400.0
     kfe_oxid = kfe_red/86400.0
     mnsetvel = mnsetvel/86400.0
     kmn_red = kfe_red/86400.0
     kmn_oxid = kfe_red/86400.0
     sdpflux = 0.0
     sdnh4flux = 0.0
     sdno3flux = 0.0
 
     firsttimeinbubbles = .TRUE.
     firstbubblesrelease = .TRUE.
     bubblesatsurface = .FALSE.
 
     bodtestout = .TRUE.
                       ! BOD testing
!    special input for CEMA code testing - start
 
!open(1399,file= 'seddatin.txt',status='old', form='formatted ' )  !Data input
 
!read(1399,*)
!read(1399,*) ninp
 
!do  i = 1,ninp
!read(1399,*)
!read(1399,*) monzz,dayzz,yearzz,sodzz,jnh4zz,jno3zz,jpo4zz,jsizz,o2zz,nh4zz,no3zz,po4zz,sizz,tempzz   !(B) ! NOTE: sodzz, jno3zz, jpo4zz, and jsizz NOT USED
!end do
 
 
!open(1199 ,file='sedin.txt' ,status='old',form='formatted')   !Parameter input
 
!read(1199,*)
!read(1199,*)
!read(1199,*) m1,m2,Dp,w2,Dd ,thtaDp,thtaDd                     !(A)
!read(1199,*)
!read(1199,*)
!read(1199,*) kappnh4 ,pienh4 ,thtanh4 ,kmnh4,thtakmnh4 ,kmnh4o2
!read(1199,*)
!read(1199,*)
!read(1199,*) kapp1no3,k2no3,thtano3
!read(1199,*)
!read(1199,*)
!read(1199,*) kappd1,kappp1,pie1s,pie2s,thtapd1,kmhso2
!read(1199,*)
!read(1199,*)
!read(1199,*) ksi,csisat,dpie1si,pie2si
!read(1199,*)
!read(1199,*)
!read(1199,*) h2ss,thtasi,kmpsi,o2critsi
!read(1199,*)
!read(1199,*)
!read(1199,*) dpie1po4,pie2po4,o2crit,kmo2Dp
!read(1199,*)
!read(1199,*)
!read(1199,*) frpon1,kpon1,thtapon1
!read(1199,*)
!read(1199,*)
!read(1199,*) frpon2,kpon2,thtapon2
!read(1199,*)
!read(1199,*)
!read(1199,*) frpon3,kpon3,thtapon3
!read(1199,*)
!read(1199,*)
!read(1199,*) frpoc1,kpoc1,thtapoc1
!read(1199,*)
!read(1199,*)
!read(1199,*) frpoc2,kpoc2,thtapoc2
!read(1199,*)
!read(1199,*)
!read(1199,*) frpoc3,kpoc3,thtapoc3
!read(1199,*)
!read(1199,*)
!read(1199,*) frpop1,kpop1,thtapop1
!read(1199,*)
!read(1199,*)
!read(1199,*) frpop2,kpop2,thtapop2
!read(1199,*)
!read(1199,*)
!read(1199,*) frpop3,kpop3,thtapop3
!read(1199,*)
!read(1199,*)
!read(1199,*) ratiocn,ratiocp,ratiocsi
 
!    CEMA testing end
 
 
     end subroutine CEMA_W2_INPUT
