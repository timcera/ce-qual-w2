!*==cemasedimentdiagenesis.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine CEMASEDIMENTDIAGENESIS
!!
 
    ! Type declarations
     use MAIN
     use GLOBAL
     use GEOMC
     use SCREENC
     use RSTART
     use EDDY
     use LOGICC
     use TVDC
     use KINETIC
     use SURFHE
     use CEMAVARS
     use TRANS
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(8) :: Alksed, Nh4sed, Phsed, Po4sed, Pocsed, T1sed, Tdssed, Ticsed
     intent (in) Alksed, Nh4sed, Po4sed, Pocsed, T1sed, Tdssed, Ticsed
     intent (inout) Phsed
!
! Local variables
!
     real(8), save :: acn, alkt, ammoniad_sd1, ammoniad_sd2, ammoniag_sd1,     &
                    & ammoniag_sd2, ammt, ao2n, bicart, cart, cellthickness,   &
                    & co2producedcon1l1, co2producedcon2l2, co2producedsrc1l1, &
                    & co2producedsrc1l2, co2producedsrc2l2, co2sed, co3sed,    &
                    & co3t, coef, coef1, coef2, coef3, coef4, c_bottom,        &
                    & c_bottom2, dh1, dh2, dh3, dhh, dissolved_alk_src,        &
                    & dissolved_ch4_src, dissolved_co2_src, dissolved_fe2_src, &
                    & dissolved_h2s_src, dissolved_mn2_src, dissolved_nh3_src, &
                    & dissolved_no3_src, dissolved_o2_snk, dissolved_po4_src,  &
                    & dissolved_so4_src, epsilon, f1, fastfluxcarbon, fetchw,  &
                    & h2co3t, h2po4t, h3po4t, hco3sed, hco3t, hion, hpo4t, hs, &
                    & ht, incr1, jnin, jpin, k1, k2, kamm, kp1, kp2, kp3, kw,  &
                    & lpomn_resuspension, lpomp_resuspension,                  &
                    & lpom_resuspension, lw, lw0, lw1, maxit, molvisc_h2o,     &
                    & mw_constituent, nh3t, nh4t, o2ss, oh, oht, omct, phost
     real(8), allocatable, dimension(:), save :: cellarea
     real(8), dimension(numgenbodconstituents, 2), save :: dissolved_bod_src,  &
          & sd_bodd
     integer, save :: genbodnum, iter1, n, sd_ch4compmethod, sd_pomresuspmethod
     integer(2), save :: initconregn, itemp, iter, regnnum
     integer(I1KIND), save :: nfcle, nflog
     real(8), save :: pht, po4t, pw_relrate, pw_relrate1, pw_relrate2,         &
                    & reyn_resusp, rpomn_resuspension, rpomp_resuspension,     &
                    & rpom_resuspension, s2, sd1_ammonia, sd1_ammonium,        &
                    & sd1_sulfide, sd1_sulfiminus, sd2_ammonia, sd2_ammonium,  &
                    & sd2_sulfide, sd2_sulfiminus, sdpoct1, sd_a11, sd_a12,    &
                    & sd_a21, sd_a22, sd_ae_ch4_co2, sd_ae_h2s_so4,            &
                    & sd_ae_hs_nh4_nit, sd_ae_hs_o2_nit, sd_ae_nh3_no3,        &
                    & sd_ae_nh3_no3_h, sd_ae_nh3_no3_l, sd_ae_no3_n2,          &
                    & sd_ae_no3_n2_h, sd_ae_no3_n2_l, sd_alk0, sd_ammonia,     &
                    & sd_ammonium, sd_an_no3_n2, sd_b1, sd_b2, sd_ben_str,     &
                    & sd_ben_strp, sd_ben_strp2, sd_bod, sd_ch40, sd_ch4sat,   &
                    & sd_ch4toco2, sd_csod, sd_csodmax, sd_depth, sd_e, sd_ea, &
                    & sd_es, sd_f12, sd_f21, sd_fd1, sd_fd2, sd_fe20,          &
                    & sd_fe2tofeooh, sd_feooh0, sd_feoohtofe2, sd_fp1, sd_fp2, &
                    & sd_h1, sd_h2, sd_jalk, sd_jc, sd_jch4, sd_jch4g, sd_jcin,&
                    & sd_jctest, sd_jc_o2equiv, sd_jdenitt, sd_jfe2, sd_jfeooh,&
                    & sd_jfeoohin, sd_jhs, sd_jinfeooh, sd_jmn2, sd_jmno2in,   &
                    & sd_jn
     real(8), save :: sd_jnh4, sd_jnin, sd_jno3, sd_jo2no3t, sd_jp, sd_jpin,   &
                    & sd_jpo4, sd_jso4, sd_jsod, sd_jt, sd_jtic, sd_k1h1d,     &
                    & sd_k1h1p, sd_k2h2d, sd_k2h2p, sd_kappah2sp1, sd_kdh2s1,  &
                    & sd_kdh2s2, sd_kdnh3, sd_kl12, sd_ksw, sd_m1, sd_m2,      &
                    & sd_minrate_poc_ine, sd_minrate_poc_lab,                  &
                    & sd_minrate_poc_ref, sd_minrate_pon_ine,                  &
                    & sd_minrate_pon_lab, sd_minrate_pon_ref,                  &
                    & sd_minrate_pop_ine, sd_minrate_pop_lab,                  &
                    & sd_minrate_pop_ref, sd_mn20, sd_mn2tomno2, sd_mno20,     &
                    & sd_mno2tomn2, sd_nh30, sd_nh3conv, sd_nh3tono3, sd_no30, &
                    & sd_normconst_h2s_so4, sd_nsod, sd_o20, sd_ox_threshol,   &
                    & sd_po40, sd_poct1, sd_poct2, sd_poc_i_fr, sd_poc_l_fr,   &
                    & sd_poc_r_fr, sd_pom, sd_pont2, sd_pon_i_fr, sd_pon_l_fr, &
                    & sd_pon_r_fr, sd_popt2, sd_pop_i_fr, sd_pop_l_fr,         &
                    & sd_pop_r_fr, sd_porosity, sd_pw_diffcoeff, sd_rho,       &
                    & sd_rhowcp, sd_s, sd_sech_arg, sd_sjch4, sd_so4, sd_so40, &
                    & sd_sod, sd_sodold, sd_sulfate_ch4_h2s, sd_sulfide,       &
                    & sd_sulfiminus, sd_taubot, sd_tc, sd_theta_ch4_co2,       &
                    & sd_theta_h2s_so4, sd_theta_nh3_no3, sd_theta_no3_n2,     &
                    & sd_theta_poc_ine
     real(8), save :: sd_theta_poc_lab, sd_theta_poc_ref, sd_theta_pon_ine,    &
                    & sd_theta_pon_lab, sd_theta_pon_ref, sd_theta_pop_ine,    &
                    & sd_theta_pop_lab, sd_theta_pop_ref, sd_theta_pw, sd_tic0,&
                    & sd_tsed, sd_tw, sd_w12, sd_w2, sd_xappd1, sd_xappp1,     &
                    & sd_xk1, sd_xk2, sedimentheat_src, sediment_heat_src,     &
                    & sedtemp, sedtemp1, sedtemp2, shields, slowfluxcarbon,    &
                    & so4consumedsnk1l2, so4producedsrc1l1, sqrs2,             &
                    & sulfided_sd1, sulfided_sd2, sulfideg_sd1, sulfideg_sd2,  &
                    & t1k, tau, ts, u2, uorb, volwater, vscour, xjn, xnh4,     &
                    & z1ss, z2ss, z3ss, z4ss, z5ass, z5bss, z5ss, z6ss
     real(8), dimension(numgenbodconstituents), save :: sedgenbodconc0
!
!*** End of declarations rewritten by SPAG
!
!!Local Variables
 
 
    ! CEMA testing variables start
    ! CEMA testing end
 
    ! sediment pH start
    ! sediment pH end
 
    !
    !  INPUTS
    !  SD_Jcin = flux to sediments from settling organic carbon
    !         from phytoplankton and detritus in oxygen equivalent units (gO2/m2/d)
    !         (NOTE: gO2/m2/d = gC/m2/d * 2.67 gO2/gC)                            [SDINC] in gC/m2/d from WATER QUALITY SEDIMENTC Subroutine
    !  SD_Jnin = nitrogen flux in settling phytoplankton and detritus (gN/m2/d)   [SDINN]
    !  SD_Jpin = phosphorus flux in settling phytoplankton and detritus (gP/m2/d) [SDINP]
    !  SD_O20 = dissolved oxygen in water overlying the sediment (mgO2/L)
    !  SD_depth = total water depth overlying the sediment (m) (used to calculate methane saturation concentration at in situ pressure)
 
    !  SD_Tw = temperature in water overlying the sediment (deg C)
    !  SD_NH30 = ammonia N in water overlying the sediment (mgN/L)
    !  SD_NO30 = nitrate N in water overlying the sediment (mgN/L)
    !  SD_PO40 = soluble reactive P in water overlying the sediment (mgP/L)
    !  SD_CH40 = fast reacting dissolved organic carbon and CBODu in the water overlying the sediment
    !         in oxygen equivalent units (mgO2/L)
    !         (NOTE: mgO2/L = mC/L * 2.67 mgO2/mgC)
    !  SD_SALw = salinity in the water overlying the sediment (ppt)
    !
    !  OUTPUTS
    !  SOD = sediment oxygen demand flux of dissolved oxygen between the water and sediment (gO2/m2/d)
    !        (positive is loss of O2 from water column)
    !  Jnh4 = flux of ammonia N between the water and sediment (gN/m2/d)
    !        (positive is source of NH4-N to water column)
    !  Jno3 = flux of nitrate N between the water and sediment (gN/m2/d)
    !        (positive is source of NO3-N to water column)
    !  Jch4 = flux of dissolved methane, fast reacting C, and CBODu between water and sediment in O2 equivalent units (gO2/m2/d)
    !
    !        (positive is source of CBOD to water column)
    !        (NOTE: gO2/m2/d = gC/m2/d * 2.67 gO2/gC)
    !        (methane is not produced in salt water)
    !  Jch4g = flux of methane gas bubbles between the water and sediment in O2 equivalent units (gO2/m2/d)
    !        (positive is source of CH4 bubbles to water column)
    !        (NOTE: gO2/m2/d = gC/m2/d * 2.67 gO2/gC)
    !        (methane is not produced in salt water)
    !  Jhs = flux of dissolved hydrogen sulfide (COD) between water and sediment in O2 equivalent units (gO2/m2/d)
    !        (positive is source of COD to water column)
    !        (hydrogen sulfide is not produced in freshwater)
    !  Jpo4 = flux of soluble reactive P between the water and sedmiment (gP/m2/d)
    !        (positive is source of PO4-P to water column)
    !  NH3(1) and NH3(2) = ammonia N in the sediment layers 1 and 2 (mgN/L)
    !  NO3(1) and NO3(2) = nitrate N in the sediment layers 1 and 2 (mgN/L)
    !  CH4(1) = dissolved methane in the aerobic sediment layer 1 (O2 equivalent units mgO2/L)
    !  HS(1) and HS(2) = dissolved sulfide in the sediment layers 1 and 2 (O2 equivalent units mgO2/L)
    !  PO4(1) and PO4(2) = soluble reactive P in the sediment layers 1 and 2 (mgP/L)
!!
!!loop start for each control volume
!!Set local default waterquality rates and constants
!!
!!Set initial time for averaging
 
     if(firsttimeincemamftseddiag)then                                !CellArea,   SW 10/3/2017
!
         allocate(cellarea(imx))
         gasreleasech4 = 0.0
                            ! SW 10/10/2017
        !GasReleaseCO2=0.0   ! SW 10/19/2017
        ! determining width of sediment cells  ! cb 6/11/16
        !Do JW=1, NWB            SW 9/18/2017
        !  Do JB=BS(JW),BE(JW)
        !    !IU = CUS(JB)
        !    IU = us(jb)
        !    ID = DS(JB)
        !    do SegNumI=iu,id
        !        sedcellwidth(SegNumI)=B(2,SegNumI)
        !    end do
        !  end do
        !end do
         call CEMAOUTPUTROUTINES  ! FIRST TIME CALL TO INITIALIZE VARIABLES  SW 3/9/16
!        !Get Rates and constants Region Numbers for each cell
         cemamft_randc_regn = 1
         do regnnum = 1, numregnssedimentdiagenesis
             do segnumi = SEDBEDDIARCREGSEGST(regnnum),                        &
               & SEDBEDDIARCREGSEGEN(regnnum)
                 cemamft_randc_regn(segnumi) = regnnum
             enddo     !SegNumI
         enddo     !RegnNum
!
!        !Get initial condition Region Numbers for each cell
         cemamft_incond_regn = 1
         do regnnum = 1, numregnssedimentbedcomposition
             do segnumi = SEDBEDINITREGSEGST(regnnum),                         &
               & SEDBEDINITREGSEGEN(regnnum)
                 cemamft_incond_regn(segnumi) = regnnum
             enddo     !SegNumI
         enddo     !RegnNum
!
!        !Initialize variables
         mftsedflxvars = 0.D00
         do jw = 1, nwb
             kt = KTWB(jw)
             do jb = BS(jw), BE(jw)
                !IU = CUS(JB)
                 iu = US(jb)
                 id = DS(jb)
                 do segnumi = iu, id
                     cellarea(segnumi) = B(KB(segnumi), segnumi)*DLX(segnumi)   ! SW 10/3/2017
                     regnnum = cemamft_randc_regn(segnumi)
                     initconregn = cemamft_incond_regn(segnumi)
 
                     call CEMAMFTRATESANDCONSTANTS
 
                     mftsedflxvars(segnumi, 16)                                &
                       & = sd_poc_l_fr*SDREGNPOC_T(initconregn)
                     mftsedflxvars(segnumi, 17)                                &
                       & = sd_poc_r_fr*SDREGNPOC_T(initconregn)
                     mftsedflxvars(segnumi, 18)                                &
                       & = sd_poc_i_fr*SDREGNPOC_T(initconregn)
                     mftsedflxvars(segnumi, 19)                                &
                       & = sd_pon_l_fr*SDREGNPON_T(initconregn)
                     mftsedflxvars(segnumi, 20)                                &
                       & = sd_pon_r_fr*SDREGNPON_T(initconregn)
                     mftsedflxvars(segnumi, 21)                                &
                       & = sd_pon_i_fr*SDREGNPON_T(initconregn)
                     mftsedflxvars(segnumi, 22)                                &
                       & = sd_pop_l_fr*SDREGNPOP_T(initconregn)
                     mftsedflxvars(segnumi, 23)                                &
                       & = sd_pop_r_fr*SDREGNPOP_T(initconregn)
                     mftsedflxvars(segnumi, 24)                                &
                       & = sd_pop_i_fr*SDREGNPOP_T(initconregn)
                     mftsedflxvars(segnumi, 25) = SDREGNSUL_T(initconregn)
                     SD_SO4CONC(segnumi) = SDREGNSUL_T(initconregn)
!                    !SD_pHValue(SegNumI)         =   CEMAPWpH
                     mftsedflxvars(segnumi, 7) = SDREGNCH4_T(initconregn)
!                    !MFTSedFlxVars(SegNumI,5)    =   SDRegnH2S_T(InitConRegn)
!                    !MFTSedFlxVars(SegNumI,6)    =   SDRegnH2S_T(InitConRegn)
!                    !MFTSedFlxVars(SegNumI,14)   =   SDRegnNH3_T(InitConRegn)
!                    !MFTSedFlxVars(SegNumI,15)   =   SDRegnNH3_T(InitConRegn)
                     mftsedflxvars(segnumi, 5) = SDREGNNH3_T(initconregn)          ! cb 2/18/13
                     mftsedflxvars(segnumi, 6) = SDREGNNH3_T(initconregn)
                     mftsedflxvars(segnumi, 14) = SDREGNH2S_T(initconregn)
                     mftsedflxvars(segnumi, 15) = SDREGNH2S_T(initconregn)
                     mftsedflxvars(segnumi, 31) = SDREGNTIC_T(initconregn)
                     mftsedflxvars(segnumi, 32) = SDREGNTIC_T(initconregn)
                     mftsedflxvars(segnumi, 33) = SDREGNALK_T(initconregn)
                     mftsedflxvars(segnumi, 34) = SDREGNALK_T(initconregn)
                     mftsedflxvars(segnumi, 10) = SDREGNPO4_T(initconregn)
                     mftsedflxvars(segnumi, 11) = SDREGNPO4_T(initconregn)
                     mftsedflxvars(segnumi, 35) = cemapwph
                     mftsedflxvars(segnumi, 36) = cemapwph
 
                     mftsedflxvars(segnumi, 37) = SDREGNFE2_T(initconregn)
                     mftsedflxvars(segnumi, 38) = SDREGNFE2_T(initconregn)
                     mftsedflxvars(segnumi, 39) = SDREGNFEOOH_T(initconregn)
                     mftsedflxvars(segnumi, 40) = SDREGNFEOOH_T(initconregn)
                     mftsedflxvars(segnumi, 41) = SDREGNMN2_T(initconregn)
                     mftsedflxvars(segnumi, 42) = SDREGNMN2_T(initconregn)
                     mftsedflxvars(segnumi, 43) = SDREGNMNO2_T(initconregn)
                     mftsedflxvars(segnumi, 44) = SDREGNMNO2_T(initconregn)
                     mftsedflxvars(segnumi, 45) = SDREGNT_T(initconregn)
                     mftsedflxvars(segnumi, 46) = SDREGNT_T(initconregn)
!
                     call UPDATECEMASEDIMENTFLUXVARIABLES   ! cb 2/13/15
 
                 enddo     !SegNumI
             enddo     !JB
         enddo !JW
 
         sd_w12 = partmixvel
                            ! Particle mixing velocity
         sd_w2 = burialvel  ! Burial Velocity
!
!        !Generic BOD code
         if(includecemagenbodconstituents)then
!            !Apply BOD constituent initial concentrations
             do i = 1, numgenbodconstituents
                 do regnnum = 1, numgenbodinitregns
                     do segnumi = SEDGENBODREGSEGST(regnnum),                  &
                       & SEDGENBODREGSEGEN(regnnum)
                         sedgenbodconc(i, segnumi, :)                          &
                           & = SEDGENBODINIT(i, regnnum)
                     enddo     !SegNum
                 enddo     !RegnNum
             enddo     !i
!
!            !Apply BOD constituent consumption rate
             do i = 1, numgenbodconstituents
                 do regnnum = 1, numgenbodconsumptionregions
                     do segnumi = SEDGENBODCONSREGSEGST(regnnum),              &
                       & SEDGENBODCONSREGSEGEN(regnnum)
                         SEDGENBODCONSRATE(i, segnumi)                         &
                           & = SEDGENBODREGNRATE(i, regnnum)
                         SEDGENBODCONSTCOEFF(i, segnumi)                       &
                           & = SEDGENBODREGNTCOEFF(i, regnnum)
                     enddo     !SegNum
                 enddo     !RegnNum
             enddo     !i
         endif
!        !End generic BOD code
!
!        !Open Output Snapshot File
         if(writecemamftsedflx)then
             open(cemasedflxfiln1, file = "CEMADiagenesisConstituent.opt")
             open(cemasedflxfiln2, file = "CEMADiagenesisSOD.opt")
             open(cemasedflxfiln3, file = "CEMADiagenesisAerobicLayer.opt")
             open(cemasedflxfiln4, file = "SOD.opt")
             write(cemasedflxfiln4, '(" SOD (gO2/m2/d)")')
             write(cemasedflxfiln4, '("    JDAY",<IMX>i10)')                   &
                 & (segnumi, segnumi = 1, imx)
             open(cemasedflxfiln5, file = "POC_Sediments.opt")
             write(cemasedflxfiln5, '(" POC (mg/l)")')
             write(cemasedflxfiln5, '("    JDAY",<IMX>i10)')                   &
                 & (segnumi, segnumi = 1, imx)
             open(cemasedflxfiln6, file = "PON_Sediments.opt")
             write(cemasedflxfiln6, '(" PON (mg/l)")')
             write(cemasedflxfiln6, '("    JDAY",<IMX>i10)')                   &
                 & (segnumi, segnumi = 1, imx)
             open(cemasedflxfiln7, file = "POP_Sediments.opt")
             write(cemasedflxfiln7, '(" POP (mg/l)")')
             write(cemasedflxfiln7, '("    JDAY",<IMX>i10)')                   &
                 & (segnumi, segnumi = 1, imx)
             open(cemasedflxfiln8, file = "CH4_Sediments.opt")
             write(cemasedflxfiln8, '(" Dissolved CH4 (mg/l)")')
             write(cemasedflxfiln8, '("    JDAY",<IMX>i10)')                   &
                 & (segnumi, segnumi = 1, imx)
             open(cemasedflxfiln9, file = "H2S_Sediments.opt")
             write(cemasedflxfiln9, '(" Dissolved H2S (mg/l)")')
             write(cemasedflxfiln9, '("    JDAY",<IMX>i10)')                   &
                 & (segnumi, segnumi = 1, imx)
            ! testing start
            !open(4478,file='CEMA_testing_output.opt',status='unknown')
            !write(4478,5598)
!5598        format("    JDAY","      SOD2      SOD3      SOD4     CSOD2    
!            CSOD3     CSOD4     NSOD2     NSOD3     NSOD4      POC2      POC3
!            POC4      PON2      PON3      PON4      SO42      SO43      SO44 
!5599        TIC2_L1   TIC2_L2   TIC3_L1   TIC3_L2   TIC4_L1   TIC4_L2  
!            ALK2_L1   ALK2_L2   ALK3_L1   ALK3_L2   ALK4_L1   ALK4_L2 
!            PO4T2_L1  PO4T2_L2  PO4T3_L1  PO4T3_L2  PO4T4_L1  PO4T4_L2  
!5601        PO42_L1   PO42_L2   PO43_L1   PO43_L2   PO44_L1   PO44_L2   
!            PH2_L1    PH2_L2    PH3_L1    PH3_L2    PH4_L1    PH4_L2
!            open(4479,file='CEMA_testing_output_metals.opt',status='unknown')
!5991        write(4479,5599) format("    JDAY","  FeII2_L1  FeII2_L2 
!            FeII3_L1  FeII3_L2  FeII4_L1  FeII4_L2 FeOOH2_L1 FeOOH2_L2
!            FeOOH3_L1 FeOOH3_L2 FeOOH4_L1 FeOOH4_L2  MnII2_L1  MnII2_L2 
!5992        MnII3_L1  MnII3_L2  MnII4_L1  MnII4_L2  MnO22_L1  MnO22_L2 
!            MnO23_L1  MnO23_L2  MnO24_L1  MnO24_L2")
 
! open(4480,file='CEMA_testing_output_temperature.opt',status='unknown') write(
         endif
         if(includecemagenbodconstituents)                                     &
           &open(cemaoutfiln5, file = "CEMAGenBOD.opt")
 
!        moved from SetupCEMASedimentModel
!Open    Output Snapshot File
         if(writebesnp .OR. writebesnp)then
             open(cemasnpoutfiln, file = "CEMAOutput.opt")
            !Open(CEMATSR1OutFilN, File = "CEMATSR1Output.opt")
             write(cemasnpoutfiln, '(a)')"$The file contains following output:"
            !Write(CEMATSR1OutFilN,'(a)')"$JDAY, KTWB(1), KB(15), z(15), BedElevation(15), BedElevationLayer(15), PorewaterRelRate(15), BedPorosity(15), BedConsolidRate(15)"
             if(writebesnp)write(cemasnpoutfiln, '(a)')"$Bed elevation (m)"
             if(writepwsnp)write(cemasnpoutfiln, '(a)')"$Bed porosity (%)"
         endif
 
        !Bed Elevation
         if(writebesnp)then
             write(cemasnpoutfiln, '("Bed elevation(m)")')
             write(cemasnpoutfiln, '("JDAY = ",f8.2)')jday
             write(cemasnpoutfiln, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
             write(cemasnpoutfiln, '(<IMX>f10.5)')                             &
                 & (BEDELEVATION(segnumi), segnumi = 1, imx)
 
             write(cemasnpoutfiln, '("Bed elevation layer(m)")')
             write(cemasnpoutfiln, '("JDAY = ",f8.2)')jday
             write(cemasnpoutfiln, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
             write(cemasnpoutfiln, '(<IMX>f10.5)')                             &
                 & (BEDELEVATIONLAYER(segnumi), segnumi = 1, imx)
         endif
 
        !Bed porosity
         if(writepwsnp)then
             write(cemasnpoutfiln, '("Bed porosity (%)")')
             write(cemasnpoutfiln, '("JDAY = ",f8.2)')jday
             write(cemasnpoutfiln, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
             write(cemasnpoutfiln, '(<IMX>f10.5)')                             &
                 & (BEDPOROSITY(segnumi)*100, segnumi = 1, imx)
             write(cemasnpoutfiln,                                             &
                  &'("Total Volume of Sediments = ",e14.6," m3")')             &
                 & totalsedimentsinbed
             write(cemasnpoutfiln, '("Total Porewater Volume = ",e14.6," m3")')&
                 & totalporewatvolume
             write(cemasnpoutfiln, '("Total Porewater Removed = ",e14.6," m3")'&
                 & )totalporewatremoved
         endif
 
        !Bottom Layer
         open(cemabtmlayfiln, file = "CEMABottomLayer.opt")
         write(cemabtmlayfiln, '(<IMX>i5)')(segnumi, segnumi = 1, imx)
         tempcntr1 = 0
         endbedconsolidation = .FALSE.
         cemalayeradded = .FALSE.
         cemassapplied = .FALSE.
 
        !Other files
         open(cemaoutfiln1, file = "CEMABubbles.opt")
         open(cemaoutfiln2, file = "CEMASedimentGas.opt")
         open(cemaoutfiln3, file = "CEMABubblesAtmosphereRelease.opt")
         open(cemaoutfiln6, file = "CEMASedDissGas.opt")
         open(cemaoutfilbub, file = "CEMABubbleReleaseSummary.csv")
         write(cemaoutfiln1, '(a)')"CEMA Sediment Bubbles Output"
         write(cemaoutfiln2, '(a)')"CEMA Sediment Bubbles Gas Output"
         write(cemaoutfiln3, '(a)')"CEMA Bubbles Release to Atmosphere Output"
         write(cemaoutfiln6, '(a)')"CEMA Sediment Dissolved Gas Output"
         write(cemaoutfilbub, *)                                               &
                     &'JDAY,WB,H2S(kg),CH4(kg),NH3(kg),CO2(kg),CH4(kg)-Method2'
        !Open(CEMAOutFilN4, File = "CEMABottomLayer.opt")
 
 
 
!
     endif
!
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
 
                 if(endbedconsolidation(segnumi))call ENDSEDIMENTFLUXVARIABLES
 
             enddo !SegNumI
         enddo !JB
     enddo !JW
!
     sd_tc = dlt/86400.0        !Time step in days
!
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
 
                 if(endbedconsolidation(segnumi))cycle
!
!	               !CellArea = B(KB(SegNumI),SegNumI)*DLX(SegNumI)    ! SW
                !CellArea =sedcellwidth(SegNumI)*DLX(SegNumI)
                !CellArea(SegNumI) = B(KB(SegNumI),SegNumI)*DLX(SegNumI)
!                10/3/2017 MOVED TO INITIAL LOOP
                 regnnum = cemamft_randc_regn(segnumi)
                 sd_h2 = BEDELEVATION(segnumi)
                 cellthickness = H1(KB(segnumi), segnumi)
 
                 call CEMAMFTRATESANDCONSTANTS
!
                 call CEMAMFTSEDIMENTFLUXMODEL
                !if(SegNumI == 3)then  ! writing out resuspension variables for testing
                  !if(nit == 1 .or. nit/10*10 == nit)then
                 !     if(SD_POMResuspMethod == 0)then
                 !       write(5387,'(f12.4,i12,30e12.4)')JDAY,nit,SD_E,FetchW,WIND(JW)*WSC(SegNumI),U2,HS,TS,LW,UORB,EPSILON,TAU,TAUCRPOM,DEPTHB(Kb(SegNumI),SegNumI),coef1,coef2,coef3,coef4
                 !     else
                 !       write(5388,'(f12.4,i12,30e12.4)')JDAY,nit,SD_E,c_bottom,c_bottom2,vscour,shields,crshield,smolvisc_h20,reyn_resusp
                 !     end if
                 ! end if
                !end if
 
             enddo     !SegNumI
         enddo !JB
     enddo !JW
 
     if(firsttimeincemamftseddiag)firsttimeincemamftseddiag = .FALSE.
     if(firsttimeinbubbles)firsttimeinbubbles = .FALSE.
 
     return
 
     entry CEMAMFTRATESANDCONSTANTS
 
     sd_poc_l_fr = SDREGNPOC_L_FR(regnnum)
     sd_poc_r_fr = SDREGNPOC_R_FR(regnnum)
     sd_poc_i_fr = 1 - sd_poc_l_fr - sd_poc_r_fr
     sd_pon_l_fr = SDREGNPON_L_FR(regnnum)
     sd_pon_r_fr = SDREGNPON_R_FR(regnnum)
     sd_pon_i_fr = 1 - sd_pon_l_fr - sd_pon_r_fr
     sd_pop_l_fr = SDREGNPOP_L_FR(regnnum)
     sd_pop_r_fr = SDREGNPOP_R_FR(regnnum)
     sd_pop_i_fr = 1 - sd_pop_l_fr - sd_pop_r_fr
     sd_pw_diffcoeff = SDREGNPW_DIFFCOEFF(regnnum)                   ! m^2/d
     sd_ox_threshol = SDREGNOX_THRESHOLD(regnnum)
     sd_ae_nh3_no3_l = SDREGNAE_NH3_NO3_L(regnnum)
     sd_ae_nh3_no3_h = SDREGNAE_NH3_NO3_H(regnnum)
     sd_ae_no3_n2_l = SDREGNAE_NO3_N2_L(regnnum)
     sd_ae_no3_n2_h = SDREGNAE_NO3_N2_H(regnnum)
     sd_an_no3_n2 = SDREGNAN_NO3_N2(regnnum)
     sd_ae_ch4_co2 = SDREGNAE_CH4_CO2(regnnum)
     sd_ae_hs_nh4_nit = SDREGNAE_HS_NH4_NIT(regnnum)
     sd_ae_hs_o2_nit = SDREGNAE_HS_O2_NIT(regnnum)
     sd_theta_pw = SDREGN_THETA_PW(regnnum)
     sd_theta_nh3_no3 = SDREGN_THETA_NH3_NO3(regnnum)
     sd_theta_no3_n2 = SDREGN_THETA_NO3_N2(regnnum)
     sd_theta_ch4_co2 = SDREGN_THETA_CH4_CO2(regnnum)
     sd_sulfate_ch4_h2s = SDREGN_SULFATE_CH4_H2S(regnnum)
     sd_ae_h2s_so4 = SDREGNAE_H2S_SO4(regnnum)
     sd_theta_h2s_so4 = SDREGN_THETA_H2S_SO4(regnnum)
     sd_normconst_h2s_so4 = SDREGN_NORMCONST_H2S_SO4(regnnum)
     sd_minrate_pon_lab = SDREGN_MINRATE_PON_LAB(regnnum)
     sd_minrate_pon_ref = SDREGN_MINRATE_PON_REF(regnnum)
     sd_minrate_pon_ine = SDREGN_MINRATE_PON_INE(regnnum)
     sd_minrate_poc_lab = SDREGN_MINRATE_POC_LAB(regnnum)
     sd_minrate_poc_ref = SDREGN_MINRATE_POC_REF(regnnum)
     sd_minrate_poc_ine = SDREGN_MINRATE_POC_INE(regnnum)
     sd_minrate_pop_lab = SDREGN_MINRATE_POP_LAB(regnnum)
     sd_minrate_pop_ref = SDREGN_MINRATE_POP_REF(regnnum)
     sd_minrate_pop_ine = SDREGN_MINRATE_POP_INE(regnnum)
     sd_theta_pon_lab = SDREGN_THETA_PON_LAB(regnnum)
     sd_theta_pon_ref = SDREGN_THETA_PON_REF(regnnum)
     sd_theta_pon_ine = SDREGN_THETA_PON_INE(regnnum)
     sd_theta_poc_lab = SDREGN_THETA_POC_LAB(regnnum)
     sd_theta_poc_ref = SDREGN_THETA_POC_REF(regnnum)
     sd_theta_poc_ine = SDREGN_THETA_POC_INE(regnnum)
     sd_theta_pop_lab = SDREGN_THETA_POP_LAB(regnnum)
     sd_theta_pop_ref = SDREGN_THETA_POP_REF(regnnum)
     sd_theta_pop_ine = SDREGN_THETA_POP_INE(regnnum)
     sd_ch4compmethod = SDREGN_CH4COMPMETHOD(regnnum)
     sd_pomresuspmethod = SDREGN_POMRESUSPMETHOD(regnnum)
 
     return
 
     entry CEMAMFTSEDIMENTFLUXMODEL
!
     call UPDATECEMASEDIMENTFLUXVARIABLES
!
!	!All fluxes from water column are zero
!
     fastfluxcarbon = 0.D0
     slowfluxcarbon = 0.D0
 
 
 ! CEMA testing start
 
     !--Compute jn !(C)
 !    xknh4= kappnh4 *thtanh4**(tempzz-20.)
 
    !--making the conversion factor for jn to sod a function of temperature
 !   ao2n = 106./16.*12./14.*32./12.*1.068**(tempzz-20.)
 !   o2ss = o2zz
 
    !--convert from mg/m2-d to g/m2-d
 !   xjnh4 = jnh4zz/1000.      ! jnh4zz is the ammonia flux to the sediments - this is the only term used to set all the other diagenesis terms
 !   xnh4= nh4zz/1000.         ! nh4zz is overlying water ammonia conc in mg/m3, this converts it to g/m3
 
    !--compute intermediate variables
 !   z1ss =xjnh4**2                              ! xjnh4: rate of ammonia flux to overlying water
 !   z2ss =xknh4**2
 !   z3ss =ao2n**2
 !   z4ss = 1./z3ss
 !   z5ss =o2ss**2
 !   z5bss = (-4.*ao2n*o2ss*xknh4**4*xnh4**3 &
 !      -z3ss*z1ss*z2ss*xnh4**2+18.*ao2n*o2ss*z1ss*z2ss*xnh4+ 27.*z5ss*z1ss*z2ss+4.*z3ss*xjnh4**4) !
 
    !--check for negative root. If this occurs, assume that the olw nh4=0
  !  if(z5bss .le. 0.0) xnh4=0.0
  !    z5ass= (z4ss*o2ss*xknh4*sqrt(-4.0*ao2n*o2ss*xknh4**4*xnh4**3-z3ss*z1ss*z2ss*xnh4**2+18.*ao2n*o2ss*z1ss*z2ss*xnh4+ &
  !    27.*z5ss*z1ss*z2ss+4.0*z3ss*xjnh4**4)/sqrt(3.)/6.0+z4ss*(9.*ao2n*o2ss*xjnh4*z2ss*xnh4+27.*z5ss*xjnh4*z2ss  &
  !    + 2.*z3ss*xjnh4**3)/54.0)                                ! Eq 16.5b
  !   z6ss = z5ass/abs(z5ass)*(abs(z5ass))**(1.0/3.0)                      ! Eq 16.5a
  !   xjn = z6ss+ (3.*o2ss*z2ss*xnh4+ao2n*z1ss)/(ao2n*z6ss)/9.0+xjnh4/3.0    ! xjn: rate at which organic matter is mineralized and ammonia is released to pore water
  !   acn= 106./16.*12./14.
  !   jcinzz= acn*xjn*1000.             !(D)   Convert ammonia flux to carbon diagenesis flux
 
  !   jnin= jcinzz/(frpon1+frpon2)/ratiocn
  !   jpin= jcinzz/(frpop1+frpop2)/ratiocp
 
 ! CEMA testing end
 
!	!
!	!Flux of caron
!	!SD_Jcin         = (FastFluxCarbon + SlowFluxCarbon)*2.67 !(gC/m2/sec to gO2/m2/sec)
        !SD_Jcin         =  sdinc(kb(segnumi),segnumi)*2.67*h(kb(segnumi),jw)          !  cb 2/28/13 (gC/m3/sec to gO2/m2/sec)
     sd_jcin = 0.0
        !do k=ktwb(jw),kb(segnumi)
     sd_jcin = SDINC(KB(segnumi), segnumi)*2.67*H(KB(segnumi), jw)                 !+SD_Jcin     !  cb 2/28/13 (gC/m3/sec to gO2/m2/sec)
        !end do
        !SD_Jcin=jcinzz  ! CEMA testing
!	!
!	!Flux of nitrogen
     sd_jnin = 0.D0
        !do k=ktwb(jw),kb(segnumi)
     sd_jnin = SDINN(KB(segnumi), segnumi)*H(KB(segnumi), jw)                          !+SD_Jnin
        !end do
!!SD_Jnin			= sdinn(kb(segnumi),segnumi)*h(kb(segnumi),jw)
    !    SD_Jnin=jnin   ! CEMA testing
!	!
!	!Flux of phosphorous
     sd_jpin = 0.D0
        !do k=ktwb(jw),kb(segnumi)  ! SW 8/31/2017
     sd_jpin = SDINP(KB(segnumi), segnumi)*H(KB(segnumi), jw)                          !+ SD_Jpin
        !end do
                !SD_Jpin			= sdinp(kb(segnumi),segnumi)*h(kb(segnumi),jw)
!    SD_Jpin=jpin
     sd_jfeoohin = SDINFEOOH(KB(segnumi), segnumi)*H(KB(segnumi), jw)
     sd_jmno2in = SDINMNO2(KB(segnumi), segnumi)*H(KB(segnumi), jw)
!
     sd_jcin = sd_jcin*86400.0D+00                            ! conversions not needed for testing...
     sd_jnin = sd_jnin*86400.0D+00
     sd_jpin = sd_jpin*86400.0D+00
     sd_jinfeooh = sd_jinfeooh*86400.0D+00
 
!	!Obtain properties from the water column
     sd_o20 = MAX(O2(KB(segnumi), segnumi), 0.0)                         !c(m1(ii,jj),k,I_DO)
        !SD_O20=o2zz  ! testing
     sd_depth = DEPTHB(KB(segnumi), segnumi)
     sd_tw = MAX(T1(KB(segnumi), segnumi), 0.0)                          !c(m1(ii,jj),k,I_Temp)
!    SD_Tw       = tempzz
     sd_nh30 = MAX(NH4(KB(segnumi), segnumi), 0.0)                       !c(m1(ii,jj),k,I_NH3)
!    SD_NH30     =nh4zz/1000.0   ! testing
     sd_no30 = MAX(NO3(KB(segnumi), segnumi), 0.0)                       !c(m1(ii,jj),k,I_NO3)
!    SD_NO30     =no3zz/1000.0
     sd_po40 = MAX(PO4(KB(segnumi), segnumi), 0.0)                       !c(m1(ii,jj),k,I_PO4)
!    SD_PO40     =po4zz/1000.0
!
     sd_ch40 = MAX(C1(KB(segnumi), segnumi, nch4), 0.0)
     sd_so40 = MAX(C1(KB(segnumi), segnumi, nso4), 0.0)
     sd_tic0 = MAX(C1(KB(segnumi), segnumi, ntic), 0.0)
     sd_alk0 = MAX(C1(KB(segnumi), segnumi, nalk), 0.0)
     sd_fe20 = MAX(C1(KB(segnumi), segnumi, nfe2), 0.0)
     sd_feooh0 = MAX(C1(KB(segnumi), segnumi, nfeooh), 0.0)
     sd_mn20 = MAX(C1(KB(segnumi), segnumi, nfe2), 0.0)
     sd_mno20 = MAX(C1(KB(segnumi), segnumi, nmno2), 0.0)
     sd_tsed = TSED(jw)
     sd_ksw = CBHE(jw)
     sd_taubot = SB(KB(segnumi), segnumi)/B(KB(segnumi), segnumi)
!
     if(includecemagenbodconstituents)then
         do genbodnum = 1, numgenbodconstituents
             sedgenbodconc0(genbodnum)                                         &
               & = MAX(CBOD(KB(segnumi), segnumi, genbodnum), 0.0)
         enddo   !GenBODNum
     endif
!
!	!SD_SO4      =   SD_SO4Conc(SegNumI)             !GET SULFATE VALUE IN SEDIMENT
!
     if(sd_o20<=1.0E-10)sd_o20 = 1.0E-10
!
     sd_rho = cemasedimentdensity
     sd_porosity = BEDPOROSITY(segnumi)
 
     do itemp = 1, 3
         sd_poct2 = sd_poct2 + SD_POC22(itemp)
         sd_pont2 = sd_pont2 + SD_PON22(itemp)
         sd_popt2 = sd_popt2 + SD_POP22(itemp)
     enddo
 
!	!Use steady state approach to initialize the values
     call DOCEMAMFTSEDIMENTDIAGEN
!
!	!
!	!Update source/sink terms
     DOSS(KB(segnumi), segnumi) = DOSS(KB(segnumi), segnumi) - dissolved_o2_snk                             !DO
     DOSEDIA(KB(segnumi), segnumi) = -dissolved_o2_snk    ! cb 3/19/13 storing for flux output
     NH4SS(KB(segnumi), segnumi) = NH4SS(KB(segnumi), segnumi)                 &
                                 & + dissolved_nh3_src                                                      !NH4
     NO3SS(KB(segnumi), segnumi) = NO3SS(KB(segnumi), segnumi)                 &
                                 & + dissolved_no3_src                                                      !NO3
     CGSS(KB(segnumi), segnumi, ngh2s) = CGSS(KB(segnumi), segnumi, ngh2s)     &
       & + dissolved_h2s_src                                                                                        !H2S
     CGSS(KB(segnumi), segnumi, ngch4) = CGSS(KB(segnumi), segnumi, ngch4)     &
       & + dissolved_ch4_src                                                                                        !CH4
     CGSS(KB(segnumi), segnumi, ngso4) = CGSS(KB(segnumi), segnumi, ngso4)     &
       & + dissolved_so4_src                                                                                !SO4
     CGSS(KB(segnumi), segnumi, ngfe2) = CGSS(KB(segnumi), segnumi, ngfe2)     &
       & + dissolved_fe2_src                                                                                !Fe(II)
     CGSS(KB(segnumi), segnumi, ngmn2) = CGSS(KB(segnumi), segnumi, ngmn2)     &
       & + dissolved_mn2_src                                                                                !Fe(II)
     TICSS(KB(segnumi), segnumi) = TICSS(KB(segnumi), segnumi)                 &
                                 & + dissolved_co2_src                                            !CO2
     ALKSS(KB(segnumi), segnumi) = ALKSS(KB(segnumi), segnumi)                 &
                                 & + dissolved_alk_src                                            !ALkalinity
     PO4SS(KB(segnumi), segnumi) = PO4SS(KB(segnumi), segnumi)                 &
                                 & + dissolved_po4_src                                            !PO4
     SDPFLUX(jw) = dissolved_po4_src*dlt*VOL(KB(segnumi), segnumi)             &
                 & /1000. + SDPFLUX(jw)                                                                                              ! SW 8/31/2017 kg
     SDNH4FLUX(jw) = dissolved_nh3_src*dlt*VOL(KB(segnumi), segnumi)           &
                   & /1000. + SDNH4FLUX(jw)
     SDNO3FLUX(jw) = dissolved_no3_src*dlt*VOL(KB(segnumi), segnumi)           &
                   & /1000. + SDNO3FLUX(jw)
     if(includecemagenbodconstituents)then
         do genbodnum = 1, numgenbodconstituents
             CBODSS(KB(segnumi), segnumi, genbodnum)                           &
               & = CBODSS(KB(segnumi), segnumi, genbodnum)                     &
               & + SEDGENBODSS(genbodnum)                                                                                   !BOD
         enddo   !GenBODNum
     endif
     LPOMSS(KB(segnumi), segnumi) = LPOMSS(KB(segnumi), segnumi)               &
                                  & + lpom_resuspension
     RPOMSS(KB(segnumi), segnumi) = RPOMSS(KB(segnumi), segnumi)               &
                                  & + rpom_resuspension
     LPOMPSS(KB(segnumi), segnumi) = LPOMPSS(KB(segnumi), segnumi)             &
                                   & + lpomp_resuspension
     RPOMPSS(KB(segnumi), segnumi) = RPOMPSS(KB(segnumi), segnumi)             &
                                   & + rpomp_resuspension
     LPOMNSS(KB(segnumi), segnumi) = LPOMNSS(KB(segnumi), segnumi)             &
                                   & + lpomp_resuspension
     RPOMNSS(KB(segnumi), segnumi) = RPOMNSS(KB(segnumi), segnumi)             &
                                   & + rpomp_resuspension
 
     TSS(KB(segnumi), segnumi) = TSS(KB(segnumi), segnumi) + sedimentheat_src                !Heat
 
     call WRITECEMASEDIMENTFLUXVARIABLES
!
     return
 
     entry UPDATECEMASEDIMENTFLUXVARIABLES
 
     SD_NO3P2(1) = mftsedflxvars(segnumi, 1)
     SD_NO3P2(2) = mftsedflxvars(segnumi, 2)
     SD_NH3P2(1) = mftsedflxvars(segnumi, 3)
     SD_NH3P2(2) = mftsedflxvars(segnumi, 4)
     SD_NH3TP2(1) = mftsedflxvars(segnumi, 5)
     SD_NH3TP2(2) = mftsedflxvars(segnumi, 6)
     SD_CH4P2(1) = mftsedflxvars(segnumi, 7)
!	!SD_PO4p2(1)		=	MFTSedFlxVars(SegNumI,8)
!	!SD_PO4p2(2)		=	MFTSedFlxVars(SegNumI,9)
     SD_PO4(1) = mftsedflxvars(segnumi, 8)
     SD_PO4(2) = mftsedflxvars(segnumi, 9)
!	!SD_PO4Tp2(1)	=	MFTSedFlxVars(SegNumI,10)
!	!SD_PO4Tp2(2)	=	MFTSedFlxVars(SegNumI,11)
     SD_PO4T(1) = mftsedflxvars(segnumi, 10)
     SD_PO4T(2) = mftsedflxvars(segnumi, 11)
     SD_HSP2(1) = mftsedflxvars(segnumi, 12)
     SD_HSP2(2) = mftsedflxvars(segnumi, 13)
     SD_HSTP2(1) = mftsedflxvars(segnumi, 14)
     SD_HSTP2(2) = mftsedflxvars(segnumi, 15)
 
 
     SD_POC22(1) = mftsedflxvars(segnumi, 16)
     SD_POC22(2) = mftsedflxvars(segnumi, 17)
     SD_POC22(3) = mftsedflxvars(segnumi, 18)
     SD_PON22(1) = mftsedflxvars(segnumi, 19)
     SD_PON22(2) = mftsedflxvars(segnumi, 20)
     SD_PON22(3) = mftsedflxvars(segnumi, 21)
     SD_POP22(1) = mftsedflxvars(segnumi, 22)
     SD_POP22(2) = mftsedflxvars(segnumi, 23)
     SD_POP22(3) = mftsedflxvars(segnumi, 24)
     sd_so4 = mftsedflxvars(segnumi, 25)
!
     sd_jsod = mftsedflxvars(segnumi, 26)
     sd_csod = mftsedflxvars(segnumi, 27)
     sd_nsod = mftsedflxvars(segnumi, 28)
 
     SD_TIC(1) = mftsedflxvars(segnumi, 31)
     SD_TIC(2) = mftsedflxvars(segnumi, 32)
     SD_ALK(1) = mftsedflxvars(segnumi, 33)
     SD_ALK(2) = mftsedflxvars(segnumi, 34)
     SD_PH(1) = mftsedflxvars(segnumi, 35)
     SD_PH(2) = mftsedflxvars(segnumi, 36)
     SD_FE2T(1) = mftsedflxvars(segnumi, 37)
     SD_FE2T(2) = mftsedflxvars(segnumi, 38)
     SD_FEOOH(1) = mftsedflxvars(segnumi, 39)
     SD_FEOOH(2) = mftsedflxvars(segnumi, 40)
     SD_MN2T(1) = mftsedflxvars(segnumi, 41)
     SD_MN2T(2) = mftsedflxvars(segnumi, 42)
     SD_MNO2(1) = mftsedflxvars(segnumi, 43)
     SD_MNO2(2) = mftsedflxvars(segnumi, 44)
     SD_T(1) = mftsedflxvars(segnumi, 45)
     SD_T(2) = mftsedflxvars(segnumi, 46)
 
     return
!
     entry WRITECEMASEDIMENTFLUXVARIABLES
!	!
!	!  NH3(1) and NH3(2) = ammonia N in the sediment layers 1 and 2 (mgN/L)
!	!  NO3(1) and NO3(2) = nitrate N in the sediment layers 1 and 2 (mgN/L)
!	!  CH4(1) = dissolved methane in the aerobic sediment layer 1 (O2 equivalent
!	!  units mgO2/L) HS(1) and HS(2) = dissolved sulfide in the sediment layers
!	!  1 and 2 (O2 equivalent units mgO2/L) PO4(1) and PO4(2) = soluble reactive
!    P in the sediment layers 1 and 2 (mgP/L)
     mftsedflxvars(segnumi, 1) = SD_NO3P2(1)
     mftsedflxvars(segnumi, 2) = SD_NO3P2(2)
     mftsedflxvars(segnumi, 3) = SD_NH3P2(1)
     mftsedflxvars(segnumi, 4) = SD_NH3P2(2)
     mftsedflxvars(segnumi, 5) = SD_NH3TP2(1)
     mftsedflxvars(segnumi, 6) = SD_NH3TP2(2)
     mftsedflxvars(segnumi, 7) = SD_CH4P2(1)
     mftsedflxvars(segnumi, 8) = SD_PO4(1)
     mftsedflxvars(segnumi, 9) = SD_PO4(2)
     mftsedflxvars(segnumi, 10) = SD_PO4T(1)
     mftsedflxvars(segnumi, 11) = SD_PO4T(2)
     mftsedflxvars(segnumi, 12) = SD_HSP2(1)
     mftsedflxvars(segnumi, 13) = SD_HSP2(2)
     mftsedflxvars(segnumi, 14) = SD_HSTP2(1)
     mftsedflxvars(segnumi, 15) = SD_HSTP2(2)
 
     mftsedflxvars(segnumi, 16) = SD_POC22(1)
     mftsedflxvars(segnumi, 17) = SD_POC22(2)
     mftsedflxvars(segnumi, 18) = SD_POC22(3)
     mftsedflxvars(segnumi, 19) = SD_PON22(1)
     mftsedflxvars(segnumi, 20) = SD_PON22(2)
     mftsedflxvars(segnumi, 21) = SD_PON22(3)
     mftsedflxvars(segnumi, 22) = SD_POP22(1)
     mftsedflxvars(segnumi, 23) = SD_POP22(2)
     mftsedflxvars(segnumi, 24) = SD_POP22(3)
     mftsedflxvars(segnumi, 25) = sd_so4
!
     mftsedflxvars(segnumi, 26) = sd_jsod
     mftsedflxvars(segnumi, 27) = sd_csod
     mftsedflxvars(segnumi, 28) = sd_nsod
!
     mftsedflxvars(segnumi, 29) = SD_POC22(1) + SD_POC22(2) + SD_POC22(3)
     mftsedflxvars(segnumi, 30) = SD_PON22(1) + SD_PON22(2) + SD_PON22(3)
     mftsedflxvars(segnumi, 31) = SD_TIC(1)
     mftsedflxvars(segnumi, 32) = SD_TIC(2)
     mftsedflxvars(segnumi, 33) = SD_ALK(1)
     mftsedflxvars(segnumi, 34) = SD_ALK(2)
     mftsedflxvars(segnumi, 35) = SD_PH(1)
     mftsedflxvars(segnumi, 36) = SD_PH(2)
     mftsedflxvars(segnumi, 37) = SD_FE2T(1)
     mftsedflxvars(segnumi, 38) = SD_FE2T(2)
     mftsedflxvars(segnumi, 39) = SD_FEOOH(1)
     mftsedflxvars(segnumi, 40) = SD_FEOOH(2)
     mftsedflxvars(segnumi, 41) = SD_MN2T(1)
     mftsedflxvars(segnumi, 42) = SD_MN2T(2)
     mftsedflxvars(segnumi, 43) = SD_MNO2(1)
     mftsedflxvars(segnumi, 44) = SD_MNO2(2)
     mftsedflxvars(segnumi, 45) = SD_T(1)
     mftsedflxvars(segnumi, 46) = SD_T(2)
     mftsedflxvars(segnumi, 47) = SD_POP22(1) + SD_POP22(2) + SD_POP22(3)
!
!	!For output
!
     CEMA_SD_VARS(segnumi, 1) = dissolved_nh3_src*cellthickness*86400.0
     CEMA_SD_VARS(segnumi, 2) = dissolved_no3_src*cellthickness*86400.0
     CEMA_SD_VARS(segnumi, 3) = dissolved_ch4_src*cellthickness*86400.0
     CEMA_SD_VARS(segnumi, 4) = dissolved_so4_src*cellthickness*86400.0
     CEMA_SD_VARS(segnumi, 5) = dissolved_co2_src*cellthickness*86400.0
     CEMA_SD_VARS(segnumi, 6) = dissolved_h2s_src*cellthickness*86400.0
     CEMA_SD_VARS(segnumi, 7) = dissolved_o2_snk*cellthickness*86400.0
!
     CEMA_SD_VARS(segnumi, 8) = sd_poct2
     CEMA_SD_VARS(segnumi, 9) = sd_pont2
     CEMA_SD_VARS(segnumi, 10) = sd_so4
!
     return
!
     entry DOCEMAMFTSEDIMENTDIAGEN
!!Adapted from original work by Steve Chapra in November 2003 version of QUAL2K.
    !Equation references are to DiToro 2001. Sediment Flux Modeling. Wiley-Interscience.
    !Greg Pelletier re-organized the specification of the rate constants and
    !added salinity-dependent nitrification/denitrification, inorganic P partitioning,
    !and sulfide production/flux from carbon diagenesis based on the work of DiToro 2001.
    !Significant modifications by Pelletier are marked by 'gp
    !
    !  OUTPUTS
    !  SOD = sediment oxygen demand flux of dissolved oxygen between the water and sediment (gO2/m2/d)
    !        (positive is loss of O2 from water column)
    !  Jnh4 = flux of ammonia N between the water and sediment (gN/m2/d)
    !        (positive is source of NH4-N to water column)
    !  Jno3 = flux of nitrate N between the water and sediment (gN/m2/d)
    !        (positive is source of NO3-N to water column)
    !  Jch4 = flux of dissolved methane, fast reacting C, and CBODu between water and sediment in O2 equivalent units (gO2/m2/d)
    !
    !        (positive is source of CBOD to water column)
    !        (NOTE: gO2/m2/d = gC/m2/d * 2.67 gO2/gC)
    !        (methane is not produced in salt water)
    !  Jch4g = flux of methane gas bubbles between the water and sediment in O2 equivalent units (gO2/m2/d)
    !        (positive is source of CH4 bubbles to water column)
    !        (NOTE: gO2/m2/d = gC/m2/d * 2.67 gO2/gC)
    !        (methane is not produced in salt water)
    !  Jhs = flux of dissolved hydrogen sulfide (COD) between water and sediment in O2 equivalent units (gO2/m2/d)
    !        (positive is source of COD to water column)
    !        (hydrogen sulfide is not produced in freshwater)
    !  Jpo4 = flux of soluble reactive P between the water and sedmiment (gP/m2/d)
    !        (positive is source of PO4-P to water column)
    !  NH3(1) and NH3(2) = ammonia N in the sediment layers 1 and 2 (mgN/L)
    !  NO3(1) and NO3(2) = nitrate N in the sediment layers 1 and 2 (mgN/L)
    !  CH4(1) = dissolved methane in the aerobic sediment layer 1 (O2 equivalent units mgO2/L)
    !  HS(1) and HS(2) = dissolved sulfide in the sediment layers 1 and 2 (O2 equivalent units mgO2/L)
    !  PO4(1) and PO4(2) = soluble reactive P in the sediment layers 1 and 2 (mgP/L)
    !
    !  REFERENCES:
    !
    !  Cerco, C.F. and T. Cole. 1995. User's guide to the CE-QUAL-ICM three-dimensional eutrophication model.
    !  Release version 1.0. U.S. Army Corps of Engineers. Waterways Experiment Station. Vicksburg, MS.
    !
    !  Chapra, S.C. and Pelletier, G.J. 2003. QUAL2K: A Modeling Framework for Simulating River and Stream
    !  Water Quality (Beta Version): Documentation and Users Manual.
    !  Civil and Environmental Engineering Dept., Tufts University, Medford, MA.
    !
    !  DiToro, D.M. 2001. Sediment Flux Modeling. Wiley-Interscience. New York, NY.
    !***** Assign rate constants (NOTE: these could be available for adjusting during model calibration)
    !reaction velocities (m/d)
!	!
!	!gp use salinity threshold to assign constants for the
!	!  nitrification/denitrification velocities for aerobic layer 1
!	!  and factor to increase aerobic layer partition for inorganic P
!
     if(sd_o20>sd_ox_threshol)then
         sd_ae_nh3_no3 = sd_ae_nh3_no3_h
         sd_ae_no3_n2 = sd_ae_no3_n2_h
     else
         sd_ae_nh3_no3 = sd_ae_nh3_no3_l
         sd_ae_no3_n2 = sd_ae_no3_n2_l
     endif
 
!	!
!	!gp
     do itemp = 1, 3
         SD_POC2(itemp) = SD_POC22(itemp)
         SD_PON2(itemp) = SD_PON22(itemp)
         SD_POP2(itemp) = SD_POP22(itemp)
     enddo
     do itemp = 1, 2
         SD_NH3TP(itemp) = SD_NH3TP2(itemp)
         SD_NO3P(itemp) = SD_NO3P2(itemp)
!	!	SD_PO4Tp(iTemp) = SD_PO4Tp2(iTemp)
         SD_HSTP(itemp) = SD_HSTP2(itemp)
     enddo
     sd_ben_strp = sd_ben_strp2
!
!	!gp assign constants for G class 1 and 2 PON, POC, and POP and calculate G class 3 as 1-fpox1-fpox2
     SD_FPON(1) = sd_pon_l_fr
     SD_FPON(2) = sd_pon_r_fr
     SD_FPON(3) = sd_pon_i_fr
     SD_FPOC(1) = sd_poc_l_fr
     SD_FPOC(2) = sd_poc_r_fr
     SD_FPOC(3) = sd_poc_i_fr
     SD_FPOP(1) = sd_pop_l_fr
     SD_FPOP(2) = sd_pop_r_fr
     SD_FPOP(3) = sd_pop_i_fr
!
!	!
!	!gp assign constants for G class 1, 2, and 3 mineralization of PON, POC, POP
     SD_KDIAPON(1) = sd_minrate_pon_lab
     SD_THTAPON(1) = sd_theta_pon_lab
     SD_KDIAPON(2) = sd_minrate_pon_ref
     SD_THTAPON(2) = sd_theta_pon_ref
     SD_KDIAPON(3) = sd_minrate_pon_ine
     SD_THTAPON(3) = sd_theta_pon_ine
     SD_KDIAPOC(1) = sd_minrate_poc_lab
     SD_THTAPOC(1) = sd_theta_poc_lab
     SD_KDIAPOC(2) = sd_minrate_poc_ref
     SD_THTAPOC(2) = sd_theta_poc_ref
     SD_KDIAPOC(3) = sd_minrate_poc_ine
     SD_THTAPOC(3) = sd_theta_poc_ine
     SD_KDIAPOP(1) = sd_minrate_pop_lab
     SD_THTAPOP(1) = sd_theta_pop_lab
     SD_KDIAPOP(2) = sd_minrate_pop_ref
     SD_THTAPOP(2) = sd_theta_pop_ref
     SD_KDIAPOP(3) = sd_minrate_pop_ine
     SD_THTAPOP(3) = sd_theta_pop_ine
!	!
!	!Compute input fluxes
     do itemp = 1, 3
         SD_JPOC(itemp) = sd_jcin*SD_FPOC(itemp)
         SD_JPON(itemp) = sd_jnin*SD_FPON(itemp)
         SD_JPOP(itemp) = sd_jpin*SD_FPOP(itemp)
     enddo
!	!
        ! computing resuspended POM
     sd_pom = sd_poct2/ORGC(jw)
     if(cema_pom_resuspension_processes)then
         if(sd_pomresuspmethod==0)then
             call CEMAWINDINDUCEDSEDIMENTRESUSPENSION
         else
             call CEMABOTTOMSCOURRESUSPENSION
         endif
 
         sd_epoc = 0.0
         sd_epon = 0.0
         sd_epop = 0.0
         do itemp = 1, 3
             if(sd_pom>nonzero .AND. sd_poct2>nonzero)sd_epoc(itemp)           &
              & = sd_poct2/sd_pom*SD_POC2(itemp)/sd_poct2*sd_e
             if(sd_pom>nonzero .AND. sd_pont2>nonzero)sd_epon(itemp)           &
              & = sd_pont2/sd_pom*SD_PON2(itemp)/sd_pont2*sd_e
             if(sd_pom>nonzero .AND. sd_popt2>nonzero)sd_epop(itemp)           &
              & = sd_popt2/sd_pom*SD_POP2(itemp)/sd_popt2*sd_e
         enddo
     endif
 
!	!Compute particulate organic forms
     sd_poct2 = 0.0
     sd_pont2 = 0.0
     sd_popt2 = 0.0
!
!	!Equation 13.32 DiToro. See also Equation 12.2 - also includes resuspension
     do itemp = 1, 3
         SD_POC2(itemp) = (SD_POC2(itemp) + SD_JPOC(itemp)*sd_tc/sd_h2 -       &
                        & sd_epoc(itemp)*sd_tc/sd_h2)                          &
                        & /(1. + SD_KDIAPOC(itemp)*SD_THTAPOC(itemp)           &
                        & **(SD_T(2) - 20.)*sd_tc)
         SD_PON2(itemp) = (SD_PON2(itemp) + SD_JPON(itemp)*sd_tc/sd_h2 -       &
                        & sd_epon(itemp)*sd_tc/sd_h2)                          &
                        & /(1. + SD_KDIAPON(itemp)*SD_THTAPON(itemp)           &
                        & **(SD_T(2) - 20.)*sd_tc)
         SD_POP2(itemp) = (SD_POP2(itemp) + SD_JPOP(itemp)*sd_tc/sd_h2 -       &
                        & sd_epop(itemp)*sd_tc/sd_h2)                          &
                        & /(1. + SD_KDIAPOP(itemp)*SD_THTAPOP(itemp)           &
                        & **(SD_T(2) - 20.)*sd_tc)
         sd_poct2 = sd_poct2 + SD_POC2(itemp)
         sd_pont2 = sd_pont2 + SD_PON2(itemp)
         sd_popt2 = sd_popt2 + SD_POP2(itemp)
     enddo
!
!	!
!	!Compute diagenesis fluxes
!	!Equation 13.31 diagenesis term only. See also Equation 12.2, 12.5 and 12.6
     sd_jc = 0
     sd_jn = 0
     sd_jp = 0
     do itemp = 1, 3
         sd_jc = sd_jc + SD_KDIAPOC(itemp)*SD_THTAPOC(itemp)**(SD_T(2) - 20.)  &
               & *SD_POC2(itemp)*sd_h2
         sd_jn = sd_jn + SD_KDIAPON(itemp)*SD_THTAPON(itemp)**(SD_T(2) - 20.)  &
               & *SD_PON2(itemp)*sd_h2
         sd_jp = sd_jp + SD_KDIAPOP(itemp)*SD_THTAPOP(itemp)**(SD_T(2) - 20.)  &
               & *SD_POP2(itemp)*sd_h2
     enddo
     sd_jctest = sd_jc
 
     maxit = 1000
     sd_es = 0.001
!
!	!Calculation of final benthic particle mixing velocity Equation 13.6 in DoToro (2001)
!	!See also Equation 2.56 and 4.48
     sd_kl12 = sd_pw_diffcoeff*sd_theta_pw**(SD_T(2) - 20.)/(sd_h2/2.)
!
     if(firsttimeincemamftseddiag)then
         sd_sodold = sd_jc + 1.714*sd_jn
     else
         sd_sodold = sd_jsod
     endif
!
     iter = 0
!	!Saturation conc. of methane in oxygen equivalent units (Equation 10.51)
     sd_ch4sat = 100.0D+00*(1.0D+00 + sd_depth/10.0D+00)                       &
               & *(1.024**(20.0D+00 - SD_T(2)))                                                          ![gmO*/m3]
 
     do while (.TRUE.)
!		!CEMA Additional Code
!		!Initialize variables to zero
         co2producedsrc1l1 = 0.D0
         so4producedsrc1l1 = 0.D0
         so4consumedsnk1l2 = 0.D0
         co2producedsrc2l2 = 0.D0
!		!End  CEMA Additional Code
!
!		!Equation 3.15 or 10.36
         if(sd_sodold>nonzero)then
             sd_s = sd_sodold/sd_o20
         else
             sd_sodold = sd_jc + 1.714*sd_jn      ! needed when segments are added
             sd_s = sd_sodold/sd_o20
         endif
         if(sd_s<nonzero)sd_s = nonzero
!		!Equation 4.51
         sd_nh3tono3 = sd_ae_nh3_no3**2.0*sd_theta_nh3_no3**(SD_T(1) - 20.)    &
                     & /sd_s*sd_ae_hs_nh4_nit/(sd_ae_hs_nh4_nit + SD_NH3(1))   &
                     & *sd_o20/(2.0*sd_ae_hs_o2_nit + sd_o20)                                                                                                                        ! [m/d] = [m2/d2] / [m/d] * [-] * [-]
!		!     cb note - SD_Ae_HS_NH4_Nit in above line may need to be temperature
!		!Calculate dissolved and particulate (sorbed) fractions
!        corrected
         sd_fd1 = (1.0D+00/(1.0D+00 + sd_m1*sd_kdnh3))                               ! cb note -SD_KdNH3 is not defined anywhere
         sd_fp1 = 1.0D+00 - sd_fd1                              != ((m1*KdNH3)/(1 + m1*KdNH3))
         sd_fd2 = (1.0D+00/(1.0D+00 + sd_m2*sd_kdnh3))                              ! cb note -SD_KdNH3 is not defined anywhere
         sd_fp2 = 1.0D+00 - sd_fd2                              != ((m2*KdNH3)/(1 + m2*KdNH3))
!		!     cb note -  since SD_KdNH3 not defined (=0), no nh3 is sorbed to
!		!Write linear system of equations around NH3T
!		!Equation 5.1
!		!     particualtes
!		!Layer 1
         sd_a11 = -sd_fd1*sd_kl12 - sd_fp1*sd_w12 - sd_fd1*sd_nh3tono3 -       &
                & sd_fd1*sd_s - sd_w2
            !SD_a11 = -SD_H1 / SD_tc - SD_fd1*SD_KL12 - SD_fp1*SD_w12 - SD_fd1*SD_NH3toNO3 - SD_fd1*SD_s - SD_w2
         sd_a12 = sd_fd2*sd_kl12 + sd_fp2*sd_w12
         sd_b1 = -sd_s*sd_nh30                                  ![m/d]*[mg/m3]
!		!
!		!Layer 2
         sd_a21 = sd_fd1*sd_kl12 + sd_fp1*sd_w12 + sd_w2
         sd_a22 = -sd_fd2*sd_kl12 - sd_fp2*sd_w12 - sd_w2 - sd_h2/sd_tc
         sd_b2 = -sd_jn - sd_h2/sd_tc*SD_NH3TP(2)
!
         call LIN_SYS(sd_a11, sd_a12, sd_a21, sd_a22, sd_b1, sd_b2, SD_NH3T(1),&
                    & SD_NH3T(2), nflog, nfcle)
!		!
!		!Dissolved Concentrations
         SD_NH3(1) = sd_fd1*SD_NH3T(1)
         SD_NH3(2) = sd_fd2*SD_NH3T(2)
!
!		!CEMA Additional Code
!		!Calculate NH4+ and NH3 concentrations
!		!
!		!SD1_Ammonia = SD_NH3(1)/(1. + 10.**(-SD_pHValue(SegNumI))/10.**(-NH4_NH3_Eqb_Const))
         sd1_ammonia = SD_NH3(1)                                               &
                     & /(1. + 10.**( - SD_PH(1))/10.**( - nh4_nh3_eqb_const))
         sd1_ammonium = SD_NH3(1) - sd1_ammonia
!		!SD2_Ammonia = SD_NH3(2)/(1. + 10.**(-SD_pHValue(SegNumI))/10.**(-NH4_NH3_Eqb_Const))
         sd2_ammonia = SD_NH3(2)                                               &
                     & /(1. + 10.**( - SD_PH(2))/10.**( - nh4_nh3_eqb_const))
         sd2_ammonium = SD_NH3(2) - sd2_ammonia
!		!Calculate dissolved <--> gaseous phase distribution
         sedtemp1 = SD_T(1) + 273.15
         sedtemp2 = SD_T(2) + 273.15
         mw_constituent = 17                !N = 14, H = 1
            !gp estimated thickness of the aerobic sediment layer 1 (DiToro Appendix B Page 576)  MOVED from below, cb 9/6/13
         sd_h1 = sd_kl12*sd_h2/sd_s
         if(sd_h1>sd_h2)sd_h1 = sd_h2
         SD_AERLAYERTHICK(segnumi) = sd_h1
!		!Layer 1
         volwater = cellarea(segnumi)*sd_h1
         call CEMADISGASPHASEDISTRIBUTION(sd1_ammonia, mw_constituent,         &
           & sedtemp1, henryconst_nh3, volwater, ammoniag_sd1, ammoniad_sd1)
!		!Layer 2
         volwater = cellarea(segnumi)*sd_h2
         call CEMADISGASPHASEDISTRIBUTION(sd2_ammonia, mw_constituent,         &
           & sedtemp2, henryconst_nh3, volwater, ammoniag_sd2, ammoniad_sd2)
!		!End  CEMA Additional Code
!
!		!
!		!Oxygen Flux due to NH3->NO2, see Equation 23.3 in Chapra (1997)
         sd_nh3conv = sd_nh3tono3*SD_NH3(1)
         sd_nsod = 2*(32.0/14.0)*sd_nh3tono3*SD_NH3(1)
!		![gmO/m2-d] = [mol O2/mol N]*[gm O2/mol O2]/[gm N/mol N]*
!		![gm/1000mg] * [m/day] * [mgN/m3]
!		!':::::::::::::::::::::::::::: BEGIN Nitrate::::::::::::::::::::::::::::
!		!
!		!Denitrification in layers 1 and 2 (Equation 4.55)
         SD_DENIT(1) = (sd_ae_no3_n2**2*sd_theta_no3_n2**(SD_T(1) - 20.)/sd_s)
         SD_DENIT(2) = sd_an_no3_n2*sd_theta_no3_n2**(SD_T(2) - 20.)
!		!
!		!Layer 1
         sd_a11 = -sd_kl12 - SD_DENIT(1) - sd_s - sd_w2
            !SD_a11 = -SD_H1 / SD_tc - SD_KL12 - SD_Denit(1) - SD_s - SD_w2
         sd_a12 = sd_kl12
         sd_b1 = -sd_s*sd_no30 - sd_nh3tono3*SD_NH3(1)
!		!Layer 2
         sd_a21 = sd_kl12 + sd_w2
         sd_a22 = -sd_kl12 - SD_DENIT(2) - sd_w2 - sd_h2/sd_tc
         sd_b2 = -sd_h2/sd_tc*SD_NO3P(2)
!
!		!
         call LIN_SYS(sd_a11, sd_a12, sd_a21, sd_a22, sd_b1, sd_b2, SD_NO3(1), &
                    & SD_NO3(2), nflog, nfcle)
!		!
!		!Nitrate Flux to water column
         sd_jno3 = sd_s*(SD_NO3(1) - sd_no30)
!		!
!		!Denitrification Flux [mgN/m2d]
         SD_JDENIT(1) = SD_DENIT(1)*SD_NO3(1)
         SD_JDENIT(2) = SD_DENIT(2)*SD_NO3(2)
         sd_jdenitt = SD_JDENIT(1) + SD_JDENIT(2)
!		!
!		!Methane consumption due to denitrification (Equation 9.16)
!		!
!		!Layer 1
!		!2    for 1/2 N2, 10/8 for eqn balance, 16/12 for O2 --> C and 12/14 for N
!        --> C
         SD_JO2NO3(1) = 2.0D+00*(16.0D+00/12.0D+00)*(10.0D+00/8.0D+00)         &
                      & *(12.0D+00/14.0D+00)*SD_JDENIT(1)
!		![gmO*/m2-d] = [molO/molC]*[(gmO/molO)/(gmC/molC)]*[molC/MolN]*
!		![(gmC/molC)/(gmN/molN)] * [mg/m2d] * [gm/1000mg]
!		!where 2*16/12 is the ubiquitous 32/12 (= 2.67) for oxidation of carbon
!		!Layer 2
         SD_JO2NO3(2) = (32.0D+00/12.0D+00)*(10.0D+00/8.0D+00)                 &
                      & *(12.0D+00/14.0D+00)*SD_JDENIT(2)
!		!
!		!Sum
         sd_jo2no3t = SD_JO2NO3(1) + SD_JO2NO3(2)
!		!
!		!Calculate methane flux in oxygen equivalent units, adjusted for
!		!the  methane consumed in denitrification
!		!gp   also used if sulfide is produced
         sd_jc_o2equiv = sd_jc - sd_jo2no3t
         if(sd_jc_o2equiv<0)sd_jc_o2equiv = 0.0D+00
!		!
!		!***** Methane/sulfide in O2 equivalents
!		!:::::::::::::::::::::::::::: BEGIN methane/sulfide::::::::::::::::::::::::::::
!		!
!		!gp   select for freshwater/saltwater methane/sulfide production from C
!        diagenesis
         if(sd_so4<=sd_sulfate_ch4_h2s)then
!			!gp      freshwater methane production, no changes to original code
!			!CSODMAX Equations 10.28 and 10.30
!			!SD_CSODmax = DMin1((2.0D+00 * SD_KL12 * SD_CH4SAT * SD_JC_O2equiv)**2.0D+00, SD_JC_O2equiv)    ![gmO*/m2-d] = sqr([m/d] * [gmO*/m3] * [gmO*/m2-d])   ! SW 10/10/2017 MAJOR ERROR
             sd_csodmax = DMIN1((2.0D+00*sd_kl12*sd_ch4sat*sd_jc_o2equiv)      &
                        & **0.5D+00, sd_jc_o2equiv)                                                            ![gmO*/m2-d] = sqr([m/d] * [gmO*/m3] * [gmO*/m2-d])   ! SW 10/10/2017
             if(sd_ch4compmethod==0)then
!				!***********************************************************************
!				!Analytical solution for methane
                 sd_sech_arg = (sd_ae_ch4_co2*sd_theta_ch4_co2**((SD_T(1) - 20.&
                             & )/2.0))/sd_s
!				!CSOD       Equation 10.35
!				!The        hyperbolic secant is defined as HSec(X) = 2 / (Exp(X) +
!                Exp(-X))
                 if(sd_sech_arg<400.0)then                            !This is the usual case
                     sd_csod = sd_csodmax*                                     &
                             & (1.0 - (2.0/(EXP(sd_sech_arg) + EXP(-sd_sech_arg&
                             & ))))
                 else                        !HSec(SECH_ARG) < 3.8E-174 ~ 0
                     sd_csod = sd_csodmax
                 endif
!				!***********************************************************************
             else
!				!
!				!NumericalSolution for methane
                 sd_ch4toco2 = (sd_ae_ch4_co2**2.0D+00*sd_theta_ch4_co2**      &
                             & ((SD_T(1) - 20.0D+00)/2.0D+00))/sd_s
                 SD_CH4(1) = (sd_csodmax + sd_s*sd_ch40)/(sd_ch4toco2 + sd_s)
                 sd_csod = sd_ch4toco2*SD_CH4(1)
             endif
!
!			!0.5CH4  + O2 --> 0.5 CO2 + H2O
!			!CO2     Produced = 0.5*(12+32)/32 = 0.6875
             co2producedsrc1l1 = sd_csod*0.6875                       !g CO2/m²/d
             co2producedcon1l1 = co2producedsrc1l1*sd_tc/sd_h1                      !g CO2/m²/d*d/m = g CO2/m³
!
         else
!			!
!			!gp      saltwater sulfide production by C diagenesis based on DiToro
!			!*****   (2001) Appendix B Calculate dissolved and particulate (sorbed)
!            fractions for sulfide
             sd_fd1 = (1.0D+00/(1.0D+00 + sd_m1*sd_kdh2s1))
             sd_fp1 = 1.0D+00 - sd_fd1                                                       != ((m1*KdH2S1)/(1 + m1*KdH2S1))
             sd_fd2 = (1.0D+00/(1.0D+00 + sd_m2*sd_kdh2s2))
             sd_fp2 = 1.0D+00 - sd_fd2                                                       != ((m2*KdH2S2)/(1 + m2*KdH2S2))
!			!
!			!*****   Temperature adjusted reaction velocities
             sd_xappd1 = sd_ae_h2s_so4*sd_theta_h2s_so4**((SD_T(1) - 20.)      &
                       & /2.0D+00)
             sd_xappp1 = sd_kappah2sp1*sd_theta_h2s_so4**((SD_T(1) - 20.)      &
                       & /2.0D+00)
!			!
!			!*****   Transport and Decay terms
!			!Equation B.19
             sd_k1h1d = sd_xappd1**2.0D+00/sd_s*(sd_o20/sd_normconst_h2s_so4)  &
                      & + sd_s
             sd_k1h1p = sd_xappp1**2.0D+00/sd_s*(sd_o20/sd_normconst_h2s_so4)
             sd_k2h2d = 0.0D+00
             sd_k2h2p = 0.0D+00
             sd_f12 = sd_w12*sd_fp1 + sd_kl12*sd_fd1
             sd_f21 = sd_w12*sd_fp2 + sd_kl12*sd_fd2
             sd_xk1 = sd_k1h1d*sd_fd1 + sd_k1h1p*sd_fp1
             sd_xk2 = sd_k2h2d*sd_fd2 + sd_k2h2p*sd_fp2
!			!
!			!*****   Matrix and forcing function
             sd_a11 = -sd_f12 - sd_xk1 - sd_w2                                                !note: -fd1 * s is included in DiToro's -xk1 term
                !SD_a11 = -SD_H1 / SD_tc - SD_F12 - SD_xk1 - SD_w2                             !note: -fd1 * s is included in DiToro's -xk1 term
             sd_a21 = sd_f12 + sd_w2
             sd_a12 = sd_f21
             sd_b1 = 0.0D+00
!
             sd_a22 = -sd_f21 - sd_xk2 - sd_w2 - sd_h2/sd_tc
             sd_b2 = -sd_jc_o2equiv - sd_h2/sd_tc*SD_HSTP(2)
!
             call LIN_SYS(sd_a11, sd_a12, sd_a21, sd_a22, sd_b1, sd_b2,        &
                        & SD_HST(1), SD_HST(2), nflog, nfcle)
!			!
!			!*****   dissolved concentrations
             SD_HS(1) = sd_fd1*SD_HST(1)
             SD_HS(2) = sd_fd2*SD_HST(2)
             sd_csod = (sd_xappd1**2.0D+00/sd_s*sd_fd1 +                       &
                     & sd_xappp1**2.0D+00/sd_s*sd_fp1)                         &
                     & *(sd_o20/sd_normconst_h2s_so4)*SD_HST(1)
!
!			!H2S     + 2 O2 --> 2 H+ + SO42-
!			!SO42-   Produced = (32+16*4)/(2*32) = 1.5
             so4producedsrc1l1 = sd_csod*1.5                       !g SO42-/m²/d
!
!			!CH2O    + 2 H+ + SO42- --> 2 CO2 + H2S + 2 H2O    !CH2O is represented in
!			!SO42-   O2 equivalent in SD_JC_O2equiv Consumed = (32+16*4)/(2*16) = 3.0
             so4consumedsnk1l2 = sd_jc_o2equiv*3.0                       !g SO42-/m²/d
!
!			!CH2O    + 2 H+ + SO42- --> 2 CO2 + H2S + 2 H2O    !CH2O is represented in
!			!CO2     O2 equivalent in SD_JC_O2equiv Produced = 2*(12+16*2)/(2*16) =
!            2.75
             co2producedsrc2l2 = sd_jc_o2equiv*2.75                       !g CO2/m²/d
             co2producedcon2l2 = co2producedsrc2l2*sd_tc/sd_h2                      !g CO2/m²/d*d/m = g CO2/m³
         endif
 
            ! Metals Start
         if(includeiron)then
            !calculating dissolved and particulate forms of ferrous iron Fe(II) (Chapra, eqn. 25.89)
             sd_fd1 = 1.0/(sd_porosity + kdfe1*sd_rho*(1.0 - sd_porosity))
             sd_fp1 = (kdfe1*sd_rho*(1.0 - sd_porosity))                       &
                    & /(sd_porosity + kdfe1*sd_rho*(1.0 - sd_porosity))
             sd_fd2 = 1.0/(sd_porosity + kdfe2*sd_rho*(1.0 - sd_porosity))
             sd_fp2 = (kdfe2*sd_rho*(1.0 - sd_porosity))                       &
                    & /(sd_porosity + kdfe2*sd_rho*(1.0 - sd_porosity))
 
             sd_fe2tofeooh = kfe_oxid*sd_o20*10**(2.0*(SD_PH(1) - 7.0))        &
                           & *sd_fd1*SD_FE2T(1)
             sd_csod = sd_csod + (0.25*2.0*16.0/55.845)*sd_fe2tofeooh ! lumping in DO consumed by Fe(II)>FeOOH into CSOD
             sd_feoohtofe2 = kfe_red*SD_FEOOH(1)
            !FeOOH + 0.25CH2O + 2H+ > Fe(II) + 0.25CO2 + 1.75H20
            ! CO2 Produced = 0.25 * (12 + 2*16) / 55.845 = 0.197         !g CO2/m²/d
             co2producedsrc1l2 = co2producedsrc1l2 + sd_feoohtofe2*0.197
 
            !Write linear system of equations around total ferrous iron SD_Fe2T
!		!Equation 5.1
!		!
!		!Layer    1
             sd_a11 = -sd_fd1*sd_kl12 - sd_fp1*sd_w12 - sd_fd1*sd_s - sd_w2 -  &
                    & sd_fe2tofeooh*sd_h1
             sd_a12 = sd_fd2*sd_kl12 + sd_fp2*sd_w12
             sd_b1 = -sd_s*sd_fe20
!		!
!		!Layer    2
             sd_a21 = sd_fd1*sd_kl12 + sd_fp1*sd_w12 + sd_w2
             sd_a22 = -sd_fd2*sd_kl12 - sd_fp2*sd_w12 - sd_w2 - sd_h2/sd_tc
             sd_b2 = -sd_h2/sd_tc*SD_FE2T(2) - sd_h2*sd_feoohtofe2
!
             call LIN_SYS(sd_a11, sd_a12, sd_a21, sd_a22, sd_b1, sd_b2,        &
                        & SD_FE2T(1), SD_FE2T(2), nflog, nfcle)
 
            ! dissolved Fe(2)
             SD_FE2(1) = sd_fd1*SD_FE2T(1)
             SD_FE2(2) = sd_fd2*SD_FE2T(2)
!
             sd_jfe2 = sd_s*(SD_FE2(1) - sd_fe20)
 
            !Write linear system of equations around total ferrous iron SD_FeOOH
!		!Equation 5.1
!		!
!		!Layer    1
             sd_a11 = -sd_w12 - sd_w2
             sd_a12 = sd_w12
             sd_b1 = -sd_jfeoohin - sd_fe2tofeooh*sd_h1
!		!
!		!Layer    2
             sd_a21 = sd_w12 + sd_w2
             sd_a22 = -sd_w12 - sd_w2 - sd_h2*sd_feoohtofe2 - sd_h2/sd_tc
             sd_b2 = -sd_h2/sd_tc*SD_FEOOH(2)
!
             call LIN_SYS(sd_a11, sd_a12, sd_a21, sd_a22, sd_b1, sd_b2,        &
                        & SD_FEOOH(1), SD_FEOOH(2), nflog, nfcle)
         endif
 
         if(includemanganese)then
            !calculating dissolved and particulate forms of Mn(II) (Chapra, eqn. 25.89)
             sd_fd1 = 1.0/(sd_porosity + kdmn1*sd_rho*(1.0 - sd_porosity))
             sd_fp1 = (kdmn1*sd_rho*(1.0 - sd_porosity))                       &
                    & /(sd_porosity + kdmn1*sd_rho*(1.0 - sd_porosity))
             sd_fd2 = 1.0/(sd_porosity + kdmn2*sd_rho*(1.0 - sd_porosity))
             sd_fp2 = (kdmn2*sd_rho*(1.0 - sd_porosity))                       &
                    & /(sd_porosity + kdmn2*sd_rho*(1.0 - sd_porosity))
 
             sd_mn2tomno2 = kmn_oxid*sd_o20*10**(2.0*(SD_PH(1) - 7.0))         &
                          & *sd_fd1*SD_MN2T(1)
             sd_csod = sd_csod + (16.0/54.94)*sd_mn2tomno2 ! lumping in DO consumed by Mn(II)>MnO2 into CSOD
             sd_mno2tomn2 = kmn_red*SD_MNO2(1)
            !MnO2 + 0.5CH2O + 2H+ > Mn(II) + 0.5CO2 + 1.5H20
            ! CO2 Produced = 0.5 * (12 + 2*16) / 54.94 = 0.400         !g CO2/m²/d
             co2producedsrc1l2 = co2producedsrc1l2 + sd_mno2tomn2*0.400
 
            !Write linear system of equations around total Mn(II) SD_Mn2T
!		!Equation 5.1
!		!
!		!Layer    1
             sd_a11 = -sd_fd1*sd_kl12 - sd_fp1*sd_w12 - sd_fd1*sd_s - sd_w2 -  &
                    & sd_mn2tomno2*sd_h1
             sd_a12 = sd_fd2*sd_kl12 + sd_fp2*sd_w12
             sd_b1 = -sd_s*sd_mn20
!		!
!		!Layer    2
             sd_a21 = sd_fd1*sd_kl12 + sd_fp1*sd_w12 + sd_w2
             sd_a22 = -sd_fd2*sd_kl12 - sd_fp2*sd_w12 - sd_w2 - sd_h2/sd_tc
             sd_b2 = -sd_h2/sd_tc*SD_MN2T(2) - sd_h2*sd_mno2tomn2
!
             call LIN_SYS(sd_a11, sd_a12, sd_a21, sd_a22, sd_b1, sd_b2,        &
                        & SD_MN2T(1), SD_MN2T(2), nflog, nfcle)
 
            ! dissolved Mn(II)
             SD_MN2(1) = sd_fd1*SD_MN2T(1)
             SD_MN2(2) = sd_fd2*SD_MN2T(2)
!
             sd_jmn2 = sd_s*(SD_MN2(1) - sd_mn20)
 
            !Write linear system of equations around manganese dioxide SD_MnO2
!		!Equation 5.1
!		!
!		!Layer    1
             sd_a11 = -sd_w12 - sd_w2
             sd_a12 = sd_w12
             sd_b1 = -sd_jmno2in - sd_mn2tomno2*sd_h1
!		!
!		!Layer    2
             sd_a21 = sd_w12 + sd_w2
             sd_a22 = -sd_w12 - sd_w2 - sd_h2*sd_mno2tomn2 - sd_h2/sd_tc
             sd_b2 = -sd_h2/sd_tc*SD_MNO2(2)
!
             call LIN_SYS(sd_a11, sd_a12, sd_a21, sd_a22, sd_b1, sd_b2,        &
                        & SD_MNO2(1), SD_MNO2(2), nflog, nfcle)
         endif
            ! Metals End
 
         sd_sod = (sd_sodold + sd_csod + sd_nsod)/2
         iter = iter + 1
         sd_ea = ABS((sd_sod - sd_sodold)/sd_sod)*100.0D+00
!		!SOD  = (SODold + CSOD + NSOD) / 2
!		!iter = iter + 1
!		!ea   = Abs((SOD - SODold) / SOD) * 100
         if(sd_ea<=sd_es)exit
         if(iter>=maxit)then
!			!MsgBox  "SOD iterations exceeded"
!			!Write(NFlog,*) 'SOD iterations exceeded'
             write(*, *)'SOD iterations exceeded'
             exit
         endif
         sd_sodold = sd_sod
     enddo
     sd_s = sd_sod/sd_o20
!	!
!	!Flux to water column
     sd_jsod = sd_sod
     sd_jnh4 = sd_s*(SD_NH3(1) - sd_nh30)
     sd_jno3 = sd_s*(SD_NO3(1) - sd_no30)
 
        !CEMA Additional Code
     if(sd_so4>sd_sulfate_ch4_h2s)then
 
            !Calculate H2S and HS- concentrations
!		!
!		!SD1_SulfiMinus = SD_HS(1)/(1 + 10**(-SD_pHValue(SegNumI))/10**(-HS_H2S_Eqb_Const))
         sd1_sulfiminus = SD_HS(1)                                             &
                        & /(1 + 10**( - SD_PH(1))/10**( - hs_h2s_eqb_const))
         sd1_sulfide = SD_HS(1) - sd1_sulfiminus
!		!SD2_SulfiMinus = SD_HS(2)/(1 + 10**(-SD_pHValue(SegNumI))/10**(-HS_H2S_Eqb_Const))
         sd2_sulfiminus = SD_HS(2)                                             &
                        & /(1 + 10**( - SD_PH(2))/10**( - hs_h2s_eqb_const))
         sd2_sulfide = SD_HS(2) - sd2_sulfiminus
!
     endif
        !End CEMA Additional Code
 
!	!
!	!gp   methane or sulfide fluxes produced from C diagenesis
     if(sd_so4<sd_sulfate_ch4_h2s)then
!		!gp   freshwater sediment fluxes - methane
         if(sd_ch4compmethod==0)then
!			!***********************************************************************
!			!Aqueous methane flux to water column
             sd_sjch4 = sd_csodmax - sd_csod
             sd_jch4 = sd_sjch4
!			!Gaseous methane flux to water column
             sd_jch4g = sd_jc_o2equiv - sd_jch4 - sd_csod                              ! flux in gO2/m2/day
!			!***********************************************************************
             SD_CH4(1) = sd_jch4/sd_s + sd_ch40
         else
!			!***********************************************************************
!			!numerical solution for methane
             sd_jch4 = sd_s*(SD_CH4(1) - sd_ch40)
             sd_jch4g = sd_jc_o2equiv - sd_jch4 - sd_csod                              ! flux in gO2/m2/day
!			!***********************************************************************
         endif
         sd_jhs = 0                    !gp
     else            !gp
!		!gp   marine sediment fluxes - sulfide
         sd_jch4 = 0
         sd_jch4g = 0
         sd_jhs = sd_s*SD_HS(1)
     endif              !gp
!
!	!CEMA Additional Code
     sd_jso4 = sd_s*sd_so4
!	!End CEMA Additional Code
!	!
!	!gp estimated thickness of the aerobic sediment layer 1 (DiToro Appendix B Page 576)  MOVED ABOVE
!	!SD_H1 = SD_KL12 * SD_H2 / SD_s
!	!If(SD_H1 > SD_H2)Then
        !		    SD_H1 = SD_H2
!	!End If
!	!SD_AerLayerThick(SegNumI) = SD_H1
!
!	!Bubbles formation
!	!Aerobic layer is thin and so ignore gas formation from aerobic layer
!	   !H2S
!	!Calculate dissolved <--> gaseous phase distribution
     sedtemp1 = SD_T(1) + 273.15
     mw_constituent = 36.            !S = 34, H = 1
!	!Calculate dissolved <--> gaseous phase distribution
!	!VolWater = CellArea(SegNumI)*SD_H2
     sedtemp2 = SD_T(2) + 273.15
!	!Layer 1
     volwater = cellarea(segnumi)*sd_h1
     call CEMADISGASPHASEDISTRIBUTION(sd1_sulfide, mw_constituent, sedtemp1,   &
                                    & henryconst_h2s, volwater, sulfideg_sd1,  &
                                    & sulfided_sd1)
!	!Layer 2
     volwater = cellarea(segnumi)*sd_h2
     call CEMADISGASPHASEDISTRIBUTION(sd2_sulfide, mw_constituent, sedtemp2,   &
                                    & henryconst_h2s, volwater, sulfideg_sd2,  &
                                    & sulfided_sd2)
!	!Old TConc(1,SegNumI) = SD2_Sulfide
!	!Old SConc(1,SegNumI) = (TConc(1,SegNumI) - TConcP(1,SegNumI))/dlt    !gm/m³/s
     TCONC(1, segnumi) = sulfideg_sd2*bubbaccfraction + TCONCP(1, segnumi)
     SCONC(1, segnumi) = (sulfideg_sd2*bubbaccfraction)/dlt              !gm/m³/s
!	!Old TConcP(1,SegNumI) = SD2_Sulfide
     TCONCP(1, segnumi) = TCONC(1, segnumi)
!
!	!CH4
     TCONC(2, segnumi) = SD_CH4(1)*bubbaccfraction + TCONCP(2, segnumi)
     SCONC(2, segnumi) = (SD_CH4(1)*bubbaccfraction)/dlt              !gm/m³/s
!	!Old TConcP(1,SegNumI) = SD2_Sulfide
     TCONCP(2, segnumi) = TCONC(2, segnumi)
!
!	!NH3
!	!Old TConc(3,SegNumI) = SD2_Ammonia
!	!Old SConc(3,SegNumI) = (TConc(3,SegNumI) - TConcP(3,SegNumI))/dlt    !gm/m³/s
     TCONC(3, segnumi) = SD_NH3(2)           !SD2_Ammonia*BubbAccFraction + TConcP(3,SegNumI)
     SCONC(3, segnumi) = (SD_NH3(2) - TCONCP(3, segnumi))/dlt           !gm/m³/s
!	!Old TConcP(3,SegNumI) = SD2_Ammonia
     TCONCP(3, segnumi) = TCONC(3, segnumi)
!
!	!CO2
!	!Old TConc(4,SegNumI) = CO2ProducedCon2L2
!	!Old SConc(4,SegNumI) = (TConc(4,SegNumI) - TConcP(4,SegNumI))/dlt    !gm/m³/s
     TCONC(4, segnumi) = co2producedcon2l2*bubbaccfraction + TCONCP(4, segnumi)
     SCONC(4, segnumi) = (co2producedcon2l2*bubbaccfraction)/dlt              !gm/m³/s
!	!Old TConcP(4,SegNumI) = CO2ProducedCon2L2
     TCONCP(4, segnumi) = TCONC(4, segnumi)
!
!	!H2S
     DISSOLVEDGASSEDIMENTS(1, segnumi) = sulfided_sd2
!    !CH4
     DISSOLVEDGASSEDIMENTS(2, segnumi) = SD_CH4(1)
!	!NH4
     DISSOLVEDGASSEDIMENTS(3, segnumi) = ammoniad_sd2
!	!CO2
     DISSOLVEDGASSEDIMENTS(4, segnumi) = 0.D0
!
     if(bubbles_calculation)call GASBUBBLESFORMATION(BUBBLERADIUSSED(segnumi), &
      & dlt, volwater)
!	!Update Sediment Bed Concentrations - additional mass balance (due to PW release)
!	!Update concentrations of SO4, NH3, NO3, H2S
!	!Sulfate
     sd_so4 = sd_so4 + (so4producedsrc1l1 - so4consumedsnk1l2)                 &
            & *sd_tc/(sd_h1 + sd_h2)
     if(sd_so4<0.0)sd_so4 = 0.0
!	!Remaining source/sink (diffusion to water column and PW release included later)
!
        !Write linear system of equations around Total inorganic carbon
!		!Equation 5.1
!		!
!		!Layer 1
     sd_a11 = -sd_kl12 - sd_s
     sd_a12 = sd_kl12
          ! COCO2ProducedSrc1L1 converted from  g CO2/m²/d to g C/m²/d ; 12/(2*16+12)=0.272
     sd_b1 = -sd_s*sd_tic0 - co2producedsrc1l1*0.272
!		!
!		!Layer 2
     sd_a21 = sd_kl12
     sd_a22 = -sd_kl12 - sd_h2/sd_tc
          ! COCO2ProducedSrc1L2 converted from  g CO2/m²/d to g C/m²/d ; 12/(2*16+12)=0.272
     sd_b2 = -sd_h2/sd_tc*SD_TIC(2) - co2producedsrc2l2*0.272
!
     call LIN_SYS(sd_a11, sd_a12, sd_a21, sd_a22, sd_b1, sd_b2, SD_TIC(1),     &
                & SD_TIC(2), nflog, nfcle)
 
            ! calculating diffusive flux between layer 1 and water column
     sd_jtic = sd_s*(SD_TIC(1) - sd_tic0)
 
     if(includealkalinity)then
        !Write linear system of equations around Total Alkalinity
        ! Nitrification of ammonium results in an alkalinity decrease: 2 eq. alk per 1 mole ammonium
        ! Denitrification of nitrate (to nitrogen gas) results in an alkalinity increase: 1 eq. alk per 1 mole nitrate
!		!Equation 5.1
!		!
!		!Layer 1
         sd_a11 = -sd_kl12 - sd_s
            !SD_a11 = -SD_H1/SD_tc -SD_KL12 -SD_s
         sd_a12 = sd_kl12
         sd_b1 = -sd_s*sd_alk0 + 2.0*sd_nh3tono3*SD_NH3(1) - SD_NO3(1)         &
               & *SD_DENIT(1)                                                                                      ![m/d]*[mg/m3]
!		!
!		!Layer 2
         sd_a21 = sd_kl12
         sd_a22 = -sd_kl12 - sd_h2/sd_tc
         sd_b2 = -sd_h2/sd_tc*SD_ALK(2) - SD_NO3(2)*SD_DENIT(2)
!
         call LIN_SYS(sd_a11, sd_a12, sd_a21, sd_a22, sd_b1, sd_b2, SD_ALK(1), &
                    & SD_ALK(2), nflog, nfcle)
 
            ! calculating diffusive flux between layer 1 and water column
         sd_jalk = sd_s*(SD_ALK(1) - sd_alk0)
 
     endif
 
            !calculating dissolved and particulate forms of phosphorus (Chapra, eqn. 25.89)
     sd_fd1 = 1.0/(sd_porosity + kdp1*sd_rho*(1.0 - sd_porosity))
     sd_fp1 = (kdp1*sd_rho*(1.0 - sd_porosity))                                &
            & /(sd_porosity + kdp1*sd_rho*(1.0 - sd_porosity))
     sd_fd2 = 1.0/(sd_porosity + kdp2*sd_rho*(1.0 - sd_porosity))
     sd_fp2 = (kdp2*sd_rho*(1.0 - sd_porosity))                                &
            & /(sd_porosity + kdp2*sd_rho*(1.0 - sd_porosity))
 
            !Write linear system of equations around total phosphate SD_PO4T
!		!Equation 5.1
!		!
!		!Layer 1
     sd_a11 = -sd_fd1*sd_kl12 - sd_fp1*sd_w12 - sd_fd1*sd_s - sd_w2
            !SD_a11 = -SD_H1/SD_tc - SD_fd1*SD_KL12 - SD_fp1*SD_w12 - SD_fd1*SD_s - SD_w2
     sd_a12 = sd_fd2*sd_kl12 + sd_fp2*sd_w12
     sd_b1 = -sd_s*sd_po40                                      ![m/d]*[mg/m3]
!		!
!		!Layer 2
     sd_a21 = sd_fd1*sd_kl12 + sd_fp1*sd_w12 + sd_w2
     sd_a22 = -sd_fd2*sd_kl12 - sd_fp2*sd_w12 - sd_w2 - sd_h2/sd_tc
     sd_b2 = -sd_jp - sd_h2/sd_tc*SD_PO4T(2)
!
     call LIN_SYS(sd_a11, sd_a12, sd_a21, sd_a22, sd_b1, sd_b2, SD_PO4T(1),    &
                & SD_PO4T(2), nflog, nfcle)
 
            ! dissolved PO4
     SD_PO4(1) = sd_fd1*SD_PO4T(1)
     SD_PO4(2) = sd_fd2*SD_PO4T(2)
!
     sd_jp = sd_s*(SD_PO4(1) - sd_po40)
 
 
        !Write linear system of equations for temperature
     sd_rhowcp = 4.186*1.0E6        ! units J g-1 C-1 * g m-3= J C-1 m-3
!		!Equation 5.1
!		!
!		!Layer 1
     sd_a11 = -sd_kl12 - sd_s
     sd_a12 = sd_kl12
     sd_b1 = -sd_s*sd_tw
!		!
!		!Layer 2
     sd_a21 = sd_kl12*sd_rhowcp
     sd_a22 = -sd_kl12*sd_rhowcp - sd_h2/sd_tc*sd_rhowcp + sd_ksw
     sd_b2 = -sd_h2*sd_rhowcp/sd_tc*SD_T(2) - sd_ksw*sd_tsed
!
     call LIN_SYS(sd_a11, sd_a12, sd_a21, sd_a22, sd_b1, sd_b2, SD_T(1),       &
                & SD_T(2), nflog, nfcle)
 
            ! calculating diffusive heat flux between layer 1 and water column
     sd_jt = sd_s*sd_rhowcp*(SD_T(1) - sd_tw)
 
     if(includedynamicph)then
        ! Sediment pH for layers 1 and 2
         sd_tds = 0.0
                    !sediments not simulating tds at the moment
         sd_poct1 = 0.0
                     ! POC not predicted for aerobic layer.
         call PH_SEDIMENTS(SD_T(1), SD_TIC(1), SD_ALK(1), SD_NH3(1), SD_PO4(1),&
                         & sd_poct1, sd_tds(1), SD_PH(1))                                               !layer 1
 
         call PH_SEDIMENTS(SD_T(2), SD_TIC(2), SD_ALK(2), SD_NH3(2), SD_PO4(2),&
                         & sd_poct2, sd_tds(2), SD_PH(2))                                               !layer 2
     endif
 
!	!Transfer to variables that will be used to hook up the model to W2 water quality model
 
        ! SW compute CO2 in layer 1 and layer 2 10/20/2017
            !HION = 10.0**SD_PH(1)
            !SD_CO2(1) = SD_TIC(1)/(1.0+K1/HION+K1*K2/(HION*HION))    ! K1 AND K2 FROM PH_SEDIMENTS SUBROUTINE
            !
            !
            !HION = 10.0**SD_PH(2)
            !SD_CO2(2) = SD_TIC(2)/(1.0+K1/HION+K1*K2/(HION*HION))
            !
 
!
!	!1
!	!Bubbles release
!	!Bubbles of CH4g, NH3g and H2Sg
!
!
!	!2
!	!Porewater release
!	!Flux of CH4d, NH3d + NH4d, H2Sd + HSd, SO42-d, NO3d, CO2d
!	!Volume of porewater
     volwater = cellarea(segnumi)*sd_h2
 
        !Aerobic Layer
     pw_relrate1 = POREWATERRELRATE(segnumi)*sd_h1/(sd_h1 + sd_h2)*86400.0     !m³/d
        !Anaerobic Layer
     pw_relrate2 = POREWATERRELRATE(segnumi)*sd_h2/(sd_h1 + sd_h2)*86400.0     !m³/d
!
!	!Dissolved_CH4_Src = SD_JCH4*SD_tc/(SD_H1 + SD_H2)   !CHECK CODE This is incorrect original code since SD_tc is in units of s*d/s or d ==> g/m3 Why didvide by SD_H1 and SD_H2 inncorrect
        !Dissolved_CH4_Src = SD_JCH4*SD_tc/(SD_H1 + SD_H2)   !SW 10/10/2017    gO2/m2/d ==> gO2/m3/d
!
!	!Dissolved_NH3_Src = (AmmoniaD_SD1*SD_H1 + AmmoniaD_SD2*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
!	!Dissolved_NH3_Src = Dissolved_NH3_Src + (SD1_Ammonium*SD_H1 + SD2_Ammonium*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))  !gm/m³/s
     dissolved_nh3_src = (ammoniad_sd1*pw_relrate1 + ammoniad_sd2*pw_relrate2) &
                       & /(cellthickness*cellarea(segnumi))                                                                              !gm/m³/d
     dissolved_nh3_src = dissolved_nh3_src +                                   &
                       & (sd1_ammonium*pw_relrate1 + sd2_ammonium*pw_relrate2) &
                       & /(cellthickness*cellarea(segnumi))                                                                                      !gm/m³/d
!
        !Dissolved_H2S_Src = (SulfideD_SD1*SD_H1 + SulfideD_SD2*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
!	!Dissolved_H2S_Src = Dissolved_H2S_Src + (SD1_SulfiMinus*SD_H1 + SD2_SulfiMinus*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))  !gm/m³/s
     dissolved_h2s_src = (sulfided_sd1*pw_relrate1 + sulfided_sd2*pw_relrate2) &
                       & /(cellthickness*cellarea(segnumi))                                                                              !gm/m³/d
     dissolved_h2s_src = dissolved_h2s_src +                                   &
                       & (sd1_sulfiminus*pw_relrate1 + sd2_sulfiminus*         &
                       & pw_relrate2)/(cellthickness*cellarea(segnumi))                                                                              !gm/m³/d
!
!	!Dissolved_NO3_Src = (SD_NO3(1)*SD_H1 + SD_NO3(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
     dissolved_no3_src = (SD_NO3(1)*pw_relrate1 + SD_NO3(2)*pw_relrate2)       &
                       & /(cellthickness*cellarea(segnumi))                                                                        !gm/m³/d
!
!	!Dissolved_SO4_Src = SD_SO4 * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
     dissolved_so4_src = sd_so4*(pw_relrate1 + pw_relrate2)                    &
                       & /(cellthickness*cellarea(segnumi))                                                           !gm/m³/s
 
!	!Dissolved_CO2_Src = (CO2ProducedCon1L1*SD_H1 + CO2ProducedCon2L2*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
     dissolved_co2_src = (co2producedcon1l1*pw_relrate1 +                      &
                       & co2producedcon2l2*pw_relrate2)                        &
                       & /(cellthickness*cellarea(segnumi))                                                                                        !gm/m³/d
        !Dissolved_ALK_Src = (SD_ALK(1)*SD_H1 + SD_ALK(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
     if(includealkalinity)then
         dissolved_alk_src = (SD_ALK(1)*pw_relrate1 + SD_ALK(2)*pw_relrate2)   &
                           & /(cellthickness*cellarea(segnumi))                                                                    !gm/m³/d
     else
         dissolved_alk_src = 0.0
     endif
        !Dissolved_PO4_Src = (SD_PO4(1)*SD_H1 + SD_PO4(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
     dissolved_po4_src = (SD_PO4(1)*pw_relrate1 + SD_PO4(2)*pw_relrate2)       &
                       & /(cellthickness*cellarea(segnumi))                                                                        !gm/m³/d
        !Dissolved_Fe2_Src = (SD_Fe2(1)*SD_H1 + SD_Fe2(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
     if(includeiron)then
         dissolved_fe2_src = (SD_FE2(1)*pw_relrate1 + SD_FE2(2)*pw_relrate2)   &
                           & /(cellthickness*cellarea(segnumi))                                                                    !gm/m³/d
     else
         dissolved_fe2_src = 0.0
     endif
        !Dissolved_Mn2_Src = (SD_Mn2(1)*SD_H1 + SD_Mn2(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
     if(includemanganese)then
         dissolved_mn2_src = (SD_MN2(1)*pw_relrate1 + SD_MN2(2)*pw_relrate2)   &
                           & /(cellthickness*cellarea(segnumi))                                                                    !gm/m³/d
     else
         dissolved_mn2_src = 0.0
     endif
        !Sediment_Heat_Src = SD_rhowcp*(SD_T(1)*SD_H1 + SD_T(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !J/m³/s
     sediment_heat_src = sd_rhowcp*(SD_T(1)*pw_relrate1 + SD_T(2)*pw_relrate2) &
                       & /(cellthickness*cellarea(segnumi))                                                                              !J/m³/d
!
        !3
        !Diffusive flux of NH3, NO3, CH4, SO4, H2S, CO2
     dissolved_nh3_src = dissolved_nh3_src + sd_jnh4/cellthickness      !g/m³/d SD_JNH4/CellThickness = g/m²/d/m = g/m³/d
     dissolved_no3_src = dissolved_no3_src + sd_jno3/cellthickness      !g/m³/d SD_JNO3/CellThickness = g/m²/d/m = g/m³/d
        !Dissolved_CH4_Src = Dissolved_CH4_Src + SD_JCH4/CellThickness   !g/m³/d SD_JCH4/CellThickness = g/m²/d/m = g/m³/d commented out because double counting...
     dissolved_ch4_src = dissolved_ch4_src + sd_jch4/cellthickness      !g/m³/d SD_JCH4/CellThickness = g/m²/d/m = g/m³/d commented out because double counting... ! SW 10/10/2017 added back because eliminated the problem above
     dissolved_so4_src = dissolved_so4_src + sd_jso4/cellthickness      !g/m³/d SD_JSO4/CellThickness = g/m²/d/m = g/m³/d
     dissolved_h2s_src = dissolved_h2s_src + sd_jhs/cellthickness
     dissolved_co2_src = dissolved_co2_src + sd_jtic/cellthickness
     dissolved_alk_src = dissolved_alk_src + sd_jalk/cellthickness
     dissolved_po4_src = dissolved_po4_src + sd_jp/cellthickness
     dissolved_fe2_src = dissolved_fe2_src + sd_jfe2/cellthickness
     dissolved_mn2_src = dissolved_mn2_src + sd_jmn2/cellthickness
 
     gasreleasech4 = gasreleasech4 + sd_jch4g*B(KB(segnumi), segnumi)          &
                   & *DLX(segnumi)*sd_tc/2.67                                                                  ! Flux of gas CH4 in gC   SW 10/10/2017
        !GasReleaseCO2=GasReleaseCO2+
 
     sediment_heat_src = sediment_heat_src + sd_jt/cellthickness     ! Heat J/m3/d  SD_JT/CellThickness = J/m3/d/m= J/m3/d
 
!	!Flux to CSOD and NSOD
     dissolved_o2_snk = sd_jsod/cellthickness              !g/m³/d SD_JSOD/CellThickness = g/m²/d/m = g/m³/d
        !BOD sediment feature
     if(includecemagenbodconstituents)call COMPUTEGENBODFATE
!
!	!Convert all source/sink to g/m³/s from g/m³/d
!	!NH3, NO3, CH4, SO4, DO, CO2, H2S
     dissolved_nh3_src = dissolved_nh3_src/86400.0                  !g/m³/d --> g/m³/s
     dissolved_no3_src = dissolved_no3_src/86400.0          !g/m³/d --> g/m³/s
     dissolved_ch4_src = dissolved_ch4_src/86400.0          !g/m³/d --> g/m³/s
     dissolved_so4_src = dissolved_so4_src/86400.0          !g/m³/d --> g/m³/s
     dissolved_co2_src = dissolved_co2_src/86400.0                  !g/m³/d --> g/m³/s
     dissolved_h2s_src = dissolved_h2s_src/86400.0                  !g/m³/d --> g/m³/s
     dissolved_o2_snk = dissolved_o2_snk/86400.0                    !g/m³/d --> g/m³/s
     dissolved_alk_src = dissolved_alk_src/86400.0           !g/m³/d --> g/m³/s
     dissolved_po4_src = dissolved_po4_src/86400.0           !g/m³/d --> g/m³/s
     dissolved_fe2_src = dissolved_fe2_src/86400.0           !g/m³/d --> g/m³/s
     dissolved_mn2_src = dissolved_mn2_src/86400.0           !g/m³/d --> g/m³/s
     if(includecemagenbodconstituents)then
         do genbodnum = 1, numgenbodconstituents
             SEDGENBODSS(genbodnum) = SEDGENBODSS(genbodnum)/cellthickness          !g/m³/d SD_JSO4/CellThickness = g/m²/d/m = g/m³/d
             SEDGENBODSS(genbodnum) = SEDGENBODSS(genbodnum)/86400.0                !g/m³/d --> g/m³/s
         enddo   !GenBODNum
     endif
     sediment_heat_src = sediment_heat_src/86400.0          !J/m³/d --> J/m³/s
 
       ! resuspension of POM, POP, and PON
     do itemp = 1, 2    ! only labile and refractory, not including inert for now
         lpom_resuspension = lpom_resuspension + sd_epoc(itemp)                &
                           & /cellthickness/ORGC(jw)
         rpom_resuspension = rpom_resuspension + sd_epoc(itemp)                &
                           & /cellthickness/ORGC(jw)
         lpomn_resuspension = lpomn_resuspension + sd_epon(itemp)/cellthickness
         rpomn_resuspension = rpomn_resuspension + sd_epon(itemp)/cellthickness
         lpomp_resuspension = lpomp_resuspension + sd_epop(itemp)/cellthickness
         rpomp_resuspension = rpomp_resuspension + sd_epop(itemp)/cellthickness
     enddo
 
!
!	!Update Sediment Bed Concentrations - additional mass balance for SO4
!	!NH3, NO3, H2S don't need to be updated as the diffusive flux is included in the mass balance equation
!	!Sulfate
     sd_so4 = sd_so4 - sd_jso4*sd_tc/(sd_h1 + sd_h2)            !Only SO4 has this term as NH3, NO3 and H2S equations include diffusion
     if(sd_so4<0.0)sd_so4 = 0.0
!	!PW release loss only reduces total mass but not the concentration. Lost mass is lost along with water so the concentration remains the same
!
     do itemp = 1, 3
         SD_POC22(itemp) = SD_POC2(itemp)
         SD_PON22(itemp) = SD_PON2(itemp)
         SD_POP22(itemp) = SD_POP2(itemp)
     enddo
 
     do itemp = 1, 2
         SD_NH3TP2(itemp) = SD_NH3T(itemp)
         SD_NO3P2(itemp) = SD_NO3(itemp)
!		!SD_PO4Tp2(iTemp) = SD_PO4T(iTemp)
         SD_NH3P2(itemp) = SD_NH3(itemp)
         SD_PO4P2(itemp) = SD_PO4(itemp)
         SD_CH4P2(itemp) = SD_CH4(itemp)                     !gp
         SD_HSTP2(itemp) = SD_HST(itemp)                     !gp
         SD_HSP2(itemp) = SD_HS(itemp)                        !gp
     enddo
     sd_ben_strp2 = sd_ben_str
!
     return
!
     entry COMPUTEGENBODFATE
!
 
    ! BOD debug
    !if(BODtestout)then
    !  open(8765,file='Sed_bod_debug.opt',status='unknown')
    !  write(8765,'("        JDAY   BOD_S2_L1   BOD_S3_L1   BOD_S4_L1   BOD_S2_L2   BOD_S3_L2   BOD_S4_L2")')
    !  BODtestout=.false.
    !end if
 
 
     sd_bodd = 0.D0
     dissolved_bod_src = 0.D0
     SEDGENBODSS = 0.D00
 
    !Decay and oxygen consumption
     if(sd_o20>sd_ox_threshol)then  !Only if aerobic conditions exist
         do genbodnum = 1, numgenbodconstituents
             do itemp = 1, 1    !Only Aerobic layer consumption and decay
 
                 SEDGENBODDECAYRATE(genbodnum, segnumi)                        &
                   & = SEDGENBODCONSRATE(genbodnum, segnumi)                   &
                   & *SEDGENBODCONSTCOEFF(genbodnum, segnumi)**(sd_tw - 20.)
                 sedgenbodconc(genbodnum, segnumi, itemp)                      &
                   & = (1 - SEDGENBODDECAYRATE(genbodnum, segnumi)*sd_tc)      &
                   & *sedgenbodconc(genbodnum, segnumi, itemp)
                 dissolved_o2_snk = dissolved_o2_snk +                         &
                                  & SEDGENBODDECAYRATE(genbodnum, segnumi)     &
                                  & *sedgenbodconc(genbodnum, segnumi, itemp)
 
             enddo     !iTemp
         enddo !GenBODNum
     endif
 
    !Limit the concentrations to minimum of zero
     do genbodnum = 1, numgenbodconstituents
         do itemp = 1, 1
                        !Only Aerobic layer consumption and decay
             if(sedgenbodconc(genbodnum, segnumi, itemp)<0.D00)                &
              & sedgenbodconc(genbodnum, segnumi, itemp) = 0.D00
         enddo !iTemp
     enddo !GenBODNum
 
    !Diffusion
     do genbodnum = 1, numgenbodconstituents
 
        ! Diffusion between layers
         sd_bod = sd_kl12*(sedgenbodconc(genbodnum, segnumi, 2)                &
                & - sedgenbodconc(genbodnum, segnumi, 1))
         sd_bodd(genbodnum, 1) = 1.D0*sd_bod
                                            !If SD_BOD > 0 then Layer 1 receives    gm/m²/d
         sd_bodd(genbodnum, 2) = -1.D0*sd_bod
                                            !If SD_BOD < 0 then Layer 2 receives    gm/m²/d
 
        !Diffusion to overlying water column
         sd_bod = 0.D0
        !Only if there are aerobic conditions
        !If(SD_O20 > 1.0e-10)SD_BOD = SD_s * (SedGenBODConc(GenBODNum,SegNumI,1) - SedGenBODConc0(GenBODNum))
         if(sd_o20>sd_ox_threshol)                                             &
          & sd_bod = sd_s*(sedgenbodconc(genbodnum, segnumi, 1)                &
          & - sedgenbodconc0(genbodnum))
         sd_bodd(genbodnum, 1) = sd_bodd(genbodnum, 1) - 1.D0*sd_bod
                                                                   !If SD_BOD < 0 then Layer 1 receives    gm/m²/d
         SEDGENBODSS(genbodnum) = 1.D0*sd_bod
 
     enddo !GenBODNum
 
    !Porewater release
     do genbodnum = 1, numgenbodconstituents
 
        !Aerobic Layer
         pw_relrate = POREWATERRELRATE(segnumi)*sd_h1/(sd_h1 + sd_h2)*86400. !m³/d
         dissolved_bod_src(genbodnum, 1)                                       &
           & = -1.D0*sedgenbodconc(genbodnum, segnumi, 1)                      &
           & *pw_relrate/cellarea(segnumi)                                                                         !gm/m²/d
 
        !Anaerobic Layer
         pw_relrate = POREWATERRELRATE(segnumi)*sd_h2/(sd_h1 + sd_h2)*86400. !m³/d
         dissolved_bod_src(genbodnum, 2)                                       &
           & = -1.D0*sedgenbodconc(genbodnum, segnumi, 2)                      &
           & *pw_relrate/cellarea(segnumi)                                                                         !gm/m²/d
 
         SEDGENBODSS(genbodnum) = SEDGENBODSS(genbodnum)                       &
                                & - dissolved_bod_src(genbodnum, 1)            &
                                & - dissolved_bod_src(genbodnum, 2)
 
     enddo !GenBODNum
 
    !Resulting BOD concentrations in each layer
     do genbodnum = 1, numgenbodconstituents
        !Aerobic Layer
         sedgenbodconc(genbodnum, segnumi, 1)                                  &
           & = sedgenbodconc(genbodnum, segnumi, 1)                            &
           & + (sd_bodd(genbodnum, 1) + dissolved_bod_src(genbodnum, 1))       &
           & /sd_h1*sd_tc
        !Anaerobic Layer
         sedgenbodconc(genbodnum, segnumi, 2)                                  &
           & = sedgenbodconc(genbodnum, segnumi, 2)                            &
           & + (sd_bodd(genbodnum, 2) + dissolved_bod_src(genbodnum, 2))       &
           & /sd_h2*sd_tc
     enddo !GenBODNum
 
    !Limit the concentrations to minimum of zero
     do genbodnum = 1, numgenbodconstituents
         do itemp = 1, 1
                        !Only Aerobic layer consumption and decay
             if(sedgenbodconc(genbodnum, segnumi, itemp)<0.D00)                &
              & sedgenbodconc(genbodnum, segnumi, itemp) = 0.D00
         enddo !iTemp
     enddo !GenBODNum
 
    !if(segnumi == 4)then
    !  write(8765,'(f12.6,100f10.6)')jday,(SedGenBODConc(1,ii,1), ii=2,4),(SedGenBODConc(1,ii,2), ii=2,4)
    !end if
 
     return
!
!
     entry ENDSEDIMENTFLUXVARIABLES
!
     mftsedflxvars(segnumi, 1) = 0.D00
     mftsedflxvars(segnumi, 2) = 0.D00
     mftsedflxvars(segnumi, 3) = 0.D00
     mftsedflxvars(segnumi, 4) = 0.D00
     mftsedflxvars(segnumi, 5) = 0.D00
     mftsedflxvars(segnumi, 6) = 0.D00
     mftsedflxvars(segnumi, 7) = 0.D00
     mftsedflxvars(segnumi, 8) = 0.D00
     mftsedflxvars(segnumi, 9) = 0.D00
     mftsedflxvars(segnumi, 10) = 0.D00
     mftsedflxvars(segnumi, 11) = 0.D00
     mftsedflxvars(segnumi, 12) = 0.D00
     mftsedflxvars(segnumi, 13) = 0.D00
     mftsedflxvars(segnumi, 14) = 0.D00
     mftsedflxvars(segnumi, 15) = 0.D00
 
     mftsedflxvars(segnumi, 16) = 0.D00
     mftsedflxvars(segnumi, 17) = 0.D00
     mftsedflxvars(segnumi, 18) = 0.D00
     mftsedflxvars(segnumi, 19) = 0.D00
     mftsedflxvars(segnumi, 20) = 0.D00
     mftsedflxvars(segnumi, 21) = 0.D00
     mftsedflxvars(segnumi, 22) = 0.D00
     mftsedflxvars(segnumi, 23) = 0.D00
     mftsedflxvars(segnumi, 24) = 0.D00
     mftsedflxvars(segnumi, 25) = 0.D00
 
     mftsedflxvars(segnumi, 26) = 0.D00
     mftsedflxvars(segnumi, 27) = 0.D00
     mftsedflxvars(segnumi, 28) = 0.D00
 
     mftsedflxvars(segnumi, 29) = 0.D00
     mftsedflxvars(segnumi, 30) = 0.D00
 
     mftsedflxvars(segnumi, 31) = 0.D00
     mftsedflxvars(segnumi, 32) = 0.D00
     mftsedflxvars(segnumi, 33) = 0.D00
     mftsedflxvars(segnumi, 34) = 0.D00
     mftsedflxvars(segnumi, 35) = 0.D00
     mftsedflxvars(segnumi, 36) = 0.D00
     mftsedflxvars(segnumi, 37) = 0.D00
     mftsedflxvars(segnumi, 38) = 0.D00
     mftsedflxvars(segnumi, 39) = 0.D00
     mftsedflxvars(segnumi, 40) = 0.D00
     mftsedflxvars(segnumi, 41) = 0.D00
     mftsedflxvars(segnumi, 42) = 0.D00
     mftsedflxvars(segnumi, 43) = 0.D00
     mftsedflxvars(segnumi, 44) = 0.D00
     mftsedflxvars(segnumi, 45) = 0.D00
     mftsedflxvars(segnumi, 46) = 0.D00
     mftsedflxvars(segnumi, 47) = 0.D00
 
!
     return
!
!!Entry ComputeCEMADiagenesisSourceSinks  SW 6/27/2017
 !
 !       Do JW=1, NWB
 !           KT = KTWB(JW)
 !           Do JB=BS(JW),BE(JW)
 !               IU = CUS(JB)
 !               ID = DS(JB)
 !               Do SegNumI = IU, ID
 !
 !
 !               End Do !SegNumI
 !           End Do !JB
 !       End Do !JW
 !  Return
 
     entry PH_SEDIMENTS(T1sed, Ticsed, Alksed, Nh4sed, Po4sed, Pocsed, Tdssed, &
                      & Phsed)                                                ! Enhancements added for buffering by ammonia, phosphate, and OM ! SR 01/01/12
!    pH and carbonate species
 
     t1k = T1sed + 273.15
     cart = Ticsed/12011. ! SR 01/01/12
     alkt = Alksed/50044. ! SR 01/01/12
     ammt = Nh4sed/14006.74 ! SR 01/01/12
     phost = Po4sed/30973.762 ! SR 01/01/12
     !OMCT = (LDOM(K,I)+RDOM(K,I))*ORGC(JW)/12011. ! moles carbon per liter from DOM ! SR 01/01/12
     omct = 0.0
               ! DOM is not simulated in the sediments yet...
     !IF (POM_BUFFERING) OMCT = OMCT + (LPOM(K,I)+RPOM(K,I))*ORGC(JW)/12011. ! SR 01/01/12
     if(pom_buffering)omct = omct + Pocsed/12011.   ! SR 01/01/12
     omct = 0.0  !
 !**** Ionic strength
     if(FRESH_WATER(jw))s2 = 2.5E-05*Tdssed
     if(SALT_WATER(jw))s2 = 1.47E-3 + 1.9885E-2*Tdssed + 3.8E-5*Tdssed*Tdssed
!**** Debye-Huckel terms and activity coefficients
     sqrs2 = SQRT(s2)
     dh1 = -0.5085*sqrs2/(1.0 + 1.3124*sqrs2) + 4.745694E-03 +                 &
         & 4.160762E-02*s2 - 9.284843E-03*s2*s2
     dh2 = -2.0340*sqrs2/(1.0 + 1.4765*sqrs2) + 1.205665E-02 +                 &
         & 9.715745E-02*s2 - 2.067746E-02*s2*s2
     dh3 = -4.5765*sqrs2/(1.0 + 1.3124*sqrs2)
                                            ! extended Debye-Huckel for PO4 ! SR 01/01/12
     dhh = -0.5085*sqrs2/(1.0 + 2.9529*sqrs2)
                                            ! extended Debye-Huckel for H+ ion ! SR 01/01/12
     h2co3t = 10.0**(0.0755*s2)
     hco3t = 10.0**dh1
     co3t = 10.0**dh2
     po4t = 10.0**dh3 ! SR 01/01/12
     ht = 10.0**dhh ! activity coefficient for H+ ! SR 01/01/12
     hpo4t = co3t ! tabled values similar to those for carbonate ! SR 01/01/12
     oht = hco3t ! tabled values similar to those for bicarbonate ! SR 01/01/12
     h2po4t = hco3t ! tabled values similar to those for bicarbonate ! SR 01/01/12
     nh4t = hco3t ! tabled values similar to those for bicarbonate ! SR 01/01/12
     nh3t = h2co3t ! neutral species, set coefficient to same as that for carbonic acid ! SR 01/01/12
     h3po4t = h2co3t ! neutral species, set coefficient to same as that for carbonic acid ! SR 01/01/12
!**** Temperature adjustment
     kw = 10.0**( - 283.971 - 0.05069842*t1k + 13323.0/t1k +                   &
        & 102.24447*LOG10(t1k) - 1119669.0/(t1k*t1k))/oht
     k1 = 10.0**( - 356.3094 - 0.06091964*t1k + 21834.37/t1k +                 &
        & 126.8339*LOG10(t1k) - 1684915/(t1k*t1k))*h2co3t/hco3t
     k2 = 10.0**( - 107.8871 - 0.03252849*t1k + 5151.79/t1k +                  &
        & 38.92561*LOG10(t1k) - 563713.9/(t1k*t1k))*hco3t/co3t
     kamm = 10.0**( - 0.09018 - 2729.92/t1k)*nh4t/nh3t
                                                    ! SR 01/01/12
     kp1 = 10.0**(4.5535 - 0.013486*t1k - 799.31/t1k)*h3po4t/h2po4t
                                                                  ! Bates (1951) ! SR 01/21/12
     kp2 = 10.0**(5.3541 - 0.019840*t1k - 1979.5/t1k)*h2po4t/hpo4t
                                                                 ! Bates and Acree (1943) ! SR 01/21/12
     kp3 = 10.0**( - 12.38)*hpo4t/po4t
                                      ! Dean (1985 )! SR 01/01/12
!**** pH evaluation
     !PHT = -PH(K,I)-2.1
     !IF (PH(K,I) <= 0.0) PHT = -14.0
     pht = -Phsed - 2.1
     if(Phsed<=0.0)pht = -14.0
     incr1 = 10.0
     do n = 1, 3
         f1 = 1.0
         incr1 = incr1/10.0
         iter1 = 0
         do while (f1>0.0 .AND. iter1<12)
             pht = pht + incr1
             hion = 10.0**pht
             f1 = cart*k1*(hion + 2.0*k2)/(hion*hion + k1*hion + k1*k2)        &
                & + kw/hion - alkt - hion/ht                                        ! SR 01/01/12
                                      ! SR 01/01/12
                                            ! SR 01/01/12
             if(ammonia_buffering)f1 = f1 + ammt*kamm/(hion + kamm)
                 ! SR 01/01/12
                                        ! SR 01/01/12
                                                                                ! SR 01/01/12
             if(phosphate_buffering)                                           &
              & f1 = f1 + phost*(kp1*kp2*hion + 2*kp1*kp2*kp3 - hion*hion*hion)&
              & /(hion*hion*hion + kp1*hion*hion + kp1*kp2*hion + kp1*kp2*kp3)
                 ! SR 01/01/12
             if(om_buffering)then
                                 ! SR 01/01/12
                 do ja = 1, nag
                        ! SR 01/01/12
                     f1 = f1 + omct*SDEN(ja)                                   &
                        & *(1.0/(1.0 + hion*(10.0**PK(ja))) - 1.0/             &
                        & (1.0 + (10.0**(PK(ja)-4.5))))                                                  ! SR 01/01/12
                 enddo
                   ! SR 01/01/12
             endif
                 ! SR 01/01/12
             iter1 = iter1 + 1
         enddo
         pht = pht - incr1
     enddo
!**** pH, carbon dioxide, bicarbonate, and carbonate concentrations
     hion = 10.0**pht
     Phsed = -pht
     co2sed = Ticsed/(1.0 + k1/hion + k1*k2/(hion*hion))
     hco3sed = Ticsed/(1.0 + hion/k1 + k2/hion)
     co3sed = Ticsed/((hion*hion)/(k1*k2) + hion/k2 + 1.0)
 
     return
 
 
     entry CEMAWINDINDUCEDSEDIMENTRESUSPENSION
 
     epsilon = 0.0
     fetchw = FETCHD(segnumi, jb)
     if(COS(PHI(jw) - PHI0(segnumi))<0.0)fetchw = FETCHU(segnumi, jb)
     fetchw = MAX(fetchw, BI(kt, segnumi), DLX(segnumi))
     u2 = WIND(jw)*WSC(segnumi)*WIND(jw)*WSC(segnumi) + nonzero
     coef1 = 0.53*(g*DEPTHB(KB(segnumi), segnumi)/u2)**0.75
     coef2 = 0.0125*(g*fetchw/u2)**0.42
     coef3 = 0.833*(g*DEPTHB(KB(segnumi), segnumi)/u2)**0.375
     coef4 = 0.077*(g*fetchw/u2)**0.25
     hs = 0.283*u2/g*0.283*TANH(coef1)*TANH(coef2/TANH(coef1))
  !TS    = 2.0*PI*U2/G*1.2*  TANH(COEF3)*TANH(COEF4/TANH(COEF3))
     ts = 2.0*pi*SQRT(u2)/g*1.2*TANH(coef3)*TANH(coef4/TANH(coef3))     ! cb 5/9/14
     lw0 = g*ts*ts/(2.0*pi)
 
     lw1 = lw0
     lw = lw0*TANH(2.0*pi*DEPTHB(KB(segnumi), segnumi)/lw1)
     do while (ABS(lw - lw1)>0.001)
         lw1 = lw
         lw = lw0*TANH(2.0*pi*DEPTHB(KB(segnumi), segnumi)/lw1)
     enddo
     coef = MIN(710.0, 2.0*pi*DEPTHB(KB(segnumi), segnumi)/lw)
     uorb = pi*hs/ts*100.0/SINH(coef)
     tau = 0.003*uorb*uorb
     if(tau - taucrpom>0.0)epsilon = MAX(0.0, 0.008/49.0*(tau - taucrpom)      &
                                   & **3*10000.0/dlt)
     sd_e = epsilon*DLX(segnumi)*BI(KB(segnumi), segnumi)                      &
          & /VOL(KB(segnumi), segnumi)                                          ! SD_E: g/m^2/s
 
     return
 
     entry CEMABOTTOMSCOURRESUSPENSION
 
     if(cao_method)then
         reyn_resusp = dia_pom*SQRT(spgrav_pom*g*dia_pom)
         if(reyn_resusp<6.0)then
             crshields = 0.1414*reyn_resusp**( - 0.2306)
         elseif(reyn_resusp>=6.0 .AND. reyn_resusp<=282.8)then
             crshields = (1.0 + (0.0223*reyn_resusp)**2.8358)                  &
                       & **0.3542/(3.0946*reyn_resusp**0.6769)
         elseif(reyn_resusp>282.8)then
             crshields = 0.045
         endif
     endif
     molvisc_h2o = 1.79E-6*EXP(0.0266*SD_T(1))
     shields = sd_taubot/(g*(spgrav_pom - 1.0)*dia_pom)
 
     vscour = 0.00033*(shields/crshields - 1.0)*(spgrav_pom - 1.0)             &
            & **0.6*g**0.6*dia_pom**0.8/molvisc_h2o
     c_bottom = (C2(KB(segnumi), segnumi, nlpom)                               &
              & + C2(KB(segnumi), segnumi, nrpom))                             &
              & *DEXP(POMS(jw)*H(KB(segnumi), jw)/DZ(KB(segnumi) - 1, segnumi))
 
     if(spgrav_pom<1.2)then
         c_bottom2 = 1.0
     elseif(spgrav_pom>=1.2 .AND. spgrav_pom<1.8)then
         c_bottom2 = 1.0*(1.8 - spgrav_pom)/0.6 + 3.0*(spgrav_pom - 1.2)/0.6
     elseif(spgrav_pom>=1.8 .AND. spgrav_pom<=2.2)then
         c_bottom2 = 3.0*(2.2 - spgrav_pom)/0.4 + 5.0*(spgrav_pom - 1.8)/0.4
     elseif(spgrav_pom>2.2)then
         c_bottom2 = 5.0
     endif
 
     c_bottom = DMIN1(c_bottom, c_bottom2)
 
     if(vscour>0.0)then
         sd_e = c_bottom*vscour
     else
         sd_e = 0.0
     endif
 
 
 
     end subroutine CEMASEDIMENTDIAGENESIS
