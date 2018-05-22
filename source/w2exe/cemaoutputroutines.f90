!*==cemaoutputroutines.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     recursive subroutine cemaoutputroutines
 
    ! Type declarations
     use MAIN
     use GLOBAL
     use SCREENC
     use GEOMC
     use KINETIC
     use NAMESC
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     character(3), dimension(imx), save :: crackstatus
     integer, save :: genbodnum, ngas
!
!*** End of declarations rewritten by SPAG
!
 
     return
 
     entry WRITECEMASEDIMENTMODELOUTPUT
 
        !Bed Elevation
     if(writebesnp)then
         write(cemasnpoutfiln, '("Bed elevation(m)")')
         write(cemasnpoutfiln, '("JDAY = ",f14.6)')jday
         write(cemasnpoutfiln, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
         write(cemasnpoutfiln, '(<IMX>f10.4)')                                 &
             & (BEDELEVATION(segnumi), segnumi = 1, imx)
 
         write(cemasnpoutfiln, '("Bed elevation layer(m)")')
         write(cemasnpoutfiln, '("JDAY = ",f14.6)')jday
         write(cemasnpoutfiln, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
         write(cemasnpoutfiln, '(<IMX>f10.4)')                                 &
             & (BEDELEVATIONLAYER(segnumi), segnumi = 1, imx)
 
            !Sediment concentration near bottom
            !Write(CEMASNPOutFilN,'("Sediment Concentration (gm/m³)")')
            !Write(CEMASNPOutFilN,'("JDAY = ",f8.2)')JDAY
            !Write(CEMASNPOutFilN,'(<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            !Write(CEMASNPOutFilN,'(<IMX>f10.3)')(CEMASedConc(SegNumI,KB(SegNumI)), SegNumI = 1, IMX)
            !Write(CEMASNPOutFilN,'("Density (gm/m³)")')
            !Write(CEMASNPOutFilN,'("JDAY = ",f8.2)')JDAY
            !Write(CEMASNPOutFilN,'(<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            !Write(CEMASNPOutFilN,'(<IMX>f10.3)')(rho(KB(SegNumI),SegNumI), SegNumI = 1, IMX)
            !Write(CEMASNPOutFilN,'("Density (gm/m³)")')
            !Write(CEMASNPOutFilN,'("JDAY = ",f8.2)')JDAY
            !Write(CEMASNPOutFilN,'(<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            !Write(CEMASNPOutFilN,'(<IMX>f10.3)')(rho(KB(SegNumI)-1,SegNumI), SegNumI = 1, IMX)
     endif
 
        !Bed porosity
     if(writepwsnp)then
         write(cemasnpoutfiln, '("Bed porosity (%)")')
         write(cemasnpoutfiln, '("JDAY = ",f8.2)')jday
         write(cemasnpoutfiln, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
         write(cemasnpoutfiln, '(<IMX>f10.4)')                                 &
             & (BEDPOROSITY(segnumi)*100, segnumi = 1, imx)
         write(cemasnpoutfiln, '("Total Volume of Sediments = ",e14.6," m3")') &
             & totalsedimentsinbed
         write(cemasnpoutfiln, '("Total Porewater Volume = ",e14.6," m3")')    &
             & totalporewatvolume
         write(cemasnpoutfiln, '("Total Porewater Removed = ",e14.6," m3")')   &
             & totalporewatremoved
!        Write(CEMATSR1OutFilN,'(f14.6 , "," ,i5,",",i5, ",", f14.6, ",",
!        f14.6, ",", f14.6, ",", f14.6, ",", f14.6, ",", f14.6, ",", f14.6,
!        ",", f14.6)')JDAY, KTWB(1), KB(4), z(4), BedElevation(4),
!        BedElevationLayer(4), PorewaterRelRate(4)*86400.d0,
!        BedPorosity(4)*100, BedConsolidRate(4)*86400.d0, EL(KTWB(1),4)-Z(4),
!        CEMACumPWRelease(4)
     endif
 
     write(cemabtmlayfiln, '("Bottom Layer (KB)")')
     write(cemabtmlayfiln, '("JDAY = ",f8.2)')jday
     tempcntr1 = tempcntr1 + 1
     if(tempcntr1>10)then
         tempcntr1 = 0
         write(cemabtmlayfiln, '(a)')                                          &
     &"________________________________________________________________________&
     &________"
         write(cemabtmlayfiln, '(<IMX>i5)')(segnumi, segnumi = 1, imx)
     endif
     write(cemabtmlayfiln, '(<IMX>i5)')(KB(segnumi), segnumi = 1, imx)
 
     return
 
     entry WRITECEMASEDIMENTFLUXOUTPUT
 
     if(writecemamftsedflx)then
            ! testing start
           ! write(4478,'(f8.3,(1000f10.3)))')jday,(MFTSedFlxVars(SegNumI,26), SegNumI = 2, IMX-1), (MFTSedFlxVars(SegNumI,27), SegNumI = 2, IMX-1),  &
           !   (MFTSedFlxVars(SegNumI,28), SegNumI = 2, IMX-1),(MFTSedFlxVars(SegNumI,29), SegNumI = 2, IMX-1),  &
           ! (MFTSedFlxVars(SegNumI,30), SegNumI = 2, IMX-1),(MFTSedFlxVars(SegNumI,25), SegNumI = 2, IMX-1),(MFTSedFlxVars(SegNumI,31),MFTSedFlxVars(SegNumI,32), SegNumI = 2, IMX-1),  &
           ! (MFTSedFlxVars(SegNumI,33),MFTSedFlxVars(SegNumI,34), SegNumI = 2, IMX-1),(MFTSedFlxVars(SegNumI,10),MFTSedFlxVars(SegNumI,11), SegNumI = 2, IMX-1),  &
           ! (MFTSedFlxVars(SegNumI,8),MFTSedFlxVars(SegNumI,9), SegNumI = 2, IMX-1),(MFTSedFlxVars(SegNumI,35),MFTSedFlxVars(SegNumI,36), SegNumI = 2, IMX-1)
 
           ! write(4479,'(f8.3,(1000f10.3)))')jday,(MFTSedFlxVars(SegNumI,37),MFTSedFlxVars(SegNumI,38), SegNumI = 2, IMX-1),(MFTSedFlxVars(SegNumI,39),MFTSedFlxVars(SegNumI,40), SegNumI = 2, IMX-1), &
           !   (MFTSedFlxVars(SegNumI,41),MFTSedFlxVars(SegNumI,42), SegNumI = 2, IMX-1),(MFTSedFlxVars(SegNumI,43),MFTSedFlxVars(SegNumI,44), SegNumI = 2, IMX-1)
           ! write(4480,'(f8.3,(1000f10.3)))')jday,(MFTSedFlxVars(SegNumI,45),MFTSedFlxVars(SegNumI,46), SegNumI = 2, IMX-1)
 
         write(cemasedflxfiln2, '("JDAY = ",f8.2)')jday
         write(cemasedflxfiln2, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
         write(cemasedflxfiln2, '("SOD (gO2/m2/d)")')
         write(cemasedflxfiln2, '(<IMX>f10.4)')                                &
             & (MFTSEDFLXVARS(segnumi, 26), segnumi = 1, imx)
 
         write(cemasedflxfiln2, '("JDAY = ",f8.2)')jday
         write(cemasedflxfiln2, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
         write(cemasedflxfiln2, '("CSOD (gO2/m2/d)")')
         write(cemasedflxfiln2, '(<IMX>f10.4)')                                &
             & (MFTSEDFLXVARS(segnumi, 27), segnumi = 1, imx)
 
         write(cemasedflxfiln2, '("JDAY = ",f8.2)')jday
         write(cemasedflxfiln2, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
         write(cemasedflxfiln2, '("NSOD (gO2/m2/d)")')
         write(cemasedflxfiln2, '(<IMX>f10.4)')                                &
             & (MFTSEDFLXVARS(segnumi, 28), segnumi = 1, imx)
 
         write(cemasedflxfiln4, '(f8.2,<IMX>f10.4)')jday,                      &
             & (MFTSEDFLXVARS(segnumi, 26), segnumi = 1, imx)
 
         write(cemasedflxfiln5, '(f8.2,<IMX>f10.1)')jday,                      &
             & (MFTSEDFLXVARS(segnumi, 29), segnumi = 1, imx)                                              ! POC
         write(cemasedflxfiln6, '(f8.2,<IMX>f10.2)')jday,                      &
             & (MFTSEDFLXVARS(segnumi, 30), segnumi = 1, imx)                                              ! PON
         write(cemasedflxfiln7, '(f8.2,<IMX>f10.3)')jday,                      &
             & (MFTSEDFLXVARS(segnumi, 47), segnumi = 1, imx)                                              ! POP
 
         write(cemasedflxfiln1, '("JDAY = ",f8.2)')jday
         write(cemasedflxfiln1, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
         write(cemasedflxfiln1, '("POC (mg/l)")')
         write(cemasedflxfiln1, '(<IMX>f10.4)')                                &
             & (MFTSEDFLXVARS(segnumi, 29), segnumi = 1, imx)
 
         write(cemasedflxfiln1, '("JDAY = ",f8.2)')jday
         write(cemasedflxfiln1, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
         write(cemasedflxfiln1, '("PON (mg/l)")')
         write(cemasedflxfiln1, '(<IMX>f10.4)')                                &
             & (MFTSEDFLXVARS(segnumi, 30), segnumi = 1, imx)
 
         write(cemasedflxfiln1, '("JDAY = ",f8.2)')jday
         write(cemasedflxfiln1, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
         write(cemasedflxfiln1, '("SO4 (mg/l)")')
         write(cemasedflxfiln1, '(<IMX>f10.4)')                                &
             & (MFTSEDFLXVARS(segnumi, 25), segnumi = 1, imx)
 
         write(cemasedflxfiln3, '("JDAY = ",f8.2)')jday
         write(cemasedflxfiln3, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
         write(cemasedflxfiln3, '("H1 (m)")')
         write(cemasedflxfiln3, '(<IMX>f10.4)')                                &
             & (SD_AERLAYERTHICK(segnumi), segnumi = 1, imx)
 
         call CEMATEMPOUTPUT
 
     endif
 
     if(includecemagenbodconstituents)then
 
         write(cemaoutfiln5,                                                   &
             &'("##################### JDAY = ",f8.2," #####################")'&
            & )jday
         write(cemaoutfiln5, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
         write(cemaoutfiln5, '("Aerobic Layer")')
         do genbodnum = 1, numgenbodconstituents
             write(cemaoutfiln5, '(a)')SEDGENBODNAME(genbodnum)
             write(cemaoutfiln5, '(<IMX>f10.4)')                               &
                 & (SEDGENBODCONC(genbodnum, segnumi, 1), segnumi = 1, imx)
         enddo         !j
         write(cemaoutfiln5, '("Anaerobic Layer")')
         do genbodnum = 1, numgenbodconstituents
             write(cemaoutfiln5, '(a)')SEDGENBODNAME(genbodnum)
             write(cemaoutfiln5, '(<IMX>f10.4)')                               &
                 & (SEDGENBODCONC(genbodnum, segnumi, 2), segnumi = 1, imx)
         enddo         !j
!
!	       Write(1715,'(f8.4,",",f14.6,",",f14.6,",",f14.6)')jday,
!        SedGenBODConc(1,ds(1),1),SedGenBODConc(1,ds(1),2), O2(KB(ds(1)),ds(1))
     endif
 
     return
 
     entry CEMATEMPOUTPUT
 
     crackstatus = " NO"
     do segnumi = 1, imx
         if(CRACKOPEN(segnumi))crackstatus(segnumi) = "YES"
     enddo     !SegNumI
 
     write(cemaoutfiln1, '("************************  JDAY = ",f8.2)')jday
     write(cemaoutfiln1, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
     write(cemaoutfiln1, '("Radius (mm)")')
     write(cemaoutfiln1, '(<IMX>f10.4)')(BUBBLERADIUSSED(segnumi)*1000, segnumi&
                                      & = 1, imx)
        !Write(CEMAOutFilN1,'("Pressure Bubble")')
        !Write(CEMAOutFilN1,'(<IMX>f10.1)')(PresBubbSed(SegNumI), SegNumI = 1, IMX)
        !Write(CEMAOutFilN1,'("Pressure Critical")')
        !Write(CEMAOutFilN1,'(<IMX>f10.1)')(PresCritSed(SegNumI), SegNumI = 1, IMX)
     write(cemaoutfiln1, '("Cg (gm/m³)")')
     write(cemaoutfiln1, '(<IMX>f10.1)')(CGSED(segnumi), segnumi = 1, imx)
     write(cemaoutfiln1, '("C0 (gm/m³)")')
     write(cemaoutfiln1, '(<IMX>f10.1)')(C0SED(segnumi), segnumi = 1, imx)
     write(cemaoutfiln1, '("Ct (gm/m³)")')
     write(cemaoutfiln1, '(<IMX>f10.1)')(CTSED(segnumi), segnumi = 1, imx)
     write(cemaoutfiln1, '("Crack Open")')
     write(cemaoutfiln1, '(<IMX>a10)')(crackstatus(segnumi), segnumi = 1, imx)
        !Write(CEMAOutFilN1,'("Bubbles released")')
        !Write(CEMAOutFilN1,'(<IMX>i10)')(MFTBubbReleased(SegNumI), SegNumI = 1, IMX)
 
     write(cemaoutfiln2, '("Gas concentration (gm/m³) at JDAY = ",f8.2)')jday
     write(cemaoutfiln3, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
     do ngas = 1, 4
         if(ngas==1)write(cemaoutfiln2, '("H2S Concentration (gm/m³)")')
         if(ngas==2)write(cemaoutfiln2, '("CH4  Concentration (gm/m³)")')
         if(ngas==3)write(cemaoutfiln2, '("NH3  Concentration (gm/m³)")')
         if(ngas==4)write(cemaoutfiln2, '("CO2  Concentration (gm/m³)")')
         write(cemaoutfiln2, '(<IMX>f10.3)')(TCONC(ngas, segnumi), segnumi = 1,&
             & imx)
     enddo     !nGas
 
     ngas = 2  ! CH4 write
        !Write(CEMASedFlxFilN8,'(f8.2,<IMX>f10.2)')JDAY,(TConc(nGas,SegNumI), SegNumI = 1, IMX)
     write(cemasedflxfiln8, '(f8.2,<IMX>f10.2)')jday,                          &
         & (DISSOLVEDGASSEDIMENTS(ngas, segnumi), segnumi = 1, imx)
     ngas = 1  ! H2S write
        !Write(CEMASedFlxFilN9,'(f8.2,<IMX>f10.2)')JDAY,(TConc(nGas,SegNumI), SegNumI = 1, IMX)
     write(cemasedflxfiln9, '(f8.2,<IMX>f10.2)')jday,                          &
         & (DISSOLVEDGASSEDIMENTS(ngas, segnumi), segnumi = 1, imx)
     do ngas = 1, 4
         if(ngas==1)write(cemaoutfiln6, '("H2S Concentration (gm/m³)")')
         if(ngas==2)write(cemaoutfiln6, '("CH4  Concentration (gm/m³)")')
         if(ngas==3)write(cemaoutfiln6, '("NH3  Concentration (gm/m³)")')
         if(ngas==4)write(cemaoutfiln6, '("CO2  Concentration (gm/m³)")')
         write(cemaoutfiln6, '(<IMX>f10.3)')                                   &
             & (DISSOLVEDGASSEDIMENTS(ngas, segnumi), segnumi = 1, imx)
     enddo     !nGas
 
     write(cemaoutfiln3, '("Gas Release to Atmosphere at JDAY = ",f8.2)')jday
     write(cemaoutfiln3, '(<IMX>i10)')(segnumi, segnumi = 1, imx)
     do ngas = 1, 4
         if(ngas==1)write(cemaoutfiln3, '("H2S Release(gm/s)")')
         if(ngas==2)write(cemaoutfiln3, '("CH4 Release (gm/s)")')
         if(ngas==3)write(cemaoutfiln3, '("NH3 Release (gm/s)")')
         if(ngas==4)write(cemaoutfiln3, '("CO2 Release (gm/s)")')
         write(cemaoutfiln3, '(<IMX>f10.3)')                                   &
             & (BRRATEAGASNET(segnumi, ngas), segnumi = 1, imx)
     enddo     !nGas
 
     do jw = 1, nwb
         write(cemaoutfilbub, '(F12.4,",",I5,",",6(E12.5,","))')jday, jw,      &
             & (BUBBLERELWB(jw, j), j = 1, 4), gasreleasech4/1000.                                                      !,GasReleaseCO2/1000.   ! SW 7/1/2017 BubbleRelWB kg  GasReleaseCH4 kg C
     enddo
     BUBBLERELWB = 0.0
 
 
        !Write(CEMAOutFilN4,'(<5>f10.3)')JDAY,(BRRateAGasNet(5, nGas), nGas = 1, 4)
        !Write(227,'(f10.5, ",",f10.5)')jday, c1(kb(6),6,7)/1000.d0
 
 
 
     end subroutine cemaoutputroutines
