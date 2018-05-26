!*==cemasedimentmodelw2.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine CEMASEDIMENTMODELW2
 
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
     use TRANS
     use CEMAVARS
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(8) :: betad, caks, capda, cemashieldsnumber, consolidratetemp,       &
              & critshldpar, d1, d2, depth, dks, dsks, nutemp, rgha, rghapda,  &
              & sdfz, sechange, sedecoeff, sgr, snetflx, srho, sscour, ssettle,&
              & sum1, surfaceareawb, tcrit_e, tcrit_s, tempvariable, tflow_b,  &
              & tsks, uctavg, uks, usr, value1, veltemp, volumeincreased,      &
              & volumenew, volumeofporewater1, volumeofporewater2,             &
              & volumeofsedimentbed1, volumeofsedimentbed2, volumepresent,     &
              & volumeupdated, wf, yalinp
     logical :: bottomupdated
     real :: dummy
     integer(4) :: genbodnum, regnnum
!
!*** End of declarations rewritten by SPAG
!
 
 
!    Common /CEMACmBl1/ CEMAShieldsNumber, SScour, SNetFlx
 
     entry SETUPCEMASEDIMENTMODEL
 
        ! determining width of sediment cells  ! cb 3/4/13
        !Do JW=1, NWB                SW 9/18/2017
        !  Do JB=BS(JW),BE(JW)
        !    !IU = CUS(JB)
        !    IU = us(jb)
        !    ID = DS(JB)
        !    do SegNumI=iu,id
        !        sedcellwidth(SegNumI)=B(2,SegNumI)
        !    end do
        !  end do
        !end do
 
        !Call to initialize region consolidation data files
     do regnnum = 1, numconsolidregns
         call INITIALIZEBEDCONSOLIDATIONFILES(1210 + regnnum,                  &
           & CONSOLIDRATEREGNFIL(regnnum))
     enddo     !RegnNum
 
     do regnnum = 1, numconsolidregns
         do segnumi = CONSREGSEGST(regnnum), CONSREGSEGEN(regnnum)
             CONSOLIDREGNNUM(segnumi) = regnnum
         enddo     !SegNumI
     enddo     !RegnNum
 
 
        !Compute total porewater content
     totalporewatvolume = 0.D00
     totalsedimentsinbed = 0.D00
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
                 totalporewatvolume = totalporewatvolume + BEDPOROSITY(segnumi)&
                                    & *BEDELEVATION(segnumi)                   &
                                    & *B(KB(segnumi), segnumi)*DLX(segnumi)
                 totalsedimentsinbed = totalsedimentsinbed +                   &
                                     & (1 - BEDPOROSITY(segnumi))              &
                                     & *BEDELEVATION(segnumi)                  &
                                     & *B(KB(segnumi), segnumi)*DLX(segnumi)
                   ! TotalPoreWatVolume = TotalPoreWatVolume + BedPorosity(SegNumI)*BedElevation(SegNumI)* sedcellwidth(SegNumI)*DLX(SegNumI)
                   ! TotalSedimentsInBed = TotalSedimentsInBed + (1-BedPorosity(SegNumI))*BedElevation(SegNumI)* sedcellwidth(SegNumI)*DLX(SegNumI)
             enddo     !SegNumI
         enddo     !JB
     enddo     !JW
 
        ! moved file opens and writes to CEMASedimentDiagenesis
 
     return
 
 
     entry CEMASEDIMENTMODEL
 
     porewaterrelrate = 0.D00
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
                 if(CEMASSAPPLIED(segnumi))then
                     CEMALAYERADDED(segnumi) = .FALSE.
                     CEMACUMPWRELEASE(segnumi) = 0.D00
                     CEMASSAPPLIED(segnumi) = .FALSE.
                 endif
             enddo     !SegNumO
         enddo     !JB
     enddo     !JW
 
        !Get Suspended sediment concentration
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
                 do k = kt, KB(segnumi)
 
                        !CEMASedConc(SegNumI, K) = C1(K,SegNumI,NSSS)
                     CEMASEDCONC(segnumi, k) = C1(k, segnumi, nmft)       ! cb 2/18/13
 
                 enddo     !k
             enddo     !SegNumI
         enddo     !JB
     enddo     !JW
 
        !If(EndBedConsolidation)Return
 
        !Get bed elevation change rate
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
 
             iu = CUS(jb)
             id = DS(jb)
 
             do segnumi = iu, id
 
                 if(ENDBEDCONSOLIDATION(segnumi))cycle
 
                 regnnum = CONSOLIDREGNNUM(segnumi)
                                                               !Constant value
                 if(CONSOLIDATIONTYPE(regnnum)==0)BEDCONSOLIDRATE(segnumi)     &
                  & = CONSTCONSOLIDRATE(regnnum)/86400.D00                              !m/d to m/s
                 if(CONSOLIDATIONTYPE(regnnum)==1)then         !Time varying value
                     call READBEDCONSOLIDATIONFILES(1210 + regnnum,            &
                       & consolidratetemp)
                     BEDCONSOLIDRATE(segnumi) = consolidratetemp/86400.D00    !m/d to m/s
                 endif
 
             enddo     !SegNumI
 
         enddo     !JB
     enddo     !JW
 
        !Compute net sediment bed elevation
     totalporewatremoved = 0.D00
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
 
                 if(ENDBEDCONSOLIDATION(segnumi))cycle
 
                 volumeofsedimentbed1 = BEDELEVATION(segnumi)                  &
                                      & *B(KB(segnumi), segnumi)*DLX(segnumi)
                    !VolumeofSedimentBed1 = BedElevation(SegNumI)*B(KT,SegNumI)*DLX(SegNumI)
                 volumeofporewater1 = volumeofsedimentbed1*BEDPOROSITY(segnumi)
                 BEDELEVATION(segnumi) = BEDELEVATION(segnumi)                 &
                   & - BEDCONSOLIDRATE(segnumi)*dlt
                 BEDELEVATIONLAYER(segnumi) = BEDELEVATIONLAYER(segnumi)       &
                   & - BEDCONSOLIDRATE(segnumi)*dlt
                 volumeofsedimentbed2 = BEDELEVATION(segnumi)                  &
                                      & *B(KB(segnumi), segnumi)*DLX(segnumi)
                    !VolumeofSedimentBed2 = BedElevation(SegNumI)*B(KT,SegNumI)*DLX(SegNumI)
                 BEDPOROSITY(segnumi) = 1.D0 - volumeofsedimentbed1/           &
                                      & volumeofsedimentbed2*                  &
                                      & (1 - BEDPOROSITY(segnumi))
                 volumeofporewater2 = volumeofsedimentbed2*BEDPOROSITY(segnumi)
                 totalporewatremoved = volumeofporewater1 - volumeofporewater2
                 porewaterrelrate(segnumi) = totalporewatremoved/dlt    !m³/s
                 CEMACUMPWRELEASE(segnumi) = CEMACUMPWRELEASE(segnumi)         &
                   & + totalporewatremoved                                                          !m³
                 CEMACUMPWRELEASERATE(segnumi) = CEMACUMPWRELEASE(segnumi)/dlt    !m³/s
 
                 if(BEDPOROSITY(segnumi)<0.D00)then
                     totalporewatremoved = 0.D00
                     BEDPOROSITY(segnumi) = 0.D00
                     ENDBEDCONSOLIDATION(segnumi) = .TRUE.
 
                     write(w2err, '(a)')                                       &
                          &"Bed reached maximum compaction level. Porosity = 0"
                     write(w2err, '(a,i4)')                                    &
                          &"No more bed compaction will occur at segment = ",  &
                         & segnumi
                     error_open = .TRUE.     ! SR 7/27/2017
                 endif
 
             enddo     !SegNumI
         enddo     !JB
     enddo     !JW
 
        !Compute net sediment and porewater volume
     totalporewatvolume = 0.D00
     totalsedimentsinbed = 0.D00
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
 
                 if(ENDBEDCONSOLIDATION(segnumi))cycle
                 totalporewatvolume = totalporewatvolume + BEDPOROSITY(segnumi)&
                                    & *BEDELEVATION(segnumi)                   &
                                    & *B(KB(segnumi), segnumi)*DLX(segnumi)
                 totalsedimentsinbed = totalsedimentsinbed +                   &
                                     & (1 - BEDPOROSITY(segnumi))              &
                                     & *BEDELEVATION(segnumi)                  &
                                     & *B(KB(segnumi), segnumi)*DLX(segnumi)
                    !TotalPoreWatVolume = TotalPoreWatVolume + BedPorosity(SegNumI)*BedElevation(SegNumI)*B(KT,SegNumI)*DLX(SegNumI)
                    !TotalSedimentsInBed = TotalSedimentsInBed + (1-BedPorosity(SegNumI))*BedElevation(SegNumI)*B(KT,SegNumI)*DLX(SegNumI)
 
             enddo     !SegNumI
         enddo     !JB
     enddo     !JW
 
        !Get Erosion/Deposition
     if(cemasedimentprocessesinc)call CEMASEDIMENTPROCESSES
 
        !Update Suspended sediment concentration
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
                 do k = kt, KB(segnumi)
 
                        !C1(K,SegNumI,NSSS) = CEMASedConc(SegNumI, K)
                     C1(k, segnumi, nmft) = CEMASEDCONC(segnumi, k)      ! cb 2/18/13
 
                 enddo     !k
             enddo     !SegNumI
         enddo     !JB
     enddo     !JW
 
 
     return
 
     entry CEMASEDIMENTPROCESSES
 
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
                 if(ENDBEDCONSOLIDATION(segnumi))cycle
 
                    !Calculate scouring
                 call CEMACALCULATESCOUR
 
                    !Calculate settling
                 if(includefftlayer .AND. fftactive)                           &
                  & ssettle = fftlayersettvel*CEMASEDCONC(segnumi, KB(segnumi))
                 if(includefftlayer .AND. .NOT.fftactive)ssettle = 0.D00
                 if(.NOT.includefftlayer)                                      &
                  & ssettle = cemasedimentsvelocity*CEMASEDCONC(segnumi,       &
                  & KB(segnumi))
 
                    !Net Flux
                 snetflx = sscour - ssettle     !gm/m2/s
 
                    !If(SNetFlx > 1)SNetFlx = 1
                    !If(SNetFlx < -1)SNetFlx = -1
 
                    !New sediment concentration
                 CEMASEDCONC(segnumi, KB(segnumi))                             &
                   & = CEMASEDCONC(segnumi, KB(segnumi))                       &
                   & + snetflx/H(KB(segnumi), jw)*dlt
 
                 BEDELEVATION(segnumi) = BEDELEVATION(segnumi)                 &
                   & - snetflx*dlt/(cemasedimentdensity*1000.D0)                                                !Density from kg/m³ to gm/m³
                 BEDELEVATIONLAYER(segnumi) = BEDELEVATIONLAYER(segnumi)       &
                   & - snetflx*dlt/(cemasedimentdensity*1000.D0)                                                          !Density from kg/m³ to gm/m³
 
             enddo     !SegNumI
         enddo     !JB
     enddo     !JW
 
     return
 
     entry CEMACALCULATESCOUR
 
     sgr = cemasedimentdensity/1000.D0
     dks = cemaparticlesize
     sscour = 0.D00
 
        !Calculate Erosion
     if(cemasedimenttype==1)then        !Cohesive Sediments
        	! Based on formulation in EFDC
!		!     Sediment Scour for Cohesive Sediments
         critshldpar = 0.1
!
!		!     Estimate critical shear stress for erosion
         tcrit_e = 0.D0
         tcrit_s = 0.D0
!		!Mass Erosion
         if(cemasedimentdensity>1065)then
             tcrit_s = (cemasedimentdensity/1000 - 1.065)**0.2
             tcrit_s = 0.001*(0.883*tcrit_s + 0.05)                    !Bulk Density divided by 1000 to convert unit to gm/cm3
         endif
!		!Surface Erosion
         if(cemasedimentdensity>1013)                                          &
          & tcrit_e = 0.001*(9.808*cemasedimentdensity/1000.0 - 9.934)                     !Bulk Density divided by 1000 to convert unit to gm/cm3
 
         veltemp = 0.5*(U(KB(segnumi), segnumi)                                &
                 & + U(KB(segnumi - 1), segnumi - 1))
         veltemp = SQRT(veltemp**2)
!		!Limit bottom velocity
         veltemp = MIN(veltemp, 0.01)
         veltemp = SQRT(g)*veltemp/FRIC(segnumi)                   !###################################CHECK
 
         tflow_b = RHO(KB(segnumi), segnumi)*(veltemp**2)
 
         if(tflow_b>tcrit_s)then
!			!Mass    Erosion
             sscour = 0.62
             sscour = MAX(sscour, 0.0)
         elseif(tflow_b>tcrit_e)then
!				!Surface Erosion
             tempvariable = 0.23*DEXP                                          &
                          & (0.198/(cemasedimentdensity/1000.0 - 1.0023))
             tempvariable = (1/360.0)*(10**tempvariable)                                !1/360 to convert from mg/hr-cm2 to gm/s-m2
             sscour = (tflow_b - tcrit_e)/tcrit_e
             sscour = (sscour)**critshldpar
             sscour = tempvariable*sscour
             sscour = DMAX1(sscour, 0.0)
         else
!				!No     Erosion
             sscour = 0.D0
         endif
     endif
 
     if(cemasedimenttype==2)then        !Non-Cohesive Sediments
         call COMPUTESHIELDSNUMBER
         uks = SQRT(cemashieldsnumber*(sgr - 1.0)*g*dks)                !Critical shear velocity (Shield's diagram)
         nutemp = 1.79E-6*DEXP( - 0.0266*0.5*(T1(KB(segnumi), segnumi)))
!
         sum1 = 0.D00
         kt = KTWB(jw)
         d1 = H(k, jw) - Z(segnumi)
         depth = 0.D00
         do k = kt, KB(segnumi)
             veltemp = 0.5*(U(k, segnumi) + U(k, segnumi - 1))
             veltemp = DSQRT(veltemp**2)
!
             sum1 = sum1 + veltemp*d1
             depth = depth + d1
             d1 = H(k + 1, jw)
         enddo                 !k
!
         rgha = 0.01*depth
         rghapda = rgha + 0.4*rgha
!
         uctavg = sum1/depth                                            !Depth Averaged velocity Pg. 34
!		!Limit bottom velocity
         uctavg = MIN(uctavg, 0.01)
         if(MANNINGS_N(jw))then
             hrad = BH1(kt, i)/(B(KTI(i), i) - B(kt + 1, i) + 2.*AVH1(kt, i))
             gc2 = g*FRIC(i)*FRIC(i)/hrad**0.33333333
         else
             gc2 = 0.0
             if(FRIC(i)/=0.0)gc2 = g/(FRIC(i)*FRIC(i))
         endif
         usr = SQRT(gc2)*uctavg
!		!
!		!Erosion process
         dsks = dks*((sgr - 1.0)*g/nutemp**2.0)**(1.0/3.0)                              !Dimensionless particle diameter Pg. 34
         tsks = (usr**2.0 - uks**2.0)/(uks**2.0)                                        !Transport stage parameter Pg. 34
         wf = cemasedimentsvelocity
         srho = cemasedimentdensity
!
         if(wf/usr>=0.1 .AND. wf/usr<1.0)then
             betad = 1.0 + 2.0*(cemasedimentsvelocity/usr)**2.0                                 !Equation 68 Page 39
         else
             betad = 1.0
         endif
         sdfz = betad*AZ(KB(segnumi) - 1, segnumi)                          !Vertical mass diffusion coefficient Equation 67 Page 39
         caks = 0.0                    !1.0e+20				! Near Bed Concentration Equation 54 Page 34
         if(tsks>0.0)caks = 0.015*(dks/rgha)*(tsks**1.5/dsks**0.3)*srho*1000.0                   !gm/m^3
!		!
!		!Near Bed concentration extrapolated from suspended sediment Equation 31
!		!Linear interpolation Ref (
!        Page 21
         d1 = H(KB(segnumi), jw)
         d2 = H(KB(segnumi) - 1, jw)
         sedecoeff = (d1 + 0.5*d2 - rghapda)/(0.5*(d1 + d2))
!
         capda = (1.0 - sedecoeff)*CEMASEDCONC(segnumi, KB(segnumi))           &
               & + sedecoeff*CEMASEDCONC(segnumi, KB(segnumi) - 1)
 
         if(capda<0)capda = 0.D0
!		!
!		!Erosion Equation 12 Page 9
         if(caks>capda)sscour = -betad*sdfz*((capda - caks)/(rghapda - rgha))
 
     endif
 
     return
 
     entry COMPUTESHIELDSNUMBER
!
     nutemp = 1.79E-6*DEXP( - 0.0266*0.5*(T1(KB(segnumi), segnumi)))
     yalinp = ((cemasedimentdensity - 1000.0)                                  &
            & *g*cemaparticlesize**3.0/(1000.0*nutemp**2.0))**0.5
 
     selectcase(INT(yalinp))
!		!
!		!Y    < 100
     case(:100)
         value1 = 0.041*(LOG10(yalinp))**2.0 - 0.356*LOG10(yalinp) - 0.977
         cemashieldsnumber = (10.0)**value1
!		!
!		!100  < Y < 3000
     case(101:3000)
         value1 = 0.132*LOG10(yalinp) - 1.804
         cemashieldsnumber = (10.0)**value1
!		!
!		!Y    > 3000
     case(3001:)
         cemashieldsnumber = 0.045
     endselect
 
     return
 
     entry CEMAUPDATEVERTICALLAYERING
 
        !Compute current volume
     volumepresent = 0.D00
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
                 do k = kt, KB(segnumi)
                     if(k==kt)volumepresent = volumepresent + B(k, segnumi)    &
                      & *DLX(segnumi)*(H(k, jw) - Z(segnumi))
                     if(k/=kt)volumepresent = volumepresent + B(k, segnumi)    &
                      & *DLX(segnumi)*H(k, jw)
                 enddo     !K
             enddo     !SegNumI
         enddo     !JB
     enddo     !JW
 
     bottomupdated = .FALSE.
 
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
 
                 if(ENDBEDCONSOLIDATION(segnumi))cycle
 
                 if(DABS(BEDELEVATIONLAYER(segnumi))                           &
                  & >layeraddthkfrac*H(KB(segnumi), jw))then                                      !Add water layer
 
                     BEDELEVATIONLAYER(segnumi) = BEDELEVATIONLAYER(segnumi)   &
                       & + H(KB(segnumi), jw)                                                       !BedElevationLayer is negative due to bed consolidation
                     KB(segnumi) = KB(segnumi) + 1       !Add layer
 
                     if(KB(segnumi)>kmx - 1)then
                         KB(segnumi) = KB(segnumi) - 1
                         write(w2err, '(a,i4)')                                &
                              &"Maximum K layer reached for Segment = ",       &
                             & segnumi
                         error_open = .TRUE.     ! SR 7/27/2017
                     else
                         k = KB(segnumi)
                            !Copy properties from old bottom layer
                         B(k, segnumi) = B(k - 1, segnumi)
                         RHO(k, segnumi) = RHO(k - 1, segnumi)
                         H1(k, segnumi) = H(k, jw)
                         H2(k, segnumi) = H(k, jw)
                         AVH1(k, segnumi)                                      &
                           & = (H1(k, segnumi) + H1(k - 1, segnumi))*0.5
                         if(.NOT.TRAPEZOIDAL(jw))then
                             BH1(k, segnumi) = BH1(k - 1, segnumi)             &
                               & - BNEW(k - 1, segnumi)*H1(k - 1, segnumi)                                                                  ! SW 1/23/06
                         else
                             call GRID_AREA1(EL(k, segnumi) - Z(segnumi),      &
                               & EL(k - 1, segnumi), BH1(k, segnumi), dummy)                                                                                            !SW 08/03/04
                             BH1(k, segnumi) = 0.25*H1(k - 1, jw)              &
                               & *(BB(k, segnumi) + 2.*B(k - 1, segnumi)       &
                               & + BB(k - 1, segnumi))
                         endif
                         VOL(k, segnumi) = BH1(kt, segnumi)*DLX(segnumi)
                         BKT(segnumi) = BH1(kt, segnumi)/H1(kt, segnumi)
                         bi(kt:KB(segnumi), segnumi)                           &
                           & = B(kt:KB(segnumi), segnumi)                               ! SW 8/26/05
                            !T1(K,SegNumI)               =   T1(K-1,SegNumI)
                            !T2(K,SegNumI)               =   T2(K-1,SegNumI)
                            !C1(K,SegNumI,CN(1:NAC))     =   C1(K-1,SegNumI,CN(1:NAC))
                            !C2(K,SegNumI,CN(1:NAC))     =   C2(K-1,SegNumI,CN(1:NAC))
                         AZ(k - 1, segnumi) = AZ(k - 2, segnumi)
 
                            !Copy Constituent properties from old bottom layer
                            !NH4(K,SegNumI)                =   NH4(K-1,SegNumI)
                            !NO3(K,SegNumI)                =   NO3(K-1,SegNumI)
                            !O2(K,SegNumI)                =   O2(K-1,SegNumI)
                            !C2(K,SegNumI,NDO)     =   C2(K-1,SegNumI,NDO)
 
                            !SP added to test 04/01/2010
                         DX(k, segnumi) = DX(k - 1, segnumi)
                            !End SP added to test 04/01/2010
 
                         bottomupdated = .TRUE.
                         CEMALAYERADDED(segnumi) = .TRUE.
 
                         write(cemalogfiln, *)"Layer added at Segment Number ",&
                             & segnumi, "on ", jday
 
                     endif
 
                 endif
 
             enddo     !SegNumI
         enddo     !JB
     enddo     !JW
 
        !Compute new volume
     volumenew = 0.D00
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
                 do k = kt, KB(segnumi)
                     if(k==kt)volumenew = volumenew + B(k, segnumi)            &
                      & *DLX(segnumi)*(H(k, jw) - Z(segnumi))
                     if(k/=kt)volumenew = volumenew + B(k, segnumi)            &
                      & *DLX(segnumi)*H(k, jw)
                 enddo     !K
             enddo     !SegNumI
         enddo     !JB
     enddo     !JW
 
     if(bottomupdated)then      !Update surface elevation to accomodate for layer addition
 
         volumeincreased = volumenew - volumepresent
 
            !Calculate surface area
         surfaceareawb = 0.D00
         do jw = 1, nwb
             kt = KTWB(jw)
             do jb = BS(jw), BE(jw)
                 iu = CUS(jb)
                 id = DS(jb)
                 do segnumi = iu, id
                     surfaceareawb = surfaceareawb + B(kt, segnumi)            &
                                   & *DLX(segnumi)
                 enddo     !SegNumI
             enddo     !JB
         enddo     !JW
 
         sechange = volumeincreased/surfaceareawb
         do jw = 1, nwb
             kt = KTWB(jw)
             do jb = BS(jw), BE(jw)
                 iu = CUS(jb)
                 id = DS(jb)
                 do segnumi = iu, id
                     Z(segnumi) = Z(segnumi) + sechange
                 enddo     !SegNumI
             enddo     !JB
         enddo     !JW
 
            !Update other variables
         do jw = 1, nwb
             kt = KTWB(jw)
             do jb = BS(jw), BE(jw)
                 iu = CUS(jb)
                 id = DS(jb)
                 do segnumi = iu, id
                     SZ(segnumi) = Z(segnumi)
 
                     do k = kt, KB(segnumi)
                         if(k==kt)H1(k, segnumi) = H(k, jw) - Z(segnumi)
                         if(k/=kt)H1(k, segnumi) = H(k, jw)
                         BH1(k, segnumi) = B(k, segnumi)*H1(k, segnumi)
                         BH2(k, segnumi) = BH1(k, segnumi)
                         BHR1(k, segnumi) = B(k, segnumi)*H1(k, segnumi)
                         BHR2(k, segnumi) = BHR1(k, segnumi)
                     enddo     !K
                 enddo     !SegNumI
             enddo     !JB
         enddo     !JW
 
            !Move KB layer to padded cells
         do jw = 1, nwb
             kt = KTWB(jw)
             do jb = BS(jw), BE(jw)
                 iu = US(jb)
                 id = DS(jb)
                 KB(iu - 1) = KB(iu)
                 KB(id + 1) = KB(id)
             enddo      !JB
         enddo     !JW
         kbi = KB
 
         do jw = 1, nwb
             kt = KTWB(jw)
             do jb = BS(jw), BE(jw)
                 iu = US(jb)
                 id = DS(jb)
                 do i = iu - 1, id
                     KBMIN(i) = MIN(KB(i), KB(i + 1))
                 enddo      !I
                 KBMIN(id + 1) = KBMIN(id)
             enddo      !JB
         enddo     !JW
 
            !Compute updated volume
         volumeupdated = 0.D00
         do jw = 1, nwb
             kt = KTWB(jw)
             do jb = BS(jw), BE(jw)
                 iu = CUS(jb)
                 id = DS(jb)
                 do segnumi = iu, id
                     do k = kt, KB(segnumi)
                         if(k==kt)volumeupdated = volumeupdated + B(k, segnumi)&
                          & *DLX(segnumi)*(H(k, jw) - Z(segnumi))
                         if(k/=kt)volumeupdated = volumeupdated + B(k, segnumi)&
                          & *DLX(segnumi)*H(k, jw)
                     enddo     !K
                 enddo     !SegNumI
             enddo     !JB
         enddo     !JW
 
         do jw = 1, nwb
             do jb = BS(jw), BE(jw)
                 do i = CUS(jb) - 1, DS(jb)
                     DEPTHB(KTWB(jw), i) = H1(KTWB(jw), i)
                     if(kbi(i)<KB(i))DEPTHB(KTWB(jw), i)                       &
                      & = (H1(KTWB(jw), i) - (EL(kbi(i) + 1, i)                &
                      & - EL(KB(i) + 1, i)))                                                                    ! SW 1/23/06
                     DEPTHM(KTWB(jw), i) = H1(KTWB(jw), i)*0.5
                     if(kbi(i)<KB(i))DEPTHM(KTWB(jw), i)                       &
                      & = (H1(KTWB(jw), i) - (EL(kbi(i) + 1, i)                &
                      & - EL(KB(i) + 1, i)))*0.5                                                                    ! SW 1/23/06
                     do k = KTWB(jw) + 1, kmx
                         DEPTHB(k, i) = DEPTHB(k - 1, i) + H1(k, i)
                         DEPTHM(k, i) = DEPTHM(k - 1, i)                       &
                                      & + (H1(k - 1, i) + H1(k, i))*0.5
                     enddo
                 enddo
             enddo
         enddo
 
         do jw = 1, nwb
             kt = KTWB(jw)
             do jb = BS(jw), BE(jw)
                 iu = CUS(jb)
                 id = DS(jb)
                 do segnumi = iu, id
                        !KTI(SegNumI) = KT
                 enddo     !SegNumI
             enddo     !JB
         enddo     !JW
 
         do jw = 1, nwb
             do i = 1, NISNP(jw)
                 KBR(jw) = MAX(KB(ISNP(i, jw)), KBR(jw))
             enddo
         enddo
 
            !SP added to test 04/01/2010
         do jw = 1, nwb
             call INTERPOLATION_MULTIPLIERS
         enddo
 
     endif
 
     if(bottomupdated .AND. movefftlayerdown)call MOVEFFTLAYERCONSOLID
 
     return
 
     entry COMPUTECEMARELATEDSOURCESINKS
 
        !VOLCEMA = 0.d00
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
 
                 if(NUMCEMAPWINST(segnumi)>5)then
                     APPLYCEMAPWRELEASE(segnumi) = .FALSE.
                     NUMCEMAPWINST(segnumi) = 0
                 endif
 
                 if(CEMALAYERADDED(segnumi))then
                     CEMACUMPWTORELEASE(segnumi) = CEMACUMPWRELEASE(segnumi)
                     CEMACUMPWRELEASED(segnumi) = CEMACUMPWTORELEASE(segnumi)  &
                       & /5.D0                                                    !Release over 100 installments
                     APPLYCEMAPWRELEASE(segnumi) = .TRUE.
                     CEMASSAPPLIED(segnumi) = .TRUE.
                 endif
 
                 if(APPLYCEMAPWRELEASE(segnumi))then
                     NUMCEMAPWINST(segnumi) = NUMCEMAPWINST(segnumi) + 1
                     QSS(KB(segnumi), segnumi) = QSS(KB(segnumi), segnumi)     &
                       & + CEMACUMPWRELEASED(segnumi)/dlt                         !m³/s
                     TSS(KB(segnumi), segnumi) = TSS(KB(segnumi), segnumi)     &
                       & + CEMACUMPWRELEASED(segnumi)/dlt*SD_T(1)                 !Cm³/s      ! cb 5/22/15
                     VOLCEMA(jb) = VOLCEMA(jb) + CEMACUMPWRELEASE(segnumi)
                 endif
 
             enddo     !SegNumI
         enddo     !JB
     enddo     !JW
 
 
 
     end subroutine CEMASEDIMENTMODELW2
