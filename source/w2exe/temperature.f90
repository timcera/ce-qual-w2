!*==temperature.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     subroutine TEMPERATURE
 
     use MAIN
     use GLOBAL
     use NAMESC
     use GEOMC
     use LOGICC
     use SURFHE
     use KINETIC
     use SHADEC
     use EDDY
     use STRUCTURES
     use TRANS
     use TVDC
     use SELWC
     use GDAYC
     use SCREENC
     use TDGAS
     use RSTART
     use MACROPHYTEC
     use POROSITYC
     use ZOOPLANKTONC
     use CEMAVARS
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND), dimension(1000) :: bta1, gma1
     real :: rn1
     external RESTART_OUTPUT
!
!*** End of declarations rewritten by SPAG
!
 
 
     do jw = 1, nwb                   ! places a limit of 1000 vertical layers
         if(READ_EXTINCTION(jw))gamma(:, US(BS(jw)):DS(BE(jw))) = EXH2O(jw)   ! SW 1/28/13
         kt = KTWB(jw)
         if(.NOT.NO_HEAT(jw))then
             if(.NOT.READ_RADIATION(jw))call SHORT_WAVE_RADIATION(jday)
             if(TERM_BY_TERM(jw))then                                  ! SW 1/25/05
                 if(TAIR(jw)>=5.0)then
           !RANLW(JW) = 5.31D-13*(273.15D0+TAIR(JW))**6*(1.0D0+0.0017D0*CLOUD(JW)**2)*0.97D0
                     RANLW(jw) = 5.31D-13*(273.15D0 + TAIR(jw))                &
                               & **6*(1.0D0 + 0.0017D0*CLOUD(jw)*CLOUD(jw))    &
                               & *0.97D0                                                              ! SW 4/20/16 SPEED
                 else
           !RANLW(JW) = 5.62D-8*(273.15D0+TAIR(JW))**4*(1.D0-0.261D0*DEXP(-7.77D-4*TAIR(JW)**2))*(1.0D0+0.0017D0*CLOUD(JW)**2)*0.97D0
                     RANLW(jw) = 5.62D-8*(273.15D0 + TAIR(jw))                 &
                               & **4*(1.D0 - 0.261D0*DEXP( - 7.77D-4*TAIR(jw)  &
                               & *TAIR(jw)))                                   &
                               & *(1.0D0 + 0.0017D0*CLOUD(jw)*CLOUD(jw))*0.97D0                                                                       ! SW 4/20/16 SPEED
                 endif
             endif
         endif
         do jb = BS(jw), BE(jw)
             if(BR_INACTIVE(jb))cycle
             iu = CUS(jb)
             id = DS(jb)
 
!******      Heat exchange
 
             if(.NOT.NO_HEAT(jw))then
                 do i = iu, id
                     if(DYNAMIC_SHADE(i))call SHADING
 
!**********          Surface
 
                     if(.NOT.ICE(i))then
                         if(TERM_BY_TERM(jw))then
                             call SURFACE_TERMS(T2(kt, i))
                             RS(i) = SRON(jw)*SHADE(i)
                             RN(i) = RS(i) + RANLW(jw) - RB(i) - RE(i) - RC(i)
                             heatex = RN(i)/rhowcp*BI(kt, i)*DLX(i)
                         else
                             call EQUILIBRIUM_TEMPERATURE
                             heatex = (ET(i) - T2(kt, i))*CSHE(i)*BI(kt, i)    &
                                    & *DLX(i)
                         endif
                         TSS(kt, i) = TSS(kt, i) + heatex
                         TSSS(jb) = TSSS(jb) + heatex*dlt
                         sroout = (1.0D0 - BETA(jw))*(SRON(jw)*SHADE(i)/rhowcp)&
                                & *BI(kt, i)*DLX(i)                            &
                                & *DEXP( - gamma(kt, i)*DEPTHB(kt, i))
                         TSS(kt, i) = TSS(kt, i) - sroout
                         TSSS(jb) = TSSS(jb) - sroout*dlt
                         if(kt==KB(i))then
                                     ! SW 4/18/07
                             srosed = sroout*TSEDF(jw)
                         else
                             srosed = sroout*(1.0D0 - BI(kt + 1, i)/BI(kt, i)) &
                                    & *TSEDF(jw)
                         endif
                         TSS(kt, i) = TSS(kt, i) + srosed
                         TSSS(jb) = TSSS(jb) + srosed*dlt
                         sroin = sroout*B(kt + 1, i)/BI(kt, i)
                         do k = kt + 1, KB(i)
                             sroout = sroin*DEXP( - gamma(k, i)*(H1(k, i)))
                             sronet = sroin - sroout
                             if(k/=KB(i))then                              ! SW 1/18/08
                                 srosed = sroout*(1.0D0 - BI(k + 1, i)/BI(k, i)&
                                   & )*TSEDF(jw)
                             else
                                 srosed = sroout*TSEDF(jw)
                             endif
                             TSS(k, i) = TSS(k, i) + sronet + srosed
                             TSSS(jb) = TSSS(jb) + (sronet + srosed)*dlt
                             sroin = sroout*B(k + 1, i)/B(k, i)
                         enddo
                     endif
 
!**********          Sediment/water
 
                     do k = kt, KB(i)
                         if(k==KB(i))then     ! SW 4/18/07
                             tflux = CBHE(jw)/rhowcp*(TSED(jw) - T2(k, i))     &
                                   & *BI(k, i)*DLX(i)
                         else
                             tflux = CBHE(jw)/rhowcp*(TSED(jw) - T2(k, i))     &
                                   & *(BI(k, i) - BI(k + 1, i))*DLX(i)
                         endif
                         TSS(k, i) = TSS(k, i) + tflux
                         TSSB(jb) = TSSB(jb) + tflux*dlt
                     enddo
                 enddo
 
!********        Ice cover
 
                 if(ICE_CALC(jw))then
 !           HIA = 0.2367*CSHE(I)/5.65E-8    ! SW 10/20/09 Duplicate line of code
                     do i = iu, id
                         ALLOW_ICE(i) = .TRUE.
                         if(T2(kt, i)>ICET2(jw))ALLOW_ICE(i) = .FALSE.                ! RC/SW 4/28/11
 !             DO K=KT,KB(I)                                                          ! RC 4/28/11: no reason to loop over all layers
 !               IF (T2(K,I) > ICET2(JW)) ALLOW_ICE(I) = .FALSE.
 !             END DO
                     enddo
  !          ICE_IN(JB) = .TRUE.                                                      ! RC/SW 4/28/11 eliminate ICE_IN
  !          DO I=IU,ID
  !            IF (ICETH(I) < ICEMIN(JW)) ICE_IN(JB) = .FALSE.
  !          END DO
                     do i = iu, id
                         if(.NOT.(SALT_WATER(jw)))then
                                            ! SW/RC 4/28/11
                             rimt = 0.0
                         elseif(TDS(kt, i)<35.)then
                             rimt = -0.0545*TDS(kt, i)                           ! REGRESSION FOR TDS BETWEEN 0 AND 35 PPT
                         else
                             rimt = -0.31462 - 0.04177*TDS(kt, i)              &
                                  & - 0.000166*TDS(kt, i)*TDS(kt, i)             ! REGRESSION EQN FOR TDS>35 PPT
                         endif
                         if(DETAILED_ICE(jw))then
                             if(T2(kt, i)<0.0)then
                                 if(.NOT.ICE(i))then
                                     iceth2 = -T2(kt, i)*RHO(kt, i)            &
                                       & *cp*H2(kt, i)/rhoirl1
                                     if(iceth2<ice_tol)then
                                         iceth2 = 0.0D0
                                     else
                                         tflux = T2(kt, i)*RHO(kt, i)          &
                                           & *cp*H2(kt, i)*BI(kt, i)           &
                                           & /(rhowcp*dlt)*DLX(i)
                                         TSS(kt, i) = TSS(kt, i) - tflux
                                         TSSICE(jb) = TSSICE(jb) - tflux*dlt
                                     endif
                                 endif
                             endif
 
!**************              Ice balance
 
                             if(ICE(i))then
                                 tice = TAIR(jw)
                                 del = 2.0D0
                                 j = 1
                                 if(TAIR(jw)>=5.0)then
                                     RANLW(jw) = 5.31D-13*(273.15D0 + TAIR(jw))&
                                       & **6*(1.0D0 + 0.0017D0*CLOUD(jw)**2)   &
                                       & *0.97D0
                                 else
                                     RANLW(jw) = 5.62D-8*(273.15D0 + TAIR(jw)) &
                                       & **4*                                  &
                                       & (1.D0 - 0.261D0*DEXP( - 7.77D-4*TAIR  &
                                       & (jw)**2))                             &
                                       & *(1.0D0 + 0.0017D0*CLOUD(jw)**2)      &
                                       & *0.97D0
                                 endif
                                 rn1 = SRON(jw)/refl*SHADE(i)                  &
                                     & *(1.0D0 - ALBEDO(jw))*BETAI(jw)         &
                                     & + RANLW(jw)                                                  ! SW 4/19/10 eliminate spurious divsion of SRO by RHOCP
                                 do while (ABS(del)>1.0 .AND. j<500)                              ! SW 4/21/10 Should have been ABS of DEL
                                     call SURFACE_TERMS(tice)
                                     RN(i) = rn1 - RB(i) - RE(i) - RC(i)
                                                     ! 4/19/10
!                                    RN(I) =
! SRON(JW)/(REFL*RHOWCP)*SHADE(I)*(1.0-ALBEDO(JW))*BETAI(JW)+RANLW(JW)-RB(I)-RE
                                     del = RN(i) + rk1*(rimt - tice)/ICETH(i)
                                                               ! RK1 is ice conductivity 2.12 W/m/oC
                                     if(ABS(del)>1.0)tice = tice + del/500.0D0
                                     j = j + 1
                                 enddo
 
!****************                Solar radiation attenuation
 
                                 tflux = DLX(i)*SRON(jw)/(rhowcp*refl)*SHADE(i)&
                                   & *(1.0D0 - ALBEDO(jw))*(1.0D0 - BETAI(jw)) &
                                   & *DEXP( - GAMMAI(jw)*ICETH(i))*BI(kt, i)                                                           !   ! SW 4/21/10 Eliminate spurious divide by RHOCP
                                 TSS(kt, i) = TSS(kt, i) + tflux
                                 TSSICE(jb) = TSSICE(jb) + tflux*dlt
                                 if(tice>0.0)then
                                     hice = rhoicp*0.5D0*tice*0.5D0*ICETH(i)   &
                                       & *BI(kt, i)/(rhowcp*dlt)
                                     icethu = -dlt*hice/B(KTI(i), i)           &
                                       & *rhowcp/rhoirl1
                                     tice = 0.0D0
                                 endif
 
!****************                Ice growth
 
                                 if(tice<0.0)                                  &
                                  & iceth1 = dlt*(rk1*(rimt - tice)/ICETH(i))  &
                                  & /rhoirl1
 
!****************                Ice melt from water-ice interface
 
                                 if(T2(kt, i)>0.0)then
                                     iceth2 = -dlt*HWI(jw)*(T2(kt, i) - rimt)  &
                                       & /rhoirl1
                                     tflux = 2.392D-7*HWI(jw)                  &
                                       & *(rimt - T2(kt, i))*BI(kt, i)*DLX(i)
                                     TSS(kt, i) = TSS(kt, i) + tflux
                                     TSSICE(jb) = TSSICE(jb) + tflux*dlt
                                 endif
                             endif
 
!**************              Ice thickness
 
                             ICETH(i) = ICETH(i) + icethu + iceth1 + iceth2
                ! SW 9/4/15
                             if(ICEC(jw)=='    ONWB')then
                                 icethicknesschange = icethu + iceth1 + iceth2
                                 if(icethicknesschange>=0.0 .AND. ICETH(i)     &
                                  & >ice_tol)then
                                     ICEQSS(i) = -icethicknesschange*BI(kt, i) &
                                       & *DLX(i)*0.917/dlt                    ! CEMA: removal of 0.917 m3 of water for every 1 m3 of ice formed
                                     ICEBANK(i) = ICEBANK(i)                   &
                                       & + icethicknesschange*BI(kt, i)*DLX(i)
                                                                              ! Icebank - always + - volume of ice in m3
                                 elseif(ICETH(i)<ice_tol)then
                                     ICETH(i) = 0.0D0
                                     ICEQSS(i) = ICEBANK(i)*0.917/dlt
                                     ICEBANK(i) = 0.0
                                     if(i==iu .AND. US(jb)<iu)then
                                                           ! CHECK SUBTRACTED SEGMENTS THAT MAY HAVE ICE, MELT THEM AND PUT WATER IN SEGMENT IU
                                         do ii = US(jb), iu - 1
                                         ICETH(ii) = 0.0D0
                                         ICEQSS(iu) = ICEQSS(iu) + ICEBANK(ii) &
                                           & *0.917/dlt
                                         ICEBANK(ii) = 0.0
                                         ICE(ii) = .FALSE.
                                         enddo
                                     endif
 
                                 else
                                     ICEQSS(i) = -(icethicknesschange/ICETH(i))&
                                       & *ICEBANK(i)*0.917/dlt                            ! SW 9/29/15
                                     ICEBANK(i)                                &
                                       & = (1. + icethicknesschange/ICETH(i))  &
                                       & *ICEBANK(i)                                      ! Note: ICETHICKNESSCHANGE is negative
                                     if(i==iu .AND. US(jb)<iu)then
                                                           ! CHECK SUBTRACTED SEGMENTS THAT MAY HAVE ICE, MELT THEM AND PUT WATER IN SEGMENT IU
                                         do ii = US(jb), iu - 1
                                         ICEQSS(iu) = ICEQSS(iu)               &
                                           & - (icethicknesschange/ICETH(iu))  &
                                           & *ICEBANK(ii)*0.917/dlt
                                         ICEBANK(ii)                           &
                                           & = (1. + icethicknesschange/ICETH  &
                                           & (iu))*ICEBANK(ii)
                                         enddo
                                     endif
                                 endif
                                 VOLICE(jb) = VOLICE(jb) + ICEQSS(i)*dlt
                             endif
 
                             if(ICETH(i)<ice_tol)ICETH(i) = 0.0D0
     !           IF (WINTER .AND. (.NOT. ICE_IN(JB))) THEN            ! RC 4/28/11 No reason for this
     !             IF (.NOT. ALLOW_ICE(I)) ICETH(I) = 0.0
     !           END IF
                             ICE(i) = ICETH(i)>0.0
                             if(ICE(i))then
                                   ! 3/27/08 SW
                                 ICESW(i) = 0.0
                             else
                                 ICESW(i) = 1.0
                             endif
                             icethu = 0.0
                             iceth1 = 0.0
                             iceth2 = 0.0
                             if(ICETH(i)<ice_tol .AND. ICETH(i)>0.0)ICETH(i)   &
                              & = ice_tol
                         else                                              ! IF no ice the preceding time step
                             if(TERM_BY_TERM(jw))call EQUILIBRIUM_TEMPERATURE
                                                                           ! SW 10/20/09 Must call this first otherwise ET and CSHE are 0
                             hia = 0.2367D0*CSHE(i)/5.65D-8                  ! JM 11/08 convert SI units of m/s to English (btu/ft2/d/F) and then back to SI W/m2/C
!                            ICETH(I) =
! MAX(0.0,ICETH(I)+DLT*((RIMT-ET(I))/(ICETH(I)/RK1+1.0/HIA)-(T2(KT,I)-RIMT))/RH
                             ICETH(i) = MAX(0.0, ICETH(i) + dlt*((rimt - ET(i))&
                                      & /(ICETH(i)/rk1 + 1.0D0/hia) - HWI(jw)  &
                                      & *(T2(kt,i) - rimt))/rhoirl1)                                                       ! SW 10/20/09 Revised missing HWI(JW)
                             ICE(i) = ICETH(i)>0.0
                             ICESW(i) = 1.0
                             if(ICE(i))then
!                                TFLUX      =
!                                2.392E-7*(RIMT-T2(KT,I))*BI(KT,I)*DLX(I)
                                 tflux = 2.392D-7*HWI(jw)*(rimt - T2(kt, i))   &
                                   & *BI(kt, i)*DLX(i)                               ! SW 10/20/09 Revised missing HWI(JW)
                                 TSS(kt, i) = TSS(kt, i) + tflux
                                 TSSICE(jb) = TSSICE(jb) + tflux*dlt
                                 ICESW(i) = 0.0
                             endif
                         endif
                     enddo
                 endif
             endif
          !DO I=IU,ID
          !  icebank_all=icebank(i)+icebank_all
          !END DO
          !write(25600,*)jday,icebank_all
!******      Heat sources/sinks and total inflow/outflow
 
             if(EVAPORATION(jw))then
                 do i = iu, id
                     TSS(kt, i) = TSS(kt, i) - EV(i)*T2(kt, i)
                     TSSEV(jb) = TSSEV(jb) - EV(i)*T2(kt, i)*dlt
                     VOLEV(jb) = VOLEV(jb) - EV(i)*dlt
                 enddo
             endif
             if(PRECIPITATION(jw))then
                 do i = iu, id
                     TSS(kt, i) = TSS(kt, i) + QPR(i)*TPR(jb)
                     TSSPR(jb) = TSSPR(jb) + QPR(i)*TPR(jb)*dlt
                     VOLPR(jb) = VOLPR(jb) + QPR(i)*dlt
                 enddo
             endif
             if(tributaries)then
                 do jt = 1, jtt
                     if(jb==JBTR(jt))then
                         i = ITR(jt)
                         if(i<CUS(jb))i = CUS(jb)
                         do k = KTTR(jt), KBTR(jt)
                             if(QTR(jt)<0)then
                                 TSS(k, i) = TSS(k, i) + T2(k, i)*QTR(jt)      &
                                   & *QTRF(k, jt)
                                 TSSTR(jb) = TSSTR(jb) + T2(k, i)*QTR(jt)      &
                                   & *QTRF(k, jt)*dlt
                             else
                                 TSS(k, i) = TSS(k, i) + TTR(jt)*QTR(jt)       &
                                   & *QTRF(k, jt)
                                 TSSTR(jb) = TSSTR(jb) + TTR(jt)*QTR(jt)       &
                                   & *QTRF(k, jt)*dlt
                             endif
                         enddo
                         VOLTRB(jb) = VOLTRB(jb) + QTR(jt)*dlt
                     endif
                 enddo
             endif
             if(DIST_TRIBS(jb))then
                 do i = iu, id
                     if(QDT(i)<0)then
                         TSS(kt, i) = TSS(kt, i) + T2(kt, i)*QDT(i)
                         TSSDT(jb) = TSSDT(jb) + T2(kt, i)*QDT(i)*dlt
                     else
                         TSS(kt, i) = TSS(kt, i) + TDTR(jb)*QDT(i)
                         TSSDT(jb) = TSSDT(jb) + TDTR(jb)*QDT(i)*dlt
                     endif
                     VOLDT(jb) = VOLDT(jb) + QDT(i)*dlt
                 enddo
             endif
             if(withdrawals)then
                 do jwd = 1, jww
                     if(QWD(jwd)/=0.0)then
                         if(jb==JBWD(jwd))then
                             i = MAX(CUS(JBWD(jwd)), IWD(jwd))
                             do k = KTW(jwd), KBW(jwd)
                                 TSS(k, i) = TSS(k, i) - T2(k, i)*QSW(k, jwd)
                                 TSSWD(jb) = TSSWD(jb) - T2(k, i)*QSW(k, jwd)  &
                                   & *dlt
                             enddo
                             VOLWD(jb) = VOLWD(jb) - QWD(jwd)*dlt
                         endif
                     endif
                 enddo
             endif
             if(UP_FLOW(jb))then
                 do k = kt, KB(iu)
                     if(.NOT.HEAD_FLOW(jb))then
                         TSS(k, iu) = TSS(k, iu) + QINF(k, jb)*QIN(jb)*TIN(jb)
                         TSSIN(jb) = TSSIN(jb) + QINF(k, jb)*QIN(jb)*TIN(jb)   &
                                   & *dlt
                     elseif(U(k, iu - 1)>=0.0)then
                         TSS(k, iu) = TSS(k, iu) + U(k, iu - 1)*BHR1(k, iu - 1)&
                                    & *T1(k, iu - 1)
                         TSSIN(jb) = TSSIN(jb) + U(k, iu - 1)*BHR1(k, iu - 1)  &
                                   & *T1(k, iu - 1)*dlt
                     else
                         TSS(k, iu) = TSS(k, iu) + U(k, iu - 1)*BHR1(k, iu - 1)&
                                    & *T1(k, iu)
                         TSSIN(jb) = TSSIN(jb) + U(k, iu - 1)*BHR1(k, iu - 1)  &
                                   & *T1(k, iu)*dlt
                     endif
                 enddo
                 VOLIN(jb) = VOLIN(jb) + QIN(jb)*dlt
             endif
             if(DN_FLOW(jb))then
                 do k = kt, KB(id)
                     TSS(k, id) = TSS(k, id) - QOUT(k, jb)*T2(k, id + 1)
                     TSSOUT(jb) = TSSOUT(jb) - QOUT(k, jb)*T2(k, id + 1)*dlt
                     VOLOUT(jb) = VOLOUT(jb) - QOUT(k, jb)*dlt
                 enddo
             endif
             if(UP_HEAD(jb))then
                 do k = kt, KB(iu)
                     iut = iu
                     if(QUH1(k, jb)>=0.0)iut = iu - 1
                     TSSUH1(k, jb) = T2(k, iut)*QUH1(k, jb)
                     TSS(k, iu) = TSS(k, iu) + TSSUH1(k, jb)
                     TSSUH(jb) = TSSUH(jb) + TSSUH1(k, jb)*dlt
                     VOLUH(jb) = VOLUH(jb) + QUH1(k, jb)*dlt
                 enddo
             endif
             if(UH_INTERNAL(jb))then
                 if(UHS(jb)/=DS(JBUH(jb)) .OR. DHS(JBUH(jb))/=US(jb))then
                     if(JBUH(jb)>=BS(jw) .AND. JBUH(jb)<=BE(jw))then
                         do k = kt, KB(iu - 1)
                             TSS(k, UHS(jb)) = TSS(k, UHS(jb)) - TSSUH2(k, jb) &
                               & /dlt
                             TSSUH(JBUH(jb)) = TSSUH(JBUH(jb)) - TSSUH2(k, jb)
                             VOLUH(JBUH(jb)) = VOLUH(JBUH(jb)) - VOLUH2(k, jb)
                         enddo
                     else
                         call UPSTREAM_CONSTITUENT(T2, TSS)
                         do k = kt, KB(iu - 1)
                             TSSUH(JBUH(jb)) = TSSUH(JBUH(jb)) - TSSUH2(k, jb)
                             VOLUH(JBUH(jb)) = VOLUH(JBUH(jb)) - VOLUH2(k, jb)
                         enddo
                     endif
                 endif
             endif
             if(DN_HEAD(jb))then
                 do k = kt, KB(id + 1)
                     idt = id + 1
                     if(QDH1(k, jb)>=0.0)idt = id
                     TSSDH1(k, jb) = T2(k, idt)*QDH1(k, jb)
                     TSS(k, id) = TSS(k, id) - TSSDH1(k, jb)
                     TSSDH(jb) = TSSDH(jb) - TSSDH1(k, jb)*dlt
                     VOLDH(jb) = VOLDH(jb) - QDH1(k, jb)*dlt
                 enddo
             endif
             if(DH_INTERNAL(jb))then
                 if(DHS(jb)/=US(JBDH(jb)) .OR. UHS(JBDH(jb))/=DS(jb))then
                     if(JBDH(jb)>=BS(jw) .AND. JBDH(jb)<=BE(jw))then
                         do k = kt, KB(id + 1)
                             TSS(k, CDHS(jb)) = TSS(k, CDHS(jb))               &
                               & + TSSDH2(k, jb)/dlt
                             TSSDH(JBDH(jb)) = TSSDH(JBDH(jb)) + TSSDH2(k, jb)
                             VOLDH(JBDH(jb)) = VOLDH(JBDH(jb)) + VOLDH2(k, jb)
                         enddo
                     else
                         call DOWNSTREAM_CONSTITUENT(T2, TSS)
                         do k = kt, KB(id + 1)
                             TSSDH(JBDH(jb)) = TSSDH(JBDH(jb)) + TSSDH2(k, jb)
                             VOLDH(JBDH(jb)) = VOLDH(JBDH(jb)) + VOLDH2(k, jb)
                         enddo
                     endif
                 endif
             endif
         enddo
     enddo
 
!**  Temperature transport
 
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             if(BR_INACTIVE(jb))cycle
             iu = CUS(jb)
             id = DS(jb)
             cold => hyd(:, :, 4)
             call HORIZONTAL_MULTIPLIERS1
             call VERTICAL_MULTIPLIERS1
             call HORIZONTAL_MULTIPLIERS
             call VERTICAL_MULTIPLIERS
     !   CNEW => T1(:,:)
     !   SSB  => TSS(:,:)
     !   SSK  => CSSB(:,:,1)
     !   CALL HORIZONTAL_TRANSPORT
 
             do i = iu, id
                      !CONCURRENT(I=IU:ID)   !
                 do k = kt, KB(i)
                          !CONCURRENT(K=KT:KB(I))      !FORALL                                         !DO K=KT,KB(I)
                     DT(k, i) = (cold(k, i)*BH2(k, i)/dlt + (ADX(k, i)*BHR1(k, &
                              & i) - ADX(k, i - 1)*BHR1(k, i - 1))/DLX(i)      &
                              & + (1.0D0 - THETA(jw))                          &
                              & *(ADZ(k, i)*BB(k, i) - ADZ(k - 1, i)           &
                              & *BB(k - 1, i)) + TSS(k, i)/DLX(i))             &
                              & *dlt/BH1(k, i)
                 enddo
             enddo
 
 
             do i = iu, id
                     !CONCURRENT(I=IU:ID)   !
                 do k = kt, KB(i)
                     AT(k, i) = 0.0D0
                     CT(k, i) = 0.0D0
                     VT(k, i) = 0.0D0                            !; DT(:,I) = 0.0D0    SW CODE SPEEDUP 6/15/13
                 enddo
                 do k = kt, KB(i)
                          !CONCURRENT(K=KT:KB(I))          !FORALL(K=KT:KB(I))                                                 !DO K=KT,KB(I)
                     AT(k, i) = -dlt/BH1(k, i)                                 &
                              & *(BB(k - 1, i)*(DZ(k - 1, i)/AVH1(k - 1, i)    &
                              & + THETA(jw)*0.5D0*W(k - 1, i)))
                     CT(k, i) = dlt/BH1(k, i)                                  &
                              & *(BB(k, i)*(THETA(jw)*0.5D0*W(k, i) - DZ(k, i) &
                              & /AVH1(k, i)))
                     VT(k, i) = 1.0D0 + dlt/BH1(k, i)                          &
                              & *(BB(k, i)*(DZ(k, i)/AVH1(k, i) + THETA(jw)    &
                              & *0.5D0*W(k, i)) + BB(k - 1, i)                 &
                              & *(DZ(k - 1, i)/AVH1(k - 1, i) - THETA(jw)      &
                              & *0.5D0*W(k - 1, i)))
       !     DT(K,I) =  CNEW(K,I)
                 enddo
         ! CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,CNEW(:,I))
                 bta1(kt) = VT(kt, i)
                 gma1(kt) = DT(kt, i)
                 do k = kt + 1, KB(i)
                     bta1(k) = VT(k, i) - AT(k, i)/bta1(k - 1)*CT(k - 1, i)
                     gma1(k) = DT(k, i) - AT(k, i)/bta1(k - 1)*gma1(k - 1)
                 enddo
                 T1(KB(i), i) = gma1(KB(i))/bta1(KB(i))
                 do k = KB(i) - 1, kt, -1
                     T1(k, i) = (gma1(k) - CT(k, i)*T1(k + 1, i))/bta1(k)
                 enddo
             enddo
         enddo
     enddo
 
 
     end subroutine TEMPERATURE
