!*==wqconstituents.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     subroutine WQCONSTITUENTS
 
 
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
     use TRIDIAG_V
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: tnalg, tnbod, tnz, tpalg, tpbod, tpz
     external RESTART_OUTPUT
!
!*** End of declarations rewritten by SPAG
!
 ! REAL(R8) :: BTA1(1000),GMA1(1000)   ! places a limit of 1000 vertical layers
 
     if(macrophyte_on .AND. update_kinetics)call POROSITY
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             if(BR_INACTIVE(jb))cycle
                                    ! SW 6/12/2017
             iu = CUS(jb)
             id = DS(jb)
 
!********    Kinetic sources/sinks
 
             if(SEDIMENT_CALC(jw))then
                 call SEDIMENT
                 call SEDIMENTP
                 call SEDIMENTN
                 call SEDIMENTC
                 if(DYNSEDK(jw)=='      ON')call SEDIMENT_DECAY_RATE
             endif
!            Amaila start
             if(sedcomp_exist)then
                               ! SW 5/26/15
                 if(SEDIMENT_CALC1(jw))call SEDIMENT1
                 if(SEDIMENT_CALC2(jw))call SEDIMENT2
             endif
!            Amaila end
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))call MACROPHYTE(m)
             enddo
 
             if(update_kinetics)then
                 if(update_rates)then
                     call TEMPERATURE_RATES
                     call KINETIC_RATES
                 endif
                 do jac = 1, nac
                     jc = CN(jac)
                     if(jc==npo4)call PHOSPHORUS
                     if(jc==nnh4)call AMMONIUM
                     if(jc==nno3)call NITRATE
                     if(jc==ndsi)call DISSOLVED_SILICA
                     if(jc==npsi)call PARTICULATE_SILICA
                     if(jc==nfe)call IRON
                     if(jc==nldom)call LABILE_DOM
                     if(jc==nrdom)call REFRACTORY_DOM
                     if(jc==nlpom)call LABILE_POM
                     if(jc==nrpom)call REFRACTORY_POM
                     if(jc==ndo)call DISSOLVED_OXYGEN
                     if(jc>=ngcs .AND. jc<=ngce)                               &
                      & call GENERIC_CONST(jc - ngcs + 1)
                     if(jc>=nsss .AND. jc<=nsse)                               &
                      & call SUSPENDED_SOLIDS(jc - nsss + 1)
                     if(jc>=nas .AND. jc<=nae)then
                         if(ALG_CALC(jc - nas + 1))call ALGAE(jc - nas + 1)
                     endif
                     if(jc>=nbods .AND. jc<=nbode)then
                         do jcb = 1, nbod
                                    ! VARIABLE STOICHIOMETRY FOR CBOD, CB 6/6/10
                             if(BOD_CALC(jcb))then
                                 if(jc==NBODC(jcb))                            &
                                  & call BIOCHEMICAL_O2_DEMAND(jcb)
                                 if(jc==NBODP(jcb) .AND. BOD_CALCP(jcb))       &
                                  & call BIOCHEMICAL_O2_DEMAND_P(jcb)                                  ! CB 5/19/2011
                                 if(jc==NBODN(jcb) .AND. BOD_CALCN(jcb))       &
                                  & call BIOCHEMICAL_O2_DEMAND_N(jcb)                                  ! CB 5/19/2011
                             endif
                         enddo
                     endif
                     if(jc>=nzoos .AND. jc<=nzooe .AND. zooplankton_calc)      &
                      & call ZOOPLANKTON
                     if(jc==nldomp)call LABILE_DOM_P
                     if(jc==nrdomp)call REFRACTORY_DOM_P
                     if(jc==nlpomp)call LABILE_POM_P
                     if(jc==nrpomp)call REFRACTORY_POM_P
                     if(jc==nldomn)call LABILE_DOM_N
                     if(jc==nrdomn)call REFRACTORY_DOM_N
                     if(jc==nlpomn)call LABILE_POM_N
                     if(jc==nrpomn)call REFRACTORY_POM_N
              !IF (JC == NALK .and. NONCON_ALKALINITY)                  CALL alkalinity
                     if(jc==nalk)call ALKALINITY                  ! NW 2/11/16
                 enddo
                 if(PH_CALC(jw))call INORGANIC_CARBON
                 if(PH_CALC(jw))then
                     if(ph_buffering)then
                                    ! enhanced pH buffering
                         call PH_CO2_NEW
                     else
                         call PH_CO2
                     endif
                 endif
             endif
             do je = 1, nep
                        ! sw 5/16/06
                 if(EPIPHYTON_CALC(jw, je))call EPIPHYTON(je)
             enddo
 
!********    External sources/sinks
 
             if(aeratec=="      ON")call AERATEMASS
             if(EVAPORATION(jw) .AND. water_age_active)then   ! CORRECT WATER AGE FOR EVAPORATION SR 7/27/2017
                 do i = iu, id
                     jc = ngcs + jg_age - 1
                     CSSB(kt, i, jc) = CSSB(kt, i, jc) - EV(i)*CG(kt, i, jc)
                 enddo
             endif
 
             do jac = 1, nac
                 jc = CN(jac)
                 if(tributaries)then
                     do jt = 1, jtt
                         if(jb==JBTR(jt))then
                             i = ITR(jt)
                             if(i<CUS(jb))i = CUS(jb)
                             do k = KTTR(jt), KBTR(jt)
                                 if(QTR(jt)<0.0)then
                                     CSSB(k, i, jc) = CSSB(k, i, jc)           &
                                       & + C1(k, i, jc)*QTR(jt)*QTRF(k, jt)
                                 else
                                     CSSB(k, i, jc) = CSSB(k, i, jc)           &
                                       & + CTR(jc, jt)*QTR(jt)*QTRF(k, jt)
                                 endif
                             enddo
                         endif
                     enddo
                 endif
                 if(DIST_TRIBS(jb))then
                     do i = iu, id
                         if(QDT(i)<0.0)then
                             CSSB(kt, i, jc) = CSSB(kt, i, jc) + C1(kt, i, jc) &
                               & *QDT(i)
                         else
                             CSSB(kt, i, jc) = CSSB(kt, i, jc) + CDTR(jc, jb)  &
                               & *QDT(i)
                         endif
                     enddo
                 endif
                 if(withdrawals)then
                     do jwd = 1, jww
                         if(QWD(jwd)/=0.0)then
                             if(jb==JBWD(jwd))then
                                 i = MAX(CUS(JBWD(jwd)), IWD(jwd))
                                 do k = KTW(jwd), KBW(jwd)
                                                !CONCURRENT(K=KTW(JWD):KBW(JWD))                 ! FORALL
                                     CSSB(k, i, jc) = CSSB(k, i, jc)           &
                                       & - C1S(k, i, jc)*QSW(k, jwd)
                                 enddo
                             endif
                         endif
                     enddo
                 endif
                 if(PRECIPITATION(jw))then
                     do i = iu, id
                            !CONCURRENT (I=IU:ID)                                  !FORALL
                         CSSB(kt, i, jc) = CSSB(kt, i, jc) + CPR(jc, jb)*QPR(i)
                     enddo
                 endif
                 if(UP_FLOW(jb))then
                     do k = kt, KB(iu)
                         if(.NOT.HEAD_FLOW(jb))then
                             CSSB(k, iu, jc) = CSSB(k, iu, jc) + QINF(k, jb)   &
                               & *QIN(jb)*CIN(jc, jb)
                         elseif(U(k, iu - 1)>=0.0)then
                             CSSB(k, iu, jc) = CSSB(k, iu, jc) + U(k, iu - 1)  &
                               & *BHR1(k, iu - 1)*C1S(k, iu - 1, jc)
                         else
                             CSSB(k, iu, jc) = CSSB(k, iu, jc) + U(k, iu - 1)  &
                               & *BHR1(k, iu - 1)*C1S(k, iu, jc)
                         endif
                     enddo
                 endif
                 if(DN_FLOW(jb))CSSB(kt:KB(id), id, jc)                        &
                  & = CSSB(kt:KB(id), id, jc) - qout(kt:KB(id), jb)            &
                  & *C1S(kt:KB(id), id, jc)
                 if(UP_HEAD(jb))then
                     do k = kt, KB(iu)
                         iut = iu
                         if(QUH1(k, jb)>=0.0)iut = iu - 1
                         CSSUH1(k, jc, jb) = C1S(k, iut, jc)*QUH1(k, jb)
                         CSSB(k, iu, jc) = CSSB(k, iu, jc) + CSSUH1(k, jc, jb)
                     enddo
                     if(UH_INTERNAL(jb))then
                         if(UHS(jb)/=DS(JBUH(jb)) .OR. DHS(JBUH(jb))/=US(jb))  &
                          & then
                             if(JBUH(jb)>=BS(jw) .AND. JBUH(jb)<=BE(jw))then
                                 i = UHS(jb)
                    !dir$ ivdep
                                 do k = kt, KB(iu)
                                     CSSB(k, i, jc) = CSSB(k, i, jc)           &
                                       & - CSSUH2(k, jc, jb)/dlt
                                 enddo
                             else
                                 call UPSTREAM_CONSTITUENT(c2(:, :, jc),       &
                                   & CSSB(:, :, jc))
                             endif
                         endif
                     endif
                 endif
                 if(DN_HEAD(jb))then
                     do k = kt, KB(id + 1)
                         idt = id + 1
                         if(QDH1(k, jb)>=0.0)idt = id
                         CSSDH1(k, jc, jb) = C1S(k, idt, jc)*QDH1(k, jb)
                         CSSB(k, id, jc) = CSSB(k, id, jc) - CSSDH1(k, jc, jb)
                     enddo
                     if(DH_INTERNAL(jb))then
                         if(DHS(jb)/=US(JBDH(jb)) .OR. UHS(JBDH(jb))/=DS(jb))  &
                          & then
                             if(JBDH(jb)>=BS(jw) .AND. JBDH(jb)<=BE(jw))then
                                 i = DHS(jb)
                                 do k = kt, KB(id + 1)
                                     CSSB(k, i, jc) = CSSB(k, i, jc)           &
                                       & + CSSDH2(k, jc, jb)/dlt
                                 enddo
                             else
                                 call DOWNSTREAM_CONSTITUENT(c2(:, :, jc),     &
                                   & CSSB(:, :, jc))
                             endif
                         endif
                     endif
                 endif
             enddo
 
             if(MASS_BALANCE(jw) .AND. CONTOUR(jw) .AND. derived_calc)then
                                                                       ! TO COMPUTE TP AND TN INFLOWS AND OUTFLOWS FOR MASSBAL.OPT FILE
                 if(tributaries)then
                     do jt = 1, jtt
                         if(jb==JBTR(jt))then
                             i = ITR(jt)
                             if(i<CUS(jb))i = CUS(jb)
                             do k = KTTR(jt), KBTR(jt)
                                 tpalg = 0.0
                                 tnalg = 0.0
                                 tpz = 0.0
                                 tnz = 0.0
                                 tpbod = 0.0
                                 tnbod = 0.0
                                 if(QTR(jt)<0.0)then
                                     do ja = 1, nal
                                         tpalg = tpalg + ALG(k, i, ja)*AP(ja)
                                         tnalg = tnalg + ALG(k, i, ja)*AN(ja)
                                     enddo
                                     do jcb = 1, nbod
                                         tpbod = tpbod + CBODP(k, i, jcb)
                                         tnbod = tnbod + CBODN(k, i, jcb)
                                     enddo
                                     do jz = 1, nzp
                                         tpz = tpz + ZOO(k, i, jz)*ZP(jz)
                                         tnz = tnz + ZOO(k, i, jz)*ZN(jz)
                                     enddo
                                     tpout = tpout -                           &
                                       & (tpalg + tpbod + tpz + PO4(k, i)      &
                                       & + LDOMP(k, i) + LPOMP(k, i)           &
                                       & + RDOMP(k, i) + RPOMP(k, i))*QTR(jt)  &
                                       & *QTRF(k, jt)*dlt/1000.
                                     tnout = tnout -                           &
                                       & (tnalg + tnbod + tnz + NO3(k, i)      &
                                       & + NH4(k, i) + LDOMN(k, i)             &
                                       & + LPOMN(k, i) + RDOMN(k, i)           &
                                       & + RPOMN(k, i))*QTR(jt)*QTRF(k, jt)    &
                                       & *dlt/1000.
                                 else
                                     do jc = nas, nae
                                         tpalg = tpalg + CTR(jc, jt)           &
                                           & *AP(jc - nas + 1)
                                         tnalg = tnalg + CTR(jc, jt)           &
                                           & *AN(jc - nas + 1)
                                     enddo
                                     do jc = nbods, nbode, 3
                                         tpbod = tpbod + CTR(jc + 1, jt)
                                         tnbod = tnbod + CTR(jc + 2, jt)
                                     enddo
                                     do jc = nzoos, nzooe
                                         tpz = tpz + CTR(jc, jt)               &
                                           & *ZP(jc - nzoos + 1)
                                         tnz = tnz + CTR(jc, jt)               &
                                           & *ZN(jc - nzoos + 1)
                                     enddo
                                     TPTRIB(jw) = TPTRIB(jw)                   &
                                       & + (tpalg + tpbod + tpz + CTR(npo4, jt)&
                                       & + CTR(nldomp, jt) + CTR(nrdomp, jt)   &
                                       & + CTR(nlpomp, jt) + CTR(nrpomp, jt))  &
                                       & *QTR(jt)*QTRF(k, jt)*dlt/1000.
                                     TNTRIB(jw) = TNTRIB(jw)                   &
                                       & + (tnalg + tnbod + tnz + CTR(nnh4, jt)&
                                       & + CTR(nno3, jt) + CTR(nldomn, jt)     &
                                       & + CTR(nrdomn, jt) + CTR(nlpomn, jt)   &
                                       & + CTR(nrpomn, jt))*QTR(jt)*QTRF(k, jt)&
                                       & *dlt/1000.
                                 endif
                             enddo
                         endif
                     enddo
                 endif
                 if(DIST_TRIBS(jb))then
                     do i = iu, id
                         tpalg = 0.0
                         tnalg = 0.0
                         tpz = 0.0
                         tnz = 0.0
                         tpbod = 0.0
                         tnbod = 0.0
                         if(QDT(i)<0.0)then
                             do ja = 1, nal
                                 tpalg = tpalg + ALG(kt, i, ja)*AP(ja)
                                 tnalg = tnalg + ALG(kt, i, ja)*AN(ja)
                             enddo
                             do jcb = 1, nbod
                                 tpbod = tpbod + CBODP(kt, i, jcb)
                                 tnbod = tnbod + CBODN(kt, i, jcb)
                             enddo
                             do jz = 1, nzp
                                 tpz = tpz + ZOO(kt, i, jz)*ZP(jz)
                                 tnz = tnz + ZOO(kt, i, jz)*ZN(jz)
                             enddo
                             tpout = tpout - (tpalg + tpbod + tpz + PO4(kt, i) &
                                   & + LDOMP(kt, i) + LPOMP(kt, i)             &
                                   & + RDOMP(kt, i) + RPOMP(kt, i))*QDT(i)     &
                                   & *dlt/1000.
                             tnout = tnout - (tnalg + tnbod + tnz + NO3(kt, i) &
                                   & + NH4(kt, i) + LDOMN(kt, i) + LPOMN(kt, i)&
                                   & + RDOMN(kt, i) + RPOMN(kt, i))*QDT(i)     &
                                   & *dlt/1000.
                         else
                             do jc = nas, nae
                                 tpalg = tpalg + CDTR(jc, jb)*AP(jc - nas + 1)
                                 tnalg = tnalg + CDTR(jc, jb)*AN(jc - nas + 1)
                             enddo
                             do jc = nbods, nbode, 3
                                 tpbod = tpbod + CDTR(jc + 1, jb)
                                 tnbod = tnbod + CDTR(jc + 2, jb)
                             enddo
                             do jc = nzoos, nzooe
                                 tpz = tpz + CDTR(jc, jb)*ZP(jc - nzoos + 1)
                                 tnz = tnz + CDTR(jc, jb)*ZN(jc - nzoos + 1)
                             enddo
                             TPDTRIB(jw) = TPDTRIB(jw)                         &
                               & + (tpalg + tpbod + tpz + CDTR(npo4, jb)       &
                               & + CDTR(nldomp, jb) + CDTR(nrdomp, jb)         &
                               & + CDTR(nlpomp, jb) + CDTR(nrpomp, jb))*QDT(i) &
                               & *dlt/1000.
                             TNDTRIB(jw) = TNDTRIB(jw)                         &
                               & + (tnalg + tnbod + tnz + CDTR(nnh4, jb)       &
                               & + CDTR(nno3, jb) + CDTR(nldomn, jb)           &
                               & + CDTR(nrdomn, jb) + CDTR(nlpomn, jb)         &
                               & + CDTR(nrpomn, jb))*QDT(i)*dlt/1000.
                         endif
                     enddo
                 endif
                 if(withdrawals)then
                     do jwd = 1, jww
                         if(QWD(jwd)/=0.0)then
                             if(jb==JBWD(jwd))then
                                 i = MAX(CUS(JBWD(jwd)), IWD(jwd))
                                 do k = KTW(jwd), KBW(jwd)
                                                !CONCURRENT(K=KTW(JWD):KBW(JWD))                 ! FORALL
                                     tpalg = 0.0
                                     tnalg = 0.0
                                     tpz = 0.0
                                     tnz = 0.0
                                     tpbod = 0.0
                                     tnbod = 0.0
                                     do jc = nas, nae
                                         tpalg = tpalg + C1S(k, i, jc)         &
                                           & *AP(jc - nas + 1)
                                         tnalg = tnalg + C1S(k, i, jc)         &
                                           & *AN(jc - nas + 1)
                                     enddo
                                     do jc = nbods, nbode, 3
                                         tpbod = tpbod + C1S(k, i, jc)
                                         tnbod = tnbod + C1S(k, i, jc)
                                     enddo
                                     do jc = nzoos, nzooe
                                         tpz = tpz + C1S(k, i, jc)             &
                                           & *ZP(jc - nzoos + 1)
                                         tnz = tnz + C1S(k, i, jc)             &
                                           & *ZN(jc - nzoos + 1)
                                     enddo
                                     TPWD(jw) = TPWD(jw)                       &
                                       & + (tpalg + tpbod + tpz +              &
                                       & C1S(k, i, npo4) + C1S(k, i, nldomp)   &
                                       & + C1S(k, i, nrdomp)                   &
                                       & + C1S(k, i, nlpomp)                   &
                                       & + C1S(k, i, nrpomp))*QSW(k, jwd)      &
                                       & *dlt/1000.
                                     TNWD(jw) = TNWD(jw)                       &
                                       & + (tnalg + tnbod + tnz +              &
                                       & C1S(k, i, nnh4) + C1S(k, i, nno3)     &
                                       & + C1S(k, i, nldomn)                   &
                                       & + C1S(k, i, nrdomn)                   &
                                       & + C1S(k, i, nlpomn)                   &
                                       & + C1S(k, i, nrpomn))*QSW(k, jwd)      &
                                       & *dlt/1000.
                                 enddo
                             endif
                         endif
                     enddo
                 endif
                 if(PRECIPITATION(jw))then
                     do i = iu, id
                            !CONCURRENT (I=IU:ID)                                  !FORALL
                         tpalg = 0.0
                         tnalg = 0.0
                         tpz = 0.0
                         tnz = 0.0
                         tpbod = 0.0
                         tnbod = 0.0
                         do jc = nas, nae
                             tpalg = tpalg + CPR(jc, jb)*AP(jc - nas + 1)
                             tnalg = tnalg + CPR(jc, jb)*AN(jc - nas + 1)
                         enddo
                         do jc = nbods, nbode, 3
                             tpbod = tpbod + CPR(jc, jb)
                             tnbod = tnbod + CPR(jc, jb)
                         enddo
                         do jc = nzoos, nzooe
                             tpz = tpz + CPR(jc, jb)*ZP(jc - nzoos + 1)
                             tnz = tnz + CPR(jc, jb)*ZN(jc - nzoos + 1)
                         enddo
                         TPPR(jw) = TPPR(jw)                                   &
                                  & + (tpalg + tpbod + tpz + CPR(npo4, jb)     &
                                  & + CPR(nldomp, jb) + CPR(nrdomp, jb)        &
                                  & + CPR(nlpomp, jb) + CPR(nrpomp, jb))*QPR(i)&
                                  & *dlt/1000.
                         TNPR(jw) = TNPR(jw)                                   &
                                  & + (tnalg + tnbod + tnz + CPR(nnh4, jb)     &
                                  & + CPR(nno3, jb) + CPR(nldomn, jb)          &
                                  & + CPR(nrdomn, jb) + CPR(nlpomn, jb)        &
                                  & + CPR(nrpomn, jb))*QPR(i)*dlt/1000.
                     enddo
                 endif
                 if(UP_FLOW(jb))then
                     do k = kt, KB(iu)
                         if(.NOT.HEAD_FLOW(jb))then
                             do jc = nas, nae
                                 tpalg = tpalg + CIN(jc, jb)*AP(jc - nas + 1)
                                 tnalg = tnalg + CIN(jc, jb)*AN(jc - nas + 1)
                             enddo
                             do jc = nbods, nbode, 3
                                 tpbod = tpbod + CIN(jc, jb)
                                 tnbod = tnbod + CIN(jc, jb)
                             enddo
                             do jc = nzoos, nzooe
                                 tpz = tpz + CIN(jc, jb)*ZP(jc - nzoos + 1)
                                 tnz = tnz + CIN(jc, jb)*ZN(jc - nzoos + 1)
                             enddo
                             TPIN(jw) = TPIN(jw)                               &
                                      & + (tpalg + tpbod + tpz + CIN(npo4, jb) &
                                      & + CIN(nldomp, jb) + CIN(nrdomp, jb)    &
                                      & + CIN(nlpomp, jb) + CIN(nrpomp, jb))   &
                                      & *QINF(k, jb)*QIN(jb)*dlt/1000.
                             TNIN(jw) = TNIN(jw)                               &
                                      & + (tnalg + tnbod + tnz + CIN(nnh4, jb) &
                                      & + CIN(nno3, jb) + CIN(nldomn, jb)      &
                                      & + CIN(nrdomn, jb) + CIN(nlpomn, jb)    &
                                      & + CIN(nrpomn, jb))*QINF(k, jb)*QIN(jb) &
                                      & *dlt/1000.
 
                         elseif(U(k, iu - 1)>=0.0)then
                    !CSSB(K,IU,JC) = CSSB(K,IU,JC)+U(K,IU-1)*BHR1(K,IU-1)*C1S(K,IU-1,JC)
                    !CSSB(K,IU,JC) = CSSB(K,IU,JC)+U(K,IU-1)*BHR1(K,IU-1)*C1S(K,IU,JC)
                         endif
                     enddo
                 endif
                 if(DN_FLOW(jb))then
                     do k = kt, KB(id)
                         tpalg = 0.0
                         tnalg = 0.0
                         tpz = 0.0
                         tnz = 0.0
                         tpbod = 0.0
                         tnbod = 0.0
                         do jc = nas, nae
                             tpalg = tpalg + C1S(k, id, jc)*AP(jc - nas + 1)
                             tnalg = tnalg + C1S(k, id, jc)*AN(jc - nas + 1)
                         enddo
                         do jc = nbods, nbode, 3
                             tpbod = tpbod + C1S(k, id, jc)
                             tnbod = tnbod + C1S(k, id, jc)
                         enddo
                         do jc = nzoos, nzooe
                             tpz = tpz + C1S(k, id, jc)*ZP(jc - nzoos + 1)
                             tnz = tnz + C1S(k, id, jc)*ZN(jc - nzoos + 1)
                         enddo
                         tpout(jw) = tpout(jw)                                 &
                                   & + (tpalg + tpbod + tpz + C1S(k, id, npo4) &
                                   & + C1S(k, id, nldomp) + C1S(k, id, nrdomp) &
                                   & + C1S(k, id, nlpomp) + C1S(k, id, nrpomp))&
                                   & *qout(k, jb)*dlt/1000.                                                                                                               ! C1S(KT:KB(ID),ID,JC)
                         tnout(jw) = tnout(jw)                                 &
                                   & + (tnalg + tnbod + tnz + C1S(k, id, nnh4) &
                                   & + C1S(k, id, nno3) + C1S(k, id, nldomn)   &
                                   & + C1S(k, id, nrdomn) + C1S(k, id, nlpomn) &
                                   & + C1S(k, id, nrpomn))*qout(k, jb)         &
                                   & *dlt/1000.
                     enddo
                 endif
 
                 if(UP_HEAD(jb))then
                     do k = kt, KB(iu)
                         iut = iu
                         if(QUH1(k, jb)>=0.0)iut = iu - 1
                !CSSUH1(K,JC,JB) = C1S(K,IUT,JC)*QUH1(K,JB)
 
                     enddo
                     if(UH_INTERNAL(jb))then
                         if(UHS(jb)/=DS(JBUH(jb)) .OR. DHS(JBUH(jb))/=US(jb))  &
                          & then
                             if(JBUH(jb)>=BS(jw) .AND. JBUH(jb)<=BE(jw))then
                                 i = UHS(jb)
                                 do k = kt, KB(iu)
                      !CSSB(K,I,JC) = CSSB(K,I,JC)-CSSUH2(K,JC,JB)/DLT
                                 enddo
                    !CALL UPSTREAM_CONSTITUENT(C2(:,:,JC),CSSB(:,:,JC))
                             endif
                         endif
                     endif
                 endif
                 if(DN_HEAD(jb))then
                     do k = kt, KB(id + 1)
                         idt = id + 1
                         if(QDH1(k, jb)>=0.0)idt = id
                !CSSDH1(K,JC,JB) = C1S(K,IDT,JC)*QDH1(K,JB)
 
                     enddo
                     if(DH_INTERNAL(jb))then
                         if(DHS(jb)/=US(JBDH(jb)) .OR. UHS(JBDH(jb))/=DS(jb))  &
                          & then
                             if(JBDH(jb)>=BS(jw) .AND. JBDH(jb)<=BE(jw))then
                                 i = DHS(jb)
                                 do k = kt, KB(id + 1)
                      !CSSB(K,I,JC) = CSSB(K,I,JC)+CSSDH2(K,JC,JB)/DLT
                                 enddo
                    !CALL DOWNSTREAM_CONSTITUENT(C2(:,:,JC),CSSB(:,:,JC))
                             endif
                         endif
                     endif
                 endif
             endif
              ! END OF TP AND TN MASS BALANCES
 
 
 
         enddo
     enddo
 
!**** Kinetic fluxes
 
     do jw = 1, nwb
         kt = KTWB(jw)   ! SW 10/25/2017
         if(FLUX(jw))call KINETIC_FLUXES
     enddo
 
    !SP CEMA
     if(sediment_diagenesis)then
         if(cemarelatedcode .AND. includebedconsolidation)                     &
          & call CEMASEDIMENTMODEL
         if(cemarelatedcode .AND. includecemaseddiagenesis)                    &
          & call CEMASEDIMENTDIAGENESIS
         if(cemarelatedcode .AND. includecemaseddiagenesis .AND.               &
          & bubbles_calculation)call CEMACALCULATERISEVELOCITY
         if(cemarelatedcode .AND. includecemaseddiagenesis .AND.               &
          & applybubbturb .AND. bubbles_calculation)                           &
          & call CEMABUBBLESRELEASETURBULENCE
         if(cemarelatedcode .AND. includecemaseddiagenesis .AND.               &
          & bubbles_calculation)call CEMABUBBLESTRANSPORT
         if(cemarelatedcode .AND. includecemaseddiagenesis .AND.               &
          & bubbles_calculation)call CEMABUBBWATTRANSFER
         if(cemarelatedcode .AND. includecemaseddiagenesis .AND.               &
          & bubbles_calculation)call CEMABUBBLESRELEASE
         if(includefftlayer)call CEMAFFTLAYERCODE
     endif
    !End SP CEMA
 
!**** Constituent transport
 
!!$OMP PARALLEL DO !!PRIVATE(I,JC,KT,JB,JW,DT)    !I,JC,KT,JW,JB,CNEW,SSB,SSK,COLD,AT,VT,CT,DT)
     do jac = 1, nac
                    !CONCURRENT(JAC=1:NAC)              !JAC=1,NAC
         jc = CN(jac)
         do jw = 1, nwb
             kt = KTWB(jw)
             do jb = BS(jw), BE(jw)
                 if(BR_INACTIVE(jb))cycle
                                   ! SW 6/12/2017
                 iu = CUS(jb)
                 id = DS(jb)
      !    DO JAC=1,NAC
      !      JC   =  CN(JAC)
                 cold => C1S(:, :, jc)
                 call HORIZONTAL_MULTIPLIERS
                 call VERTICAL_MULTIPLIERS
 
     !       CNEW => C1(:,:,JC)
     !       SSB  => CSSB(:,:,JC)
     !       SSK  => CSSK(:,:,JC)
     !       CALL HORIZONTAL_TRANSPORT
                 do i = iu, id
                     do k = kt, KB(i)
                         DT(k, i) = (C1S(k, i, jc)*BH2(k, i)/dlt + (ADX(k, i)* &
                                  & BHR1(k, i) - ADX(k, i - 1)*BHR1(k, i - 1)) &
                                  & /DLX(i) + (1.0D0 - THETA(jw))              &
                                  & *(ADZ(k, i)*BB(k, i) - ADZ(k - 1, i)       &
                                  & *BB(k - 1, i)) + CSSB(k, i, jc)/DLX(i))    &
                                  & *dlt/BH1(k, i) + CSSK(k, i, jc)*dlt
                     enddo
                 enddo
                 do i = iu, id
        !      CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,CNEW(:,I))
                     BTA1(kt) = VT(kt, i)
                     GMA1(kt) = DT(kt, i)
                     do k = kt + 1, KB(i)
                         BTA1(k) = VT(k, i) - AT(k, i)/BTA1(k - 1)*CT(k - 1, i)
                         GMA1(k) = DT(k, i) - AT(k, i)/BTA1(k - 1)*GMA1(k - 1)
                     enddo
                     C1(KB(i), i, jc) = GMA1(KB(i))/BTA1(KB(i))
                     do k = KB(i) - 1, kt, -1
                         C1(k, i, jc) = (GMA1(k) - CT(k, i)*C1(k + 1, i, jc))  &
                                      & /BTA1(k)
                     enddo
                 enddo
             enddo
         enddo
     enddo
!!$OMP END PARALLEL DO
     if(derived_calc)call DERIVED_CONSTITUENTS
 
     end subroutine WQCONSTITUENTS
