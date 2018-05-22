!*==balances.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
 
     subroutine BALANCES
 
     use MAIN
     use GLOBAL
     use NAMESC
     use GEOMC
     use LOGICC
     use PREC
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
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: tnplant, tnsed, tnwb, tpplant, tpsed, tpwb, voldtjw, volevjw,     &
           & volicejw, volinjw, voloutjw, volprjw, voltrbjw, volwdjw
     external RESTART_OUTPUT
!
!*** End of declarations rewritten by SPAG
!
 
!***********************************************************************************************************************************
!*   TASK 2.6: BALANCES                                                       
!***********************************************************************************************************************************
 
    !QINT  = 0.0
    !QOUTT = 0.0
!    **
     volsr = 0.0
     voltr = 0.0
 
     do jw = 1, nwb
         volinjw = 0.0
         volprjw = 0.0
         voloutjw = 0.0
         volwdjw = 0.0
         volevjw = 0.0
         voldtjw = 0.0
         voltrbjw = 0.0
         volicejw = 0.0
 
         tpwb = 0.0
         tpsed = 0.0
         tnwb = 0.0
         tnsed = 0.0
         tnplant = 0.0
         tpplant = 0.0
 
         kt = KTWB(jw)
         if(VOLUME_BALANCE(jw))then
             do jb = BS(jw), BE(jw)
                 if(BR_INACTIVE(jb))cycle
                 VOLSBR(jb) = VOLSBR(jb) + DLVOL(jb)
                 VOLTBR(jb) = VOLEV(jb) + VOLPR(jb) + VOLTRB(jb) + VOLDT(jb)   &
                            & + VOLWD(jb) + VOLUH(jb) + VOLDH(jb) + VOLIN(jb)  &
                            & + VOLOUT(jb) + VOLICE(jb)
          !SP CEMA
                 if(sediment_diagenesis)then
                     if(cemarelatedcode .AND. includebedconsolidation)         &
                      & VOLTBR(jb) = VOLTBR(jb) + VOLCEMA(jb)
                 endif
          !End SP CEMA
                 volsr(jw) = volsr(jw) + VOLSBR(jb)
                 voltr(jw) = voltr(jw) + VOLTBR(jb)
          !QINT(JW)   = QINT(JW) +VOLIN(JB)+VOLTRB(JB)+VOLDT(JB)+VOLPR(JB)
          !QOUTT(JW)  = QOUTT(JW)-VOLEV(JB)-VOLWD(JB) -VOLOUT(JB)
                 if(ABS(VOLSBR(jb) - VOLTBR(jb))>vtol .AND. VOLTBR(jb)         &
                  & >100.0*vtol)then
                     if(volume_warning)then
                         write(wrn, '(A,F0.3,3(:/A,E15.8,A))')                 &
                              &'COMPUTATIONAL WARNING AT JULIAN DAY = ', jday, &
                              &'SPATIAL CHANGE  =', VOLSBR(jb), ' M^3',        &
                              &'TEMPORAL CHANGE =', VOLTBR(jb), ' M^3',        &
                              &'VOLUME ERROR    =', VOLSBR(jb) - VOLTBR(jb),   &
                              &' M^3'
                         write(wrn, *)'LAYER CHANGE:', LAYERCHANGE(jw)
                         write(wrn, *)'SZ', sz, 'Z', z, 'H2KT', h2(kt, 1:imx), &
                                     &'H1KT', h1(kt, 1:imx), 'WSE', elws, 'Q', &
                                    & q, 'QC', qc, 'T1', t1(kt, 1:imx), 'T2',  &
                                    & t2(kt, 1:imx), 'SUKT', su(kt, 1:imx),    &
                                     &'UKT', u(kt, 1:imx), 'QIN', qinsum,      &
                                    & 'QTR', qtr, 'QWD', qwd
                         warning_open = .TRUE.
                         volume_warning = .FALSE.
                     endif
                 endif
                 if(volsr(jw)/=0.0)DLVR(jw) = (voltr(jw) - volsr(jw))/volsr(jw)&
                  & *100.0
                 volinjw = volinjw + VOLIN(jb)
                 volprjw = volprjw + VOLPR(jb)
                 voloutjw = voloutjw + VOLOUT(jb)
                 volwdjw = volwdjw + VOLWD(jb)
                 volevjw = volevjw + VOLEV(jb)
                 voldtjw = voldtjw + VOLDT(jb)
                 voltrbjw = voltrbjw + VOLTRB(jb)
                 volicejw = volicejw + VOLICE(jb)
             enddo
 
             if(CONTOUR(jw))then
                 if(jday + (dlt/day)>=NXTMCP(jw) .OR. jday + (dlt/day)         &
                  & >=CPLD(CPLDP(jw) + 1, jw))then
                     if(VOLUME_BALANCE(jw))then
                         write(flowbfn,                                        &
                              &'(F10.3,",",1X,I3,",",11(E16.8,",",1X))')jday,  &
                             & jw, volinjw, volprjw, voloutjw, volwdjw,        &
                             & volevjw, voldtjw, voltrbjw, volicejw, DLVR(jw)
                     else
                         write(flowbfn,                                        &
                              &'(F10.3,",",1X,I3,",",10(E16.8,",",1X))')jday,  &
                             & jw, volinjw, volprjw, voloutjw, volwdjw,        &
                             & volevjw, voldtjw, voltrbjw, volicejw
                     endif
                 endif
 
             endif
                ! CONTOUR INTERVAL FOR WRITING OUT FLOW BALANCE
         endif  ! VOLUME BALANCE
         if(ENERGY_BALANCE(jw))then
             ESR(jw) = 0.0
             ETR(jw) = 0.0
             do jb = BS(jw), BE(jw)
                 ETBR(jb) = EBRI(jb) + TSSEV(jb) + TSSPR(jb) + TSSTR(jb)       &
                          & + TSSDT(jb) + TSSWD(jb) + TSSUH(jb) + TSSDH(jb)    &
                          & + TSSIN(jb) + TSSOUT(jb) + TSSS(jb) + TSSB(jb)     &
                          & + TSSICE(jb)
                 ESBR(jb) = 0.0
                 do i = CUS(jb), DS(jb)
                     do k = kt, KB(i)
                         ESBR(jb) = ESBR(jb) + t1(k, i)*DLX(i)*BH1(k, i)
                     enddo
                 enddo
                 ETR(jw) = ETR(jw) + ETBR(jb)
                 ESR(jw) = ESR(jw) + ESBR(jb)
             enddo
         endif
         if(MASS_BALANCE(jw))then
             do jb = BS(jw), BE(jw)
                 do jc = 1, nac
                     CMBRS(CN(jc), jb) = 0.0
                     do i = CUS(jb), DS(jb)
                         do k = kt, KB(i)
                             CMBRS(CN(jc), jb) = CMBRS(CN(jc), jb)             &
                               & + C1(k, i, CN(jc))*DLX(i)*BH1(k, i)
                             CMBRT(CN(jc), jb) = CMBRT(CN(jc), jb)             &
                               & + (CSSB(k, i, CN(jc)) + CSSK(k, i, CN(jc))    &
                               & *BH1(k, i)*DLX(i))*dlt
                         enddo
                     enddo
                 enddo
                 if(derived_calc)then
                     do i = CUS(jb), DS(jb)
                         do k = kt, KB(i)
                             tpwb = tpwb + TP(k, i)*VOL(k, i)/1000.
                                                              ! kg
                             tpsed = tpsed + SEDP(k, i)*VOL(k, i)/1000.
                                                              ! kg
                             tnwb = tnwb + TN(k, i)*VOL(k, i)/1000.
                                                              ! kg
                             tnsed = tnsed + SEDN(k, i)*VOL(k, i)/1000.
                                                              ! kg
                             PFLUXIN(jw) = PFLUXIN(jw) + SEDPINFLUX(k, i)      &
                               & *VOL(k, i)/1000.                            ! kg
                             NFLUXIN(jw) = NFLUXIN(jw) + SEDNINFLUX(k, i)      &
                               & *VOL(k, i)/1000.                            ! kg
                             do m = 1, nmc
                                 tnplant = tnplant + MAC(k, i, m)*VOL(k, i)    &
                                   & *MN(m)/1000.
                                 tpplant = tpplant + MAC(k, i, m)*VOL(k, i)    &
                                   & *MP(m)/1000.
                             enddo
                             do m = 1, nep
                                 tnplant = tnplant + EPM(k, i, m)*EN(m)/1000.
                                 tpplant = tpplant + EPM(k, i, m)*EP(m)/1000.
                             enddo
                         enddo
                     enddo
                 endif
 
!                MACROPHYTES
                 do m = 1, nmc
                     if(MACROPHYTE_CALC(jw, m))then
                         MACMBRS(jb, m) = 0.0
                         do i = CUS(jb), DS(jb)
                             if(KTICOL(i))then
                                 jt = KTI(i)
                             else
                                 jt = KTI(i) + 1
                             endif
                             je = KB(i)
                             do j = jt, je
                                 if(j<kt)then
                                     colb = EL(j + 1, i)
                                 else
                                     colb = EL(kt + 1, i)
                                 endif
                  !COLDEP=ELWS(I)-COLB
                                 coldep = EL(kt, i) - z(i)*COSA(jb) - colb
                                                      ! cb 3/7/16
                                 MACMBRS(jb, m) = MACMBRS(jb, m)               &
                                   & + MACRM(j, kt, i, m)
                                 MACMBRT(jb, m) = MACMBRT(jb, m)               &
                                   & + (MACSS(j, kt, i, m)*coldep*CW(j, i)     &
                                   & *DLX(i))*dlt
                             enddo
                             do k = kt + 1, KB(i)
                                 jt = k
                                 je = KB(i)
                                 do j = jt, je
                                     MACMBRS(jb, m) = MACMBRS(jb, m)           &
                                       & + MACRM(j, k, i, m)
!                                    MACMBRT(JB,M) =
! MACMBRT(JB,M)+(MACSS(J,K,I,M)*H2(K,I)*CW(J,I)*DLX(I))*DLT
                                     MACMBRT(jb, m) = MACMBRT(jb, m)           &
                                       & + (MACSS(j, k, i, m)                  &
                                       & *(CW(j, i)/B(k, i))*BH1(k, i)*DLX(i)) &
                                       & *dlt
                                 enddo
                             enddo
                         enddo
                     endif
                 enddo
!                END MACROPHYTES
             enddo
 
             if(CONTOUR(jw) .AND. derived_calc)then
                 if(jday + (dlt/day)>=NXTMCP(jw) .OR. jday + (dlt/day)         &
                  & >=CPLD(CPLDP(jw) + 1, jw))then
                     if(sediment_diagenesis)then
                         write(massbfn,                                        &
                              &'(F10.3,",",1X,I3,",",40(E16.8,",",1X))')jday,  &
                             & jw, tpwb, tpsed, tpplant, TPOUT(jw), TPTRIB(jw),&
                             & TPDTRIB(jw), TPWD(jw), TPPR(jw), TPIN(jw),      &
                             & TP_SEDSOD_PO4(jw), PFLUXIN(jw), SDPFLUX(jw),    &
                             & tnwb, tnsed, tnplant, TNOUT(jw), TNTRIB(jw),    &
                             & TNDTRIB(jw), TNWD(jw), TNPR(jw), TNIN(jw),      &
                             & TN_SEDSOD_NH4(jw), NFLUXIN(jw), SDNH4FLUX(jw),  &
                             & SDNO3FLUX(jw)
 
                     else
                         write(massbfn,                                        &
                              &'(F10.3,",",1X,I3,",",30(E16.8,",",1X))')jday,  &
                             & jw, tpwb, tpsed, tpplant, TPOUT(jw), TPTRIB(jw),&
                             & TPDTRIB(jw), TPWD(jw), TPPR(jw), TPIN(jw),      &
                             & TP_SEDSOD_PO4(jw), PFLUXIN(jw), tnwb, tnsed,    &
                             & tnplant, TNOUT(jw), TNTRIB(jw), TNDTRIB(jw),    &
                             & TNWD(jw), TNPR(jw), TNIN(jw), TN_SEDSOD_NH4(jw),&
                             & NFLUXIN(jw)
 
                     endif
                 endif
 
             endif
                ! CONTOUR INTERVAL FOR WRITING OUT FLOW BALANCE
 
 
         endif
              ! MASS BALANCE
     enddo
 
     end subroutine BALANCES
