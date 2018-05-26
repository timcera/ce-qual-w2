!*==hydroinout.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine HYDROINOUT
 
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
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND) :: akbr, cgas, dtvl, elw, fw, qsumm, qtrfr, rhotr, tm, tsum,     &
               & vptg, vqtr, vqtri
     integer :: jbd, jbu, jlat, jwu
     external RESTART_OUTPUT
!
!*** End of declarations rewritten by SPAG
!
!***********************************************************************************************************************************
!**  Task 2.1: Hydrodynamic sources/sinks                                     
!***********************************************************************************************************************************
 
!    **
     qinsum = 0.0
!***********************************************************************************************************************************
!**  Task 2.1: Hydrodynamic sources/sinks                                     
!***********************************************************************************************************************************
 
!    **
     tinsum = 0.0
!***********************************************************************************************************************************
!**  Task 2.1: Hydrodynamic sources/sinks                                     
!***********************************************************************************************************************************
 
!    **
     cinsum = 0.0
!***********************************************************************************************************************************
!**  Task 2.1: Hydrodynamic sources/sinks                                     
!***********************************************************************************************************************************
 
!    **
     uxbr = 0.0
!***********************************************************************************************************************************
!**  Task 2.1: Hydrodynamic sources/sinks                                     
!***********************************************************************************************************************************
 
!    **
     uybr = 0.0
!***********************************************************************************************************************************
!**  Task 2.1: Hydrodynamic sources/sinks                                     
!***********************************************************************************************************************************
 
!    **
     tdgon = .FALSE.                                                                        ! cb 1/16/13
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             if(BR_INACTIVE(jb))cycle
             iu = CUS(jb)
             id = DS(jb)
             tsum = 0.0
             csum = 0.0
             QSUM(jb) = 0.0
             qout(:, jb) = 0.0
             TOUT(jb) = 0.0
             cout(:, jb) = 0.0
 
!******      Densities
 
             do i = iu - 1, id + 1
                 do k = kt, KB(i)
                     TISS(k, i) = 0.0
                     do js = 1, nss
                         TISS(k, i) = TISS(k, i) + SS(k, i, js)
                     enddo
                     RHO(k, i) = DENSITY(T2(k, i), DMAX1(TDS(k, i), 0.0D0),    &
                               & DMAX1(TISS(k, i), 0.0D0))
                 enddo
             enddo
!            v3.5 deleted pumpback code from v3.2
             do js = 1, NSTR(jb)
                 if(QSTR(js, jb)/=0.0)call DOWNSTREAM_WITHDRAWAL(js)
             enddo
             do k = kt, KB(id)
                 QSUM(jb) = QSUM(jb) + qout(k, jb)
                 tsum = tsum + qout(k, jb)*T2(k, id)
                 csum(cn(1:nac)) = csum(cn(1:nac)) + qout(k, jb)               &
                                 & *C2(k, id, cn(1:nac))
             enddo
             if(QSUM(jb)/=0.0)then
                 TOUT(jb) = tsum/QSUM(jb)
                 cout(cn(1:nac), jb) = csum(cn(1:nac))/QSUM(jb)
             endif
             if(QSUM(jb)/=0.0 .AND. DAM_OUTFLOW(jb))then                                                               !TC 08/03/04
                 tinsum(JBDAM(jb)) = (tsum + qinsum(JBDAM(jb))*tinsum(JBDAM(jb)&
                                   & ))/(QSUM(jb) + qinsum(JBDAM(jb)))
                 cinsum(cn(1:nac), JBDAM(jb))                                  &
                   & = (csum(cn(1:nac)) + qinsum(JBDAM(jb))                    &
                   & *cinsum(cn(1:nac), JBDAM(jb)))                            &
                   & /(QSUM(jb) + qinsum(JBDAM(jb)))
                 qinsum(JBDAM(jb)) = qinsum(JBDAM(jb)) + QSUM(jb)
             endif
         enddo
     enddo
     ilat = 0
     jww = nwd
     withdrawals = jww>0
     if(nwdt>nwd)qwd(nwd + 1:nwdt) = 0.0
                                       ! SW 10/30/2017
     jtt = ntr
     tributaries = jtt>0
     if(ntrt>ntr)qtr(ntr + 1:ntrt) = 0.0
                                       ! SW 10/30/2017
     jss = NSTR
     if(spillway)then
         call SPILLWAY_FLOW
         do js = 1, nsp
 
!******      Positive flows
 
             jlat = 0
             jbu = JBUSP(js)
             jbd = JBDSP(js)
             tdgon = .FALSE. ! cb 1/16/13
             jsg = js
             nnsg = 0
             if(CAC(ndo)=='      ON' .AND. GASSPC(js)=='      ON')             &
              & tdgon = .TRUE.
             if(QSP(js)>=0.0)then
                 if(LATERAL_SPILLWAY(js))then
                     jww = jww + 1
                     IWD(jww) = IUSP(js)
                     qwd(jww) = QSP(js)
                     KTWD(jww) = KTUSP(js)
                     KBWD(jww) = KBUSP(js)
                     EWD(jww) = ESP(js)
                     JBWD(jww) = jbu
                     i = MAX(CUS(JBWD(jww)), IWD(jww))
                     jb = JBWD(jww)
                     jw = JWUSP(js)
                     kt = KTWB(jw)
                     jwd = jww
                     call LATERAL_WITHDRAWAL
                                    !(JWW)
                     do k = KTW(jww), KBW(jww)
                         QSS(k, i) = QSS(k, i) - QSW(k, jww)
                     enddo
                     if(IDSP(js)/=0)then ! cb 9/11/13
                         jtt = jtt + 1
                         qtr(jtt) = QSP(js)
                         ITR(jtt) = IDSP(js)
                         PLACE_QTR(jtt) = PDSPC(js)==' DENSITY'
                         SPECIFY_QTR(jtt) = PDSPC(js)==' SPECIFY'
                         if(SPECIFY_QTR(jtt))then
                             ELTRT(jtt) = ETDSP(js)
                             ELTRB(jtt) = EBDSP(js)
                         endif
                         JBTR(jtt) = jbd
                     endif
                        ! cb 9/11/13
                     if(IDSP(js)/=0 .AND. QSP(js)>0.0)then
                         tsum = 0.0
                         qsumm = 0.0
                         csum = 0.0
                         do k = KTW(jww), KBW(jww)
                             qsumm = qsumm + QSW(k, jww)
                             tsum = tsum + QSW(k, jww)*T2(k, IWD(jww))
                             csum(cn(1:nac)) = csum(cn(1:nac)) + QSW(k, jww)   &
                               & *C2(k, IWD(jww), cn(1:nac))
                         enddo
                         TTR(jtt) = tsum/qsumm
                         do jc = 1, nac
                             CTR(cn(jc), jtt) = csum(cn(jc))/qsumm
                             if(cn(jc)==ndo .AND. GASSPC(js)=='      ON' .AND. &
                              & QSP(js)>0.0)then
                                 TDG_SPILLWAY(jww, js) = .TRUE.
                                 call TOTAL_DISSOLVED_GAS(0, PALT(i), 0, js,   &
                                   & TTR(jtt), CTR(cn(jc), jtt))                     ! DO
                             endif
                             if(cn(jc)==ngn2 .AND. GASSPC(js)=='      ON' .AND.&
                              & QSP(js)>0.0)then
                                 TDG_SPILLWAY(jww, js) = .TRUE.
                                 call TOTAL_DISSOLVED_GAS(1, PALT(i), 0, js,   &
                                   & TTR(jtt), CTR(cn(jc), jtt))                     ! N2
                             endif
                         enddo
                     elseif(CAC(ndo)=='      ON' .AND. GASSPC(js)              &
                          & =='      ON' .AND. QSP(js)>0.0)then
                         TDG_SPILLWAY(jww, js) = .TRUE.
                     endif
 
                 else
                     jss(jbu) = jss(jbu) + 1
                     KTSW(jss(jbu), jbu) = KTUSP(js)
                     KBSW(jss(jbu), jbu) = KBUSP(js)
                     jb = jbu
                     POINT_SINK(jss(jbu), jbu) = .TRUE.
                     id = IUSP(js)
                     QSTR(jss(jbu), jbu) = QSP(js)
                     ESTR(jss(jbu), jbu) = ESP(js)
                     kt = KTWB(JWUSP(js))
                     jw = JWUSP(js)
                     call DOWNSTREAM_WITHDRAWAL(jss(jbu))
                     QSUM(jb) = 0.0
                     tsum = 0.0
                     csum = 0.0
                     do k = kt, KB(id)
                         QSUM(jb) = QSUM(jb) + qout(k, jb)
                         tsum = tsum + qout(k, jb)*T2(k, id)
                         do jc = 1, nac
                             if(cn(jc)==ndo .AND. CAC(ndo)=='      ON' .AND.   &
                              & GASSPC(js)=='      ON' .AND. QSP(js)>0.0)then                                             ! MM 5/21/2009
                                 t2r4 = T2(k, id)
                                 cgas = C2(k, id, cn(jc))                                                                   ! MM 5/21/2009
                                 call TOTAL_DISSOLVED_GAS(0, PALT(id), 0, js,  &
                                   & t2r4, cgas)                          ! O2
                                 csum(cn(jc)) = csum(cn(jc)) + qout(k, jb)*cgas
                !ELSEIF (CN(JC)==NGN2 .AND. CAC(NGN2) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN     ! SW 10/27/15
                             elseif(cn(jc)==ngn2 .AND. GASSPC(js)              &
                                  & =='      ON' .AND. QSP(js)>0.0)then                           ! SW 10/27/15
                                 if(CAC(ngn2)=='      ON')then                                                                   ! cb 1/13/16
                                     t2r4 = T2(k, id)
                                     cgas = C2(k, id, cn(jc))                                                                 !
                                     call TOTAL_DISSOLVED_GAS(1, PALT(id), 0,  &
                                       & js, t2r4, cgas)                   ! N2
                                     csum(cn(jc)) = csum(cn(jc)) + qout(k, jb) &
                                       & *cgas
                                 endif
                             else
                                 csum(cn(jc)) = csum(cn(jc)) + qout(k, jb)     &
                                   & *C2(k, id, cn(jc))
                             endif
                         enddo
                     enddo
                     if(QSUM(jb)/=0.0)then
                         TOUT(jb) = tsum/QSUM(jb)
                         cout(cn(1:nac), jb) = csum(cn(1:nac))/QSUM(jb)
                     endif
                     if(IDSP(js)/=0 .AND. US(jbd)==IDSP(js))then
                         qsumm = 0.0
                         tsum = 0.0
                         csum = 0.0
                         do k = kt, KB(id)
                             qsumm = qsumm + QNEW(k)
                             tsum = tsum + QNEW(k)*T2(k, id)
                             do jc = 1, nac
                                 if(cn(jc)==ndo .AND. CAC(ndo)                 &
                                  & =='      ON' .AND. GASSPC(js)              &
                                  & =='      ON' .AND. QSP(js)>0.0)then                                                             ! MM 5/21/2009
                                     t2r4 = T2(k, id)
                                     cgas = C2(k, id, cn(jc))                                                                       ! MM 5/21/2009
                                     call TOTAL_DISSOLVED_GAS(0, PALT(id), 0,  &
                                       & js, t2r4, cgas)
                                     csum(cn(jc)) = csum(cn(jc)) + QNEW(k)*cgas
                  !ELSEIF (CN(JC)==NGN2 .AND. CAC(NGN2) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN             ! SW 10/27/15
                                 elseif(cn(jc)==ngn2 .AND. GASSPC(js)          &
                                       &=='      ON' .AND. QSP(js)>0.0)then                                 ! SW 10/27/15
                                     if(CAC(ngn2)=='      ON')then                                                                 ! cb 1/13/16
                                         t2r4 = T2(k, id)
                                         cgas = C2(k, id, cn(jc))                                                                     !
                                         call TOTAL_DISSOLVED_GAS(1, PALT(id), &
                                           & 0, js, t2r4, cgas)
                                         csum(cn(jc)) = csum(cn(jc)) + QNEW(k) &
                                           & *cgas
                                     endif
                                 else
                                     csum(cn(jc)) = csum(cn(jc)) + QNEW(k)     &
                                       & *C2(k, id, cn(jc))
                                 endif
                             enddo
                         enddo
                         if(qsumm/=0.0)then
                             tinsum(jbd) = (tsum + qinsum(jbd)*tinsum(jbd))    &
                               & /(qsumm + qinsum(jbd))
                             cinsum(cn(1:nac), jbd)                            &
                               & = (csum(cn(1:nac)) + qinsum(jbd)              &
                               & *cinsum(cn(1:nac), jbd))/(qsumm + qinsum(jbd))
                             qinsum(jbd) = qinsum(jbd) + qsumm
                         endif
                     elseif(IDSP(js)/=0)then
                         jtt = jtt + 1
                         qtr(jtt) = QSP(js)
                         ITR(jtt) = IDSP(js)
                         PLACE_QTR(jtt) = PDSPC(js)==' DENSITY'
                         SPECIFY_QTR(jtt) = PDSPC(js)==' SPECIFY'
                         if(SPECIFY_QTR(jtt))then
                             ELTRT(jtt) = ETDSP(js)
                             ELTRB(jtt) = EBDSP(js)
                         endif
                         JBTR(jtt) = jbd
                         TTR(jtt) = TOUT(jb)
                         do jc = 1, nac
                             CTR(cn(jc), jtt) = cout(cn(jc), jb)
                             if(cn(jc)==ndo .AND. GASSPC(js)=='      ON' .AND. &
                              & QSP(js)>0.0)                                   &
                              & call TOTAL_DISSOLVED_GAS(0, PALT(ITR(jtt)), 0, &
                              & js, TTR(jtt), CTR(cn(jc), jtt))
                             if(cn(jc)==ngn2 .AND. GASSPC(js)=='      ON' .AND.&
                              & QSP(js)>0.0)                                   &
                              & call TOTAL_DISSOLVED_GAS(1, PALT(ITR(jtt)), 0, &
                              & js, TTR(jtt), CTR(cn(jc), jtt))                                   ! SW 10/27/15
                         enddo
                     endif
                 endif
             elseif(QSP(js)<0.0)then
                 jtt = jtt + 1
                 jww = jww + 1
                 IWD(jww) = IDSP(js)
                 ITR(jtt) = IUSP(js)
                 qtr(jtt) = -QSP(js)
                 qwd(jww) = -QSP(js)
                 KTWD(jww) = KTDSP(js)
                 KBWD(jww) = KBDSP(js)
                 EWD(jww) = ESP(js)
                 PLACE_QTR(jtt) = PUSPC(js)==' DENSITY'
                 SPECIFY_QTR(jtt) = PUSPC(js)==' SPECIFY'
                 if(SPECIFY_QTR(jtt))then
                     ELTRT(jtt) = ETUSP(js)
                     ELTRB(jtt) = EBUSP(js)
                 endif
                 JBTR(jtt) = jbu
                 JBWD(jww) = jbd
                 i = MAX(CUS(JBWD(jww)), IWD(jww))
                 jb = JBWD(jww)
                 jw = JWDSP(js)
                 kt = KTWB(jw)
                 jwd = jww
                 call LATERAL_WITHDRAWAL
                                  !(JWW)
                 do k = KTW(jww), KBW(jww)
                     QSS(k, i) = QSS(k, i) - QSW(k, jww)
                 enddo
                 if(IDSP(js)/=0)then
                     tsum = 0.0
                     qsumm = 0.0
                     csum = 0.0
                     do k = KTW(jww), KBW(jww)
                         qsumm = qsumm + QSW(k, jww)
                         tsum = tsum + QSW(k, jww)*T2(k, IWD(jww))
                         csum(cn(1:nac)) = csum(cn(1:nac)) + QSW(k, jww)       &
                           & *C2(k, IWD(jww), cn(1:nac))
                     enddo
                     TTR(jtt) = tsum/qsumm
                     do jc = 1, nac
                         CTR(cn(jc), jtt) = csum(cn(jc))/qsumm
                         if(cn(jc)==ndo .AND. GASSPC(js)=='      ON' .AND.     &
                          & QSP(js)>0.0)then
                             TDG_SPILLWAY(jww, js) = .TRUE.
                             call TOTAL_DISSOLVED_GAS(0, PALT(i), 0, js,       &
                               & TTR(jtt), CTR(cn(jc), jtt))                          ! O2
                         endif
                         if(cn(jc)==ngn2 .AND. GASSPC(js)=='      ON' .AND.    &
                          & QSP(js)>0.0)then
                             TDG_SPILLWAY(jww, js) = .TRUE.
                             call TOTAL_DISSOLVED_GAS(1, PALT(i), 0, js,       &
                               & TTR(jtt), CTR(cn(jc), jtt))                          ! N2
                         endif
                     enddo
                 elseif(CAC(ndo)=='      ON' .AND. GASSPC(js)=='      ON' .AND.&
                      & QSP(js)>0.0)then
                     TDG_SPILLWAY(jww, js) = .TRUE.
                 endif
             endif
         enddo
     endif
     if(pumps)then
         do jp = 1, npu
             jlat = 0
             jwu = JWUPU(jp)
             jbu = JBUPU(jp)
             jbd = JBDPU(jp)
             tdgon = .FALSE. ! cb 1/16/13
             if(LATERAL_PUMP(jp))then
                 elw = EL(KTWB(jwu), IUPU(jp)) - Z(IUPU(jp))*COSA(jbu)
        !    JWW       = JWW+1      ! SW 9/25/13
        !    JBWD(JWW) = JBU
        !    IWD(JWW)  = IUPU(JP)
                 jww = jww + 1     ! SW 10/30/2017
                 JBWD(jww) = jbu
                 IWD(jww) = IUPU(jp)
             else
                 elw = EL(KTWB(jwu), IUPU(jp)) - Z(IUPU(jp))*COSA(jbu)         &
                     & - SINA(jbu)*DLX(IUPU(jp))*0.5
        !    JSS(JBU)                 =  JSS(JBU)+1     ! SW 9/25/13
                 jss(jbu) = jss(jbu) + 1               ! SW 10/30/2017
             endif
             if(jday>=ENDPU(jp))PUMPON(jp) = .FALSE.                                                       !  CB 1/13/06
             if(jday>=STRTPU(jp) .AND. jday<ENDPU(jp))then
 
                 if(elw<=EOFFPU(jp))PUMPON(jp) = .FALSE.                                                    ! CB 1/13/06
                 if(elw>EOFFPU(jp) .AND. QPU(jp)>0.0)then
                     if(elw>=EONPU(jp))PUMPON(jp) = .TRUE.
                     if(PUMPON(jp))then
                         if(LATERAL_PUMP(jp))then
                             jlat = 1
                !JWW       = JWW+1               ! SW 9/25/13
                !JBWD(JWW) = JBU
                !IWD(JWW)  = IUPU(JP)
                             qwd(jww) = QPU(jp)
                             KTWD(jww) = KTPU(jp)
                             KBWD(jww) = KBPU(jp)
                             EWD(jww) = EPU(jp)
                             i = MAX(CUS(JBWD(jww)), IWD(jww))
                             jb = JBWD(jww)
                             jw = jwu
                             kt = KTWB(jw)
                             jwd = jww
                             call LATERAL_WITHDRAWAL
                                                ! (JWW)
                             do k = KTW(jww), KBW(jww)
                                 QSS(k, i) = QSS(k, i) - QSW(k, jww)
                             enddo
                             if(IDPU(jp)/=0)then  ! MOVED CODE SW 9/25/13
                                 jtt = jtt + 1
                                 qtr(jtt) = QPU(jp)
                                 ITR(jtt) = IDPU(jp)
                                 PLACE_QTR(jtt) = PPUC(jp)==' DENSITY'
                                 SPECIFY_QTR(jtt) = PPUC(jp)==' SPECIFY'
                                 if(SPECIFY_QTR(jtt))then
                                     ELTRT(jtt) = ETPU(jp)
                                     ELTRB(jtt) = EBPU(jp)
                                 endif
                                 JBTR(jtt) = jbd
                                 tsum = 0.0
                                 qsumm = 0.0
                                 csum(cn(1:nac)) = 0.0
                                 do k = KTW(jww), KBW(jww)
                                     qsumm = qsumm + QSW(k, jww)
                                     tsum = tsum + QSW(k, jww)*T2(k, IWD(jww))
                                     csum(cn(1:nac)) = csum(cn(1:nac))         &
                                       & + QSW(k, jww)                         &
                                       & *C2(k, IWD(jww), cn(1:nac))
                                 enddo
                                 if(qsumm>0.0)then
                                     TTR(jtt) = tsum/qsumm
                                     CTR(cn(1:nac), jtt) = csum(cn(1:nac))     &
                                       & /qsumm
                                 endif
                             endif
                              ! SW 9/25/13 END MOVED CODE
                         else
                !JSS(JBU)                 =  JSS(JBU)+1     ! SW 9/25/13
                             KTSW(jss(jbu), jbu) = KTPU(jp)
                             KBSW(jss(jbu), jbu) = KBPU(jp)
                             jb = jbu
                             POINT_SINK(jss(jbu), jbu) = .TRUE.
                             id = IUPU(jp)
                             QSTR(jss(jbu), jbu) = QPU(jp)
                             ESTR(jss(jbu), jbu) = EPU(jp)
                             kt = KTWB(jwu)
                             jw = jwu
                             call DOWNSTREAM_WITHDRAWAL(jss(jbu))
                             if(IDPU(jp)/=0 .AND. US(jbd)==IDPU(jp))then
                                 qsumm = 0.0
                                 tsum = 0.0
                                 csum = 0.0
                                 do k = kt, KB(id)
                                     qsumm = qsumm + QNEW(k)
                                     tsum = tsum + QNEW(k)*T2(k, id)
                                     csum(cn(1:nac)) = csum(cn(1:nac))         &
                                       & + QNEW(k)*C2(k, id, cn(1:nac))
                                 enddo
                                 if(qsumm/=0.0)then
                                     tinsum(jbd)                               &
                                       & = (tsum + tinsum(jbd)*qinsum(jbd))    &
                                       & /(qsumm + qinsum(jbd))
                                     cinsum(cn(1:nac), jbd)                    &
                                       & = (csum(cn(1:nac)) +                  &
                                       & cinsum(cn(1:nac), jbd)*qinsum(jbd))   &
                                       & /(qsumm + qinsum(jbd))
                                     qinsum(jbd) = qinsum(jbd) + qsumm
                                 endif
                             endif
                             QSUM(jb) = 0.0
                             tsum = 0.0
                             csum = 0.0
                             do k = kt, KB(id)
                                 QSUM(jb) = QSUM(jb) + qout(k, jb)
                                 tsum = tsum + qout(k, jb)*T2(k, id)
                                 csum(cn(1:nac)) = csum(cn(1:nac))             &
                                   & + qout(k, jb)*C2(k, id, cn(1:nac))
                             enddo
                             if(QSUM(jb)/=0.0)then
                                 TOUT(jb) = tsum/QSUM(jb)
                                 cout(cn(1:nac), jb) = csum(cn(1:nac))/QSUM(jb)
                             endif
                             if(IDPU(jp)/=0)then
                                             ! SW 9/25/13 Moved code start
                                 if(US(jbd)/=IDPU(jp) .OR. HEAD_FLOW(jbd) .OR. &
                                  & UP_HEAD(jbd))then
                                     jtt = jtt + 1
                                     qtr(jtt) = QPU(jp)
                                     ITR(jtt) = IDPU(jp)
                                     PLACE_QTR(jtt) = PPUC(jp)==' DENSITY'
                                     SPECIFY_QTR(jtt) = PPUC(jp)==' SPECIFY'
                                     if(SPECIFY_QTR(jtt))then
                                         ELTRT(jtt) = ETPU(jp)
                                         ELTRB(jtt) = EBPU(jp)
                                     endif
                                     JBTR(jtt) = jbd
                                     TTR(jtt) = TOUT(jb)
                                     CTR(cn(1:nac), jtt) = cout(cn(1:nac), jb)
                                 endif
                             endif        ! Moved code end SW 9/25/13
                         endif
            !  IF (IDPU(JP) /= 0) THEN
                !IF (US(JBD) /= IDPU(JP) .OR. HEAD_FLOW(JBD) .OR. UP_HEAD(JBD)) THEN      ! SW 9/25/13 MOVED CODE ABOVE
                !  JTT              = JTT+1
                !  QTR(JTT)         = QPU(JP)
                !  ITR(JTT)         = IDPU(JP)
                !  PLACE_QTR(JTT)   = PPUC(JP) == ' DENSITY'
                !  SPECIFY_QTR(JTT) = PPUC(JP) == ' SPECIFY'
                !  IF (SPECIFY_QTR(JTT)) THEN
                !    ELTRT(JTT) = ETPU(JP)
                !    ELTRB(JTT) = EBPU(JP)
                !  END IF
                !  JBTR(JTT) = JBD
                !  IF (JLAT == 1) THEN
                !    TSUM = 0.0; QSUMM = 0.0; CSUM(CN(1:NAC)) = 0.0
                !    DO K=KTW(JWW),KBW(JWW)
                !      QSUMM           = QSUMM          +QSW(K,JWW)
                !      TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                !      CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
                !    END DO
                !    TTR(JTT)           = TSUM           /QSUMM
                !    CTR(CN(1:NAC),JTT) = CSUM(CN(1:NAC))/QSUMM
                !  ELSE
                !    TTR(JTT)          = TOUT(JB)
                !    CTR(CN(1:NAC),JTT)= COUT(CN(1:NAC),JB)
                !  END IF
                !ELSE IF (LATERAL_PUMP(JP)) THEN
                !  TSUM      = 0.0; QSUMM = 0.0; CSUM = 0.0
                !  ILAT(JWW) = 1
                !  DO K=KTW(JWW),KBW(JWW)
                !    QSUMM           = QSUMM          +QSW(K,JWW)
                !    TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                !    CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
                !  END DO
                !  TINSUM(JBD)           = (TSUM           +TINSUM(JBD)          *QINSUM(JBD))/(QSUMM+QINSUM(JBD))
                !  CINSUM(CN(1:NAC),JBD) = (CSUM(CN(1:NAC))+CINSUM(CN(1:NAC),JBD)*QINSUM(JBD))/(QSUMM+QINSUM(JBD))
                !  QINSUM(JBD)           =  QINSUM(JBD)+QSUMM
                !END IF
            !  END IF
                     endif
                 endif
             endif
         enddo
     endif
     if(pipes)then
         yss = ys
         vss = vs
         vsts = vst
         ysts = yst
         dtps = dtp
         qolds = qold
         call PIPE_FLOW     ! (NIT)
         do jp = 1, npi
 
             if(DYNPIPE(jp)=='      ON')QPI(jp) = QPI(jp)*BP(jp)
                                                              ! SW 5/10/10
 
 
!******      Positive flows
 
             jlat = 0
             jbu = JBUPI(jp)
             jbd = JBDPI(jp)
             tdgon = .FALSE. ! cb 1/16/13
             if(QPI(jp)>=0.0)then
                 if(LATERAL_PIPE(jp))then
                     jlat = 1
                     jww = jww + 1
                     IWD(jww) = IUPI(jp)
                     qwd(jww) = QPI(jp)
                     KTWD(jww) = KTUPI(jp)
                     KBWD(jww) = KBUPI(jp)
                     EWD(jww) = EUPI(jp)
                     JBWD(jww) = jbu
                     i = MAX(CUS(JBWD(jww)), IWD(jww))
                     jb = JBWD(jww)
                     jw = JWUPI(jp)
                     kt = KTWB(jw)
                     jwd = jww
                     call LATERAL_WITHDRAWAL
                                       !(JWW)
                     do k = KTW(jww), KBW(jww)
                         QSS(k, i) = QSS(k, i) - QSW(k, jww)
                     enddo
                 else
                     jss(jbu) = jss(jbu) + 1
                     KTSW(jss(jbu), jbu) = KTDPI(jp)
                     KBSW(jss(jbu), jbu) = KBDPI(jp)
                     jb = jbu
                     POINT_SINK(jss(jbu), jbu) = .TRUE.
                     id = IUPI(jp)
                     QSTR(jss(jbu), jbu) = QPI(jp)
                     ESTR(jss(jbu), jbu) = EUPI(jp)
                     kt = KTWB(JWUPI(jp))
                     jw = JWUPI(jp)
                     call DOWNSTREAM_WITHDRAWAL(jss(jbu))
                     if(IDPI(jp)/=0 .AND. US(jbd)==IDPI(jp))then
                         qsumm = 0.0
                         tsum = 0.0
                         csum = 0.0
                         do k = kt, KB(id)
                             qsumm = qsumm + QNEW(k)
                             tsum = tsum + QNEW(k)*T2(k, id)
                             csum(cn(1:nac)) = csum(cn(1:nac)) + QNEW(k)       &
                               & *C2(k, id, cn(1:nac))
                         enddo
                         if(qsumm/=0.0)then
                             tinsum(jbd) = (tsum + qinsum(jbd)*tinsum(jbd))    &
                               & /(qsumm + qinsum(jbd))
                             cinsum(cn(1:nac), jbd)                            &
                               & = (csum(cn(1:nac)) + qinsum(jbd)              &
                               & *cinsum(cn(1:nac), jbd))/(qsumm + qinsum(jbd))
                             qinsum(jbd) = qinsum(jbd) + qsumm
                         endif
                     endif
                     QSUM(jb) = 0.0
                     tsum = 0.0
                     csum = 0.0
                     do k = kt, KB(id)
                         QSUM(jb) = QSUM(jb) + qout(k, jb)
                         tsum = tsum + qout(k, jb)*T2(k, id)
                         csum(cn(1:nac)) = csum(cn(1:nac)) + qout(k, jb)       &
                           & *C2(k, id, cn(1:nac))
                     enddo
                     if(QSUM(jb)/=0.0)then
                         TOUT(jb) = tsum/QSUM(jb)
                         cout(cn(1:nac), jb) = csum(cn(1:nac))/QSUM(jb)
                     endif
                 endif
                 if(IDPI(jp)/=0)then
                     if(US(jbd)/=IDPI(jp) .OR. HEAD_FLOW(jbd) .OR. UP_HEAD(jbd)&
                      & )then
                         jtt = jtt + 1
                         qtr(jtt) = QPI(jp)
                         ITR(jtt) = IDPI(jp)
                         PLACE_QTR(jtt) = PDPIC(jp)==' DENSITY'
                         SPECIFY_QTR(jtt) = PDPIC(jp)==' SPECIFY'
                         if(SPECIFY_QTR(jtt))then
                             ELTRT(jtt) = ETDPI(jp)
                             ELTRB(jtt) = EBDPI(jp)
                         endif
                         JBTR(jtt) = jbd
                         if(jlat==1)then
                             tsum = 0.0
                             qsumm = 0.0
                             csum(cn(1:nac)) = 0.0
                             do k = KTW(jww), KBW(jww)
                                 qsumm = qsumm + QSW(k, jww)
                                 tsum = tsum + QSW(k, jww)*T2(k, IWD(jww))
                                 csum(cn(1:nac)) = csum(cn(1:nac))             &
                                   & + QSW(k, jww)*C2(k, IWD(jww), cn(1:nac))
                             enddo
                             TTR(jtt) = tsum/qsumm
                             CTR(cn(1:nac), jtt) = csum(cn(1:nac))/qsumm
                         else
                             TTR(jtt) = TOUT(jb)
                             CTR(cn(1:nac), jtt) = cout(cn(1:nac), jb)
                         endif
                     elseif(LATERAL_PIPE(jp))then
                         tsum = 0.0
                         qsumm = 0.0
                         csum = 0.0
                         ilat(jww) = 1
                         jb = jbd
                         do k = KTW(jww), KBW(jww)
                             qsumm = qsumm + QSW(k, jww)
                             tsum = tsum + QSW(k, jww)*T2(k, IWD(jww))
                             csum(cn(1:nac)) = csum(cn(1:nac)) + QSW(k, jww)   &
                               & *C2(k, IWD(jww), cn(1:nac))
                         enddo
                         tinsum(jb) = (tinsum(jb)*qinsum(jb) + tsum)           &
                                    & /(qsumm + qinsum(jb))
                         cinsum(cn(1:nac), jb)                                 &
                           & = (cinsum(cn(1:nac), jb)*qinsum(jb)               &
                           & + csum(cn(1:nac)))/(qsumm + qinsum(jb))
                         qinsum(jb) = qsumm + qinsum(jb)
                     endif
                 endif
             else
                 jtt = jtt + 1
                 jww = jww + 1
                 IWD(jww) = IDPI(jp)
                 ITR(jtt) = IUPI(jp)
                 qtr(jtt) = -QPI(jp)
                 qwd(jww) = -QPI(jp)
                 KTWD(jww) = KTDPI(jp)
                 KBWD(jww) = KBDPI(jp)
                 EWD(jww) = EDPI(jp)
                 PLACE_QTR(jtt) = PUPIC(jp)==' DENSITY'
                 SPECIFY_QTR(jtt) = PUPIC(jp)==' SPECIFY'
                 if(SPECIFY_QTR(jtt))then
                     ELTRT(jtt) = ETUPI(jp)
                     ELTRB(jtt) = EBUPI(jp)
                 endif
                 JBTR(jtt) = jbu
                 JBWD(jww) = jbd
                 i = MAX(CUS(JBWD(jww)), IWD(jww))
                 jb = JBWD(jww)
                 jw = JWDPI(jp)
                 kt = KTWB(jw)
                 jwd = jww
                 call LATERAL_WITHDRAWAL
                                    !(JWW)
                 do k = KTW(jww), KBW(jww)
                     QSS(k, i) = QSS(k, i) - QSW(k, jww)
                 enddo
                 if(IDPI(jp)/=0)then
                     tsum = 0.0
                     qsumm = 0.0
                     csum = 0.0
                     do k = KTW(jww), KBW(jww)
                         qsumm = qsumm + QSW(k, jww)
                         tsum = tsum + QSW(k, jww)*T2(k, IWD(jww))
                         csum(cn(1:nac)) = csum(cn(1:nac)) + QSW(k, jww)       &
                           & *C2(k, IWD(jww), cn(1:nac))
                     enddo
                     TTR(jtt) = tsum/qsumm
                     CTR(cn(1:nac), jtt) = csum(cn(1:nac))/qsumm
                 endif
             endif
         enddo
     endif
     if(gates)then
         call GATE_FLOW
         do jg = 1, ngt
 
!******      Positive flows
 
             jlat = 0
             jbu = JBUGT(jg)
             jbd = JBDGT(jg)
             tdgon = .FALSE. ! cb 1/16/13
             jsg = jg
             nnsg = 1
             if(CAC(ndo)=='      ON' .AND. GASGTC(jg)=='      ON')             &
              & tdgon = .TRUE.
             if(QGT(jg)>=0.0)then
                 if(LATERAL_GATE(jg))then
 !           JLAT      = 1
                     jww = jww + 1
                     IWD(jww) = IUGT(jg)
                     qwd(jww) = QGT(jg)
                     KTWD(jww) = KTUGT(jg)
                     KBWD(jww) = KBUGT(jg)
                     EWD(jww) = EGT(jg)
                     if(DYNGTC(jg)=='     ZGT' .AND. gt2char=='EGT2ELEV')then               ! SW 2/25/11
                         if(EGT2(jg)/=0.0)EWD(jww) = EGT2(jg)
                     endif
                     JBWD(jww) = jbu
                     i = MAX(CUS(JBWD(jww)), IWD(jww))
                     jw = JWUGT(jg)
                     jb = JBWD(jww)
                     kt = KTWB(jw)
                     jwd = jww
                     call LATERAL_WITHDRAWAL
                                       !(JWW)
                     do k = KTW(jww), KBW(jww)
                         QSS(k, i) = QSS(k, i) - QSW(k, jww)
                     enddo
                     if(IDGT(jg)/=0)then
                         csum(cn(1:nac)) = 0.0
                         tsum = 0.0
                         qsumm = 0.0
                         jtt = jtt + 1                     !  SW 4/1/09
                         qtr(jtt) = QGT(jg)                !  SW 4/1/09
                         ITR(jtt) = IDGT(jg)               !  SW 4/1/09
                         PLACE_QTR(jtt) = PDGTC(jg)==' DENSITY'
                                                           !  SW 4/1/09
                         SPECIFY_QTR(jtt) = PDGTC(jg)==' SPECIFY'
                                                           !  SW 4/1/09
                         if(SPECIFY_QTR(jtt))then          !  SW 4/1/09
                             ELTRT(jtt) = ETDGT(jg)        !  SW 4/1/09
                             ELTRB(jtt) = EBDGT(jg)        !  SW 4/1/09
                         endif
                         JBTR(jtt) = jbd                   !  SW 4/1/09
                         do k = KTW(jww), KBW(jww)
                             qsumm = qsumm + QSW(k, jww)
                             tsum = tsum + QSW(k, jww)*T2(k, IWD(jww))
                             csum(cn(1:nac)) = csum(cn(1:nac)) + QSW(k, jww)   &
                               & *C2(k, IWD(jww), cn(1:nac))
                         enddo
                         if(qsumm==0.0)then
                             TTR(jtt) = 0.0
                             CTR(:, jtt) = 0.0
                         else
                             TTR(jtt) = tsum/qsumm
                             do jc = 1, nac
                                 CTR(cn(jc), jtt) = csum(cn(jc))/qsumm
                                 if(cn(jc)==ndo .AND. GASGTC(jg)               &
                                  & =='      ON' .AND. QGT(jg)>0.0)then
                                     TDG_GATE(jww, jg) = .TRUE.
                                     call TOTAL_DISSOLVED_GAS(0, PALT(id), 1,  &
                                       & jg, TTR(jtt), CTR(cn(jc), jtt))                ! O2
                                 endif
                                 if(cn(jc)==ngn2 .AND. GASGTC(jg)              &
                                  & =='      ON' .AND. QGT(jg)>0.0)then
                                     TDG_GATE(jww, jg) = .TRUE.
                                     call TOTAL_DISSOLVED_GAS(1, PALT(id), 1,  &
                                       & jg, TTR(jtt), CTR(cn(jc), jtt))                ! N2
                                 endif
                             enddo
                         endif
                     elseif(CAC(ndo)=='      ON' .AND. GASGTC(jg)              &
                          & =='      ON' .AND. QGT(jg)>0.0)then
                         TDG_GATE(jww, jg) = .TRUE.
                     endif
                 else
                     jss(jbu) = jss(jbu) + 1
                     KTSW(jss(jbu), jbu) = KTUGT(jg)
                     KBSW(jss(jbu), jbu) = KBUGT(jg)
                     jb = jbu
                     POINT_SINK(jss(jbu), jbu) = .TRUE.
                     id = IUGT(jg)
                     ESTR(jss(jbu), jbu) = EGT(jg)
                     if(DYNGTC(jg)=='     ZGT' .AND. gt2char=='EGT2ELEV')then               ! SW 2/25/11
                         if(EGT2(jg)/=0.0)ESTR(jss(jbu), jbu) = EGT2(jg)
                     endif
                     QSTR(jss(jbu), jbu) = QGT(jg)
                     kt = KTWB(JWUGT(jg))
                     jw = JWUGT(jg)
                     call DOWNSTREAM_WITHDRAWAL(jss(jbu))
                     QSUM(jb) = 0.0
                     tsum = 0.0
                     csum = 0.0
                     do k = kt, KB(id)
                         QSUM(jb) = QSUM(jb) + qout(k, jb)
                         tsum = tsum + qout(k, jb)*T2(k, id)
                         do jc = 1, nac
                             if(cn(jc)==ndo .AND. CAC(ndo)=='      ON' .AND.   &
                              & GASGTC(jg)=='      ON' .AND. QGT(jg)>0.0)then                                             ! MM 5/21/2009
                                 t2r4 = T2(k, id)
                                 cgas = C2(k, id, cn(jc))                                                                 ! MM 5/21/2009
                                 call TOTAL_DISSOLVED_GAS(0, PALT(id), 1, jg,  &
                                   & t2r4, cgas)
                                 csum(cn(jc)) = csum(cn(jc)) + qout(k, jb)*cgas
                 !ELSEIF (CN(JC) == NGN2 .AND. CAC(NGN2) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN   ! SW 10/27/15
                             elseif(cn(jc)==ngn2 .AND. GASGTC(jg)              &
                                  & =='      ON' .AND. QGT(jg)>0.0)then                           ! SW 10/27/15
                                 if(CAC(ngn2)=='      ON')then                                                                   ! cb 1/13/16
                                     t2r4 = T2(k, id)
                                     cgas = C2(k, id, cn(jc))                                                                 !
                                     call TOTAL_DISSOLVED_GAS(1, PALT(id), 1,  &
                                       & jg, t2r4, cgas)
                                     csum(cn(jc)) = csum(cn(jc)) + qout(k, jb) &
                                       & *cgas
                                 endif
                             else
                                 csum(cn(jc)) = csum(cn(jc)) + qout(k, jb)     &
                                   & *C2(k, id, cn(jc))
                             endif
                         enddo
                     enddo
                     if(QSUM(jb)/=0.0)then
                         TOUT(jb) = tsum/QSUM(jb)
                         cout(cn(1:nac), jb) = csum(cn(1:nac))/QSUM(jb)
                     endif
                     if(IDGT(jg)/=0 .AND. US(jbd)==IDGT(jg))then
                         qsumm = 0.0
                         tsum = 0.0
                         csum = 0.0
                         do k = kt, KB(id)
                             qsumm = qsumm + QNEW(k)
                             tsum = tsum + QNEW(k)*T2(k, id)
                             do jc = 1, nac
                                 if(cn(jc)==ndo .AND. CAC(ndo)                 &
                                  & =='      ON' .AND. GASGTC(jg)              &
                                  & =='      ON' .AND. QGT(jg)>0.0)then                                                     ! MM 5/21/2009
                                     t2r4 = T2(k, id)
                                     cgas = C2(k, id, cn(jc))                                                               ! MM 5/21/2009
                                     call TOTAL_DISSOLVED_GAS(0, PALT(id), 1,  &
                                       & jg, t2r4, cgas)                  ! O2
                                     csum(cn(jc)) = csum(cn(jc)) + QNEW(k)*cgas
                  !ELSEIF (CN(JC) == NGN2 .AND. CAC(NGN2) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN   ! SW 10/27/15
                                 elseif(cn(jc)==ngn2 .AND. GASGTC(jg)          &
                                       &=='      ON' .AND. QGT(jg)>0.0)then                         ! SW 10/27/15
                                     if(CAC(ngn2)=='      ON')then
                                         t2r4 = T2(k, id)
                                         cgas = C2(k, id, cn(jc))                                                             !
                                         call TOTAL_DISSOLVED_GAS(1, PALT(id), &
                                           & 1, jg, t2r4, cgas)             ! N2
                                         csum(cn(jc)) = csum(cn(jc)) + QNEW(k) &
                                           & *cgas
                                     endif
                                 else
                                     csum(cn(jc)) = csum(cn(jc)) + QNEW(k)     &
                                       & *C2(k, id, cn(jc))
                                 endif
                             enddo
                         enddo
                         if(qsumm/=0.0)then
                             tinsum(jbd) = (tsum + qinsum(jbd)*tinsum(jbd))    &
                               & /(qsumm + qinsum(jbd))
                             cinsum(cn(1:nac), jbd)                            &
                               & = (csum(cn(1:nac)) + qinsum(jbd)              &
                               & *cinsum(cn(1:nac), jbd))/(qsumm + qinsum(jbd))
                             qinsum(jbd) = qinsum(jbd) + qsumm
                         endif
                     elseif(IDGT(jg)/=0)then
                         jtt = jtt + 1
                         qtr(jtt) = QGT(jg)
                         ITR(jtt) = IDGT(jg)
                         PLACE_QTR(jtt) = PDGTC(jg)==' DENSITY'
                         SPECIFY_QTR(jtt) = PDGTC(jg)==' SPECIFY'
                         if(SPECIFY_QTR(jtt))then
                             ELTRT(jtt) = ETDGT(jg)
                             ELTRB(jtt) = EBDGT(jg)
                         endif
                         JBTR(jtt) = jbd
                         TTR(jtt) = TOUT(jb)
                         do jc = 1, nac
                             CTR(cn(jc), jtt) = cout(cn(jc), jb)
                             if(cn(jc)==ndo .AND. GASGTC(jg)=='      ON' .AND. &
                              & QGT(jg)>0.0)                                   &
                              & call TOTAL_DISSOLVED_GAS(0, PALT(id), 0, js,   &
                              & TTR(jtt), CTR(cn(jc), jtt))                                ! O2
                             if(cn(jc)==ngn2 .AND. GASGTC(jg)=='      ON' .AND.&
                              & QGT(jg)>0.0)                                   &
                              & call TOTAL_DISSOLVED_GAS(1, PALT(id), 0, js,   &
                              & TTR(jtt), CTR(cn(jc), jtt))                                ! N2
                         enddo
                     endif
                 endif
             elseif(QGT(jg)<0.0)then
                 jtt = jtt + 1
                 jww = jww + 1
                 IWD(jww) = IDGT(jg)
                 ITR(jtt) = IUGT(jg)
                 qtr(jtt) = -QGT(jg)
                 qwd(jww) = -QGT(jg)
                 KTWD(jww) = KTDGT(jg)
                 KBWD(jww) = KBDGT(jg)
                 EWD(jww) = EGT(jg)
                 PLACE_QTR(jtt) = PUGTC(jg)==' DENSITY'
                 SPECIFY_QTR(jtt) = PUGTC(jg)==' SPECIFY'
                 if(SPECIFY_QTR(jtt))then
                     ELTRT(jtt) = ETUGT(jg)
                     ELTRB(jtt) = EBUGT(jg)
                 endif
                 JBTR(jtt) = jbu
                 JBWD(jww) = jbd
                 i = MAX(CUS(JBWD(jww)), IWD(jww))
                 jw = JWDGT(jg)
                 jb = JBWD(jww)
                 kt = KTWB(jw)
                 jwd = jww
                 call LATERAL_WITHDRAWAL
                                  !(JWW)
                 do k = KTW(jww), KBW(jww)
                     QSS(k, i) = QSS(k, i) - QSW(k, jww)
                 enddo
                 if(IDGT(jg)/=0)then
                     csum(cn(1:nac)) = 0.0
                     tsum = 0.0
                     qsumm = 0.0
                     do k = KTW(jww), KBW(jww)
                         qsumm = qsumm + QSW(k, jww)
                         tsum = tsum + QSW(k, jww)*T2(k, IWD(jww))
                         csum(cn(1:nac)) = csum(cn(1:nac)) + QSW(k, jww)       &
                           & *C2(k, IWD(jww), cn(1:nac))
                     enddo
                     TTR(jtt) = tsum/qsumm
                     do jc = 1, nac
                         CTR(cn(jc), jtt) = csum(cn(jc))/qsumm
                         if(cn(jc)==ndo .AND. GASGTC(jg)=='      ON' .AND.     &
                          & QGT(jg)>0.0)then
                             TDG_GATE(jww, jg) = .TRUE.
                             call TOTAL_DISSOLVED_GAS(0, PALT(id), 1, jg,      &
                               & TTR(jtt), CTR(cn(jc), jtt))
                         endif
                         if(cn(jc)==ngn2 .AND. GASGTC(jg)=='      ON' .AND.    &
                          & QGT(jg)>0.0)then
                             TDG_GATE(jww, jg) = .TRUE.
                             call TOTAL_DISSOLVED_GAS(1, PALT(id), 1, jg,      &
                               & TTR(jtt), CTR(cn(jc), jtt))
                         endif
                     enddo
                 elseif(CAC(ndo)=='      ON' .AND. GASGTC(jg)=='      ON' .AND.&
                      & QGT(jg)>0.0)then
                     TDG_GATE(jww, jg) = .TRUE.
                 endif
             endif
         enddo
     endif
 
     tdgon = .FALSE.                      ! cb 1/17/13
     tributaries = jtt>0
     withdrawals = jww>0
 
 
     do jw = 1, nwb
         do jb = BS(jw), BE(jw)
             if(BR_INACTIVE(jb))then
                ! CONVERT INFLOWS TO TRIBS SET TO THE CUS(1) LOCATION
                 jtt = jtt + 1
                 ITR(jtt) = CUS(1)           ! HARDWIRED TO FIRST BRANCH
                 qtr(jtt) = QIN(jb)
                 TTR(jtt) = TIN(jb)
                 do jc = 1, nac
                     CTR(cn(jc), jtt) = CIN(cn(jc), jb)
                 enddo
               !  PLACE_QTR(JTT)   =  '   DISTR'
                 PLACE_QTR(jtt) = .FALSE.             !SR 01/22/2018
                 JBTR(jtt) = 1
             endif
         enddo
     enddo
 
 
 
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             if(BR_INACTIVE(jb))cycle
                                ! SW 6/12/2017
             iu = CUS(jb)
             id = DS(jb)
             if(EVAPORATION(jw))then
                 EVBR(jb) = 0.0
                 do i = iu, id
                     fw = AFW(jw) + BFW(jw)*WIND2(i)**CFW(jw)
                     if(RH_EVAP(jw))then
                         ea = EXP                                              &
                            & (2.3026*(7.5*TDEW(jw)/(TDEW(jw) + 237.3) + 0.6609&
                            & ))
                         es = EXP(2.3026*(7.5*T2(kt, i)/(T2(kt,i) + 237.3) +   &
                            & 0.6609))
                         if(TDEW(jw)<0.0)                                      &
                          & ea = EXP(2.3026*(9.5*TDEW(jw)/(TDEW(jw) + 265.5)   &
                          & + 0.6609))
                         if(T2(kt, i)<0.0)                                     &
                          & es = EXP(2.3026*(9.5*T2(kt, i)/(T2(kt,i) + 265.5)  &
                          & + 0.6609))
                         tairv = (TAIR(jw) + 273.0)/(1.0 - 0.378*ea/760.0)
                         dtv = (T2(kt, i) + 273.0)/(1.0 - 0.378*es/760.0)      &
                             & - tairv
                         dtvl = 0.0084*WIND2(i)**3
                         if(dtv<dtvl)dtv = dtvl
                         fw = (3.59*dtv**0.3333333 + 4.26*WIND2(i))
                     endif
                     tm = (T2(kt, i) + TDEW(jw))*0.5
                     vptg = 0.35 + 0.015*tm + 0.0012*tm*tm
                     EV(i) = vptg*(T2(kt, i) - TDEW(jw))*fw*BI(kt, i)*DLX(i)   &
                           & /2.45E9
                     if(EV(i)<0.0 .OR. ICE(i))EV(i) = 0.0
                     QSS(kt, i) = QSS(kt, i) - EV(i)
                     EVBR(jb) = EVBR(jb) + EV(i)
                 enddo
             endif
             if(PRECIPITATION(jw))then
                 QPRBR(jb) = 0.0
                 do i = iu, id
                     QPR(i) = PR(jb)*BI(kt, i)*DLX(i)
                     QPRBR(jb) = QPRBR(jb) + QPR(i)
                     QSS(kt, i) = QSS(kt, i) + QPR(i)
                 enddo
             endif
             if(tributaries)then
                 do jt = 1, jtt
 
!**********          Inflow fractions
 
                     if(jb==JBTR(jt))then
                         i = MAX(ITR(jt), iu)
                         qtrf(kt:KB(i), jt) = 0.0
                         if(PLACE_QTR(jt))then
 
!**************              Inflow layer
 
                             sstot = 0.0
                             do j = nsss, nsse
                                 sstot = sstot + CTR(j, jt)
                             enddo
                             rhotr = DENSITY(TTR(jt), CTR(ntds, jt), sstot)
                             k = kt
                             do while (rhotr>RHO(k, i) .AND. k<KB(i))
                                 k = k + 1
                             enddo
                             KTTR(jt) = k
                             KBTR(jt) = k
 
!**************              Layer inflows
 
                             vqtr = qtr(jt)*dlt
                             vqtri = vqtr
                             qtrfr = 1.0
                             incr = -1
                             do while (qtrfr>0.0)
                                 if(k<=KB(i))then
                                     v1 = VOL(k, i)
                                     if(vqtr>0.5*v1)then
                                         qtrf(k, jt) = 0.5*v1/vqtri
                                         qtrfr = qtrfr - qtrf(k, jt)
                                         vqtr = vqtr - qtrf(k, jt)*vqtri
                                         if(k==kt)then
                                         k = KBTR(jt)
                                         incr = 1
                                         endif
                                     else
                                         qtrf(k, jt) = qtrfr
                                         qtrfr = 0.0
                                     endif
                                     if(incr<0)KTTR(jt) = k
                                     if(incr>0)KBTR(jt) = MIN(KB(i), k)
                                     k = k + incr
                                 else
                                     qtrf(kt, jt) = qtrf(kt, jt) + qtrfr
                                     qtrfr = 0.0
                                 endif
                             enddo
                         else
                             if(SPECIFY_QTR(jt))then
                                 KTTR(jt) = 2
     !             DO WHILE (EL(KTTR(JT),I) > ELTRT(JT))
                                 do while (EL(KTTR(jt), i)>ELTRT(jt) .AND.     &
                                   & EL(KTTR(jt) + 1, i)>ELTRT(jt))                            ! SW 10/3/13
                                     KTTR(jt) = KTTR(jt) + 1
                                 enddo
                                 KBTR(jt) = kmx - 1
                                 do while (EL(KBTR(jt), i)<ELTRB(jt))
                                     KBTR(jt) = KBTR(jt) - 1
                                 enddo
                             else
                                 KTTR(jt) = kt
                                 KBTR(jt) = KB(i)
                             endif
                             KTTR(jt) = MAX(kt, KTTR(jt))
                             KBTR(jt) = MIN(KB(i), KBTR(jt))
                             if(KBTR(jt)<KTTR(jt))KBTR(jt) = KTTR(jt)
                             bhsum = 0.0
                             do k = KTTR(jt), KBTR(jt)
                                 bhsum = bhsum + BH2(k, i)
                             enddo
                             do k = KTTR(jt), KBTR(jt)
                                 qtrf(k, jt) = BH2(k, i)/bhsum
                             enddo
                         endif
                         do k = KTTR(jt), KBTR(jt)
                             QSS(k, i) = QSS(k, i) + qtr(jt)*qtrf(k, jt)
                         enddo
                     endif
                 enddo
             endif
             if(DIST_TRIBS(jb))then
                 akbr = 0.0
                 do i = iu, id
                     akbr = akbr + BI(kt, i)*DLX(i)
                 enddo
                 do i = iu, id
                     QDT(i) = QDTR(jb)*BI(kt, i)*DLX(i)/akbr
                     QSS(kt, i) = QSS(kt, i) + QDT(i)
                 enddo
             endif
             if(withdrawals)then
                 do jwd = 1, nwd
                     if(jb==JBWD(jwd))then
                         i = MAX(CUS(JBWD(jwd)), IWD(jwd))
                         call LATERAL_WITHDRAWAL
                                         !(JWD)
                         do k = KTW(jwd), KBW(jwd)
                             QSS(k, i) = QSS(k, i) - QSW(k, jwd)
                         enddo
                     endif
                 enddo
             endif
             if(UH_INTERNAL(jb))then
                 if(UHS(jb)/=DS(JBUH(jb)) .OR. DHS(JBUH(jb))/=US(jb))then
                     if(JBUH(jb)>=BS(jw) .AND. JBUH(jb)<=BE(jw))then
                         do k = kt, KB(iu - 1)
                             QSS(k, UHS(jb)) = QSS(k, UHS(jb)) - VOLUH2(k, jb) &
                               & /dlt
                         enddo
                     else
                         call UPSTREAM_FLOW
                     endif
                 endif
             endif
             if(DH_INTERNAL(jb))then
                 if(DHS(jb)/=US(JBDH(jb)) .OR. UHS(JBDH(jb))/=DS(jb))then
                     if(JBDH(jb)>=BS(jw) .AND. JBDH(jb)<=BE(jw))then
                         do k = kt, KB(id + 1)
                             QSS(k, CDHS(jb)) = QSS(k, CDHS(jb))               &
                               & + VOLDH2(k, jb)/dlt
                         enddo
                     else
                         call DOWNSTREAM_FLOW
                     endif
                 endif
             endif
         enddo
     enddo
 
!**  Compute tributary contribution to cross-shear
 
     if(tributaries)then
         do jw = 1, nwb
             do jb = BS(jw), BE(jw)
                 do jt = 1, jtt
                     if(jb==JBTR(jt))then
                         i = MAX(CUS(jb), ITR(jt))
                         do k = KTWB(jw), KBMIN(i)
                             uybr(k, i) = uybr(k, i) + ABS(qtr(jt))*qtrf(k, jt)
                         enddo
                     endif
                 enddo
             enddo
         enddo
     endif
 
     end subroutine HYDROINOUT
