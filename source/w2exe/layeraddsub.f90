!*==layeraddsub.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     subroutine LAYERADDSUB
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
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND) :: dummy, tmac, w1, w2, w3
     integer :: i_br_num, jbiz, kkb, ktmax
!
!*** End of declarations rewritten by SPAG
!
 
 
!***********************************************************************************************************************************
!**  Task 2.5: Layer - Segment Additions and Subtractions                     
!***********************************************************************************************************************************
 
!**  ** Water surface minimum thickness
 
     do jw = 1, nwb
         kt = KTWB(jw)
         ZMIN(jw) = -1000.0
         ktmax = 2                                                                                                  ! SR 10/17/05
         do jb = BS(jw), BE(jw)
             if(BR_INACTIVE(jb))then
                 if(DS(jb) - US(jb) + 1>=3)then
                     i_br_num = DS(jb) - 2
                 else
                     i_br_num = DS(jb) - 1
                                   ! FOR BRANCHES WITH LESS THAN 3 SEGMENTS
                 endif
 
                 if(CUS(JBDH(jb))<=DHS(jb) .AND. ELWS(DHS(jb))                 &
                  & >EL(KB(i_br_num), i_br_num))then                                       ! ***
                     BR_INACTIVE(jb) = .FALSE.
                     if(SNAPSHOT(jw))write(SNP(jw),                            &
                       &'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))')   &
                       &'   Branch Active: ', jb, ' at Julian day = ', jday,   &
                       &'   NIT = ', nit
                     if(SNAPSHOT(jw))write(SNP(jw), '(/17X,2(A,I0))')          &
                       &' Add segments ', DS(jb) - 1, ' through ', DS(jb)
                     CUS(jb) = DS(jb) - 1
                     do i = DS(jb) - 1, DS(jb)
                         Z(i) = Z(DHS(jb))
                         KTI(i) = KTI(DHS(jb))
                         H1(kt + 1, i) = H(kt + 1, jw)
                         AVH1(kt + 1, i) = (H1(kt + 1, i) + H1(kt + 2, i))*0.5
                         AVHR(kt + 1, i) = H1(kt + 1, i)                       &
                           & + (H1(kt + 1, i + 1) - H1(kt + 1, i))             &
                           & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                                               !SW 07/29/04
                         if(.NOT.TRAPEZOIDAL(jw))then
                             BH1(kt + 1, i) = B(kt + 1, i)*H(kt + 1, jw)
                             H1(kt, i) = H(kt, jw) - Z(i)
                             bi(kt:KB(i), i) = B(kt:KB(i), i)                                                                      ! SW 4/18/07
                             bi(kt, i) = B(KTI(i), i)
                             AVH1(kt, i) = (H1(kt, i) + H1(kt + 1, i))*0.5
                             AVHR(kt, i) = H1(kt, i)                           &
                               & + (H1(kt, i + 1) - H1(kt, i))                 &
                               & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                                           !SW 07/29/04
                             BH1(kt, i) = bi(kt, i)                            &
                               & *(EL(kt, i) - Z(i)*COSA(jb) -                 &
                               & EL(KTI(i) + 1, i))/COSA(jb)
                             if(KTI(i)>=KB(i))BH1(kt, i) = B(kt, i)*H1(kt, i)
                             do k = KTI(i) + 1, kt
                                 BH1(kt, i) = BH1(kt, i) + BH(k, i)
                             enddo
                         else
                             call GRID_AREA1(EL(kt, i) - Z(i), EL(kt + 1, i),  &
                               & BH1(kt, i), bi(kt, i))                                                                                      !SW 08/03/04
                             BH1(kt + 1, i) = 0.25*H(kt + 1, jw)               &
                               & *(BB(kt, i) + 2.*B(kt + 1, i) + BB(kt + 1, i))
                             H1(kt, i) = H(kt, jw) - Z(i)
                             AVH1(kt, i) = (H1(kt, i) + H1(kt + 1, i))*0.5
                             AVHR(kt, i) = H1(kt, i)                           &
                               & + (H1(kt, i + 1) - H1(kt, i))                 &
                               & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                                           !SW 07/29/04
                         endif
                         BKT(i) = BH1(kt, i)/H1(kt, i)
                         DEPTHB(kt, i) = H1(kt, i)
                         DEPTHM(kt, i) = H1(kt, i)*0.5
                         do k = kt + 1, KB(i)
                             DEPTHB(k, i) = DEPTHB(k - 1, i) + H1(k, i)
                             DEPTHM(k, i) = DEPTHM(k - 1, i)                   &
                               & + (H1(k - 1, i) + H1(k, i))*0.5
                         enddo
                     enddo
                     do i = DS(jb) - 1, DS(jb)
                         BHR1(kt + 1, i) = BH1(kt + 1, i)                      &
                           & + (BH1(kt + 1, i + 1) - BH1(kt + 1, i))           &
                           & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                                               !SW 07/29/04
                         BHR1(kt, i) = BH1(kt, i)                              &
                                     & + (BH1(kt, i + 1) - BH1(kt, i))         &
                                     & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                                     !SW 07/29/04
                     enddo
                     do i = DS(jb) - 1, DS(jb)
                         WIND2(i) = WIND2(DHS(jb))
                         if(DYNAMIC_SHADE(i))call SHADING
                         do k = kt, KB(i)
                             U(k, i) = 0.0
                             SDKV(k, i) = SDK(jw)
                             if(DXI(jw)>=0.0)then
                                 DX(k, i) = DXI(jw)
                             else
                                 DX(k, i) = ABS(U(k, i))*ABS(DXI(jw))*H(k, jw)
                                                                              ! SW 8/2/2017
                             endif
                             if(INTERNAL_WEIR(k, i))then
                                 DX(k, i) = 0.0
                                 U(k, i) = 0.0
                             endif
                             T1(k, i) = T1(k, DHS(jb))
                             T2(k, i) = T1(k, DHS(jb))
                             SU(k, i) = 0.0
                             C1(k, i, cn(1:nac)) = C1(k, DHS(jb), cn(1:nac))
                             C2(k, i, cn(1:nac)) = C1(k, DHS(jb), cn(1:nac))
                             do je = 1, nep
                                 EPD(k, i, je) = 0.01
                                 EPC(k, i, je) = 0.01/H1(k, i)
                             enddo
                             CMBRT(cn(1:nac), jb) = CMBRT(cn(1:nac), jb)       &
                               & + C1(k, DHS(jb), cn(1:nac))*DLX(i)*BH1(k, i)
                             EBRI(jb) = EBRI(jb) + T1(k, DHS(jb))*DLX(i)       &
                                      & *BH1(k, i)
                         enddo
                         do k = kt, KB(i) - 1
                             AZ(k, i) = AZ(k, DHS(jb))
                             TKE(k, i, 1) = TKE(k, DHS(jb), 1)    !sg 10/4/07
                             TKE(k, i, 2) = TKE(k, DHS(jb), 2)
                                                          !sg 10/4/07
                             SAZ(k, i) = AZ(k, DHS(jb))
                             if(INTERNAL_WEIR(k, i))then
                                 AZ(k, i) = 0.0
                                 TKE(k, i, 1) = 0.0
                                                   !sg 10/4/07
                                 TKE(k, i, 2) = 0.0
                                               !sg 10/4/07
                                 SAZ(k, i) = 0.0
                             endif
                         enddo
                     enddo
                 endif
             endif
             if(.NOT.BR_INACTIVE(jb))then
                 do i = CUS(jb), DS(jb)
                     if(KB(i)>ktmax)ktmax = KB(i)                                                                   ! SR 10/17/05
                     if(Z(i)>ZMIN(jw))then
                         IZMIN(jw) = i
                         jbiz = jb
                     endif
                     ZMIN(jw) = MAX(ZMIN(jw), Z(i))
                 enddo
             endif
         enddo
         add_layer = ZMIN(jw)< - 0.85*H(kt - 1, jw) .AND. kt/=2
         sub_layer = ZMIN(jw)>0.60*H(kt, jw) .AND. kt<ktmax                                                          ! SR 10/17/05
         if(KTWB(jw)==kmx - 1 .AND. SLOPE(jbiz)>0.0 .AND. sub_layer .AND.      &
          & ONE_LAYER(IZMIN(jw)))then
             if(ZMIN(jw)>0.99*H(kt, jw))then
                 write(wrn, '(A,I0,2(A,F0.3))')'Low water in segment ',        &
                     & IZMIN(jw), ' water surface deviation' // ' = ', ZMIN(jw)&
                     & , ' at day ', jday
                 warning_open = .TRUE.
             endif
             sub_layer = .FALSE.
         endif
 
         if(add_layer==.TRUE. .OR. sub_layer==.TRUE.)then
             LAYERCHANGE(jw) = .TRUE.
         else
             LAYERCHANGE(jw) = .FALSE.
         endif
 
!****    Add layers
 
         do while (add_layer)
             if(SNAPSHOT(jw))write(SNP(jw),                                    &
                           &'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))'&
                          & )'   Add layer ', kt - 1, ' at Julian day = ',     &
                           & jday, '   NIT = ', nit, ' IZMIN =', IZMIN(jw)                                              ! SW 1/23/06
 
!******      Variable initialization
 
             KTWB(jw) = KTWB(jw) - 1
             kt = KTWB(jw)
             ilayer = 0
 
  ! RECOMPUTE INTERNAL WEIR FOR FLOATING WEIR
             if(weir_calc)then
                          !  SW 3/16/18
                 do jwr = 1, niw
                     if(IWR(jwr)>=US(BS(jw)) .AND. IWR(jwr)<=DS(BE(jw)))then
                         if(EKTWR(jwr)==0.0)then
                             KTWR(jwr) = KTWB(jw)
                         else
                             KTWR(jwr) = INT(EKTWR(jwr))
                         endif
                         if(EKBWR(jwr)<=0.0)then
                             do k = KTWR(jwr), KB(IWR(jwr))
                                 if(DEPTHB(k, IWR(jwr))>ABS(EKBWR(jwr)))then
                                     KBWR(jwr) = k
                                     exit
                                 endif
                             enddo
                         else
                             KBWR(jwr) = INT(EKBWR(jwr))
                         endif
 
                         do k = 2, kmx - 1
                             if((k>=KTWR(jwr) .AND. k<=KBWR(jwr)))             &
                              & INTERNAL_WEIR(k, IWR(jwr)) = .TRUE.
                         enddo
                     endif
                 enddo
             endif
 
 
             do jb = BS(jw), BE(jw)
                 if(BR_INACTIVE(jb))cycle
                                    ! SW 6/12/2017
                 iu = CUS(jb)
                 id = DS(jb)
                 do i = iu - 1, id + 1
                     Z(i) = H(kt, jw) + Z(i)
                     H1(kt, i) = H(kt, jw) - Z(i)
                     H1(kt + 1, i) = H(kt + 1, jw)
                     H2(kt + 1, i) = H(kt + 1, jw)
                     AVH1(kt, i) = (H1(kt, i) + H1(kt + 1, i))*0.5D0
                     AVH1(kt + 1, i) = (H1(kt + 1, i) + H1(kt + 2, i))*0.5D0
                     if(.NOT.TRAPEZOIDAL(jw))then
                         BH1(kt, i) = BH1(kt + 1, i) - BNEW(kt + 1, i)         &
                                    & *H1(kt + 1, i)                                         ! SW 1/23/06
                         BH1(kt + 1, i) = BNEW(kt + 1, i)*H1(kt + 1, i)                      ! SW 1/23/06
                     else
                         call GRID_AREA1(EL(kt, i) - Z(i), EL(kt + 1, i),      &
                           & BH1(kt, i), dummy)                                                                                  !SW 08/03/04
                         BH1(kt + 1, i) = 0.25D0*H1(kt + 1, jw)                &
                           & *(BB(kt, i) + 2.D0*B(kt + 1, i) + BB(kt + 1, i))
                     endif
                     VOL(kt, i) = BH1(kt, i)*DLX(i)
                     VOL(kt + 1, i) = BH1(kt + 1, i)*DLX(i)
                     BKT(i) = BH1(kt, i)/H1(kt, i)
                     DEPTHB(kt, i) = H1(kt, i)
                     DEPTHM(kt, i) = H1(kt, i)*0.5
                     bi(kt:KB(i), i) = B(kt:KB(i), i)
                                              ! SW 8/26/05
                     bi(kt, i) = B(KTI(i), i)
                     T1(kt, i) = T1(kt + 1, i)
                     SDKV(kt, i) = SDKV(kt + 1, i)
                                              ! SW 1/18/08
                     SED(kt, i) = SED(kt + 1, i)
                                              ! SW 1/18/08
                     SEDN(kt, i) = SEDN(kt + 1, i)
                                              ! SW 1/18/08
                     SEDP(kt, i) = SEDP(kt + 1, i)
                                              ! SW 1/18/08
                     SEDC(kt, i) = SEDC(kt + 1, i)
                                              ! SW 1/18/08
!                    RHO(KT,I)     =
!                    DENSITY(T1(KT,I),MAX(TDS(KT,I),0.0),MAX(TISS(KT,I),0.0)) 
!                    ! SR 5/15/06
                     if(SDFIRSTADD(kt, i))then
                         SED1(kt, i) = SED1IC(kt, i)
                                       ! cb 6/17/17
                         SED2(kt, i) = SED2IC(kt, i)
                                       ! cb 6/17/17
                         SDFIRSTADD(kt, i) = .FALSE.
                     endif
                     do k = kt + 1, kmx
                         DEPTHB(k, i) = DEPTHB(k - 1, i) + H1(k, i)
                         DEPTHM(k, i) = DEPTHM(k - 1, i)                       &
                                      & + (H1(k - 1, i) + H1(k, i))*0.5D0
                     enddo
                     C1(kt, i, cn(1:nac)) = C1(kt + 1, i, cn(1:nac))
                     CSSK(kt, i, cn(1:nac)) = CSSK(kt + 1, i, cn(1:nac))
                     KF(kt, i, kfcn(1:NAF(jw), jw))                            &
                       & = KF(kt + 1, i, kfcn(1:NAF(jw), jw))
                     KFS(kt, i, kfcn(1:NAF(jw), jw))                           &
                       & = KFS(kt + 1, i, kfcn(1:NAF(jw), jw))                              !KF(KT+1,I,KFCN(1:NAF(JW),JW))   CODE ERROR FIX SW 10/24/2017
                     KF(kt + 1, i, kfcn(1:NAF(jw), jw)) = 0.0
                     KFS(kt + 1, i, kfcn(1:NAF(jw), jw)) = 0.0
                     if(kt>=KBI(i))then  ! CB 5/24/06
                         ADX(kt + 1, i) = 0.0
                                        ! CB 5/15/06
                         C1(kt + 1, i, cn(1:nac)) = 0.0
                                        ! CB 5/15/06
                         CSSK(kt + 1, i, cn(1:nac)) = 0.0
                                        ! CB 5/15/06
                     endif              ! CB 5/15/06
                     do je = 1, nep
                         if(kt<KBI(i))then
                                     ! CB 4/28/06
                             EPD(kt, i, je) = EPD(kt + 1, i, je)
                             EPM(kt, i, je) = EPD(kt, i, je)                   &
                               & *(bi(kt, i) - B(kt + 1, i) + 2.0*H1(kt, i))   &
                               & *DLX(i)                                                         ! SR 5/15/06
                             EPM(kt + 1, i, je) = EPM(kt + 1, i, je)           &
                               & - EPM(kt, i, je)
                             EPC(kt, i, je) = EPM(kt, i, je)/VOL(kt, i)
                             EPC(kt + 1, i, je) = EPM(kt + 1, i, je)           &
                               & /VOL(kt + 1, i)
                         else
                             EPD(kt, i, je) = EPD(kt + 1, i, je)
                                              ! SW 5/15/06
                             EPM(kt, i, je) = EPM(kt + 1, i, je)
                             EPC(kt, i, je) = EPC(kt + 1, i, je)
                             EPD(kt + 1, i, je) = 0.0
                             EPM(kt + 1, i, je) = 0.0
                             EPC(kt + 1, i, je) = 0.0
                         endif      ! CB 4/28/06
                     enddo
                 enddo
 
!********macrophytes...
                 do i = iu, id
                     jt = KTI(i)
                     je = KB(i)
                     do j = jt, je
                         if(j<kt)then
                             colb = EL(j + 1, i)
                         else
                             colb = EL(kt + 1, i)
                         endif
         !     COLDEP=ELWS(I)-COLB
                         coldep = EL(kt, i) - Z(i)*COSA(jb) - colb
                                                   ! cb 3/7/16
         !     MACT(J,KT,I)=MACT(J,KT+1,I)
                         if(macrophyte_on)MACT(j, kt, i) = MACT(j, kt + 1, i)
                                                              ! SW 9/28/13
                         do m = 1, nmc
                             if(MACROPHYTE_CALC(jw, m))then
                                 MACRC(j, kt, i, m) = MACRC(j, kt + 1, i, m)
                                 MACRM(j, kt, i, m) = MACRC(j, kt, i, m)       &
                                   & *CW(j, i)*coldep*DLX(i)
                             endif
                         enddo
                     enddo
 
                     jt = kt + 1
                     je = KB(i)
                     do j = jt, je
                         do m = 1, nmc
                             if(MACROPHYTE_CALC(jw, m))MACRM(j, kt + 1, i, m)  &
                              & = MACRC(j, kt + 1, i, m)*CW(j, i)*H(kt + 1, jw)&
                              & *DLX(i)
                         enddo
                     enddo
                 enddo
                 do i = iu - 1, id
                     AVHR(kt + 1, i) = H1(kt + 1, i)                           &
                                     & + (H1(kt + 1, i + 1) - H1(kt + 1, i))   &
                                     & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                           !SW 07/29/04 (H1(KT+1,I+1) +H1(KT+1,I))*0.5
                     AVHR(kt, i) = H1(kt, i) + (H1(kt, i + 1) - H1(kt, i))     &
                                 & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                               !SW 07/29/04 (H1(KT,I+1)   +H1(KT,I))*0.5
                     BHR1(kt, i) = BH1(kt, i) + (BH1(kt, i + 1) - BH1(kt, i))  &
                                 & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                               !SW 07/29/04  (BH1(KT,I+1)  +BH1(KT,I))*0.5
                     BHR1(kt + 1, i) = BH1(kt + 1, i)                          &
                                     & + (BH1(kt + 1, i + 1) - BH1(kt + 1, i)) &
                                     & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                           !SW 07/29/04  (BH1(KT+1,I+1)+BH1(KT+1,I))*0.5
                     U(kt, i) = U(kt + 1, i)
                 enddo
                 do i = iu, id
                     if(ONE_LAYER(i))then
                         W(kt, i) = 0.0
                     else
                         w1 = W(kt + 1, i)*BB(kt + 1, i)
                         w2 = (BHR1(kt + 1, i)*U(kt + 1, i)                    &
                            & - BHR1(kt + 1, i - 1)*U(kt + 1, i - 1))/DLX(i)
                         w3 = ( - QSS(kt + 1, i)*BH1(kt + 1, i)                &
                            & /(BH1(kt + 1, i) + BH1(kt, i)))/DLX(i)
                         W(kt, i) = (w1 + w2 + w3)/BB(kt, i)
                     endif
                 enddo
                 if(UP_HEAD(jb))then
                     bhsum = BHR1(kt, iu - 1) + BHR1(kt + 1, iu - 1)
                     QUH1(kt, jb) = QUH1(kt + 1, jb)*BHR1(kt, iu - 1)/bhsum
                     QUH1(kt + 1, jb) = QUH1(kt + 1, jb)*BHR1(kt + 1, iu - 1)  &
                                      & /bhsum
                     TSSUH1(kt, jb) = TSSUH1(kt + 1, jb)*BHR1(kt, iu - 1)/bhsum
                     TSSUH1(kt + 1, jb) = TSSUH1(kt + 1, jb)                   &
                       & *BHR1(kt + 1, iu - 1)/bhsum
                     CSSUH1(kt, cn(1:nac), jb) = CSSUH1(kt + 1, cn(1:nac), jb) &
                       & *BHR1(kt, iu - 1)/bhsum
                     CSSUH1(kt + 1, cn(1:nac), jb)                             &
                       & = CSSUH1(kt + 1, cn(1:nac), jb)*BHR1(kt + 1, iu - 1)  &
                       & /bhsum
                 endif
                 if(DN_HEAD(jb))then
                     bhsum = BHR1(kt, id) + BHR1(kt + 1, id)
                     QDH1(kt, jb) = QDH1(kt + 1, jb)*BHR1(kt, id)/bhsum
                     QDH1(kt + 1, jb) = QDH1(kt + 1, jb)*BHR1(kt + 1, id)/bhsum
                     TSSDH1(kt, jb) = TSSDH1(kt + 1, jb)*BHR1(kt, id)/bhsum
                     TSSDH1(kt + 1, jb) = TSSDH1(kt + 1, jb)*BHR1(kt + 1, id)  &
                       & /bhsum
                     CSSDH1(kt, cn(1:nac), jb) = CSSDH1(kt + 1, cn(1:nac), jb) &
                       & *BHR1(kt, id)/bhsum
                     CSSDH1(kt + 1, cn(1:nac), jb)                             &
                       & = CSSDH1(kt + 1, cn(1:nac), jb)*BHR1(kt + 1, id)/bhsum
                 endif
                 do i = iu, id - 1
                     if(DXI(jw)>=0.0)then
                         DX(kt, i) = DXI(jw)
                     else
                         DX(kt, i) = ABS(U(kt, i))*ABS(DXI(jw))*H(k, jw)
                                                                        ! SW 8/2/2017
                     endif
                     if(INTERNAL_WEIR(kt, i))DX(kt, i) = 0.0
                 enddo
                 iut = iu
                 idt = id - 1
                 if(UP_HEAD(jb))iut = iu - 1
                 if(DN_HEAD(jb))idt = id
                 do i = iut, idt
                     AZ(kt, i) = azmin
                     TKE(kt, i, 1) = 1.25E-7
                                      !sg 10/4/07
                     TKE(kt, i, 2) = 1.0E-9
                                      !sg 10/4/07
                     SAZ(kt, i) = azmin
                     if(INTERNAL_WEIR(kt, i))then
                         AZ(kt, i) = 0.0
                         TKE(kt, i, 1) = 0.0
                                      !sg  10/4/07
                         TKE(kt, i, 2) = 0.0
                                      !sg  10/4/07
                         SAZ(kt, i) = 0.0
                     endif
                 enddo
                 if(constituents)then
                     call TEMPERATURE_RATES
                     call KINETIC_RATES
                 endif
 
!********        Upstream active segment
 
                 iut = US(jb)
                 if(SLOPE(jb)==0.0)then
                     do i = US(jb), DS(jb)
                         if(KB(i) - kt<NL(jb) - 1)iut = i + 1
                     enddo
                 else
                     do i = US(jb) - 1, DS(jb) + 1
                         if(KB(i)>KBI(i))then
                             BNEW(KB(i), i) = B(KB(i), i)                                      ! SW 1/23/06   ! SW 3/2/05
                             DX(KB(i), i) = 0.0
                             KB(i) = KB(i) - 1
                             ilayer(i) = 1
                             U(KB(i) + 1, i) = 0.0                                             ! SW 1/23/06
                             write(wrn, '(2(A,I8),A,F0.3)')                    &
                                  &'Raising bottom layer at segment ', i,      &
                                  &' at iteration ', nit, ' at Julian day ',   &
                                 & jday
                             warning_open = .TRUE.
                         endif
                     enddo           ! SW 1/23/06
                     do i = US(jb) - 1, DS(jb) + 1
                                     ! SW 1/23/06
!                        IF (KB(I)-KT < NL(JB)-1) IUT = I+1    ! SW 1/23/06
!                        IF (I /= DS(JB)+1) KBMIN(I)   = MIN(KB(I),KB(I+1))   
!                        ! SW 1/23/06                                 ! SW
!                        3/2/05
                         if(i/=US(jb) - 1)KBMIN(i - 1) = MIN(KB(i - 1), KB(i))         ! SW 1/23/06
                         if(KBI(i)<KB(i))then
                             BKT(i) = BH1(kt, i)                               &
                                    & /(H1(kt, i) - (EL(KBI(i) + 1, i)         &
                                    & - EL(KB(i) + 1, i)))                    ! SW 1/23/06
                             DEPTHB(KTWB(jw), i)                               &
                               & = (H1(KTWB(jw), i) - (EL(KBI(i) + 1, i)       &
                               & - EL(KB(i) + 1, i)))                                 ! SW 1/23/06
                             DEPTHM(KTWB(jw), i)                               &
                               & = (H1(KTWB(jw), i) - (EL(KBI(i) + 1, i)       &
                               & - EL(KB(i) + 1, i)))*0.5                                 ! SW 1/23/06
                             AVHR(kt, i)                                       &
                               & = (H1(kt, i) - (EL(KBI(i) + 1, i) - EL(KB(i)  &
                               & + 1, i)))                                     &
                               & + (H1(kt, i + 1) - (EL(KBI(i) + 1, i + 1)     &
                               & - EL(KB(i) + 1, i + 1)) - H1(kt, i)           &
                               & + (EL(KBI(i) + 1, i) - EL(KB(i) + 1, i)))     &
                               & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                                                                                ! SW 1/23/06
                         endif
                     enddo
                     do i = US(jb) - 1, DS(jb) + 1
                                 ! SW 1/23/06
                         do k = KBMIN(i) + 1, KB(i)
                             U(k, i) = 0.0
                         enddo
                     enddo
                  ! SW 11/9/07
 
                     do i = US(jb), DS(jb)
                               ! SW 11/9/07
                         if(ilayer(i)==1 .AND. ilayer(i + 1)==0)then
                             bhrsum = 0.0
                             Q(i) = 0.0
                             do k = kt, KBMIN(i)
                                 if(.NOT.INTERNAL_WEIR(k, i))then
                                     bhrsum = bhrsum + BHR1(k, i)
                                     Q(i) = Q(i) + U(k, i)*BHR1(k, i)
                                 endif
                             enddo
                             do k = kt, KBMIN(i)
                                 if(INTERNAL_WEIR(k, i))then
                                     U(k, i) = 0.0
                                 else
                                     U(k, i) = U(k, i) + (QC(i) - Q(i))/bhrsum
                                 endif
                             enddo
                         elseif(ilayer(i)==1 .AND. ilayer(i - 1)==0)then
                             bhrsum = 0.0
                             Q(i - 1) = 0.0
                             do k = kt, KBMIN(i - 1)
                                 if(.NOT.INTERNAL_WEIR(k, i - 1))then
                                     bhrsum = bhrsum + BHR1(k, i - 1)
                                     Q(i - 1) = Q(i - 1) + U(k, i - 1)         &
                                       & *BHR1(k, i - 1)
                                 endif
                             enddo
                             do k = kt, KBMIN(i - 1)
                                 if(INTERNAL_WEIR(k, i - 1))then
                                     U(k, i - 1) = 0.0
                                 else
                                     U(k, i - 1) = U(k, i - 1)                 &
                                       & + (QC(i - 1) - Q(i - 1))/bhrsum
                                 endif
                             enddo
                         endif
 
                     enddo
                 endif
 
!********        Segment addition
 
                 if(iut/=iu)then
                     if(SNAPSHOT(jw))write(SNP(jw), '(/17X,2(A,I0))')          &
                       &' Add segments ', iut, ' through ', iu - 1
                     do i = iut - 1, iu - 1
                         Z(i) = Z(iu)
                         KTI(i) = KTI(iu)
                         H1(kt + 1, i) = H(kt + 1, jw)
                         AVH1(kt + 1, i) = (H1(kt + 1, i) + H1(kt + 2, i))*0.5
                         AVHR(kt + 1, i) = H1(kt + 1, i)                       &
                           & + (H1(kt + 1, i + 1) - H1(kt + 1, i))             &
                           & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                                   !SW 07/29/04
                         if(.NOT.TRAPEZOIDAL(jw))then
                             BH1(kt + 1, i) = B(kt + 1, i)*H(kt + 1, jw)
                             H1(kt, i) = H(kt, jw) - Z(i)
                             bi(kt:KB(i), i) = B(kt:KB(i), i)                                                          ! SW 4/18/07
                             bi(kt, i) = B(KTI(i), i)
                             AVH1(kt, i) = (H1(kt, i) + H1(kt + 1, i))*0.5
                             AVHR(kt, i) = H1(kt, i)                           &
                               & + (H1(kt, i + 1) - H1(kt, i))                 &
                               & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                               !SW 07/29/04
                             BH1(kt, i) = bi(kt, i)                            &
                               & *(EL(kt, i) - Z(i)*COSA(jb) -                 &
                               & EL(KTI(i) + 1, i))/COSA(jb)
                             if(KTI(i)>=KB(i))BH1(kt, i) = B(kt, i)*H1(kt, i)
                             do k = KTI(i) + 1, kt
                                 BH1(kt, i) = BH1(kt, i) + BH(k, i)
                             enddo
                         else
                             call GRID_AREA1(EL(kt, i) - Z(i), EL(kt + 1, i),  &
                               & BH1(kt, i), bi(kt, i))                                                                          !SW 08/03/04
                             BH1(kt + 1, i) = 0.25*H(kt + 1, jw)               &
                               & *(BB(kt, i) + 2.*B(kt + 1, i) + BB(kt + 1, i))
                             H1(kt, i) = H(kt, jw) - Z(i)
                             AVH1(kt, i) = (H1(kt, i) + H1(kt + 1, i))*0.5
                             AVHR(kt, i) = H1(kt, i)                           &
                               & + (H1(kt, i + 1) - H1(kt, i))                 &
                               & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                               !SW 07/29/04
                         endif
                         BKT(i) = BH1(kt, i)/H1(kt, i)
                         DEPTHB(kt, i) = H1(kt, i)
                         DEPTHM(kt, i) = H1(kt, i)*0.5
                         do k = kt + 1, KB(i)
                             DEPTHB(k, i) = DEPTHB(k - 1, i) + H1(k, i)
                             DEPTHM(k, i) = DEPTHM(k - 1, i)                   &
                               & + (H1(k - 1, i) + H1(k, i))*0.5
                         enddo
                     enddo
                     do i = iut - 1, iu - 1
                         BHR1(kt + 1, i) = BH1(kt + 1, i)                      &
                           & + (BH1(kt + 1, i + 1) - BH1(kt + 1, i))           &
                           & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                                   !SW 07/29/04
                         BHR1(kt, i) = BH1(kt, i)                              &
                                     & + (BH1(kt, i + 1) - BH1(kt, i))         &
                                     & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                         !SW 07/29/04
                     enddo
                     do i = iut, iu - 1
              !ICE(I)   = ICE(IU)   ! SW 9/29/15
              !ICETH(I) = ICETH(IU)
                         WIND2(i) = WIND2(iu)
                         if(DYNAMIC_SHADE(i))call SHADING
                         do k = kt, KB(i)
                             U(k, i) = U(k, iu)
                             if(DXI(jw)>=0.0)then
                                 DX(k, i) = DXI(jw)
                             else
                                 DX(k, i) = ABS(U(k, i))*ABS(DXI(jw))*H(k, jw) ! SW 8/2/2017
                             endif
                             SDKV(k, i) = SDK(jw)
                                      ! SW 1/18/08
                             if(INTERNAL_WEIR(k, i))then
                                 DX(k, i) = 0.0
                                 U(k, i) = 0.0
                             endif
                             T1(k, i) = T1(k, iu)
                             T2(k, i) = T1(k, iu)
                             SU(k, i) = U(k, iu)
                             C1(k, i, cn(1:nac)) = C1(k, iu, cn(1:nac))
                             C2(k, i, cn(1:nac)) = C1(k, iu, cn(1:nac))
                             do je = 1, nep
                                 EPD(k, i, je) = 0.01
                                 EPC(k, i, je) = 0.01/H1(k, i)
                             enddo
                             CMBRT(cn(1:nac), jb) = CMBRT(cn(1:nac), jb)       &
                               & + C1(k, iu, cn(1:nac))*DLX(i)*BH1(k, i)
                             EBRI(jb) = EBRI(jb) + T1(k, iu)*DLX(i)*BH1(k, i)
                         enddo
                         do k = kt, KB(i) - 1
                             AZ(k, i) = AZ(k, iu)
                             TKE(k, i, 1) = TKE(k, iu, 1)
                                                 !sg 10/4/07
                             TKE(k, i, 2) = TKE(k, iu, 2)
                                         !sg 10/4/07
                             SAZ(k, i) = AZ(k, iu)
                             if(INTERNAL_WEIR(k, i))then
                                 AZ(k, i) = 0.0
                                 TKE(k, i, 1) = 0.0
                                       !sg 10/4/07
                                 TKE(k, i, 2) = 0.0
                                   !sg 10/4/07
                                 SAZ(k, i) = 0.0
                             endif
                         enddo
            !END DO
!*********macrophytes
                         do m = 1, nmc
                             if(MACROPHYTE_CALC(jw, m))then
                                 jt = KTI(i)
                                 je = KB(i)
                                 do j = jt, je
                                     if(j<kt)then
                                         colb = EL(j + 1, i)
                                     else
                                         colb = EL(kt + 1, i)
                                     endif
                    !COLDEP=ELWS(I)-COLB
                                     coldep = EL(kt, i) - Z(i)*COSA(jb) - colb
                                                        ! cb 3/7/16
                    !MACRC(J,KT,I,M)=MACWBCI(JW,M)
                                     if(ISO_MACROPHYTE(jw, m))                 &
                                      & MACRC(j, kt, i, m) = MACWBCI(jw, m)        ! cb 3/7/16
                                     if(VERT_MACROPHYTE(jw, m))                &
                                      & MACRC(j, kt, i, m) = 0.1
                                     if(LONG_MACROPHYTE(jw, m))                &
                                      & MACRC(j, kt, i, m) = 0.1
                                     MACRM(j, kt, i, m) = MACRC(j, kt, i, m)   &
                                       & *CW(j, i)*coldep*DLX(i)
                                     MACMBRT(jb, m) = MACMBRT(jb, m)           &
                                       & + MACRM(j, kt, i, m)
                                 enddo
                                 do k = kt + 1, KB(i)
                                     jt = k
                                     je = KB(i)
                                     do j = jt, je
                      !MACRC(J,K,I,M)=MACWBCI(JW,M)
                                         if(ISO_MACROPHYTE(jw, m))             &
                                           & MACRC(j, k, i, m) = MACWBCI(jw, m)     ! cb 3/7/16
                                         if(VERT_MACROPHYTE(jw, m))            &
                                           & MACRC(j, k, i, m) = 0.1
                                         if(LONG_MACROPHYTE(jw, m))            &
                                           & MACRC(j, k, i, m) = 0.1
                                         MACRM(j, k, i, m) = MACRC(j, k, i, m) &
                                           & *CW(j, i)*H2(k, i)*DLX(i)
                                         MACMBRT(jb, m) = MACMBRT(jb, m)       &
                                           & + MACRM(j, k, i, m)
                                     enddo
                                 enddo
                             endif
                         enddo
                     enddo
                    ! cb 3/7/16  moved enddo to include macrophytes
                     U(KB(iut):KB(iu), iu - 1) = 0.0
                     SU(KB(iut):KB(iu), iu - 1) = 0.0
                     ADX(KB(iut):KB(iu), iu) = 0.0
                     iu = iut
                     CUS(jb) = iu
                     if(UH_EXTERNAL(jb))KB(iu - 1) = KB(iu)
                     if(UH_INTERNAL(jb))then
                         if(JBUH(jb)>=BS(jw) .AND. JBUH(jb)<=BE(jw))then
                             KB(iu - 1) = MIN(KB(UHS(jb)), KB(iu))
                         else
                             do kkb = kt, kmx
                                 if(EL(kkb, iu)<=EL(KB(UHS(jb)), UHS(jb)))exit
                             enddo
                             KB(iu - 1) = MIN(kkb, KB(iu))
                         endif
                     endif
                     if(UP_HEAD(jb))then
                         AZ(kt:KB(iu - 1) - 1, iu - 1) = azmin
                         TKE(kt:KB(iu - 1) - 1, iu - 1, 1) = 1.25E-7
                                                      !SG 10/4/07
                         TKE(kt:KB(iu - 1) - 1, iu - 1, 2) = 1.0E-9
                                                      !SG 10/4/07
                         SAZ(kt:KB(iu - 1) - 1, iu - 1) = azmin
                     endif
                 endif
                 if(constituents)then
                     call TEMPERATURE_RATES
                     call KINETIC_RATES
                 endif
 
!********        Total active cells and single layers
 
                 do i = iu, id
                     ntac = ntac + 1
                     ONE_LAYER(i) = KTWB(jw)==KB(i)
                 enddo
                 ntacmx = MAX(ntac, ntacmx)
             enddo
             call INTERPOLATION_MULTIPLIERS
 
!******      Additional layers
 
             ZMIN(jw) = -1000.0
             do jb = BS(jw), BE(jw)
                 if(BR_INACTIVE(jb))cycle
                                    ! SW 6/12/2017
                 do i = CUS(jb), DS(jb)
                     ZMIN(jw) = MAX(ZMIN(jw), Z(i))
                 enddo
             enddo
             add_layer = ZMIN(jw)< - 0.80*H(kt - 1, jw) .AND. kt/=2
         enddo
 
!****    Subtract layers
 
         do while (sub_layer)
             if(SNAPSHOT(jw))write(SNP(jw),                                    &
                        &'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,1x,13("*"))'&
                       & )'Subtract layer ', kt, ' at Julian day = ', jday,    &
                         &' NIT = ', nit, ' IZMIN =', IZMIN(jw)                                                                            ! SW 1/23/06
 
!******      Variable initialization
 
             KTWB(jw) = KTWB(jw) + 1
             kt = KTWB(jw)
             ilayer = 0
                       ! SW 1/23/06  11/7/07
 
!            RECOMPUTE INTERNAL WEIR FOR FLOATING WEIR
             if(weir_calc)then
                          !  SW 3/16/18
                 do jwr = 1, niw
                     if(IWR(jwr)>=US(BS(jw)) .AND. IWR(jwr)<=DS(BE(jw)))then
                         if(EKTWR(jwr)==0.0)then
                             KTWR(jwr) = KTWB(jw)
                         else
                             KTWR(jwr) = INT(EKTWR(jwr))
                         endif
                         if(EKBWR(jwr)<=0.0)then
                             do k = KTWR(jwr), KB(IWR(jwr))
                                 if(DEPTHB(k, IWR(jwr))>ABS(EKBWR(jwr)))then
                                     KBWR(jwr) = k
                                     exit
                                 endif
                             enddo
                         else
                             KBWR(jwr) = INT(EKBWR(jwr))
                         endif
 
                         do k = 2, kmx - 1
                             if((k>=KTWR(jwr) .AND. k<=KBWR(jwr)))             &
                              & INTERNAL_WEIR(k, IWR(jwr)) = .TRUE.
                         enddo
                     endif
                 enddo
             endif
 
 
 
             do jb = BS(jw), BE(jw)
                 if(BR_INACTIVE(jb))cycle
                                    ! SW 6/12/2017
                 iu = CUS(jb)
                 id = DS(jb)
                 if(constituents)do1(kt - 1, iu - 1:id + 1) = 0.0
                 do i = iu - 1, id + 1
                     Z(i) = Z(i) - H(kt - 1, jw)
                     H1(kt - 1, i) = H(kt - 1, jw)
                     H1(kt, i) = H(kt, jw) - Z(i)
                     bi(kt, i) = B(KTI(i), i)
                     bi(kt - 1, i) = B(kt - 1, i)
                     AVH1(kt - 1, i) = (H(kt - 1, jw) + H(kt, jw))*0.5
                     AVH1(kt, i) = (H1(kt, i) + H1(kt + 1, i))*0.5
                     if(TRAPEZOIDAL(jw))then
                         call GRID_AREA1(EL(kt, i) - Z(i), EL(kt + 1, i),      &
                           & BH1(kt, i), dummy)                                                                                  !SW 08/03/04
                         BH1(kt - 1, i) = 0.25*H1(kt - 1, jw)                  &
                           & *(BB(kt - 2, i) + 2.*B(kt - 1, i) + BB(kt - 1, i))
                     elseif(KB(i)>=kt)then                   ! SW 1/23/06
                         BH1(kt, i) = BH1(kt - 1, i) + BH1(kt, i)
                                                             ! SW 1/23/06
                         BH1(kt - 1, i) = B(kt - 1, i)*H(kt - 1, jw)
                     else
                    ! SW 1/23/06
                         BH1(kt, i) = BH1(kt - 1, i)
                                                 ! SW 1/23/06
                     endif
                     VOL(kt, i) = BH1(kt, i)*DLX(i)
                     VOL(kt - 1, i) = BH1(kt - 1, i)*DLX(i)
                     BKT(i) = BH1(kt, i)/H1(kt, i)
                     if(KB(i)>=kt)then                      ! SW 1/23/06
                         U(kt, i) = (U(kt - 1, i)*BHR1(kt - 1, i) + U(kt, i)   &
                                  & *BHR(kt, i))/(BHR1(kt - 1, i) + BHR(kt, i))
                         T1(kt, i) = (T1(kt - 1, i)*(BH1(kt, i) - BH(kt, i))   &
                                   & + T1(kt, i)*BH(kt, i))/BH1(kt, i)
                     else
!                        EBRI(JB) = EBRI(JB)-T1(KT,I)*VOL(KT,I)   1/23/06
                         U(kt, i) = U(kt - 1, i)
                                     ! SW 1/23/06
                         T1(kt, i) = T1(kt - 1, i)
                                     ! SW 1/23/06
                     endif
                     if(KB(i)>=kt)then
                                      ! SW 1/23/06
                         C1(kt, i, cn(1:nac))                                  &
                           & = (C1(kt - 1, i, cn(1:nac))*(BH1(kt, i)           &
                           & - BH(kt, i)) + C1(kt, i, cn(1:nac))*BH(kt, i))    &
                           & /BH1(kt, i)
                         CSSK(kt, i, cn(1:nac))                                &
                           & = (CSSK(kt - 1, i, cn(1:nac))*(BH1(kt, i)         &
                           & - BH(kt, i)) + CSSK(kt, i, cn(1:nac))*BH(kt, i))  &
                           & /BH1(kt, i)
                     else
                         C1(kt, i, cn(1:nac)) = C1(kt - 1, i, cn(1:nac))
                                                                     ! SW 1/23/06
                         CSSK(kt, i, cn(1:nac)) = CSSK(kt - 1, i, cn(1:nac))
                                                                      ! SW 1/23/06
                     endif
                    ! SW 1/23/06
                     CSSB(kt, i, cn(1:nac)) = CSSB(kt - 1, i, cn(1:nac))       &
                       & + CSSB(kt, i, cn(1:nac))
                     KF(kt, i, kfcn(1:NAF(jw), jw))                            &
                       & = KF(kt - 1, i, kfcn(1:NAF(jw), jw))
                     KFS(kt, i, kfcn(1:NAF(jw), jw))                           &
                       & = KFS(kt - 1, i, kfcn(1:NAF(jw), jw))
                     C1(kt - 1, i, cn(1:nac)) = 0.0
                     C2(kt - 1, i, cn(1:nac)) = 0.0
                     CSSB(kt - 1, i, cn(1:nac)) = 0.0
                     CSSK(kt - 1, i, cn(1:nac)) = 0.0
                     KF(kt - 1, i, kfcn(1:NAF(jw), jw)) = 0.0
                     KFS(kt - 1, i, kfcn(1:NAF(jw), jw)) = 0.0
                     do je = 1, nep
                         if(kt<=KBI(i))then
                                      ! CB 4/28/06
                             EPM(kt, i, je) = EPM(kt - 1, i, je)               &
                               & + EPM(kt, i, je)
                             EPD(kt, i, je) = EPM(kt, i, je)                   &
                               & /((bi(kt, i) - bi(kt + 1, i) + 2.0*H1(kt, i)) &
                               & *DLX(i))
                             EPC(kt, i, je) = EPM(kt, i, je)/VOL(kt, i)
                             EPM(kt - 1, i, je) = 0.0
                             EPD(kt - 1, i, je) = 0.0
                             EPC(kt - 1, i, je) = 0.0
                         else        ! SW 5/15/06
                             EPM(kt, i, je) = EPM(kt - 1, i, je)
                             EPD(kt, i, je) = EPD(kt - 1, i, je)
                             EPC(kt, i, je) = EPC(kt - 1, i, je)
                             EPM(kt - 1, i, je) = 0.0
                             EPD(kt - 1, i, je) = 0.0
                             EPC(kt - 1, i, je) = 0.0
                         endif        ! CB 4/28/06
                     enddo
                 enddo
                 do i = iu, id
                     do m = 1, nmc
                         if(MACROPHYTE_CALC(jw, m))then
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
                                 coldep = EL(kt, i) - Z(i)*COSA(jb) - colb
                                                       ! cb 3/7/16
                                 if(j<kt)then
                                     MACRM(j, kt, i, m)                        &
                                       & = MACRM(j, kt - 1, i, m)
                                 else
                                     MACRM(j, kt, i, m)                        &
                                       & = MACRM(j, kt - 1, i, m)              &
                                       & + MACRM(j, kt, i, m)
                                 endif
                                 if(CW(j, i)>0.0)then
                                     MACRC(j, kt, i, m) = MACRM(j, kt, i, m)   &
                                       & /(CW(j, i)*coldep*DLX(i))
                                 else
                                     MACRC(j, kt, i, m) = 0.0
                                 endif
                                 MACRM(j, kt - 1, i, m) = 0.0
                                 MACRC(j, kt - 1, i, m) = 0.0
 
                             enddo
                         endif
                     enddo
                     jt = KTI(i)
                     je = KB(i)
                     do j = jt, je
                         MACT(j, kt, i) = 0.0
                         MACT(j, kt - 1, i) = 0.0
                     enddo
                     do m = 1, nmc
                         if(MACROPHYTE_CALC(jw, m))then
                             do j = jt, je
                                 MACT(j, kt, i) = MACRC(j, kt, i, m)           &
                                   & + MACT(j, kt, i)
                             enddo
                         endif
                     enddo
                     do m = 1, nmc
                         tmac = 0.0
                         if(MACROPHYTE_CALC(jw, m))then
                             jt = KTI(i)
                             je = KB(i)
                             do j = jt, je
                                 tmac = tmac + MACRM(j, kt, i, m)
                             enddo
                         endif
                         MAC(kt, i, m) = tmac/(BH1(kt, i)*DLX(i))
                     enddo
                     do m = 1, nmc
                         if(MACROPHYTE_CALC(jw, m))MAC(kt - 1, i, m) = 0.0
                     enddo
                 enddo
 
                 do i = iu - 1, id
                     AVHR(kt - 1, i) = H1(kt - 1, i)                           &
                                     & + (H1(kt - 1, i + 1) - H1(kt - 1, i))   &
                                     & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                         !SW 07/29/04
                     AVHR(kt, i) = H1(kt, i) + (H1(kt, i + 1) - H1(kt, i))     &
                                 & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                             !SW 07/29/04
                     BHR1(kt, i) = BH1(kt, i) + (BH1(kt, i + 1) - BH1(kt, i))  &
                                 & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                             !SW 07/29/04
                     BHR1(kt - 1, i) = BH1(kt - 1, i)                          &
                                     & + (BH1(kt - 1, i + 1) - BH1(kt - 1, i)) &
                                     & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                         !SW 07/29/04
                 enddo
                 U(kt - 1, iu - 1:id + 1) = 0.0
                 W(kt - 1, iu - 1:id + 1) = 0.0
                 sw(kt - 1, iu - 1:id + 1) = 0.0                                                                       !TC 3/9/05
                 p(kt - 1, iu - 1:id + 1) = 0.0
                 AZ(kt - 1, iu - 1:id + 1) = 0.0
                 TKE(kt - 1, iu - 1:id + 1, 1) = 0.0
                                        !sg 10/4/07
                 TKE(kt - 1, iu - 1:id + 1, 2) = 0.0
                                        !sg 10/4/07
                 dz(kt - 1, iu - 1:id + 1) = 0.0
                 admz(kt - 1, iu - 1:id + 1) = 0.0
                 adz(kt - 1, iu - 1:id + 1) = 0.0
                 decay(kt - 1, iu - 1:id + 1) = 0.0
                 if(UP_HEAD(jb))then
                     QUH1(kt, jb) = QUH1(kt, jb) + QUH1(kt - 1, jb)
                     TSSUH1(kt, jb) = TSSUH1(kt - 1, jb) + TSSUH1(kt, jb)
                     CSSUH1(kt, cn(1:nac), jb) = CSSUH1(kt - 1, cn(1:nac), jb) &
                       & + CSSUH1(kt, cn(1:nac), jb)
                 endif
                 if(DN_HEAD(jb))then
                     QDH1(kt, jb) = QDH1(kt, jb) + QDH1(kt - 1, jb)
                     TSSDH1(kt, jb) = TSSDH1(kt - 1, jb) + TSSDH1(kt, jb)
                     CSSDH1(kt, cn(1:nac), jb) = CSSDH1(kt - 1, cn(1:nac), jb) &
                       & + CSSDH1(kt, cn(1:nac), jb)
                 endif
 
!********        Upstream active segment
 
                 iut = US(jb)
                 if(SLOPE(jb)/=0.0)then
                     do i = US(jb) - 1, DS(jb) + 1
                         if(KB(i)<kt)then                                                                            ! SR 10/17/05
                             KB(i) = kt
                             BNEW(KB(i), i) = 0.000001
                                                   ! sw 1/23/06
                             if(DXI(jw)>=0.0)then
                                 DX(KB(i), i) = DXI(jw)
                             else
                                 DX(KB(i), i) = ABS(U(KB(i), i))*ABS(DXI(jw))  &
                                   & *H(k, jw)                             ! SW 8/2/2017
                             endif
                             ilayer(i) = 1
                             T1(KB(i), i) = T1(kt, i)              !    SW 5/15/06    T1(KB(I)-1,I)
                             C1(KB(i), i, cn(1:nac)) = C1(kt, i, cn(1:nac))
                                                                   !    SW 5/15/06    C1(KB(I)-1,I,CN(1:NAC))
                             write(wrn, '(2(A,I8),A,F0.3,A,F0.3)')             &
                                  &'Lowering bottom segment ', i,              &
                                  &' at iteration ', nit, ' at Julian day ',   &
                                 & jday, ' Z(I)=', Z(i)
                             warning_open = .TRUE.
                         endif
                     enddo           ! SW 1/23/06
                     do i = US(jb) - 1, DS(jb) + 1
                                     ! SW 1/23/06
!                        IF (I /= DS(JB)+1) KBMIN(I)   = MIN(KB(I),KB(I+1))  !
!                        SW 1/23/06
                         if(i/=US(jb) - 1)KBMIN(i - 1) = MIN(KB(i - 1), KB(i))
                                                                    ! SW 1/23/06
                         if(KBI(i)<KB(i))then
                             BKT(i) = BH1(kt, i)                               &
                                    & /(H1(kt, i) - (EL(KBI(i) + 1, i)         &
                                    & - EL(KB(i) + 1, i)))                    ! SW 1/23/06
                             DEPTHB(KTWB(jw), i)                               &
                               & = (H1(KTWB(jw), i) - (EL(KBI(i) + 1, i)       &
                               & - EL(KB(i) + 1, i)))                                 ! SW 1/23/06
                             DEPTHM(KTWB(jw), i)                               &
                               & = (H1(KTWB(jw), i) - (EL(KBI(i) + 1, i)       &
                               & - EL(KB(i) + 1, i)))*0.5                                 ! SW 1/23/06
                             AVHR(kt, i)                                       &
                               & = (H1(kt, i) - (EL(KBI(i) + 1, i) - EL(KB(i)  &
                               & + 1, i)))                                     &
                               & + (H1(kt, i + 1) - (EL(KBI(i) + 1, i + 1)     &
                               & - EL(KB(i) + 1, i + 1)) - H1(kt, i)           &
                               & + (EL(KBI(i) + 1, i) - EL(KB(i) + 1, i)))     &
                               & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                                                                                                          ! SW 1/23/06
                         endif
                     enddo
 
                     do i = US(jb), DS(jb)
                             ! SW 1/23/06   11/13/07   US(JB)-1,DS(JB)+1
 
                         if(ilayer(i)==1 .AND. ilayer(i + 1)==0)then
                                                      ! SW 1/23/06
                             bhrsum = 0.0
                             Q(i) = 0.0
                             do k = kt, KBMIN(i)
                                 if(.NOT.INTERNAL_WEIR(k, i))then
                                     bhrsum = bhrsum + BHR1(k, i)
                                     Q(i) = Q(i) + U(k, i)*BHR1(k, i)
                                 endif
                             enddo
                             do k = kt, KBMIN(i)
                                 if(INTERNAL_WEIR(k, i))then
                                     U(k, i) = 0.0
                                 else
                                     U(k, i) = U(k, i) + (QC(i) - Q(i))/bhrsum
                                 endif
                             enddo
                         elseif(ilayer(i)==1 .AND. ilayer(i - 1)==0)then
                             bhrsum = 0.0
                             Q(i - 1) = 0.0
                             do k = kt, KBMIN(i - 1)
                                 if(.NOT.INTERNAL_WEIR(k, i - 1))then
                                     bhrsum = bhrsum + BHR1(k, i - 1)
                                     Q(i - 1) = Q(i - 1) + U(k, i - 1)         &
                                       & *BHR1(k, i - 1)
                                 endif
                             enddo
                             do k = kt, KBMIN(i - 1)
                                 if(INTERNAL_WEIR(k, i - 1))then
                                     U(k, i - 1) = 0.0
                                 else
                                     U(k, i - 1) = U(k, i - 1)                 &
                                       & + (QC(i - 1) - Q(i - 1))/bhrsum
                                 endif
                             enddo
                         endif
                   ! SW 1/23/06
 
                     enddo
                 endif
                 do i = US(jb), DS(jb)
                     if(KB(i) - kt<NL(jb) - 1)iut = i + 1
                     ONE_LAYER(i) = KTWB(jw)==KB(i)
                 enddo
                 if(iut>DS(jb))then
                     if(jb==1)then
                         write(w2err, '(A,I0/A,F0.2,2(A,I0))')                 &
                              &'Fatal error - insufficient segments in branch '&
                             & , jb, 'Julian day = ', jday, ' at iteration ',  &
                             & nit, ' with water surface layer = ', kt                                                           !    ! SEE NEW CODE BELOW
                         write(w2err, '(2(A,I0))')                             &
                              &'Minimum water surface located at segment ',    &
                             & IZMIN(jw), ' with bottom layer at ',            &
                             & KB(IZMIN(jw))
                         text = 'Runtime error - see w2.err'
                         error_open = .TRUE.
                         return
                     else
                         BR_INACTIVE(jb) = .TRUE.
                                    ! SW 6/12/2017
                         if(SNAPSHOT(jw))write(SNP(jw),                        &
                           &'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))'&
                          & )'   Branch Inactive: ', jb, ' at Julian day = ',  &
                           & jday, '   NIT = ', nit
                         write(wrn,                                            &
                           &'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))'&
                          & )'   Branch Inactive: ', jb, ' at Julian day = ',  &
                           & jday, '   NIT = ', nit                                                                                               ! SW 11/16/2018
                         cycle
                     endif
 
       ! Go to 230
                 endif
 
!********        Segment subtraction
 
                 if(iut/=iu)then
                     if(SNAPSHOT(jw))write(SNP(jw), '(/17X,A,I0,A,I0)')        &
                       &' Subtract segments ', iu, ' through ', iut - 1
                     do i = iu, iut - 1
                         do k = kt, KB(i)
                             EBRI(jb) = EBRI(jb) - T1(k, i)*VOL(k, i)
                             CMBRT(cn(1:nac), jb) = CMBRT(cn(1:nac), jb)       &
                               & - C1(k, i, cn(1:nac))*VOL(k, i)               &
                               & + (CSSB(k, i, cn(1:nac))                      &
                               & + CSSK(k, i, cn(1:nac))*VOL(k, i))*dlt
                         enddo
                     enddo
 
                     do i = iu, iut - 1
                         do m = 1, nmc
                             if(MACROPHYTE_CALC(jw, m))then
                                 jt = KTI(i)
                                 je = KB(i)
                                 do j = jt, je
                                     if(j<kt)then
                                         colb = EL(j + 1, i)
                                     else
                                         colb = EL(kt + 1, i)
                                     endif
                    !COLDEP=ELWS(I)-COLB
                                     coldep = EL(kt, i) - Z(i)*COSA(jb) - colb
                                                         ! cb 3/7/16
!                                    MACMBRT(JB,M) =
! MACMBRT(JB,M)-MACRM(J,KT,I,M)+(MACSS(J,KT,I,M)*COLDEP*CW(J,I)*DLX(I))*DLT
                                     MACMBRT(jb, m) = MACMBRT(jb, m)           &
                                       & - MACRM(j, kt, i, m)
                                 enddo
                                 do k = kt + 1, KB(i)
                                     jt = k
                                     je = KB(i)
                                     do j = jt, je
!                                        MACMBRT(JB,M) =
! MACMBRT(JB,M)-MACRM(J,K,I,M)+(MACSS(J,K,I,M)*H2(K,I)*CW(J,I)*DMX(I))*DLT
                                         MACMBRT(jb, m) = MACMBRT(jb, m)       &
                                           & - MACRM(j, k, i, m)
                                     enddo
                                 enddo
                             endif
                         enddo
                     enddo
 
                     f(iu - 1:iut - 1) = 0.0
                     Z(iu - 1:iut - 1) = 0.0
            !ICETH(IU-1:IUT-1) =  0.0    ! SW 9/29/15
                     bhrho(iu - 1:iut - 1) = 0.0
            !ICE(IU-1:IUT-1)   = .FALSE.  ! SW 9/29/15
                     do k = kt, KB(iut)
                         ADX(k, iu - 1:iut - 1) = 0.0
                         DX(k, iu - 1:iut - 1) = 0.0
                         AZ(k, iu - 1:iut - 1) = 0.0
                         TKE(k, iu - 1:iut - 1, 1) = 0.0
                                                 !SG  10/4/07
                         TKE(k, iu - 1:iut - 1, 2) = 0.0
                                                 !SG  10/4/07
                         SAZ(k, iu - 1:iut - 1) = 0.0
                         U(k, iu - 1:iut - 1) = 0.0
                         SU(k, iu - 1:iut - 1) = 0.0
                         T1(k, iu - 1:iut - 1) = 0.0
                         tss(k, iu - 1:iut - 1) = 0.0
                         QSS(k, iu - 1:iut - 1) = 0.0
                         C1(k, iu - 1:iut - 1, cn(1:nac)) = 0.0
                         C2(k, iu - 1:iut - 1, cn(1:nac)) = 0.0
                         c1s(k, iu - 1:iut - 1, cn(1:nac)) = 0.0
                         CSSB(k, iu - 1:iut - 1, cn(1:nac)) = 0.0
                         CSSK(k, iu - 1:iut - 1, cn(1:nac)) = 0.0
                     enddo
 
                     do m = 1, nmc
                         if(MACROPHYTE_CALC(jw, m))then
                             MAC(k, i, m) = 0.0
                             MACT(j, k, i) = 0.0
                         endif
                     enddo
                     jt = KTI(i)
                     je = KB(i)
                     do j = jt, je
                         do m = 1, nmc
                             if(MACROPHYTE_CALC(jw, m))MACRC(j, k, i, m) = 0.0
                         enddo
                     enddo
 
                     iu = iut
                     CUS(jb) = iu
                     Z(iu - 1) = (EL(kt, iu - 1) - (EL(kt, iu) - Z(iu)*COSA(jb)&
                               & ))/COSA(jb)
                     SZ(iu - 1) = Z(iu)
                     KTI(iu - 1) = KTI(iu)
                     if(.NOT.TRAPEZOIDAL(jw))then
                         bi(kt, iu - 1) = B(KTI(iu - 1), i)
                         H1(kt, iu - 1) = H(kt, jw) - Z(iu - 1)
                         BH1(kt, iu - 1) = BNEW(KTI(iu - 1), iu - 1)           &
                           & *(EL(kt, iu - 1) - EL(KTI(iu - 1) + 1, iu - 1)    &
                           & - Z(iu - 1)*COSA(jb))/COSA(jb)                                                      ! sw 1/23/06  Bnew(KTI(IU-1),IU-1)*(EL(KT,IU-1)-EL(KTI(IU-1)+1,IU-1)-Z(IU-1)*COSA(JB))/COSA(JB)     ! SR 10/17/05
                         if(kt>=KB(iu - 1))BH1(kt, iu - 1) = BNEW(kt, iu - 1)  &
                          & *H1(kt, iu - 1)                                  ! sw 1/23/06
                         do k = KTI(iu - 1) + 1, kt
                             BH1(kt, iu - 1) = BH1(kt, iu - 1) + BH1(k, iu - 1)
                         enddo
                     else
                         call GRID_AREA1(EL(kt, i) - Z(i), EL(kt + 1, iu - 1), &
                           & BH1(kt, iu - 1), bi(kt, iu - 1))                                                                    !SW 08/03/04
                         BH1(kt, i) = 0.25*H(kt, jw)                           &
                                    & *(BB(kt - 1, i) + 2.*B(kt, i) + BB(kt, i)&
                                    & )
                         H1(kt, i) = H(kt, jw) - Z(i)
                     endif
                     BKT(iu - 1) = BH1(kt, iu - 1)/H1(kt, iu - 1)
                     BHR1(kt, iu - 1) = BH1(kt, iu - 1)                        &
                                      & + (BH1(kt, iu) - BH1(kt, iu - 1))      &
                                      & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                        !SW 07/29/04
                     if(UH_EXTERNAL(jb))KB(iu - 1) = KB(iu)
                     if(UH_INTERNAL(jb))then
                         if(JBUH(jb)>=BS(jw) .AND. JBUH(jb)<=BE(jw))then
                             KB(iu - 1) = MIN(KB(UHS(jb)), KB(iu))
                         else
                             do kkb = kt, kmx
                                 if(EL(kkb, iu)<=EL(KB(UHS(jb)), UHS(jb)))exit
                             enddo
                             KB(iu - 1) = MIN(kkb, KB(iu))
                         endif
                     endif
                 endif
                 if(constituents)then ! SW 5/15/06
                     call TEMPERATURE_RATES
                     call KINETIC_RATES
                 endif
 
 
!********        Total active cells
 
                 do i = iu, id
                     ntac = ntac - 1
                 enddo
                 ntacmn = MIN(ntac, ntacmn)
             enddo
             call INTERPOLATION_MULTIPLIERS
 
!******      Additional layer subtractions
 
             ZMIN(jw) = -1000.0
             do jb = BS(jw), BE(jw)
                 if(BR_INACTIVE(jb))cycle
                                    ! SW 6/12/2017
                 do i = CUS(jb), DS(jb)
                     ZMIN(jw) = MAX(ZMIN(jw), Z(i))
                 enddo
             enddo
             sub_layer = ZMIN(jw)>0.60*H(kt, jw) .AND. kt<ktmax                                                       ! SR 10/17/05
         enddo
     enddo
 
!**  Temporary downstream head segment
 
     do jb = 1, nbr
         if(BR_INACTIVE(jb))cycle
                                ! SW 6/12/2017
         if(DHS(jb)>0)then
             do jjb = 1, nbr
                 if(DHS(jb)>=US(jjb) .AND. DHS(jb)<=DS(jjb))exit
             enddo
             if(CUS(jjb)>DHS(jb))CDHS(jb) = CUS(jjb)
         endif
     enddo
     end subroutine LAYERADDSUB
