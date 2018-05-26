!*==initcond.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
!***********************************************************************************************************************************
!**                                              Task 1.4.4: Initial conditions                                                   **
!***********************************************************************************************************************************
     subroutine INITCOND
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
     character(8) :: iblank
     character(1) :: ichar
     real :: tmac, xsar
     external RESTART_OUTPUT
!
!*** End of declarations rewritten by SPAG
!
 
 
     bic = b
 
     do jw = 1, nwb
         kt = KTWB(jw)
         if(VERT_PROFILE(jw))then
 
!****        Temperature and water quality
 
             open(VPR(jw), file = VPRFN(jw), status = 'OLD')
             read(VPR(jw), '(A1)')ichar
 
             if(ichar=='$')then
                 read(VPR(jw), '(/)')
                 if(VERT_TEMP(jw))read(VPR(jw), *)iblank,                      &
                                     & (tvp(k, jw), k = kt, KBMAX(jw))
                 if(constituents)then
                     do jc = 1, nct
                         if(VERT_CONC(jc, jw))read(VPR(jw), *)iblank,          &
                          & (cvp(k, jc, jw), k = kt, KBMAX(jw))
                     enddo
                     do je = 1, nep
                         if(VERT_EPIPHYTON(jw, je))read(VPR(jw), *)iblank,     &
                          & (epivp(k, jw, je), k = kt, KBMAX(jw))
                     enddo
                     do m = 1, nmc                               ! cb 8/21/15
                         if(VERT_MACROPHYTE(jw, m))read(VPR(jw), *)iblank,     &
                          & (macrcvp(k, jw, m), k = kt, KBMAX(jw))
                     enddo
                     if(VERT_SEDIMENT(jw))read(VPR(jw), *)iblank,              &
                      & (sedvp(k, jw), k = kt, KBMAX(jw))
                 endif
             else
                 if(VERT_TEMP(jw))read(VPR(jw), '(//(8X,9F8.0))')              &
                                     & (tvp(k, jw), k = kt, KBMAX(jw))
                 if(constituents)then
                     do jc = 1, nct
                         if(VERT_CONC(jc, jw))read(VPR(jw), '(//(8X,9F8.0))')  &
                          & (cvp(k, jc, jw), k = kt, KBMAX(jw))
                     enddo
                     do je = 1, nep
                         if(VERT_EPIPHYTON(jw, je))                            &
                           &read(VPR(jw), '(//(8X,9F8.0))')                    &
                          & (epivp(k, jw, je), k = kt, KBMAX(jw))
                     enddo
                     do m = 1, nmc                                 ! cb 8/21/15
                         if(VERT_MACROPHYTE(jw, m))                            &
                           &read(VPR(jw), '(//(8X,9F8.0))')                    &
                          & (macrcvp(k, jw, m), k = kt, KBMAX(jw))
                     enddo
                     if(VERT_SEDIMENT(jw))read(VPR(jw), '(//(8X,9F8.0))')      &
                      & (sedvp(k, jw), k = kt, KBMAX(jw))
                 endif
             endif
         endif
 
!**      Longitudinal/vertical initial profiles
 
         if(LONG_PROFILE(jw))then
             open(LPR(jw), file = LPRFN(jw), status = 'OLD')
             read(LPR(jw), '(A1)')ichar
             if(ichar=='$')read(LPR(jw), *)
         endif
 
!**      Branch related variables
 
         if(.NOT.restart_in)then
             if(LONG_TEMP(jw) .AND. ichar=='$')read(LPR(jw), *)
             do jb = BS(jw), BE(jw)
 
!******          Temperature
 
                 do i = CUS(jb), DS(jb)
                     if(LONG_TEMP(jw))then
                         if(ichar=='$')then
                             read(LPR(jw), *)iblank, (t1(k, i), k = kt, KB(i))
                         else
                             read(LPR(jw), '(//(8X,9F8.0))')                   &
                                & (t1(k, i), k = kt, KB(i))
                         endif
                     endif
                     do k = kt, KB(i)
                         if(ISO_TEMP(jw))t1(k, i) = T2I(jw)
                         if(VERT_TEMP(jw))t1(k, i) = tvp(k, jw)
                         T2(k, i) = t1(k, i)
                     enddo
                 enddo
             enddo
 
!****        Constituents
 
             do jc = 1, nac
                 if(LONG_CONC(CN(jc), jw) .AND. ichar=='$')read(LPR(jw), *)
                 do jb = BS(jw), BE(jw)
                     do i = CUS(jb), DS(jb)
                         jac = CN(jc)
                         if(LONG_CONC(jac, jw))then
                             if(ichar=='$')then
                                 read(LPR(jw), *)iblank,                       &
                                    & (c2(k, i, jac), k = kt, KB(i))
                             else
                                 read(LPR(jw), '(//(8X,9F8.0))')               &
                                    & (c2(k, i, jac), k = kt, KB(i))
                             endif
                         endif
                         do k = kt, KB(i)
                             if(ISO_CONC(jac, jw))c2(k, i, jac) = C2I(jac, jw)
                             if(VERT_CONC(jac, jw))c2(k, i, jac)               &
                              & = cvp(k, jac, jw)
                             C1(k, i, jac) = c2(k, i, jac)
                             C1S(k, i, jac) = C1(k, i, jac)
                         enddo
                     enddo
                 enddo
             enddo
 
!****        Epiphyton
 
 
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))then
                     if(LONG_EPIPHYTON(jw, je) .AND. ichar=='$')               &
                      & read(LPR(jw), *)
                     do jb = BS(jw), BE(jw)
                         do i = CUS(jb), DS(jb)
                             if(LONG_EPIPHYTON(jw, je))then
                                 if(ichar=='$')then
                                     read(LPR(jw), *)iblank,                   &
                                       & (epd(k, i, je), k = kt, KB(i))
                                 else
                                     read(LPR(jw), '(//(8X,9F8.0))')           &
                                       & (epd(k, i, je), k = kt, KB(i))
                                 endif
                             endif
                             if(ISO_EPIPHYTON(jw, je))epd(:, i, je)            &
                              & = EPICI(jw, je)
                             if(VERT_EPIPHYTON(jw, je))epd(:, i, je)           &
                              & = epivp(:, jw, je)                                                       ! CB 5/16/2009
                         enddo
                     enddo
                 endif
             enddo
 
!****        macrophytes - added 8/21/15
 
 
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     if(LONG_MACROPHYTE(jw, m) .AND. ichar=='$')               &
                      & read(LPR(jw), *)
                     do jb = BS(jw), BE(jw)
                         do i = CUS(jb), DS(jb)
                             if(LONG_MACROPHYTE(jw, m))then
                                 if(ichar=='$')then
                                     read(LPR(jw), *)iblank,                   &
                                       & (macrclp(k, i, m), k = kt, KB(i))
                                 else
                                     read(LPR(jw), '(//(8X,9F8.0))')           &
                                       & (macrclp(k, i, m), k = kt, KB(i))
                                 endif
                             endif
              !IF (ISO_macrophyte(JW,m))  macrc(:,I,m) = macwbci(JW,m)
              !IF (VERT_macrophyte(JW,m)) macrc(:,I,m) = macrcvp(:,JW,m)
                         enddo
                     enddo
                 endif
             enddo
 
!****        Sediments
 
             do jb = BS(jw), BE(jw)
 !     SDKV(:,US(JB):DS(JB))=SDK(JW)
                 sdkv(:, US(jb) - 1:DS(jb) + 1) = SDK(jw)
                                             ! SW 9/28/13
                 if(SEDIMENT_CALC(jw))then
                     if(LONG_SEDIMENT(jw) .AND. jb==BS(jw))read(LPR(jw), *)
                     do i = CUS(jb), DS(jb)
                         if(LONG_SEDIMENT(jw))then
                             if(ichar=='$')then
                                 read(LPR(jw), *)iblank,                       &
                                    & (sed(k, i), k = kt, KB(i))
                             else
                                 read(LPR(jw), '(//(8X,9F8.0))')               &
                                    & (sed(k, i), k = kt, KB(i))
                             endif
                         endif
                         do k = kt, KB(i)
                             if(ISO_SEDIMENT(jw))sed(k, i) = SEDCI(jw)
                             if(VERT_SEDIMENT(jw))sed(k, i) = sedvp(k, jw)
                         enddo
                         sed(kt, i) = sed(kt, i)/H2(kt, i)
                         sed(kt + 1:KB(i), i) = sed(kt + 1:KB(i), i)           &
                           & /H2(kt + 1:KB(i), i)
                     enddo
                 endif
             enddo
             do jb = BS(jw), BE(jw)
                 if(SEDIMENT_CALC(jw))then
                     do i = CUS(jb), DS(jb)
                         do k = kt, KB(i)
                             if(ISO_SEDIMENT(jw))SEDP(k, i) = ORGP(jw)         &
                              & *SEDCI(jw)
                             if(VERT_SEDIMENT(jw))SEDP(k, i) = sedvp(k, jw)    &
                              & *ORGP(jw)
                             if(LONG_SEDIMENT(jw))SEDP(k, i) = ORGP(jw)        &
                              & *sed(k, i)
                         enddo
                         SEDP(kt, i) = SEDP(kt, i)/H2(kt, i)
                         SEDP(kt + 1:KB(i), i) = SEDP(kt + 1:KB(i), i)         &
                           & /H2(kt + 1:KB(i), i)
                     enddo
                 endif
             enddo
             do jb = BS(jw), BE(jw)
                 if(SEDIMENT_CALC(jw))then
                     do i = CUS(jb), DS(jb)
                         do k = kt, KB(i)
                             if(ISO_SEDIMENT(jw))SEDN(k, i) = ORGN(jw)         &
                              & *SEDCI(jw)
                             if(VERT_SEDIMENT(jw))SEDN(k, i) = sedvp(k, jw)    &
                              & *ORGN(jw)
                             if(LONG_SEDIMENT(jw))SEDN(k, i) = ORGN(jw)        &
                              & *sed(k, i)
                         enddo
                         SEDN(kt, i) = SEDN(kt, i)/H2(kt, i)
                         SEDN(kt + 1:KB(i), i) = SEDN(kt + 1:KB(i), i)         &
                           & /H2(kt + 1:KB(i), i)
                     enddo
                 endif
             enddo
             do jb = BS(jw), BE(jw)
                 if(SEDIMENT_CALC(jw))then
                     do i = CUS(jb), DS(jb)
                         do k = kt, KB(i)
                             if(ISO_SEDIMENT(jw))SEDC(k, i) = SEDCI(jw)        &
                              & *ORGC(jw)
                             if(VERT_SEDIMENT(jw))SEDC(k, i) = sedvp(k, jw)    &
                              & *ORGC(jw)
                             if(LONG_SEDIMENT(jw))SEDC(k, i) = ORGC(jw)        &
                              & *sed(k, i)
                         enddo
                         SEDC(kt, i) = SEDC(kt, i)/H2(kt, i)
                         SEDC(kt + 1:KB(i), i) = SEDC(kt + 1:KB(i), i)         &
                           & /H2(kt + 1:KB(i), i)
                     enddo
                 endif
             enddo
 
             sed(:, US(BS(jw)):DS(BE(jw))) = sed(:, US(BS(jw)):DS(BE(jw)))     &
               & *FSED(jw)
             SEDP(:, US(BS(jw)):DS(BE(jw))) = SEDP(:, US(BS(jw)):DS(BE(jw)))   &
               & *FSED(jw)
             SEDN(:, US(BS(jw)):DS(BE(jw))) = SEDN(:, US(BS(jw)):DS(BE(jw)))   &
               & *FSED(jw)
             SEDC(:, US(BS(jw)):DS(BE(jw))) = SEDC(:, US(BS(jw)):DS(BE(jw)))   &
               & *FSED(jw)
 
!            Amaila start Additional sediment compartments
             do jb = BS(jw), BE(jw)
                 if(SEDIMENT_CALC1(jw))then
            !IF(LONG_SEDIMENT(JW).AND.JB==BS(JW))READ (LPR(JW),*)
                     if(LONG_SEDIMENT(jw) .AND. ichar=='$')read(LPR(jw), *)
                                                                      ! cb 6/10/13
          !DO I=CUS(JB),DS(JB)
                     do i = US(jb), DS(jb) ! cb 6/17/17
                         if(LONG_SEDIMENT1(jw))then
                             if(ichar=='$')then
                    !READ (LPR(JW),*)IBLANK, (SED1(K,I),K=KT,KB(I))
                                 read(LPR(jw), *)iblank,                       &
                                    & (sed1(k, i), k = 2, KB(i))    ! cb 6/17/17
                             else
                    !READ (LPR(JW),'(//(8X,9F8.0))') (SED1(K,I),K=KT,KB(I))
                                 read(LPR(jw), '(//(8X,9F8.0))')               &
                                    & (sed1(k, i), k = 2, KB(i))           ! cb 6/17/17
                             endif
                         endif
                         do k = kt, KB(i)
                             if(ISO_SEDIMENT1(jw))sed1(k, i) = SEDCI1(jw)
                             if(VERT_SEDIMENT1(jw))sed1(k, i) = SEDVP1(k, jw)
                         enddo
            !SED1(KT,I)         = SED1(KT,I)/H2(KT,I)   ! intial conditions for "tree" sediment compartments are given in g/m^3
            !SED1(KT+1:KB(I),I) = SED1(KT+1:KB(I),I)/H2(KT+1:KB(I),I)
                     enddo
                 endif
             enddo
 
             do jb = BS(jw), BE(jw)
                 if(SEDIMENT_CALC2(jw))then
            !IF(LONG_SEDIMENT(JW).AND.JB==BS(JW))READ (LPR(JW),*)
                     if(LONG_SEDIMENT(jw) .AND. ichar=='$')read(LPR(jw), *)
                                                                      ! cb 6/10/13
          !DO I=CUS(JB),DS(JB)
                     do i = US(jb), DS(jb)
                                          ! cb 6/17/17
                         if(LONG_SEDIMENT2(jw))then
                             if(ichar=='$')then
                    !READ (LPR(JW),*)IBLANK, (SED2(K,I),K=KT,KB(I))
                                 read(LPR(jw), *)iblank,                       &
                                    & (sed2(k, i), k = 2, KB(i))    ! cb 6/17/17
                             else
                    !READ (LPR(JW),'(//(8X,9F8.0))') (SED2(K,I),K=KT,KB(I))
                                 read(LPR(jw), '(//(8X,9F8.0))')               &
                                    & (sed2(k, i), k = 2, KB(i))           ! cb 6/17/17
                             endif
                         endif
                         do k = kt, KB(i)
                             if(ISO_SEDIMENT2(jw))sed2(k, i) = SEDCI2(jw)
                             if(VERT_SEDIMENT2(jw))sed2(k, i) = SEDVP2(k, jw)
                         enddo
            !SED2(KT,I)         = SED2(KT,I)/H2(KT,I)      ! intial conditions for "tree" sediment compartments are given in g/m^3
            !SED2(KT+1:KB(I),I) = SED2(KT+1:KB(I),I)/H2(KT+1:KB(I),I)
                     enddo
                 endif
             enddo
 
             do jb = BS(jw), BE(jw)
                              ! 9/3/17
                 do i = US(jb), DS(jb)
                     do k = kt, KB(i)
                         SDFIRSTADD(k, i) = .FALSE.
                     enddo
                 enddo
             enddo
 
             sed1(:, US(BS(jw)):DS(BE(jw))) = sed1(:, US(BS(jw)):DS(BE(jw)))   &
               & *FSEDC1(jw)                                                     ! cb 6/7/17
             sed2(:, US(BS(jw)):DS(BE(jw))) = sed2(:, US(BS(jw)):DS(BE(jw)))   &
               & *FSEDC2(jw)
             sed1ic = sed1
                       ! cb 6/17/17
             sed2ic = sed2
                       ! cb 6/17/17
 
 
!            Amaila end
 
             do jb = BS(jw), BE(jw)
                 do m = 1, nmc
                     if(MACROPHYTE_CALC(jw, m))then
 
!C                       DISTRIBUTING INITIAL MACROPHYTE CONC TO BOTTOM COLUMN
!                        CELLS; MACWBCI = G/M^3
                         do i = CUS(jb), DS(jb)
 
                             depkti = ELWS(i) - EL(KTI(i) + 1, i)
 
                             if(depkti>=thrkti)then
                                 KTICOL(i) = .TRUE.
                                 jt = KTI(i)
                             else
                                 KTICOL(i) = .FALSE.
                                 jt = KTI(i) + 1
                             endif
 
                             je = KB(i)
                             do j = jt, je
                                 if(j<=kt)then
                                     k = kt
                                 else
                                     k = j
                                 endif
                !MACRC(J,K,I,M) = MACWBCI(JW,M)
                !SMACRC(J,K,I,M) = MACWBCI(JW,M)
                                 if(ISO_MACROPHYTE(jw, m))MACRC(j, k, i, m)    &
                                  & = MACWBCI(jw, m)                          ! cb 8/24/15
                                 if(VERT_MACROPHYTE(jw, m))MACRC(j, k, i, m)   &
                                  & = macrcvp(k, jw, m)
                                 if(LONG_MACROPHYTE(jw, m))MACRC(j, k, i, m)   &
                                  & = macrclp(k, i, m)
                                 SMACRC(j, k, i, m) = MACRC(j, k, i, m)
                             enddo
                         enddo
 
                         do i = CUS(jb), DS(jb)
                             tmac = 0.0
                             xsar = 0.0
                             do k = KTI(i), kt
                                 jt = k
                                 je = KB(i)
                                 colb = EL(k + 1, i)
                                 coldep = ELWS(i) - colb
                                 do j = jt, je
                                     tmac = tmac + MACRC(j, kt, i, m)*CW(j, i) &
                                       & *coldep
                                     xsar = xsar + CW(j, i)*coldep
                                 enddo
                             enddo
                             MAC(kt, i, m) = tmac/xsar
                             SMAC(kt, i, m) = MAC(kt, i, m)
 
                             do k = kt + 1, KB(i)
                                 jt = k
                                 je = KB(i)
                                 tmac = 0.0
                                 do j = jt, je
                                     tmac = tmac + MACRC(j, k, i, m)*CW(j, i)
                                 enddo
                                 MAC(k, i, m) = tmac/b(k, i)
                                 SMAC(k, i, m) = MAC(k, i, m)
                             enddo
                         enddo
 
                         do i = CUS(jb), DS(jb)
                             jt = KTI(i)
                             je = KB(i)
                             do j = jt, je
                                 if(j<kt)then
                                     colb = EL(j + 1, i)
                                 else
                                     colb = EL(kt + 1, i)
                                 endif
                                 coldep = ELWS(i) - colb
                                 MACRM(j, kt, i, m) = MACRC(j, kt, i, m)       &
                                   & *coldep*CW(j, i)*DLX(i)
                                 SMACRM(j, kt, i, m) = MACRM(j, kt, i, m)
                             enddo
 
                             do k = kt + 1, KB(i)
 
                                 jt = k
                                 je = KB(i)
 
                                 do j = jt, je
 
                                     MACRM(j, k, i, m) = MACRC(j, k, i, m)     &
                                       & *H2(k, i)*CW(j, i)*DLX(i)
                                     SMACRM(j, k, i, m) = MACRM(j, k, i, m)
                                 enddo
 
                             enddo
                         enddo
 
                     endif
                 enddo
             enddo
!            V3.5 END
 
!****        ENERGY
 
             do jb = BS(jw), BE(jw)
                 do i = CUS(jb), DS(jb)
                     if(ENERGY_BALANCE(jw))then
                         do k = kt, KB(i)
                             EBRI(jb) = EBRI(jb) + T2(k, i)*DLX(i)*BH2(k, i)
                         enddo
                     endif
                     do k = kt, KB(i)
                         CMBRT(CN(1:nac), jb) = CMBRT(CN(1:nac), jb)           &
                           & + c2(k, i, CN(1:nac))*DLX(i)*BH2(k, i)
                     enddo
                 enddo
 
!                V3.5 START
!C               INITIALIZING MACROPHYTE TEMPORAL MASS BALANCE TERM....
                 do m = 1, nmc
                     if(MACROPHYTE_CALC(jw, m))then
                         do i = CUS(jb), DS(jb)
                             if(KTICOL(i))then
                                 jt = KTI(i)
                             else
                                 jt = KTI(i) + 1
                             endif
                             je = KB(i)
                             do j = jt, je
                                 MACMBRT(jb, m) = MACMBRT(jb, m)               &
                                   & + MACRM(j, kt, i, m)
                             enddo
                             do k = kt + 1, KB(i)
                                 jt = k
                                 je = KB(i)
                                 do j = jt, je
                                     MACMBRT(jb, m) = MACMBRT(jb, m)           &
                                       & + MACRM(j, k, i, m)
                                 enddo
                             enddo
                         enddo
                     endif
                 enddo
 
!******          Ice cover
 
                 if(ICE_CALC(jw))then
                     iceth(CUS(jb):DS(jb)) = ICETHI(jw)    ! SW 9/29/15 only initialize CUS to DS
                     ice(US(jb):DS(jb)) = iceth(US(jb):DS(jb))>0.0
                     do i = CUS(jb), DS(jb)                ! SW 9/29/15
                         ICEBANK(i) = iceth(i)*BI(kt, i)*DLX(i)
                                                           ! Initial volume of ice in m3
                     enddo
                 endif
 
!******          Vertical eddy viscosity
 
                 iut = CUS(jb)
                 idt = DS(jb) - 1
                 if(UP_HEAD(jb))iut = iu - 1
                 if(DN_HEAD(jb))idt = id
                 do i = iut, idt
                     do k = kt, KB(i) - 1
                         AZ(k, i) = azmin
                         TKE(k, i, 1) = 1.25E-7
                         TKE(k, i, 2) = 1.0E-9
                     enddo
                 enddo
                 do jwr = 1, niw
                     if(weir_calc)AZ(MAX(kt, KTWR(jwr) - 1):KBWR(jwr), IWR(jwr)&
                      & ) = 0.0
                 enddo
             enddo
         endif
 
!**      Horizontal diffusivities
 
         do jb = BS(jw), BE(jw)
             do i = CUS(jb), DS(jb) - 1
                 do k = kt, KBMIN(i)
                     DX(k, i) = ABS(DXI(jw))
                                    ! SW 8/2/2017 FIRST TIME STEP EVEN IF NEGATIVE USE AS ABS OF DX SINCE IT WILL ALWAYS BE LESS THAN 1
                     if(INTERNAL_WEIR(k, i))DX(k, i) = 0.0
                 enddo
             enddo
         enddo
         if(VERT_PROFILE(jw))close(VPR(jw))
         if(LONG_PROFILE(jw))close(LPR(jw))
     enddo
 
!    Atmospheric pressure
     if(constituents)palt(:) = (1.0 - ELWS(:)/1000.0/44.3)**5.25     ! SW 2/3/08
 
 
     end subroutine INITCOND
