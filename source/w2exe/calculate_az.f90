!*==calculate_az.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                        S U B R O U T I N E   C A L C U L A T E  A Z                                           **
!***********************************************************************************************************************************
 
     recursive subroutine calculate_az
     use GEOMC
     use GLOBAL
     use TRANS
     use EDDY
     use KINETIC
     use MACROPHYTEC
     use LOGICC
     use MAIN, ONLY:warning_open
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND), save :: ab, ata, az0, bouk, buoy, convex, convez, convkx,       &
                     & convkz, depth, depthl, depthr, effric, expaz, gc2, hrad,&
                     & prdk, prhe, prhk, ri, riaz0, riaz1, slm, udl, udr, unse,&
                     & unst, ustar, ustarb, ustarbkt, viscf, visck, zd, zdlr
     integer, save :: idt, iut, j, jj, k, kl
     real(R8KIND), dimension(2), save :: sig
     real(R8KIND), save :: tkemin1, tkemin2
     real(R8KIND), external :: WALLFUNCTION
!
!*** End of declarations rewritten by SPAG
!
     data sig/1.0D0, 1.3D0/, tkemin1/1.25D-7/, tkemin2/1.0D-9/
 
     if(AZC(jw)=='     TKE')then
         call CALCULATE_TKE
     elseif(AZC(jw)=='    TKE1')then
         call CALCULATE_TKE1
     else
         do k = kt, KBMIN(i) - 1
             VSH(k, i) = ((U(k + 1, i) - U(k, i))/((AVH2(k,i) + AVH2(k,i+1))*  &
                       & 0.5))**2
         enddo
 
         do k = kt, KBMIN(i) - 1
             call CALCULATE_AZ0
             buoy = (RHO(k + 1, i) - RHO(k, i) + RHO(k + 1, i + 1)             &
                  & - RHO(k, i + 1))/(AVH2(k, i) + AVH2(k, i + 1))                ! 2.0*AVH2(K,I)
             riaz0 = DLOG(az0/AZMAX(jw))*0.666666667D0                            ! /1.5
             ri = g*buoy/(rhow*VSH(k, i) + nonzero)
             riaz1 = DMAX1(ri, riaz0)
             riaz1 = DMIN1(riaz1, 10.0D0)
             expaz = DEXP( - 1.5D0*riaz1)
             AZ(k, i) = DMAX1(azmin, az0*expaz + azmin*(1.0D0 - expaz)) ! AZ computed at lower edge of cell
             DZT(k, i) = DMAX1(dzmin, frazdz*(az0*expaz + dzmin*(1.0 - expaz)))
                                                                       ! DZ computed at lower edge of cell - later averaged to cell center lower edge
         enddo
     endif
     return
 
     entry CALCULATE_AZ0
     if(k==kt)then
         depthl = (ELWS(i) - EL(KB(i), i) + H2(KB(i), i)*COSA(jb))/COSA(jb) !(EL(KT,I)  -Z(I)  *COSA(JB)-EL(KB(I),I)    +H2(KB(I),I)  *COSA(JB))/COSA(JB)
         depthr = (ELWS(i + 1) - EL(KB(i + 1), i + 1) + H2(KB(i + 1), i)       &
                & *COSA(jb))/COSA(jb)                                       !(EL(KT,I+1)-Z(I+1)*COSA(JB)-EL(KB(I+1),I+1)+H2(KB(I+1),I)*COSA(JB))/COSA(JB)
         zdlr = depthl - H2(kt, i) + depthr - H2(kt, i + 1)
         if(AZC(jw)=='     RNG' .OR. AZC(jw)=='   PARAB')then
             depth = (depthr + depthl)*0.5
             ustar = SQRT(g*depth*SLOPEC(jb))
             if(SLOPEC(jb)==0.0)then
                 ustar = 0.0
                 do kl = kt, KBMIN(i)
                     ustar = ustar + SQRT(AZ(kl - 1, i)*SQRT(VSH(kl - 1, i))   &
                           & /RHO(kl, i))
                 enddo
                 ustar = ustar/(KBMIN(i) - kt + 1)
             endif
         endif
     else
         zdlr = (EL(k - 1, i) - EL(KB(i), i) + H2(KB(i), i)*COSA(jb))/COSA(jb) &
              & + (EL(k - 1, i + 1) - EL(KB(i), i + 1) + H2(KB(i + 1), i)      &
              & *COSA(jb))/COSA(jb)
     endif
     zd = zdlr/(depthl + depthr)
     if(AZC(jw)=='    NICK')then
         slm = (depth*(0.14 - 0.08*(1.0 - zd)**2 - 0.06*(1.0 - zd)**4))**2
         az0 = MAX(azmin, slm*SQRT(VSH(k, i)))
     elseif(AZC(jw)=='     RNG')then
         visck = DEXP((T2(k, i) + 495.691)/( - 37.3877))
         if(T2(k, i)>30.0)visck = DEXP((T2(k, i) + 782.190)/( - 57.7600))
         viscf = MAX(0.0, 0.08477*((zdlr*0.5*ustar/visck)**3)                  &
               & *((1.0 - zdlr*0.5/depth)**3) - 100.0)
         viscf = (1.0 + viscf)**0.33333333333333
         az0 = MAX(azmin, visck*viscf)
     elseif(AZC(jw)=='   PARAB')then
         az0 = MAX(azmin, 0.41*ustar*zdlr*0.5*(1.0 - zd))
     else
         slm = hmax2
         if(AZC(jw)=='     W2N')slm = ((depthr + depthl)*0.5*(0.14 - 0.08*(1.0 &
                                    & - zd)**2 - 0.06*(1.0 - zd)**4))**2
         az0 = 0.4*slm*SQRT(VSH(k, i)                                          &
             & + ((FRICBR(k,i) + WSHY(i)*DECAY(k,i))/(AZ(k,i) + nonzero))**2)  &
             & + azmin
     endif
     return
 
     entry CALCULATE_TKE
     ustar = SQRT(1.25*CZ(i)*WIND10(i)**2/RHO(kt, i))
     if(MANNINGS_N(jw))then
         hrad = BH1(kt, i)/(B(KTI(i), i) - B(kt + 1, i) + 2.*AVH1(kt, i))
                                                               ! HRAD = BHR1(KT,I)/(BR(KTI(I),I)-BR(KT+1,I)+2.*AVH1(KT,I))  SW 10/5/07  These calculations are at the segment centers and vertical center of a layer
         if(macrophyte_on .AND. MANNINGS_N(jw))then
             call MACROPHYTE_FRICTION(hrad, FRIC(i), effric, kt, i)
             gc2 = g*effric*effric/hrad**0.33333333
         elseif(.NOT.macrophyte_on .AND. MANNINGS_N(jw))then
             gc2 = g*FRIC(i)*FRIC(i)/hrad**0.33333333
         endif
     else
         gc2 = 0.0
         if(FRIC(i)/=0.0)gc2 = g/(FRIC(i)*FRIC(i))
     endif
     ustarbkt = SQRT(gc2)*ABS(0.5*(U(kt, i) + U(kt, i - 1)))                      ! SG 10/4/07
     TKE(kt, i, 1) = (3.33*(ustar*ustar + ustarbkt*ustarbkt))                  &
                   & *(BH2(kt, i)/BH1(kt, i))                                     ! SG 10/4/07
     TKE(kt, i, 2) = (ustar*ustar*ustar + ustarbkt*ustarbkt*ustarbkt)          &
                   & *5.0/H1(kt, i)*(BH2(kt, i)/BH1(kt, i))                       ! SG 10/4/07
     do k = kt + 1, KB(i) - 1
         bouk = MAX(AZ(k, i)*g*(RHO(k + 1, i) - RHO(k, i))/(H(k, jw)*rhow),    &
              & 0.0)
         prdk = AZ(k, i)                                                       &
              & *(0.5*(U(k, i) + U(k, i - 1) - U(k + 1, i) - U(k + 1, i - 1))  &
              & /(H(k, jw)*0.5 + H(k + 1, jw)*0.5))**2.0                          ! SG 10/4/07
         prhe = 10.0*gc2**1.25*ABS(0.5*(U(k, i) + U(k, i - 1)))                &
              & **4.0/(0.5*B(k, i))**2.0
         if(MANNINGS_N(jw))then
             hrad = BH(k, i)/(B(k, i) - B(k + 1, i) + 2.0*H(k, jw))
                                                            ! HRAD = BHR(K,I)/(BR(K,I)-BR(K+1,I)+2.0*H(K,JW))  SW 10/5/07
             if(macrophyte_on .AND. MANNINGS_N(jw))then
                 call MACROPHYTE_FRICTION(hrad, FRIC(i), effric, k, i)
                 gc2 = g*effric*effric/hrad**0.33333333
             elseif(.NOT.macrophyte_on .AND. MANNINGS_N(jw))then
                 gc2 = g*FRIC(i)*FRIC(i)/hrad**0.33333333
             endif
         endif
         prhk = gc2/(0.5*B(k, i))*ABS(0.5*(U(k, i) + U(k, i - 1)))**3.0
         unst = prdk - TKE(k, i, 2)
         unse = 1.44*TKE(k, i, 2)/TKE(k, i, 1)                                 &
              & *prdk - 1.92*(TKE(k, i, 2)/TKE(k, i, 1)*TKE(k, i, 2))
         TKE(k, i, 1) = TKE(k, i, 1) + dlt*(unst + prhk - bouk)
         TKE(k, i, 2) = TKE(k, i, 2) + dlt*(unse + prhe)
     enddo
     ustarb = SQRT(gc2)*ABS(0.5*(U(KB(i), i) + U(KB(i), i - 1)))
     TKE(KB(i), i, 1) = 0.5*(3.33*ustarb*ustarb + TKE(KB(i), i, 1))
     TKE(KB(i), i, 2) = 0.5*(ustarb*ustarb*ustarb*5.0/H(KB(i), jw)             &
                      & + TKE(KB(i), i, 2))
 
     do j = 1, 2                                              ! SG 10/4/07 Series of bug fixes for TKE
         k = kt
         AT(k, i) = 0.0
         CT(k, i) = 0.0
         VT(k, i) = 1.0
         DT(k, i) = TKE(k, i, j)
         do k = kt + 1, KB(i) - 1
             AT(k, i) = -dlt/BH1(k, i)*BB(k - 1, i)/sig(j)*AZ(k - 1, i)        &
                      & /AVH1(k - 1, i)
             CT(k, i) = -dlt/BH1(k, i)*BB(k, i)/sig(j)*AZ(k, i)/AVH1(k, i)
             VT(k, i) = 1.0 - AT(k, i) - CT(k, i)
             DT(k, i) = TKE(k, i, j)
         enddo
         k = KB(i)
         AT(k, i) = 0.0
         CT(k, i) = 0.0
         VT(k, i) = 1.0
         DT(k, i) = TKE(k, i, j)
         call TRIDIAG(AT(:, i), VT(:, i), CT(:, i), DT(:, i), kt, KB(i), kmx,  &
                    & TKE(:, i, j))
     enddo
 
     do k = kt, KB(i)
         TKE(k, i, 1) = MAX(TKE(k, i, 1), tkemin1)
         TKE(k, i, 2) = MAX(TKE(k, i, 2), tkemin2)
         AZT(k, i) = 0.09*TKE(k, i, 1)*TKE(k, i, 1)/TKE(k, i, 2)
     enddo
 
     do k = kt, KB(i) - 1
         AZ(k, i) = 0.5*(AZT(k, i) + AZT(k + 1, i))
         AZ(k, i) = MAX(azmin, AZ(k, i))
         AZ(k, i) = MIN(AZMAX(jw), AZ(k, i))
         DZ(k, i) = MAX(dzmin, frazdz*AZ(k, i))    ! No need to average DZ further since defined at cell, center bottom. AZ needs to be averaged to RHS of cell.
     enddo
     AZ(KB(i), i) = azmin
     AZT(KB(i), i) = azmin
     return
 
     entry CALCULATE_TKE1
!    Subroutine based on work of Gould(2006) RODI WITHOUT WIND
!    SW 6-28-04 DEADSEA  TKE
!    BASED ON SUBROUTINE FIRST DEVELOPED BY CHAPMAN AND COLE (1995)
!    BUG FIXES, CORRECTIONS, ADDITIONS SCOTT WELLS (2001)
 
!    BOTTOM LAYER BC
     if(ABS((U(KB(i),i) + U(KB(i),i-1))*.5)>nonzero)then
         ustarb = WALLFUNCTION()
!        ERROR TRAPPING
         if(ustarb==0.0)then
!	           IF (MANNINGS_N(JW)) THEN
!	           HRAD = BH1(KT,I)/(B(KTI(I),I)-B(KT+1,I)+2.*AVH1(KT,I))
!            if(macrophyte_on.and.mannings_n(jw))then
!            call macrophyte_friction(hrad,fric(i),effric,kt,i)
!            gc2=g*effric*effric/hrad**0.33333333
!            else if(.not.macrophyte_on.and.mannings_n(jw))then
!            gc2=g*fric(i)*fric(i)/hrad**0.33333333
!            end if
!            ELSE
!            GC2 = 0.0
!            IF (FRIC(I) /= 0.0) GC2 = G/(FRIC(I)*FRIC(I))
!            END IF
!	           USTARB = SQRT(GC2)*ABS(0.5*(U(KB(I),I)+U(KB(I),I-1)))
!            END ERROR TRAPPING
         endif
         TKE(KB(i), i, 1) = (3.33*ustarb*ustarb)
         TKE(KB(i), i, 2) = (ustarb*ustarb*ustarb)/(0.41*(H(KB(i), jw))*0.5)                  !/2.0
         USTARBTKE(i) = ustarb
         if(STRICKON(jw))call CALCFRIC()
     else
         ustarb = 0
         TKE(KB(i), i, 1) = tkemin1
         TKE(KB(i), i, 2) = tkemin2
         USTARBTKE(i) = ustarb
     endif
 
!INLET SET TO MOLECULAR CONDITIONS
     iut = FIRSTI(jw)
     if(i==iut)then
         do k = kt, KB(i)
             TKE(k, iut - 1, 2) = tkemin2
             TKE(k, iut - 1, 1) = tkemin1
         enddo
     endif
!OUTLET SET TO EQUAL TO I=I_OUTLET - 1
     idt = LASTI(jw)
     if(i==idt - 1)then
         do k = kt, KB(i)
             TKE(k, idt, 2) = TKE(k, idt - 1, 2)
             TKE(k, idt, 1) = TKE(k, idt - 1, 1)
         enddo
     endif
 
!CALCULATE HORIZONTAL PORTION OF THE SPLIT K-E EQUATIONS FOR INTERIOR LAYERS
     do k = kt + 1, KB(i) - 1
         if(MANNINGS_N(jw))then
             hrad = BH1(kt, i)/(B(KTI(i), i) - B(kt + 1, i) + 2.D0*AVH1(kt, i))   !HRAD = BHR(K,I)/(BR(K,I)-BR(K+1,I)+2.0*H(K,JW))
             gc2 = g*FRIC(i)*FRIC(i)/hrad**0.333
         else
             gc2 = 0.0
             if(FRIC(i)/=0.0)gc2 = g/(FRIC(i)*FRIC(i))
         endif
 
         bouk = DMAX1(AZ(k, i)*g*(RHO(k + 1, i) - RHO(k, i))/(H(k, jw)*rhow),  &
              & 0.0D0)
         prdk = AZ(k, i)                                                       &
              & *(0.5D0*(U(k, i) + U(k, i - 1) - U(k + 1, i) - U(k + 1, i - 1))&
              & /(H(k, jw)*0.5D0 + H(k + 1, jw)*0.5D0))**2.0                      !/2.0  /2.0
         if(.NOT.TKELATPRD(jw))then
             prhe = 0.0
             prhk = 0.0
         else
             prhe = TKELATPRDCONST(jw)                                         &
                  & *gc2**1.25*DABS(0.5D0*(U(k, i) + U(k, i - 1)))             &
                  & **4.0/(0.5D0*B(k, i))**2.0
             prhk = gc2/(0.5D0*B(k, i))*DABS(0.5D0*(U(k, i) + U(k, i - 1)))    &
                  & **3.0
         endif
         unst = prdk - TKE(k, i, 2)
         unse = 1.44D0*TKE(k, i, 2)/TKE(k, i, 1)                               &
              & *prdk - 1.92D0*(TKE(k, i, 2)/TKE(k, i, 1)*TKE(k, i, 2))
 
         udr = (1.0D0 + DSIGN(1.0D0, U(k, i)))*0.5D0
         udl = (1.0D0 + DSIGN(1.0D0, U(k, i - 1)))*0.5D0
         convkx = dlt/(DLX(i)*BH1(k, i))                                       &
                & *((U(k, i)*(udr*TKE(k,i,1) + (1.0-udr)*TKE(k,i+1,1))*BR(k, i)&
                & *H2(k, i) - U(k, i - 1)                                      &
                & *(udl*TKE(k,i-1,1) + (1.0-udl)*TKE(k,i,1))*BR(k, i - 1)      &
                & *H2(k, i - 1)))
         convex = dlt/(DLX(i)*BH1(k, i))                                       &
                & *((U(k, i)*(udr*TKE(k,i,2) + (1.0-udr)*TKE(k,i+1,2))*BR(k, i)&
                & *H2(k, i) - U(k, i - 1)                                      &
                & *(udl*TKE(k,i-1,2) + (1.0-udl)*TKE(k,i,2))*BR(k, i - 1)      &
                & *H2(k, i - 1)))
         if(IMPTKE(jw)=='     IMP')then
             TKE(k, i, 1) = TKE(k, i, 1) - convkx + dlt*(unst + prhk - bouk)
             TKE(k, i, 2) = TKE(k, i, 2) - convex + dlt*(unse + prhe)
         else
             ata = (1.D00 + DSIGN(1.0D0, W(k - 1, i)))*0.5D0
             ab = (1.0D0 + DSIGN(1.0D0, W(k, i)))*0.5D0
!
             convkz = dlt/(H2(k, i)*BH1(k, i))                                 &
                    & *((ab*TKE(k, i, 1) + (1.0 - ab)*TKE(k + 1, i, 1))*W(k, i)&
                    & *BB(k, i)*AVH2(k, i)                                     &
                    & - (ata*(TKE(k-1,i,1) + (1.0-ata)*TKE(k,i,1))*W(k - 1, i) &
                    & *BB(k - 1, i)*AVH2(k - 1, i)))
             convez = dlt/(H2(k, i)*BH1(k, i))                                 &
                    & *((ab*TKE(k, i, 2) + (1.0 - ab)*TKE(k + 1, i, 2))*W(k, i)&
                    & *BB(k, i)*AVH2(k, i)                                     &
                    & - (ata*(TKE(k-1,i,2) + (1.0-ata)*TKE(k,i,2))*W(k - 1, i) &
                    & *BB(k - 1, i)*AVH2(k - 1, i)))
 
             TKE(k, i, 1) = TKE(k, i, 1) - convkx - convkz +                   &
                          & dlt*(unst + prhk - bouk)
             TKE(k, i, 2) = TKE(k, i, 2) - convex - convez + dlt*(unse + prhe)
         endif
     enddo
     do
 
!CALCULATE TOP LAYER BOUNDARY CONDITION AND SOLVE THE K-E EQUATION USING TRIDIAGONAL SOLVER
         if(TKEBC(jw)==1)then
!RODI        NOWIND
             depthl = (ELWS(i) - EL(KB(i), i) + H2(KB(i), i)*COSA(jb))/COSA(jb)   !EL(KT,I)  -Z(I)  *COSA(JB)
             depthr = (ELWS(i + 1) - EL(KB(i + 1), i + 1) + H2(KB(i + 1), i)   &
                    & *COSA(jb))/COSA(jb)                                         !EL(KT,I+1)-Z(I+1)*COSA(JB)
             depth = (depthr + depthl)*0.5
             do j = 1, 2
                 k = kt
                 if(j==1)then
                     AT(k, i) = 0.0
                     CT(k, i) = -1.0
                     VT(k, i) = 1.0
                     DT(k, i) = 0.0
                 else
                     TKE(kt, i, 2) = TKE(k, i, 1)**(1.5)/(ARODI(jw)*depth)     &
                                   & *BH2(kt, i)/BH1(kt, i)                       !3.0/2.0
                     AT(k, i) = 0.0
                     CT(k, i) = 0.0
                     VT(k, i) = 1.0
                     DT(k, i) = TKE(kt, i, 2)
                 endif
                 do k = kt + 1, KB(i) - 1
                     if(IMPTKE(jw)=='     IMP')then
                         AT(k, i) = -dlt/BH1(k, i)                             &
                                  & *(BB(k - 1, i)/sig(j)*AZ(k - 1, i)         &
                                  & /AVH1(k - 1, i) + 0.5*W(k - 1, i)          &
                                  & *BB(k - 1, i))
                         CT(k, i) = -dlt/BH1(k, i)                             &
                                  & *(BB(k, i)/sig(j)*AZ(k, i)/AVH1(k, i)      &
                                  & - 0.5*W(k, i)*BB(k, i))
                         VT(k, i) = 1 + dlt/BH1(k, i)                          &
                                  & *(BB(k, i)/sig(j)*AZ(k, i)/AVH1(k, i)      &
                                  & + BB(k - 1, i)/sig(j)*AZ(k - 1, i)         &
                                  & /AVH1(k - 1, i)                            &
                                  & + 0.5*(W(k, i)*BB(k, i) - W(k - 1, i)      &
                                  & *BB(k - 1, i)))
                         DT(k, i) = TKE(k, i, j)
                     else
                         AT(k, i) = -dlt/BH1(k, i)*BB(k - 1, i)/sig(j)         &
                                  & *AZ(k - 1, i)/AVH1(k - 1, i)
                         CT(k, i) = -dlt/BH1(k, i)*BB(k, i)/sig(j)*AZ(k, i)    &
                                  & /AVH1(k, i)
                         VT(k, i) = 1.0 - AT(k, i) - CT(k, i)
                         DT(k, i) = TKE(k, i, j)
                     endif
                 enddo
                 k = KB(i)
                 AT(k, i) = 0.0
                 CT(k, i) = 0.0
                 VT(k, i) = 1.0
                 DT(k, i) = TKE(k, i, j)
                 call TRIDIAG(AT(:, i), VT(:, i), CT(:, i), DT(:, i), kt,      &
                            & KB(i), kmx, TKE(:, i, j))
             enddo
         elseif(TKEBC(jw)==2)then
!RODI        WIND
             depthl = (EL(kt, i) - Z(i)*COSA(jb) - EL(KB(i), i) + H2(KB(i), i) &
                    & *COSA(jb))/COSA(jb)
             depthr = (EL(kt, i + 1) - Z(i + 1)*COSA(jb) - EL(KB(i + 1), i + 1)&
                    & + H2(KB(i + 1), i)*COSA(jb))/COSA(jb)
             depth = (depthr + depthl)*0.5
             do jj = 1, 3
                 if(jj==1)then
                     j = 1
                     do k = kt + 1, KB(i) - 1
                         TKE(k, i, 3) = TKE(k, i, 1)
                     enddo
                 elseif(jj==2)then
                     j = 1
                 elseif(jj==3)then
                     j = 2
                 endif
!                ! WIND SHEAR AT SURFACE
                 ustar = SQRT(1.25*CZ(i)*WIND10(i)**2/RHO(kt, i))
                 if(jj==1)then
                     AT(kt, i) = 0.0
                     CT(kt, i) = -1.0
                     VT(kt, i) = 1.0
                     DT(kt, i) = 0.0
                 elseif(jj==2)then
                     if(TKE(kt, i, 1)*0.3<(ustar*ustar))then
                         TKE(kt, i, 1) = (ustar*ustar)*3.33333*BH2(kt, i)      &
                           & /BH1(kt, i)
                         AT(kt, i) = 0.0
                         CT(kt, i) = 0.0
                         VT(kt, i) = 1.0
                         DT(kt, i) = TKE(kt, i, j)
                         do k = kt + 1, KB(i) - 1
                             TKE(k, i, 1) = TKE(k, i, 3)
                         enddo
                     else
                         AT(kt, i) = 0.0
                         CT(kt, i) = -1.0
                         VT(kt, i) = 1.0
                         DT(kt, i) = 0.0
                     endif
                 endif
                 if(jj==3)then
                     TKE(kt, i, 2) = (TKE(kt, i, 1)*0.3)**(1.5)                &
                                   & /(0.41*(H2(kt, i)*0.5 + ARODI(jw)         &
                                   & *depth*(1 - (ustar*ustar)                 &
                                   & /(TKE(kt,i,1)*0.3))))*BH2(kt, i)          &
                                   & /BH1(kt, i)                                  !3.0/2.0    /2.0
                     AT(kt, i) = 0.0
                     CT(kt, i) = 0.0
                     VT(kt, i) = 1.0
                     DT(kt, i) = TKE(kt, i, 2)
                 endif
                 do k = kt + 1, KB(i) - 1
                     if(IMPTKE(jw)=='     IMP')then
                         AT(k, i) = -dlt/BH1(k, i)                             &
                                  & *(BB(k - 1, i)/sig(j)*AZ(k - 1, i)         &
                                  & /AVH1(k - 1, i) + 0.5*W(k - 1, i)          &
                                  & *BB(k - 1, i))
                         CT(k, i) = -dlt/BH1(k, i)                             &
                                  & *(BB(k, i)/sig(j)*AZ(k, i)/AVH1(k, i)      &
                                  & - 0.5*W(k, i)*BB(k, i))
                         VT(k, i) = 1 + dlt/BH1(k, i)                          &
                                  & *(BB(k, i)/sig(j)*AZ(k, i)/AVH1(k, i)      &
                                  & + BB(k - 1, i)/sig(j)*AZ(k - 1, i)         &
                                  & /AVH1(k - 1, i)                            &
                                  & + 0.5*(W(k, i)*BB(k, i) - W(k - 1, i)      &
                                  & *BB(k - 1, i)))
                         DT(k, i) = TKE(k, i, j)
                     else
                         AT(k, i) = -dlt/BH1(k, i)*BB(k - 1, i)/sig(j)         &
                                  & *AZ(k - 1, i)/AVH1(k - 1, i)
                         CT(k, i) = -dlt/BH1(k, i)*BB(k, i)/sig(j)*AZ(k, i)    &
                                  & /AVH1(k, i)
                         VT(k, i) = 1.0 - AT(k, i) - CT(k, i)
                         DT(k, i) = TKE(k, i, j)
                     endif
                 enddo
                 k = KB(i)
                 AT(k, i) = 0.0
                 CT(k, i) = 0.0
                 VT(k, i) = 1.0
                 DT(k, i) = TKE(k, i, j)
                 call TRIDIAG(AT(:, i), VT(:, i), CT(:, i), DT(:, i), kt,      &
                            & KB(i), kmx, TKE(:, i, j))
             enddo
         elseif(TKEBC(jw)==3)then
!CEQUALW2    WIND
!            TOP LAYER BC
             ustar = SQRT(1.25*CZ(i)*WIND10(i)**2/RHO(kt, i))
             if(MANNINGS_N(jw))then
                 hrad = BH1(kt, i)                                             &
                      & /(B(KTI(i), i) - B(kt + 1, i) + 2.*AVH1(kt, i))          !HRAD = BHR1(KT,I)/(BR(KTI(I),I)-BR(KT+1,I)+2.*AVH1(KT,I))
                 gc2 = g*FRIC(i)*FRIC(i)/hrad**0.333
             else
                 gc2 = 0.0
                 if(FRIC(i)/=0.0)gc2 = g/(FRIC(i)*FRIC(i))
             endif
             ustarbkt = SQRT(gc2)*ABS(0.5*(U(kt, i) + U(kt, i - 1)))
             TKE(kt, i, 1) = (3.33*(ustar*ustar + ustarbkt*ustarbkt))          &
                           & *(BH2(kt, i)/BH1(kt, i))
             TKE(kt, i, 2) = (ustar*ustar*ustar + ustarbkt*ustarbkt*ustarbkt)  &
                           & *5.0/H1(kt, i)*(BH2(kt, i)/BH1(kt, i))
             do j = 1, 2
                 k = kt
                 AT(k, i) = 0.0
                 CT(k, i) = 0.0
                 VT(k, i) = 1.0
                 DT(k, i) = TKE(k, i, j)
                 do k = kt + 1, KB(i) - 1
                     if(IMPTKE(jw)=='     IMP')then
                         AT(k, i) = -dlt/BH1(k, i)                             &
                                  & *(BB(k - 1, i)/sig(j)*AZ(k - 1, i)         &
                                  & /AVH1(k - 1, i) + 0.5*W(k - 1, i)          &
                                  & *BB(k - 1, i))
                         CT(k, i) = -dlt/BH1(k, i)                             &
                                  & *(BB(k, i)/sig(j)*AZ(k, i)/AVH1(k, i)      &
                                  & - 0.5*W(k, i)*BB(k, i))
                         VT(k, i) = 1 + dlt/BH1(k, i)                          &
                                  & *(BB(k, i)/sig(j)*AZ(k, i)/AVH1(k, i)      &
                                  & + BB(k - 1, i)/sig(j)*AZ(k - 1, i)         &
                                  & /AVH1(k - 1, i)                            &
                                  & + 0.5*(W(k, i)*BB(k, i) - W(k - 1, i)      &
                                  & *BB(k - 1, i)))
                         DT(k, i) = TKE(k, i, j)
                     else
                         AT(k, i) = -dlt/BH1(k, i)*BB(k - 1, i)/sig(j)         &
                                  & *AZ(k - 1, i)/AVH1(k - 1, i)
                         CT(k, i) = -dlt/BH1(k, i)*BB(k, i)/sig(j)*AZ(k, i)    &
                                  & /AVH1(k, i)
                         VT(k, i) = 1.0 - AT(k, i) - CT(k, i)
                         DT(k, i) = TKE(k, i, j)
                     endif
                 enddo
                 k = KB(i)
                 AT(k, i) = 0.0
                 CT(k, i) = 0.0
                 VT(k, i) = 1.0
                 DT(k, i) = TKE(k, i, j)
                 call TRIDIAG(AT(:, i), VT(:, i), CT(:, i), DT(:, i), kt,      &
                            & KB(i), kmx, TKE(:, i, j))
             enddo
         else
             warning_open = .TRUE.
             write(wrn, *)'UNKNOWN TKEBC OPTION CHECK W2_CON.NPT'
             write(wrn, *)'Setting TKE calculation to TKEBC(JW)=3'
             TKEBC(jw) = 3
             cycle
         endif
 
         do k = kt, KB(i)
             TKE(k, i, 1) = MAX(TKE(k, i, 1), tkemin1)
             TKE(k, i, 2) = MAX(TKE(k, i, 2), tkemin2)
             AZT(k, i) = 0.09*TKE(k, i, 1)*TKE(k, i, 1)/TKE(k, i, 2)
         enddo
 
         do k = kt, KB(i) - 1
             AZ(k, i) = 0.5*(AZT(k, i) + AZT(k + 1, i))
             AZ(k, i) = MAX(azmin, AZ(k, i))
             AZ(k, i) = MIN(AZMAX(jw), AZ(k, i))
             DZ(k, i) = MAX(dzmin, frazdz*AZ(k, i))
                                                 ! No need to average further, DZ is now defined at cell bottom, AZ will need to be averaged to the RHS of cell
         enddo
         AZ(KB(i), i) = azmin
         AZT(KB(i), i) = azmin
         exit
     enddo
     end subroutine calculate_az
