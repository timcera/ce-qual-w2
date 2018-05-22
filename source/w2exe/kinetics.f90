!*==kinetics.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                            S U B R O U T I N E   K I N E T I C S                                              **
!***********************************************************************************************************************************
 
     subroutine KINETICS
  !SP CEMA
  !End SP CEMA
 
! Type declarations
     use SCREENC
     use GLOBAL
     use KINETIC
     use GEOMC
     use TVDC
     use LOGICC
     use SURFHE
     use MACROPHYTEC
     use ZOOPLANKTONC
     use MAIN, ONLY:ngctdg, EPIPHYTON_CALC, BOD_CALC, ALG_CALC, BOD_CALCN,     &
       & BOD_CALCP, po4_calc, n_calc, dsi_calc, sedcomp_exist, jg_age,         &
       & water_age_active
     use CEMAVARS
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     integer :: J, Jbod, Jg, Llm
     intent (in) Llm
!
! Local variables
!
     real, save :: algex, algn, algp, alkt, ammt, bicart, blim, bmass,         &
                 & bmasstest, bodtot, cart, cbodct, cbodnt, cbodpt, cbodset,   &
                 & ch4ex, co2ex, co3t, coef, coef1, coef2, coef3, coef4, colb, &
                 & coldep, cvol, dh1, dh2, dh3, dhh, dosat, ea, epsilon, f,    &
                 & fdpo4, fetch, h2co3t, h2po4t, h2sex, h3po4t, hco3t, hion,   &
                 & hpo4t, hs, ht, incr, k1, k2, kamm, kp1, kp2, kp3, kw, l, l0,&
                 & l1, lam1, lam2, lavg, light, limit, ltcoefm, macext,        &
                 & macext1, n2sat, nh3t, nh4pr, nh4t, no3pr, o2ex, oh, oht,    &
                 & omct, phost, pht, po4t, s2, sdalgc, sdalgn, sdalgp
     real, allocatable, dimension(:, :, :), save :: atrm, atrmf, atrmr, etrm,  &
       & etrmf, etrmr
     real, allocatable, dimension(:, :), save :: bibh2, dom, lam2m, nh4bod,    &
       & nh4trm, no3trm, omtrm, po4bod, pom, sodtrm, ticbod
     logical, save :: femn
     integer, save :: ibod, iter, ja, jaf, jcb, jcg, jd, je, jj, jjz, js, jt,  &
                    & k, m, mi, n
     real, save :: sdbodc, sdbodn, sdbodp, sdepc, sdepn, sdepp, sedem, sedsi,  &
                 & sedsic, sedsidk, sedsin, sedsip, sedso, sedsoc, sedson,     &
                 & sedsop, sedsum, sedsumk, sqrs2, ssext, ssr, t1k, tau,       &
                 & ticch4, tmac, totmac, totss0, tpss, ts, u2, uorb, xdum, xx, &
                 & zgztot, zminfac, zooext, zoon, zoop
!
!*** End of declarations rewritten by SPAG
!
 
                                                                                                               ! SW 10/17/15
                                                                                      ! CB 4/20/11
                                                                             ! CB 6/6/10
  ! enhanced pH buffering start
  ! enhanced pH buffering end
!    CEMA start
!    CEMA end
 
!    Allocation declarations
 
     allocate(omtrm(kmx, imx), sodtrm(kmx, imx), nh4trm(kmx, imx),             &
            & no3trm(kmx, imx), dom(kmx, imx), pom(kmx, imx))          ! SW 4/5/09
     allocate(po4bod(kmx, imx), nh4bod(kmx, imx), ticbod(kmx, imx))
     allocate(atrm(kmx, imx, nal), atrmr(kmx, imx, nal), atrmf(kmx, imx, nal))
     allocate(etrm(kmx, imx, nep), etrmr(kmx, imx, nep), etrmf(kmx, imx, nep))
     allocate(lam2m(kmx, kmx), bibh2(kmx, imx))
     return
 
!***********************************************************************************************************************************
!**  T E M P E R A T U R E  R A T E  M U L T I P L I E R S                    
!***********************************************************************************************************************************
 
!    **
     entry TEMPERATURE_RATES
     do i = iu, id
         do k = kt, KB(i)
             lam1 = FR(T1(k, i), NH4T1(jw), NH4T2(jw), NH4K1(jw), NH4K2(jw))
             nh4trm(k, i) = lam1/(1.0 + lam1 - NH4K1(jw))
             lam1 = FR(T1(k, i), NO3T1(jw), NO3T2(jw), NO3K1(jw), NO3K2(jw))
             no3trm(k, i) = lam1/(1.0 + lam1 - NO3K1(jw))
             lam1 = FR(T1(k, i), OMT1(jw), OMT2(jw), OMK1(jw), OMK2(jw))
             omtrm(k, i) = lam1/(1.0 + lam1 - OMK1(jw))
             lam1 = FR(T1(k, i), SODT1(jw), SODT2(jw), SODK1(jw), SODK2(jw))
             sodtrm(k, i) = lam1/(1.0 + lam1 - SODK1(jw))
             do ja = 1, nal
                 if(ALG_CALC(ja))then
                     lam1 = FR(T1(k, i), AT1(ja), AT2(ja), AK1(ja), AK2(ja))
                     lam2 = FF(T1(k, i), AT3(ja), AT4(ja), AK3(ja), AK4(ja))
                     atrmr(k, i, ja) = lam1/(1.0 + lam1 - AK1(ja))
                     atrmf(k, i, ja) = lam2/(1.0 + lam2 - AK4(ja))
                     atrm(k, i, ja) = atrmr(k, i, ja)*atrmf(k, i, ja)
                 endif
             enddo
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))then
                     lam1 = FR(T1(k, i), ET1(je), ET2(je), EK1(je), EK2(je))
                     lam2 = FF(T1(k, i), ET3(je), ET4(je), EK3(je), EK4(je))
                     etrmr(k, i, je) = lam1/(1.0 + lam1 - EK1(je))
                     etrmf(k, i, je) = lam2/(1.0 + lam2 - EK4(je))
                     etrm(k, i, je) = etrmr(k, i, je)*etrmf(k, i, je)
                 endif
             enddo
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     lam1 = FR(T1(k, i), MT1(m), MT2(m), MK1(m), MK2(m))
                     lam2 = FF(T1(k, i), MT3(m), MT4(m), MK3(m), MK4(m))
                     MACTRMR(k, i, m) = lam1/(1.0 + lam1 - MK1(m))
                     MACTRMF(k, i, m) = lam2/(1.0 + lam2 - MK4(m))
                     MACTRM(k, i, m) = MACTRMR(k, i, m)*MACTRMF(k, i, m)
                 endif
             enddo
             if(zooplankton_calc)then
                 do jz = 1, nzp
                     lam1 = FR(T1(k, i), ZT1(jz), ZT2(jz), ZK1(jz), ZK2(jz))
                     lam2 = FF(T1(k, i), ZT3(jz), ZT4(jz), ZK3(jz), ZK4(jz))
                     ZOORMR(k, i, jz) = lam1/(1. + lam1 - ZK1(jz))
                     ZOORMF(k, i, jz) = lam2/(1. + lam2 - ZK4(jz))
                     ZOORM(k, i, jz) = ZOORMR(k, i, jz)*ZOORMF(k, i, jz)
                 enddo
             endif
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  K I N E T I C   R A T E S                                                
!***********************************************************************************************************************************
 
!    **
     entry KINETIC_RATES
 
!    Decay rates
!!$OMP PARALLEL DO
     do i = iu, id
         do k = kt, KB(i)
             DO1(k, i) = O2(k, i)/(O2(k, i) + kdo)
             DO2(k, i) = 1.0 - DO1(k, i)                         !O2(K,I)/(O2(K,I)+KDO)
             DO3(k, i) = (1.0 + SIGN(1.0, O2(k, i) - 1.E-10))*0.5
             SEDD(k, i) = sodtrm(k, i)*SDKV(k, i)*SED(k, i)*DO3(k, i)      !CB 10/22/06
             SEDDP(k, i) = sodtrm(k, i)*SDKV(k, i)*SEDP(k, i)*DO3(k, i)
             SEDDN(k, i) = sodtrm(k, i)*SDKV(k, i)*SEDN(k, i)*DO3(k, i)
             SEDDC(k, i) = sodtrm(k, i)*SDKV(k, i)*SEDC(k, i)*DO3(k, i)
!            Amaila start
             if(sedcomp_exist)then
                 SEDD1(k, i) = sodtrm(k, i)*SDK1(jw)*SED1(k, i)*DO3(k, i)
                 SEDD2(k, i) = sodtrm(k, i)*SDK2(jw)*SED2(k, i)*DO3(k, i)
             endif
!            amaila end
             SEDBR(k, i) = SEDB(jw)*SED(k, i)                               !CB 11/30/06
             SEDBRP(k, i) = SEDB(jw)*SEDP(k, i)                             !CB 11/30/06
             SEDBRN(k, i) = SEDB(jw)*SEDN(k, i)                             !CB 11/30/06
             SEDBRC(k, i) = SEDB(jw)*SEDC(k, i)                             !CB 11/30/06
             NH4D(k, i) = nh4trm(k, i)*NH4DK(jw)*NH4(k, i)*DO1(k, i)
             NO3D(k, i) = no3trm(k, i)*NO3DK(jw)*NO3(k, i)*DO2(k, i)
             LDOMD(k, i) = omtrm(k, i)*LDOMDK(jw)*LDOM(k, i)*DO3(k, i)
             RDOMD(k, i) = omtrm(k, i)*RDOMDK(jw)*RDOM(k, i)*DO3(k, i)
             LPOMD(k, i) = omtrm(k, i)*LPOMDK(jw)*LPOM(k, i)*DO3(k, i)
             RPOMD(k, i) = omtrm(k, i)*RPOMDK(jw)*RPOM(k, i)*DO3(k, i)
             LRDOMD(k, i) = omtrm(k, i)*LRDDK(jw)*LDOM(k, i)*DO3(k, i)
             LRPOMD(k, i) = omtrm(k, i)*LRPDK(jw)*LPOM(k, i)*DO3(k, i)
             cbodd(k, i, 1:nbod) = kbod(1:nbod)*tbod(1:nbod)**(T1(k, i) - 20.0)&
                                 & *DO3(k, i)
             if(k==KB(i))then  ! SW 4/18/07
                 SODD(k, i) = SOD(i)/BH2(k, i)*sodtrm(k, i)*BI(k, i)
             else
                 SODD(k, i) = SOD(i)/BH2(k, i)*sodtrm(k, i)                    &
                            & *(BI(k, i) - BI(k + 1, i))
             endif
 
!            Inorganic suspended solids settling rates - P adsorption onto SS
!            and Fe
             FPSS(k, i) = PARTP(jw)/(PARTP(jw)*TISS(k, i) + PARTP(jw)*FE(k, i) &
                        & *DO1(k, i) + 1.0)
             FPFE(k, i) = PARTP(jw)*FE(k, i)                                   &
                        & /(PARTP(jw)*TISS(k, i) + PARTP(jw)*FE(k, i)*DO1(k, i)&
                        & + 1.0)
             SSSI(k, i) = SSSO(k - 1, i)
             totss0 = 0.0
             do js = 1, nss
                 totss0 = totss0 + SSS(js)*FPSS(k, i)*SS(k, i, js)
             enddo
             SSSO(k, i) = (totss0 + FES(jw)*FPFE(k, i))*BI(k, i)/BH2(k, i)     &
                        & *DO1(k, i)                                                  ! SW 11/7/07
             FPSS(k, i) = FPSS(k, i)*TISS(k, i)
 
!            OM stoichiometry
             ORGPLD(k, i) = 0.0
             ORGPRD(k, i) = 0.0
             ORGPLP(k, i) = 0.0
             ORGPRP(k, i) = 0.0
             ORGNLD(k, i) = 0.0
             ORGNRD(k, i) = 0.0
             ORGNLP(k, i) = 0.0
             ORGNRP(k, i) = 0.0
             if(CAC(nldomp)/='      ON')then
                 ORGPLD(k, i) = ORGP(jw)
             elseif(LDOM(k, i)>0.0)then
                 ORGPLD(k, i) = LDOMP(k, i)/LDOM(k, i)
             else
                 ORGPLD(k, i) = ORGP(jw)
             endif
             if(CAC(nrdomp)/='      ON')then
                 ORGPRD(k, i) = ORGP(jw)
             elseif(RDOM(k, i)>0.0)then
                 ORGPRD(k, i) = RDOMP(k, i)/RDOM(k, i)
             else
                 ORGPRD(k, i) = ORGP(jw)
             endif
             if(CAC(nlpomp)/='      ON')then
                 ORGPLP(k, i) = ORGP(jw)
             elseif(LPOM(k, i)>0.0)then
                 ORGPLP(k, i) = LPOMP(k, i)/LPOM(k, i)
             else
                 ORGPLP(k, i) = ORGP(jw)
             endif
             if(CAC(nrpomp)/='      ON')then
                 ORGPRP(k, i) = ORGP(jw)
             elseif(RPOM(k, i)>0.0)then
                 ORGPRP(k, i) = RPOMP(k, i)/RPOM(k, i)
             else
                 ORGPRP(k, i) = ORGP(jw)
             endif
             if(CAC(nldomn)/='      ON')then
                 ORGNLD(k, i) = ORGN(jw)
             elseif(LDOM(k, i)>0.0)then
                 ORGNLD(k, i) = LDOMN(k, i)/LDOM(k, i)
             else
                 ORGNLD(k, i) = ORGN(jw)
             endif
             if(CAC(nrdomn)/='      ON')then
                 ORGNRD(k, i) = ORGN(jw)
             elseif(RDOM(k, i)>0.0)then
                 ORGNRD(k, i) = RDOMN(k, i)/RDOM(k, i)
             else
                 ORGNRD(k, i) = ORGN(jw)
             endif
             if(CAC(nlpomn)/='      ON')then
                 ORGNLP(k, i) = ORGN(jw)
             elseif(LPOM(k, i)>0.0)then
                 ORGNLP(k, i) = LPOMN(k, i)/LPOM(k, i)
             else
                 ORGNLP(k, i) = ORGN(jw)
             endif
             if(CAC(nrpomn)/='      ON')then  ! SR 8/2/2017
                 ORGNRP(k, i) = ORGN(jw)
             elseif(RPOM(k, i)>0.0)then
                 ORGNRP(k, i) = RPOMN(k, i)/RPOM(k, i)
             else
                 ORGNRP(k, i) = ORGN(jw)
             endif
 
!            Light Extinction Coefficient
             if(.NOT.READ_EXTINCTION(jw))then
                 algex = 0.0
                 ssext = 0.0
                 zooext = 0.0                                                                    ! SW 11/8/07
                 do ja = 1, nal
                     if(ALG_CALC(ja))algex = algex + EXA(ja)*ALG(k, i, ja)
                 enddo
                 do js = 1, nss
                     ssext = ssext + EXSS(jw)*SS(k, i, js)
                 enddo
 !       TOTMAC=0.0                                                                ! SW 4/20/11 Delete this section?
 !       DO M=1,NMC
 !         IF(MACROPHYTE_CALC(JW,M))THEN
 !           JT=KTI(I)
 !           JE=KB(I)
 !           DO JJ=JT,JE
 !             TOTMAC = EXM(M)*MACRM(JJ,K,I,M)+TOTMAC
 !           END DO
 !         END IF
 !       END DO
 !       MACEXT=TOTMAC/(BH2(K,I)*DLX(I))
 
                 if(zooplankton_calc)then
                     do jz = 1, nzp
                         zooext = zooext + ZOO(k, i, jz)*EXZ(jz)
                     enddo
                 endif
!
                 GAMMA(k, i) = EXH2O(jw) + ssext + EXOM(jw)                    &
                             & *(LPOM(k, i) + RPOM(k, i)) + algex + zooext                       ! sw 4/21/11
!
                 if(nmc>0)then
                             ! cb 4/20/11
                     macext1 = 0.0
                             ! cb 4/20/11
                     if(KTICOL(i))then
                         jt = KTI(i)
                     else
                         jt = KTI(i) + 1
                     endif
                     je = KB(i)
                     do jj = jt, je
                         totmac = 0.0
                         do m = 1, nmc
                             if(MACROPHYTE_CALC(jw, m))totmac = EXM(m)         &
                              & *MACRM(jj, k, i, m) + totmac
                         enddo
                         if(CW(jj, i)>0.0)then
                             macext = totmac/(CW(jj, i)*DLX(i)*H2(k, i))
                         else
                             macext = 0.0
                         endif
                         GAMMAJ(jj, k, i) = GAMMA(k, i) + macext ! SW 4/20/11
                         macext1 = macext*CW(jj, i) + macext1
                                                 ! cb 4/20/11
                     enddo
                     GAMMA(k, i) = GAMMA(k, i) + macext1/B(jt, i)                        ! SW 4/21/11
                 endif
             else
                 GAMMA(k, i) = EXH2O(jw)
             endif
 
!            Zooplankton Rates
             if(zooplankton_calc)then
                 do jz = 1, nzp
                     TGRAZE(k, i, jz) = PREFP(jz)*LPOM(k, i)
                     do jjz = 1, nzp
                         TGRAZE(k, i, jz) = TGRAZE(k, i, jz) + PREFZ(jjz, jz)  &
                           & *ZOO(k, i, jjz)                                    !CB 5/17/2007
                     enddo
                     do ja = 1, nal
                         if(ALG_CALC(ja))TGRAZE(k, i, jz) = PREFA(ja, jz)      &
                          & *ALG(k, i, ja) + TGRAZE(k, i, jz)
                     enddo
                     zminfac = (1.0 + SIGN(1.0, ZOO(k, i, jz) - ZOOMIN(jz)))   &
                             & *0.5
                     ZRT(k, i, jz) = ZOORMR(k, i, jz)*ZR(jz)*zminfac*DO3(k, i)
                     if(TGRAZE(k, i, jz)<=0.0 .OR. O2(k, i)<2.0)then
                         ZMU(k, i, jz) = 0.0
                         agz(k, i, 1:nal, jz) = 0.0
                         zgz(k, i, jz, :) = 0.0
                         if(O2(k, i)<2.0)zminfac = 2*zminfac
                     else
                         ZMU(k, i, jz)                                         &
                           & = MAX(ZOORM(k, i, jz)*ZG(jz)*(TGRAZE(k, i, jz)    &
                           & - ZOOMIN(jz))/(TGRAZE(k, i, jz) + ZS2P(jz)), 0.0)
                         do ja = 1, nal
                             if(ALG_CALC(ja))agz(k, i, ja, jz) = ZMU(k, i, jz) &
                              & *ZOO(k, i, jz)                                 &
                              & *(ALG(k, i, ja)*PREFA(ja, jz)/TGRAZE(k, i, jz))                                                   !  KV 5/26/2007
                         enddo
                         do jjz = 1, nzp
                         ! OMNIVOROUS ZOOPLANKTON
                             zgz(k, i, jjz, jz) = ZMU(k, i, jz)*ZOO(k, i, jz)  &
                               & *(ZOO(k, i, jjz)*PREFZ(jjz, jz)               &
                               & /TGRAZE(k, i, jz))                                                      !KV 5/26/2007
                         enddo
                     endif
                     ZMT(k, i, jz) = MAX(1.0 - ZOORMF(k, i, jz), 0.02)*ZM(jz)  &
                                   & *zminfac
                 enddo
               ! ZOOP LOOP
             endif
 
         enddo
           ! K LOOP
     enddo ! I LOOP
!!$OMP END PARALLEL DO
 
!    Algal rates
     do ja = 1, nal
         if(ALG_CALC(ja))then
             do i = iu, id
!****            Limiting factor
                 light = (1.0 - BETA(jw))*SRON(jw)*SHADE(i)/ASAT(ja)
                 lam1 = light
                 lam2 = light
                 do k = kt, KB(i)
 
!******              Limiting factor
                     lam1 = lam2
                     lam2 = lam1*EXP( - GAMMA(k, i)*H2(k, i))
                     fdpo4 = 1.0 - FPSS(k, i) - FPFE(k, i)
                     ALLIM(k, i, ja) = 2.718282*(EXP( - lam2) - EXP( - lam1))  &
                                     & /(GAMMA(k, i)*H2(k, i))
                     if(AHSP(ja)/=0.0 .AND. po4_calc)APLIM(k, i, ja)           &
                      & = fdpo4*PO4(k, i)                                      &
                      & /(fdpo4*PO4(k, i) + AHSP(ja) + nonzero)                                                      ! cb 10/12/11
                     if(AHSN(ja)/=0.0 .AND. n_calc)ANLIM(k, i, ja)             &
                      & = (NH4(k, i) + NO3(k, i))                              &
                      & /(NH4(k, i) + NO3(k, i) + AHSN(ja) + nonzero)                                                ! cb 10/12/11
                     if(AHSSI(ja)/=0.0 .AND. dsi_calc)ASLIM(k, i, ja)          &
                      & = DSI(k, i)/(DSI(k, i) + AHSSI(ja) + nonzero)                                                ! cb 10/12/11
                     limit = MIN(APLIM(k, i, ja), ANLIM(k, i, ja),             &
                           & ASLIM(k, i, ja), ALLIM(k, i, ja))
 
!******              Algal rates
                     AGR(k, i, ja) = atrm(k, i, ja)*AG(ja)*limit
                     ARR(k, i, ja) = atrm(k, i, ja)*AR(ja)*DO3(k, i)
                     AMR(k, i, ja) = (atrmr(k, i, ja) + 1.0 - atrmf(k, i, ja)) &
                                   & *AM(ja)
                     AER(k, i, ja) = MIN((1.0 - ALLIM(k, i, ja))*AE(ja)*atrm(k,&
                                   & i, ja), AGR(k, i, ja))
                     if(AS(ja)>=0.0)then
                         if(k==kt)then
                             ASR(k, i, ja) = AS(ja)*( - ALG(k, i, ja))*BI(k, i)&
                               & /BH2(k, i)
                         else
                             ASR(k, i, ja) = AS(ja)                            &
                               & *(ALG(k - 1, i, ja) - ALG(k, i, ja))*BI(k, i) &
                               & /BH2(k, i)
                         endif
                     elseif(k==KB(i))then
                         ASR(k, i, ja) = -AS(ja)                               &
                           & *( - ALG(k, i, ja)*BI(k, i)/BH2(k, i))                                           !SW 11/8/07
                     elseif(k==kt)then
                         ASR(k, i, ja) = -AS(ja)*ALG(k + 1, i, ja)*BI(k + 1, i)&
                           & *DLX(i)/VOL(k, i)                                                               !SW 11/8/07
                     else
                         ASR(k, i, ja) = -AS(ja)                               &
                           & *(ALG(k + 1, i, ja)*BI(k + 1, i)/BH2(k, i)        &
                           & - ALG(k, i, ja)*BI(k, i)/BH2(k, i))                                              !SP 8/27/07
                     endif
                 enddo
             enddo
         endif
     enddo  ! ALGAE LOOP
 
!    Macrophyte Light/Nutrient Limitation and kinetic rates
     do m = 1, nmc
         mgr(:, :, iu:id, m) = 0.0
         mrr(:, iu:id, m) = 0.0
         mmr(:, iu:id, m) = 0.0                                 ! cb 3/8/16
         if(MACROPHYTE_CALC(jw, m))then
             do i = iu, id
                 ltcoefm = (1.0 - BETA(jw))*SRON(jw)*SHADE(i)
                 if(KTICOL(i))then
                     jt = KTI(i)
                 else
                     jt = KTI(i) + 1
                 endif
                 je = KB(i)
                 do jj = jt, je
                     lam1 = ltcoefm
                     lam2m(jj, kt) = lam1*EXP( - GAMMAJ(jj, kt, i)*H2(kt, i))
                     lavg = (lam1 - lam2m(jj, kt))                             &
                          & /(GAMMAJ(jj, kt, i)*H2(kt, i))
                     MLLIM(jj, kt, i, m) = lavg/(lavg + MSAT(m))
                     if(MHSP(m)/=0.0 .AND. PSED(m)<1.0)then
                         MPLIM(kt, i, m) = fdpo4*PO4(kt, i)                    &
                           & /(fdpo4*PO4(kt, i) + MHSP(m) + nonzero)
                     else
                         MPLIM(kt, i, m) = 1.0
                     endif
                     if(MHSN(m)/=0.0 .AND. NSED(m)<1.0)then
                         MNLIM(kt, i, m) = NH4(kt, i)                          &
                           & /(NH4(kt, i) + MHSN(m) + nonzero)
                     else
                         MNLIM(kt, i, m) = 1.0
                     endif
                     if(MHSC(m)/=0.0)MCLIM(kt, i, m) = CO2(kt, i)              &
                      & /(CO2(kt, i) + MHSC(m) + nonzero)
                     limit = MIN(MPLIM(kt, i, m), MNLIM(kt, i, m),             &
                           & MCLIM(kt, i, m), MLLIM(jj, kt, i, m))
 
!*************       sources/sinks
 
                     mgr(jj, kt, i, m) = MACTRM(kt, i, m)*MG(m)*limit
 
                 enddo
 
                 mrr(kt, i, m) = MACTRM(kt, i, m)*MR(m)*DO3(kt, i)
                 mmr(kt, i, m) = (MACTRMR(kt, i, m) + 1.0 - MACTRMF(kt, i, m)) &
                               & *MM(m)
 
                 do k = kt + 1, KB(i)
                     jt = k
                     je = KB(i)
                     do jj = jt, je
                         lam1 = lam2m(jj, k - 1)
                         lam2m(jj, k) = lam1*EXP( - GAMMAJ(jj, k, i)*H2(k, i))
                         lavg = (lam1 - lam2m(jj, k))                          &
                              & /(GAMMAJ(jj, k, i)*H2(k, i))
                         MLLIM(jj, k, i, m) = lavg/(lavg + MSAT(m))
                         if(MHSP(m)/=0.0 .AND. PSED(m)<1.0)then
                             MPLIM(k, i, m) = fdpo4*PO4(k, i)                  &
                               & /(fdpo4*PO4(k, i) + MHSP(m) + nonzero)
                         else
                             MPLIM(k, i, m) = 1.0
                         endif
                         if(MHSN(m)/=0.0 .AND. NSED(m)<1.0)then
                             MNLIM(k, i, m) = NH4(k, i)                        &
                               & /(NH4(k, i) + MHSN(m) + nonzero)
                         else
                             MNLIM(k, i, m) = 1.0
                         endif
                         if(MHSC(m)/=0.0)MCLIM(k, i, m) = CO2(k, i)            &
                          & /(CO2(k, i) + MHSC(m) + nonzero)
                         limit = MIN(MPLIM(k, i, m), MNLIM(k, i, m),           &
                               & MCLIM(k, i, m), MLLIM(jj, k, i, m))
 
!*************           sources/sinks
 
                         mgr(jj, k, i, m) = MACTRM(k, i, m)*MG(m)*limit
 
                     enddo
 
                     mrr(k, i, m) = MACTRM(k, i, m)*MR(m)*DO3(k, i)
                     mmr(k, i, m) = (MACTRMR(k, i, m) + 1.0 - MACTRMF(k, i, m))&
                                  & *MM(m)
                 enddo
             enddo
         endif
     enddo
 
     return
 
!***********************************************************************************************************************************
!**  G E N E R I C   C O N S T I T U E N T                                    
!***********************************************************************************************************************************
 
!    **
     entry GENERIC_CONST(Jg)
 
     if(water_age_active .AND. Jg==jg_age)then
                                             ! SW 7/27/2017  Speed of computation
         do i = iu, id
             do k = kt, KB(i)
                 CGSS(k, i, Jg) = -CG0DK(Jg)
             enddo
         enddo
     else
 
         xx = 0.0
         femn = .FALSE.
         if(sd_global)then
             sdinfeooh(:, iu:id) = 0.0
             sdinmno2(:, iu:id) = 0.0
             if(Jg==ngfe2 .OR. Jg==ngfeooh .OR. Jg==ngmn2 .OR. Jg==ngmno2)     &
              & femn = .TRUE.
             if(Jg==ngh2s)then
                 h2sd(:, iu:id) = 0.0
                 h2sreaer(:, iu:id) = 0.0
             endif
             if(Jg==ngch4)then
                 ch4d(:, iu:id) = 0.0
                 ch4reaer(:, iu:id) = 0.0
             endif
             if(Jg==ngfe2)fe2d(:, iu:id) = 0.0
             if(Jg==ngmn2)mn2d(:, iu:id) = 0.0
         endif
!        CEMA end
         do i = iu, id
             light = (1.0 - BETA(jw))*SRON(jw)*SHADE(i) !LCJ 2/26/15
             lam1 = light
             lam2 = light
 
             do k = kt, KB(i)
                 lam1 = lam2
                 lam2 = lam1*EXP( - GAMMA(k, i)*H2(k, i))
                 light = lam1*(1. - EXP( - GAMMA(k, i)*H2(k, i)))              &
                       & /(GAMMA(k, i)*H2(k, i))                                     ! SW 10/17/15
 
!                CEMA start
                 if(femn)then
                     if(Jg==ngfeooh)then
                         if(k==kt)then
                             xx = fesetvel*(CG(k - 1, i, Jg) - CG(k, i, Jg))   &
                                & *BI(k, i)/BH2(k, i)
                         else
                             xx = fesetvel*(CG(k - 1, i, Jg) - CG(k, i, Jg))   &
                                & *BI(k, i)/BH2(k, i)
                         endif
                         sdinfeooh(k, i) = xx
                     endif
                     if(Jg==ngmno2)then
                         if(k==kt)then
                             xx = mnsetvel*(CG(k - 1, i, Jg) - CG(k, i, Jg))   &
                                & *BI(k, i)/BH2(k, i)
                         else
                             xx = mnsetvel*(CG(k - 1, i, Jg) - CG(k, i, Jg))   &
                                & *BI(k, i)/BH2(k, i)
                         endif
                         sdinmno2(k, i) = xx
                     endif
!                    CEMA end
                 elseif(CGS(Jg)>0.0)then
                     if(k==kt)then
                         xx = CGS(Jg)*(CG(k - 1, i, Jg) - CG(k, i, Jg))        &
                            & *BI(k, i)/BH2(k, i)                     ! AS(JA)*(-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
                     else
                         xx = CGS(Jg)*(CG(k - 1, i, Jg) - CG(k, i, Jg))        &
                            & *BI(k, i)/BH2(k, i)                      !AS(JA)*(ALG(K-1,I,JA)-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
                     endif
                 elseif(CGS(Jg)<0.0)then
                     if(k==KB(i))then
                         xx = -CGS(Jg)*( - CG(k, i, Jg))*BI(k, i)/BH2(k, i)
                                                            !-AS(JA)*(-ALG(K,I,JA)  *BI(K,I)/BH2(K,I))                                           !SW 11/8/07
                     elseif(k==kt)then
                         xx = -CGS(Jg)*CG(k + 1, i, Jg)*BI(k + 1, i)*DLX(i)    &
                            & /VOL(k, i)                            !-AS(JA)* ALG(K+1,I,JA)*BI(K+1,I)*DLX(I)/VOL(K,I)                                   !SW 11/8/07
                     else
                         xx = -CGS(Jg)*(CG(k + 1, i, Jg)*BI(k + 1, i)/BH2(k, i)&
                            & - CG(k, i, Jg)*BI(k, i)/BH2(k, i))                           !-AS(JA)*(ALG(K+1,I,JA)*BI(K+1,I)/BH2(K,I)-ALG(K,I,JA)*BI(K,I)/BH2(K,I))             !SP 8/27/07
                     endif
                 endif
                ! CEMA
 
 
!                CEMA start
                 if(femn)then
                     if(Jg==ngfeooh)CGSS(k, i, Jg) = kfe_oxid*O2(k, i)         &
                      & *10**(2.0*(PH(k, i) - 7.0))*CG(k, i, ngfe2)            &
                      & - kfe_red*(kfeooh_halfsat/(O2(k, i) + kfeooh_halfsat)) &
                      & *CG(k, i, Jg) + xx
                     if(Jg==ngfe2)then
                         fe2d(k, i) = -kfe_oxid*O2(k, i)                       &
                                    & *10**(2.0*(PH(k, i) - 7.0))*CG(k, i, Jg)
                         CGSS(k, i, Jg) = fe2d(k, i)                           &
                           & + kfe_red*(kfeooh_halfsat/(O2(k, i)               &
                           & + kfeooh_halfsat))*CG(k, i, ngfeooh) + xx
                     endif
                     if(Jg==ngmno2)CGSS(k, i, Jg) = kmn_oxid*O2(k, i)          &
                      & *10**(2.0*(PH(k, i) - 7.0))*CG(k, i, ngmn2)            &
                      & - kmn_red*(kmno2_halfsat/(O2(k, i) + kmno2_halfsat))   &
                      & *CG(k, i, Jg) + xx
                     if(Jg==ngmn2)then
                         mn2d(k, i) = -kmn_oxid*O2(k, i)                       &
                                    & *10**(2.0*(PH(k, i) - 7.0))*CG(k, i, Jg)
                         CGSS(k, i, Jg) = mn2d(k, i)                           &
                           & + kmn_red*(kmno2_halfsat/(O2(k, i)                &
                           & + kmno2_halfsat))*CG(k, i, ngmno2) + xx
                     endif
!                    CEMA end
 
                 elseif(CGQ10(Jg)/=0.0)then
!                    CGSS(K,I,JG) =
! -CG0DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)-CG1DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)*CG(K
                     if(Jg==ngh2s)then
                         h2sd(k, i) = ( - CG0DK(Jg)*CGQ10(Jg)**(T1(k, i) - 20.0&
                                    & ) - CG1DK(Jg)*CGQ10(Jg)                  &
                                    & **(T1(k, i) - 20.0)*CG(k, i, Jg))        &
                                    & *DO3(k, i)
                         CGSS(k, i, Jg) = h2sd(k, i)
                     elseif(Jg==ngch4)then
                         ch4d(k, i) = ( - CG0DK(Jg)*CGQ10(Jg)**(T1(k, i) - 20.0&
                                    & ) - CG1DK(Jg)*CGQ10(Jg)                  &
                                    & **(T1(k, i) - 20.0)*CG(k, i, Jg))        &
                                    & *DO3(k, i)
                         CGSS(k, i, Jg) = ch4d(k, i)
                     else
                         CGSS(k, i, Jg) = -CG0DK(Jg)*CGQ10(Jg)                 &
                           & **(T1(k, i) - 20.0) - CG1DK(Jg)*CGQ10(Jg)         &
                           & **(T1(k, i) - 20.0)*CG(k, i, Jg) - CGLDK(Jg)      &
                           & *light*CG(k, i, Jg) + xx
                     endif
!                    CEMA end
!                    CGSS(K,I,JG) = -CG0DK(JG)-CG1DK(JG)*CG(K,I,JG)+xx        
!                    ! SW 4/5/09
!                    CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I) CEMA
!                    start
                 elseif(Jg==ngh2s)then
                     h2sd(k, i) = -CG0DK(Jg) - CG1DK(Jg)*CG(k, i, Jg)*DO3(k, i)
                     CGSS(k, i, Jg) = h2sd(k, i) + xx
                 elseif(Jg==ngch4)then
                     ch4d(k, i) = -CG0DK(Jg) - CG1DK(Jg)*CG(k, i, Jg)*DO3(k, i)
                     CGSS(k, i, Jg) = ch4d(k, i) + xx
                 else
                     CGSS(k, i, Jg) = -CG0DK(Jg) - CG1DK(Jg)*CG(k, i, Jg)      &
                                    & - CGLDK(Jg)*light*CG(k, i, Jg) + xx
!                    CEMA end
                 endif
                ! CEMA
                 if(Jg==ngso4)CGSS(k, i, Jg) = CGSS(k, i, Jg) - h2sd(k, i)
                                                                 ! sulfate production from sulfide decay, negative sign because 'h2sd' is negative
             enddo
!            CEMA start
 
             if(Jg==ngh2s)then
                 if(.NOT.ICE(i))then
                     if(REAER(i)==0.0)call GAS_TRANSFER
                     h2sex = REAER(i)*0.984
                     h2sreaer(kt, i) = h2sex*( - CG(kt, i, Jg))*BI(kt, i)      &
                                     & /BH2(kt, i)
                     CGSS(kt, i, Jg) = CGSS(kt, i, Jg) + h2sreaer(kt, i)
                 endif
             endif
             if(Jg==ngch4)then
                 if(.NOT.ICE(i))then
                     if(REAER(i)==0.0)call GAS_TRANSFER
                     ch4ex = REAER(i)*1.188
                     ch4reaer(kt, i) = ch4ex*( - CG(kt, i, Jg))*BI(kt, i)      &
                                     & /BH2(kt, i)
                     CGSS(kt, i, Jg) = CGSS(kt, i, Jg) + ch4reaer(kt, i)
                 endif
             endif
             if(CGKLF(Jg)/=0.0)then
                 if(.NOT.ICE(i))then
                     if(REAER(i)==0.0)call GAS_TRANSFER
                     if(CGCS(Jg)== - 1.0)then
                                      ! THIS IS FOR N2 GAS
                         ea = DEXP                                             &
                            & (2.3026D0*(7.5D0*TDEW(jw)/(TDEW(jw) + 237.3D0)   &
                            & + 0.6609D0))*0.001316                                      ! in mm Hg   0.0098692atm=7.5006151mmHg
            ! PN2=0.79*(PALT(I)-EA)   ! atm with water vapor correction since 0.79 atm is for dry air
                         n2sat = 1.5568D06*0.79*(PALT(i) - ea)                 &
                               & *(1.8816D-5 - 4.116D-7*T1(kt, i)              &
                               & + 4.6D-9*T1(kt, i)**2)
     ! N2SS(KT,I) = (N2SAT-N2(KT,I))*REAER(I)*1.304*BI(KT,I)/BH2(KT,I)
                         CGSS(kt, i, Jg) = CGSS(kt, i, Jg) + REAER(i)*CGKLF(Jg)&
                           & *(n2sat - CG(kt, i, Jg))*BI(kt, i)/BH2(kt, i)                            ! fixed value of KLN2=1.034*KLO2
                         do k = kt, KB(i)        ! NOTE THERE IS 1 TIME STEP LAG WITH TDG SINCE PLACED HERE
                             dosat = SATO(T1(k, i), TDS(k, i), PALT(i),        &
                                   & SALT_WATER(jw))
                             TDG(k, i)                                         &
                               & = 100.*((0.79*CG(k, i, ngctdg)/n2sat) + O2(k, &
                               & i)/dosat*0.21)                                     !(1.5568D06*0.79*(PALT(I)-EA))*(1.8816D-5 - 4.116D-7 * T1(KT,I) + 4.6D-9 * T1(KT,I)**2))
                         enddo
 
                     else
 
                         CGSS(kt, i, Jg) = CGSS(kt, i, Jg) + REAER(i)*CGKLF(Jg)&
                           & *(CGCS(Jg) - CG(kt, i, Jg))*BI(kt, i)/BH2(kt, i)
                     endif
                 endif
             endif
 
         enddo
     endif
 
     return
 
 
 
!    IF (CGQ10(JG) /= 0.0) THEN
!    DO I=IU,ID
!    DO K=KT,KB(I)
!    CGSS(K,I,JG) =
! -CG0DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)-CG1DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)*CG(K
!    END DO
!    ELSE
!    DO I=IU,ID
!    DO K=KT,KB(I)
!    CGSS(K,I,JG) = -CG0DK(JG)-CG1DK(JG)*CG(K,I,JG)+                          
!    ! SW 4/5/09 CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I) END DO
!    END DO
!    END IF
!RETURN
 
!***********************************************************************************************************************************
!**  S U S P E N D E D   S O L I D S                                          
!***********************************************************************************************************************************
 
!    **
     entry SUSPENDED_SOLIDS(J)
 
    !SP CEMA
     if(sediment_diagenesis)then
        !All resuspension done in CEMA code
         if(includebedconsolidation)SEDIMENT_RESUSPENSION(J) = .FALSE.
     endif
    !End SP CEMA
 
     do i = iu, id
         ssr = 0.0
         if(SEDIMENT_RESUSPENSION(J))then
             fetch = FETCHD(i, jb)
             if(COS(PHI(jw) - PHI0(i))<0.0)fetch = FETCHU(i, jb)
             fetch = MAX(fetch, BI(kt, i), DLX(i))
             u2 = WIND(jw)*WSC(i)*WIND(jw)*WSC(i) + nonzero
             coef1 = 0.53*(g*DEPTHB(kt, i)/u2)**0.75
             coef2 = 0.0125*(g*fetch/u2)**0.42
             coef3 = 0.833*(g*DEPTHB(kt, i)/u2)**0.375
             coef4 = 0.077*(g*fetch/u2)**0.25
             hs = 0.283*u2/g*0.283*TANH(coef1)*TANH(coef2/TANH(coef1))
      !TS    = 2.0*PI*U2/G*1.2*  TANH(COEF3)*TANH(COEF4/TANH(COEF3))
             ts = 2.0*pi*SQRT(u2)/g*1.2*TANH(coef3)*TANH(coef4/TANH(coef3)) ! cb 7/15/14
             l0 = g*ts*ts/(2.0*pi)
         endif
         SSSS(kt, i, J) = -SSS(J)*SS(kt, i, J)*BI(kt, i)/BH2(kt, i) + ssr
 !   DO K=KT-1,KB(I)-1                                             ! SW 4/3/09   KT,KB
         do k = kt + 1, KB(i) - 1      ! cb 9/29/14
             if(SEDIMENT_RESUSPENSION(J))then
                 l1 = l0
                 l = l0*TANH(2.0*pi*DEPTHB(k, i)/l1)
                 do while (ABS(l - l1)>0.001)
                     l1 = l
                     l = l0*TANH(2.0*pi*DEPTHB(k, i)/l1)
                 enddo
                 coef = MIN(710.0, 2.0*pi*DEPTHB(k, i)/l)
                 uorb = pi*hs/ts*100.0/SINH(coef)
                 tau = 0.003*uorb*uorb
                 if(tau - TAUCR(J)>0.0)                                        &
                  & epsilon = MAX(0.0, 0.008/49.0*(tau - TAUCR(J))             &
                  & **3*10000.0/dlt)
                 if(k==KB(i))then    ! SW 4/18/07
                     ssr = epsilon*DLX(i)*BI(k, i)/VOL(k, i)
                 else
                     ssr = epsilon*DLX(i)*(BI(k, i) - BI(k + 1, i))/VOL(k, i)
                 endif
             endif
             SSSS(k, i, J) = SSS(J)*(SS(k - 1, i, J) - SS(k, i, J))*BI(k, i)   &
                           & /BH2(k, i) + ssr
         enddo
         if(SEDIMENT_RESUSPENSION(J))ssr = epsilon*DLX(i)*BI(KB(i), i)         &
          & /VOL(KB(i), i)
!CEMA    SP
         if(sediment_diagenesis)then
             if(includefftlayer .AND. fftactive)SSSS(KB(i), i, J)              &
              & = (SSS(J)*SS(KB(i) - 1, i, J) - fftlayersettvel*SS(KB(i), i, J)&
              & )/H(KB(i), jw) + ssr
             if(includefftlayer .AND. .NOT.fftactive)SSSS(KB(i), i, J) = 0.D0
  ! Flocculation              !SR                                                      !New section on flocculation          !SR 04/21/13
    !DO K=KT,KB(I)
    !  SSF = 0.0
    !  IF (J > 1 .AND. SSFLOC(J-1) > 0.0) THEN
    !    IF (FLOCEQN(J-1) == 0) THEN
    !      SSF = MIN(SSFLOC(J-1), SS(K,I,J-1)/DLT)
    !    ELSE IF (FLOCEQN(J-1) == 1) THEN
    !      SSF = SSFLOC(J-1)*SS(K,I,J-1)
    !    ELSE IF (FLOCEQN(J-1) == 2) THEN
    !      SSF = SSFLOC(J-1)*SS(K,I,J-1)*SS(K,I,J-1)
    !    END IF
    !  END IF
    !  IF (J < NSS .AND. SSFLOC(J) > 0.0) THEN
    !    IF (FLOCEQN(J) == 0) THEN
    !      SSF = SSF - MIN(SSFLOC(J), SS(K,I,J)/DLT)
    !    ELSE IF (FLOCEQN(J) == 1) THEN
    !      SSF = SSF - SSFLOC(J)*SS(K,I,J)
    !    ELSE IF (FLOCEQN(J) == 2) THEN
    !      SSF = SSF - SSFLOC(J)*SS(K,I,J)*SS(K,I,J)
    !    END IF
    !  END IF
    !  SSSS(K,I,J) = SSSS(K,I,J) + SSF
    !END DO                                                                        !End new section on flocculation      !SR 04/21/13
             if(.NOT.includefftlayer)SSSS(KB(i), i, J) = SSS(J)                &
              & *(SS(KB(i) - 1, i, J) - SS(KB(i), i, J))/H(KB(i), jw) + ssr
         endif
    !End CEMA SP
     enddo
     return
 
!***********************************************************************************************************************************
!**  P H O S P H O R U S                                                     
!***********************************************************************************************************************************
 
!    **
     entry PHOSPHORUS
     po4ar(:, iu:id) = 0.0
     po4ag(:, iu:id) = 0.0
     po4er(:, iu:id) = 0.0
     po4eg(:, iu:id) = 0.0
     po4bod(:, iu:id) = 0.0
     po4mr(:, iu:id) = 0.0
     po4mg(:, iu:id) = 0.0
     po4zr(:, iu:id) = 0.0
 
     do i = iu, id
         do k = kt, KB(i)
             do jcb = 1, nbod
!                IF(BOD_CALC(JCB))PO4BOD(K,I) =
!                PO4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODP(JCB)
                 if(BOD_CALCP(jcb))then                                        ! cb 5/19/11
                     po4bod(k, i) = po4bod(k, i) + cbodd(k, i, jcb)            &
                                  & *CBODP(k, i, jcb)
                 else
                     po4bod(k, i) = po4bod(k, i) + cbodd(k, i, jcb)            &
                                  & *CBOD(k, i, jcb)*BODP(jcb)
                 endif
             enddo
             do ja = 1, nal
                 if(ALG_CALC(ja))then
                     po4ag(k, i) = po4ag(k, i) + AGR(k, i, ja)*ALG(k, i, ja)   &
                                 & *AP(ja)
                     po4ar(k, i) = po4ar(k, i) + ARR(k, i, ja)*ALG(k, i, ja)   &
                                 & *AP(ja)
                 endif
             enddo
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))then
                     po4eg(k, i) = po4eg(k, i) + EGR(k, i, je)*EPC(k, i, je)   &
                                 & *EP(je)
                     po4er(k, i) = po4er(k, i) + ERR(k, i, je)*EPC(k, i, je)   &
                                 & *EP(je)
                 endif
             enddo
             PO4EP(k, i) = po4er(k, i) - po4eg(k, i)
             PO4AP(k, i) = po4ar(k, i) - po4ag(k, i)
             PO4POM(k, i) = ORGPLP(k, i)*LPOMD(k, i) + ORGPRP(k, i)*RPOMD(k, i)
             PO4DOM(k, i) = ORGPLD(k, i)*LDOMD(k, i) + ORGPRD(k, i)*RDOMD(k, i)
             PO4OM(k, i) = PO4POM(k, i) + PO4DOM(k, i)
             if(sedcomp_exist)then ! SW 5/26/15
           ! PO4SD(K,I)  = SEDDp(K,I)+sedd1(k,i)*orgp(jw) + sedd2(k,i)*orgp(jw)   ! Amaila
                 PO4SD(k, i) = SEDDP(k, i) + SEDD1(k, i)*PBIOM(jw)             &
                             & + SEDD2(k, i)*PBIOM(jw)                             ! Amaila, cb 6/7/17
             else
                 PO4SD(k, i) = SEDDP(k, i)
             endif
             PO4SR(k, i) = PO4R(jw)*SODD(k, i)*DO2(k, i)
             PO4NS(k, i) = SSSI(k, i)*PO4(k - 1, i) - SSSO(k, i)*PO4(k, i)
 
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     if(k==kt)then
                         jt = KTI(i)
                     else
                         jt = k
                     endif
                     je = KB(i)
                     do jj = jt, je
                         po4mg(k, i) = po4mg(k, i) + mgr(jj, k, i, m)          &
                                     & *MACRM(jj, k, i, m)*MP(m)               &
                                     & *(1.0 - PSED(m))
                         po4mr(k, i) = po4mr(k, i) + mrr(k, i, m)              &
                                     & *MACRM(jj, k, i, m)*MP(m)
                     enddo
                 endif
             enddo
             po4mr(k, i) = po4mr(k, i)/(DLX(i)*BH(k, i))
             po4mg(k, i) = po4mg(k, i)/(DLX(i)*BH(k, i))
             if(zooplankton_calc)then
                 do jz = 1, nzp
                     po4zr(k, i) = po4zr(k, i) + ZRT(k, i, jz)*ZOO(k, i, jz)   &
                                 & *ZP(jz)
                 enddo
             endif
 
 
             PO4SS(k, i) = PO4AP(k, i) + PO4EP(k, i) + PO4OM(k, i)             &
                         & + PO4SD(k, i) + PO4SR(k, i) + PO4NS(k, i)           &
                         & + po4bod(k, i) + po4mr(k, i) - po4mg(k, i)          &
                         & + po4zr(k, i)
 
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  A M M O N I U M                                                        **
!***********************************************************************************************************************************
 
     entry AMMONIUM
     nh4ag(:, iu:id) = 0.0
     nh4ar(:, iu:id) = 0.0
     nh4er(:, iu:id) = 0.0
     nh4eg(:, iu:id) = 0.0
     nh4bod(:, iu:id) = 0.0
     nh4mg(:, iu:id) = 0.0
     nh4mr(:, iu:id) = 0.0
     nh4zr(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             do jcb = 1, nbod
!                IF(BOD_CALC(JCB))NH4BOD(K,I) = 
!                NH4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODN(JCB)
                 if(BOD_CALCN(jcb))then                                        ! cb 5/19/11
                     nh4bod(k, i) = nh4bod(k, i) + cbodd(k, i, jcb)            &
                                  & *CBODN(k, i, jcb)
                 else
                     nh4bod(k, i) = nh4bod(k, i) + cbodd(k, i, jcb)            &
                                  & *CBOD(k, i, jcb)*BODN(jcb)
                 endif
             enddo
             do ja = 1, nal
                 if(ALG_CALC(ja))then
                     if(ANEQN(ja)==2)then
                         nh4pr = NH4(k, i)*NO3(k, i)                           &
                               & /((ANPR(ja) + NH4(k, i))*(ANPR(ja) + NO3(k, i)&
                               & )) + NH4(k, i)*ANPR(ja)                       &
                               & /((NO3(k, i) + NH4(k, i) + nonzero)           &
                               & *(ANPR(ja) + NO3(k, i)))
                     else
                         nh4pr = NH4(k, i)/(NH4(k, i) + NO3(k, i) + nonzero)
                     endif
                     if(AHSN(ja)>0.0)nh4ag(k, i) = nh4ag(k, i) + AGR(k, i, ja) &
                      & *ALG(k, i, ja)*AN(ja)*nh4pr
                     nh4ar(k, i) = nh4ar(k, i) + ARR(k, i, ja)*ALG(k, i, ja)   &
                                 & *AN(ja)
                 endif
             enddo
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))then
                     if(ENEQN(je)==2)then
                         nh4pr = NH4(k, i)*NO3(k, i)                           &
                               & /((ENPR(je) + NH4(k, i))*(ENPR(je) + NO3(k, i)&
                               & )) + NH4(k, i)*ENPR(je)                       &
                               & /((NO3(k, i) + NH4(k, i) + nonzero)           &
                               & *(ENPR(je) + NO3(k, i)))
                     else
                         nh4pr = NH4(k, i)/(NH4(k, i) + NO3(k, i) + nonzero)
                     endif
                     if(EHSN(je)>0.0)nh4eg(k, i) = nh4eg(k, i) + EGR(k, i, je) &
                      & *EPC(k, i, je)*EN(je)*nh4pr
                     nh4er(k, i) = nh4er(k, i) + ERR(k, i, je)*EPC(k, i, je)   &
                                 & *EN(je)
                 endif
             enddo
             NH4EP(k, i) = nh4er(k, i) - nh4eg(k, i)
             NH4AP(k, i) = nh4ar(k, i) - nh4ag(k, i)
 
             NH4DOM(k, i) = LDOMD(k, i)*ORGNLD(k, i) + RDOMD(k, i)*ORGNRD(k, i)
             NH4POM(k, i) = LPOMD(k, i)*ORGNLP(k, i) + RPOMD(k, i)*ORGNRP(k, i)
 
             NH4OM(k, i) = NH4DOM(k, i) + NH4POM(k, i)
 
             if(sedcomp_exist)then ! SW 5/26/15
            !NH4SD(K,I)  =  SEDDn(K,I) +sedd1(k,i)*orgn(jw) + sedd2(k,i)*orgn(jw)   ! Amaila
                 NH4SD(k, i) = SEDDN(k, i) + SEDD1(k, i)*NBIOM(jw)             &
                             & + SEDD2(k, i)*NBIOM(jw)                               ! Amaila, cb 6/7/17
             else
                 NH4SD(k, i) = SEDDN(k, i)
             endif
 
             NH4SR(k, i) = NH4R(jw)*SODD(k, i)*DO2(k, i)
 
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     if(k==kt)then
                         jt = KTI(i)
                     else
                         jt = k
                     endif
                     je = KB(i)
                     do jj = jt, je
                         nh4mr(k, i) = nh4mr(k, i) + mrr(k, i, m)              &
                                     & *MACRM(jj, k, i, m)*MN(m)
                         nh4mg(k, i) = nh4mg(k, i) + mgr(jj, k, i, m)          &
                                     & *MACRM(jj, k, i, m)*MN(m)               &
                                     & *(1.0 - NSED(m))
                     enddo
                 endif
             enddo
             nh4mr(k, i) = nh4mr(k, i)/(DLX(i)*BH(k, i))
             nh4mg(k, i) = nh4mg(k, i)/(DLX(i)*BH(k, i))
             if(zooplankton_calc)then
                 do jz = 1, nzp
                     nh4zr(k, i) = nh4zr(k, i) + ZRT(k, i, jz)*ZOO(k, i, jz)   &
                                 & *ZN(jz)
                 enddo
             endif
             NH4SS(k, i) = NH4AP(k, i) + NH4EP(k, i) + NH4OM(k, i)             &
                         & + NH4SD(k, i) + NH4SR(k, i) + nh4bod(k, i)          &
                         & - NH4D(k, i) + nh4mr(k, i) - nh4mg(k, i)            &
                         & + nh4zr(k, i)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  N I T R A T E                                                        **
!***********************************************************************************************************************************
 
     entry NITRATE
     no3ag(:, iu:id) = 0.0
     no3eg(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             do ja = 1, nal
                 if(ALG_CALC(ja))then
                     no3pr = 1.0 - NH4(k, i)/(NH4(k, i) + NO3(k, i) + nonzero)
                     if(ANEQN(ja)==2)no3pr = 1.0 -                             &
                      & (NH4(k, i)*NO3(k, i)/((ANPR(ja) + NH4(k,i))            &
                      & *(ANPR(ja) + NO3(k,i))) + NH4(k, i)*ANPR(ja)           &
                      & /((NO3(k,i) + NH4(k,i) + nonzero)*(ANPR(ja) + NO3(k,i))&
                      & ))
                     if(AHSN(ja)>0.0)no3ag(k, i) = no3ag(k, i) + AGR(k, i, ja) &
                      & *ALG(k, i, ja)*no3pr*AN(ja)
                 endif
             enddo
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))then
                     no3pr = 1.0 - NH4(k, i)/(NH4(k, i) + NO3(k, i) + nonzero)
                     if(ENEQN(je)==2)no3pr = 1.0 -                             &
                      & (NH4(k, i)*NO3(k, i)/((ENPR(je) + NH4(k,i))            &
                      & *(ENPR(je) + NO3(k,i))) + NH4(k, i)*ENPR(je)           &
                      & /((NO3(k,i) + NH4(k,i) + nonzero)*(ENPR(je) + NO3(k,i))&
                      & ))
                     if(EHSN(je)>0.0)no3eg(k, i) = no3eg(k, i) + EGR(k, i, je) &
                      & *EPC(k, i, je)*no3pr*EN(je)
                 endif
             enddo
             if(k==KB(i))then  ! SW 4/18/07
                 NO3SED(k, i) = NO3(k, i)*NO3S(jw)*no3trm(k, i)*(BI(k, i))     &
                              & /BH2(k, i)
             else
                 NO3SED(k, i) = NO3(k, i)*NO3S(jw)*no3trm(k, i)                &
                              & *(BI(k, i) - BI(k + 1, i))/BH2(k, i)
             endif
             NO3SS(k, i) = NH4D(k, i) - NO3D(k, i) - no3ag(k, i) - no3eg(k, i) &
                         & - NO3SED(k, i)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  D I S S O L V E D   S I L I C A                                          
!***********************************************************************************************************************************
 
!    **
     entry DISSOLVED_SILICA
     dsiag(:, iu:id) = 0.0
     dsieg(:, iu:id) = 0.0                                            !; DSIBOD = 0.0
     do i = iu, id
         do k = kt, KB(i)
             do ja = 1, nal
                 if(ALG_CALC(ja))dsiag(k, i) = dsiag(k, i) + AGR(k, i, ja)     &
                  & *ALG(k, i, ja)*ASI(ja)
             enddo
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))dsieg(k, i) = dsieg(k, i)           &
                  & + EGR(k, i, je)*EPC(k, i, je)*ESI(je)
             enddo
             DSID(k, i) = PSIDK(jw)*PSI(k, i)
             DSISD(k, i) = SEDD(k, i)*ORGSI(jw)
             DSISR(k, i) = DSIR(jw)*SODD(k, i)*DO2(k, i)
             DSIS(k, i) = (SSSI(k, i)*DSI(k - 1, i) - SSSO(k, i)*DSI(k, i))    &
                        & *PARTSI(jw)
             DSISS(k, i) = DSID(k, i) + DSISD(k, i) + DSISR(k, i) + DSIS(k, i) &
                         & - dsiag(k, i) - dsieg(k, i)                                 !+DSIBOD
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  P A R T I C U L A T E   S I L I C A                                      
!***********************************************************************************************************************************
 
!    **
     entry PARTICULATE_SILICA
     psiam(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             do ja = 1, nal
                 if(ALG_CALC(ja))psiam(k, i) = psiam(k, i) + AMR(k, i, ja)     &
                  & *ALG(k, i, ja)*ASI(ja)                          !   PSI(K,I)   HA-Z  12/2016
             enddo
             PSID(k, i) = PSIDK(jw)*PSI(k, i)
             PSINS(k, i) = PSIS(jw)*(PSI(k - 1, i)*DO1(k - 1, i) - PSI(k, i)   &
                         & *DO1(k, i))*BI(k, i)/BH2(k, i)
             PSISS(k, i) = psiam(k, i) - PSID(k, i) + PSINS(k, i)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  I R O N                                                            **
!***********************************************************************************************************************************
 
     entry IRON
     do i = iu, id
         do k = kt, KB(i)
             FENS(k, i) = FES(jw)                                              &
                        & *(FE(k - 1, i)*DO1(k - 1, i) - FE(k, i)*DO1(k, i))   &
                        & *BI(k, i)/BH2(k, i)
             FESR(k, i) = FER(jw)*SODD(k, i)*DO2(k, i)
             FESS(k, i) = FESR(k, i) + FENS(k, i)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  L A B I L E   D O M                                                     **
!***********************************************************************************************************************************
 
     entry LABILE_DOM
     ldomap(:, iu:id) = 0.0
     ldomep(:, iu:id) = 0.0
     ldommac(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             do ja = 1, nal
                 if(ALG_CALC(ja))ldomap(k, i) = ldomap(k, i)                   &
                  & + (AER(k, i, ja) + (1.0 - APOM(ja))*AMR(k, i, ja))         &
                  & *ALG(k, i, ja)
             enddo
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))ldomep(k, i) = ldomep(k, i)         &
                  & + (EER(k, i, je) + (1.0 - EPOM(je))*EMR(k, i, je))         &
                  & *EPC(k, i, je)
             enddo
 
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     if(k==kt)then
                         jt = KTI(i)
                     else
                         jt = k
                     endif
                     je = KB(i)
                     do jj = jt, je
                         ldommac(k, i) = ldommac(k, i) + (1.0 - MPOM(m))       &
                           & *mmr(k, i, m)*MACRM(jj, k, i, m)
                     enddo
                 endif
             enddo
             ldommac(k, i) = ldommac(k, i)/(DLX(i)*BH(k, i))
             LDOMSS(k, i) = ldomap(k, i) + ldomep(k, i) - LDOMD(k, i)          &
                          & - LRDOMD(k, i) + ldommac(k, i)
 
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  R E F R A C T O R Y   D O M                                              
!***********************************************************************************************************************************
 
!    **
     entry REFRACTORY_DOM
     do i = iu, id
         do k = kt, KB(i)
             RDOMSS(k, i) = LRDOMD(k, i) - RDOMD(k, i)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  L A B I L E   P O M                                                     
!***********************************************************************************************************************************
 
!    **
     entry LABILE_POM
     lpomap(:, iu:id) = 0.0
     lpomep(:, iu:id) = 0.0
     lpommac(:, iu:id) = 0.0
     lpzooin(:, iu:id) = 0.0
     lpzooout(:, iu:id) = 0.0                                                                                          ! cb 5/19/06
     do i = iu, id
         do k = kt, KB(i)
             do ja = 1, nal
                 if(ALG_CALC(ja))lpomap(k, i) = lpomap(k, i) + APOM(ja)        &
                  & *(AMR(k, i, ja)*ALG(k, i, ja))
             enddo
             do je = 1, nep                                                ! cb 5/19/06
                 if(EPIPHYTON_CALC(jw, je))lpomep(k, i) = lpomep(k, i)         &
                  & + EPOM(je)*(EMR(k, i, je)*EPC(k, i, je))                                         ! cb 5/19/06
             enddo                                                         ! cb 5/19/06
             LPOMNS(k, i) = POMS(jw)*(LPOM(k - 1, i) - LPOM(k, i))*BI(k, i)    &
                          & /BH2(k, i)
 
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     jt = k
                     je = KB(i)
                     do jj = jt, je
                         lpommac(k, i) = lpommac(k, i) + MPOM(m)*LRPMAC(m)     &
                           & *mmr(k, i, m)*MACRM(jj, k, i, m)
                     enddo
                 endif
             enddo
             lpommac(k, i) = lpommac(k, i)/(DLX(i)*BH(k, i))
             if(zooplankton_calc)then
                 do jz = 1, nzp
                     if(TGRAZE(k, i, jz)>0.0)then
                         lpzooout(k, i) = lpzooout(k, i) + ZOO(k, i, jz)       &
                           & *(ZMT(k, i, jz)                                   &
                           & + (ZMU(k, i, jz) - (ZMU(k,i,jz)*ZEFF(jz))))
                         lpzooin(k, i) = lpzooin(k, i) + ZOO(k, i, jz)         &
                           & *PREFP(jz)*ZMU(k, i, jz)*LPOM(k, i)               &
                           & /TGRAZE(k, i, jz)
                     else
                         lpzooout(k, i) = lpzooout(k, i) + ZOO(k, i, jz)       &
                           & *(ZMT(k, i, jz)                                   &
                           & + (ZMU(k, i, jz) - (ZMU(k,i,jz)*ZEFF(jz))))
                         lpzooin(k, i) = 0.0
                     endif
                 enddo
             endif
             LPOMSS(k, i) = lpomap(k, i) + lpomep(k, i) - LPOMD(k, i)          &
                          & + LPOMNS(k, i) - LRPOMD(k, i) + lpommac(k, i)      &
                          & + lpzooout(k, i) - lpzooin(k, i)                                                                 ! cb 5/19/06
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  R E F R A C T O R Y   P O M                                              
!***********************************************************************************************************************************
 
!    **
     entry REFRACTORY_POM
     rpommac(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             RPOMNS(k, i) = POMS(jw)*(RPOM(k - 1, i) - RPOM(k, i))*BI(k, i)    &
                          & /BH2(k, i)
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     jt = k
                     je = KB(i)
                     do jj = jt, je
                         rpommac(k, i) = rpommac(k, i) + MPOM(m)               &
                           & *(1.0 - LRPMAC(m))*mmr(k, i, m)*MACRM(jj, k, i, m)
                     enddo
                 endif
             enddo
             rpommac(k, i) = rpommac(k, i)/(DLX(i)*BH(k, i))
             RPOMSS(k, i) = LRPOMD(k, i) + RPOMNS(k, i) - RPOMD(k, i)          &
                          & + rpommac(k, i)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  A L G A E                                                             **
!***********************************************************************************************************************************
 
     entry ALGAE(J)
     agzt(:, iu:id, J) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             if(zooplankton_calc)then
                 do jz = 1, nzp
                     agzt(k, i, J) = agzt(k, i, J) + agz(k, i, J, jz)     ! CB 5/26/07
                 enddo
             endif
             ASS(k, i, J) = ASR(k, i, J)                                       &
                          & + (AGR(k, i, J) - AER(k, i, J) - AMR(k, i, J)      &
                          & - ARR(k, i, J))*ALG(k, i, J) - agzt(k, i, J)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  B I O C H E M I C A L   O 2   D E M A N D                                
!***********************************************************************************************************************************
 
!    **
     entry BIOCHEMICAL_O2_DEMAND(Jbod)
     if(Jbod==1)cbodns(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             cbodset = CBODS(Jbod)*(CBOD(k - 1, i, Jbod) - CBOD(k, i, Jbod))   &
                     & *BI(k, i)/BH2(k, i)
             cbodns(k, i) = cbodns(k, i) + cbodset
             CBODSS(k, i, Jbod) = -cbodd(k, i, Jbod)*CBOD(k, i, Jbod) + cbodset
         enddo
     enddo
     return
 
!    VARIABLE STOCHIOMETRY FOR CBOD SECTION ! CB 6/6/10
!***********************************************************************************************************************************
!**  B I O C H E M I C A L   O 2   D E M A N D   P H O S P H O R U S          
!***********************************************************************************************************************************
 
!    **
     entry BIOCHEMICAL_O2_DEMAND_P(Jbod)
     if(Jbod==1)cbodnsp(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             cbodset = CBODS(Jbod)*(CBODP(k - 1, i, Jbod) - CBODP(k, i, Jbod)) &
                     & *BI(k, i)/BH2(k, i)
             cbodnsp(k, i) = cbodnsp(k, i) + cbodset
             CBODPSS(k, i, Jbod) = -cbodd(k, i, Jbod)*CBODP(k, i, Jbod)        &
                                 & + cbodset
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  B I O C H E M I C A L   O 2   D E M A N D   N I T R O G E N              
!***********************************************************************************************************************************
 
!    **
     entry BIOCHEMICAL_O2_DEMAND_N(Jbod)
     if(Jbod==1)cbodnsn(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             cbodset = CBODS(Jbod)*(CBODN(k - 1, i, Jbod) - CBODN(k, i, Jbod)) &
                     & *BI(k, i)/BH2(k, i)
             cbodnsn(k, i) = cbodnsn(k, i) + cbodset
             CBODNSS(k, i, Jbod) = -cbodd(k, i, Jbod)*CBODN(k, i, Jbod)        &
                                 & + cbodset
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  D I S S O L V E D   O X Y G E N                                          
!***********************************************************************************************************************************
 
!    **
     entry DISSOLVED_OXYGEN
     doap(:, iu:id) = 0.0
     doar(:, iu:id) = 0.0
     doep(:, iu:id) = 0.0
     doer(:, iu:id) = 0.0
     dobod(:, iu:id) = 0.0
     domp(:, iu:id) = 0.0
     domr(:, iu:id) = 0.0
     dozr(:, iu:id) = 0.0
     doh2s(:, iu:id) = 0.0
     doch4(:, iu:id) = 0.0
     dofe2(:, iu:id) = 0.0
     domn2(:, iu:id) = 0.0                                                                  ! CEMA
 
     do i = iu, id
         DOSS(kt, i) = 0.0
         do k = kt, KB(i)
             do jcb = 1, nbod
                 if(BOD_CALC(jcb))dobod(k, i) = dobod(k, i) + RBOD(jcb)        &
                  & *cbodd(k, i, jcb)*CBOD(k, i, jcb)
             enddo
             do ja = 1, nal
                 if(ALG_CALC(ja))then
                     doap(k, i) = doap(k, i) + AGR(k, i, ja)*ALG(k, i, ja)     &
                                & *O2AG(ja)
                     doar(k, i) = doar(k, i) + ARR(k, i, ja)*ALG(k, i, ja)     &
                                & *O2AR(ja)
                 endif
             enddo
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))then
                     doep(k, i) = doep(k, i) + EGR(k, i, je)*EPC(k, i, je)     &
                                & *O2EG(je)
                     doer(k, i) = doer(k, i) + ERR(k, i, je)*EPC(k, i, je)     &
                                & *O2ER(je)
                 endif
             enddo
 
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     if(k==kt)then
                         jt = KTI(i)
                     else
                         jt = k
                     endif
                     je = KB(i)
                     do jj = jt, je
                         domp(k, i) = domp(k, i) + mgr(jj, k, i, m)            &
                                    & *MACRM(jj, k, i, m)*O2MG(m)
                         domr(k, i) = domr(k, i) + mrr(k, i, m)                &
                                    & *MACRM(jj, k, i, m)*O2MR(m)
                     enddo
                 endif
             enddo
             domp(k, i) = domp(k, i)/(DLX(i)*BH(k, i))
             domr(k, i) = domr(k, i)/(DLX(i)*BH(k, i))
             DOPOM(k, i) = (LPOMD(k, i) + RPOMD(k, i))*O2OM(jw)
             DODOM(k, i) = (LDOMD(k, i) + RDOMD(k, i))*O2OM(jw)
             DOOM(k, i) = DOPOM(k, i) + DODOM(k, i) + dobod(k, i)
             DONIT(k, i) = NH4D(k, i)*O2NH4(jw)
             if(sedcomp_exist)then ! SW 5/26/15
                 DOSED(k, i) = SEDD(k, i)*O2OM(jw) + SEDD1(k, i)*O2OM(jw)      &
                             & + SEDD2(k, i)*O2OM(jw)                                   !Amaila
             else
                 DOSED(k, i) = SEDD(k, i)*O2OM(jw)
             endif
             DOSOD(k, i) = SODD(k, i)*DO3(k, i)
!            CEMA start
             do jcg = ngcs, ngce
                 if(jcg==nch4)doch4(k, i) = ch4d(k, i)*o2ch4
                 if(jcg==nh2s)doh2s(k, i) = h2sd(k, i)*o2h2s
                 if(jcg==nfe2)dofe2(k, i) = fe2d(k, i)*o2fe2
                 if(jcg==nmn2)domn2(k, i) = mn2d(k, i)*o2mn2
             enddo
!            CEMA end
             if(zooplankton_calc)then
                 do jz = 1, nzp
                     dozr(k, i) = dozr(k, i) + ZRT(k, i, jz)*ZOO(k, i, jz)     &
                                & *O2ZR(jz)
                 enddo
             endif
    !DOSS(K,I)  =  DOAP(K,I)+DOEP(K,I)-DOAR(K,I)-DOER(K,I)-DOOM(K,I)-DONIT(K,I)-DOSOD(K,I)-DOSED(K,I)  &
    !                +DOMP(K,I)-DOMR(K,I)-DOZR(K,I)
             DOSS(k, i) = doap(k, i) + doep(k, i) - doar(k, i) - doer(k, i)    &
                        & - DOOM(k, i) - DONIT(k, i) - DOSOD(k, i)             &
                        & - DOSED(k, i) + domp(k, i) - domr(k, i) - dozr(k, i) &
                        & + doch4(k, i) + doh2s(k, i) + dofe2(k, i)            &
                        & + domn2(k, i)                                                                                                                                         !&     ! CEMA
                    !+DOMP(K,I)-DOMR(K,I)-DOZR(K,I)+doch4(k,i)+doh2s(k,i)   ! doch4, doh2s,dofe2 already negative...
         enddo
         dosat = SATO(T1(kt, i), TDS(kt, i), PALT(i), SALT_WATER(jw))
         if(.NOT.ICE(i))then
             call GAS_TRANSFER
             o2ex = REAER(i)
             DOAE(kt, i) = (dosat - O2(kt, i))*o2ex*BI(kt, i)/BH2(kt, i)
             DOSS(kt, i) = DOSS(kt, i) + DOAE(kt, i)
         endif
     enddo
     return
 
!***********************************************************************************************************************************
!**  I N O R G A N I C   C A R B O N                                          
!***********************************************************************************************************************************
 
!    **
     entry INORGANIC_CARBON
     ticap(:, iu:id) = 0.0
     ticep(:, iu:id) = 0.0
     ticbod(:, iu:id) = 0.0
     ticmc(:, iu:id) = 0.0
     ticzr(:, iu:id) = 0.0                  !v3.5
     ticch4 = 0.0
              ! CEMA
     do i = iu, id
         do k = kt, KB(i)
             do jcb = 1, nbod
                 if(BOD_CALC(jcb))ticbod(k, i) = ticbod(k, i)                  &
                  & + cbodd(k, i, jcb)*CBOD(k, i, jcb)*BODC(jcb)
             enddo
             do ja = 1, nal
                 if(ALG_CALC(ja))ticap(k, i) = ticap(k, i) + AC(ja)            &
                  & *(ARR(k, i, ja) - AGR(k, i, ja))*ALG(k, i, ja)
             enddo
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))ticep(k, i) = ticep(k, i) + EC(je)  &
                  & *(ERR(k, i, je) - EGR(k, i, je))*EPC(k, i, je)
             enddo
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     if(k==kt)then
                         jt = KTI(i)
                     else
                         jt = k
                     endif
                     je = KB(i)
                     do jj = jt, je
                         ticmc(k, i) = ticmc(k, i)                             &
                                     & + (mrr(k, i, m) - mgr(jj, k, i, m))     &
                                     & *MACRM(jj, k, i, m)*MC(m)
                     enddo
                 endif
             enddo
             ticmc(k, i) = ticmc(k, i)/(DLX(i)*BH(k, i))
             if(zooplankton_calc)then
                 do jz = 1, nzp
                     ticzr(k, i) = ticzr(k, i) + ZRT(k, i, jz)*ZOO(k, i, jz)   &
                                 & *ZC(jz)                   !MLM
                 enddo
      ! CEMA start
                 do jcg = ngcs, ngce
                     if(jcg==nch4)ticch4 = ch4d(k, i)
                 enddo
      ! CEMA end
             endif
 
             if(sedcomp_exist)then ! SW 5/26/15
               !TICSS(K,I) = TICAP(K,I)+TICEP(K,I)+SEDDC(K,I)+ORGC(JW)*(LPOMD(K,I)+RPOMD(K,I)+LDOMD(K,I)+RDOMD(K,I))                          &
               !    +CO2R(JW)*SODD(K,I)*DO3(K,I)+TICBOD(K,I)+TICMC(K,I)+TICZR(K,I)  + ticch4               &
               !    +sedd1(k,i)*orgc(jw) + sedd2(k,i)*orgc(jw)             ! Amaila
                 TICSS(k, i) = ticap(k, i) + ticep(k, i) + SEDDC(k, i)         &
                             & + ORGC(jw)                                      &
                             & *(LPOMD(k, i) + RPOMD(k, i) + LDOMD(k, i)       &
                             & + RDOMD(k, i)) + CO2R(jw)*SODD(k, i)*DO3(k, i)  &
                             & + ticbod(k, i) + ticmc(k, i) + ticzr(k, i)      &
                             & + ticch4 + SEDD1(k, i)*CBIOM(jw) + SEDD2(k, i)  &
                             & *CBIOM(jw)                                   ! Amaila, cb 6/7/17
 
             else
                 TICSS(k, i) = ticap(k, i) + ticep(k, i) + SEDDC(k, i)         &
                             & + ORGC(jw)                                      &
                             & *(LPOMD(k, i) + RPOMD(k, i) + LDOMD(k, i)       &
                             & + RDOMD(k, i)) + CO2R(jw)*SODD(k, i)*DO3(k, i)  &
                             & + ticbod(k, i) + ticmc(k, i) + ticzr(k, i)
             endif
 
         enddo
         if(.NOT.ICE(i))then
             if(REAER(i)==0.0)call GAS_TRANSFER
             co2ex = REAER(i)*0.923
             CO2REAER(kt, i) = co2ex*(0.286*EXP( - 0.0314*(T2(kt,i))*PALT(i))  &
                             & - CO2(kt, i))*BI(kt, i)/BH2(kt, i)
             TICSS(kt, i) = TICSS(kt, i) + CO2REAER(kt, i)
         endif
     enddo
     return
 
!***********************************************************************************************************************************
!**  S E D I M E N T                                                         
!***********************************************************************************************************************************
 
!    **
     entry SEDIMENT
     sedas(:, iu:id) = 0.0
     lpomep(:, iu:id) = 0.0
     sedcb(:, iu:id) = 0.0
     do i = iu, id
         sedsi = 0.0
         do k = kt, KB(i)
             if(k==KB(i))then
                 bibh2(k, i) = BI(k, i)/BH2(k, i)
             else
                 bibh2(k, i) = BI(k, i)/BH2(k, i)*(1.0 - BI(k + 1, i)/BI(k, i))
             endif
             do ja = 1, nal
                 if(ALG_CALC(ja))sedas(k, i) = sedas(k, i) + MAX(AS(ja), 0.0)  &
                  & *ALG(k, i, ja)*bibh2(k, i)                                                        !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
             enddo
             sedem = 0.0
                    ! CB 5/19/06
             do je = 1, nep
!                LPOMEP(K,I) = LPOMEP(K,I)+EPOM(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
                 if(EPIPHYTON_CALC(jw, je))sedem = sedem + EBR(k, i, je)       &
                  & /H1(k, i)*EPC(k, i, je)                                        ! cb 5/19/06
             enddo
             do jd = 1, nbod
                 if(BOD_CALC(jd))sedcb(k, i) = sedcb(k, i)                     &
                  & + MAX(CBODS(jd), 0.0)*CBOD(k, i, jd)*bibh2(k, i)/O2OM(jw)                                 !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
             enddo
             SEDOMS(k, i) = POMS(jw)*(LPOM(k, i) + RPOM(k, i))*bibh2(k, i)           !cb 10/22/06
             if(k==KB(i))then
                 sedso = 0.0
             else
                 sedso = SEDS(jw)*SED(k, i)*BI(k + 1, i)/BH2(k, i)             &
                       & *(1.0 - BI(k + 1, i)/BI(k, i))
             endif
             SEDNS(k, i) = sedsi - sedso
             sedsi = sedso
             if(k<KB(i))then
                          ! CEMA sediment in kb layer goes to sediment diagenesis model
                 SED(k, i) = MAX(SED(k, i) + (sedem + sedas(k, i) + sedcb(k, i)&
                           & + SEDOMS(k, i) + SEDNS(k, i) - SEDD(k, i)         &
                           & - SEDBR(k, i))*dlt, 0.0)                                                                   ! cb 11/30/06
             elseif(k==KB(i) .AND. .NOT.sediment_diagenesis)then
                 SED(k, i) = MAX(SED(k, i) + (sedem + sedas(k, i) + sedcb(k, i)&
                           & + SEDOMS(k, i) + SEDNS(k, i) - SEDD(k, i)         &
                           & - SEDBR(k, i))*dlt, 0.0)
             endif
         enddo
     enddo
     return
 
 
!***********************************************************************************************************************************
!**  S E D I M E N T   P H O S P H O R U S                                   
!***********************************************************************************************************************************
 
!    **
     entry SEDIMENTP
     sedasp(:, iu:id) = 0.0
     lpomepp(:, iu:id) = 0.0
     sedcbp(:, iu:id) = 0.0
     sdinp(:, iu:id) = 0.0
                       ! CEMA
     do i = iu, id
         sedsip = 0.0
         do k = kt, KB(i)
             sdalgp = 0.0
             do ja = 1, nal
                 if(ALG_CALC(ja))then
                     sedasp(k, i) = sedasp(k, i) + MAX(AS(ja), 0.0)*AP(ja)     &
                                  & *ALG(k, i, ja)*bibh2(k, i)                             !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
                     if(k==KB(i))sdalgp = sdalgp + MAX(AS(ja), 0.0)*AP(ja)     &
                      & *ALG(k, i, ja)*bibh2(k, i)                                  ! CEMA
                 endif
             enddo
             sdepp = 0.0
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))then
                     lpomepp(k, i) = lpomepp(k, i) + EPOM(je)*EP(je)           &
                                   & *(EMR(k, i, je)*EPC(k, i, je))
                     if(k==KB(i))sdepp = sdepp + EPOM(je)*EP(je)               &
                      & *(EMR(k, i, je)*EPC(k, i, je))                       ! CEMA
                 endif
             enddo
             sdbodp = 0.0
                   ! CEMA
             do jd = 1, nbod
! IF(BOD_CALC(JD))SEDCBP(K,I)=SEDCBP(K,I)+MAX(CBODS(JD),0.0)*BODP(JD)*CBOD(K,I,
                 if(BOD_CALC(jd))then
                     sedcbp(k, i) = sedcbp(k, i) + MAX(CBODS(jd), 0.0)         &
                                  & *CBODP(k, i, jd)*bibh2(k, i)                 ! CB 6/6/10
                     if(k==KB(i))sdbodp = sdbodp + MAX(CBODS(jd), 0.0)         &
                      & *CBODP(k, i, jd)*bibh2(k, i)                              ! CEMA
                 endif
             enddo
             SEDOMSP(k, i) = POMS(jw)*(LPOMP(k, i) + RPOMP(k, i))*bibh2(k, i)            !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  !CB 10/22/06
             if(k==KB(i))then
                 sedsop = 0.0
             else
                 sedsop = SEDS(jw)*SEDP(k, i)*BI(k + 1, i)/BH2(k, i)           &
                        & *(1.0 - BI(k + 1, i)/BI(k, i))
             endif
             SEDNSP(k, i) = sedsip - sedsop
             sedsip = sedsop
!            CEMA start
             SEDPINFLUX(k, i) = (lpomepp(k, i) + sedasp(k, i) + SEDOMSP(k, i)  &
                              & + sedcbp(k, i))*dlt
             if(k<KB(i))then
        !SEDP(K,I)    = MAX(SEDP(K,I)+(LPOMEPP(K,I)+SEDASP(K,I)+SEDOMSP(K,I)+SEDCBP(K,I)+SEDNSP(K,I)-SEDDP(K,I)   &   ! SW 4/8/16
                 SEDP(k, i) = MAX(SEDP(k, i) + SEDPINFLUX(k, i) + (SEDNSP(k, i)&
                            & - SEDDP(k, i) - SEDBRP(k, i))*dlt, 0.0)                                       !cb 11/30/06
             elseif(k==KB(i) .AND. .NOT.sediment_diagenesis)then
!                SEDP(K,I)    =
! MAX(SEDP(K,I)+(LPOMEPP(K,I)+SEDASP(K,I)+SEDOMSP(K,I)+SEDCBP(K,I)+SEDNSP(K,I)-
                 SEDP(k, i) = MAX(SEDP(k, i) + SEDPINFLUX(k, i) + (SEDNSP(k, i)&
                            & - SEDDP(k, i) - SEDBRP(k, i))*dlt, 0.0)                                       !cb 11/30/06
             elseif(k==KB(i) .AND. sediment_diagenesis)then
                 sdinp(k, i) = sdepp + sdalgp + sdbodp + SEDOMSP(k, i)         &
                             & + SEDNSP(k, i)                       ! CEMA calculating P flux to sediment diagnesis model
             endif
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  S E D I M E N T   N I T R O G E N                                       
!***********************************************************************************************************************************
 
!    **
     entry SEDIMENTN
     sedasn(:, iu:id) = 0.0
     lpomepn(:, iu:id) = 0.0
     sedcbn(:, iu:id) = 0.0
     sdinn(:, iu:id) = 0.0       ! CEMA
     do i = iu, id
         sedsin = 0.0
         do k = kt, KB(i)
             sdalgn = 0.0
             do ja = 1, nal
                 if(ALG_CALC(ja))then
                     sedasn(k, i) = sedasn(k, i) + MAX(AS(ja), 0.0)*AN(ja)     &
                                  & *ALG(k, i, ja)*bibh2(k, i)                               !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
                     if(k==KB(i))sdalgn = sdalgn + MAX(AS(ja), 0.0)*AN(ja)     &
                      & *ALG(k, i, ja)*bibh2(k, i)
                 endif
             enddo
             sdepn = 0.0
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))then
                     lpomepn(k, i) = lpomepn(k, i) + EPOM(je)*EN(je)           &
                                   & *(EMR(k, i, je)*EPC(k, i, je))
                     if(k==KB(i))sdepn = sdepn + EPOM(je)*EN(je)               &
                      & *(EMR(k, i, je)*EPC(k, i, je))
                 endif
             enddo
             sdbodn = 0.0
             do jd = 1, nbod
! IF(BOD_CALC(JD))SEDCBN(K,I)=SEDCBN(K,I)+MAX(CBODS(JD),0.0)*BODN(JD)*CBOD(K,I,
                 if(BOD_CALC(jd))then
                     sedcbn(k, i) = sedcbn(k, i) + MAX(CBODS(jd), 0.0)         &
                                  & *CBODN(k, i, jd)*bibh2(k, i)                 ! CB 6/6/10
                     if(k==KB(i))sdbodn = sdbodn + MAX(CBODS(jd), 0.0)         &
                      & *CBODN(k, i, jd)*bibh2(k, i)                              ! CEMA
                 endif
             enddo
             SEDOMSN(k, i) = POMS(jw)*(LPOMN(k, i) + RPOMN(k, i))*bibh2(k, i)              !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  !CB 10/22/06
 
             if(k==KB(i))then  ! SW 12/16/07
      !SEDNO3(K,I)  = FNO3SED(JW)*NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I))/BH2(K,I)      ! CEMA - KB layer N goes to sediment diagensis
                 sedson = 0.0
             else
                 SEDNO3(k, i) = FNO3SED(jw)*NO3(k, i)*NO3S(jw)*no3trm(k, i)    &
                              & *(BI(k, i) - BI(k + 1, i))/BH2(k, i)
                 sedson = SEDS(jw)*SEDN(k, i)*BI(k + 1, i)/BH2(k, i)           &
                        & *(1.0 - BI(k + 1, i)/BI(k, i))
             endif
             SEDNSN(k, i) = sedsin - sedson
             sedsin = sedson
!            CEMA start
             SEDNINFLUX(k, i) = (lpomepn(k, i) + sedasn(k, i) + SEDOMSN(k, i)  &
                              & + sedcbn(k, i))*dlt
             if(k<KB(i))then
!                SEDN(K,I)    =
! MAX(SEDN(K,I)+(LPOMEPN(K,I)+SEDASN(K,I)+SEDOMSN(K,I)+SEDCBN(K,I)+SEDNSN(K,I)+
                 SEDN(k, i) = MAX(SEDN(k, i) + SEDNINFLUX(k, i) + (SEDNSN(k, i)&
                            & + SEDNO3(k, i) - SEDDN(k, i) - SEDBRN(k, i))*dlt,&
                            & 0.0)                      !CB 11/30/06
             elseif(k==KB(i) .AND. .NOT.sediment_diagenesis)then
!                SEDN(K,I)    =
! MAX(SEDN(K,I)+(LPOMEPN(K,I)+SEDASN(K,I)+SEDOMSN(K,I)+SEDCBN(K,I)+SEDNSN(K,I)+
                 SEDN(k, i) = MAX(SEDN(k, i) + SEDNINFLUX(k, i) + (SEDNSN(k, i)&
                            & + SEDNO3(k, i) - SEDDN(k, i) - SEDBRN(k, i))*dlt,&
                            & 0.0)                      !CB 11/30/06
             elseif(k==KB(i) .AND. sediment_diagenesis)then
                 sdinn(k, i) = sdepn + sdalgn + sdbodn + +SEDOMSN(k, i)        &
                             & + SEDNSN(k, i)
             endif
!            CEMA end
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  S E D I M E N T   C A R B O N                                           
!***********************************************************************************************************************************
 
!    **
     entry SEDIMENTC
     sedasc(:, iu:id) = 0.0
     lpomepc(:, iu:id) = 0.0
     sedcbc(:, iu:id) = 0.0
     sdinc(:, iu:id) = 0.0
                       ! CEMA
     do i = iu, id
         do k = kt, KB(i)
             sedsip = 0.0
             sdalgc = 0.0
             do ja = 1, nal
                 if(ALG_CALC(ja))then
                     sedasc(k, i) = sedasc(k, i) + MAX(AS(ja), 0.0)*AC(ja)     &
                                  & *ALG(k, i, ja)*bibh2(k, i)                                !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
                     if(k==KB(i))sdalgc = sdalgc + MAX(AS(ja), 0.0)*AC(ja)     &
                      & *ALG(k, i, ja)*bibh2(k, i)                                  ! CEMA
                 endif
             enddo
             sdepc = 0.0
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))then
                     lpomepc(k, i) = lpomepc(k, i) + EPOM(je)*EC(je)           &
                                   & *(EMR(k, i, je)*EPC(k, i, je))
                     if(k==KB(i))sdepc = sdepc + EPOM(je)*EC(je)               &
                      & *(EMR(k, i, je)*EPC(k, i, je))                       ! CEMA
                 endif
             enddo
             sdbodc = 0.0
                    ! CEMA
             do jd = 1, nbod
                 if(BOD_CALC(jd))then
                     sedcbc(k, i) = sedcbc(k, i) + MAX(CBODS(jd), 0.0)*BODC(jd)&
                                  & *CBOD(k, i, jd)*bibh2(k, i)                               !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
                     if(k==KB(i))sdbodc = sdbodc + MAX(CBODS(jd), 0.0)*BODC(jd)&
                      & *CBOD(k, i, jd)*bibh2(k, i)                                       ! CEMA
                 endif
             enddo
             SEDOMSC(k, i) = POMS(jw)*ORGC(jw)*(LPOM(k, i) + RPOM(k, i))       &
                           & *bibh2(k, i)                                                   !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))   !CB 10/22/06
             if(k==KB(i))then
                 sedsoc = 0.0
             else
                 sedsoc = SEDS(jw)*SEDC(k, i)*BI(k + 1, i)/BH2(k, i)           &
                        & *(1.0 - BI(k + 1, i)/BI(k, i))
             endif
             SEDNSC(k, i) = sedsic - sedsoc
             sedsic = sedsoc
!            CEMA start
             if(k<KB(i))then
                 SEDC(k, i) = MAX(SEDC(k, i) + (lpomepc(k, i) + sedasc(k, i) + &
                            & SEDOMSC(k, i) + sedcbc(k, i) + SEDNSC(k, i)      &
                            & - SEDDC(k, i) - SEDBRC(k, i))*dlt, 0.0)
             elseif(k==KB(i) .AND. .NOT.sediment_diagenesis)then
                 SEDC(k, i) = MAX(SEDC(k, i) + (lpomepc(k, i) + sedasc(k, i) + &
                            & SEDOMSC(k, i) + sedcbc(k, i) + SEDNSC(k, i)      &
                            & - SEDDC(k, i) - SEDBRC(k, i))*dlt, 0.0)
             elseif(k==KB(i) .AND. sediment_diagenesis)then
                 sdinc(k, i) = sdepc + sdalgc + sdbodc + SEDOMSC(k, i)         &
                             & + SEDNSC(k, i)                       ! CEMA calculating C flux to sediment diagnesis model
             endif
!            CEMA end
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  S E D I M E N T   D E C A Y    R A T E                                  
!***********************************************************************************************************************************
 
!    **
     entry SEDIMENT_DECAY_RATE
     do i = iu, id
         sedsidk = 0.0
         do k = kt, KB(i)
             sedsum = 0.0
             sedsumk = 0.0
 
             do ja = 1, nal
                 if(ALG_CALC(ja))then
                     xdum = MAX(AS(ja), 0.0)*ALG(k, i, ja)*bibh2(k, i)
                     sedsumk = sedsumk + xdum*LPOMDK(jw)
                     sedsum = sedsum + xdum
                 endif
             enddo
 
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))then
                     xdum = EPOM(je)*(EMR(k, i, je)*EPC(k, i, je))
                     sedsumk = sedsumk + xdum*LPOMDK(jw)
                     sedsum = sedsum + xdum
                 endif
             enddo
 
             do jd = 1, nbod
                 if(BOD_CALC(jd))then
                     xdum = MAX(CBODS(jd), 0.0)*CBOD(k, i, jd)*bibh2(k, i)     &
                          & *RBOD(jd)/O2OM(jw)
                     sedsumk = sedsumk + xdum*cbodd(k, i, jd)
                     sedsum = sedsum + xdum
                 endif
             enddo
 
             sedsumk = sedsumk + POMS(jw)                                      &
                     & *(LPOM(k, i)*LPOMDK(jw) + RPOM(k, i)*RPOMDK(jw))        &
                     & *bibh2(k, i)                                                              !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  ! CB 10/22/06
             sedsum = sedsum + POMS(jw)*(LPOM(k, i) + RPOM(k, i))*bibh2(k, i)
 
             sedsumk = sedsumk*dlt
             sedsum = sedsum*dlt
 
             if((sedsum + SED(k, i))>0.0)then
                 SDKV(k, i) = (sedsumk + SED(k, i)*SDKV(k, i))                 &
                            & /(sedsum + SED(k, i))
             else
                 SDKV(k, i) = 0.0
             endif
 
         enddo
     enddo
     return
 
!    Amaila start
!    additional sediment compartments simulate slow and fast decaying OM left
!***********************************************************************************************************************************
!**  in standing trees S E D I M E N T  1                                     
!***********************************************************************************************************************************
 
!    **
     entry SEDIMENT1
     do i = iu, id
         do k = kt, KB(i)
             SED1(k, i) = MAX(SED1(k, i) + ( - SEDD1(k, i))*dlt, 0.0)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  S E D I M E N T  2                                                      
!***********************************************************************************************************************************
 
!    **
     entry SEDIMENT2
     do i = iu, id
         do k = kt, KB(i)
             SED2(k, i) = MAX(SED2(k, i) + ( - SEDD2(k, i))*dlt, 0.0)
         enddo
     enddo
     return
 
!    Amaila end
 
!***********************************************************************************************************************************
!*   E P I P H Y T O N                                                      **
!***********************************************************************************************************************************
 
     entry EPIPHYTON(J)
     do i = iu, id
 
!**      Limiting factor
 
         light = (1.0 - BETA(jw))*SRON(jw)*SHADE(i)/ESAT(J)
         lam2 = light
         lam1 = light
         do k = kt, KB(i)
 
!****        Limiting factor
 
             lam1 = lam2
             lam2 = lam1*EXP( - GAMMA(k, i)*H1(k, i))
             fdpo4 = 1.0 - FPSS(k, i) - FPFE(k, i)
             ELLIM(k, i, J) = 2.718282*(EXP( - lam2) - EXP( - lam1))           &
                            & /(GAMMA(k, i)*H1(k, i))
             if(EHSP(J)/=0.0)EPLIM(k, i, J) = fdpo4*PO4(k, i)                  &
              & /(fdpo4*PO4(k, i) + EHSP(J) + nonzero)
             if(EHSN(J)/=0.0)ENLIM(k, i, J) = (NH4(k, i) + NO3(k, i))          &
              & /(NH4(k, i) + NO3(k, i) + EHSN(J) + nonzero)
             if(EHSSI(J)/=0.0)ESLIM(k, i, J) = DSI(k, i)                       &
              & /(DSI(k, i) + EHSSI(J) + nonzero)
             limit = MIN(EPLIM(k, i, J), ENLIM(k, i, J), ESLIM(k, i, J),       &
                   & ELLIM(k, i, J))
             blim = 1.0 - (EPD(k, i, J)/(EPD(k, i, J) + EHS(J)))
 
!****        Sources/sinks
 
             EGR(k, i, J) = MIN(etrm(k, i, J)*EG(J)*limit*blim, PO4(k, i)      &
                          & /(EP(J)*dlt*EPD(k, i, J)/H1(kt, i) + nonzero),     &
                          & (NH4(k, i) + NO3(k, i))                            &
                          & /(EN(J)*dlt*EPD(k, i, J)/H1(k, i) + nonzero))
             ERR(k, i, J) = etrm(k, i, J)*ER(J)*DO3(k, i)
             EMR(k, i, J) = (etrmr(k, i, J) + 1.0 - etrmf(k, i, J))*EM(J)
             EER(k, i, J) = MIN((1.0 - ELLIM(k, i, J))*EE(J)*etrm(k, i, J),    &
                          & EGR(k, i, J))
!            EPD(K,I,J) = 
! MAX(EPD(K,I,J)+EPD(K,I,J)*(EGR(K,I,J)-ERR(K,I,J)-EMR(K,I,J)-EER(K,I,J)-EBR(K,
             EPD(k, i, J) = MAX(EPD(k, i, J) + EPD(k, i, J)*(EGR(k, i, J) - ERR&
                          & (k, i, J) - EMR(k, i, J) - EER(k, i, J)            &
                          & - EBR(k, i, J)/H1(k, i))*dlt, 0.00)                                                            ! cb 5/18/06
             if(k==KB(i))then  ! SW 12/16/07
                 EPM(k, i, J) = EPD(k, i, J)*(BI(k, i) + 2.0*H1(k, i))*DLX(i)
             else
                 EPM(k, i, J) = EPD(k, i, J)                                   &
                              & *(BI(k, i) - BI(k + 1, i) + 2.0*H1(k, i))      &
                              & *DLX(i)
             endif
             EPC(k, i, J) = EPM(k, i, J)/VOL(k, i)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  L A B I L E   D O M   P H O S P H O R U S                               **
!***********************************************************************************************************************************
 
     entry LABILE_DOM_P
     ldompap(:, iu:id) = 0.0
     ldompep(:, iu:id) = 0.0
     ldompmp(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             do ja = 1, nal
                 if(ALG_CALC(ja))ldompap(k, i) = ldompap(k, i)                 &
                  & + (AER(k, i, ja) + (1.0 - APOM(ja))*AMR(k, i, ja))         &
                  & *ALG(k, i, ja)*AP(ja)
             enddo
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))ldompep(k, i) = ldompep(k, i)       &
                  & + (EER(k, i, je) + (1.0 - EPOM(je))*EMR(k, i, je))         &
                  & *EPC(k, i, je)*EP(je)
             enddo
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     if(k==kt)then
                         jt = KTI(i)
                     else
                         jt = k
                     endif
                     je = KB(i)
                     do jj = jt, je
                         ldompmp(k, i) = ldompmp(k, i) + (1.0 - MPOM(m))       &
                           & *mmr(k, i, m)*MACRM(jj, k, i, m)*MP(m)
                     enddo
                 endif
             enddo
             ldompmp(k, i) = ldompmp(k, i)/(DLX(i)*BH(k, i))
             LDOMPSS(k, i) = ldompap(k, i) + ldompep(k, i) + ldompmp(k, i)     &
                           & - (LDOMD(k, i) + LRDOMD(k, i))*ORGPLD(k, i)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  R E F R A C T O R Y   D O M   P H O S P H O R U S                        
!***********************************************************************************************************************************
 
!    **
     entry REFRACTORY_DOM_P
     do i = iu, id
         do k = kt, KB(i)
             RDOMPSS(k, i) = LRDOMD(k, i)*ORGPLD(k, i) - RDOMD(k, i)           &
                           & *ORGPRD(k, i)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  L A B I L E   P O M   P H O S P H O R U S                               
!***********************************************************************************************************************************
 
!    **
     entry LABILE_POM_P
     lpompap(:, iu:id) = 0.0
     lpompmp(:, iu:id) = 0.0
     lpzooinp(:, iu:id) = 0.0
     lpzoooutp(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             do ja = 1, nal
                 if(ALG_CALC(ja))lpompap(k, i) = lpompap(k, i) + APOM(ja)      &
                  & *(AMR(k, i, ja)*ALG(k, i, ja))*AP(ja)
             enddo
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     jt = k
                     je = KB(i)
                     do jj = jt, je
                         lpompmp(k, i) = lpompmp(k, i) + MPOM(m)*LRPMAC(m)     &
                           & *mmr(k, i, m)*MACRM(jj, k, i, m)*MP(m)
                     enddo
                 endif
             enddo
             lpompmp(k, i) = lpompmp(k, i)/(DLX(i)*BH(k, i))
             if(zooplankton_calc)then
                 do jz = 1, nzp
                     if(TGRAZE(k, i, jz)>0.0)then
                         lpzoooutp(k, i) = lpzoooutp(k, i) + ZOO(k, i, jz)     &
                           & *(ZMT(k, i, jz)                                   &
                           & + (ZMU(k, i, jz) - (ZMU(k,i,jz)*ZEFF(jz))))*ZP(jz)
                         lpzooinp(k, i) = lpzooinp(k, i) + ZOO(k, i, jz)       &
                           & *ZMU(k, i, jz)*PREFP(jz)*LPOM(k, i)               &
                           & /TGRAZE(k, i, jz)*ZP(jz)
                     else
                         lpzoooutp(k, i) = lpzoooutp(k, i) + ZOO(k, i, jz)     &
                           & *(ZMT(k, i, jz)                                   &
                           & + (ZMU(k, i, jz) - (ZMU(k,i,jz)*ZEFF(jz))))*ZP(jz)
                         lpzooinp(k, i) = 0.0
                     endif
                 enddo
             endif
             LPOMPNS(k, i) = POMS(jw)                                          &
                           & *(LPOM(k - 1, i)*ORGPLP(k - 1, i) - LPOM(k, i)    &
                           & *ORGPLP(k, i))*BI(k, i)/BH2(k, i)
             LPOMPSS(k, i) = lpompap(k, i) + lpompmp(k, i) - LPOMD(k, i)       &
                           & *ORGPLP(k, i) + LPOMPNS(k, i) - LRPOMD(k, i)      &
                           & *ORGPLP(k, i)
!!   DO JZ = 1,NZP                                           ! KV 4/24/12
!!   END DO                                                  ! KV 4/24/12
             if(zooplankton_calc)LPOMPSS(k, i) = LPOMPSS(k, i)                 &
              & + lpzoooutp(k, i) - lpzooinp(k, i)
 
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  R E F R A C T O R Y   P O M   P H O S P H O R U S                        
!***********************************************************************************************************************************
 
!    **
     entry REFRACTORY_POM_P
     rpompmp(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     jt = k
                     je = KB(i)
                     do jj = jt, je
                         rpompmp(k, i) = rpompmp(k, i) + MPOM(m)               &
                           & *(1.0 - LRPMAC(m))*mmr(k, i, m)*MACRM(jj, k, i, m)&
                           & *MP(m)
                     enddo
                 endif
             enddo
             rpompmp(k, i) = rpompmp(k, i)/(DLX(i)*BH(k, i))
             RPOMPNS(k, i) = POMS(jw)                                          &
                           & *(RPOM(k - 1, i)*ORGPRP(k - 1, i) - RPOM(k, i)    &
                           & *ORGPRP(k, i))*BI(k, i)/BH2(k, i)
             RPOMPSS(k, i) = LRPOMD(k, i)*ORGPLP(k, i) + RPOMPNS(k, i)         &
                           & - RPOMD(k, i)*ORGPRP(k, i) + rpompmp(k, i)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  L A B I L E   D O M   N I T R O G E N                                   **
!***********************************************************************************************************************************
 
     entry LABILE_DOM_N
     ldomnap(:, iu:id) = 0.0
     ldomnep(:, iu:id) = 0.0
     ldomnmp(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             do ja = 1, nal
                 if(ALG_CALC(ja))ldomnap(k, i) = ldomnap(k, i)                 &
                  & + (AER(k, i, ja) + (1.0 - APOM(ja))*AMR(k, i, ja))         &
                  & *ALG(k, i, ja)*AN(ja)
             enddo
             do je = 1, nep
                 if(EPIPHYTON_CALC(jw, je))ldomnep(k, i) = ldomnep(k, i)       &
                  & + (EER(k, i, je) + (1.0 - EPOM(je))*EMR(k, i, je))         &
                  & *EPC(k, i, je)*EN(je)
             enddo
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     if(k==kt)then
                         jt = KTI(i)
                     else
                         jt = k
                     endif
                     je = KB(i)
                     do jj = jt, je
                         ldomnmp(k, i) = ldomnmp(k, i) + (1.0 - MPOM(m))       &
                           & *mmr(k, i, m)*MACRM(jj, k, i, m)*MN(m)
                     enddo
                 endif
             enddo
             ldomnmp(k, i) = ldomnmp(k, i)/(DLX(i)*BH(k, i))
             LDOMNSS(k, i) = ldomnap(k, i) + ldomnep(k, i) + ldomnmp(k, i)     &
                           & - (LDOMD(k, i) + LRDOMD(k, i))*ORGNLD(k, i)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  R E F R A C T O R Y   D O M   N I T R O G E N                            
!***********************************************************************************************************************************
 
!    **
     entry REFRACTORY_DOM_N
     do i = iu, id
         do k = kt, KB(i)
             RDOMNSS(k, i) = LRDOMD(k, i)*ORGNLD(k, i) - RDOMD(k, i)           &
                           & *ORGNRD(k, i)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  L A B I L E   P O M   N I T R O G E N                                   
!***********************************************************************************************************************************
 
!    **
     entry LABILE_POM_N
     lpomnap(:, iu:id) = 0.0
     lpomnmp(:, iu:id) = 0.0
     lpzooinn(:, iu:id) = 0.0
     lpzoooutn(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             do ja = 1, nal
                 if(ALG_CALC(ja))lpomnap(k, i) = lpomnap(k, i) + APOM(ja)      &
                  & *(AMR(k, i, ja)*ALG(k, i, ja))*AN(ja)
             enddo
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     jt = k
                     je = KB(i)
                     do jj = jt, je
                         lpomnmp(k, i) = lpomnmp(k, i) + MPOM(m)*LRPMAC(m)     &
                           & *mmr(k, i, m)*MACRM(jj, k, i, m)*MN(m)
                     enddo
                 endif
             enddo
             lpomnmp(k, i) = lpomnmp(k, i)/(DLX(i)*BH(k, i))
             if(zooplankton_calc)then
                 do jz = 1, nzp
                     if(TGRAZE(k, i, jz)>0.0)then
                         lpzoooutn(k, i) = lpzoooutn(k, i) + ZOO(k, i, jz)     &
                           & *(ZMT(k, i, jz)                                   &
                           & + (ZMU(k, i, jz) - (ZMU(k,i,jz)*ZEFF(jz))))*ZN(jz)
                         lpzooinn(k, i) = lpzooinn(k, i) + ZOO(k, i, jz)       &
                           & *PREFP(jz)*ZMU(k, i, jz)*LPOM(k, i)               &
                           & /TGRAZE(k, i, jz)*ZN(jz)
                     else
                         lpzoooutn(k, i) = lpzoooutn(k, i) + ZOO(k, i, jz)     &
                           & *(ZMT(k, i, jz)                                   &
                           & + (ZMU(k, i, jz) - (ZMU(k,i,jz)*ZEFF(jz))))*ZN(jz)
                         lpzooinn(k, i) = 0.0
                     endif
                 enddo
             endif
             LPOMNNS(k, i) = POMS(jw)                                          &
                           & *(LPOM(k - 1, i)*ORGNLP(k - 1, i) - LPOM(k, i)    &
                           & *ORGNLP(k, i))*BI(k, i)/BH2(k, i)
             LPOMNSS(k, i) = lpomnap(k, i) + lpomnmp(k, i) - LPOMD(k, i)       &
                           & *ORGNLP(k, i) + LPOMNNS(k, i) - LRPOMD(k, i)      &
                           & *ORGNLP(k, i) + lpzoooutn(k, i) - lpzooinn(k, i)
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  R E F R A C T O R Y   P O M   N I T R O G E N                            
!***********************************************************************************************************************************
 
!    **
     entry REFRACTORY_POM_N
     rpomnmp(:, iu:id) = 0.0
     do i = iu, id
         do k = kt, KB(i)
             do m = 1, nmc
                 if(MACROPHYTE_CALC(jw, m))then
                     jt = k
                     je = KB(i)
                     do jj = jt, je
                         rpomnmp(k, i) = rpomnmp(k, i) + MPOM(m)               &
                           & *(1.0 - LRPMAC(m))*mmr(k, i, m)*MACRM(jj, k, i, m)&
                           & *MN(m)
                     enddo
                 endif
             enddo
             rpomnmp(k, i) = rpomnmp(k, i)/(DLX(i)*BH(k, i))
             RPOMNNS(k, i) = POMS(jw)                                          &
                           & *(RPOM(k - 1, i)*ORGNRP(k - 1, i) - RPOM(k, i)    &
                           & *ORGNRP(k, i))*BI(k, i)/BH2(k, i)
             RPOMNSS(k, i) = LRPOMD(k, i)*ORGNLP(k, i) + RPOMNNS(k, i)         &
                           & - RPOMD(k, i)*ORGNRP(k, i) + rpomnmp(k, i)
         enddo
     enddo
     return
 
 
!************************************************************************
!**  M A C R O P H Y T E                       **
!************************************************************************
 
     entry MACROPHYTE(Llm)
     m = Llm
     do i = iu, id
         if(KTICOL(i))then
             jt = KTI(i)
         else
             jt = KTI(i) + 1
         endif
         je = KB(i)
         do jj = jt, je
             if(jj<kt)then
                 colb = EL(jj + 1, i)
             else
                 colb = EL(kt + 1, i)
             endif
      !COLDEP=ELWS(I)-COLB
             coldep = EL(kt, i) - Z(i)*COSA(jb) - colb
                                           ! cb 3/7/16
             if(MACRC(jj, kt, i, m)>MMAX(m))mgr(jj, kt, i, m) = 0.0
             MACSS(jj, kt, i, m) = (mgr(jj, kt, i, m) - mmr(kt, i, m) - mrr(kt,&
                                 & i, m))*MACRC(jj, kt, i, m)
             MACRM(jj, kt, i, m) = MACRM(jj, kt, i, m) + MACSS(jj, kt, i, m)   &
                                 & *dlt*coldep*CW(jj, i)*DLX(i)
         enddo
 
         do k = kt + 1, KB(i)
             jt = k
             je = KB(i)
             do jj = jt, je
                 if(MACRC(jj, k, i, m)>MMAX(m))mgr(jj, k, i, m) = 0.0
                 MACSS(jj, k, i, m) = (mgr(jj, k, i, m) - mmr(k, i, m) - mrr(k,&
                                    & i, m))*MACRC(jj, k, i, m)
                 if(MACT(jj, k, i)>MBMP(m) .AND. MACT(jj, k - 1, i)<MBMP(m)    &
                  & .AND. MACSS(jj, k, i, m)>0.0)then
                     if(k - 1==kt)then
                         bmass = MACSS(jj, k, i, m)*dlt*H2(k, i)*CW(jj, i)     &
                               & *DLX(i)
                         MACRM(jj, k - 1, i, m) = MACRM(jj, k - 1, i, m)       &
                           & + bmass
                         colb = EL(kt + 1, i)
            !COLDEP=ELWS(I)-COLB
                         coldep = EL(kt, i) - Z(i)*COSA(jb) - colb
                                                 ! cb 3/7/16
                         MACSS(jj, k - 1, i, m)                                &
                           & = bmass/dlt/(coldep*CW(jj, i)*DLX(i))             &
                           & + MACSS(jj, k - 1, i, m)
                     else
                         bmass = MACSS(jj, k, i, m)*dlt*H2(k, i)*CW(jj, i)     &
                               & *DLX(i)
                         MACRM(jj, k - 1, i, m) = MACRM(jj, k - 1, i, m)       &
                           & + bmass
                         MACSS(jj, k - 1, i, m)                                &
                           & = bmass/dlt/(H2(k - 1, i)*CW(jj, i)*DLX(i))       &
                           & + MACSS(jj, k - 1, i, m)
                     endif
                     MACSS(jj, k, i, m) = 0.0
                 else
                     bmasstest = MACRM(jj, k, i, m) + MACSS(jj, k, i, m)       &
                               & *dlt*H2(k, i)*CW(jj, i)*DLX(i)
                     if(bmasstest>=0.0)then
                         MACRM(jj, k, i, m) = bmasstest
                     else
                         MACSS(jj, k, i, m) = -MACRM(jj, k, i, m)              &
                           & /dlt/(H2(k, i)*CW(jj, i)*DLX(i))
                         MACRM(jj, k, i, m) = 0.0
                     endif
                 endif
             enddo
         enddo
     enddo
     do i = iu, id
         tmac = 0.0
         cvol = 0.0
         if(KTICOL(i))then
             jt = KTI(i)
         else
             jt = KTI(i) + 1
         endif
         je = KB(i)
 
         do jj = jt, je
             if(jj<kt)then
                 colb = EL(jj + 1, i)
             else
                 colb = EL(kt + 1, i)
             endif
      !COLDEP=ELWS(I)-COLB
             coldep = EL(kt, i) - Z(i)*COSA(jb) - colb
                                           ! cb 3/7/16
             if(CW(jj, i)>0.0)then
                 MACRC(jj, kt, i, m) = MACRM(jj, kt, i, m)                     &
                                     & /(CW(jj, i)*coldep*DLX(i))
             else
                 MACRC(jj, kt, i, m) = 0.0
             endif
             tmac = tmac + MACRM(jj, kt, i, m)
             cvol = cvol + CW(jj, i)*coldep*DLX(i)
         enddo
 
         MAC(kt, i, m) = tmac/cvol
 
         do k = kt + 1, KB(i)
             jt = k
             je = KB(i)
             tmac = 0.0
             cvol = 0.0
             do jj = jt, je
                 if(CW(jj, i)>0.0)then
                     MACRC(jj, k, i, m) = MACRM(jj, k, i, m)                   &
                       & /(CW(jj, i)*H2(k, i)*DLX(i))
                 else
                     MACRC(jj, k, i, m) = 0.0
                 endif
                 tmac = tmac + MACRM(jj, k, i, m)
                 cvol = cvol + CW(jj, i)*H2(k, i)*DLX(i)
             enddo
             MAC(k, i, m) = tmac/cvol
         enddo
     enddo
 
     do i = iu, id
         tmac = 0.0
         cvol = 0.0
         do k = kt, KB(i)
             if(k==kt)then
                 jt = KTI(i)
             else
                 jt = k
             endif
             je = KB(i)
             do jj = jt, je
                 MACT(jj, k, i) = 0.0
                 do mi = 1, nmc
                     if(MACROPHYTE_CALC(jw, mi))MACT(jj, k, i)                 &
                      & = MACRC(jj, k, i, mi) + MACT(jj, k, i)
                 enddo
             enddo
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!*   K I N E T I C   F L U X E S                                              
!***********************************************************************************************************************************
 
!    **
     entry KINETIC_FLUXES
     do jaf = 1, NAF(jw)
         do jb = BS(jw), BE(jw)        ! SW 3/9/16
             do i = CUS(jb), DS(jb)
                 do k = kt, KB(i)
                     KFS(k, i, KFCN(jaf, jw)) = KFS(k, i, KFCN(jaf, jw))       &
                       & + KF(k, i, KFCN(jaf, jw))*VOL(k, i)*dlt                             ! KF IN G/M3/S x VOL M3 x DT S == G
                 enddo
             enddo
         enddo
     enddo
 
     return
 
!***********************************************************************************************************************************
!**  p H   C O 2                                                             **
!***********************************************************************************************************************************
 
     entry PH_CO2
 
!    pH and carbonate species
 
     do i = iu, id
         do k = kt, KB(i)
             cart = TIC(k, i)/12000.0        ! CART=equivalents/liter of C    TIC=mg/l C (MW=12g/mole)
             alkt = ALK(k, i)/5.0E+04        ! ALK=mg/l as CaCO3 (MW=50 g/mole; EQ=50g/eq))      ALKT=equivalents/l
             t1k = T1(k, i) + 273.15
 
!****        Ionic strength
 
             if(FRESH_WATER(jw))s2 = 2.5E-05*TDS(k, i)
             if(SALT_WATER(jw))s2 = 1.47E-3 + 1.9885E-2*TDS(k, i)              &
                                  & + 3.8E-5*TDS(k, i)*TDS(k, i)
 
!****        Debye-Huckel terms and activity coefficients
 
             sqrs2 = SQRT(s2)
             dh1 = -0.5085*sqrs2/(1.0 + 1.3124*sqrs2) + 4.745694E-03 +         &
                 & 4.160762E-02*s2 - 9.284843E-03*s2*s2
             dh2 = -2.0340*sqrs2/(1.0 + 1.4765*sqrs2) + 1.205665E-02 +         &
                 & 9.715745E-02*s2 - 2.067746E-02*s2*s2
             h2co3t = 10.0**(0.0755*s2)
             hco3t = 10.0**dh1
             co3t = 10.0**dh2
             oh = hco3t
 
!****        Temperature adjustment
 
             kw = 10.0**( - 283.971 - 0.05069842*t1k + 13323.0/t1k +           &
                & 102.24447*LOG10(t1k) - 1119669.0/(t1k*t1k))/oh
             k1 = 10.0**( - 3404.71/t1k + 14.8435 - 0.032786*t1k)*h2co3t/hco3t
             k2 = 10.0**( - 2902.39/t1k + 6.4980 - 0.023790*t1k)*hco3t/co3t
 
!****        pH evaluation
 
             pht = -PH(k, i) - 2.1
             if(PH(k, i)<=0.0)pht = -14.0
             incr = 10.0
             do n = 1, 3
                 f = 1.0
                 incr = incr/10.0
                 iter = 0
                 do while (f>0.0 .AND. iter<12)
                     pht = pht + incr
                     hion = 10.0**pht
                     bicart = cart*k1*hion/(k1*hion + k1*k2 + hion*hion)
                     f = bicart*(hion + 2.0*k2)/hion + kw/hion - alkt - hion/oh
                     iter = iter + 1
                 enddo
                 pht = pht - incr
             enddo
 
!****        pH, carbon dioxide, bicarbonate, and carbonate concentrations
 
             hion = 10.0**pht
             PH(k, i) = -pht
             CO2(k, i) = TIC(k, i)/(1.0 + k1/hion + k1*k2/(hion*hion))
                                                                     ! mg/l as C
             HCO3(k, i) = TIC(k, i)/(1.0 + hion/k1 + k2/hion)        ! mg/l as C
             CO3(k, i) = TIC(k, i)/((hion*hion)/(k1*k2) + hion/k2 + 1.0)
                                                                     ! mg/l as C
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  p H   C O 2   N E W                                                     **
!***********************************************************************************************************************************
 
     entry PH_CO2_NEW
                 ! Enhancements added for buffering by ammonia, phosphate, and OM ! SR 01/01/12
!    pH and carbonate species
     do i = iu, id
         do k = kt, KB(i)
             t1k = T1(k, i) + 273.15
             cart = TIC(k, i)/12011.
                            ! SR 01/01/12
             alkt = ALK(k, i)/50044.
                            ! SR 01/01/12
             ammt = NH4(k, i)/14006.74
                              ! SR 01/01/12
             phost = PO4(k, i)/30973.762
                                ! SR 01/01/12
             omct = (LDOM(k, i) + RDOM(k, i))*ORGC(jw)/12011.
                                                  ! moles carbon per liter from DOM ! SR 01/01/12
             if(pom_buffering)omct = omct + (LPOM(k, i) + RPOM(k, i))*ORGC(jw) &
                                   & /12011.                                ! SR 01/01/12
 !**** Ionic strength
             if(FRESH_WATER(jw))s2 = 2.5E-05*TDS(k, i)
             if(SALT_WATER(jw))s2 = 1.47E-3 + 1.9885E-2*TDS(k, i)              &
                                  & + 3.8E-5*TDS(k, i)*TDS(k, i)
!****        Debye-Huckel terms and activity coefficients
             sqrs2 = SQRT(s2)
             dh1 = -0.5085*sqrs2/(1.0 + 1.3124*sqrs2) + 4.745694E-03 +         &
                 & 4.160762E-02*s2 - 9.284843E-03*s2*s2
             dh2 = -2.0340*sqrs2/(1.0 + 1.4765*sqrs2) + 1.205665E-02 +         &
                 & 9.715745E-02*s2 - 2.067746E-02*s2*s2
             dh3 = -4.5765*sqrs2/(1.0 + 1.3124*sqrs2)
                                            ! extended Debye-Huckel for PO4 ! SR 01/01/12
             dhh = -0.5085*sqrs2/(1.0 + 2.9529*sqrs2)
                                            ! extended Debye-Huckel for H+ ion ! SR 01/01/12
             h2co3t = 10.0**(0.0755*s2)
             hco3t = 10.0**dh1
             co3t = 10.0**dh2
             po4t = 10.0**dh3
                      ! SR 01/01/12
             ht = 10.0**dhh
                    ! activity coefficient for H+ ! SR 01/01/12
             hpo4t = co3t
                  ! tabled values similar to those for carbonate ! SR 01/01/12
             oht = hco3t
                 ! tabled values similar to those for bicarbonate ! SR 01/01/12
             h2po4t = hco3t
                    ! tabled values similar to those for bicarbonate ! SR 01/01/12
             nh4t = hco3t
                  ! tabled values similar to those for bicarbonate ! SR 01/01/12
             nh3t = h2co3t
                   ! neutral species, set coefficient to same as that for carbonic acid ! SR 01/01/12
             h3po4t = h2co3t
                     ! neutral species, set coefficient to same as that for carbonic acid ! SR 01/01/12
!****        Temperature adjustment
             kw = 10.0**( - 283.971 - 0.05069842*t1k + 13323.0/t1k +           &
                & 102.24447*LOG10(t1k) - 1119669.0/(t1k*t1k))/oht
             k1 = 10.0**( - 356.3094 - 0.06091964*t1k + 21834.37/t1k +         &
                & 126.8339*LOG10(t1k) - 1684915/(t1k*t1k))*h2co3t/hco3t
             k2 = 10.0**( - 107.8871 - 0.03252849*t1k + 5151.79/t1k +          &
                & 38.92561*LOG10(t1k) - 563713.9/(t1k*t1k))*hco3t/co3t
             kamm = 10.0**( - 0.09018 - 2729.92/t1k)*nh4t/nh3t
                                                    ! SR 01/01/12
             kp1 = 10.0**(4.5535 - 0.013486*t1k - 799.31/t1k)*h3po4t/h2po4t
                                                                  ! Bates (1951) ! SR 01/21/12
             kp2 = 10.0**(5.3541 - 0.019840*t1k - 1979.5/t1k)*h2po4t/hpo4t
                                                                 ! Bates and Acree (1943) ! SR 01/21/12
             kp3 = 10.0**( - 12.38)*hpo4t/po4t
                                      ! Dean (1985) ! SR 01/01/12
!****        pH evaluation
             pht = -PH(k, i) - 2.1
             if(PH(k, i)<=0.0)pht = -14.0
             incr = 10.0
             do n = 1, 3
                 f = 1.0
                 incr = incr/10.0
                 iter = 0
                 do while (f>0.0 .AND. iter<12)
                     pht = pht + incr
                     hion = 10.0**pht
                     f = cart*k1*(hion + 2.0*k2)/(hion*hion + k1*hion + k1*k2) &
                       & + kw/hion - alkt - hion/ht                                ! SR 01/01/12
                                      ! SR 01/01/12
                                          ! SR 01/01/12
                     if(ammonia_buffering)f = f + ammt*kamm/(hion + kamm)
                 ! SR 01/01/12
                                        ! SR 01/01/12
                                                                                ! SR 01/01/12
                     if(phosphate_buffering)                                   &
                      & f = f + phost*(kp1*kp2*hion + 2*kp1*kp2*kp3 -          &
                      & hion*hion*hion)                                        &
                      & /(hion*hion*hion + kp1*hion*hion + kp1*kp2*hion +      &
                      & kp1*kp2*kp3)
                 ! SR 01/01/12
                     if(om_buffering)then
                                 ! SR 01/01/12
                         do ja = 1, nag
                        ! SR 01/01/12
                             f = f + omct*SDEN(ja)                             &
                               & *(1.0/(1.0 + hion*(10.0**PK(ja)))             &
                               & - 1.0/(1.0 + (10.0**(PK(ja)-4.5))))                                   ! SR 01/01/12
                         enddo
                   ! SR 01/01/12
                     endif
                 ! SR 01/01/12
                     iter = iter + 1
                 enddo
                 pht = pht - incr
             enddo
!****        pH, carbon dioxide, bicarbonate, and carbonate concentrations
             hion = 10.0**pht
             PH(k, i) = -pht
             CO2(k, i) = TIC(k, i)/(1.0 + k1/hion + k1*k2/(hion*hion))
             HCO3(k, i) = TIC(k, i)/(1.0 + hion/k1 + k2/hion)
             CO3(k, i) = TIC(k, i)/((hion*hion)/(k1*k2) + hion/k2 + 1.0)
         enddo
     enddo
     return
 
 
!**********************************************************
!**  SUBROUTINE ZOOPLANKTON                     **
!**********************************************************
 
     entry ZOOPLANKTON
     do i = iu, id
         do k = kt, KB(i)
             do jz = 1, nzp
                 zgztot = 0.0                                                                                            ! KV 5/9/2007
                 do jjz = 1, nzp
!                    ZGZTOT=ZGZTOT+ZGZ(K,I,JZ,JJZ)*ZOO(K,I,JZ)                
!                    ! KV 5/9/2007
                     zgztot = zgztot + zgz(k, i, jz, jjz)                                                             ! CB 5/26/07
                 enddo
                 ZOOSS(k, i, jz) = (ZMU(k, i, jz)*ZEFF(jz) - ZRT(k, i, jz)     &
                                 & - ZMT(k, i, jz))*ZOO(k, i, jz) - zgztot                   ! OMNIVOROUS ZOOPLANKTON    ! KV 5/9/2007
             enddo
         enddo
     enddo
     return
 
!***********************************************************************************************************************************
!**  D E R I V E D   C O N S T I T U E N T S                                  
!***********************************************************************************************************************************
 
!    **
     entry DERIVED_CONSTITUENTS
     apr = 0.0
     atot = 0.0
     totss = 0.0
     chla = 0.0
     cbodu = 0.0
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             do i = CUS(jb), DS(jb)
                 do k = kt, KB(i)
                     do ja = 1, nal
                         if(ALG_CALC(ja))apr(k, i) = apr(k, i)                 &
                          & + (AGR(k, i, ja) - ARR(k, i, ja))*ALG(k, i, ja)    &
                          & *H2(k, i)*day
                     enddo
                 enddo
                 do k = kt, KB(i)
                     cbodct = 0.0
                     cbodnt = 0.0
                     cbodpt = 0.0
                     bodtot = 0.0
                     algp = 0.0
                     algn = 0.0                                                           ! cb 6/6/10
                     do ja = 1, nal
                         if(ALG_CALC(ja))atot(k, i) = atot(k, i)               &
                          & + ALG(k, i, ja)
                     enddo
                     do ibod = 1, nbod
                         if(BOD_CALC(ibod))then
                             cbodct = cbodct + CBOD(k, i, ibod)*BODC(ibod)
                                                          ! cb 6/6/10
                             cbodnt = cbodnt + CBODN(k, i, ibod)
                                                          ! cb 6/6/10
                             cbodpt = cbodpt + CBODP(k, i, ibod)
                                                          ! cb 6/6/10
                             bodtot = bodtot + CBOD(k, i, ibod)
                             if(CBODS(ibod)>0.0)totss(k, i) = totss(k, i)      &
                              & + CBOD(k, i, ibod)/O2OM(jw)                                  ! SW 9/5/13  Added particulate CBOD to TSS computation
                         endif
                     enddo
                     dom(k, i) = LDOM(k, i) + RDOM(k, i)
                     pom(k, i) = LPOM(k, i) + RPOM(k, i)
                     DOC(k, i) = dom(k, i)*ORGC(jw) + cbodct
                                                          ! cb 6/6/10
                     POC(k, i) = pom(k, i)*ORGC(jw)
                     do ja = 1, nal
                         if(ALG_CALC(ja))then
                             POC(k, i) = POC(k, i) + ALG(k, i, ja)*AC(ja)
                             algp = algp + ALG(k, i, ja)*AP(ja)
                             algn = algn + ALG(k, i, ja)*AN(ja)
                         endif
                     enddo
                     if(zooplankton_calc)then
                         do jz = 1, nzp
                             POC(k, i) = POC(k, i) + ZC(jz)*ZOO(k, i, jz)
                                                     !MLM BAULK
                             zoop = ZOO(k, i, jz)*ZP(jz)
                                        !MLM BAULK
                             zoon = ZOO(k, i, jz)*ZN(jz)
                                        !MLM BAULK
                             cbodu(k, i) = cbodu(k, i) + O2OM(jw)*ZOO(k, i, jz)
                             totss(k, i) = totss(k, i) + ZOO(k, i, jz)
                                                                  ! SW 9/5/13  Added zooplankton to TSS computation
                         enddo
                     endif
                     TOC(k, i) = DOC(k, i) + POC(k, i)
                     DOP(k, i) = LDOM(k, i)*ORGPLD(k, i) + RDOM(k, i)          &
                               & *ORGPRD(k, i) + cbodpt                      ! CB 6/6/10
                     DON(k, i) = LDOM(k, i)*ORGNLD(k, i) + RDOM(k, i)          &
                               & *ORGNRD(k, i) + cbodnt                      ! CB 6/6/10
                     POP(k, i) = LPOM(k, i)*ORGPLP(k, i) + RPOM(k, i)          &
                               & *ORGPRP(k, i) + algp + zoop
                     PON(k, i) = LPOM(k, i)*ORGNLP(k, i) + RPOM(k, i)          &
                               & *ORGNRP(k, i) + algn + zoop
                     TOP(k, i) = DOP(k, i) + POP(k, i)
                     TON(k, i) = DON(k, i) + PON(k, i)
                     TKN(k, i) = TON(k, i) + NH4(k, i)
                     cbodu(k, i) = cbodu(k, i) + O2OM(jw)                      &
                                 & *(dom(k, i) + pom(k, i) + atot(k, i))       &
                                 & + bodtot
                     tpss = 0.0
                     do js = 1, nss
                         tpss = tpss + SS(k, i, js)*PARTP(jw)
                     enddo
                     TP(k, i) = TOP(k, i) + PO4(k, i) + tpss
                     TN(k, i) = TON(k, i) + NH4(k, i) + NO3(k, i)
                     O2DG(k, i) = (O2(k, i)/SATO(T1(k, i), TDS(k, i), PALT(i), &
                                & SALT_WATER(jw)))*100.0
                     do ja = 1, nal
                         if(ALG_CALC(ja))then
                             chla(k, i) = chla(k, i) + ALG(k, i, ja)/ACHLA(ja)
                             totss(k, i) = totss(k, i) + ALG(k, i, ja)
                         endif
                     enddo
                     totss(k, i) = totss(k, i) + TISS(k, i) + pom(k, i)
                 enddo
             enddo
         enddo
     enddo
     return
 
!********************************************************************************************************************
!**  A L K A L I N I T Y                                                **
!********************************************************************************************************************
     entry ALKALINITY
                 ! entire subroutine added ! SR 01/01/12
!    According to Stumm and Morgan (1996), table 4.5 on page 173:
!    Utilization of ammonium during photosynthesis results in an alkalinity
!    decrease: 14 eq. alk per 16 moles ammonium Utilization of nitrate during
!    photosynthesis results in an alkalinity increase: 18 eq. alk per 16 moles
!    nitrate Production of ammonium during respiration results in an
!    alkalinity increase: 14 eq. alk per 16 moles ammonium Nitrification of
!    ammonium results in an alkalinity decrease: 2 eq. alk per 1 mole ammonium
!    Denitrification of nitrate (to nitrogen gas) results in an alkalinity
 
 
!    increase: 1 eq. alk per 1 mole nitrate Alkalinity is represented as mg/L
!    CaCO3 (MW=100.088). CaCO3 has 2 equivalents of alk per mole. Nitrogen has
!    an atomic mass of 14.00674. These numbers account for the factor of
!    50.044/14.00674 used below.
     do i = iu, id
         do k = kt, KB(i)
             if(noncon_alkalinity)then
                 ALKSS(k, i) = (50.044/14.00674)                               &
                             & *(14./16.*(NH4AP(k, i) + NH4EP(k, i)            &
                             & + nh4zr(k, i) + nh4mr(k, i) - nh4mg(k, i))      &
                             & + 18./16.*(no3ag(k, i) + no3eg(k, i))           &
                             & - 2.*NH4D(k, i) + NO3D(k, i) + NO3SED(k, i)     &
                             & *(1 - FNO3SED(jw)))
             else
                 ALKSS(k, i) = 0.0
                          ! NW 2/11/16
             endif
         enddo
     enddo
 
     return
 
     entry DEALLOCATE_KINETICS
     deallocate(omtrm, sodtrm, nh4trm, no3trm, dom, pom, po4bod, nh4bod,       &
              & ticbod, atrm, atrmr, atrmf, etrm, etrmr, etrmf, bibh2)
     deallocate(lam2m)
     end subroutine KINETICS
