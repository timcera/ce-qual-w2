!*==withdrawal.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************************************************************************
!**                                            S U B R O U T I N E   W I T H D R A W A L                                          **
!***********************************************************************************************************************************
 
     subroutine WITHDRAWAL
 
     use GLOBAL                                                                                                                             ! cb 1/16/13
     use GEOMC
     use TVDC
     use SELWC
     use LOGICC
     use MAIN, ONLY:derived_calc, cdn, tdgon, jsg, nnsg, ndo, jwd, ngn2,       &
       & ngctdg, ea
     use SCREENC, ONLY:JDAY
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real :: Estrtest, Tempest
     integer :: Jjwd, Js
     intent (in) Estrtest
     intent (inout) Tempest
!
! Local variables
!
     real :: coef, dlrhob, dlrhomax, dlrhot, elr, elstr, elwd, fracv, hb, hswb,&
           & hswt, ht, hwdb, hwdt, qsumjs, qsumwd, ratio, rhofb, rhoft, vsum,  &
           & wsel
     real(r8) :: dosat, n2sat
     integer :: k, kbot, kstr, ktop, kwd
!
!*** End of declarations rewritten by SPAG
!
                              ! cb 11/7/17
 
     return                                   ! jwd
 
!***********************************************************************************************************************************
!**  D O W N S T R E A M   W I T H D R A W A L                                
!***********************************************************************************************************************************
 
!    **
     entry DOWNSTREAM_WITHDRAWAL(Js)
 
!    Variable initialization
 
     hswt = 0.0
 
!    Variable initialization
 
     hswb = 0.0
 
!    Variable initialization
 
     vnorm = 0.0
 
!    Variable initialization
 
     qnew = 0.0
 
!    Water surface elevation
 
     elr = SINA(jb)*DLX(id)*0.5
     wsel = ELWS(id) - elr              !EL(KT,ID)-Z(ID)*COSA(JB)
 
!    Structure layer
 
     do k = kt, KB(id)
         if(EL(k, id) - elr<ESTR(Js, jb))exit
     enddo
     kstr = MAX(k - 1, kt)
     kstr = MIN(kstr, KB(id))
 
!    Initial withdrawal limits
 
     ktop = MAX(KTSW(Js, jb), kt)
     if(kstr<ktop)ktop = kstr
     kbot = MIN(KBSW(Js, jb), KB(id))
     if(kbot<=kt .AND. kbot/=KB(id))kbot = kt + 1
     if(kbot>KB(id))kbot = KB(id)
     elstr = ESTR(Js, jb)
     if(ESTR(Js, jb)<=EL(KB(id) + 1, id + 1) - elr)then
         kstr = KB(id)
         elstr = EL(KB(id), id) - elr
     endif
     if(ESTR(Js, jb)>EL(kt, id) - elr)elstr = wsel
     if(KBSW(Js, jb)<kstr)then
         kstr = kt
         elstr = wsel
     endif
 
!    Boundary interference
 
     coef = 1.0
     if((wsel - EL(kbot, id) - elr)/=0.0)then
         ratio = (elstr - (EL(kbot, id) - elr))/(wsel - (EL(kbot, id) - elr))
         if(ratio<0.1 .OR. ratio>0.9)coef = 2.0
     endif
 
!    Withdrawal zone above structure
 
     do k = kstr - 1, ktop, -1
 
!**      Density frequency
 
         ht = (EL(k, id) - elr) - elstr
         rhoft = MAX                                                           &
               & (SQRT((ABS(RHO(k,id)-RHO(kstr,id)))/(ht*RHO(kstr,id) + nonzero&
               & )*g), nonzero)
 
!**      Thickness
 
         if(POINT_SINK(Js, jb))then
             hswt = (coef*QSTR(Js, jb)/rhoft)**0.333333
         else
             hswt = SQRT(2.0*coef*QSTR(Js, jb)/(WSTR(Js, jb)*rhoft))
         endif
         if(ht>=hswt)then
             ktop = k
             exit
         endif
     enddo
 
!    Reference density
 
     if((elstr + hswt)<wsel)then
         dlrhot = ABS(RHO(kstr, id) - RHO(ktop, id))
     elseif(wsel==elstr)then
         dlrhot = nonzero
     else
         dlrhot = ABS(RHO(kstr, id) - RHO(kt, id))*hswt/(wsel - elstr)
     endif
     dlrhot = MAX(dlrhot, nonzero)
 
!    Withdrawal zone below structure
 
     do k = kstr + 1, kbot
 
!**      Density frequency
 
         hb = elstr - (EL(k, id) - elr)
         rhofb = MAX                                                           &
               & (SQRT((ABS(RHO(k,id)-RHO(kstr,id)))/(hb*RHO(kstr,id) + nonzero&
               & )*g), nonzero)
 
!**      Thickness
 
         if(POINT_SINK(Js, jb))then
             hswb = (coef*QSTR(Js, jb)/rhofb)**0.333333
         else
             hswb = SQRT(2.0*coef*QSTR(Js, jb)/(WSTR(Js, jb)*rhofb))
         endif
         if(hb>=hswb)then
             kbot = k
             exit
         endif
     enddo
 
!    Reference density
 
     if((elstr - hswb)>EL(kbot + 1, id))then
         dlrhob = ABS(RHO(kstr, id) - RHO(kbot, id))
     elseif((EL(kbot + 1, id) - elr)==elstr)then                                                                       !SR 03/24/13
         dlrhob = nonzero                                                                                              !SR 03/24/13
     else
         dlrhob = ABS(RHO(kstr, id) - RHO(kbot, id))                           &
                & *hswb/(elstr - (EL(kbot + 1, id) - elr))
     endif
     dlrhob = MAX(dlrhob, nonzero)
 
!    Velocity profile
 
     vsum = 0.0
!    DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)                      ! GH 1/31/08
     do k = ktop, kbot
!        VNORM(K) = ABS(1.0-((RHO(K,ID)-RHO(KSTR,ID))/DLRHOMAX)**2)*BHR2(K,ID)
         if(k>kstr)then
             dlrhomax = MAX(dlrhob, 1.0E-10)                   !GH 1/31/08
         else
             dlrhomax = MAX(dlrhot, 1.0E-10)                   !GH 1/31/08
         endif
         vnorm(k) = 1.0 - ((RHO(k, id) - RHO(kstr, id))/dlrhomax)**2
         if(vnorm(k)>1.0)vnorm(k) = 1.0                           !GH 1/31/08
         if(vnorm(k)<0.0)vnorm(k) = 0.0                           !GH 1/31/08
         vnorm(k) = vnorm(k)*BHR2(k, id)
         vsum = vsum + vnorm(k)
     enddo
 
!    OUTFLOWS
     qsumjs = 0.0                                             ! SW 7/30/09
     TAVG(Js, jb) = 0.0                                              ! CB 5/12/10
     if(constituents)CAVG(Js, jb, cn(1:nac)) = 0.0
     if(derived_calc)CDAVG(Js, jb, cdn(1:NACD(jw), jw)) = 0.0
     do k = ktop, kbot
         qnew(k) = (vnorm(k)/vsum)*QSTR(Js, jb)
         QOUT(k, jb) = QOUT(k, jb) + qnew(k)
         TAVG(Js, jb) = TAVG(Js, jb) + qnew(k)*T2(k, id)      ! SW 7/30/09
         if(constituents)CAVG(Js, jb, cn(1:nac)) = CAVG(Js, jb, cn(1:nac))     &
          & + qnew(k)*C2(k, id, cn(1:nac))
         if(derived_calc)CDAVG(Js, jb, cdn(1:NACD(jw), jw))                    &
          & = CDAVG(Js, jb, cdn(1:NACD(jw), jw)) + qnew(k)                     &
          & *CD(k, id, cdn(1:NACD(jw), jw))
         qsumjs = qsumjs + qnew(k)
     enddo
     if(qsumjs>0.0)then
         TAVG(Js, jb) = TAVG(Js, jb)/qsumjs
         if(constituents)then             ! cb 1/16/13
             CAVG(Js, jb, cn(1:nac)) = CAVG(Js, jb, cn(1:nac))/qsumjs
             if(tdgon)then
                 call TOTAL_DISSOLVED_GAS(0, PALT(id), nnsg, jsg, TAVG(Js, jb),&
                   & CAVG(Js, jb, ndo))
                 if(ngctdg/=0)call TOTAL_DISSOLVED_GAS(1, PALT(id), nnsg, jsg, &
                  & TAVG(Js, jb), CAVG(Js, jb, ngn2))                                                ! n2 GAS
             endif
         endif
         if(derived_calc)then             ! cb 1/16/13
             CDAVG(Js, jb, cdn(1:NACD(jw), jw))                                &
               & = CDAVG(Js, jb, cdn(1:NACD(jw), jw))/qsumjs
    !if(tdgon)then                  ! cb 11/6/17
      !cdavg(js,jb,16)  = (cavg(js,jb,ndo)/exp(7.7117-1.31403*(log(tavg(js,jb)+45.93)))*palt(id))*100.0
             dosat = EXP(7.7117 - 1.31403*(LOG(TAVG(Js,jb) + 45.93)))*PALT(id)
             CDAVG(Js, jb, 16) = (CAVG(Js, jb, ndo)/dosat)*100.0
             if(ngctdg/=0)then
                 ea = DEXP(2.3026D0*(7.5D0*TDEW(jw)/(TDEW(jw) + 237.3D0)       &
                    & + 0.6609D0))*0.001316                                           ! in mm Hg   0.0098692atm=7.5006151mmHg
          !cdavg(js,jb,NDC)  = (cavg(js,jb,NGN2)/(1.5568D06*0.79*(PALT(ID)-EA)*(1.8816D-5 - 4.116D-7 * Tavg(js,jb) + 4.6D-9 * Tavg(js,jb)**2)))*100.0    ! SW 10/27/15
                 n2sat = 1.5568D06*0.79*(PALT(id) - ea)                        &
                       & *(1.8816D-5 - 4.116D-7*TAVG(Js, jb)                   &
                       & + 4.6D-9*TAVG(Js, jb)**2)
                 CDAVG(Js, jb, ndc) = 100.*(0.79*(CAVG(Js, jb, ngn2)/n2sat) +  &
                                    & 0.21*(CAVG(Js, jb, ndo)/dosat))
             endif
    !end if
         endif
     else
         TAVG(Js, jb) = -99.0
         if(constituents)CAVG(Js, jb, cn(1:nac)) = -99.0
         if(derived_calc)CDAVG(Js, jb, cdn(1:NACD(jw), jw)) = -99.0
     endif
 
!    Inactive layers and total outflow
 
     if(Js==nst)WHERE(QOUT(:, jb)==0.0)u(:, id) = 0.0
     return
!***********************************************************************************************************************************
!**  D O W N S T R E A M   W I T H D R A W A L  ESTIMATE                      
!***********************************************************************************************************************************
 
!    **
     entry DOWNSTREAM_WITHDRAWAL_ESTIMATE(Js, Tempest, Estrtest)
 
!    VARIABLE INITIALIZATION
 
     hswt = 0.0
 
!    VARIABLE INITIALIZATION
 
     hswb = 0.0
 
!    VARIABLE INITIALIZATION
 
     vnorm = 0.0
 
!    VARIABLE INITIALIZATION
 
     qnew = 0.0
 
!    Water surface elevation
 
     elr = SINA(jb)*DLX(id)*0.5
     wsel = EL(kt, id) - Z(id)*COSA(jb) - elr
 
!    Structure layer
 
     do k = kt, KB(id)
         if(EL(k, id) - elr<Estrtest)exit
     enddo
     kstr = MAX(k - 1, kt)
     kstr = MIN(kstr, KB(id))
 
!    Initial withdrawal limits
 
     ktop = MAX(KTSW(Js, jb), kt)
     if(kstr<ktop)ktop = kstr
     kbot = MIN(KBSW(Js, jb), KB(id))
     if(kbot<=kt .AND. kbot/=KB(id))kbot = kt + 1
     if(kbot>KB(id))kbot = KB(id)                                                                                      !SW 06/03/02
     elstr = Estrtest
     if(Estrtest<=EL(KB(id) + 1, id + 1) - elr)then                                                                 !SW 10/17/01
         kstr = KB(id)
         elstr = EL(KB(id), id) - elr                                                                                  !SW 10/17/01
     endif
     if(Estrtest>EL(kt, id) - elr)elstr = wsel
     if(KBSW(Js, jb)<kstr)then
         kstr = kt
         elstr = wsel                                                                                                  !SW 10/05/00
     endif
 
!    Boundary interference
 
     coef = 1.0
     if((wsel - EL(kbot, id) - elr)/=0.0)then
         ratio = (elstr - (EL(kbot, id) - elr))/(wsel - (EL(kbot, id) - elr))                                          !SW 10/17/01
         if(ratio<0.1 .OR. ratio>0.9)coef = 2.0
     endif
 
!    Withdrawal zone above structure
 
     do k = kstr - 1, ktop, -1
 
!**      Density frequency
 
         ht = (EL(k, id) - elr) - elstr
         rhoft = MAX                                                           &
               & (SQRT((ABS(RHO(k,id)-RHO(kstr,id)))/(ht*RHO(kstr,id) + nonzero&
               & )*g), nonzero)
 
!**      Thickness
 
         if(POINT_SINK(Js, jb))then
             hswt = (coef*QSTR(Js, jb)/rhoft)**0.333333
         else
             hswt = SQRT(2.0*coef*QSTR(Js, jb)/(WSTR(Js, jb)*rhoft))
         endif
         if(ht>=hswt)then
             ktop = k
             exit
         endif
     enddo
 
!    Reference density
 
     if((elstr + hswt)<wsel)then
         dlrhot = ABS(RHO(kstr, id) - RHO(ktop, id))
     elseif(wsel==elstr)then
         dlrhot = nonzero
     else
         dlrhot = ABS(RHO(kstr, id) - RHO(kt, id))*hswt/(wsel - elstr)
     endif
     dlrhot = MAX(dlrhot, nonzero)
 
!    Withdrawal zone below structure
 
     do k = kstr + 1, kbot
 
!**      Density frequency
 
         hb = elstr - (EL(k, id) - elr)                                                                                !SW 10/17/01
         rhofb = MAX                                                           &
               & (SQRT((ABS(RHO(k,id)-RHO(kstr,id)))/(hb*RHO(kstr,id) + nonzero&
               & )*g), nonzero)
 
!**      Thickness
 
         if(POINT_SINK(Js, jb))then
             hswb = (coef*QSTR(Js, jb)/rhofb)**0.333333
         else
             hswb = SQRT(2.0*coef*QSTR(Js, jb)/(WSTR(Js, jb)*rhofb))
         endif
         if(hb>=hswb)then
             kbot = k
             exit
         endif
     enddo
 
!    Reference density
 
     if((elstr - hswb)>EL(kbot + 1, id))then
         dlrhob = ABS(RHO(kstr, id) - RHO(kbot, id))
     elseif((EL(kbot + 1, id) - elr)==elstr)then                                                                       !SR 03/24/13
         dlrhob = nonzero                                                                                              !SR 03/24/13
     else
         dlrhob = ABS(RHO(kstr, id) - RHO(kbot, id))                           &
                & *hswb/(elstr - (EL(kbot + 1, id) - elr))                                                             !SW 10/17/01
     endif
     dlrhob = MAX(dlrhob, nonzero)
 
!    Velocity profile
 
     vsum = 0.0
!    DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)
     do k = ktop, kbot
         if(k>kstr)then
             dlrhomax = MAX(dlrhob, 1.0E-10)                   !GH 1/31/08
         else
             dlrhomax = MAX(dlrhot, 1.0E-10)                   !GH 1/31/08
         endif
         vnorm(k) = 1.0 - ((RHO(k, id) - RHO(kstr, id))/dlrhomax)**2
         if(vnorm(k)>1.0)vnorm(k) = 1.0                        !GH 1/31/08
         if(vnorm(k)<0.0)vnorm(k) = 0.0                            !GH 1/31/08
         vnorm(k) = vnorm(k)*BHR2(k, id)
         vsum = vsum + vnorm(k)
     enddo
 
!    Outflows
 
     Tempest = 0.0
     do k = ktop, kbot
         Tempest = Tempest + T2(k, id)*(vnorm(k)/vsum)*QSTR(Js, jb)
     enddo
 
     if(QSTR(Js, jb)>0.0)Tempest = Tempest/QSTR(Js, jb)
 
     return
!***********************************************************************************************************************************
!**  L A T E R A L   W I T H D R A W A L                                      
!***********************************************************************************************************************************
 
!    **
     entry LATERAL_WITHDRAWAL
 
!    Variable initialization
 
     vnorm = 0.0
 
!    Variable initialization
 
     qsw(:, jwd) = 0.0
 
!    Variable initialization
 
     hwdt = 0.0
 
!    Variable initialization
 
     hwdb = 0.0
 
!    Structure layer
 
     k = kt
     do k = kt, KB(i)
         if(EL(k, i)<EWD(jwd))exit
     enddo
     kwd = MAX(k - 1, kt)
     kwd = MIN(kwd, KB(i))
 
!    Initial withdrawal limits
 
     ktop = MAX(KTWD(jwd), kt)
     if(kwd<ktop)ktop = kwd
     kbot = MIN(KBWD(jwd), KB(i))
     if(kbot<=kt .AND. KB(i)/=kbot)kbot = kt + 1
     if(kbot>KB(i))kbot = KB(i)
     elwd = EWD(jwd)
     if(EWD(jwd)<=EL(KB(i) + 1, i))then
         kwd = KB(i)
         elwd = EL(KB(i), i)
     endif
     if(EWD(jwd)>EL(kt, i))elwd = EL(kt, i)
     if(KBWD(jwd)<kwd)then
         kwd = kt
         elwd = EL(kt, i)
     endif
 
!    Boundary interference
 
     coef = 1.0
     if(kt/=kbot)then
         ratio = (elwd - EL(kbot, i))/(EL(kt, i) - EL(kbot, i))
         if(ratio<0.1 .OR. ratio>0.9)coef = 2.0
     endif
 
!    Withdrawal zone above structure
 
     do k = kwd - 1, ktop, -1
 
!**      Density frequency
 
         ht = EL(k, i) - elwd
         rhoft = MAX(SQRT((ABS(RHO(k,i)-RHO(kwd,i)))/(ht*RHO(kwd,i) + nonzero) &
               & *g), nonzero)
 
!**      Thickness
 
         hwdt = (coef*QWD(jwd)/rhoft)**0.333333
         if(ht>=hwdt)then
             ktop = k
             exit
         endif
     enddo
 
!    Reference density
 
     if((elwd + hwdt)<EL(kt, i))then
         dlrhot = ABS(RHO(kwd, i) - RHO(ktop, i))
     elseif(EL(kt, i)==elwd)then
         dlrhot = nonzero
     else
         dlrhot = ABS(RHO(kwd, i) - RHO(kt, i))*hwdt/(EL(kt, i) - elwd)
     endif
     dlrhot = MAX(dlrhot, nonzero)
 
!    Withdrawal zone below structure
 
     do k = kwd + 1, kbot
 
!**      Density frequency
 
         hb = elwd - EL(k, i)
         rhofb = MAX(SQRT((ABS(RHO(k,i)-RHO(kwd,i)))/(hb*RHO(kwd,i) + nonzero) &
               & *g), nonzero)
 
!**      Thickness
 
         hwdb = (coef*QWD(jwd)/rhofb)**0.333333
         if(hb>=hwdb)then
             kbot = k
             exit
         endif
     enddo
 
!    Reference density
 
     if((elwd - hwdb)>EL(kbot + 1, i))then
         dlrhob = ABS(RHO(kwd, i) - RHO(kbot, i))
     elseif(EL(kbot + 1, i)==elwd)then                                                                                 !SR 03/24/13
         dlrhob = nonzero                                                                                              !SR 03/24/13
     else
         dlrhob = ABS(RHO(kwd, i) - RHO(kbot, i))*hwdb/(elwd - EL(kbot + 1, i))
     endif
     dlrhob = MAX(dlrhob, nonzero)
 
!    Velocity profile
 
     vsum = 0.0
!    DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)                                    
!    ! SW 1/24/05
     do k = ktop, kbot
!        VNORM(K) = ABS(1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2)*BHR2(K,I)
         if(k>kwd)then
             dlrhomax = MAX(dlrhob, 1.0E-10)                   !GH 1/31/08
         else
             dlrhomax = MAX(dlrhot, 1.0E-10)                   !GH 1/31/08
         endif
         vnorm(k) = 1.0 - ((RHO(k, i) - RHO(kwd, i))/dlrhomax)**2
         if(vnorm(k)>1.0)vnorm(k) = 1.0                           !GH 1/31/08
         if(vnorm(k)<0.0)vnorm(k) = 0.0                           !GH 1/31/08
         vnorm(k) = vnorm(k)*BHR2(k, i)
         vsum = vsum + vnorm(k)
     enddo
 
!    Outflows
     qsumwd = 0.0                                             ! SW 7/30/09
     TAVGW(jwd) = 0.0
     if(constituents)CAVGW(jwd, cn(1:nac)) = 0.0
     if(derived_calc)CDAVGW(jwd, cdn(1:NACD(jw), jw)) = 0.0
 
     do k = ktop, kbot
         fracv = (vnorm(k)/vsum)
         qsw(k, jwd) = qsw(k, jwd) + fracv*QWD(jwd)
         TAVGW(jwd) = TAVGW(jwd) + fracv*QWD(jwd)*T2(k, i)        ! SW 7/30/09
         if(constituents)CAVGW(jwd, cn(1:nac)) = CAVGW(jwd, cn(1:nac))         &
          & + fracv*QWD(jwd)*C2(k, i, cn(1:nac))
         if(derived_calc)CDAVGW(jwd, cdn(1:NACD(jw), jw))                      &
          & = CDAVGW(jwd, cdn(1:NACD(jw), jw)) + fracv*QWD(jwd)                &
          & *CD(k, i, cdn(1:NACD(jw), jw))
         qsumwd = qsumwd + fracv*QWD(jwd)
     enddo
  ! Debug
  !if(qwd(jwd)>0.0 .and. qsumwd <= 0.0)then
  !    write(9575,'(A,f8.3,1x,i5,1x,i5,1x,i5,f8.4,1x,f8.4,1x,f10.2)')'JDAY, ktop, kbot, kwd, qwd, qsumwd, elwd:',JDAY, ktop, kbot, kwd, qwd, qsumwd, elwd
  !endif
  ! Debug
     if(qsumwd>0.0)then
         TAVGW(jwd) = TAVGW(jwd)/qsumwd        ! SW 7/30/09
         if(constituents)then                  ! cb 1/16/13
             CAVGW(jwd, cn(1:nac)) = CAVGW(jwd, cn(1:nac))/qsumwd
             if(tdgon)then
                 call TOTAL_DISSOLVED_GAS(0, PALT(i), nnsg, jsg, TAVGW(jwd),   &
                   & CAVGW(jwd, ndo))
                 if(ngctdg/=0)call TOTAL_DISSOLVED_GAS(1, PALT(i), nnsg, jsg,  &
                  & TAVGW(jwd), CAVGW(jwd, ngn2))                                                   ! n2 GAS
             endif
         endif
         if(derived_calc)then
             CDAVGW(jwd, cdn(1:NACD(jw), jw))                                  &
               & = CDAVGW(jwd, cdn(1:NACD(jw), jw))/qsumwd
      !if(tdgon)then                ! cb 11/6/17
        !cdavgw(jwd,16)  = (cavgw(jwd,ndo)/exp(7.7117-1.31403*(log(tavgw(jwd)+45.93)))*palt(i))*100.0
             dosat = EXP(7.7117 - 1.31403*(LOG(TAVGW(jwd) + 45.93)))*PALT(i)
             CDAVGW(jwd, 16) = (CAVGW(jwd, ndo)/dosat)*100.0
             if(ngctdg/=0)then
                 ea = DEXP(2.3026D0*(7.5D0*TDEW(jw)/(TDEW(jw) + 237.3D0)       &
                    & + 0.6609D0))*0.001316                                           ! in mm Hg   0.0098692atm=7.5006151mmHg
          !cdavgw(jwd,NDC)  = (cavgw(jwd,NGN2)/(1.5568D06*0.79*(PALT(I)-EA)*(1.8816D-5 - 4.116D-7 * Tavgw(jwd) + 4.6D-9 * Tavgw(jwd)**2)))*100.0    ! SW 10/27/15
                 n2sat = 1.5568D06*0.79*(PALT(i) - ea)                         &
                       & *(1.8816D-5 - 4.116D-7*TAVGW(jwd) + 4.6D-9*TAVGW(jwd) &
                       & **2)
                 CDAVGW(jwd, ndc) = 100.*(0.79*(CAVGW(jwd, ngn2)/n2sat) + 0.21*&
                                  & (CAVGW(jwd, ndo)/dosat))
             endif
      !end if
         endif
     else
         TAVGW(jwd) = -99.0
         if(constituents)CAVGW(jwd, cn(1:nac)) = -99.0
         if(derived_calc)CDAVGW(jwd, cdn(1:NACD(jw), jw)) = -99.0
     endif
     KTW(jwd) = ktop
     KBW(jwd) = kbot
     return
!***********************************************************************************************************************************
!**  L A T E R A L   W I T H D R A W A L ESTIMATE                             
!***********************************************************************************************************************************
 
!    **
     entry LATERAL_WITHDRAWAL_ESTIMATE(Jjwd, Tempest, Estrtest)
 
!    VARIABLE INITIALIZATION
 
     vnorm = 0.0
 
!    VARIABLE INITIALIZATION
 
     qsw(:, Jjwd) = 0.0
 
!    VARIABLE INITIALIZATION
 
     hwdt = 0.0
 
!    VARIABLE INITIALIZATION
 
     hwdb = 0.0
 
!    Structure layer
 
     k = kt
     do k = kt, KB(i)
         if(EL(k, i)<Estrtest)exit
     enddo
     kwd = MAX(k - 1, kt)
     kwd = MIN(kwd, KB(i))
 
!    Initial withdrawal limits
 
     ktop = MAX(KTWD(Jjwd), kt)
     if(kwd<ktop)ktop = kwd
     kbot = MIN(KBWD(Jjwd), KB(i))
     if(kbot<=kt .AND. KB(i)/=kbot)kbot = kt + 1
     if(kbot>KB(i))kbot = KB(i)
     elwd = Estrtest
     if(Estrtest<=EL(KB(i) + 1, i))then
         kwd = KB(i)
         elwd = EL(KB(i), i)
     endif
     if(Estrtest>EL(kt, i))elwd = EL(kt, i)
     if(KBWD(Jjwd)<kwd)then
         kwd = kt
         elwd = EL(kt, i)
     endif
 
!    Boundary interference
 
     coef = 1.0
     if(kt/=kbot)then
         ratio = (elwd - EL(kbot, i))/(EL(kt, i) - EL(kbot, i))
         if(ratio<0.1 .OR. ratio>0.9)coef = 2.0
     endif
 
!    Withdrawal zone above structure
 
     do k = kwd - 1, ktop, -1
 
!**      Density frequency
 
         ht = EL(k, i) - elwd
         rhoft = MAX(SQRT((ABS(RHO(k,i)-RHO(kwd,i)))/(ht*RHO(kwd,i) + nonzero) &
               & *g), nonzero)
 
!**      Thickness
 
         hwdt = (coef*QWD(Jjwd)/rhoft)**0.333333
         if(ht>=hwdt)then
             ktop = k
             exit
         endif
     enddo
 
!    Reference density
 
     if((elwd + hwdt)<EL(kt, i))then
         dlrhot = ABS(RHO(kwd, i) - RHO(ktop, i))
     elseif(EL(kt, i)==elwd)then
         dlrhot = nonzero
     else
         dlrhot = ABS(RHO(kwd, i) - RHO(kt, i))*hwdt/(EL(kt, i) - elwd)
     endif
     dlrhot = MAX(dlrhot, nonzero)
 
!    Withdrawal zone below structure
 
     do k = kwd + 1, kbot
 
!**      Density frequency
 
         hb = elwd - EL(k, i)
         rhofb = MAX(SQRT((ABS(RHO(k,i)-RHO(kwd,i)))/(hb*RHO(kwd,i) + nonzero) &
               & *g), nonzero)
 
!**      Thickness
 
         hwdb = (coef*QWD(Jjwd)/rhofb)**0.333333
         if(hb>=hwdb)then
             kbot = k
             exit
         endif
     enddo
 
!    Reference density
 
     if((elwd - hwdb)>EL(kbot + 1, i))then
         dlrhob = ABS(RHO(kwd, i) - RHO(kbot, i))
     elseif(EL(kbot + 1, i)==elwd)then                                                                                 !SR 03/24/13
         dlrhob = nonzero                                                                                              !SR 03/24/13
     else
         dlrhob = ABS(RHO(kwd, i) - RHO(kbot, i))*hwdb/(elwd - EL(kbot + 1, i))
     endif
     dlrhob = MAX(dlrhob, nonzero)
 
!    Velocity profile
 
     vsum = 0.0
!    DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)                                    
!    ! SW 1/24/05
     do k = ktop, kbot
!        VNORM(K) = ABS(1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2)*BHR2(K,I)
         if(k>kwd)then
             dlrhomax = MAX(dlrhob, 1.0E-10)                   !GH 1/31/08
         else
             dlrhomax = MAX(dlrhot, 1.0E-10)                   !GH 1/31/08
         endif
         vnorm(k) = 1.0 - ((RHO(k, i) - RHO(kwd, i))/dlrhomax)**2
         if(vnorm(k)>1.0)vnorm(k) = 1.0                            !GH 1/31/08
         if(vnorm(k)<0.0)vnorm(k) = 0.0                            !GH 1/31/08
         vnorm(k) = vnorm(k)*BHR2(k, i)
         vsum = vsum + vnorm(k)
     enddo
 
!    Outflows
 
     do k = ktop, kbot
         Tempest = Tempest + T2(k, i)*(vnorm(k)/vsum)*QWD(Jjwd)
     enddo
     if(QWD(Jjwd)>0.0)Tempest = Tempest/QWD(Jjwd)
     KTW(Jjwd) = ktop
     KBW(Jjwd) = kbot
 
 
     end subroutine WITHDRAWAL
