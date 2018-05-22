!*==waterbody.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                            S U B R O U T I N E    W A T E R B O D Y                                           **
!***********************************************************************************************************************************
 
     subroutine WATERBODY
 
! Type declarations
     use GLOBAL
     use GEOMC
     use TVDC
     use LOGICC
     use PREC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(r8), dimension(kmx, imx) :: C, Ss
     intent (in) C
     intent (inout) Ss
!
! Local variables
!
     real :: b1, b2, brtot, el1, elw, frac, ht, q1, t1l, t2l, u1
     real, allocatable, dimension(:), save :: cl, ell, elr
     integer :: idt, iut, jjw, k, kl, kr
     real, allocatable, dimension(:, :), save :: qd, qu
!
!*** End of declarations rewritten by SPAG
!
 
 
!    Allocation declarations
 
     allocate(ell(kmx), elr(kmx), cl(nct), qu(kmx, imx), qd(kmx, imx))
 
!    Variable initialization
 
     ell = 0.0
 
!    Variable initialization
 
     elr = 0.0
 
!    Variable initialization
 
     cl = 0.0
 
!    Variable initialization
 
     qu = 0.0
 
!    Variable initialization
 
     qd = 0.0
 
!    Debug variable
!    ncount=0
!    End debug
 
     return
 
!***********************************************************************************************************************************
!**  U P S T R E A M   V E L O C I T Y                                        
!***********************************************************************************************************************************
 
!    **
     entry UPSTREAM_VELOCITY
     do jjb = 1, nbr
         if(UHS(jb)>=US(jjb) .AND. UHS(jb)<=DS(jjb))exit
     enddo
     do jjw = 1, nwb
         if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
     enddo
     do k = KTWB(jjw), KB(UHS(jb)) + 1
         ell(k) = EL(k, UHS(jb)) - SINA(jjb)*DLX(UHS(jb))*0.5
     enddo
     do k = kt, KB(iu) + 1
         elr(k) = EL(k, iu) + SINA(jb)*DLX(iu)*0.5
     enddo
     elw = ELWS(UHS(jb))                  !EL(KTWB(JJW),UHS(JB))-Z(UHS(JB))*COSA(JJB)
     el1 = elw - SINA(jjb)*DLX(UHS(jb))*0.5
!    IF(SLOPE(JB) /= 0.0)THEN
!    EL1  =
!    ELW-(ELWS(UHS(JB)+1)-ELW)/(0.5*(DLX(UHS(JB)+DLX(JB)+1))*DLX(UHS(JB))*0.5 
!    ! SW 7/17/09 ELSE EL1=ELW
!    ENDIF
     elw = el1
     kl = KTWB(jjw) + 1
     do k = kt + 1, KB(iu) + 1
         if(elr(k)>=ell(kl))then
             if(kl==KTWB(jjw) + 1)then
                 q1 = U(kl - 1, UHS(jb))*BHR1(KTWB(jjw), UHS(jb))
                 if(kl==KB(UHS(jb)) + 1 .AND. ell(kl)<elr(KB(iu) + 1))then
                     ht = elw - elr(KB(iu) + 1)
                 else
                     ht = H1(KTWB(jjw), UHS(jb))
                 endif
             else
                 q1 = U(kl - 1, UHS(jb))*BHR1(kl - 1, UHS(jb))
                 ht = H1(kl - 1, UHS(jb))
             endif
             if(k==kt + 1)then
                 el1 = EL(KTWB(jjw), UHS(jb)) - Z(UHS(jb))*COSA(jjb)           &
                     & - SINA(jjb)*DLX(UHS(jb))*0.5
                 U(k - 1, iu - 1) = q1*((el1 - elr(k))/ht)/BHR1(kt, iu - 1)
             else
                 U(k - 1, iu - 1) = q1*((el1 - elr(k))/ht)/BHR1(k - 1, iu - 1)
             endif
             el1 = elr(k)
             if(elr(k)==ell(kl))kl = kl + 1
         else
             q1 = 0.0
             do while (elr(k)<=ell(kl))
                 if(kl==KTWB(jjw) + 1 .AND. k==kt + 1)then
                     q1 = q1 + U(kl - 1, UHS(jb))*BHR1(KTWB(jjw), UHS(jb))
                 elseif(kl==KTWB(jjw) + 1)then
                     q1 = q1 + U(kl - 1, UHS(jb))*BHR1(KTWB(jjw), UHS(jb))     &
                        & *(el1 - ell(kl))/H1(KTWB(jjw), UHS(jb))
                 else
                     frac = (el1 - ell(kl))/H1(kl - 1, UHS(jb))
                     q1 = q1 + U(kl - 1, UHS(jb))*BHR1(kl - 1, UHS(jb))*frac
                 endif
                 el1 = ell(kl)
                 kl = kl + 1
                 if(kl>KB(UHS(jb)))exit
             enddo
             if(k==kt + 1)then
                 brtot = BHR1(kt, iu - 1)
             else
                 brtot = BHR1(k - 1, iu - 1)
             endif
             frac = 0.0
!            IF (KL < KMX) THEN                             ! SW 6/29/06
             ht = H1(kl - 1, UHS(jb))
             frac = (el1 - elr(k))/ht
             if(KB(UHS(jb))>=KB(iu - 1) .AND. k>KB(iu - 1))                    &
              & frac = (el1 - ell(kl))/ht                                           ! SW 6/29/06
             if(kl>=kmx .AND. elr(k)<ell(kl))frac = 1.0         ! SW 6/29/06
             q1 = q1 + U(kl - 1, UHS(jb))*BHR1(kl - 1, UHS(jb))*frac
!            ELSE                                           ! SW 6/29/06
!            Q1 = Q1+U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))   ! SW 6/29/06
!            END IF                                         ! SW 6/29/06
             U(k - 1, iu - 1) = q1/brtot
             if(kl>KB(UHS(jb)))then
                 if(frac<1.0 .AND. frac/=0.0)then
                     if(k==KB(iu) + 1)then
                         U(k - 1, iu - 1)                                      &
                           & = (q1 + U(kl - 1, UHS(jb))*BHR1(kl - 1, UHS(jb))  &
                           & *(1.0 - frac))/brtot
                     else
                         U(k, iu - 1) = U(kl - 1, UHS(jb))                     &
                                      & *BHR1(kl - 1, UHS(jb))*(1.0 - frac)    &
                                      & /BHR(k, iu - 1)
                     endif
                 endif
                 exit
             endif
             el1 = elr(k)
         endif
     enddo
 
!    Debug
!    QL=0.0; QR=0.0
!    do k=ktwb(jjw),kb(uhs(jb))
!    QL=QL+u(k,uhs(jb))*BHR1(k,uhs(jb))
!    enddo
!    do k=kt,kb(iu-1)
!    QR=QR+u(k,iu-1)*BHR1(k,iu-1)
!    enddo
!    if((QR-QL)/QL > 0.02)then
!    if(ncount.eq.0)open(1299,file='debug_out.txt',status='unknown')
!    ncount=ncount+1
!    write(1299,*)'***QL=',QL,' QR=',QR
!    write(1299,*)'UHS(JB)=',uhs(jb),'  iu-1=',iu-1,' ELWS=', elw
!    write(1299,*)'LEFT: K      U     BHR1    ELEV      H1     BHR    ELL'
!    do k=ktwb(jjw),kb(uhs(jb))
! write(1299,'(i8,f8.3,f8.3,5f8.3)')k,u(k,uhs(jb)),BHR1(k,uhs(jb)),el(k,uhs(jb)
!    enddo
!    write(1299,*)'RIGHT: K     U     BHR1    ELEV      H1     BHR    ELR'
!    do k=kt,kb(iu-1)
! write(1299,'(i8,f8.3,f8.3,5f8.3)')k,u(k,iu-1),BHR1(k,iu-1),el(k,iu-1),h1(k,iu
!    enddo
!    end if
!    End debug
 
     return
 
!***********************************************************************************************************************************
!**  U P S T R E A M   W A T E R B O D Y                                      
!***********************************************************************************************************************************
 
!    **
     entry UPSTREAM_WATERBODY
     do jjb = 1, nbr
         if(UHS(jb)>=US(jjb) .AND. UHS(jb)<=DS(jjb))exit
     enddo
     do jjw = 1, nwb
         if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
     enddo
     do k = KTWB(jjw), KB(UHS(jb)) + 1
         ell(k) = EL(k, UHS(jb)) - SINA(jjb)*DLX(UHS(jb))*0.5
     enddo
     do k = kt, KB(iu) + 1
         elr(k) = EL(k, iu) + SINA(jb)*DLX(iu)*0.5
     enddo
     elw = ELWS(UHS(jb))                              !EL(KTWB(JJW),UHS(JB))-Z(UHS(JB))*COSA(JJB)
     el1 = elw - SINA(jjb)*DLX(UHS(jb))*0.5
 ! EL1  = ELW-(ELWS(UHS(JB)+1)-ELW)/(0.5*(DLX(UHS(JB)+DLX(JB)+1))*DLX(UHS(JB))*0.5   ! SW 7/17/09
     kl = KTWB(jjw) + 1
     do k = kt + 1, KB(iu) + 1
         if(elr(k)>=ell(kl))then
             T1(k - 1, iu - 1) = T1(kl - 1, UHS(jb))
             T2(k - 1, iu - 1) = T2(kl - 1, UHS(jb))
             C1S(k - 1, iu - 1, cn(1:nac)) = C1S(kl - 1, UHS(jb), cn(1:nac))
             C1(k - 1, iu - 1, cn(1:nac)) = C1S(kl - 1, UHS(jb), cn(1:nac))
             C2(k - 1, iu - 1, cn(1:nac)) = C1S(kl - 1, UHS(jb), cn(1:nac))
             el1 = elr(k)
             if(elr(k)==ell(kl))kl = kl + 1
         else
             brtot = 0.0
             cl = 0.0
             t1l = 0.0
             t2l = 0.0
             do while (elr(k)<=ell(kl))
                 if(kl==KTWB(jjw) + 1 .AND. k==kt + 1)then
                     b1 = BH2(KTWB(jjw), UHS(jb))
                 else
                     b1 = B(kl - 1, UHS(jb))*(el1 - ell(kl))
                 endif
                 brtot = brtot + b1
                 t1l = t1l + b1*T1(kl - 1, UHS(jb))
                 t2l = t2l + b1*T2(kl - 1, UHS(jb))
                 cl(cn(1:nac)) = cl(cn(1:nac))                                 &
                               & + b1*C1S(kl - 1, UHS(jb), cn(1:nac))
                 el1 = ell(kl)
                 kl = kl + 1
                 if(kl>KB(UHS(jb)) + 1)exit
             enddo
             if(kl<=KB(UHS(jb)) + 1)then
                 b1 = B(kl - 1, UHS(jb))*(el1 - elr(k))
                 brtot = brtot + b1
                 if(brtot>0.0)then
                     T1(k - 1, iu - 1) = (t1l + b1*T1(kl - 1, UHS(jb)))/brtot
                     T2(k - 1, iu - 1) = (t2l + b1*T2(kl - 1, UHS(jb)))/brtot
                     C1S(k - 1, iu - 1, cn(1:nac))                             &
                       & = (cl(cn(1:nac)) + b1*C1S(kl - 1, UHS(jb), cn(1:nac)))&
                       & /brtot
                     C1(k - 1, iu - 1, cn(1:nac))                              &
                       & = C1S(k - 1, iu - 1, cn(1:nac))
                     C2(k - 1, iu - 1, cn(1:nac))                              &
                       & = C1S(k - 1, iu - 1, cn(1:nac))
                 else
                     T1(k - 1, iu - 1) = T1(kl - 1, UHS(jb))
                     T2(k - 1, iu - 1) = T2(kl - 1, UHS(jb))
                     C1S(k - 1, iu - 1, cn(1:nac))                             &
                       & = C1S(kl - 1, UHS(jb), cn(1:nac))
                     C1(k - 1, iu - 1, cn(1:nac))                              &
                       & = C1S(k - 1, iu - 1, cn(1:nac))
                     C2(k - 1, iu - 1, cn(1:nac))                              &
                       & = C1S(k - 1, iu - 1, cn(1:nac))
                 endif
             else
                 if(brtot>0.0)then
                     T1(k - 1, iu - 1) = t1l/brtot
                     T2(k - 1, iu - 1) = t2l/brtot
                     C1S(k - 1, iu - 1, cn(1:nac)) = cl(cn(1:nac))/brtot
                     C1(k - 1, iu - 1, cn(1:nac))                              &
                       & = C1S(k - 1, iu - 1, cn(1:nac))
                     C2(k - 1, iu - 1, cn(1:nac))                              &
                       & = C1S(k - 1, iu - 1, cn(1:nac))
                 else
                     T1(k - 1, iu - 1) = T1(kl - 1, UHS(jb))
                     T2(k - 1, iu - 1) = T2(kl - 1, UHS(jb))
                     C1S(k - 1, iu - 1, cn(1:nac))                             &
                       & = C1S(kl - 1, UHS(jb), cn(1:nac))
                     C1(k - 1, iu - 1, cn(1:nac))                              &
                       & = C1S(k - 1, iu - 1, cn(1:nac))
                     C2(k - 1, iu - 1, cn(1:nac))                              &
                       & = C1S(k - 1, iu - 1, cn(1:nac))
                 endif
                 exit
             endif
             el1 = elr(k)
         endif
     enddo
     return
 
!***********************************************************************************************************************************
!**  D O W N S T R E A M   W A T E R B O D Y                                  
!***********************************************************************************************************************************
 
!    **
     entry DOWNSTREAM_WATERBODY
     do jjb = 1, nbr
         if(CDHS(jb)>=CUS(jjb) .AND. CDHS(jb)<=DS(jjb))exit
     enddo
     do jjw = 1, nwb
         if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
     enddo
     kr = KTWB(jjw) + 1
     do k = kt, KB(id) + 1
         ell(k) = EL(k, id) - SINA(jb)*DLX(id)*0.5
     enddo
     do k = KTWB(jjw), KB(CDHS(jb)) + 1
         elr(k) = EL(k, CDHS(jb)) + SINA(jjb)*DLX(CDHS(jb))*0.5
     enddo
     elw = ELWS(id)                                          !EL(KTWB(JW),ID)-Z(ID)*COSA(JB)
     el1 = elw - SINA(jb)*DLX(id)*0.5
!    EL1  = ELW+(ELW-ELWS(ID-1))/(0.5*(DLX(ID)+DLX(ID-1))*DLX(ID)*0.5   ! SW
!    7/17/09
     idt = id + 1
     do k = kt + 1, KB(id) + 1
         if(ell(k)>=elr(kr))then
             T1(k - 1, idt) = T1(kr - 1, CDHS(jb))
             T2(k - 1, idt) = T2(kr - 1, CDHS(jb))
             C1S(k - 1, idt, cn(1:nac)) = C1S(kr - 1, CDHS(jb), cn(1:nac))
             C1(k - 1, idt, cn(1:nac)) = C1S(kr - 1, CDHS(jb), cn(1:nac))
             C2(k - 1, idt, cn(1:nac)) = C1S(kr - 1, CDHS(jb), cn(1:nac))
             el1 = ell(k)
             if(ell(k)==elr(kr))kr = kr + 1
             if(kr>KB(CDHS(jb)) + 1)exit
         else
             brtot = 0.0
             cl = 0.0
             t1l = 0.0
             t2l = 0.0
             do while (ell(k)<=elr(kr))
                 if(kr==KTWB(jjw) + 1 .AND. k==kt + 1)then
                     b1 = BH2(KTWB(jjw), CDHS(jb))
                     brtot = brtot + b1
                 else
                     b1 = B(kr - 1, CDHS(jb))*(el1 - elr(kr))
                     brtot = brtot + b1
                 endif
                 t1l = t1l + b1*T1(kr - 1, CDHS(jb))
                 t2l = t2l + b1*T2(kr - 1, CDHS(jb))
                 cl(cn(1:nac)) = cl(cn(1:nac))                                 &
                               & + b1*C1S(kr - 1, CDHS(jb), cn(1:nac))
                 el1 = elr(kr)
                 kr = kr + 1
                 if(kr>KB(CDHS(jb)) + 1)exit
             enddo
             if(kr<=KB(CDHS(jb) + 1))then
                 b1 = B(kr - 1, CDHS(jb))*(el1 - ell(k))
             else
                 b1 = 0.0
             endif
             brtot = brtot + b1
             if(brtot==0.0)exit
             T1(k - 1, idt) = (t1l + b1*T1(kr - 1, CDHS(jb)))/brtot
             T2(k - 1, idt) = (t2l + b1*T2(kr - 1, CDHS(jb)))/brtot
             C1S(k - 1, idt, cn(1:nac))                                        &
               & = (cl(cn(1:nac)) + b1*C1S(kr - 1, CDHS(jb), cn(1:nac)))/brtot
             C1(k - 1, idt, cn(1:nac)) = C1S(k - 1, idt, cn(1:nac))
             C2(k - 1, idt, cn(1:nac)) = C1S(k - 1, idt, cn(1:nac))
             el1 = ell(k)
         endif
     enddo
     return
 
!***********************************************************************************************************************************
!**  U P S T R E A M   B R A N C H                                            
!***********************************************************************************************************************************
 
!    **
     entry UPSTREAM_BRANCH
     do jjw = 1, nwb
         if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
     enddo
     do k = kt, KB(i) + 1
         ell(k) = EL(k, i)
     enddo
     do k = KTWB(JWUH(jb)), KB(CUS(jjb)) + 1
         elr(k) = EL(k, CUS(jjb)) + SINA(jjb)*DLX(CUS(jjb))*0.5
     enddo
     elw = ELWS(CUS(jjb))                             !EL(KTWB(JWUH(JB)),CUS(JJB))-Z(CUS(JJB))*COSA(JJB)
     el1 = elw + SINA(jjb)*DLX(CUS(jjb))*0.5
!    EL1  =
! ELW-(ELWS(CUS(JJB)+1)-ELW)/(0.5*(DLX(CUS(JJB)+DLX(CUS(JJB)+1))*DLX(CUS(JJB))*
     kr = KTWB(JWUH(jb)) + 1
     do k = kt + 1, KB(i) + 1
         if(ell(k)>=elr(kr))then
             q1 = VOLUH2(kr - 1, jjb)/dlt
             b1 = BHR(kr - 1, CUS(jjb) - 1)
             if(kr==KTWB(JWUH(jb)) + 1)b1 = BHR2(kt, CUS(jjb) - 1)
             u1 = U(kr - 1, CUS(jjb) - 1)*b1
             ht = H2(kr - 1, JWUH(jb))
             if(kr==KTWB(JWUH(jb)) + 1)ht = H2(kt, CUS(jjb))
             q1 = q1*((el1 - ell(k))/ht)
             b2 = BHR(k - 1, i)
             if(k==KTWB(jw) + 1)b2 = BHR2(kt, i)
             UXBR(k - 1, i) = UXBR(k - 1, i)                                   &
                            & + (ABS((u1/b2)*COS(betabr)*q1)/DLX(i))
             UYBR(k - 1, i) = UYBR(k - 1, i) + ABS(q1*SIN(betabr))
             el1 = ell(k)
             if(ell(k)==elr(kr))kr = kr + 1
         else
             u1 = 0.0
             q1 = 0.0
             brtot = 0.0
             do while (ell(k)<=elr(kr))
                 if(kr/=KTWB(JWUH(jb)) + 1)then
                     frac = (el1 - elr(kr))/(H(kr - 1, JWUH(jb)))
                     b1 = BHR(kr - 1, CUS(jjb) - 1)
                 else
                     frac = (el1 - elr(kr))/H2(kt, CUS(jjb))
                     b1 = BHR2(kt, CUS(jjb) - 1)
                 endif
                 u1 = u1 + U(kr - 1, CUS(jjb) - 1)*b1*frac
                 q1 = q1 + VOLUH2(kr - 1, jjb)/dlt*frac
                 el1 = elr(kr)
                 kr = kr + 1
                 if(kr>KB(CUS(jjb) + 1))exit
             enddo
             if(k==KTWB(jw) + 1)then
                 b2 = BHR2(kt, i)
             else
                 b2 = BHR(k - 1, i)
             endif
             if(H(kr - 1, JWUH(jb))/=0.0)then
                 if(kr - 1==KTWB(JWUH(jb)))then
                     ht = H2(kt, CUS(jjb))
                 else
                     ht = H2(kr - 1, JWUH(jb))
                 endif
                 frac = (el1 - ell(kr - 1))/ht
                 q1 = q1 + frac*VOLUH2(kr - 1, jjb)/dlt
                 UXBR(k - 1, i) = UXBR(k - 1, i)                               &
                                & + (ABS((u1/b2)*COS(betabr)*q1)/DLX(i))
                 UYBR(k - 1, i) = UYBR(k - 1, i) + ABS(q1*SIN(betabr))
             endif
             if(kr>KB(CUS(jjb) + 1))exit
             el1 = ell(k)
         endif
     enddo
     return
 
!***********************************************************************************************************************************
!**  D O W N S T R E A M   B R A N C H                                        
!***********************************************************************************************************************************
 
!    **
     entry DOWNSTREAM_BRANCH
     do jjw = 1, nwb
         if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
     enddo
     do k = kt, KB(i) + 1
         elr(k) = EL(k, i)
     enddo
     do k = KTWB(jjw), KB(DS(jjb)) + 1
         ell(k) = EL(k, DS(jjb)) + SINA(jjb)*DLX(DS(jjb))*0.5
     enddo
     elw = ELWS(DS(jjb))                                   !EL(KTWB(JJW),DS(JJB))-Z(DS(JJB))*COSA(JJB)
     el1 = elw - SINA(jjb)*DLX(DS(jjb))*0.5
!    IF(SLOPE(JJB) /= 0.0)THEN
!    EL1  =
! ELW+(ELW-ELWS(DS(JJB)-1))/(0.5*(DLX(DS(JJB)+DLX(DS(JJB)+1))*DLX(DS(JJB))*0.5 
!    EL1=ELW
!    ENDIF
     kl = KTWB(jjw) + 1
     do k = kt + 1, KB(i) + 1
         if(elr(k)>=ell(kl))then
             q1 = VOLDH2(kl - 1, jjb)/dlt
             b1 = BHR(kl - 1, DS(jjb) - 1)
             if(kl==KTWB(jjw) + 1)b1 = BHR2(kt, DS(jjb))
             u1 = U(kl - 1, DS(jjb) - 1)*b1
             ht = H2(kl - 1, jjw)
             if(kl==KTWB(jjw) + 1)ht = H2(kt, DS(jjb))
             q1 = q1*((el1 - elr(k))/ht)
             b2 = BHR(k - 1, i)
             if(k==KTWB(jw) + 1)b2 = BHR2(kt, i)
             UXBR(k - 1, i) = UXBR(k - 1, i)                                   &
                            & + (ABS((u1/b2)*COS(betabr)*q1)/DLX(i))
             UYBR(k - 1, i) = UYBR(k - 1, i) + ABS(q1*SIN(betabr))
             el1 = elr(k)
             if(elr(k)==ell(kl))kl = kl + 1
         else
             u1 = 0.0
             q1 = 0.0
             brtot = 0.0
             do while (elr(k)<=ell(kl))
                 if(kl/=KTWB(jjw) + 1)then
                     frac = (el1 - ell(kl))/H(kl - 1, jjw)
                     b1 = BHR(kl - 1, DS(jjb) - 1)
                 else
                     frac = (el1 - ell(kl))/H2(kt, DS(jjb))
                     b1 = BHR2(kt, DS(jjb))
                 endif
                 u1 = u1 + U(kl - 1, DS(jjb))*b1*frac
                 q1 = q1 + (VOLDH2(kl - 1, jjb)/dlt)*frac
                 el1 = ell(kl)
                 kl = kl + 1
                 if(kl>KB(DS(jjb) + 1))exit
             enddo
             if(k==KTWB(jw) + 1)then
                 b2 = BHR2(kt, i)
             else
                 b2 = BHR(k - 1, i)
             endif
             if(H(kl - 1, jjw)/=0.0)then
                 if(kl - 1==KTWB(jjw))then
                     ht = H2(kt, DS(jjb))
                 else
                     ht = H2(kl - 1, jjw)
                 endif
                 frac = (el1 - ell(kl - 1))/ht
                 q1 = q1 + frac*(VOLDH2(kl - 1, jjb)/dlt)
                 UXBR(k - 1, i) = UXBR(k - 1, i)                               &
                                & + (ABS((u1/b2)*COS(betabr)*q1)/DLX(i))
                 UYBR(k - 1, i) = UYBR(k - 1, i) + ABS(q1*SIN(betabr))
             endif
             if(kl>KB(DS(jjb) + 1))exit
             el1 = elr(k)
         endif
     enddo
     return
 
!***********************************************************************************************************************************
!**  U P S T R E A M   F L O W                                                
!***********************************************************************************************************************************
 
!    **
     entry UPSTREAM_FLOW
     do jjb = 1, nbr
         if(UHS(jb)>=US(jjb) .AND. UHS(jb)<=DS(jjb))exit
     enddo
     do jjw = 1, nwb
         if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
     enddo
     do k = KTWB(JWUH(jb)), KB(UHS(jb)) + 1
         ell(k) = EL(k, UHS(jb))
     enddo
     do k = kt, KB(iu) + 1
         elr(k) = EL(k, iu) + SINA(jb)*DLX(iu)*0.5
     enddo
     elw = ELWS(iu)                               !EL(KTWB(JW),IU)-Z(IU)*COSA(JB)
     el1 = elw + SINA(jb)*DLX(iu)*0.5
!    IF(SLOPE(JB) /= 0.0)THEN
!    EL1  = ELW-(ELWS(IU+1)-ELW)/(0.5*(DLX(IU)+DLX(IU+1))*DLX(IU)*0.5   ! SW
!    7/17/09 ELSE
!    EL1=ELW
!    ENDIF
     kr = kt + 1
     do k = KTWB(JWUH(jb)) + 1, KB(UHS(jb)) + 1
         if(ell(k)>=elr(kr))then
             q1 = VOLUH2(kr - 1, jb)/dlt
             ht = H2(kr - 1, jw)
             if(kr==KTWB(jw) + 1)then
                 qu(k - 1, UHS(jb)) = q1
                 QSS(k - 1, UHS(jb)) = QSS(k - 1, UHS(jb)) - q1
             else
                 qu(k - 1, UHS(jb)) = q1*((el1 - ell(k))/ht)
                 QSS(k - 1, UHS(jb)) = QSS(k - 1, UHS(jb)) - qu(k - 1, UHS(jb))
             endif
             el1 = ell(k)
             if(ell(k)==elr(kr))kr = kr + 1
         else
             q1 = 0.0
             do while (ell(k)<=elr(kr))
                 if(kr/=KTWB(jw) + 1)then
                     frac = (el1 - elr(kr))/(H(kr - 1, jw))
                 else
                     frac = (el1 - elr(kr))/(H2(kt, iu))
                 endif
                 q1 = q1 + VOLUH2(kr - 1, jb)/dlt*frac
                 el1 = elr(kr)
                 kr = kr + 1
                 if(kr>KB(iu) + 1)exit
             enddo
             if(H(kr - 1, jw)/=0.0)then
                 frac = (el1 - ell(k))/H(kr - 1, jw)
                 qu(k - 1, UHS(jb)) = q1 + frac*VOLUH2(kr - 1, jb)/dlt
             else
                 qu(k - 1, UHS(jb)) = q1
             endif
             QSS(k - 1, UHS(jb)) = QSS(k - 1, UHS(jb)) - qu(k - 1, UHS(jb))
             if(kr>KB(iu) + 1)exit
             el1 = ell(k)
         endif
     enddo
     return
 
!***********************************************************************************************************************************
!**  D O W N S T R E A M   F L O W                                            
!***********************************************************************************************************************************
 
!    **
     entry DOWNSTREAM_FLOW
     do jjb = 1, nbr
         if(CDHS(jb)>=US(jjb) .AND. CDHS(jb)<=DS(jjb))exit
     enddo
     do jjw = 1, nwb
         if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
     enddo
     do k = KTWB(jjw), KB(CDHS(jb)) + 1
         elr(k) = EL(k, CDHS(jb))
     enddo
     do k = kt, KB(id) + 1
         ell(k) = EL(k, id) - SINA(jb)*DLX(id)*0.5
     enddo
     elw = ELWS(id)                          !EL(KTWB(JW),ID)-Z(ID)*COSA(JB)
     el1 = elw - SINA(jb)*DLX(id)*0.5
!    IF(SLOPE(JB) /= 0.0)THEN
!    EL1  = ELW+(ELW-ELWS(ID-1))/(0.5*(DLX(ID)+DLX(ID-1))*DLX(ID)*0.5   ! SW
!    7/17/09 ELSE
!    EL1=ELW
!    ENDIF
     kl = kt + 1
     do k = KTWB(jjw) + 1, KB(CDHS(jb)) + 1
         if(elr(k)>=ell(kl))then
             q1 = VOLDH2(kl - 1, jb)/dlt
             ht = H2(kl - 1, jw)
             if(kl==KTWB(jw) + 1)ht = H2(kt, id)
             qd(k - 1, CDHS(jb)) = q1*((el1 - elr(k))/ht)
             QSS(k - 1, CDHS(jb)) = QSS(k - 1, CDHS(jb)) + qd(k - 1, DHS(jb))
             el1 = elr(k)
             if(elr(k)==ell(kl))kl = kl + 1
         else
             q1 = 0.0
             do while (elr(k)<=ell(kl))
                 if(kl/=KTWB(jw) + 1)then
                     frac = (el1 - ell(kl))/(H(kl - 1, jw))
                 else
                     frac = (el1 - ell(kl))/(H2(kt, id))
                 endif
                 q1 = q1 + (VOLDH2(kl - 1, jb)/dlt)*frac
                 el1 = ell(kl)
                 kl = kl + 1
                 if(kl>KB(id) + 1)exit
             enddo
             if(H(kl - 1, jw)/=0.0)then
                 frac = (el1 - elr(k))/H(kl - 1, jw)
                 qd(k - 1, CDHS(jb)) = q1 + frac*(VOLDH2(kl - 1, jb)/dlt)
             else
                 qd(k - 1, CDHS(jb)) = q1
             endif
             QSS(k - 1, CDHS(jb)) = QSS(k - 1, CDHS(jb)) + qd(k - 1, CDHS(jb))
             if(kl>KB(id) + 1)exit
             el1 = elr(k)
         endif
     enddo
     return
 
!***********************************************************************************************************************************
!**  U P S T R E A M   C O N S T I T U E N T                                  
!***********************************************************************************************************************************
 
!    **
     entry UPSTREAM_CONSTITUENT(C, Ss)
     do jjb = 1, nbr
         if(UHS(jb)>=US(jjb) .AND. UHS(jb)<=DS(jjb))exit
     enddo
     do jjw = 1, nwb
         if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
     enddo
     do k = KTWB(JWUH(jb)), KB(UHS(jb)) + 1
         ell(k) = EL(k, UHS(jb))
     enddo
     do k = kt, KB(iu) + 1
         elr(k) = EL(k, iu) + SINA(jb)*DLX(iu)*0.5
     enddo
     elw = ELWS(iu)                          !EL(KTWB(JW),IU)-Z(IU)*COSA(JB)
     el1 = elw + SINA(jb)*DLX(iu)*0.5
!    IF(SLOPE(JB) /= 0.0)THEN
!    EL1  = ELW-(ELWS(IU+1)-ELW)/(0.5*(DLX(IU)+DLX(IU+1))*DLX(IU)*0.5   ! SW
!    7/17/09 ELSE
!    EL1=ELW
!    ENDIF
     kr = kt + 1
     do k = KTWB(JWUH(jb)) + 1, KB(UHS(jb)) + 1
         iut = iu
         if(qu(k - 1, UHS(jb))>=0.0)iut = iu - 1
         if(ell(k)>=elr(kr))then
             t1l = C(kr - 1, iut)
             Ss(k - 1, UHS(jb)) = Ss(k - 1, UHS(jb)) - t1l*qu(k - 1, UHS(jb))
             el1 = ell(k)
             if(ell(k)==elr(kr))kr = kr + 1
         else
             t1l = 0.0
             brtot = 0.0
             do while (ell(k)<=elr(kr))
                 if(kr==kt + 1 .AND. k==KTWB(JWUH(jb)) + 1)then
                     b1 = BH2(kt, iu)
                     brtot = brtot + b1
                 else
                     b1 = B(kr - 1, iu)*(el1 - elr(kr))
                     brtot = brtot + b1
                 endif
                 iut = iu
                 if(qu(k - 1, UHS(jb))>=0.0)iut = iu - 1
                 t1l = t1l + b1*C(kr - 1, iut)
                 el1 = elr(kr)
                 kr = kr + 1
                 if(kr>KB(iu) + 1)exit
             enddo
             iut = iu
             if(qu(k - 1, UHS(jb))>=0.0)iut = iu - 1
             b1 = B(kr - 1, iu)*(el1 - ell(k))
             brtot = brtot + b1
             t1l = (t1l + b1*C(kr - 1, iut))/brtot
             Ss(k - 1, UHS(jb)) = TSS(k - 1, UHS(jb)) - t1l*qu(k - 1, UHS(jb))
             if(kr>KB(iu) + 1)exit
             el1 = ell(k)
         endif
     enddo
     return
 
!***********************************************************************************************************************************
!**  D O W N S T R E A M   C O N S T I T U E N T                              
!***********************************************************************************************************************************
 
!    **
     entry DOWNSTREAM_CONSTITUENT(C, Ss)
     do jjb = 1, nbr
         if(CDHS(jb)>=US(jjb) .AND. CDHS(jb)<=DS(jjb))exit
     enddo
     do jjw = 1, nwb
         if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
     enddo
     do k = KTWB(jjw), KB(CDHS(jb)) + 1
         elr(k) = EL(k, CDHS(jb))
     enddo
     do k = kt, KB(id) + 1
         ell(k) = EL(k, id) + SINA(jb)*DLX(id)*0.5
     enddo
     elw = ELWS(id)                                   !EL(KTWB(JW),ID)-Z(ID)*COSA(JB)
     el1 = elw + SINA(jb)*DLX(id)*0.5
!    IF(SLOPE(JB) /= 0.0)THEN
!    EL1  = ELW+(ELW-ELWS(ID-1))/(0.5*(DLX(ID)+DLX(ID-1))*DLX(ID)*0.5   ! SW
!    7/17/09 ELSE
!    EL1=ELW
!    ENDIF
     kl = kt + 1
     do k = KTWB(jjw) + 1, KB(CDHS(jb)) + 1
         idt = id + 1
         if(qd(k - 1, CDHS(jb))>=0.0)idt = id
         if(elr(k)>=ell(kl))then
             t1l = C(kl - 1, idt)
             Ss(k - 1, CDHS(jb)) = Ss(k - 1, CDHS(jb))                         &
                                 & + t1l*qd(k - 1, CDHS(jb))
             el1 = elr(k)
             if(elr(k)==ell(kl))kl = kl + 1
         else
             t1l = 0.0
             brtot = 0.0
             do while (elr(k)<=ell(kl))
                 if(kl==KTWB(jw) + 1 .AND. k==KTWB(jjw) + 1)then
                     b1 = BH2(kt, id)
                     brtot = brtot + b1
                 else
                     b1 = B(kl - 1, id)*(el1 - ell(kl))
                     brtot = brtot + b1
                 endif
                 t1l = t1l + b1*C(kl - 1, idt)
                 el1 = ell(kl)
                 kl = kl + 1
                 if(kl>KB(id))then
                     b1 = B(kl - 1, id)*(el1 - elr(k))
                     brtot = brtot + b1
                     t1l = (t1l + b1*C(kl - 1, idt))/brtot
                     Ss(k - 1, CDHS(jb)) = Ss(k - 1, CDHS(jb))                 &
                       & + t1l*qd(k - 1, CDHS(jb))
                     goto 100
                 endif
             enddo
             b1 = B(kl - 1, id)*(el1 - elr(k))
             brtot = brtot + b1
             t1l = (t1l + b1*C(kl - 1, idt))/brtot
             Ss(k - 1, CDHS(jb)) = Ss(k - 1, CDHS(jb))                         &
                                 & + t1l*qd(k - 1, CDHS(jb))
             el1 = ell(kl)
         endif
     enddo
100   return
     entry DEALLOCATE_WATERBODY
     deallocate(ell, elr, cl, qu, qd)
     end subroutine WATERBODY
