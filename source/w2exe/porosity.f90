!*==porosity.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
!************************************************************************
!**               S U B R O U T I N E    porosity                      **
!************************************************************************
 
     subroutine POROSITY
 
 
 
 
 
     use GEOMC
     use GLOBAL
     use MACROPHYTEC
     use POROSITYC
     use LOGICC
     use SCREENC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: b11, el1, el2, ell, ell1, ell2, elr, elr2, err1, err2, vstot
     integer :: iexit, iut, k, k1, kdn, kup, m
!
!*** End of declarations rewritten by SPAG
!
 
 
     if(nit==0)then
 
9001     format((8x, i8, 3F8.0))
 
         do jw = 1, nwb
             kt = KTWB(jw)
             do jb = BS(jw), BE(jw)
                 iu = CUS(jb)
                 id = DS(jb)
                 do i = iu, id
                     do k = 2, KB(i)
                         VOLI(k, i) = BH(k, i)*DLX(i)
                     enddo
                     VOLI(kt, i) = BH2(kt, i)*DLX(i)
                 enddo
             enddo
         enddo
 
     endif
 
     do jb = 1, nbr
         COSA(jb) = COS(ALPHA(jb))
     enddo
 
!C   CALCULATING # OF MACROPHYTE STEMS IN EACH CELL
 
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do i = iu, id
!                HKTI  = H(KT,JW)-Z(I)    REPLACED BY H1(KT,I)
                 if(kt==KTI(i))then
                     VOLKTI(i) = H1(kt, i)*BIC(kt, i)*DLX(i)
                 else
                     VOLKTI(i) = BIC(KTI(i), i)                                &
                               & *(EL(kt, i) - EL(KTI(i) + 1, i) - Z(i)        &
                               & *COSA(jb))/COSA(jb)*DLX(i)
                 endif
                 do k = KTI(i) + 1, kt
                     VOLKTI(i) = VOLKTI(i) + VOLI(k, i)
                 enddo
 
                 do m = 1, nmc
                     VSTEMKT(i, m) = (MAC(kt, i, m)*VOLKTI(i))/DWV(m)
                                                           !CB 6/29/06
                 enddo
 
                 do k = kt + 1, KB(i)
                     do m = 1, nmc
                         VSTEM(k, i, m) = (MAC(k, i, m)*VOLI(k, i))/DWV(m)
                                                           !CB 6/29/06
                     enddo
                 enddo
             enddo
         enddo
     enddo
 
     por = 1.0
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
 
             do i = iu, id
                 do k = kt, KB(i)
                     if(k==kt)then
                         vstot = 0.0
                         do m = 1, nmc
                             vstot = vstot + VSTEMKT(i, m)
                         enddo
                         por(kt, i) = (VOLKTI(i) - vstot)/VOLKTI(i)
                     else
                         vstot = 0.0
                         do m = 1, nmc
                             vstot = vstot + VSTEM(k, i, m)
                         enddo
                         por(k, i) = (VOLI(k, i) - vstot)/VOLI(k, i)
                     endif
                 enddo
             enddo
 
             do i = iu, id
                 do k = KTI(i), KB(i)
                     if(k<=kt)then
                         B(k, i) = por(kt, i)*BIC(k, i)
                     else
                         B(k, i) = por(k, i)*BIC(k, i)
                     endif
 
                 enddo
             enddo
 
         enddo
     enddo
 
 
 
!    BOUNDARY WIDTHS
 
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = US(jb)
             id = DS(jb)
             do i = iu - 1, id + 1
                 B(1, i) = B(2, i)
                 do k = KB(i) + 1, kmx
                     B(k, i) = B(KB(i), i)
                 enddo
             enddo
         enddo
     enddo
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = US(jb)
             id = DS(jb)
             iexit = 0
             do k = 1, kmx - 1
                 B(k, iu - 1) = B(k, iu)
                 if(UH_INTERNAL(jb) .OR. HEAD_FLOW(jb))then
                     if(JBUH(jb)>=BS(jw) .AND. JBUH(jb)<=BE(jw))then
                         B(k, iu - 1) = B(k, UHS(jb))
                     else
                         elr = EL(k, iu) + SINA(jb)*DLX(iu)*0.5
                         ell = EL(2, UHS(jb)) - SINA(JBUH(jb))*DLX(UHS(jb))*0.5
                         if(elr>=ell)then
                             B(k, iu - 1) = B(2, UHS(jb))
                         else
                             do kup = 2, kmx - 1
                                 ell1 = EL(kup, UHS(jb)) - SINA(JBUH(jb))      &
                                      & *DLX(UHS(jb))*0.5
                                 ell2 = EL(kup + 1, UHS(jb)) - SINA(JBUH(jb))  &
                                      & *DLX(UHS(jb))*0.5
                                 if(ell1>elr .AND. ell2<=elr)then
                                     if(kup>KB(UHS(jb)))then
                                         KB(iu - 1) = k - 1
                                         KBMIN(iu - 1)                         &
                                           & = MIN(KB(iu), KB(iu - 1))
                                         iexit = 1
                                         exit
                                     endif
                                     elr2 = EL(k + 1, iu) + SINA(jb)*DLX(iu)   &
                                       & *0.5
                                     if(elr2>=ell2)then
                                         B(k, iu - 1) = B(kup, UHS(jb))
                                     else
                                         k1 = kup + 1
                                         if(k1>kmx)exit
                                         b11 = 0.0
                                         el1 = elr
                                         el2 = EL(k1, UHS(jb)) - SINA(JBUH(jb))&
                                           & *DLX(iu)*0.5
                                         do while (elr2<=el2)
                                         b11 = b11 + (el1 - el2)               &
                                           & *B(k1 - 1, UHS(jb))
                                         el1 = el2
                                         k1 = k1 + 1
                                         if(k1>=kmx + 1 .OR. el2==elr2)exit
                                         el2 = EL(k1, UHS(jb)) - SINA(JBUH(jb))&
                                           & *DLX(UHS(jb))*0.5
                                         if(el2<=elr2)el2 = elr2
                                         enddo
                                         B(k, iu - 1) = b11/H(k, jw)
                                     endif
                                     exit
                                 endif
                             enddo
                             if(EL(kmx, UHS(jb))>EL(k, iu))B(k, iu - 1)        &
                              & = B(k - 1, iu - 1)
                             if(B(k, iu - 1)==0.0)B(k, iu - 1)                 &
                              & = B(k - 1, iu - 1)
                             if(iexit==1)exit
                         endif
                     endif
                 endif
             enddo
             iexit = 0
             do k = 1, kmx - 1
                 B(k, id + 1) = B(k, id)
                 if(DH_INTERNAL(jb))then
                     if(JBDH(jb)>=BS(jw) .AND. JBDH(jb)<=BE(jw))then
                         B(k, id + 1) = B(k, DHS(jb))
                     else
                         ell = EL(k, id) - SINA(jb)*DLX(id)*0.5
                         elr = EL(2, DHS(jb)) + SINA(JBDH(jb))*DLX(DHS(jb))*0.5
                         if(ell>=elr)then
                             B(k, id + 1) = B(2, DHS(jb))
                         else
                             do kdn = 2, kmx - 1
                                 err1 = EL(kdn, DHS(jb)) + SINA(JBDH(jb))      &
                                      & *DLX(DHS(jb))*0.5
                                 err2 = EL(kdn + 1, DHS(jb)) + SINA(JBDH(jb))  &
                                      & *DLX(DHS(jb))*0.5
                                 if(err1>=ell .AND. err2<ell)then
                                     if(kdn>KB(DHS(jb)))then
                                         KB(id + 1) = k - 1
                                         KBMIN(id) = MIN(KB(id), KB(id + 1))
                                         iexit = 1
                                         exit
                                     endif
                                     ell2 = EL(k + 1, id) - SINA(jb)*DLX(id)   &
                                       & *0.5
                                     if(ell2>=err2)then
                                         B(k, id + 1) = B(kdn, DHS(jb))
                                     else
                                         k1 = kdn + 1
                                         if(k1>kmx)exit
                                         b11 = 0.0
                                         el2 = ell
                                         el1 = EL(k1, DHS(jb)) + SINA(JBDH(jb))&
                                           & *DLX(DHS(jb))*0.5
                                         do while (ell2<=el1)
                                         b11 = b11 + (el2 - el1)               &
                                           & *B(k1 - 1, DHS(jb))
                                         el2 = el1
                                         k1 = k1 + 1
                                         if(k1>=kmx + 1 .OR. el1==ell2)exit
                                         el1 = EL(k1, DHS(jb)) + SINA(JBDH(jb))&
                                           & *DLX(DHS(jb))*0.5
                                         if(el1<=ell2)el1 = ell2
                                         enddo
                                         B(k, id + 1) = b11/H(k, jw)
                                     endif
                                     exit
                                 endif
                             enddo
                             if(EL(kmx, DHS(jb))>EL(k, id))B(k, id + 1)        &
                              & = B(k - 1, id + 1)
                             if(B(k, id + 1)==0.0)B(k, id + 1)                 &
                              & = B(k - 1, id + 1)
                             if(iexit==1)exit
                         endif
                     endif
                 endif
             enddo
 
!****        AREAS AND BOTTOM WIDTHS
 
             if(.NOT.TRAPEZOIDAL(jw))then                                                                              !SW 07/16/04
                 do i = iu - 1, id + 1
                     do k = 1, kmx - 1
                         BH2(k, i) = B(k, i)*H(k, jw)
                         BH(k, i) = B(k, i)*H(k, jw)
                         BB(k, i) = B(k, i) - (B(k, i) - B(k + 1, i))          &
                                  & /(0.5*(H(k, jw) + H(k + 1, jw)))*H(k, jw)  &
                                  & *0.5                                                                               !SW 08/02/04
                     enddo
                     BH(kmx, i) = BH(kmx - 1, i)
                 enddo
!******          DERIVED GEOMETRY
 
                 do i = iu - 1, id + 1
                     BH2(kt, i) = B(KTI(i), i)                                 &
                                & *(EL(kt, i) - EL(KTI(i) + 1, i) - Z(i)       &
                                & *COSA(jb))/COSA(jb)
                     if(kt==KTI(i))BH2(kt, i) = H2(kt, i)*B(kt, i)
                     do k = KTI(i) + 1, kt
                         BH2(kt, i) = BH2(kt, i) + BH(k, i)
                     enddo
                     BKT(i) = BH2(kt, i)/H2(kt, i)
                     BI(kt, i) = B(KTI(i), i)
                 enddo
             else                                                                                                      !SW 07/16/04
                 do i = iu - 1, id + 1
                     do k = 1, kmx - 1
                         BB(k, i) = B(k, i) - (B(k, i) - B(k + 1, i))          &
                                  & /(0.5*(H(k, jw) + H(k + 1, jw)))*H(k, jw)  &
                                  & *0.5
                     enddo
                     BB(KB(i), i) = B(KB(i), i)*0.5
                     BH2(1, i) = B(1, i)*H(1, jw)
                     BH(1, i) = BH2(1, i)
                     do k = 2, kmx - 1
                         BH2(k, i) = 0.25*H(k, jw)                             &
                                   & *(BB(k - 1, i) + 2.*B(k, i) + BB(k, i))
                         BH(k, i) = BH2(k, i)
                     enddo
                     BH(kmx, i) = BH(kmx - 1, i)
                 enddo
                 do i = iu - 1, id + 1
                     call GRID_AREA1(EL(kt, i) - Z(i), EL(kt + 1, i),          &
                                   & BH2(kt, i), BI(kt, i))
                     BKT(i) = BH2(kt, i)/H2(kt, i)
                 enddo
             endif
             do i = iu - 1, id
                 do k = 1, kmx - 1
                     AVH2(k, i) = (H2(k, i) + H2(k + 1, i))*0.5
                     AVHR(k, i) = H2(k, i) + (H2(k, i + 1) - H2(k, i))         &
                                & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                              !SW 07/29/04
                 enddo
                 AVH2(kmx, i) = H2(kmx, i)
                 do k = 1, kmx
                     BR(k, i) = B(k, i) + (B(k, i + 1) - B(k, i))              &
                              & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                                !SW 07/29/04
                     BHR(k, i) = BH(k, i) + (BH(k, i + 1) - BH(k, i))          &
                               & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                               !SW 07/29/04
                     BHR2(k, i) = BH2(k, i) + (BH2(k, i + 1) - BH2(k, i))      &
                                & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                              !SW 07/29/04
                 enddo
             enddo
             do k = 1, kmx - 1
                 AVH2(k, id + 1) = (H2(k, id + 1) + H2(k + 1, id + 1))*0.5
                 BR(k, id + 1) = B(k, id + 1)
                 BHR(k, id + 1) = BH(k, id + 1)
             enddo
             AVH2(kmx, id + 1) = H2(kmx, id + 1)
             AVHR(kt, id + 1) = H2(kt, id + 1)
             BHR2(kt, id + 1) = BH2(kt, id + 1)
             iut = iu
             if(UP_HEAD(jb))iut = iu - 1
             do i = iut, id
                 do k = 1, kmx - 1
                     VOL(k, i) = B(k, i)*H2(k, i)*DLX(i)
                 enddo
                 VOL(kt, i) = BH2(kt, i)*DLX(i)
                 DEPTHB(kt, i) = H2(kt, i)
                 DEPTHM(kt, i) = H2(kt, i)*0.5
                 do k = kt + 1, kmx
                     DEPTHB(k, i) = DEPTHB(k - 1, i) + H2(k, i)
                     DEPTHM(k, i) = DEPTHM(k - 1, i)                           &
                                  & + (H2(k - 1, i) + H2(k, i))*0.5
                 enddo
             enddo
         enddo
     enddo
 
 
     end subroutine POROSITY
