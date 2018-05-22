!*==initgeom.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
!***********************************************************************************************************************************
!**                                                    Task 1.1.3: Geometry                                                       **
!***********************************************************************************************************************************
 
! Layer elevations
     subroutine INITGEOM
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
     use POROSITYC
     use MACROPHYTEC
     use RSTART
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: b11, el1, el2, ell, ell1, ell2, elr, elr2, err1, err2
     integer :: iexit, jjw, k1, kdn, ktmax, kup, ncbp, ninternal, nnbp, nup
!
!*** End of declarations rewritten by SPAG
!
 
 
     npoint = 0
     do jw = 1, nwb
         if(ZERO_SLOPE(jw))then
             do i = US(BS(jw)) - 1, DS(BE(jw)) + 1
                 EL(kmx, i) = ELBOT(jw)
                 do k = kmx - 1, 1, -1
                     EL(k, i) = EL(k + 1, i) + H(k, jw)
                 enddo
             enddo
         else
             EL(kmx, DS(JBDN(jw)) + 1) = ELBOT(jw)
             jb = JBDN(jw)
             npoint(jb) = 1
             nnbp = 1
             ncbp = 0
             ninternal = 0
             nup = 0
             do while (nnbp<=(BE(jw) - BS(jw) + 1))
                 ncbp = ncbp + 1
                 if(ninternal==0)then
                     if(nup==0)then
                         do i = DS(jb), US(jb), -1
                             if(i/=DS(jb))then
                                 EL(kmx, i) = EL(kmx, i + 1) + SINA(jb)        &
                                   & *(DLX(i) + DLX(i + 1))*0.5
                             else
                                 EL(kmx, i) = EL(kmx, i + 1)
                             endif
                             do k = kmx - 1, 1, -1
                                 EL(k, i) = EL(k + 1, i) + H(k, jw)*COSA(jb)
                             enddo
                         enddo
                     else
                         do i = US(jb), DS(jb)
                             if(i/=US(jb))then
                                 EL(kmx, i) = EL(kmx, i - 1) - SINA(jb)        &
                                   & *(DLX(i) + DLX(i - 1))*0.5
                             else
                                 EL(kmx, i) = EL(kmx, i - 1)
                             endif
                             do k = kmx - 1, 1, -1
                                 EL(k, i) = EL(k + 1, i) + H(k, jw)*COSA(jb)
                             enddo
                         enddo
                         nup = 0
                     endif
                     do k = kmx, 1, -1
                         if(UP_HEAD(jb))then
                             EL(k, US(jb) - 1) = EL(k, US(jb)) + SINA(jb)      &
                               & *DLX(US(jb))
                         else
                             EL(k, US(jb) - 1) = EL(k, US(jb))
                         endif
                         if(DN_HEAD(jb))then
                             EL(k, DS(jb) + 1) = EL(k, DS(jb)) - SINA(jb)      &
                               & *DLX(DS(jb))
                         else
                             EL(k, DS(jb) + 1) = EL(k, DS(jb))
                         endif
                     enddo
                 else
                     do k = kmx - 1, 1, -1
                         EL(k, UHS(jjb)) = EL(k + 1, UHS(jjb)) + H(k, jw)      &
                           & *COSA(jb)
                     enddo
                     do i = UHS(jjb) + 1, DS(jb)
                         EL(kmx, i) = EL(kmx, i - 1) - SINA(jb)                &
                                    & *(DLX(i) + DLX(i - 1))*0.5
                         do k = kmx - 1, 1, -1
                             EL(k, i) = EL(k + 1, i) + H(k, jw)*COSA(jb)
                         enddo
                     enddo
                     do i = UHS(jjb) - 1, US(jb), -1
                         EL(kmx, i) = EL(kmx, i + 1) + SINA(jb)                &
                                    & *(DLX(i) + DLX(i + 1))*0.5
                         do k = kmx - 1, 1, -1
                             EL(k, i) = EL(k + 1, i) + H(k, jw)*COSA(jb)
                         enddo
                     enddo
                     ninternal = 0
                 endif
                 if(nnbp==(BE(jw) - BS(jw) + 1))exit
 
!******          Find next branch connected to furthest downstream branch
 
                 do jb = BS(jw), BE(jw)
                     if(npoint(jb)/=1)then
                         do jjb = BS(jw), BE(jw)
                             if(DHS(jb)>=US(jjb) .AND. DHS(jb)<=DS(jjb) .AND.  &
                              & npoint(jjb)==1)then
                                 npoint(jb) = 1
                                 EL(kmx, DS(jb) + 1) = EL(kmx, DHS(jb))        &
                                   & + SINA(jb)*(DLX(DS(jb)) + DLX(DHS(jb)))   &
                                   & *0.5
                                 nnbp = nnbp + 1
                                 exit
                             endif
                             if(UHS(jjb)==DS(jb) .AND. npoint(jjb)==1)then
                                 npoint(jb) = 1
                                 EL(kmx, DS(jb) + 1) = EL(kmx, US(jjb))        &
                                   & + (SINA(jjb)*DLX(US(jjb)) + SINA(jb)      &
                                   & *DLX(DS(jb)))*0.5
                                 nnbp = nnbp + 1
                                 exit
                             endif
                             if(UHS(jjb)>=US(jb) .AND. UHS(jjb)<=DS(jb) .AND.  &
                              & npoint(jjb)==1)then
                                 npoint(jb) = 1
                                 EL(kmx, UHS(jjb)) = EL(kmx, US(jjb))          &
                                   & + SINA(jjb)*DLX(US(jjb))*0.5
                                 nnbp = nnbp + 1
                                 ninternal = 1
                                 exit
                             endif
                             if(UHS(jb)>=US(jjb) .AND. UHS(jb)<=DS(jjb) .AND.  &
                              & npoint(jjb)==1)then
                                 npoint(jb) = 1
                                 EL(kmx, US(jb) - 1) = EL(kmx, UHS(jb))        &
                                   & - SINA(jb)*DLX(US(jb))*0.5
                                 nnbp = nnbp + 1
                                 nup = 1
                                 exit
                             endif
                         enddo
                         if(npoint(jb)==1)exit
                     endif
                 enddo
             enddo
         endif
     enddo
 
!    Minimum/maximum layer heights
 
     do jw = 1, nwb
         do k = kmx - 1, 1, -1
             hmin = DMIN1(H(k, jw), hmin)
             hmax = DMAX1(H(k, jw), hmax)
         enddo
     enddo
     hmax2 = hmax**2
 
!    Water surface and bottom layers
 
     do jw = 1, nwb
         do jb = BS(jw), BE(jw)
             do i = US(jb) - 1, DS(jb) + 1
                 if(.NOT.restart_in)then
                     KTI(i) = 2
                     do while (EL(KTI(i), i)>ELWS(i))
                         KTI(i) = KTI(i) + 1
                         if(KTI(i)==kmx)then
                                   ! cb 7/7/2010 if elws below grid, setting to elws to just within grid ! SW 5/27/17 DELETED
                             KTI(i) = kmx - 1
                            ! 2
                             ELWS(i) = EL(KTI(i), i)
                             exit
                         endif
                     enddo
 
                     Z(i) = (EL(KTI(i), i) - ELWS(i))/COSA(jb)
 
 
                     ZMIN(jw) = DMAX1(ZMIN(jw), Z(i))
          !KTI(I)   =  MAX(KTI(I)-1,2)   ! MOVED SW 5/27/17
                     ktmax = MAX(2, KTI(i))
                     KTWB(jw) = MAX(ktmax, KTWB(jw))
                     KTI(i) = MAX(KTI(i) - 1, 2)
                                         ! original    IF(KTI(I) /= KMX)
                     if(Z(i)>ZMIN(jw))IZMIN(jw) = i
                 endif
                 k = 2
                 do while (B(k, i)>0.0)
                     KB(i) = k
                     k = k + 1
                 enddo
                 KBMAX(jw) = MAX(KBMAX(jw), KB(i))
             enddo
             KB(US(jb) - 1) = KB(US(jb))
             KB(DS(jb) + 1) = KB(DS(jb))
         enddo
 
 
!**      Correct for water surface going over several layers
 
         if(.NOT.restart_in)then
             kt = KTWB(jw)
             do jb = BS(jw), BE(jw)
                 do i = US(jb) - 1, DS(jb) + 1
                     H2(kt, i) = H(kt, jw) - Z(i)
                     k = KTI(i) + 1
                     do while (kt>k)
                         Z(i) = Z(i) - H(k, jw)
                         H2(kt, i) = H(kt, jw) - Z(i)
                         k = k + 1
                     enddo
                 enddo
             enddo
         endif
         ELKT(jw) = EL(KTWB(jw), DS(BE(jw))) - Z(DS(BE(jw)))*COSA(BE(jw))
     enddo
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = US(jb)
             id = DS(jb)
 
!****        Boundary bottom layers
 
             if(UH_EXTERNAL(jb))KB(iu - 1) = KB(iu)                                 ! CB 6/12/07
             if(DH_EXTERNAL(jb))KB(id + 1) = KB(id)
 
!****        Branch numbers corresponding to tributaries, withdrawals, and head
 
             if(tributaries)then
                 do jt = 1, ntr
                     if(ITR(jt)>=US(jb) .AND. ITR(jt)<=DS(jb))JBTR(jt) = jb
                 enddo
             endif
             if(withdrawals)then
                 do jwd = 1, nwd
                     if(IWD(jwd)>=US(jb) .AND. IWD(jwd)<=DS(jb))JBWD(jwd) = jb
                 enddo
             endif
             if(UH_INTERNAL(jb))then
                 do jjb = 1, nbr
                     JBUH(jb) = jjb
                     if(UHS(jb)>=US(jjb) .AND. UHS(jb)<=DS(jjb))exit
                 enddo
                 do jjw = 1, nwb
                     JWUH(jb) = jjw
                     if(JBUH(jb)>=BS(jjw) .AND. JBUH(jb)<=BE(jjw))exit
                 enddo
             endif
             if(INTERNAL_FLOW(jb))then
                 do jjb = 1, nbr
                     JBUH(jb) = jjb
                     if(UHS(jb)>=US(jjb) .AND. UHS(jb)<=DS(jjb))exit
                 enddo
                 do jjw = 1, nwb
                     JWUH(jb) = jjw
                     if(JBUH(jb)>=BS(jjw) .AND. JBUH(jb)<=BE(jjw))exit
                 enddo
             endif
             if(DH_INTERNAL(jb))then
                 do jjb = 1, nbr
                     JBDH(jb) = jjb
                     if(DHS(jb)>=US(jjb) .AND. DHS(jb)<=DS(jjb))exit
                 enddo
                 do jjw = 1, nwb
                     JWDH(jb) = jjw
                     if(JBDH(jb)>=BS(jjw) .AND. JBDH(jb)<=BE(jjw))exit
                 enddo
             endif
 
!****        Bottom boundary cells
 
             if(UH_INTERNAL(jb))then
                 if(JBUH(jb)>=BS(jw) .AND. JBUH(jb)<=BE(jw))then
                     KB(iu - 1) = MIN(KB(UHS(jb)), KB(iu))                   ! CB 6/12/07
                 elseif(EL(KB(iu), iu)>=EL(KB(UHS(jb)), UHS(jb)))then        ! CB 6/12/07
                     KB(iu - 1) = KB(iu)                                     ! CB 6/12/07
                 else
                     do k = kt, KB(iu)                                       ! CB 6/12/07
                         if(EL(KB(UHS(jb)), UHS(jb))>=EL(k, iu))then         ! CB 6/12/07
                             KB(iu - 1) = k
                             exit                                            ! CB 6/12/07
                         endif
                     enddo
                 endif
             endif
             if(DH_INTERNAL(jb))then
                 if(JBDH(jb)>=BS(jw) .AND. JBDH(jb)<=BE(jw))then
                     KB(id + 1) = MIN(KB(DHS(jb)), KB(id))
                 elseif(EL(KB(id), id)>=EL(KB(DHS(jb)), DHS(jb)))then
                     KB(id + 1) = KB(id)
                 else
                     do k = kt, KB(id)
                         if(EL(KB(DHS(jb)), DHS(jb))>=EL(k, id))then
                             KB(id + 1) = k
                             exit
                         endif
                     enddo
                 endif
             endif
 
!****        Boundary segment lengths
 
             DLX(iu - 1) = DLX(iu)
             DLX(id + 1) = DLX(id)
 
!****        Minimum bottom layers and average segment lengths
 
             do i = iu - 1, id
                 KBMIN(i) = MIN(KB(i), KB(i + 1))
                 DLXR(i) = (DLX(i) + DLX(i + 1))*0.5
             enddo
             KBMIN(id + 1) = KBMIN(id)
             DLXR(id + 1) = DLX(id)
 
!****        Minimum/maximum segment lengths
 
             do i = iu, id
                 dlxmin = DMIN1(dlxmin, DLX(i))
                 dlxmax = DMAX1(dlxmax, DLX(i))
             enddo
         enddo
     enddo
 
!    Boundary widths
 
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
         enddo
             ! SW 1/23/06
     enddo   ! SW 1/23/06
     bnew = B
             ! SW 1/23/06
     kbi = KB
             ! SW 10/29/2010
 
!**** Upstream active segment and single layer  ! 1/23/06 entire section moved SW
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = US(jb)
             id = DS(jb)
             iexit = 0
             if(SLOPE(jb)/=0.0)then
                 do i = US(jb) - 1, DS(jb) + 1
                     if(kbi(i)<kt)then
                                     ! SW 1/23/06
                         do k = kbi(i) + 1, kt
                             bnew(k, i) = 0.000001
                                      ! SW 1/23/06
                         enddo
                         KB(i) = kt
                     endif
                 enddo
             endif
             iut = iu
             do i = iu, id
                 if(KB(i) - kt<NL(jb) - 1)iut = i + 1
                 ONE_LAYER(i) = kt==KB(i)
             enddo
             do i = iu - 1, id
                     ! recompute kbmin after one_layer computation   11/12/07
                 KBMIN(i) = MIN(KB(i), KB(i + 1))
             enddo
             KBMIN(id + 1) = KBMIN(id)
 
             CUS(jb) = iut
             if(iut>=DS(jb))BR_INACTIVE(jb) = .TRUE.
                                                ! SW 6/12/2017
 
 
!****        Areas and bottom widths
 
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
 
!                column widths
                 do i = iu - 1, id + 1
                     CW(KB(i), i) = B(KB(i), i)
                     do k = 1, KB(i) - 1
                         CW(k, i) = B(k, i) - B(k + 1, i)
                     enddo
                 enddo
 
!******          Derived geometry
 
                 do i = iu - 1, id + 1
                     BH2(kt, i) = B(KTI(i), i)                                 &
                                & *(EL(kt, i) - EL(KTI(i) + 1, i) - Z(i)       &
                                & *COSA(jb))/COSA(jb)
                     if(kt==KTI(i))BH2(kt, i) = H2(kt, i)*B(kt, i)
                     do k = KTI(i) + 1, kt
                         BH2(kt, i) = BH2(kt, i) + bnew(k, i)*H(k, jw)                   ! sw 1/23/06    BH(K,I)
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
     h1 = H2
     bh1 = BH2
     bhr1 = BHR2
     avh1 = AVH2
 
!    Temporary downstream head segment
 
     do jb = 1, nbr
         if(DHS(jb)>0)then
             do jjb = 1, nbr
                 if(DHS(jb)>=US(jjb) .AND. DHS(jb)<=DS(jjb))exit
             enddo
             if(CUS(jjb)>DHS(jb))CDHS(jb) = CUS(jjb)
         endif
     enddo
 
!    Total active cells
 
     do jw = 1, nwb
         do jb = BS(jw), BE(jw)
             if(BR_INACTIVE(jb))cycle
                                 ! SW 6/12/2017
             do i = CUS(jb), DS(jb)
                 do k = KTWB(jw), KB(i)
                     ntac = ntac + 1
                 enddo
             enddo
             ntacmx = ntac
             ntacmn = ntac
 
!****        Wind fetch lengths
 
             do i = US(jb), DS(jb)
                 FETCHD(i, jb) = FETCHD(i - 1, jb) + DLX(i)
             enddo
             do i = DS(jb), US(jb), -1
                 FETCHU(i, jb) = FETCHU(i + 1, jb) + DLX(i)
             enddo
         enddo
     enddo
 
!    Segment heights
 
     do jw = 1, nwb
         do jb = BS(jw), BE(jw)
             do i = US(jb) - 1, DS(jb) + 1
                 do k = MIN(kmx - 1, KB(i)), 2, -1
                     HSEG(k, i) = HSEG(k + 1, i) + H2(k, i)
                 enddo
             enddo
         enddo
     enddo
 
!    Beginning and ending segments/layers for snapshots
 
     do jw = 1, nwb
         do i = 1, NISNP(jw)
             KBR(jw) = MAX(KB(ISNP(i, jw)), KBR(jw))
         enddo
     enddo
 
     end subroutine INITGEOM
