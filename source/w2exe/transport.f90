!*==transport.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                           S U B R O U T I N E   T R A N S P O R T                                             **
!***********************************************************************************************************************************
 
     subroutine TRANSPORT
     use GLOBAL
     use GEOMC
     use TVDC
     use TRANS
     use LOGICC
     use STRUCTURES
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND) :: acurz, adelc, c1x, c1z, c2x, c2z, c3x, c3z, calf, cart, cmax1,&
               & cmin1, cour, cref, delc, dltdlxr, dlxmin, dlxt, flux, ftemp,  &
               & hb, hm, hmin, ht, ratdi, ratzi
     real(R8KIND), allocatable, dimension(:, :), save :: ad1x, ad1z, ad2x, ad2z,   &
           & ad3x, ad3z, alfa, alfaz, curs1, curs1z, curs2, curs2z, curs3,     &
           & curs3z, dx1, dx2, dx3, rats, ratsz, sf1x, sf1z
     real(R8KIND), allocatable, dimension(:), save :: curx1, curx2, curx3, ratd
     integer :: k
     real(R8KIND), allocatable, dimension(:, :, :), save :: sf10x, sf10z, sf11x,   &
           & sf12x, sf13x, sf2x, sf2z, sf3x, sf3z, sf4x, sf4z, sf5x, sf5z,     &
           & sf6x, sf6z, sf7x, sf7z, sf8x, sf8z, sf9x, sf9z
!
!*** End of declarations rewritten by SPAG
!
 
!    Type declarations
 
!    REAL,     SAVE, ALLOCATABLE, DIMENSION(:,:)   :: RATZ,   CURZ1,  CURZ2, 
!    CURZ3
 
 
!    Allocation declarations
 
     allocate(ratd(imx), curx1(imx), curx2(imx), curx3(imx))
     allocate(sf1x(kmx, imx), sf1z(kmx, nwb))
!    ALLOCATE (RATZ(KMX,NWB),   CURZ1(KMX,NWB),  CURZ2(KMX,NWB),  
!    CURZ3(KMX,NWB))
     allocate(dx1(kmx, imx), dx2(kmx, imx), dx3(kmx, imx))
     allocate(ad1x(kmx, imx), ad2x(kmx, imx), ad3x(kmx, imx))
     allocate(ad1z(kmx, imx), ad2z(kmx, imx), ad3z(kmx, imx))
     allocate(rats(kmx, imx), curs1(kmx, imx), curs2(kmx, imx), curs3(kmx, imx)&
            & , alfa(kmx, imx))
     allocate(ratsz(kmx, imx), curs1z(kmx, imx), curs2z(kmx, imx),             &
            & curs3z(kmx, imx), alfaz(kmx, imx))
     allocate(sf2x(kmx, imx, 2), sf3x(kmx, imx, 2), sf4x(kmx, imx, 2),         &
            & sf5x(kmx, imx, 2), sf6x(kmx, imx, 2), sf7x(kmx, imx, 2))
     allocate(sf8x(kmx, imx, 2), sf9x(kmx, imx, 2), sf10x(kmx, imx, 2),        &
            & sf11x(kmx, imx, 2), sf12x(kmx, imx, 2), sf13x(kmx, imx, 2))
     allocate(sf2z(kmx, 2, nwb), sf3z(kmx, 2, nwb), sf4z(kmx, 2, nwb),         &
            & sf5z(kmx, 2, nwb), sf6z(kmx, 2, nwb), sf7z(kmx, 2, nwb))
     allocate(sf8z(kmx, 2, nwb), sf9z(kmx, 2, nwb), sf10z(kmx, 2, nwb))
 
!    Variable initialization
 
     ct = 0.0
 
!    Variable initialization
 
     at = 0.0
 
!    Variable initialization
 
     vt = 0.0
 
!    Variable initialization
 
     dt = 0.0
 
!    Variable initialization
 
     dx1 = 0.0
 
!    Variable initialization
 
     dx2 = 0.0
 
!    Variable initialization
 
     dx3 = 0.0
 
!    Variable initialization
 
     adz = 0.0
 
!    Variable initialization
 
     adx = 0.0
 
!    Variable initialization
 
     ad1x = 0.0
     ad2x = 0.0
     ad3x = 0.0
     ad1z = 0.0
     ad2z = 0.0
     ad3z = 0.0
     return
 
!***********************************************************************************************************************************
!**  I N T E R P O L A T I O N  M U L T I P L I E R S                         
!***********************************************************************************************************************************
 
!    **
     entry INTERPOLATION_MULTIPLIERS
 
!    Positive horizontal flows
 
     do i = 2, imx - 1
         do k = 2, kmx - 1
             dlxt = DLX(i - 1)
             if(k>KB(i - 1) .OR. INTERNAL_WEIR(k, i))dlxt = DLX(i)
             dlxmin = DMIN1(DLX(i + 1), DLX(i))
             sf1x(k, i) = (DLX(i + 1) + DLX(i))*0.5D0
             sf2x(k, i, 1) = DLX(i)/(DLX(i) + DLX(i + 1))
      !SF3X(K,I,1)  =  DLX(I)**2
             sf3x(k, i, 1) = DLX(i)*DLX(i)   ! SW 4/20/16 SPEED
             sf4x(k, i, 1) = DLX(i + 1)/(DLX(i) + DLX(i + 1))
             sf5x(k, i, 1) = 0.25D0*(dlxt + 2.0D0*DLX(i) + DLX(i + 1))         &
                           & *(dlxt + DLX(i))
             sf6x(k, i, 1) = -0.25D0*(DLX(i) + DLX(i + 1))*(dlxt + DLX(i))
             sf7x(k, i, 1) = 0.25D0*(DLX(i) + DLX(i + 1))                      &
                           & *(dlxt + 2.0D0*DLX(i) + DLX(i + 1))
             sf8x(k, i, 1) = 0.50D0*(DLX(i) - DLX(i + 1))*dlxmin
             sf9x(k, i, 1) = 0.50D0*(dlxt + 2.0D0*DLX(i) - DLX(i + 1))*dlxmin
             sf10x(k, i, 1) = 0.50D0*(dlxt + 3.0D0*DLX(i))*dlxmin
      !SF11X(K,I,1) =  SF8X(K,I,1) /SF5X(K,I,1)/SF1X(K,I)
      !SF12X(K,I,1) =  SF9X(K,I,1) /SF6X(K,I,1)/SF1X(K,I)
      !SF13X(K,I,1) =  SF10X(K,I,1)/SF7X(K,I,1)/SF1X(K,I)
             sf11x(k, i, 1) = sf8x(k, i, 1)/(sf5x(k, i, 1)*sf1x(k, i)) ! SW 4/20/16 SPEED
             sf12x(k, i, 1) = sf9x(k, i, 1)/(sf6x(k, i, 1)*sf1x(k, i))
             sf13x(k, i, 1) = sf10x(k, i, 1)/(sf7x(k, i, 1)*sf1x(k, i))
         enddo
     enddo
 
!    Negative horizontal flows
 
     do i = 2, imx - 2
         do k = 2, kmx - 1
             dlxt = DLX(i + 2)
             if(k>KB(i + 2))dlxt = DLX(i + 1)
             dlxmin = DMIN1(DLX(i), DLX(i + 1))
             sf1x(k, i) = (DLX(i + 1) + DLX(i))*0.5D0
             sf2x(k, i, 2) = DLX(i + 1)/(DLX(i) + DLX(i + 1))
      !SF3X(K,I,2)  =  DLX(I+1)**2
             sf3x(k, i, 2) = DLX(i + 1)*DLX(i + 1) ! SW 4/20/16 SPEED
             sf4x(k, i, 2) = DLX(i)/(DLX(i) + DLX(i + 1))
             sf5x(k, i, 2) = 0.25D0*(DLX(i) + 2.0D0*DLX(i + 1) + dlxt)         &
                           & *(DLX(i) + DLX(i + 1))
             sf6x(k, i, 2) = -0.25D0*(DLX(i + 1) + dlxt)*(DLX(i) + DLX(i + 1))
             sf7x(k, i, 2) = 0.25D0*(DLX(i) + 2.0D0*DLX(i + 1) + dlxt)         &
                           & *(DLX(i + 1) + dlxt)
             sf8x(k, i, 2) = -0.50D0*(3.0D0*DLX(i + 1) + dlxt)*dlxmin
             sf9x(k, i, 2) = 0.50D0*(DLX(i) - 2.0D0*DLX(i + 1) - dlxt)*dlxmin
             sf10x(k, i, 2) = 0.50D0*(DLX(i) - DLX(i + 1))*dlxmin
      !SF11X(K,I,2) =  SF8X(K,I,2) /SF5X(K,I,2)/SF1X(K,I)
      !SF12X(K,I,2) =  SF9X(K,I,2) /SF6X(K,I,2)/SF1X(K,I)
      !SF13X(K,I,2) =  SF10X(K,I,2)/SF7X(K,I,2)/SF1X(K,I)
             sf11x(k, i, 2) = sf8x(k, i, 2)/(sf5x(k, i, 2)*sf1x(k, i))
                                                                   ! SW 4/20/16 SPEED
             sf12x(k, i, 2) = sf9x(k, i, 2)/(sf6x(k, i, 2)*sf1x(k, i))
             sf13x(k, i, 2) = sf10x(k, i, 2)/(sf7x(k, i, 2)*sf1x(k, i))
         enddo
     enddo
 
!    Ultimate multipliers
 
     if(ULTIMATE(jw))then
         do jb = BS(jw), BE(jw)
             do i = US(jb), DS(jb)
                            !CONCURRENT(I=US(JB):DS(JB))      !FORALL                                                        !DO I=US(JB),DS(JB)
                 ratd(i) = DLXR(i - 1)/DLXR(i)
        !CURX1(I) =  2.0D0*DLX(I)**2/(DLXR(I)+DLXR(I-1))/DLXR(I-1)   ! code speed improvement SW 4/20/16
        !CURX2(I) = -2.0D0*DLX(I)**2/(DLXR(I)*DLXR(I-1))
        !CURX3(I) =  2.0D0*DLX(I)**2/(DLXR(I)+DLXR(I-1))/DLXR(I)
                 curx1(i) = 2.0D0*DLX(i)*DLX(i)                                &
                          & /((DLXR(i) + DLXR(i - 1))*DLXR(i - 1))
                 curx2(i) = -2.0D0*DLX(i)*DLX(i)/(DLXR(i)*DLXR(i - 1))
                 curx3(i) = 2.0D0*DLX(i)*DLX(i)                                &
                          & /((DLXR(i) + DLXR(i - 1))*DLXR(i))
             enddo
         enddo
     endif
 
!    Vertical positive flows
 
     do k = 2, kmx - 1
         ht = H(k - 1, jw)
         hm = H(k, jw)
         hb = H(k + 1, jw)
         hmin = DMIN1(hb, hm)
         sf1z(k, jw) = (hb + hm)*0.5D0
    !SF2Z(K,1,JW)  =  HM**2
         sf2z(k, 1, jw) = hm*hm  ! SW 4/20/16
         sf3z(k, 1, jw) = hm/(hm + hb)
         sf4z(k, 1, jw) = hb/(hm + hb)
         sf5z(k, 1, jw) = 0.25D0*(ht + 2.0D0*hm + hb)*(ht + hm)
         sf6z(k, 1, jw) = -0.25D0*(hm + hb)*(ht + hm)
         sf7z(k, 1, jw) = 0.25D0*(hm + hb)*(ht + 2.0D0*hm + hb)
         sf8z(k, 1, jw) = 0.50D0*(hm - hb)*hmin
         sf9z(k, 1, jw) = 0.50D0*(ht + 2.0D0*hm - hb)*hmin
         sf10z(k, 1, jw) = 0.50D0*(ht + 3.0D0*hm)*hmin
     enddo
 
!    Vertical negative flows
 
     do k = 2, kmx - 2
         ht = H(k, jw)
         hm = H(k + 1, jw)
         hb = H(k + 2, jw)
         hmin = DMIN1(ht, hm)
         sf1z(k, jw) = (hm + ht)*0.5D0
    !SF2Z(K,2,JW)  =  HM**2
         sf2z(k, 2, jw) = hm*hm ! SW 4/20/16
         sf3z(k, 2, jw) = hm/(ht + hm)
         sf4z(k, 2, jw) = ht/(ht + hm)
         sf5z(k, 2, jw) = 0.25D0*(ht + 2.0D0*hm + hb)*(ht + hm)
         sf6z(k, 2, jw) = -0.25D0*(hm + hb)*(ht + hm)
         sf7z(k, 2, jw) = 0.25D0*(ht + 2.0D0*hm + hb)*(hm + hb)
         sf8z(k, 2, jw) = -0.50D0*(3.0D0*hm + hb)*hmin
         sf9z(k, 2, jw) = 0.50D0*(ht - 2.0D0*hm - hb)*hmin
         sf10z(k, 2, jw) = 0.50D0*(ht - hm)*hmin
     enddo
 
!    Ultimate multipliers
 
     if(ULTIMATE(jw))then   !also called during UPDATE since surface layer properties change
         do k = 2, kmx
             RATZ(k, jw) = AVH2(k - 1, DS(BE(jw)))/AVH2(k, DS(BE(jw)))                               ! SW 5/20/05
      !CURZ1(K,JW) =  2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))/AVH2(K-1,DS(BE(JW)))   ! SW 5/20/05
      !CURZ2(K,JW) = -2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))*AVH2(K,DS(BE(JW))))                        ! SW 5/20/05
      !CURZ3(K,JW) =  2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))/AVH2(K,DS(BE(JW)))     ! SW 5/20/05
             CURZ1(k, jw) = 2.0*H(k, jw)*H(k, jw)                              &
                          & /((AVH2(k - 1, DS(BE(jw))) + AVH2(k, DS(BE(jw))))  &
                          & *AVH2(k - 1, DS(BE(jw))))                                                      ! SW 5/20/05   SPEED 4/20/16
             CURZ2(k, jw) = -2.0*H(k, jw)*H(k, jw)                             &
                          & /(AVH2(k - 1, DS(BE(jw)))*AVH2(k, DS(BE(jw))))                                 ! SW 5/20/05
             CURZ3(k, jw) = 2.0*H(k, jw)*H(k, jw)                              &
                          & /((AVH2(k - 1, DS(BE(jw))) + AVH2(k, DS(BE(jw))))  &
                          & *AVH2(k, DS(BE(jw))))                                                          ! SW 5/20/05
         enddo
     endif
     return
 
!***********************************************************************************************************************************
!**  H O R I Z O N T A L  M U L T I P L I E R S                               
!***********************************************************************************************************************************
 
!    **
     entry HORIZONTAL_MULTIPLIERS1
 
!    Horizontal advection and diffusion multipliers FIRST PASS
 
     if(UPWIND(jw))then
         do i = iu, id - 1
             do k = kt, KB(i)
                 if(U(k, i)>=0.0)then
                     dx2(k, i) = -DX(k, i)/sf1x(k, i)
                     dx3(k, i) = DX(k, i)/sf1x(k, i)
                 else
                     dx1(k, i) = -DX(k, i)/sf1x(k, i)
                     dx2(k, i) = DX(k, i)/sf1x(k, i)
                 endif
             enddo
         enddo
     else
!!$OMP   PARALLEL DO
         do i = iu, id - 1
             dltdlxr = dlt/DLXR(i)
             do k = kt, KB(i)
                 cour = U(k, i)*dltdlxr
                 if(U(k, i)>=0.0)then
                     rats(k, i) = ratd(i)
                     curs1(k, i) = curx1(i)
                     curs2(k, i) = curx2(i)
                     curs3(k, i) = curx3(i)
                     dx1(k, i) = DX(k, i)*sf11x(k, i, 1)
                     dx2(k, i) = DX(k, i)*sf12x(k, i, 1)
                     dx3(k, i) = DX(k, i)*sf13x(k, i, 1)
                     alfa(k, i) = 2.0D0*(DX(k, i)*dlt/(sf1x(k, i)*sf1x(k, i))  &
                                & - (1.0D0 - cour*cour)*0.1666667D0)           &
                                & *sf3x(k, i, 1)                                                                    !/6.0
                     ad1x(k, i) = (alfa(k, i) - cour*sf8x(k, i, 1)*0.5D0)      &
                                & /sf5x(k, i, 1)
                     ad2x(k, i) = sf4x(k, i, 1)                                &
                                & + (alfa(k, i) - cour*sf9x(k, i, 1)*0.5D0)    &
                                & /sf6x(k, i, 1)
                     ad3x(k, i) = sf2x(k, i, 1)                                &
                                & + (alfa(k, i) - cour*sf10x(k, i, 1)*0.5D0)   &
                                & /sf7x(k, i, 1)
                 else
                     rats(k, i) = ratd(i + 1)
                     curs1(k, i) = curx1(i + 1)
                     curs2(k, i) = curx2(i + 1)
                     curs3(k, i) = curx3(i + 1)
                     dx1(k, i) = DX(k, i)*sf11x(k, i, 2)
                     dx2(k, i) = DX(k, i)*sf12x(k, i, 2)
                     dx3(k, i) = DX(k, i)*sf13x(k, i, 2)
                     alfa(k, i) = 2.0D0*(DX(k, i)*dlt/(sf1x(k, i)*sf1x(k, i))  &
                                & - (1.0D0 - cour*cour)*0.1666667D0)           &
                                & *sf3x(k, i, 2)                                                                   !/6.0
                     ad1x(k, i) = sf2x(k, i, 2)                                &
                                & + (alfa(k, i) - cour*sf8x(k, i, 2)*0.5D0)    &
                                & /sf5x(k, i, 2)
                     ad2x(k, i) = sf4x(k, i, 2)                                &
                                & + (alfa(k, i) - cour*sf9x(k, i, 2)*0.5D0)    &
                                & /sf6x(k, i, 2)
                     ad3x(k, i) = (alfa(k, i) - cour*sf10x(k, i, 2)*0.5D0)     &
                                & /sf7x(k, i, 2)
                 endif
             enddo
         enddo
!!$OMP   END PARALLEL DO
     endif
     return
 
     entry HORIZONTAL_MULTIPLIERS
 
!    Horizontal advection and diffusion multipliers
 
     if(UPWIND(jw))then
         do i = iu, id - 1
             do k = kt, KB(i)
                 if(U(k, i)>=0.0)then
                     c2x = COLD(k, i)
                     c3x = COLD(k, i + 1)
                     adx(k, i) = (dx2(k, i) - U(k, i))*c2x + dx3(k, i)*c3x
                 else
                     c1x = COLD(k, i)
                     c2x = COLD(k, i + 1)
                     adx(k, i) = dx1(k, i)*c1x + (dx2(k, i) - U(k, i))*c2x
                 endif
             enddo
         enddo
     else
!!$OMP   PARALLEL DO
         do i = iu, id - 1
             dltdlxr = dlt/DLXR(i)
             do k = kt, KB(i)
                 cour = U(k, i)*dltdlxr
                 if(U(k, i)>=0.0)then
                     c1x = COLD(k, i - 1)
                     c2x = COLD(k, i)
                     c3x = COLD(k, i + 1)
                     if(U(k, i - 1)<=0.0 .OR. k>KB(i - 1) .OR.                 &
                      & INTERNAL_WEIR(k, i - 1))c1x = COLD(k, i)
                     if(INTERNAL_WEIR(k, i))c3x = COLD(k, i)
                     cart = c3x
                     calf = c1x
                 else
                     c1x = COLD(k, i)
                     c2x = COLD(k, i + 1)
                     c3x = COLD(k, i + 2)
                     if(U(k, i + 2)>=0.0 .OR. k>KB(i + 2) .OR. i==id - 1 .OR.  &
                      & INTERNAL_WEIR(k, i + 1))c3x = COLD(k, i + 1)
                     if(INTERNAL_WEIR(k, i))then
                         c2x = COLD(k, i)
                         c3x = COLD(k, i)
                     endif
                     cart = c1x
                     calf = c3x
                 endif
                 if(.NOT.ULTIMATE(jw))then
                     adx(k, i) = (dx1(k, i) - U(k, i)*ad1x(k, i))              &
                               & *c1x + (dx2(k, i) - U(k, i)*ad2x(k, i))       &
                               & *c2x + (dx3(k, i) - U(k, i)*ad3x(k, i))*c3x
      !  IF (ULTIMATE(JW)) THEN    ! SW Code speedup 6/16/13
                 else
                     ratdi = 1.0/rats(k, i)
                     delc = rats(k, i)*c3x + (ratdi - rats(k, i))              &
                          & *c2x - ratdi*c1x
                     delc = DSIGN(1.0, U(k, i))*delc
                     adelc = DABS(delc)
                     acurz = DABS(curs3(k, i)*c3x + curs2(k, i)                &
                           & *c2x + curs1(k, i)*c1x)
                     if(acurz<=0.6*adelc)then
                         flux = ad1x(k, i)*c1x + ad2x(k, i)*c2x + ad3x(k, i)   &
                              & *c3x
                     elseif(acurz>=adelc)then
                         flux = c2x
                     elseif(ABS(cour)>0.0)then
                         ftemp = ad1x(k, i)*c1x + ad2x(k, i)*c2x + ad3x(k, i)  &
                               & *c3x
                         cref = calf + (c2x - calf)/ABS(cour)
                         if(delc>0.0)then
                             cmax1 = DMIN1(cref, cart)
                             if(cref<c2x)cmax1 = cart
                             flux = 0.5D0*(c2x + cmax1)
                             if(ftemp<=cmax1 .AND. ftemp>=c2x)flux = ftemp
                         else
                             cmin1 = DMAX1(cref, cart)
                             if(cref>c2x)cmin1 = cart
                             if(ftemp>=cmin1 .AND. ftemp<=c2x)then
                                 flux = ftemp
                             elseif(ftemp>0.0)then
                                 flux = 0.5D0*(c2x + cmin1)
                             else
                                 flux = 0.0D0
                             endif
                         endif
                     else
                         flux = 0.0D0
                     endif
                     adx(k, i) = (dx1(k, i)*c1x + dx2(k, i)*c2x + dx3(k, i)    &
                               & *c3x) - U(k, i)*flux
                 endif
             enddo
         enddo
!!$OMP   END PARALLEL DO
     endif
     return
 
!***********************************************************************************************************************************
!**  V E R T I C A L  M U L T I P L I E R S                                   
!***********************************************************************************************************************************
 
!    **
     entry VERTICAL_MULTIPLIERS1
                               ! FIRST PASS
 
!    Vertical advection multipliers
 
     if(.NOT.UPWIND(jw))then
!!$OMP   PARALLEL DO
         do i = iu, id
             do k = kt, KB(i) - 1
                 if(W(k, i)>=0.0)then
                     ratsz(k, i) = RATZ(k, jw)
                     curs1z(k, i) = CURZ1(k, jw)
                     curs2z(k, i) = CURZ2(k, jw)
                     curs3z(k, i) = CURZ3(k, jw)
                     if(k<=kt + 1)then
                         ht = H1(kt, i)
                         hm = H1(k, i)
                         hb = H1(k + 1, i)
                         ratsz(k, i) = AVH1(kt, i)/AVH1(k, i)
            !CURS1Z(K,I) =  2.0D0*HM*HM/(AVH1(KT,I)+AVH1(K,I))/AVH1(KT,I)
                         curs1z(k, i) = 2.0D0*hm*hm/((AVH1(kt, i) + AVH1(k, i))&
                                      & *AVH1(kt, i))                           ! SW 4/20/16 SPEED
                         curs2z(k, i) = -2.0D0*hm*hm/(AVH1(kt, i)*AVH1(k, i))
            !CURS3Z(K,I) =  2.0D0*HM*HM/(AVH1(KT,I)+AVH1(K,I))/AVH1(K,I)
                         curs3z(k, i) = 2.0D0*hm*hm/((AVH1(kt, i) + AVH1(k, i))&
                                      & *AVH1(k, i))                            ! SW 4/20/16 SPEED
                         if(k==kt)then
                             hm = H1(kt, i)
                             ratsz(k, i) = 1.0D0
                             curs3z(k, i) = 1.0D0
                             curs2z(k, i) = -2.0D0
                             curs1z(k, i) = 1.0D0
                         endif
                         hmin = DMIN1(hb, hm)
                         sf1z(k, jw) = (hb + hm)*0.5D0
            !SF2Z(K,1,JW)  =  HM**2
                         sf2z(k, 1, jw) = hm*hm
                                          ! SW 4/20/16
                         sf3z(k, 1, jw) = hm/(hm + hb)
                         sf4z(k, 1, jw) = hb/(hm + hb)
                         sf5z(k, 1, jw) = 0.25D0*(ht + 2.0D0*hm + hb)*(ht + hm)
                         sf6z(k, 1, jw) = -0.25D0*(hm + hb)*(ht + hm)
                         sf7z(k, 1, jw) = 0.25D0*(hm + hb)*(ht + 2.0D0*hm + hb)
                         sf8z(k, 1, jw) = 0.5D0*(hm - hb)*hmin
                         sf9z(k, 1, jw) = 0.5D0*(ht + 2.0D0*hm - hb)*hmin
                         sf10z(k, 1, jw) = 0.5D0*(ht + 3.0D0*hm)*hmin
                     endif
                     cour = W(k, i)*dlt/sf1z(k, jw)
                     alfaz(k, i) = 2.0D0*(DZQ(k, i)*dlt/(sf1z(k, jw)*sf1z(k, jw&
                                 & )) - (1.0D0 - cour*cour)*0.16666667D0)      &
                                 & *sf2z(k, 1, jw)                                                                   !/6.0
                     ad1z(k, i) = (alfaz(k, i) - cour*sf8z(k, 1, jw)*0.5D0)    &
                                & /sf5z(k, 1, jw)
                     ad2z(k, i) = sf4z(k, 1, jw)                               &
                                & + (alfaz(k, i) - cour*sf9z(k, 1, jw)*0.5D0)  &
                                & /sf6z(k, 1, jw)
                     ad3z(k, i) = sf3z(k, 1, jw)                               &
                                & + (alfaz(k, i) - cour*sf10z(k, 1, jw)*0.5D0) &
                                & /sf7z(k, 1, jw)
                 else
                     curs3z(k, i) = CURZ3(k + 1, jw)
                     curs2z(k, i) = CURZ2(k + 1, jw)
                     curs1z(k, i) = CURZ1(k + 1, jw)
                     ratsz(k, i) = AVH1(k, i)/AVH1(k + 1, i)
                     if(k==kt)then
                         ht = H1(kt, i)
                         hm = H1(kt + 1, i)
                         hb = H1(kt + 2, i)
                         hmin = DMIN1(ht, hm)
                         ratsz(k, i) = AVH1(kt, i)/AVH1(k, i)
            !CURS1Z(K,I)         =  2.0D0*HM*HM/(AVH1(KT,I)+AVH1(K,I))/AVH1(KT,I)        ! SW 4/20/16 SPEED
                         curs1z(k, i) = 2.0D0*hm*hm/((AVH1(kt, i) + AVH1(k, i))&
                                      & *AVH1(kt, i))
                         curs2z(k, i) = -2.0D0*hm*hm/(AVH1(kt, i)*AVH1(k, i))
            !CURS3Z(K,I)         =  2.0D0*HM*HM/(AVH1(KT,I)+AVH1(K,I))/AVH1(K,I)          ! SW 4/20/16 SPEED
                         curs3z(k, i) = 2.0D0*hm*hm/((AVH1(kt, i) + AVH1(k, i))&
                                      & *AVH1(k, i))
                         sf1z(k, jw) = (hm + ht)*0.5D0
            !SF2Z(K,2,JW)  =  HM**2
                         sf2z(k, 2, jw) = hm*hm
                                        ! SW 4/20/16
                         sf3z(k, 2, jw) = hm/(ht + hm)
                         sf4z(k, 2, jw) = ht/(ht + hm)
                         sf5z(k, 2, jw) = 0.25D0*(ht + 2.0D0*hm + hb)*(ht + hm)
                         sf6z(k, 2, jw) = -0.25D0*(hm + hb)*(ht + hm)
                         sf7z(k, 2, jw) = 0.25D0*(ht + 2.0D0*hm + hb)*(hm + hb)
                         sf8z(k, 2, jw) = -0.5D0*(3.0D0*hm + hb)*hmin
                         sf9z(k, 2, jw) = 0.5D0*(ht - 2.0D0*hm - hb)*hmin
                         sf10z(k, 2, jw) = 0.5D0*(ht - hm)*hmin
                     endif
                     cour = W(k, i)*dlt/sf1z(k, jw)
                     alfaz(k, i) = 2.0D0*(DZQ(k, i)*dlt/(sf1z(k, jw)*sf1z(k, jw&
                                 & )) - (1.0D0 - cour*cour)*0.16666667D0)      &
                                 & *sf2z(k, 2, jw)                                                                  !/6.0
                     ad1z(k, i) = sf3z(k, 2, jw)                               &
                                & + (alfaz(k, i) - cour*sf8z(k, 2, jw)*0.5D0)  &
                                & /sf5z(k, 2, jw)
                     ad2z(k, i) = sf4z(k, 2, jw)                               &
                                & + (alfaz(k, i) - cour*sf9z(k, 2, jw)*0.5D0)  &
                                & /sf6z(k, 2, jw)
                     ad3z(k, i) = (alfaz(k, i) - cour*sf10z(k, 2, jw)*0.5D0)   &
                                & /sf7z(k, 2, jw)
                 endif
             enddo
         enddo
!$END    PARALLEL DO
     endif
     return
 
     entry VERTICAL_MULTIPLIERS
 
!    Vertical advection multipliers
 
     if(UPWIND(jw))then
         do i = iu, id
             do k = kt, KB(i) - 1
                 c2z = COLD(k + 1, i)
                 if(W(k, i)>=0.0)c2z = COLD(k, i)
                 adz(k, i) = -W(k, i)*c2z
             enddo
         enddo
     else
!!$OMP   PARALLEL DO
         do i = iu, id
             do k = kt, KB(i) - 1
                 if(W(k, i)>=0.0)then
                     c1z = COLD(k - 1, i)
                     c2z = COLD(k, i)
                     c3z = COLD(k + 1, i)
                     cart = c3z
                     calf = c1z
                     if(k<=kt + 1)then
                         c1z = COLD(kt, i)
                         calf = c1z
                     endif
                 else
                     c1z = COLD(k, i)
                     c2z = COLD(k + 1, i)
                     c3z = COLD(k + 2, i)
                     if(k==KB(i) - 1)c3z = COLD(k + 1, i)
                     cart = c1z
                     calf = c3z
                 endif
                 if(.NOT.ULTIMATE(jw))then
                     adz(k, i) = -W(k, i)                                      &
                               & *(ad1z(k, i)*c1z + ad2z(k, i)*c2z + ad3z(k, i)&
                               & *c3z)
                 else
    !    IF (ULTIMATE(JW)) THEN    ! SW code speedup 6/16/13
                     cour = W(k, i)*dlt/sf1z(k, jw)
                     ratzi = 1.0D0/ratsz(k, i)
                     delc = ratsz(k, i)*c3z + (ratzi - ratsz(k, i))            &
                          & *c2z - ratzi*c1z
                     delc = DSIGN(1.0, W(k, i))*delc
                     adelc = DABS(delc)
                     acurz = DABS(curs3z(k, i)*c3z + curs2z(k, i)              &
                           & *c2z + curs1z(k, i)*c1z)
                     if(acurz<=0.6*adelc)then
                         flux = ad1z(k, i)*c1z + ad2z(k, i)*c2z + ad3z(k, i)   &
                              & *c3z
                     elseif(acurz>=adelc)then
                         flux = c2z
                     elseif(DABS(cour)>0.0)then
                         ftemp = ad1z(k, i)*c1z + ad2z(k, i)*c2z + ad3z(k, i)  &
                               & *c3z
                         cref = calf + (c2z - calf)/DABS(cour)
                         if(delc>0.0)then
                             cmax1 = cart
                             if(cref>=c2z)cmax1 = DMIN1(cref, cart)
                             flux = 0.5*(c2z + cmax1)
                             if(ftemp<=cmax1 .AND. ftemp>=c2z)flux = ftemp
                         else
                             cmin1 = DMAX1(cref, cart)
                             if(cref>c2z)cmin1 = cart
                             if(ftemp>=cmin1 .AND. ftemp<=c2z)then
                                 flux = ftemp
                             elseif(ftemp>0.0)then
                                 flux = 0.5D0*(c2z + cmin1)
                             else
                                 flux = 0.0D0
                             endif
                         endif
                     else
                         flux = 0.0D0
                     endif
                     adz(k, i) = -W(k, i)*flux
                 endif
             enddo
         enddo
!$END    PARALLEL DO
     endif
     return
 
!***********************************************************************************************************************************
!**  H O R I Z O N T A L  T R A N S P O R T                                   
!***********************************************************************************************************************************
 
!    **
     entry HORIZONTAL_TRANSPORT
     if(constituents)then
         do i = iu, id
             do k = kt, KB(i)
                        !CONCURRENT(K=KT:KB(I))    !FORALL                                                        !DO K=KT,KB(I)
                 CNEW(k, i) = (COLD(k, i)*BH2(k, i)/dlt + (adx(k, i)*BHR1(k, i)&
                            & - adx(k, i - 1)*BHR1(k, i - 1))/DLX(i)           &
                            & + (1.0D0 - THETA(jw))                            &
                            & *(adz(k, i)*BB(k, i) - adz(k - 1, i)*BB(k - 1, i)&
                            & ) + SSB(k, i)/DLX(i))*dlt/BH1(k, i) + SSK(k, i)  &
                            & *dlt
             enddo
         enddo
     else
         do i = iu, id
             do k = kt, KB(i)
                        !CONCURRENT(K=KT:KB(I))      !FORALL                                         !DO K=KT,KB(I)
                 CNEW(k, i) = (COLD(k, i)*BH2(k, i)/dlt + (adx(k, i)*BHR1(k, i)&
                            & - adx(k, i - 1)*BHR1(k, i - 1))/DLX(i)           &
                            & + (1.0D0 - THETA(jw))                            &
                            & *(adz(k, i)*BB(k, i) - adz(k - 1, i)*BB(k - 1, i)&
                            & ) + SSB(k, i)/DLX(i))*dlt/BH1(k, i)
             enddo
         enddo
     endif
     return
     entry DEALLOCATE_TRANSPORT
     deallocate(ratd, curx1, curx2, curx3, sf1x, sf1z, RATZ, CURZ1, CURZ2,     &
              & CURZ3, ct, at, vt, dt, rats, curs1, curs2, curs3, alfa)
     deallocate(dx1, dx2, dx3, ad1x, ad2x, ad3x, ad1z, ad2z, ad3z, sf2x, sf3x, &
              & sf4x, sf5x, sf6x, sf7x)
     deallocate(ratsz, curs1z, curs2z, curs3z, alfaz)
     deallocate(sf8x, sf9x, sf10x, sf11x, sf12x, sf13x, sf2z, sf3z, sf4z, sf5z,&
              & sf6z, sf7z, sf8z, sf9z, sf10z)
     end subroutine TRANSPORT
