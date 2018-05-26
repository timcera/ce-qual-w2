!*==outputa.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     subroutine OUTPUTA
 
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
     use BIOENERGETICS
     use CEMAVARS
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND) :: cgas, dle, dlvbr, tgate, tspill
     real :: cgasd, qsumm, tvolavg, voltot, xdum
     integer :: ic, iflag, itot, jad, jaf, jj, jjc, jwwd, nlines, numoutlets
     integer, dimension(100) :: jsss
     real, dimension(100) :: qoutlet, toutlet
     real(R4KIND), allocatable, dimension(:, :), save :: wdsi
     real(R4KIND), allocatable, dimension(:), save :: wsel
!
!*** End of declarations rewritten by SPAG
!
 
 
 
  ! *** DSI W2_TOOL LINKAGE
 
     if(.NOT.ALLOCATED(wsel))then
         allocate(wsel(imx))
         allocate(wdsi(kmx, imx))
     endif
 
  !***********************************************************************************************************************************
  !*                                                    Task 2.8: Output Results                                                    **
  !***********************************************************************************************************************************
 
!    code for error checking
  !   if(macrophyte_on)then
  !     if(nit.eq.1)open(689,file='totmac.dat',status='unknown')
  !     if((nit/200)*200.eq.nit)then
  !       bmass=0.0
  !       totar=0.0
  !       DO Jw=1,Nwb
  !         KT = KTwb(Jw)
  !         DO JB=BS(Jw),BE(Jw)
  !           IU = CUS(JB)
  !           ID = DS(JB)
  !           DO I=IU,ID
  !             smass=0.0
  !             do k=kt,kb(i)
  !               do m=1,nmc
  !                 bmass=bmass+mac(k,i,m)*bh2(k,i)*dlx(i)
  !                 smass=smass+mac(k,i,m)*bh2(k,i)*dlx(i)
  !               end do
  !             end do
  !             totar=totar+dlx(i)*b(kti(i),i)
  !             armac(i)=smass/(dlx(i)*b(kti(i),i))
  !           end do
  !         end do
  !       end do
  !       bmassar=bmass/totar
  !       write(689,'(f10.4,2f10.3)')jday,bmass/1000.0,bmassar
  !     end if
  !   end if
 
  ! BIOEXP MLM
     if(bioexp)then
         if(jday>=nxtbio)then
             nxtbio = nxtbio + BIOF(1)
                               ! mlm3
             do j = 1, nibio
 
                 i = IBIO(j)
                 do jw = 1, nwb
                     if(i>=US(BS(jw)) - 1 .AND. i<=DS(BE(jw)) + 1)exit
                 enddo
        ! bioexp
                 do k = KTWB(jw), KBI(i)
                     do jc = nzoos, nzooe
                         C2ZOO(k, i, jc) = C2(k, i, jc)*CMULT(jc)
                     enddo
                     write(BIOEXPFN(j),                                        &
     &'(F8.2,",",I8,",",3(F8.2,","),<NZOOE-NZOOS+1>(F8.3,","),I8,",",2(F8.2,","&
     &),A,",",I0,",",I0)')jday, IBIO(j), DEPTHM(k, i), T1(k, i), GAMMA(k, i),  &
                        & (C2ZOO(k, i, jc), jc = nzoos, nzooe), k, BH(k, i),   &
                        & EL(k, i), month, gday, year                                                                                                               !          ! CB 1/6/17
 
                 enddo
           ! VOLUME WEIGHTING OF ACTIVE CONSTITUENTS     !MLM 18.07.06
                 klim = MIN(KB(i), KTWB(jw) + 5)
                 do k = KTWB(jw), klim
                     do jjc = 1, nac
                         C2W(i, jjc) = C2W(i, jjc) + C2(k, i, CN(jjc))         &
                                     & *CMULT(CN(jjc))*BH(k, i)
                     enddo
                     C2W(i, nac + 1) = C2W(i, nac + 1) + CD(k, i, 20)*BH(k, i)
                                                                            !PH
                     C2W(i, nac + 2) = C2W(i, nac + 2) + CD(k, i, 12)*BH(k, i)  ! TOTAL PHOSPHOROUS
                     VOLROOS(i) = VOLROOS(i) + BH(k, i)                         ! VOLUME WEIGHTED
                 enddo
                 do jjc = 1, nac + 2
                     C2W(i, jjc) = C2W(i, jjc)/VOLROOS(i)
                 enddo
 
                 write(WEIGHTNUM(j), '(F10.3,",",<NAC>(F9.3,","),2(F9.3,","))')&
                     & jday, (C2W(i, jjc), jjc = 1, nac + 2)
                 VOLROOS = 0.0
 
             enddo
         endif
     endif
 
  !** Time series
 
     if(time_series)then
         if(jday>=nxtmts .OR. jday>=TSRD(tsrdp + 1))then
             if(jday>=TSRD(tsrdp + 1))then
                 tsrdp = tsrdp + 1
                 nxtmts = TSRD(tsrdp)
             endif
             nxtmts = nxtmts + TSRF(tsrdp)
 
      ! write out water level File
             write(wlfn, '(f10.3,",",1000(f8.3,","))')jday,                    &
                 & ((ELWS(i), i = US(jb), DS(jb)), jb = 1, nbr)
 
             do j = 1, niktsr
                 i = ITSR(j)
                 do jw = 1, nwb
                     if(i>=US(BS(jw)) - 1 .AND. i<=DS(BE(jw)) + 1)exit
                 enddo
        ! TEMP VOL WEIGHTED AVERAGE  SW 6/1/2015
                 tvolavg = 0.0
                            ! SW 6/1/2015
                 voltot = 0.0
                 do k = KTWB(jw), KB(i)
                     voltot = voltot + VOL(k, i)
                     tvolavg = tvolavg + T1(k, i)*VOL(k, i)
                 enddo
                 if(KB(i)>=KTWB(jw))tvolavg = tvolavg/voltot
 
                 if(ETSR(j)<0)then
                     k = INT(ABS(ETSR(j)))
                 else
                     do k = KTWB(jw), KB(i)
                         if(DEPTHB(k, i)>ETSR(j))exit
                     enddo
                     if(k>KB(i))cycle
                 endif
 
                 if(k>=KTWB(jw))then
                               ! SW 4/4/2018 ADDED TO ELIMINATE THE LAST VALUE OF A VARIABLE BEING USED FOR CASE WHEN ESTR IS NEGATIVE
 
                     do jac = 1, nac
                         l = LEN_TRIM(FMTC(CN(jac)))
                         write(C2CH(jac), FMTC(CN(jac))(1:l))C2(k, i, CN(jac)) &
                             & *CMULT(CN(jac))
                     enddo
                     do jf = 1, NAF(jw)
                         write(KFCH(jf), '(E10.3)')KF(k, i, (KFCN(jf, jw)))    &
                             & *VOL(k, i)/1000./day
                     enddo
                     do jad = 1, NACD(jw)
                         l = LEN_TRIM(FMTCD(CDN(jad, jw)))
                         write(CDCH(jad), FMTCD(CDN(jad, jw))(1:l))            &
                             & CD(k, i, CDN(jad, jw))*CDMULT(CDN(jad, jw))
                     enddo
                     do je = 1, nep
                         write(EPCH(je), '(F10.3)')EPD(k, i, je)                                    ! SW 8/13/06
                     enddo
                     do ja = 1, nal
                         write(APCH(ja), '(F10.3)')APLIM(k, i, ja)                                    ! SW 8/13/06
                         write(ANCH(ja), '(F10.3)')ANLIM(k, i, ja)                                    ! SW 8/13/06
                         write(ALCH(ja), '(F10.3)')ALLIM(k, i, ja)                                    ! SW 8/13/06
                     enddo
                     do jm = 1, nmc
                         write(MACCH(jm), '(F10.3)')MAC(k, i, jm)                                   ! SW 8/13/06
                     enddo
                     if(SEDIMENT_CALC(jw))then
                         write(sedch, '(F10.3)')SED(k, i)                                           ! SW 8/13/06
                         write(sedpch, '(F10.3)')SEDP(k, i)
                         write(sednch, '(F10.3)')SEDN(k, i)
                         write(sedcch, '(F10.3)')SEDC(k, i)
                     endif
                     if(ice_computation)then
                         if(SEDIMENT_CALC(jw))then
                                                                                                                                                      !     ! SW 8/13/06
                                                                                                                                                  !  ! CB 7/26/07
                             write(TSR(j),                                     &
                                  &'(f10.3,",",18(F10.3,","),1000(A,","))')    &
                                 & jday, dlt, ELWS(i), T1(k, i), U(k, i),      &
                                 & QC(i), SRON(jw)*1.06, GAMMA(k, i),          &
                                 & DEPTHB(KB(i), i), BI(KTWB(jw), i), SHADE(i),&
                                 & ICETH(i), tvolavg, RN(i), RS(i), RANLW(jw), &
                                 & RB(i), RE(i), RC(i),                        &
                                 & (ADJUSTR(C2CH(jac)), jac = 1, nac),         &
                                 & (ADJUSTR(EPCH(je)), je = 1, nep),           &
                                 & (ADJUSTR(MACCH(jm)), jm = 1, nmc), sedch,   &
                                 & sedpch, sednch, sedcch,                     &
                                 & (ADJUSTR(CDCH(jad)), jad = 1, NACD(jw)),    &
                                 & (ADJUSTR(KFCH(jf)), jf = 1, NAF(jw)),       &
                                 & (ADJUSTR(APCH(ja)), ja = 1, nal),           &
                                 & (ADJUSTR(ANCH(ja)), ja = 1, nal),           &
                                 & (ADJUSTR(ALCH(ja)), ja = 1, nal)                                                                                                        ! SW 10/20/15
                         else
                                                                                                                                                        !     ! SW 8/13/06
                                                                                                                                                  !  ! CB 7/26/07
                             write(TSR(j),                                     &
                                  &'(f10.3,",",18(F10.3,","),1000(A,","))')    &
                                 & jday, dlt, ELWS(i), T1(k, i), U(k, i),      &
                                 & QC(i), SRON(jw)*1.06, GAMMA(k, i),          &
                                 & DEPTHB(KB(i), i), BI(KTWB(jw), i), SHADE(i),&
                                 & ICETH(i), tvolavg, RN(i), RS(i), RANLW(jw), &
                                 & RB(i), RE(i), RC(i),                        &
                                 & (ADJUSTR(C2CH(jac)), jac = 1, nac),         &
                                 & (ADJUSTR(EPCH(je)), je = 1, nep),           &
                                 & (ADJUSTR(MACCH(jm)), jm = 1, nmc),          &
                                 & (ADJUSTR(CDCH(jad)), jad = 1, NACD(jw)),    &
                                 & (ADJUSTR(KFCH(jf)), jf = 1, NAF(jw)),       &
                                 & (ADJUSTR(APCH(ja)), ja = 1, nal),           &
                                 & (ADJUSTR(ANCH(ja)), ja = 1, nal),           &
                                 & (ADJUSTR(ALCH(ja)), ja = 1, nal)                                                                                                        ! SW 10/20/15
                         endif
                     elseif(SEDIMENT_CALC(jw))then
                                                                                                                                                      !     ! SW 8/13/06
                                                                                                                                                  !  ! CB 7/26/07
                         write(TSR(j), '(f10.3,",",17(F10.3,","),1000(A,","))')&
                             & jday, dlt, ELWS(i), T1(k, i), U(k, i), QC(i),   &
                             & SRON(jw)*1.06, GAMMA(k, i), DEPTHB(KB(i), i),   &
                             & BI(KTWB(jw), i), SHADE(i), tvolavg, RN(i),      &
                             & RS(i), RANLW(jw), RB(i), RE(i), RC(i),          &
                             & (ADJUSTR(C2CH(jac)), jac = 1, nac),             &
                             & (ADJUSTR(EPCH(je)), je = 1, nep),               &
                             & (ADJUSTR(MACCH(jm)), jm = 1, nmc), sedch,       &
                             & sedpch, sednch, sedcch,                         &
                             & (ADJUSTR(CDCH(jad)), jad = 1, NACD(jw)),        &
                             & (ADJUSTR(KFCH(jf)), jf = 1, NAF(jw)),           &
                             & (ADJUSTR(APCH(ja)), ja = 1, nal),               &
                             & (ADJUSTR(ANCH(ja)), ja = 1, nal),               &
                             & (ADJUSTR(ALCH(ja)), ja = 1, nal)                                                                                                            ! SW 10/20/15
                     else
                                                                                                                                                        !      ! SW 8/13/06
                                                                                                                                                  !  ! CB 7/26/07
                         write(TSR(j), '(f10.3,",",17(F10.3,","),1000(A,","))')&
                             & jday, dlt, ELWS(i), T1(k, i), U(k, i), QC(i),   &
                             & SRON(jw)*1.06, GAMMA(k, i), DEPTHB(KB(i), i),   &
                             & BI(KTWB(jw), i), SHADE(i), tvolavg, RN(i),      &
                             & RS(i), RANLW(jw), RB(i), RE(i), RC(i),          &
                             & (ADJUSTR(C2CH(jac)), jac = 1, nac),             &
                             & (ADJUSTR(EPCH(je)), je = 1, nep),               &
                             & (ADJUSTR(MACCH(jm)), jm = 1, nmc),              &
                             & (ADJUSTR(CDCH(jad)), jad = 1, NACD(jw)),        &
                             & (ADJUSTR(KFCH(jf)), jf = 1, NAF(jw)),           &
                             & (ADJUSTR(APCH(ja)), ja = 1, nal),               &
                             & (ADJUSTR(ANCH(ja)), ja = 1, nal),               &
                             & (ADJUSTR(ALCH(ja)), ja = 1, nal)                                                                        ! SW 10/20/15
                     endif
                 elseif(ice_computation)then
                     if(SEDIMENT_CALC(jw))then
                                                                                                                                           !     ! SW 8/13/06
                                                                                                                                    !  ! CB 7/26/07
                         write(TSR(j), '(f10.3,",",18(F10.3,","),1000(A,","))')&
                             & jday, dlt, ELWS(i), -99., -99., QC(i), SRON(jw) &
                             & *1.06, -99., DEPTHB(KB(i), i), BI(KTWB(jw), i), &
                             & SHADE(i), ICETH(i), tvolavg, RN(i), RS(i),      &
                             & RANLW(jw), RB(i), RE(i), RC(i),                 &
                             & ( - 99., jac = 1, nac), ( - 99., je = 1, nep),  &
                             & ( - 99., jm = 1, nmc), sedch, sedpch, sednch,   &
                             & sedcch, ( - 99., jad = 1, NACD(jw)),            &
                             & ( - 99., jf = 1, NAF(jw)), ( - 99., ja = 1, nal)&
                             & , ( - 99., ja = 1, nal), ( - 99., ja = 1, nal)                             ! SW 10/20/15
                     else
                                                                                                                                             !     ! SW 8/13/06
                                                                                                                                    !  ! CB 7/26/07
                         write(TSR(j), '(f10.3,",",18(F10.3,","),1000(A,","))')&
                             & jday, dlt, ELWS(i), -99., -99., QC(i), SRON(jw) &
                             & *1.06, -99., DEPTHB(KB(i), i), BI(KTWB(jw), i), &
                             & SHADE(i), ICETH(i), tvolavg, RN(i), RS(i),      &
                             & RANLW(jw), RB(i), RE(i), RC(i),                 &
                             & ( - 99., jac = 1, nac), ( - 99., je = 1, nep),  &
                             & ( - 99., jm = 1, nmc),                          &
                             & ( - 99., jad = 1, NACD(jw)),                    &
                             & ( - 99., jf = 1, NAF(jw)), ( - 99., ja = 1, nal)&
                             & , ( - 99., ja = 1, nal), ( - 99., ja = 1, nal)                             ! SW 10/20/15
                     endif
                 elseif(SEDIMENT_CALC(jw))then
                                                                                                                                           !     ! SW 8/13/06
                                                                                                                           !  ! CB 7/26/07
                     write(TSR(j), '(f10.3,",",18(F10.3,","),1000(A,","))')    &
                         & jday, dlt, ELWS(i), -99., -99., QC(i), SRON(jw)     &
                         & *1.06, -99., DEPTHB(KB(i), i), BI(KTWB(jw), i),     &
                         & SHADE(i), tvolavg, RN(i), RS(i), RANLW(jw), RB(i),  &
                         & RE(i), RC(i), ( - 99., jac = 1, nac),               &
                         & ( - 99., je = 1, nep), ( - 99., jm = 1, nmc), sedch,&
                         & sedpch, sednch, sedcch, ( - 99., jad = 1, NACD(jw)),&
                         & ( - 99., jf = 1, NAF(jw)), ( - 99., ja = 1, nal),   &
                         & ( - 99., ja = 1, nal), ( - 99., ja = 1, nal)                                  ! SW 10/20/15
                 else
                                                                                                                                             !     ! SW 8/13/06
                                                                                                                           !  ! CB 7/26/07
                     write(TSR(j), '(f10.3,",",18(F10.3,","),1000(A,","))')    &
                         & jday, dlt, ELWS(i), -99., -99., QC(i), SRON(jw)     &
                         & *1.06, -99., DEPTHB(KB(i), i), BI(KTWB(jw), i),     &
                         & SHADE(i), tvolavg, RN(i), RS(i), RANLW(jw), RB(i),  &
                         & RE(i), RC(i), ( - 99., jac = 1, nac),               &
                         & ( - 99., je = 1, nep), ( - 99., jm = 1, nmc),       &
                         & ( - 99., jad = 1, NACD(jw)),                        &
                         & ( - 99., jf = 1, NAF(jw)), ( - 99., ja = 1, nal),   &
                         & ( - 99., ja = 1, nal), ( - 99., ja = 1, nal)                                  ! SW 10/20/15
                 endif
             ! SW 4/4/2018
             enddo
         endif
     endif
 
     ! SEDIMENT DIAGENESIS FREQUENCY OUTPUT
     if(constituents .AND. sediment_diagenesis .AND. jday>=nxtsediag)then
         nxtsediag = nxtsediag + sediagfreq
         call WRITECEMASEDIMENTMODELOUTPUT
         call WRITECEMASEDIMENTFLUXOUTPUT
     endif
 
     do jw = 1, nwb
 
    !**** Inactive segments
 
         jb = BS(jw)
         NBL(jw) = 1
         IBPR(jw) = 1
         do i = 1, NISNP(jw) - 1
             if(CUS(jb)>ISNP(i, jw))then
                 BL(NBL(jw), jw) = i
                 NBL(jw) = NBL(jw) + 1
                 IBPR(jw) = i + 1
             endif
             if(ISNP(i + 1, jw)>DS(jb))jb = jb + 1
         enddo
         NBL(jw) = NBL(jw) - 1
 
    !**** Snapshots
 
         if(SNAPSHOT(jw))then
             if(jday>=NXTMSN(jw) .OR. jday>=SNPD(SNPDP(jw) + 1, jw))then
                 if(jday>=SNPD(SNPDP(jw) + 1, jw))then
                     SNPDP(jw) = SNPDP(jw) + 1
                     NXTMSN(jw) = SNPD(SNPDP(jw), jw)
                 endif
                 NXTMSN(jw) = NXTMSN(jw) + SNPF(SNPDP(jw), jw)
                 write(SNP(jw), 9001)w2ver, (TITLE(j), j = 1, 10)
 
  ! Snapshot formats
 
9001             format('CE-QUAL-W2 VERSION', f4.2/(1x, a72))
                 write(SNP(jw), 9002)'Time Parameters', month, gday, year,     &
                                   & INT(jday), (jday - INT(jday))*24.0,       &
                                   & INT(eltmjd), (eltmjd - INT(eltmjd))*24.0, &
                                   & INT(dlts1), kloc, iloc, INT(mindlt),      &
                                   & INT(jdmin), (jdmin - INT(jdmin))*24.0,    &
                                   & kmin, imin
9002             format(/1x, a/3x, 'Gregorian date      [GDAY] =', a19, 1x, i0,&
                       &', ', i0/3x, 'Julian date         [JDAY] =', i10,      &
                       &' days', f6.2, ' hours'/3x,                            &
                       &'Elapsed time      [ELTMJD] =', i10, ' days', f6.2,    &
                       &' hours'/3x, 'Timestep             [DLT] =', i10,      &
                       &' sec'/3x, '  at location  [KLOC,ILOC] = (', i0, ',',  &
                      & i0, ')'/3x, 'Minimum timestep  [MINDLT] =', i10,       &
                       &' sec '/3x, '  at Julian day    [JDMIN] =', i10,       &
                       &' days', f6.2, ' hours'/3x,                            &
                       &'  at location  [KMIN,IMIN] = (', i0, ',', i0, ')')
                 if(LIMITING_DLT(jw))write(SNP(jw), 9003)kmin, imin
9003             format(3x, 'Limiting timestep'/3x,                            &
                       &'  at location  [KMIN,IMIN] = (', i0, ',', i0, ')')
                 write(SNP(jw), 9004)INT(dltav), nit, nv
9004             format(3x, 'Average timestep   [DLTAV] =', i10, ' sec'/3x,    &
                       &'Number of iterations [NIT] =', i10/3x,                &
                       &'Number of violations  [NV] =', i10/)
                 write(SNP(jw), 9005)'Meteorological Parameters'
9005             format(1x, a)
                 write(SNP(jw), 9006)TAIR(jw), deg, TDEW(jw), deg, PHI(jw),    &
                                   & CLOUD(jw), ET(DS(1)), deg, CSHE(DS(1)),   &
                                   & SRON(jw), deg
9006             format(3x, 'Input'/3x, '  Air temperature          [TAIR] =', &
                      & f9.2, 1x, a/3x, '  Dewpoint temperature     [TDEW] =', &
                      & f9.2, 1x, a/3x, '  Wind direction            [PHI] =', &
                      & f9.2, ' rad'/3x, '  Cloud cover             [CLOUD] =',&
                      & f9.2/3x, '  Calculated'/5x,                            &
                       &'  Equilibrium temperature    [ET] =', f9.2, 1x, a/5x, &
                       &'  Surface heat exchange    [CSHE] =', e9.2,           &
                      & ' m/sec'/5x, '  Net short wave radiation [SRON] =',    &
                      & e9.2, 1x, a, ' W/m^2'/)
                 write(SNP(jw), 9007)'Inflows', 'Upstream inflows'
9007             format(1x, a/3x, a)
                 do jb = BS(jw), BE(jw)
                     if(UP_FLOW(jb))write(SNP(jw), 9008)jb, KTQIN(jb),         &
                      & KBQIN(jb), QIN(jb), TIN(jb), deg
9008                 format(5x, 'Branch ', i0/5x, '  Layer       [KQIN] = ',   &
                          & i0, '-', i0/5x, '  Inflow       [QIN] =', f8.2,    &
                           &' m^3/sec'/5x, '  Temperature  [TIN] =', f8.2, 1x, &
                          & a)
                 enddo
                 do jb = BS(jw), BE(jw)
                     if(DIST_TRIBS(jb))then
                         write(SNP(jw), 9009)
9009                     format(/3x, 'Distributed Tributaries')
                         write(SNP(jw), 9010)jb, QDTR(jb), TDTR(jb), deg
9010                     format(5x, 'Branch ', i0/5x, '  Inflow      [QDTR] =',&
                              & f8.2, ' m^3/sec'/5x, '  Temperature [TDTR] =', &
                              & f8.2, 1x, a)
                     endif
                 enddo
                 if(tributaries)then
                     write(SNP(jw), 9011)(ITR(jt), jt = 1, jtt)
9011                 format(:/3x, 'Tributaries'/5x, 'Segment     [ITR] =',     &
                          & 11I8:/(t25, 11I8))
                     write(SNP(jw), 9012)(KTTR(jt), KBTR(jt), jt = 1, jtt)
9012                 format(:5x, 'Layer      [KTWB] = ', 11(i0, '-', i0, 2x)   &
                          & :/(t25, 11(i0, '-', i0)))
                     write(SNP(jw), 9013)(QTR(jt), jt = 1, jtt)
9013                 format(:5x, 'Inflow      [QTR] =', 11F8.2:/(t25, 11F8.1))
                     write(SNP(jw), 9014)(TTR(jt), jt = 1, jtt)
9014                 format(:5x, 'Temperature [TTR] =', 11F8.2:/(t25, 11F8.1))
                 endif
                 write(SNP(jw), 9015)
9015             format(/1x, 'Outflows')
                 do jb = BS(jw), BE(jw)
                     if(DN_FLOW(jb))then
                         write(SNP(jw), 9016)jb,                               &
                             & (QSTR(js, jb), js = 1, JSS(jb))
9016                     format(3x, 'Structure outflows [QSTR]'/3x,            &
                              & '  Branch ', i0, ' = ', 11F8.2:/(t16, 11F8.2))
                         write(SNP(jw), 9017)QSUM(jb),                         &
                             & (k, k = KTWB(jw), KB(DS(jb)))
9017                     format(:/3x, 'Total outflow [QOUT] =', f8.2,          &
                              & ' m^3/s'/5x, 'Outlets'/5x,                     &
                               &'  Layer             [KOUT] =',                &
                              & 12I7:/(33x, 12I7))
                         write(SNP(jw), 9018)                                  &
                             & (QOUT(k, jb), k = KTWB(jw), KB(DS(jb)))
9018                     format(:7x, 'Outflow (m^3/sec) [QOUT] =',             &
                              & 12F7.2:/(33x, 12F7.2))
                     endif
                 enddo
                 if(withdrawals)then
                     do jwd = 1, jww
                         write(SNP(jw), 9019)MAX(CUS(JBWD(jwd)), IWD(jwd)),    &
                             & QWD(jwd)
9019                     format(:5x, 'Withdrawals'/5x,                         &
                               &'  Segment            [IWD] =', i7/5x,         &
                               &'  Outflow (m^3/sec)  [QWD] =', f7.2)
                         if(QWD(jwd)/=0.0)then
                             write(SNP(jw), 9046)(k, k = KTW(jwd), KBW(jwd))
                             write(SNP(jw), 9047)                              &
                                 & (QSW(k, jwd), k = KTW(jwd), KBW(jwd))
                         else
                             write(SNP(jw), 9046)
                             write(SNP(jw), 9047)QWD(jwd)
                         endif
                     enddo
                 endif
                 if(constituents)then
                     write(SNP(jw), 9020)'Constituent Inflow Concentrations'
9020                 format(/'1', a)
                     do jb = BS(jw), BE(jw)
                         if(UP_FLOW(jb) .AND. NACIN(jb)>0)write(SNP(jw), 9021) &
                          & jb,                                                &
                          & (cname1(INCN(jc, jb))(1:18), CIN(INCN(jc, jb), jb),&
                          & CUNIT2(INCN(jc, jb)), jc = 1, NACIN(jb))
9021                     format(3x, 'Branch ', i0,                             &
                              & ' [CIN]'/(5x, a, t25, '=', f9.3, 1x, a))
                         if(DIST_TRIBS(jb) .AND. NACDT(jb)>0)                  &
                          & write(SNP(jw), 9022)jb,                            &
                          & (cname1(DTCN(jc, jb))(1:18), CDTR(DTCN(jc, jb), jb)&
                          & , CUNIT2(DTCN(jc, jb)), jc = 1, NACDT(jb))
9022                     format(3x, 'Distributed tributary ', i0,              &
                              & ' [CDT]'/(5x, a, t25, '=', f9.3, 1x, a))
                     enddo
                     do jt = 1, ntr
                         if(NACTR(jt)>0)write(SNP(jw), 9023)jt,                &
                          & (cname1(TRCN(jc, jt))(1:18), CTR(TRCN(jc, jt), jt),&
                          & CUNIT2(TRCN(jc, jt)), jc = 1, NACTR(jt))
9023                     format(3x, 'Tributary ', i0,                          &
                              & ' [CTR]'/(5x, a, t25, '=', f9.3, 1x, a))
                     enddo
                 endif
                 if(EVAPORATION(jw) .OR. PRECIPITATION(jw))write(SNP(jw), 9024)
9024             format(/'Surface calculations')
                 if(EVAPORATION(jw))then
                     write(SNP(jw), 9025)(jb, EVBR(jb), jb = BS(jw), BE(jw))
                                                                ! SW 9/15/05
9025                 format(3x, 'Evaporation rate [EV]'/(:3x, '  Branch ', i0, &
                           &' = ', e10.3, ' m^3/s'))      ! SW 9/15/05 4/21/10
                     write(SNP(jw), 9026)(jb, -VOLEV(jb), jb = BS(jw), BE(jw))
9026                 format(3x,                                                &
                           &'Cumulative evaporation [VOLEV]'/(:3x, '  Branch ',&
                          & i0, ' = ', f0.1, ' m^3'))
                 endif
                 if(PRECIPITATION(jw))write(SNP(jw), 9027)                     &
                  & (jb, PR(jb), jb = BS(jw), BE(jw))
9027             format(3x, 'Precipitation [PR]'/(3x, '  Branch ', i0, ' = ',  &
                      & f8.6), ' m/s')
                 if(HEAD_BOUNDARY(jw))then
                     write(SNP(jw), 9028)
9028                 format(/1x, 'External head boundary elevations'/)
                     do jb = BS(jw), BE(jw)
                         if(UH_EXTERNAL(jb))write(SNP(jw), 9029)jb, ELUH(jb)
9029                     format(3x, 'Branch ', i0/5x,                          &
                               &'Upstream elevation   [ELUH] =', f8.3, ' m')
                         if(DH_EXTERNAL(jb))write(SNP(jw), 9030)jb, ELDH(jb)
9030                     format(3x, 'Branch ', i0/5x,                          &
                               &'Downstream elevation [ELDH] =', f8.3, ' m')
                     enddo
                 endif
                 if(VOLUME_BALANCE(jw))then
                     write(SNP(jw), 9031)
9031                 format(/'Water Balance')
                     write(SNP(jw), 9032)jw, VOLSR(jw), VOLTR(jw), VOLTR(jw)   &
                         & - VOLSR(jw), DLVR(jw)
9032                 format(3x, 'Waterbody ', i0/3x,                           &
                           &'  Spatial change  [VOLSR]  = ', e15.8, ' m^3'/3x, &
                           &'  Temporal change [VOLTR]  = ', e15.8, ' m^3'/3x, &
                           &'  Volume error             = ', e15.8, ' m^3'/3x, &
                           &'  Percent error            = ', e15.8, ' %')
                     do jb = BS(jw), BE(jw)
                         if(VOLSBR(jb)/=0.0)dlvbr = (VOLTBR(jb) - VOLSBR(jb))  &
                          & /VOLSBR(jb)
                         write(SNP(jw), 9033)jb, VOLSBR(jb), VOLTBR(jb),       &
                             & VOLTBR(jb) - VOLSBR(jb), dlvbr*100.0
9033                     format(3x, 'Branch ', i0/3x,                          &
                               &'  Spatial change  [VOLSBR] = ', e15.8,        &
                              & ' m^3'/3x, '  Temporal change [VOLTBR] = ',    &
                              & e15.8, ' m^3'/3x,                              &
                               &'  Volume error             = ', e15.8,        &
                              & ' m^3'/3x, '  Percent error            = ',    &
                              & e15.8, ' %')
                     enddo
                 endif
                 if(ENERGY_BALANCE(jw))then
                     write(SNP(jw), 9034)
9034                 format(/1x, 'Energy Balance')
                     if(ESR(jw)/=0.0)dle = (ESR(jw) - ETR(jw))/ESR(jw)
                     write(SNP(jw), 9035)jw, ESR(jw)*4.184E3, ETR(jw)*4.184E3, &
                         & (ESR(jw) - ETR(jw))*4.184E3, dle*100.0
9035                 format(3x, 'Waterbody ', i0/3x,                           &
                           &'  Spatially integrated energy   [ESR] = ', e15.8, &
                           &' kJ'/3x,                                          &
                           &'  Temporally integrated energy  [ETR] = ', e15.8, &
                           &' kJ'/3x,                                          &
                           &'  Energy error                        = ', e15.8, &
                           &' kJ'/3x,                                          &
                           &'  Percent error                       = ', e15.8, &
                           &' %')
                     do jb = BS(jw), BE(jw)
                         write(SNP(jw), 9048)jb
                         if(ESBR(jb)/=0.0)dle = (ESBR(jb) - ETBR(jb))/ESBR(jb)
                         write(SNP(jw), 9036)ESBR(jb)*4.184E3, ETBR(jb)        &
                             & *4.1843E3, (ESBR(jb) - ETBR(jb))*4.1843E3,      &
                             & dle*100.0
9036                     format(3x, '  Spatially integrated energy  [ESBR] = ',&
                              & e15.8, ' kJ'/3x,                               &
                               &'  Temporally integrated energy [ETBR] = ',    &
                              & e15.8, ' kJ'/3x,                               &
                               &'  Energy error                        = ',    &
                              & e15.8, ' kJ'/3x,                               &
                               &'  Percent error                       = ',    &
                              & e15.8, ' %')
                     enddo
                 endif
                 if(MASS_BALANCE(jw))then
                     write(SNP(jw), 9037)
9037                 format(/1x, 'Mass Balance')
                     do jb = BS(jw), BE(jw)
                         write(SNP(jw), 9048)jb
                         do jc = 1, nac
                             if(CMBRS(CN(jc), jb)/=0.0)                        &
                              & dlmr = (CMBRT(CN(jc), jb) - CMBRS(CN(jc), jb)) &
                              & /(CMBRS(CN(jc), jb) + nonzero)*100.0
                             write(SNP(jw), 9038)cname1(CN(jc)),               &
                                 & CMBRS(CN(jc), jb), CUNIT1(CN(jc)),          &
                                 & CMBRT(CN(jc), jb), CUNIT1(CN(jc)),          &
                                 & (CMBRT(CN(jc), jb) - CMBRS(CN(jc), jb)),    &
                                 & CUNIT1(CN(jc)), dlmr
9038                         format(5x, a/5x,                                  &
                                   &'  Spatially integrated mass  [CMBRS] = ', &
                                  & e15.8, 1x, a/5x,                           &
                                   &'  Temporally integrated mass [CMBRT] = ', &
                                  & e15.8, 1x, a/5x,                           &
                                   &'  Mass error                         = ', &
                                  & e15.8, 1x, a/5x,                           &
                                   &'  Percent error                      = ', &
                                  & e15.8, ' %')
                         enddo
 
                         do m = 1, nmc
                             if(MACROPHYTE_CALC(jw, m))then
                                 if(MACMBRS(jb, m)/=0.0)                       &
                                  & dlmr = (MACMBRT(jb, m) - MACMBRS(jb, m))   &
                                  & /(MACMBRS(jb, m) + nonzero)
                                 write(SNP(jw), 9039)m, MACMBRS(jb, m),        &
                                     & MACMBRT(jb, m),                         &
                                     & (MACMBRT(jb, m) - MACMBRS(jb, m)),      &
                                     & dlmr*100.0
9039                             format(5x, 'Macrophyte spec ', i2, /7x,       &
                                      &'Spatially integrated mass [MACMBRS] = '&
                                     & , 1pe15.8E2, 1x, 'g ', /7x,             &
                                     &'Temporally integrated mass [MACMBRT] = '&
                                    & , 1pe15.8E2, 1x, 'g ', /7x,              &
                                     &'Mass error                         = ', &
                                    & 1pe15.8E2, 1x, 'g ', /7x,                &
                                     &'Percent error                      = ', &
                                    & 1pe15.8E2, ' %')
                             endif
                         enddo
 
                     enddo
                 endif
                 write(SNP(jw), 9040)'Geometry', KTWB(jw), ELKT(jw)
9040             format(/1x, a/3x, 'Surface layer [KT] = ', i0/3x,             &
                       &'Elevation   [ELKT] =', f10.3, ' m')
                 write(SNP(jw), 9041)(jb, CUS(jb), jb = BS(jw), BE(jw))
9041             format(/3x, 'Current upstream segment [CUS]'/(3x, '  Branch ',&
                      & i0, ' = ', i0))
                 call OUTPUT(jday, IBPR(jw), NISNP(jw), KBR(jw), ISNP,         &
                           & BL(1, jw), NBL(jw))
          !!SP CEMA    Code moved on 5/25/2017 SW
          !if(sediment_diagenesis)then
          !Call WriteCEMASedimentModelOutput
          !Call WriteCEMASedimentFluxOutput
          !end if
          !!End SP CEMA
             endif
         endif
 
    !**** Vertical profiles
 
         if(PROFILE(jw))then
             if(jday>=NXTMPR(jw) .OR. jday>=PRFD(PRFDP(jw) + 1, jw))then
                 if(jday>=PRFD(PRFDP(jw) + 1, jw))then
                     PRFDP(jw) = PRFDP(jw) + 1
                     NXTMPR(jw) = PRFD(PRFDP(jw), jw)
                 endif
                 NXTMPR(jw) = NXTMPR(jw) + PRFF(PRFDP(jw), jw)
                 NSPRF(jw) = NSPRF(jw) + 1
                 if(IPRF(1, 1)/= - 1)then
                                    ! SW 4/1/2016
                     write(PRF(jw), '(F8.3,1X,A3,I3,A,2I4,F8.4,I8)')jday,      &
                         & ADJUSTL(month), gday, ', ', year, KTWB(jw),         &
                         & SNGL(Z(DS(BS(jw)))), NSPRF(jw)
                     do jp = 1, NIPRF(jw)
                         nrs = KB(IPRF(jp, jw)) - KTWB(jw) + 1
                         write(PRF(jw), '(A8,I4/(8F10.2))')'TEMP    ', nrs,    &
                             & (T2(k, IPRF(jp, jw)), k = KTWB(jw),             &
                             & KB(IPRF(jp, jw)))
                     enddo
                     do jc = 1, nac
                         if(PRINT_CONST(CN(jc), jw))then
                             do jp = 1, NIPRF(jw)
                                 nrs = KB(IPRF(jp, jw)) - KTWB(jw) + 1
                                 write(PRF(jw), '(A8,I4/(8(E13.6,2X)))')       &
                                     & ADJUSTL(CNAME2(CN(jc))), nrs,           &
                                     & (C2(k, IPRF(jp, jw), CN(jc))            &
                                     & *CMULT(CN(jc)), k = KTWB(jw),           &
                                     & KB(IPRF(jp, jw)))                                                                         !     ! CB 1/24/05
                             enddo
                         endif
                     enddo
                     if(constituents)then
                         do jd = 1, NACD(jw)
                             do jp = 1, NIPRF(jw)
                                 nrs = KB(IPRF(jp, jw)) - KTWB(jw) + 1
                                 write(PRF(jw), '(A8,I4/(8(E13.6,2X)))')       &
                                     & ADJUSTL(CDNAME2(CDN(jd, jw))), nrs,     &
                                     & (CD(k, IPRF(jp, jw), CDN(jd, jw))       &
                                     & *CDMULT(CDN(jd, jw)), k = KTWB(jw),     &
                                     & KB(IPRF(jp, jw)))                                                                         !      ! CB 1/24/05
                             enddo
                         enddo
                     endif
 
                 elseif(jw==1)then
                              ! write out individual files on these days
                     write(segnum, '(F8.2)')jday    !   '(I0)'    INT(JDAY)
                     segnum = ADJUSTL(segnum)
                     l = LEN_TRIM(segnum)
                     open(nunit, file = 'ProfLongJD' // segnum(1:l) // '.csv', &
                        & status = 'UNKNOWN')
                     if(constituents)then
                         write(nunit, '(1000(A,","))')'Seg#',                  &
                              &'ElevWaterSurf(m)', 'Q(m3/s)', 'Temp(oC)',      &
                              &'Depth(m)', 'Width(m)',                         &
                             & (CNAME2(CN(jc)), jc = 1, nac),                  &
                             & (CDNAME2(CDN(jd, jw)), jd = 1, NACD(jw))
                     else
                         write(nunit, '(1000(A,","))')'Seg#',                  &
                              &'ElevWaterSurf(m)', 'Q(m3/s)', 'Temp(oC)',      &
                              &'Depth(m)', 'Width(m)'
                     endif
 
                     do jj = 1, nwb
                         k = KTWB(jj)
                         do jb = BS(jj), BE(jj)
                             do i = CUS(jb), DS(jb)
                                 if(constituents)then
                                     write(nunit, '(I5,",",100(F12.3,","))')i, &
                                       & ELWS(i), QC(i), T2(k, i),             &
                                       & DEPTHB(KB(i), i), B(KTI(i), i),       &
                                       & (C2(k, i, CN(jac))*CMULT(CN(jac)),    &
                                       & jac = 1, nac),                        &
                                       & (CD(k, i, CDN(jd, jw))*CDMULT         &
                                       & (CDN(jd, jw)), jd = 1, NACD(jw))
                                 else
                                     write(nunit, '(I5,",",100(F12.3,","))')i, &
                                       & ELWS(i), QC(i), T2(k, i),             &
                                       & DEPTHB(KB(i), i), B(KTI(i), i)
                                 endif
                             enddo
                         enddo
                     enddo
                     close(nunit)
                 endif
                    ! end of iprf /= -1
 
 
 
 
 
 
             endif
         endif
 
    !**** Spreadsheet
 
         if(SPREADSHEET(jw))then
             if(jday>=NXTMSP(jw) .OR. jday>=SPRD(SPRDP(jw) + 1, jw))then
                 if(jday>=SPRD(SPRDP(jw) + 1, jw))then
                     SPRDP(jw) = SPRDP(jw) + 1
                     NXTMSP(jw) = SPRD(SPRDP(jw), jw)
                 endif
                 conv1 = blank1
                 NXTMSP(jw) = NXTMSP(jw) + SPRF(SPRDP(jw), jw)
                 do j = 1, NISPR(jw)
                     KBMAX(jw) = MAX(KB(ISPR(j, jw)), KBMAX(jw))
                     do k = KTWB(jw), KB(ISPR(j, jw))
                         write(conv1(k, j), '(F10.2)')T2(k, ISPR(j, jw))
                     enddo
                 enddo
                 do k = KTWB(jw), KBMAX(jw)
                     write(SPR(jw),                                            &
                          &'(A,",",2(F10.3,","),1000(F10.3,",",A,","))')       &
                          &'Temperature', jday, DEPTHM(k, DS(BS(jw))),         &
                         & (ELWS(ISPR(j, jw)) - DEPTHM(k, ISPR(j, jw)),        &
                         & conv1(k, j), j = 1, NISPR(jw))
                 enddo
                 do jc = 1, nac
                     if(PRINT_CONST(CN(jc), jw))then
                         do j = 1, NISPR(jw)
                             do k = KTWB(jw), KB(ISPR(j, jw))
                                 write(conv1(k, j), '(F10.3)')                 &
                                     & C2(k, ISPR(j, jw), CN(jc))*CMULT(CN(jc))                                                   ! SW 8/13/06
                             enddo
                         enddo
                         do k = KTWB(jw), KBMAX(jw)
                             write(SPR(jw),                                    &
                                &'(A38,",",2(F10.3,","),1000(F10.3,",",A,","))'&
                               & )CNAME3(CN(jc)), jday, DEPTHM(k, DS(BS(jw))), &
                                & (ELWS(ISPR(j, jw)) - DEPTHM(k, ISPR(j, jw)), &
                                & conv1(k, j), j = 1, NISPR(jw))
                         enddo
                     endif
                 enddo
                 if(constituents)then
                     do jd = 1, NACD(jw)
                         if(PRINT_DERIVED(CDN(jd, jw), jw))then
                             do j = 1, NISPR(jw)
                                 do k = KTWB(jw), KB(ISPR(j, jw))
                                     write(conv1(k, j), '(F10.3)')             &
                                       & CD(k, ISPR(j, jw), CDN(jd, jw))       &
                                       & *CDMULT(CDN(jd, jw))                                                                     ! SW 8/13/06
                                 enddo
                             enddo
                             do k = KTWB(jw), KBMAX(jw)
                                 write(SPR(jw),                                &
                                &'(A38,",",2(F10.3,","),1000(F10.3,",",A,","))'&
                               & )CDNAME3(CDN(jd, jw)), jday,                  &
                                & DEPTHM(k, DS(BS(jw))),                       &
                                & (ELWS(ISPR(j, jw)) - DEPTHM(k, ISPR(j, jw)), &
                                & conv1(k, j), j = 1, NISPR(jw))
                             enddo
                         endif
                     enddo
                 endif
             endif
         endif
 
    !**** Contours
 
         if(CONTOUR(jw))then
             if(jday>=NXTMCP(jw) .OR. jday>=CPLD(CPLDP(jw) + 1, jw))then
                 if(jday>=CPLD(CPLDP(jw) + 1, jw))then
                     CPLDP(jw) = CPLDP(jw) + 1
                     NXTMCP(jw) = CPLD(CPLDP(jw), jw)
                 endif
                 NXTMCP(jw) = NXTMCP(jw) + CPLF(CPLDP(jw), jw)
 
                 if(TECPLOT(jw)/='      ON')then
                     write(CPL(jw), '(A,F12.4,5X,A9,5X,I2,5X,I4)')'New date ', &
                         & jday, month, gday, year
                     write(CPL(jw), '(9(I8,2X))')KTWB(jw)
                     write(CPL(jw), '(9(E13.6,2X))')(QTR(jt), jt = 1, ntr)
                     write(CPL(jw), '(9(E13.6,2X))')(TTR(jt), jt = 1, ntr)
                     do jt = 1, ntr
                         do jac = 1, NACTR(jt)
                             if(PRINT_CONST(TRCN(jac, jt), jw))                &
                               &write(CPL(jw), '(9(E13.6,2X))')                &
                              & CTR(TRCN(jac, jt), jt)
                         enddo
                     enddo
                     do jb = BS(jw), BE(jw)
                         write(CPL(jw), '(9(I8,2X))')CUS(jb)
                         write(CPL(jw), '(9(E13.6,2X))')QIN(jb), QSUM(jb)
                         do i = CUS(jb), DS(jb)
                             write(CPL(jw), '(A38/(9(E13.6,2X)))')'BHR',       &
                                 & (BHR1(k, i), k = KTWB(jw) + 1, KB(i))
                         enddo
                         do i = CUS(jb), DS(jb)
                             write(CPL(jw), '(A38/(9(E13.6,2X)))')'U',         &
                                 & (U(k, i), k = KTWB(jw), KB(i))
                         enddo
                         write(CPL(jw), '(A38/(9(E13.6,2X)))')'QC',            &
                             & (QC(i), i = CUS(jb), DS(jb))
                         write(CPL(jw), '(A38/(9(E13.6,2X)))')'Z',             &
                             & (Z(i), i = CUS(jb), DS(jb))
                         write(CPL(jw), '(A38/(9(I8,2X)))')'KTI',              &
                             & (KTI(i), i = CUS(jb), DS(jb))                              ! v3.5
                         do i = CUS(jb), DS(jb)
                             write(CPL(jw), '(A38/(9(E13.6,2X)))')             &
                                 & 'Temperature',                              &
                                 & (T2(k, i), k = KTWB(jw), KB(i))
                         enddo
                         do jc = 1, nac
                             if(PRINT_CONST(CN(jc), jw))then
                                 do i = CUS(jb), DS(jb)
                                     write(CPL(jw), '(A38/(9(E13.6,2X)))')     &
                                       & CNAME(CN(jc)),                        &
                                       & (C2(k, i, CN(jc))*CMULT(CN(jc)),      &
                                       & k = KTWB(jw), KB(i))
                                 enddo
                             endif
                         enddo
                         do je = 1, nep
                             do i = CUS(jb), DS(jb)
                                 if(PRINT_EPIPHYTON(jw, je))                   &
                                   &write(CPL(jw), '(A38/(9(E13.6,2X)))')      &
                                   &'Epiphyton',                               &
                                  & (EPD(k, i, je), k = KTWB(jw), KB(i))
                             enddo
                         enddo
                         if(PRINT_SEDIMENT(jw))then
                             do i = CUS(jb), DS(jb)
                                 write(CPL(jw), '(A38/(9(E13.6,2X)))')         &
                                      &'Sediment',                             &
                                     & (SED(k, i), k = KTWB(jw), KB(i))
                             enddo
                             do i = CUS(jb), DS(jb)
                                 write(CPL(jw), '(A38/(9(E13.6,2X)))')         &
                                      &'Sediment P',                           &
                                     & (SEDP(k, i), k = KTWB(jw), KB(i))
                             enddo
                             do i = CUS(jb), DS(jb)
                                 write(CPL(jw), '(A38/(9(E13.6,2X)))')         &
                                      &'Sediment N',                           &
                                     & (SEDN(k, i), k = KTWB(jw), KB(i))
                             enddo
                             do i = CUS(jb), DS(jb)
                                 write(CPL(jw), '(A38/(9(E13.6,2X)))')         &
                                      &'Sediment C',                           &
                                     & (SEDC(k, i), k = KTWB(jw), KB(i))
                             enddo
                         endif
                         do m = 1, nmc
                             if(PRINT_MACROPHYTE(jw, m))then
                                 do i = CUS(jb), DS(jb)
                                     write(CPL(jw), '(A38/(9(E13.6,2X)))')     &
                                        &'Macrophytes',                        &
                                       & ((MACRC(j, k, i, m), j = KTI(i),      &
                                       & KB(i)), k = KTWB(jw), KB(i))
                                 enddo
                             endif
                         enddo
                         if(constituents)then
                             do jd = 1, NACD(jw)
                                 if(PRINT_DERIVED(CDN(jd, jw), jw))then
                                     do i = CUS(jb), DS(jb)
                                         write(CPL(jw), '(A38/(9(F10.3,2X)))') &
                                           & CDNAME(CDN(jd, jw)),              &
                                           & (CD(k, i, CDN(jd, jw))            &
                                           & *CDMULT(CDN(jd, jw)), k = KTWB(jw)&
                                           & , KB(i))                                                                                       ! cb 6/28/13
                                     enddo
                  !WRITE (CPL(JW),'(A38/(9(F10.3,2X)))') CDNAME(CDN(JD,JW)),((CD(K,I,CDN(JD,JW))*CDMULT(CDN(JD,JW)),             &        ! SW 8/12/06
                  !K=KTWB(JW),KB(I)),I=CUS(JB),DS(JB))  ! CB 1/03/05
                                 endif
                             enddo
                         endif
                     enddo
                 else
          !         ICPL=ICPL+1
                     ICPL(jw) = ICPL(jw) + 1
                               ! cb 1/26/09
                     itot = 0
          !         do jb=1,nbr
                     do jb = BS(jw), BE(jw)
                               ! cb 1/26/09
                         itot = itot + DS(jb) - CUS(jb) + 2
                     enddo
                     write(CPL(jw), 9042)jday, kmx - KTWB(jw) + 2, itot
9042                 format('ZONE T="', f9.3, '"', ' I=', i3, ' J=', i3,       &
                           &' F=POINT')
                     do jb = BS(jw), BE(jw)
                         do i = CUS(jb), DS(jb) + 1
                             k = KTWB(jw)
                          ! PRINT AN EXTRA LINE FOR THE SURFACE
                             if(i==DS(jb) + 1)then
                                 xdum = -99.0
                                 if(habtatc=='      ON')then                                   ! SW 7/15/14
                                     write(CPL(jw), 9045)X1(i), ELWS(i), xdum, &
                                       & xdum, xdum, xdum, xdum,               &
                                       & (xdum, jj = 1, nac),                  &
                                       & (xdum, jj = 1, NACD(jw))
                                 else                                                          ! SW 7/15/14
                                     write(CPL(jw), 9045)X1(i), ELWS(i), xdum, &
                                       & xdum, xdum, xdum, (xdum, jj = 1, nac),&
                                       & (xdum, jj = 1, NACD(jw))                                                      ! SW 1/17/17 7/15/14
                                 endif                                                         ! SW 7/15/14
                             elseif(habtatc=='      ON')then
                                 write(CPL(jw), 9045)X1(i), ELWS(i), U(k, i),  &
                                     & -W(k, i), T1(k, i), RHO(k, i), HAB(k, i)&
                                     & , (C2(k, i, CN(jc)), jc = 1, nac),      &
                                     & (CD(k, i, CDN(jd, jw)), jd = 1, NACD(jw)&
                                     & )                                                                                                                     ! SW 1/17/17
                             else
                                 write(CPL(jw), 9045)X1(i), ELWS(i), U(k, i),  &
                                     & -W(k, i), T1(k, i), RHO(k, i),          &
                                     & (C2(k, i, CN(jc)), jc = 1, nac),        &
                                     & (CD(k, i, CDN(jd, jw)), jd = 1, NACD(jw)&
                                     & )                                                                                                                     ! SW 1/17/17
 
                             endif
                             do k = KTWB(jw), kmx - 1
                                 if(i/=DS(jb) + 1 .AND. k<=KB(i))then
                                     if(habtatc=='      ON')then
                                         write(CPL(jw), 9045)X1(i), ELWS(i)    &
                                           & - DEPTHM(k, i), U(k, i), -W(k, i),&
                                           & T1(k, i), RHO(k, i), HAB(k, i),   &
                                           & (C2(k, i, CN(jc)), jc = 1, nac),  &
                                           & (CD(k, i, CDN(jd, jw)), jd = 1,   &
                                           & NACD(jw))                                                                                                                   ! SW 1/17/17
                                     else
                                         write(CPL(jw), 9045)X1(i), ELWS(i)    &
                                           & - DEPTHM(k, i), U(k, i), -W(k, i),&
                                           & T1(k, i), RHO(k, i),              &
                                           & (C2(k, i, CN(jc)), jc = 1, nac),  &
                                           & (CD(k, i, CDN(jd, jw)), jd = 1,   &
                                           & NACD(jw))                                                                                                          ! SW 1/17/17
                                     endif
 
                                     if(k==KB(i))then
                                         if(habtatc=='      ON')then
                                         write(CPL(jw), 9045)X1(i), ELWS(i)    &
                                           & - DEPTHB(k, i), U(k, i), -W(k, i),&
                                           & T1(k, i), RHO(k, i), HAB(k, i),   &
                                           & (C2(k, i, CN(jc)), jc = 1, nac),  &
                                           & (CD(k, i, CDN(jd, jw)), jd = 1,   &
                                           & NACD(jw))                                                                                                                      ! SW 1/17/17
                                         else
                                         write(CPL(jw), 9045)X1(i), ELWS(i)    &
                                           & - DEPTHB(k, i), U(k, i), -W(k, i),&
                                           & T1(k, i), RHO(k, i),              &
                                           & (C2(k, i, CN(jc)), jc = 1, nac),  &
                                           & (CD(k, i, CDN(jd, jw)), jd = 1,   &
                                           & NACD(jw))                                                                                                             ! SW 1/17/17
                                         endif
 
                                     endif
                                 else
                                     xdum = -99.0
                                     if(habtatc=='      ON')then
                                         write(CPL(jw), 9045)X1(i), ELWS(i - 1)&
                                           & - DEPTHM(k, i - 1), xdum, xdum,   &
                                           & xdum, xdum, xdum,                 &
                                           & (xdum, jj = 1, nac),              &
                                           & (xdum, jj = 1, NACD(jw))
                                     else
                                         write(CPL(jw), 9045)X1(i), ELWS(i - 1)&
                                           & - DEPTHM(k, i - 1), xdum, xdum,   &
                                           & xdum, xdum, (xdum, jj = 1, nac),  &
                                           & (xdum, jj = 1, NACD(jw))                                                                ! SW 7/15/14
                                     endif
 
                                     if(k==KB(i))then
                                         if(habtatc=='      ON')then
                                         write(CPL(jw), 9045)X1(i), ELWS(i - 1)&
                                           & - DEPTHB(k, i - 1), xdum, xdum,   &
                                           & xdum, xdum, xdum,                 &
                                           & (xdum, jj = 1, nac),              &
                                           & (xdum, jj = 1, NACD(jw))
                                         else
                                         write(CPL(jw), 9045)X1(i), ELWS(i - 1)&
                                           & - DEPTHB(k, i - 1), xdum, xdum,   &
                                           & xdum, xdum, (xdum, jj = 1, nac),  &
                                           & (xdum, jj = 1, NACD(jw))                                                                     ! SW 7/15/14
                                         endif
 
                                     endif
                                 endif
                             enddo
                         enddo
                     enddo
                     write(CPL(jw), 9043)ICPL(jw), imon, gday, year
                                                      ! cb 1/26/09
9043                 format('TEXT X=0.75, y=0.85, H=2.8,ZN=', i4, ',',         &
                           &' C=BLACK,', 'T= "', i2, '/', i2, '/', i4, '"')
                     write(CPL(jw), 9044)ICPL(jw), jday
                                            ! cb 1/26/09
9044                 format('TEXT X=0.75, y=0.90, H=2.8,ZN=', i4, ',',         &
                           &' C=BLACK,', 'T= "Julian Day ', f9.3, '"')
                 endif
             endif
         endif
 
    !**** Fluxes   KF is the instantaneous flux in g/m3/s, KFS is the summed flux in g eventually divided by elapsed time between calls to FLUX output and converted to kg below
 
 
         if(FLUX(jw))then
             if(jday>=NXTMFL(jw) .OR. jday>=FLXD(FLXDP(jw) + 1, jw))then
                 if(jday>=FLXD(FLXDP(jw) + 1, jw))then
                     FLXDP(jw) = FLXDP(jw) + 1
                     NXTMFL(jw) = FLXD(FLXDP(jw), jw)
                 endif
 
                 if(CONTOUR(jw) .AND. derived_calc)then
                     do jb = BS(jw), BE(jw)
                         do i = CUS(jb), DS(jb)
                             do k = KTWB(jw), KB(i)
                      !   TP_SEDBURIAL(JW) = TP_SEDBURIAL(JW) +KFS(K,I,KF_SED_PBURIAL)/1000.
                      !   TN_SEDBURIAL(JW) = TN_SEDBURIAL(JW) +KFS(K,I,KF_SED_NBURIAL)/1000.
                                 TN_SEDSOD_NH4(jw) = TN_SEDSOD_NH4(jw)         &
                                   & + KFS(k, i, kf_nh4_sod)                   &
                                   & /1000. + KFS(k, i, kf_nh4_sed)/1000.
                                 TP_SEDSOD_PO4(jw) = TP_SEDSOD_PO4(jw)         &
                                   & + KFS(k, i, kf_po4_sod)                   &
                                   & /1000. + KFS(k, i, kf_po4_sed)/1000.
                             enddo
                         enddo
                     enddo
                 endif
 
 
                 nlines = 0
                    ! SW 3/8/16
                 NXTMFL(jw) = NXTMFL(jw) + FLXF(FLXDP(jw), jw)
                 conv = blank
                 do jaf = 1, NAF(jw)
                     do jb = BS(jw), BE(jw)
                         do i = CUS(jb), DS(jb)
                             do k = KTWB(jw), KB(i)
                                 KFS(k, i, KFCN(jaf, jw))                      &
                                   & = day*KFS(k, i, KFCN(jaf, jw))            &
                                   & /(1000.*ELTMF(jw))                              ! KFS IN G, 86400 S/D * G /ELAPSED TIME IN S/1000 G/KG == KG/D
                                 KFJW(jw, KFCN(jaf, jw))                       &
                                   & = KFJW(jw, KFCN(jaf, jw))                 &
                                   & + KFS(k, i, KFCN(jaf, jw))                      ! SUM UP FOR ENTIRE WATERBODY
                             enddo
                         enddo
                     enddo
                     do i = 1, NISNP(jw)
                         do k = KTWB(jw), KB(ISNP(i, jw))
                             write(conv(k, i), '(E10.3)')                      &
                                 & KFS(k, ISNP(i, jw), KFCN(jaf, jw))         ! KG/D
                         enddo
                     enddo
                     if(new_page)then
                         write(FLX(jw), '(/(A72))')(TITLE(j), j = 1, 11)
                         nlines = kmx - KTWB(jw) + 14
                         new_page = .FALSE.
                     endif
                     nlines = nlines + kmx - KTWB(jw) + 11
                     new_page = nlines>72
                     write(FLX(jw), '(/A,F10.3,X,3(A,I0),A,F0.2,A/)')          &
                          &'New date ', jday, month // ' ', gday, ', ', year,  &
                          &'   Julian Date = ', INT(jday), ' days ',           &
                         & (jday - INT(jday))*24.0, ' hours           ' //     &
                         & KFNAME(KFCN(jaf, jw))
                     write(FLX(jw), '(3X,2000I10)')                            &
                         & (ISNP(i, jw), i = 1, NISNP(jw))
                     do k = KTWB(jw), KBR(jw)
                         write(FLX(jw), '(1X,I2,200A)')k,                      &
                             & (conv(k, i), i = 1, NISNP(jw))
                     enddo
                 enddo
                 write(FLX2(jw), '(F10.3,",",f8.3,",",1000(E12.4,","))')jday,  &
                     & ELTMF(jw)/day, (KFJW(jw, KFCN(k, jw)), k = 1, NAF(jw))
                 ELTMF(jw) = 0.0
                 KF(:, CUS(BS(jw)):DS(BE(jw)), KFCN(1:NAF(jw), jw)) = 0.0
                 KFS(:, CUS(BS(jw)):DS(BE(jw)), KFCN(1:NAF(jw), jw)) = 0.0
                 KFJW(jw, KFCN(1:NAF(jw), jw)) = 0.0
             endif
         endif
 
     enddo
 
  !** Downstream flow, temperature, and constituent files
 
     if(downstream_outflow)then
         if(jday>=nxtmwd .OR. jday>=WDOD(wdodp + 1))then
             if(jday>=WDOD(wdodp + 1))then
                 wdodp = wdodp + 1
                 nxtmwd = WDOD(wdodp)
             endif
             if(wdoc=='      ON')then
                 nxtmwd = nxtmwd + WDOF(wdodp)
             elseif(wdoc=='     ONS')then
                 nxtmwd_sec = nxtmwd_sec + WDOF(wdodp)
                                             ! cb 4/6/18 frequency test seconds
                 nxtmwd = nxtmwd_sec/86400.0
             elseif(wdoc=='     ONH')then
                 nxtmwd_sec = nxtmwd_sec + WDOF(wdodp)*3600.
                                                   ! cb 4/6/18 frequency test hourly
                 nxtmwd = nxtmwd_sec/86400.0
             endif
 
             jfile = 0
 
             do j = 1, niwdo
                 QWDO(j) = 0.0
                 TWDO(j) = 0.0
                 cwdo(:, j) = 0.0
                 cdwdo(:, j) = 0.0
                 cdtot = 0.0
                 numoutlets = 0
                 jwd = 0           ! SW 5/17/13
 
                 do jw = 1, nwb
                     do jb = BS(jw), BE(jw)
                         if(DS(jb)==IWDO(j))then
 
                             do js = 1, JSS(jb)
                               !NSTR(JB)
                                 numoutlets = numoutlets + 1
                                 if(QSTR(js, jb)==0.0)then
                                     TAVG(js, jb) = -99.0
                                     cavg(js, jb, :) = -99.0
                                     cdavg(js, jb, :) = -99.0
                                 endif
                                 qoutlet(numoutlets) = QSTR(js, jb)
                                 toutlet(numoutlets) = TAVG(js, jb)
                             enddo
 
              ! cb 1/16/13 removed old code
 
              ! OUTPUT INDIVIDUAL FILES
                             do js = 1, NSTR(jb)
                                 jfile = jfile + 1
                                 QWDO(j) = QWDO(j) + QSTR(js, jb)                   ! cb 1/16/13
                                 TWDO(j) = TWDO(j) + QSTR(js, jb)*TAVG(js, jb)
                                 write(WDO2(jfile, 1), '(F10.3,",",F10.4)')    &
                                     & jday, QSTR(js, jb)
                                 write(WDO2(jfile, 2), '(F10.3,",",F8.2)')jday,&
                                     & TAVG(js, jb)
                                 if(constituents)then
                                     cwdo(CN(1:nac), j) = cwdo(CN(1:nac), j)   &
                                       & + QSTR(js, jb)*cavg(js, jb, CN(1:nac))                          ! cb 1/16/13
                                     write(WDO2(jfile, 3),                     &
                                        &'(F10.3,",",1000(F10.4,","))')jday,   &
                                       & (cavg(js, jb, CN(jc)), jc = 1, nac)
                                 endif
                                 if(derived_calc)then
                                     cdwdo(CDN(1:NACD(jw), jw), j)             &
                                       & = cdwdo(CDN(1:NACD(jw), jw), j)       &
                                       & + QSTR(js, jb)                        &
                                       & *cdavg(js, jb, CDN(1:NACD(jw), jw))                                                            ! cb 1/16/13
                                     write(WDO2(jfile, 4),                     &
                                        &'(F10.3,",",1000(F10.4,","))')jday,   &
                                       & (cdavg(js, jb, CDN(jd, jw)), jd = 1,  &
                                       & NACD(jw))
                                 endif
                             enddo
                             jsss(jb) = NSTR(jb)
                         endif
                     enddo
                 enddo
 
        ! OUTPUT INDIVIDUAL FILES
        ! Order: spillways NSP, pumps NPU, gates NGT, pipes NPI
   !     JWD=0                        ! sw 5/17/13
                 do jw = 1, nwb       ! cb 1/16/13
                     if(IWDO(j)>=US(BS(jw)) .AND. IWDO(j)<=DS(BE(jw)))exit
                 enddo
                 do jj = 1, nwd
                     jwd = jwd + 1
                     if(QWD(jwd)==0.0)then
                         TAVGW(jwd) = -99.0
                         cavgw(jwd, :) = -99.0
                         cdavgw(jwd, :) = -99.0
                     endif
                     if(IWD(jwd)==IWDO(j))then
                         jfile = jfile + 1
                         QWDO(j) = QWDO(j) + QWD(jwd)                        ! cb 1/16/13
                         TWDO(j) = TWDO(j) + QWD(jwd)*TAVGW(jwd)                        ! cb 1/16/13
                         write(WDO2(jfile, 1), '(F10.3,",",F10.4)')jday,       &
                             & QWD(jwd)
                         write(WDO2(jfile, 2), '(F10.3,",",F8.2)')jday,        &
                             & TAVGW(jwd)
                         if(constituents)then
                             cwdo(CN(1:nac), j) = cwdo(CN(1:nac), j) + QWD(jwd)&
                               & *cavgw(jwd, CN(1:nac))                                          ! cb 1/16/13
                             write(WDO2(jfile, 3),                             &
                                  &'(F10.3,",",1000(F10.4,","))')jday,         &
                                 & (cavgw(jwd, CN(jc)), jc = 1, nac)
                         endif
                         if(derived_calc)then
              !DO JW=1,NWB
              !  IF (IWDO(J) >= US(BS(JW)) .AND. IWDO(J) <= DS(BE(JW))) EXIT
              !END DO
                             cdwdo(CDN(1:NACD(jw), jw), j)                     &
                               & = cdwdo(CDN(1:NACD(jw), jw), j) + QWD(jwd)    &
                               & *cdavgw(jwd, CDN(1:NACD(jw), jw))                                                              ! cb 1/16/13
                             write(WDO2(jfile, 4),                             &
                                  &'(F10.3,",",1000(F10.4,","))')jday,         &
                                 & (cdavgw(jwd, CDN(jd, jw)), jd = 1, NACD(jw))
                         endif
                     endif
                 enddo
                 do js = 1, nsp
                     ! spillways
                     if(LATERAL_SPILLWAY(js))then
                         jwd = jwd + 1
                     else
                         jsss(JBUSP(js)) = jsss(JBUSP(js)) + 1
                     endif
 
                     if(IWDO(j)==IUSP(js))then
                         jfile = jfile + 1
                         write(WDO2(jfile, 1), '(F10.3,",",F10.4)')jday,       &
                             & QSP(js)
                         if(LATERAL_SPILLWAY(js))then
            !  JWD=JWD+1
                             QWDO(j) = QWDO(j) + QSP(js)                      ! cb 1/16/13
                             TWDO(j) = TWDO(j) + QSP(js)*TAVGW(jwd)                      ! cb 1/16/13
                             if(QSP(js)==0.0)then
                                 TAVGW(jwd) = -99.0
                                 cavgw(jwd, :) = -99.0
                                 cdavgw(jwd, :) = -99.0
                             endif
                             write(WDO2(jfile, 2), '(F10.3,",",F8.2)')jday,    &
                                 & TAVGW(jwd)
                             if(constituents)then
                                 cwdo(CN(1:nac), j) = cwdo(CN(1:nac), j)       &
                                   & + QSP(js)*cavgw(jwd, CN(1:nac))                              ! cb 1/16/13
                                 write(WDO2(jfile, 3),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cavgw(jwd, CN(jc)), jc = 1, nac)
                             endif
                             if(derived_calc)then
                                 cdwdo(CDN(1:NACD(jw), jw), j)                 &
                                   & = cdwdo(CDN(1:NACD(jw), jw), j) + QSP(js) &
                                   & *cdavgw(jwd,                              &
                                   & CDN(1:NACD(JWUSP(js)), JWUSP(js)))                                                               ! cb 1/16/13
                                 write(WDO2(jfile, 4),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cdavgw(jwd, CDN(jd, JWUSP(js))),       &
                                     & jd = 1, NACD(JWUSP(js)))
                             endif
                         else
            !  JSSS(JBUSP(JS))=JSSS(JBUSP(JS))+1
                             QWDO(j) = QWDO(j) + QSP(js)                        ! cb 1/16/13
                             TWDO(j) = TWDO(j) + QSP(js)                       &
                                     & *TAVG(jsss(JBUSP(js)), JBUSP(js))                           ! cb 1/16/13
                             if(QSP(js)==0.0)then
                                 TAVG(jsss(JBUSP(js)), JBUSP(js)) = -99.0
                                 cavg(jsss(JBUSP(js)), JBUSP(js), :) = -99.0
                                 cdavg(jsss(JBUSP(js)), JBUSP(js), :) = -99.0
                             endif
                             write(WDO2(jfile, 2), '(F10.3,",",F8.2)')jday,    &
                                 & TAVG(jsss(JBUSP(js)), JBUSP(js))
                             if(constituents)then
                                 cwdo(CN(1:nac), j) = cwdo(CN(1:nac), j)       &
                                   & + QSP(js)                                 &
                                   & *cavg(jsss(JBUSP(js)), JBUSP(js), CN      &
                                   & (1:nac))                                                                     ! cb 1/16/13
                                 write(WDO2(jfile, 3),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cavg(jsss(JBUSP(js)), JBUSP(js), CN(jc)&
                                     & ), jc = 1, nac)
                             endif
                             if(derived_calc)then
                                 cdwdo(CDN(1:NACD(jw), jw), j)                 &
                                   & = cdwdo(CDN(1:NACD(jw), jw), j) + QSP(js) &
                                   & *cdavg(jsss(JBUSP(js)), JBUSP(js),        &
                                   & CDN(1:NACD(JWUSP(js)), JWUSP(js)))                                                                                ! cb 1/16/13
                                 write(WDO2(jfile, 4),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cdavg(jsss(JBUSP(js)), JBUSP(js),      &
                                     & CDN(jd, JWUSP(js))), jd = 1,            &
                                     & NACD(JWUSP(js)))
                             endif
                         endif
                     endif
                 enddo
                 do js = 1, npu
                     ! PUMP
                     if(LATERAL_PUMP(js))then
                         jwd = jwd + 1
                     else
                         jsss(JBUPU(js)) = jsss(JBUPU(js)) + 1
                     endif
                     if(IWDO(j)==IUPU(js))then
                         jfile = jfile + 1
                         if(PUMPON(js))then
                             write(WDO2(jfile, 1), '(F10.3,",",F10.4)')jday,   &
                                 & QPU(js)
                         else
                             write(WDO2(jfile, 1), '(F10.3,",",F8.3)')jday, 0.0
                         endif
                         if(LATERAL_PUMP(js))then
            !  JWD=JWD+1
                             if(QPU(js)==0.0)then
                                 TAVGW(jwd) = -99.0
                                 cavgw(jwd, :) = -99.0
                                 cdavgw(jwd, :) = -99.0
                             endif
                             if(PUMPON(js))then    ! cb 1/16/13
                                 QWDO(j) = QWDO(j) + QPU(js)
                                 TWDO(j) = TWDO(j) + QPU(js)*TAVGW(jwd)
                             endif
                             write(WDO2(jfile, 2), '(F10.3,",",F8.2)')jday,    &
                                 & TAVGW(jwd)
              ! Debug
              !WRITE(WDO2(JFILE,2),'(F10.3,",",F8.2,",",f8.3,",",i5,",",i5)')JDAY,TAVGW(JWD),qpu(js),js,jwd      ! Debug
                             if(constituents)then
                                 if(PUMPON(js))cwdo(CN(1:nac), j)              &
                                  & = cwdo(CN(1:nac), j) + QPU(js)             &
                                  & *cavgw(jwd, CN(1:nac))
                                                   ! cb 1/16/13
                                 write(WDO2(jfile, 3),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cavgw(jwd, CN(jc)), jc = 1, nac)
                             endif
                             if(derived_calc)then
                                 if(PUMPON(js))cdwdo(CDN(1:NACD(jw), jw), j)   &
                                  & = cdwdo(CDN(1:NACD(jw), jw), j) + QPU(js)  &
                                  & *cdavgw(jwd,                               &
                                  & CDN(1:NACD(JWUPU(js)), JWUPU(js)))
                                     ! cb 1/16/13
                                 write(WDO2(jfile, 4),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cdavgw(jwd, CDN(jd, JWUPU(js))),       &
                                     & jd = 1, NACD(JWUPU(js)))
                             endif
                         else
            !  JSSS(JBUPU(JS))=JSSS(JBUPU(JS))+1
                             if(QPU(js)==0.0)then
                                 TAVG(jsss(JBUPU(js)), JBUPU(js)) = -99.0
                                 cavg(jsss(JBUPU(js)), JBUPU(js), :) = -99.0
                                 cdavg(jsss(JBUPU(js)), JBUPU(js), :) = -99.0
                             endif
                             if(PUMPON(js))then  ! cb 1/16/13
                                 QWDO(j) = QWDO(j) + QPU(js)
                                 TWDO(j) = TWDO(j) + QPU(js)                   &
                                   & *TAVG(jsss(JBUPU(js)), JBUPU(js))
                             endif
                             write(WDO2(jfile, 2), '(F10.3,",",F8.2)')jday,    &
                                 & TAVG(jsss(JBUPU(js)), JBUPU(js))
                             if(constituents)then
                                 if(PUMPON(js))cwdo(CN(1:nac), j)              &
                                  & = cwdo(CN(1:nac), j) + QPU(js)             &
                                  & *cavg(jsss(JBUPU(js)), JBUPU(js), CN(1:nac)&
                                  & )             ! cb 1/16/13
                                 write(WDO2(jfile, 3),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cavg(jsss(JBUPU(js)), JBUPU(js), CN(jc)&
                                     & ), jc = 1, nac)
                             endif
                             if(derived_calc)then
                                 if(PUMPON(js))cdwdo(CDN(1:NACD(jw), jw), j)   &
                                  & = cdwdo(CDN(1:NACD(jw), jw), j) + QPU(js)  &
                                  & *cdavg(jsss(JBUPU(js)), JBUPU(js),         &
                                  & CDN(1:NACD(JWUPU(js)), JWUPU(js)))
                                     ! cb 1/16/13
                                 write(WDO2(jfile, 4),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cdavg(jsss(JBUPU(js)), JBUPU(js),      &
                                     & CDN(jd, JWUPU(js))), jd = 1,            &
                                     & NACD(JWUPU(js)))
                             endif
                         endif
                     endif
                 enddo
                 do js = 1, npi
                     ! pipes
                     if(LATERAL_PIPE(js))then
                         jwd = jwd + 1
                     else
                         jsss(JBUPI(js)) = jsss(JBUPI(js)) + 1
                     endif
                     if(IWDO(j)==IUPI(js))then
                         jfile = jfile + 1
                         write(WDO2(jfile, 1), '(F10.3,",",F10.4)')jday,       &
                             & QPI(js)
                         if(LATERAL_PIPE(js))then
           !   JWD=JWD+1
                             QWDO(j) = QWDO(j) + QPI(js)            ! cb 1/16/13
                             TWDO(j) = TWDO(j) + QPI(js)*TAVGW(jwd)
                             if(QPI(js)==0.0)then
                                 TAVGW(jwd) = -99.0
                                 cavgw(jwd, :) = -99.0
                                 cdavgw(jwd, :) = -99.0
                             endif
                             write(WDO2(jfile, 2), '(F10.3,",",F8.2)')jday,    &
                                 & TAVGW(jwd)
                             if(constituents)then
                                 cwdo(CN(1:nac), j) = cwdo(CN(1:nac), j)       &
                                   & + QPI(js)*cavgw(jwd, CN(1:nac))                  ! cb 1/16/13
                                 write(WDO2(jfile, 3),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cavgw(jwd, CN(jc)), jc = 1, nac)
                             endif
                             if(derived_calc)then
                                 cdwdo(CDN(1:NACD(jw), jw), j)                 &
                                   & = cdwdo(CDN(1:NACD(jw), jw), j) + QPI(js) &
                                   & *cdavgw(jwd,                              &
                                   & CDN(1:NACD(JWUPI(js)), JWUPI(js)))                                                            ! cb 1/16/13
                                 write(WDO2(jfile, 4),                         &
                                      &'(F10.3,1000(F10.4,","))')jday,         &
                                     & (cdavgw(jwd, CDN(jd, JWUPI(js))),       &
                                     & jd = 1, NACD(JWUPI(js)))
                             endif
                         else
            !  JSSS(JBUPI(JS))=JSSS(JBUPI(JS))+1
                             QWDO(j) = QWDO(j) + QPI(js)            ! cb 1/16/13
                             TWDO(j) = TWDO(j) + QPI(js)                       &
                                     & *TAVG(jsss(JBUPI(js)), JBUPI(js))
                             if(QPI(js)==0.0)then
                                 TAVG(jsss(JBUPI(js)), JBUPI(js)) = -99.0
                                 cavg(jsss(JBUPI(js)), JBUPI(js), :) = -99.0
                                 cdavg(jsss(JBUPI(js)), JBUPI(js), :) = -99.0
                             endif
                             write(WDO2(jfile, 2), '(F10.3,",",F8.2)')jday,    &
                                 & TAVG(jsss(JBUPI(js)), JBUPI(js))
                             if(constituents)then
                                 cwdo(CN(1:nac), j) = cwdo(CN(1:nac), j)       &
                                   & + QPI(js)                                 &
                                   & *cavg(jsss(JBUPI(js)), JBUPI(js), CN      &
                                   & (1:nac))                                                              ! cb 1/16/13
                                 write(WDO2(jfile, 3),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cavg(jsss(JBUPI(js)), JBUPI(js), CN(jc)&
                                     & ), jc = 1, nac)
                             endif
                             if(derived_calc)then
                                 cdwdo(CDN(1:NACD(jw), jw), j)                 &
                                   & = cdwdo(CDN(1:NACD(jw), jw), j) + QPI(js) &
                                   & *cdavg(jsss(JBUPI(js)), JBUPI(js),        &
                                   & CDN(1:NACD(JWUPI(js)), JWUPI(js)))                                                                                 ! cb 1/16/13
                                 write(WDO2(jfile, 4),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cdavg(jsss(JBUPI(js)), JBUPI(js),      &
                                     & CDN(jd, JWUPI(js))), jd = 1,            &
                                     & NACD(JWUPI(js)))
                             endif
                         endif
                     endif
                 enddo
                 do js = 1, ngt
                     ! gates
                     if(LATERAL_GATE(js))then
                         jwd = jwd + 1
                     else
                         jsss(JBUGT(js)) = jsss(JBUGT(js)) + 1
                     endif
 
                     if(IWDO(j)==IUGT(js))then
                         jfile = jfile + 1
                         write(WDO2(jfile, 1), '(F10.3,",",F10.4)')jday,       &
                             & QGT(js)
                         if(LATERAL_GATE(js))then
          !    JWD=JWD+1
                             QWDO(j) = QWDO(j) + QGT(js)          ! cb 1/16/13
                             TWDO(j) = TWDO(j) + QGT(js)*TAVGW(jwd)
                             if(QGT(js)==0.0)then
                                 TAVGW(jwd) = -99.0
                                 cavgw(jwd, :) = -99.0
                                 cdavgw(jwd, :) = -99.0
                             endif
                             write(WDO2(jfile, 2), '(F10.3,",",F8.2)')jday,    &
                                 & TAVGW(jwd)
                             if(constituents)then
                                 cwdo(CN(1:nac), j) = cwdo(CN(1:nac), j)       &
                                   & + QGT(js)*cavgw(jwd, CN(1:nac))                  ! cb 1/16/13
                                 write(WDO2(jfile, 3),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cavgw(jwd, CN(jc)), jc = 1, nac)
                             endif
                             if(derived_calc)then
                                 cdwdo(CDN(1:NACD(jw), jw), j)                 &
                                   & = cdwdo(CDN(1:NACD(jw), jw), j) + QGT(js) &
                                   & *cdavgw(jwd,                              &
                                   & CDN(1:NACD(JWUGT(js)), JWUGT(js)))                                                            ! cb 1/16/13
                                 write(WDO2(jfile, 4),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cdavgw(jwd, CDN(jd, JWUGT(js))),       &
                                     & jd = 1, NACD(JWUGT(js)))
                             endif
                         else
           !   JSSS(JBUGT(JS))=JSSS(JBUGT(JS))+1
                             QWDO(j) = QWDO(j) + QGT(js)           ! cb 1/16/13
                             TWDO(j) = TWDO(j) + QGT(js)                       &
                                     & *TAVG(jsss(JBUGT(js)), JBUGT(js))
                             if(QGT(js)==0.0)then
                                 TAVG(jsss(JBUGT(js)), JBUGT(js)) = -99.0
                                 cavg(jsss(JBUGT(js)), JBUGT(js), :) = -99.0
                                 cdavg(jsss(JBUGT(js)), JBUGT(js), :) = -99.0
                             endif
                             write(WDO2(jfile, 2), '(F10.3,",",F8.2)')jday,    &
                                 & TAVG(jsss(JBUGT(js)), JBUGT(js))
                             if(constituents)then
                                 cwdo(CN(1:nac), j) = cwdo(CN(1:nac), j)       &
                                   & + QGT(js)                                 &
                                   & *cavg(jsss(JBUGT(js)), JBUGT(js), CN      &
                                   & (1:nac))                                                              ! cb 1/16/13
                                 write(WDO2(jfile, 3),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cavg(jsss(JBUGT(js)), JBUGT(js), CN(jc)&
                                     & ), jc = 1, nac)
                             endif
                             if(derived_calc)then
                                 cdwdo(CDN(1:NACD(jw), jw), j)                 &
                                   & = cdwdo(CDN(1:NACD(jw), jw), j) + QGT(js) &
                                   & *cdavg(jsss(JBUGT(js)), JBUGT(js),        &
                                   & CDN(1:NACD(JWUGT(js)), JWUGT(js)))                                                                                 ! cb 1/16/13
                                 write(WDO2(jfile, 4),                         &
                                      &'(F10.3,",",1000(F10.4,","))')jday,     &
                                     & (cdavg(jsss(JBUGT(js)), JBUGT(js),      &
                                     & CDN(jd, JWUGT(js))), jd = 1,            &
                                     & NACD(JWUGT(js)))
                             endif
                         endif
                     endif
                 enddo
 
        ! cb 1/16/13 deleted old withdrawal output code
 
                 if(QWDO(j)/=0.0)TWDO(j) = TWDO(j)/QWDO(j)
                 do jc = 1, nac
                     if(QWDO(j)/=0.0)cwdo(CN(jc), j) = cwdo(CN(jc), j)/QWDO(j)
                     write(CWDOC(CN(jc)), '(F10.4)')cwdo(CN(jc), j)             ! SW 9/23/13 Changed format from G8.3 to F8.3 to avoid format overflow
                     CWDOC(CN(jc)) = ADJUSTR(CWDOC(CN(jc)))
                 enddo
                 do jw = 1, nwb
                     if(IWDO(j)>=US(BS(jw)) .AND. IWDO(j)<=DS(BE(jw)))exit
                 enddo
                 do jd = 1, NACD(jw)
                     if(QWDO(j)/=0.0)cdwdo(CDN(jd, jw), j)                     &
                      & = cdwdo(CDN(jd, jw), j)/QWDO(j)
                     write(CDWDOC(CDN(jd, jw)), '(F10.4)')cdwdo(CDN(jd, jw), j) ! SW 9/23/13 Changed format from G8.3 to F8.3 to avoid format overflow
                     CDWDOC(CDN(jd, jw)) = ADJUSTR(CDWDOC(CDN(jd, jw)))
                 enddo
                 write(WDO(j, 1), '(F10.3,",",F9.3,",",8X,100(F9.3,","))')jday,&
                     & QWDO(j), (qoutlet(i), i = 1, numoutlets)
                 write(WDO(j, 2), '(F10.3,",",F8.2,",",8X,100(F8.2,","))')jday,&
                     & TWDO(j), (toutlet(i), i = 1, numoutlets)
                 if(constituents)write(WDO(j, 3), '(F10.3,",",1000(A10,","))') &
                                     & jday, (CWDOC(CN(jc)), jc = 1, nac)
                 if(derived_calc)write(WDO(j, 4), '(F10.3,",",1000(A10,","))') &
                                     & jday,                                   &
                                     & (CDWDOC(CDN(jd, jw)), jd = 1, NACD(jw))
             enddo
         endif
     endif
 
  !**** DSI W2 Linkage File (W2L) (Supercedes Old Velocity vectors)
     if(VECTOR(1))then
    ! *** Apply the same linkage settings for all waterbodies
         if(jday>=NXTMVP(1) .OR. jday>=VPLD(VPLDP(1) + 1, 1))then
             if(jday>=VPLD(VPLDP(1) + 1, 1))then
                 VPLDP(1) = VPLDP(1) + 1
                 NXTMVP(1) = VPLD(VPLDP(1), 1)
             endif
             NXTMVP(1) = NXTMVP(1) + VPLF(VPLDP(1), 1)
 
      ! *** Write the W2L Snapshot (DSI)
             write(VPL(1))REAL(jday, 4)
 
      ! *** Compute the elevation, with ELWS zeroed for segments < CUS
             do jw = 1, nwb
                 do jb = BS(jw), BE(jw)
                     do i = US(jb) - 1, DS(jb) + 1
                         if(i<CUS(jb))then
                             wsel(i) = -9999
                         else
                             wsel(i) = ELWS(i)
                         endif
                     enddo
                 enddo
             enddo
 
             write(VPL(1))((wsel(i)), i = 1, imx)
             write(VPL(1))((REAL(U(k,i), 4), k = 1, kmx), i = 1, imx)
             write(VPL(1))((REAL(W(k,i), 4), k = 1, kmx), i = 1, imx)
 
             do i = 1, imx
        ! *** ABOVE ACTIVE LAYER
                 do k = 1, KTI(i) - 1
                     wdsi(k, i) = -9999.
                 enddo
        ! *** ACTIVE LAYERS
                 do k = KTI(i), KB(i)
                     wdsi(k, i) = T2(k, i)
                 enddo
        ! *** BELOW ACTIVE LAYER
                 do k = KB(i) + 1, kmx
                     wdsi(k, i) = -9999.
                 enddo
             enddo
             write(VPL(1))((wdsi(k, i), k = 1, kmx), i = 1, imx)
 
      ! *** Constituent data
             do jc = 1, nac
                 ic = CN(jc)
                 do i = 1, imx
          ! *** ABOVE ACTIVE LAYER
                     do k = 1, KTI(i) - 1
                         wdsi(k, i) = -9999.
                     enddo
          ! *** ACTIVE LAYERS
                     do k = KTI(i), KB(i)
                         wdsi(k, i) = C2(k, i, ic)
                     enddo
          ! *** BELOW ACTIVE LAYER
                     do k = KB(i) + 1, kmx
                         wdsi(k, i) = -9999.
                     enddo
                 enddo
                 write(VPL(1))((wdsi(k, i), k = 1, kmx), i = 1, imx)
             enddo
 
      !          WRITE (VPL(JW),*)  'New date ',JDAY,MONTH//GDCH//',',YEAR,KTWB(JW),(US(JB),JB=BS(JW),BE(JW))
      !          WRITE (VPL(JW),*) ((Z(I)*COSA(BS(JW))),     I=US(BS(JW)),DS(BE(JW)))
      !          WRITE (VPL(JW),*) ((EL(K,I),K=KTWB(JW),KMX),I=US(BS(JW)),DS(BE(JW)))
      !          WRITE (VPL(JW),*) ((U(K,I), K=KTWB(JW),KMX),I=US(BS(JW)),DS(BE(JW)))
      !          WRITE (VPL(JW),*) ((W(K,I), K=KTWB(JW),KMX),I=US(BS(JW)),DS(BE(JW)))
         endif
     endif
 
  !** Restart
 
     if(restart_out)then
         if(jday>=nxtmrs .OR. jday>=RSOD(rsodp + 1))then
             if(jday>=RSOD(rsodp + 1))then
                 rsodp = rsodp + 1
                 nxtmrs = RSOD(rsodp)
             endif
             nxtmrs = nxtmrs + RSOF(rsodp)
             if(RSOF(rsodp)>=1.0)then
                 write(ext, '(I0)')INT(jday)
             else
                 write(ext, '(I0,"_",I2)')INT(jday),                           &
                     & INT(100.*(jday - INT(jday)))                    ! Allows for writing out file names with JDAY FREQ less than 1 day
             endif
             ext = ADJUSTL(ext)
             l = LEN_TRIM(ext)
             rsofn = 'rso' // ext(1:l) // '.opt'
             call RESTART_OUTPUT(rsofn)
         endif
     endif
9045  format(200(e13.6, 1x))
9046  format(5x, '  Layer              [KWD] =', 12I7/(33x, 12I7))
9047  format(:5x, '  Outflow (m^3/sec)  [QSW] =', 12F7.2/(33x, 12F7.2))
9048  format(3x, 'Branch ', i0)
 
 
     end subroutine OUTPUTA
