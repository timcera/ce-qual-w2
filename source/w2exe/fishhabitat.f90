!*==fishhabitat.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine FISHHABITAT(Iopenfish)
     use GLOBAL
     use MAIN
     use SCREENC
     use KINETIC, ONLY:O2, CHLA, NO3, NH4, PO4, TP, GAMMA, SED, SATO
     use TVDC, ONLY:CONSTITUENTS
     use NAMESC, ONLY:CNAME2, CDNAME2
     use LOGICC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     integer :: Iopenfish
     intent (in) Iopenfish
!
! Local variables
!
     real, allocatable, dimension(:), save :: cchla, cchlas, cdo, cdos, cgamma,&
       & cnh4, cnh4s, cno3, cno3s, cpo4, cpo4s, ctotp, ctotps, fishdo,         &
       & fishtemph, fishtempl, habvol, phabvol, ssedd, voltotbr, voltotwb
     character(80), save :: conavg, conhab, consod, consurf
     real, save :: dosat, o2corr, voltot
     character(80), allocatable, dimension(:), save :: fishname
     real, allocatable, dimension(:, :), save :: habvolbr, habvolwb, phabvolbr,&
       & phabvolwb
     integer, save :: ifish, jbfile, jbfile1, jjw, jwfile, jwfile1, kkmax,     &
                    & kseg, n, nseg
     integer, allocatable, dimension(:), save :: isegvol
!
!*** End of declarations rewritten by SPAG
!
 
     jwfile1 = 9549
     jbfile1 = 9749
 
     if(Iopenfish/=3)then
!        read input file
 
         if(nit==1 .OR. Iopenfish==0)then
             open(fishhabfn, file = 'w2_habitat.npt', status = 'old')
!            skip 1st 2 lines
             read(fishhabfn, *)
             read(fishhabfn, *)
             read(fishhabfn, *)ifish, conhab
             allocate(fishname(ifish), fishtempl(ifish), fishtemph(ifish),     &
                    & fishdo(ifish), habvol(ifish), phabvol(ifish),            &
                    & habvolbr(nbr, ifish), habvolwb(nwb, ifish),              &
                    & phabvolbr(nbr, ifish), phabvolwb(nwb, ifish),            &
                    & voltotbr(nbr), voltotwb(nwb))
             read(fishhabfn, *)
             do i = 1, ifish
                 read(fishhabfn, *)fishname(i), fishtempl(i), fishtemph(i),    &
                                 & fishdo(i)
             enddo
             read(fishhabfn, *)
             read(fishhabfn, *)nseg, conavg
                                 ! volume weighted averages of critical WQ parameters
             read(fishhabfn, *)
             allocate(isegvol(nseg), cdo(nseg), cpo4(nseg), cno3(nseg),        &
                    & cnh4(nseg), cchla(nseg), ctotp(nseg), cdos(nseg),        &
                    & cpo4s(nseg), cno3s(nseg), cnh4s(nseg), cchlas(nseg),     &
                    & ctotps(nseg), cgamma(nseg))
             allocate(ssedd(imx))
             read(fishhabfn, *)(isegvol(i), i = 1, nseg)
             read(fishhabfn, *)
             read(fishhabfn, *)kseg, consurf
                                  ! # of layers for surface average
             read(fishhabfn, *)
             read(fishhabfn, *)consod
             close(fishhabfn)
 
 
             if(restart_in)then
                 open(fishhabfn, file = conhab, position = 'APPEND')
                 jday1 = 0.0
                 rewind(fishhabfn)
                 read(fishhabfn, '(///)')
                 do i = 1, ifish
                     read(fishhabfn, *)
                 enddo
                 read(fishhabfn, '(//)', end = 10)
 
                 do while (jday1<jday)
                     read(fishhabfn, '(F10.0)', end = 10)jday1
                 enddo
                 backspace(fishhabfn)
10               jday1 = 0.0
                 jbfile = jbfile1
                 jwfile = jwfile1
                 do jw = 1, nwb
                     jwfile = jwfile + 1
                     write(segnum, '(I0)')jw
                     segnum = ADJUSTL(segnum)
                     l = LEN_TRIM(segnum)
                     open(jwfile, file = 'fish_habitat_wb' // segnum(1:l)      &
                         & // '.opt', position = 'APPEND')
                     jday1 = 0.0
                     rewind(jwfile)
                     read(jwfile, '(///)')
                     do i = 1, ifish
                         read(jwfile, *)
                     enddo
                     read(jwfile, '(//)', end = 15)
                     do while (jday1<jday)
                         read(jwfile, '(F10.0)', end = 15)jday1
                     enddo
                     backspace(jwfile)
15                   jday1 = 0.0
 
                     do jb = BS(jw), BE(jw)
                         jbfile = jbfile + 1
                         write(segnum, '(I0)')jb
                         segnum = ADJUSTL(segnum)
                         l = LEN_TRIM(segnum)
                         open(jbfile, file = 'fish_habitat_br' // segnum(1:l)  &
                            & // '.opt', position = 'APPEND')
                         jday1 = 0.0
                         rewind(jbfile)
                         read(jbfile, '(///)')
                         do i = 1, ifish
                             read(jbfile, *)
                         enddo
                         read(jbfile, '(//)', end = 16)
                         do while (jday1<jday)
                             read(jbfile, '(F10.0)', end = 16)jday1
                         enddo
                         backspace(jbfile)
16                       jday1 = 0.0
 
                     enddo
                 enddo
 
                 if(oxygen_demand)then
                     open(fishhabfn + 1, file = conavg, position = 'APPEND')
                     rewind(fishhabfn + 1)
                     read(fishhabfn + 1, '(//)')
                     do while (jday1<jday)
                         read(fishhabfn + 1, '(F10.0)', end = 20)jday1
                     enddo
                     backspace(fishhabfn + 1)
20                   jday1 = 0.0
                     open(fishhabfn + 2, file = consurf, position = 'APPEND')
                     rewind(fishhabfn + 2)
                     read(fishhabfn + 2, '(//)')
                     do while (jday1<jday)
                         read(fishhabfn + 2, '(F10.0)', end = 25)jday1
                     enddo
                     backspace(fishhabfn + 2)
25                   jday1 = 0.0
 
                     do jjw = 1, nwb
                         if(SEDIMENT_CALC(jjw))then
                             open(fishhabfn + 3, file = consod,                &
                                 &position = 'APPEND')
                             rewind(fishhabfn + 3)
                             read(fishhabfn + 3, '(/)')
                             do while (jday1<jday)
                                 read(fishhabfn + 3, '(F10.0)', end = 26)jday1
                             enddo
                             backspace(fishhabfn + 3)
26                           jday1 = 0.0
                             exit
                         endif
                     enddo
                 endif
             else
 
                 open(fishhabfn, file = conhab, status = 'unknown')
                 write(fishhabfn, *)                                           &
                             &'Fish habitat analysis: CE-QUAL-W2 model results'
                 write(fishhabfn, *)
                 write(fishhabfn, *)                                           &
     &'Species, Temperature minimum, Temperature maximum, Dissolved oxygen mini&
     &mum'
                 do i = 1, ifish
                     write(fishhabfn, "(a,',',t25,f8.2,',',f8.2,',',f8.2)")    &
                         & TRIM(fishname(i)), fishtempl(i), fishtemph(i),      &
                         & fishdo(i)
                 enddo
                 write(fishhabfn, *)
                 write(fishhabfn, 9003)(TRIM(fishname(i)), TRIM(fishname(i)),  &
                                     & i = 1, ifish)
 
                 jbfile = jbfile1
 
                 jwfile = jwfile1
                 do jw = 1, nwb
                     jwfile = jwfile + 1
                     write(segnum, '(I0)')jw
                     segnum = ADJUSTL(segnum)
                     l = LEN_TRIM(segnum)
                     open(jwfile, file = 'fish_habitat_wb' // segnum(1:l)      &
                         & // '.opt', status = 'UNKNOWN')
 
                     write(jwfile, *)                                          &
                             &'Fish habitat analysis: CE-QUAL-W2 model results'
                     write(jwfile, *)'FOR WATERBODY:', jw
                     write(jwfile, *)                                          &
     &'Species, Temperature minimum, Temperature maximum, Dissolved oxygen mini&
     &mum'
                     do i = 1, ifish
                         write(jwfile, "(a,',',t25,f8.2,',',f8.2,',',f8.2)")   &
                             & TRIM(fishname(i)), fishtempl(i), fishtemph(i),  &
                             & fishdo(i)
                     enddo
                     write(jwfile, *)
                     write(jwfile, 9003)(TRIM(fishname(i)), TRIM(fishname(i)), &
                                      & i = 1, ifish)
 
                     do jb = BS(jw), BE(jw)
                         jbfile = jbfile + 1
                         write(segnum, '(I0)')jb
                         segnum = ADJUSTL(segnum)
                         l = LEN_TRIM(segnum)
                         open(jbfile, file = 'fish_habitat_br' // segnum(1:l)  &
                            & // '.opt', status = 'UNKNOWN')
                         write(jbfile, *)                                      &
                             &'Fish habitat analysis: CE-QUAL-W2 model results'
                         write(jbfile, *)'FOR BRANCH:', jb
                         write(jbfile, *)                                      &
     &'Species, Temperature minimum, Temperature maximum, Dissolved oxygen mini&
     &mum'
                         do i = 1, ifish
                             write(jbfile,                                     &
                                 & "(a,',',t25,f8.2,',',f8.2,',',f8.2)")       &
                                 & TRIM(fishname(i)), fishtempl(i),            &
                                 & fishtemph(i), fishdo(i)
                         enddo
                         write(jbfile, *)
                         write(jbfile, 9003)                                   &
                             & (TRIM(fishname(i)), TRIM(fishname(i)), i = 1,   &
                             & ifish)
                     enddo
                 enddo
 
                 if(oxygen_demand)then
                     open(fishhabfn + 1, file = conavg, status = 'unknown')
                     write(fishhabfn + 1, '(a,80(1x,i4))')                     &
                          &'Volume weighted WQ parameters at segments:',       &
                         & (isegvol(i), i = 1, nseg)
                     write(fishhabfn + 1, 9001)                                &
                         & (TRIM(CNAME2(npo4)), isegvol(i), TRIM(CNAME2(nnh4)),&
                         & isegvol(i), TRIM(CNAME2(nno3)), isegvol(i),         &
                         & TRIM(CNAME2(ndo)), isegvol(i), TRIM(CDNAME2(12)),   &
                         & isegvol(i), TRIM(CDNAME2(14)), isegvol(i), i = 1,   &
                         & nseg)                                                                                                                                                                                        ! Chlor a and TP
9001                 format('JDAY,', <nseg>(6((a,'-',i3,','))))
 
                     open(fishhabfn + 2, file = consurf, status = 'unknown')
                     write(fishhabfn + 2, '(a,i4,a,80(1x,i4))')                &
                         & 'Surface (upper', kseg,                             &
                   &' model layers) Volume weighted WQ parameters at segments:'&
                  & , (isegvol(i), i = 1, nseg)
                     write(fishhabfn + 2, 9002)                                &
                         & (TRIM(CNAME2(npo4)), isegvol(i), TRIM(CNAME2(nnh4)),&
                         & isegvol(i), TRIM(CNAME2(nno3)), isegvol(i),         &
                         & TRIM(CNAME2(ndo)), isegvol(i), TRIM(CDNAME2(12)),   &
                         & isegvol(i), TRIM(CDNAME2(14)), isegvol(i),          &
                          &'Gamma(m-1)', isegvol(i), i = 1, nseg)
9002                 format('JDAY,', <nseg>(7((a,'-',i3,','))))
 
                     do jjw = 1, nwb
                         if(SEDIMENT_CALC(jjw))then
                             open(fishhabfn + 3, file = consod,                &
                                 &status = 'unknown')
                             write(fishhabfn + 3, "('JDAY,',1000(i3,','))")    &
                                 & (((i, i = US(jb), DS(jb)), jb = BS(jw),     &
                                 & BE(jw)), jw = 1, nwb)
                             exit
                         endif
                     enddo
                 endif
 
             endif
         endif
 
!        compute total volume and habitat volume
         habvol = 0.0
         voltot = 0.0
         hab = 100.0
         voltotbr = 0.0
         voltotwb = 0.0
         habvolbr = 0.0
         habvolwb = 0.0
         do jw = 1, nwb
             do jb = BS(jw), BE(jw)
                 do i = CUS(jb), DS(jb)
                     do k = KTWB(jw), KB(i)
                         voltot = voltot + VOL(k, i)
                         voltotbr(jb) = voltotbr(jb) + VOL(k, i)
                         voltotwb(jw) = voltotwb(jw) + VOL(k, i)
                         do ii = ifish, 1, -1
                             if(oxygen_demand)then
                                 if(T2(k, i)<=fishtemph(ii) .AND. T2(k, i)     &
                                  & >fishtempl(ii) .AND. O2(k, i)>=fishdo(ii)) &
                                  & then
                                     habvol(ii) = habvol(ii) + VOL(k, i)
                                     habvolbr(jb, ii) = habvolbr(jb, ii)       &
                                       & + VOL(k, i)
                                     habvolwb(jw, ii) = habvolwb(jw, ii)       &
                                       & + VOL(k, i)
                                     hab(k, i) = ii
                                 endif
                             elseif(T2(k, i)<=fishtemph(ii) .AND. T2(k, i)     &
                                  & >fishtempl(ii))then
                                 habvol(ii) = habvol(ii) + VOL(k, i)
                                 habvolbr(jb, ii) = habvolbr(jb, ii)           &
                                   & + VOL(k, i)
                                 habvolwb(jw, ii) = habvolwb(jw, ii)           &
                                   & + VOL(k, i)
                                 hab(k, i) = ii
                             endif
                         enddo
                     enddo
                 enddo
             enddo
         enddo
 
         do ii = 1, ifish
             phabvol(ii) = habvol(ii)/voltot
             do jw = 1, nwb
                 phabvolwb(jw, ii) = habvolwb(jw, ii)/voltotwb(jw)
                 do jb = BS(jw), BE(jw)
                     phabvolbr(jb, ii) = habvolbr(jb, ii)/voltotbr(jb)
                 enddo
             enddo
         enddo
 
!        write out results
 
         write(fishhabfn, 9004)jday, (100.*phabvol(i), habvol(i), i = 1, ifish)
         jbfile = jbfile1
         jwfile = jwfile1
         do jw = 1, nwb
             jwfile = jwfile + 1
             write(jwfile, 9004)jday, (100.*phabvolwb(jw, i), habvolwb(jw, i), &
                              & i = 1, ifish)
             do jb = BS(jw), BE(jw)
                 jbfile = jbfile + 1
                 write(jbfile, 9004)jday,                                      &
                                  & (100.*phabvolbr(jb, i), habvolbr(jb, i),   &
                                  & i = 1, ifish)
             enddo
         enddo
 
         if(oxygen_demand)then
 
             cno3 = 0.0
             cdo = 0.0
             cpo4 = 0.0
             cchla = 0.0
             cnh4 = 0.0
             cgamma = 0.0
             ctotp = 0.0
 
             do n = 1, nseg
 
                 i = isegvol(n)
 
    ! Find waterbody associated with this segment
                 do jjw = 1, nwb
                     if(i>=US(BS(jjw)) .AND. i<=DS(BE(jjw)))exit
                 enddo
 
                 voltot = 0.0
                 kkmax = MIN(kseg, KB(i) - KTWB(jjw))
                                         ! kseg is the # of layers
                 if(kkmax<0)cycle
                 do k = KTWB(jjw), KB(i)
                     voltot = voltot + VOL(k, i)
                     cpo4(n) = cpo4(n) + PO4(k, i)*VOL(k, i)
                     if(k<=KTWB(jjw) + kkmax)cgamma(n) = cgamma(n)             &
                      & + GAMMA(k, i)*VOL(k, i)
        ! NOTE*********** No credit for superstauration - if DO > saturation, then set DO=100% saturation
                     dosat = SATO(T2(k, i), 0.D0, PALT(i), SALT_WATER(jjw))
                     if(O2(k, i)>dosat)then
                         o2corr = dosat
                     else
                         o2corr = O2(k, i)
                     endif
 
                     cdo(n) = cdo(n) + o2corr*VOL(k, i)
                     cno3(n) = cno3(n) + NO3(k, i)*VOL(k, i)
                     cchla(n) = cchla(n) + CHLA(k, i)*VOL(k, i)
                     ctotp(n) = ctotp(n) + TP(k, i)*VOL(k, i)
                     cnh4(n) = cnh4(n) + NH4(k, i)*VOL(k, i)
 
                     if(k==KTWB(jjw) + kkmax)then
                         cdos(n) = cdo(n)/voltot
                         cpo4s(n) = cpo4(n)/voltot
                         cno3s(n) = cno3(n)/voltot
                         cnh4s(n) = cnh4(n)/voltot
                         cchlas(n) = cchla(n)/voltot
                         ctotps(n) = ctotp(n)/voltot
                         cgamma(n) = cgamma(n)/voltot
                     endif
 
 
                 enddo
                 cpo4(n) = cpo4(n)/voltot
                 cdo(n) = cdo(n)/voltot
                 cno3(n) = cno3(n)/voltot
                 cnh4(n) = cnh4(n)/voltot
                 cchla(n) = cchla(n)/voltot
                 ctotp(n) = ctotp(n)/voltot
             enddo
 
             write(fishhabfn + 2, 9005)jday,                                   &
                                     & (cpo4s(n), cnh4s(n), cno3s(n), cdos(n), &
                                     & ctotps(n), cchlas(n), cgamma(n), n = 1, &
                                     & nseg)
             write(fishhabfn + 1, 9005)jday,                                   &
                                     & (cpo4(n), cnh4(n), cno3(n), cdo(n),     &
                                     & ctotp(n), cchla(n), n = 1, nseg)
 
 
 
!            write out sed for each segment
             ssedd = 0.0
             do jw = 1, nwb
                 if(SEDIMENT_CALC(jw))then
                     do jb = BS(jw), BE(jw)
                         do i = US(jb), DS(jb)
                             if(KTWB(jw)<=KB(i))then
                                 do k = KTWB(jw), KB(i)
                                     ssedd(i) = ssedd(i) + SED(k, i)*VOL(k, i)
                                 enddo
                             else
                                 ssedd(i) = -99.
                             endif
                         enddo
                     enddo
                 endif
             enddo
 
             do jjw = 1, nwb
                 if(SEDIMENT_CALC(jjw))then
                     write(fishhabfn + 3, "(f10.3,',',1000(e13.4,','))")jday,  &
                         & (((ssedd(i), i = US(jb), DS(jb)), jb = BS(jw),      &
                         & BE(jw)), jw = 1, nwb)                                                                 ! OUTPUT IS IN GRAMS
                     exit
                 endif
             enddo
 
         endif
 
 
     else
 
         deallocate(isegvol, cdo, cpo4, cno3, cnh4, cchla, ctotp, cdos, cpo4s, &
                  & cno3s, cnh4s, cchlas, ctotps, cgamma, ssedd, fishname,     &
                  & fishtempl, fishtemph, fishdo, habvol, phabvol, habvolbr,   &
                  & habvolwb, phabvolbr, phabvolwb, voltotbr, voltotwb)
         close(fishhabfn)
         close(fishhabfn + 1)
         close(fishhabfn + 2)
         close(fishhabfn + 3)
         jbfile = jbfile1
         jwfile = jwfile1
         do jw = 1, nwb
             jwfile = jwfile + 1
             close(jwfile)
             do jb = BS(jw), BE(jw)
                 jbfile = jbfile + 1
                 close(jbfile)
             enddo
         enddo
 
     endif
9003  format('JDAY,', <ifish>('%VOL-', a, ',', 'HAB-VOL(m3)-', a, ','))
9004  format(f10.3, ',', <ifish>(f8.2, ',', e12.4, ','))
9005  format(f10.3, ',', <nseg>(7(f10.4, ',')))
 
 
     end subroutine FISHHABITAT
