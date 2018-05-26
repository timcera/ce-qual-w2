!*==selective.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine SELECTIVE
 
     use SELECTIVE1
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
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: daytest, elr, qall, tcomp, tempbot, tempest, temptop, tmod, wsel
     integer :: ifile, jj, jjw, kk, ks, kstr
!
!*** End of declarations rewritten by SPAG
!
  !** Timestep violation entry point  210 CONTINUE
 
     if(tspltc=='      ON')then
         do j = 1, numtsplt
             if(TSYEARLY(j)=='     OFF')then
                 daytest = jday
             else
                 daytest = REAL(jdayg) + jday - INT(jday)
             endif
             if(nxtsplit>TSTSRT(j) .AND. daytest<=TSTSRT(j))                   &
              & nxtsplit = TSTSRT(j)
         enddo
     endif
 
     if(tspltc=='      ON' .AND. jday>=nxtsplit)then
 
         do j = 1, numtsplt
             if(TSYEARLY(j)=='     OFF')then
                 daytest = jday
             else
                 daytest = REAL(jdayg) + jday - INT(jday)
             endif
             if(daytest>=TSTSRT(j) .AND. daytest<TSTEND(j))then
    ! DO STRUCTURES FIRST
                 do jw = 1, nwb
                     do jb = BS(jw), BE(jw)
                         if(TSPLTJB(j)==jb .AND. TSPLTCNTR(j)=='      ST')then
                             qall = 0.0
                             do jj = 1, NOUTS(j)
                                 qall = qall + QSTR(JSTSPLT(j, jj), TSPLTJB(j))
                                                           ! SUM UP ALL THE FLOWS
                                 elr = SINA(jb)*DLX(DS(jb))*0.5
                                 do k = KTWB(jw), KB(DS(jb))
                                     if(EL(k, DS(jb))                          &
                                      & - elr<ESTR(JSTSPLT(j, jj), TSPLTJB(j)))&
                                      & exit                                                                                                                  !SW 10/17/01
                                 enddo
                                 kstr = k - 1
                                 KSTRSPLT(jj) = MIN(kstr, KB(DS(jb)))
                             enddo
                             do jj = 1, NOUTS(j)
                                             ! cb 11/11/12 dividing total flow between outlets for temperature test - if no flow there is no temperature test
                                 QSTR(JSTSPLT(j, jj), TSPLTJB(j))              &
                                   & = qall/REAL(NOUTS(j))
                             enddo
                             id = DS(jb)
                             elr = SINA(jb)*DLX(id)*0.5
                                                   ! CB 10/14/11
                             wsel = ELWS(id) - elr ! CB 10/14/11
                             call DOWNSTREAM_WITHDRAWAL_ESTIMATE(JSTSPLT(j, 1),&
                               & temptop, ESTR(JSTSPLT(j, 1), TSPLTJB(j)))
                             call DOWNSTREAM_WITHDRAWAL_ESTIMATE(JSTSPLT(j, 2),&
                               & tempbot, ESTR(JSTSPLT(j, 2), TSPLTJB(j)))
                             if(ESTR(JSTSPLT(j, 1), TSPLTJB(j))>wsel .AND.     &
                              & ELCONTSPL(j)=='     OFF')then                                ! NO FLOWS THROUG THIS OUTLET IF WSEL BELOW LEVEL OF OUTLET  ! CB 10/14/11
                                 QSTR(JSTSPLT(j, 1), TSPLTJB(j)) = 0.0
                                 QSTRFRAC(JSTSPLT(j, 1), TSPLTJB(j)) = 0.0
 
                             elseif(temptop>TSPLTT(j) .AND. tempbot>TSPLTT(j)) &
                                  & then                                     ! NO FLOWS THROUG THIS OUTLET IF T1 AND T2 > TCRITERIA
                                 QSTR(JSTSPLT(j, 1), TSPLTJB(j)) = 0.0
                                 QSTRFRAC(JSTSPLT(j, 1), TSPLTJB(j)) = 0.0
 
              !ELSEIF(T2(KSTRSPLT(1),DS(JB)) < TSPLTT(J)) THEN   ! ALL FLOWS FROM TOP IF TCRITERIA < TOUTLET
                             elseif(temptop<TSPLTT(j))then
                                                 ! ALL FLOWS FROM TOP IF TCRITERIA < TOUTLET
                                 QSTR(JSTSPLT(j, 1), TSPLTJB(j)) = qall
                                 QSTRFRAC(JSTSPLT(j, 1), TSPLTJB(j)) = 1.0
 
                !QSTR(JSTSPLT(J,1),TSPLTJB(J))=QALL*(TSPLTT(J)-T2(KSTRSPLT(2),DS(JB)))/(T2(KSTRSPLT(1),DS(JB))-T2(KSTRSPLT(2),DS(JB)))
                             elseif(ABS(temptop - tempbot)<0.0001)then
                                 QSTR(JSTSPLT(j, 1), TSPLTJB(j)) = qall
                                 QSTRFRAC(JSTSPLT(j, 1), TSPLTJB(j)) = 1.0
                             else
                                 QSTR(JSTSPLT(j, 1), TSPLTJB(j))               &
                                   & = qall*(TSPLTT(j) - tempbot)              &
                                   & /(temptop - tempbot)
                                 QSTRFRAC(JSTSPLT(j, 1), TSPLTJB(j))           &
                                   & = QSTR(JSTSPLT(j, 1), TSPLTJB(j))/qall
                             endif
 
                             QSTR(JSTSPLT(j, 2), TSPLTJB(j))                   &
                               & = qall - QSTR(JSTSPLT(j, 1), TSPLTJB(j))
                             QSTRFRAC(JSTSPLT(j, 2), TSPLTJB(j))               &
                               & = QSTR(JSTSPLT(j, 2), TSPLTJB(j))/qall
                             exit
                         endif
                     enddo
                 enddo
    ! DO WITHDRAWALS NEXT
                 do jwd = 1, nwd
                     if(TSPLTCNTR(j)=='      WD')then
                         qall = 0.0
                         do jjb = 1, nbr
                             if(IWD(jwd)>=US(jjb) .AND. IWD(jwd)<=DS(jjb))exit
                         enddo
                         do jjw = 1, nwb
                             if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
                         enddo
                         do jj = 1, NOUTS(j)
                             qall = qall + QWD(JSTSPLT(j, jj))
                                               ! SUM UP ALL THE FLOWS
                             elr = SINA(jjb)*DLX(IWD(jwd))*0.5
                             do k = KTWB(jjw), KB(IWD(jwd))
                                 if(EL(k, IWD(jwd)) - elr<EWD(JSTSPLT(j, jj))) &
                                  & exit                                                                                                            !SW 10/17/01
                             enddo
                             kstr = k - 1
                             KSTRSPLT(jj) = MIN(kstr, KB(IWD(jwd)))
                         enddo
                         jj = 1   ! ASSIGN FLOW TO FIRST OUTLET
                         wsel = ELWS(IWD(jwd)) - elr      ! CB 10/14/11
                         i = IWD(jwd)
                         call LATERAL_WITHDRAWAL_ESTIMATE(JSTSPLT(j, 1),       &
                           & temptop, EWD(JSTSPLT(j, 1)))
                         call LATERAL_WITHDRAWAL_ESTIMATE(JSTSPLT(j, 2),       &
                           & tempbot, EWD(JSTSPLT(j, 2)))
                         if(EWD(JSTSPLT(j, 1))>wsel .AND. TCELEVCON(j)         &
                           &=='     OFF')then
                             QWD(JSTSPLT(j, 1)) = 0.0
                             QWDFRAC(JSTSPLT(j, 1)) = 0.0
                         elseif(temptop>TSPLTT(j) .AND. tempbot>TSPLTT(j))then
                                                                              ! NO FLOWS THROUG THIS OUTLET IF T1 AND T2 > TCRITERIA
                             QWD(JSTSPLT(j, 1)) = 0.0
                             QWDFRAC(JSTSPLT(j, 1)) = 0.0
 
                         elseif(temptop<TSPLTT(j))then
                                                 ! ALL FLOWS FROM TOP IF TCRITERIA < TOUTLET
                             QWD(JSTSPLT(j, 1)) = qall
                             QWDFRAC(JSTSPLT(j, 1)) = 1.0
 
                !QWD(JSTSPLT(J,1))=QALL*(TSPLTT(J)-T2(KSTRSPLT(2),IWD(JWD)))/(T2(KSTRSPLT(1),IWD(JWD))-T2(KSTRSPLT(2),IWD(JWD)))
                         elseif(ABS(temptop - tempbot)<0.0001)then
                             QWD(JSTSPLT(j, 1)) = qall
                             QWDFRAC(JSTSPLT(j, 1)) = 1.0
                         else
                             QWD(JSTSPLT(j, 1)) = qall*(TSPLTT(j) - tempbot)   &
                               & /(temptop - tempbot)
                             QWDFRAC(JSTSPLT(j, 1)) = QWD(JSTSPLT(j, 1))/qall
                         endif
 
                         QWD(JSTSPLT(j, 2)) = qall - QWD(JSTSPLT(j, 1))
                         QWDFRAC(JSTSPLT(j, 2)) = QWD(JSTSPLT(j, 2))/qall
                         exit
                     endif
                 enddo
             endif
         enddo
 
         nxtsplit = nxtsplit + tcdfreq
     endif
     if(tspltc=='      ON')then
 
         do j = 1, numtsplt
             if(TSYEARLY(j)=='     OFF')then
                 daytest = jday
             else
                 daytest = REAL(jdayg) + jday - INT(jday)
             endif
             if(daytest>=TSTSRT(j) .AND. daytest<TSTEND(j))then
    ! DO STRUCTURES FIRST
                 do jw = 1, nwb
                     do jb = BS(jw), BE(jw)
                         if(TSPLTJB(j)==jb .AND. TSPLTCNTR(j)=='      ST')then
                             qall = 0.0
                             do jj = 1, NOUTS(j)
                                 qall = qall + QSTR(JSTSPLT(j, jj), TSPLTJB(j))
                                                           ! SUM UP ALL THE FLOWS
                                 elr = SINA(jb)*DLX(DS(jb))*0.5
                                 do k = KTWB(jw), KB(DS(jb))
                                     if(EL(k, DS(jb))                          &
                                      & - elr<ESTR(JSTSPLT(j, jj), TSPLTJB(j)))&
                                      & exit                                                                                                                  !SW 10/17/01
                                 enddo
                                 kstr = k - 1
                                 KSTRSPLT(jj) = MIN(kstr, KB(DS(jb)))
                             enddo
                             QSTR(JSTSPLT(j, 1), TSPLTJB(j))                   &
                               & = QSTRFRAC(JSTSPLT(j, 1), TSPLTJB(j))*qall
                             QSTR(JSTSPLT(j, 2), TSPLTJB(j))                   &
                               & = QSTRFRAC(JSTSPLT(j, 2), TSPLTJB(j))*qall
                             exit
                         endif
                     enddo
                 enddo
    ! DO WITHDRAWALS NEXT
                 do jwd = 1, nwd
                     if(TSPLTCNTR(j)=='      WD')then
                         qall = 0.0
                         do jjb = 1, nbr
                             if(IWD(jwd)>=US(jjb) .AND. IWD(jwd)<=DS(jjb))exit
                         enddo
                         do jjw = 1, nwb
                             if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
                         enddo
                         do jj = 1, NOUTS(j)
                             qall = qall + QWD(JSTSPLT(j, jj))
                                               ! SUM UP ALL THE FLOWS
                             elr = SINA(jjb)*DLX(IWD(jwd))*0.5
                             do k = KTWB(jjw), KB(IWD(jwd))
                                 if(EL(k, IWD(jwd)) - elr<EWD(JSTSPLT(j, jj))) &
                                  & exit                                                                                                            !SW 10/17/01
                             enddo
                             kstr = k - 1
                             KSTRSPLT(jj) = MIN(kstr, KB(IWD(jwd)))
                         enddo
                         QWD(JSTSPLT(j, 1)) = QWDFRAC(JSTSPLT(j, 1))*qall
                         QWD(JSTSPLT(j, 2)) = QWDFRAC(JSTSPLT(j, 2))*qall
                         exit
                     endif
                 enddo
             endif
         enddo
 
     endif
 
 
     if(jday>=nxtstr)then
         nxtstr = nxtstr + tfrqtmp
         ifile = 1949
         do jb = 1, nbr
             if(NSTR(jb)>0)then
                 ifile = ifile + 1
                 write(ifile,                                                  &
                     &'(F10.4,<NSTR(JB)>F10.2,<NSTR(JB)>F10.2,<NSTR(JB)>F10.2)'&
                    & )jday, (TAVG(i, jb), i = 1, NSTR(jb)),                   &
                     & (QSTR(i, jb), i = 1, NSTR(jb)),                         &
                     & (ESTR(i, jb), i = 1, NSTR(jb))
             endif
         enddo
         if(nwd>0)then
             ifile = ifile + 1
             write(ifile, '(F10.4,<NWD>F10.2,<NWD>F10.2,<NWD>F10.2)')jday,     &
                 & (TAVGW(i), i = 1, nwd), (QWD(i), i = 1, nwd),               &
                 & (EWD(i), i = 1, nwd)
         endif
         ! TEMPERATURE CONTROL LOGIC
 
         ! COMPUTING RESERVOIR VOLUME AND VOLUME BELOW 'TEMPCRIT'
         volmc = 0.0
         volm = 0.0
         do jw = 1, nwb
             kt = KTWB(jw)
             do jb = BS(jw), BE(jw)
                 do i = CUS(jb), DS(jb)
                     volm(jw) = volm(jw) + BH2(kt, i)*DLX(i)
                     do k = kt + 1, KB(i)
                         volm(jw) = volm(jw) + BH(k, i)*DLX(i)
                     enddo
                     do kk = 1, tempn
                         if(T2(kt, i)<=TEMPCRIT(jw, kk))volmc(jw, kk)          &
                          & = volmc(jw, kk) + BH2(kt, i)*DLX(i)
                         do k = kt + 1, KB(i)
                             if(T2(k, i)<=TEMPCRIT(jw, kk))volmc(jw, kk)       &
                              & = volmc(jw, kk) + BH(k, i)*DLX(i)
                         enddo
                     enddo
                 enddo
             enddo
 
             ifile = ifile + 1
             write(ifile, 9001)jday, volm(jw), (volmc(jw, kk), kk = 1, tempn)
9001         format(f8.2, 100(g12.4, g12.4))
         enddo
 
     endif
 
 
 
     if(tempc=='      ON' .AND. jday>=nxttcd)then
 
!        IF DYNAMIC SELECTIVE CHANGE TEMPERATURE
 
         do j = 1, numtempc
             if(DYNSEL(j)=='      ON')then
                 do while (jday>=NXSEL(j))
                     TCTEMP(j) = TEMP2(j)
                     read(SELD(j), '(1000F8.0)')NXSEL(j), TEMP2(j)
                 enddo
             endif
         enddo
 
 
!        STRUCTURES
 
         do jw = 1, nwb
             do jb = BS(jw), BE(jw)
                 do js = 1, nst
                     do j = 1, numtempc
 
 
                         if(TCJB(j)==jb .AND. TCJS(j)==js .AND. TCNTR(j)       &
                           &=='      ST')then
                             if(TCISEG(j)==0)then
                                 tcomp = TAVG(TCJS(j), TCJB(j))
                                          !CB 9/8/06   TAVG(JSMON(J),JBMON(J))
                             elseif(TCISEG(j)<0)then
                                 tcomp = TWDO(ABS(TCISEG(j)))
                                            ! SW 11/26/10
                             else
 
!                                CHECKING TO SEE IF THE MONITORING SEGMENT
!                                TCISEG IS IN THE SAME BRANCH AND WATER BODY
!                                AS THE STRUCTURE
                                 do jjb = 1, nbr
                                     if(TCISEG(j)>=US(jjb) .AND. TCISEG(j)     &
                                      & <=DS(jjb))exit
                                 enddo
                                 do jjw = 1, nwb
                                     if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
                                 enddo
 
                                 if(TCKLAY(j)<0)then
                                     k = INT(ABS(TCKLAY(j)))
                                 else
                                     do k = KTWB(jjw), KB(TCISEG(j))
                                         if(DEPTHB(k, TCISEG(j))>TCKLAY(j))exit
                                     enddo
                                     k = MIN(k, KB(TCISEG(j)))
                                 endif
                                 tcomp = T2(k, TCISEG(j))
                             endif
                             if(TCYEARLY(j)=='     OFF')then
                                 daytest = jday
                             else
                                 daytest = REAL(jdayg) + jday - INT(jday)
                             endif
                             if(daytest>=TCTSRT(j) .AND. daytest<TCTEND(j))then
                                 if(tcomp>TCTEMP(j) .AND. TCNELEV(j)           &
                                  & >NCOUNTC(js, jb))then
               ! MAKING SURE THAT THE NEXT LOWER STRUCTURE FOR A PARTICULAR 'J' IS FOUND
                                     do nn = NCOUNTC(js, jb) + 1, TCNELEV(j)
                                         if(TCELEV(j, nn)<ESTR(js, jb))then
                                         NCOUNTC(js, jb) = nn
                                         ESTR(js, jb)                          &
                                           & = TCELEV(j, NCOUNTC(js, jb))
                                         exit
                                         endif
                                     enddo
                                 elseif(tcomp<TCTEMP(j) .AND. NCOUNTC(js, jb)  &
                                      & >1)then
                 ! TO PREVENT THIS HAPPENING AT EACH TIME IT CHECKS IT AND HENCE OSCIALLTING BACK AND FORTH - CHECK THE TEMP AT THE UPPER OUTLET ALSO
                                     if(TCISEG(j)>0)then
                                         if(jb==jjb)then
                                         do ks = KTWB(jw), KB(DS(jb))
                                         if(DEPTHB(ks, TCISEG(j))              &
                                           & >TCELEV(j, NCOUNTC(js, jb) - 1))  &
                                           & exit                                                                                                            !TC 01/03/02
                                         enddo
                                         ks = MIN(ks, KB(TCISEG(j)))
                                         tmod = T2(ks, DS(jb))
                                         else
                                         tmod = T2(k, TCISEG(j))
                                         endif
                                         if(tmod<TCTEMP(j) .AND.               &
                                           & TCELEV(j, NCOUNTC(js, jb) - 1)    &
                                           & <ELWS(DS(jb)))then
                      ! MAKING SURE THAT THE NEXT UPPER STRUCTURE FOR A PARTICULAR 'J' IS FOUND
                                         do nn = NCOUNTC(js, jb) - 1, 1, -1
                                         if(TCELEV(j, nn)>ESTR(js, jb))then
                                         NCOUNTC(js, jb) = nn
                                         ESTR(js, jb)                          &
                                           & = TCELEV(j, NCOUNTC(js, jb))
                                         exit
                                         endif
                                         enddo
                                         endif
                                     endif
                         ! CB 9/8/06
                                     if(TCISEG(j)==0)then
!                                        CALCULATE THE ESTIMATED OUTFLOW
!                                        TEMPERATURE AT HIGHER PORTS WHEN
!                                        TCOMP<TCTEMP(J), AND MOVE UP IF
!                                        HIGHER PORT STILL MEETS TO CRITERIA -
!                                        THIS DOESN'T HAPPEN WHEN TCISEG < 0
                                         do nn = 1, NCOUNTC(js, jb) - 1
                                         id = DS(jb)
                                         kt = KTWB(jw)
                                         call DOWNSTREAM_WITHDRAWAL_ESTIMATE   &
                                           & (js, tempest, TCELEV(j, nn))
                                         if(tempest<TCTEMP(j) .AND.            &
                                           & TCELEV(j, nn)<ELWS(DS(jb)))then
                                         NCOUNTC(js, jb) = nn
                                         ESTR(js, jb)                          &
                                           & = TCELEV(j, NCOUNTC(js, jb))
                                         exit
                                         endif
                                         enddo
                                     endif
                                 endif
                                 if(TCELEVCON(j)=='      ON' .AND. TCNELEV(j)  &
                                  & >NCOUNTC(js, jb) .AND. ESTR(js, jb)        &
                                  & >ELWS(DS(jb)))then
                                     NCOUNTC(js, jb) = NCOUNTC(js, jb) + 1
                                     ESTR(js, jb) = TCELEV(j, NCOUNTC(js, jb))
                                 endif
                             endif
                         endif
                     enddo
                 enddo
             enddo
         enddo
 
 ! Withdrawals
 
 
         do jwd = 1, nwd
             do j = 1, numtempc
                 if(TCJS(j)==jwd .AND. TCNTR(j)=='      WD')then
                     if(TCISEG(j)==0)then
!                        TCOMP=TOUT(JB)
                         tcomp = TAVGW(TCJS(j))
                                   !CB 9/8/06   TAVGW(JSMON(J))
                     elseif(TCISEG(j)<0)then
                         tcomp = TWDO(ABS(TCISEG(j)))
                     else
 
!                        CHECKING TO SEE IF THE MONITORING SEGMENT TCISEG IS
!                        IN THE SAME BRANCH AND WATER BODY AS THE WITHDRAWAL
                         do jjb = 1, nbr
                             if(TCISEG(j)>=US(jjb) .AND. TCISEG(j)<=DS(jjb))   &
                              & exit
                         enddo
                         do jjw = 1, nwb
                             if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
                         enddo
 
                         if(TCKLAY(j)<0)then
                             k = INT(ABS(TCKLAY(j)))
                         else
                             do k = KTWB(jjw), KB(TCISEG(j))
                                 if(DEPTHB(k, TCISEG(j))>TCKLAY(j))exit
                             enddo
                             k = MIN(k, KB(TCISEG(j)))
                         endif
                         tcomp = T2(k, TCISEG(j))
                     endif
                     if(TCYEARLY(j)=='     OFF')then
                         daytest = jday
                     else
                         daytest = REAL(jdayg) + jday - INT(jday)
                     endif
                     if(daytest>=TCTSRT(j) .AND. daytest<TCTEND(j))then
                         if(tcomp>TCTEMP(j) .AND. TCNELEV(j)>NCOUNTCW(jwd))then
               ! MAKING SURE THAT THE NEXT LOWER STRUCTURE FOR A PARTICULAR 'J' IS FOUND
                             do nn = NCOUNTCW(jwd) + 1, TCNELEV(j)
                                 if(TCELEV(j, nn)<EWD(jwd))then
                                     NCOUNTCW(jwd) = nn
                                     EWD(jwd) = TCELEV(j, NCOUNTCW(jwd))
                                     exit
                                 endif
                             enddo
                         elseif(tcomp<TCTEMP(j) .AND. NCOUNTCW(jwd)>1)then
                 ! TO PREVENT THIS HAPPENING AT EACH TIME IT CHECKS IT AND HENCE OSCIALLTING BACK AND FORTH - CHECK THE TEMP AT THE UPPER OUTLET ALSO
                             if(TCISEG(j)>0)then
                                 tmod = T2(k, TCISEG(j))
                                 if(tmod<TCTEMP(j) .AND.                       &
                                  & TCELEV(j, NCOUNTCW(jwd) - 1)<ELWS(IWD(jwd))&
                                  & )then
                      ! MAKING SURE THAT THE NEXT UPPER STRUCTURE FOR A PARTICULAR 'J' IS FOUND
                                     do nn = NCOUNTCW(jwd) - 1, 1, -1
                                         if(TCELEV(j, nn)>EWD(jwd))then
                                         NCOUNTCW(jwd) = nn
                                         EWD(jwd) = TCELEV(j, NCOUNTCW(jwd))
                                         exit
                                         endif
                                     enddo
                                 endif
                             endif
                         ! CB 9/8/06
                             if(TCISEG(j)==0)then
!                                CALCULATE THE ESTIMATED OUTFLOW TEMPERATURE
!                                AT HIGHER PORTS WHEN TCOMP<TCTEMP(J), AND
!                                MOVE UP IF HIGHER PORT STILL MEETS TO CRITERIA
                                 i = MAX(CUS(JBWD(jwd)), IWD(jwd))
                                 do nn = 1, NCOUNTCW(jwd) - 1
                                     call LATERAL_WITHDRAWAL_ESTIMATE(jwd,     &
                                       & tempest, TCELEV(j, nn))
                                     if(tempest<TCTEMP(j) .AND. TCELEV(j, nn)  &
                                      & <ELWS(IWD(jwd)))then
                                         NCOUNTCW(jwd) = nn
                                         EWD(jwd) = TCELEV(j, NCOUNTCW(jwd))
                                         exit
                                     endif
                                 enddo
                             endif
                         endif
                         if(TCELEVCON(j)=='      ON' .AND. TCNELEV(j)          &
                          & >NCOUNTCW(jwd) .AND. EWD(jwd)>ELWS(IWD(jwd)))then
                             NCOUNTCW(jwd) = NCOUNTCW(jwd) + 1
                             EWD(jwd) = TCELEV(j, NCOUNTCW(jwd))
                         endif
                     endif
                 endif
             enddo
         enddo
 
         nxttcd = nxttcd + tcdfreq
     endif
 
     return
     entry DEALLOCATE_SELECTIVE
     deallocate(TCNELEV, TCJB, TCJS, TCELEV, TCTEMP, TCTEND, TCTSRT, NCOUNTC,  &
              & TCISEG, TCKLAY, TCELEVCON, ELCONTSPL)
     deallocate(TSPLTJB, TSPLTT, NOUTS, JSTSPLT, KSTRSPLT, TCYEARLY, jbmon,    &
              & jsmon, TCNTR, TSPLTCNTR, jstspltt)
     deallocate(volm, monctr, NCOUNTCW, QWDFRAC, QSTRFRAC)
     deallocate(TEMPCRIT, volmc, DYNSEL, SELD, NXSEL, TEMP2, TSYEARLY, TSTEND, &
              & TSTSRT)
 
     end subroutine SELECTIVE
