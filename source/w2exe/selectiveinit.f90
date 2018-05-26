!*==selectiveinit.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine SELECTIVEINIT
 
 
 
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
     real :: daytest
     integer :: ifile, n
!
!*** End of declarations rewritten by SPAG
!
 
 
!**  Task 2: Calculations                                                     
!***********************************************************************************************************************************
 
!    **
     ifile = 1949
     tavg = 0.0
     tavgw = 0.0
     do jb = 1, nbr
         if(NSTR(jb)>0)then
             ifile = ifile + 1
             write(segnum, '(I0)')jb
             segnum = ADJUSTL(segnum)
             l = LEN_TRIM(segnum)
 
             if(restart_in)then
                 open(ifile, file = 'str_br' // segnum(1:l) // '.opt',         &
                    & position = 'APPEND')
                 jday1 = 0.0
                 rewind(ifile)
                 read(ifile, '(/)', end = 10)
                 do while (jday1<jday)
                     read(ifile, '(F10.0)', end = 10)jday1
                 enddo
                 backspace(ifile)
10               jday1 = 0.0
             else
                 open(ifile, file = 'str_br' // segnum(1:l) // '.opt',         &
                     &status = 'unknown')
                 write(ifile, *)'Branch:', jb, ' # of structures:', NSTR(jb),  &
                               &' outlet temperatures'
                 write(ifile,                                                  &
     &'("      JDAY",<nstr(jb)>(6x,"T(C)"),<nstr(jb)>(3x,"Q(m3/s)"),<nstr(jb)>(&
     &4x,"ELEVCL"))')
             endif
         endif
     enddo
 
     if(nwd>0)then
         ifile = ifile + 1
         if(restart_in)then
             open(ifile, file = 'wd_out.opt', position = 'APPEND')
             jday1 = 0.0
             rewind(ifile)
             read(ifile, '(/)', end = 20)
             do while (jday1<jday)
                 read(ifile, '(F10.0)', end = 20)jday1
             enddo
             backspace(ifile)
20           jday1 = 0.0
         else
             open(ifile, file = 'wd_out.opt', status = 'unknown')
             write(ifile, *)'Withdrawals: # of withdrawals:', nwd,             &
                           &' outlet temperatures'
             write(ifile,                                                      &
      &'("      JDAY",<nwd>(6x,"T(C)"),<nwd>(3x,"Q(m3/s)"),<nwd>(4x,"ELEVCL"))'&
     & )
         endif
     endif
 
 
     open(nunit, file = 'w2_selective.npt', status = 'old')
     read(nunit, '(///8X,F8.0)')tfrqtmp
     nxtstr = tmstrt
     read(nunit, '(//8X,A8,I8,F8.0)')tempc, numtempc, tcdfreq
     nxttcd = tmstrt
     nxtsplit = tmstrt
 
     allocate(tcnelev(numtempc), tcjb(numtempc), tcjs(numtempc),               &
            & tcelev(numtempc, 11), tctemp(numtempc), tctend(numtempc),        &
            & tctsrt(numtempc), ncountc(nst, nbr), tciseg(numtempc),           &
            & tcklay(numtempc), tcelevcon(numtempc))
     allocate(tcyearly(numtempc), jbmon(numtempc), jsmon(numtempc),            &
            & tcntr(numtempc))
     allocate(volm(nwb), monctr(numtempc), ncountcw(nwd), qwdfrac(nwd),        &
            & qstrfrac(nst, nbr), dynsel(numtempc), seld(numtempc),            &
            & nxsel(numtempc), temp2(numtempc))
 
     do j = 1, 2
         read(nunit, *)
     enddo
     ncountc = 0
     do j = 1, numtempc
         read(nunit, '(8X,A8,I8,I8,A8,F8.0,F8.0,F8.0,I8,10(F8.0))')tcntr(j),   &
            & tcjb(j), tcjs(j), tcyearly(j), tctsrt(j), tctend(j), tctemp(j),  &
            & tcnelev(j), (tcelev(j, n), n = 1, tcnelev(j))
         if(tcntr(j)=='      ST')then
             tcelev(j, tcnelev(j) + 1) = ESTR(tcjs(j), tcjb(j))
                                                       ! ALWAYS PUT THE ORIGINAL ELEVATION AS THE LAST ELEVATION
         else
             tcelev(j, tcnelev(j) + 1) = EWD(tcjs(j))
                                              ! ALWAYS PUT THE ORIGINAL ELEVATION AS THE LAST ELEVATION
         endif
     enddo
     do j = 1, 2
         read(nunit, *)
     enddo
     do j = 1, numtempc
         read(nunit, '(8X,I8,F8.0,A8)')tciseg(j), tcklay(j), dynsel(j)
     enddo
     do j = 1, 2
         read(nunit, *)
     enddo
     do j = 1, numtempc
         read(nunit, '(8X,A8)')tcelevcon(j)
     enddo
     do j = 1, 2
         read(nunit, *)
     enddo
     read(nunit, '(8X,A8,I8)')tspltc, numtsplt
 
     allocate(tsyearly(numtsplt), tstsrt(numtsplt), tstend(numtsplt),          &
            & tspltjb(numtsplt), tspltt(numtsplt), nouts(numtsplt),            &
            & jstsplt(numtsplt, 10), kstrsplt(numtsplt), tspltcntr(numtsplt))
     allocate(jstspltt(numtsplt, 10), elcontspl(numtsplt))
 
     do j = 1, 2
         read(nunit, *)
     enddo
     do j = 1, numtsplt
      !READ(NUNIT,'(8X,A8,I8,A8,F8.0,F8.0,F8.0,I8,10I8)')TSPLTCNTR(J),TSPLTJB(J),TSYEARLY(J),TSTSRT(J),TSTEND(J),TSPLTT(J),NOUTS(J),(JSTSPLTT(J,N),N=1,NOUTS(J))
         read(nunit, '(8X,A8,I8,A8,F8.0,F8.0,F8.0,I8,2I8,A8)')tspltcntr(j),    &
            & tspltjb(j), tsyearly(j), tstsrt(j), tstend(j), tspltt(j),        &
            & nouts(j), (jstspltt(j, n), n = 1, 2), elcontspl(j)
         nouts(j) = 2           ! NUMBER OF OUTLETS FOR EACH SPLIT FLOW PERIOD LIMITED TO 2
      !IF(NOUTS(J).GT.2)WRITE(*,*)'TCD NOUTS > 2 - ONLY FIRST 2 WILL BE USED'
     enddo
     jstsplt = jstspltt                                                                                            ! CB 10/14/11 START
     do j = 1, numtsplt
                       !REODERING OUTLETS SO THAT HIGHEST ELEVATION STRUCTURE ON TOP (ASSUMING 2 SPLIT OUTLETS)
!        IF(TCNTR(J) == '      ST')THEN
         if(tspltcntr(j)=='      ST')then                                                                         ! cb 11/11/12
             if(ESTR(jstspltt(j, 1), tspltjb(j))                               &
              & <ESTR(jstspltt(j, 2), tspltjb(j)))then
                 jstsplt(j, 1) = jstspltt(j, 2)
                 jstsplt(j, 2) = jstspltt(j, 1)
             endif
!            ELSE IF(TCNTR(J) == '      WD')THEN
         elseif(tspltcntr(j)=='      WD')then                                                                          ! cb 11/11/12
             if(EWD(jstspltt(j, 1))<EWD(jstspltt(j, 2)))then
                 jstsplt(j, 1) = jstspltt(j, 2)
                 jstsplt(j, 2) = jstspltt(j, 1)
             endif
         endif
     enddo                                                                                                         ! CB 10/14/11 END
     do j = 1, 2
         read(nunit, *)
     enddo
     read(nunit, '(8X,I8)')tempn
     do j = 1, 2
         read(nunit, *)
     enddo
     allocate(tempcrit(nwb, tempn), volmc(nwb, tempn))
     do j = 1, tempn
         read(nunit, '(8X,100F8.0)')(tempcrit(jw, j), jw = 1, nwb)
                                                              ! NOTE MAX OF 100 WATERBODIES   sw 4/20/15
     enddo
     close(nunit)
 
 
     do jw = 1, nwb
         ifile = ifile + 1
         write(segnum, '(I0)')jw
         segnum = ADJUSTL(segnum)
         l = LEN_TRIM(segnum)
         if(restart_in)then
             open(ifile, file = 'VOLUME_WB' // segnum(1:l) // '.OPT',          &
                & position = 'APPEND')
             jday1 = 0.0
             rewind(ifile)
             read(ifile, '(/)', end = 40)
             do while (jday1<jday)
                 read(ifile, '(F10.0)', end = 40)jday1
             enddo
             backspace(ifile)
40           jday1 = 0.0
         else
             open(ifile, file = 'VOLUME_WB' // segnum(1:l) // '.OPT',          &
                 &status = 'UNKNOWN')
             write(ifile, 9001)
 
9001         format("JDAY    VOLUME    ", <tempn>("VOLCRIT      "))
         endif
     enddo
 
 
!    INITIALIZING STRUCTURE ELEVATION IF STRUCTURE
     if(tempc=='      ON')then
         do jw = 1, nwb
             do jb = BS(jw), BE(jw)
                 do js = 1, nst
                     do j = 1, numtempc
                         if(tcjb(j)==jb .AND. tcjs(j)==js .AND. tcntr(j)       &
                           &=='      ST')then
                             if(tcyearly(j)=='     OFF')then
                                 daytest = jday
                             else
                                 daytest = REAL(jdayg) + jday - INT(jday)
                             endif
                             if(daytest>=tctsrt(j) .AND. daytest<tctend(j))then
               ! MAKING SURE THAT STRUCTURE IS BELOW WATER SURFACE
                                 do nn = 1, tcnelev(j)
                                     if(tcelev(j, nn)<ELWS(DS(jb)))then
                                         ncountc(js, jb) = nn
                                         ESTR(js, jb)                          &
                                           & = tcelev(j, ncountc(js, jb))
                                         exit
                                     endif
                                 enddo
                             endif
                         endif
                     enddo
                 enddo
             enddo
         enddo
 
 
   ! INITIALIZING STRUCTURE ELEVATION IF WITHDRAWAL
 
         do jwd = 1, nwd
 
             do j = 1, numtempc
                 if(tcjs(j)==jwd .AND. tcntr(j)=='      WD')then
                     if(tcyearly(j)=='     OFF')then
                         daytest = jday
                     else
                         daytest = REAL(jdayg) + jday - INT(jday)
                     endif
                     if(daytest>=tctsrt(j) .AND. daytest<tctend(j))then
               ! MAKING SURE THAT STRUCTURE IS BELOW WATER SURFACE
                         do nn = 1, tcnelev(j)
                             if(tcelev(j, nn)<ELWS(IWD(jwd)))then
                                 ncountcw(jwd) = nn
                                 EWD(jwd) = tcelev(j, ncountcw(jwd))
                                 exit
                             endif
                         enddo
                     endif
                 endif
             enddo
 
         enddo
     endif
 
  ! OPEN DYNAMIC SELECTIVE WITHDRAWAL FILES
 
     do j = 1, numtempc
         if(dynsel(j)=='      ON')then
             write(segnum, '(I0)')j
             segnum = ADJUSTL(segnum)
             l = LEN_TRIM(segnum)
             seld(j) = 1009 + j
             open(seld(j), file = 'dynselective' // segnum(1:l) // '.npt',     &
                & status = 'OLD')
             read(seld(j), '(///1000F8.0)')nxsel(j), temp2(j)
             tctemp(j) = temp2(j)
             read(seld(j), '(1000F8.0)')nxsel(j), temp2(j)
         endif
     enddo
 
 
     end subroutine SELECTIVEINIT
