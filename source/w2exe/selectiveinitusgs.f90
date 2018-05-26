!*==selectiveinitusgs.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
     subroutine SELECTIVEINITUSGS
 
 
     use SELECTIVE1USGS
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
     integer :: ifile, jj, n, nj
     character(8) :: tsshare
     character(20) :: fmtstr
!
!*** End of declarations rewritten by SPAG
!
 
 
     ifile = 1949
     tavg = 0.0
     tavgw = 0.0
     do jb = 1, nbr
         if(NSTR(jb)>0)then
             ifile = ifile + 1
             write(segnum, '(i0)')jb
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
                 write(fmtstr, '(I20)') nstr(jb)
                 write(ifile,                                                  &
     &'("      JDAY",'//trim(fmtstr)//'(6x,"T(C)"),'//trim(fmtstr)//'(3x,"Q(m3/s)"),'  &
     &//trim(fmtstr)//'(4x,"ELEVCL"))')
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
             write(fmtstr, '(I20)') nwd
             write(ifile,                                                      &
      &'("      JDAY",'//trim(fmtstr)//'(6x,"T(C)"),'//trim(fmtstr)//'(3x,"Q(m3/s)"),'           &
      &//trim(fmtstr)//'(4x,"ELEVCL"))')
         endif
     endif
 
     open(nunit, file = 'w2_selective.npt', status = 'old')
     read(nunit, '(///8x,f8.0)')tfrqtmp
     read(nunit, '(//8x,a8,i8,f8.0)')tempc, numtempc, tcdfreq
     nxtstr = tmstrt
     nxttcd = tmstrt
     nxtsplit = tmstrt
 
     allocate(tcnelev(numtempc), tcjb(numtempc), tcjs(numtempc),               &
            & tcelev(numtempc, 11), tctemp(numtempc), tctend(numtempc),        &
            & tctsrt(numtempc))
     allocate(ncountc(nst, nbr), tciseg(numtempc), tcklay(numtempc),           &
            & tcelevcon(numtempc))
     allocate(tcyearly(numtempc), tcntr(numtempc))
     allocate(volm(nwb), ncountcw(nwd), qwdfrac(nwd), qstrfrac(nst, nbr),      &
            & dynsel(numtempc), seld(numtempc), nxsel(numtempc),               &
            & temp2(numtempc))
 
     ncountc = 0
     read(nunit, '(/)')
     do j = 1, numtempc
         read(nunit, '(8x,a8,i8,i8,a8,f8.0,f8.0,f8.0,i8,10(f8.0))')tcntr(j),   &
            & tcjb(j), tcjs(j), tcyearly(j), tctsrt(j), tctend(j), tctemp(j),  &
            & tcnelev(j), (tcelev(j, n), n = 1, tcnelev(j))
         if(tcntr(j)=='      ST')then
             tcelev(j, tcnelev(j) + 1) = ESTR(tcjs(j), tcjb(j))
                                                     ! always put the original elevation as the last elevation
         else
             tcelev(j, tcnelev(j) + 1) = EWD(tcjs(j))
                                            ! always put the original elevation as the last elevation
         endif
     enddo
     read(nunit, '(/)')
     do j = 1, numtempc
         read(nunit, '(8x,i8,f8.0,A8)')tciseg(j), tcklay(j), dynsel(j)
     enddo
     read(nunit, '(/)')
     do j = 1, numtempc
         read(nunit, '(8x,a8)')tcelevcon(j)
     enddo
     read(nunit, '(//8x,a8,i8,2f8.0)')tspltc, numtsplt, tspltfreq, tsconv
 
     allocate(tsyearly(numtsplt), tstsrt(numtsplt), tstend(numtsplt),          &
            & tspltjb(numtsplt), tspltt(numtsplt), nouts(numtsplt))
     allocate(jstsplt(numtsplt, 10), kstrsplt(10), tspltcntr(numtsplt),        &
            & elcontspl(numtsplt))
     allocate(tsdepth(numtsplt, 10), tstype(numtsplt, 10),                     &
            & tsminfrac(numtsplt, 10), tsprior(numtsplt, 10))
     allocate(tsminhead(numtsplt, 10), tsmaxhead(numtsplt, 10),                &
            & tsmaxflow(numtsplt, 10), no_flow(numtsplt, 10),                  &
            & share_flow(numtsplt))
     allocate(tsdynsel(numtsplt), tsseld(numtsplt), nxtssel(numtsplt),         &
            & tstemp2(numtsplt))
     allocate(nout0(10), nout1(10), nout2(10), minfrac1(10), maxfrac1(10),     &
            & minfrac2(10), maxfrac2(10), splt2t(10), splt2e(10))
     allocate(ewdsav(nwd), wd_active(nwd), estrsav(nst, nbr),                  &
            & str_active(nst, nbr))
 
     read(nunit, '(/)')
     do j = 1, numtsplt
         read(nunit, '(8x,a8,i8,a8,3f8.0,2a8,i8,a8)')tspltcntr(j), tspltjb(j), &
            & tsyearly(j), tstsrt(j), tstend(j), tspltt(j), tsdynsel(j),       &
            & elcontspl(j), nouts(j), tsshare
         if(tspltc=='      ON')then
             if(nouts(j)<2)then
                 write(w2err, '(A,I0)')                                        &
                 &'ERROR-- Less than two outlets specified for blending group '&
                & , j
                 error_open = .TRUE.
                                ! will trigger the program to end when this subroutine is completed
                 return
             elseif(nouts(j)>10)then
                 write(w2err, '(A,I0)')                                        &
                 &'ERROR-- More than ten outlets specified for blending group '&
                & , j
                 error_open = .TRUE.
                 return
             endif
         endif
         share_flow(j) = tsshare=='      ON'
     enddo
     read(nunit, '(/)')
     do j = 1, numtsplt
         read(nunit, '(8x,10i8)')(jstsplt(j, n), n = 1, nouts(j))
     enddo
     read(nunit, '(/)')
     do j = 1, numtsplt
         read(nunit, '(8x,10f8.0)')(tsdepth(j, n), n = 1, nouts(j))
     enddo
     read(nunit, '(/)')
     do j = 1, numtsplt
         read(nunit, '(8x,10f8.0)')(tsminfrac(j, n), n = 1, nouts(j))
     enddo
     read(nunit, '(/)')
     do j = 1, numtsplt
         read(nunit, '(8x,10i8)')(tsprior(j, n), n = 1, nouts(j))
     enddo
     read(nunit, '(/)')
     do j = 1, numtsplt
         read(nunit, '(8x,10f8.0)')(tsminhead(j, n), n = 1, nouts(j))
     enddo
     read(nunit, '(/)')
     do j = 1, numtsplt
         read(nunit, '(8x,10f8.0)')(tsmaxhead(j, n), n = 1, nouts(j))
     enddo
     read(nunit, '(/)')
     do j = 1, numtsplt
         read(nunit, '(8x,10f8.0)')(tsmaxflow(j, n), n = 1, nouts(j))
     enddo
 
     estrsav = ESTR ! Save the original structure elevations
     ewdsav = EWD   ! Save the original withdrawal elevations
     do j = 1, numtsplt
         do n = 1, nouts(j)
             tstype(j, n) = "FIXED"
             if(tsdepth(j, n)>0.0)tstype(j, n) = "FLOAT"
             if(tsminfrac(j, n)>1.0)tsminfrac(j, n) = 1.0
                                                        ! remove unrealistic input value
             if(tsminhead(j, n)<0.0)tsminhead(j, n) = 0.0
                                                        ! remove unrealistic input value
             if(tsmaxhead(j, n)<0.0)tsmaxhead(j, n) = 0.0
                                                        ! remove unrealistic input value
             if(tsmaxflow(j, n)<0.0)tsmaxflow(j, n) = 0.0
                                                        ! remove unrealistic input value
         enddo
     enddo
     if(tsconv<=0.0)tsconv = 0.005    ! constrain the convergence criterion to be > 0.0 and <= 0.1
     if(tsconv>0.1)tsconv = 0.1
 
     read(nunit, '(//8x,i8)')tempn
     allocate(tempcrit(nwb, tempn), volmc(nwb, tempn))
     read(nunit, '(/)')
     do j = 1, tempn
         read(nunit, '(8x,10f8.0)')(tempcrit(jw, j), jw = 1, nwb)
                                                            ! Note max of 10 waterbodies
     enddo
     close(nunit)
 
     do jw = 1, nwb
         ifile = ifile + 1
         write(segnum, '(i0)')jw
         segnum = ADJUSTL(segnum)
         l = LEN_TRIM(segnum)
         if(restart_in)then
             open(ifile, file = 'Volume_wb' // segnum(1:l) // '.opt',          &
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
             open(ifile, file = 'Volume_wb' // segnum(1:l) // '.opt',          &
                 &status = 'unknown')
             write(ifile, 9001)
 
9001         format("jday    Volume    ", ("Volcrit      "))
         endif
     enddo
 
     if(tempc=='      ON')then
         do j = 1, numtempc
             if(tcyearly(j)=='     OFF')then
                 daytest = jday
             else
                 daytest = REAL(jdayg) + jday - INT(jday)
             endif
             if(daytest>=tctsrt(j) .AND. daytest<tctend(j))then
 
      ! initializing structure elevation
                 if(tcntr(j)=='      ST')then
                     jb = tcjb(j)                              ! set branch index
                     js = tcjs(j)                              ! set structure index
                     do nj = 1, tcnelev(j)                     ! making sure that structure is below water surface
                         if(tcelev(j, nj)<ELWS(DS(jb)))then
                             ncountc(js, jb) = nj
            ! ESTR(js,jb)=tcelev(j,ncountc(js,jb))    ! don't alter the elevation at this point.  Set it later.
                             exit
                         endif
                     enddo
 
      ! initializing withdrawal elevation
                 elseif(tcntr(j)=='      WD')then
                     jwd = tcjs(j)                             ! set withdrawal index
                     do nj = 1, tcnelev(j)                     ! making sure that structure is below water surface
                         if(tcelev(j, nj)<ELWS(IWD(jwd)))then
                             ncountcw(jwd) = nj
            ! EWD(jwd)=tcelev(j,ncountcw(jwd))        ! don't alter the elevation at this point.  Set it later.
                             exit
                         endif
                     enddo
                 endif
             endif
 
    ! Open dynamic selective withdrawal files
             if(dynsel(j)=='      ON')then
                 write(segnum, '(I0)')j
                 segnum = ADJUSTL(segnum)
                 l = LEN_TRIM(segnum)
                 seld(j) = 1009 + j
                 open(seld(j), file = 'dynselective' // segnum(1:l) // '.npt', &
                    & status = 'OLD')
                 read(seld(j), '(///1000F8.0)')nxsel(j), temp2(j)
                 tctemp(j) = temp2(j)
                 read(seld(j), '(1000F8.0)')nxsel(j), temp2(j)
             endif
         enddo
     endif
 
!    Open dynamic temperature target files for blending outlets
     if(tspltc=='      ON')then
         do j = 1, numtsplt
             if(tsdynsel(j)=='      ON')then
                 write(segnum, '(i0)')j
                 segnum = ADJUSTL(segnum)
                 l = LEN_TRIM(segnum)
                 tsseld(j) = 1009 + numtempc + j
                 open(tsseld(j), file = 'dynsplit_selective' // segnum(1:l)    &
                     & // '.npt', status = 'old')
                 read(tsseld(j), '(///2F8.0)')nxtssel(j), tstemp2(j)
                 tspltt(j) = tstemp2(j)
                 read(tsseld(j), '(2F8.0)')nxtssel(j), tstemp2(j)
             endif
         enddo
     endif
 
!    Test to see if the user specified inconsistent inputs. If so, stop with
!    an error message.
     if(tspltc=='      ON')then
         do j = 1, numtsplt
             do n = 1, nouts(j) - 1
                 do nj = n + 1, nouts(j)
                     if(jstsplt(j, n)==jstsplt(j, nj))then
                         write(w2err, '(A,I0)')                                &
      &'w2_selective.npt USGS ERROR-- Duplicate split outlet numbers in group '&
     & , j
                         error_open = .TRUE.
                                    ! will trigger the program to end when this subroutine is completed
                     endif
                 enddo
             enddo
         enddo
         do j = 1, numtsplt - 1
             do jj = j + 1, numtsplt
                 if((tstsrt(jj)>=tstsrt(j) .AND. tstsrt(jj)<tstend(j)) .OR.    &
                  & (tstend(jj)>tstsrt(j) .AND. tstend(jj)<=tstend(j)))then
                     if(tspltcntr(j)==tspltcntr(jj) .AND.                      &
                       &(tspltcntr(j)=='      WD' .OR. tspltjb(j)==tspltjb(jj))&
                      & )then
                         do n = 1, nouts(j)
                             do nj = 1, nouts(jj)
                                 if(jstsplt(j, n)==jstsplt(jj, nj))then
                                     write(w2err, '(A,I0,A)')                  &
                          &'w2_selective.npt USGS ERROR-- Split outlet number '&
                         & , jstsplt(j, n),                                    &
                          &' used in more than one group at a time.'
                                     error_open = .TRUE.
                                            ! will trigger the program to end when this subroutine is completed
                                 endif
                             enddo
                         enddo
                     endif
                 endif
             enddo
         enddo
         do j = 1, numtsplt
             do n = 1, nouts(j)
                 if(tsprior(j, n)< - 1)then
                     write(w2err, '(A,I0,A,I0,A)')                             &
                    &'w2_selective.npt USGS ERROR-- Priority input for outlet '&
                   & , jstsplt(j, n), ' in group ', j, ' is less than -1.'
                     error_open = .TRUE. ! will trigger the program to end when this subroutine is completed
                 endif
                 if(tsminhead(j, n)>0.0 .AND. tsmaxhead(j, n)>0.0 .AND.        &
                  & tsminhead(j, n)>tsmaxhead(j, n))then
                     write(wrn, '(A,I0,A,I0,A)')                               &
     &'w2_selective.npt USGS WARNING-- Minimum and maximum head constraints for&
     & outlet ', jstsplt(j, n), ' in group ', j,                               &
     &' are such that the outlet cannot ever be used.'
                     warning_open = .TRUE.
                 endif
                 if(tsdepth(j, n)>0.0 .AND. tsminhead(j, n)>0.0 .AND.          &
                  & tsdepth(j, n)<tsminhead(j, n))then
                     write(wrn, '(A,I0,A,I0,A)')                               &
                   &'w2_selective.npt USGS WARNING-- Depth of floating outlet '&
                  & , jstsplt(j, n), ' in group ', j,                          &
     &' is shallower than the minimum head constraint.  To honor the head const&
     &raint, no flow is possible for that outlet.'
                     warning_open = .TRUE.
                 endif
             enddo
         enddo
 
         if(tempc=='      ON')then
             do j = 1, numtsplt
                 do jj = 1, numtempc
                     if((tctsrt(jj)>=tstsrt(j) .AND. tctsrt(jj)<tstend(j)) .OR.&
                      & (tctend(jj)>tstsrt(j) .AND. tctend(jj)<=tstend(j)))then
                         if(tspltcntr(j)==tcntr(jj) .AND.                      &
                           &(tspltcntr(j)=='      WD' .OR. tspltjb(j)==tcjb(jj)&
                          & ))then
                             do n = 1, nouts(j)
                                 if(jstsplt(j, n)==tcjs(jj))then
                                     write(w2err, '(A,I0,A)')                  &
                                &'w2_selective.npt USGS ERROR-- Outlet number '&
                               & , tcjs(jj),                                   &
                             &' used in tower and blending group at same time.'
                                     error_open = .TRUE.
                                            ! will trigger the program to end when this subroutine is completed
                                 endif
                             enddo
                         endif
                     endif
                 enddo
             enddo
         endif
     endif
 
     if(tempc=='      ON')then
         do j = 1, numtempc - 1
             do jj = j + 1, numtempc
                 if((tctsrt(jj)>=tctsrt(j) .AND. tctsrt(jj)<tctend(j)) .OR.    &
                  & (tctend(jj)>tctsrt(j) .AND. tctend(jj)<=tctend(j)))then
                     if(tcntr(j)==tcntr(jj) .AND.                              &
                      & (tcntr(j)=='      WD' .OR. tcjb(j)==tcjb(jj)))then
                         if(tcjs(j)==tcjs(jj))then
                             write(w2err, '(A,I0,A)')                          &
                          &'w2_selective.npt USGS ERROR-- Tower outlet number '&
                         & , tcjs(j),                                          &
                          &' used more than once for overlapping dates.'
                             error_open = .TRUE.
                                        ! will trigger the program to end when this subroutine is completed
                         endif
                     endif
                 endif
             enddo
         enddo
     endif
 
 
     end subroutine SELECTIVEINITUSGS
