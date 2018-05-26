!*==selectiveusgs.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************************************************************************
!**                                                   S E L E C T I V E                                                           **
!***********************************************************************************************************************************
 
     subroutine SELECTIVEUSGS
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
     real :: addfrac, blendfrac, daytest, elev1, elev2, elr, etemp, etemp1,    &
           & etemp2, excess_frac, lastfrac, lastfrac2, maxelev, maxtemp,       &
           & minelev, mintemp, qall, q_notblended, sumelev, sumfrac, sumtemp,  &
           & sum_maxfrac1, sum_maxfrac2, sum_minfrac0, tcomp, tempest, tmod,   &
           & ttarg, wsel
     integer :: ifile, j2hi, j2lo, j2max, j2min, j2pref, jj, jjw, jst, kk, ks, &
              & kstr, n, ng0, ng1max, ng1min, nj, num_left, num_noflow, prior1,&
              & prior2
     character(20) :: fmtstr
!
!*** End of declarations rewritten by SPAG
!
 
 
 
!    qstr = qstrsav   ! xxx not sure how to do this yet -- need to reset QSTR
!    when control periods expire qwd  = qwdsav    ! xxx not sure how to do
 
!    this yet -- need to reset QWD when control periods expire
     str_active = .FALSE.
     wd_active = .FALSE.
     j2pref = 1
 
!    Update some date variables and determine which outlets are being actively
!    blended or adjusted.
     if(tspltc=='      ON')then
         do j = 1, numtsplt
             if(TSYEARLY(j)=='     OFF')then
                 daytest = jday
             else
                 daytest = REAL(jdayg) + jday - INT(jday)
             endif
             if(nxtsplit>TSTSRT(j) .AND. daytest<=TSTSRT(j))                   &
              & nxtsplit = TSTSRT(j)
             if(daytest>=TSTSRT(j) .AND. daytest<TSTEND(j))then
                 do jj = 1, NOUTS(j)
                     if(TSPLTCNTR(j)=='      ST')then
                         str_active(JSTSPLT(j, jj), TSPLTJB(j)) = .TRUE.
                     elseif(TSPLTCNTR(j)=='      WD')then
                         wd_active(JSTSPLT(j, jj)) = .TRUE.
                     endif
                 enddo
             endif
         enddo
     endif
     if(tempc=='      ON')then
         do j = 1, numtempc
             if(TCYEARLY(j)=='     OFF')then
                 daytest = jday
             else
                 daytest = REAL(jdayg) + jday - INT(jday)
             endif
             if(daytest>=TCTSRT(j) .AND. daytest<TCTEND(j))then
                 if(TCNTR(j)=='      ST')then
                     str_active(TCJS(j), TCJB(j)) = .TRUE.
                 elseif(TCNTR(j)=='      WD')then
                     wd_active(TCJS(j)) = .TRUE.
                 endif
             endif
         enddo
     endif
 
!    Reset elevations of outlets back to original values outside of control
!    periods
     do jst = 1, nst
         do jb = 1, nbr
             if(.NOT.str_active(jst, jb))ESTR(jst, jb) = ESTRSAV(jst, jb)
         enddo
     enddo
     do jwd = 1, nwd
         if(.NOT.wd_active(jwd))EWD(jwd) = EWDSAV(jwd)
     enddo
 
 
!    Check to see if it's time to update temperature targets and flow
!    fractions for blended groups.
     if(tspltc=='      ON' .AND. jday>=nxtsplit)then
 
  ! Update the temperature targets
         do j = 1, numtsplt
             if(TSDYNSEL(j)=='      ON')then
                 do while (jday>=NXTSSEL(j))
                     TSPLTT(j) = TSTEMP2(j)
                     read(TSSELD(j), '(2F8.0)')NXTSSEL(j), TSTEMP2(j)
                 enddo
             endif
         enddo
 
         do j = 1, numtsplt
             qall = 0.0                                                      ! sum up all the flows
             sumfrac = 0.0                                                   ! sum of flow fraction multipliers
             do jj = 1, NOUTS(j)
                 if(TSPLTCNTR(j)=='      ST')then
                     qall = qall + QSTR(JSTSPLT(j, jj), TSPLTJB(j))
                     sumfrac = sumfrac + QSTRFRAC(JSTSPLT(j, jj), TSPLTJB(j))
                 elseif(TSPLTCNTR(j)=='      WD')then
                     qall = qall + QWD(JSTSPLT(j, jj))
                     sumfrac = sumfrac + QWDFRAC(JSTSPLT(j, jj))
                 endif
             enddo
             if(TSYEARLY(j)=='     OFF')then
                 daytest = jday
             else
                 daytest = REAL(jdayg) + jday - INT(jday)
             endif
 
    ! Do blending calculations if date is in correct window or the window was just entered.
    ! If flows are zero and flow fractions are already initialized, then leave everything alone.
             if(daytest>=TSTSRT(j) .AND. daytest<TSTEND(j) .AND.               &
              & (qall>0.0 .OR. sumfrac<0.001 .OR. daytest<TSTSRT(j)            &
              & + tspltfreq))then
 
      ! Do structures first
                 if(TSPLTCNTR(j)=='      ST')then
                     jb = TSPLTJB(j)                                         ! set branch index
                     do jw = 1, nwb                                          ! set waterbody index
                         if(jb>=BS(jw) .AND. jb<=BE(jw))exit
                     enddo
                     no_flow(j, :) = .FALSE.
                     num_noflow = 0
                     do jj = 1, NOUTS(j)
                         jst = JSTSPLT(j, jj)
                         elr = SINA(jb)*DLX(DS(jb))*0.5
                         wsel = ELWS(DS(jb)) - elr                           ! compute water-surface elevation
                         ESTR(jst, jb) = ESTRSAV(jst, jb)                    ! reset outlet elevation to original
                         if(TSTYPE(j, jj)=="FLOAT")then
                             ESTR(jst, jb) = wsel - TSDEPTH(j, jj)
                         elseif(ESTR(jst, jb)>wsel)then
                             if(ELCONTSPL(j)=='     OFF')then
                                 no_flow(j, jj) = .TRUE.                     ! no flow-- high and dry
                                 num_noflow = num_noflow + 1
                             else
                                 ESTR(jst, jb) = wsel                        ! poor man's floating outlet
                             endif
                         endif
                         if(.NOT.no_flow(j, jj) .AND. TSMINHEAD(j, jj)         &
                          & >0.0 .AND. wsel - ESTR(jst, jb)<TSMINHEAD(j, jj))  &
                          & then
                             no_flow(j, jj) = .TRUE.                         ! minimum head criterion not met -- no flow
                             num_noflow = num_noflow + 1
                         endif
                         if(.NOT.no_flow(j, jj) .AND. TSMAXHEAD(j, jj)         &
                          & >0.0 .AND. wsel - ESTR(jst, jb)>TSMAXHEAD(j, jj))  &
                          & then
                             no_flow(j, jj) = .TRUE.                         ! maximum head criterion exceeded -- no flow
                             num_noflow = num_noflow + 1
                         endif
                         do k = KTWB(jw), KB(DS(jb))
                             if(EL(k, DS(jb)) - elr<ESTR(jst, jb))exit
                         enddo
                         KSTRSPLT(jj) = MIN(k - 1, KB(DS(jb)))
                         QSTRFRAC(jst, jb) = 0.0                             ! initialize flow fractions
                     enddo
 
        ! Use priority inputs to determine which outlets to use
                     prior1 = -999
                     do jj = 1, NOUTS(j)
                         if(.NOT.no_flow(j, jj) .AND. TSPRIOR(j, jj)>=0)then
                             if(prior1== - 999 .OR. TSPRIOR(j, jj)<prior1)     &
                              & prior1 = TSPRIOR(j, jj)
                         endif
                     enddo
                     prior2 = -999
                     do jj = 1, NOUTS(j)
                         if(.NOT.no_flow(j, jj) .AND. TSPRIOR(j, jj)>=0 .AND.  &
                          & TSPRIOR(j, jj)>prior1)then
                             if(prior2== - 999 .OR. TSPRIOR(j, jj)<prior2)     &
                              & prior2 = TSPRIOR(j, jj)
                         endif
                     enddo
 
        ! Outlets with a priority of -1 get used, but are not blended
                     ng0 = 0
                     q_notblended = 0.0
                     do jj = 1, NOUTS(j)
                         jst = JSTSPLT(j, jj)
                         if(.NOT.no_flow(j, jj) .AND. TSPRIOR(j, jj)== - 1)then
                             ng0 = ng0 + 1
                             NOUT0(ng0) = jj
                             if(QSTR(jst, jb)>TSMAXFLOW(j, jj) .AND.           &
                              & TSMAXFLOW(j, jj)>0.0)then
                                 q_notblended = q_notblended + TSMAXFLOW(j, jj)
                                 QSTRFRAC(jst, jb) = TSMAXFLOW(j, jj)/qall
                             elseif(qall>0.0)then
                                 q_notblended = q_notblended + QSTR(jst, jb)
                                 QSTRFRAC(jst, jb) = QSTR(jst, jb)/qall
                             endif
                         endif
                     enddo
                     sum_minfrac0 = 0.0
                     if(qall>0.0)sum_minfrac0 = q_notblended/qall
 
        ! Outlets with priority 1 and 2 may be used and blended.
                     ng1 = 0
                     ng2 = 0
                     sum_minfrac1 = 0.0
                     sum_minfrac2 = 0.0
                     sum_maxfrac1 = 0.0
                     sum_maxfrac2 = 0.0
                     do jj = 1, NOUTS(j)
                         if(.NOT.no_flow(j, jj))then
                             if(TSPRIOR(j, jj)==prior1)then
                                 ng1 = ng1 + 1
                                 NOUT1(ng1) = jj
                                 MAXFRAC1(ng1) = 1.0
                                 if(qall>0.0 .AND. TSMAXFLOW(j, jj)>0.0)       &
                                  & MAXFRAC1(ng1)                              &
                                  & = MIN(1.0, TSMAXFLOW(j, jj)/qall)
                                 MINFRAC1(ng1) = TSMINFRAC(j, jj)
                                 if(TSMINFRAC(j, jj)<0.0)then
                                     MINFRAC1(ng1) = 0.0
                                     if(qall>0.0)MINFRAC1(ng1)                 &
                                      & = MIN(1.0, ABS(TSMINFRAC(j, jj))/qall)
                                 endif
                                 if(MINFRAC1(ng1)>MAXFRAC1(ng1))MINFRAC1(ng1)  &
                                  & = MAXFRAC1(ng1)
                                 sum_minfrac1 = sum_minfrac1 + MINFRAC1(ng1)
                                 sum_maxfrac1 = sum_maxfrac1 + MAXFRAC1(ng1)
 
                             elseif(TSPRIOR(j, jj)==prior2)then
                                 ng2 = ng2 + 1
                                 NOUT2(ng2) = jj
                                 MAXFRAC2(ng2) = 1.0
                                 if(qall>0.0 .AND. TSMAXFLOW(j, jj)>0.0)       &
                                  & MAXFRAC2(ng2)                              &
                                  & = MIN(1.0, TSMAXFLOW(j, jj)/qall)
                                 MINFRAC2(ng2) = TSMINFRAC(j, jj)
                                 if(TSMINFRAC(j, jj)<0.0)then
                                     MINFRAC2(ng2) = 0.0
                                     if(qall>0.0)MINFRAC2(ng2)                 &
                                      & = MIN(1.0, ABS(TSMINFRAC(j, jj))/qall)
                                 endif
                                 if(MINFRAC2(ng2)>MAXFRAC2(ng2))MINFRAC2(ng2)  &
                                  & = MAXFRAC2(ng2)
                                 sum_minfrac2 = sum_minfrac2 + MINFRAC2(ng2)
                                 sum_maxfrac2 = sum_maxfrac2 + MAXFRAC2(ng2)
                             endif
                         endif
                     enddo
 
        ! If minimum flows are overspecified, then the priority 2 minimum flow fractions are decreased.
                     if(ng2>0 .AND. sum_minfrac0 + sum_minfrac1 + sum_minfrac2>&
                      & 1.0)then
                         if(sum_minfrac0 + sum_minfrac1>=1.0)then
                             ng2 = 0
                             sum_minfrac2 = 0.0
                         else
                             do n = 1, ng2
                                 MINFRAC2(n) = MINFRAC2(n)                     &
                                   & *(1.0 - sum_minfrac0 - sum_minfrac1)      &
                                   & /sum_minfrac2
                             enddo
                             sum_minfrac2 = 1.0 - sum_minfrac0 - sum_minfrac1
                         endif
                     endif
 
        ! If minimum flows are still overspecified, then the priority 1 minimum flow fractions are decreased.
                     if(ng1>0 .AND. sum_minfrac0 + sum_minfrac1>1.0)then
                         if(sum_minfrac0>=1.0)then
                             ng1 = 0
                             sum_minfrac1 = 0.0
                         else
                             do n = 1, ng1
                                 MINFRAC1(n) = MINFRAC1(n)*(1.0 - sum_minfrac0)&
                                   & /sum_minfrac1
                             enddo
                             sum_minfrac1 = 1.0 - sum_minfrac0
                         endif
                     endif
 
        ! If group 1 has 3 or more outlets and group 2 has no outlets, then redistribute priorities based on elevation.
        ! Keep the highest and lowest elevation outlets in group 1, and put other active outlets into nonblended category
        ! with their minimum flows.  If ties in elevation exist, go with the first in the list.
                     if(ng1>2 .AND. ng2==0)then
                         ng1max = 1
                         ng1min = 1
                         jst = JSTSPLT(j, NOUT1(1))
                         maxelev = ESTR(jst, jb)
                         minelev = ESTR(jst, jb)
                         do n = 2, ng1
                             jst = JSTSPLT(j, NOUT1(n))
                             if(ESTR(jst, jb)>maxelev)then
                                 maxelev = ESTR(jst, jb)
                                 ng1max = n
                             elseif(ESTR(jst, jb)<minelev)then
                                 minelev = ESTR(jst, jb)
                                 ng1min = n
                             endif
                         enddo
                         blendfrac = 1.0 - sum_minfrac0 - sum_minfrac1 +       &
                                   & MINFRAC1(ng1max) + MINFRAC1(ng1min)
                         if(MAXFRAC1(ng1max) + MAXFRAC1(ng1min)<blendfrac)then
                             if(sum_maxfrac1<1.0 - sum_minfrac0)then
                                 write(wrn, '(A,I0,A,F0.3)')                   &
                     &'Warning-- Maximum flows for outlets exceeded for group '&
                    & , j, ' at day ', jday
                                 warning_open = .TRUE.
                                 do n = 1, ng1
                                     if(n/=ng1max .AND. n/=ng1min)MINFRAC1(n)  &
                                      & = MAXFRAC1(n)
                                 enddo
                             else
                                 excess_frac = blendfrac - MAXFRAC1(ng1max)    &
                                   & - MAXFRAC1(ng1min)
                                 num_left = ng1 - 2
                                 do nj = 1, ng1                    ! iterative process to redistribute excess flows
                                     if(num_left>0 .AND. excess_frac>0.0)then
                                         addfrac = excess_frac/num_left
                                         do n = 1, ng1
                                         if(n/=ng1max .AND. n/=ng1min .AND.    &
                                           & MAXFRAC1(n) - MINFRAC1(n)>0.00001)&
                                           & then
                                         if(MINFRAC1(n) + addfrac>MAXFRAC1(n)) &
                                           & then
                                         num_left = num_left - 1
                                         excess_frac = excess_frac -           &
                                           & (MAXFRAC1(n) - MINFRAC1(n))
                                         MINFRAC1(n) = MAXFRAC1(n)
                                         else
                                         excess_frac = excess_frac - addfrac
                                         MINFRAC1(n) = MINFRAC1(n) + addfrac
                                         endif
                                         endif
                                         enddo
                                     endif
                                 enddo
                             endif
                         endif
                         do n = 1, ng1                       ! assign the other priority 1 outlets to nonblended status
                             if(n/=ng1max .AND. n/=ng1min)then
                                 ng0 = ng0 + 1
                                 NOUT0(ng0) = NOUT1(n)
                                 jst = JSTSPLT(j, NOUT1(n))
                                 sum_minfrac0 = sum_minfrac0 + MINFRAC1(n)
                                 q_notblended = q_notblended + qall*MINFRAC1(n)
                                 QSTRFRAC(jst, jb) = MINFRAC1(n)
                             endif
                         enddo
                         ng1 = 1                   ! rearrange outlets-- one in each priority group, but same priority
                         ng2 = 1
                         NOUT2(1) = NOUT1(ng1min)
                         MINFRAC2(1) = MINFRAC1(ng1min)
                         MAXFRAC2(1) = MAXFRAC1(ng1min)
                         sum_minfrac2 = MINFRAC1(ng1min)
                         sum_maxfrac2 = MAXFRAC1(ng1min)
                         NOUT1(1) = NOUT1(ng1max)
                         MINFRAC1(1) = MINFRAC1(ng1max)
                         MAXFRAC1(1) = MAXFRAC1(ng1max)
                         sum_minfrac1 = MINFRAC1(ng1max)
                         sum_maxfrac1 = MAXFRAC1(ng1max)
                         prior2 = prior1
                     endif
 
        ! If only two blended outlets, ensure that they are in separate groups.
                     if(ng1==2 .AND. ng2==0)then
                         ng1 = 1
                         ng2 = 1
                         NOUT2(1) = NOUT1(2)
                         MINFRAC2(1) = MINFRAC1(2)
                         MAXFRAC2(1) = MAXFRAC1(2)
                         sum_minfrac2 = MINFRAC1(2)
                         sum_maxfrac2 = MAXFRAC1(2)
                         sum_minfrac1 = MINFRAC1(1)
                         sum_maxfrac1 = MAXFRAC1(1)
                         prior2 = prior1
                     endif
 
 
        ! Begin the blending decisions.
        ! No usable outlets.  All flow fractions remain at zero.
                     if(NOUTS(j)==num_noflow)then
                         write(wrn, '(A,I0,A,F0.3)')                           &
                            &'Warning-- All outlets dry or unusable for group '&
                           & , j, ' at day ', jday
                         warning_open = .TRUE.
 
        ! Only nonblended outlets.
                     elseif(NOUTS(j)==ng0)then
                         write(wrn, '(A,I0,A,F0.3)')                           &
                         &'Warning-- Only nonblended outlets present in group '&
                        & , j, ' at day ', jday
                         warning_open = .TRUE.
 
        ! Only one blended outlet.  It gets all of the blended flow, but must not exceed its maximum flow criterion.
                     elseif(ng1 + ng2==1)then
                         jst = JSTSPLT(j, NOUT1(1))
                         QSTRFRAC(jst, jb) = 1.0 - sum_minfrac0
                         if(qall - q_notblended>TSMAXFLOW(j, NOUT1(1)) .AND.   &
                          & TSMAXFLOW(j, NOUT1(1))>0.0)then
                             QSTRFRAC(jst, jb) = TSMAXFLOW(j, NOUT1(1))/qall
                             write(wrn, '(A,A,I0,A,I0,A,F0.3)')                &
     &'Warning-- Total release flow rate decreased to comply with maximum flow &
     &', 'criterion for structure ', jst, ' in group ', j, ' at day ', jday
                             warning_open = .TRUE.
                         endif
 
        ! Minimum flows comprise entire release.  No blending calculations required.
                     elseif(ABS(1.0 - sum_minfrac0 - sum_minfrac1 -            &
                          & sum_minfrac2)<=0.000001)then
                         do n = 1, ng1
                             jst = JSTSPLT(j, NOUT1(n))
                             QSTRFRAC(jst, jb) = MINFRAC1(n)
                         enddo
                         do n = 1, ng2
                             jst = JSTSPLT(j, NOUT2(n))
                             QSTRFRAC(jst, jb) = MINFRAC2(n)
                         enddo
 
        ! More than one usable outlet, and blending among priority 1 outlet(s) and priority 2 outlet(s) required.
                     else
                         id = DS(jb)                                         ! needed for downstream_withdrawal_estimate
                         kt = KTWB(jw)                                       ! needed for downstream_withdrawal_estimate
 
          ! Warn the user if maximum flow criteria are likely to decrease the specified outflows.
                         if(sum_minfrac0 + sum_maxfrac1 + sum_maxfrac2<1.0)then
                             write(wrn, '(A,A,I0,A,F0.3)')                     &
     &'Warning-- Total release flow rate may be decreased to comply with maximu&
     &m flow ', 'criteria for structures in group ', j, ' at day ', jday
                             warning_open = .TRUE.
                         endif
 
          ! Set the initial release fractions.  Ensure that maximum flows are not exceeded.
                         qfrac1 = sum_minfrac1 +                               &
                                & 0.5*(1.0 - sum_minfrac0 - sum_minfrac1 -     &
                                & sum_minfrac2)
                         qfrac1 = MIN(sum_maxfrac1, qfrac1)
                         qfrac2 = 1.0 - sum_minfrac0 - qfrac1
                         if(qfrac2>sum_maxfrac2)then
                             excess_frac = qfrac2 - sum_maxfrac2
                             if(qfrac1 + excess_frac<=sum_maxfrac1)then
                                 qfrac1 = qfrac1 + excess_frac
                             else
                                 qfrac1 = sum_maxfrac1
                             endif
                             qfrac2 = sum_maxfrac2
                         endif
                         j2pref = 1
                         call SET_FLOW_FRACS(j, jb, j2pref)                  ! set flow fractions; redistribute if maxfrac exceeded
 
          ! If priority 2 outlets are not sharing flows, identify the two with the highest and lowest elevations
                         if(.NOT.SHARE_FLOW(j) .AND. ng2>1)then
                             j2hi = 1
                             j2lo = 1
                             jst = JSTSPLT(j, NOUT2(1))
                             maxelev = ESTR(jst, jb)
                             minelev = ESTR(jst, jb)
                             do n = 2, ng2
                                 jst = JSTSPLT(j, NOUT2(n))
                                 if(ESTR(jst, jb)>maxelev)then
                                     maxelev = ESTR(jst, jb)
                                     j2hi = n
                                 elseif(ESTR(jst, jb)<minelev)then
                                     minelev = ESTR(jst, jb)
                                     j2lo = n
                                 endif
                             enddo
                         endif
 
          ! Get weighted blend of all nonblended release temperatures
                         ttarg = TSPLTT(j)
                         if(sum_minfrac0>0.0)then
                             sumtemp = 0.0
                             do n = 1, ng0
                                 jst = JSTSPLT(j, NOUT0(n))
                                 QSTR(jst, jb) = qall*QSTRFRAC(jst, jb)
                                 if(QSTR(jst, jb)>0.0)then
                                     call DOWNSTREAM_WITHDRAWAL_ESTIMATE(jst,  &
                                       & etemp, ESTR(jst, jb))
                                 else
                                     etemp = T2(KSTRSPLT(NOUT0(n)), DS(jb))  ! Use temperature at outlet elevation if no flow
                                 endif
                                 sumtemp = sumtemp + QSTRFRAC(jst, jb)*etemp
                             enddo
                             etemp = sumtemp/sum_minfrac0
                             ttarg = (ttarg - sum_minfrac0*etemp)              &
                                   & /(1.0 - sum_minfrac0)                   ! New temperature target for blended releases
                         endif
 
          ! Need an iterative approach because released T depends on Q
                         lastfrac = qfrac1
                         do jj = 1, 8                                        ! Maximum of eight iterations
                             lastfrac2 = lastfrac
                             lastfrac = qfrac1
 
                             sumtemp = 0.0
                             sumelev = 0.0
                             do n = 1, ng1                                   ! Get weighted temp and elevation for group 1
                                 jst = JSTSPLT(j, NOUT1(n))
                                 QSTR(jst, jb) = qall*QSTRFRAC(jst, jb)
                                 if(QSTR(jst, jb)>0.0)then
                                     call DOWNSTREAM_WITHDRAWAL_ESTIMATE(jst,  &
                                       & etemp, ESTR(jst, jb))
                                 else
                                     etemp = T2(KSTRSPLT(NOUT1(n)), DS(jb))  ! Use temperature at outlet elevation if no flow
                                 endif
                                 if(qfrac1>0.0)then
                                     sumtemp = sumtemp + QSTRFRAC(jst, jb)     &
                                       & *etemp
                                     sumelev = sumelev + QSTRFRAC(jst, jb)     &
                                       & *ESTR(jst, jb)
                                 else
                                     sumtemp = sumtemp + etemp
                                     sumelev = sumelev + ESTR(jst, jb)
                                 endif
                             enddo
                             if(qfrac1>0.0)then
                                 etemp1 = sumtemp/qfrac1                     ! Weighted temperature from group 1 outlets
                                 elev1 = sumelev/qfrac1                      ! Weighted elevation of group 1 outlets
                             else
                                 etemp1 = sumtemp/ng1
                                 elev1 = sumelev/ng1
                             endif
 
                             if(SHARE_FLOW(j) .OR. ng2<2)then                ! Get weighted temp and elevation for group 2
                                 sumtemp = 0.0                               ! ...when flows are shared among outlets
                                 sumelev = 0.0
                                 do n = 1, ng2
                                     jst = JSTSPLT(j, NOUT2(n))
                                     QSTR(jst, jb) = qall*QSTRFRAC(jst, jb)
                                     if(QSTR(jst, jb)>0.0)then
                                         call DOWNSTREAM_WITHDRAWAL_ESTIMATE   &
                                           & (jst, etemp, ESTR(jst, jb))
                                     else
                                         etemp = T2(KSTRSPLT(NOUT2(n)), DS(jb))
                                                                             ! Use temperature at outlet elevation if no flow
                                     endif
                                     if(qfrac2>0.0)then
                                         sumtemp = sumtemp + QSTRFRAC(jst, jb) &
                                           & *etemp
                                         sumelev = sumelev + QSTRFRAC(jst, jb) &
                                           & *ESTR(jst, jb)
                                     else
                                         sumtemp = sumtemp + etemp
                                         sumelev = sumelev + ESTR(jst, jb)
                                     endif
                                 enddo
                                 if(qfrac2>0.0)then
                                     etemp2 = sumtemp/qfrac2                 ! Weighted temperature from group 2 outlets
                                     elev2 = sumelev/qfrac2                  ! Weighted elevation of group 2 outlets
                                 else
                                     etemp2 = sumtemp/ng2
                                     elev2 = sumelev/ng2
                                 endif
 
                             else                                            ! ...and when flows are not shared
                                 if(qfrac2==0.0)then
                                     do n = 1, ng2
                                         jst = JSTSPLT(j, NOUT2(n))
                                         SPLT2T(n)                             &
                                           & = T2(KSTRSPLT(NOUT2(n)), DS(jb))
                                         SPLT2E(n) = ESTR(jst, jb)
                                     enddo
                                 else
                                     do nj = 1, ng2                          ! Find the temperatures produced in group 2
                                         sumtemp = 0.0                       ! by testing when each outlet is preferred
                                         sumelev = 0.0
                                         call SET_FLOW_FRACS2(j, jb, nj)
                                         do n = 1, ng2
                                         jst = JSTSPLT(j, NOUT2(n))
                                         QSTR(jst, jb) = qall*QSTRFRAC(jst, jb)
                                         if(QSTR(jst, jb)>0.0)then
                                         call DOWNSTREAM_WITHDRAWAL_ESTIMATE   &
                                           & (jst, etemp, ESTR(jst, jb))
                                         else
                                         etemp = T2(KSTRSPLT(NOUT2(n)), DS(jb))
                                         endif
                                         sumtemp = sumtemp + QSTRFRAC(jst, jb) &
                                           & *etemp
                                         sumelev = sumelev + QSTRFRAC(jst, jb) &
                                           & *ESTR(jst, jb)
                                         enddo
                                         SPLT2T(nj) = sumtemp/qfrac2
                                         SPLT2E(nj) = sumelev/qfrac2
                                     enddo
                                 endif
                                 j2max = 1
                                 j2min = 1
                                 maxtemp = SPLT2T(1)
                                 mintemp = SPLT2T(1)
                                 do n = 2, ng2
                                     if(SPLT2T(n)>maxtemp)then
                                         maxtemp = SPLT2T(n)
                                         j2max = n
                                     elseif(SPLT2T(n)<mintemp)then
                                         mintemp = SPLT2T(n)
                                         j2min = n
                                     endif
                                 enddo
                                 if(ttarg<etemp1 - 0.001)then                ! need a colder temp from group 2
                                     if(maxtemp - mintemp>0.001)then
                                         etemp2 = SPLT2T(j2min)              ! preferred outlet is the coldest one
                                         elev2 = SPLT2E(j2min)
                                         j2pref = j2min
                                     else
                                         etemp2 = SPLT2T(j2lo)               ! preferred outlet is the lowest one
                                         elev2 = SPLT2E(j2lo)
                                         j2pref = j2lo
                                     endif
                                 elseif(ttarg>etemp1 + 0.001)then            ! need a warmer temp from group 2
                                     if(maxtemp - mintemp>0.001)then
                                         etemp2 = SPLT2T(j2max)              ! preferred outlet is the warmest one
                                         elev2 = SPLT2E(j2max)
                                         j2pref = j2max
                                     else
                                         etemp2 = SPLT2T(j2hi)               ! preferred outlet is the highest one
                                         elev2 = SPLT2E(j2hi)
                                         j2pref = j2hi
                                     endif
                                 else
                                     etemp2 = SPLT2T(1)                      ! if temp is close to target, choose first outlet
                                     elev2 = SPLT2E(1)
                                     j2pref = 1
                                 endif
                             endif
 
            ! Target temperature is less than either outlet temperature.
                             if(ttarg<etemp1 .AND. ttarg<etemp2)then
                                 qfrac1 = sum_minfrac1                       ! default for if/then cases
                                 if(ABS(etemp1 - etemp2)<0.001)then
                                     if(prior1==prior2)then                  ! If each outlet has the same priority level, then...
                                         if(elev1<=elev2)qfrac1 = 1.0 -        &
                                           & sum_minfrac0 - sum_minfrac2     ! Choose lower outlet if both have same temperature.
                                     elseif(prior1<prior2)then
                                         qfrac1 = 1.0 - sum_minfrac0 -         &
                                           & sum_minfrac2                    ! Choose higher priority outlet if temps are the same.
                                     endif
                                 elseif(etemp1<etemp2)then                   ! If temps are different, choose the one closer
                                     qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                                                                             ! to target temperature
                                 endif
 
            ! Target temperature is greater than either outlet temperature.
                             elseif(ttarg>etemp1 .AND. ttarg>etemp2)then
                                 qfrac1 = sum_minfrac1                       ! default for if/then cases
                                 if(ABS(etemp1 - etemp2)<0.001)then
                                     if(prior1==prior2)then                  ! If each outlet has the same priority level, then...
                                         if(elev1>=elev2)qfrac1 = 1.0 -        &
                                           & sum_minfrac0 - sum_minfrac2     ! Choose upper outlet if both have same temperature.
                                     elseif(prior1<prior2)then
                                         qfrac1 = 1.0 - sum_minfrac0 -         &
                                           & sum_minfrac2                    ! Choose higher priority outlet if temps are the same.
                                     endif
                                 elseif(etemp1>etemp2)then                   ! If temps are different, choose the one closer
                                     qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                                                                             ! to target temperature
                                 endif
 
            ! Target temperature is essentially the same as the two outlet temperatures.
                             elseif(ABS(etemp1 - etemp2)<0.001)then
                                 qfrac1 = sum_minfrac1                       ! default for if/then cases
                                 if(prior1==prior2)then                      ! If each outlet has the same priority level, then...
                                     qfrac1 = sum_minfrac1 +                   &
                                       & 0.5*(1.0 - sum_minfrac0 -             &
                                       & sum_minfrac1 - sum_minfrac2)                                 ! Split the flow equally.
                                 elseif(prior1<prior2)then
                                     qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                                                                             ! Choose higher priority outlet if temps are the same.
                                 endif
 
            ! Target temperature is between the two outlet temperatures.
                             else
                                 qfrac1 = (1.0 - sum_minfrac0)                 &
                                   & *ABS((ttarg - etemp2)                     &
                                   & /(etemp1 - etemp2 + nonzero))
                                 qfrac1 = MAX(sum_minfrac1, qfrac1)
                                 qfrac1 = MIN(1.0 - sum_minfrac0 -             &
                                   & sum_minfrac2, qfrac1)
                             endif
                             qfrac1 = MIN(sum_maxfrac1, qfrac1)
                             qfrac2 = 1.0 - sum_minfrac0 - qfrac1
                             if(qfrac2>sum_maxfrac2)then
                                 excess_frac = qfrac2 - sum_maxfrac2
                                 if(qfrac1 + excess_frac<=sum_maxfrac1)then
                                     qfrac1 = qfrac1 + excess_frac
                                 else
                                     qfrac1 = sum_maxfrac1
                                 endif
                                 qfrac2 = sum_maxfrac2
                             endif
 
            ! Set flow fractions for individual outlets and redistribute flows if maximum flow fractions exceeded.
                             call SET_FLOW_FRACS(j, jb, j2pref)
 
            ! Exit the loop if the latest flow fraction calculation agrees with the previous one.
            ! Exit if no flow, because no iteration requried in that case.
                             if(ABS(lastfrac - qfrac1)<tsconv .OR. qall==0.0)  &
                              & exit
                         enddo
 
          ! Check to see if iterative solution did not converge.
                         if(ABS(lastfrac - qfrac1)>=tsconv .AND. qall>0.0)then
                             write(wrn, '(A,F0.3,3(A,F0.4))')                  &
                           &'Flow fraction calculations not converging at day '&
                          & , jday, '  Current: ', qfrac1, ' Last: ', lastfrac,&
                           &' Next-to-last: ', lastfrac2
                             warning_open = .TRUE.
 
            ! Check to see if the iterative solution is unstable.  If so, use priorities to assign releases.
            ! Criteria:  change is at least 0.1 and most recent change is in opposite direction from previous change.
                             if(ABS(lastfrac - qfrac1)>=0.1 .AND.              &
                              & (qfrac1 - lastfrac)*(lastfrac - lastfrac2)<0.0)&
                              & then
                                 qfrac1 = sum_minfrac1                       ! default for if/then cases
                                 if(prior1<prior2)then                       ! group 1 is higher priority
                                     qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                                 else                                        ! else, fulfill minima and split the rest
                                     qfrac1 = sum_minfrac1 +                   &
                                       & 0.5*(1.0 - sum_minfrac0 -             &
                                       & sum_minfrac1 - sum_minfrac2)
                                 endif
                                 qfrac1 = MIN(sum_maxfrac1, qfrac1)
                                 qfrac2 = 1.0 - sum_minfrac0 - qfrac1
                                 if(qfrac2>sum_maxfrac2)then
                                     excess_frac = qfrac2 - sum_maxfrac2
                                     if(qfrac1 + excess_frac<=sum_maxfrac1)then
                                         qfrac1 = qfrac1 + excess_frac
                                     else
                                         qfrac1 = sum_maxfrac1
                                     endif
                                     qfrac2 = sum_maxfrac2
                                 endif
                                 call SET_FLOW_FRACS(j, jb, j2pref)          ! set flow fractions; redistribute if maxfrac exceeded
                             endif
                         endif
                     endif
 
        ! Set final flows to go with the flow fractions.  May not be necessary, but do it anyway.
                     do jj = 1, NOUTS(j)
                         QSTR(JSTSPLT(j, jj), jb)                              &
                           & = qall*QSTRFRAC(JSTSPLT(j, jj), jb)
                     enddo
 
 
      ! Do Withdrawals next
                 elseif(TSPLTCNTR(j)=='      WD')then
                     jwd = JSTSPLT(j, 1)                                     ! assume withdrawals are from same branch and waterbody
                     do jb = 1, nbr
                         if(IWD(jwd)>=US(jb) .AND. IWD(jwd)<=DS(jb))exit
                     enddo
                     do jw = 1, nwb
                         if(jb>=BS(jw) .AND. jb<=BE(jw))exit
                     enddo
                     no_flow(j, :) = .FALSE.
                     num_noflow = 0
                     do jj = 1, NOUTS(j)
                         jwd = JSTSPLT(j, jj)
                         elr = SINA(jb)*DLX(IWD(jwd))*0.5
                         wsel = ELWS(IWD(jwd)) - elr                         ! compute water-surface elevation
                         EWD(jwd) = EWDSAV(jwd)                              ! reset outlet elevation to original
                         if(TSTYPE(j, jj)=="FLOAT")then
                             EWD(jwd) = wsel - TSDEPTH(j, jj)
                         elseif(EWD(jwd)>wsel)then
                             if(ELCONTSPL(j)=='     OFF')then
                                 no_flow(j, jj) = .TRUE.                     ! no flow-- high and dry
                                 num_noflow = num_noflow + 1
                             else
                                 EWD(jwd) = wsel                             ! poor man's floating outlet
                             endif
                         endif
                         if(.NOT.no_flow(j, jj) .AND. TSMINHEAD(j, jj)         &
                          & >0.0 .AND. wsel - EWD(jwd)<TSMINHEAD(j, jj))then
                             no_flow(j, jj) = .TRUE.                         ! minimum head criterion not met -- no flow
                             num_noflow = num_noflow + 1
                         endif
                         if(.NOT.no_flow(j, jj) .AND. TSMAXHEAD(j, jj)         &
                          & >0.0 .AND. wsel - EWD(jwd)>TSMAXHEAD(j, jj))then
                             no_flow(j, jj) = .TRUE.                         ! maximum head criterion exceeded -- no flow
                             num_noflow = num_noflow + 1
                         endif
                         do k = KTWB(jw), KB(IWD(jwd))
                             if(EL(k, IWD(jwd)) - elr<EWD(jwd))exit
                         enddo
                         KSTRSPLT(jj) = MIN(k - 1, KB(IWD(jwd)))
                         QWDFRAC(jwd) = 0.0                                  ! initialize flow fractions
                     enddo
 
        ! Use priority inputs to determine which outlets to use
                     prior1 = -999
                     do jj = 1, NOUTS(j)
                         if(.NOT.no_flow(j, jj) .AND. TSPRIOR(j, jj)>=0)then
                             if(prior1== - 999 .OR. TSPRIOR(j, jj)<prior1)     &
                              & prior1 = TSPRIOR(j, jj)
                         endif
                     enddo
                     prior2 = -999
                     do jj = 1, NOUTS(j)
                         if(.NOT.no_flow(j, jj) .AND. TSPRIOR(j, jj)>=0 .AND.  &
                          & TSPRIOR(j, jj)>prior1)then
                             if(prior2== - 999 .OR. TSPRIOR(j, jj)<prior2)     &
                              & prior2 = TSPRIOR(j, jj)
                         endif
                     enddo
 
        ! Outlets with a priority of -1 get used, but are not blended
                     ng0 = 0
                     q_notblended = 0.0
                     do jj = 1, NOUTS(j)
                         jwd = JSTSPLT(j, jj)
                         if(.NOT.no_flow(j, jj) .AND. TSPRIOR(j, jj)== - 1)then
                             ng0 = ng0 + 1
                             NOUT0(ng0) = jj
                             if(QWD(jwd)>TSMAXFLOW(j, jj) .AND.                &
                              & TSMAXFLOW(j, jj)>0.0)then
                                 q_notblended = q_notblended + TSMAXFLOW(j, jj)
                                 QWDFRAC(jwd) = TSMAXFLOW(j, jj)/qall
                             elseif(qall>0.0)then
                                 q_notblended = q_notblended + QWD(jwd)
                                 QWDFRAC(jwd) = QWD(jwd)/qall
                             endif
                         endif
                     enddo
                     sum_minfrac0 = 0.0
                     if(qall>0.0)sum_minfrac0 = q_notblended/qall
 
        ! Outlets with priority 1 and 2 may be used and blended.
                     ng1 = 0
                     ng2 = 0
                     sum_minfrac1 = 0.0
                     sum_minfrac2 = 0.0
                     sum_maxfrac1 = 0.0
                     sum_maxfrac2 = 0.0
                     do jj = 1, NOUTS(j)
                         if(.NOT.no_flow(j, jj))then
                             if(TSPRIOR(j, jj)==prior1)then
                                 ng1 = ng1 + 1
                                 NOUT1(ng1) = jj
                                 MAXFRAC1(ng1) = 1.0
                                 if(qall>0.0 .AND. TSMAXFLOW(j, jj)>0.0)       &
                                  & MAXFRAC1(ng1)                              &
                                  & = MIN(1.0, TSMAXFLOW(j, jj)/qall)
                                 MINFRAC1(ng1) = TSMINFRAC(j, jj)
                                 if(TSMINFRAC(j, jj)<0.0)then
                                     MINFRAC1(ng1) = 0.0
                                     if(qall>0.0)MINFRAC1(ng1)                 &
                                      & = MIN(1.0, ABS(TSMINFRAC(j, jj))/qall)
                                 endif
                                 if(MINFRAC1(ng1)>MAXFRAC1(ng1))MINFRAC1(ng1)  &
                                  & = MAXFRAC1(ng1)
                                 sum_minfrac1 = sum_minfrac1 + MINFRAC1(ng1)
                                 sum_maxfrac1 = sum_maxfrac1 + MAXFRAC1(ng1)
 
                             elseif(TSPRIOR(j, jj)==prior2)then
                                 ng2 = ng2 + 1
                                 NOUT2(ng2) = jj
                                 MAXFRAC2(ng2) = 1.0
                                 if(qall>0.0 .AND. TSMAXFLOW(j, jj)>0.0)       &
                                  & MAXFRAC2(ng2)                              &
                                  & = MIN(1.0, TSMAXFLOW(j, jj)/qall)
                                 MINFRAC2(ng2) = TSMINFRAC(j, jj)
                                 if(TSMINFRAC(j, jj)<0.0)then
                                     MINFRAC2(ng2) = 0.0
                                     if(qall>0.0)MINFRAC2(ng2)                 &
                                      & = MIN(1.0, ABS(TSMINFRAC(j, jj))/qall)
                                 endif
                                 if(MINFRAC2(ng2)>MAXFRAC2(ng2))MINFRAC2(ng2)  &
                                  & = MAXFRAC2(ng2)
                                 sum_minfrac2 = sum_minfrac2 + MINFRAC2(ng2)
                                 sum_maxfrac2 = sum_maxfrac2 + MAXFRAC2(ng2)
                             endif
                         endif
                     enddo
 
        ! If minimum flows are overspecified, then the priority 2 minimum flow fractions are decreased.
                     if(ng2>0 .AND. sum_minfrac0 + sum_minfrac1 + sum_minfrac2>&
                      & 1.0)then
                         if(sum_minfrac0 + sum_minfrac1>=1.0)then
                             ng2 = 0
                             sum_minfrac2 = 0.0
                         else
                             do n = 1, ng2
                                 MINFRAC2(n) = MINFRAC2(n)                     &
                                   & *(1.0 - sum_minfrac0 - sum_minfrac1)      &
                                   & /sum_minfrac2
                             enddo
                             sum_minfrac2 = 1.0 - sum_minfrac0 - sum_minfrac1
                         endif
                     endif
 
        ! If minimum flows are still overspecified, then the priority 1 minimum flow fractions are decreased.
                     if(ng1>0 .AND. sum_minfrac0 + sum_minfrac1>1.0)then
                         if(sum_minfrac0>=1.0)then
                             ng1 = 0
                             sum_minfrac1 = 0.0
                         else
                             do n = 1, ng1
                                 MINFRAC1(n) = MINFRAC1(n)*(1.0 - sum_minfrac0)&
                                   & /sum_minfrac1
                             enddo
                             sum_minfrac1 = 1.0 - sum_minfrac0
                         endif
                     endif
 
        ! If group 1 has 3 or more outlets and group 2 has no outlets, then redistribute priorities based on elevation.
        ! Keep the highest and lowest elevation outlets in group 1, and put other active outlets into nonblended category
        ! with their minimum flows.  If ties in elevation exist, go with the first in the list.
                     if(ng1>2 .AND. ng2==0)then
                         ng1max = 1
                         ng1min = 1
                         jwd = JSTSPLT(j, NOUT1(1))
                         maxelev = EWD(jwd)
                         minelev = EWD(jwd)
                         do n = 2, ng1
                             jwd = JSTSPLT(j, NOUT1(n))
                             if(EWD(jwd)>maxelev)then
                                 maxelev = EWD(jwd)
                                 ng1max = n
                             elseif(EWD(jwd)<minelev)then
                                 minelev = EWD(jwd)
                                 ng1min = n
                             endif
                         enddo
                         blendfrac = 1.0 - sum_minfrac0 - sum_minfrac1 +       &
                                   & MINFRAC1(ng1max) + MINFRAC1(ng1min)
                         if(MAXFRAC1(ng1max) + MAXFRAC1(ng1min)<blendfrac)then
                             if(sum_maxfrac1<1.0 - sum_minfrac0)then
                                 write(wrn, '(A,I0,A,F0.3)')                   &
                     &'Warning-- Maximum flows for outlets exceeded for group '&
                    & , j, ' at day ', jday
                                 warning_open = .TRUE.
                                 do n = 1, ng1
                                     if(n/=ng1max .AND. n/=ng1min)MINFRAC1(n)  &
                                      & = MAXFRAC1(n)
                                 enddo
                             else
                                 excess_frac = blendfrac - MAXFRAC1(ng1max)    &
                                   & - MAXFRAC1(ng1min)
                                 num_left = ng1 - 2
                                 do nj = 1, ng1                    ! iterative process to redistribute excess flows
                                     if(num_left>0 .AND. excess_frac>0.0)then
                                         addfrac = excess_frac/num_left
                                         do n = 1, ng1
                                         if(n/=ng1max .AND. n/=ng1min .AND.    &
                                           & MAXFRAC1(n) - MINFRAC1(n)>0.00001)&
                                           & then
                                         if(MINFRAC1(n) + addfrac>MAXFRAC1(n)) &
                                           & then
                                         num_left = num_left - 1
                                         excess_frac = excess_frac -           &
                                           & (MAXFRAC1(n) - MINFRAC1(n))
                                         MINFRAC1(n) = MAXFRAC1(n)
                                         else
                                         excess_frac = excess_frac - addfrac
                                         MINFRAC1(n) = MINFRAC1(n) + addfrac
                                         endif
                                         endif
                                         enddo
                                     endif
                                 enddo
                             endif
                         endif
                         do n = 1, ng1                       ! assign the other priority 1 outlets to nonblended status
                             if(n/=ng1max .AND. n/=ng1min)then
                                 ng0 = ng0 + 1
                                 NOUT0(ng0) = NOUT1(n)
                                 jwd = JSTSPLT(j, NOUT1(n))
                                 sum_minfrac0 = sum_minfrac0 + MINFRAC1(n)
                                 q_notblended = q_notblended + qall*MINFRAC1(n)
                                 QWDFRAC(jwd) = MINFRAC1(n)
                             endif
                         enddo
                         ng1 = 1                   ! rearrange outlets-- one in each priority group, but same priority
                         ng2 = 1
                         NOUT2(1) = NOUT1(ng1min)
                         MINFRAC2(1) = MINFRAC1(ng1min)
                         MAXFRAC2(1) = MAXFRAC1(ng1min)
                         sum_minfrac2 = MINFRAC1(ng1min)
                         sum_maxfrac2 = MAXFRAC1(ng1min)
                         NOUT1(1) = NOUT1(ng1max)
                         MINFRAC1(1) = MINFRAC1(ng1max)
                         MAXFRAC1(1) = MAXFRAC1(ng1max)
                         sum_minfrac1 = MINFRAC1(ng1max)
                         sum_maxfrac1 = MAXFRAC1(ng1max)
                         prior2 = prior1
                     endif
 
        ! If only two blended outlets, ensure that they are in separate groups.
                     if(ng1==2 .AND. ng2==0)then
                         ng1 = 1
                         ng2 = 1
                         NOUT2(1) = NOUT1(2)
                         MINFRAC2(1) = MINFRAC1(2)
                         MAXFRAC2(1) = MAXFRAC1(2)
                         sum_minfrac2 = MINFRAC1(2)
                         sum_maxfrac2 = MAXFRAC1(2)
                         sum_minfrac1 = MINFRAC1(1)
                         sum_maxfrac1 = MAXFRAC1(1)
                         prior2 = prior1
                     endif
 
 
        ! Begin the blending decisions.
        ! No usable outlets.  All flow fractions remain at zero.
                     if(NOUTS(j)==num_noflow)then
                         write(wrn, '(A,I0,A,F0.3)')                           &
                            &'Warning-- All outlets dry or unusable for group '&
                           & , j, ' at day ', jday
                         warning_open = .TRUE.
 
        ! Only nonblended outlets.
                     elseif(NOUTS(j)==ng0)then
                         write(wrn, '(A,I0,A,F0.3)')                           &
                         &'Warning-- Only nonblended outlets present in group '&
                        & , j, ' at day ', jday
                         warning_open = .TRUE.
 
        ! Only one blended outlet.  It gets all of the blended flow, but must not exceed its maximum flow criterion.
                     elseif(ng1 + ng2==1)then
                         jwd = JSTSPLT(j, NOUT1(1))
                         QWDFRAC(jwd) = 1.0 - sum_minfrac0
                         if(qall - q_notblended>TSMAXFLOW(j, NOUT1(1)) .AND.   &
                          & TSMAXFLOW(j, NOUT1(1))>0.0)then
                             QWDFRAC(jwd) = TSMAXFLOW(j, NOUT1(1))/qall
                             write(wrn, '(A,A,I0,A,I0,A,F0.3)')                &
     &'Warning-- Total release flow rate decreased to comply with maximum flow &
     &', 'criterion for withdrawal ', jwd, ' in group ', j, ' at day ', jday
                             warning_open = .TRUE.
                         endif
 
        ! Minimum flows comprise entire release.  No blending calculations required.
                     elseif(ABS(1.0 - sum_minfrac0 - sum_minfrac1 -            &
                          & sum_minfrac2)<=0.000001)then
                         do n = 1, ng1
                             jwd = JSTSPLT(j, NOUT1(n))
                             QWDFRAC(jwd) = MINFRAC1(n)
                         enddo
                         do n = 1, ng2
                             jwd = JSTSPLT(j, NOUT2(n))
                             QWDFRAC(jwd) = MINFRAC2(n)
                         enddo
 
        ! More than one usable outlet, and blending among priority 1 outlet(s) and priority 2 outlet(s) required.
                     else
                         i = IWD(JSTSPLT(j, NOUT1(1)))                       ! needed for lateral_withdrawal_estimate
                         kt = KTWB(jw)                                       ! needed for lateral_withdrawal_estimate
 
          ! Warn the user if maximum flow criteria are likely to decrease the specified outflows.
                         if(sum_minfrac0 + sum_maxfrac1 + sum_maxfrac2<1.0)then
                             write(wrn, '(A,A,I0,A,F0.3)')                     &
     &'Warning-- Total release flow rate may be decreased to comply with maximu&
     &m flow ', 'criteria for withdrawals in group ', j, ' at day ', jday
                             warning_open = .TRUE.
                         endif
 
          ! Set the initial release fractions.  Ensure that maximum flows are not exceeded.
                         qfrac1 = sum_minfrac1 +                               &
                                & 0.5*(1.0 - sum_minfrac0 - sum_minfrac1 -     &
                                & sum_minfrac2)
                         qfrac1 = MIN(sum_maxfrac1, qfrac1)
                         qfrac2 = 1.0 - sum_minfrac0 - qfrac1
                         if(qfrac2>sum_maxfrac2)then
                             excess_frac = qfrac2 - sum_maxfrac2
                             if(qfrac1 + excess_frac<=sum_maxfrac1)then
                                 qfrac1 = qfrac1 + excess_frac
                             else
                                 qfrac1 = sum_maxfrac1
                             endif
                             qfrac2 = sum_maxfrac2
                         endif
                         j2pref = 1
                         call SET_FLOW_FRACS(j, jb, j2pref)                  ! set flow fractions; redistribute if maxfrac exceeded
 
          ! If priority 2 outlets are not sharing flows, identify the two with the highest and lowest elevations
                         if(.NOT.SHARE_FLOW(j) .AND. ng2>1)then
                             j2hi = 1
                             j2lo = 1
                             jwd = JSTSPLT(j, NOUT2(1))
                             maxelev = EWD(jwd)
                             minelev = EWD(jwd)
                             do n = 2, ng2
                                 jwd = JSTSPLT(j, NOUT2(n))
                                 if(EWD(jwd)>maxelev)then
                                     maxelev = EWD(jwd)
                                     j2hi = n
                                 elseif(EWD(jwd)<minelev)then
                                     minelev = EWD(jwd)
                                     j2lo = n
                                 endif
                             enddo
                         endif
 
          ! Get weighted blend of all nonblended release temperatures
                         ttarg = TSPLTT(j)
                         if(sum_minfrac0>0.0)then
                             sumtemp = 0.0
                             do n = 1, ng0
                                 jwd = JSTSPLT(j, NOUT0(n))
                                 QWD(jwd) = qall*QWDFRAC(jwd)
                                 if(QWD(jwd)>0.0)then
                                     call LATERAL_WITHDRAWAL_ESTIMATE(jwd,     &
                                       & etemp, EWD(jwd))                    ! Get an estimate of the temperature of the outflow
                                 else
                                     etemp = T2(KSTRSPLT(NOUT0(n)), IWD(jwd))
                                                                             ! Use temperature at outlet elevation if no flow
                                 endif
                                 sumtemp = sumtemp + QWDFRAC(jwd)*etemp
                             enddo
                             etemp = sumtemp/sum_minfrac0
                             ttarg = (ttarg - sum_minfrac0*etemp)              &
                                   & /(1.0 - sum_minfrac0)                   ! New temperature target for blended releases
                         endif
 
          ! Need an iterative approach because released T depends on Q
                         lastfrac = qfrac1
                         do jj = 1, 8                                        ! Maximum of eight iterations
                             lastfrac2 = lastfrac
                             lastfrac = qfrac1
 
                             sumtemp = 0.0
                             sumelev = 0.0
                             do n = 1, ng1                                   ! Get weighted temp and elevation for group 1
                                 jwd = JSTSPLT(j, NOUT1(n))
                                 QWD(jwd) = qall*QWDFRAC(jwd)
                                 if(QWD(jwd)>0.0)then
                                     call LATERAL_WITHDRAWAL_ESTIMATE(jwd,     &
                                       & etemp, EWD(jwd))
                                 else
                                     etemp = T2(KSTRSPLT(NOUT1(n)), IWD(jwd))
                                                                             ! Use temperature at outlet elevation if no flow
                                 endif
                                 if(qfrac1>0.0)then
                                     sumtemp = sumtemp + QWDFRAC(jwd)*etemp
                                     sumelev = sumelev + QWDFRAC(jwd)*EWD(jwd)
                                 else
                                     sumtemp = sumtemp + etemp
                                     sumelev = sumelev + EWD(jwd)
                                 endif
                             enddo
                             if(qfrac1>0.0)then
                                 etemp1 = sumtemp/qfrac1                     ! Weighted temperature from group 1 outlets
                                 elev1 = sumelev/qfrac1                      ! Weighted elevation of group 1 outlets
                             else
                                 etemp1 = sumtemp/ng1
                                 elev1 = sumelev/ng1
                             endif
 
                             if(SHARE_FLOW(j) .OR. ng2<2)then                ! Get weighted temp and elevation for group 2
                                 sumtemp = 0.0                               ! ...when flows are shared among outlets
                                 sumelev = 0.0
                                 do n = 1, ng2
                                     jwd = JSTSPLT(j, NOUT2(n))
                                     QWD(jwd) = qall*QWDFRAC(jwd)
                                     if(QWD(jwd)>0.0)then
                                         call LATERAL_WITHDRAWAL_ESTIMATE(jwd, &
                                           & etemp, EWD(jwd))
                                     else
                                         etemp = T2(KSTRSPLT(NOUT2(n)),        &
                                           & IWD(jwd))                       ! Use temperature at outlet elevation if no flow
                                     endif
                                     if(qfrac2>0.0)then
                                         sumtemp = sumtemp + QWDFRAC(jwd)*etemp
                                         sumelev = sumelev + QWDFRAC(jwd)      &
                                           & *EWD(jwd)
                                     else
                                         sumtemp = sumtemp + etemp
                                         sumelev = sumelev + EWD(jwd)
                                     endif
                                 enddo
                                 if(qfrac2>0.0)then
                                     etemp2 = sumtemp/qfrac2                 ! Weighted temperature from group 2 outlets
                                     elev2 = sumelev/qfrac2                  ! Weighted elevation of group 2 outlets
                                 else
                                     etemp2 = sumtemp/ng2
                                     elev2 = sumelev/ng2
                                 endif
 
                             else                                            ! ...and when flows are not shared
                                 if(qfrac2==0.0)then
                                     do n = 1, ng2
                                         jwd = JSTSPLT(j, NOUT2(n))
                                         SPLT2T(n)                             &
                                           & = T2(KSTRSPLT(NOUT2(n)), IWD(jwd))
                                         SPLT2E(n) = EWD(jwd)
                                     enddo
                                 else
                                     do nj = 1, ng2                          ! Find the temperatures produced in group 2
                                         sumtemp = 0.0                       ! by testing when each outlet is preferred
                                         sumelev = 0.0
                                         call SET_FLOW_FRACS2(j, jb, nj)
                                         do n = 1, ng2
                                         jwd = JSTSPLT(j, NOUT2(n))
                                         QWD(jwd) = qall*QWDFRAC(jwd)
                                         if(QWD(jwd)>0.0)then
                                         call LATERAL_WITHDRAWAL_ESTIMATE(jwd, &
                                           & etemp, EWD(jwd))
                                         else
                                         etemp = T2(KSTRSPLT(NOUT2(n)),        &
                                           & IWD(jwd))
                                         endif
                                         sumtemp = sumtemp + QWDFRAC(jwd)*etemp
                                         sumelev = sumelev + QWDFRAC(jwd)      &
                                           & *EWD(jwd)
                                         enddo
                                         SPLT2T(nj) = sumtemp/qfrac2
                                         SPLT2E(nj) = sumelev/qfrac2
                                     enddo
                                 endif
                                 j2max = 1
                                 j2min = 1
                                 maxtemp = SPLT2T(1)
                                 mintemp = SPLT2T(1)
                                 do n = 2, ng2
                                     if(SPLT2T(n)>maxtemp)then
                                         maxtemp = SPLT2T(n)
                                         j2max = n
                                     elseif(SPLT2T(n)<mintemp)then
                                         mintemp = SPLT2T(n)
                                         j2min = n
                                     endif
                                 enddo
                                 if(ttarg<etemp1 - 0.001)then                ! need a colder temp from group 2
                                     if(maxtemp - mintemp>0.001)then
                                         etemp2 = SPLT2T(j2min)              ! preferred outlet is the coldest one
                                         elev2 = SPLT2E(j2min)
                                         j2pref = j2min
                                     else
                                         etemp2 = SPLT2T(j2lo)               ! preferred outlet is the lowest one
                                         elev2 = SPLT2E(j2lo)
                                         j2pref = j2lo
                                     endif
                                 elseif(ttarg>etemp1 + 0.001)then            ! need a warmer temp from group 2
                                     if(maxtemp - mintemp>0.001)then
                                         etemp2 = SPLT2T(j2max)              ! preferred outlet is the warmest one
                                         elev2 = SPLT2E(j2max)
                                         j2pref = j2max
                                     else
                                         etemp2 = SPLT2T(j2hi)               ! preferred outlet is the highest one
                                         elev2 = SPLT2E(j2hi)
                                         j2pref = j2hi
                                     endif
                                 else
                                     etemp2 = SPLT2T(1)                      ! if temp is close to target, choose first outlet
                                     elev2 = SPLT2E(1)
                                     j2pref = 1
                                 endif
                             endif
 
            ! Target temperature is less than either outlet temperature.
                             if(ttarg<etemp1 .AND. ttarg<etemp2)then
                                 qfrac1 = sum_minfrac1                       ! default for if/then cases
                                 if(ABS(etemp1 - etemp2)<0.001)then
                                     if(prior1==prior2)then                  ! If each outlet has the same priority level, then...
                                         if(elev1<=elev2)qfrac1 = 1.0 -        &
                                           & sum_minfrac0 - sum_minfrac2     ! Choose lower outlet if both have same temperature.
                                     elseif(prior1<prior2)then
                                         qfrac1 = 1.0 - sum_minfrac0 -         &
                                           & sum_minfrac2                    ! Choose higher priority outlet if temps are the same.
                                     endif
                                 elseif(etemp1<etemp2)then                   ! If temps are different, choose the one closer
                                     qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                                                                             ! to target temperature
                                 endif
 
            ! Target temperature is greater than either outlet temperature.
                             elseif(ttarg>etemp1 .AND. ttarg>etemp2)then
                                 qfrac1 = sum_minfrac1                       ! default for if/then cases
                                 if(ABS(etemp1 - etemp2)<0.001)then
                                     if(prior1==prior2)then                  ! If each outlet has the same priority level, then...
                                         if(elev1>=elev2)qfrac1 = 1.0 -        &
                                           & sum_minfrac0 - sum_minfrac2     ! Choose upper outlet if both have same temperature.
                                     elseif(prior1<prior2)then
                                         qfrac1 = 1.0 - sum_minfrac0 -         &
                                           & sum_minfrac2                    ! Choose higher priority outlet if temps are the same.
                                     endif
                                 elseif(etemp1>etemp2)then                   ! If temps are different, choose the one closer
                                     qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                                                                             ! to target temperature
                                 endif
 
            ! Target temperature is essentially the same as the two outlet temperatures.
                             elseif(ABS(etemp1 - etemp2)<0.001)then
                                 qfrac1 = sum_minfrac1                       ! default for if/then cases
                                 if(prior1==prior2)then                      ! If each outlet has the same priority level, then...
                                     qfrac1 = sum_minfrac1 +                   &
                                       & 0.5*(1.0 - sum_minfrac0 -             &
                                       & sum_minfrac1 - sum_minfrac2)                                 ! Split the flow equally.
                                 elseif(prior1<prior2)then
                                     qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                                                                             ! Choose higher priority outlet if temps are the same.
                                 endif
 
            ! Target temperature is between the two outlet temperatures.
                             else
                                 qfrac1 = (1.0 - sum_minfrac0)                 &
                                   & *ABS((ttarg - etemp2)                     &
                                   & /(etemp1 - etemp2 + nonzero))
                                 qfrac1 = MAX(sum_minfrac1, qfrac1)
                                 qfrac1 = MIN(1.0 - sum_minfrac0 -             &
                                   & sum_minfrac2, qfrac1)
                             endif
                             qfrac1 = MIN(sum_maxfrac1, qfrac1)
                             qfrac2 = 1.0 - sum_minfrac0 - qfrac1
                             if(qfrac2>sum_maxfrac2)then
                                 excess_frac = qfrac2 - sum_maxfrac2
                                 if(qfrac1 + excess_frac<=sum_maxfrac1)then
                                     qfrac1 = qfrac1 + excess_frac
                                 else
                                     qfrac1 = sum_maxfrac1
                                 endif
                                 qfrac2 = sum_maxfrac2
                             endif
 
            ! Set flow fractions for individual outlets and redistribute flows if maximum flow fractions exceeded.
                             call SET_FLOW_FRACS(j, jb, j2pref)
 
            ! Exit the loop if the latest flow fraction calculation agrees with the previous one.
            ! Exit if no flow, because no iteration requried in that case.
                             if(ABS(lastfrac - qfrac1)<tsconv .OR. qall==0.0)  &
                              & exit
                         enddo
 
          ! Check to see if iterative solution did not converge.
                         if(ABS(lastfrac - qfrac1)>=tsconv .AND. qall>0.0)then
                             write(wrn, '(A,F0.3,3(A,F0.4))')                  &
                           &'Flow fraction calculations not converging at day '&
                          & , jday, '  Current: ', qfrac1, ' Last: ', lastfrac,&
                           &' Next-to-last: ', lastfrac2
                             warning_open = .TRUE.
 
            ! Check to see if the iterative solution is unstable.  If so, use priorities to assign releases.
            ! Criteria:  change is at least 0.1 and most recent change is in opposite direction from previous change.
                             if(ABS(lastfrac - qfrac1)>=0.1 .AND.              &
                              & (qfrac1 - lastfrac)*(lastfrac - lastfrac2)<0.0)&
                              & then
                                 qfrac1 = sum_minfrac1                       ! default for if/then cases
                                 if(prior1<prior2)then                       ! group 1 is higher priority
                                     qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                                 else                                        ! else, fulfill minima and split the rest
                                     qfrac1 = sum_minfrac1 +                   &
                                       & 0.5*(1.0 - sum_minfrac0 -             &
                                       & sum_minfrac1 - sum_minfrac2)
                                 endif
                                 qfrac1 = MIN(sum_maxfrac1, qfrac1)
                                 qfrac2 = 1.0 - sum_minfrac0 - qfrac1
                                 if(qfrac2>sum_maxfrac2)then
                                     excess_frac = qfrac2 - sum_maxfrac2
                                     if(qfrac1 + excess_frac<=sum_maxfrac1)then
                                         qfrac1 = qfrac1 + excess_frac
                                     else
                                         qfrac1 = sum_maxfrac1
                                     endif
                                     qfrac2 = sum_maxfrac2
                                 endif
                                 call SET_FLOW_FRACS(j, jb, j2pref)          ! set flow fractions; redistribute if maxfrac exceeded
                             endif
                         endif
                     endif
 
        ! Set final flows to go with the flow fractions.  May not be necessary, but do it anyway.
                     do jj = 1, NOUTS(j)
                         QWD(JSTSPLT(j, jj)) = qall*QWDFRAC(JSTSPLT(j, jj))
                     enddo
                 endif
             endif
         enddo
 
         nxtsplit = nxtsplit + tspltfreq
     endif
 
!    Use the flow fractions to set flows in blended groups.
     if(tspltc=='      ON')then
         do j = 1, numtsplt
             if(TSYEARLY(j)=='     OFF')then
                 daytest = jday
             else
                 daytest = REAL(jdayg) + jday - INT(jday)
             endif
             if(daytest>=TSTSRT(j) .AND. daytest<TSTEND(j))then
                 qall = 0.0
 
      ! Do structures first
                 if(TSPLTCNTR(j)=='      ST')then
                     do jj = 1, NOUTS(j)
                         qall = qall + QSTR(JSTSPLT(j, jj), TSPLTJB(j))             ! sum up all the flows
                     enddo
                     do jj = 1, NOUTS(j)                                            ! set the flows and honor the maximum flow
                         jst = JSTSPLT(j, jj)
                         QSTR(jst, TSPLTJB(j)) = QSTRFRAC(jst, TSPLTJB(j))*qall
                         if(TSMAXFLOW(j, jj)>0.0 .AND. QSTR(jst, TSPLTJB(j))   &
                          & >TSMAXFLOW(j, jj))QSTR(jst, TSPLTJB(j))            &
                          & = TSMAXFLOW(j, jj)
                     enddo
 
      ! Do Withdrawals next
                 elseif(TSPLTCNTR(j)=='      WD')then
                     do jj = 1, NOUTS(j)
                         qall = qall + QWD(JSTSPLT(j, jj))                          ! sum up all the flows
                     enddo
                     do jj = 1, NOUTS(j)                                            ! set the flows and honor the maximum flow
                         jwd = JSTSPLT(j, jj)
                         QWD(jwd) = QWDFRAC(jwd)*qall
                         if(TSMAXFLOW(j, jj)>0.0 .AND. QWD(jwd)                &
                          & >TSMAXFLOW(j, jj))QWD(jwd) = TSMAXFLOW(j, jj)
                     enddo
                 endif
             endif
         enddo
     endif
 
!    Output some results.
     if(jday>=nxtstr)then
         nxtstr = nxtstr + tfrqtmp
         ifile = 1949
         do jb = 1, nbr
             if(NSTR(jb)>0)then
                 ifile = ifile + 1
                 write(fmtstr, '(I20)') nstr(jb)
                 write(ifile,                                                  &
                     &'(f10.4,'//trim(fmtstr)//'f10.2,'//trim(fmtstr)//'f10.2,'//trim(fmtstr)//'f10.2)'&
                    & )jday, (TAVG(i, jb), i = 1, NSTR(jb)),                   &
                     & (QSTR(i, jb), i = 1, NSTR(jb)),                         &
                     & (ESTR(i, jb), i = 1, NSTR(jb))
             endif
         enddo
         if(nwd>0)then
             ifile = ifile + 1
             write(fmtstr, '(I20)') nwd
             write(ifile, '(f10.4,'//trim(fmtstr)//'f10.2,'//trim(fmtstr)//'f10.2,'//trim(fmtstr)//'f10.2)')jday,     &
                 & (TAVGW(i), i = 1, nwd), (QWD(i), i = 1, nwd),               &
                 & (EWD(i), i = 1, nwd)
         endif
 
  ! computing reservoir volume and volume below 'tempcrit'
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
 
!    Check elevations and status of temperature control towers.
     if(tempc=='      ON' .AND. jday>=nxttcd)then
 
  ! Update the temperature targets
         do j = 1, numtempc
             if(DYNSEL(j)=='      ON')then
                 do while (jday>=NXSEL(j))
                     TCTEMP(j) = TEMP2(j)
                     read(SELD(j), '(1000F8.0)')NXSEL(j), TEMP2(j)
                 enddo
             endif
         enddo
 
         do j = 1, numtempc
 
    ! Structures
             if(TCNTR(j)=='      ST')then
                 js = TCJS(j)                        ! set structure index
                 jb = TCJB(j)                        ! set branch index
                 do jw = 1, nwb                      ! set waterbody index
                     if(jb>=BS(jw) .AND. jb<=BE(jw))exit
                 enddo
 
                 if(TCISEG(j)==0)then
                     tcomp = TAVG(js, jb)   !cb 9/8/06
                 elseif(TCISEG(j)<0)then
                     tcomp = TWDO(ABS(TCISEG(j)))
                                            ! sw 11/26/10
                 else
 
        ! Check to see if the monitoring segment tciseg is in the same branch and waterbody as the structure
                     do jjb = 1, nbr
                         if(TCISEG(j)>=US(jjb) .AND. TCISEG(j)<=DS(jjb))exit
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
                     ESTR(js, jb) = TCELEV(j, NCOUNTC(js, jb))    ! initialize the structure elevation
 
                     if(tcomp>TCTEMP(j) .AND. TCNELEV(j)>NCOUNTC(js, jb))then
          ! making sure that the next lower structure for a particular 'j' is found
                         do nj = NCOUNTC(js, jb) + 1, TCNELEV(j)
                             if(TCELEV(j, nj)<ESTR(js, jb))then
                                 NCOUNTC(js, jb) = nj
                                 ESTR(js, jb) = TCELEV(j, NCOUNTC(js, jb))
                                 exit
                             endif
                         enddo
 
                     elseif(tcomp<TCTEMP(j) .AND. NCOUNTC(js, jb)>1)then
          ! to prevent this happening at each time it checks it and hence oscillating back and forth - check the temp at the upper outlet also
                         if(TCISEG(j)>0)then
                             if(jb==jjb)then
                                 wsel = ELWS(DS(jb)) - SINA(jb)*DLX(DS(jb))*0.5           ! compute water-surface elevation !SR 03/24/13
                                 do ks = KTWB(jw), KB(DS(jb))
                ! if (depthb(ks,tciseg(j)) > tcelev(j,ncountc(js,jb)-1)) exit             !??can't be right-- SR 03/24/13
                                     if(wsel - DEPTHB(ks, TCISEG(j))           &
                                      & <TCELEV(j, NCOUNTC(js, jb) - 1))exit              !SR 03/24/13
                                 enddo
                                 ks = MIN(ks, KB(TCISEG(j)))
                                 tmod = T2(ks, DS(jb))
                             else
                                 tmod = T2(k, TCISEG(j))
                             endif
                             if(tmod<TCTEMP(j) .AND.                           &
                              & TCELEV(j, NCOUNTC(js, jb) - 1)<ELWS(DS(jb)))   &
                              & then
              ! making sure that the next upper structure for a particular 'j' is found
                                 do nj = NCOUNTC(js, jb) - 1, 1, -1
                                     if(TCELEV(j, nj)>ESTR(js, jb))then
                                         NCOUNTC(js, jb) = nj
                                         ESTR(js, jb)                          &
                                           & = TCELEV(j, NCOUNTC(js, jb))
                                         exit
                                     endif
                                 enddo
                             endif
 
                         elseif(TCISEG(j)==0)then
            ! calculate the estimated outflow temperature at higher ports when tcomp < tctemp(j),
            ! and move up if higher port still meets criteria - this doesn't happen when tciseg < 0
                             do nj = 1, NCOUNTC(js, jb) - 1
                                 id = DS(jb)
                                 kt = KTWB(jw)
                                 call DOWNSTREAM_WITHDRAWAL_ESTIMATE(js,       &
                                   & tempest, TCELEV(j, nj))
                                 if(tempest<TCTEMP(j) .AND. TCELEV(j, nj)      &
                                  & <ELWS(DS(jb)))then
                                     NCOUNTC(js, jb) = nj
                                     ESTR(js, jb) = TCELEV(j, NCOUNTC(js, jb))
                                     exit
                                 endif
                             enddo
                         endif
                     endif
                     if(TCELEVCON(j)=='      ON' .AND. TCNELEV(j)              &
                      & >NCOUNTC(js, jb) .AND. ESTR(js, jb)>ELWS(DS(jb)))then
                         NCOUNTC(js, jb) = NCOUNTC(js, jb) + 1
                         ESTR(js, jb) = TCELEV(j, NCOUNTC(js, jb))
                     endif
                 endif
 
    ! Withdrawals
             elseif(TCNTR(j)=='      WD')then
                 jwd = TCJS(j)
                 if(TCISEG(j)==0)then
        ! tcomp = tout(jb)
                     tcomp = TAVGW(TCJS(j))
                                   !cb 9/8/06
                 elseif(TCISEG(j)<0)then
                     tcomp = TWDO(ABS(TCISEG(j)))
                 else
 
        ! checking to see if the monitoring segment tciseg is in the same branch and water body as the withdrawal
                     do jjb = 1, nbr
                         if(TCISEG(j)>=US(jjb) .AND. TCISEG(j)<=DS(jjb))exit
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
                     EWD(jwd) = TCELEV(j, NCOUNTCW(jwd))      ! initialize the withdrawal elevation
 
                     if(tcomp>TCTEMP(j) .AND. TCNELEV(j)>NCOUNTCW(jwd))then
          ! making sure that the next lower structure for a particular 'j' is found
                         do nj = NCOUNTCW(jwd) + 1, TCNELEV(j)
                             if(TCELEV(j, nj)<EWD(jwd))then
                                 NCOUNTCW(jwd) = nj
                                 EWD(jwd) = TCELEV(j, NCOUNTCW(jwd))
                                 exit
                             endif
                         enddo
 
                     elseif(tcomp<TCTEMP(j) .AND. NCOUNTCW(jwd)>1)then
          ! to prevent this happening at each time it checks it and hence oscillating back and forth - check the temp at the upper outlet also
                         if(TCISEG(j)>0)then
                             tmod = T2(k, TCISEG(j))
                             if(tmod<TCTEMP(j) .AND.                           &
                              & TCELEV(j, NCOUNTCW(jwd) - 1)<ELWS(IWD(jwd)))   &
                              & then
              ! making sure that the next upper structure for a particular 'j' is found
                                 do nj = NCOUNTCW(jwd) - 1, 1, -1
                                     if(TCELEV(j, nj)>EWD(jwd))then
                                         NCOUNTCW(jwd) = nj
                                         EWD(jwd) = TCELEV(j, NCOUNTCW(jwd))
                                         exit
                                     endif
                                 enddo
                             endif
 
                         elseif(TCISEG(j)==0)then
            ! calculate estimated outflow temperature at higher ports when tcomp < tctemp(j),
            ! and move up if higher port still meets criteria
                             i = MAX(CUS(JBWD(jwd)), IWD(jwd))
                             do jjb = 1, nbr
                                 if(i>=US(jjb) .AND. i<=DS(jjb))exit
                             enddo
                             do jjw = 1, nwb
                                 if(jjb>=BS(jjw) .AND. jjb<=BE(jjw))exit
                             enddo
                             kt = KTWB(jjw)
                             do nj = 1, NCOUNTCW(jwd) - 1
                                 call LATERAL_WITHDRAWAL_ESTIMATE(jwd, tempest,&
                                   & TCELEV(j, nj))
                                 if(tempest<TCTEMP(j) .AND. TCELEV(j, nj)      &
                                  & <ELWS(IWD(jwd)))then
                                     NCOUNTCW(jwd) = nj
                                     EWD(jwd) = TCELEV(j, NCOUNTCW(jwd))
                                     exit
                                 endif
                             enddo
                         endif
                     endif
                     if(TCELEVCON(j)=='      ON' .AND. TCNELEV(j)>NCOUNTCW(jwd)&
                      & .AND. EWD(jwd)>ELWS(IWD(jwd)))then
                         NCOUNTCW(jwd) = NCOUNTCW(jwd) + 1
                         EWD(jwd) = TCELEV(j, NCOUNTCW(jwd))
                     endif
                 endif
             endif
         enddo
 
         nxttcd = nxttcd + tcdfreq
     endif
     return
!
!
     entry DEALLOCATE_SELECTIVEUSGS
     deallocate(TCNELEV, TCJB, TCJS, TCELEV, TCTEMP, TCTEND, TCTSRT, NCOUNTC,  &
              & TCISEG, TCKLAY, TCELEVCON, ELCONTSPL)
     deallocate(TSPLTJB, TSPLTT, NOUTS, JSTSPLT, KSTRSPLT, TCYEARLY, TCNTR,    &
              & TSPLTCNTR)
     deallocate(volm, NCOUNTCW, QWDFRAC, QSTRFRAC)
     deallocate(TEMPCRIT, volmc, DYNSEL, SELD, NXSEL, TEMP2, TSYEARLY, TSTEND, &
              & TSTSRT)
     deallocate(TSDEPTH, TSTYPE, TSMINFRAC, TSPRIOR, TSMINHEAD, TSMAXHEAD,     &
              & TSMAXFLOW, no_flow)
     deallocate(TSDYNSEL, TSSELD, NXTSSEL, TSTEMP2, EWDSAV, ESTRSAV,           &
              & SHARE_FLOW, wd_active, str_active)
     deallocate(NOUT0, NOUT1, NOUT2, MINFRAC1, MINFRAC2, MAXFRAC1, MAXFRAC2,   &
              & SPLT2T, SPLT2E)
 
     end subroutine SELECTIVEUSGS
