!*==set_flow_fractions.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************************************************************************
!**                                                S E T _ F L O W _ F R A C T I O N S                                            **
!***********************************************************************************************************************************
                                                                              ! Entire routine added/modified by S. Rounds, 06/26/13
     subroutine SET_FLOW_FRACTIONS
                                                                              ! This routine sets the flow fractions, and then
     use SELECTIVE1USGS                                                       ! redistributes flow to other outlets in the same
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     integer :: J, J2pref, Jb
     intent (inout) J2pref
!
! Local variables
!
     real :: addfrac, excess_frac
     integer :: jst, jwd, n, nexcess, nj
!
!*** End of declarations rewritten by SPAG
!
                                                                              ! group when one or more outlets exceed their maximum
                                                                              ! flow rates.  Excess flow that cannot be accommodated
     return                                                                   ! within each group will be discarded.
 
 
     entry SET_FLOW_FRACS(J, Jb, J2pref)
     excess_frac = 0.0                                                        ! Excess flow above maximum release rates
     nexcess = 0                                                              ! Number of group 1 outlets exceeding maximum rates
 
!    Set and rebalance release fractions for group 1 structures
     if(TSPLTCNTR(J)=='      ST')then
         do n = 1, ng1                                                        ! Find maxed-out outlets in group 1; set flows to max
             jst = JSTSPLT(J, NOUT1(n))
             QSTRFRAC(jst, Jb) = MINFRAC1(n) + (qfrac1 - sum_minfrac1)/ng1
             if(QSTRFRAC(jst, Jb)>MAXFRAC1(n))then
                 nexcess = nexcess + 1
                 excess_frac = excess_frac + QSTRFRAC(jst, Jb) - MAXFRAC1(n)
                 QSTRFRAC(jst, Jb) = MAXFRAC1(n)
             endif
         enddo
         if(excess_frac>0.0 .AND. ng1>nexcess)then                            ! Redistribute excess flow to other outlets in group
             do nj = 1, ng1                                                   ! Iterative process, in case others get maxed-out
                 if(ng1==nexcess .OR. excess_frac<=0.00001)exit
                 addfrac = excess_frac/(ng1 - nexcess)
                 do n = 1, ng1
                     jst = JSTSPLT(J, NOUT1(n))
                     if(MAXFRAC1(n) - QSTRFRAC(jst, Jb)>0.00001)then
                         if(QSTRFRAC(jst, Jb) + addfrac>MAXFRAC1(n))then
                             nexcess = nexcess + 1
                             excess_frac = excess_frac -                       &
                               & (MAXFRAC1(n) - QSTRFRAC(jst, Jb))
                             QSTRFRAC(jst, Jb) = MAXFRAC1(n)
                         else
                             excess_frac = excess_frac - addfrac
                             QSTRFRAC(jst, Jb) = QSTRFRAC(jst, Jb) + addfrac
                         endif
                     endif
                 enddo
             enddo
         endif
 
!        Set and rebalance release fractions for group 1 withdrawals
     else
         do n = 1, ng1                                                        ! Find maxed-out outlets in group 1; set flows to max
             jwd = JSTSPLT(J, NOUT1(n))
             QWDFRAC(jwd) = MINFRAC1(n) + (qfrac1 - sum_minfrac1)/ng1
             if(QWDFRAC(jwd)>MAXFRAC1(n))then
                 nexcess = nexcess + 1
                 excess_frac = excess_frac + QWDFRAC(jwd) - MAXFRAC1(n)
                 QWDFRAC(jwd) = MAXFRAC1(n)
             endif
         enddo
         if(excess_frac>0.0 .AND. ng1>nexcess)then                            ! Redistribute excess flow to other outlets in group
             do nj = 1, ng1                                                   ! Iterative process, in case others get maxed-out
                 if(ng1==nexcess .OR. excess_frac<=0.00001)exit
                 addfrac = excess_frac/(ng1 - nexcess)
                 do n = 1, ng1
                     jwd = JSTSPLT(J, NOUT1(n))
                     if(MAXFRAC1(n) - QWDFRAC(jwd)>0.00001)then
                         if(QWDFRAC(jwd) + addfrac>MAXFRAC1(n))then
                             nexcess = nexcess + 1
                             excess_frac = excess_frac -                       &
                               & (MAXFRAC1(n) - QWDFRAC(jwd))
                             QWDFRAC(jwd) = MAXFRAC1(n)
                         else
                             excess_frac = excess_frac - addfrac
                             QWDFRAC(jwd) = QWDFRAC(jwd) + addfrac
                         endif
                     endif
                 enddo
             enddo
         endif
     endif
 
     entry SET_FLOW_FRACS2(J, Jb, J2pref)                                     ! Separate entry just for group 2 outlets
 
     if(J2pref==0)J2pref = 1                                                  ! Preferred outlet number, if not sharing
     excess_frac = 0.0                                                        ! Excess flow above maximum release rates
     nexcess = 0                                                              ! Number of group 2 outlets exceeding maximum rates
 
!    Set and rebalance release fractions for group 2 structures
     if(TSPLTCNTR(J)=='      ST')then
         do n = 1, ng2                                                        ! Find maxed-out outlets in group 2; set flows to max
             jst = JSTSPLT(J, NOUT2(n))
             if(.NOT.(.NOT.SHARE_FLOW(J) .AND. ng2>1))then                    ! Direct flow to preferred outlet if not shared
                 QSTRFRAC(jst, Jb) = MINFRAC2(n) + (qfrac2 - sum_minfrac2)/ng2
             elseif(n==J2pref)then
                 QSTRFRAC(jst, Jb) = qfrac2 - sum_minfrac2 + MINFRAC2(n)
             else
                 QSTRFRAC(jst, Jb) = MINFRAC2(n)
             endif
             if(QSTRFRAC(jst, Jb)>MAXFRAC2(n))then
                 nexcess = nexcess + 1
                 excess_frac = excess_frac + QSTRFRAC(jst, Jb) - MAXFRAC2(n)
                 QSTRFRAC(jst, Jb) = MAXFRAC2(n)
             endif
         enddo
         if(excess_frac>0.0 .AND. ng2>nexcess)then                            ! Redistribute excess flow to other outlets in group
             do nj = 1, ng2                                                   ! Iterative process, in case others get maxed-out
                 if(ng2==nexcess .OR. excess_frac<=0.00001)exit
                 addfrac = excess_frac/(ng2 - nexcess)
                 do n = 1, ng2
                     jst = JSTSPLT(J, NOUT2(n))
                     if(MAXFRAC2(n) - QSTRFRAC(jst, Jb)>0.00001)then
                         if(QSTRFRAC(jst, Jb) + addfrac>MAXFRAC2(n))then
                             nexcess = nexcess + 1
                             excess_frac = excess_frac -                       &
                               & (MAXFRAC2(n) - QSTRFRAC(jst, Jb))
                             QSTRFRAC(jst, Jb) = MAXFRAC2(n)
                         else
                             excess_frac = excess_frac - addfrac
                             QSTRFRAC(jst, Jb) = QSTRFRAC(jst, Jb) + addfrac
                         endif
                     endif
                 enddo
             enddo
         endif
 
!        Set and rebalance release fractions for group 2 withdrawals
     elseif(TSPLTCNTR(J)=='      WD')then
         do n = 1, ng2                                                        ! Find maxed-out outlets in group 2; set flows to max
             jwd = JSTSPLT(J, NOUT2(n))
             if(.NOT.(.NOT.SHARE_FLOW(J) .AND. ng2>1))then                    ! Direct flow to preferred outlet if not shared
                 QWDFRAC(jwd) = MINFRAC2(n) + (qfrac2 - sum_minfrac2)/ng2
             elseif(n==J2pref)then
                 QWDFRAC(jwd) = qfrac2 - sum_minfrac2 + MINFRAC2(n)
             else
                 QWDFRAC(jwd) = MINFRAC2(n)
             endif
             if(QWDFRAC(jwd)>MAXFRAC2(n))then
                 nexcess = nexcess + 1
                 excess_frac = excess_frac + QWDFRAC(jwd) - MAXFRAC2(n)
                 QWDFRAC(jwd) = MAXFRAC2(n)
             endif
         enddo
         if(excess_frac>0.0 .AND. ng2>nexcess)then                            ! Redistribute excess flow to other outlets in group
             do nj = 1, ng2                                                   ! Iterative process, in case others get maxed-out
                 if(ng2==nexcess .OR. excess_frac<=0.00001)exit
                 addfrac = excess_frac/(ng2 - nexcess)
                 do n = 1, ng2
                     jwd = JSTSPLT(J, NOUT2(n))
                     if(MAXFRAC2(n) - QWDFRAC(jwd)>0.00001)then
                         if(QWDFRAC(jwd) + addfrac>MAXFRAC2(n))then
                             nexcess = nexcess + 1
                             excess_frac = excess_frac -                       &
                               & (MAXFRAC2(n) - QWDFRAC(jwd))
                             QWDFRAC(jwd) = MAXFRAC2(n)
                         else
                             excess_frac = excess_frac - addfrac
                             QWDFRAC(jwd) = QWDFRAC(jwd) + addfrac
                         endif
                     endif
                 enddo
             enddo
         endif
     endif
     end subroutine SET_FLOW_FRACTIONS
