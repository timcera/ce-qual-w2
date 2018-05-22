!*==normal_depth.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**        S U B R O U T I N E    N O R M A L    D E P T H                                                                        **
!***********************************************************************************************************************************
 
     subroutine NORMAL_DEPTH(Flow)
     use GLOBAL
     use GEOMC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
     integer, parameter :: JMAX = 40
!
! Dummy arguments
!
     real(r8) :: Flow
!
! Local variables
!
     real :: dx, fmid, func1, funcval1, funcval2, r8, rtbis, x1, x2, xacc, xmid
     real :: EL
     real :: ELWS
     real :: GEOMC, GLOBAL
     integer :: i, j, jj
     integer :: KB
!
!*** End of declarations rewritten by SPAG
!
 
 
 
 
!    FIRST, BRACKETING ROOT
     x1 = 0.001
     x2 = 1.0
     call MANNINGS_EQN(Flow, x1, funcval1)
     call MANNINGS_EQN(Flow, x2, funcval2)
 
     do jj = 1, JMAX
         if(funcval1*funcval2<=0.0)exit
         if(ABS(funcval1)<ABS(funcval2))then
             x1 = x1/2.0
             call MANNINGS_EQN(Flow, x1, funcval1)
         else
             x2 = x2 + 1.5*(x2 - x1)
             call MANNINGS_EQN(Flow, x2, funcval2)
         endif
     enddo
 
!    FINDING ROOT BY BISECTION
     xacc = 0.01
     call MANNINGS_EQN(Flow, x2, fmid)
     call MANNINGS_EQN(Flow, x1, func1)
  !    IF(FUNC1*FMID.GE.0.) PAUSE 'ROOT MUST BE BRACKETED IN RTBIS'
     if(func1<0.)then
         rtbis = x1
         dx = x2 - x1
     else
         rtbis = x2
         dx = x1 - x2
     endif
     do j = 1, JMAX
         dx = dx*.5
         xmid = rtbis + dx
         call MANNINGS_EQN(Flow, xmid, fmid)
         if(fmid<=0.)rtbis = xmid
         if(ABS(dx)<xacc .OR. fmid==0.)then
             ELWS(i) = rtbis + EL(KB(i) + 1, i)                            ! SW 4/5/13
             return
         endif
     enddo
  !    PAUSE 'TOO MANY BISECTIONS IN RTBIS'
     end subroutine NORMAL_DEPTH
