!*==mannings_eqn.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************************************************************************
!**        S U B R O U T I N E    M A N N I N G S    E Q U A T I O N                                                              **
!***********************************************************************************************************************************
 
     subroutine MANNINGS_EQN(Flow, Depth, Funcvalue)
     use GLOBAL
     use GEOMC
     use EDDY
     use LOGICC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real :: Depth, Funcvalue
     real(r8) :: Flow
     intent (in) Depth, Flow
     intent (out) Funcvalue
!
! Local variables
!
     real :: B, EL, FRIC, SLOPEC
     real :: EDDY, GEOMC, GLOBAL
     real :: fmann, hrad, r8, wper, wsurf, xarea
     integer :: i, jb, jw
     integer :: KB, KTI, MANNINGS_N
     integer :: LOGICC
!
!*** End of declarations rewritten by SPAG
!
 
 
 
 
 
 
 !     WSURF=EL(KB(I)-1,I)+DEPTH
     wsurf = EL(KB(i) + 1, i) + Depth
                                  ! CB 7/7/10
     call XSECTIONAL_AREA(wsurf, xarea)
     wper = B(KTI(i), i) + 2.0*Depth
     hrad = xarea/wper
     if(MANNINGS_N(jw))then
         fmann = FRIC(i)
     else
         fmann = hrad**0.166666667/FRIC(i)
     endif
     Funcvalue = Flow - xarea*hrad**0.6667*SLOPEC(jb)**0.5/fmann         ! SW 4/5/2013
 
     end subroutine MANNINGS_EQN
