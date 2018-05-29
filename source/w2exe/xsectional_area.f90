!*==xsectional_area.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**        S U B R O U T I N E    C R O S S    S E C T I O N A L    A R E A                                                       **
!***********************************************************************************************************************************
 
     subroutine XSECTIONAL_AREA(Wsurf, Xarea)
     use GLOBAL
     use GEOMC
     use MAIN
     use INITIALVELOCITY
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real :: Wsurf, Xarea
     intent (in) Wsurf
     intent (inout) Xarea
!
! Local variables
!
     integer :: kttop
!
!*** End of declarations rewritten by SPAG
!
 
 
 
 
 
                               ! 4/5/13 SW
 
 !     KTTOP = 2
     do k = 2, kmx - 1         ! 4/5/13 SW
         if(EL(k, i)<Wsurf)then
             kttop = k - 1
             exit
         endif
         kttop = k  ! CB 8/10/10
     enddo
 !     DO WHILE (EL(KTTOP,I) > WSURF)
 !        KTTOP = KTTOP+1
 !     END DO
     Xarea = (Wsurf - EL(kttop + 1, i))*BSAVE(kttop, i)
     do k = kttop + 1, KBI(i)
         Xarea = Xarea + BSAVE(k, i)*H(k, jw)
     enddo
 
     end subroutine XSECTIONAL_AREA
