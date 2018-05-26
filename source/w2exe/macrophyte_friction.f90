!*==macrophyte_friction.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!************************************************************************
!**               S U B R O U T I N E    MACROPHYTE_FRICTION           **
!************************************************************************
 
     subroutine MACROPHYTE_FRICTION(Hrad, Bedfr, Effric, K, Ii)
 
 
 
 
     use GEOMC
     use GLOBAL
     use MACROPHYTEC
     use POROSITYC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(R8KIND) :: Bedfr, Effric, Hrad
     integer :: Ii, K
     intent (in) Bedfr, Hrad
     intent (out) Effric
!
! Local variables
!
     real :: artot, cdavg, frin, savolrat, sctot, tsarea, xsarea
     integer :: m
!
!*** End of declarations rewritten by SPAG
!
 
     do m = 1, nmc
         savolrat = DWV(m)/DWSA(m)
                                !CB 6/29/2006
         if(K==kt)then
!            SAREA(M)=VSTEMKT(II,M)*SAVOLRAT/PI
             SAREA(m) = VSTEMKT(Ii, m)*savolrat*ANORM(m)
                                                   !CB 6/29/2006
         else
!            SAREA(M)=VSTEM(K,II,M)*SAVOLRAT/PI
             SAREA(m) = VSTEM(K, Ii, m)*savolrat*ANORM(m)
                                                   !CB 6/29/2006
         endif
     enddo
     xsarea = BH2(K, Ii)
 
     tsarea = 0.0
     artot = 0.0
     sctot = 0.0
     do m = 1, nmc
         artot = artot + SAREA(m)
         sctot = sctot + CDDRAG(m)*SAREA(m)
         tsarea = tsarea + SAREA(m)
     enddo
 
     if(artot>0.0)then
         cdavg = sctot/artot
         frin = cdavg*tsarea*Hrad**(4./3.)/(2.0*g*xsarea*DLX(Ii)*Bedfr**2)
         Effric = Bedfr*SQRT(1.0 + frin)
     else
         Effric = Bedfr
     endif
 
     end subroutine MACROPHYTE_FRICTION
