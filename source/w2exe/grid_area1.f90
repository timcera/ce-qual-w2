!*==grid_area1.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                             S U B R O U T I N E   G R I D  A R E A 1                                          **
!***********************************************************************************************************************************
 
     subroutine GRID_AREA1(El1, El2, Diff, Btop)
     use GLOBAL
     use GEOMC
     use PREC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(r8) :: Btop, Diff, El1, El2
     intent (in) El1, El2
     intent (out) Btop, Diff
!
! Local variables
!
     real :: B, BB, BH, EL, H
     real(r8) :: barea1, barea2, dist, dist1, dist2, slpe
     real :: GEOMC, GLOBAL, PREC
     integer :: i, jw, k, k1, k2
     integer :: KB
     real :: r8
!
!*** End of declarations rewritten by SPAG
!
 
 
!    Difference in areas for trapezoidal geometry
 
     do k = 2, KB(i)
         if(EL(k, i)<=El1)then
             k1 = k
             exit
         endif
     enddo
     do k = 2, KB(i)
         if(EL(k, i)<=El2)then
             k2 = k
             exit
         endif
     enddo
     barea1 = 0.0
     barea2 = 0.0
     do k = KB(i), k1, -1
         barea1 = barea1 + BH(k, i)
     enddo
     dist = El1 - EL(k1, i)
     if(H(k1 - 1, jw)/2.0<dist)then
         dist1 = H(k1 - 1, jw)*0.5
         slpe = (B(k1 - 1, i) - BB(k1 - 1, i))/(0.5*H(k1 - 1, jw))
         barea1 = barea1 + BB(k1 - 1, i)*dist1 + 0.5*slpe*dist1*dist1
         dist2 = dist - H(k1 - 1, jw)*0.5
         slpe = (BB(k1 - 2, i) - B(k1 - 1, i))/(0.5*H(k1 - 1, jw))
         barea1 = barea1 + B(k1 - 1, i)*dist2 + 0.5*slpe*dist2*dist2
         Btop = B(k1 - 1, i) + dist2*slpe
     else
         slpe = (B(k1 - 1, i) - BB(k1 - 1, i))/(0.5*H(k1 - 1, jw))
         barea1 = barea1 + BB(k1 - 1, i)*dist + 0.5*slpe*dist*dist
         Btop = BB(k1 - 1, i) + dist*slpe
     endif
     do k = KB(i), k2, -1
         barea2 = barea2 + BH(k, i)
     enddo
     dist = El2 - EL(k2, i)
     if(H(k2 - 1, jw)/2.<dist)then
         dist1 = H(k2 - 1, jw)*0.5
         slpe = (B(k2 - 1, i) - BB(k2 - 1, i))/(0.5*H(k2 - 1, jw))
         barea2 = barea2 + BB(k2 - 1, i)*dist1 + 0.5*slpe*dist1*dist1
         dist2 = dist - H(k2 - 1, jw)*0.5
         slpe = (BB(k2 - 2, i) - B(k2 - 1, i))/(0.5*H(k2 - 1, jw))
         barea2 = barea2 + B(k2 - 1, i)*dist2 + 0.5*slpe*dist2*dist2
     else
         slpe = (B(k2 - 1, i) - BB(k2 - 1, i))/(0.5*H(k2 - 1, jw))
         barea2 = barea2 + BB(k2 - 1, i)*dist + 0.5*slpe*dist*dist
     endif
     Diff = barea1 - barea2
     end subroutine GRID_AREA1
