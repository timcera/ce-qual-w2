!*==calcfric.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine CALCFRIC()
     use GLOBAL
     use GEOMC
     use EDDY
     use KINETIC
     use LOGICC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(r8) :: area, depth, depthl, depthr, perimeter, uavg
     integer :: k
!
!*** End of declarations rewritten by SPAG
!
     depthl = (ELWS(i) - EL(KB(i), i) + H2(KB(i), i)*COSA(jb))/COSA(jb)                                           !EL(KT,I)  -Z(I)  *COSA(JB)
     depthr = (ELWS(i + 1) - EL(KB(i + 1), i + 1) + H2(KB(i + 1), i)*COSA(jb)) &
            & /COSA(jb)                                                                                                    !EL(KT,I+1)-Z(I+1)*COSA(JB)
     depth = (depthr + depthl)*0.5
     perimeter = 0.0
     perimeter = BR(KB(i), i) + 2.0*H(KB(i), jw)
     do k = kt, KB(i) - 1
         perimeter = perimeter + (B(k, i) - B(k + 1, i) + 2.0*H(k, jw))
     enddo
     area = 0.0
     do k = kt, KB(i)
         area = area + BH1(k, i)
     enddo
     uavg = 0.5*(QC(i) + QC(i - 1))/area
     if(ABS(uavg)>nonzero)then
         if(MANNINGS_N(jw))then
             FRIC(i) = (area/perimeter)**(0.16666666667)                       &
                     & /SQRT(9.81*(uavg/USTARBTKE(i))**2)                                                    !1.0/6.0   (9.81/(USTARBTKE(I)/UAVG)**2)
         else
             FRIC(i) = SQRT(9.81*(uavg/USTARBTKE(i))**2)
         endif
     endif
     end subroutine CALCFRIC
