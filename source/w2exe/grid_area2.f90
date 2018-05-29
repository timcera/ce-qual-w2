!*==grid_area2.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                             S U B R O U T I N E   G R I D  A R E A 2                                          **
!***********************************************************************************************************************************
 
     subroutine GRID_AREA2
     use F77KINDS
     use GLOBAL
     use GEOMC
     use RSTART
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND) :: area, a_coef, b_coef, c_coef, sl
     integer :: k
!
!*** End of declarations rewritten by SPAG
!
 
     area = (EL(kt, i) - SZ(i) - (EL(kt, i) - Z(i)))*BI(kt, i)
     sl = (B(kt, i) - BB(kt, i))/(0.5*H(kt, jw))
     a_coef = -1.0
     b_coef = SZ(i)*2. + BI(kt, i)/(0.5*sl)
     c_coef = -area/(0.5*sl) - SZ(i)**2 - BI(kt, i)*2.*SZ(i)/sl
     Z(i) = ( - b_coef + SQRT(b_coef**2 - 4.*a_coef*c_coef))/(2.0*a_coef)
     KTI(i) = 2
     do k = 2, KB(i)
         if(EL(k, i)<=EL(kt, i) - Z(i))then
             KTI(i) = k - 1
             exit
         endif
     enddo
     end subroutine GRID_AREA2
