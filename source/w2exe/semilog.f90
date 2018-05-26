!*==semilog.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     function SEMILOG(Ust, Ures, Y, V)
     use EDDY
     use GLOBAL
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(R8KIND) :: Ures, Ust, V, Y
     real(R8KIND) :: SEMILOG
     intent (in) Ures, Ust, V, Y
!
!*** End of declarations rewritten by SPAG
!
     SEMILOG = 0.41*Ures/DLOG(E(i)*Y*Ust/V) - Ust
     end function SEMILOG
