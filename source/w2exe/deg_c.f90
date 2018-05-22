!*==deg_c.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     function DEG_C(X)
     use PREC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(r8) :: X
     real :: DEG_C
     intent (in) X
!
! Local variables
!
     real :: PREC
     real :: r8
!
!*** End of declarations rewritten by SPAG
!
     DEG_C = (X - 32.0)*5.0/9.0
     end function DEG_C
