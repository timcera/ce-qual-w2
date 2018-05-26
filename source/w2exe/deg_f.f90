!*==deg_f.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
! Function declaration
 
     function DEG_F(X)
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(R8KIND) :: X
     real :: DEG_F
     intent (in) X
!
! Local variables
!
!
!*** End of declarations rewritten by SPAG
!
     DEG_F = X*1.8 + 32.0
     end function DEG_F
