!*==cdfunc.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                                  F U N C T I O N   C D F U N C                                                **
!***********************************************************************************************************************************
 
     function CDFUNC(Depth, Flow)
     use F77KINDS
     use STRUCTURES
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(R8KIND) :: Depth, Flow
     real :: CDFUNC
     intent (in) Flow
!
! Local variables
!
     real(R8KIND) :: BAREA, TWIDTH
!
!*** End of declarations rewritten by SPAG
!
     CDFUNC = (Flow**2*TWIDTH(Depth, dia))/(BAREA(Depth, dia)**3*9.81D0)       &
            & - 1.0D0
     end function CDFUNC
