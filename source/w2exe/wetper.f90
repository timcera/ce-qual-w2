!*==wetper.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                                  F U N C T I O N   W E T P E R                                                **
!***********************************************************************************************************************************
 
     function WETPER(Depth, Dia)
     use PREC
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
     real(r8), parameter :: PI = 3.14159265359D0
!
! Dummy arguments
!
     real(r8) :: Depth, Dia
     real :: WETPER
     intent (in) Depth, Dia
!
! Local variables
!
     real :: PREC
     real :: r8
!
!*** End of declarations rewritten by SPAG
!
     if(Depth<Dia)then
         WETPER = Dia*(DASIN((2.0D0/Dia)*(Depth - Dia*0.5D0)) + PI*0.5D0)
     else
         WETPER = PI*Dia
     endif
     end function WETPER
