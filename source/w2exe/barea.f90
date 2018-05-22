!*==barea.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                                  F U N C T I O N   B A R E A                                                  **
!***********************************************************************************************************************************
 
     function BAREA(Depth, Dia)
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
     real :: BAREA
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
         BAREA = (Depth - Dia*0.5D0)*DSQRT(Depth*Dia - Depth**2)               &
               & + (Dia**2*0.25D0)*DASIN((2.0D0/Dia)*(Depth - Dia*0.5D0))      &
               & + (PI*Dia**2)/8.0D0
     else
         BAREA = (PI*Dia**2)*0.25D0
     endif
     end function BAREA
