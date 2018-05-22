!*==twidth.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                                  F U N C T I O N   T W I D T H                                                **
!***********************************************************************************************************************************
 
     function TWIDTH(Depth, Dia)
     use PREC
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(r8) :: Depth, Dia
     real :: TWIDTH
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
         TWIDTH = 2.0D0*DSQRT((Dia*Depth) - Depth**2)
     else
         TWIDTH = 0.005D0*Dia
     endif
     end function TWIDTH
