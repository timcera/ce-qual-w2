!*==random.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************
!***********************************************************************
!*                                                                    **
!*           R A N D O M   N U M B E R   G E N E R A T O R            **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls No Other Subroutines
 
 
     subroutine RANDOM(Seed, Rndx)             ! Random Number Generator Subroutine
     implicit none                             !   obtained from "STRUCTURED FORTRAN
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real :: Rndx
     integer :: Seed
     intent (out) Rndx
     intent (inout) Seed
!
!*** End of declarations rewritten by SPAG
!
                                               !   77 for Engineers and Scientists, 5th
                                               !   Edition - Author: Delores M. Etter"
 
 
     Seed = 2045*Seed + 1                      !   The random # is between 0.0 and 1.0
     Seed = Seed - (Seed/1048576)*1048576
     Rndx = REAL(Seed + 1)/1048577.0
 
     end subroutine RANDOM
