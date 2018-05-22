!*==density.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
!***********************************************************************************************************************************
!**                                            F U N C T I O N   D E N S I T Y                                                    **
!***********************************************************************************************************************************
 
     function DENSITY(T, Tds, Ss)
     use PREC
     use LOGICC, ONLY:susp_solids, FRESH_WATER, SALT_WATER
     use GLOBAL, ONLY:jw
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(r8) :: Ss, T, Tds
     real :: DENSITY
     intent (in) Ss, T, Tds
!
!*** End of declarations rewritten by SPAG
!
     DENSITY = ((((6.536332D-9*T-1.120083D-6)*T + 1.001685D-4)*T - 9.09529D-3) &
             & *T + 6.793952D-2)*T + 0.842594D0
     if(susp_solids)DENSITY = DENSITY + 6.2D-4*Ss
     if(FRESH_WATER(jw))DENSITY = DENSITY +                                    &
                                & Tds*((4.99D-8*T - 3.87D-6)*T + 8.221D-4)
     if(SALT_WATER(jw))DENSITY = DENSITY +                                     &
                               & Tds*((((5.3875D-9*T-8.2467D-7)*T + 7.6438D-5) &
                               & *T - 4.0899D-3)*T + 0.824493D0)               &
                               & + (( - 1.6546D-6*T + 1.0227D-4)               &
                               & *T - 5.72466D-3)*Tds**1.5D0 +                 &
                               & 4.8314D-4*Tds*Tds
     DENSITY = DENSITY + 999.0D0
     end function DENSITY
