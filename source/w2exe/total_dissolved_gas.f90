!*==total_dissolved_gas.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                       S U B R O U T I N E   T O T A L  D I S S O L V E D  G A S                               **
!***********************************************************************************************************************************
 
     subroutine TOTAL_DISSOLVED_GAS(Nsat, P, Nsg, N, T, C)
     use TDGAS
     use STRUCTURES
     use GLOBAL
     use MAIN, ONLY:ea
     use TVDC, ONLY:TDEW
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(R8KIND) :: C, P, T
     integer :: N, Nsat, Nsg
     intent (in) Nsat, Nsg, P, T
     intent (inout) C
!
! Local variables
!
     real :: da, db, sat, tdg
!
!*** End of declarations rewritten by SPAG
!
 
     if(Nsat==0)then ! DISSOLVED OXYGEN
         sat = EXP(7.7117 - 1.31403*(LOG(T + 45.93)))*P
     else            ! N2 GAS
         ea = DEXP(2.3026D0*(7.5D0*TDEW(jw)/(TDEW(jw) + 237.3D0) + 0.6609D0))  &
            & *0.001316                                                       ! in mm Hg   0.0098692atm=7.5006151mmHg
  !SAT=(1.5568D06*0.79*(P-EA)*(1.8816D-5 - 4.116D-7 * T + 4.6D-9 * T**2))   ! SW 10/27/15
         sat = (1.5568D06*0.79*(P - ea)*(1.8816D-5 - 4.116D-7*T + 4.6D-9*T*T))
                                                                           ! SW 10/27/15    4/20/16 SPEED
     endif
 
     if(Nsg==0)then
         if(EQSP(N)==1)then
             tdg = AGASSP(N)*.035313*QSP(N) + BGASSP(N)
             if(tdg>145.0)tdg = 145.0
             C = sat
             if(tdg>=100.0)C = tdg*sat/100.0
         elseif(EQSP(N)==2)then
             tdg = AGASSP(N) + BGASSP(N)*EXP(0.035313*QSP(N)*CGASSP(N))
             if(tdg>145.0)tdg = 145.0
             C = sat
             if(tdg>=100.0)C = tdg*sat/100.0
         else
             da = sat - C                                                                      ! MM 5/21/2009 DA: Deficit upstream
             db = da/(1.0 + 0.38*AGASSP(N)*BGASSP(N)*CGASSP(N)                 &
                & *(1.0 - 0.11*CGASSP(N))*(1.0 + 0.046*T))                                     ! DB: deficit downstream
             C = sat - db
         endif
     elseif(EQGT(N)==1)then
         tdg = AGASGT(N)*0.035313*QGT(N) + BGASGT(N)
         if(tdg>145.0)tdg = 145.0
         C = sat
         if(tdg>=100.0)C = tdg*sat/100.0
     elseif(EQGT(N)==2)then
         tdg = AGASGT(N) + BGASGT(N)*EXP(.035313*QGT(N)*CGASGT(N))
         if(tdg>145.0)tdg = 145.0
         C = sat
         if(tdg>=100.0)C = tdg*sat/100.0
     else
         da = sat - C                                                                        ! MM 5/21/2009 DA: Deficit upstream
         db = da/(1.0 + 0.38*AGASGT(N)*BGASGT(N)*CGASGT(N)                     &
            & *(1.0 - 0.11*CGASGT(N))*(1.0 + 0.046*T))                                       ! DB: deficit downstream
         C = sat - db
     endif
     end subroutine TOTAL_DISSOLVED_GAS
