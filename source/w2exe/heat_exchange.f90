!*==heat_exchange.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                        S U B R O U T I N E   H E A T  E X C H A N G E                                         **
!***********************************************************************************************************************************
 
     subroutine HEAT_EXCHANGE
     use GLOBAL
     use GDAYC
     use SURFHE
     use TVDC
     use SHADEC
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real :: Jday
     real(R8KIND) :: Tsur
     intent (in) Jday, Tsur
!
! Local variables
!
     real(R8KIND) :: a0, aconv, bconv, beta, dtv, dtvl, ea, es, etp, fw, hour, ra, &
               & sinal, sro_br, standard, tairv, tair_f, taud, tdew_f, tstar,  &
               & wind2m, wind_mph
     real(R8KIND), save :: bowen_constant, btu_ft2_day_to_w_m2, flux_br_to_flux_si,&
                     & mps_to_mph, w_m2_to_btu_ft2_day
     real(R8KIND) :: DEG_C, DEG_F
     integer :: iday, j
     real :: local
!
!*** End of declarations rewritten by SPAG
!
 
!    Type declaration
 
!    Data declaration
 
     data mps_to_mph/2.23714D0/, w_m2_to_btu_ft2_day/7.60796D0/,               &
        & flux_br_to_flux_si/0.23659D0/
     data btu_ft2_day_to_w_m2/0.1314D0/
     data bowen_constant/0.47D0/
 
     return
 
!***********************************************************************************************************************************
!**  S H O R T  W A V E  R A D I A T I O N                                    
!***********************************************************************************************************************************
 
!    **
     entry SHORT_WAVE_RADIATION(Jday)
     local = LONGIT(jw)
     standard = 15.0*INT(LONGIT(jw)/15.0)
     hour = (Jday - INT(Jday))*24.0
     iday = Jday - ((INT(Jday/365))*365)
     iday = iday + INT(INT(Jday/365)/4)
     taud = (2.*pi*(iday - 1))/365.
     eqtnew = 0.170*SIN(4.*pi*(iday - 80)/373.)                                &
            & - 0.129*SIN(2.*pi*(iday - 8)/355.)
     HH(jw) = 0.261799*(hour - (local - standard)*0.0666667 + eqtnew - 12.0)
     DECL(jw) = 0.006918 - 0.399912*COS(taud) + 0.070257*SIN(taud)             &
              & - 0.006758*COS(2.*taud) + 0.000907*SIN(2.*taud)                &
              & - 0.002697*COS(3.*taud) + 0.001480*SIN(3.*taud)
     sinal = SIN(LAT(jw)*.0174533)*SIN(DECL(jw)) + COS(LAT(jw)*.0174533)       &
           & *COS(DECL(jw))*COS(HH(jw))
     A00(jw) = 57.2957795*ASIN(sinal)
     a0 = A00(jw)
     if(a0>0.0)then
         SRON(jw) = (1.0 - 0.0065*CLOUD(jw)**2)                                &
                  & *24.0*(2.044*a0 + 0.1296*a0**2 - 1.941E-3*a0**3 +          &
                  & 7.591E-6*a0**4)*btu_ft2_day_to_w_m2
     else
         SRON(jw) = 0.0
     endif
     return
 
!***********************************************************************************************************************************
!**  E Q U I L I B R I U M  T E M P E R A T U R E                             
!***********************************************************************************************************************************
 
!    **
     entry EQUILIBRIUM_TEMPERATURE
 
!    British units
 
     tdew_f = DEG_F(TDEW(jw))
     tair_f = DEG_F(TAIR(jw))
     sro_br = SRON(jw)*w_m2_to_btu_ft2_day*SHADE(i)
     wind_mph = WIND(jw)*WSC(i)*mps_to_mph
     wind2m = wind_mph*DLOG(2.0D0/Z0(jw))/DLOG(WINDH(jw)/Z0(jw)) + nonzero  ! SW 11/28/07  old version z0=0.003
     aconv = w_m2_to_btu_ft2_day
     if(CFW(jw)==1.0)bconv = 3.401062
     if(CFW(jw)==2.0)bconv = 1.520411
 
!    Equilibrium temperature and heat exchange coefficient
 
     ET(i) = tdew_f
     tstar = (ET(i) + tdew_f)*0.5
     beta = 0.255 - (8.5E-3*tstar) + (2.04E-4*tstar*tstar)
     fw = aconv*AFW(jw) + bconv*BFW(jw)*wind2m**CFW(jw)
     CSHE(i) = 15.7 + (0.26 + beta)*fw
     ra = 3.1872E-08*(tair_f + 459.67)**4
     etp = (sro_br + ra - 1801.0)/CSHE(i) + (CSHE(i) - 15.7)                   &
         & *(0.26*tair_f + beta*tdew_f)/(CSHE(i)*(0.26 + beta))
     j = 0
     do while (ABS(etp - ET(i))>0.05 .AND. j<10)
         ET(i) = etp
         tstar = (ET(i) + tdew_f)*0.5
         beta = 0.255 - (8.5E-3*tstar) + (2.04E-4*tstar*tstar)
         CSHE(i) = 15.7 + (0.26 + beta)*fw
         etp = (sro_br + ra - 1801.0)/CSHE(i) + (CSHE(i) - 15.7)               &
             & *(0.26*tair_f + beta*tdew_f)/(CSHE(i)*(0.26 + beta))
         j = j + 1
     enddo
 
!    SI units
 
     ET(i) = DEG_C(ET(i))
     CSHE(i) = CSHE(i)*flux_br_to_flux_si/rhowcp
     return
 
!***********************************************************************************************************************************
!**  S U R F A C E   T E R M S                                                
!***********************************************************************************************************************************
 
!    **
     entry SURFACE_TERMS(Tsur)
 
!    Partial water vapor pressure of air (mm hg)
 
!    EA = EXP(2.3026*(9.5*TDEW(JW)/(TDEW(JW)+265.5)+0.6609))         ! SW
!    6/10/2011 IF (TDEW(JW) > 0.0) EA =
!    EXP(2.3026*(7.5*TDEW(JW)/(TDEW(JW)+237.3)+0.6609))
     ea = DEXP(2.3026D0*(7.5D0*TDEW(jw)/(TDEW(jw) + 237.3D0) + 0.6609D0))
 
!    Partial water vapor pressure at the water surface
 
     if(Tsur<0.0)then
         es = DEXP(2.3026D0*(9.5D0*Tsur/(Tsur + 265.5D0) + 0.6609D0))
     else
         es = DEXP(2.3026D0*(7.5D0*Tsur/(Tsur + 237.3D0) + 0.6609D0))
     endif
 
!    Wind function
 
     if(RH_EVAP(jw))then
         tairv = (TAIR(jw) + 273.0D0)/(1.0D0 - 0.378D0*ea/760.0D0)
         dtv = (Tsur + 273.0D0)/(1.0D0 - 0.378D0*es/760.0D0) - tairv
         dtvl = 0.0084D0*WIND2(i)**3
         if(dtv<dtvl)dtv = dtvl
         fw = (3.59D0*dtv**0.3333D0 + 4.26D0*WIND2(i))
     else
         fw = AFW(jw) + BFW(jw)*WIND2(i)**CFW(jw)
     endif
 
!    Evaporative flux
 
     RE(i) = fw*(es - ea)
  !IF(RE(I) < 0.0)RE(I)=0.0     ! SW 6/22/2016  SHOULD WE USE THIS AS AN ANALOG FOR CONDENSATION WHEN LESS THAN 0? TVA(1972) SUGGESTS YOU CAN - SEE SECTION 4 AND 5 but Ryan, Harleman suggest setting RE=0 if less than zero since the condensation coefficients are unknown...
 
!    Conductive flux
 
     RC(i) = fw*bowen_constant*(Tsur - TAIR(jw))
 
!    Back radiation flux
 
     RB(i) = 5.51D-8*(Tsur + 273.15D0)**4
     end subroutine HEAT_EXCHANGE
