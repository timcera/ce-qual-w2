!*==shading.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************************************************************************
!**                                                S U B R O U T I N E   S H A D I N G                                            **
!***********************************************************************************************************************************
 
     subroutine SHADING
     use SHADEC
     use GLOBAL
     use GDAYC
     use SURFHE
     use GEOMC
     use SCREENC
     use LOGICC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: a0, a02, ang1, ang2, ax, az00, azt, cline, edaz, edge, hour, ht,  &
           & local, sfact, sinal, sn, sred, standard, stlen, taud, topoang
     character(1) :: bank
     integer :: iday, j
!
!*** End of declarations rewritten by SPAG
!
 
!    Calculate solar altitude, declination, and local hour angle when
 
!    short-wave solar radiation is provided as input
     if(READ_RADIATION(jw))then
         local = LONGIT(jw)
         standard = 15.0*INT(LONGIT(jw)/15.0)
         hour = (jday - INT(jday))*24.0
         iday = jday - ((INT(jday/365))*365)
         iday = iday + INT(INT(jday/365)/4)
         taud = (2*pi*(iday - 1))/365
         eqtnew = 0.170*SIN(4*pi*(iday - 80)/373)                              &
                & - 0.129*SIN(2*pi*(iday - 8)/355)
         HH(jw) = 0.261799*(hour - (local - standard)*0.0666667 + eqtnew -     &
                & 12.0)
         DECL(jw) = 0.006918 - 0.399912*COS(taud) + 0.070257*SIN(taud)         &
                  & - 0.006758*COS(2*taud) + 0.000907*SIN(2*taud)              &
                  & - 0.002697*COS(3*taud) + 0.001480*SIN(3*taud)
         sinal = SIN(LAT(jw)*.0174533)*SIN(DECL(jw)) + COS(LAT(jw)*.0174533)   &
               & *COS(DECL(jw))*COS(HH(jw))
         A00(jw) = 57.2957795*ASIN(sinal)
     endif
 
!    If the sun is below the horizon, set SHADE(I) to 0
 
     if(A00(jw)<0.0)then
         SHADE(i) = 0.0
     else
 
!**      Calculate solar azimuth angle
 
         a02 = A00(jw)/57.2957795
         ax = (SIN(DECL(jw))*COS(LAT(jw)*0.017453) - COS(DECL(jw))*COS(HH(jw)) &
            & *SIN(LAT(jw)*0.017453))/COS(a02)
         if(ax>1.0)ax = 1.0
         if(ax< - 1.0)ax = -1.0
         azt = ACOS(ax)
         if(HH(jw)<0.0)then
             az00 = azt
         else
             az00 = 2.0*pi - azt
         endif
         a0 = a02
 
!**      Interpolate the topographic shade angle
 
         do j = 1, iang - 1
             if(az00>ANG(j) .AND. az00<=ANG(j + 1))then
                 ang1 = az00 - ANG(j)
                 ang2 = (TOPO(i, j + 1) - TOPO(i, j))/gama     ! SW 10/17/05
                 topoang = TOPO(i, j) + ang2*ang1
             endif
         enddo
         if(az00>ANG(iang) .AND. az00<=2*pi)then
             ang1 = az00 - ANG(iang)
             ang2 = (TOPO(i, 1) - TOPO(i, iang))/gama          ! SW 10/17/05
             topoang = TOPO(i, iang) + ang2*ang1
         endif
 
!**      Complete topographic shading if solar altitude less than topo angle
 
         if(a0<=topoang)then
             sfact = 0.90
             goto 50
         endif
 
!**      No vegetative shading if azimuth angle is oriented parallel to stream
 
         if(az00==PHI0(i) .OR. az00==PHI0(i) + pi .OR. az00 + pi==PHI0(i))then
             sfact = 0.0
             goto 50
         endif
 
!**      Bank with the controlling vegetation
 
         if(PHI0(i)>0.0 .AND. PHI0(i)<=pi)then
             if(az00>PHI0(i) .AND. az00<=PHI0(i) + pi)bank = 'L'
             if(az00>0.0 .AND. az00<=PHI0(i))bank = 'R'
             if(az00>PHI0(i) + pi .AND. az00<2.0*pi)bank = 'R'
         elseif(PHI0(i)>pi .AND. PHI0(i)<=2.0*pi)then
             if(az00>=PHI0(i) .AND. az00<2.0*pi)bank = 'L'
             if(az00>=0.0 .AND. az00<PHI0(i) - pi)bank = 'L'
             if(az00>=PHI0(i) - pi .AND. az00<PHI0(i))bank = 'R'
         endif
 
!**      No topographic shading
 
         if(bank=='L')then
             if(TTLB(i)<ELWS(i))then
                 sfact = 0.0
                 goto 50
             else
                 ht = TTLB(i) - ELWS(i)
                 cline = CLLB(i)
                 sred = SRLB2(i)
                 if(jdayg>SRFJD1(i) .AND. jdayg<=SRFJD2(i))sred = SRLB1(i)
             endif
         elseif(TTRB(i)<ELWS(i))then
             sfact = 0.0
             goto 50
         else
             ht = TTRB(i) - ELWS(i)
             cline = CLRB(i)
             sred = SRRB2(i)
             if(jdayg>SRFJD1(i) .AND. jdayg<=SRFJD2(i))sred = SRRB1(i)
         endif
         stlen = ht/TAN(a0)
         edge = MAX(0.0, cline - BI(kt, i)/2.0)
 
!**      Distance from vegetation to water edge on line parallel to azimuth
 
         edaz = edge/ABS(SIN(PHI0(i) - az00))
         if(stlen<=edaz)then
             sfact = 0.0
             goto 50
         endif
 
!**      Distance shadow extends over water (perpendicular to segment
 
!        orientation)
         sn = MIN(ht*ABS(SIN(ABS(PHI0(i)-az00)))/TAN(a0) - edge, BI(kt, i))
         sfact = sred*sn/BI(kt, i)
50       SHADE(i) = MAX(0.0, 1 - sfact)
         SHADE(i) = MIN(ABS(SHADEI(i)), SHADE(i))        ! SW 10/2/2017 Allows for fixed canopy cover over top of channel - only used if shade is less than shadei only valid for -0.99 and 0.0
     endif
     end subroutine SHADING
