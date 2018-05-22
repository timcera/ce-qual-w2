!*==gas_transfer.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
!************************************
 
!***********************************************************************************************************************************
!**                                          S U B R O U T I N E   G A S   T R A N S F E R                                        **
!***********************************************************************************************************************************
 
     subroutine GAS_TRANSFER
     use GLOBAL
     use GEOMC
     use KINETIC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
     real, parameter :: THETA_REAERATION = 1.024, M_TO_FT = 3.2808
!
! Local variables
!
     real :: a, adepth, area, bcoef, dmo2, hdepth, s, uavg, ustar
     integer :: k
!
!*** End of declarations rewritten by SPAG
!
 
     if(REAERC(jw)=='   RIVER')then
 
!**      Average depth in ft
 
         area = 0.0
         do k = kt, KBMIN(i)
             area = area + BHR1(k, i)
         enddo
         adepth = area/BR(KTI(i), i)*M_TO_FT
 
!**      Average velocity in feet/second
 
         uavg = ABS(QC(i))/area*M_TO_FT
 
!**      Reaeration factor
 
         if(NEQN(jw)==0)then
             if(adepth<=2.0)then
                 REAER(i) = 21.64*uavg**0.67/adepth**1.85
             elseif(uavg<=1.8)then
                 REAER(i) = 12.96*SQRT(uavg)/adepth**1.5
             else
                 hdepth = -11.875*uavg + 23.375
                 if(hdepth>=adepth)then
                     REAER(i) = 12.96*SQRT(uavg)/adepth**1.5
                 else
                     REAER(i) = 11.57*uavg**0.969/adepth**1.673
                 endif
             endif
         elseif(NEQN(jw)==1)then                                                                              !O'connor-Dobbins
             REAER(i) = 12.96*SQRT(uavg)/adepth**1.5
         elseif(NEQN(jw)==2)then                                                                              !Churchill
             REAER(i) = 11.57*uavg**0.969/adepth**1.673
         elseif(NEQN(jw)==3)then                                                                              !Tsivoglou
             s = SLOPEC(jb)*5280.0
             if(ABS(QC(i))*35.5>=10.0)then
                 REAER(i) = 0.88*s*uavg
             else
                 REAER(i) = 1.8*s*uavg
             endif
         elseif(NEQN(jw)==4)then                                                                              !Owens
             REAER(i) = 21.64*uavg**0.67/adepth**1.85
         elseif(NEQN(jw)==5)then                                                                              !Thackston and Krenkel
             ustar = SQRT(adepth*SLOPEC(jb)*32.2)                                                              ! SR 5/10/05
             REAER(i) = 24.88*(1.0 + SQRT(0.176*uavg/SQRT(adepth)))            &
                      & *ustar/adepth                                                                         ! SR 5/10/05
         elseif(NEQN(jw)==6)then                                                                              !Langbien and Durum
             REAER(i) = 7.60*uavg/adepth**1.33
         elseif(NEQN(jw)==7)then                                                                              !Melching and Flores
             uavg = uavg/M_TO_FT
             if(QC(i)==0.0)then
                 REAER(i) = 0.0
             elseif(ABS(QC(i))<0.556)then
                 REAER(i) = 517.0*((uavg*SLOPEC(jb))**0.524)*ABS(QC(i))        &
                          & **( - 0.242)
             else
                 REAER(i) = 596.0*((uavg*SLOPEC(jb))**0.528)*ABS(QC(i))        &
                          & **( - 0.136)
             endif
         elseif(NEQN(jw)==8)then                                                                              !Melching and Flores
             uavg = uavg/M_TO_FT
             adepth = adepth/M_TO_FT
             if(ABS(QC(i))<0.556)then
                 REAER(i) = 88.0*((uavg*SLOPEC(jb))**0.313)*adepth**( - 0.353)
             else
                 REAER(i) = 142.0*((uavg*SLOPEC(jb))**0.333)*adepth**( - 0.66) &
                          & *BI(kt, i)**( - 0.243)
             endif
         elseif(NEQN(jw)==9)then                                                                              !User defined SI units
             uavg = uavg/M_TO_FT
             adepth = adepth/M_TO_FT
             REAER(i) = RCOEF1(jw)*(uavg**RCOEF2(jw))*(adepth**RCOEF3(jw))     &
                      & *(SLOPEC(jb)**RCOEF4(jw))
         elseif(NEQN(jw)==10)then                                                                             ! Thackston and Krenkel - updated
             ustar = SQRT(adepth*SLOPEC(jb)*32.2)                                                              ! SR 5/10/05
             REAER(i) = 4.99*(1.0 + 9.0*(0.176*uavg/SQRT(adepth))**0.25)       &
                      & *ustar/adepth                                                                         ! SR 5/10/05
         endif
         REAER(i) = REAER(i)*adepth/M_TO_FT
     elseif(REAERC(jw)=='    LAKE')then
         if(NEQN(jw)==1)then                                                                                  !Broecker
             REAER(i) = 0.864*WIND10(i)
         elseif(NEQN(jw)==2)then
             if(WIND10(i)<=3.5)then                                                                           !Gelda
                 a = 0.2
                 bcoef = 1.0
             else
                 a = 0.057
                 bcoef = 2.0
             endif
             REAER(i) = a*WIND10(i)**bcoef
         elseif(NEQN(jw)==3)then                                                                              !Banks & Herrera
             REAER(i) = (0.728*SQRT(WIND10(i)) - 0.317*WIND10(i)               &
                      & + 0.0372*WIND10(i)**2)
         elseif(NEQN(jw)==4)then                                                                              !Wanninkhof
             REAER(i) = 0.0986*WIND10(i)**1.64
         elseif(NEQN(jw)==5)then                                                                              !Chen & Kanwisher
             dmo2 = 2.04E-9
             REAER(i) = day*dmo2/((200.0 - 60.0*SQRT(MIN(WIND10(i),11.0)))     &
                      & *1.E-6)
         elseif(NEQN(jw)==6)then                                                                              !Cole & Buchak
             REAER(i) = (0.5 + 0.05*WIND10(i)*WIND10(i))
         elseif(NEQN(jw)==7)then                                                                              !Banks
             if(WIND10(i)<=5.5)then
                 REAER(i) = 0.362*SQRT(WIND10(i))
             else
                 REAER(i) = 0.0277*WIND10(i)**2
             endif
         elseif(NEQN(jw)==8)then                                                                              !Smith
             REAER(i) = 0.64 + 0.128*WIND10(i)**2
         elseif(NEQN(jw)==9)then                                                                              !Liss
             if(WIND10(i)<=4.1)then
                 REAER(i) = 0.156*WIND10(i)**0.63
             else
                 REAER(i) = 0.0269*WIND10(i)**1.9
             endif
         elseif(NEQN(jw)==10)then                                                                             !Downing and Truesdale
             REAER(i) = 0.0276*WIND10(i)**2
         elseif(NEQN(jw)==11)then                                                                             !Kanwisher
             REAER(i) = 0.0432*WIND10(i)**2
         elseif(NEQN(jw)==12)then                                                                             !Yu, et al
             REAER(i) = 0.319*WIND10(i)
         elseif(NEQN(jw)==13)then                                                                             !Weiler
             if(WIND10(i)<=1.6)then
                 REAER(i) = 0.398
             else
                 REAER(i) = 0.155*WIND10(i)**2
             endif
         elseif(NEQN(jw)==14)then                                                                             !User defined
             REAER(i) = RCOEF1(jw) + RCOEF2(jw)*WIND10(i)**RCOEF3(jw)
         endif
     elseif(REAERC(jw)==' ESTUARY')then
         area = 0.0
         do k = kt, KBMIN(i)
             area = area + BHR1(k, i)
         enddo
 
         adepth = area/BR(KTI(i), i)
         uavg = ABS(QC(i))/area
 
!**      Reaeration factor
 
         if(NEQN(jw)==0)then
             adepth = adepth*M_TO_FT
             uavg = uavg*M_TO_FT
             if(adepth<=2.0)then
                 REAER(i) = 21.64*uavg**0.67/adepth**1.85
             elseif(uavg<=1.8)then
                 REAER(i) = 12.96*SQRT(uavg)/adepth**1.5
             else
                 hdepth = -11.875*uavg + 23.375
                 if(hdepth>=adepth)then
                     REAER(i) = 12.96*SQRT(uavg)/adepth**1.5
                 else
                     REAER(i) = 11.57*uavg**0.969/adepth**1.673
                 endif
             endif
         elseif(NEQN(jw)==1)then                                                                           !Thomann and Fitzpatrick
             REAER(i) = (0.728*SQRT(WIND10(i)) - 0.317*WIND10(i)               &
                      & + 0.0372*WIND10(i)**2) + 3.93*SQRT(uavg)/(adepth)**0.5
         elseif(NEQN(jw)==2)then                                                                           !User defined
             REAER(i) = RCOEF1(jw)*(uavg**RCOEF2(jw))*(adepth**RCOEF3(jw))     &
                      & + (0.5 + RCOEF4(jw)*WIND10(i)*WIND10(i))
         endif
 
  !  ADEPTH = AREA/BR(KTI(I),I)*M_TO_FT
  !  UAVG   = ABS(QC(I))/AREA*M_TO_FT
 
!**      Reaeration factor
 
   ! IF (NEQN(JW) == 0) THEN
   !   IF (ADEPTH <= 2.0) THEN
   !     REAER(I) = 21.64*UAVG**0.67/ADEPTH**1.85
   !   ELSE IF (UAVG <= 1.8) THEN
   !     REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
   !   ELSE
   !     HDEPTH = -11.875*UAVG+23.375
   !     IF (HDEPTH >= ADEPTH) THEN
   !       REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
   !     ELSE
   !       REAER(I) = 11.57*UAVG**0.969/ADEPTH**1.673
   !     END IF
   !   END IF
   ! ELSE IF (NEQN(JW) == 1) THEN                                                                           !Thomann and Fitzpatrick
   !   REAER(I) = (0.728*SQRT(WIND10(I))-0.317*WIND10(I)+0.0372*WIND10(I)**2)+3.93*SQRT(UAVG/M_TO_FT)/(ADEPTH/M_TO_FT)**0.5
   ! END IF
 
     endif
     if(REAER(i)<=0.6)REAER(i) = 0.6
     REAER(i) = REAER(i)*THETA_REAERATION**(T1(kt, i) - 20.0)
     REAER(i) = REAER(i)/day
     end subroutine GAS_TRANSFER
