!*==wallfunction.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     function WALLFUNCTION()
     use GLOBAL
     use GEOMC
     use EDDY
     use LOGICC
     use MAIN, ONLY:warning_open
     use SCREENC, ONLY:jday
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real :: WALLFUNCTION
!
! Local variables
!
     real(r8) :: area, bs1, ks, lbound, perimeter, ubound, ust, v, visck
     integer :: k, multiplier
     real(r8), external :: RTBIS, SEMILOG
!
!*** End of declarations rewritten by SPAG
!
                                                                                                          !,USTARB
     visck = DEXP( - (T2(KB(i), i) + 495.691)*0.026746764)                                               !/(-37.3877)
     if(T2(KB(i), i)>30.0)visck = DEXP( - (T2(KB(i), i) + 782.190)*            &
                                & 0.01731301939)                                                         !/(-57.7600)
     ubound = ABS((U(KB(i), i) + U(KB(i), i - 1))*.5)
     if(SEMILOG(ubound, ABS((U(KB(i),i)+U(KB(i),i-1))*.5), (H(KB(i),jw))/2.0,  &
      & visck)>0)then                                                                                  !FRIC(I),
         multiplier = 2
         do while (SEMILOG(ubound, ABS((U(KB(i),i)+U(KB(i),i-1))*.5),          &
                 & (H(KB(i),jw))/2.0, visck)>0)                                                            !,FRIC(I)
             ubound = ABS((U(KB(i), i) + U(KB(i), i - 1))*.5)*REAL(multiplier)
             multiplier = multiplier + 1
             if(multiplier>30)then
                 warning_open = .TRUE.
                 write(wrn, *)                                                 &
                &'WALLFUNCTION: UPPER BOUND NOT FOUND IN TKE1 ROUTINE ON JDAY:'&
               & , jday
                 write(wrn, *)'Setting TKE calculation to TKE from TKE1'
                 AZC(jw) = '     TKE'
                 WALLFUNCTION = 0.0
                 return
             endif
         enddo
     endif
     lbound = visck/(E(i)*(H(KB(i), jw))/2.0D0) + nonzero
     multiplier = 1
     do while (SEMILOG(lbound, DABS((U(KB(i),i)+U(KB(i),i-1))*.5D0),           &
             & (H(KB(i),jw))/2.0D0, visck)<0)
         lbound = visck/(E(i)*(H(KB(i), jw))/2.0) + 10**multiplier*nonzero
         multiplier = multiplier + 1
         if(multiplier>30)then
             warning_open = .TRUE.
             write(wrn, *)                                                     &
                &'WALLFUNCTION: LOWER BOUND NOT FOUND IN TKE1 ROUTINE ON JDAY:'&
               & , jday
             write(wrn, *)'Setting TKE calculation to TKE from TKE1'
             AZC(jw) = '     TKE'
             WALLFUNCTION = 0.0
             return
         endif
     enddo
     WALLFUNCTION = RTBIS(lbound, ubound, 1D-15,                               &
                  & DABS((U(KB(i),i) + U(KB(i),i-1))*.5D0), (H(KB(i), jw))     &
                  & /2.0D0, visck)                                                                                     !,FRIC(I)
     if(STRICKON(jw))then
         v = visck
         ust = WALLFUNCTION
         if(MANNINGS_N(jw))then
             ks = ((STRICK(jw)*FRIC(i))**6)
         else
             perimeter = B(KB(i), i) + 2.0D0*H(KB(i), jw)
             do k = kt, KB(i) - 1
                 perimeter = perimeter +                                       &
                           & (B(k, i) - B(k + 1, i) + 2.0D0*H(k, jw))
             enddo
             area = 0
             do k = kt, KB(i)
                 area = area + BH1(k, i)
             enddo
             ks = ((STRICK(jw)*1.0D0/FRIC(i)*(area/perimeter)**(0.16666666667D0&
                & ))**6)                                                                    !1.0/6.0
         endif
         if(ust*ks/v>1.0)then
             bs1 = (5.5D0 + 2.5D0*DLOG(ust*ks/v))                              &
                 & *DEXP( - 0.217D0*DLOG(ust*ks/v)**2)                         &
                 & + 8.5D0*(1 - DEXP( - 0.217D0*DLOG(ust*ks/v)**2))
             E(i) = DEXP((.41*bs1)/(ust*ks/v))
         else
             E(i) = 9.535D0
         endif
     endif
     end function WALLFUNCTION
