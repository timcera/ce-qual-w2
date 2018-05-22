!*==pipe_flow_initialize.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                            S U B R O U T I N E   P I P E  F L O W                                             **
!***********************************************************************************************************************************
 
     subroutine PIPE_FLOW_INITIALIZE
     use GLOBAL
     use GEOMC
     use STRUCTURES
     use SCREENC, ONLY:nit
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real, save :: begin, clen, closs, dia, dlt, dnie, eps2, fman, omega, qold,&
                 & r8, upie, vmax, wlflag
     real(r8), save :: d1, d2, dcheck, dcrit, dltx, dtest, dtq, el1, el2, eps, &
                     & hie, tott, vtot
     real(r8) :: DEPTHCRIT
     real :: DLX, DLXPI, EDPI, ELWS, EUPI, FMINPI, FPI, US, WPI
     real :: GEOMC, GLOBAL, SCREENC, STRUCTURES
     integer :: IDPI, IUPI, JBDPI, LATERAL_PIPE
     integer, save :: jp, nc, nit, npi
     real :: QPI
!
!*** End of declarations rewritten by SPAG
!
 
 
     allocate(begin(npi), wlflag(npi), vmax(npi))
                     !,NIT
     qold = 0.01
     vmax = 0.01
     begin = .TRUE.
     wlflag = .TRUE.
     return
 
     entry PIPE_FLOW !(NIT)
     dtq = dlt/10.0
     do jp = 1, npi
         dia = WPI(jp)
         clen = DLXPI(jp)
         fman = FPI(jp)
         closs = FMINPI(jp)
         upie = EUPI(jp)
         dnie = EDPI(jp)
         dltx = clen/(REAL(nc - 1)*0.5)
         if(LATERAL_PIPE(jp))then
             el1 = ELWS(IUPI(jp))                                                            !EL(KTWB(JWUPI(JP)),IUPI(JP))-Z(IUPI(JP))*COSA(JBUPI(JP))
         else
     ! EL1   = ELWS(IUPI(JP))-SINA(JBDPI(JP))*DLX(IUPI(JP))*0.5                               !EL(KTWB(JWUPI(JP)),IUPI(JP))-Z(IUPI(JP))*COSA(JBUPI(JP))-SINA(JBDPI(JP))*DLX(IUPI(JP))*0.5
             el1 = ELWS(IUPI(jp)) + (ELWS(IUPI(jp)) - ELWS(IUPI(jp) - 1))      &
                 & /(0.5*(DLX(IUPI(jp)) + DLX(IUPI(jp) - 1)))*DLX(IUPI(jp))*0.5                                             ! LINEAR INTERPOLATION OF THE WATER LEVEL TO THE EDGE
         endif
         if(IDPI(jp)==0)then
             el2 = -1.0
         elseif(US(JBDPI(jp))/=IDPI(jp))then
             el2 = ELWS(IDPI(jp))                                                            !EL(KTWB(JWDPI(JP)),IDPI(JP))-Z(IDPI(JP))*COSA(JBDPI(JP))
         else
       ! EL2   = ELWS(IDPI(JP))+SINA(JBDPI(JP))*DLX(IDPI(JP))*0.5                             !EL(KTWB(JWDPI(JP)),IDPI(JP))-Z(IDPI(JP))*COSA(JBDPI(JP))+SINA(JBDPI(JP))*DLX(IDPI(JP))*0.5
             el2 = ELWS(IDPI(jp)) - (ELWS(IDPI(jp) + 1) - ELWS(IDPI(jp)))      &
                 & /(0.5*(DLX(IDPI(jp)) + DLX(IDPI(jp) + 1)))*DLX(IDPI(jp))*0.5
         endif
         hie = MAX(upie, dnie)
         if(dia==0.0)then
             QPI(jp) = 0.0
             wlflag(jp) = .TRUE.
             goto 100
         endif
         eps = 0.001
         if((hie + eps)>=el1 .AND. (hie + eps)>=el2)then
             QPI(jp) = 0.0
             wlflag(jp) = .TRUE.
             goto 100
         endif
         if(el1>el2)then
             dcheck = el1 - upie
         else
             dcheck = el2 - dnie
         endif
         if(dcheck<0.02)then
             QPI(jp) = 0.0
             wlflag(jp) = .TRUE.
             goto 100
         endif
         if(ABS(qold(jp))<0.001)qold(jp) = 0.001
         if(el1>=(upie + dia) .AND. el2>=(dnie + dia))then
             d1 = el1
             d2 = el2
             goto 50
         endif
         if(el1>el2)then
             dtest = el2 - dnie
         else
             dtest = el1 - upie
         endif
         dcrit = DEPTHCRIT(ABS(qold(jp)))
         if(dtest<=dcrit)then
             if(el1<=el2)then
                 d1 = upie + dcrit
                 d2 = el2
             else
                 d1 = el1
                 d2 = dnie + dcrit
             endif
             vtot = 0.0
             tott = 0.0
             do
                 if(nit/=0)then
                     dtq = omega*dltx/vmax(jp)
                     if(dtq>(dlt - tott))then
                         dtq = dlt - tott
                     elseif((2.0*dtq)>(dlt - tott))then
                         dtq = (dlt - tott)*0.5
                     endif
                 endif
                 call OPEN_CHANNEL(d1, d2, QPI(jp), jp, dtq)
                 dcrit = DEPTHCRIT(ABS(QPI(jp)))
                 if(el1<=el2)then
                     d1 = upie + dcrit
                 else
                     d2 = dnie + dcrit
                 endif
                 vtot = vtot + dtq*QPI(jp)
                 tott = dtq + tott
                 if(tott>=(dlt - eps2))then
                     QPI(jp) = vtot/dlt
                     goto 100
                 endif
             enddo
         endif
         d1 = el1
         d2 = el2
50       tott = 0.0
         vtot = 0.0
         do
             if(nit/=0)then
                 dtq = omega*dltx/vmax(jp)
                 if(dtq>(dlt - tott))then
                     dtq = dlt - tott
                 elseif((2.0*dtq)>(dlt - tott))then
                     dtq = (dlt - tott)*0.5
                 endif
             endif
             call OPEN_CHANNEL(d1, d2, QPI(jp), jp, dtq)
             vtot = vtot + dtq*QPI(jp)
             tott = dtq + tott
             if(tott>=(dlt - eps2))then
                 QPI(jp) = vtot/dlt
                 exit
             endif
         enddo
100      qold(jp) = QPI(jp)
         if(QPI(jp)==0.0)wlflag(jp) = .TRUE.
     enddo
     return
     entry DEALLOCATE_PIPE_FLOW
     deallocate(begin, wlflag, vmax)
     end subroutine PIPE_FLOW_INITIALIZE
