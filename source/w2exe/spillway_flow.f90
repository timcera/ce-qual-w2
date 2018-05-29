!*==spillway_flow.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                         S U B R O U T I N E   S P I L L W A Y  F L O W                                        **
!***********************************************************************************************************************************
 
     subroutine SPILLWAY_FLOW
     use F77KINDS
     use STRUCTURES
     use GLOBAL
     use GEOMC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND) :: dlel, elid, eliu, henergy, htail
     integer :: isub, js
!
!*** End of declarations rewritten by SPAG
!
 
     do js = 1, nsp
         if(LATERAL_SPILLWAY(js))then
             eliu = ELWS(IUSP(js))                                   !EL(KTWB(JWUSP(JS)),IUSP(JS))-Z(IUSP(JS))*COSA(BS(JWUSP(JS)))
         else
      ! ELIU =  ELWS(IUSP(JS))-SINA(JBUSP(JS))*DLX(IUSP(JS))*0.5      !EL(KTWB(JWUSP(JS)),IUSP(JS))-Z(IUSP(JS))*COSA(BS(JWUSP(JS)))-SINA(JBUSP(JS))*DLX(IUSP(JS))*0.5
             eliu = ELWS(IUSP(js)) + (ELWS(IUSP(js)) - ELWS(IUSP(js) - 1))     &
                  & /(0.5*(DLX(IUSP(js)) + DLX(IUSP(js) - 1)))*DLX(IUSP(js))   &
                  & *0.5                                                                                                     ! LINEAR INTERPOLATION OF THE WATER LEVEL TO THE EDGE
         endif
         if(IDSP(js)==0)then
             elid = -1.0
         elseif(US(JBDSP(js))/=IDSP(js))then
             elid = ELWS(IDSP(js))                                   !EL(KTWB(JWDSP(JS)),IDSP(JS))-Z(IDSP(JS))*COSA(BS(JWDSP(JS)))
         else
      !   ELID = ELWS(IDSP(JS))+SINA(JBDSP(JS))*DLX(IDSP(JS))*0.5     !EL(KTWB(JWDSP(JS)),IDSP(JS))-Z(IDSP(JS))*COSA(BS(JWDSP(JS)))+SINA(JBDSP(JS))*DLX(IDSP(JS))*0.5
             elid = ELWS(IDSP(js)) - (ELWS(IDSP(js) + 1) - ELWS(IDSP(js)))     &
                  & /(0.5*(DLX(IDSP(js)) + DLX(IDSP(js) + 1)))*DLX(IDSP(js))   &
                  & *0.5
         endif
         if(elid>=ESP(js) .OR. eliu>=ESP(js))then
             isub = 0
             if(A2SP(js)/=0.0 .AND. IDSP(js)/=0)then
                 htail = elid - ESP(js)                                         ! SW 5/10/05
                 if(htail>0)then
                     henergy = (U(KTWB(JWUSP(js)), IUSP(js))**2)/(2.0*g)       &
                             & + eliu - ESP(js)                                 ! SW 5/10/05
                     if(htail/henergy>0.67)isub = 1
                 endif
             endif
             if(isub==0)then
                 dlel = eliu - ESP(js)
                 if(dlel<0.0)then
                     dlel = -dlel
                     QSP(js) = -A1SP(js)*dlel**B1SP(js)
                 else
                     QSP(js) = A1SP(js)*dlel**B1SP(js)
                 endif
             elseif(elid>eliu)then
                 dlel = elid - eliu
                 QSP(js) = -A2SP(js)*dlel**B2SP(js)
             else
                 dlel = eliu - elid
                 QSP(js) = A2SP(js)*dlel**B2SP(js)
             endif
         else
             QSP(js) = 0.0
         endif
     enddo
     end subroutine SPILLWAY_FLOW
