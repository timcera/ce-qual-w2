!*==gate_flow.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                          S U B R O U T I N E   G A T E  F L O W                                               **
!***********************************************************************************************************************************
 
     subroutine GATE_FLOW
     use STRUCTURES
     use GLOBAL
     use GEOMC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(r8) :: dlel, elid, eliu, henergy, htail
     integer :: igt, isub, jg
!
!*** End of declarations rewritten by SPAG
!
 
 
     do jg = 1, ngt
         if(DYNGTC(jg)=='    FLOW')then
             QGT(jg) = BGT(jg)
         else
 
 !   ELIU =  ELWS(IUGT(JG))       !EL(KTWB(JWUGT(JG)),IUGT(JG))-Z(IUGT(JG))*COSA(BS(JWUGT(JG)))
             if(LATERAL_GATE(jg))then
                 eliu = ELWS(IUGT(jg))
                                  !EL(KTWB(JWUGT(JG)),IUGT(JG))-Z(IUGT(JG))*COSA(BS(JWUGT(JG)))
             else
      !ELIU=ELWS(IUGT(JG))-SINA(JBUGT(JG))*DLX(IUGT(JG))*0.5     !EL(KTWB(JWUGT(JG)),IUGT(JG))-Z(IUGT(JG))*COSA(BS(JWUGT(JG)))-SINA(JBUGT(JG))*DLX(IUGT(JG))*0.5
                 eliu = ELWS(IUGT(jg)) + (ELWS(IUGT(jg)) - ELWS(IUGT(jg) - 1)) &
                      & /(0.5*(DLX(IUGT(jg)) + DLX(IUGT(jg) - 1)))             &
                      & *DLX(IUGT(jg))*0.5                                                                                 ! LINEAR INTERPOLATION OF THE WATER LEVEL TO THE EDGE
             endif
             if(IDGT(jg)==0)then
                 elid = -100.0
             elseif(US(JBDGT(jg))/=IDGT(jg))then
                 elid = ELWS(IDGT(jg))
                                  !EL(KTWB(JWDGT(JG)),IDGT(JG))-Z(IDGT(JG))*COSA(BS(JWDGT(JG)))
             else
        !ELID = ELWS(IDGT(JG))+SINA(JBDGT(JG))*DLX(IDGT(JG))*0.5          !EL(KTWB(JWDGT(JG)),IDGT(JG))-Z(IDGT(JG))*COSA(BS(JWDGT(JG)))+SINA(JBDGT(JG))*DLX(IDGT(JG))*0.5
                 elid = ELWS(IDGT(jg)) - (ELWS(IDGT(jg) + 1) - ELWS(IDGT(jg))) &
                      & /(0.5*(DLX(IDGT(jg)) + DLX(IDGT(jg) + 1)))             &
                      & *DLX(IDGT(jg))*0.5
             endif
             if(BGT(jg)==0.0)then
                 QGT(jg) = 0.0
             elseif(elid>EGT(jg) .OR. eliu>EGT(jg))then
                 isub = 0
                 if(A2GT(jg)/=0.0 .AND. IDGT(jg)/=0)then                                 ! SW 8/21/2013
                     htail = elid - EGT(jg)                                              ! SW 5/10/05
                     if(htail>0)then
                         henergy = (U(KTWB(JWUGT(jg)), IUGT(jg))**2)/(2.0*g)   &
                                 & + eliu - EGT(jg)                                      ! SW 5/10/05
                         if(htail/henergy>0.67)isub = 1
                     endif
                 endif
                 igt = 0
                 if(BGT(jg)>=0.8*(eliu - EGT(jg)) .AND. GTA1(jg)/=0.0)igt = 1
                 if(igt==0)then
                     if(isub==0)then
                         dlel = eliu - EGT(jg)
                         if(A2GT(jg)==0.0 .AND. G2GT(jg)/=0.0)dlel = eliu -    &
                          & G2GT(jg)
                         if(dlel<0.0)then
                             dlel = -dlel
                             QGT(jg) = -A1GT(jg)*(dlel**B1GT(jg))*BGT(jg)      &
                                     & **G1GT(jg)
                         else
                             QGT(jg) = A1GT(jg)*(dlel**B1GT(jg))*BGT(jg)       &
                                     & **G1GT(jg)
                         endif
                     elseif(elid>eliu)then
                         dlel = elid - eliu
                         QGT(jg) = -A2GT(jg)*dlel**B2GT(jg)*BGT(jg)**G2GT(jg)
                     else
                         dlel = eliu - elid
                         QGT(jg) = A2GT(jg)*dlel**B2GT(jg)*BGT(jg)**G2GT(jg)
                     endif
                 elseif(isub==0)then
                     dlel = eliu - EGT(jg)
                     if(elid>EGT(jg))dlel = eliu - elid
                     if(dlel<0.0)then
                         dlel = -dlel
                         QGT(jg) = -GTA1(jg)*dlel**GTB1(jg)
                     else
                         QGT(jg) = GTA1(jg)*dlel**GTB1(jg)
                     endif
                 elseif(elid>eliu)then
                     dlel = elid - eliu
                     QGT(jg) = -GTA2(jg)*dlel**GTB2(jg)
                 else
                     dlel = eliu - elid
                     QGT(jg) = GTA2(jg)*dlel**GTB2(jg)
                 endif
             else
                 QGT(jg) = 0.0
             endif
         endif
     enddo
     end subroutine GATE_FLOW
