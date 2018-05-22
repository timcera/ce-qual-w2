!*==tag124578.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************
!***********************************************************************
!*                                                                    **
!*   I N F L U E N C E   O F   O N E   J O I N I N G   B R A N C H    **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls the Following Subroutines:
!      FINDNEWBR
 
 
     subroutine TAG124578
     use FISHY
     use GLOBAL
     use GEOMC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: avebrvel, cbid, cbidl, cbidr, flowhere, totbrarea, usbrflow, width
     integer :: fimptmp, fkmptmp, usbrimp
!
!*** End of declarations rewritten by SPAG
!
 
 
 
     if(tag==1)usbrimp = LIMPBR(fimp + 1, 1)       ! USBRIMP = UpStream left BRanch IMP
     if(tag==2)usbrimp = RIMPBR(fimp + 1, 1)       ! USBRIMP = UpStream right BRanch IMP
     if(tag==4)usbrimp = LIMPBR(fimp - 1, 1)       ! USBRIMP = UpStream left BRanch IMP
     if(tag==5)usbrimp = RIMPBR(fimp - 1, 1)       ! USBRIMP = UpStream right BRanch IMP
     if(tag==7)usbrimp = LIMPBR(fimp, 1)           ! USBRIMP = UpStream left BRanch IMP
     if(tag==8)usbrimp = RIMPBR(fimp, 1)           ! USBRIMP = UpStream right BRanch IMP
 
     if(debug)then
         write(datadebugfn, 9001)tag, fimp, usbrimp
9001     format(' TAG=', i4, ' INITIAL FIMP=', i6, ' USBRIMP=', i6)
     endif
     fkmptmp = KTI(usbrimp)
     usbrflow = 0.0                                ! USBRFLOW = UpStream (incoming) BRanch FLOW
     totbrarea = 0.0                               ! X-Sectional Area of incoming branch
     do
         if(fkmptmp<=KB(usbrimp))then
             usbrflow = usbrflow + U(fkmptmp, usbrimp)                         &
                      & *(B(fkmptmp, usbrimp)*H(fkmptmp, fjr))
             totbrarea = totbrarea + B(fkmptmp, usbrimp)*H(fkmptmp, fjr)
             fkmptmp = fkmptmp + 1
             cycle
         endif
         if((tag==1) .OR. (tag==2))then
             fkmptmp = KTI(fimp + 1)
             fimptmp = fimp + 1
         elseif((tag==4) .OR. (tag==5))then
             fkmptmp = KTI(fimp - 1)
             fimptmp = fimp - 1
         elseif((tag==7) .OR. (tag==8))then
             fkmptmp = KTI(fimp)
             fimptmp = fimp
         endif
         flowhere = 0.0                                           ! FLOWHERE = FLOW in current branch
         exit
     enddo
     do
         if(fkmptmp<=KB(fimptmp))then
             flowhere = flowhere + ABS(U(fkmptmp, fimptmp))                    &
                      & *(B(fkmptmp, fimptmp)*H(fkmptmp, fjr))
             fkmptmp = fkmptmp + 1
             cycle
         endif
         fkmptmp = fkmp                                           ! FKMPTMP=where fish is. FIMPTMP=where fish will be.
         exit
     enddo
     do
         if((fkmptmp<KTI(fimptmp)) .OR. (fkmptmp>KB(fimptmp)))then    ! Check to make sure there is an active cell
             if(debug)write(datadebugfn, *)' B EQ 0 CHECK ACTIVATED (1)' !   where calculating CBID
             if(fkmptmp<KTI(fimptmp))then
                 fkmptmp = fkmptmp + 1
             elseif(fkmptmp>KB(fimptmp))then
                 fkmptmp = fkmptmp - 1
             else
                 write(*, *)'ERROR: Problem finding where to calculate CBID'
             endif
             cycle
         endif
                                                    !   ! CBID = Combining Branch Influence Distance
         cbid = ABS(usbrflow)/(ABS(usbrflow) + flowhere)*B(fkmptmp, fimptmp)
                                                       !        into current branch
         avebrvel = usbrflow/totbrarea                 ! AVEBRVEL = Average BRanch VELocity, = Flow/X-sect Area
         if(debug)then
             write(datadebugfn, 9002)B(fkmptmp, fimptmp), cbid, avebrvel, fyloc
9002         format(' B=', f9.2, ' CBID=', f9.2, ' AVEBRVEL=', e9.3,           &
                   &' INITIAL FYLOC=', f9.3)
         endif
         if((tag==1) .OR. (tag==4) .OR. (tag==7))then                  !Incoming left branch
             cbidl = cbid                                              ! CBID from the left bank
             if((fyloc>0) .AND. (fyloc<cbidl))then
                 fyloc = fyloc + (avebrvel + fyvel + vfish)*dlt + rdy  ! New updated Fish Y-LOCation
                 if(debug)write(datadebugfn, *)                                &
                               &' FYLOC INFLUENCED BY LEFT BRANCH'
             endif
         elseif((tag==2) .OR. (tag==5) .OR. (tag==8))then              !Incoming right branch
             cbidr = B(fkmptmp, fimptmp) - cbid                        ! CBID from the right bank
             if((fyloc>cbidr) .AND. (fyloc<B(fkmptmp, fimptmp)))then
                 fyloc = fyloc - (avebrvel - fyvel - vfish)*dlt + rdy  ! New updated Fish Y-LOCation
                 if(debug)write(datadebugfn, *)                                &
                               &' FYLOC INFLUENCED BY RIGHT BRANCH'
             endif
         endif
         if(debug)write(datadebugfn, 9003)fyloc
9003     format(' THE NEW FYLOC=', f9.3)
 
!        Calculating New Fish Location: FIMP, FXLOC, FYLOC, and FNBP (and
 
!        FKMP, if necessary)
         width = B(fkmptmp, fimptmp)
         if((fyloc<0) .AND. ((tag==1) .OR. (tag==4) .OR. (tag==7)))then            ! ((CASE 1))
             if(debug)write(datadebugfn, *)' ((CASE 1)) ACTIVATED'
             fyloctmp = fyloc
             fimp = usbrimp
             do
                 if((fkmptmp<KTI(fimp)) .OR. (fkmptmp>KB(fimp)))then
                                                                    ! Check to make sure there is an active cell
                     if(debug)write(datadebugfn, *)                            &
                                   &' B EQ 0 CHECK ACTIVATED (2)'          !   where calculating new fish location
                     if(fkmptmp<KTI(fimp))then
                         fkmptmp = fkmptmp + 1
                     elseif(fkmptmp>KB(fimp))then
                         fkmptmp = fkmptmp - 1
                     else
                         write(*, *)                                           &
        &'ERROR: Problem finding upstream KMP for fish in left upstream branch'
                         stop
                     endif
                     cycle
                 endif
                 fyloc = B(fkmptmp, fimp)/2                         ! When fish moves into upstream branch I
                 call FINDNEWBR                                     !   arbitrarily set its new FYLOC as the middle
                 fxloc = DLX(fimp) + fyloctmp                       ! A negative FXLOC will be caught by the Subroutine FISH
                 exit
             enddo
         elseif(fyloc<0)then                                        ! ((CASE 2))
             if(debug)write(datadebugfn, *)' ((CASE 2)) ACTIVATED'
             fyloc = B(fkmptmp, fimptmp)*yrefl                      ! Logic for rebounding off of left bank
             if(tag==2)then                                         ! Moving downstream to FIMP+1
                 fxloc = fxloc - DLX(fimp)
                 fimp = fimptmp
             elseif(tag==5)then                                     ! Moving upstream to FIMP-1
                 fimp = fimptmp
                 fxloc = DLX(fimp) + fxloc
             elseif(tag/=8)then
            ! NOTHING
                 write(*, *)                                                   &
                   &'ERROR: Unable to calculate fish movement off of left bank'
                 stop
             endif
         elseif((fyloc>width) .AND. ((tag==2) .OR. (tag==5) .OR. (tag==8)))then
                                                                    ! ((CASE 3))
             if(debug)write(datadebugfn, *)' ((CASE 3)) ACTIVATED'
             fyloctmp = fyloc
             fimp = usbrimp
             do
                 if((fkmptmp<KTI(fimp)) .OR. (fkmptmp>KB(fimp)))then
                                                                    ! Check to make sure there is an active cell
                     if(debug)write(datadebugfn, *)                            &
                                   &' B EQ 0 CHECK ACTIVATED (3)'          !   where calculating new fish location
                     if(fkmptmp<KTI(fimp))then
                         fkmptmp = fkmptmp + 1
                     elseif(fkmptmp>KB(fimp))then
                         fkmptmp = fkmptmp - 1
                     else
                         write(*, *)                                           &
       &'ERROR: Problem finding upstream KMP for fish in right upstream branch'
                         stop
                     endif
                     cycle
                 endif
                 fyloc = B(fkmptmp, fimp)/2                         ! When fish moves into upstream branch I
                 call FINDNEWBR                                     !   arbitrarily set its new FYLOC as the middle
                 fxloc = DLX(fimp) - (fyloctmp - width)             ! A negative FXLOC will be caught by the Subroutine FISH
                 exit
             enddo
         elseif(fyloc>width)then                                    ! ((CASE 4))
             if(debug)write(datadebugfn, *)' ((CASE 4)) ACTIVATED'
             fyloc = B(fkmptmp, fimptmp)*(1 - yrefl)                ! Logic for rebounding off of right bank
             if(tag==1)then                                         ! Moving downstream to FIMP+1
                 fxloc = fxloc - DLX(fimp)
                 fimp = fimptmp
             elseif(tag==4)then                                     ! Moving upstream to FIMP-1
                 fimp = fimptmp
                 fxloc = DLX(fimp) + fxloc
             elseif(tag/=7)then
            ! NOTHING
                 write(*, *)                                                   &
                  &'ERROR: Unable to calculate fish movement off of right bank'
                 stop
             endif
         else                                                       ! ((CASE 5))
             if(debug)write(datadebugfn, *)' ((CASE 5)) ACTIVATED'
             if((tag==1) .OR. (tag==2))then                         ! Moving downstream to FIMP+1
                 fxloc = fxloc - DLX(fimp)
                 fimp = fimptmp
             elseif((tag==4) .OR. (tag==5))then                     ! Moving upstream to FIMP-1
                 fimp = fimptmp
                 fxloc = DLX(fimp) + fxloc
             elseif((tag/=7) .AND. (tag/=8))then
            ! NOTHING
                 write(*, *)                                                   &
                          &'ERROR: Unable to calculate fish FXLOC in TAG124578'
                 stop
             endif
         endif
         if(debug)then
             write(datadebugfn, 9004)fnbp, fimp, fxloc, fyloc
9004         format(' NEW: FNBP=', i6, ' FIMP=', i6, ' FXLOC=', f9.2,          &
                  & ' FYLOC=', f9.3)
         endif
         exit
     enddo
 
     end subroutine TAG124578
