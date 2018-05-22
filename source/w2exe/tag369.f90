!*==tag369.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************
!***********************************************************************
!*                                                                    **
!*  I N F L U E N C E   O F   T W O   J O I N I N G   B R A N C H E S **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls the Following Subroutines:
!      FINDNEWBR
 
 
     subroutine TAG369
 
 
 
 
     use FISHY
     use GLOBAL
     use GEOMC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real, dimension(2) :: avebrvel, cbid, totbrarea, usbrflow
     real :: cbidr, flowhere, width
     integer :: fimptmp, jbr
     integer, dimension(2) :: fkmptmp, usbrimp
!
!*** End of declarations rewritten by SPAG
!
 
 
 
 
     if(tag==3)then
         usbrimp(1) = LIMPBR(fimp + 1, 1)          ! USBRIMP(1) = UpStream left BRanch IMP
         usbrimp(2) = RIMPBR(fimp + 1, 1)          ! USBRIMP(2) = UpStream right BRanch IMP
     elseif(tag==6)then
         usbrimp(1) = LIMPBR(fimp - 1, 1)          ! USBRIMP(1) = UpStream left BRanch IMP
         usbrimp(2) = RIMPBR(fimp - 1, 1)          ! USBRIMP(2) = UpStream right BRanch IMP
     elseif(tag==9)then
         usbrimp(1) = LIMPBR(fimp, 1)              ! USBRIMP(1) = UpStream left BRanch IMP
         usbrimp(2) = RIMPBR(fimp, 1)              ! USBRIMP(2) = UpStream right BRanch IMP
     endif
     fkmptmp(1) = KTI(usbrimp(1))                  ! Surface Layer of incoming left branch
     fkmptmp(2) = KTI(usbrimp(2))                  ! Surface Layer of incoming right branch
     do jbr = 1, 2                                  ! BR: =1 incoming left branch, =2 incoming right branch
         usbrflow(jbr) = 0.0                        ! USBRFLOW = UpStream (incoming) BRanch FLOW
         totbrarea(jbr) = 0.0                       ! X-Sectional Area of incoming branch
         do while (fkmptmp(jbr)<=KB(usbrimp(jbr)))
             usbrflow(jbr) = usbrflow(jbr) + U(fkmptmp(jbr), usbrimp(jbr))     &
                           & *(B(fkmptmp(jbr), usbrimp(jbr))                   &
                           & *H(fkmptmp(jbr), fjr))
             totbrarea(jbr) = totbrarea(jbr) + B(fkmptmp(jbr), usbrimp(jbr))   &
                            & *H(fkmptmp(jbr), fjr)
             fkmptmp(jbr) = fkmptmp(jbr) + 1
         enddo
     enddo
     if(tag==3)then
         fkmptemp = KTI(fimp + 1)
         fimptmp = fimp + 1
     elseif(tag==6)then
         fkmptemp = KTI(fimp - 1)
         fimptmp = fimp - 1
     elseif(tag==9)then
         fkmptemp = KTI(fimp)
         fimptmp = fimp
     endif
     flowhere = 0.0                                       ! FLOWHERE = FLOW in current branch
     do
         if(fkmptemp<=KB(fimptmp))then
             flowhere = flowhere + ABS(U(fkmptemp, fimptmp))                   &
                      & *(B(fkmptemp, fimptmp)*H(fkmptemp, fjr))
             fkmptemp = fkmptemp + 1
             cycle
         endif
         fkmptemp = fkmp                                  ! FKMPTEMP=where fish is. FIMPTMP=where fish will be.
         exit
     enddo
     do
         if((fkmptemp<KTI(fimptmp)) .OR. (fkmptemp>KB(fimptmp)))then
                                                          ! Check to make sure there is an active cell
             if(fkmptemp<KTI(fimptmp))then                !   where calculating CBID
                 fkmptemp = fkmptemp + 1
             elseif(fkmptemp>KB(fimptmp))then
                 fkmptemp = fkmptemp - 1
             else
                 write(*, *)                                                   &
     &'ERROR: Problem finding where to calculate CBID with 2 incoming branches'
             endif
             cycle
         endif
         do jbr = 1, 2
             cbid(jbr) = ABS(usbrflow(jbr))                                    &
                       & /((ABS(usbrflow(1)) + ABS(usbrflow(2)) + flowhere)    &
                       & *B(fkmptemp, fimptmp))           ! CBID = Combining Branch Influence Distance
             avebrvel(jbr) = usbrflow(jbr)/totbrarea(jbr)    ! AVEBRVEL = Average BRanch VELocity, = Flow/X-sect Area
         enddo
         width = B(fkmptemp, fimptmp)
         if((fyloc>0) .AND. (fyloc<cbid(1)))then
             fyloc = fyloc + (avebrvel(1) + fyvel + vfish)*dlt + rdy
                                                          ! New updated Fish Y-LOCation
         else
             cbidr = width - cbid(2)                      ! CBIDR = CBID from the right bank
             if((fyloc>cbidr) .AND. (fyloc<width))fyloc = fyloc -              &
              & (avebrvel(2) - fyvel - vfish)*dlt + rdy   ! New updated Fish Y-LOCation
         endif
 
!        Calculating New Fish Location: FIMP, FXLOC, FYLOC, and FNBP (and
 
!        FKMP, if necessary)
         if(fyloc<0)then
             fyloctmp = fyloc
             fimp = usbrimp(1)
             do
                 if((fkmptemp<KTI(fimp)) .OR. (fkmptemp>KB(fimp)))then
                                                          ! Check to make sure there is an active cell
                     if(fkmptemp<KTI(fimp))then           !   where calculating new fish location
                         fkmptemp = fkmptemp + 1
                     elseif(fkmptemp>KB(fimp))then
                         fkmptemp = fkmptemp - 1
                     else
                         write(*, *)                                           &
     &'ERROR: Problem finding upstream KMP for fish in left upstream branch wit&
     &h 2 incoming branches'
                         stop
                     endif
                     cycle
                 endif
                 fyloc = B(fkmptemp, fimp)/2              ! When fish moves into upstream branch I
                 call FINDNEWBR                           !   arbitrarily set its new FYLOC as the middle
                 fxloc = DLX(fimp) + fyloctmp             ! A negative FXLOC will be caught by the Subroutine FISH
                 exit
             enddo
         elseif(fyloc>width)then
             fyloctmp = fyloc
             fimp = usbrimp(2)
             do
                 if((fkmptemp<KTI(fimp)) .OR. (fkmptemp>KB(fimp)))then
                                                          ! Check to make sure there is an active cell
                     if(fkmptemp<KTI(fimp))then           !   where calculating new fish location
                         fkmptemp = fkmptemp + 1
                     elseif(fkmptemp>KB(fimp))then
                         fkmptemp = fkmptemp - 1
                     else
                         write(*, *)                                           &
     &'ERROR: Problem finding upstream KMP for fish in right upstream branch wi&
     &th 2 incoming branches'
                         stop
                     endif
                     cycle
                 endif
                 fyloc = B(fkmptemp, fimp)/2              ! When fish moves into upstream branch I
                 call FINDNEWBR                           !   arbitrarily set its new FYLOC as the middle
                 fxloc = DLX(fimp) - (fyloctmp - width)   ! A negative FXLOC will be caught by the Subroutine FISH
                 exit
             enddo
         elseif(tag==3)then                               ! Moving downstream to FIMP+1
             fxloc = fxloc - DLX(fimp)
             fimp = fimptmp
         elseif(tag==6)then                               ! Moving upstream to FIMP-1
             fimp = fimptmp
             fxloc = DLX(fimp) + fxloc
         elseif(tag/=9)then
            ! NOTHING
             write(*, *)'ERROR: Unable to calculate fish FXLOC in TAG369'
             stop
         endif
         exit
     enddo
 
     end subroutine TAG369
