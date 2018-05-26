!*==cemabubblesturbulence.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine CEMABUBBLESTURBULENCE
     use MAIN
     use GLOBAL
     use GEOMC
     use SCREENC
     use KINETIC
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: bubbcntr, bubblnumber, nrelarr
     real(8) :: tempbubbdiam, temprelvelocity
!
!*** End of declarations rewritten by SPAG
!
 
     segnumi = i
     do k = kt, KBMIN(segnumi) - 1
         bubbcntr = 0
         tempbubbdiam = 0
         temprelvelocity = 0
         do nrelarr = 1, numbubrelarr
             if(BUBBLESSTATUS(segnumi, nrelarr)==1)then
                 bubblnumber = BUBBLESLNUMBER(segnumi, nrelarr)
 
                 if(bubblnumber==k)then
                     bubbcntr = bubbcntr + 1
                     tempbubbdiam = tempbubbdiam + 2.0*BUBBLESRADIUS(segnumi,  &
                                  & nrelarr)
                     temprelvelocity = temprelvelocity +                       &
                                     & BUBBLESRISEV(segnumi, nrelarr)          &
                                     & + W(k - 1, segnumi)                        !Rise velocity is +ve upwards and W is +ve downwards
                 endif
 
             endif
         enddo !nRelArr
 
         if(bubbcntr>0)then
             tempbubbdiam = tempbubbdiam/bubbcntr
             temprelvelocity = temprelvelocity/bubbcntr
             if(k==KB(segnumi))then
                 AZ(k, segnumi) = AZ(k, segnumi) + BOTTOMTURBULENCE(segnumi)
             else
                !AZ(K,SegNumI) = AZ(K,SegNumI)  + TempBubbDiam/TempRelVelocity
                 AZ(k, segnumi) = AZ(k, segnumi) + tempbubbdiam*temprelvelocity  ! cb 2/7/13
             endif
         endif
 
         AZ(KB(segnumi) - 1, segnumi) = AZ(KB(segnumi) - 1, segnumi)           &
                                      & + BOTTOMTURBULENCE(segnumi)
 
     enddo
 
     end subroutine CEMABUBBLESTURBULENCE
