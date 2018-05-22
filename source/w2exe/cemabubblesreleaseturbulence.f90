!*==cemabubblesreleaseturbulence.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine CEMABUBBLESRELEASETURBULENCE
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
 
     bottomturbulence = 0.D00
     do jw = 1, nwb
         k = KB(segnumi)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
 
                 bubbcntr = 0
                 tempbubbdiam = 0
                 temprelvelocity = 0
                 do nrelarr = 1, numbubrelarr
                     if(BUBBLESSTATUS(segnumi, nrelarr)==1)then
                         bubblnumber = BUBBLESLNUMBER(segnumi, nrelarr)
 
                         if(bubblnumber==k)then
                             bubbcntr = bubbcntr + 1
                             tempbubbdiam = tempbubbdiam +                     &
                               & 2.0*BUBBLESRADIUS(segnumi, nrelarr)
                             temprelvelocity = temprelvelocity +               &
                               & BUBBLESRISEV(segnumi, nrelarr)                &
                               & + W(k - 1, segnumi)                                                            !Rise velocity is +ve upwards and W is +ve downwards
                         endif
 
                     endif
                 enddo !nRelArr
 
                 if(bubbcntr>0)then
                     tempbubbdiam = tempbubbdiam/bubbcntr
                     temprelvelocity = temprelvelocity/bubbcntr
                     bottomturbulence(segnumi)                                 &
                       & = cematurbulencescaling*tempbubbdiam/temprelvelocity
                 endif
 
             enddo !SegNumI
         enddo  !JB
     enddo !JW
 
     do jw = 1, nwb
         k = KB(segnumi)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
 
                 if(bottomturbulence(segnumi)>0)then
                 endif
 
             enddo
         enddo
     enddo
 
 
     end subroutine CEMABUBBLESRELEASETURBULENCE
