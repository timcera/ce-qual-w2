!*==cemabubblestransport.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
     subroutine CEMABUBBLESTRANSPORT
 
     use MAIN
     use GLOBAL
     use GEOMC
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: bubblayer, nrelarr
     real(8) :: vdisttravbubble, vlocationbubble
!
!*** End of declarations rewritten by SPAG
!
 
 
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
                 do nrelarr = 1, numbubrelarr
                     if(BUBBLESSTATUS(segnumi, nrelarr)==1 .AND.               &
                      & .NOT.BUBBLESATSURFACE(segnumi, nrelarr))then
 
                         bubblayer = BUBBLESLNUMBER(segnumi, nrelarr)
                         vlocationbubble = 0.5*(EL(bubblayer + 1, segnumi)     &
                           & + EL(bubblayer, segnumi))
                         vdisttravbubble = BUBBLESRISEV(segnumi, nrelarr)*dlt
                         vlocationbubble = vlocationbubble + vdisttravbubble
 
                        !Locate vertical location
                         BUBBLESLNUMBER(segnumi, nrelarr) = kt
                         do k = kt, KB(segnumi)
                             if(vlocationbubble<EL(k, segnumi))                &
                              & BUBBLESLNUMBER(segnumi, nrelarr) = k
                         enddo !K
 
                         if(BUBBLESLNUMBER(segnumi, nrelarr)==kt)then
                             BUBBLESATSURFACE(segnumi, nrelarr) = .TRUE.
                             FIRSTBUBBLESRELEASE(segnumi, nrelarr) = .TRUE.
                             BUBBLESRELEASEALLVALUE(segnumi, nrelarr)          &
                               & = bubbrelfractionatm*BUBBLESCARRIED(segnumi,  &
                               & nrelarr)
                         endif
 
                     endif
                 enddo
             enddo
         enddo
     enddo
 
     end subroutine CEMABUBBLESTRANSPORT
