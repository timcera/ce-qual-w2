!*==cemacalculaterisevelocity.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine CEMACALCULATERISEVELOCITY
 
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
     integer :: ngas, nrelarr
     real(8) :: rhog
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
                     if(BUBBLESSTATUS(segnumi, nrelarr)==0)then
 
                         do ngas = 1, numgas
                             BRVOLUAGAS(segnumi, nrelarr, ngas) = 0.D00
                             BRRATEAGAS(segnumi, nrelarr, ngas) = 0.D00
                             BRRATEAGASNET(segnumi, ngas) = 0.D00
                         enddo !nGas
 
                     endif
                 enddo
             enddo
         enddo
     enddo
 
 
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
                 do nrelarr = 1, numbubrelarr
                     if(BUBBLESSTATUS(segnumi, nrelarr)==1)then
 
                         rhog = 0.D0
                         do ngas = 1, numgas
                             rhog = rhog + BUBBLESGASCONC(segnumi, nrelarr,    &
                                  & ngas)/1000                                    !kg/m³
                         enddo !nGas
 
                         call CEMABUBBLESRISEVELOCITY                          &
                           & (BUBBLESRADIUS(segnumi, nrelarr), rhog,           &
                           & BUBBLESRISEV(segnumi, nrelarr))
 
                     endif
                 enddo
             enddo
         enddo
     enddo
 
     end subroutine CEMACALCULATERISEVELOCITY
