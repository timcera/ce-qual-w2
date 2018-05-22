!*==cemabubblesrelease.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine CEMABUBBLESRELEASE
     use MAIN
     use GLOBAL
     use GEOMC
     use CEMAVARS
     use SCREENC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: BE, BS, BUBBLESATSURFACE, BUBBLESGASCONC, BUBBLESRADIUS,          &
           & BUBBLESRELEASEALLVALUE, BUBBLESSTATUS, CUS, DS
     real :: BRRATEAGAS, BRVOLUAGAS, BUBBLERELWB, BUBBLESCARRIED
     real :: brrateagasnet, dlt, segnumi
     real :: CEMAVARS, GEOMC, GLOBAL, SCREENC
     integer :: ICE, KTWB
     integer :: id, iu, jb, jw, kt, ngas, nrelarr, numbubrelarr, numgas, nwb
     integer :: MAIN
     real(8) :: tempbubblesrelvolume
!
!*** End of declarations rewritten by SPAG
!
 
    !IMPLICIT NONE
 
 
     brrateagasnet = 0.D00
 
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
                 do nrelarr = 1, numbubrelarr
                     if(BUBBLESSTATUS(segnumi, nrelarr)==1 .AND.               &
                      & BUBBLESATSURFACE(segnumi, nrelarr) .AND.               &
                      & .NOT.ICE(segnumi))then
 
                         do ngas = 1, numgas
                            !TempBubblesRelVolume = 4/3*3.14*BubblesRadius(SegNumI, nRelArr)**3
                             tempbubblesrelvolume = 4./3.*3.14*BUBBLESRADIUS   &
                               & (segnumi, nrelarr)**3                                              ! SW 10/10/2017
                             BRVOLUAGAS(segnumi, nrelarr, ngas)                &
                               & = BUBBLESRELEASEALLVALUE(segnumi, nrelarr)    &
                               & *tempbubblesrelvolume*BUBBLESGASCONC(segnumi, &
                               & nrelarr, ngas)                                                                                                                         !gm
                             BRRATEAGAS(segnumi, nrelarr, ngas)                &
                               & = BRVOLUAGAS(segnumi, nrelarr, ngas)/dlt                               !gm/s
                             brrateagasnet(segnumi, ngas)                      &
                               & = brrateagasnet(segnumi, ngas)                &
                               & + BRRATEAGAS(segnumi, nrelarr, ngas)                                                        !gm/s
                             BUBBLERELWB(jw, ngas) = BUBBLERELWB(jw, ngas)     &
                               & + dlt*brrateagasnet(segnumi, ngas)/1000.                                         ! SW 7/1/2017 Convert from gm/s to kg
                         enddo !nGas
 
                         BUBBLESCARRIED(segnumi, nrelarr)                      &
                           & = BUBBLESCARRIED(segnumi, nrelarr)                &
                           & - BUBBLESRELEASEALLVALUE(segnumi, nrelarr)
                         if(BUBBLESCARRIED(segnumi, nrelarr)<=0.D00)then
                             BUBBLESRELEASEALLVALUE(segnumi, nrelarr) = 0.D00
                             BUBBLESCARRIED(segnumi, nrelarr) = 0.D00
                             BUBBLESGASCONC(segnumi, nrelarr, :) = 0.D00
                             BUBBLESATSURFACE(segnumi, nrelarr) = .FALSE.
                             BUBBLESSTATUS(segnumi, nrelarr) = 0
                         endif
 
                     endif
                 enddo
             enddo
         enddo
     enddo
 
     end subroutine CEMABUBBLESRELEASE
