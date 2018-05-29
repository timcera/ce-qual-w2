!*==acoustics.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************
!***********************************************************************
!*                                                                    **
!*  V I R T U A L   S A M P L I N G  -  H Y D R O A C O U S T I C S   **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls No Other Subroutines
 
 
     subroutine ACOUSTICS
     use FISHY
     use GEOMC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: anum
     real :: newxdist, newydistfromcenter, newzdist, oldxdist,                 &
           & oldydistfromcenter, oldzdist, sndylft, sndyrgt, sndzbot, sndztop, &
           & tagdepth, ydist
     real :: HAOPERAT(2, 13), OLDHALOC(2, 13)
!
!*** End of declarations rewritten by SPAG
!
 
 
                                                                      ! (3*NFISH) = Only a Guess; Impossible to Know How
                                                                      !   Many Fish will be Detected during an HA Survey
 
 
     oldxdist = NODES(fkmp, fimp, 1) + oldfxloc                       ! Previous Distance Downstream from Dam (m) to Location of Fish
     oldydistfromcenter = oldfyloc - B(fkmp, fimp)/2                  ! Previous Width Distance (m) from Lake Center
     oldzdist = (NODES(fkmp, fimp, 2) - NODES(ktwbf, fimp, 2)) - oldfzloc
!                                                                     ! Previous Depth (m) of Fish
     newxdist = NODES(fkmp, fimp, 1) + fxloc                          ! New Distance Downstream from Dam (m) to Location of Fish
     newydistfromcenter = fyloc - B(fkmp, fimp)/2                     ! New Width Distance (m) from Lake Center
     newzdist = (NODES(fkmp, fimp, 2) - NODES(ktwbf, fimp, 2)) - fzloc
                                                                      ! New Depth (m) of Fish
 
     do anum = 1, numacoustics                                        ! NUMACOUSTICS = Total # of Hydroacoustic Surveys
         if(HAOPERAT(anum, 13)==1)then                                ! Hydroacoustic Survey is Active
             sndylft = HAOPERAT(anum, 8)                              ! Meters from Center of Lake HA Beam will extend to Left Bank
             sndyrgt = HAOPERAT(anum, 9)                              ! Meters from Center of Lake HA Beam will extend to Right Bank
             sndztop = HAOPERAT(anum, 10)* - 1                        ! Depth Below Water Surface (m) to Top of HA Beam: (+) Number
             sndzbot = HAOPERAT(anum, 11)* - 1                        ! Depth Below Water Surface (m) to Bottom of HA Beam: (+) Number
             tagdepth = 0
                                                  !                    ! Fish May Have Crossed HA Beam Approaching from Upstream
             if(((oldxdist<=OLDHALOC(anum,1)) .AND. (newxdist>=HAOPERAT(anum,5)&
              & )) .OR.                                                        &
              & ((oldxdist>=OLDHALOC(anum,2)) .AND. (newxdist<=HAOPERAT(anum,6)&
              & )))then                                               ! Fish May Have Crossed HA Beam Approaching from Downstream
                 tagdepth = (newzdist + oldzdist)/2                   ! Approx. Depth of Fish when it crosses HA Beam
                 ydist = (newydistfromcenter + oldydistfromcenter)/2  ! Approx. Location of Fish in Width Dimension when crossing HA Beam
             endif
             if((tagdepth/=0) .OR.                                             &
              & ((newxdist>=HAOPERAT(anum,5)) .AND. (newxdist<=HAOPERAT(anum,6)&
              & )))then
                 if(tagdepth==0)then
                     tagdepth = newzdist
                     ydist = newydistfromcenter
                 endif
                 if(((ydist>=sndylft) .AND. (ydist<=sndyrgt)) .AND.            &
                  & ((tagdepth<=sndztop) .AND. (tagdepth>=sndzbot)))then
                                                                        ! Fish is Detected by Hydroacoustic Beam
                     SNDCOUNT(anum) = SNDCOUNT(anum) + 1                ! Tally # of Fish Detected by Hydroacoustic Beam
                     SNDCATCH(anum, SNDCOUNT(anum), 1) = SNDCOUNT(anum) ! Tally # of Fish Detected by Hydroacoustic Beam
                     SNDCATCH(anum, SNDCOUNT(anum), 2) = tagdepth       ! Depth (m) at Which Fish is Detected by Hydroacoustic Beam
                 endif
             endif
         endif
     enddo
 
     end subroutine ACOUSTICS
