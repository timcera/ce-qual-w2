!*==vgillnets.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************
!***********************************************************************
!*                                                                    **
!*       V I R T U A L   S A M P L I N G   -   G I L L N E T S        **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls the Following Subroutines:
!      SPLINE
 
     subroutine VGILLNETS
     use GLOBAL
     use FISHY
     use GEOMC
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND) :: dovalue, tempvalue
     real :: fractionx, netxput, netylft, netyrgt, netzbot, netztop, newxdist, &
           & newydistfromcenter, newzdist, oldxdist, oldydistfromcenter,       &
           & oldzdist, xlocatnet, ydistatnet, zdistatnet, zlocatnet
     integer :: nnum
!
!*** End of declarations rewritten by SPAG
!
 
 
 
 
     oldxdist = NODES(fkmp, fimp, 1) + oldfxloc                       ! Previous Distance Downstream from Dam (m) to Location of Fish
     oldydistfromcenter = oldfyloc - B(fkmp, fimp)/2                  ! Previous Width Distance (m) from Lake Center
     oldzdist = (NODES(fkmp, fimp, 2) - NODES(ktwbf, fimp, 2)) - oldfzloc
                                                                      ! Previous Depth (m) of Fish
     newxdist = NODES(fkmp, fimp, 1) + fxloc                          ! New Distance Downstream from Dam (m) to Location of Fish
     newydistfromcenter = fyloc - B(fkmp, fimp)/2                     ! New Width Distance (m) from Lake Center
     newzdist = (NODES(fkmp, fimp, 2) - NODES(ktwbf, fimp, 2)) - fzloc
                                                                      ! New Depth (m) of Fish
 
     do nnum = 1, numgillnets                                         ! NUMGILLNETS = Total Number of Gillnets
         if((GILLNETS(nnum, 11)==0) .AND. (fsnag==nnum))then          ! Release all Fish caught in a gillnet no longer active
             fsnag = 0
         elseif((GILLNETS(nnum, 11)==1) .AND. (fsnag==nnum))then      ! All Fish caught in a gillnet do not move
             fxloc = oldfxloc
             fyloc = oldfyloc
             fzloc = oldfzloc
         elseif((GILLNETS(nnum, 11)==1) .AND. (fsnag==0))then         ! Test to see if the free Fish will get snagged in a gillnet
             netxput = GILLNETS(nnum, 5)                              ! X-Location (m) Downstream Where Gillnet is Placed
             netylft = GILLNETS(nnum, 6)                              ! Meters from Center of Lake Gillnet will extend to Left Bank
             netyrgt = GILLNETS(nnum, 7)                              ! Meters from Center of Lake Gillnet will extend to Right Bank
             netztop = GILLNETS(nnum, 8)* - 1                         ! Depth Below Water Surface (m) to Top of Gillnet: (+) Number
             netzbot = GILLNETS(nnum, 9)* - 1                         ! Depth Below Water Surface (m) to Bottom of Gillnet: (+) Number
             fractionx = 0
             if((oldxdist<=netxput) .AND. (newxdist>=netxput))then    ! Fish May be Caught in Gillnet Approaching from Upstream
                 fractionx = (netxput - oldxdist)/(newxdist - oldxdist)
             elseif((oldxdist>=netxput) .AND. (newxdist<=netxput))then
                                                                      ! Fish May be Caught in Gillnet Approaching from Downstream
                 fractionx = (oldxdist - netxput)/(oldxdist - newxdist)
             endif
             if(fractionx/=0)then                                     ! Fish May be Caught in Gillnet
                 ydistatnet = oldydistfromcenter +                             &
                            & (newydistfromcenter - oldydistfromcenter)        &
                            & *fractionx                              ! Location of Fish in Width Dimension when it crosses Gillnet
                 zdistatnet = oldzdist + (newzdist - oldzdist)*fractionx
                                                                      ! Depth of Fish when it crosses Gillnet
                 if(((ydistatnet>=netylft) .AND. (ydistatnet<=netyrgt)) .AND.  &
                  & ((zdistatnet<=netztop) .AND. (zdistatnet>=netzbot)))then
                                                                      ! Fish has just been caught in Gillnet
                     fsnag = nnum                                     ! Fish caught in Gillnet # NNUM
                     SNAGCOUNT(nnum) = SNAGCOUNT(nnum) + 1            ! Tally # of Fish caught in Gillnet # NNUM
                     fxloc = netxput - NODES(fkmp, fimp, 1)           ! Location of Fish in Gillnet
                     fyloc = ydistatnet + B(fkmp, fimp)/2             ! Location of Fish in Gillnet
                     fzloc = zdistatnet -                                      &
                           & (NODES(fkmp, fimp, 2) - NODES(ktwbf, fimp, 2))
                                                                       ! Location of Fish in Gillnet
                     xlocatnet = netxput - NODES(fkmp, fimp, 1)
                     zlocatnet = -1*(zdistatnet - NODES(fkmp, fimp, 2))
                     ilok = fimp                                      ! Segment I where Variable Value is desired
                     xlok = xlocatnet                                 ! Location within Segment I where Variable Value is desired
                     klok = fkmp                                      ! Layer K where Variable Value is desired
                     zlok = zlocatnet                                 ! Location within Layer K where Variable Value is desired
                     NETCATCH(nnum, SNAGCOUNT(nnum), 1) = SNAGCOUNT(nnum)
                                                                      ! Tally # of Fish Caught
                     NETCATCH(nnum, SNAGCOUNT(nnum), 2) = netxput     ! X-Location (m) of Gillnet Downstream from RBR Dam
                     NETCATCH(nnum, SNAGCOUNT(nnum), 3) = zdistatnet  ! Depth (m) at Which the Fish was Caught
                     variable = 3
                         ! SW 1/15/01
                     call SPLINE                                      ! VARIABLE: 1=Horz Vel, 2=Vert Vel, 3=Temp, 4=Diss Oxy
                     NETCATCH(nnum, SNAGCOUNT(nnum), 4) = value   ! Temperature Where the Fish was Caught
                     variable = 4
                         ! SW 1/15/01
                     call SPLINE                                      ! VARIABLE: 1=Horz Vel, 2=Vert Vel, 3=Temp, 4=Diss Oxy
                     NETCATCH(nnum, SNAGCOUNT(nnum), 5) = value     ! Dissolved Oxygen Where the Fish was Caught
                 endif
             endif
         endif
     enddo
 
     end subroutine VGILLNETS
