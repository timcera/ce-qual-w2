!*==fishoutput.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!*****************************************
     subroutine FISHOUTPUT
!*****************************************
! Numerical Fish Surrogate Output
!*****************************************
 
!*****************************************
! Numerical Fish Surrogate Output
!*****************************************
 
     use FISHY
     use GLOBAL
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: avedepth, aveintdo, aveinttemp, datend, datstrt, dintcount,       &
           & haspeed, hrsend, hrsstrt, maxdo, maxtemp, mindo, mintemp,         &
           & netdropday, netdrophr, netpullday, netpullhr, netx, netylft,      &
           & netyrgt, netzbot, netztop, sndylft, sndyrgt, sndzbot, sndztop
     integer :: catch, catchdepth, dint, dintlast, dotally, fij, fvar2, fvar3, &
              & jf, tagged, taggeddepth, temptally
!
!*** End of declarations rewritten by SPAG
!
 
 
 
!    final fish output
     write(finalfn, '(a291)')                                                  &
     &'Part#,Seg#,XLocationwithinSegmentfromUpstreamSide(m),Layer#,VerticalDist&
     &fromTop(m),LateralDistfromLeftBank,Branch#,ParticleInModel(=0),JDAYleftsy&
     &stem,DetentionTime(days),RemovalMechanism,SedVelocity(m/d),DateStart'
     do jf = 1, nfish
         write(finalfn, '(i7,",",1x,10(f10.3,","),f8.4,",",f12.4)')jf,         &
             & FISHES(jf, 1), FISHES(jf, 2), FISHES(jf, 3), FISHES(jf, 4),     &
             & FISHES(jf, 5), FISHES(jf, 6), FISHES(jf, 12), FISHES(jf, 14),   &
             & FISHES(jf, 15), FISHES(jf, 16), SEDVEL(jf), DELAYDATE(jf)
     enddo
 
     close(finalfn)
!    End fish output section
                                                             !FISH
 
     end subroutine FISHOUTPUT
