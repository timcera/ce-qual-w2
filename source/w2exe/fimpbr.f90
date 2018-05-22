!*==fimpbr.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************
!***********************************************************************
!*                                                                    **
!* S E G M E N T S   W I T H   C O N N E C T I N G   B R A N C H E S  **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls No Other Subroutines
 
     subroutine FIMPBR
 
 
 
 
 
     use FISHY
     use GEOMC
     use GLOBAL
     use SURFHE
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: brangle
!
!*** End of declarations rewritten by SPAG
!
 
!    *** ONLY HERE IS PHI0 used SW 1/14/01, DHS is used here also
 
     do i = 1, imx                                  !For all segments in the system
         nr = 0                                     ! Number of branches connecting to the right bank
         nl = 0                                     ! Number of branches connecting to the left bank
         RIMPBR(i, 1) = 0                           ! Right bank joining segment
         RIMPBR(i, 2) = 0                           ! Right bank joining branch
         LIMPBR(i, 1) = 0                           ! Left bank joining segment
         LIMPBR(i, 2) = 0                           ! Left bank joining branch
         do jb = 1, nbr                             !Examine all branches in the system
       ! IF(BR_INACTIVE(JB))CYCLE  ! SW 6/12/2017  No Need to do this for an inactive branch since this is called only once to set up geometery
             if(i==DHS(jb))then                     ! If find a connecting branch
                 if(PHI0(DS(jb))>PHI0(i))then
                     brangle = PHI0(DS(jb)) - PHI0(i)
                                                    ! BRANGLE = Angle of Incoming Branch relative to
                 else                               !           Branch I
                     brangle = PHI0(DS(jb)) + (2.*(22./7.) - PHI0(i))
                 endif
                 if((brangle>0) .AND. (brangle<(22./7.)))then   ! From BRANGLE, one can determine if the
                     nr = nr + 1                              !   Incoming Branch is on the right or
                     RIMPBR(i, 1) = DS(jb)                    !   left of Branch I
                     RIMPBR(i, 2) = jb
                 else
                     nl = nl + 1
                     LIMPBR(i, 1) = DS(jb)
                     LIMPBR(i, 2) = jb
                 endif
             endif
         enddo
         if((nr>1) .OR. (nl>1))then                             ! CHECK: Can only have one Branch join a segment from each side from the same side
             write(datadebugfn, *)'ERROR: More than 1 branch joining segment'
             write(datadebugfn, 9001)i
9001         format('ERROR occurred at segment ', i6)
             stop
         endif
     enddo
 
     end subroutine FIMPBR
