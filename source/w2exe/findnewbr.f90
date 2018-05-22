!*==findnewbr.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************
!***********************************************************************
!*                                                                    **
!*   F I N D   C U R R E N T   B R A N C H   O F   L O C A T I O N    **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls No Other Subroutines
 
 
     subroutine FINDNEWBR
 
 
 
 
 
     use FISHY
     use GEOMC
     use GLOBAL
     use SCREENC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: trib
!
!*** End of declarations rewritten by SPAG
!
 
100   trib = 1                                              ! TRIB = Tributary (or Branch)
     do while ((fimp<US(trib)) .OR. (fimp>DS(trib)))
         trib = trib + 1
         if(trib>nbr)then                                     ! CHECK: TRIB can't be more
             write(datadebugfn, *)'ERROR: Finding new FNBP in Subroutine'
                                                                        !        than NBR # BRANCHES
             write(datadebugfn, *)'ERROR: FIMP=', fimp, ' JDAY=', jday
             write(datadebugfn, *)'FIMP=FIMP-1'
             fimp = fimp - 1
             goto 100
          !STOP
         endif
     enddo
     fnbp = trib
 
     end subroutine FINDNEWBR
