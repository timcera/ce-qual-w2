!*==whatjr.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************
!***********************************************************************
!*                                                                    **
!*        F I N D   W A T E R B O D Y   O F   L O C A T I O N         **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls No Other Subroutines
 
     subroutine WHATJR
 
     use FISHY
     use GLOBAL
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: wbdy
!
!*** End of declarations rewritten by SPAG
!
 
 
     wbdy = 0                                               ! WBDY = Water Body
     do
         wbdy = wbdy + 1
         if(wbdy>nwb)then                                   ! CHECK: WBDY can't be more
             write(datadebugfn, *)'ERROR: Finding current WB in Subroutine'
                                                                      !        than NWB
             stop
         endif
         do jb = BS(wbdy), BE(wbdy)                         ! Branches within a particular
             if((fimp>=US(jb)) .AND. (fimp<=DS(jb)))then    !   Water Body
                 fjr = wbdy
                 goto 99999
             endif
         enddo
     enddo
 
 
99999 end subroutine WHATJR
