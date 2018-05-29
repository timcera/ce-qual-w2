!*==initial_u_velocity.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************************************************************************
!**        S U B R O U T I N E    I N I T I A L    H O R I Z O N T A L    V E L O C I T Y                                         **
!***********************************************************************************************************************************
 
     subroutine INITIAL_U_VELOCITY
     use GLOBAL
     use GEOMC
     use INITIALVELOCITY
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: k
     real :: wsurf, xarea
!
!*** End of declarations rewritten by SPAG
!
 
 
 
 
     do jw = 1, nwb
         kt = KTWB(jw)
 
         do jb = BS(jw), BE(jw)
             if(SLOPE(jb)>0.0 .AND. .NOT.LOOP_BRANCH(jb))then
                 iu = CUS(jb)
                 id = DS(jb)
                 do i = iu, id
                     wsurf = ELWS(i)
                     call XSECTIONAL_AREA(wsurf, xarea)
                     UAVG(i) = QSSI(i)/xarea
                     do k = kt, KB(i)
                         U(k, i) = UAVG(i)
                     enddo
                 enddo
             endif
         enddo
     enddo
 
     end subroutine INITIAL_U_VELOCITY
