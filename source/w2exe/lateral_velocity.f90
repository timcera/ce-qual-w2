!*==lateral_velocity.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!**********************************
     subroutine LATERAL_VELOCITY
! SW 2/01/02   7/31/2017
 
 
     use FISHY
     use GLOBAL
     use GEOMC
     use MAIN, ONLY:EV, QPR, IWD
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: iw
     logical :: with_at_particle_location
     real :: xarea
!
!*** End of declarations rewritten by SPAG
!
 
!    Concept all QSS from each cell will be treated as a lateral withdrawal -
!    each withdrawal will be assigned a RHS or LHS looking downstream
!    location; lateral velocity origin is the segment/cell center. Velocities
!    to the RHS are + and those to the LHS are - IF RHS LOOKING DOWNSTREAM
!    THEN VELOCITY IS POSITIVE, RIMPBR(I,1)>0 IF LHS LOOKING DOWNSTREAM THEN
   ! TRIBS AND INFLOWS IN QSS ARE ASSIGNED POSITIVE VALUES, WITHDRAWALS ARE NEGATIVE
   !
!    VELOCITY IS NEGATIVE, LIMPBR(I,1)>0
     with_at_particle_location = .FALSE.
!    CHECK FOR WITHDRAWALS AT FIMP
     do iw = 1, nwd
         if(fimp==IWD(iw))then
             with_at_particle_location = .TRUE.
             exit
         endif
     enddo
 
     if(RIMPBR(fimp, 1)==0 .AND. LIMPBR(fimp, 1)==0 .AND.                      &
      & .NOT.with_at_particle_location)then
         fyvel = 0.0
 
     elseif(fkmp<=KTWB(fjr))then
         xarea = H1(KTWB(fjr), fimp)*DLX(fimp)
         fyvel = (QSS(KTWB(fjr), fimp) + EV(fimp) - QPR(fimp))/xarea    ! ADD EVAPORATION AND REMOVE PRECIP BACK TO QSS...be careful about evaporation which is also a -QSS flow
 
         if(RIMPBR(fimp, 1)>0)then
             if(BR_INACTIVE(RIMPBR(fimp, 2)))then
                 fyvel = 0.0
             else
                 fyvel = -fyvel
             endif
         elseif(LIMPBR(fimp, 1)>0)then
             if(BR_INACTIVE(LIMPBR(fimp, 2)))then
                 fyvel = 0.0
             else
                 fyvel = fyvel
             endif
         endif
 
     else
         xarea = H1(fkmp, fimp)*DLX(fimp)
         fyvel = QSS(fkmp, fimp)/xarea
 
         if(RIMPBR(fimp, 1)>0)then
             if(.NOT.BR_INACTIVE(RIMPBR(fimp, 2)) .AND.                        &
              & fkmp<=KB(DS(RIMPBR(fimp,2))))then
                 fyvel = -fyvel
             else
                 fyvel = 0.0
             endif
         elseif(LIMPBR(fimp, 1)>0)then
             if(.NOT.BR_INACTIVE(LIMPBR(fimp, 2)) .AND.                        &
              & fkmp<=KB(DS(LIMPBR(fimp,2))))then
                 fyvel = fyvel
             else
                 fyvel = 0.0
             endif
         endif
 
     endif
 
     end subroutine LATERAL_VELOCITY
