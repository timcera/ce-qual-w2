!*==open_channel_initialize.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                        S U B R O U T I N E   O P E N  C H A N N E L                                           **
!***********************************************************************************************************************************
 
     subroutine OPEN_CHANNEL_INITIALIZE
     use GLOBAL
     use STRUCTURES
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
     real(R8KIND), parameter :: THETA = 0.55
!
! Dummy arguments
!
     real(R8KIND) :: Dt, El1, El2, Qout
     integer :: Ic
     intent (in) Dt, El1, El2
     intent (inout) Qout
!
! Local variables
!
     real(R8KIND), allocatable, dimension(:, :), save :: al, daa
     real(R8KIND), allocatable, dimension(:), save :: b, belev, carea, q, rt,      &
           & tarea, topw, topwt, v, vold, vpr, vt, y, yold, ypr, yt
     real(R8KIND), save :: bar1, bar2, bc1, bc2, bepr1, bepr2, d, dist, dltx,      &
                     & dltx2, phi, qavg, qsum, rad1, rad2, slope, vavg, vtot,  &
                     & wlslope
     real(R8KIND) :: BAREA, TWIDTH, WETPER
     integer, save :: j, n, np, nqcnt
     integer, allocatable, dimension(:), save :: indx
     logical, save :: smooth_water_levels
!
!*** End of declarations rewritten by SPAG
!
 
 
!    Type declarations
 
                                                                                                    ! CB 10/4/07
 
!    Allocation declarations
 
                                                                   !, OPENWRN
     allocate(y(nn), v(nn), carea(nn), topw(nn), belev(nn), q(nn), vold(nn),   &
            & yold(nn), b(nn))                                                                                     ! CB 10/4/07
     allocate(yt(nn), vt(nn), vpr(nn), ypr(nn), tarea(nn), topwt(nn), rt(nn),  &
            & indx(nn))
     allocate(al(nn, 2), daa(nn, nn))
     return
 
     entry OPEN_CHANNEL(El1, El2, Qout, Ic, Dt)
 
!    Variable initializtion
 
     b = 0.0
 
!    Variable initializtion
 
     y = 0.0
 
!    Variable initializtion
 
     v = 0.0
 
!    Variable initializtion
 
     vt = 0.0
 
!    Variable initializtion
 
     yt = 0.0
 
!    Variable initializtion
 
     rt = 0.0
 
!    Variable initializtion
 
     daa = 0.0
 
!    Variable initializtion
 
     ypr = 0.0
 
!    Variable initializtion
 
     vpr = 0.0
 
!    Variable initializtion
 
     topw = 0.0
 
!    Variable initializtion
 
     topwt = 0.0
     carea = 0.0
     tarea = 0.0
     belev(1) = upie
     belev(nc) = dnie
     phi = ASIN((upie - dnie)/clen)
     dltx = clen/(REAL(nc - 1)*0.5)
     do j = 2, nc - 1
         dltx2 = dltx*0.5
         slope = (upie - dnie)/clen
         dist = (REAL(j - 1)*dltx2)
         belev(j) = upie - slope*dist
     enddo
     bepr1 = upie + slope*dltx2
     bepr2 = dnie - slope*dltx2
     bc1 = (El1 - bepr1)*COS(phi)
     if(bc1<=0.0)bc1 = El1 - upie
     bc2 = (El2 - bepr2)*COS(phi)
     if(bc2<=0.0)bc2 = El2 - dnie
     if(.NOT.BEGIN(Ic))then
         if(WLFLAG(Ic))then
             do j = 2, nc - 1, 2
                 wlslope = ((bc1 - bc2)/(clen + dltx))*DCOS(phi)
                 dist = (REAL(j - 1)*0.5*dltx) + dltx2
                 y(j) = bc1 - wlslope*dist
                 yt(j) = y(j)
                 DTP(Ic) = Dt
             enddo
         else
             do i = 2, nc - 1, 2
                 y(i) = YS(i, Ic)
                 yt(i) = YST(i, Ic)
             enddo
         endif
     endif
     do i = 1, nc, 2
         v(i) = VS(i, Ic)
         vt(i) = VST(i, Ic)
     enddo
     if(BEGIN(Ic))then
         BEGIN(Ic) = .FALSE.
         do j = 2, nc - 1, 2
             wlslope = ((bc1 - bc2)/(clen + dltx))*DCOS(phi)
             dist = (REAL(j - 1)*0.5*dltx) + dltx2
             y(j) = bc1 - wlslope*dist
             yt(j) = y(j)
             DTP(Ic) = Dt
         enddo
         do j = 1, nc, 2
             v(j) = 0.0
             vt(j) = v(j)
         enddo
!        OPENWRN = .TRUE.
     endif
     smooth_water_levels = .FALSE.
     do n = 1, nc, 2
         if(n==nc)then
             bar1 = BAREA(bc2, dia)
             rad1 = bar1/WETPER(bc2, dia)
         else
             bar1 = BAREA(y(n + 1), dia)
             rad1 = bar1/WETPER(y(n + 1), dia)
         endif
         if(n==1)then
             bar2 = BAREA(bc1, dia)
             rad2 = bar2/WETPER(bc1, dia)
         else
             bar2 = BAREA(y(n - 1), dia)
             rad2 = bar2/WETPER(y(n - 1), dia)
         endif
         rt(n) = (rad1 + rad2)*0.5
     enddo
     do n = 2, nc - 1, 2
         tarea(n) = BAREA(y(n), dia)
         topwt(n) = TWIDTH(y(n), dia)
         carea(n) = BAREA(y(n), dia)
     enddo
 
!    Projected water levels and velocities
 
     do j = 1, nc, 2
         vpr(j) = v(j) + Dt*(v(j) - vt(j))/DTP(Ic)
     enddo
     do j = 2, nc - 1, 2
         ypr(j) = y(j) + Dt*(y(j) - yt(j))/DTP(Ic)
     enddo
 
!    Matrix setup
 
     vtot = 0.0
     do j = 1, nc, 2
         vtot = vtot + v(j)
     enddo
     vavg = vtot/(REAL(nc - 1)*0.5)
 
!    Continuity
 
     do n = 2, nc - 1, 2
         vpr(n) = (vpr(n - 1) + vpr(n + 1))*0.5D0
         v(n) = (v(n - 1) + v(n + 1))*0.5D0
         if(n/=2)daa(n, n - 2) = -THETA*(Dt/dltx)*(vpr(n)*0.5)
         daa(n, n - 1) = -THETA*(Dt/dltx)*(tarea(n)/topwt(n))
         daa(n, n) = 1.0D0
         daa(n, n + 1) = THETA*(Dt/dltx)*(tarea(n)/topwt(n))
         if(n/=nc - 1)daa(n, n + 2) = THETA*(Dt/dltx)*(vpr(n)*0.5D0)
         if(n==2)then
             b(n) = y(n) - (1.0D0 - THETA)*(Dt/dltx)*(tarea(n)/topwt(n))       &
                  & *(v(n + 1) - v(n - 1)) - (1.0D0 - THETA)*(Dt/dltx)         &
                  & *(v(n)*0.5D0)*(y(n + 2) - bc1) + THETA*(Dt/dltx)           &
                  & *(vpr(n)*0.5D0)*bc1
         elseif(n==nc - 1)then
             b(n) = y(n) - (1.0D0 - THETA)*(Dt/dltx)*(tarea(n)/topwt(n))       &
                  & *(v(n + 1) - v(n - 1)) - (1.0D0 - THETA)*(Dt/dltx)         &
                  & *(v(n)*0.5D0)*(bc2 - y(n - 2)) - THETA*(Dt/dltx)           &
                  & *(vpr(n)*0.5D0)*bc2
         else
             b(n) = y(n) - (1.0D0 - THETA)*(Dt/dltx)*(tarea(n)/topwt(n))       &
                  & *(v(n + 1) - v(n - 1)) - (1.0D0 - THETA)*(Dt/dltx)         &
                  & *(v(n)*0.5D0)*(y(n + 2) - y(n - 2))
         endif
     enddo
     if(vavg>0.0 .OR. (vavg==0.0 .AND. El1>El2))then
 
!**      Momentum
 
         do n = 1, nc, 2
             if(n/=1)then
                 daa(n, n - 2) = -THETA*(Dt/dltx)*vpr(n)
                 daa(n, n - 1) = -THETA*(Dt/dltx)*g*DCOS(phi)
             endif
             daa(n, n) = 1.0 + THETA*Dt*g*(fman**2)*DABS(vpr(n))               &
                       & /(rt(n)**(4.0/3.0)) + THETA*(Dt/dltx)*vpr(n)          &
                       & + THETA*(closs*0.5D0)*(Dt/clen)*DABS(vpr(n))
             if(n/=nc)daa(n, n + 1) = THETA*(Dt/dltx)*g*DCOS(phi)
             if(n==1)then
                 b(n) = v(n) - (1.0D0 - THETA)*(Dt/dltx)*g*(y(n + 1) - bc1)    &
                      & *DCOS(phi) - (1.0D0 - THETA)*v(n)*(Dt/dltx)*v(n)       &
                      & - (1.0D0 - THETA)*Dt*g*(fman**2)/(rt(n)**(4.0/3.0))    &
                      & *v(n)*DABS(v(n)) + Dt*g*DSIN(phi) - (1.0D0 - THETA)    &
                      & *(Dt/clen)*(closs*0.5D0)*v(n)*DABS(v(n))               &
                      & + THETA*(Dt/dltx)*g*DCOS(phi)*bc1
             elseif(n==nc)then
                 b(n) = v(n) - (1.0D0 - THETA)*(Dt/dltx)*g*(bc2 - y(n - 1))    &
                      & *DCOS(phi) - (1.0D0 - THETA)*v(n)*(Dt/dltx)            &
                      & *(v(n) - v(n - 2)) - (1.0D0 - THETA)*Dt*g*(fman**2)    &
                      & /(rt(n)**(4.0/3.0))*v(n)*DABS(v(n)) + Dt*g*DSIN(phi)   &
                      & - (1.0D0 - THETA)*(Dt/clen)*(closs*0.5D0)*v(n)         &
                      & *DABS(v(n)) - THETA*(Dt/dltx)*g*DCOS(phi)*bc2
             else
                 b(n) = v(n) - (1.0D0 - THETA)*(Dt/dltx)                       &
                      & *g*(y(n + 1) - y(n - 1))*COS(phi) - (1.0D0 - THETA)    &
                      & *v(n)*(Dt/dltx)*(v(n) - v(n - 2)) - (1.0D0 - THETA)    &
                      & *Dt*g*(fman**2)/(rt(n)**(4.0/3.0))*v(n)*DABS(v(n))     &
                      & + Dt*g*DSIN(phi) - (1.0D0 - THETA)*(Dt/clen)           &
                      & *(closs*0.5D0)*v(n)*DABS(v(n))
             endif
         enddo
     else
         do n = 1, nc, 2
             if(n/=nc)then
                 daa(n, n + 2) = THETA*(Dt/dltx)*vpr(n)
                 daa(n, n + 1) = THETA*(Dt/dltx)*g*DCOS(phi)
             endif
             daa(n, n) = 1.0 + THETA*Dt*g*(fman**2)*DABS(vpr(n))               &
                       & /(rt(n)**(4.0/3.0)) - THETA*(Dt/dltx)*vpr(n)          &
                       & + THETA*(closs*0.5D0)*(Dt/clen)*DABS(vpr(n))
             if(n/=1)daa(n, n - 1) = -THETA*(Dt/dltx)*g*DCOS(phi)
             if(n==nc)then
                 b(n) = v(n) - (1.0D0 - THETA)*(Dt/dltx)*g*(bc2 - y(n - 1))    &
                      & *DCOS(phi) - (1.0 - THETA)*v(n)*(Dt/dltx)*( - v(n))    &
                      & - (1.0D0 - THETA)*Dt*g*(fman**2)/(rt(n)**(4.0/3.0))    &
                      & *v(n)*DABS(v(n)) + Dt*g*DSIN(phi) - (1.0 - THETA)      &
                      & *(Dt/clen)*(closs*0.5)*v(n)*DABS(v(n))                 &
                      & - THETA*(Dt/dltx)*g*DCOS(phi)*bc2
             elseif(n==1)then
                 b(n) = v(n) - (1.0D0 - THETA)*(Dt/dltx)*g*(y(n + 1) - bc1)    &
                      & *DCOS(phi) - (1.0 - THETA)*v(n)*(Dt/dltx)              &
                      & *(v(n + 2) - v(n)) - (1.0D0 - THETA)*Dt*g*(fman**2)    &
                      & /(rt(n)**(4.0/3.0))*v(n)*ABS(v(n)) + Dt*g*SIN(phi)     &
                      & - (1.0 - THETA)*(Dt/clen)*(closs*0.5D0)*v(n)*DABS(v(n))&
                      & + THETA*(Dt/dltx)*g*DCOS(phi)*bc1
             else
                 b(n) = v(n) - (1.0D0 - THETA)*(Dt/dltx)                       &
                      & *g*(y(n + 1) - y(n - 1))*DCOS(phi) - (1.0D0 - THETA)   &
                      & *v(n)*(Dt/dltx)*(v(n + 2) - v(n)) - (1.0D0 - THETA)    &
                      & *Dt*g*(fman**2)/(rt(n)**(4.0/3.0))*v(n)*DABS(v(n))     &
                      & + Dt*g*DSIN(phi) - (1.0D0 - THETA)*(Dt/clen)           &
                      & *(closs*0.5D0)*v(n)*DABS(v(n))
             endif
         enddo
     endif
     np = nn
     call LUDCMP(daa, nc, np, indx, d)
     call LUBKSB(daa, nc, np, indx, b)
     do i = 2, nc - 1, 2
         yold(i) = y(i)
         YST(i, Ic) = y(i)
     enddo
     do i = 2, nc - 1, 2
         y(i) = b(i)
     enddo
 
!    Smooth water levels
 
     do i = 2, nc - 1, 2
 !     IF (OPENWRN) THEN
 !       OPEN (391,FILE='culvert.wrn',STATUS='unknown')
 !       OPENWRN = .FALSE.
 !     END IF
         if(y(i)<=0.0)smooth_water_levels = .TRUE.
     enddo
     if(smooth_water_levels)then
         do j = 2, nc - 1, 2
             wlslope = ((bc1 - bc2)/(clen + dltx))*DCOS(phi)
             dist = (REAL(j - 1)*0.5D0*dltx) + dltx2
             y(j) = bc1 - wlslope*dist
         enddo
  !  WRITE (391,10010) IC, JDAY
         smooth_water_levels = .FALSE.
     endif
 
!    Flows
 
     nqcnt = 0
     qsum = 0.0
     do i = 1, nc, 2
         vold(i) = v(i)
         VST(i, Ic) = v(i)
         v(i) = b(i)
         if(i==nc)then
             bar1 = BAREA(bc2, dia)
         else
             bar1 = BAREA(y(i + 1), dia)
         endif
         if(i==1)then
             bar2 = BAREA(bc1, dia)
         else
             bar2 = BAREA(y(i - 1), dia)
         endif
         carea(i) = (bar1 + bar2)*0.5D0
         q(i) = v(i)*carea(i)
         nqcnt = nqcnt + 1
         qsum = qsum + q(i)
     enddo
     qavg = qsum/REAL(nqcnt)
     do i = 2, nc - 1, 2
         YS(i, Ic) = y(i)
     enddo
     VMAX(Ic) = 0.0
     do i = 1, nc, 2
         VS(i, Ic) = v(i)
         VMAX(Ic) = MAX(ABS(v(i)), VMAX(Ic))
     enddo
     DTP(Ic) = Dt
     Qout = qavg
     QOLD(Ic) = Qout
     WLFLAG(Ic) = .FALSE.
9001  format('water levels for culvert ', i3, ' on Julian Day ', f10.3,        &
            &' are <= 0 - predictions have been smoothed')
     return
     entry DEALLOCATE_OPEN_CHANNEL
     deallocate(y, v, carea, topw, belev, q, vold, yold, b, yt, vt, vpr, ypr,  &
              & tarea, topwt, rt, indx, al, daa)                                                                  ! CB 10/4/07
     end subroutine OPEN_CHANNEL_INITIALIZE
