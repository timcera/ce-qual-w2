!*==part_transport.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine PART_TRANSPORT
! Compute passive particle transport including sedimentation
 
     use FISHY
     use TRANS
     use GLOBAL
     use GEOMC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: costheta, dispx, dispz, dx1, dx2, dxavg, dz1, dz2, dzavg, r1, r2, &
           & rx, rz, sintheta, sk, vel, wpart, xarea
     real, save :: dxmin, dzmax1, dzmin1
!
!*** End of declarations rewritten by SPAG
!
 
     data dxmin, dzmin1, dzmax1/0.0005, 0.0001, 0.5/
                                                   ! SW 2/01/01   in units of m2/s
 
     fxloc = fxloc + (FXVEL(5))*(nfsfreq)*24.*3600.         ! New updated Part X-Location
     fyloc = fyloc + (fyvel)*(nfsfreq)*24.*3600.         ! New updated Part Y-Location
     fzloc = fzloc + (FZVEL(5) + SEDVEL(fn))*(nfsfreq)*24.*3600.         ! New updated Part Z-Location
 
 
     vel = SQRT(FXVEL(5)*FXVEL(5) + FZVEL(5)*FZVEL(5))
 
    ! Compute "Dispersion Processes"
     if(vel>1.0E-09)then
         costheta = FXVEL(5)/vel
         sintheta = FZVEL(5)/vel
!        Dispx=alphax*abs(fxvel(5))
!        Dispz=alphaz*abs(fzvel(5))
 
!        average to segment centers - should be interpolated though   ! SW
!        2/01/01
         dx1 = DX(fkmp, fimp)
         dx2 = DX(fkmp, fimp - 1)
         if(dx1<=0.0)dx1 = dxmin
         if(dx2<=0.0)dx2 = dxmin
         dxavg = 0.5*(dx1 + dx2)
 
         dz1 = DZ(fkmp, fimp)
         if(fkmp<=ktwbf)then
             dz2 = DZ(fkmp, fimp)
         else
             dz2 = DZ(fkmp - 1, fimp)
         endif
!        Note during density inversions DZ is set to DZMAX1=1000. This leads
!        to incredible variations in RZ - constrain to DZ=0.5
         if(dz1>10.)dz1 = dzmax1
         if(dz2>10.)dz2 = dzmax1
 
         if(dz1<=0.0)dz1 = dzmin1
         if(dz2<=0.0)dz2 = dzmin1
         dzavg = 0.5*(dz1 + dz2)
 
!if(dzavg.gt.10)then
!write(DIAGFN,*)'DZAVG>10 JDAY:',jday,'dz1,dz2,fkmp,fimp,dz(fkmp,fimp),dz(fkmp-1,fimp)'
!write(DIAGFN,*)dz1,dz2,fkmp,fimp,dz(fkmp,fimp),dz(fkmp-1,fimp)
!end     if
!if(dxavg.gt.20)then
!write(DIAGFN,*)'DXAVG>10 JDAY:',jday,'dx1,dx2,fkmp,fimp,dx(fkmp,fimp),dx(fkmp-1,fimp)'
!write(DIAGFN,*)dx1,dx2,fkmp,fimp,dx(fkmp,fimp),dx(fkmp,fimp-1)
!end     if
 
         if(dxtheory==' ON')then
             dispx = 5.84E-4*DLX(fimp)**1.1
             dispz = dzavg
         else
             dispx = dxavg
                         !*alphax            ! Note these should be interpolated rather than using nearest cell #
             dispz = dzavg
                         !*alphaz
         endif
 
!        COMPUTE DZ and DX by interpolating
 
         call RANDOM(seed, r1)
         call RANDOM(seed, r2)
 
!        r1=gasdev(seed)
!        r2=gasdev(seed)
 
         r1 = (r1 - 0.5)*2.
         r2 = (r2 - 0.5)*2.
 
         rx = SQRT(6.0*dispx*nfsfreq*86400.)*(r1*costheta - r2*sintheta)
                                                                        ! DX is m2/s    nfsfreq is days
         rz = SQRT(6.0*dispz*nfsfreq*86400.)*(r1*sintheta + r2*costheta)
 
        !write(9500,*)(r1*costheta-r2*sintheta),(r1*sintheta+r2*costheta)
 
!        constrain random component to segment length and cell layer height
         if(ABS(rz)>H(ktwbf, fjr))rz = SIGN(1.0, rz)*H(ktwbf, fjr)
                                                                  !.gt.h(fkmp,KTWBF))rz=sign(1.0,rz)*h(fkmp,KTWBF)    SW 7/1/2017
         if(ABS(rx)>DLX(fimp))rx = SIGN(1.0, rx)*DLX(fimp)
 
!if(rz.gt.2..or.rx.gt.1000.)then
!write(DIAGFN,*)'rx,rz,dispz,nfsfreq,r1,sintheta,r2,costheta,dzavg,alphaz,dxavg,alphax,dz1,dz2,dx1,dx2,vel,fxvel(5),fzvel(5)'
!write(DIAGFN,*)rx,rz,dispz,nfsfreq,r1,sintheta,r2,costheta,dzavg,alphaz,dxavg,alphax,dz1,dz2,dx1,dx2,vel,fxvel(5),fzvel(5)
!end     if
 
        !RX=SQRT(24.0*Dispx*nfsfreq)*(r1-0.5)
        !RZ=SQRT(24.0*Dispz*nfsfreq)*(r2-0.5)
 
     else
         rx = 0.0
         rz = 0.0
 
     endif
 
     fxloc = fxloc + rx      ! New updated Part X-Location
     fyloc = fyloc + rx      ! New updated Part Y-Location
     fzloc = fzloc + rz       ! New updated Part Z-Location
 
!if(rx.eq.0.0)then
!write(DIAGFN,*)'RX=0.0:rx,rz,dispx,r1,sintheta,r2,costheta,vel,fxvel(5)'  ! debug
!write(DIAGFN,*)rx,rz,dispx,r1,sintheta,r2,costheta,vel,fxvel(5)  ! debug
!end if
 
 
     end subroutine PART_TRANSPORT
