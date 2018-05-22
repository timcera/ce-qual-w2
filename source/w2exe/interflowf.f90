!*==interflowf.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************
!***********************************************************************
!*                                                                    **
!*  I N T E R P O L A T    A L G O R I T H M  -  F L O W   F I E L D  **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!    1                   2                                        3    *
!     O--------[VV]-------O------------------[VV]------------------O   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!    [HV]                [HV]                                     [HV] *
!     |                   |                                        |   *
!     |                   |                                        |   *
!    4|                  5|                                       6|   *
!     O--------[VV]-------O------------------[VV]------------------O   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!    [HV]                [HV]                                     [HV] *
!     |                   |                                        |   *
!    7|                  8|                                       9|   *
!     O--------[VV]-------O------------------[VV]------------------O   *
!                                                                      *
!     [HV] = Horizontal Velocity Given                                 *
!     [VV] = Vertical Velocity Given                                   *
!      O   = Node where Information is Needed for TecPlot Purposes     *
!                                                                      *
!***********************************************************************
 
!  This Subroutine Calls No Other Subroutines
 
 
     subroutine INTERFLOWF
 
 
 
 
     use FISHY
     use GLOBAL
     use GEOMC
     use SCREENC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
     integer, parameter :: N = 3
!
! Local variables
!
     real :: ad, at, bd, bt, ct, cv, lfti, rgti, xrange, xv, yrange
     real, dimension(N + 1) :: dlxho, hvert, uvert, whorz
     real, dimension(N + 1, N + 1) :: fddu, fddw
     integer :: ii, j, jjz, k, kk, nn, totnn, wb
     real, save :: lastjday, topk
!
!*** End of declarations rewritten by SPAG
!
 
 
 
 
 
 
                                                ! The order of the Newton Interpolating Polynomial
     nn = 0                                                        ! NN = Global Node # (UpperLeft Corner of cell(K,I))
     do wb = 1, nwb                                                ! WB = Water Body # / NWB = Total # of Water Bodies
         do jb = BS(wb), BE(wb)                                  ! JB = Branch #
             do k = 2, kmx                                         ! K = Layer #
                 do i = US(jb), (DS(jb) + 1)                   ! I = Segment #  SW 5/2015
                     if((k - 1<=KB(i)) .OR. (k - 1<=KB(i - 1)))then
                                                                   ! All Nodes at or above Water Body Bottom
                         nn = nn + 1                               ! PRE-CHECK CALCULATION
                         kk = NDINFO(nn, 3)                        ! PRE-CHECK CALCULATION
                         ii = NDINFO(nn, 4)                        ! PRE-CHECK CALCULATION
                         if((i/=ii) .OR. (k/=kk))then              ! CHECK
                             write(datadebugfn, *)                             &
                    &'ERROR: Nodes not matching up for Flowfield Interpolation'
                             write(datadebugfn, *)'WB=', wb, ' JB=', jb
                             write(datadebugfn, *)'NN=', nn, ' kk=', kk,       &
                                 & ' ii=', ii, ' k=', k, ' i=', i
 
                             stop
                         endif
 
!LINEAR                  INTERPOLATION of the Flow Field (i.e. Velocities)
                         if(.NOT.(linear))then
                                              ! LINEAR INTERPOLATION
 
!3rd                         ORDER NEWTON INTERPOLATING POLYNOMIAL used for
!                            interpolation of the Flow Field (i.e. Velocities)
!                            1
! O-----[VV]----O------[VV]------O-----[VV]----O------[VV]------O |            
!                            |             |                |             |   
!                            | [HV]          [HV]             [HV]         
!                            [HV]             [HV] |             |            
!                            |             |                | |             | 
!                            2|             |                |
! O-----[VV]----O------[VV]------O-----[VV]----O------[VV]------O |            
!                            |             |                |             |   
!                            | [HV]          [HV]             [HV]         
!                            [HV]             [HV] |             |            
!                            |             |                | 1|            2|
!                            |            3|               4|
! O-----[VV]----O------[VV]------X-----[VV]----O------[VV]------O |            
!                            |             |                |             |   
!                            | [HV]          [HV]             [HV]  (K,I)  
!                            [HV]             [HV] |             |            
!                            |             |                | |             | 
!                            3|             |                |
! O-----[VV]----O------[VV]------O-----[VV]----O------[VV]------O |            
!                            |             |                |             |   
!                            | [HV]          [HV]             [HV]         
!                            [HV]             [HV] |             |            
!                            |             |                | |             | 
!                            4|             |                |
! O-----[VV]----O------[VV]------O-----[VV]----O------[VV]------O
                             xv = 0.0               ! Location (X ABOVE) of Interpolated Velocity (both U & W)
                             if(k - 1<KTWB(wb))then  !Both Top Nodes are above Water Surface
                                 uvert(1) = U(k, i - 1)
                                                    ! This Preserves the Horizontal Velocity at the Water Surface
                                 uvert(2) = U(k, i - 1)
                                                    ! This Preserves the Horizontal Velocity at the Water Surface
                                 hvert(1) = -H(k, wb)/10000.
                                                        ! Arbitrary value
                                 hvert(2) = -H(k, wb)/10001.
                                                        ! Arbitrary value
                             elseif(k - 2<KTWB(wb))then
                                                     !Only Top-most Node above Water Surface
                                 uvert(1) = U(k - 1, i - 1)
                                                    ! This Preserves the Horizontal Velocity at the Water Surface
                                 uvert(2) = U(k - 1, i - 1)
                                 hvert(1) = -H(k - 1, wb)
                                                       ! Arbitrary value
                                 hvert(2) = -0.5*H(k - 1, wb)
                             else                   !Both Top Nodes exist
                                 uvert(1) = U(k - 2, i - 1)
                                 uvert(2) = U(k - 1, i - 1)
                                 hvert(1) = -H(k - 1, wb) - 0.5*H(k - 2, wb)
                                 hvert(2) = -0.5*H(k - 1, wb)
                             endif
                             if(k>KB(i))then        !Both Bottom Nodes are below Water Body Bottom
                                 uvert(3) = 0.0       ! Helps satisfy no-slip condition on bottom (i.e., U=0)
                                 uvert(4) = 0.0       ! Helps satisfy no-slip condition on bottom (i.e., U=0)
                                 hvert(3) = H(k - 1, wb)/10001.
                                                        ! Arbitrary value
                                 hvert(4) = H(k - 1, wb)/10000.
                                                        ! Arbitrary value
                             elseif(k + 1>KB(i))then
                                                    !Only Bottom-most Node below Water Body Bottom
                                 uvert(3) = U(k, i - 1)
                                 uvert(4) = 0.0       ! Helps satisfy no-slip condition on bottom (i.e., U=0)
                                 hvert(3) = 0.5*H(k, wb)
                                 hvert(4) = H(k, wb)   ! Helps satisfy no-slip condition on bottom (i.e., U=0)
                             else                   !Both Bottom Nodes exist
                                 uvert(3) = U(k, i - 1)
                                 uvert(4) = U(k + 1, i - 1)
                                 hvert(3) = 0.5*H(k, wb)
                                 hvert(4) = H(k, wb) + 0.5*H(k + 1, wb)
                             endif
                             if(i - 1<US(jb))then !Both Left Nodes are Upstream of Water Body
                                 whorz(1) = 0.0       ! Helps satisfy no-slip condition at Upstream Boundary
                                 whorz(2) = 0.0       ! Helps satisfy no-slip condition at Upstream Boundary
                                 dlxho(1) = -DLX(i)/10000.
                                                     ! Arbitrary value
                                 dlxho(2) = -DLX(i)/10001.
                                                     ! Arbitrary value
                             elseif(i - 2<US(jb))then
                                                  !Only Left-most Node Upstream of Water Body
                                 whorz(1) = 0       ! Helps satisfy no-slip condition at Upstream Boundary
                                 whorz(2) = W(k - 1, i - 1)
                                 dlxho(1) = -DLX(i - 1)
                                                    ! Arbitrary value
                                 dlxho(2) = -0.5*DLX(i - 1)
                             else                   !Both Upstream Nodes exist
                                 whorz(1) = W(k - 1, i - 2)
                                 whorz(2) = W(k - 1, i - 1)
                                 dlxho(1) = -DLX(i - 1) - 0.5*DLX(i - 2)
                                 dlxho(2) = -0.5*DLX(i - 1)
                             endif
                             if(i>DS(jb))then     !Both Right Nodes are Downstream of Water Body
                                 whorz(3) = 0       ! Helps satisfy no-slip condition at Downstream Boundary
                                 whorz(4) = 0       ! Helps satisfy no-slip condition at Downstream Boundary
                                 dlxho(3) = DLX(i - 1)/10001
                                                    ! Arbitrary value
                                 dlxho(4) = DLX(i - 1)/10000
                                                    ! Arbitrary value
                             elseif(i + 1>DS(jb))then
                                                  !Only Right-most Node Downstream of Water Body
                                 whorz(3) = W(k - 1, i)
                                 whorz(4) = 0       ! Helps satisfy no-slip condition at Downstream Boundary
                                 dlxho(3) = 0.5*DLX(i)
                                 dlxho(4) = DLX(i)  ! Arbitrary value
                             else                   !Both Downstream Nodes exist
                                 whorz(3) = W(k - 1, i)
                                 whorz(4) = W(k - 1, i + 1)
                                 dlxho(3) = 0.5*DLX(i)
                                 dlxho(4) = DLX(i) + 0.5*DLX(i + 1)
                             endif
!                            Begin Newton Interpolating Polynomial
!                            Algorithm(Algorithm obtained from "NUMERICAL
!                            METHODS for Engineers - Authors: Steven Chapra &
!                            Raymond Canale")
                             do jjz = 1, N + 1
                                 fddu(jjz, 1) = uvert(jjz)
                                 fddw(jjz, 1) = whorz(jjz)
                             enddo
                             do j = 2, N + 1
                                 do jjz = 1, (N + 1 - j) + 1
                                     fddu(jjz, j)                              &
                                       & = (fddu(jjz + 1, j - 1) - fddu(jjz,   &
                                       & j - 1))                               &
                                       & /(hvert(jjz + j - 1) - hvert(jjz))
                                     fddw(jjz, j)                              &
                                       & = (fddw(jjz + 1, j - 1) - fddw(jjz,   &
                                       & j - 1))                               &
                                       & /(dlxho(jjz + j - 1) - dlxho(jjz))
                                 enddo
                             enddo
                             if(k<KTWB(wb))then                     !Above Water Surface
                                 FLOWFIELD(k, i, 1) = 0.0          ! Horiz. Vel.
                                 FLOWFIELD(k, i, 2) = 0.0          ! Vert.  Vel. - Original Vert Vel from W2 are in millimeters
                             else
                                 FLOWFIELD(k, i, 1) = fddu(1, 1) + fddu(1, 2)  &
                                   & *(xv - hvert(1)) + fddu(1, 3)             &
                                   & *(xv - hvert(1))*(xv - hvert(2))          &
                                   & + fddu(1, 4)*(xv - hvert(1))              &
                                   & *(xv - hvert(2))*(xv - hvert(3))
                                 FLOWFIELD(k, i, 2)                            &
                                   & = (fddw(1, 1) + fddw(1, 2)*(xv - dlxho(1))&
                                   & + fddw(1, 3)*(xv - dlxho(1))              &
                                   & *(xv - dlxho(2)) + fddw(1, 4)             &
                                   & *(xv - dlxho(1))*(xv - dlxho(2))          &
                                   & *(xv - dlxho(3)))             ! Original Vertical Vel from W2 are in millimeters
                             endif
                         elseif(k<KTWB(wb))then                     !Above Water Surface
                             FLOWFIELD(k, i, 1) = 0.0              !   Horiz. Vel.
                             FLOWFIELD(k, i, 2) = 0.0              !   Vert.  Vel.
                         elseif(k==KTWB(wb))then                    !At Water Surface
                             FLOWFIELD(k, i, 1) = U(k, i - 1)      !   Horiz. Vel.
                             FLOWFIELD(k, i, 2) = 0.0              !   Vert.  Vel.
                         elseif(k==KB(i) + 1)then                  !Water Body Bottom
                             FLOWFIELD(k, i, 1) = U(k - 1, i - 1)
                             FLOWFIELD(k, i, 2) = 0.0
                         elseif(i==US(jb))then                   !Upstream-most Nodes
                             yrange = .5*H(k - 1, wb) + .5*H(k, wb)
                             FLOWFIELD(k, i, 1)                                &
                               & = (U(k - 1, i - 1)*(yrange - .5*H(k - 1, wb)) &
                               & + U(k, i - 1)*(yrange - .5*H(k, wb)))/yrange
                             FLOWFIELD(k, i, 2) = W(k - 1, i)
                         elseif(i==DS(jb) + 1)then               !Downstream-most Nodes
                             yrange = .5*H(k - 1, wb) + .5*H(k, wb)
                             FLOWFIELD(k, i, 1)                                &
                               & = (U(k - 1, i - 1)*(yrange - .5*H(k - 1, wb)) &
                               & + U(k, i - 1)*(yrange - .5*H(k, wb)))/yrange
                             FLOWFIELD(k, i, 2) = W(k - 1, i - 1)
                         else                                      !Analogous to Node 5 above (i.e., in the middle)
                             yrange = .5*H(k - 1, wb) + .5*H(k, wb)
                             FLOWFIELD(k, i, 1)                                &
                               & = (U(k - 1, i - 1)*(yrange - .5*H(k - 1, wb)) &
                               & + U(k, i - 1)*(yrange - .5*H(k, wb)))/yrange
                             xrange = .5*DLX(i - 1) + .5*DLX(i)
                             FLOWFIELD(k, i, 2)                                &
                               & = ((W(k - 1, i - 1)*(xrange - .5*DLX(i-1))    &
                               & + W(k - 1, i)*(xrange - .5*DLX(i)))/xrange)
!                            End Newton Interpolating Polynomial Algorithm
                         endif
 
 
!                        Begin Calculation of Accelerations ==> Time
!                        Derivative of Velocity (Backward Differencing)
                         if(nit/=0)then
                             if((jday - lastjday)==0)then
                                 FLOWFIELD(k, i, 3) = 0.0
                                 FLOWFIELD(k, i, 4) = 0.0
                             else
                                 FLOWFIELD(k, i, 3)                            &
                                   & = (FLOWFIELD(k, i, 1) - LASTFIELD(k, i, 1)&
                                   & )/((jday - lastjday)*86400.)                                      ! 24.*3600.
                                 FLOWFIELD(k, i, 4)                            &
                                   & = (FLOWFIELD(k, i, 2) - LASTFIELD(k, i, 2)&
                                   & )/((jday - lastjday)*86400.)
                             endif
                             LASTFIELD(k, i, 1) = FLOWFIELD(k, i, 1)
                             LASTFIELD(k, i, 2) = FLOWFIELD(k, i, 2)
                             lastjday = jday
                         else                     ! No Time Derivative if Time=0 (i.e., NIT=0)
                             FLOWFIELD(k, i, 3) = 0
                             FLOWFIELD(k, i, 4) = 0
                             LASTFIELD(k, i, 1) = FLOWFIELD(k, i, 1)
                             LASTFIELD(k, i, 2) = FLOWFIELD(k, i, 2)
                             lastjday = jday
                         endif
!                        End Calculation of Accelerations
 
                     endif
                 enddo
             enddo
         enddo
     enddo
 
     end subroutine INTERFLOWF
