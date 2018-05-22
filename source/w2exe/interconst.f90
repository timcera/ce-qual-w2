!*==interconst.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************
!***********************************************************************
!*                                                                    **
!* I N T E R P O L A T   A L G O R I T H M - C O N S T I T U E N T S  **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!    1                   2                                        3    *
!     O-------------------O----------------------------------------O   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!     |         +         |                    +                   |   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!    4|                  5|                                       6|   *
!     O-------------------O----------------------------------------O   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!     |         +         |                    +                   |   *
!     |                   |                                        |   *
!    7|                  8|                                       9|   *
!     O-------------------O----------------------------------------O   *
!                                                                      *
!     + = Constituent Information Given in W2                          *
!     O = Node where Information is Needed for TecPlot Purposes        *
!                                                                      *
!***********************************************************************
!  This Subroutine Calls No Other Subroutines
 
     subroutine INTERCONST
 
 
 
 
 
 
     use FISHY
     use GLOBAL
     use GEOMC
     use SCREENC
     use KINETIC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: ad, at, bd, bt, ct, cv, dv, dw, dx, dy, h1do, h1t, h2do, h2t,     &
           & h3do, h3t, h4do, h4t, lfti, rgti, xrange, yrange
     integer :: ii, k, kk, nn, totnn, wb
     real, save :: topk
!
!*** End of declarations rewritten by SPAG
!
 
 
 
 
     nn = 0                                                        ! NN = Global Node # (UpperLeft Corner of cell(K,I))
     do wb = 1, nwb                                                ! WB = Water Body # / NWB = Total # of Water Bodies
         do jb = BS(wb), BE(wb)                                    ! JB = Branch #
             do k = 2, kmx                                         ! K = Layer #
                 do i = US(jb), (DS(jb) + 1)                   ! I = Segment #  changed from us to cus sw 5/2015
                     if((k - 1<=KB(i)) .OR. (k - 1<=KB(i - 1)))then
                                                                   ! All Nodes at or above Water Body Bottom
                         nn = nn + 1                               ! PRE-CHECK CALCULATION
                         kk = NDINFO(nn, 3)                        ! PRE-CHECK CALCULATION
                         ii = NDINFO(nn, 4)                        ! PRE-CHECK CALCULATION
                         if((i/=ii) .OR. (k/=kk))then              ! CHECK
 
                             write(datadebugfn, *)                             &
                  &'ERROR: Nodes not matching up for Constituent Interpolation'
                             write(datadebugfn, *)'JDAY=', jday
                             write(datadebugfn, *)'WB=', wb, ' JR=', jb
                             write(datadebugfn, *)'NN=', nn, ' kk=', kk,       &
                                 & ' ii=', ii, ' k=', k, ' i=', i
                             stop
                         endif
 
!                        LINEAR INTERPOLATION for Boundary Nodes; those
!                        analogous to Node 5 use Bilinear Splines below
                         if(k<KTWB(wb))then                              !Above Water Surface
                             WQFIELD(k, i, 1) = 0.0
                             WQFIELD(k, i, 2) = 0.0
                         elseif((k==KTWB(wb)) .AND. (i==US(jb)))then   !Analogous to Node 1 above
                             WQFIELD(k, i, 1) = T2(k, i)
                             WQFIELD(k, i, 2) = O2(k, i)
                         elseif((k==KTWB(wb)) .AND. (i==DS(jb) + 1))then
                                                                       !Analogous to Node 3 above
                             WQFIELD(k, i, 1) = T2(k, i - 1)
                             WQFIELD(k, i, 2) = O2(k, i - 1)
                         elseif(k==KTWB(wb))then                         !Analogous to Node 2 above
                             xrange = .5*DLX(i - 1) + .5*DLX(i)
                             WQFIELD(k, i, 1)                                  &
                               & = (T2(k, i - 1)*(xrange - .5*DLX(i - 1))      &
                               & + T2(k, i)*(xrange - .5*DLX(i)))/xrange
                             WQFIELD(k, i, 2)                                  &
                               & = (O2(k, i - 1)*(xrange - .5*DLX(i - 1))      &
                               & + O2(k, i)*(xrange - .5*DLX(i)))/xrange
                         elseif((k==KB(i) + 1) .AND. (i==US(jb)))then !Analogous to Node 7 above
                             WQFIELD(k, i, 1) = T2(k - 1, i)
                             WQFIELD(k, i, 2) = O2(k - 1, i)
                         elseif((k==KB(i) + 1) .AND. (i==DS(jb) + 1))then
                                                                      !Analogous to Node 9 above
                             WQFIELD(k, i, 1) = T2(k - 1, i - 1)
                             WQFIELD(k, i, 2) = O2(k - 1, i - 1)
                         elseif(k==KB(i) + 1)then                       !Analogous to Node 8 above
                             xrange = .5*DLX(i - 1) + .5*DLX(i)
                             WQFIELD(k, i, 1)                                  &
                               & = (T2(k - 1, i - 1)*(xrange - .5*DLX(i - 1))  &
                               & + T2(k - 1, i)*(xrange - .5*DLX(i)))/xrange
                             WQFIELD(k, i, 2)                                  &
                               & = (O2(k - 1, i - 1)*(xrange - .5*DLX(i - 1))  &
                               & + O2(k - 1, i)*(xrange - .5*DLX(i)))/xrange
                         elseif(i==US(jb))then                        !Analogous to Node 4 above
                             yrange = .5*H(k - 1, wb) + .5*H(k, wb)
                             WQFIELD(k, i, 1)                                  &
                               & = (T2(k - 1, i)*(yrange - .5*H(k - 1, wb))    &
                               & + T2(k, i)*(yrange - .5*H(k, wb)))/yrange
                             WQFIELD(k, i, 2)                                  &
                               & = (O2(k - 1, i)*(yrange - .5*H(k - 1, wb))    &
                               & + O2(k, i)*(yrange - .5*H(k, wb)))/yrange
                         elseif(i==DS(jb) + 1)then                    !Analogous to Node 6 above
                             yrange = .5*H(k - 1, wb) + .5*H(k, wb)
                             WQFIELD(k, i, 1)                                  &
                               & = (T2(k - 1, i - 1)*(yrange - .5*H(k - 1, wb))&
                               & + T2(k, i - 1)*(yrange - .5*H(k, wb)))/yrange
                             WQFIELD(k, i, 2)                                  &
                               & = (O2(k - 1, i - 1)*(yrange - .5*H(k - 1, wb))&
                               & + O2(k, i - 1)*(yrange - .5*H(k, wb)))/yrange
                         else                                      ! Analogous to Node 5 above (i.e., in the middle)
!                            Begin Bilinear Spline Interpolation Algorithm
!                            (Algorithm obtained from "TWO DIMENSIONAL SPLINE
!                            INTERPOLATION ALGORITHMS - Author: Helmuth Spath")
                             dv = 0.0 - ( - 0.5*DLX(i - 1))
                                                     ! 0.0 = X-coord of the desired interpolated location
                             dw = 0.0 - ( - 0.5*H(k - 1, wb))
                                                        ! 0.0 = Y-coord of the desired interpolated location
                             dx = 0.5*DLX(i) - ( - 0.5*DLX(i - 1))
                             dy = 0.5*H(k, wb) - ( - 0.5*H(k - 1, wb))
                             h1t = T2(k - 1, i - 1)  ! H1T,H2T,H3T,H4T = Value of Temperature at 4 cell centers as
                             h2t = T2(k - 1, i)      !   established in W2
                             h3t = T2(k, i - 1)
                             h4t = T2(k, i)
                             h1do = O2(k - 1, i - 1) ! H1DO,H2DO,H3DO,H4DO = Value of Diss. Oxyg. at 4 cell centers as
                             h2do = O2(k - 1, i)     !   established in W2
                             h3do = O2(k, i - 1)
                             h4do = O2(k, i)
                             WQFIELD(k, i, 1) = h1t + (h2t - h1t)              &
                               & /dx*dv + (h3t - h1t)                          &
                               & /dy*dw + (h4t - h3t - h2t + h1t)/(dx*dy)*dv*dw
                                                                  !                ! Temperature at grid node
                             WQFIELD(k, i, 2) = h1do + (h2do - h1do)           &
                               & /dx*dv + (h3do - h1do)                        &
                               & /dy*dw + (h4do - h3do - h2do + h1do)/(dx*dy)  &
                               & *dv*dw                               !            ! Dissolved Oxygen at grid node
!                            End Bilinear Spline Interpolation Algorithm
                         endif
 
 
                     endif
                 enddo
             enddo
         enddo
     enddo
 
     end subroutine INTERCONST
