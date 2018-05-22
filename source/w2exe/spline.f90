!*==spline.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************
!***********************************************************************
!*                                                                    **
!*     B I L I N E A R   S P L I N E   I N T E R P O L A T I O N      **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!                                                                      *
!   Locations:1,2,3,4,5            4 (Above Fish)                      *
!    (LOCATE)                      ^                                   *
!                               UZSENSPH                               *
!                                  |                                   *
!                                _____Fish                             *
!  (Behind Fish)             |\_/     \_     (Forward of Fish)         *
!       3  <--- BXSENSPH --- | _   5   _> ----- FXSENSPH ----->  1     *
!                            |/ \_____/                                *
!                                  |                                   *
!                               DZSENSPH                               *
!                                  V                                   *
!                                  2 (Below Fish)                      *
!                                                                      *
!     Variables are relative to the fish, and the NFS compensates      *
!                 if fish x-axis orientation changes                   *
!                                                                      *
!***********************************************************************
 
!  This Subroutine Calls No Other Subroutines
 
 
     subroutine SPLINE              ! Bilinear Spline Algorithm obtained from
                                    !   "TWO DIMENSIONAL SPLINE INTERPOLATION
 
     use FISHY                               !    ALGORITHMS - Author: Helmuth Spath"
     use GLOBAL
     use GEOMC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: dv, dw, dx, dy, h3, h4, xh1, xh2
!
!*** End of declarations rewritten by SPAG
!
 
     do
 
 !Check for Boundary Violations: Horizontal Plane
 
         if(xlok>DLX(ilok))then                               ! Must look into the next segment downstream
             if((ilok + 1)>DS(fnbp))then                      ! Already at most downstream segment in branch
                 xlok = DLX(ilok)
             elseif(klok>KB(ilok + 1))then
                 xlok = DLX(ilok)
             else                                             ! More segments exist downstream of fish location
                 xlok = xlok - DLX(ilok)
                 ilok = ilok + 1                              ! Looking into next segment downstream
             endif
             if(xlok>DLX(ilok))cycle
         elseif(xlok<0)then                                   ! Must look into the next segment upstream
             if(ilok==CUS(fnbp))then                          ! Already at most upstream segment in branch
                 xlok = 0.0
             elseif(klok>KB(ilok - 1))then
                 xlok = 0.0
             else                                             ! More segments exist upstream of current segment
                 ilok = ilok - 1                              ! Looking into next segment upstream
                 xlok = DLX(ilok) + xlok
             endif
             if(xlok<0)cycle
         endif
         exit
     enddo
     do
 
!Check   for Boundary Violations: Vertical Plane
 
         if(zlok>H(klok, fjr))then                                ! Must look into the next layer below
             if(klok==KB(ilok))then                           ! Already at water body bottom
                 zlok = H(klok, fjr)
             else                                             ! More layers exist below fish location
                 zlok = zlok - H(klok, fjr)
                 klok = klok + 1                              ! Looking into next layer below
             endif
             if(zlok>H(klok, fjr))cycle
         elseif(zlok<0)then                                   ! Must look into the next layer above
             if(klok==ktwbf .OR. klok==1)then      ! SW 1/17/01        ! Already at water surface layer
                 zlok = 0.0
             else                                             ! More layers exist above fish location
                 klok = klok - 1                              ! Looking into next layer above
                 zlok = H(klok, fjr) + zlok
             endif
             if(zlok<0)cycle
         endif
         exit
     enddo
 
!Begin Bilinear Spline Interpolation Algorithm
 
     dv = xlok - 0.0                      ! XLOK = X-coord of interpolated location; 0.0 = X-coord of left nodes
     dw = zlok - 0.0                      ! ZLOK = Y-coord of interpolated location; 0.0 = Y-coord of upper nodes
     dx = DLX(ilok)                       ! DX = X distance between left nodes and right nodes
     dy = H(klok, fjr)                        ! DY = Y distance between upper nodes and lower nodes
 
     if(variable==1)then                  ! Horizontal Velocity Calculations
         xh1 = FLOWFIELD(klok, ilok, 1)    ! UpperLeft Node
         xh2 = FLOWFIELD(klok, ilok + 1, 1)
                                           ! UpperRight Node
         h3 = FLOWFIELD(klok + 1, ilok, 1)
                                          ! LowerLeft Node
         h4 = FLOWFIELD(klok + 1, ilok + 1, 1)
                                          ! LowerRight Node
     elseif(variable==2)then              ! Vertical Velocity Calculations
         xh1 = FLOWFIELD(klok, ilok, 2)    ! XH1,XH2,H3,H4 = Value of Vertical Velocity at 4 cell nodes
         xh2 = FLOWFIELD(klok, ilok + 1, 2)
         h3 = FLOWFIELD(klok + 1, ilok, 2)
         h4 = FLOWFIELD(klok + 1, ilok + 1, 2)
     elseif(variable==3)then              ! Temperature Calculations
         xh1 = WQFIELD(klok, ilok, 1)      ! XH1,XH2,H3,H4 = Value of Temperature at 4 cell nodes
         xh2 = WQFIELD(klok, ilok + 1, 1)
         h3 = WQFIELD(klok + 1, ilok, 1)
         h4 = WQFIELD(klok + 1, ilok + 1, 1)
     elseif(variable==4)then              ! Dissolved Oxygen Calculations
         xh1 = WQFIELD(klok, ilok, 2)      ! XH1,XH2,H3,H4 = Value of Dissolved Oxygen at 4 cell nodes
         xh2 = WQFIELD(klok, ilok + 1, 2)
         h3 = WQFIELD(klok + 1, ilok, 2)
         h4 = WQFIELD(klok + 1, ilok + 1, 2)
     endif
 
     value = xh1 + (xh2 - xh1)/dx*dv + (h3 - xh1)/dy*dw + (h4 - h3 - xh2 + xh1)&
           & /(dx*dy)*dv*dw
 
!End Bilinear Spline Interpolation Algorithm
 
     end subroutine SPLINE
