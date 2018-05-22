!*==gridplot.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************
!***********************************************************************
!*                                                                    **
!*   P R E P A R E   S Y S T E M   G R I D   F O R   D I S P L A Y    **
!*                                                                    **
!***********************************************************************
!***********************************************************************
! Number of TecPlot Zones       = # of timesteps            (Animation by Zone)
!           TecPlot Limit     ==> 32,700 per Data Set
! Number of TecPlot Data Points = # of nodes                (for Grid Display)
!                                 # of fish                 (for Fish Display)
!           TecPlot Limit     ==> 2 billion per Variable
! Number of TecPlot Variables   = # of flow & const. param  (for Grid Display)
!                                 # of fish parameters      (for Fish Display)
!           TecPlot Limit     ==> 32,700 per Data Set
! Number of TecPlot Data Sets   = 2*(# of branches)         (for Grid & Fish Displays)
!           TecPlot Limit     ==> 128
 
!  This Subroutine Calls the Following Subroutines:
!      INTERCONST,INTERFLOWF
 
 
     subroutine GRIDPLOT     ! This Subroutine preps Grid, Flow, & WQ info for TecPlot
 
 
     use FISHY
     use GLOBAL
     use GEOMC
     use SCREENC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: brchnn, ee, ele, ii, joiningtrib, k, kk, nd, nn, nnlast, totnn,&
              & water, wb
     character(1) :: dtyp1
     character(2) :: dtyp2
     logical, save :: fegrid
     character(12) :: filename
     character(6), save :: file_prefix
     character(4), save :: file_suffix
     integer, allocatable, dimension(:), save :: totbrnn, totele
     real :: vertflow, xdist, zdist
!
!*** End of declarations rewritten by SPAG
!
 
 
     data file_prefix/'Branch'/
     data file_suffix/'.dat'/
 
     if(nit==0)then
         fegrid = .FALSE.               ! The first pass is used to set up the FE grid
         allocate(totele(nbr), totbrnn(nbr))
                                            ! SW 1/14/01
     endif
     do
 
         if(fegrid)call INTERCONST                    ! Subroutine to interpolate Constituent values
         if(fegrid)call INTERFLOWF                    ! Subroutine to interpolate Flow Field values
         nn = 0                                       ! NN = Global Node # (UpperLeft Corner of cell(K,I))
         do wb = 1, nwb                               ! WB = Water Body # / NWB = Total # of Water Bodies
             do jb = BS(wb), BE(wb)                 ! JB = Branch #
                 if(.NOT.(fegrid .EQV. .FALSE.))then  ! Skip to setting up FE grid if NIT = 0
 
                     if(nit==0)then
                         if(jb<10)then                ! Creating the output file names based on Branch #
                             write(dtyp1, 9001)jb
9001                         format(i1)
                             filename = file_prefix // dtyp1 // file_suffix
                         elseif(jb<=64)then           ! TecPlot is limited to 128 Data Sets (i.e. Data Files)
                             write(dtyp2, 9002)jb   !   That leaves: 64 Data Sets for Grid Display and
9002                         format(i2)
                             filename = file_prefix // dtyp2 // file_suffix
                                                      !                64 Data Sets for Fish Display
                         else                           ! CHECK
                             write(datadebugfn, *)                             &
     &'ERROR: Exceeded Maximum Number of TecPlot Data Sets ==> # Data Sets = # &
     &of Branches'
                             stop
                         endif
 
                         open(20000 + jb, file = filename, status = 'UNKNOWN')
                                                                ! Opening the Output Files Based on Branch #
                         write(20000 + jb, 9003)jb            ! Creating Output File Header for TecPlot
9003                     format('TITLE = "GRID NODE INFO for Branch', i6, '"')
                         write(20000 + jb, 9004)                ! Creating Output File Header for TecPlot
9004                     format(                                               &
     &'VARIABLES = "X", "Z", "WATER", "HorizVel", "VertVel", "Temp", "DO", "KLa&
     &yer", "ISegment"')
                         write(20000 + jb, 9005)jday, totbrnn(jb), totele(jb)
                                                            ! Creating Output File Header for TecPlot
9005                     format('ZONE T="JDAY ', f9.2, '", N=', i8, ', E=', i8,&
                               &', F=FEPOINT, ET=QUADRILATERAL')
                     else                                         ! Subsequent Passes thru Subroutine
                         write(20000 + jb, *)' '
                         write(20000 + jb, 9006)jday, totbrnn(jb), totele(jb)
9006                     format('ZONE T="JDAY ', f9.2, '", N=', i8, ', E=', i8,&
                               &', F=FEPOINT, ET=QUADRILATERAL, D=(FECONNECT)')
                     endif
                 endif
                 zdist = H(2, wb)                                    ! Begin Calculations to set up FE grid
                 nnlast = nn
                 if(nnlast==0)nnlast = 1
                                   ! sw 1/9/01
                 brchnn = 0                         ! BRCHNN = Branch Node # (each Branch starts off with '0')
                 do k = 2, kmx                                    ! K = Layer #
         ! DO K=2,KMX
                     zdist = zdist - H(k, wb)                        ! ZDIST = Depth to Node in Branch JB
                     xdist = 0                                    ! XDIST = X Distance to Node in Branch JB
                     do i = US(jb), DS(jb) + 1 ! I = Segment #!  DO I=CUS(JB),DS(JB)+1             ! DO 118  I=US(JB),DS(JB)+1         ! I = Segment #
              !DO K=KTWB(WB),KB(I)   ! SW 5/2015
                         if((k - 1<=KB(i)) .OR. (k - 1<=KB(i - 1)))then
                                                                  ! All Nodes at or above Water Body Bottom
                             nn = nn + 1            ! NN = Global Node # (UpperLeft Corner of cell(K,I))
                             brchnn = brchnn + 1    ! BRCHNN = Branch Node # (each Branch starts off with '0')
                             if(fegrid .EQV. .FALSE.)then
                                                    ! This info is stored for later use
                                 NDINFO(nn, 1) = nn
                                 NDINFO(nn, 2) = brchnn
                                 NDINFO(nn, 3) = k
                                 NDINFO(nn, 4) = i
                                 NODES(k, i, 1) = xdist
                                 NODES(k, i, 2) = zdist
                                 totnn = nn
              !  ELSE
                             endif
                             if(fegrid)then
                                 if(k<KTWB(wb))then                     ! WATER used for value-blanking cells above 9952 surface in TecPlot
                                     if(.NOT.(showsky))then            ! Show Colors of Day Above Water Surface
                                         water = -1                    ! WATER = -1 ==> Cells above water surface are WHITE
                                     elseif((miltime>skynight) .OR.            &
                                       & (miltime<=skydawn))then       ! During the Night, the Sky is
                                         water = -100                  !   BLACK
                                     elseif((miltime>skydawn) .AND.            &
                                       & (miltime<=skyday))then        ! During the Morning, the Sky is
                                         water = -75                   !   YELLOW
                                     elseif((miltime>skyday) .AND.             &
                                       & (miltime<=skydusk))then       ! During the Day, the Sky is
                                         water = -50                   !   BLUE
                                     elseif((miltime>skydusk) .AND.            &
                                       & (miltime<=skynight))then      ! During the Evening, the Sky is
                                         water = -25                   !   ORANGE
                                     else                              ! CHECK
                                         write(datadebugfn, *)                 &
     &'ERROR: MILTIME not Corresponding with Intervals Set for Determining Sky &
     &Color for TecPlot'
                                         stop
                                     endif
                                     WQFIELD(k, i, 1) = WQFIELD(KTWB(wb), i, 1)
                                                           ! This is done to prevent odd-looking contours in TecPlot near
                                     WQFIELD(k, i, 2) = WQFIELD(KTWB(wb), i, 2)
                                                           !   the water surface since WQ does not exist above water surface
                                 else
                                     water = 1                         ! WATER = 1  ==> cell at or below current water surface
                                 endif
                                 if(.NOT.(wbskip .AND. (wb/=wbrun)))then
                                                                       ! Option:Skip output for Water Bodies, except WBRUN
                                     if(ABS(FLOWFIELD(k, i, 2))<vvelcap)then
                                                                       ! This keeps relatively large Vert Vel vectors in TecPlot
                                         vertflow = FLOWFIELD(k, i, 2)
                                     elseif(FLOWFIELD(k, i, 2)>=0)then !   from blacking out the screen since these vectors
                                         vertflow = vvelcap            !   would be so big
                                     else
                                         vertflow = -vvelcap
                                     endif
!                                    JOININGTRIB = 0
!                                    IF (RIMPBR(I,1).GT.0) JOININGTRIB = 1    
!                                    ! Branch joining right bank at Segment I
!                                    IF (LIMPBR(I,1).GT.0) JOININGTRIB = 1    
!                                    ! Branch joining left bank at Segment I
!                                    WRITE(20000+JB,9260) XDIST,ZDIST,WATER,& 
!                                    ! Outputting node information
                                     if(k==KTWB(wb))then
                                         zdist = ELWS(i)
                                     else
                                         zdist = ELWS(i) - DEPTHM(k, i)
                                     endif
 
                                                       !            ! Outputting node information  SW output elevation in m
                                     write(20000 + jb, 9007)xdist, zdist,      &
                                       & water, FLOWFIELD(k, i, 1), -vertflow, &
                                       & WQFIELD(k, i, 1), WQFIELD(k, i, 2), k,&
                                       & i         !                !   to a file for TecPlot to display
 
9007                                 format(f15.1, f10.2, i7, e15.2, e15.2,    &
                                       & f8.2, f8.2, i5, i5)
                                 endif
!                                .                   F10.2,I5)
                                                 ! XDIST                 = X Distance to Node in Branch JB
                                                 ! ZDIST                 = Depth to Node in Branch JB
                                                 ! NN                    = Global Node #
                                                 ! BRCHNN                = Branch Node #
                                                 ! K                     = Layer #
                                                 ! I                     = Segment #
                                                 ! JB                  = Branch #
                                                 ! WATER                 = Water in cell? [-1=No , 1=Yes]
                                                 ! FLOWFIELD(K,I,1)      = Velocity: Horizontal (m/s)
                                                 ! FLOWFIELD(K,I,2)                  Vertical   (m/s):W(K,I) is (+) Downward
                                                 ! FLOWFIELD(K,I,3)      = Acceleration: Horizontal (m/s^2)
                                                 ! FLOWFIELD(K,I,4)                      Vertical   (m/s^2)
                                                 ! WQFIELD(K,I,1)        = Temperature              (deg C)
                                                 ! WQFIELD(K,I,2)        = Dissolved Oxygen         (g/m^3)
                                                 ! B(K,I)                = Width of Cell (K,I)
                                                 ! JOININGTRIB           = Incoming Tributary at Segment I? [0=No , 1=Yes]
 
                             endif
 
                         endif
                         xdist = xdist + DLX(i)
    !ENDDO
    !ENDDO
 
                     enddo
                 enddo
                 if(.NOT.(fegrid))then                     ! If FEGRID = .TRUE., then FE connectivity already set up
                     ele = 0                               ! ELE = Element #
                     do kk = 2, kmx                        ! KK = Layer #
                         do ii = US(jb), DS(jb)        ! II = Segment #
                             if(kk<=KB(ii))then            ! Any Layer above Water Body Bottom
                                 ele = ele + 1
                                 do nd = nnlast, nn        ! NNLAST prevents having to scan the entire matrix for
                                     if((NDINFO(nd, 3)==kk) .AND.              &
                                      & (NDINFO(nd, 4)==ii))then
                                                 !          !   the desired information needed to determine connectivity
                                         CORNERS(jb, ele, 1) = NDINFO(nd, 2)
                                     elseif((NDINFO(nd, 3)==kk) .AND.          &
                                       & (NDINFO(nd, 4)==ii + 1))then
                                         CORNERS(jb, ele, 2) = NDINFO(nd, 2)
                                     elseif((NDINFO(nd, 3)==kk + 1) .AND.      &
                                       & (NDINFO(nd, 4)==ii + 1))then
                                         CORNERS(jb, ele, 3) = NDINFO(nd, 2)
                                     elseif((NDINFO(nd, 3)==kk + 1) .AND.      &
                                       & (NDINFO(nd, 4)==ii))then
                                         CORNERS(jb, ele, 4) = NDINFO(nd, 2)
                                     endif
                                 enddo
                             endif
                         enddo
                     enddo
                     totele(jb) = ele                    ! TOTELE(JB) = Total # of elements in Branch JB
                     totbrnn(jb) = brchnn                ! TOTBRNN(JB) = Total # of nodes in Branch JB
                 elseif(nit==0)then                        ! Once FE connectivity is established (below), it still
                     write(20000 + jb, *)' '             !   must be outputted the first time node info is outputted
                     if(.NOT.(wbskip .AND. (wb/=wbrun)))then
                                                           ! Option:Skip output for all Water Bodies, except WBRUN
                         do ee = 1, totele(jb)
                             write(20000 + jb, 9008)CORNERS(jb, ee, 1),        &
                                 & CORNERS(jb, ee, 2), CORNERS(jb, ee, 3),     &
                                 & CORNERS(jb, ee, 4) ! ! Writing connectivity information to output file
9008                         format(i8, i8, i8, i8)
                         enddo
                     endif
                 endif
             enddo
         enddo
 
         if(.NOT.(fegrid))then
             fegrid = .TRUE.                               ! Set FEGRID = .TRUE. once connectivity has been
             cycle                                         !   determined
         endif
         exit
     enddo
 
     end subroutine GRIDPLOT
