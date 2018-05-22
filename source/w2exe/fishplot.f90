!*==fishplot.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!***********************************************************************
!***********************************************************************
!*                                                                    **
!*         P R E P A R E   F I S H   F O R   D I S P L A Y            **
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
 
!  This Subroutine Calls No Other Subroutines
 
 
     subroutine FISHPLOT                 ! This Subroutine preps fish information for TecPlot
 
 
 
 
 
 
 
 
 
 
 
     use FISHY
     use GEOMC
     use GLOBAL
     use GDAYC
     use SCREENC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: brf, fbrch, fishnbp, haimp, k, netimp, water, wb
     real :: fag, fsz, netxput, netzbot, netztop, sndxdwn, sndxups, sndztop,   &
           & uf, vf, wbbottom, wf, xloc, yloc, zloc
     character(10) :: filename
     character(4), save :: file_prefix, file_suffix
     character(1) :: ftyp1
     character(2) :: ftyp2
!
!*** End of declarations rewritten by SPAG
!
 
 
 
     data file_prefix/'Part'/
     data file_suffix/'.dat'/
 
     if(nit==0)nzones = 0                             ! This # is inputted into a TecPlot Macro for animation
     nzones = nzones + 1                              ! This # is inputted into a TecPlot Macro for animation
     do wb = 1, nwb                                   ! WB = Water Body # / NWB = Total # of Water Bodies
         do jb = BS(wb), BE(wb)                     ! JB = Branch #
             NBRF(jb) = 0                           ! NBRF(JB) = # of Fish in Branch JB
             do fn = 1, nfish                         ! NFISH = Total # of Fish in System
                 if(INT(FISHES(fn, 12))/=1 .AND. DELAYDATE(fn)<=jday)then      ! If Fish left system do not output its info to TecPlot
                     fishnbp = INT(FISHES(fn, 6))     ! FISHNBP = Branch where current fish resides
                     if(fishnbp==jb)then
                         NBRF(jb) = NBRF(jb) + 1          ! NBRF    = # of Fish in Branch JB
                         BRCHFISH(jb, NBRF(jb), 1) = FISHES(fn, 1)
                                                          ! FIMP    = Segment IMP where fish is located
                         BRCHFISH(jb, NBRF(jb), 2) = FISHES(fn, 2)
                                                          ! FXLOC   = Fish location within segment IMP
                         BRCHFISH(jb, NBRF(jb), 3) = FISHES(fn, 3)
                                                          ! FKMP    = Layer KMP where fish is located
                         BRCHFISH(jb, NBRF(jb), 4) = FISHES(fn, 4)
                                                          ! FZLOC   = Fish location within layer KMP
                         BRCHFISH(jb, NBRF(jb), 5) = FISHES(fn, 5)
                                                          ! FYLOC   = Lateral fish location
                         BRCHFISH(jb, NBRF(jb), 6) = FISHES(fn, 6)
                                                          ! FNBP    = Branch where fish is located
                         BRCHFISH(jb, NBRF(jb), 7) = FISHES(fn, 7)
                                                          ! FSIZE   = Size (i.e., length) of fish in meters
                         BRCHFISH(jb, NBRF(jb), 8) = FISHES(fn, 8)
                                                          ! FAGE    = Age of the fish
                         BRCHFISH(jb, NBRF(jb), 9) = FISHES(fn, 9)
                                                          ! UFISH   = Longitudinal velocity of the fish
                         BRCHFISH(jb, NBRF(jb), 10) = FISHES(fn, 10)
                                                          ! VFISH   = Lateral velocity of the fish
                         BRCHFISH(jb, NBRF(jb), 11) = FISHES(fn, 11)
                                                          ! WFISH   = Vertical velocity of the fish
                     endif
                 endif
             enddo
 
             if(nit==0)then
 
                 if(jb<10)then                      ! Creating the output file names based on Branch #
                     write(ftyp1, 9001)jb
9001                 format(i1)
                     filename = file_prefix // ftyp1 // file_suffix
                 elseif(jb<=64)then                 ! TecPlot is limited to 128 Data Sets (i.e. Data Files)
                     write(ftyp2, 9002)jb           !   That leaves: 64 Data Sets for Grid Display and
9002                 format(i2)
                     filename = file_prefix // ftyp2 // file_suffix
                                                      !                64 Data Sets for Fish Display
                 else                                 ! CHECK
                     write(datadebugfn, *)                                     &
     &'ERROR: Exceeded Maximum Number of TecPlot Data Sets ==># Data Sets = # o&
     &f Branches'
                     stop
                 endif
 
                 open(30000 + jb, file = filename, status = 'UNKNOWN')
                                                                      ! Opening the Output Files Based on Branch #
                 write(30000 + jb, 9003)jb                          ! Creating Output File Header for TecPlot
9003             format('TITLE = "PARTICLE INFO for Branch', i6, '"')
                 write(30000 + jb, 9004)                              ! Creating Output File Header for TecPlot
9004             format(                                                       &
                 &'VARIABLES = "X", "Z", "FSIZE", "FAGE", "WATER", "ISEG", "K"'&
                & )
             else
                 write(30000 + jb, *)' '
             endif
             write(30000 + jb, 9005)jday, (1 + NBRF(jb))            ! Creating Output File Header for TecPlot
9005         format('ZONE T="JDAY ', f9.2, '", I=', i7, ', F=POINT')
             write(30000 + jb, 9009)0.01, 0.01, 0.0, 0.0, -1, 0, 0
                                                                ! This is a dumby fish so TecPlot can at least
                                                              !   open the file if no fish are in JB
             do brf = 1, NBRF(jb)                           ! NBRF = # of Fish in Branch JB
                 k = INT(BRCHFISH(jb, brf, 3))              ! K = Layer # where Fish is located
                 i = INT(BRCHFISH(jb, brf, 1))              ! I = Segment # where Fish is located
                 xloc = NODES(k, i, 1) + BRCHFISH(jb, brf, 2)
                                                            ! XLOC = X Distance (within Branch) where Fish is located
                 if(k<KTWB(wb))then
         !     ZLOC = NODES(KTWB(WB),I,2) - 0                   ! ZLOC(Adjusted): Due to inability of TecPlot to show
                     zloc = ELWS(i)
                           ! SW 5/2015 CHANGE TO ELEVATION RATHER THAN DEPTH
                 else                                         !                 actual Water Surface
         !     ZLOC = NODES(K,I,2) - BRCHFISH(JB,BRF,4)      ! ZLOC = Z Distance where Fish is located
                     zloc = EL(k, i) - BRCHFISH(jb, brf, 4)
                                               ! SW 5/2015
                 endif
                 yloc = BRCHFISH(jb, brf, 5)                ! YLOC = Y Distance where Fish is located
                 fbrch = INT(BRCHFISH(jb, brf, 6))          ! FBRCH = Branch # where Fish is located
                 fsz = BRCHFISH(jb, brf, 7)                 ! FSZ = Fish Size
                 fag = BRCHFISH(jb, brf, 8)                 ! FAG = Fish Age
                 water = 1                                    ! WATER = 1 (for value-blanking purposes in TecPlot)
                 uf = BRCHFISH(jb, brf, 9)                  ! UF = Fish Horizontal Velocity (m/s)
                 vf = BRCHFISH(jb, brf, 10)                 ! VF = Fish Lateral Velocity (m/s)
                 wf = BRCHFISH(jb, brf, 11)                 ! WF = Fish Vertical Velocity (m/s)
                 if(.NOT.(wbskip .AND. (wb/=wbrun)))write(30000 + jb, 9009)    &
                  & xloc, zloc, fsz, fag, water, i, k         ! Option:Skip output for all Water Bodies, except WBRUN
             enddo
             write(30000 + jb, 9006)month, gday, year, INT(hr), ampm, jday,    &
                                  & nzones                    ! Output info used as dynamic text in TecPlot Animation
9006         format(                                                           &
     &'TEXT X=47, Y=90, F=HELV-BOLD, HU=FRAME, AN=MIDCENTER, C=RED, H=2.5, T="'&
    & , a10, i3, ',', i5, 4x, i3, a2, '    (JDAY', f8.3, ')", ZN=', i6)
             write(30000 + jb, 9007)NBRF(jb), nzones      ! Output info used as dynamic text in TecPlot Animation
9007         format(                                                           &
     &'TEXT X=30.0, Y=21.0, F=HELV-BOLD, HU=FRAME, AN=MIDRIGHT, C=BLACK, H=2.1,&
     & T="', i6, '", ZN=', i6)
 
          !  The Following is Used for Creating and Positioning Static Text in the TecPlot Animation
 
             if(nit==0)then
                 write(30000 + jb, 9008)
9008             format(                                                       &
     &'TEXT X=25.0, Y=25.0, F=HELV-BOLD, HU=FRAME, AN=MIDCENTER, C=BLACK, H=2.1&
     &, T="# of Particles"')
             endif
         enddo
     enddo
9009  format(f15.2, f10.2, f9.2, f9.2, i4, 1x, i3, 1x, i3)
 
     end subroutine FISHPLOT
