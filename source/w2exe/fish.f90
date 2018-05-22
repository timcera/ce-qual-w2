!*==fish.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
!  This Subroutine Calls the Following Subroutines:
!      FIMPBR,GRIDPLOT,WHATJR,SPLINE,RANDOM
!      FINDNEWBR,FISHPLOT,INTERCONST,INTERFLOWF
 
 
     subroutine FISH
 
 
 
 
 
 
 
     use SURFHE
     use FISHY
     use GDAYC
     use SCREENC
     use GEOMC
     use GLOBAL
     use MAIN, ONLY:IWD, fish_particle_exist
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: dz
     integer :: grouplast, iw, jf, kk, n
!
!*** End of declarations rewritten by SPAG
!
 
 
     if(nit==0)then              ! The very first time the Subroutine FISH is called
 
!*****   Open 'Numerical Fish Surrogate' files
         call READ_FISH_DATA
         if(.NOT.fish_particle_exist)return
                                          !STOP ALL PROCESSING
         open(diagfn, file = 'DIAGNOSTICS.OUT', status = 'UNKNOWN')                 !FISH
         open(datadebugfn, file = 'DATADEBUG.OUT', status = 'UNKNOWN')                   !FISH
         open(finalfn, file = 'finalparticle.csv', status = 'unknown')                     !FISH
         open(initialfn, file = 'initialparticle.csv', status = 'unknown')
 
!        Allocate arrays
!        Read input data files
 
         allocate(fishes(nfish, fpara), wqfield(kmx, imx, wqp),                &
                & flowfield(kmx, imx, ffp), netcatch(numgillnets, nfish, 5),   &
                & sndcatch(numacoustics, 3*nfish, 2), rimpbr(imx, 2),          &
                & limpbr(imx, 2), nodes(kmx, imx, 2),                          &
                & ndinfo(imx*kmx + imx + kmx - 1, 4),                          &
                & corners(nbr, imx*kmx + imx + kmx - 1, 4),                    &
                & lastfield(kmx, imx, 2), brchfish(nbr, nfish, fpara),         &
                & nbrf(nbr))
         flowfield = 0.0
         wqfield = 0.0
         corners = 0
         ndinfo = 0
         nodes = 0.0
         netcatch = 0
         sndcatch = 0.0
         lastfield = 0.0
         brchfish = 0.0
         nbrf = 0
         fishes = 0.0
         wbskip = .FALSE.
         showsky = .FALSE.
         call FIMPBR                 ! Establishes data set of what segments are in each Branch
         seed = 92                   ! Seed for the Random Number Generator
 
 
!------------------------------------------ PARTICLE DATA INITIALIZATION ------------------------------------------
 
!Multiple Particles
         ngrpfh = 0                                  ! NGRPFH  = The number of particles/fish being initiated
         do jj = 1, nfishseg
             do fk = IFISHT(jj), IFISHB(jj)           !KTWB(2),(KB(FI)-2),1   ! Fish will be placed at each one of these Layer nodes
                 do fatpt = 1, nfishpcel             ! The Number of Fish that will be placed at each drop location
                     ngrpfh = ngrpfh + 1                 ! NGRPFH  = The number of particles/fish being initiated
                     fishes(ngrpfh, 1) = IFISH(jj)
                                              !FI        ! FIMP    = Initial segment IMP where particle/fish is released
                     fishes(ngrpfh, 3) = fk              ! FKMP    = Initial layer KMP where particle/fish is released
                 enddo
             enddo
         enddo
!_____________________________________________________
!Data    for All Fish
         grouplast = 1
         n = nfishpcel
         do fn = 1, nfish
             if(.NOT.line)then
                 fishes(fn, 4) = FZLOCI(fn)
                                       !0.0   ! FZLOC   = Location of fish within layer KMP from top side
             elseif(n==nfishpcel .OR. grouplast/=GROUP(fn))then
                 fishes(fn, 4) = 0.      ! FZLOC   = assumed to be zero for LINE
                 fimp = fishes(fn, 1)
                 call WHATJR
                 dz = H(fishes(fn, 3), fjr)/nfishpcel
                 grouplast = GROUP(fn)
                 n = 1
             else
                 fishes(fn, 4) = fishes(fn - 1, 4) + dz
                 n = n + 1
             endif
             fishes(fn, 2) = FXLOCI(fn)
                                       !0.0   ! FXLOC   = Location of fish within segment IMP from upstream side
             fishes(fn, 5) = B(INT(fishes(fn, 3)), INT(fishes(fn, 1)))*0.5
                                     ! FYLOC   = Lateral fish release location (from left bank in plan view)
          !                     Looking down on a segment
          !                                downstream
          !        y=0                        y=B/2                       y=B
          !         +---------------------------+--------------------------+
          !         |                                                      |
          !         |                                                      |
          !  L      |                           X                          |       R     ! location of particle in lateral
          !         |                                                      |
          !         |                                                      |
          !         |                                                      |
          !         +----------------------------+-------------------------+
          !                                 upstream
 
             call FINDBRANCH(fishes(fn, 1))
             fishes(fn, 6) = REAL(fnbp)
                                       !6     ! FNBP    = Branch where fish is released
             fishes(fn, 7) = fsize
                                  !0.178 ! FSIZE   = Size (i.e., length) of fish in meters (1 inch = 0.0254 meters)
             fishes(fn, 8) = fage
                                 !2.0   ! FAGE    = Age of the fish
             fishes(fn, 9) = ufish
                                  !0.0   ! UFISH   = Initial longitudinal velocity of the fish relative to water
             fishes(fn, 10) = vfish
                                  !0.0   ! VFISH   = Initial lateral velocity of the fish relative to water
             fishes(fn, 11) = wfish
                                  !0.0  ! WFISH   = Initial vertical velocity of the fish relative to water
             fishes(fn, 12) = 0      ! ISWITCH = Is fish still in system: Yes=0  No=1
             fishes(fn, 13) = 0      ! FSNAG   = Gillnet # fish is snagged in; = 0 if fish not in a gillnet
         enddo
     endif
 
!*********************************************************************************************
!******************************  NO MORE INPUT BELOW THIS POINT  *****************************
!*********************************************************************************************
 
 
!Initialize Parameters Used to Evaluate Implementation of the Stimuli-Response Rules,
!    Prepare Text for Output Purposes, and Output Input Data to
 
 
!    DIAGNOSTICS.OUT File
     if(nit==0)then
         vvelcap = 4.E-4            ! Any Vert Vel that exceeds this value will have its TECPLOT VECTOR truncated back to this value
                                   !   Suggested Value = 4E-4
!Set     Counting and Averaging Variables to Zero
 
         fx00count = 0                           ! Counts # of times No Gradient is Dominant in the X-direction
         fxhvcount = 0                           ! Counts # of times Horizontal Velocity is Dominant in the X-direction
         fxvvcount = 0                           ! Counts # of times Vertical Velocity is Dominant in the X-direction
         fxdocount = 0                           ! Counts # of times Dissolved Oxygen is Dominant in the X-direction
         fxtpcount = 0                           ! Counts # of times Temperature is Dominant in the X-direction
         fxrdcount = 0                           ! Counts # of times Randomization is Dominant in the X-direction
         totfxwgt = 0                            ! Sums All Weights Used in making X-Movement Decisions
         fz00count = 0                           ! Counts # of times No Gradient is Dominant in the Z-direction
         fzhvcount = 0                           ! Counts # of times Horizontal Velocity is Dominant in the Z-direction
         fzvvcount = 0                           ! Counts # of times Vertical Velocity is Dominant in the Z-direction
         fzdocount = 0                           ! Counts # of times Dissolved Oxygen is Dominant in the Z-direction
         fztpcount = 0                           ! Counts # of times Temperature is Dominant in the Z-direction
         fzrdcount = 0                           ! Counts # of times Randomization is Dominant in the Z-direction
         totfzwgt = 0                            ! Sums All Weights Used in making Z-Movement Decisions
         xhvtrunc = 0                            ! Counts # of times Scaled Horz Vel Gradient in X-dir must be truncated to 1.0
         xhvtrucsum = 0                          ! Sums the Values of all Scaled Horz Vel Gradients that were Truncated in X-dir
         zhvtrunc = 0                            ! Counts # of times Scaled Horz Vel Gradient in Z-dir must be truncated to 1.0
         zhvtrucsum = 0                          ! Sums the Values of all Scaled Horz Vel Gradients that were Truncated in Z-dir
         xvvtrunc = 0                            ! Counts # of times Scaled Vert Vel Gradient in X-dir must be truncated to 1.0
         xvvtrucsum = 0                          ! Sums the Values of all Scaled Vert Vel Gradients that were Truncated in X-dir
         zvvtrunc = 0                            ! Counts # of times Scaled Vert Vel Gradient in Z-dir must be truncated to 1.0
         zvvtrucsum = 0                          ! Sums the Values of all Scaled Vert Vel Gradients that were Truncated in Z-dir
         xdotrunc = 0                            ! Counts # of times Scaled Diss Oxyg Gradient in X-dir must be truncated to 1.0
         xdotrucsum = 0                          ! Sums the Values of all Scaled Diss Oxyg Gradients that were Truncated in X-dir
         zdotrunc = 0                            ! Counts # of times Scaled Diss Oxyg Gradient in Z-dir must be truncated to 1.0
         zdotrucsum = 0                          ! Sums the Values of all Scaled Diss Oxyg Gradients that were Truncated in Z-dir
         xtptrunc = 0                            ! Counts # of times Scaled Temperature Gradient in X-dir must be truncated to 1.0
         xtptrucsum = 0                          ! Sums the Values of all Scaled Temperature Gradients that were Truncated in X-dir
         ztptrunc = 0                            ! Counts # of times Scaled Temperature Gradient in Z-dir must be truncated to 1.0
         ztptrucsum = 0                          ! Sums the Values of all Scaled Temperature Gradients that were Truncated in Z-dir
         countfortruc = 0                        ! Counts total # of times Stimuli-Response Rules are implemented
 
!Preparing Text for TecPlot Animation and Output Files
 
 
         if(particle)then
             txtparticle = '  ON'
         else
             txtparticle = ' OFF'
         endif
 
 
!        initial fish output
         write(initialfn, '(a291)')                                            &
     &'Part#,Seg#,XLocationwithinSegmentfromUpstreamSide(m),Layer#,VerticalDist&
     &fromTop(m),LateralDistfromLeftBank,Branch#,ParticleInModel(=0),JDAYleftsy&
     &stem,DetentionTime(days),RemovalMechanism,SedVelocity(m/d),DateStart'
         do jf = 1, nfish
             write(initialfn, '(i7,",",1x,10(f10.3,","),f8.4,",",f12.4)')jf,   &
                 & fishes(jf, 1), fishes(jf, 2), fishes(jf, 3), fishes(jf, 4), &
                 & fishes(jf, 5), fishes(jf, 6), fishes(jf, 12), fishes(jf, 14)&
                 & , fishes(jf, 15), fishes(jf, 16), SEDVEL(jf), DELAYDATE(jf)
         enddo
 
         close(initialfn)
 
     endif
 
 
!Specification of the Time Stepping used in Running the NFS (W2 may, at times, run at timesteps of down to 2 minutes!)
 
!    COMPUTE DYNAMIC NFSFREQ SW 1/15/01
 
     nfsfreq = dlt/86400.
 
     if(nit==0)lrunday = jday                                 ! LRUNDAY is used in determining NFS run frequency
     rundiff = jday - lrunday                                 ! RUNDIFF is used in determining NFS run frequency
 
 
!Specification of when and how often to output information/data
 
     if(nit==0)ljday = jday                                   ! LJDAY is used for TecPlot output frequency purposes
     jdaydiff = jday - ljday                                  ! JDAYDIFF is used for TecPlot output frequency purposes
     outfreq = outfreqp
 
     if(jdaydiff>=outfreq)then          ! Set FCOUNT = 1, so information will be outputted
         fcount = 1                                           ! Information outputted iff FCOUNT = 1
         ljday = jday
     else
         fcount = 0                                           ! Set FCOUNT = 0, so information won't be outputted
     endif
 
 
!Calculate Time of Day/Night
 
     miltime = (jday - REAL(INT(jday)))*24.       !(JDAY-REAL(JDAYG))*24                         ! Military Time - MILTIME used for NFS Calculations
     hr = miltime                                             ! Convert MILTIME to hour of day - HR used for TecPlot
     ampm = 'am'                                              ! Is HR am or pm?
     if(hr>=12)ampm = 'pm'                                    ! Is HR am or pm?
     if(hr>=13)hr = hr - 12                                   ! Value of HR after 12 noon
     if(hr<1)hr = hr + 12                                     ! Value of HR if time is between 12 midnight and 1am
 
 
!Call the Subroutine GRIDPLOT to set up the FE Grid for TecPlot and to Calc node information for Traffic Rules
 
     if((fcount==1) .OR. (nit==0))call GRIDPLOT               ! This prepares the FE Grid for output display and
                                                              !   calculates the Flow and WQ information at nodes
 
!Load Fish Information and Begin NFS Logic
 
      !IF (JDAY .LT. DELAYDATE) GOTO 28
 
     call INTERCONST                                          ! Subroutine to interpolate Constituent values
     call INTERFLOWF                                          ! Subroutine to interpolate Flow Field values
 
     do fn = 1, nfish                ! This is the only loop where statements are not indented
         if(jday<DELAYDATE(fn))cycle
         fimp = INT(fishes(fn, 1))   ! Segment IMP where fish is located
         fxloc = fishes(fn, 2)       ! Location of fish within segment IMP from upstream side
         fkmp = INT(fishes(fn, 3))   ! Layer KMP where fish is located
         fzloc = fishes(fn, 4)       ! Location of fish within layer KMP from top side
         fyloc = fishes(fn, 5)       ! Location of fish laterally from left bank (in plan view)
         fnbp = INT(fishes(fn, 6))   ! Branch where fish is located
         fsize = fishes(fn, 7)       ! Size (i.e., length) of fish in meters
         fage = fishes(fn, 8)        ! Age of the fish
         ufish = fishes(fn, 9)       ! Longitudinal velocity of fish relative to water
         vfish = fishes(fn, 10)      ! Lateral velocity of fish relative to water
         wfish = fishes(fn, 11)      ! Vertical velocity of fish relative to water
         iswitch = INT(fishes(fn, 12))
                                     ! Is fish still in system?: Yes=0  No=1
         fsnag = INT(fishes(fn, 13)) ! Gillnet # fish is snagged in; = 0 if fish not in a gillnet
 
         if(iswitch==0)then          ! If fish has already left the system, then skip to the next fish
 
!Determine   Horizontal & Vertical Velocity, Temperature, and DO Gradients in
!            Fish Sensory Sphere: How far to seek out gradients (i.e. edges of
 
!            the fish sensory sphere given timestep used)?
             call RANDOM(seed, rrr)                                     ! Subroutine which calculates a random number
 
             fxsensph = (fbdysearch*fsize)*rrr                          ! FXSENSPH = Leading edge of the sensory sphere
             bxsensph = (bbdysearch*fsize)*rrr                          ! BXSENSPH = Trailing edge of the sensory sphere
             uzsensph = (ubdysearch*fsize)*rrr                          ! UZSENSPH = Upper edge of the sensory sphere
             dzsensph = (dbdysearch*fsize)*rrr                          ! DZSENSPH = Lower edge of the sensory sphere
 
             call WHATJR                                                ! What waterbody,JR, is the fish in?
             ktwbf = KTWB(fjr)                                            ! Water Surface Layer of waterbody fish is located in
 
!            Find Flow and Water Quality Constituent Values at places of
 
!            interest:
             dir = 1                                          ! DIR = Direction Fish is swimming: 1=downstream, 2=upstream
             if(ufish<0)dir = -1                              ! UFISH is (+) when Fish is swimming downstream
             do variable = 1, 4                               ! There are 4 Variables: 1=Horz Vel, 2=Vert Vel, 3=Temp, 4=Diss Oxy
                 do locate = 1, 5                             ! There are 5 Locations around Fish of interest-See Subroutine SPLINE
                     ilok = fimp                              ! Segment I where Variable Value is desired
                     xlok = fxloc                             ! Location within Segment I where Variable Value is desired
                     klok = fkmp                              ! Layer K where Variable Value is desired
                     zlok = fzloc                             ! Location within Layer K where Variable Value is desired
                     if((locate==1) .AND. (dir==1))xlok = fxloc + fxsensph
                                                                        ! Location forward of Fish when Fish swimming DOWNstream
                     if(locate==2)zlok = fzloc + dzsensph               ! Location below Fish
                     if((locate==3) .AND. (dir==1))xlok = fxloc - bxsensph
                                                                        ! Location behind Fish when Fish swimming DOWNstream
                     if(locate==4)zlok = fzloc - uzsensph               ! Location above Fish
                     if((locate==1) .AND. (dir== - 1))xlok = fxloc - fxsensph
                                                                        ! Location forward of Fish when Fish swimming UPstream
                     if((locate==3) .AND. (dir== - 1))xlok = fxloc + bxsensph
                                                                        ! Location behind Fish when Fish swimming UPstream
                     call SPLINE                                        ! Bilinear Spline Interpolation Subroutine
                     if(variable==1)FXVEL(locate) = value               ! Horiz Vel at Location of Interest: (+) is Downstream !!!
                     if(variable==2)FZVEL(locate) = value               ! Vert Vel at Location of Interest:  (+) is Downward !!!
                     if(variable==3)FISHTEMP(locate) = value            ! Temperature at Location of Interest
                     if(variable==4)FISHDO(locate) = value              ! Dissolved Oxygen at Location of Interest
                 enddo
             enddo
 
 
!            Determine lateral velocity based on lateral withdrawals or
!            connected branches  SW 2/01/01
             call LATERAL_VELOCITY
!            FYVEL = 0.0                                                     !
 
 
!            Lateral Velocity Calculate (Linear) Gradients within Fish Sensory
!            Sphere:
             call PART_TRANSPORT
 
 
 
!***************************************************************
!            C H E C K   F O R   B O U N D A R Y   V I O L A T I O N S   *
!***************************************************************
 
 
!Check       for Boundary Violations: Lateral Direction Check (1 of 2)
 
             if(fyvel==0.0)then
                         ! SW 2/01/01 logic for removing particles laterally also
                 if(fyloc<0)then
                                ! reflect particles
                     if(hitstickside)then
                         iswitch = 1
                         fishes(fn, 14) = jday   ! Time fish left system
                         fishes(fn, 15) = jday - DELAYDATE(fn)
                                                     ! Detention time of fish in system
                         fishes(fn, 12) = 1.0    ! FISH LEAVES SYSTEM
                         fishes(fn, 16) = 6.0    ! Particle leaves by hitting side wall and sticking-LHS
                     else
                         fyloc = B(fkmp, fimp)*yrefl
                     endif
                 elseif(fyloc>B(fkmp, fimp))then
                     if(hitstickside)then
                         iswitch = 1
                         fishes(fn, 14) = jday   ! Time fish left system
                         fishes(fn, 15) = jday - DELAYDATE(fn)
                                                     ! Detention time of fish in system
                         fishes(fn, 12) = 1.0    ! FISH LEAVES SYSTEM
                         fishes(fn, 16) = 7.0    ! Particle leaves by hitting side wall and sticking-RHS
                     else
                         fyloc = B(fkmp, fimp)*(1 - yrefl)
                     endif
                 endif
             elseif(fyloc<0 .AND. fyvel<0.0 .AND. limpbr(fimp, 1)==0)then
                                                                     !
            ! check for a withdrawal
                 do iw = 1, nwd
                     if(fimp==IWD(iw))then
                                       ! remove particle through withdrawal
                         iswitch = 1
                         fishes(fn, 14) = jday   ! Time fish left system
                         fishes(fn, 15) = jday - DELAYDATE(fn)
                                                     ! Detention time of fish in system
                         fishes(fn, 12) = 1.0    ! FISH LEAVES SYSTEM
                         fishes(fn, 16) = 1.0    ! Signifies the lateral removal by withdrawal
                         write(diagfn, *)                                      &
                              &'Withdrawal: Particle Leaves System on JDAY:',  &
                             & jday, 'Data:', fishes(fn, :)
                         goto 50
                     endif
                 enddo
                 if(hitstickside)then
                     iswitch = 1
                     fishes(fn, 14) = jday       ! Time fish left system
                     fishes(fn, 15) = jday - DELAYDATE(fn)
                                                     ! Detention time of fish in system
                     fishes(fn, 12) = 1.0        ! FISH LEAVES SYSTEM
                     fishes(fn, 16) = 6.0        ! Particle leaves by hitting side wall and sticking-LHS
                 else
                     fyloc = B(fkmp, fimp)*yrefl
                                           ! reflect
                 endif
             elseif(fyloc>B(fkmp, fimp) .AND. fyvel>0.0 .AND. rimpbr(fimp, 1)  &
                  & >0)then
                 write(diagfn, *)                                              &
                          &'LateralRIGHT: Particle move to new branch on JDAY:'&
                         & , jday, 'FYVEL=', fyvel, 'OLD I:', fimp, 'NEW I:',  &
                         & rimpbr(fimp, 1), ' OLD branch:', fnbp,              &
                          &' NEW branch:', rimpbr(fimp, 2)
                 fimp = rimpbr(fimp, 1)
                 fyloc = B(fkmp, fimp)*0.5
                                      ! PLACE IN MIDDLE LATERALLY
                 fxloc = DLX(fimp)*0.9
                                      ! PLACE IN NEAR END OF SEGMENT
 
           !ISWITCH = 1
           !! SW 2/01/01 Track timing of fish movement from system
           !fishes(fn,14)=JDAY               ! Time fish left system
           !fishes(fn,15)=JDAY-DELAYDATE     ! Detention time of fish in system
           !fishes(fn,16)=1.0                ! Signifies the lateral removal
           !fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
           !WRITE(DIAGFN,*) 'LateralRIGHT: Particle Leaves System on JDAY:',jday,'FYVEL=',fyvel, 'Data:',fishes(fn,:)
           ! End Section SW 2/01/01
           !GOTO 20
           ! Particle is removed       ***need to transfer to another branch - it is only removed****
             elseif(fyloc<0.0 .AND. fyvel<0.0 .AND. limpbr(fimp, 1)>0)then
                 write(diagfn, *)                                              &
                           &'LateralLEFT: Particle move to new branch on JDAY:'&
                          & , jday, 'FYVEL=', fyvel, 'OLD I:', fimp, 'NEW I:', &
                          & limpbr(fimp, 1), ' OLD branch:', fnbp,             &
                           &' NEW branch:', limpbr(fimp, 2)
                 fimp = limpbr(fimp, 1)
                 fyloc = B(fkmp, fimp)*0.5
                                      ! PLACE IN MIDDLE LATERALLY
                 fxloc = DLX(fimp)*0.9
                                      ! PLACE NEAR END OF SEGMENT
           !ISWITCH = 1
           !! SW 2/01/01 Track timing of fish movement from system
           !fishes(fn,14)=JDAY               ! Time fish left system
           !fishes(fn,15)=JDAY-DELAYDATE     ! Detention time of fish in system
           !fishes(fn,16)=1.0                ! Signifies the lateral removal
           !fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
           !WRITE(DIAGFN,*) 'LateralLEFT: Particle Leaves System on JDAY:',jday,'FYVEL=',fyvel, 'Data:',fishes(fn,:)
           ! End Section SW 2/01/01
           !GOTO 20
 
             elseif(fyloc>B(fkmp, fimp) .AND. rimpbr(fimp, 1)==0.0)then
                                                                      ! reflect
                 if(hitstickside)then
                     iswitch = 1
                     fishes(fn, 14) = jday       ! Time fish left system
                     fishes(fn, 15) = jday - DELAYDATE(fn)
                                                     ! Detention time of fish in system
                     fishes(fn, 12) = 1.0        ! FISH LEAVES SYSTEM
                     fishes(fn, 16) = 7.0        ! Particle leaves by hitting side wall and sticking-RHS
                 else
 
                     fyloc = B(fkmp, fimp)*(1 - yrefl)
                 endif
 
             elseif(fyloc<0.0 .AND. limpbr(fimp, 1)==0)then
                                                           ! reflect
                 if(hitstickside)then
                     iswitch = 1
                     fishes(fn, 14) = jday       ! Time fish left system
                     fishes(fn, 15) = jday - DELAYDATE(fn)
                                                     ! Detention time of fish in system
                     fishes(fn, 12) = 1.0        ! FISH LEAVES SYSTEM
                     fishes(fn, 16) = 6.0        ! Particle leaves by hitting side wall and sticking-LHS
                 else
                     fyloc = B(fkmp, fimp)*yrefl
                 endif
             endif
 
!Check       for Boundary Violations: Horizontal Direction
 
             if(fxloc>DLX(fimp))then
 
                 if(fimp<DS(fnbp))then
                     fimp = fimp + 1
                     fxloc = fxloc - DLX(fimp)
                 elseif(FXVEL(5)>0.0 .AND. DHS(fnbp)==0)then
               ! LOSE PARTCLE THROUGH DAM
                     iswitch = 1
                     fishes(fn, 14) = jday       ! Time fish left system
                     fishes(fn, 15) = jday - DELAYDATE(fn)
                                                     ! Detention time of fish in system
                     fishes(fn, 12) = 1.0        ! FISH LEAVES SYSTEM
                     fishes(fn, 16) = 5.0
                                   ! Particle leaves at downstream structure/dam/hydraulic structure
                     write(diagfn, *)'At Dam: Particle Leaves System on JDAY:',&
                                   & jday, 'Data:', fishes(fn, :)
                     goto 50
                 elseif(FXVEL(5)>0.0 .AND. DHS(fnbp)== - 1)then
               ! LOSE PARTCLE to external head BC
                     iswitch = 1
                     fishes(fn, 14) = jday       ! Time fish left system
                     fishes(fn, 15) = jday - DELAYDATE(fn)
                                                     ! Detention time of fish in system
                     fishes(fn, 12) = 1.0        ! FISH LEAVES SYSTEM
                     fishes(fn, 16) = 4.0
                                   ! Particle leaves at external head BC Downstream
                     write(diagfn, *)                                          &
       &'At External Downstream Head Boundary: Particle Leaves System on JDAY:'&
      & , jday, 'Data:', fishes(fn, :)
                     goto 50
                 elseif(FXVEL(5)>0.0 .AND. DHS(fnbp)>0)then
                     write(diagfn, *)                                          &
                                   &'DHS: Particle move to new branch on JDAY:'&
                                  & , jday, 'FXVEL=', FXVEL(5), 'OLD I:', fimp,&
                                   &'NEW I:', DHS(fnbp), ' OLD branch:', fnbp
                     fimp = DHS(fnbp)
                     fyloc = B(fkmp, fimp)*0.5
                                          ! PLACE IN MIDDLE LATERALLY
                     fxloc = DLX(fimp)*0.5
                                          ! PLACE IN MIDDLE OF SEGMENT
               ! MOVE PARTICLE TO NEW BRANCH
                !ISWITCH = 1
                !fishes(fn,14)=JDAY               ! Time fish left system
                !fishes(fn,15)=JDAY-DELAYDATE     ! Detention time of fish in system
                !fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
                !WRITE(DIAGFN,*) 'At Dam: Particle Leaves System on JDAY:',jday, 'Data:',fishes(fn,:)
                 else
               ! REFLECT OFF DAM
                     fxloc = DLX(fimp)*(1.0 - xrefl)
                 endif
 
             elseif(fxloc<0)then
 
                 if(fimp==CUS(fnbp) .AND. UHS(fnbp)==0)then   ! The UpStream end of Branch JBP or UNBP
                     fxloc = DLX(fimp)*xrefl
                                         ! REflect partilce off UPSTREAM boundary
 
                 elseif(fimp==CUS(fnbp) .AND. UHS(fnbp)== - 1 .AND. FXVEL(5)   &
                      & <0.0)then
                            ! LOSE PARTCLE to external head BC
                     iswitch = 1
                     fishes(fn, 14) = jday       ! Time fish left system
                     fishes(fn, 15) = jday - DELAYDATE(fn)
                                                     ! Detention time of fish in system
                     fishes(fn, 12) = 1.0        ! FISH LEAVES SYSTEM
                     fishes(fn, 16) = 3.0
                                   ! Particle leaves at external head BC Upstream
                     write(diagfn, *)                                          &
         &'At External Upstream Head Boundary: Particle Leaves System on JDAY:'&
        & , jday, 'Data:', fishes(fn, :)
                     goto 50
                 elseif(fimp==CUS(fnbp) .AND. UHS(fnbp)== - 1 .AND. FXVEL(5)   &
                      & >=0.0)then
                     fxloc = DLX(fimp)*xrefl
                                        ! REflect partilce off UPSTREAM boundary
                 elseif(fimp==CUS(fnbp) .AND. UHS(fnbp)>0 .AND. FXVEL(5)>=0.0) &
                      & then
                     fxloc = DLX(fimp)*xrefl   ! REflect partilce off UPSTREAM boundary
                 elseif(fimp==CUS(fnbp) .AND. UHS(fnbp)>0 .AND. FXVEL(5)<0.0)  &
                      & then
                     write(diagfn, *)                                          &
                                   &'UHS: Particle move to new branch on JDAY:'&
                                  & , jday, 'FXVEL=', FXVEL(5), 'OLD I:', fimp,&
                                   &'NEW I:', DHS(fnbp), ' OLD branch:', fnbp
                     fimp = UHS(fnbp)
                     fyloc = B(fkmp, fimp)*0.5
                                          ! PLACE IN MIDDLE LATERALLY
                     fxloc = DLX(fimp)*0.5
                                          ! PLACE IN MIDDLE OF SEGMENT
                 elseif(fimp==CUS(fnbp) .AND. UHS(fnbp)< - 1)then      ! at a DAM
                     fxloc = DLX(fimp)*xrefl   ! REflect partilce off UPSTREAM boundary
                 elseif(fimp>CUS(fnbp))then
                     fimp = fimp - 1                           ! There are still segments upstream for the fish to
                     fxloc = DLX(fimp) + fxloc                 !    move to
                 else
                     write(diagfn, *)                                          &
                &'***NO RESOLUTION OF FXLOC<0*** JDAY,I,K,BR#,FN#,FYLOC,FXLOC:'&
               & , jday, fimp, fkmp, fnbp, fn, fyloc, fxloc
                 endif
 
        !ELSEIF(FXLOC.EQ.0.0.and.cus(fnbp).eq.fimp.and.fxvel(5).eq.0.)then
        !     FXLOC = DLX(FIMP)*XREFL    ! REflect partilce off UPSTREAM boundary
        !
        !ELSEIF(FXLOC == DLX(FIMP) .AND. DS(FNBP)==FIMP   .AND.  FXVEL(5) <= 0.0)THEN
        !     FXLOC = DLX(FIMP)*(1.0-XREFL)  ! REFLECT OFF DOWNSTREAM BOUNDARY
        !ELSE
 
             endif
 
 
!Check       for Boundary Violations: Vertical Direction Check (1 of 2)
 
             ktwbf = KTWB(fjr)                       ! Water Surface Layer of waterbody fish is located in
!            This first section below checks and corrects fish vertically
!            above or below the bottom
             do while ((fkmp<=KTI(fimp)) .OR. (fkmp>KB(fimp)))
                                   !   ! SW 2/01/01 Change so that do a check if in "air" or not
                 if(fkmp>KB(fimp))then                 !  FKMP is below KB
                     fkmp = fkmp - 1
                     fzloc = H(fkmp, fjr)*(1 - zbotrefl)
                                                       !    Logic for rebounding off of the bottom
                     cycle
                 elseif(fkmp<=KTI(fimp))then          ! Don't leave fish in the air - move to surface layer
                     fkmp = KTI(fimp)    ! check and make sure the fish/particles are below the water surface   ! SW 2/01/01
                     if(KTI(fimp)==ktwbf)then
                         if(fzloc<Z(fimp))fzloc = Z(fimp) + H(fkmp, fjr)       &
                          & *(1 - zsurrefl)
                     elseif(fzloc<=H(KTI(fimp), fjr) + Z(fimp))then
                         if(Z(fimp)>0.0)write(datadebugfn, *)                  &
                           &'Debug Error in vertical reflection:FN,I,K,JDAY',  &
                          & fn, fimp, fkmp, jday
                         fzloc = (H(KTI(fimp), fjr) + Z(fimp)) + H(fkmp, fjr)  &
                               & *(1 - zsurrefl)
                     endif
                 endif
                 exit
             enddo
             do
 
!Check           for Boundary Violations: Vertical Direction Check (2 of 2)
 
                 if(fzloc>=H(fkmp, fjr))then               ! Did fish move down below current layer?  -- YES
                     if(debug)write(datadebugfn, *)'FZLOC GT H'
                     if(fkmp + 1<=KB(fimp))then        ! Did fish move below bottom layer? -- YES
                         fzloc = fzloc - H(fkmp, fjr)
                         fkmp = fkmp + 1
                         if(fzloc>H(fkmp, fjr))cycle       ! Did fish move down more than one layer
                     elseif(hitstickbottom)then
                  ! PARTICLE LOST FROM SYSTEM
                                              ! LOSE PARTCLE
                         iswitch = 1
                         fishes(fn, 14) = jday   ! Time fish left system
                         fishes(fn, 15) = jday - DELAYDATE(fn)
                                                     ! Detention time of fish in system
                         fishes(fn, 12) = 1.0    ! FISH LEAVES SYSTEM
                         fishes(fn, 16) = 2.0
                                   ! Particle hits bottom and sticks
                         write(diagfn, *)                                      &
              &'Particle Hit and Stick Bottom: Particle Leaves System on JDAY:'&
             & , jday, 'Data:', fishes(fn, :)
                     else
                         fzloc = H(fkmp, fjr)*(1 - zbotrefl)       !    Logic for rebounding off of the bottom
 
                     endif
                 elseif(fzloc<0)then                   ! Did fish move up?  --  YES
                     if(debug)write(datadebugfn, *)'FZLOC LT 0'
                     if(fkmp==(ktwbf + 1))then          ! Is fish just below surface layer KTWBF?  --  YES
                         if(debug)write(datadebugfn, *)'FKMP EQ KTWBF+1'
                         if(ABS(fzloc)>(H(ktwbf, fjr) - Z(fimp)))then
                                                            ! Did fish move above water surface?  --  YES
                             fromktbot = (H(ktwbf, fjr) - Z(fimp))             &
                               & *(1 - zsurrefl)
                         else                          ! Did fish move above water surface?  --  NO
                             fromktbot = ABS(fzloc)
                         endif
                         surfcalc = 1
                     elseif(fkmp==ktwbf)then
                         if(debug)write(datadebugfn, *)'FKMP EQ KTWBF'
                         if((ABS(fzloc) + H(ktwbf, fjr))                       &
                          & >(H(ktwbf, fjr) - Z(fimp)))then
                             fromktbot = (H(ktwbf, fjr) - Z(fimp))             &
                               & *(1 - zsurrefl)
                         else
                             fromktbot = ABS(fzloc) + H(ktwbf, fjr)
                         endif
                         surfcalc = 1
                     elseif(fkmp<=(ktwbf - 1))then
                         if(debug)write(datadebugfn, *)'FKMP LE KTWBF-1'
                         fromktbot = 0.0
                         do kk = ktwbf, (fkmp + 1), -1
                             fromktbot = fromktbot + H(kk, fjr)
                         enddo
                         fromktbot = fromktbot + ABS(fzloc)
                         if(fromktbot>(H(ktwbf, fjr) - Z(fimp)))then
                             fromktbot = (H(ktwbf, fjr) - Z(fimp))             &
                               & *(1 - zsurrefl)
                         else
                             fromktbot = fromktbot
                         endif
                         surfcalc = 1
                     else
                         if(debug)write(datadebugfn, *)'ORDINARY FZLOC LT 0'
                         fkmp = fkmp - 1
                         fzloc = H(fkmp, fjr) + fzloc
                         if(fzloc<0)cycle              ! Did fish move up more than one layer - must improve this.
                         surfcalc = 0
                     endif
                     if(surfcalc==1)then
                         fkmptemp = ktwbf
                         do while (fromktbot>H(fkmptemp, fjr))
                             fromktbot = fromktbot - H(fkmptemp, fjr)
                             fkmptemp = fkmptemp - 1
                         enddo
                         fzloc = H(fkmptemp, fjr) - fromktbot
                         fkmp = fkmptemp
                     endif
                 endif
                 exit
             enddo
             if(fkmp<KTI(fimp))then
                 write(datadebugfn, 9001)fkmp, fzloc, KTI(fimp),               &
                     & (H(ktwbf, fjr) - Z(fimp))
9001             format('FKMP.LT.KTI(FIMP) ==> FKMP=', i6, ' FZLOC=', f10.3,   &
                       &' KTI(FIMP)=', i6, ' and (H(KTWBF)-Z(FIMP))=', f10.3)
             endif
 
!Check       for Boundary Violations: Lateral Direction Check (2 of 2)
 
             if(fyloc<0)then
                 fyloc = B(fkmp, fimp)*yrefl
             elseif(fyloc>B(fkmp, fimp))then
                 fyloc = B(fkmp, fimp)*(1 - yrefl)
             endif
 
         endif
 
 
50       call FINDNEWBR
         fishes(fn, 1) = REAL(fimp)   ! Segment IMP where fish is located
         fishes(fn, 2) = fxloc        ! Location of fish within segment IMP from upstream side
         fishes(fn, 3) = REAL(fkmp)   ! Layer KMP where fish is located
         fishes(fn, 4) = fzloc        ! Location of fish within layer KMP from top side
         fishes(fn, 5) = fyloc        ! Location of fish laterally from left bank (in plan view)
         fishes(fn, 6) = REAL(fnbp)   ! Branch where fish is located
         fishes(fn, 7) = fsize        ! Size (i.e., length) of fish in meters
         fishes(fn, 8) = fage + nfsfreq    ! Age of the fish  ! SW 2/01/01
         fishes(fn, 9) = ufish        ! Longitudinal velocity of fish relative to water
         fishes(fn, 10) = vfish       ! Lateral velocity of fish relative to water
         fishes(fn, 11) = wfish       ! Vertical velocity of fish relative to water
         fishes(fn, 12) = REAL(iswitch)
                                      ! Is fish still in system?: Yes=0  No=1
         fishes(fn, 13) = REAL(fsnag) ! Gillnet # fish is snagged in; = 0 if fish not in a gillnet
 
     enddo       !end of NFISH loop
 
 
     lrunday = jday
 
 
!*******************************************************
!    O U T P U T   A N D   E N D   S U B R O U T I N E   *
!*******************************************************
 
     if((fcount==1) .OR. (nit==0))call FISHPLOT                         ! This plots the location of virtual fish
 
     end subroutine FISH
