!*==initial_water_level.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
!***********************************************************************************************************************************
!**    S U B R O U T I N E    I N I T I A L    W A T E R    L E V E L                                                            **
!**********************************************************************************************************************************
 
     subroutine INITIAL_WATER_LEVEL
 
     use MAIN
     use GLOBAL
     use NAMESC
     use GEOMC
     use LOGICC
     use PREC
     use SURFHE
     use KINETIC
     use SHADEC
     use EDDY
     use STRUCTURES
     use TRANS
     use TVDC
     use SELWC
     use GDAYC
     use SCREENC
     use TDGAS
     use RSTART
     use MACROPHYTEC
     use POROSITYC
     use ZOOPLANKTONC
     use INITIALVELOCITY
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: brlen, dist, qgate, wldiff, wlslope, wsup
     integer :: jbd, jbu, jjw
     external RESTART_OUTPUT
!
!*** End of declarations rewritten by SPAG
!
 
     loop_branch = .FALSE.
 
!    estimating initial flows in each segment
     qssi = 0.0
 
!    first considering specified flows: upstream inflows, tributaries,
!    distributed tribs, structural withdrawals and withdrawals
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do i = iu, id
                 if(i==iu .AND. UP_FLOW(jb))qssi(i) = QIN(jb) + qssi(i)
!                if downstream of structure
          !IF (i == iu .and. UHS(JB) < 0) then
                 if(i==iu .AND. DAM_INFLOW(jb))then ! CB 4/27/2011
                     do jjw = 1, nwb
                         do jjb = BS(jjw), BE(jjw)
                             if(DS(jjb)==ABS(UHS(jb)))then
                                 do js = 1, NSTR(jjb)
                                     qssi(i) = qssi(i) + QSTR(js, jjb)
                                 enddo
                             endif
                         enddo
                     enddo
                 endif
 
                 if(tributaries)then
                     do jt = 1, ntr
                         if(ITR(jt)==i)qssi(i) = qssi(i) + QTR(jt)
                     enddo
                 endif
                 if(DIST_TRIBS(jb))qssi(i) = qssi(i) + QDTR(jb)                &
                  & /REAL(id - iu + 1)                ! SINCE INITIAL WL UNKNOWN, DISTRIBUTING FLOW EVENLY BTW. SEGS.
                 if(withdrawals)then
                     do jwd = 1, nwd
                         if(IWD(jwd)==i)qssi(i) = qssi(i) - QWD(jwd)
                     enddo
                 endif
                 if(i==id)then
                     do js = 1, NSTR(jb)
                         qssi(i) = qssi(i) - QSTR(js, jb)
                     enddo
                 endif
             enddo
         enddo
     enddo
 
!    INCLUDING FLOWS WITHIN BRANCH UPSTEAM OF SEGMENT
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do i = iu + 1, id
                 qssi(i) = qssi(i) + qssi(i - 1)
             enddo
         enddo
     enddo
 
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do i = iu, id
!                DETERMINING IF SEGMENT IS DOWNSTREAM INTERNAL HEAD BOUNDARY
!                OF ANOTHER *UPSTREAM* BRANCH, AND ADDING FLOW TO SEGMENT AND
!                SEGMENTS DOWNSTREAM
                 do jjw = 1, nwb
                     do jjb = BS(jjw), BE(jjw)
                         if(DHS(jjb)==i)then
                             do ii = i, id
                                 qssi(ii) = qssi(ii) + qssi(DS(jjb))
                             enddo
                         endif
                     enddo
                 enddo
!                DETERMINING IF SEGMENT IS DOWNSTREAM OF SPILLWAY BELOW
!                ANOTHER BRANCH
                 do js = 1, nsp
                     if(ESP(js)<EL(2, i))then
                                       ! DISREGARDING IF CREST ABOVE GRID
                         if(i==IDSP(js))then
                             do ii = i, id
                                 qssi(ii) = qssi(ii) + qssi(IUSP(js))
                             enddo
                         endif
                     endif
                 enddo
!                DETERMINING IF SEGMENT IS DOWNSTREAM OF GATE BELOW ANOTHER
!                BRANCH
                 do jg = 1, ngt
                     if(EGT(jg)<EL(2, i))then
                                       ! DISREGARDING IF CREST ABOVE GRID
                         if(i==IDGT(jg))then
                             if(DYNGTC(jg)=='    FLOW')then
                                 qgate = BGT(jg)
                                 do ii = i, id
                                     qssi(ii) = qssi(ii) + qgate
                                 enddo
                             else
                                 do ii = i, id
                                     qssi(ii) = qssi(ii) + qssi(IUGT(jg))
                                 enddo
                             endif
                         endif
                     endif
                 enddo
             enddo
         enddo
     enddo
 
!    DETERMINING IF BRANCH IS A SECONDARY BRANCH THAT LOOPS AROUND AN ISLAND
!    WITH INTERNAL HEAD BOUNDARIES AT UPSTREAM AND DOWNSTREAM BC ATTACHED TO A
!    SINGLE BRANCH
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             if(UH_INTERNAL(jb) .AND. DH_INTERNAL(jb))then
                 do jjw = 1, nwb
                     do jjb = BS(jjw), BE(jjw)
                         if(UHS(jb)>US(jjb) .AND. UHS(jb)<DS(jjb))jbu = jjb
                     enddo
                 enddo
                 do jjw = 1, nwb
                     do jjb = BS(jjw), BE(jjw)
                         if(DHS(jb)>US(jjb) .AND. DHS(jb)<DS(jjb))jbd = jjb
                                                                        ! WW 8/19/2013
                     enddo
                 enddo
                 loop_branch(jb) = jbu==jbd
             endif
         enddo
     enddo
 
 
!    GIVEN ESTIMATED FLOWS FOR EACH SEGMENT, ESTIMATING WL WITH NORMAL DEPTH
!    EQUATION
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
!            ONLY CONSIDERING BRANCHES WITH SLOPES > 0
             if(SLOPE(jb)>0.0 .AND. .NOT.loop_branch(jb))then
                 iu = CUS(jb)
                 id = DS(jb)
                 do i = iu, id
                     call NORMAL_DEPTH(qssi(i))
                 enddo
             endif
         enddo
     enddo
 
!    SMOOTHING WATER SURFACE ELEVATIONS, FIRST WITHIN BRANCHES
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             if(SLOPE(jb)>0.0 .AND. .NOT.loop_branch(jb))then
                 iu = CUS(jb)
                 id = DS(jb)
                 do i = id - 1, iu, -1
                     if(ELWS(i + 1)>ELWS(i))ELWS(i) = ELWS(i + 1)
                 enddo
             endif
         enddo
     enddo
 
!    SMOOTHING WATER SURFACE ELEVATIONS AT INTERNAL HEAD BOUNDARIES
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             if(SLOPE(jb)>0.0 .AND. .NOT.loop_branch(jb) .AND. DH_INTERNAL(jb))&
              & then
                 iu = CUS(jb)
                 id = DS(jb)
                 if(ELWS(DHS(jb))>ELWS(id))then
                     ELWS(id) = ELWS(DHS(jb))
                     do i = id - 1, iu, -1
                         if(ELWS(i + 1)>ELWS(i))ELWS(i) = ELWS(i + 1)
                     enddo
                 endif
             endif
         enddo
     enddo
 
!    IF SPILLWAY OR GATE AT DOWNSTREAM END OF BRANCH, MAKING SURE WATER LEVEL
!    IS ABOVE CREST ELEVATION
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             if(SLOPE(jb)>0.0 .AND. .NOT.loop_branch(jb) .AND. NSTR(jb)==0)then
                 iu = CUS(jb)
                 id = DS(jb)
!                SPILLWAYS
                 do js = 1, nsp
                     if(ESP(js)<EL(2, i))then
                                       ! DISREGARDING IF CREST ABOVE GRID
                         if(id==IUSP(js))then
                             wsup = ESP(js) + (qssi(id)/A1SP(js))              &
                                  & **(1.0/B1SP(js))               ! ESTIMATING UPSTREAM WS ELEV
                             if(ELWS(id)<wsup)then
                                 if(IDSP(js)/=0)then
                                              ! CB 8/10/10
                                     if(ELWS(IDSP(js))>wsup)                   &
                                      & wsup = ELWS(IDSP(js))        ! CHECKING TO SEE IF DOWNSTREAM WS ELEVATION ISN'T ALREADY 'HIGH'
                                     ELWS(id) = wsup
                                 endif        ! CB 8/10/10
                                 do i = id - 1, iu, -1
                                     if(ELWS(i + 1)>ELWS(i))ELWS(i)            &
                                      & = ELWS(i + 1)
                                 enddo
                             endif
                         endif
                     endif
                 enddo
!                GATES
                 do jg = 1, ngt
                     if(EGT(jg)<EL(2, i) .AND. DYNGTC(jg)/='    FLOW')then
                                                                        ! DISREGARDING IF CREST ABOVE GRID
                         if(id==IUGT(jg))then
                             if(DYNGTC(jg)=='       B')wsup = EGT(jg)          &
                              & + (qssi(id)/(A1GT(jg)*BGT(jg)**G1GT(jg)))      &
                              & **(1.0/B1GT(jg))
                             if(DYNGTC(jg)=='     ZGT')wsup = EGT(jg)          &
                              & + (qssi(id)/A1GT(jg))**(1.0/B1GT(jg))
                             if(ELWS(id)<wsup)then
                                 if(ELWS(IDGT(jg))>wsup)wsup = ELWS(IDGT(jg))
                                                                     ! CHECKING TO SEE IF DOWNSTREAM WS ELEVATION ISN'T ALREADY 'HIGH'  WX 8/21/13
                                 ELWS(id) = wsup
                                 do i = id - 1, iu, -1
                                     if(ELWS(i + 1)>ELWS(i))ELWS(i)            &
                                      & = ELWS(i + 1)
                                 enddo
                             endif
                         endif
                     endif
                 enddo
             endif
         enddo
     enddo
 
!    SMOOTHING WATER LEVEL AROUND LOOP BRANCHES
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             if(loop_branch(jb))then
                 wldiff = ELWS(UHS(jb)) - ELWS(UHS(jb))
                 brlen = 0.0
                 do i = iu, id
                     brlen = brlen + DLX(i)
                 enddo
                 wlslope = wldiff/brlen
                 dist = DLX(iu)/2.0
                 ELWS(iu) = ELWS(UHS(jb)) - wlslope*dist
                 do i = iu + 1, id
                     dist = dist + (DLX(i - 1) + DLX(i))/2.0
                     ELWS(iu) = ELWS(UHS(jb)) - wlslope*dist
                 enddo
             endif
         enddo
     enddo
 
     end subroutine INITIAL_WATER_LEVEL
