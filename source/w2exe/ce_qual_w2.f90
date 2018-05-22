!*==ce_qual_w2.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
! CE-QUAL-W2 computations
     function CE_QUAL_W2(Dlg)
 
! IVF/CVF specific code
 
 
! IVF/CVF specific code
 
 
! IVF/CVF specific code
 
 
 !DEC$ATTRIBUTES STDCALL   :: ce_qual_w2
 !DEC$ATTRIBUTES REFERENCE :: Dlg
 !SP CEMA
 !End SP CEMA
 !End SP CEMA
!  include "omp_lib.h"      ! OPENMP directive to adjust the # of processors TOGGLE FOR DEBUG
     use DFLOGM              ! to get current working directory
     use MSCLIB
     use DFWIN, RENAMED => dlt
     use IFPORT
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
     use CEMAVARS
     use INITIALVELOCITY
     use ENVIRPMOD
     use BIOENERGETICS
     use BUILDVERSION
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     type(DIALOG) :: Dlg
     integer :: CE_QUAL_W2
!
! Local variables
!
     real :: depth
     character(255) :: dirc
     integer(4) :: istatus, length
     integer :: result
!
!*** End of declarations rewritten by SPAG
!
 
 
!***********************************************************************************************************************************
!**  Task 1: Inputs                                                          **
!***********************************************************************************************************************************
 
!    call omp_set_num_threads(4)   ! set # of processors to NPROC  Moved to
 
!    INPUT subroutine
     if(end_run .OR. error_open)stop
                                 ! SW 6/26/15 3/18/16 Added code to prevent a thread from reinitializing output files as dialog box is closing...intermittant error Updated 8/23/2017
 
     call GET_COMMAND_ARGUMENT(1, dirc, length, istatus)
     dirc = TRIM(dirc)
 
!    IF(ISTATUS.NE.0)WRITE(*,*)'GET_COMMAND_ARGUMENT FAILED: STATUS=',ISTATUS
 
     if(length/=0)then
         istatus = CHDIR(dirc)
         selectcase(istatus)
         case(2)
               ! ENOENT
             write(w2err, *)'The directory does not exist:', dirc
             write(w2err, *)'Run stopped'
         case(20)
                 ! ENOTDIR
             write(w2err, *)'This is not a directory:', dirc
             write(w2err, *)'Run stopped'
         case(0) ! NO ERROR
         endselect
     endif
 
     moddir = file$curdrive         !  GET CURRENT DIRECTORY
     length = GETDRIVEDIRQQ(moddir)
 
     open(10, file = 'W2CodeCompilerVersion.opt', status = 'unknown')
     write(10, '(A,F5.2)')' CE-QUAL-W2 Version #:', w2ver
     write(10, *)'Compiler Version and Code Compile Date'
     write(10, *)'INTEL_COMPILER_VERSION:', intel_compiler_version
     write(10, *)'INTEL_COMPILER_BUILD_DATE:', intel_compiler_build_date
     write(10, *)'CE-QUAL-W2 Version compile date:', buildtime
     close(10)
 
!    open ancillary control file for fish bioenergetics !mlm bioexp
     bioexp = .FALSE.           ! INITIALIZE LOGICAL VARIABLE THIS IS READ IN THE CONTROL FILE
     fishbio = .FALSE.
     inquire(file = 'W2_con_anc.npt', exist = fishbio)
                                                   ! SW 5/26/15
     if(fishbio)then
 
         open(fishbiofn, file = 'W2_con_anc.npt', status = 'old')
         do ii = 1, 16
             read(fishbiofn, '(A8)')bioc
                                ! DUMMY VARIABLE AT THIS POINT FIX THIS
         enddo
     endif
 
     fish_particle_exist = .FALSE.
     inquire(file = 'particle.csv', exist = fish_particle_exist)
                                                            ! SW 4/30/15
 
!    Open control file
     iopenfish = 0
     open(con, file = confn, status = 'OLD', iostat = i)
     if(i/=0)then
         text = 'Could not open w2_con.npt'
         goto 300
     endif
 
     call INPUT
!SP  CEMA
     call CEMA_W2_INPUT
!End SP CEMA
     restart_in = rsic=='      ON' .OR. restart_pushed
!    Restart data
     if(restart_pushed)rsifn = 'rso.opt'
 
     jday = tmstrt
     if(restart_in)then
         rsi = 10
             ! SW 5/26/15
         vert_profile = .FALSE.
         long_profile = .FALSE.
         open(rsi, file = rsifn, form = 'UNFORMATTED', status = 'OLD')
         read(rsi)nit, nv, kmin, imin, nsprf, cmbrt, zmin, izmin, start,       &
                & current
         read(rsi)dltdp, snpdp, tsrdp, vpldp, prfdp, cpldp, sprdp, rsodp,      &
                & scrdp, flxdp, wdodp
         read(rsi)jday, eltm, eltmf, dlt, dltav, dlts, mindlt, jdmin, curmax
         read(rsi)nxtmsn, nxtmts, nxtmpr, nxtmcp, nxtmvp, nxtmrs, nxtmsc,      &
                & nxtmsp, nxtmfl, nxtmwd
         read(rsi)volin, volout, voluh, voldh, volpr, voltrb, voldt, volwd,    &
                & volev, volsbr, voltr, volsr, volice, icebank
         read(rsi)tssev, tsspr, tsstr, tssdt, tsswd, tssin, tssout, tsss, tssb,&
                & tssice
         read(rsi)tssuh, tssdh, tssuh2, tssdh2, cssuh2, cssdh2, voluh2, voldh2,&
                & quh1
         read(rsi)esbr, etbr, ebri
         read(rsi)z, sz, elws, savh2, savhr, h2
         read(rsi)ktwb, kti, skti, sbkt
         read(rsi)ice, iceth, cuf, qsum
         read(rsi)u, w, su, sw, az, saz, dltlim
         read(rsi)t1, t2, c1, c2, c1s, sed, kfs, cssk
         read(rsi)epd, epm
         read(rsi)macmbrt, macrc, smacrc, mac, smac, macrm, macss
         read(rsi)sedc, sedn, sedp, zoo, cd
                                           ! mlm 10/06
         read(rsi)sdkv                     ! cb 11/30/06
         read(rsi)tke                      ! sw 10/4/07
         if(envirpc=='      ON')then
 
             allocate(cc_e(nct), c_int(nct), c_top(nct), cd_e(ndc), cd_int(ndc)&
                    & , cd_top(ndc), c_avg(nct), cd_avg(ndc), cn_e(nct),       &
                    & cdn_e(ndc))
             cc_e = '   '
             c_int = 0.0
             c_top = 0.0
             cd_e = '   '
             cd_int = 0.0
             cd_top = 0.0
             c_avg = 0.0
             cd_avg = 0.0
             cn_e = 0.0
             cdn_e = 0.0
             nac_e = 0
             nacd_e = 0
             open(cone, file = 'w2_envirprf.npt', status = 'old')
             read(cone, '(//I1,7X,I8,5x,a3,f8.0,f8.0,9(i8,i8))')i_segint,      &
                & numclass, selectivec, sjday1, sjday2, istart(1), iend(1),    &
                & (istart(i), iend(i), i = 2, i_segint)
             if(i_segint==0)i_segint = 1
             read(cone, '(//8x,3(5x,a3,f8.3,f8.3))')vel_vpr, vel_int, vel_top, &
                & temp_vpr, temp_int, temp_top, depth_vpr, d_int, d_top
             read(cone, '(//(8X,(5X,A3,F8.0,F8.0)))')                          &
                & (cc_e(jc), c_int(jc), c_top(jc), jc = 1, nct)
             read(cone, '(//(8X,(5X,A3,F8.0,F8.0)))')                          &
                & (cd_e(jd), cd_int(jd), cd_top(jd), jd = 1, ndc)
             close(cone)
 
             do jc = 1, nct
                 if(cc_e(jc)==' ON')then
                     nac_e = nac_e + 1
                     cn_e(nac_e) = jc
                 endif
             enddo
             do jd = 1, ndc
                 if(cd_e(jd)==' ON')then
                     nacd_e = nacd_e + 1
                     cdn_e(nacd_e) = jd
                 endif
             enddo
 
             allocate(c_cnt(i_segint, nct), cd_cnt(i_segint, ndc),             &
                    & c_class(i_segint, nct, numclass),                        &
                    & cd_class(i_segint, ndc, numclass), c_tot(i_segint, nct), &
                    & cd_tot(i_segint, ndc), t_class(i_segint, numclass),      &
                    & v_class(i_segint, numclass), c_sum(nct), cd_sum(ndc))
             allocate(conc_c(nct, numclass), conc_cd(ndc, numclass))
             allocate(d_class(i_segint, numclass))
             allocate(d_tot(i_segint), d_cnt(i_segint), t_tot(i_segint),       &
                    & t_cnt(i_segint))
             allocate(v_tot(i_segint), v_cnt(i_segint), volgl(i_segint),       &
                    & sumvolt(i_segint))
 
             read(rsi)t_class, v_class, c_class, cd_class, t_tot, t_cnt,       &
                    & sumvolt, v_cnt, v_tot, c_tot, c_cnt, cd_tot, cd_cnt
 
         endif
         if(npi>0)read(rsi)ys, vs, vst, yst, dtp, qold
         read(rsi)tpout, tptrib, tpdtrib, tpwd, tppr, tpin, tp_sedsod_po4,     &
                & pfluxin, tnout, tntrib, tndtrib, tnwd, tnpr, tnin,           &
                & tn_sedsod_nh4, nfluxin                                                                                            ! TP_SEDBURIAL,TN_SEDBURIAL,
         close(rsi)
     endif
     CE_QUAL_W2 = 1
     call INIT
 
!    determining initial horizontal velocities and water levels
     once_through = .TRUE.
     if(inituwl=='      ON')init_vel = .TRUE.
     if(.NOT.restart_in)then
         if(init_vel)then
             allocate(qssi(imx), loop_branch(nbr), elwss(imx), uavg(imx))
             elwss = elws
             call INITIAL_WATER_LEVEL
             b = bsave
             call INITGEOM
             call INITIAL_U_VELOCITY
             open(nunit, file = 'init_wl_u_check.dat', status = 'unknown')
             write(nunit,                                                      &
                  &'("       i elws_calc    qssi       u   depth elws_init")')
             do jw = 1, nwb
                 do jb = BS(jw), BE(jw)
                     iu = CUS(jb)
                     id = DS(jb)
                     do i = iu, id
                         depth = elws(i) - EL(KBI(i) + 1, i)
                         write(nunit, '(i8,f8.3,f8.2,f8.3,2f8.2)')i, elws(i),  &
                             & qssi(i), u(kt, i), depth, elwss(i)
                     enddo
                 enddo
             enddo
             close(nunit)
             deallocate(qssi, loop_branch, elwss)
         endif
     endif
 
     if(.NOT.restart_in)then
         line = cctime(1:2) // ':' // cctime(3:4) // ':' // cctime(5:6)
         result = DLGSET(Dlg, starting_time, TRIM(line))                                              !Display starting time
         result = DLGSET(Dlg, status, 'Executing')                                                    !Display execution status
         current = 0.0
     else
         call CPU_TIME(current)
     endif
 
     call OUTPUTINIT
     if(restart_in)then
         do jw = 1, nwb
             if(SCREEN_OUTPUT(jw))call SCREEN_UPDATE(Dlg)
         enddo
     endif
 
     if(.NOT.restart_in)call CPU_TIME(start)
 
     if(macrophyte_on .AND. constituents)call POROSITY
     if(selectc=='      ON')call SELECTIVEINIT  ! new subroutine for selecting water temperature target
     if(selectc=='    USGS')call SELECTIVEINITUSGS  ! new subroutine for selecting water temperature target
     if(aeratec=='      ON' .AND. oxygen_demand)call AERATE
 !SP CEMA
     if(sediment_diagenesis)then
         if(cemarelatedcode .AND. includebedconsolidation)                     &
          & call SETUPCEMASEDIMENTMODEL
         if(includefftlayer)call CEMAFFTLAYERCODE
     endif
 !End SP CEMA
 
 ! CEMA testing start
  ! open(1081,file='sod_ss_testing.dat',status='unknown')
  ! write(1081,'("   xjnh4      Jc     Jcs     SOD")')
 ! CEMA testing end
 
!***********************************************************************************************************************************
!**  Task 2: Calculations                                                     
!***********************************************************************************************************************************
!    **
     do while (.NOT.end_run .AND. .NOT.stop_pushed)
         if(jday>=nxtvd)call READ_INPUT_DATA(nxtvd)
         call INTERPOLATE_INPUTS
         dlttvd = (nxtvd - jday)*day
         dlt = MIN(dlt, dlttvd + 1.0)
         dlts1 = dlt
         if(dlt<=dlttvd + 0.999)then
             dlts = dlt
         else
             kloc = 1
             iloc = 1
         endif
 
    ! update wind at 2m for evaopration and evaoprative heat flux computations  ! SW 5/21/15
         do jw = 1, nwb
             do i = CUS(BS(jw)), DS(BE(jw))
                 WIND2(i) = WIND(jw)*WSC(i)*DLOG(2.0D0/Z0(jw))                 &
                          & /DLOG(WINDH(jw)/Z0(jw))
             enddo
         enddo
50       if(selectc=='      ON')call SELECTIVE
                                           ! new subroutine for selecting water temperature target
         if(selectc=='    USGS')call SELECTIVEUSGS
                                               ! new subroutine for selecting water temperature target
 
         call HYDROINOUT
 
!SP      CEMA
         if(sediment_diagenesis)then
             if(cemarelatedcode .AND. includebedconsolidation)                 &
              & call COMPUTECEMARELATEDSOURCESINKS
  !  If(CEMARelatedCode .and. IncludeCEMASedDiagenesis)Call ComputeCEMADiagenesisSourceSinks  SW 6/27/2017
         endif
!End     SP CEMA
 
!***********************************************************************************************************************************
!**      Task 2.2: Hydrodynamic calculations                                  
!***********************************************************************************************************************************
 
!        **
         do jw = 1, nwb
             kt = ktwb(jw)
             do jb = BS(jw), BE(jw)
                 if(BR_INACTIVE(jb))cycle
                                    ! SW 6/12/17
                 iu = CUS(jb)
                 id = DS(jb)
 
!***********************************************************************************************************************************
!**              Task 2.2.1: Boundary concentrations, temperatures, and
!***********************************************************************************************************************************
 
!                densities                               **
                 iut = iu
                 idt = id
                 if(UP_FLOW(jb))then
                     if(.NOT.INTERNAL_FLOW(jb))then
                         do k = kt, KB(iu)
                             if(QIND(jb) + QINSUM(jb)>0.0)then
                                 TIN(jb)                                       &
                                   & = (TINSUM(jb)*QINSUM(jb) + TIND(jb)*QIND  &
                                   & (jb))/(QIND(jb) + QINSUM(jb))
                                 CIN(cn(1:nac), jb)                            &
                                   & = MAX((CINSUM(cn(1:nac), jb)*QINSUM(jb)   &
                                   & + CIND(cn(1:nac), jb)*QIND(jb))           &
                                   & /(QIND(jb) + QINSUM(jb)), 0.0)
                                 t1(k, iu - 1) = TIN(jb)
                                 t2(k, iu - 1) = TIN(jb)
                                 c1s(k, iu - 1, cn(1:nac)) = CIN(cn(1:nac), jb)
                                 QIN(jb) = QIND(jb) + QINSUM(jb)
                             else
                                 QIN(jb) = 0.0
                                 TIN(jb) = TIND(jb)
                                 t1(k, iu - 1) = TIND(jb)
                                 t2(k, iu - 1) = TIND(jb)
                                 c1s(k, iu - 1, cn(1:nac))                     &
                                   & = CIND(cn(1:nac), jb)
                             endif
                         enddo
                     elseif(.NOT.DAM_INFLOW(jb))then                                                                   !TC 08/03/04
                         if(JBUH(jb)>=BS(jw) .AND. JBUH(jb)<=BE(jw))then
                             TIN(jb) = t1(kt, UHS(jb))
                             CIN(cn(1:nac), jb)                                &
                               & = MAX(c1s(kt, UHS(jb), cn(1:nac)), 0.0)
                             do k = kt, KB(iu)
                                !CONCURRENT(K=KT:KB(IU))   !FORALL(K=KT:KB(IU))           !DO K=KT,KB(IU)
                                 t1(k, iu - 1) = t1(k, UHS(jb))
                                 t2(k, iu - 1) = t1(k, UHS(jb))
                                 c1s(k, iu - 1, cn(1:nac))                     &
                                   & = c1s(k, UHS(jb), cn(1:nac))
                                 c1(k, iu - 1, cn(1:nac))                      &
                                   & = c1s(k, UHS(jb), cn(1:nac))
                                 c2(k, iu - 1, cn(1:nac))                      &
                                   & = c1s(k, UHS(jb), cn(1:nac))
                             enddo
                         else
                             call UPSTREAM_WATERBODY
                             TIN(jb) = t1(kt, iu - 1)
                             CIN(cn(1:nac), jb)                                &
                               & = MAX(c1(kt, iu - 1, cn(1:nac)), 0.0)
                         endif
                     else
                         TIN(jb) = TINSUM(jb)
                         QIN(jb) = QINSUM(jb)
                         CIN(cn(1:nac), jb) = MAX(CINSUM(cn(1:nac), jb), 0.0)
                         do k = kt, KB(id)
                             t1(k, iu - 1) = TIN(jb)
                             t2(k, iu - 1) = TIN(jb)
                             c1s(k, iu - 1, cn(1:nac)) = CIN(cn(1:nac), jb)
                         enddo
                     endif
                 endif
                 if(DN_FLOW(jb))then
                     do k = kt, KB(id)
                         t1(k, id + 1) = t2(k, id)
                         t2(k, id + 1) = t2(k, id)
                         c1s(k, id + 1, cn(1:nac)) = c1s(k, id, cn(1:nac))
                     enddo
                 endif
                 if(UP_HEAD(jb))then
                     iut = iu - 1
                     if(UH_INTERNAL(jb))then
                         if(JBUH(jb)>=BS(jw) .AND. JBUH(jb)<=BE(jw))then
                             do k = kt, KB(iut)
                                 RHO(k, iut) = RHO(k, UHS(jb))
                                 t1(k, iut) = t2(k, UHS(jb))
                                 t2(k, iut) = t2(k, UHS(jb))
                                 c1s(k, iut, cn(1:nac))                        &
                                   & = c1s(k, UHS(jb), cn(1:nac))
                                 c1(k, iut, cn(1:nac))                         &
                                   & = c1s(k, UHS(jb), cn(1:nac))
                                 c2(k, iut, cn(1:nac))                         &
                                   & = c1s(k, UHS(jb), cn(1:nac))
                             enddo
                         else
                             call UPSTREAM_WATERBODY
                         endif
                         do k = kt, KB(iut)
                             RHO(k, iut)                                       &
                               & = DENSITY(t2(k, iut), DMAX1(TDS(k, iut),      &
                               & 0.0D0), DMAX1(TISS(k, iut), 0.0D0))
                         enddo
                     elseif(UH_EXTERNAL(jb))then
                         do k = kt, KB(iut)
                             RHO(k, iut)                                       &
                               & = DENSITY(TUH(k, jb), DMAX1(TDS(k, iut),      &
                               & 0.0D0), DMAX1(TISS(k, iut), 0.0D0))
                             t1(k, iut) = TUH(k, jb)
                             t2(k, iut) = TUH(k, jb)
                             c1s(k, iut, cn(1:nac)) = CUH(k, cn(1:nac), jb)
                             c1(k, iut, cn(1:nac)) = CUH(k, cn(1:nac), jb)
                             c2(k, iut, cn(1:nac)) = CUH(k, cn(1:nac), jb)
                         enddo
                     endif
                 endif
                 if(DN_HEAD(jb))then
                     idt = id + 1
                     if(DH_INTERNAL(jb))then
                         if(JBDH(jb)>=BS(jw) .AND. JBDH(jb)<=BE(jw))then
                             do k = kt, KB(idt)
                                 RHO(k, idt) = RHO(k, DHS(jb))
                                 t1(k, idt) = t2(k, DHS(jb))
                                 t2(k, idt) = t2(k, DHS(jb))
                                 c1s(k, idt, cn(1:nac))                        &
                                   & = c1s(k, DHS(jb), cn(1:nac))
                                 c1(k, idt, cn(1:nac))                         &
                                   & = c1s(k, DHS(jb), cn(1:nac))
                                 c2(k, idt, cn(1:nac))                         &
                                   & = c1s(k, DHS(jb), cn(1:nac))
                             enddo
                         else
                             call DOWNSTREAM_WATERBODY
                         endif
                         do k = kt, KB(id)
                             RHO(k, idt)                                       &
                               & = DENSITY(t2(k, idt), DMAX1(TDS(k, idt),      &
                               & 0.0D0), DMAX1(TISS(k, idt), 0.0D0))
                         enddo
                     elseif(DH_EXTERNAL(jb))then
                         do k = kt, KB(idt)
                             RHO(k, idt)                                       &
                               & = DENSITY(TDH(k, jb), DMAX1(TDS(k, idt),      &
                               & 0.0D0), DMAX1(TISS(k, idt), 0.0D0))
                             t1(k, idt) = TDH(k, jb)
                             t2(k, idt) = TDH(k, jb)
                             c1s(k, idt, cn(1:nac)) = CDH(k, cn(1:nac), jb)
                             c1(k, idt, cn(1:nac)) = CDH(k, cn(1:nac), jb)
                             c2(k, idt, cn(1:nac)) = CDH(k, cn(1:nac), jb)
                         enddo
                     endif
                 endif
 
!***********************************************************************************************************************************
!**              Task 2.2.2: Momentum terms                                   
!***********************************************************************************************************************************
 
!******          ** Density pressures
 
                 do i = iut, idt
                     do k = kt, KB(i)
                         P(k, i) = P(k - 1, i) + RHO(k, i)*g*H(k, jw)*COSA(jb)
                     enddo
                 enddo
 
!******          Horizontal density gradients
 
                 do i = iut, idt - 1
                     HDG(kt, i) = DLXRHO(i)*(BKT(i) + BKT(i + 1))              &
                                & *0.5D0*H(kt, jw)*(P(kt, i + 1) - P(kt, i))
                     do k = kt + 1, KBMIN(i)
                         HDG(k, i) = DLXRHO(i)*BHR2(k, i)                      &
                                   & *((P(k - 1, i + 1) - P(k - 1, i))         &
                                   & + (P(k, i + 1) - P(k, i)))
                     enddo
                 enddo
 
!******          Adjusted wind speed and surface wind shear drag coefficient
 
                 do i = iu - 1, id + 1
                     WIND10(i) = WIND(jw)*WSC(i)*DLOG(10.0D0/Z0(jw))           &
                               & /DLOG(WINDH(jw)/Z0(jw))                             ! older  version z0=0.01                      ! SW 11/28/07
                     FETCH(i) = FETCHD(i, jb)
                     if(COS(PHI(jw) - PHI0(i))<0.0)FETCH(i) = FETCHU(i, jb)
                     if(FETCH(i)<=0.0)FETCH(i) = DLX(i)
                     if(FETCH_CALC(jw))then
                         zb = 0.8D0*DLOG(FETCH(i)*0.5D0) - 1.0718D0
                         WIND10(i) = WIND10(i)*(5.0D0*zb + 4.6052D0)           &
                                   & /(3.0D0*zb + 9.2103D0)
                     endif
 
                     if(WIND10(i)>=15.0)then            ! SW 1/19/2008
                         CZ(i) = 0.0026D0
                     elseif(WIND10(i)>=4.0)then
                         CZ(i) = 0.0005D0*DSQRT(WIND10(i))
                     elseif(WIND10(i)>=0.5)then
                         CZ(i) = 0.0044D0*WIND10(i)**( - 1.15D0)
                     else
                         CZ(i) = 0.01D0
                     endif
 
  !        CZ(I) = 0.0
  !        IF (WIND10(I) >= 1.0)  CZ(I) = 0.0005*SQRT(WIND10(I))
  !        IF (WIND10(I) >= 4.0) CZ(I) = 0.0005*SQRT(WIND10(I))
  !        IF (WIND10(I) >= 15.0) CZ(I) = 0.0026
                 enddo
 
!******          Longitudinal and lateral surface wind shear and exponential
 
!                decay
                 do i = iut, idt - 1
          !WSHX(I) = CZ(I)*WIND10(I)**2*RHOA/RHOW*    DCOS(PHI(JW)-PHI0(I))* ICESW(I)
          !WSHY(I) = CZ(I)*WIND10(I)**2*RHOA/RHOW*DABS(DSIN(PHI(JW)-PHI0(I)))*ICESW(I)
                     WSHX(i) = CZ(i)*WIND10(i)*WIND10(i)                       &
                             & *rhoa/rhow*DCOS(PHI(jw) - PHI0(i))*ICESW(i)                     ! SW 4/20/16 SPEED
                     WSHY(i) = CZ(i)*WIND10(i)*WIND10(i)                       &
                             & *rhoa/rhow*DABS(DSIN(PHI(jw) - PHI0(i)))        &
                             & *ICESW(i)
                     wwt = 0.0
                     if(WIND10(i)/=0.0)wwt = 6.95D-2*(FETCH(i)**0.233D0)       &
                      & *WIND10(i)**0.534D0
                     dfc = -8.0D0*pi*pi/(g*wwt*wwt + nonzero)
                     do k = kt, KBMIN(i)
                         DECAY(k, i) = DEXP(DMAX1(dfc*DEPTHB(k, i), -30.0D0))
                     enddo
 
!********            Branch inflow lateral shear and friction
 
                     do jjb = 1, nbr
                         if(BR_INACTIVE(jjb))cycle
                                       ! SW 6/12/2017
                         if(i==UHS(jjb) .AND. .NOT.INTERNAL_FLOW(jjb))then
                             betabr = (PHI0(i) - PHI0(US(jjb)))
                             if(jjb>=BS(jw) .AND. jjb<=BE(jw))then
                                 do k = kt, KBMIN(i)
                                     if(u(k, US(jjb))<0.0)then
                                         UXBR(k, i) = UXBR(k, i)               &
                                           & + ABS(u(k, US(jjb)))*DCOS(betabr) &
                                           & *voluh2(k, jjb)/(dlt*DLX(i))
                                         UYBR(k, i) = UYBR(k, i)               &
                                           & + ABS(DSIN(betabr))               &
                                           & *ABS(voluh2(k, jjb))/dlt
                                     endif
                                 enddo
                             else
                                 call UPSTREAM_BRANCH
                             endif
                         endif
                         if(i==DHS(jjb))then
                             betabr = (PHI0(i) - PHI0(DS(jjb)))
                             if(i==US(jb) .AND. UHS(jb)/=DS(jjb))then
                                 if(jjb>=BS(jw) .AND. jjb<=BE(jw))then
                                     do k = kt, KBMIN(i)
                                         if(u(k, DS(jjb))>=0.0)then
                                         UXBR(k, i) = UXBR(k, i)               &
                                           & + u(k, DS(jjb))*DCOS(betabr)      &
                                           & *voldh2(k, jjb)/(dlt*DLX(i))
                                         UYBR(k, i) = UYBR(k, i)               &
                                           & + ABS(DSIN(betabr))*voldh2(k, jjb)&
                                           & /dlt
                                         endif
                                     enddo
                                 else
                                     call DOWNSTREAM_BRANCH
                                 endif
                             elseif(i/=US(jb))then
                                 if(jjb>=BS(jw) .AND. jjb<=BE(jw))then
                                     do k = kt, KBMIN(i)
                                         if(u(k, DS(jjb))>=0.0)then
                                         UXBR(k, i) = UXBR(k, i)               &
                                           & + u(k, DS(jjb))*DCOS(betabr)      &
                                           & *voldh2(k, jjb)/(dlt*DLX(i))
                                         UYBR(k, i) = UYBR(k, i)               &
                                           & + ABS(DSIN(betabr))*voldh2(k, jjb)&
                                           & /dlt
                                         endif
                                     enddo
                                 else
                                     call DOWNSTREAM_BRANCH
                                 endif
                             endif
                         endif
                     enddo
                     do k = kt, KBMIN(i)
                         FRICBR(k, i) = (FI(jw)/8.0D0)*RHO(k, i)               &
                                      & *(UYBR(k, i)/(DLX(i)*h2(k, i)))**2
                     enddo
                 enddo
 
!******          Vertical eddy viscosities/diffusivities
                 FIRSTI(jw) = iut
                 LASTI(jw) = idt
                 do i = iut, idt - 1
                     call CALCULATE_AZ
          !SP CEMA
                     if(sediment_diagenesis)then
                         if(cemarelatedcode .AND.                              &
                          & includecemaseddiagenesis .AND. applybubbturb)      &
                          & call CEMABUBBLESTURBULENCE
                     endif
          !End SP CEMA
                     if(KBMIN(i)<=kt + 1 .AND. KB(i)>KBMIN(i))then
                         az(KBMIN(i), i) = azmin
                         DZ(KBMIN(i), i) = dzmin
                     endif
                 enddo
                 if(AZC(jw)=='     TKE' .OR. AZC(jw)=='    TKE1')then
                     azt(:, idt - 1) = az(:, idt - 1)
                     do i = iut, idt - 2
                         do k = kt, KBMIN(i)
                             azt(k, i) = 0.5*(az(k, i) + az(k, i + 1))
                         enddo
                         az(KBMIN(i), i) = azmin
                                              !SG 10/4/07 SW 10/4/07
                     enddo
                     az(kt:kmx - 1, iut:idt - 1) = azt(kt:kmx - 1, iut:idt - 1)
                 endif
                 do jwr = 1, niw
                     if(weir_calc)az(KTWR(jwr) - 1:KBWR(jwr), iwr(1:niw)) = 0.0
                 enddo
 
!******          Average eddy diffusivities
 
                 if(AZC(jw)/='     TKE' .AND. AZC(jw)/='    TKE1')then
                     DZ(kt:KB(idt) - 1, idt) = dzt(kt:KB(idt) - 1, idt - 1)
                                                          ! DZT is only used for non-TKE algorithms
                 else
                     DZ(kt:KB(idt) - 1, idt) = DZ(kt:KB(idt) - 1, idt - 1)
                 endif
                 do i = iut, idt - 1
                     do k = kt, KB(i) - 1
                         if(k>=KBMIN(i))then
                             if(KB(i - 1)>=KB(i) .AND. i/=iut)then
                                 DZ(k, i) = DZ(k, i - 1)
                             else
                                 DZ(k, i) = dzmin
                             endif
                         elseif(AZC(jw)/='     TKE' .AND. AZC(jw)/='    TKE1') &
                              & then
                             if(i==iut)then                   ! SW 10/20/07
                                 DZ(k, i) = dzt(k, i)
                             else
                                 DZ(k, i) = (dzt(k, i) + dzt(k, i - 1))*0.5D0
                                                                 ! SW 10/20/07  (DZT(K,I)+DZT(K+1,I))*0.5 ! FOR NON-TKE ALGORITHMS, AVERAGE DZ FROM EDGES TO CELL CENTERS
                             endif
                         endif
                     enddo
                 enddo
 
!                Hypolimnetic aeration
 
                 if(aeratec=='      ON' .AND. oxygen_demand)call DZAERATE
 
!******          Density inversions
 
                 do i = iut, idt
                     do k = kt, KB(i) - 1
                         DZQ(k, i) = MIN(1.0D-2, DZ(k, i))                    !MIN(1.0E-4,DZ(K,I)) No reason to limit DZ in rivers/estuaries-used in ULTIMATE scheme
                         if(RHO(k, i)>RHO(k + 1, i))DZ(k, i) = dzmax
                     enddo
                 enddo
 
!******          Wind, velocity, and bottom shear stresses @ top and bottom of
 
!                cell
                 sb(:, iut:idt - 1) = 0.0
                 do i = iut, idt - 1
                     ST(kt, i) = WSHX(i)*BR(kti(i), i)
                     do k = kt + 1, KBMIN(i)
                         ST(k, i) = WSHX(i)*DECAY(k - 1, i)*BR(k, i)
                         if(.NOT.IMPLICIT_VISC(jw))ST(k, i) = ST(k, i)         &
                          & + az(k - 1, i)*(BR(k - 1, i) + BR(k, i))           &
                          & *0.5D0*(u(k - 1, i) - u(k, i))                     &
                          & /((AVH2(k - 1, i) + AVH2(k - 1, i + 1))*0.5D0)
                     enddo
                     gc2 = 0.0
                     if(FRIC(i)/=0.0)gc2 = g/(FRIC(i)*FRIC(i))
 
                     hrad = BHR2(kt, i)                                        &
                          & /(BR(kti(i), i) - BR(kt + 1, i) + 2.0D0*AVHR(kt, i)&
                          & )
                     if(macrophyte_on .AND. MANNINGS_N(jw))then
                         call MACROPHYTE_FRICTION(hrad, FRIC(i), effric, kt, i)
                         gc2 = g*effric*effric/hrad**0.33333333D0
                     elseif(.NOT.macrophyte_on .AND. MANNINGS_N(jw))then
                         gc2 = g*FRIC(i)*FRIC(i)/hrad**0.33333333D0
                     endif
                     if(ONE_LAYER(i))then
                         sb(kt, i) = ST(kt + 1, i)                             &
                                   & + gc2*(BR(kti(i), i) + 2.0D0*AVHR(kt, i)) &
                                   & *u(kt, i)*DABS(u(kt, i))
                     else
                         sb(kt, i) = gc2*(BR(kti(i), i) - BR(kt + 1, i) + 2.0D0&
                                   & *AVHR(kt, i))*u(kt, i)*DABS(u(kt, i))
                         do k = kt + 1, KBMIN(i) - 1
                             hrad = (BHR2(k, i)                                &
                                  & /(BR(k, i) - BR(k + 1, i) + 2.0D0*H(k, jw))&
                                  & )
                             if(macrophyte_on .AND. MANNINGS_N(jw))then
                                 call MACROPHYTE_FRICTION(hrad, FRIC(i),       &
                                   & effric, k, i)
                                 gc2 = g*effric*effric/hrad**0.33333333D0
                             elseif(.NOT.macrophyte_on .AND. MANNINGS_N(jw))   &
                                  & then
                                 gc2 = g*FRIC(i)*FRIC(i)/hrad**0.33333333D0
                             endif
                             sb(k, i) = gc2*(BR(k, i) - BR(k + 1, i) + 2.0D0*H(&
                                      & k, jw))*u(k, i)*DABS(u(k, i))
                         enddo
                         if(kt/=KBMIN(i))then
                             hrad = (BHR2(KBMIN(i), i)                         &
                                  & /(BR(KBMIN(i), i) + 2.0D0*H(KBMIN(i), jw)))
                             if(macrophyte_on .AND. MANNINGS_N(jw))then
                                 call MACROPHYTE_FRICTION(hrad, FRIC(i),       &
                                   & effric, KBMIN(i), i)
                                 gc2 = g*effric*effric/hrad**0.33333333D0
                             elseif(.NOT.macrophyte_on .AND. MANNINGS_N(jw))   &
                                  & then
                                 gc2 = g*FRIC(i)*FRIC(i)/hrad**0.33333333D0
                             endif
 
                             if(KBMIN(i)/=KB(i))then
                                 sb(KBMIN(i), i)                               &
                                   & = gc2*(BR(KBMIN(i), i) - BR(KBMIN(i) + 1, &
                                   & i) + 2.0D0*h2(k, i))*u(KBMIN(i), i)       &
                                   & *DABS(u(KBMIN(i), i))
                             else
                                 sb(KBMIN(i), i)                               &
                                   & = gc2*(BR(KBMIN(i), i) + 2.0D0*h2(k, i))  &
                                   & *u(KBMIN(i), i)*DABS(u(KBMIN(i), i))
                             endif
                         endif
                     endif
                     do k = kt, KBMIN(i) - 1
                         sb(k, i) = sb(k, i) + ST(k + 1, i)
                     enddo
                     sb(KBMIN(i), i) = sb(KBMIN(i), i) + WSHX(i)               &
                                     & *DECAY(KBMIN(i), i)                     &
                                     & *(BR(KBMIN(i) - 1, i) + BR(KBMIN(i), i))&
                                     & *0.5D0
                 enddo
 
!******          Horizontal advection of momentum
 
                 do i = iu, id - 1
                     do k = kt, KBMIN(i)
                         udr = (1.0D0 + DSIGN(1.0D0, (u(k,i) + u(k,i+1))*0.5D0)&
                             & )*0.5D0
                         udl = (1.0D0 + DSIGN(1.0D0, (u(k,i) + u(k,i-1))*0.5D0)&
                             & )*0.5D0
                         ADMX(k, i) = (BH2(k, i + 1)*(u(k, i + 1) + u(k, i))   &
                                    & *0.5D0*(udr*u(k, i) + (1.0 - udr)        &
                                    & *u(k, i + 1)) - BH2(k, i)                &
                                    & *(u(k, i) + u(k, i - 1))                 &
                                    & *0.5D0*(udl*u(k, i - 1) + (1.0D0 - udl)  &
                                    & *u(k, i)))/DLXR(i)
                     enddo
                 enddo
 
!******          Horizontal dispersion of momentum
 
                 do i = iu, id - 1
                     do k = kt, KBMIN(i)
                         if(AX(jw)>=0.0)then
                             DM(k, i) = AX(jw)                                 &
                                      & *(BH2(k, i + 1)*(u(k, i + 1) - u(k, i))&
                                      & /DLX(i + 1) - BH2(k, i)                &
                                      & *(u(k, i) - u(k, i - 1))/DLX(i))       &
                                      & /DLXR(i)
                         else
                             DM(k, i) = ABS(u(k, i))*ABS(AX(jw))*H(k, jw)      &
                                      & *(BH2(k, i + 1)*(u(k, i + 1) - u(k, i))&
                                      & /DLX(i + 1) - BH2(k, i)                &
                                      & *(u(k, i) - u(k, i - 1))/DLX(i))       &
                                      & /DLXR(i)                                                                                            ! SW 8/2/2017 SCALE AX WITH U, FOR EXAMPLE AX=0.1U
                         endif
                     enddo
                 enddo
 
!******          Vertical advection of momentum
 
                 do i = iu, id - 1
                     do k = kt, KB(i) - 1
                         ab = (1.0D0 + DSIGN(1.0D0, (w(k,i+1) + w(k,i))*0.5D0))&
                            & *0.5D0
                         ADMZ(k, i) = (BR(k, i) + BR(k + 1, i))                &
                                    & *0.5D0*(w(k, i + 1) + w(k, i))           &
                                    & *0.5D0*(ab*u(k, i) + (1.0 - ab)          &
                                    & *u(k + 1, i))
                     enddo
                 enddo
 
!******          Gravity force due to channel slope
 
                 do i = iu - 1, id
                     GRAV(kt, i) = AVHR(kt, i)*(BKT(i) + BKT(i + 1))           &
                                 & *0.5D0*g*SINAC(jb)
                     do k = kt + 1, KB(i)
                         GRAV(k, i) = BHR2(k, i)*g*SINAC(jb)
                     enddo
                 enddo
 
                 if(ICEC(jw)=='    ONWB')then
                     do i = iu, id
                                ! water loss due to ice formation or water gain due to ice melting
                         QSS(kt, i) = QSS(kt, i) + ICEQSS(i)
                         ICEQSS(i) = 0.0D00
                     enddo
                 endif
 
 
!***********************************************************************************************************************************
!**              Task 2.2.3: Water surface elevation                          
!***********************************************************************************************************************************
 
!******          ** Tridiagonal coefficients
 
                 bhrho(iu - 1:id + 1) = 0.0D0
 
 
!***********************************************************************************************************************************
!**              Task 2.2.3: Water surface elevation                          
!***********************************************************************************************************************************
 
!******          ** Tridiagonal coefficients
 
                 d(iu - 1:id + 1) = 0.0D0
 
 
!***********************************************************************************************************************************
!**              Task 2.2.3: Water surface elevation                          
!***********************************************************************************************************************************
 
!******          ** Tridiagonal coefficients
 
                 f(iu - 1:id + 1) = 0.0D0
                 do i = iu, id - 1
                     do k = kt, KBMIN(i)
                         bhrho(i) = bhrho(i)                                   &
                                  & + (BH2(k, i + 1)/RHO(k, i + 1) + BH2(k, i) &
                                  & /RHO(k, i))
                     enddo
                     do k = kt, KB(i)
                         d(i) = d(i)                                           &
                              & + (u(k, i)*BHR2(k, i) - u(k, i - 1)*BHR2(k,    &
                              & i - 1) - QSS(k, i)                             &
                              & + (UXBR(k, i) - UXBR(k, i - 1))*dlt)
                         f(i) = f(i)                                           &
                              & + ( - sb(k, i) + ST(k, i) - ADMX(k, i) + DM(k, &
                              & i) - HDG(k, i) + GRAV(k, i))
                     enddo
                 enddo
 
!******          Boundary tridiagonal coefficients
 
                 d(iu) = 0.0D0
                 do k = kt, KB(iu)
                     d(iu) = d(iu) + (u(k, iu)*BHR2(k, iu) - QSS(k, iu))       &
                           & + UXBR(k, iu)*dlt
                 enddo
                 if(DN_FLOW(jb))then
                     do k = kt, KB(id)
                         d(id) = d(id) - u(k, id - 1)*BHR2(k, id - 1)          &
                               & - QSS(k, id) + (UXBR(k, id) - UXBR(k, id - 1))&
                               & *dlt + QOUT(k, jb)
                     enddo
                 endif
                 if(UP_HEAD(jb))then
                     do k = kt, KBMIN(iu - 1)
                         bhrho(iu - 1) = bhrho(iu - 1)                         &
                           & + (BH2(k, iu)/RHO(k, iu) + BH2(k, iu - 1)         &
                           & /RHO(k, iu - 1))
                     enddo
                     do k = kt, KB(iu)
                         d(iu) = d(iu) - u(k, iu - 1)*BHR2(k, iu - 1)
                         f(iu - 1) = f(iu - 1)                                 &
                                   & - (sb(k, iu - 1) - ST(k, iu - 1) +        &
                                   & HDG(k, iu - 1) - GRAV(k, iu - 1))
                     enddo
                 endif
                 if(DN_HEAD(jb))then
                     do k = kt, KBMIN(id)
                         bhrho(id) = bhrho(id)                                 &
                                   & + (BH2(k, id + 1)/RHO(k, id + 1) +        &
                                   & BH2(k, id)/RHO(k, id))
                     enddo
                     do k = kt, KB(id)
                         d(id) = d(id)                                         &
                               & + (u(k, id)*BHR2(k, id) - u(k, id - 1)*BHR2(k,&
                               & id - 1) - QSS(k, id))                         &
                               & + (UXBR(k, id) - UXBR(k, id - 1))*dlt
                         f(id) = f(id) + ( - sb(k, id) + ST(k, id) - HDG(k, id)&
                               & + GRAV(k, id))
                     enddo
                 endif
             enddo
         enddo
         do jw = 1, nwb
             kt = ktwb(jw)
             do jb = BS(jw), BE(jw)
                 if(BR_INACTIVE(jb))cycle
                 iu = CUS(jb)
                 id = DS(jb)
                 if(INTERNAL_FLOW(jb) .AND. .NOT.DAM_INFLOW(jb))then                                                   !TC 08/03/04
                     QIN(jb) = 0.0D0
                     do k = ktwb(JWUH(jb)), KB(UHS(jb))
                         QIN(jb) = QIN(jb) + u(k, UHS(jb))*BHR2(k, UHS(jb))
                     enddo
                 endif
                 if(UP_FLOW(jb))d(iu) = d(iu) - QIN(jb)
 
!******          Boundary surface elevations
 
                 if(UH_INTERNAL(jb))then
                     z(iu - 1) = (( - EL(ktwb(JWUH(jb)), UHS(jb)) + z(UHS(jb))*&
                               & COSA(JBUH(jb))) + EL(kt, iu - 1) + SINA(jb)   &
                               & *DLXR(iu - 1))/COSA(jb)
                     elws(iu - 1) = EL(kt, iu - 1) - z(iu - 1)*COSA(jb)
                     kti(iu - 1) = 2
                     do while (EL(kti(iu - 1), iu - 1)>elws(iu - 1))
                         kti(iu - 1) = kti(iu - 1) + 1
                     enddo
                     kti(iu - 1) = MAX(kti(iu - 1) - 1, 2)
                 endif
                 if(UH_EXTERNAL(jb))z(iu - 1)                                  &
                  & = (EL(kt, iu - 1) - (ELUH(jb) + SINA(jb)*DLX(iu)*0.5D0))   &
                  & /COSA(jb)
                 if(DH_INTERNAL(jb))then
                     z(id + 1) = (( - EL(ktwb(JWDH(jb)), DHS(jb)) + z(DHS(jb))*&
                               & COSA(JBDH(jb))) + EL(kt, id + 1))/COSA(jb)
                     elws(id + 1) = EL(kt, id + 1) - z(id + 1)*COSA(jb)
                     kti(id + 1) = 2
                     do while (EL(kti(id + 1), id + 1)>elws(id + 1))
                         kti(id + 1) = kti(id + 1) + 1
                     enddo
                     kti(id + 1) = MAX(kti(id + 1) - 1, 2)
                     if(kti(id + 1)>=KB(id))then
                         z(id + 1) = z(id) - SLOPE(jb)*DLX(id)/2.0D0
                         elws(id + 1) = EL(kt, id + 1) - z(id + 1)*COSA(jb)
                         kti(id + 1) = 2
                         do while (EL(kti(id + 1), id + 1)>elws(id + 1))
                             kti(id + 1) = kti(id + 1) + 1
                         enddo
                         kti(id + 1) = MAX(kti(id + 1) - 1, 2)
                     endif
                 endif
                 if(DH_EXTERNAL(jb))z(id + 1)                                  &
                  & = (EL(kt, id + 1) - (ELDH(jb) - SINA(jb)*DLX(id)*0.5D0))   &
                  & /COSA(jb)
 
!******          Implicit water surface elevation solution
 
                 do i = iu, id
                     !CONCURRENT(I=IU:ID)   !FORALL(I=IU:ID)                       !DO I=IU,ID
                     A(i) = -RHO(kt, i - 1)*g*COSA(jb)*dlt*dlt*bhrho(i - 1)    &
                          & *0.5D0/DLXR(i - 1)
                     C(i) = -RHO(kt, i + 1)*g*COSA(jb)*dlt*dlt*bhrho(i)        &
                          & *0.5D0/DLXR(i)
                     V(i) = RHO(kt, i)*g*COSA(jb)                              &
                          & *dlt*dlt*(bhrho(i)*0.5D0/DLXR(i) + bhrho(i - 1)    &
                          & *0.5D0/DLXR(i - 1)) + DLX(i)*BI(kt, i)
                     d(i) = dlt*(d(i) + dlt*(f(i) - f(i - 1))) + DLX(i)        &
                          & *BI(kt, i)*z(i)
                 enddo
                 if(UP_HEAD(jb))d(iu) = d(iu) - A(iu)*z(iu - 1)
                 if(DN_HEAD(jb))d(id) = d(id) - C(id)*z(id + 1)
                 BTA(iu) = V(iu)
                 GMA(iu) = d(iu)
                 do i = iu + 1, id
                     BTA(i) = V(i) - A(i)/BTA(i - 1)*C(i - 1)
                     GMA(i) = d(i) - A(i)/BTA(i - 1)*GMA(i - 1)
                 enddo
                 z(id) = GMA(id)/BTA(id)
                 do i = id - 1, iu, -1
                     z(i) = (GMA(i) - C(i)*z(i + 1))/BTA(i)
                 enddo
 
!******          Boundary water surface elevations
 
                 if(UP_FLOW(jb) .AND. .NOT.HEAD_FLOW(jb))z(iu - 1) = z(iu)
                 if(UP_FLOW(jb) .AND. HEAD_FLOW(jb))z(iu - 1)                  &
                  & = ( - EL(ktwb(JWUH(jb)), UHS(jb)) + z(UHS(jb))             &
                  & *COSA(JBUH(jb)) + EL(kt, iu - 1) + SINA(JBUH(jb))          &
                  & *DLXR(iu - 1))/COSA(JBUH(jb))
                 if(DN_FLOW(jb))z(id + 1) = z(id)
 
!******          Updated surface layer and geometry
 
                 if(.NOT.TRAPEZOIDAL(jw))then                                                                          !SW 07/16/04
                     do i = iu - 1, id + 1
                         if(EL(kt, i) - z(i)*COSA(jb)>EL(kti(i), i))then
                             do while (EL(kt, i) - z(i)*COSA(jb)>EL(kti(i), i) &
                                     & .AND. kti(i)/=2)
                                 z(i) = (EL(kt, i) - EL(kti(i), i) - (EL(kt, i)&
                                      & - EL(kti(i), i) - z(i)*COSA(jb))       &
                                      & *(b(kti(i), i)/b(kti(i) - 1, i)))      &
                                      & /COSA(jb)
 
                                 if(macrophyte_on)then
                                     ktip = kti(i)
!C                                   KEEPING TRACK IF COLUMN KTI HAS
!                                    MACROPHYTES
                                     if(ktip>2)KTICOL(i) = .FALSE.
                                 endif
 
                                 kti(i) = MAX(kti(i) - 1, 2)
                             enddo
                         elseif(EL(kt, i) - z(i)*COSA(jb)<EL(kti(i) + 1, i))   &
                              & then
                             do while (EL(kt, i) - z(i)*COSA(jb)               &
                                     & <EL(kti(i) + 1, i) .AND. kti(i)<KB(i))                           ! sw 7/18/11
                                 z(i) = (EL(kt, i) - EL(kti(i) + 1, i) - (EL(kt&
                                      & , i) - EL(kti(i) + 1, i) - z(i)        &
                                      & *COSA(jb))                             &
                                      & *(b(kti(i), i)/b(kti(i) + 1, i)))      &
                                      & /COSA(jb)
                                 kti(i) = kti(i) + 1
                                 if(macrophyte_on)KTICOL(i) = .TRUE.
                                 if(kti(i)>=KB(i))exit
                             enddo
                         endif
                         BI(kt:KB(i), i) = b(kt:KB(i), i)
                         BI(kt, i) = b(kti(i), i)
                         H1(kt, i) = H(kt, jw) - z(i)
                         AVH1(kt, i) = (H1(kt, i) + H1(kt + 1, i))*0.5D0
                         if(kt==kti(i) .OR. kti(i)>=KB(i))then
                             BH1(kt, i) = b(kt, i)*H1(kt, i)
                         else
                             BH1(kt, i) = BI(kt, i)                            &
                               & *(EL(kt, i) - z(i)*COSA(jb) -                 &
                               & EL(kti(i) + 1, i))/COSA(jb)
                         endif
                         do k = kti(i) + 1, kt
                             BH1(kt, i) = BH1(kt, i) + BNEW(k, i)*H(k, jw)
                                                      !BNEW(K,I)*H(K,JW)   ! SW 1/23/06
                         enddo
                         BKT(i) = BH1(kt, i)/H1(kt, i)
                         if(KBI(i)<KB(i))BKT(i) = BH1(kt, i)                   &
                          & /(H1(kt, i) - (EL(KBI(i) + 1, i) - EL(KB(i) + 1, i)&
                          & ))                                                              ! SW 1/23/06
                         VOL(kt, i) = BH1(kt, i)*DLX(i)
                     enddo
                     do i = iu - 1, id
                         AVHR(kt, i) = H1(kt, i) + (H1(kt, i + 1) - H1(kt, i)) &
                                     & /(0.5D0*(DLX(i) + DLX(i + 1)))          &
                                     & *0.5D0*DLX(i)                                                                       !SW 07/29/04  (H1(KT,I+1) +H1(KT,I))*0.5
                         if(KBI(i)<KB(i))AVHR(kt, i)                           &
                          & = (H1(kt, i) - (EL(KBI(i) + 1, i)                  &
                          & - EL(KB(i) + 1, i)))                               &
                          & + (H1(kt, i + 1) - (EL(KBI(i) + 1, i + 1)          &
                          & - EL(KB(i) + 1, i + 1)) - H1(kt, i)                &
                          & + (EL(KBI(i) + 1, i) - EL(KB(i) + 1, i)))          &
                          & /(0.5D0*(DLX(i) + DLX(i + 1)))*0.5D0*DLX(i)       ! SW 1/23/06
                         BHR1(kt, i) = BH1(kt, i)                              &
                                     & + (BH1(kt, i + 1) - BH1(kt, i))         &
                                     & /(0.5D0*(DLX(i) + DLX(i + 1)))          &
                                     & *0.5D0*DLX(i)                                                                        !SW 07/29/04 (BH1(KT,I+1)+BH1(KT,I))*0.5
                     enddo
                     AVHR(kt, id + 1) = H1(kt, id + 1)
                     BHR1(kt, id + 1) = BH1(kt, id + 1)
                     DLVOL(jb) = 0.0
                 else                                                                                                  !SW 07/16/04
                     do i = iu - 1, id + 1
                         BI(kt:KB(i), i) = b(kt:KB(i), i)
                         call GRID_AREA2
                         H1(kt, i) = H(kt, jw) - z(i)
                         AVH1(kt, i) = (H1(kt, i) + H1(kt + 1, i))*0.5
                         call GRID_AREA1(EL(kt, i) - z(i), EL(kt + 1, i),      &
                           & BH1(kt, i), BI(kt, i))
                         BKT(i) = BH1(kt, i)/H1(kt, i)
                         if(KBI(i)<KB(i))BKT(i) = BH1(kt, i)                   &
                          & /(H1(kt, i) - (EL(KBI(i) + 1, i) - EL(KB(i) + 1, i)&
                          & ))                                                              ! SW 1/23/06
                         VOL(kt, i) = BH1(kt, i)*DLX(i)
                     enddo
                     do i = iu - 1, id
                         AVHR(kt, i) = H1(kt, i) + (H1(kt, i + 1) - H1(kt, i)) &
                                     & /(0.5D0*(DLX(i) + DLX(i + 1)))          &
                                     & *0.5D0*DLX(i)                                                                       !SW 07/29/04
                         if(KBI(i)<KB(i))AVHR(kt, i)                           &
                          & = (H1(kt, i) - (EL(KBI(i) + 1, i)                  &
                          & - EL(KB(i) + 1, i)))                               &
                          & + (H1(kt, i + 1) - (EL(KBI(i) + 1, i + 1)          &
                          & - EL(KB(i) + 1, i + 1)) - H1(kt, i)                &
                          & + (EL(KBI(i) + 1, i) - EL(KB(i) + 1, i)))          &
                          & /(0.5D0*(DLX(i) + DLX(i + 1)))*0.5D0*DLX(i)                                                    ! SW 1/23/06
                         BHR1(kt, i) = BH1(kt, i)                              &
                                     & + (BH1(kt, i + 1) - BH1(kt, i))         &
                                     & /(0.5*(DLX(i) + DLX(i + 1)))*0.5*DLX(i)                                         !SW 07/29/04
                     enddo
                     AVHR(kt, id + 1) = H1(kt, id + 1)
                     BHR1(kt, id + 1) = BH1(kt, id + 1)
                     DLVOL(jb) = 0.0
                 endif
                 elws(CUS(jb):DS(jb) + 1) = EL(kt, CUS(jb):DS(jb) + 1)         &
                   & - z(CUS(jb):DS(jb) + 1)*COSA(jb)
                 do i = iu, id
                     DLVOL(jb) = DLVOL(jb) + (BH1(kt, i) - BH2(kt, i))*DLX(i)
                     if(kt==2 .AND. H1(kt, i)>H(2, jw) .AND.                   &
                      & .NOT.surface_warning)then
                         write(wrn, '(A,I0,A,F0.3)')                           &
                       &'Water surface is above the top of layer 2 in segment '&
                      & , i, ' at day ', jday
                         warning_open = .TRUE.
                         surface_warning = .TRUE.
                     endif
                 enddo
 
                 if(macrophyte_on)then
!C                   IF DEPTH IN KTI LAYER BECOMES GREATER THAN THRESHOLD,
!C                   SETTING MACROPHYTE CONC. IN KTI COLUMN TO INITIAL CONC.
                     do i = iu, id
                         depkti = elws(i) - EL(kti(i) + 1, i)
 
!*******                 MACROPHYTES, SETTING CONC. OF MACROPHYTES IN NEW
!*********               COLUMNS TO INITIAL CONCENTRATION IF COLUMN DEPTH IS
!                        GREATER THAN 'THRKTI'
                         if(.NOT.KTICOL(i) .AND. depkti>=thrkti)then
                             KTICOL(i) = .TRUE.
                             jt = kti(i)
                             MACT(jt, kt, i) = 0.0
                             do m = 1, nmc
                !MACRC(JT,KT,I,M)=MACWBCI(JW,M)
                                 if(ISO_MACROPHYTE(jw, m))macrc(jt, kt, i, m)  &
                                  & = MACWBCI(jw, m)                            ! cb 3/7/16
                                 if(VERT_MACROPHYTE(jw, m))macrc(jt, kt, i, m) &
                                  & = 0.1
                                 if(LONG_MACROPHYTE(jw, m))macrc(jt, kt, i, m) &
                                  & = 0.1
                                 colb = EL(kti(i) + 1, i)
                                 coldep = elws(i) - colb
                !MACRM(JT,KT,I,M)=MACWBCI(JW,M)*COLDEP*CW(JT,I)*DLX(I)
                                 macrm(jt, kt, i, m) = macrc(jt, kt, i, m)     &
                                   & *coldep*CW(jt, i)*DLX(i)                    ! cb 3/17/16
                                 MACT(jt, kt, i) = MACT(jt, kt, i)             &
                                   & + MACWBCI(jw, m)
                                 macmbrt(jb, m) = macmbrt(jb, m)               &
                                   & + macrm(jt, kt, i, m)
                             enddo
                         endif
 
!******                  MACROPHYTES, WHEN COLUMN DEPTH IS LESS THAN 'THRKTI',
!                        ZEROING OUT CONC.
                         if(KTICOL(i) .AND. depkti<thrkti)then
                             KTICOL(i) = .FALSE.
                             jt = kti(i)
                             MACT(jt, kt, i) = 0.0
                             do m = 1, nmc
                                 macmbrt(jb, m) = macmbrt(jb, m)               &
                                   & - macrm(jt, kt, i, m)
                                 macrc(jt, kt, i, m) = 0.0
                                 macrm(jt, kt, i, m) = 0.0
                             enddo
                         endif
                     enddo
                 endif
 
!***********************************************************************************************************************************
!**              Task 2.2.4: Longitudinal velocities                          
!***********************************************************************************************************************************
 
!                **
                 iut = iu
                 idt = id
                 if(UP_HEAD(jb))iut = iu - 1
                 if(DN_HEAD(jb))idt = id + 1
 
!******          Pressures
 
                 do i = iut, idt
                     do k = kt, KB(i)
                         P(k, i) = P(k - 1, i) + RHO(k, i)*g*H1(k, i)*COSA(jb)
                     enddo
                 enddo
 
!******          Horizontal pressure gradients
 
                 do i = iut, idt - 1
                     HPG(kt, i) = DLXRHO(i)*(BKT(i) + BKT(i + 1))              &
                                & *0.5D0*(H1(kt, i + 1)*P(kt, i + 1)           &
                                & - H1(kt, i)*P(kt, i))
                     do k = kt + 1, KBMIN(i)
                         HPG(k, i) = DLXRHO(i)*BHR2(k, i)                      &
                                   & *((P(k - 1, i + 1) - P(k - 1, i))         &
                                   & + (P(k, i + 1) - P(k, i)))
                     enddo
                 enddo
 
!******          Boundary horizontal velocities
 
                 if(UP_FLOW(jb))then
                     if(.NOT.HEAD_FLOW(jb))then
                         qinf(:, jb) = 0.0
                         if(PLACE_QIN(jw))then
 
!************                Inflow layer
 
                             k = kt
                             sstot = 0.0
                             do jc = nsss, nsse
                                 sstot = sstot + CIN(jc, jb)
                             enddo
                             rhoin = DENSITY(TIN(jb), DMAX1(CIN(1, jb), 0.0D0),&
                                   & DMAX1(sstot, 0.0D0))
                             do while (rhoin>RHO(k, iu) .AND. k<KB(iu))
                                 k = k + 1
                             enddo
                             KTQIN(jb) = k
                             KBQIN(jb) = k
 
!************                Layer inflows
 
                             vqin = QIN(jb)*dlt
                             vqini = vqin
                             qinfr = 1.0
                             incr = -1
                             do while (qinfr>0.0D0)
                                 v1 = VOL(k, iu)
                                 if(k<=KB(iu))then
                                     if(vqin>0.5D0*v1)then
                                         qinf(k, jb) = 0.5D0*v1/vqini
                                         qinfr = qinfr - qinf(k, jb)
                                         vqin = vqin - qinf(k, jb)*vqini
                                         if(k==kt)then
                                         k = KBQIN(jb)
                                         incr = 1
                                         endif
                                     else
                                         qinf(k, jb) = qinfr
                                         qinfr = 0.0D0
                                     endif
                                     if(incr<0)KTQIN(jb) = k
                                     if(incr>0)KBQIN(jb) = MIN(KB(iu), k)
                                     k = k + incr
                                 else
                                     qinf(kt, jb) = qinf(kt, jb) + qinfr
                                     qinfr = 0.0D0
                                 endif
                             enddo
                         else
                             KTQIN(jb) = kt
                             KBQIN(jb) = KB(iu)
                             bhsum = 0.0D0
                             do k = kt, KB(iu)
                                 bhsum = bhsum + BH1(k, iu)
                             enddo
                             do k = kt, KB(iu)
                                 qinf(k, jb) = BH1(k, iu)/bhsum
                             enddo
                         endif
                         do k = kt, KB(iu)
                             u(k, iu - 1) = qinf(k, jb)*QIN(jb)/BHR1(k, iu - 1)
                         enddo
                     else
                         KTQIN(jb) = kt
                         KBQIN(jb) = KB(iu)
                         if(JBUH(jb)<=BE(jw) .AND. JBUH(jb)>=BS(jw))then
                             do k = kt, KB(iu)
                                 u(k, iu - 1) = u(k, UHS(jb))*BHR1(k, UHS(jb)) &
                                   & /BHR1(k, iu - 1)
                             enddo
                         else
                             call UPSTREAM_VELOCITY
                         endif
                     endif
                 endif
                 if(DN_FLOW(jb))then
                     do k = kt, KB(id)
                         u(k, id) = QOUT(k, jb)/BHR1(k, id)
                     enddo
                 endif
                 if(UP_HEAD(jb))then
                     do k = kt, KB(iu - 1)
                         u(k, iu - 1) = (BHR2(k, iu - 1)*u(k, iu - 1) + dlt*( -&
                                      & sb(k, iu - 1) + ST(k, iu - 1)          &
                                      & - HPG(k, iu - 1) + GRAV(k, iu - 1)))   &
                                      & /BHR1(k, iu - 1)
                     enddo
                 endif
                 if(DN_HEAD(jb))then
                     do k = kt, KB(id + 1)
                         u(k, id) = (BHR2(k, id)*u(k, id) + dlt*( - sb(k, id) +&
                                  & ST(k, id) - HPG(k, id) + GRAV(k, id)))     &
                                  & /BHR1(k, id)
                     enddo
                 endif
 
!******          Horizontal velocities
 
                 do i = iu, id - 1
                     do k = kt, KBMIN(i)
                         u(k, i) = (BHR2(k, i)*u(k, i))/BHR1(k, i)             &
                                 & + (dlt*( - sb(k, i) + ST(k, i) - ADMZ(k, i) &
                                 & + ADMZ(k - 1, i) - ADMX(k, i) + DM(k, i)    &
                                 & - HPG(k, i) + GRAV(k, i) + UXBR(k, i)       &
                                 & /h2(k, i)))/BHR1(k, i)
                         if(INTERNAL_WEIR(k, i))u(k, i) = 0.0D0
                     enddo
                 enddo
 
!******          Implicit vertical eddy viscosity
 
                 if(IMPLICIT_VISC(jw))then
        !  AT = 0.0D0; CT = 0.0D0; VT = 0.0D0; DT = 0.0D0
                     do i = iut, idt - 1
                                      ! SW CODE SPEEDUP
                         do k = kt, KBMIN(i)
                             AT(k, i) = 0.0D0
                             CT(k, i) = 0.0D0
                             VT(k, i) = 0.0D0
                             DT(k, i) = 0.0D0
                         enddo
                     enddo
                     do i = iut, idt - 1
                         do k = kt, KBMIN(i)
                                        !KB(I)  SW 10/7/07
                             AT(k, i) = -dlt/BHR1(k, i)*az(k - 1, i)           &
                                      & *(BHR1(k - 1, i)/AVHR(k - 1, i)        &
                                      & + BR(k, i))                            &
                                      & /(AVH1(k - 1, i) + AVH1(k - 1, i + 1))
                             CT(k, i) = -dlt/BHR1(k, i)*az(k, i)               &
                                      & *(BHR1(k, i)/AVHR(k, i) + BR(k + 1, i))&
                                      & /(AVH1(k, i) + AVH1(k, i + 1))
                             VT(k, i) = 1.0D0 - AT(k, i) - CT(k, i)
                             DT(k, i) = u(k, i)
                         enddo
                         call TRIDIAG(AT(:, i), VT(:, i), CT(:, i), DT(:, i),  &
                                    & kt, KBMIN(i), kmx, u(:, i))
                     enddo
                 endif
 
!******          Corrected horizontal velocities
 
                 if(UP_HEAD(jb))then
                     is = id
                     ie = iu - 1
                     incr = -1
                     Q(is) = 0.0D0
                     do k = kt, KB(id)
                         Q(is) = Q(is) + u(k, is)*BHR1(k, is)
                     enddo
                     QSSUM(is) = 0.0D0
                     do k = kt, KB(is)
                         QSSUM(is) = QSSUM(is) + QSS(k, is)
                     enddo
                 else
                     is = iu - 1
                     ie = id
                     incr = 1
                     if(DN_FLOW(jb))ie = id - 1
                     Q(is) = 0.0D0
                     do k = kt, KB(iu)
                         Q(is) = Q(is) + u(k, is)*BHR1(k, is)
                     enddo
                 endif
                 QC(is) = Q(is)
                 do i = is + incr, ie, incr
                     QSSUM(i) = 0.0D0
                     do k = kt, KB(i)
                         QSSUM(i) = QSSUM(i) + QSS(k, i)
                     enddo
                     bhrsum = 0.0D0
                     Q(i) = 0.0D0
                     do k = kt, KBMIN(i)
                         if(.NOT.INTERNAL_WEIR(k, i))then
                             bhrsum = bhrsum + BHR1(k, i)
                             Q(i) = Q(i) + u(k, i)*BHR1(k, i)
                         endif
                     enddo
                     if(UP_HEAD(jb))then
                         QC(i) = QC(i + 1) + (BH1(kt, i + 1) - BH2(kt, i + 1)) &
                               & *DLX(i + 1)/dlt - QSSUM(i + 1)
                     else
                         QC(i) = QC(i - 1) - (BH1(kt, i) - BH2(kt, i))*DLX(i)  &
                               & /dlt + QSSUM(i)
                     endif
                     do k = kt, KBMIN(i)
                         if(INTERNAL_WEIR(k, i))then
                             u(k, i) = 0.0D0
                         else
                             u(k, i) = u(k, i) + (QC(i) - Q(i))/bhrsum
                             if(Q(i)/=0.0)QERR(i) = (Q(i) - QC(i))/Q(i)*100.0
                         endif
                     enddo
                 enddo
 
!******          Head boundary flows
 
                 if(UP_HEAD(jb))quh1(kt:KB(iu - 1), jb)                        &
                  & = u(kt:KB(iu - 1), iu - 1)*BHR1(kt:KB(iu - 1), iu - 1)
                 if(DN_HEAD(jb))qdh1(kt:KB(id + 1), jb) = u(kt:KB(id + 1), id) &
                  & *BHR1(kt:KB(id + 1), id)
 
!***********************************************************************************************************************************
!**              Task 2.2.5: Vertical velocities                              
!***********************************************************************************************************************************
 
!                **
                 do i = iu, id
                     do k = KB(i) - 1, kt, -1
                         wt1 = w(k + 1, i)*BB(k + 1, i)
                         wt2 = (BHR(k + 1, i)*u(k + 1, i) - BHR(k + 1, i - 1)  &
                             & *u(k + 1, i - 1) - QSS(k + 1, i))/DLX(i)
                         w(k, i) = (wt1 + wt2)/BB(k, i)
                     enddo
                 enddo
             enddo
         enddo
 
!***********************************************************************************************************************************
!**      Task 2.2.6: Autostepping                                             
!***********************************************************************************************************************************
 
!        **
         do jw = 1, nwb
             kt = ktwb(jw)
             do jb = BS(jw), BE(jw)
                 do i = CUS(jb), DS(jb)
 
                     if(DLTADD(jw)=='      ON' .AND. ABS(H1(kt, i) - h2(kt, i))&
                      & /h2(kt, i)>0.5)then
                         write(wrn, '(A,F0.3,A,I0,A,F0.3/3(A,F0.3),a,i10,a)')  &
                      &'Computational warning |h1-h2|/h2>0.25 on Julian day = '&
                     & , jday, ' at segment ', i, ' timestep DLT= ', dlt,      &
                      &'   Water surface deviation [Z,m] = ', z(i),            &
                      &' H1 layer thickness(m) = ', H1(kt, i),                 &
                      &' H2 layer thickenss(m)=', h2(kt, i),                   &
                     & ' Iteration[NIT]=', nit, ' DLT reduced'
                         warning_open = .TRUE.
                         if(H1(kt, i)>0.0)then
                             curmax = dlt*0.5
                         else
                             curmax = dltmin
                         endif
                         goto 100
 
                     elseif(H1(kt, i)<0.0)then
                         write(wrn, '(A,F0.3,A,I0/4(A,F0.3))')                 &
                              &'Computational warning at Julian day = ', jday, &
                              &' at segment ', i, 'timestep = ', dlt,          &
                              &' water surface deviation [Z] = ', z(i),        &
                              &' m  layer thickness = ', H1(kt, i), ' m'
                         warning_open = .TRUE.
                         if(dlt>dltmin)then
                             write(wrn, '(A,I0/2(A,F0.3),A,I0)')               &
                                &'Negative surface layer thickness in segment '&
                               & , i, '  time step reduced to ', dltmin,       &
                                &' s on day ', jday, ' at iteration ', nit
                             warning_open = .TRUE.
                             curmax = dltmin
                             goto 100
                         else
                             write(w2err, '(A,F0.3/A,I0)')                     &
                                  &'Unstable water surface elevation on day ', &
                                 & jday, 'negative surface layer thickness ' //&
                                  &'using minimum timestep at iteration ', nit
                             write(w2err, *)'Branch #:', jb, ' in Waterbody:', &
                                 & jw
                             write(w2err, '(A)')                               &
        &'Segment, Surface layer thickness, m, Flow m3/s, U(KT,I) m/s, ELWS, m'
                             do ii = MAX(CUS(jb), i - 3), MIN(DS(jb), i + 3)
                                 write(w2err,                                  &
                               &'(T4,I3,T19,F10.2,t37,f10.2,1x,f10.2,2x,f10.2)'&
                              & )ii, H1(kt, ii), QC(ii), u(kt, ii), elws(ii)                                                                  ! SW 7/13/10
                             enddo
                             text = 'Runtime error - see w2.err'
                             error_open = .TRUE.
                             goto 200
                         endif
                     endif
                 enddo
                 do i = CUS(jb), DS(jb)
                     if(VISCOSITY_LIMIT(jw))then
                         if(AX(jw)>=0.0)tau1 = 2.0*AX(jw)/(DLX(i)*DLX(i))
                     endif
                     if(CELERITY_LIMIT(jw))                                    &
                      & celrty = SQRT((ABS(RHO(KB(i),i) - RHO(kt,i)))          &
                      & /1000.0*g*DEPTHB(KBI(i), i)*0.5)                                                                    ! SW 1/23/06
                     do k = kt, KB(i)
                         if(VISCOSITY_LIMIT(jw) .AND. .NOT.IMPLICIT_VISC(jw))  &
                          & tau2 = 2.0*az(k, i)/(H1(k, jw)*H1(k, jw))
                         QTOT(k, i) = (ABS(u(k, i))*BHR1(k, i) + ABS(u(k, i - 1&
                                    & ))*BHR1(k, i - 1)                        &
                                    & + (ABS(w(k,i))*BB(k, i) + ABS(w(k-1,i))  &
                                    & *BB(k - 1, i))*DLX(i) + DLX(i)           &
                                    & *ABS(BH2(k, i) - BH1(k, i))              &
                                    & /dlt + ABS(QSS(k, i)))*0.5
                         if(VISCOSITY_LIMIT(jw) .AND. AX(jw)<0.0)              &
                          & tau1 = 2.0*ABS(u(k, i))*ABS(AX(jw))*H(k, jw)       &
                          & /(DLX(i)*DLX(i))
                         dltcal = 1.0/((QTOT(k, i)/BH1(k, i) + celrty)/DLX(i)  &
                                & + tau1 + tau2 + nonzero)
                         if(dltcal<curmax)then
                             kloc = k
                             iloc = i
                             curmax = dltcal
                             if(dltff*curmax<mindlt)then
                                 kmin = k
                                 imin = i
                             endif
                         endif
                     enddo
                 enddo
             enddo
         enddo
 
!**      Restore timestep dependent variables and restart calculations
 
100      if(curmax<dlt .AND. dlt>dltmin)then
             dlt = dltff*curmax
             if(dlt<=dltmin)then
                 write(wrn, '(A,F0.3/A,F0.3,A)')                               &
                      &'Computational warning at Julian day = ', jday,         &
                      &' timestep = ', dlt, ' sec: DLT<DLTMIN set DLT=DLTMIN'
                 warning_open = .TRUE.
                 dlt = dltmin
             endif
             nv = nv + 1
             z = sz
             u = su
             w = sw
             az = saz
             AVH2 = savh2
             AVHR = savhr
             kti = skti
             BKT = sbkt
             QSS = 0.0
      !SP CEMA
      !if(sediment_diagenesis)then
      !  If(CEMARelatedCode .and. IncludeBedConsolidation)TSS       = 0.0  ! SW 7/27/2017
      !end if
      !End SP CEMA
             sb = 0.0
             dlts = dlt
 
             do jw = 1, nwb                                                                    ! SW 8/25/05
                 do jb = BS(jw), BE(jw)
                     do i = US(jb) - 1, DS(jb) + 1
                         VOL(ktwb(jw), i) = BH2(ktwb(jw), i)*DLX(i)
                         BI(ktwb(jw), i) = b(kti(i), i)
                     enddo
                 enddo
             enddo
 
 
             curmax = dltmaxx/dltff
             if(pipes)then
                 ys = yss
                 vs = vss
                 vst = vsts
                 yst = ysts
                 dtp = dtps
                 qold = qolds
             endif
 
!**********  Macrophytes
             do jw = 1, nwb
                 do m = 1, nmc
                     if(MACROPHYTE_CALC(jw, m))then
                         kt = ktwb(jw)
                         do jb = BS(jw), BE(jw)
                             do i = CUS(jb), DS(jb)
                                 do k = kt, KB(i)
                                     mac(k, i, m) = smac(k, i, m)
                                     if(KTICOL(i))then
                                         jt = kti(i)
                                     else
                                         jt = kti(i) + 1
                                     endif
                                     je = KB(i)
                                     do j = jt, je
                                         macrc(j, k, i, m) = smacrc(j, k, i, m)
                                         macrm(j, k, i, m) = SMACRM(j, k, i, m)
                                     enddo
                                 enddo
                             enddo
                         enddo
                     endif
                 enddo
             enddo
 
             goto 50
         endif
         dltlim(kmin, imin) = dltlim(kmin, imin) + 1.0
 
!**      Layer bottom and middle depths
 
         do jw = 1, nwb
             do jb = BS(jw), BE(jw)
                 do i = CUS(jb) - 1, DS(jb)
                     DEPTHB(ktwb(jw), i) = H1(ktwb(jw), i)
                     DEPTHM(ktwb(jw), i) = H1(ktwb(jw), i)*0.5
                     if(KBI(i)<KB(i) .AND.                                     &
                      & (EL(KBI(i) + 1, i) - EL(KB(i) + 1, i))<H1(ktwb(jw), i))&
                      & then                                                                  ! SW 7/22/10 if h1 < elev diff this means depth is below the bottom - if we ignore that the run will continue but if dpethb is negative it will bomb in computing DECAY
                         DEPTHB(ktwb(jw), i)                                   &
                           & = (H1(ktwb(jw), i) - (EL(KBI(i) + 1, i)           &
                           & - EL(KB(i) + 1, i)))                                  ! SW 1/23/06
                         DEPTHM(ktwb(jw), i)                                   &
                           & = (H1(ktwb(jw), i) - (EL(KBI(i) + 1, i)           &
                           & - EL(KB(i) + 1, i)))*0.5
                     endif
                     do k = ktwb(jw) + 1, kmx
                         DEPTHB(k, i) = DEPTHB(k - 1, i) + H1(k, i)
                         DEPTHM(k, i) = DEPTHM(k - 1, i)                       &
                                      & + (H1(k - 1, i) + H1(k, i))*0.5
                     enddo
                 enddo
             enddo
         enddo
 
         call TEMPERATURE
 
         if(constituents)call WQCONSTITUENTS
 
         if(fish_particle_exist)call FISH
                                 ! SW 4/30/15
 
!SP      CEMA
         if(sediment_diagenesis)then
             if(cemarelatedcode .AND. includebedconsolidation)                 &
              & call CEMAUPDATEVERTICALLAYERING
             if(cemarelatedcode .AND. includebedconsolidation)                 &
              & call CEMACOMPUTETURBIDITY
         endif
!End     SP CEMA
 
         call LAYERADDSUB
         if(error_open)exit
 
         call BALANCES
 
!SP      CEMA
!if(sediment_diagenesis)then
!        CEMATSSCopy = TSS
!end     if
!End     SP CEMA
 
         call UPDATE
 
         if(restart_in)then
             if(iopenfish==0)nxtmts = jday
         endif
 
         if(jday>=nxtmts .OR. jday>=TSRD(tsrdp + 1) .OR. nit==1)then
                                                                  ! OUTPUT AT FREQUENCY OF TSR FILES
             if(habtatc=='      ON')call FISHHABITAT(iopenfish)                                        ! OUTPUT AT FREQUENCY OF TSR FILES
             if(aeratec=='      ON' .AND. oxygen_demand)call AERATEOUTPUT
             if(restart_in)iopenfish = 1
             if(envirpc=='      ON')call ENVIRP
             iopenfish = 1
         endif                                                    ! OUTPUT AT FREQUENCY OF TSR FILES
 
 
         call OUTPUTA
!****    Screen output
         do jw = 1, nwb
             if(SCREEN_OUTPUT(jw))then
                 if(jday>=nxtmsc(jw) .OR. jday>=SCRD(scrdp(jw) + 1, jw))then
                     if(jday>=SCRD(scrdp(jw) + 1, jw))then
                         scrdp(jw) = scrdp(jw) + 1
                         nxtmsc(jw) = SCRD(scrdp(jw), jw)
                     endif
                     kt = ktwb(jw)
                     nxtmsc(jw) = nxtmsc(jw) + SCRF(scrdp(jw), jw)
                     call SCREEN_UPDATE(Dlg)
                     call DATE_AND_TIME(cdate, cctime)
 !         DO JH=1,NHY
 !           IF (HYDRO_PLOT(JH))       CALL GRAPH_UPDATE (JH,HYD(:,:,JH),     HNAME(JH), HYMIN(JH),1.0,       LNAME(JH))
 !         END DO
 !         DO JC=1,NCT
 !           IF (CONSTITUENT_PLOT(JC)) CALL GRAPH_UPDATE (JH+JC,C2(:,:,JC),   CNAME(JC), CMIN(JC), CMULT(JC), LNAME(JH+JC))
 !         END DO
 !         DO JD=1,NDC
 !           IF (DERIVED_PLOT(JD))     CALL GRAPH_UPDATE (JH+JC+JD,CD(:,:,JD),CDNAME(JD),CDMIN(JD),CDMULT(JD),LNAME(JH+JC+JD))
 !         END DO
                 endif
             endif
         enddo
     enddo
          ! END OF MAIN DO WHILE LOOP
 
200   if(stop_pushed)then
         text = 'Execution stopped at ' // cctime(1:2) // ':' // cctime(3:4)   &
               & // ':' // cctime(5:6) // ' on ' // cdate(5:6) // '/' //       &
              & cdate(7:8) // '/' // cdate(3:4)
         call RESTART_OUTPUT('rso.opt')
     endif
 
     if(end_run .AND. restart_out)call RESTART_OUTPUT('rso.opt')
                                                                ! cb 4/9/15 writing restart output at end of simulation if RSOC='ON'
     iopenfish = 3
     if(envirpc=="      ON")call ENVIRP
     if(habtatc=="      ON")call FISHHABITAT(iopenfish)                                       ! FINAL OUTPUT FOR ENVIR PERFORMANCE
     call ENDSIMULATION
!    FISH OUTPUT SW 4/30/15  *************
     if(fish_particle_exist)call FISHOUTPUT
 
!    CALL DEALLOCATE_GRAPH
 
300   if(closec=='      ON' .AND. end_run)then
         call EXITDIALOG(Dlg, text)
     else
         call STOP_W2(Dlg, text)
     endif
     end function CE_QUAL_W2
