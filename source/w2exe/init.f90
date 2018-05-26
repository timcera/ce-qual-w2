!*==init.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     subroutine INIT
 
     use MAIN
     use GLOBAL
     use NAMESC
     use GEOMC
     use LOGICC
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
     use BIOENERGETICS
     use CEMAVARS, ONLY:cemarelatedcode, nxtsediag
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: nb
     external RESTART_OUTPUT
!
!*** End of declarations rewritten by SPAG
!
 
!***********************************************************************************************************************************
!**  Task 1.1: Variable Initialization                                        
!***********************************************************************************************************************************
 
!***********************************************************************************************************************************
!**  ** Task 1.1.1: Zero Variables                                            
!***********************************************************************************************************************************
!    **
     iceqss = 0.0D00
 
!***********************************************************************************************************************************
!**  Task 1.1: Variable Initialization                                        
!***********************************************************************************************************************************
 
!***********************************************************************************************************************************
!**  ** Task 1.1.1: Zero Variables                                            
!***********************************************************************************************************************************
!    **
     water_age_active = .FALSE.               ! SR 7/27/2017
     kb = 0
     kbr = 0
     nac = 0
     ntac = 0
     nacd = 0
     nacin = 0
     nactr = 0
     nacdt = 0
     nacpr = 0
     ndsp = 0
     hmax = 0
     kbmax = 0
     dlxmax = 0
     kbqin = 0
     ktqin = 0
     qgt = 0.0
     qsp = 0.0                                                                                                     ! SW 8/26/15 Initialize Qgt and Qsp for screen output on restart
     naf = 0
     tiss = 0.0
     cshe = 0.0
     cin = 0.0
     tin = 0.0
     ev = 0.0
     dz = 0.0D0
     et = 0.0
     cshe = 0.0
     a = 0.0D0
     f = 0.0D0
     d = 0.0D0
     c = 0.0D0
     eltmf = 0.0
     el = 0.0
     dx = 0.0D0
     st = 0.0D0
     sb = 0.0D0
     dzq = 0.0D0
     tss = 0.0
     hseg = 0.0
     qss = 0.0D0
     hpg = 0.0D0
     hdg = 0.0D0
     vsh = 0.0D0
     qdh1 = 0.0D0
     admx = 0.0D0
     decay = 0.0D0
     admz = 0.0D0
     uybr = 0.0D0
     grav = 0.0D0
     fetch = 0.0D0
     fetchu = 0.0D0
     fetchd = 0.0D0
     dlttvd = 0.0
     icethu = 0.0
     iceth1 = 0.0
     iceth2 = 0.0
     p = 0.0D0
     celrty = 0.0
     tau1 = 0.0D0
     tau2 = 0.0D0
     volsr = 0.0
     voltr = 0.0
     af = 0.0
     ef = 0.0
     eltms = 0.0
     dm = 0.0D0
     qin = 0.0D0
     reaer = 0.0
     st = 0.0D0
     sb = 0.0D0
     admx = 0.0D0
     admz = 0.0D0
     hpg = 0.0D0
     hdg = 0.0D0
     rho = 0.0D0
     jdayts = 0.0
     jday1 = 0.0
     depthb = 0.0D0
     depthm = 0.0D0
     uxbr = 0.0D0
     bhrho = 0.0D0
     dlmr = 0.0D0
     sron = 0.0
     cssb = 0.0
     q = 0.0D0
     bh1 = 0.0D0
     bh2 = 0.0D0
     bhr1 = 0.0D0
     bhr2 = 0.0D0
     avhr = 0.0D0
     grav = 0.0D0
     kbp = 0
     dzt = 0.0D0
     azt = 0.0D0
     kfjw = 0.0
     qc = 0.0D0
     yss = 0.0
     ysts = 0.0
     qwd = 0.0D0
     qdtr = 0.0D0
     ttr = 0.0
     ctr = 0.0
     tdtr = 0.0
     qolds = 0.0D0
     dtps = 0.0
     vsts = 0.0
     vss = 0.0
     egt2 = 0.0
     hab = 100.0
     sedpinflux = 0.0
     sedninflux = 0.0
     rs = 0.0
     rn = 0.0
     rb = 0.0
     re = 0.0
     rc = 0.0
     ranlw = 0.0                                            ! SW 10/20/15
     br_inactive = .FALSE.
     sdfirstadd = .TRUE.
                      ! cb 9/3/17
     if(.NOT.restart_in)then
         nsprf = 0
         izmin = 0
         ktwb = 2
         kmin = 1
         imin = 1
         t1 = 0.0D0
         t2 = 0.0D0
         c1 = 0.0D0
         c2 = 0.0D0
         cd = 0.0
         cin = 0.0
         c1s = 0.0
         kf = 0.0
         cmbrt = 0.0
         kfs = 0.0
         u = 0.0D0
         w = 0.0D0
         su = 0.0D0
         sw = 0.0D0
         saz = 0.0D0
         az = 0.0D0
         esbr = 0.0
         epd = 0.0
         etbr = 0.0
         ebri = 0.0
         dltlim = 0.0
         volev = 0.0
         volpr = 0.0
         voldt = 0.0
         volwd = 0.0
         current = 0.0
         volice = 0.0
         icebank = 0.0
         voluh = 0.0
         voldh = 0.0
         volin = 0.0
         volout = 0.0
         volsbr = 0.0
         voltrb = 0.0
         tsss = 0.0
         tssb = 0.0
         ef = 0.0
         tssev = 0.0
         tsspr = 0.0
         tsstr = 0.0
         tssdt = 0.0
         tsswd = 0.0
         tssuh = 0.0
         tssdh = 0.0
         tssin = 0.0
         cssk = 0.0
         tssout = 0.0
         tssice = 0.0
         tssuh1 = 0.0
         tssuh2 = 0.0
         cssuh1 = 0.0
         cssuh2 = 0.0
         tssdh1 = 0.0
         tssdh2 = 0.0
         cssdh1 = 0.0
         cssdh2 = 0.0
         qind = 0.0
         tind = 0.0
         cind = 0.0
         savh2 = 0.0
         savhr = 0.0
         voluh2 = 0.0
         avh1 = 0.0
         avh2 = 0.0
         voldh2 = 0.0
         z = 0.0D0
         quh1 = 0.0D0
         sed = 0.0
         sedc = 0.0
         sedn = 0.0
         vs = 0.0
         ys = 0.0
         yst = 0.0
         vst = 0.0
         dtp = 0.0
         qold = 0.0
         qsum = 0.0
         voltbr = 0.0
         dlvol = 0.0
         evbr = 0.0                                                                                                               ! SW 7/24/2017
         tpout = 0.0
         tptrib = 0.0
         tpdtrib = 0.0
         tpwd = 0.0
         tppr = 0.0
         tpin = 0.0
         tnout = 0.0
         tntrib = 0.0
         tndtrib = 0.0
         tnwd = 0.0
         tnpr = 0.0
         tnin = 0.0
         tn_sedsod_nh4 = 0.0
         tp_sedsod_po4 = 0.0                                                                                                                                      ! SW 2/19/16  TP_SEDBURIAL=0.0;TN_SEDBURIAL=0.0;
         sedp = 0.0
         iceth = 0.0
         pfluxin = 0.0
         nfluxin = 0.0                                                                  ! SW 4/19/10
         zmin = -1000.0
         tke = 0.0                                              ! SG 10/4/07
         sedp = 0.0
         sedc = 0.0
         sedn = 0.0
         macmbrt = 0.0
         macrc = 0.0
         smacrc = 0.0
         mac = 0.0
         smac = 0.0
         macrm = 0.0
         epm = 0.0
         macss = 0.0                                                               ! cb 3/8/16
         kticol = .FALSE.
     endif
     anlim = 1.0
     aplim = 1
     aslim = 1.0
     allim = 1.0
     enlim = 1.0
     eplim = 1
     eslim = 1.0
     ellim = 1.0
     kloc = 1
     iloc = 1
     mnlim = 1.0
     mplim = 1
     mclim = 1.0
     mllim = 1.0
     icesw = 1.0
     hmin = 1.0E10
     dlxmin = 1.0E10
     lfpr = blank
     conv = blank
     conv1 = blank1
     cname2 = ADJUSTR(cname2)
     cdname2 = ADJUSTR(cdname2)
     kfname2 = ADJUSTR(kfname2)
     TITLE(11) = ' '
     text = ' '
     icpl = 0
     if(.NOT.constituents)then
         if(nbod>0)deallocate(nbodc, nbodn, nbodp)
         nal = 0
         nep = 0
         nss = 0
         nbod = 0
     endif
     do jw = 1, nwb
         gamma(:, US(BS(jw)):DS(BE(jw))) = EXH2O(jw)
     enddo
  !***********************************************************************************************************************************
!**  Task 1.1.2: Miscellaneous Variables                                      
!***********************************************************************************************************************************
 
!    ** Logical controls
 
     new_page = .TRUE.
  !***********************************************************************************************************************************
!**  Task 1.1.2: Miscellaneous Variables                                      
!***********************************************************************************************************************************
 
!    ** Logical controls
 
     volume_warning = .TRUE.
  !***********************************************************************************************************************************
!**  Task 1.1.2: Miscellaneous Variables                                      
!***********************************************************************************************************************************
 
!    ** Logical controls
 
     initialize_graph = .TRUE.
  !***********************************************************************************************************************************
!**  Task 1.1.2: Miscellaneous Variables                                      
!***********************************************************************************************************************************
 
!    ** Logical controls
 
     update_graph = .TRUE.
     ice = .FALSE.
     flux = .FALSE.
     pumpon = .FALSE.
     tdg_gate = .FALSE.
     tdg_spillway = .FALSE.
     internal_weir = .FALSE.
     surface_warning = .FALSE.
     warning_open = .FALSE.
     print_const = .FALSE.
     print_derived = .FALSE.
     error_open = .FALSE.
     limiting_factor = .FALSE.
     head_boundary = .FALSE.
     print_hydro = .FALSE.
     one_layer = .FALSE.
     zero_slope = .TRUE.
     internal_flow = .FALSE.
     dam_inflow = .FALSE.
     dam_outflow = .FALSE.
     head_flow = .FALSE.                                                                                               !TC 08/03/04
     update_rates = .FALSE.                                                                                             !TC 08/03/04
     weir_calc = niw>0
     gates = ngt>0
     pipes = npi>0
     pumps = npu>0
     spillway = nsp>0
     tributaries = ntr>0
     withdrawals = nwd>0
     volume_balance = vbc=='      ON'
     place_qin = pqc=='      ON'
     evaporation = evc=='      ON'
     energy_balance = ebc=='      ON'
     rh_evap = rhevc=='      ON'
     precipitation = prc=='      ON'
     restart_out = rsoc=='      ON'
     interp_tribs = tric=='      ON'
     interp_dtribs = dtric=='      ON'
     interp_head = hdic=='      ON'
     interp_inflow = qinic=='      ON'
     interp_outflow = stric=='      ON'
     interp_withdrawal = wdic=='      ON'
     interp_gate = gtic=='      ON'                       ! cb 8/13/2010
  !INTERP_METEOROLOGY    = METIC       == '      ON'; DOWNSTREAM_OUTFLOW = WDOC   == '      ON'
     interp_meteorology = metic=='      ON'
     if(wdoc=='      ON' .OR. wdoc=='     ONH' .OR. wdoc=='     ONS')          &
      & downstream_outflow = .TRUE.                                                                          ! cb 4/11/18
     celerity_limit = celc=='      ON'
     viscosity_limit = visc=='      ON'
     hydro_plot = hpltc=='      ON'
     print_hydro = hprwbc=='      ON'
     limiting_dlt = hprwbc(1, :)=='      ON'
     fetch_calc = fetchc=='      ON'
     screen_output = scrc=='      ON'
     snapshot = snpc=='      ON'
     contour = cplc=='      ON'
     vector = vplc=='      ON'
     profile = prfc=='      ON'
     spreadsheet = sprc=='      ON'
     time_series = tsrc=='      ON'
     read_radiation = sroc=='      ON'
     ice_calc = icec=='      ON' .OR. icec=='    ONWB'
     interp_extinction = exic=='      ON'
     read_extinction = exc=='      ON'
     no_inflow = qinc=='     OFF'
     no_outflow = qoutc=='     OFF'
     no_heat = heatc=='     OFF'
     no_wind = windc=='     OFF'
     specify_qtr = trc==' SPECIFY'
     dist_tribs = dtrc=='      ON'
     implicit_visc = azslc=='     IMP'
     upwind = sltrc=='  UPWIND'
     ultimate = sltrc=='ULTIMATE'
     term_by_term = slhtc=='    TERM'
     mannings_n = fricc=='    MANN'
     place_qtr = trc==' DENSITY'
     lateral_spillway = latspc/='    DOWN'
     lateral_pump = latpuc/='    DOWN'
     lateral_gate = latgtc/='    DOWN'
     lateral_pipe = latpic/='    DOWN'
     trapezoidal = gridc=='    TRAP'                                                                                   !SW 07/16/04
     epiphyton_calc = constituents .AND. epic=='      ON'
     mass_balance = constituents .AND. mbc=='      ON'
     susp_solids = constituents .AND. CAC(nsss)=='      ON'
     oxygen_demand = constituents .AND. CAC(ndo)=='      ON'
     sediment_calc = constituents .AND. sedcc=='      ON'
     zooplankton_calc = constituents .AND. CAC(nzoos)=='      ON'
     sediment_resuspension = constituents .AND. sedrc=='      ON'
     derived_plot = constituents .AND. cdpltc=='      ON'
     derived_calc = constituents .AND. ANY(cdwbc=='      ON')
     ph_calc = constituents .AND. cdwbc(20, :)=='      ON'
     print_epiphyton = constituents .AND. epiprc=='      ON' .AND.             &
                     & epiphyton_calc
     print_sediment = constituents .AND. sedprc=='      ON' .AND. sediment_calc
     sediment_calc1 = constituents .AND. sedcc1=='      ON'
     sediment_calc2 = constituents .AND. sedcc2=='      ON'
     print_sediment1 = constituents .AND. sedprc1=='      ON' .AND.            &
                     & sediment_calc1
     print_sediment2 = constituents .AND. sedprc2=='      ON' .AND.            &
                     & sediment_calc2
     fresh_water = constituents .AND. wtypec=='   FRESH' .AND. CAC(ntds)       &
                  &=='      ON'
     salt_water = constituents .AND. wtypec=='    SALT' .AND. CAC(ntds)        &
                 &=='      ON'
     constituent_plot = constituents .AND. cpltc=='      ON' .AND.             &
                       &CAC=='      ON'
     detailed_ice = ice_calc .AND. slicec=='  DETAIL'
     leap_year = MOD(year, 4)==0
     ice_computation = ANY(ice_calc)
     end_run = jday>tmend
     update_kinetics = constituents
  ! IF (WEIR_CALC) THEN   ! SW MOVED TO AFTER GEOM DEFINITION BECAUSE OF NEW FEATURE 3/18/16
  !  DO JWR=1,NIW
  !    DO K=2,KMX-1
  !      IF ((K >= KTWR(JWR) .AND. K <= KBWR(JWR))) INTERNAL_WEIR(K,IWR(JWR)) = .TRUE.
  !    END DO
  !  END DO
  !END IF
     where(read_extinction)
         exom = 0.0
         exss = 0.0
     endwhere
     if(constituents)then
         jg_age = 0
         do jg = 1, ngc
                  !SR 7/27/2017
             if(CGQ10(jg)==0.0 .AND. CG0DK(jg)== - 1.0 .AND. CG1DK(jg)         &
              & ==0.0 .AND. CGS(jg)==0.0 .AND. CGLDK(jg)==0.0 .AND. CGKLF(jg)  &
              & ==0.0 .AND. CGCS(jg)==0.0)then
                 water_age_active = .TRUE.
                 jg_age = jg
                 exit
             endif
         enddo
 
         susp_solids = .FALSE.
         flux = flxc=='      ON'
         print_const = cprwbc=='      ON'
         print_derived = cdwbc=='      ON'
         if(ANY(CAC(nsss:nsse)=='      ON'))susp_solids = .TRUE.
         if(ANY(CAC(nsss:nct)=='      ON'))update_rates = .TRUE.
         do ja = 1, nal
             limiting_factor(ja) = constituents .AND. CAC(nas - 1 + ja)        &
                                  &=='      ON' .AND. limc=='      ON'
             ALG_CALC(ja) = CAC(nas - 1 + ja)=='      ON'
         enddo
         do nb = 1, nbod
             BOD_CALC(nb) = CAC(nbods - 1 + nb)=='      ON'
             BOD_CALCP(nb) = CAC(nbods - 1 + nbod + nb)=='      ON'                  ! cb 5/19/2011
             BOD_CALCN(nb) = CAC(nbods - 1 + 2*nbod + nb)=='      ON'                  ! cb 5/19/2011
         enddo
         dsi_calc = CAC(ndsi)=='      ON'                                                              ! cb 10/12/11
         po4_calc = CAC(npo4)=='      ON'                                                              ! cb 10/12/11
         n_calc = CAC(nnh4)=='      ON' .OR. CAC(nno3)=='      ON'                                     ! cb 10/12/11
 
     endif
     jbdam = 0
     cdhs = dhs
     do jb = 1, nbr
         UP_FLOW(jb) = UHS(jb)==0
         DN_FLOW(jb) = dhs(jb)==0
         UP_HEAD(jb) = UHS(jb)/=0
         UH_INTERNAL(jb) = UHS(jb)>0
         if(UP_HEAD(jb))then
             do jjb = 1, nbr
                 if(ABS(UHS(jb))>=US(jjb) .AND. ABS(UHS(jb))<=DS(jjb))then
                     if(ABS(UHS(jb))==DS(jjb))then
                         if(dhs(jjb)==US(jb))then
                             UP_FLOW(jb) = .TRUE.
                             head_flow(jb) = .TRUE.
                             internal_flow(jb) = .TRUE.
                             UP_HEAD(jb) = .FALSE.
                             UH_INTERNAL(jb) = .FALSE.
                         endif
                         if(UHS(jb)<0)then
                             do jjjb = 1, nbr
                                 if(ABS(UHS(jb))==DS(jjjb))exit                                                        ! CB 1/2/05
                             enddo
                             UP_FLOW(jb) = .TRUE.
                             dam_inflow(jb) = .TRUE.                                                                   !TC 08/03/04
                             dam_outflow(jjjb) = .TRUE.                                                                !TC 08/03/04
                             internal_flow(jb) = .TRUE.
                             UP_HEAD(jb) = .FALSE.
                             UHS(jb) = ABS(UHS(jb))
                             jbdam(jjjb) = jb
                         endif
                     endif
                     exit
                 endif
             enddo
         endif
         DH_INTERNAL(jb) = dhs(jb)>0
         DN_HEAD(jb) = dhs(jb)/=0
         UH_EXTERNAL(jb) = UHS(jb)== - 1
         DH_EXTERNAL(jb) = dhs(jb)== - 1
         UQ_EXTERNAL(jb) = UHS(jb)==0
         DQ_EXTERNAL(jb) = dhs(jb)==0
         DQ_INTERNAL(jb) = DQB(jb)>0
         UQ_INTERNAL(jb) = UQB(jb)>0 .AND. .NOT.dam_inflow(jb)                                                         !TC 08/03/04
     enddo
     do jw = 1, nwb
         if(TKELATPRDCONST(jw)>0.0)TKELATPRD(jw) = .TRUE.
         if(STRICK(jw)>0.0)STRICKON(jw) = .TRUE.
         do jb = BS(jw), BE(jw)
             if(UH_EXTERNAL(jb) .OR. DH_EXTERNAL(jb))head_boundary(jw) = .TRUE.
             if(SLOPE(jb)/=0.0)zero_slope(jw) = .FALSE.
         enddo
     enddo
     WHERE(CAC=='     OFF')cpltc = '     OFF'
 
!    Kinetic flux variables
 
     KFNAME(1) = 'TISS settling in - source, kg/day            '
 
!    Kinetic flux variables
 
     KFNAME(2) = 'TISS settling out - sink, kg/day             '
     KFNAME(3) = 'PO4 algal respiration - source, kg/day       '
     KFNAME(4) = 'PO4 algal growth - sink, kg/day              '
     KFNAME(5) = 'PO4 algal net- source/sink, kg/day           '
     KFNAME(6) = 'PO4 epiphyton respiration - source, kg/day   '
     KFNAME(7) = 'PO4 epiphyton growth - sink, kg/day          '
     KFNAME(8) = 'PO4 epiphyton net- source/sink, kg/day       '
     KFNAME(9) = 'PO4 POM decay - source, kg/day               '
     KFNAME(10) = 'PO4 DOM decay - source, kg/day               '
     kf_po4_sed = 12
     kf_po4_sod = 13
     KFNAME(11) = 'PO4 OM decay - source, kg/day                '
     KFNAME(kf_po4_sed) = 'PO4 sediment decay - source, kg/day          '
     KFNAME(kf_po4_sod) = 'PO4 SOD release - source, kg/day             '
     KFNAME(14) = 'PO4 net settling  - source/sink, kg/day      '
     KFNAME(15) = 'NH4 nitrification - sink, kg/day             '
     KFNAME(16) = 'NH4 algal respiration - source, kg/day       '
     KFNAME(17) = 'NH4 algal growth - sink, kg/day              '
     KFNAME(18) = 'NH4 algal net - source/sink, kg/day          '
     KFNAME(19) = 'NH4 epiphyton respiration - source, kg/day   '
     KFNAME(20) = 'NH4 epiphyton growth - sink, kg/day          '
     KFNAME(21) = 'NH4 epiphyton net - source/sink, kg/day      '
     KFNAME(22) = 'NH4 POM decay - source, kg/day               '
     KFNAME(23) = 'NH4 DOM decay  - source, kg/day              '
     KFNAME(24) = 'NH4 OM decay - source, kg/day                '
     kf_nh4_sed = 25
     kf_nh4_sod = 26
     KFNAME(kf_nh4_sed) = 'NH4 sediment decay - source, kg/day          '
     KFNAME(kf_nh4_sod) = 'NH4 SOD release - source, kg/day             '
     KFNAME(27) = 'NO3 denitrification - sink, kg/day           '
     KFNAME(28) = 'NO3 algal growth - sink, kg/day              '
     KFNAME(29) = 'NO3 epiphyton growth - sink, kg/day          '
     KFNAME(30) = 'NO3 sediment uptake - sink, kg/day           '
     KFNAME(31) = 'DSi algal growth - sink, kg/day              '
     KFNAME(32) = 'DSi epiphyton growth - sink, kg/day          '
     KFNAME(33) = 'DSi PBSi decay - source, kg/day              '
     KFNAME(34) = 'DSi sediment decay - source, kg/day          '
     KFNAME(35) = 'DSi SOD release  - source, kg/day            '
     KFNAME(36) = 'DSi net settling - source/sink, kg/day       '
     KFNAME(37) = 'PBSi algal mortality  - source, kg/day       '
     KFNAME(38) = 'PBSi net settling - source/sink, kg/day      '
     KFNAME(39) = 'PBSi decay - sink, kg/day                    '
     KFNAME(40) = 'Fe net settling - source/sink, kg/day        '
     KFNAME(41) = 'Fe sediment release - source, kg/day         '
     KFNAME(42) = 'LDOM decay - sink, kg/day                    '
     KFNAME(43) = 'LDOM decay to RDOM - sink, kg/day            '
     KFNAME(44) = 'RDOM decay - sink, kg/day                    '
     KFNAME(45) = 'LDOM algal mortality - source, kg/day        '
     KFNAME(46) = 'LDOM epiphyton mortality - source, kg/day    '
     KFNAME(47) = 'LPOM decay - sink, kg/day                    '
     KFNAME(48) = 'LPOM decay to RPOM - sink, kg/day            '
     KFNAME(49) = 'RPOM decay - sink, kg/day                    '
     KFNAME(50) = 'LPOM algal production - source, kg/day       '
     KFNAME(51) = 'LPOM epiphyton production - source, kg/day   '
     KFNAME(52) = 'LPOM net settling - source/sink, kg/day      '
     KFNAME(53) = 'RPOM net settling - source/sink, kg/day      '
     KFNAME(54) = 'CBOD decay - sink, kg/day                    '
     KFNAME(55) = 'DO algal production  - source, kg/day        '
     KFNAME(57) = 'DO algal respiration - sink, kg/day          '                                                              ! cb 6/2/2009
     KFNAME(56) = 'DO epiphyton production  - source, kg/day    '
     KFNAME(58) = 'DO epiphyton respiration - sink, kg/day      '                                                              ! cb 6/2/2009
     KFNAME(59) = 'DO POM decay - sink, kg/day                  '
     KFNAME(60) = 'DO DOM decay - sink, kg/day                  '
     KFNAME(61) = 'DO OM decay - sink, kg/day                   '
     KFNAME(62) = 'DO nitrification - sink, kg/day              '
     KFNAME(63) = 'DO CBOD uptake - sink, kg/day                '
     KFNAME(64) = 'DO reaeration - source/sink, kg/day          '
     kf_do_sed = 65
     kf_do_sod = 66
     KFNAME(kf_do_sed) = 'DO sediment uptake - sink, kg/day            '
     KFNAME(kf_do_sod) = 'DO SOD uptake - sink, kg/day                 '
     KFNAME(67) = 'TIC algal uptake - sink, kg/day              '
     KFNAME(68) = 'TIC epiphyton uptake - sink, kg/day          '
     KFNAME(69) = 'Sediment decay - sink, kg/day                '
     KFNAME(70) = 'Sediment algal settling - sink, kg/day       '
     KFNAME(71) = 'Sediment LPOM settling - source,kg/day       '
     KFNAME(72) = 'Sediment net settling - source/sink, kg/day  '
     KFNAME(73) = 'SOD decay - sink, kg/day                     '
 
     KFNAME(74) = 'LDOM P algal mortality - source, kg/day      '
 
     KFNAME(75) = 'LDOM P epiphyton mortality - source, kg/day  '
     KFNAME(76) = 'LPOM P algal production- source, kg/day      '
     KFNAME(77) = 'LPOM P net settling - source/sink, kg/day    '
     KFNAME(78) = 'RPOM P net settling - source/sink, kg/day    '
     KFNAME(79) = 'LDOM P algal mortality - source, kg/day      '
     KFNAME(80) = 'LDOM P epiphyton mortality - source, kg/day  '
     KFNAME(81) = 'LPOM P algal production- source, kg/day      '
     KFNAME(82) = 'LPOM P net settling - source/sink, kg/day    '
     KFNAME(83) = 'RPOM P net settling - source/sink, kg/day    '
     KFNAME(84) = 'Sediment P decay - sink, kg/day              '
     KFNAME(85) = 'Sediment algal P settling - source, kg/day   '
     KFNAME(86) = 'Sediment P LPOM settling - source,kg/day     '
     KFNAME(87) = 'Sediment net P settling - source/sink, kg/day'
     KFNAME(88) = 'Sediment epiphyton P settling - source,kg/day'
     KFNAME(89) = 'Sediment N decay - sink, kg/day              '
     KFNAME(90) = 'Sediment algal N settling - source, kg/day   '
     KFNAME(91) = 'Sediment N LPOM settling - source,kg/day     '
     KFNAME(92) = 'Sediment net N settling - source/sink, kg/day'
     KFNAME(93) = 'Sediment epiphyton N settling - source,kg/day'
     KFNAME(94) = 'Sediment C decay - sink, kg/day              '
     KFNAME(95) = 'Sediment algal C settling - source, kg/day   '
     KFNAME(96) = 'Sediment C LPOM settling - source,kg/day     '
     KFNAME(97) = 'Sediment net C settling - source/sink, kg/day'
     KFNAME(98) = 'Sediment epiphyton C settling - source,kg/day'
     KFNAME(99) = 'Sediment N denitrification - source, kg/day  '
     KFNAME(100) = 'PO4 macrophyte resp - source, kg/day         '
     KFNAME(101) = 'PO4 macrophyte growth - sink, kg/day         '
     KFNAME(102) = 'NH4 macrophyte resp - source, kg/day         '
     KFNAME(103) = 'NH4 macrophyte growth - sink, kg/day         '
     KFNAME(104) = 'LDOM macrophyte mort  - source, kg/day       '
     KFNAME(105) = 'LPOM macrophyte mort  - source, kg/day       '
     KFNAME(106) = 'RPOM macrophyte mort  - source, kg/day       '
     KFNAME(107) = 'DO  macrophyte production  - source, kg/day  '
     KFNAME(108) = 'DO  macrophyte respiration - sink, kg/day    '
     KFNAME(109) = 'TIC macrophyte growth/resp  - S/S, kg/day    '
     KFNAME(110) = 'CBOD settling - sink, kg/day                 '
     KFNAME(111) = 'Sediment CBOD settling - source, kg/day      '
     KFNAME(112) = 'Sediment CBOD P settling - source, kg/day    '
     KFNAME(113) = 'Sediment CBOD N settling - source, kg/day    '
     KFNAME(114) = 'Sediment CBOD C settling - source, kg/day    '
     KFNAME(115) = 'Sediment Burial - sink, kg/day               '
     kf_sed_pburial = 116
     kf_sed_nburial = 117
     KFNAME(kf_sed_pburial) = 'Sediment P Burial - sink, kg/day             '
     KFNAME(kf_sed_nburial) = 'Sediment N Burial - sink, kg/day             '
     KFNAME(118) = 'Sediment C Burial - sink, kg/day             '
     KFNAME(119) = 'CBOD P settling - sink, kg/day               '
     KFNAME(120) = 'CBOD N settling - sink, kg/day               '
     KFNAME(121) = 'CO2 gas exchange air/water interface, kg/day '
     KFNAME(122) = 'DO H2S decay - sink, kg/day                  '
     KFNAME(123) = 'DO CH4 decay - sink, kg/day                  '
     KFNAME(124) = 'H2S gas exchange air/water interface, kg/day '
     KFNAME(125) = 'CH4 gas exchange air/water interface, kg/day '
     KFNAME(126) = 'H2S decay - sink, kg/day                     '
     KFNAME(127) = 'CH4 decay - sink, kg/day                     '
     KFNAME(128) = 'C to Sed. Diagenesis module - source, kg/day '
     KFNAME(129) = 'N to Sed. Diagenesis module - source, kg/day '
     KFNAME(130) = 'P to Sed. Diagenesis module - source, kg/day '
     KFNAME(131) = 'DO sediment diagenesis uptake - sink, kg/day '
     KFNAME(132) = 'Fe(II) oxidation water column - sink, kg/day '
     KFNAME(133) = 'DO Fe(II) oxidation water col.- sink, kg/day '
     KFNAME(134) = 'FeOOH settling from water col. - sink, kg/day'
     KFNAME(135) = 'MnO2 settling from water col. - sink, kg/day '
     KFNAME(136) = 'Mn(II) oxidation water column - sink, kg/day '
     KFNAME(137) = 'DO Mn(II) oxidation water col.- sink, kg/day '
     KFNAME(138) = 'Labile standing biomass decay- sink, kg/day  '
     KFNAME(139) = 'Refract. stand. biomass decay- sink, kg/day  '
 
 
 
!    Convert rates from per-day to per-second
 
     if(constituents)then
         ae = ae/day
         am = am/day
         ar = ar/day
         ag = ag/day
         as = as/day
         ee = ee/day
         em = em/day
         er = er/day
         eg = eg/day
         eb = eb/day
         CGS = CGS/day
         CG0DK = CG0DK/day
         CG1DK = CG1DK/day
         sss = sss/day
         fes = fes/day
         psis = psis/day
         poms = poms/day
         sdk = sdk/day
         nh4dk = nh4dk/day
         no3dk = no3dk/day
         no3s = no3s/day
         psidk = psidk/day
         lrddk = lrddk/day
         lrpdk = lrpdk/day
         ldomdk = ldomdk/day
         lpomdk = lpomdk/day
         rdomdk = rdomdk/day
         rpomdk = rpomdk/day
         kbod = kbod/day
         seds = seds/day                                                                                   !v3.5
         sdk1 = sdk1/day
         sdk2 = sdk2/day                          ! Amaila
         sedb = sedb/day
                    !CB 11/27/06
         cbods = cbods/day
                         !CB 7/23/07
!        SSFLOC = SSFLOC/DAY                                                  
!        !SR 04/21/13
         do jw = 1, nwb
             sod(US(BS(jw)) - 1:DS(BE(jw)) + 1)                                &
               & = (sod(US(BS(jw)) - 1:DS(BE(jw)) + 1)/day)*FSOD(jw)
         enddo
         do j = 1, nep
             ebr(:, :, j) = eb(j)
         enddo
 
         mg = mg/day
         mr = mr/day
         mm = mm/day
         zg = zg/day
         zr = zr/day
         zm = zm/day
 
 
     endif
 
!    Convert slope to angle alpha in radians
 
     alpha = ATAN(SLOPE)
     sina = SIN(alpha)
     sinac = SIN(ATAN(slopec))
     cosa = COS(alpha)
 
!    Time and printout control variables
 
     if(.NOT.restart_in)then
         jday = tmstrt
         eltm = tmstrt*day
         dlt = DLTMAX(1)
         dlts = dlt
         mindlt = dlt
         nit = 0
         nv = 0
         dltdp = 1
         rsodp = 1
         tsrdp = 1
         snpdp = 1
         vpldp = 1
         prfdp = 1
         sprdp = 1
         cpldp = 1
         scrdp = 1
         flxdp = 1
         wdodp = 1
         nxtsediag = tmstrt
         do jw = 1, nwb
             do j = 1, nod
                 if(tmstrt>SNPD(j, jw))SNPD(j, jw) = tmstrt
                 if(tmstrt>PRFD(j, jw))PRFD(j, jw) = tmstrt
                 if(tmstrt>SPRD(j, jw))SPRD(j, jw) = tmstrt
                 if(tmstrt>CPLD(j, jw))CPLD(j, jw) = tmstrt
                 if(tmstrt>VPLD(j, jw))VPLD(j, jw) = tmstrt
                 if(tmstrt>SCRD(j, jw))SCRD(j, jw) = tmstrt
                 if(tmstrt>FLXD(j, jw))FLXD(j, jw) = tmstrt
             enddo
             NXTMSN(jw) = SNPD(snpdp(jw), jw)
             NXTMPR(jw) = PRFD(prfdp(jw), jw)
             NXTMSP(jw) = SPRD(sprdp(jw), jw)
             NXTMCP(jw) = CPLD(cpldp(jw), jw)
             NXTMVP(jw) = VPLD(vpldp(jw), jw)
             NXTMSC(jw) = SCRD(scrdp(jw), jw)
             NXTMFL(jw) = FLXD(flxdp(jw), jw)
         enddo
         do j = 1, nod
             if(tmstrt>TSRD(j))TSRD(j) = tmstrt
             if(tmstrt>WDOD(j))WDOD(j) = tmstrt
             if(tmstrt>RSOD(j))RSOD(j) = tmstrt
             if(tmstrt>DLTD(j))DLTD(j) = tmstrt
         enddo
         nxtmts = TSRD(tsrdp)
         nxtmrs = RSOD(rsodp)
         nxtmwd = WDOD(wdodp)
 
         nxtmwd_sec = WDOD(wdodp)*86400.
                                       ! cb 4/6/18 frequency test seconds
 
         if(bioexp)then
             biodp = 1
                  ! MLM BIOEXP -BIOENERGETICS
             nxbio = BIOD(1)
                        !BIOD(BIODP) ! MLM BIOEXP  BIOENERGETICS
             nxtbio = BIOD(1)
                         ! MLM INITIALIZE THE OUTPUT VARIABLES
             BIOD(nbio + 1:nod) = tmend + 1.0
                                     ! MLM BIOEXP BIOENERGETICS
         endif
     endif
     TSRD(ntsr + 1:nod) = tmend + 1.0
     WDOD(nwdo + 1:nod) = tmend + 1.0
     RSOD(nrso + 1:nod) = tmend + 1.0
     DLTD(ndlt + 1:nod) = tmend + 1.0
 
     do jw = 1, nwb
         SNPD(NSNP(jw) + 1:nod, jw) = tmend + 1.0
         PRFD(NPRF(jw) + 1:nod, jw) = tmend + 1.0
         SPRD(NSPR(jw) + 1:nod, jw) = tmend + 1.0
         VPLD(NVPL(jw) + 1:nod, jw) = tmend + 1.0
         CPLD(NCPL(jw) + 1:nod, jw) = tmend + 1.0
         SCRD(NSCR(jw) + 1:nod, jw) = tmend + 1.0
         FLXD(NFLX(jw) + 1:nod, jw) = tmend + 1.0
     enddo
     jdayg = jday
     jdaynx = jdayg + 1
     nxtvd = jday
     do j = 1, nod
                  ! SW 12/9/2016
         if(DLTD(j)/=DLTD(j + 1))exit
         dltdp = dltdp + 1
     enddo
 
     dltmaxx = DLTMAX(dltdp)
     dltff = DLTF(dltdp)
     curmax = DLTMAX(dltdp)/DLTF(dltdp)
 
!    Hydraulic structures
 
     if(spillway)then
         do js = 1, nsp
             if(lateral_spillway(js))then
                 if(IDSP(js)/=0)then
                     tributaries = .TRUE.
                     withdrawals = .TRUE.
                 else
                     withdrawals = .TRUE.
                 endif
             endif
             do jb = 1, nbr
                 if(IUSP(js)>=US(jb) .AND. IUSP(js)<=DS(jb))exit
             enddo
             JBUSP(js) = jb
             if(IUSP(js)==DS(JBUSP(js)) .AND. .NOT.lateral_spillway(js))       &
              & nst = nst + 1
             do jw = 1, nwb
                 if(jb>=BS(jw) .AND. jb<=BE(jw))exit
             enddo
             JWUSP(js) = jw
             if(IDSP(js)>0)then
                 do jb = 1, nbr
                     if(IDSP(js)>=US(jb) .AND. IDSP(js)<=DS(jb))exit
                 enddo
                 JBDSP(js) = jb
                 do jw = 1, nwb
                     if(jb>=BS(jw) .AND. jb<=BE(jw))exit
                 enddo
                 JWDSP(js) = jw
             else
                 JBDSP(js) = 1
                 JWDSP(js) = 1
             endif
         enddo
     endif
     if(pipes)then
         do jp = 1, npi
             if(lateral_pipe(jp))then
                 if(IDPI(jp)/=0)then
                     tributaries = .TRUE.
                     withdrawals = .TRUE.
                 else
                     withdrawals = .TRUE.
                 endif
             endif
             do jb = 1, nbr
                 if(IUPI(jp)>=US(jb) .AND. IUPI(jp)<=DS(jb))exit
             enddo
             JBUPI(jp) = jb
             if(IUPI(jp)==DS(JBUPI(jp)) .AND. .NOT.lateral_pipe(jp))nst = nst +&
              & 1
             do jw = 1, nwb
                 if(jb>=BS(jw) .AND. jb<=BE(jw))exit
             enddo
             JWUPI(jp) = jw
             if(IDPI(jp)>0)then
                 do jb = 1, nbr
                     if(IDPI(jp)>=US(jb) .AND. IDPI(jp)<=DS(jb))exit
                 enddo
                 JBDPI(jp) = jb
                 do jw = 1, nwb
                     if(jb>=BS(jw) .AND. jb<=BE(jw))exit
                 enddo
                 JWDPI(jp) = jw
             else
                 JBDPI(jp) = 1
                 JWDPI(jp) = 1
             endif
         enddo
     endif
     if(gates)then
         do jg = 1, ngt
             if(lateral_gate(jg))then
                 if(IDGT(jg)/=0)then
                     tributaries = .TRUE.
                     withdrawals = .TRUE.
                 else
                     withdrawals = .TRUE.
                 endif
             endif
             do jb = 1, nbr
                 if(IUGT(jg)>=US(jb) .AND. IUGT(jg)<=DS(jb))exit
             enddo
             JBUGT(jg) = jb
             if(IUGT(jg)==DS(JBUGT(jg)) .AND. .NOT.lateral_gate(jg))nst = nst +&
              & 1
             do jw = 1, nwb
                 if(jb>=BS(jw) .AND. jb<=BE(jw))exit
             enddo
             JWUGT(jg) = jw
             if(IDGT(jg)>0)then
                 do jb = 1, nbr
                     if(IDGT(jg)>=US(jb) .AND. IDGT(jg)<=DS(jb))exit
                 enddo
                 JBDGT(jg) = jb
                 do jw = 1, nwb
                     if(jb>=BS(jw) .AND. jb<=BE(jw))exit
                 enddo
                 JWDGT(jg) = jw
             else
                 JBDGT(jg) = 1  ! SW 3/24/10
                 JWDGT(jg) = 1  ! SW 3/24/10
             endif
         enddo
     endif
     if(pumps)then
         do jp = 1, npu
             if(lateral_pump(jp))then
                 if(IDPU(jp)/=0)then
                     tributaries = .TRUE.
                     withdrawals = .TRUE.
                 else
                     withdrawals = .TRUE.
                 endif
             endif
             do jb = 1, nbr
                 if(IUPU(jp)>=US(jb) .AND. IUPU(jp)<=DS(jb))exit
             enddo
             JBUPU(jp) = jb
             if(IUPU(jp)==DS(JBUPU(jp)) .AND. .NOT.lateral_pump(jp))nst = nst +&
              & 1
             do jw = 1, nwb
                 if(jb>=BS(jw) .AND. jb<=BE(jw))exit
             enddo
             JWUPU(jp) = jw
             if(IDPU(jp)>0)then
                 do jb = 1, nbr
                     if(IDPU(jp)>=US(jb) .AND. IDPU(jp)<=DS(jb))exit
                 enddo
                 JBDPU(jp) = jb
                 do jw = 1, nwb
                     if(jb>=BS(jw) .AND. jb<=BE(jw))exit
                 enddo
                 JWDPU(jp) = jw
             else
                 JBDPU(jp) = 1
                 JWDPU(jp) = 1
             endif
         enddo
     endif
 
     allocate(estr(nst, nbr), wstr(nst, nbr), qstr(nst, nbr), ktsw(nst, nbr),  &
            & kbsw(nst, nbr), sinkc(nst, nbr), point_sink(nst, nbr), qnew(kmx),&
            & tavg(nst, nbr), tavgw(nwd + nsp + ngt + npi + npu),              &
            & cavg(nst, nbr, nct), cdavg(nst, nbr, ndc),                       &
            & cavgw(nwd + nsp + ngt + npi + npu, nct),                         &
            & cdavgw(nwd + nsp + ngt + npi + npu, ndc))
     tavgw = 0.0
     tavg = 0.0
     cavg = 0.0
     cdavg = 0.0
     cavgw = 0.0
     cdavgw = 0.0
 
     qstr = 0.0
     do jb = 1, nbr
         estr(1:NSTR(jb), jb) = estrt(1:NSTR(jb), jb)
         ktsw(1:NSTR(jb), jb) = ktswt(1:NSTR(jb), jb)
         kbsw(1:NSTR(jb), jb) = kbswt(1:NSTR(jb), jb)
         wstr(1:NSTR(jb), jb) = wstrt(1:NSTR(jb), jb)
         sinkc(1:NSTR(jb), jb) = sinkct(1:NSTR(jb), jb)
         point_sink(1:NSTR(jb), jb) = sinkc(1:NSTR(jb), jb)=='   POINT'
                                                                     ! SW 9/27/13
     enddo
     deallocate(estrt, kbswt, ktswt, wstrt, sinkct)
 
!    Active constituents, derived constituents, and fluxes
 
     if(constituents)then
         do jc = 1, nct
             if(CAC(jc)=='      ON')then
                 nac = nac + 1
                 CN(nac) = jc
             endif
             do jb = 1, nbr
                 if(CINBRC(jc, jb)=='      ON')then
                     nacin(jb) = nacin(jb) + 1
                     INCN(nacin(jb), jb) = jc
                 endif
                 if(CDTBRC(jc, jb)=='      ON')then
                     nacdt(jb) = nacdt(jb) + 1
                     DTCN(nacdt(jb), jb) = jc
                 endif
                 if(CPRBRC(jc, jb)=='      ON')then
                     nacpr(jb) = nacpr(jb) + 1
                     PRCN(nacpr(jb), jb) = jc
                 endif
             enddo
             do jt = 1, ntr
                 if(CTRTRC(jc, jt)=='      ON')then
                     nactr(jt) = nactr(jt) + 1
                     TRCN(nactr(jt), jt) = jc
                 endif
             enddo
         enddo
         do jw = 1, nwb
             do jd = 1, ndc
                 if(cdwbc(jd, jw)=='      ON')then
                     nacd(jw) = nacd(jw) + 1
                     CDN(nacd(jw), jw) = jd
                 endif
             enddo
             do jf = 1, nfl
                 if(KFWBC(jf, jw)=='      ON')then
                     naf(jw) = naf(jw) + 1
                     KFCN(naf(jw), jw) = jf
                 elseif(ph_calc(jw) .AND. jf==121)then
                     naf(jw) = naf(jw) + 1
                     KFCN(naf(jw), jw) = jf
                 elseif(jf>=121 .AND. cemarelatedcode)then
                                                      ! CEMA turning on flux output  ! SW 4/14/2017 ONLY TURN ON IF sediment diagenesis IS ON
                     naf(jw) = naf(jw) + 1
                     KFCN(naf(jw), jw) = jf
                 endif
             enddo
         enddo
     endif
 
!    Starting time
 
     deg = CHAR(248) // 'C'
     esc = CHAR(027)
     call DATE_AND_TIME(cdate, cctime)
     do jw = 1, nwb
         TITLE(11) = 'Model run at ' // cctime(1:2) // ':' // cctime(3:4)      &
                    & // ':' // cctime(5:6) // ' on ' // cdate(5:6) // '/' //  &
                   & cdate(7:8) // '/' // cdate(3:4)
         if(restart_in)TITLE(11) = 'Model restarted at ' // cctime(1:2)        &
                                 & // ':' // cctime(3:4) // ':' // cctime(5:6) &
                                 & // ' on ' // cdate(5:6) // '/' // cdate(7:8)&
                                 & // '/' // cdate(3:4)
     enddo
 
     call INITGEOM
                 ! Call initial geometry
 
!    Density related derived constants
 
     rhowcp = rhow*cp
     rhoicp = rhoi*cp
     rhoirl1 = rhoi*rl1
     dlxrho = 0.0D0
     do jw = 1, nwb
         do jb = BS(jw), BE(jw)
             dlxrho(US(jb):DS(jb)) = 0.5D0/(dlxr(US(jb):DS(jb))*rhow)
             if(UP_HEAD(jb))dlxrho(US(jb) - 1) = 0.5D0/(dlxr(US(jb))*rhow)
         enddo
     enddo
 
!    Transport interpolation multipliers
 
     do jw = 1, nwb
         call INTERPOLATION_MULTIPLIERS
     enddo
 
     call INITCOND
                 ! CALL ROUTINE TO SET UP IC
 
     if(weir_calc)then  ! MOVED FROM ABOVE AFTER GEOMETRY SETUP  SW 3/16/18
         do jwr = 1, niw
             if(EKTWR(jwr)==0.0)then
                 do jw = 1, nwb
                     if(IWR(jwr)>=US(BS(jw)) .AND. IWR(jwr)<=DS(BE(jw)))then
                         KTWR(jwr) = ktwb(jw)
                         exit
                     endif
                 enddo
             else
                 KTWR(jwr) = INT(EKTWR(jwr))
             endif
             if(EKBWR(jwr)<=0.0)then
                 do k = KTWR(jwr), kb(IWR(jwr))
                     if(depthb(k, IWR(jwr))>=ABS(EKBWR(jwr)))then
                         KBWR(jwr) = k
                         exit
                     endif
                 enddo
             else
                 KBWR(jwr) = INT(EKBWR(jwr))
             endif
 
             do k = 2, kmx - 1
                 if((k>=KTWR(jwr) .AND. k<=KBWR(jwr)))                         &
                  & internal_weir(k, IWR(jwr)) = .TRUE.
             enddo
         enddo
     endif
 
 
!    Saved variables for autostepping
 
     if(.NOT.restart_in)then
         sz = z
         su = u
         sw = w
         saz = az
         skti = kti
         sbkt = bkt
         savh2 = avh2
         savhr = avhr
     endif
     call GREGORIAN_DATE
!    CALL TIME_VARYING_DATA
!    if(.not. once_through)
     call TIME_VARYING_DATA              ! w2-ressim
!    CALL READ_INPUT_DATA (NXTVD)
!    if(.not. once_through)
     call READ_INPUT_DATA(nxtvd)         ! w2-ressim
     if(constituents)then
         do jw = 1, nwb
             kt = ktwb(jw)
             do jb = BS(jw), BE(jw)
                 iu = US(jb)
                 id = DS(jb)
                 call TEMPERATURE_RATES
                 call KINETIC_RATES
             enddo
         enddo
     endif
 
     end subroutine INIT
