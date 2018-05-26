!*==endsimulation.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     subroutine ENDSIMULATION
 
  ! CEMA testing start
  ! CEMA testing end
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
     use TRIDIAG_V
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: ifile
     external RESTART_OUTPUT
!
!*** End of declarations rewritten by SPAG
!
!***********************************************************************************************************************************
!*   Task 3: End Simulation                                                   
!***********************************************************************************************************************************
 
 ! CEMA testing start
   !real jcs
   !!jcs=jcinzz*2.67/1000.0  ! converting C flux to DO flux, assuming 2.67 gO/gC
   !jcs=SD_jctest*2.67  ! converting C flux to DO flux, assuming 2.67 gO/gC
   ! !write(1081,'(5g12.5)')xjnh4,Jcinzz/1000.0,Jcs,MFTSedFlxVars(2,26)
   ! write(1081,'(5g12.5)')xjnh4,SD_Jctest,Jcs,MFTSedFlxVars(2,26)
 ! CEMA testing end
 
!    **
     call DATE_AND_TIME(cdate, cctime)
     if(.NOT.error_open)text = 'Normal termination at ' // cctime(1:2)         &
                             & // ':' // cctime(3:4) // ':' // cctime(5:6)     &
                              & // ' on ' // cdate(5:6) // '/' // cdate(7:8)   &
                              & // '/' // cdate(3:4)
     text = ADJUSTL(TRIM(text))
     call CPU_TIME(current)
     do jw = 1, nwb
         if(SNAPSHOT(jw))then
             write(SNP(jw), '(/A/)')ADJUSTL(TRIM(text))
             write(SNP(jw), '(A)')'Runtime statistics'
             write(SNP(jw), '(2(A,I0))')'  Grid                 = ', imx,      &
                                      & ' x ', kmx
             write(SNP(jw), '(A,I0)')'  Maximum active cells = ', ntacmx,      &
                                    &'  Minimum active cells = ', ntacmn
             write(SNP(jw), '(3(A,F0.1))')'  Segment lengths, m   = ', dlxmin, &
                  &'-', dlxmax
             write(SNP(jw), '(3(A,F0.1))')'  Layer heights, m     = ', hmin,   &
                  &'-', hmax
             write(SNP(jw), '(A)')'  Timestep'
             write(SNP(jw), '(A,I0)')'    Total iterations   = ', nit
             write(SNP(jw), '(A,I0)')'    # of violations    = ', nv
             write(SNP(jw), '(A,F0.2)')'    % violations       = ', FLOAT(nv)  &
                                     & /FLOAT(nit)*100.0
             write(SNP(jw), '(A,I0,A)')'    Average timestep   = ', INT(dltav),&
                                      &' sec'
             write(SNP(jw), '(A,I0,A,F0.2,A)')'  Simulation time      = ',     &
                 & INT(eltmjd), ' days ', (eltmjd - INT(eltmjd))*24.0, ' hours'
             write(SNP(jw), '(A,F0.2,A)')'  Total CPU runtime    = ',          &
                 & (current - start)/60.0, ' min'
             close(SNP(jw))
         endif
    !IF (VECTOR(JW))      CLOSE (VPL(JW))
         if(PROFILE(jw))close(PRF(jw))
         if(SPREADSHEET(jw))close(SPR(jw))
         if(CONTOUR(jw))close(CPL(jw))
     enddo
 
  ! *** DSI W2_TOOL LINKAGE
     if(VECTOR(1))close(VPL(1))
 
     if(time_series)then
         do j = 1, niktsr
             close(TSR(j))
         enddo
         close(wlfn)
                  ! WL output file  ! SW 9/25/13
     endif
     if(warning_open)then
         close(wrn)
     else
         close(wrn, status = 'DELETE')
     endif
     if(error_open)then
         close(w2err)
     else
         close(w2err, status = 'DELETE')
     endif
     do j = 40, nopen
         close(j)
     enddo
     do jw = 1, nwb
         if(VOLUME_BALANCE(jw) .AND. CONTOUR(jw))then
                                                    ! SW 2/19/16
             close(flowbfn)
                       ! flowbal file
             exit
         endif
     enddo
     do jw = 1, nwb
         if(MASS_BALANCE(jw) .AND. derived_calc .AND. CONTOUR(jw))then
                                                                  ! SW 2/19/16
             close(massbfn)
                       ! MASS BALANCE file
             exit
         endif
     enddo
 
     if(selectc=='      ON')then         ! SW 9/25/13 New Section on closing files
         ifile = 1949
         do jb = 1, nbr
             if(NSTR(jb)>0)then
                 ifile = ifile + 1
                 close(ifile)
             endif
         enddo
         if(nwd>0)then
             ifile = ifile + 1
             close(ifile)
         endif
         do jw = 1, nwb
                ! sw 4/20/15
             ifile = ifile + 1
             close(ifile)
         enddo
 
     endif
 
     if(downstream_outflow)then
         jfile = 0
         do jwd = 1, niwdo
             close(WDO(jwd, 1))
             close(WDO(jwd, 2))
             if(constituents)close(WDO(jwd, 3))
             if(derived_calc)close(WDO(jwd, 4))
 
        ! Determine the # of withdrawals at the WITH SEG
             do jb = 1, nbr
                     ! structures
                 if(IWDO(jwd)==DS(jb) .AND. NSTR(jb)/=0)then
                     do js = 1, NSTR(jb)
                         jfile = jfile + 1
                         close(WDO2(jfile, 1))
                         close(WDO2(jfile, 2))
                         if(constituents)close(WDO2(jfile, 3))
                         if(derived_calc)close(WDO2(jfile, 4))
                     enddo
                 endif
             enddo
 
             do js = 1, nwd
                     ! withdrawals
                 if(IWDO(jwd)==IWD(js))then
                     jfile = jfile + 1
                     close(WDO2(jfile, 1))
                     close(WDO2(jfile, 2))
                     if(constituents)close(WDO2(jfile, 3))
                     if(derived_calc)close(WDO2(jfile, 4))
                 endif
             enddo
 
             do js = 1, nsp
                     ! spillways
                 if(IWDO(jwd)==IUSP(js))then
                     jfile = jfile + 1
                     close(WDO2(jfile, 1))
                     close(WDO2(jfile, 2))
                     if(constituents)close(WDO2(jfile, 3))
                     if(derived_calc)close(WDO2(jfile, 4))
                 endif
             enddo
 
             do js = 1, npu
                     ! pumps
                 if(IWDO(jwd)==IUPU(js))then
                     jfile = jfile + 1
                     close(WDO2(jfile, 1))
                     close(WDO2(jfile, 2))
                     if(constituents)close(WDO2(jfile, 3))
                     if(derived_calc)close(WDO2(jfile, 4))
                 endif
             enddo
 
             do js = 1, npi
                     ! pipes
                 if(IWDO(jwd)==IUPI(js))then
                     jfile = jfile + 1
                     close(WDO2(jfile, 1))
                     close(WDO2(jfile, 2))
                     if(constituents)close(WDO2(jfile, 3))
                     if(derived_calc)close(WDO2(jfile, 4))
                 endif
             enddo
 
             do js = 1, ngt
                     ! gates
                 if(IWDO(jwd)==IUGT(js))then
                     jfile = jfile + 1
                     close(WDO2(jfile, 1))
                     close(WDO2(jfile, 2))
                     if(constituents)close(WDO2(jfile, 3))
                     if(derived_calc)close(WDO2(jfile, 4))
                 endif
             enddo
 
         enddo
     endif
 
 
 
     if(error_open)then
         open(w2err, file = 'W2Errordump.opt', status = 'unknown')
         write(w2err, *)'JDAY', jday, 'SZ', sz, 'Z', z, 'H2KT', h2(kt, 1:imx), &
                       &'H1KT', h1(kt, 1:imx), 'BHR1', bhr1(kt, 1:imx), 'BHR2',&
                      & bhr2(kt, 1:imx), 'WSE', elws, 'Q', q, 'QC', qc, 'QERR',&
                      & qerr, 'T1', t1(kt, 1:imx), 'T2', t2(kt, 1:imx), 'SUKT',&
                      & su(kt, 1:imx), 'UKT', u(kt, 1:imx), 'QIN', qin, 'QTR', &
                      & qtr, 'QWD', qwd
         close(w2err)
     endif
 
     deallocate(layerchange, tecplot, x1, hab)
     deallocate(TSR, WDO, WDO2, etsr, IWDO, itsr, title, cdac, wsc, estr, wstr,&
              & qstr, ktsw, kbsw, sinkc)
     deallocate(ebc, mbc, pqc, evc, prc, windc, qinc, qoutc, heatc, slhtc,     &
              & qinic, dtric, tric, wdic)
     deallocate(exc, exic, vbc, metic, sltrc, theta, fricc, naf, eltmf, zmin,  &
              & izmin, c2ch, cdch, epch, kfch, apch, anch, alch)
     deallocate(cpltc, hpltc, cmin, cmax, hymin, hymax, cdmin, cdmax, jbdam,   &
              & ilat, cdpltc, qinsum, tinsum, tind)
     deallocate(qold, dtp, dtps, qolds, qind, jss, hdic, qnew, yss, vss, ys,   &
              & vs, vsts, nsprf)
     deallocate(latgtc, latspc, latpic, dynpipe, latpuc, dyngtc, opt, cind,    &
              & cinsum, cdwbc, kfwbc, cprwbc, cinbrc, ctrtrc, cdtbrc, dynpump)                                                                   ! SW 5/10/10
     deallocate(ysts, yst, vst, allim, aplim, anlim, aslim, ellim, eplim,      &
              & enlim, eslim, cssk, c1, c2, z0)
     deallocate(kfs, af, ef, hyd, kf, azslc, stric, cprbrc, cd, kbmax, elkt,   &
              & wind2, visc, celc, dltadd)
     deallocate(qoavr, qimxr, qomxr, reaerc, lat, longit, elbot, bth, vpr, lpr,&
              & nisnp, niprf, nispr, decl)
     deallocate(a00, hh, t2i, ktwb, kbr, ibpr, dlvr, esr, etr, nbl, lprfn,     &
              & extfn, bthfn, metfn)
     deallocate(snpfn, prffn, sprfn, cplfn, vplfn, flxfn, flxfn2, vprfn, afw,  &
              & bfw, cfw, windh, rhevc, fetchc, jbdn)
     deallocate(kbi, macch, gridc, gma, bta, qtot, sedcip, sedcin, sedcic,     &
              & sedcis)                                                                    ! SW 9/27/2007
     deallocate(sdk, fsod, fsed, sedci, sedcc, sedprc, icec, slicec, icethi,   &
              & albedo, hwi, betai, gammai, icemin)
     deallocate(seds, sedb)    !CB 11/28/06
     deallocate(exh2o, beta, exom, exss, dxi, cbhe, tsed, tsedf, fi, icet2,    &
              & azc, azmax)                                                                                      ! QINT,   QOUTT
     deallocate(ax, wtypec, tair, tdew, wind, phi, cloud, cshe, sron, ranlw,   &
              & rb, rc, re, shade)
     deallocate(et, rs, rn, snpc, scrc, prfc, sprc, cplc, vplc, flxc, nxtmcp,  &
              & nxtmvp, nxtmfl, gamma)
     deallocate(nxtmsn, nxtmsc, nxtmpr, nxtmsp, snpdp, scrdp, prfdp, sprdp,    &
              & cpldp, vpldp, flxdp, ncpl, nvpl, nflx)
     deallocate(nsnp, nscr, nprf, nspr, neqn, po4r, partp, nh4dk, nh4r, no3dk, &
              & no3s, fer, fes, cdsum)
     deallocate(co2r, sroc, o2er, o2eg, caq10, cadk, cas, bodp, bodn, bodc,    &
              & kbod, tbod, rbod, dtrc)
     deallocate(ldomdk, rdomdk, lrddk, omt1, omt2, omk1, omk2, lpomdk, rpomdk, &
              & lrpdk, poms, orgp, orgn, orgc)
     deallocate(rcoef1, rcoef2, rcoef3, rcoef4, orgsi, nh4t1, nh4t2, nh4k1,    &
              & nh4k2, no3t1, no3t2, no3k1, no3k2, NSTR)
     deallocate(dsir, psis, psidk, partsi, sodt1, sodt2, sodk1, sodk2, o2nh4,  &
              & o2om, o2ar, o2ag, cg1dk, cgs)
     deallocate(cgq10, cg0dk, cgldk, cgklf, cgcs, cunit, cunit1, cunit2, cac,  &
              & incac, trcac, dtcac, prcac, cname, cname1, cname2, cmult)                                                                      !LCJ 2/26/15
     deallocate(cn, incn, dtcn, prcn, csum, dltmax, qwdo, twdo, sss, sedrc,    &
              & taucr, xbr, fno3sed, dynstruc)
!    DEALLOCATE (SSFLOC, FLOCEQN)
  !DEALLOCATE (SEDCC1,SEDCC2, ICEQSS,SDK1,sdk2,SEDCI1,SEDCI2,SEDPRC1,SEDPRC2,SEDVP1,SEDVP2,SED1,SED2)
  !DEALLOCATE (SEDCC1,SEDCC2, ICEQSS,SDK1,sdk2,SEDCI1,SEDCI2,SEDPRC1,SEDPRC2,SEDVP1,SEDVP2,SED1,SED2,fsedc1,fsedc2,pbiom,nbiom,cbiom)   ! Amaila, cb 6/7/17
     deallocate(sedcc1, sedcc2, iceqss, sdk1, sdk2, sedci1, sedci2, sedprc1,   &
              & sedprc2, sedvp1, sedvp2, sed1, sed2, fsedc1, fsedc2, pbiom,    &
              & nbiom, cbiom, sed1ic, sed2ic, sdfirstadd)                                                                                                       ! cb 9/3/17
     deallocate(iso_sediment1, vert_sediment1, long_sediment1, iso_sediment2,  &
              & vert_sediment2, long_sediment2, print_sediment1,               &
              & print_sediment2)
     deallocate(sediment_calc1, sediment_calc2)
     deallocate(qtavb, qtmxb, bs, be, jbuh, jbdh, tsss, tssb, tssice, esbr,    &
              & etbr, ebri, qdtr, evbr)
     deallocate(qin, pr, qprbr, tin, tout, tpr, tdtr, tpb, nacpr, nacin, nacdt,&
              & nactr, nacd, eldh)
     deallocate(qsum, nout, ktqin, kbqin, eluh, nl, npoint, slope, slopec,     &
              & alpha, cosa, sina, sinac, tdhfn, qotfn, prefn)
     deallocate(cprfn, euhfn, tuhfn, cuhfn, edhfn, qinfn, tinfn, cinfn, cdhfn, &
              & qdtfn, tdtfn, cdtfn, tprfn, volev)
     deallocate(volwd, volsbr, voltbr, dlvol, volg, volsr, voltr, volb, volpr, &
              & voltrb, voldt, voluh, voldh, volin, volice, icebank)
     deallocate(us, DS, cus, uhs, dhs, uqb, dqb, cdhs, volout, tsswd, tssuh,   &
              & tssdh, tssin, tssout)
     deallocate(tssev, tsspr, tsstr, tssdt, sod, elws, bkt, reaer, iceth, ice, &
              & icesw, q, qc, qerr)
     deallocate(kti, srosh, seg, dlxrho, qssum, dlx, dlxr, quh1, qdh1, bi,     &
              & jwuh, jwdh)
     deallocate(a, c, d, f, v, skti, kbmin, ev, qdt, qpr, sbkt, bhrho)
     deallocate(sz, wshx, wshy, wind10, cz, fetch, phi0, fric, adz, hmult,     &
              & fmtc, fmtcd, cname3, cdname3)
     deallocate(z, kb, vnorm, anpr, aneqn, apom, achla, ahsp, ahsn, ahssi)
     deallocate(ac, asi, at1, at2, at3, at4, ak1, ak2, ak3, ak4, exa, asat, ap,&
              & an)
     deallocate(ag, ar, ae, am, as, enpr, eneqn, eg, er, ee, em, eb, esat, ep)
     deallocate(ec, esi, echla, ehsp, ehsn, ehssi, epom, ehs, en, et4, ek1,    &
              & ek2, ek3, ek4)
     deallocate(et1, et2, et3, hname, fmth, kfac, kfname, kfname2, kfcn, c2i,  &
              & trcn, cdn, cdname, cdname2, cdmult)
     deallocate(cmbrs, cmbrt, fetchu, fetchd, iprf, isnp, ispr, bl, lfpr, do3, &
              & sed, tke, palt)
     deallocate(adx, do1, do2, b, conv, conv1, el, dz, dzq, dx, saz, t1, tss,  &
              & qss, bnew, ilayer)                                                                                               ! SW 1/23/06
     deallocate(p, su, sw, bb, br, bh, bhr, vol, hseg, decay, fpfe, fricbr,    &
              & uxbr, uybr)
     deallocate(depthb, depthm, fpss, tuh, tdh, tssuh1, tssuh2, tssdh1, tssdh2,&
              & sedvp, h, epc)
     deallocate(tvp, qinf, qout, kout, voluh2, voldh2, cwdo, cdwdo, cwdoc,     &
              & cdwdoc, cdtot, cpr, cpb, cout)
     deallocate(cin, cdtr, rsod, rsof, dltd, dltf, tsrd, tsrf, wdod, wdof,     &
              & snpd, snpf, sprd, sprf)
     deallocate(scrd, scrf, prfd, prff, cpld, cplf, vpld, vplf, flxd, flxf,    &
              & epic, epici, epiprc, epivp)
     deallocate(cuh, cdh, epm, epd, c1s, cssb, cvp, cssuh1, cssuh2, cssdh2,    &
              & cssdh1, lname, iwr, ktwr, ektwr, ekbwr)
     deallocate(jwusp, jwdsp, qsp, kbwr, ktwd, kbwd, jbwd, gta1, gtb1, gta2,   &
              & gtb2, bgt, IUGT, idgt)
     deallocate(qtr, ttr, kttr, kbtr, egt, egt2, agasgt, bgasgt, cgasgt,       &
              & gasgtc, pugtc, etugt, ebugt, ktugt, kbugt)
     deallocate(pdgtc, etdgt, ebdgt, ktdgt, kbdgt, a1gt, b1gt, g1gt, a2gt,     &
              & b2gt, g2gt, jwugt, jwdgt, qgt)
     deallocate(eqgt, jbugt, jbdgt, jbupi, jbdpi, jwupi, jwdpi, qpi, IUPI,     &
              & idpi, eupi, edpi, wpi, dlxpi, bp)                                                                                  ! SW 5/5/10
     deallocate(etupi, ebupi, ktupi, kbupi, pdpic, etdpi, ebdpi, ktdpi, kbdpi, &
              & fpi, fminpi, pupic, etdsp, ebdsp)
     deallocate(puspc, etusp, ebusp, ktusp, kbusp, pdspc, ktdsp, kbdsp, IUSP,  &
              & idsp, esp, a1sp, b1sp, a2sp)
     deallocate(b2sp, agassp, bgassp, cgassp, eqsp, gasspc, jbusp, jbdsp,      &
              & strtpu, endpu, eonpu, eoffpu, qpu, ppuc)
     deallocate(IUPU, idpu, epu, etpu, ebpu, ktpu, kbpu, jwupu, jwdpu, jbupu,  &
              & jbdpu, pumpon, ktw, kbw)
     deallocate(IWD, kwd, qwd, ewd, itr, qtrfn, ttrfn, ctrfn, eltrt, eltrb,    &
              & trc, jbtr, qtrf, clrb)
     deallocate(ttlb, ttrb, cllb, srlb1, srrb1, srlb2, srrb2, srfjd1, shadei,  &
              & srfjd2, topo, qsw, ctr)                                                                               ! SW 10/17/05
     deallocate(h1, h2, bh1, bh2, bhr1, bhr2, avh1, avh2, savh2, avhr, savhr,  &
              & cbodd)
     deallocate(point_sink, hprwbc, read_extinction, read_radiation)
     deallocate(dist_tribs, upwind, ultimate, fresh_water, salt_water,         &
              & limiting_factor)
     deallocate(uh_external, dh_external, uh_internal, dh_internal,            &
              & uq_internal, dq_internal)
     deallocate(uq_external, dq_external, up_flow, dn_flow, up_head, dn_head)
     deallocate(internal_flow, dam_inflow, dam_outflow, head_flow,             &
              & head_boundary)                                                                                !TC 08/03/04
     deallocate(iso_conc, vert_conc, long_conc, vert_sediment, long_sediment)
     deallocate(iso_sediment, viscosity_limit, celerity_limit, implicit_az,    &
              & one_layer, implicit_visc)
     deallocate(fetch_calc, limiting_dlt, term_by_term, mannings_n, place_qtr, &
              & specify_qtr)
     deallocate(place_qin, print_const, print_hydro, print_sediment,           &
              & energy_balance, MASS_BALANCE)
     deallocate(VOLUME_BALANCE, detailed_ice, ice_calc, allow_ice, ph_calc,    &
              & br_inactive)                                                                                            ! ICE_IN,       RC/SW 4/28/11
     deallocate(bod_calcp, bod_calcn)
     deallocate(evaporation, precipitation, rh_evap, no_inflow, no_outflow,    &
              & no_heat)
     deallocate(iso_temp, vert_temp, long_temp, vert_profile, long_profile,    &
              & no_wind)
     deallocate(SNAPSHOT, PROFILE, VECTOR, CONTOUR, SPREADSHEET, internal_weir)
     deallocate(screen_output, flux, dynamic_shade, trapezoidal, bod_calc,     &
              & alg_calc)
     deallocate(sediment_calc, epiphyton_calc, print_derived, print_epiphyton, &
              & tdg_spillway, tdg_gate, dynsedk)
     deallocate(iso_epiphyton, vert_epiphyton, long_epiphyton,                 &
              & lateral_spillway, lateral_gate, lateral_pump)
     deallocate(iso_macrophyte, vert_macrophyte, long_macrophyte, macrcvp,     &
              & macrclp)                                                                    ! cb 8/21/15
     deallocate(interp_head, interp_withdrawal, interp_extinction,             &
              & interp_dtribs, lateral_pipe, interp_tribs)
     deallocate(interp_outflow, interp_inflow, interp_meteorology,             &
              & constituent_plot, derived_plot, zero_slope)
     deallocate(hydro_plot, sediment_resuspension)
     deallocate(orgpld, orgprd, orgplp, orgprp, orgnld, orgnrd, orgnlp)
     deallocate(icpl, tavg, tavgw, cavg, cavgw, cdavg, cdavgw)
     deallocate(orgnrp)
     deallocate(print_macrophyte, macrophyte_calc, macwbc, conv2)
     deallocate(mac, macrc, mact, macrm, macss)
     deallocate(mgr, mmr, mrr)
     deallocate(smacrc, smacrm)
     deallocate(smact, smac)
     deallocate(mt1, mt2, mt3, mt4, mk1, mk2, mk3, mk4, mg, mr, mm)
     deallocate(mp, mn, mc, psed, nsed, mhsp, mhsn, mhsc, msat)
     deallocate(cddrag, kticol, armac, macwbci, anorm, dwv, dwsa)
     deallocate(mbmp, mmax, mpom, lrpmac, o2mr, o2mg)
     deallocate(macmbrs, macmbrt, ssmacmb)
     deallocate(cw, bic)
     deallocate(mactrmr, mactrmf, mactrm)
     deallocate(mlfpr)
     deallocate(mllim, mplim, mclim, mnlim)
     deallocate(gammaj)
     deallocate(por, volkti, voli, vstem, vstemkt, sarea)
     deallocate(iwind)
                     ! MLM 08/12/05
     deallocate(zg, zm, zeff, prefp, zr, zoomin, zs2p, exz, zt1, zt2, zt3, zt4,&
              & zk1, zk2)
     deallocate(ldompmp, ldomnmp, lpompmp, lpomnmp, rpompmp, rpomnmp, o2zr)
                                                                    ! MLM 06/10/06
     deallocate(mprwbc)                                             ! MLM 06/10/06
     deallocate(exm)                                                ! MLM 06/10/06
     deallocate(ustarbtke, e, erough, arodi, strick, tkelatprdconst, azt, dzt)
     deallocate(firsti, lasti, tkelatprd, strickon, wallpnt, imptke, tkebc)
     deallocate(zk3, zk4, zp, zn, zc, prefa, zmu, tgraze, zrt, zmt, zoorm,     &
              & zoormr, zoormf)                                              ! POINTERS ,ZOO,ZOOSS,
     deallocate(lpzooout, lpzooin, po4zr, nh4zr, dozr, ticzr, agz, agzt)
     deallocate(gtic, bgto, egto)
                                ! CB 8/13/2010
     deallocate(interp_gate)                   ! CB 8/13/2010
     deallocate(zgz, prefz)
                         !OMNIVOROUS ZOOPLANKTON
     deallocate(lpzooinp, lpzooinn, lpzoooutp, lpzoooutn)
     deallocate(sedc, sedn, sedp, sedninflux, sedpinflux, pfluxin, nfluxin)
     deallocate(sedvpc, sedvpp, sedvpn)
     deallocate(sdkv, seddktot)
     deallocate(cbods, kfjw)
     deallocate(bsave, gma1, bta1)
     deallocate(tn_sedsod_nh4, tp_sedsod_po4, tpout, tptrib, tpdtrib, tpwd,    &
              & tppr, tpin, tnout, tntrib, tndtrib, tnwd, tnpr, tnin)                                               ! TP_SEDBURIAL,TN_SEDBURIAL,
     if(nbod>0)deallocate(nbodc, nbodn, nbodp)
 
     if(fishbio)then
         deallocate(bioexpfn, weightnum, c2zoo, volroos, c2w, ibio)
         deallocate(biod, biof, biodp)
     endif
 
     call DEALLOCATE_TIME_VARYING_DATA
     call DEALLOCATE_TRANSPORT
     call DEALLOCATE_KINETICS
     call DEALLOCATE_WATERBODY
     call DEALLOCATE_PIPE_FLOW
     call DEALLOCATE_OPEN_CHANNEL
     if(aeratec=='      ON')call DEALLOCATE_AERATE
     if(selectc=='      ON')call DEALLOCATE_SELECTIVE
     if(selectc=='    USGS')call DEALLOCATE_SELECTIVEUSGS
 
     if(constituents .AND. phbuff_exist)then
         deallocate(sdeni, pki, pksd)
         if(om_buffering)then
             if(omtype=='    DIST')then
                 deallocate(sden, pk, fract)
             else
                 deallocate(sden, pk)
             endif
         endif
     endif
     if(cemarelatedcode)call DEALLOCATE_CEMA
 
 
     end subroutine ENDSIMULATION
