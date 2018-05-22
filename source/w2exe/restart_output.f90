!*==restart_output.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
!************************************
!***********************************************************************************************************************************
!*                                       S U B R O U T I N E    R E S T A R T   O U T P U T                                       **
!***********************************************************************************************************************************
 
     subroutine RESTART_OUTPUT(Rsofn)
     use GLOBAL
     use SCREENC
     use RSTART
     use GDAYC
     use GEOMC
     use KINETIC, ONLY:epm, epd, sedc, sedn, sedp, PH, sdkv
     use TVDC, ONLY:qsum
     use KINETIC, ONLY:sed, pfluxin, nfluxin
     use ZOOPLANKTONC, ONLY:zoo
     use EDDY, ONLY:tke
     use MAIN, ONLY:envirpc
     use ENVIRPMOD
     use LOGICC
     use STRUCTURES
     use MACROPHYTEC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     character(*) :: Rsofn
     intent (in) Rsofn
!
!*** End of declarations rewritten by SPAG
!
 
     open(rso, file = Rsofn, form = 'UNFORMATTED', status = 'UNKNOWN')
     write(rso)nit, nv, kmin, imin, nsprf, cmbrt, zmin, izmin, start, current
     write(rso)dltdp, snpdp, tsrdp, vpldp, prfdp, cpldp, sprdp, rsodp, scrdp,  &
             & flxdp, wdodp
     write(rso)jday, eltm, eltmf, dlt, dltav, dlts, mindlt, jdmin, curmax
     write(rso)nxtmsn, nxtmts, nxtmpr, nxtmcp, nxtmvp, nxtmrs, nxtmsc, nxtmsp, &
             & nxtmfl, nxtmwd
     write(rso)volin, volout, voluh, voldh, volpr, voltrb, voldt, volwd, volev,&
             & volsbr, voltr, volsr, volice, icebank
     write(rso)tssev, tsspr, tsstr, tssdt, tsswd, tssin, tssout, tsss, tssb,   &
             & tssice
     write(rso)tssuh, tssdh, tssuh2, tssdh2, cssuh2, cssdh2, voluh2, voldh2,   &
             & quh1
     write(rso)esbr, etbr, ebri
     write(rso)z, sz, elws, savh2, savhr, h2
     write(rso)ktwb, kti, skti, sbkt
     write(rso)ice, iceth, cuf, qsum
     write(rso)u, w, su, sw, az, saz, dltlim
     write(rso)t1, t2, c1, c2, c1s, sed, kfs, cssk
     write(rso)epd, epm
     write(rso)macmbrt, macrc, smacrc, mac, smac, macrm, macss
     write(rso)sedc, sedn, sedp, zoo, cd ! mlm 10/06
     write(rso)sdkv                      ! MLM 6/10/07
     write(rso)tke                       ! SW 10/4/07
     if(envirpc=='      ON')write(rso)t_class, v_class, c_class, cd_class,     &
                                    & t_tot, t_cnt, sumvolt, v_cnt, v_tot,     &
                                    & c_tot, c_cnt, cd_tot, cd_cnt
     if(pipes)write(rso)ys, vs, vst, yst, dtp, qold
     write(rso)tpout, tptrib, tpdtrib, tpwd, tppr, tpin, tp_sedsod_po4,        &
             & pfluxin, tnout, tntrib, tndtrib, tnwd, tnpr, tnin,              &
             & tn_sedsod_nh4, nfluxin                                                                                              !TP_SEDBURIAL,TN_SEDBURIAL,
     close(rso)
     end subroutine RESTART_OUTPUT
