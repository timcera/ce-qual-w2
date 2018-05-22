!*==deallocate_cema.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine DEALLOCATE_CEMA
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
!*** End of declarations rewritten by SPAG
!
 
     deallocate(consolidationtype, constconsolidrate)
     deallocate(constporewtrrate, consolidrateregnfil)
     deallocate(consregsegst, consregsegen, bubblerelwb)       ! SW 7/1/2017
 
     deallocate(bedelevation, bedelevationlayer, bedporosity)
    !DEALLOCATE (sedcellwidth)
     deallocate(consolidregnnum, bedconsolidrate, porewaterrelrate)
     deallocate(cemasedconc)
     deallocate(cemacumpwrelease, cemalayeradded, cemassapplied)
     deallocate(cemacumpwtorelease, cemacumpwreleased)
     deallocate(numcemapwinst)
     deallocate(applycemapwrelease)
     deallocate(cemacumpwreleaserate)
     deallocate(endbedconsolidation)
     deallocate(cematsscopy)
     deallocate(volcema)
 
     deallocate(fftactprdst, fftactprden)
     deallocate(fftlayconc)
 
     deallocate(sdregnpoc_t, sdregnpon_t, sdregnsul_t)
     deallocate(sdregnpop_t)
     deallocate(sdregnh2s_t, sdregnnh3_t, sdregnch4_t)
     deallocate(sdregntic_t, sdregnalk_t, sdregnpo4_t)
     deallocate(sdregnfe2_t, sdregnfeooh_t)
     deallocate(sdregnmn2_t, sdregnmno2_t)
     deallocate(sdregnt_t)
     deallocate(sedbedinitregsegst, sedbedinitregsegen)
     deallocate(sedgenbodname)
     deallocate(sedgenbodinit)
     deallocate(sedgenbodregsegst, sedgenbodregsegen)
     deallocate(sedgenbodconsregsegst, sedgenbodconsregsegen)
     deallocate(sedgenbodregnrate, sedgenbodregntcoeff)
 
     deallocate(sdregnpoc_l_fr, sdregnpoc_r_fr, sdregnpon_l_fr)
     deallocate(sdregnpon_r_fr, sdregnpw_diffcoeff, sdregnox_threshold)
     deallocate(sdregnpop_l_fr, sdregnpop_r_fr)
     deallocate(sdregnae_nh3_no3_l, sdregnae_nh3_no3_h, sdregnae_no3_n2_l)
     deallocate(sdregnae_no3_n2_h, sdregnan_no3_n2, sdregnae_ch4_co2)
     deallocate(sdregnae_hs_nh4_nit, sdregnae_hs_o2_nit, sdregn_theta_pw)
     deallocate(sdregn_theta_nh3_no3, sdregn_theta_no3_n2,                     &
              & sdregn_theta_ch4_co2)
     deallocate(sdregn_sulfate_ch4_h2s, sdregnae_h2s_so4, sdregn_theta_h2s_so4)
     deallocate(sdregn_normconst_h2s_so4, sdregn_minrate_pon_lab,              &
              & sdregn_minrate_pon_ref)
     deallocate(sdregn_minrate_pon_ine, sdregn_minrate_poc_lab,                &
              & sdregn_minrate_poc_ref)
     deallocate(sdregn_minrate_poc_ine, sdregn_theta_pon_lab,                  &
              & sdregn_theta_pon_ref)
     deallocate(sdregn_theta_pon_ine, sdregn_theta_poc_lab,                    &
              & sdregn_theta_poc_ref)
     deallocate(sdregn_theta_poc_ine, sdregn_ch4compmethod,                    &
              & sdregn_pomresuspmethod)
     deallocate(sdregn_theta_pop_lab, sdregn_theta_pop_ref,                    &
              & sdregn_theta_pop_ine)
     deallocate(sdregn_minrate_pop_lab, sdregn_minrate_pop_ref,                &
              & sdregn_minrate_pop_ine)
     deallocate(sedbeddiarcregsegst, sedbeddiarcregsegen)
 
     deallocate(cemamft_randc_regn, cemamft_incond_regn, mftsedflxvars,        &
              & cema_sd_vars)
     deallocate(sd_no3p2, sd_nh3p2, sd_nh3tp2, sd_ch4p2, sd_po4p2, sd_po4tp2)
     deallocate(sd_hsp2, sd_hstp2, sd_poc22, sd_pon22, sd_pop22)
     deallocate(sd_poc2, sd_pon2, sd_pop2, sd_nh3tp, sd_no3p, sd_po4tp,        &
              & sd_hstp)
     deallocate(sd_fpon, sd_fpoc, sd_kdiapon, sd_thtapon, sd_kdiapoc,          &
              & sd_thtapoc)
     deallocate(sd_jpoc, sd_jpon, sd_jpop, sd_nh3, sd_tic, sd_alk, sd_ph,      &
              & sd_tds)
     deallocate(sd_epoc, sd_epon, sd_epop)
     deallocate(sd_t)
     deallocate(sd_denit, sd_jdenit, sd_jo2no3, sd_ch4, sd_hs, sd_po4t)
     deallocate(sd_fe2t, sd_fe2, sd_feooh, sd_mn2t, sd_mn2, sd_mno2)
     deallocate(sd_kdiapop, sd_thtapop, sd_nh3t, sd_no3, sd_hst, sd_po4,       &
              & sd_fpop)
     deallocate(sd_so4conc, sd_phvalue)
     deallocate(sd_aerlayerthick)
     deallocate(h2sdis, h2sgas, ch4dis, ch4gas)
     deallocate(nh4dis, nh4gas, co2dis, co2gas)
     deallocate(bubbleradiussed, presbubbsed, prescritsed)
     deallocate(cgsed, c0sed, ctsed)
     deallocate(tconc, tconcp, sconc)
     deallocate(dissolvedgassediments)
     deallocate(crackopen, mftbubbreleased, lastdiffvolume)
     deallocate(bubblescarried, bubblesradius)
     deallocate(bubbleslnumber, bubblesstatus)
     deallocate(bubblesrisev)
     deallocate(bubblesgasconc)
     deallocate(brvoluagas, brrateagas)
     deallocate(firstbubblesrelease, bubblesreleaseallvalue)
     deallocate(brrateagasnet)
     deallocate(bubblesatsurface)
     deallocate(bottomturbulence)
!!Sediment Generic BOD
     deallocate(sedgenbodconc, sedgenboddecayrate)
     deallocate(sedgenbodconsrate, sedgenbodconstcoeff)
     deallocate(sedgenbodss, sdpflux, sdnh4flux, sdno3flux)
 
 
     end subroutine DEALLOCATE_CEMA
