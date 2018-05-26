!*==input.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     subroutine INPUT
 
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
     use MSCLIB, ONLY:restart_pushed
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     character(8) :: aid
     character(1) :: char1
     integer :: nndc, nproc
     real :: sum
     external RESTART_OUTPUT
!
!*** End of declarations rewritten by SPAG
!
 
           ! enhanced pH buffering
 
!    Title and array dimensions
 
     allocate(title(11))                                                         ! SW 7/13/09
     read(con, '(///(8X,A72))')(title(j), j = 1, 10)
     read(con, '(//8X,5I8,2A8)')nwb, nbr, imx, kmx, nproc, closec                   ! SW 7/31/09
     read(con, '(//8X,8I8)')ntr, nst, niw, nwd, ngt, nsp, npi, npu
     read(con, '(//8X,7I8,a8)')ngc, nss, nal, nep, nbod, nmc, nzp
     read(con, '(//8X,I8,5A8)')nod, selectc, habtatc, envirpc, aeratec, inituwl
  !READ (CON,'(//8X,I8,7A8)')  NOD,SELECTC,HABTATC,ENVIRPC,AERATEC,inituwl,PHBUFC,NCALKC  ! cb 10/25/13
 
     if(nproc==0)nproc = 1                                                              ! SW 7/31/09
!    call omp_set_num_threads(NPROC)   ! set # of processors to NPROC  Moved
!    to INPUT subroutine  TOGGLE FOR DEBUG
     if(selectc=='        ')selectc = '     OFF'
 
 
!    Constituent numbers
 
     ntds = 1
     ngcs = 2
     ngce = ngcs + ngc - 1
     nsss = ngce + 1
     nsse = nsss + nss - 1
     npo4 = nsse + 1
     nnh4 = npo4 + 1
     nno3 = nnh4 + 1
     ndsi = nno3 + 1
     npsi = ndsi + 1
     nfe = npsi + 1
     nldom = nfe + 1
     nrdom = nldom + 1
     nlpom = nrdom + 1
     nrpom = nlpom + 1
     nbods = nrpom + 1
     if(nbod>0)then  ! VARIABLE STOICHIOMETRY FOR CBOD    ! CB 6/6/10
         allocate(nbodc(nbod), nbodp(nbod), nbodn(nbod))
         ibod = nbods
         nbodcs = ibod
         do jcb = 1, nbod
             nbodc(jcb) = ibod
             ibod = ibod + 1
         enddo
         nbodce = ibod - 1
         nbodps = ibod
         do jcb = 1, nbod
             nbodp(jcb) = ibod
             ibod = ibod + 1
         enddo
         nbodpe = ibod - 1
         nbodns = ibod
         do jcb = 1, nbod
             nbodn(jcb) = ibod
             ibod = ibod + 1
         enddo
         nbodne = ibod - 1
     else
         nbodns = 1
         nbodne = 1
         nbodps = 1
         nbodpe = 1
         nbodcs = 1
         nbodce = 1
 
     endif
     nbode = nbods + nbod*3 - 1
                             ! each BOD group has C, N and P groups
     nas = nbode + 1
     nae = nas + nal - 1
     ndo = nae + 1
     ntic = ndo + 1
     nalk = ntic + 1
     nzoos = nalk + 1
     nzooe = nzoos + nzp - 1
     nldomp = nzooe + 1
     nrdomp = nldomp + 1
     nlpomp = nrdomp + 1
     nrpomp = nlpomp + 1
     nldomn = nrpomp + 1
     nrdomn = nldomn + 1
     nlpomn = nrdomn + 1
     nrpomn = nlpomn + 1
     nct = nrpomn
 
!    Constituent, tributary, and widthdrawal totals
 
     ntrt = ntr + ngt + nsp + npi + npu + (nbr - 1)
                                         ! ADDING NBR FOR RESERVOIR FILLING    SW 6/12/2017
     nwdt = nwd + ngt + nsp + npi + npu
     nept = MAX(nep, 1)
     nmct = MAX(nmc, 1)
     nzpt = MAX(nzp, 1)
     if(restart_pushed)ndc = 23
                             ! SW 4/14/2017
  !IF(.NOT.RESTART_PUSHED)NDC=NDC+1   ! SW 10/20/15 THIS ALLOWS FOR THE POSSIBILITY OF TDG AS AN ADDED DERIVED VARIABLE - THIS WILL BE REDUCED IF TDG IS NOT PRESENT
     ndc = ndc + 1
                ! SW 4/14/17  ADDING TDG
     allocate(cdac(ndc), x1(imx), tecplot(nwb))
     allocate(bta1(kmx), gma1(kmx))
     allocate(wsc(imx), kbi(imx))
     allocate(vbc(nwb), ebc(nwb), mbc(nwb), pqc(nwb), evc(nwb), prc(nwb))
     allocate(windc(nwb), qinc(nwb), qoutc(nwb), heatc(nwb), slhtc(nwb))
     allocate(qinic(nbr), dtric(nbr), tric(ntr), wdic(nwd), hdic(nbr),         &
            & metic(nwb))
     allocate(exc(nwb), exic(nwb))
     allocate(sltrc(nwb), theta(nwb), fricc(nwb), naf(nwb), eltmf(nwb), z0(nwb)&
            & )
     allocate(zmin(nwb), izmin(nwb))
     allocate(c2ch(nct), cdch(ndc), epch(nept), macch(nmct), kfch(nfl),        &
            & apch(nal), anch(nal), alch(nal))
     allocate(cpltc(nct), hpltc(nhy), cdpltc(ndc))
     allocate(cmin(nct), cmax(nct), hymin(nhy), hymax(nhy), cdmin(ndc),        &
            & cdmax(ndc))
     allocate(jbdam(nbr), ilat(nwdt))
     allocate(qinsum(nbr), tinsum(nbr), tind(nbr), jss(nbr), qind(nbr))
     allocate(qold(npi), dtp(npi), dtps(npi), qolds(npi))
     allocate(latgtc(ngt), latspc(nsp), latpic(npi), dynpipe(npi), dynpump(npu)&
            & , latpuc(npu), dyngtc(ngt))
     allocate(gtic(ngt), bgto(ngt), egto(ngt))                                              ! cb 8/13/2010
     allocate(interp_gate(ngt))                                                             ! cb 8/13/2010
     allocate(opt(nwb, 7), cind(nct, nbr), cinsum(nct, nbr))
     allocate(cdwbc(ndc, nwb), kfwbc(nfl, nwb), cprwbc(nct, nwb),              &
            & cinbrc(nct, nbr), ctrtrc(nct, ntr))
     allocate(cdtbrc(nct, nbr), cprbrc(nct, nbr))
     allocate(yss(nnpipe, npi), vss(nnpipe, npi), ys(nnpipe, npi),             &
            & vs(nnpipe, npi), vsts(nnpipe, npi))
     allocate(ysts(nnpipe, npi), yst(nnpipe, npi), vst(nnpipe, npi))
     allocate(cbodd(kmx, imx, nbod))
     allocate(allim(kmx, imx, nal), aplim(kmx, imx, nal), anlim(kmx, imx, nal),&
            & aslim(kmx, imx, nal))
     allocate(ellim(kmx, imx, nep), eplim(kmx, imx, nep), enlim(kmx, imx, nep),&
            & eslim(kmx, imx, nep))
     allocate(cssk(kmx, imx, nct), c1(kmx, imx, nct), c2(kmx, imx, nct),       &
            & cd(kmx, imx, ndc), kf(kmx, imx, nfl))
     allocate(kfs(kmx, imx, nfl), af(kmx, imx, nal, 5), ef(kmx, imx, nep, 5),  &
            & hyd(kmx, imx, nhy), kfjw(nwb, nfl))
     allocate(tke(kmx, imx, 3), azt(kmx, imx), dzt(kmx, imx))
     allocate(ustarbtke(imx), e(imx), erough(nwb), arodi(nwb), strick(nwb),    &
            & tkelatprdconst(nwb))
     allocate(firsti(nwb), lasti(nwb), tkelatprd(nwb), strickon(nwb),          &
            & wallpnt(nwb), imptke(nwb), tkebc(nwb))
     allocate(hydro_plot(nhy), constituent_plot(nct), derived_plot(ndc))
     allocate(zero_slope(nwb), dynamic_shade(imx))
     allocate(azslc(nwb))
     allocate(nsprf(nwb))
     allocate(kbmax(nwb), elkt(nwb), wind2(imx))
     allocate(visc(nwb), celc(nwb), dltadd(nwb), reaerc(nwb))
     allocate(qoavr(nwb), qimxr(nwb), qomxr(nwb))
     allocate(lat(nwb), longit(nwb), elbot(nwb))
     allocate(bth(nwb), vpr(nwb), lpr(nwb))
     allocate(nisnp(nwb), niprf(nwb), nispr(nwb))
     allocate(icpl(nwb))
     allocate(tn_sedsod_nh4(nwb), tp_sedsod_po4(nwb), tpout(nwb), tptrib(nwb), &
            & tpdtrib(nwb), tpwd(nwb), tppr(nwb), tpin(nwb), tnout(nwb),       &
            & tntrib(nwb), tndtrib(nwb), tnwd(nwb), tnpr(nwb), tnin(nwb))                                                                                                               ! TP_SEDBURIAL(NWB),TN_SEDBURIAL(NWB),
     allocate(a00(nwb), hh(nwb), decl(nwb))
     allocate(t2i(nwb), ktwb(nwb), kbr(nwb), ibpr(nwb))
     allocate(dlvr(nwb), esr(nwb), etr(nwb), nbl(nwb))
     allocate(lprfn(nwb), extfn(nwb), bthfn(nwb), metfn(nwb), vprfn(nwb))
     allocate(snpfn(nwb), prffn(nwb), sprfn(nwb), cplfn(nwb), vplfn(nwb),      &
            & flxfn(nwb), flxfn2(nwb))
     allocate(afw(nwb), bfw(nwb), cfw(nwb), windh(nwb), rhevc(nwb), fetchc(nwb)&
            & )
     allocate(sdk(nwb), fsod(nwb), fsed(nwb), sedci(nwb), sedcc(nwb),          &
            & sedprc(nwb), seds(nwb), sedb(nwb), dynsedk(nwb))                                                                  !cb 11/28/06
  !ALLOCATE (SDK1(NWB),sdk2(nwb),SEDCI1(NWB),SEDCI2(NWB),SEDPRC1(NWB),SEDPRC2(NWB),SEDCC1(NWB),SEDCC2(NWB),fsedc1(nwb),fsedc2(nwb))   ! Amaila, cb 6/8/17
     allocate(sdk1(nwb), sdk2(nwb), sedci1(nwb), sedci2(nwb), sedprc1(nwb),    &
            & sedprc2(nwb), sedcc1(nwb), sedcc2(nwb), fsedc1(nwb), fsedc2(nwb))                                                      ! cb 6/17/17
     allocate(icec(nwb), slicec(nwb), icethi(nwb), albedo(nwb), hwi(nwb),      &
            & betai(nwb), gammai(nwb), icemin(nwb), icet2(nwb))
     allocate(exh2o(nwb), beta(nwb), exom(nwb), exss(nwb), dxi(nwb), cbhe(nwb),&
            & tsed(nwb), tsedf(nwb), fi(nwb))
     allocate(ax(nwb), wtypec(nwb), jbdn(nwb), azc(nwb), azmax(nwb), gridc(nwb)&
            & )                                                                             !SW 07/14/04    !  QINT(NWB),   QOUTT(NWB),
     allocate(tair(nwb), tdew(nwb), wind(nwb), phi(nwb), cloud(nwb), cshe(imx),&
            & sron(nwb), ranlw(nwb))
     allocate(snpc(nwb), scrc(nwb), prfc(nwb), sprc(nwb), cplc(nwb), vplc(nwb),&
            & flxc(nwb))
     allocate(nxtmsn(nwb), nxtmsc(nwb), nxtmpr(nwb), nxtmsp(nwb), nxtmcp(nwb), &
            & nxtmvp(nwb), nxtmfl(nwb))
     allocate(snpdp(nwb), scrdp(nwb), prfdp(nwb), sprdp(nwb), cpldp(nwb),      &
            & vpldp(nwb), flxdp(nwb))
     allocate(nsnp(nwb), nscr(nwb), nprf(nwb), nspr(nwb), ncpl(nwb), nvpl(nwb),&
            & nflx(nwb))
     allocate(neqn(nwb), po4r(nwb), partp(nwb))
     allocate(nh4dk(nwb), nh4r(nwb))
     allocate(no3dk(nwb), no3s(nwb), fno3sed(nwb))
     allocate(fer(nwb), fes(nwb))
     allocate(co2r(nwb), sroc(nwb))
     allocate(o2er(nept), o2eg(nept))
     allocate(caq10(nwb), cadk(nwb), cas(nwb))
     allocate(bodp(nbod), bodn(nbod), bodc(nbod))
     allocate(kbod(nbod), tbod(nbod), rbod(nbod))
     allocate(ldomdk(nwb), rdomdk(nwb), lrddk(nwb))
     allocate(omt1(nwb), omt2(nwb), omk1(nwb), omk2(nwb))
     allocate(lpomdk(nwb), rpomdk(nwb), lrpdk(nwb), poms(nwb))
     allocate(orgp(nwb), orgn(nwb), orgc(nwb), orgsi(nwb))
     allocate(pbiom(nwb), nbiom(nwb), cbiom(nwb))    ! Amaila, cb 6/8/17
     allocate(rcoef1(nwb), rcoef2(nwb), rcoef3(nwb), rcoef4(nwb))
     allocate(nh4t1(nwb), nh4t2(nwb), nh4k1(nwb), nh4k2(nwb))
     allocate(no3t1(nwb), no3t2(nwb), no3k1(nwb), no3k2(nwb))
     allocate(dsir(nwb), psis(nwb), psidk(nwb), partsi(nwb))
     allocate(sodt1(nwb), sodt2(nwb), sodk1(nwb), sodk2(nwb))
     allocate(o2nh4(nwb), o2om(nwb))
     allocate(o2ar(nal), o2ag(nal))
     allocate(cgq10(ngc), cg0dk(ngc), cg1dk(ngc), cgs(ngc), cgldk(ngc),        &
            & cgklf(ngc), cgcs(ngc))                                                           !LCJ 2/26/15
     allocate(cunit(nct), cunit1(nct), cunit2(nct))
     allocate(cac(nct), incac(nct), trcac(nct), dtcac(nct), prcac(nct))
     allocate(cname(nct), cname1(nct), cname2(nct), cname3(nct), cmult(nct),   &
            & csum(nct))
     allocate(cn(nct))
     allocate(sss(nss), taucr(nss), sedrc(nss))     !,  SSFLOC(NSS), FLOCEQN(NSS))                                           !SR 04/21/13
     allocate(cdsum(ndc))
     allocate(dtrc(nbr))
     allocate(nstr(nbr), xbr(nbr), dynstruc(nbr))
     allocate(qtavb(nbr), qtmxb(nbr))
     allocate(bs(nwb), be(nwb), jbuh(nbr), jbdh(nbr), jwuh(nbr), jwdh(nbr))
     allocate(tsss(nbr), tssb(nbr), tssice(nbr))
     allocate(esbr(nbr), etbr(nbr), ebri(nbr))
     allocate(qin(nbr), pr(nbr), qprbr(nbr), qdtr(nbr), evbr(nbr))
     allocate(tin(nbr), tout(nbr), tpr(nbr), tdtr(nbr), tpb(nbr))
     allocate(nacpr(nbr), nacin(nbr), nacdt(nbr), nactr(ntr), nacd(nwb))
     allocate(qsum(nbr), nout(nbr), ktqin(nbr), kbqin(nbr), eluh(nbr),         &
            & eldh(nbr))
     allocate(nl(nbr), npoint(nbr), slope(nbr), slopec(nbr), alpha(nbr),       &
            & cosa(nbr), sina(nbr), sinac(nbr), ilayer(imx))
     allocate(cprfn(nbr), euhfn(nbr), tuhfn(nbr), cuhfn(nbr), edhfn(nbr),      &
            & tdhfn(nbr), qotfn(nbr), prefn(nbr))
     allocate(qinfn(nbr), tinfn(nbr), cinfn(nbr), cdhfn(nbr), qdtfn(nbr),      &
            & tdtfn(nbr), cdtfn(nbr), tprfn(nbr))
     allocate(volwd(nbr), volsbr(nbr), voltbr(nbr), dlvol(nbr), volg(nwb),     &
            & volsr(nwb), voltr(nwb), volev(nbr), volice(nbr), icebank(imx))
     allocate(volb(nbr), volpr(nbr), voltrb(nbr), voldt(nbr), voluh(nbr),      &
            & voldh(nbr), volin(nbr), volout(nbr))
     allocate(us(nbr), ds(nbr), cus(nbr), uhs(nbr), dhs(nbr), uqb(nbr),        &
            & dqb(nbr), cdhs(nbr))
     allocate(tssev(nbr), tsspr(nbr), tsstr(nbr), tssdt(nbr), tsswd(nbr),      &
            & tssuh(nbr), tssdh(nbr), tssin(nbr), tssout(nbr))
     allocate(et(imx), rs(imx), rn(imx), rb(imx), rc(imx), re(imx), shade(imx))
     allocate(dltmax(nod), qwdo(imx), twdo(imx))                                                                        ! SW 1/24/05
     allocate(sod(imx), elws(imx), bkt(imx), reaer(imx))
     allocate(iceth(imx), ice(imx), icesw(imx))
     allocate(q(imx), qc(imx), qerr(imx), qssum(imx))
     allocate(kti(imx), skti(imx), srosh(imx), seg(imx), dlxrho(imx))
     allocate(dlx(imx), dlxr(imx))
     allocate(a(imx), c(imx), d(imx), f(imx), v(imx), bta(imx), gma(imx))
     allocate(kbmin(imx), ev(imx), qdt(imx), qpr(imx), sbkt(imx), bhrho(imx))
     allocate(sz(imx), wshx(imx), wshy(imx), wind10(imx), cz(imx), fetch(imx), &
            & phi0(imx), fric(imx))
     allocate(z(imx), kb(imx), palt(imx))
     allocate(vnorm(kmx))
     allocate(anpr(nal), aneqn(nal), apom(nal))
     allocate(ac(nal), asi(nal), achla(nal), ahsp(nal), ahsn(nal), ahssi(nal))
     allocate(at1(nal), at2(nal), at3(nal), at4(nal), ak1(nal), ak2(nal),      &
            & ak3(nal), ak4(nal))
     allocate(ag(nal), ar(nal), ae(nal), am(nal), as(nal), exa(nal), asat(nal),&
            & ap(nal), an(nal))
     allocate(enpr(nept), eneqn(nept))
     allocate(eg(nept), er(nept), ee(nept), em(nept), eb(nept), esat(nept),    &
            & ep(nept), en(nept))
     allocate(ec(nept), esi(nept), echla(nept), ehsp(nept), ehsn(nept),        &
            & ehssi(nept), epom(nept), ehs(nept))
     allocate(et1(nept), et2(nept), et3(nept), et4(nept), ek1(nept), ek2(nept),&
            & ek3(nept), ek4(nept))
     allocate(hname(nhy), fmth(nhy), hmult(nhy), fmtc(nct), fmtcd(ndc))
     allocate(kfac(nfl), kfname(nfl), kfname2(nfl), kfcn(nfl, nwb))
     allocate(c2i(nct, nwb), trcn(nct, ntr))
     allocate(cdn(ndc, nwb), cdname(ndc), cdname2(ndc), cdname3(ndc),          &
            & cdmult(ndc))
     allocate(cmbrs(nct, nbr), cmbrt(nct, nbr), incn(nct, nbr), dtcn(nct, nbr),&
            & prcn(nct, nbr))
     allocate(fetchu(imx, nbr), fetchd(imx, nbr))
     allocate(iprf(imx, nwb), isnp(imx, nwb), ispr(imx, nwb), bl(imx, nwb))
     allocate(h1(kmx, imx), h2(kmx, imx), bh1(kmx, imx), bh2(kmx, imx),        &
            & bhr1(kmx, imx), bhr2(kmx, imx), qtot(kmx, imx))
     allocate(savh2(kmx, imx), avh1(kmx, imx), avh2(kmx, imx), avhr(kmx, imx), &
            & savhr(kmx, imx))
     allocate(lfpr(kmx, imx), bi(kmx, imx), bnew(kmx, imx))     ! SW 1/23/06
     allocate(adx(kmx, imx), adz(kmx, imx), do1(kmx, imx), do2(kmx, imx),      &
            & do3(kmx, imx), sed(kmx, imx))
     allocate(b(kmx, imx), conv(kmx, imx), conv1(kmx, imx), el(kmx, imx),      &
            & dz(kmx, imx), dzq(kmx, imx), dx(kmx, imx))
     allocate(p(kmx, imx), su(kmx, imx), sw(kmx, imx), saz(kmx, imx),          &
            & t1(kmx, imx), tss(kmx, imx), qss(kmx, imx))
     allocate(bb(kmx, imx), br(kmx, imx), bh(kmx, imx), bhr(kmx, imx),         &
            & vol(kmx, imx), hseg(kmx, imx), decay(kmx, imx))
     allocate(depthb(kmx, imx), depthm(kmx, imx), fpss(kmx, imx),              &
            & fpfe(kmx, imx), fricbr(kmx, imx), uxbr(kmx, imx), uybr(kmx, imx))
     allocate(quh1(kmx, nbr), qdh1(kmx, nbr), voluh2(kmx, nbr),                &
            & voldh2(kmx, nbr), tuh(kmx, nbr), tdh(kmx, nbr))
     allocate(tssuh1(kmx, nbr), tssuh2(kmx, nbr), tssdh1(kmx, nbr),            &
            & tssdh2(kmx, nbr))
     allocate(tvp(kmx, nwb), sedvp(kmx, nwb), h(kmx, nwb))
  !ALLOCATE (SEDVP1(KMX,NWB),  SEDVP2(KMX,NWB),SED1(KMX,IMX),SED2(KMX,IMX))   ! Amaila
     allocate(sedvp1(kmx, nwb), sedvp2(kmx, nwb), sed1(kmx, imx),              &
            & sed2(kmx, imx), sed1ic(kmx, imx), sed2ic(kmx, imx),              &
            & sdfirstadd(kmx, imx))                                                                                              !  cb 9/3/17
     allocate(qinf(kmx, nbr), qout(kmx, nbr), kout(kmx, nbr))
     allocate(ct(kmx, imx), at(kmx, imx), vt(kmx, imx), dt(kmx, imx),          &
            & gamma(kmx, imx))
     allocate(cwdo(nct, nod), cdwdo(ndc, nod), cwdoc(nct), cdwdoc(ndc),        &
            & cdtot(ndc))
     allocate(cin(nct, nbr), cdtr(nct, nbr), cpr(nct, nbr), cpb(nct, nbr),     &
            & cout(nct, nbr))
     allocate(rsod(nod), rsof(nod), dltd(nod), dltf(nod))
     allocate(tsrd(nod), tsrf(nod), wdod(nod), wdof(nod))
     allocate(snpd(nod, nwb), snpf(nod, nwb), sprd(nod, nwb), sprf(nod, nwb))
     allocate(scrd(nod, nwb), scrf(nod, nwb), prfd(nod, nwb), prff(nod, nwb))
     allocate(cpld(nod, nwb), cplf(nod, nwb), vpld(nod, nwb), vplf(nod, nwb),  &
            & flxd(nod, nwb), flxf(nod, nwb))
     allocate(epic(nwb, nept), epici(nwb, nept), epiprc(nwb, nept))
     allocate(epivp(kmx, nwb, nep), macrcvp(kmx, nwb, nmc),                    &
            & macrclp(kmx, imx, nmc))                                        ! cb 8/21/15
     allocate(cuh(kmx, nct, nbr), cdh(kmx, nct, nbr))
     allocate(epm(kmx, imx, nept), epd(kmx, imx, nept), epc(kmx, imx, nept))
     allocate(c1s(kmx, imx, nct), cssb(kmx, imx, nct), cvp(kmx, nct, nwb))
     allocate(cssuh1(kmx, nct, nbr), cssuh2(kmx, nct, nbr),                    &
            & cssdh2(kmx, nct, nbr), cssdh1(kmx, nct, nbr))
     allocate(read_extinction(nwb), read_radiation(nwb))
     allocate(dist_tribs(nbr), limiting_factor(nal))
     allocate(upwind(nwb), ultimate(nwb))
     allocate(stric(nst, nbr), estrt(nst, nbr), wstrt(nst, nbr),               &
            & ktswt(nst, nbr), kbswt(nst, nbr), sinkct(nst, nbr))
     allocate(fresh_water(nwb), salt_water(nwb), trapezoidal(nwb))                                                     !SW 07/16/04
     allocate(uh_external(nbr), dh_external(nbr), uh_internal(nbr),            &
            & dh_internal(nbr))
     allocate(uq_external(nbr), dq_external(nbr), uq_internal(nbr),            &
            & dq_internal(nbr))
     allocate(up_flow(nbr), dn_flow(nbr), up_head(nbr), dn_head(nbr))
     allocate(internal_flow(nbr), dam_inflow(nbr), dam_outflow(nbr),           &
            & head_flow(nbr), head_boundary(nwb))                                                                      !TC 08/03/04
     allocate(iso_conc(nct, nwb), vert_conc(nct, nwb), long_conc(nct, nwb))
     allocate(iso_sediment(nwb), vert_sediment(nwb), long_sediment(nwb))
     allocate(iso_sediment1(nwb), vert_sediment1(nwb), long_sediment1(nwb))     !Amaila
     allocate(iso_sediment2(nwb), vert_sediment2(nwb), long_sediment2(nwb))     !Amaila
     allocate(viscosity_limit(nwb), celerity_limit(nwb), implicit_az(nwb))
     allocate(fetch_calc(nwb), one_layer(imx), implicit_visc(nwb))
     allocate(limiting_dlt(nwb), term_by_term(nwb), mannings_n(nwb))
     allocate(place_qin(nwb), place_qtr(ntrt), specify_qtr(ntrt))
     allocate(print_const(nct, nwb), print_hydro(nhy, nwb), print_sediment(nwb)&
            & )
     allocate(print_sediment1(nwb), print_sediment2(nwb)) ! Amaila
     allocate(volume_balance(nwb), energy_balance(nwb), mass_balance(nwb))
     allocate(detailed_ice(nwb), ice_calc(nwb), allow_ice(imx),                &
            & br_inactive(nbr))                                                                            !   ICE_IN(NBR),    RC/SW 4/28/11
     allocate(evaporation(nwb), precipitation(nwb), rh_evap(nwb), ph_calc(nwb))
     allocate(no_inflow(nwb), no_outflow(nwb), no_heat(nwb), no_wind(nwb))
     allocate(iso_temp(nwb), vert_temp(nwb), long_temp(nwb), vert_profile(nwb),&
            & long_profile(nwb))
     allocate(snapshot(nwb), profile(nwb), vector(nwb), contour(nwb),          &
            & spreadsheet(nwb))
     allocate(screen_output(nwb), flux(nwb))
     allocate(print_derived(ndc, nwb), print_epiphyton(nwb, nept))
     allocate(sediment_calc(nwb), epiphyton_calc(nwb, nept),                   &
            & sediment_resuspension(nss), bod_calc(nbod), alg_calc(nal))
     allocate(sediment_calc1(nwb), sediment_calc2(nwb))
                                                       ! Amaila
     allocate(bod_calcp(nbod), bod_calcn(nbod))                                              ! cb 5/19/2011
     allocate(tdg_spillway(nwdt, nsp), tdg_gate(nwdt, ngt),                    &
            & internal_weir(kmx, imx))
     allocate(iso_epiphyton(nwb, nept), vert_epiphyton(nwb, nept),             &
            & long_epiphyton(nwb, nept))
     allocate(iso_macrophyte(nwb, nmc), vert_macrophyte(nwb, nmc),             &
            & long_macrophyte(nwb, nmc))                                                     ! cb 8/21/15
     allocate(lateral_spillway(nsp), lateral_gate(ngt), lateral_pump(npu),     &
            & lateral_pipe(npi))
     allocate(interp_head(nbr), interp_withdrawal(nwd), interp_extinction(nwb),&
            & interp_dtribs(nbr))
     allocate(interp_outflow(nst, nbr), interp_inflow(nbr),                    &
            & interp_meteorology(nwb), interp_tribs(ntr))
     allocate(lname(nct + nhy + ndc))
     allocate(iwr(niw), ktwr(niw), kbwr(niw), ektwr(niw), ekbwr(niw))      ! SW 3/18/16
     allocate(jwusp(nsp), jwdsp(nsp), qsp(nsp))
     allocate(ktwd(nwdt), kbwd(nwdt), jbwd(nwdt))
     allocate(gta1(ngt), gtb1(ngt), gta2(ngt), gtb2(ngt))
     allocate(bgt(ngt), iugt(ngt), idgt(ngt), egt(ngt), egt2(ngt))
     allocate(qtr(ntrt), ttr(ntrt), kttr(ntrt), kbtr(ntrt))
     allocate(agasgt(ngt), bgasgt(ngt), cgasgt(ngt), gasgtc(ngt))
     allocate(pugtc(ngt), etugt(ngt), ebugt(ngt), ktugt(ngt), kbugt(ngt))
     allocate(pdgtc(ngt), etdgt(ngt), ebdgt(ngt), ktdgt(ngt), kbdgt(ngt))
     allocate(a1gt(ngt), b1gt(ngt), g1gt(ngt), a2gt(ngt), b2gt(ngt), g2gt(ngt))
     allocate(eqgt(ngt), jbugt(ngt), jbdgt(ngt), jwugt(ngt), jwdgt(ngt),       &
            & qgt(ngt))
     allocate(jbupi(npi), jbdpi(npi), jwupi(npi), jwdpi(npi), qpi(npi), bp(npi)&
            & )                                                                                                 ! SW 5/10/10
     allocate(iupi(npi), idpi(npi), eupi(npi), edpi(npi), wpi(npi), dlxpi(npi),&
            & fpi(npi), fminpi(npi), pupic(npi))
     allocate(etupi(npi), ebupi(npi), ktupi(npi), kbupi(npi), pdpic(npi),      &
            & etdpi(npi), ebdpi(npi), ktdpi(npi), kbdpi(npi))
     allocate(puspc(nsp), etusp(nsp), ebusp(nsp), ktusp(nsp), kbusp(nsp),      &
            & pdspc(nsp), etdsp(nsp), ebdsp(nsp))
     allocate(ktdsp(nsp), kbdsp(nsp), iusp(nsp), idsp(nsp), esp(nsp), a1sp(nsp)&
            & , b1sp(nsp), a2sp(nsp))
     allocate(b2sp(nsp), agassp(nsp), bgassp(nsp), cgassp(nsp), eqsp(nsp),     &
            & gasspc(nsp), jbusp(nsp), jbdsp(nsp))
     allocate(iupu(npu), idpu(npu), epu(npu), strtpu(npu), endpu(npu),         &
            & eonpu(npu), eoffpu(npu), qpu(npu), ppuc(npu))
     allocate(etpu(npu), ebpu(npu), ktpu(npu), kbpu(npu), jwupu(npu),          &
            & jwdpu(npu), jbupu(npu), jbdpu(npu), pumpon(npu))
     allocate(iwd(nwdt), kwd(nwdt), qwd(nwdt), ewd(nwdt), ktw(nwdt), kbw(nwdt))
     allocate(itr(ntrt), qtrfn(ntr), ttrfn(ntr), ctrfn(ntr), eltrt(ntrt),      &
            & eltrb(ntrt), trc(ntrt), jbtr(ntrt), qtrf(kmx, ntrt))
     allocate(ttlb(imx), ttrb(imx), cllb(imx), clrb(imx))
     allocate(srlb1(imx), srrb1(imx), srlb2(imx), srrb2(imx), srfjd1(imx),     &
            & shadei(imx), srfjd2(imx))
     allocate(topo(imx, iang))                                                                                     ! SW 10/17/05
     allocate(qsw(kmx, nwdt), ctr(nct, ntrt), hprwbc(nhy, nwb))
     allocate(ratz(kmx, nwb), curz1(kmx, nwb), curz2(kmx, nwb), curz3(kmx, nwb)&
            & )                                                                   ! SW 5/15/06
     allocate(zg(nzpt), zm(nzpt), zeff(nzpt), prefp(nzpt), zr(nzpt),           &
            & zoomin(nzpt), zs2p(nzpt), exz(nzpt), prefz(nzpt, nzpt))
     allocate(zt1(nzpt), zt2(nzpt), zt3(nzpt), zt4(nzpt), zk1(nzpt), zk2(nzpt),&
            & zk3(nzpt), zk4(nzpt), o2zr(nzpt))
     allocate(zp(nzpt), zn(nzpt), zc(nzpt))
     allocate(prefa(nal, nzpt))
     allocate(po4zr(kmx, imx), nh4zr(kmx, imx))
     allocate(zmu(kmx, imx, nzp), tgraze(kmx, imx, nzp), zrt(kmx, imx, nzp),   &
            & zmt(kmx, imx, nzp))                                                   ! MLM POINTERS:,ZOO(KMX,IMX,NZP),ZOOSS(KMX,IMX,NZP))
     allocate(zoorm(kmx, imx, nzp), zoormr(kmx, imx, nzp),                     &
            & zoormf(kmx, imx, nzp))
     allocate(lpzooout(kmx, imx), lpzooin(kmx, imx), dozr(kmx, imx),           &
            & ticzr(kmx, imx))
     allocate(agz(kmx, imx, nal, nzp), zgz(kmx, imx, nzp, nzp),                &
            & agzt(kmx, imx, nal))                                        !OMNIVOROUS ZOOPLANKTON
     allocate(orgpld(kmx, imx), orgprd(kmx, imx), orgplp(kmx, imx),            &
            & orgprp(kmx, imx), orgnld(kmx, imx), orgnrd(kmx, imx),            &
            & orgnlp(kmx, imx))
     allocate(orgnrp(kmx, imx))
     allocate(ldompmp(kmx, imx), ldomnmp(kmx, imx), lpompmp(kmx, imx),         &
            & lpomnmp(kmx, imx), rpompmp(kmx, imx), rpomnmp(kmx, imx))
     allocate(lpzooinp(kmx, imx), lpzooinn(kmx, imx), lpzoooutp(kmx, imx),     &
            & lpzoooutn(kmx, imx))
     allocate(sedvpp(kmx, nwb), sedvpc(kmx, nwb), sedvpn(kmx, nwb))
     allocate(sedp(kmx, imx), sedn(kmx, imx), sedc(kmx, imx),                  &
            & sedninflux(kmx, imx), sedpinflux(kmx, imx), pfluxin(nwb),        &
            & nfluxin(nwb))
     allocate(sdkv(kmx, imx), seddktot(kmx, imx))
     allocate(sedcip(nwb), sedcin(nwb), sedcic(nwb), sedcis(nwb))
     allocate(cbods(nbod), cbodns(kmx, imx), sedcb(kmx, imx), sedcbp(kmx, imx),&
            & sedcbn(kmx, imx), sedcbc(kmx, imx))
     allocate(print_macrophyte(nwb, nmct), macrophyte_calc(nwb, nmct),         &
            & macwbc(nwb, nmct), conv2(kmx, kmx), mprwbc(nwb, nmct))
     allocate(mac(kmx, imx, nmct), macrc(kmx, kmx, imx, nmct),                 &
            & mact(kmx, kmx, imx), macrm(kmx, kmx, imx, nmct),                 &
            & macss(kmx, kmx, imx, nmct))
     allocate(mgr(kmx, kmx, imx, nmct), mmr(kmx, imx, nmct),                   &
            & mrr(kmx, imx, nmct))
     allocate(smacrc(kmx, kmx, imx, nmct), smacrm(kmx, kmx, imx, nmct))
     allocate(smact(kmx, kmx, imx), smac(kmx, imx, nmct))
     allocate(mt1(nmct), mt2(nmct), mt3(nmct), mt4(nmct), mk1(nmct), mk2(nmct),&
            & mk3(nmct), mk4(nmct), mg(nmct), mr(nmct), mm(nmct))
     allocate(mbmp(nmct), mmax(nmct), cddrag(nmct), dwv(nmct), dwsa(nmct),     &
            & anorm(nmct))
     allocate(mp(nmct), mn(nmct), mc(nmct), psed(nmct), nsed(nmct), mhsp(nmct),&
            & mhsn(nmct), mhsc(nmct), msat(nmct), exm(nmct))
     allocate(o2mg(nmct), o2mr(nmct), lrpmac(nmct), mpom(nmct))
     allocate(kticol(imx), armac(imx), macwbci(nwb, nmct))
     allocate(macmbrs(nbr, nmct), macmbrt(nbr, nmct), ssmacmb(nbr, nmct))
     allocate(cw(kmx, imx), bic(kmx, imx))
     allocate(mactrmr(kmx, imx, nmct), mactrmf(kmx, imx, nmct),                &
            & mactrm(kmx, imx, nmct))
     allocate(mlfpr(kmx, kmx, imx, nmct))
     allocate(mllim(kmx, kmx, imx, nmct), mplim(kmx, imx, nmct),               &
            & mclim(kmx, imx, nmct), mnlim(kmx, imx, nmct))
     allocate(gammaj(kmx, kmx, imx))
     allocate(por(kmx, imx), volkti(imx), voli(kmx, imx), vstem(kmx, imx, nmct)&
            & , vstemkt(imx, nmct), sarea(nmct))
     allocate(iwind(nwb))
     allocate(layerchange(nwb))
     allocate(cbodp(kmx, imx, nbod), cbodn(kmx, imx, nbod))
                                                           ! CB 6/6/10
     allocate(hab(kmx, imx))
     allocate(iceqss(imx))
                         ! CEMA
 ! BIOENERGETICS !mlm
     if(fishbio)then
         allocate(biod(nod), biof(nod), biodp(nod))
 
  ! BIOENERGETICS OUTPUT CARDS
         read(fishbiofn, '(//(8X,A8,2I8))')bioc, nbio, nibio
         allocate(bioexpfn(nibio), weightnum(nibio), c2zoo(kmx, imx, nct),     &
                & volroos(imx), c2w(imx, nct + 2), ibio(nibio))
         c2w = 0.0
         read(fishbiofn, '(//(:8X,9F8.0))')(biod(ii), ii = 1, nbio)
         read(fishbiofn, '(//(:8X,9F8.0))')(biof(ii), ii = 1, nbio)
         read(fishbiofn, '(//(:8X,9I8))')(ibio(ii), ii = 1, nibio)
         read(fishbiofn, '(//(8x,a72))')biofn
         close(fishbiofn)
         weightfn = 'weight.opt'
                            ! for output filename generation
         bioexp = bioc=='      ON'
         if(bioexp)then
             BHEAD(1) = 'Jday99'
             BHEAD(2) = 'Segmnt'
             BHEAD(3) = 'Dpth_m'
             BHEAD(4) = 'T_C'
             BHEAD(5) = 'gamma'
 
             do j = 1, (nzooe - nzoos + 1)
                 write(segnum, '(i0)')j
                 segnum = ADJUSTL(segnum)
                 l = LEN_TRIM(segnum)
                 BHEAD(5 + j) = 'Zoo' // segnum(1:l)
             enddo
!!bhead(6)   = 'Zoo1'
!!bhead(7)   = 'Zoo2'
             BHEAD(5 + j) = 'K'
             BHEAD(6 + j) = 'BH'
             BHEAD(7 + j) = 'EL'
             BHEAD(8 + j) = 'date'
         endif
  ! NVIOL CARD
 ! read(1222,'(//(8x,a8))') nviolc
 ! NVIOL_PRINT = NVIOLC        == '      ON'
 ! if(nviol_print) then
 !   open(1333,file='nviol.dat',status='unknown')
 !	allocate(nviol_loc(kmx,imx))
 !	nviol_loc = 0
 ! end if
     endif
!    Allocate subroutine variables
 
     call TRANSPORT
     call KINETICS
     call WATERBODY
     call OPEN_CHANNEL_INITIALIZE
     call PIPE_FLOW_INITIALIZE
 
!    State variables
 
     tds => c2(:, :, 1)
 
!    State variables
 
     po4 => c2(:, :, npo4)
 
!    State variables
 
     nh4 => c2(:, :, nnh4)
 
!    State variables
 
     no3 => c2(:, :, nno3)
 
!    State variables
 
     dsi => c2(:, :, ndsi)
     psi => c2(:, :, npsi)
     fe => c2(:, :, nfe)
     ldom => c2(:, :, nldom)
     rdom => c2(:, :, nrdom)
     lpom => c2(:, :, nlpom)
     rpom => c2(:, :, nrpom)
     o2 => c2(:, :, ndo)
     tic => c2(:, :, ntic)
     alk => c2(:, :, nalk)
     cg => c2(:, :, ngcs:ngce)
     ss => c2(:, :, nsss:nsse)
     alg => c2(:, :, nas:nae)
     cbod => c2(:, :, nbodcs:nbodce)
     cbodp => c2(:, :, nbodps:nbodpe)
     cbodn => c2(:, :, nbodns:nbodne)                                                                                                                   ! CB 6/6/10
     zoo => c2(:, :, nzoos:nzooe)
     ldomp => c2(:, :, nldomp)
     rdomp => c2(:, :, nrdomp)
     lpomp => c2(:, :, nlpomp)
     rpomp => c2(:, :, nrpomp)
     ldomn => c2(:, :, nldomn)
     rdomn => c2(:, :, nrdomn)
     lpomn => c2(:, :, nlpomn)
     rpomn => c2(:, :, nrpomn)
 
!    State variable source/sinks
 
     cgss => cssk(:, :, ngcs:ngce)
 
!    State variable source/sinks
 
     ssss => cssk(:, :, nsss:nsse)
 
!    State variable source/sinks
 
     po4ss => cssk(:, :, npo4)
 
!    State variable source/sinks
 
     nh4ss => cssk(:, :, nnh4)
     no3ss => cssk(:, :, nno3)
     dsiss => cssk(:, :, ndsi)
     psiss => cssk(:, :, npsi)
     fess => cssk(:, :, nfe)
     ldomss => cssk(:, :, nldom)
     rdomss => cssk(:, :, nrdom)
     lpomss => cssk(:, :, nlpom)
     rpomss => cssk(:, :, nrpom)
     ass => cssk(:, :, nas:nae)
     doss => cssk(:, :, ndo)
     ticss => cssk(:, :, ntic)
     cbodss => cssk(:, :, nbodcs:nbodce)
     cbodpss => cssk(:, :, nbodps:nbodpe)
     cbodnss => cssk(:, :, nbodns:nbodne)                                                                                     ! CB 6/6/10
     zooss => cssk(:, :, nzoos:nzooe)
     ldompss => cssk(:, :, nldomp)
     rdompss => cssk(:, :, nrdomp)
     lpompss => cssk(:, :, nlpomp)
     rpompss => cssk(:, :, nrpomp)
     ldomnss => cssk(:, :, nldomn)
     rdomnss => cssk(:, :, nrdomn)
     lpomnss => cssk(:, :, nlpomn)
     rpomnss => cssk(:, :, nrpomn)
     alkss => cssk(:, :, nalk)        ! enhanced pH buffering
 
!    Derived variables
 
     doc => cd(:, :, 1)
 
!    Derived variables
 
     poc => cd(:, :, 2)
 
!    Derived variables
 
     toc => cd(:, :, 3)
 
!    Derived variables
 
     don => cd(:, :, 4)
 
!    Derived variables
 
     pon => cd(:, :, 5)
 
!    Derived variables
 
     ton => cd(:, :, 6)
     tkn => cd(:, :, 7)
     tn => cd(:, :, 8)
     dop => cd(:, :, 9)
     pop => cd(:, :, 10)
     top => cd(:, :, 11)
     tp => cd(:, :, 12)
     apr => cd(:, :, 13)
     chla => cd(:, :, 14)
     atot => cd(:, :, 15)
     o2dg => cd(:, :, 16)
     totss => cd(:, :, 17)
     tiss => cd(:, :, 18)
     cbodu => cd(:, :, 19)
     ph => cd(:, :, 20)
     co2 => cd(:, :, 21)
     hco3 => cd(:, :, 22)
     co3 => cd(:, :, 23)
     tdg => cd(:, :, 24)                                                                                                    ! SW 10/20/15
 
!    Kinetic fluxes
 
     sssi => kf(:, :, 1)
 
!    Kinetic fluxes
 
     ssso => kf(:, :, 2)
 
!    Kinetic fluxes
 
     po4ar => kf(:, :, 3)
 
!    Kinetic fluxes
 
     po4ag => kf(:, :, 4)
 
!    Kinetic fluxes
 
     po4ap => kf(:, :, 5)
     po4er => kf(:, :, 6)
     po4eg => kf(:, :, 7)
     po4ep => kf(:, :, 8)
     po4pom => kf(:, :, 9)
     po4dom => kf(:, :, 10)
     po4om => kf(:, :, 11)
     po4sd => kf(:, :, 12)
     po4sr => kf(:, :, 13)
     po4ns => kf(:, :, 14)
     nh4d => kf(:, :, 15)
     nh4ar => kf(:, :, 16)
     nh4ag => kf(:, :, 17)
     nh4ap => kf(:, :, 18)
     nh4er => kf(:, :, 19)
     nh4eg => kf(:, :, 20)
     nh4ep => kf(:, :, 21)
     nh4pom => kf(:, :, 22)
     nh4dom => kf(:, :, 23)
     nh4om => kf(:, :, 24)
     nh4sd => kf(:, :, 25)
     nh4sr => kf(:, :, 26)
     no3d => kf(:, :, 27)
     no3ag => kf(:, :, 28)
     no3eg => kf(:, :, 29)
     no3sed => kf(:, :, 30)
     dsiag => kf(:, :, 31)
     dsieg => kf(:, :, 32)
     dsid => kf(:, :, 33)
     dsisd => kf(:, :, 34)
     dsisr => kf(:, :, 35)
     dsis => kf(:, :, 36)
     psiam => kf(:, :, 37)
     psins => kf(:, :, 38)
     psid => kf(:, :, 39)
     fens => kf(:, :, 40)
     fesr => kf(:, :, 41)
     ldomd => kf(:, :, 42)
     lrdomd => kf(:, :, 43)
     rdomd => kf(:, :, 44)
     ldomap => kf(:, :, 45)
     ldomep => kf(:, :, 46)
     lpomd => kf(:, :, 47)
     lrpomd => kf(:, :, 48)
     rpomd => kf(:, :, 49)
     lpomap => kf(:, :, 50)
     lpomep => kf(:, :, 51)
     lpomns => kf(:, :, 52)
     rpomns => kf(:, :, 53)
     cboddk => kf(:, :, 54)
     doap => kf(:, :, 55)
 ! DOAR   => KF(:,:,56); DOEP   => KF(:,:,57); DOER   => KF(:,:,58); DOPOM  => KF(:,:,59); DODOM  => KF(:,:,60)   ! cb 6/2/2009
     doep => kf(:, :, 56)
 ! DOAR   => KF(:,:,56); DOEP   => KF(:,:,57); DOER   => KF(:,:,58); DOPOM  => KF(:,:,59); DODOM  => KF(:,:,60)   ! cb 6/2/2009
     doar => kf(:, :, 57)
 ! DOAR   => KF(:,:,56); DOEP   => KF(:,:,57); DOER   => KF(:,:,58); DOPOM  => KF(:,:,59); DODOM  => KF(:,:,60)   ! cb 6/2/2009
     doer => kf(:, :, 58)
 ! DOAR   => KF(:,:,56); DOEP   => KF(:,:,57); DOER   => KF(:,:,58); DOPOM  => KF(:,:,59); DODOM  => KF(:,:,60)   ! cb 6/2/2009
     dopom => kf(:, :, 59)
 ! DOAR   => KF(:,:,56); DOEP   => KF(:,:,57); DOER   => KF(:,:,58); DOPOM  => KF(:,:,59); DODOM  => KF(:,:,60)   ! cb 6/2/2009
     dodom => kf(:, :, 60)                                                                                       ! cb 9/16/2015
     doom => kf(:, :, 61)
     donit => kf(:, :, 62)
     dobod => kf(:, :, 63)
     doae => kf(:, :, 64)
     dosed => kf(:, :, 65)
     dosod => kf(:, :, 66)
     ticap => kf(:, :, 67)
     ticep => kf(:, :, 68)
     sedd => kf(:, :, 69)
     sedas => kf(:, :, 70)
     sedoms => kf(:, :, 71)
     sedns => kf(:, :, 72)
     sodd => kf(:, :, 73)
 
     ldompap => kf(:, :, 74)
 
     ldompep => kf(:, :, 75)
 
     lpompap => kf(:, :, 76)
 
     lpompns => kf(:, :, 77)
 
     rpompns => kf(:, :, 78)
     ldomnap => kf(:, :, 79)
     ldomnep => kf(:, :, 80)
     lpomnap => kf(:, :, 81)
     lpomnns => kf(:, :, 82)
     rpomnns => kf(:, :, 83)
     seddp => kf(:, :, 84)
     sedasp => kf(:, :, 85)
     sedomsp => kf(:, :, 86)
     sednsp => kf(:, :, 87)
     lpomepp => kf(:, :, 88)
     seddn => kf(:, :, 89)
     sedasn => kf(:, :, 90)
     sedomsn => kf(:, :, 91)
     sednsn => kf(:, :, 92)
     lpomepn => kf(:, :, 93)
     seddc => kf(:, :, 94)
     sedasc => kf(:, :, 95)
     sedomsc => kf(:, :, 96)
     sednsc => kf(:, :, 97)
     lpomepc => kf(:, :, 98)
     sedno3 => kf(:, :, 99)
     po4mr => kf(:, :, 100)
     po4mg => kf(:, :, 101)
     nh4mr => kf(:, :, 102)
     nh4mg => kf(:, :, 103)
     ldommac => kf(:, :, 104)
     rpommac => kf(:, :, 105)
     lpommac => kf(:, :, 106)
     domp => kf(:, :, 107)
     domr => kf(:, :, 108)
     ticmc => kf(:, :, 109)
     cbodns => kf(:, :, 110)
     sedcb => kf(:, :, 111)
     sedcbp => kf(:, :, 112)
     sedcbn => kf(:, :, 113)
     sedcbc => kf(:, :, 114)
     sedbr => kf(:, :, 115)
     sedbrp => kf(:, :, 116)
     sedbrn => kf(:, :, 117)
     sedbrc => kf(:, :, 118)
     co2reaer => kf(:, :, 121)
     cbodnsp => kf(:, :, 119)
     cbodnsn => kf(:, :, 120)                          ! CB 6/6/10
  ! CEMA start
     doh2s => kf(:, :, 122)
  ! CEMA start
     doch4 => kf(:, :, 123)
  ! CEMA start
     h2sreaer => kf(:, :, 124)
  ! CEMA start
     ch4reaer => kf(:, :, 125)
  ! CEMA start
     h2sd => kf(:, :, 126)
     ch4d => kf(:, :, 127)
     sdinc => kf(:, :, 128)
     sdinn => kf(:, :, 129)
     sdinp => kf(:, :, 130)
     dosedia => kf(:, :, 131)
     fe2d => kf(:, :, 132)
     dofe2 => kf(:, :, 133)
     sdinfeooh => kf(:, :, 134)
     sdinmno2 => kf(:, :, 135)
     mn2d => kf(:, :, 136)
     domn2 => kf(:, :, 137)
     sedd1 => kf(:, :, 138)
     sedd2 => kf(:, :, 139)
  ! CEMA end
 
 
!    Algal rate variables
 
     agr => af(:, :, :, 1)
  ! CEMA end
 
 
!    Algal rate variables
 
     arr => af(:, :, :, 2)
  ! CEMA end
 
 
!    Algal rate variables
 
     aer => af(:, :, :, 3)
  ! CEMA end
 
 
!    Algal rate variables
 
     amr => af(:, :, :, 4)
  ! CEMA end
 
 
!    Algal rate variables
 
     asr => af(:, :, :, 5)
     egr => ef(:, :, :, 1)
     err => ef(:, :, :, 2)
     eer => ef(:, :, :, 3)
     emr => ef(:, :, :, 4)
     ebr => ef(:, :, :, 5)
 
!    Hydrodynamic variables
 
     dltlim => hyd(:, :, 1)
 
!    Hydrodynamic variables
 
     u => hyd(:, :, 2)
 
!    Hydrodynamic variables
 
     w => hyd(:, :, 3)
 
!    Hydrodynamic variables
 
     t2 => hyd(:, :, 4)
 
!    Hydrodynamic variables
 
     rho => hyd(:, :, 5)
 
!    Hydrodynamic variables
 
     az => hyd(:, :, 6)
     vsh => hyd(:, :, 7)
     st => hyd(:, :, 8)
     sb => hyd(:, :, 9)
     admx => hyd(:, :, 10)
     dm => hyd(:, :, 11)
     hdg => hyd(:, :, 12)
     admz => hyd(:, :, 13)
     hpg => hyd(:, :, 14)
     grav => hyd(:, :, 15)
 
!    I/O units
 
     snp => opt(:, 1)
 
!    I/O units
 
     prf => opt(:, 2)
 
!    I/O units
 
     vpl => opt(:, 3)
 
!    I/O units
 
     cpl => opt(:, 4)
 
!    I/O units
 
     spr => opt(:, 5)
 
!    I/O units
 
     flx => opt(:, 6)
 
!    I/O units
 
     flx2 => opt(:, 7)
 
!    Zero variables
 
     itr = 0
 
!    Zero variables
 
     jbtr = 0
 
!    Zero variables
 
     kttr = 0
 
!    Zero variables
 
     kbtr = 0
 
!    Zero variables
 
     qtr = 0.0
 
!    Zero variables
 
     ttr = 0.0
 
!    Zero variables
 
     ctr = 0.0
 
!    Zero variables
 
     qtrf = 0.0
 
!    Zero variables
 
     snpd = 0.0
 
!    Zero variables
 
     tsrd = 0.0
     prfd = 0.0
     sprd = 0.0
     cpld = 0.0
     vpld = 0.0
     scrd = 0.0
     flxd = 0.0
     wdod = 0.0
     rsod = 0.0
     eltrb = 0.0
     eltrt = 0.0
 
!    Input file unit numbers
 
     nunit = 40
     do jw = 1, nwb
         bth(jw) = nunit
         vpr(jw) = nunit + 1
         lpr(jw) = nunit + 2
         nunit = nunit + 3
     enddo
     grf = nunit
     nunit = nunit + 1
 
!    Time control cards
 
     read(con, '(//8X,2F8.0,I8)')tmstrt, tmend, year
     read(con, '(//8X,I8,F8.0,a8)')ndlt, dltmin, dltinter
     dltd = 0.0                                                                ! SW 9/28/13 INITIALIZE ARRAY TO NOD SINCE ONLY NDLT ASSIGNED
     read(con, '(//(:8X,9F8.0))')(dltd(j), j = 1, ndlt)
     read(con, '(//(:8X,9F8.0))')(dltmax(j), j = 1, ndlt)
     read(con, '(//(:8X,9F8.0))')(dltf(j), j = 1, ndlt)
     read(con, '(//(8X,3A8))')(visc(jw), celc(jw), dltadd(jw), jw = 1, nwb)
 
!    Grid definition cards
 
     read(con, '(//(8X,7I8,F8.0,F8.0))')(us(jb), ds(jb), uhs(jb), dhs(jb), uqb(&
                                      & jb), dqb(jb), nl(jb), slope(jb),       &
                                      & slopec(jb), jb = 1, nbr)
     read(con, '(//(8X,3F8.0,3I8))')(lat(jw), longit(jw), elbot(jw), bs(jw),   &
                                  & be(jw), jbdn(jw), jw = 1, nwb)
 
!    Initial condition cards
 
     read(con, '(//(8X,2F8.0,2A8))')(t2i(jw), icethi(jw), wtypec(jw), gridc(jw)&
                                  & , jw = 1, nwb)
     read(con, '(//(8X,6A8))')(vbc(jw), ebc(jw), mbc(jw), pqc(jw), evc(jw),    &
                            & prc(jw), jw = 1, nwb)
     read(con, '(//(8X,4A8))')(windc(jw), qinc(jw), qoutc(jw), heatc(jw),      &
                            & jw = 1, nwb)
     read(con, '(//(8X,3A8))')(qinic(jb), dtric(jb), hdic(jb), jb = 1, nbr)
     read(con, '(//(8X,5A8,4F8.0))')(slhtc(jw), sroc(jw), rhevc(jw), metic(jw),&
                                  & fetchc(jw), afw(jw), bfw(jw), cfw(jw),     &
                                  & windh(jw), jw = 1, nwb)
     read(con, '(//(8X,2A8,6F8.0))')(icec(jw), slicec(jw), albedo(jw), hwi(jw),&
                                  & betai(jw), gammai(jw), icemin(jw),         &
                                  & icet2(jw), jw = 1, nwb)
     read(con, '(//(8X,A8,F8.0))')(sltrc(jw), theta(jw), jw = 1, nwb)
     read(con, '(//(8X,6F8.0,A8,F8.0))')(ax(jw), dxi(jw), cbhe(jw), tsed(jw),  &
                                      & fi(jw), tsedf(jw), fricc(jw), z0(jw),  &
                                      & jw = 1, nwb)
     read(con, '(//(8X,2A8,F8.0,I8,F8.0,F8.0,F8.0,F8.0,A8))')                  &
        & (azc(jw), azslc(jw), azmax(jw), tkebc(jw), erough(jw), arodi(jw),    &
        & strick(jw), tkelatprdconst(jw), imptke(jw), jw = 1, nwb)                                           !,PHISET(JW
 
     do jw = 1, nwb
         if(z0(jw)<=0.0)z0(jw) = 0.001
                                     ! SW 11/28/07
         do jb = bs(jw), be(jw)
             do i = us(jb), ds(jb)
                 e(i) = erough(jw)
             enddo
         enddo
     enddo
 
!    Inflow-outflow cards
 
     read(con, '(//(8X,I8,A8))')(nstr(jb), dynstruc(jb), jb = 1, nbr)
     read(con, '(/)')
     do jb = 1, nbr
         read(con, '(:8X,9A8)')(stric(js, jb), js = 1, nstr(jb))
     enddo
     read(con, '(/)')
     do jb = 1, nbr
         read(con, '(:8X,9I8)')(ktswt(js, jb), js = 1, nstr(jb))
     enddo
     read(con, '(/)')
     do jb = 1, nbr
         read(con, '(:8X,9I8)')(kbswt(js, jb), js = 1, nstr(jb))
     enddo
     read(con, '(/)')
     do jb = 1, nbr
         read(con, '(:8X,9A8)')(sinkct(js, jb), js = 1, nstr(jb))
     enddo
     read(con, '(/)')
     do jb = 1, nbr
         read(con, '(:8X,9F8.0)')(estrt(js, jb), js = 1, nstr(jb))
     enddo
     read(con, '(/)')
     do jb = 1, nbr
         read(con, '(:8X,9F8.0)')(wstrt(js, jb), js = 1, nstr(jb))
     enddo
     read(con, '(//(:8X,2I8,6F8.0,A8,A8))')                                    &
        & (iupi(jp), idpi(jp), eupi(jp), edpi(jp), wpi(jp), dlxpi(jp), fpi(jp),&
        & fminpi(jp), latpic(jp), dynpipe(jp), jp = 1, npi)
     read(con, '(//(:8X,A8,2F8.0,2I8))')(pupic(jp), etupi(jp), ebupi(jp), ktupi&
                                      & (jp), kbupi(jp), jp = 1, npi)
     read(con, '(//(:8X,A8,2F8.0,2I8))')(pdpic(jp), etdpi(jp), ebdpi(jp), ktdpi&
                                      & (jp), kbdpi(jp), jp = 1, npi)
     read(con, '(//(:8X,2I8,5F8.0,A8))')(iusp(js), idsp(js), esp(js), a1sp(js),&
                                      & b1sp(js), a2sp(js), b2sp(js),          &
                                      & latspc(js), js = 1, nsp)
     read(con, '(//(:8X,A8,2F8.0,2I8))')(puspc(js), etusp(js), ebusp(js), ktusp&
                                      & (js), kbusp(js), js = 1, nsp)
     read(con, '(//(:8X,A8,2F8.0,2I8))')(pdspc(js), etdsp(js), ebdsp(js), ktdsp&
                                      & (js), kbdsp(js), js = 1, nsp)
     read(con, '(//(:8X,A8,I8,3F8.0))')(gasspc(js), eqsp(js), agassp(js),      &
                                     & bgassp(js), cgassp(js), js = 1, nsp)
     read(con, '(//(:8X,2I8,7F8.0,A8))')(iugt(jg), idgt(jg), egt(jg), a1gt(jg),&
                                      & b1gt(jg), g1gt(jg), a2gt(jg), b2gt(jg),&
                                      & g2gt(jg), latgtc(jg), jg = 1, ngt)
     read(con, '(//(:8X,4F8.0,2A8))')(gta1(jg), gtb1(jg), gta2(jg), gtb2(jg),  &
                                   & dyngtc(jg), gtic(jg), jg = 1, ngt)                                                  ! cb 8/13/2010
     read(con, '(//(:8X,A8,2F8.0,2I8))')(pugtc(jg), etugt(jg), ebugt(jg), ktugt&
                                      & (jg), kbugt(jg), jg = 1, ngt)
     read(con, '(//(:8X,A8,2F8.0,2I8))')(pdgtc(jg), etdgt(jg), ebdgt(jg), ktdgt&
                                      & (jg), kbdgt(jg), jg = 1, ngt)
     read(con, '(//(:8X,A8,I8,3F8.0))')(gasgtc(jg), eqgt(jg), agasgt(jg),      &
                                     & bgasgt(jg), cgasgt(jg), jg = 1, ngt)
     read(con, '(//(:8X,2I8,6F8.0,2A8))')                                      &
        & (iupu(jp), idpu(jp), epu(jp), strtpu(jp), endpu(jp), eonpu(jp),      &
        & eoffpu(jp), qpu(jp), latpuc(jp), dynpump(jp), jp = 1, npu)
     read(con, '(//(:8X,A8,2F8.0,2I8))')(ppuc(jp), etpu(jp), ebpu(jp), ktpu(jp)&
                                      & , kbpu(jp), jp = 1, npu)
     read(con, '(//(:8X,9I8))')(iwr(jw), jw = 1, niw)
     read(con, '(//(:8X,9F8.0))')(ektwr(jw), jw = 1, niw)                   ! SW 3/18/16
     read(con, '(//(:8X,9F8.0))')(ekbwr(jw), jw = 1, niw)                   ! SW 3/18/16
     read(con, '(//(:8X,9A8))')(wdic(jw), jw = 1, nwd)
     read(con, '(//(:8X,9I8))')(iwd(jw), jw = 1, nwd)
     read(con, '(//(:8X,9F8.0))')(ewd(jw), jw = 1, nwd)
     read(con, '(//(:8X,9I8))')(ktwd(jw), jw = 1, nwd)
     read(con, '(//(:8X,9I8))')(kbwd(jw), jw = 1, nwd)
     trc = '      '                                                          ! SW 9/27/13 INITIALIZATION SINCE ALLOCATION IS TO NTRT
     read(con, '(//(:8X,9A8))')(trc(jt), jt = 1, ntr)
     read(con, '(//(:8X,9A8))')(tric(jt), jt = 1, ntr)
     read(con, '(//(:8X,9I8))')(itr(jt), jt = 1, ntr)
     read(con, '(//(:8X,9F8.0))')(eltrt(jt), jt = 1, ntr)
     read(con, '(//(:8X,9F8.0))')(eltrb(jt), jt = 1, ntr)
     read(con, '(//(8X,A8))')(dtrc(jb), jb = 1, nbr)
 
!    Output control cards (excluding constituents)
 
     read(con, '(/)')
     do jh = 1, nhy
         read(con, '(:8X,9A8)')(hprwbc(jh, jw), jw = 1, nwb)
     enddo
     read(con, '(//(8X,A8,2I8))')(snpc(jw), nsnp(jw), nisnp(jw), jw = 1, nwb)
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(snpd(j, jw), j = 1, nsnp(jw))
     enddo
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(snpf(j, jw), j = 1, nsnp(jw))
     enddo
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9I8)')(isnp(i, jw), i = 1, nisnp(jw))
     enddo
     read(con, '(//(8X,A8,I8))')(scrc(jw), nscr(jw), jw = 1, nwb)
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(scrd(j, jw), j = 1, nscr(jw))
     enddo
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(scrf(j, jw), j = 1, nscr(jw))
     enddo
     read(con, '(//(8X,A8,2I8))')(prfc(jw), nprf(jw), niprf(jw), jw = 1, nwb)
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(prfd(j, jw), j = 1, nprf(jw))
     enddo
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(prff(j, jw), j = 1, nprf(jw))
     enddo
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9I8)')(iprf(j, jw), j = 1, niprf(jw))
     enddo
     read(con, '(//(8X,A8,2I8))')(sprc(jw), nspr(jw), nispr(jw), jw = 1, nwb)
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(sprd(j, jw), j = 1, nspr(jw))
     enddo
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(sprf(j, jw), j = 1, nspr(jw))
     enddo
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9I8)')(ispr(j, jw), j = 1, nispr(jw))
     enddo
     read(con, '(//(8X,A8,I8))')(vplc(jw), nvpl(jw), jw = 1, nwb)
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(vpld(j, jw), j = 1, nvpl(jw))
     enddo
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(vplf(j, jw), j = 1, nvpl(jw))
     enddo
     read(con, '(//(8X,A8,I8,A8))')(cplc(jw), ncpl(jw), tecplot(jw), jw = 1,   &
                                 & nwb)
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(cpld(j, jw), j = 1, ncpl(jw))
     enddo
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(cplf(j, jw), j = 1, ncpl(jw))
     enddo
     read(con, '(//(8X,A8,I8))')(flxc(jw), nflx(jw), jw = 1, nwb)
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(flxd(j, jw), j = 1, nflx(jw))
     enddo
     read(con, '(/)')
     do jw = 1, nwb
         read(con, '(:8X,9F8.0)')(flxf(j, jw), j = 1, nflx(jw))
     enddo
     read(con, '(//8X,A8,2I8)')tsrc, ntsr, niktsr
     allocate(itsr(MAX(1, niktsr)), etsr(MAX(1, niktsr)))
     read(con, '(//(:8X,9F8.0))')(tsrd(j), j = 1, ntsr)
     read(con, '(//(:8X,9F8.0))')(tsrf(j), j = 1, ntsr)
     read(con, '(//(:8X,9I8))')(itsr(j), j = 1, niktsr)
     read(con, '(//(:8X,9F8.0))')(etsr(j), j = 1, niktsr)
     read(con, '(//8X,A8,2I8)')wdoc, nwdo, niwdo
     allocate(iwdo(MAX(1, niwdo)))
     read(con, '(//(:8X,9F8.0))')(wdod(j), j = 1, nwdo)
     read(con, '(//(:8X,9F8.0))')(wdof(j), j = 1, nwdo)
     read(con, '(//(8X,9I8))')(iwdo(j), j = 1, niwdo)
     read(con, '(//8X,A8,I8,A8)')rsoc, nrso, rsic
     rsod = 0.0                                                         ! SW 9/27/13 INITIALIZE SINCE ALLOCATED AS NOD BUT ONLY NRSO USED
     read(con, '(//(:8X,9F8.0))')(rsod(j), j = 1, nrso)
     read(con, '(//(:8X,9F8.0))')(rsof(j), j = 1, nrso)
 
!    Constituent control cards
 
     read(con, '(//8X,2A8,I8)')ccc, limc, cuf
     read(con, '(//(2A8))')(cname2(jc), cac(jc), jc = 1, nct)
     read(con, '(/)')
 
  !IF(.NOT.RESTART_PUSHED)THEN
     nndc = ndc - 1    ! SW 4/14/2017
    !ELSE
    !    NNDC=NDC
    !ENDIF
 
     do jd = 1, nndc
                  !NDC-1    ! SW 10/20/15 ADDED AN EXTRA DERIVED VARIABLE INTERNALLY
         if(nwb<10)read(con, '(A8,(:9A8))')cdname2(jd),                        &
                      & (cdwbc(jd, jw), jw = 1, nwb)
         if(nwb>=10)read(con, '(A8,9A8,/(:8X,9A8))')cdname2(jd),               &
                       & (cdwbc(jd, jw), jw = 1, nwb)                                                     !cb 9/13/12  sw 2/18/13  Foramt 6/16/13 8/13/13
     enddo
 
     cdname2(ndc) = ' TDG(%)'       ! SW 10/20/15
     cdwbc(ndc, :) = '      ON'
 
     read(con, '(/)')
     kfname2 = '     '
                     ! SW 9/27/13 INITIALIZE ENTIRE ARRAY
     kfwbc = '     ' ! SW 9/27/13 INITIALIZE ENTIRE ARRAY
!    DO JF=1,NFL
     do jf = 1, 73
               ! Fix this later
         if(nwb<10)read(con, '(A8,(:9A8))')kfname2(jf),                        &
                      & (kfwbc(jf, jw), jw = 1, nwb)
         if(nwb>=10)read(con, '(A8,9A8,/(:8X,9A8))')kfname2(jf),               &
                       & (kfwbc(jf, jw), jw = 1, nwb)                                                    !cb 9/13/12  sw2/18/13  Foramt 6/16/13 8/13/13
         kfname2(jf) = kfname2(jf)(1:8) // '(kg/d)'
     enddo
     kfname2(121) = 'CO2GASX(kg/d)'
!    CEMA and Amaila start
     kfname2(122) = 'DOH2S(kg/d)'
!    CEMA and Amaila start
     kfname2(123) = 'DOCH4(kg/d)'
!    CEMA and Amaila start
     kfname2(124) = 'H2SGASX(kg/d)'
!    CEMA and Amaila start
     kfname2(125) = 'CH4GASX(kg/d)'
     kfname2(126) = 'H2SDK(kg/d)'
     kfname2(127) = 'CH4DK(kg/d)'
     kfname2(128) = 'SD_C_IN(kg/d)'
     kfname2(129) = 'SD_N_IN(kg/d)'
     kfname2(130) = 'SD_P_IN(kg/d)'
     kfname2(131) = 'DOSEDIA(kg/d)'
     kfname2(132) = 'Fe2D(kg/d)'
     kfname2(133) = 'DOFe2(kg/d)'
     kfname2(134) = 'SDINFeOOH(kg/d)'
     kfname2(135) = 'SDINMnO2(kg/d)'
     kfname2(136) = 'Mn2d(kg/d)'
     kfname2(137) = 'DOMn2(kg/d)'
     kfname2(138) = 'SEDD1(kg/d)'
     kfname2(139) = 'SEDD2(kg/d)'
!    CEMA end
     read(con, '(/)')
     do jc = 1, nct
         read(con, '(:8X,9F8.0)')(c2i(jc, jw), jw = 1, nwb)
     enddo
     read(con, '(/)')
     do jc = 1, nct
         read(con, '(:8X,9A8)')(cprwbc(jc, jw), jw = 1, nwb)
     enddo
     read(con, '(/)')
     do jc = 1, nct
         read(con, '(:8X,9A8)')(cinbrc(jc, jb), jb = 1, nbr)
     enddo
     read(con, '(/)')
     do jc = 1, nct
         read(con, '(:8X,9A8)')(ctrtrc(jc, jt), jt = 1, ntr)
     enddo
     read(con, '(/)')
     do jc = 1, nct
         read(con, '(:8X,9A8)')(cdtbrc(jc, jb), jb = 1, nbr)
     enddo
     read(con, '(/)')
     do jc = 1, nct
         read(con, '(:8X,9A8)')(cprbrc(jc, jb), jb = 1, nbr)
     enddo
 
!    Kinetics coefficients
 
     read(con, '(//(8X,4F8.0,2A8))')(exh2o(jw), exss(jw), exom(jw), beta(jw),  &
                                  & exc(jw), exic(jw), jw = 1, nwb)
     read(con, '(//(8X,9F8.0))')(exa(ja), ja = 1, nal)
     read(con, '(//(8X,9F8.0))')(exz(jz), jz = 1, nzpt)
     read(con, '(//(8X,9F8.0))')(exm(jm), jm = 1, nmct)
!    READ (CON,'(//(8X,4F8.0))')         (CGQ10(JG),  CG0DK(JG),  CG1DK(JG), 
!    CGS(JG),                            JG=1,NGC)
     read(con, '(//(8X,7F8.0))')(cgq10(jg), cg0dk(jg), cg1dk(jg), cgs(jg),     &
                              & cgldk(jg), cgklf(jg), cgcs(jg), jg = 1, ngc)                                                          !LCJ 2/26/15
     read(con, '(//(8X,F8.0,A,F8.0))')(sss(js), sedrc(js), taucr(js), js = 1,  &
                                    & nss)                                                                                   ! READ (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
     read(con, '(//(8X,9F8.0))')(ag(ja), ar(ja), ae(ja), am(ja), as(ja),       &
                              & ahsp(ja), ahsn(ja), ahssi(ja), asat(ja),       &
                              & ja = 1, nal)
     read(con, '(//(8X,8F8.0))')(at1(ja), at2(ja), at3(ja), at4(ja), ak1(ja),  &
                              & ak2(ja), ak3(ja), ak4(ja), ja = 1, nal)
     read(con, '(//(8X,6F8.0,I8,F8.0))')(ap(ja), an(ja), ac(ja), asi(ja), achla&
                                      & (ja), apom(ja), aneqn(ja), anpr(ja),   &
                                      & ja = 1, nal)
     read(con, '(//(8X,9A8))')(epic(jw, 1), jw = 1, nwb)
     do je = 2, nept
         read(con, '(8X,9A8)')(epic(jw, je), jw = 1, nwb)
     enddo
     read(con, '(//(8X,9A8))')(epiprc(jw, 1), jw = 1, nwb)
     do je = 2, nept
         read(con, '(8X,9A8)')(epiprc(jw, je), jw = 1, nwb)
     enddo
     read(con, '(//(8X,9F8.0))')(epici(jw, 1), jw = 1, nwb)
     do je = 2, nept
         read(con, '(8X,9F8.0)')(epici(jw, je), jw = 1, nwb)
     enddo
     read(con, '(//(8X,8F8.0))')(eg(je), er(je), ee(je), em(je), eb(je),       &
                              & ehsp(je), ehsn(je), ehssi(je), je = 1, nept)                                                     !JE=1,NEP)  SW 9/27/13
     read(con, '(//(8X,2F8.0,I8,F8.0))')(esat(je), ehs(je), eneqn(je), enpr(je)&
                                      & , je = 1, nept)                                                                          !JE=1,NEP)  SW 9/27/13
     read(con, '(//(8X,8F8.0))')(et1(je), et2(je), et3(je), et4(je), ek1(je),  &
                              & ek2(je), ek3(je), ek4(je), je = 1, nept)                                                         !JE=1,NEP)  SW 9/27/13
     read(con, '(//(8X,6F8.0))')(ep(je), en(je), ec(je), esi(je), echla(je),   &
                              & epom(je), je = 1, nept)                                                                          !JE=1,NEP)  SW 9/27/13
     read(con, '(//(8X,7F8.0))')(zg(jz), zr(jz), zm(jz), zeff(jz), prefp(jz),  &
                              & zoomin(jz), zs2p(jz), jz = 1, nzpt)
 
     read(con, '(//(8X,8F8.0))')(prefa(ja, 1), ja = 1, nal)                                                                       ! MM 7/13/06
     do jz = 2, nzpt
         read(con, '((8X,8F8.0))')(prefa(ja, jz), ja = 1, nal)
     enddo
     read(con, '(//(8X,8F8.0))')(prefz(jjz, 1), jjz = 1, nzpt)
     do jz = 2, nzpt
         read(con, '((8X,8F8.0))')(prefz(jjz, jz), jjz = 1, nzpt)                                                                  ! MM 7/13/06
     enddo
     read(con, '(//(8X,8F8.0))')(zt1(jz), zt2(jz), zt3(jz), zt4(jz), zk1(jz),  &
                              & zk2(jz), zk3(jz), zk4(jz), jz = 1, nzpt)
     read(con, '(//(8X,3F8.0))')(zp(jz), zn(jz), zc(jz), jz = 1, nzpt)
     read(con, '(//(8X,9A8))')(macwbc(jw, 1), jw = 1, nwb)
     do jm = 2, nmct
         read(con, '(8X,9A8)')(macwbc(jw, jm), jw = 1, nwb)
     enddo
     read(con, '(//(8X,9A8))')(mprwbc(jw, 1), jw = 1, nwb)
     do jm = 2, nmct
         read(con, '(8X,9A8)')(mprwbc(jw, jm), jw = 1, nwb)
     enddo
     read(con, '(//(8X,9F8.0))')(macwbci(jw, 1), jw = 1, nwb)
     do jm = 2, nmct
         read(con, '(8X,9F8.0)')(macwbci(jw, jm), jw = 1, nwb)
     enddo
     read(con, '(//(8X,9F8.0))')(mg(jm), mr(jm), mm(jm), msat(jm), mhsp(jm),   &
                              & mhsn(jm), mhsc(jm), mpom(jm), lrpmac(jm),      &
                              & jm = 1, nmct)
     read(con, '(//(8X,2F8.0))')(psed(jm), nsed(jm), jm = 1, nmct)
     read(con, '(//(8X,2F8.0))')(mbmp(jm), mmax(jm), jm = 1, nmct)
     read(con, '(//(8X,4F8.0))')(cddrag(jm), dwv(jm), dwsa(jm), anorm(jm),     &
                              & jm = 1, nmct)                                                                             !CB 6/29/06
     read(con, '(//(8X,8F8.0))')(mt1(jm), mt2(jm), mt3(jm), mt4(jm), mk1(jm),  &
                              & mk2(jm), mk3(jm), mk4(jm), jm = 1, nmct)
     read(con, '(//(8X,3F8.0))')(mp(jm), mn(jm), mc(jm), jm = 1, nmct)
     read(con, '(//(8X,3F8.0))')(ldomdk(jw), rdomdk(jw), lrddk(jw), jw = 1,    &
                              & nwb)
     read(con, '(//(8X,4F8.0))')(lpomdk(jw), rpomdk(jw), lrpdk(jw), poms(jw),  &
                              & jw = 1, nwb)
     read(con, '(//(8X,4F8.0))')(orgp(jw), orgn(jw), orgc(jw), orgsi(jw),      &
                              & jw = 1, nwb)
     read(con, '(//(8X,4F8.0))')(omt1(jw), omt2(jw), omk1(jw), omk2(jw),       &
                              & jw = 1, nwb)
     read(con, '(//(8X,4F8.0))')(kbod(jb), tbod(jb), rbod(jb), cbods(jb),      &
                              & jb = 1, nbod)
     read(con, '(//(8X,3F8.0))')(bodp(jb), bodn(jb), bodc(jb), jb = 1, nbod)
     read(con, '(//(8X,2F8.0))')(po4r(jw), partp(jw), jw = 1, nwb)
     read(con, '(//(8X,2F8.0))')(nh4r(jw), nh4dk(jw), jw = 1, nwb)
     read(con, '(//(8X,4F8.0))')(nh4t1(jw), nh4t2(jw), nh4k1(jw), nh4k2(jw),   &
                              & jw = 1, nwb)
     read(con, '(//(8X,3F8.0))')(no3dk(jw), no3s(jw), fno3sed(jw), jw = 1, nwb)
     read(con, '(//(8X,4F8.0))')(no3t1(jw), no3t2(jw), no3k1(jw), no3k2(jw),   &
                              & jw = 1, nwb)
     read(con, '(//(8X,4F8.0))')(dsir(jw), psis(jw), psidk(jw), partsi(jw),    &
                              & jw = 1, nwb)
     read(con, '(//(8X,2F8.0))')(fer(jw), fes(jw), jw = 1, nwb)
     read(con, '(//(8X,F8.0))')(co2r(jw), jw = 1, nwb)
     read(con, '(//(8X,2F8.0))')(o2nh4(jw), o2om(jw), jw = 1, nwb)
     read(con, '(//(8X,2F8.0))')(o2ar(ja), o2ag(ja), ja = 1, nal)
     read(con, '(//(8X,2F8.0))')(o2er(je), o2eg(je), je = 1, nept)
     read(con, '(//(8X,F8.0))')(o2zr(jz), jz = 1, nzpt)
     read(con, '(//(8X,2F8.0))')(o2mr(jm), o2mg(jm), jm = 1, nmct)
     read(con, '(//(8X,F8.0))')kdo
     if(kdo==0.0)kdo = 0.01                                  ! SW 10/24/15 ERROR TRAPPING
     read(con, '(//(8X,2A8,6F8.0,A8))')(sedcc(jw), sedprc(jw), sedci(jw), sdk( &
                                     & jw), seds(jw), fsod(jw), fsed(jw),      &
                                     & sedb(jw), dynsedk(jw), jw = 1, nwb)                                                                                   ! cb 11/28/06
     read(con, '(//(8X,4F8.0))')(sodt1(jw), sodt2(jw), sodk1(jw), sodk2(jw),   &
                              & jw = 1, nwb)
     read(con, '(//(8X,9F8.0))')(sod(i), i = 1, imx)
     read(con, '(//(8X,A8,I8,4F8.2))')(reaerc(jw), neqn(jw), rcoef1(jw), rcoef2&
                                    & (jw), rcoef3(jw), rcoef4(jw), jw = 1,    &
                                    & nwb)
 
!    Input filenames
 
     read(con, '(//(8X,A72))')rsifn
     read(con, '(//(8X,A72))')qwdfn
     read(con, '(//(8X,A72))')qgtfn
     read(con, '(//(8X,A72))')wscfn
     read(con, '(//(8X,A72))')shdfn
     read(con, '(//(8X,A72))')(bthfn(jw), jw = 1, nwb)
     read(con, '(//(8X,A72))')(metfn(jw), jw = 1, nwb)
     read(con, '(//(8X,A72))')(extfn(jw), jw = 1, nwb)
     read(con, '(//(8X,A72))')(vprfn(jw), jw = 1, nwb)
     read(con, '(//(8X,A72))')(lprfn(jw), jw = 1, nwb)
     read(con, '(//(8X,A72))')(qinfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(tinfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(cinfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(qotfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(qtrfn(jt), jt = 1, ntr)
     read(con, '(//(8X,A72))')(ttrfn(jt), jt = 1, ntr)
     read(con, '(//(8X,A72))')(ctrfn(jt), jt = 1, ntr)
     read(con, '(//(8X,A72))')(qdtfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(tdtfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(cdtfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(prefn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(tprfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(cprfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(euhfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(tuhfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(cuhfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(edhfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(tdhfn(jb), jb = 1, nbr)
     read(con, '(//(8X,A72))')(cdhfn(jb), jb = 1, nbr)
 
!    Output filenames
 
     read(con, '(//(8X,A72))')(snpfn(jw), jw = 1, nwb)
     read(con, '(//(8X,A72))')(prffn(jw), jw = 1, nwb)
     read(con, '(//(8X,A72))')(vplfn(jw), jw = 1, nwb)
     read(con, '(//(8X,A72))')(cplfn(jw), jw = 1, nwb)
     read(con, '(//(8X,A72))')(sprfn(jw), jw = 1, nwb)
     read(con, '(//(8X,A72))')(flxfn(jw), jw = 1, nwb)
     read(con, '(//(8X,A72))')tsrfn1
     read(con, '(//(8X,A72))')wdofn
     close(con)
 
!    Bathymetry file
 
     do jw = 1, nwb
         open(bth(jw), file = bthfn(jw), status = 'OLD')
         read(bth(jw), '(a1)')char1                                 ! New Bathymetry format option SW 6/22/09
         if(char1=='$')then
             read(bth(jw), *)
             read(bth(jw), *)aid, (dlx(i), i = us(bs(jw)) - 1, ds(be(jw)) + 1)
             read(bth(jw), *)aid, (elws(i), i = us(bs(jw)) - 1, ds(be(jw)) + 1)
             read(bth(jw), *)aid, (phi0(i), i = us(bs(jw)) - 1, ds(be(jw)) + 1)
             read(bth(jw), *)aid, (fric(i), i = us(bs(jw)) - 1, ds(be(jw)) + 1)
             read(bth(jw), *)
             do k = 1, kmx
                 read(bth(jw), *)h(k, jw),                                     &
                               & (b(k, i), i = us(bs(jw)) - 1, ds(be(jw)) + 1)
             enddo
             do i = us(bs(jw)) - 1, ds(be(jw)) + 1
                 h2(:, i) = h(:, jw)
             enddo
         else
             read(bth(jw), '(//(10F8.0))')                                     &
                & (dlx(i), i = us(bs(jw)) - 1, ds(be(jw)) + 1)
             read(bth(jw), '(//(10F8.0))')                                     &
                & (elws(i), i = us(bs(jw)) - 1, ds(be(jw)) + 1)
             read(bth(jw), '(//(10F8.0))')                                     &
                & (phi0(i), i = us(bs(jw)) - 1, ds(be(jw)) + 1)
             read(bth(jw), '(//(10F8.0))')                                     &
                & (fric(i), i = us(bs(jw)) - 1, ds(be(jw)) + 1)
             read(bth(jw), '(//(10F8.0))')(h(k, jw), k = 1, kmx)
             do i = us(bs(jw)) - 1, ds(be(jw)) + 1
                 read(bth(jw), '(//(10F8.0))')(b(k, i), k = 1, kmx)
                 h2(:, i) = h(:, jw)
             enddo
         endif
         close(bth(jw))
     enddo
     h1 = h2
     bi = b
 
     allocate(bsave(kmx, imx))
     bsave = 0.0
     bsave = b
 
!    Amaila start - reading additional sediment compartments coefficients
     sedcomp_exist = .FALSE.
     inquire(file = 'w2_amaila.npt', exist = sedcomp_exist)
                                                       ! SW 4/30/15
     if(sedcomp_exist)then
         open(nunit, file = 'w2_amaila.npt', status = 'old')
   !READ (NUNIT,'(//(8X,2A8,6F8.0,A8))')     (SEDCC1(JW),   SEDPRC1(JW), SEDCI1(JW),  SDK1(JW),   JW=1,NWB)
         read(nunit, '(//(8X,2A8,3F8.0))')(sedcc1(jw), sedprc1(jw), sedci1(jw),&
            & sdk1(jw), fsedc1(jw), jw = 1, nwb)                                                                      ! cb 6/7/17
   !READ (NUNIT,'(//(8X,2A8,6F8.0,A8))')     (SEDCC2(JW),   SEDPRC2(JW), SEDCI2(JW),  SDK2(JW),   JW=1,NWB)
         read(nunit, '(//(8X,2A8,3F8.0))')(sedcc2(jw), sedprc2(jw), sedci2(jw),&
            & sdk2(jw), fsedc2(jw), jw = 1, nwb)
         read(nunit, '(//(8X,3F8.0))')(pbiom(jw), nbiom(jw), cbiom(jw), jw = 1,&
                                    & nwb)                                            ! cb 6/7/17
         close(nunit)
     endif
 
 
!    Output file unit numbers
 
     allocate(tsr(niktsr))
     allocate(wdo(niwdo, 4), wdo2(nwd + nst + ngt + nsp + npu + npi, 4))
     do j = 1, 7
         do jw = 1, nwb
             opt(jw, j) = nunit
             nunit = nunit + 1
         enddo
     enddo
     do j = 1, niktsr
         tsr(j) = nunit
         nunit = nunit + 1
     enddo
     do jw = 1, niwdo
         wdo(jw, 1) = nunit
         nunit = nunit + 1
         wdo(jw, 2) = nunit
         nunit = nunit + 1
         wdo(jw, 3) = nunit
         nunit = nunit + 1
         wdo(jw, 4) = nunit
         nunit = nunit + 1
     enddo
 
!    BIOENERGETICS bioexp mlm output filenumber assigment
     if(fishbio)then
         do j = 1, nibio
             bioexpfn(j) = nunit
             nunit = nunit + 1
         enddo
         do j = 1, nibio
             weightnum(j) = nunit
             nunit = nunit + 1
         enddo
     endif
!    Variable names, formats, multipliers, and Compaq Visual FORTRAN array
 
!    viewer controls
     open(grf, file = 'graph.npt', status = 'OLD')
     read(grf, '(///(A43,1X,A9,3F8.0,A8))')                                    &
        & (hname(j), fmth(j), hmult(j), hymin(j), hymax(j), hpltc(j), j = 1,   &
        & nhy)
     read(grf, '(// (A43,1X,A9,3F8.0,A8))')                                    &
        & (cname(j), fmtc(j), cmult(j), cmin(j), cmax(j), cpltc(j), j = 1, nct)
     read(grf, '(// (A43,1X,A9,3F8.0,A8))')                                    &
        & (cdname(j), fmtcd(j), cdmult(j), cdmin(j), cdmax(j), cdpltc(j),      &
        & j = 1, nndc)                                                                                              ! SW 10/20/15 INTERNAL TDG
     close(grf)
 
     cdname(ndc) = 'TDG(%)' ! SW 10/17/15
     fmtcd(ndc) = ' (F10.3)'
     cdmult(ndc) = 1.0
 
 
     do jc = 1, nct
         l3 = 1
         l1 = SCAN(cname(jc), ',') + 2
         if(l1==2)l1 = 43       ! SW 12/3/2012   Implies no comma found
         l2 = SCAN(cname(jc)(l1:43), '  ') + l1
         if(l2>43)l2 = 43       ! SW 12/3/2012
         cunit(jc) = cname(jc)(l1:l2)
         cname1(jc) = cname(jc)(1:l1 - 3)
         cname3(jc) = cname1(jc)
         do while (l3<l1 - 3)
             if(cname(jc)(l3:l3)==' ')cname3(jc)(l3:l3) = '_'
             l3 = l3 + 1
         enddo
         cunit1(jc) = cunit(jc)(1:1)
         cunit2(jc) = cunit(jc)
         if(cunit(jc)(1:2)=='mg')then
             cunit1(jc) = 'g'
             cunit2(jc) = 'g/m^3'
         endif
         if(cunit(jc)(1:2)/='g/' .AND. cunit(jc)(1:2)/='mg')cunit1(jc) = '  '
     enddo
     do jc = 1, ndc
         l1 = 1
         l2 = MAX(4, SCAN(cdname(jc), ',') - 1)
         cdname3(jc) = cdname(jc)(1:l2)
         do while (l1<l2)
             if(cdname(jc)(l1:l1)==' ')cdname3(jc)(l1:l1) = '_'
             l1 = l1 + 1
         enddo
     enddo
     fmth(1:nhy) = ADJUSTL(fmth(1:nhy))
 
!    Initialize logical control variables
 
     vert_profile = .FALSE.
     long_profile = .FALSE.
     constituents = ccc=='      ON'
     do jw = 1, nwb
         iso_temp(jw) = t2i(jw)>=0
         vert_temp(jw) = t2i(jw)== - 1
         long_temp(jw) = t2i(jw)< - 1
         if(constituents)then                ! CB 12/04/08
             iso_sediment(jw) = sedci(jw)>=0 .AND. sedcc(jw)=='      ON'
             vert_sediment(jw) = sedci(jw)== - 1.0 .AND. sedcc(jw)=='      ON'
             long_sediment(jw) = sedci(jw)< - 1.0 .AND. sedcc(jw)=='      ON'
!            Amaila Start
             iso_sediment1(jw) = sedci1(jw)>=0 .AND. sedcc1(jw)=='      ON'
             vert_sediment1(jw) = sedci1(jw)== - 1.0 .AND. sedcc1(jw)          &
                                 &=='      ON'
             long_sediment1(jw) = sedci1(jw)< - 1.0 .AND. sedcc1(jw)           &
                                 &=='      ON'
             iso_sediment2(jw) = sedci2(jw)>=0 .AND. sedcc2(jw)=='      ON'
             vert_sediment2(jw) = sedci2(jw)== - 1.0 .AND. sedcc2(jw)          &
                                 &=='      ON'
             long_sediment2(jw) = sedci2(jw)< - 1.0 .AND. sedcc2(jw)           &
                                 &=='      ON'
!            Amaila End
             iso_epiphyton(jw, :) = epici(jw, :)>=0 .AND. epic(jw, :)          &
                                   &=='      ON'
             vert_epiphyton(jw, :) = epici(jw, :)== - 1.0 .AND. epic(jw, :)    &
                                    &=='      ON'
             long_epiphyton(jw, :) = epici(jw, :)< - 1.0 .AND. epic(jw, :)     &
                                    &=='      ON'
             iso_macrophyte(jw, :) = macwbci(jw, :)>=0 .AND. macwbc(jw, :)     &
                                    &=='      ON'                                    ! cb 8/21/15
             vert_macrophyte(jw, :) = macwbci(jw, :)== - 1.0 .AND.             &
                                    & macwbc(jw, :)=='      ON'                      ! cb 8/21/15
             long_macrophyte(jw, :) = macwbci(jw, :)< - 1.0 .AND. macwbc(jw, :)&
                                    & =='      ON'                                   ! cb 8/21/15
             do jc = 1, nct
                 iso_conc(jc, jw) = c2i(jc, jw)>=0.0
                 vert_conc(jc, jw) = c2i(jc, jw)== - 1.0 .AND. cac(jc)         &
                                    &=='      ON'
                 long_conc(jc, jw) = c2i(jc, jw)< - 1.0 .AND. cac(jc)          &
                                    &=='      ON'
                 if(vert_conc(jc, jw))vert_profile(jw) = .TRUE.
                 if(long_conc(jc, jw))long_profile(jw) = .TRUE.
             enddo
             if(vert_sediment(jw))vert_profile(jw) = .TRUE.
             if(vert_sediment1(jw))vert_profile(jw) = .TRUE.   ! amaila
             if(vert_sediment2(jw))vert_profile(jw) = .TRUE.   ! amaila
             if(long_sediment(jw))long_profile(jw) = .TRUE.
             if(long_sediment1(jw))long_profile(jw) = .TRUE.   ! amaila
             if(long_sediment2(jw))long_profile(jw) = .TRUE.   ! amaila
             if(ANY(vert_epiphyton(jw, :)))vert_profile(jw) = .TRUE.
             if(ANY(long_epiphyton(jw, :)))long_profile(jw) = .TRUE.
             if(ANY(vert_macrophyte(jw, :)))vert_profile(jw) = .TRUE.
                                                                ! cb 8/21/15
             if(ANY(long_macrophyte(jw, :)))long_profile(jw) = .TRUE.
                                                                ! cb 8/21/15
         endif                      ! cb 12/04/08
         if(vert_temp(jw))vert_profile(jw) = .TRUE.
         if(long_temp(jw))long_profile(jw) = .TRUE.
         do m = 1, nmc
      !MACROPHYTE_CALC(JW,M) = CONSTITUENTS.AND.MACWBC(JW,M).EQ.' ON'
      !PRINT_MACROPHYTE(JW,M) = MACROPHYTE_CALC(JW,M).AND.MPRWBC(JW,M).EQ.' ON'
             macrophyte_calc(jw, m) = constituents .AND. macwbc(jw, m)         &
                                     &=='      ON'                                     ! cb 8/24/15
             print_macrophyte(jw, m) = macrophyte_calc(jw, m) .AND.            &
                                     & mprwbc(jw, m)=='      ON'                      ! cb 8/24/15
             if(macrophyte_calc(jw, m))macrophyte_on = .TRUE.
         enddo
     enddo
 
!    Open error and warning files
 
     open(w2err, file = 'w2.err', status = 'UNKNOWN')
 
!    Open error and warning files
 
     open(wrn, file = 'w2.wrn', status = 'UNKNOWN')
 
!    Initialize variables for enhanced pH buffering ! entire section ! SR
!    01/01/12
     phbuff_exist = .FALSE.
     inquire(file = 'pH_buffering.npt', exist = phbuff_exist)
     if(constituents .AND. phbuff_exist)then
         open(nunit, file = 'ph_buffering.npt', status = 'OLD')
         read(nunit, '(///8X,2A8)')phbufc, ncalkc
         read(nunit, '(//8X,3A8)')nh4bufc, po4bufc, ombufc
         read(nunit, '(//8X,A8,I8,A8)')omtype, nagi, pombufc
         allocate(sdeni(nagi), pki(nagi), pksd(nagi))
         read(nunit, '(//(:8X,9F8.0))')(sdeni(j), j = 1, nagi)
         read(nunit, '(//(:8X,9F8.0))')(pki(j), j = 1, nagi)
         read(nunit, '(//(:8X,9F8.0))')(pksd(j), j = 1, nagi)
         close(nunit)
         ph_buffering = phbufc=='      ON'
         noncon_alkalinity = ncalkc=='      ON'
         ammonia_buffering = nh4bufc=='      ON'
         phosphate_buffering = po4bufc=='      ON'
         om_buffering = ombufc=='      ON'
         pom_buffering = pombufc=='      ON' .AND. om_buffering
         if(om_buffering)then
             sdeni = ABS(sdeni)
             if(omtype=='    DIST')then
                 if(ANY(pksd<=0))then
                     warning_open = .TRUE.
                     write(wrn, '(A)')                                         &
     &'WARNING -- PKSD inputs in the ph_buffering.npt file must be greater than&
     & zero.'
                     write(wrn, '(A/)')                                        &
      &'Please fix your inputs. For now, PKSD values of zero will be set to 1.'
                 endif
                 do ja = 1, nagi
                     if(pksd(ja)<=0)pksd(ja) = 1.0
                 enddo
                 nag = 27
                 allocate(sden(nag), pk(nag), fract(nag))
                 sden = 0.0
                 do j = 1, nag
                     pk(j) = 0.5*j
                 enddo
                 do ja = 1, nagi
                     sum = 0.0
                     do j = 1, nag
                         fract(j) = EXP( - 0.5*(((pk(j)-pki(ja))/pksd(ja))**2))
                         sum = sum + fract(j)
                     enddo
                     do j = 1, nag
                         sden(j) = sden(j) + sdeni(ja)*fract(j)/sum
                     enddo
                 enddo
             else
                 allocate(sden(nagi), pk(nagi))
                 nag = nagi
                 sden = sdeni
                 pk = pki
                 omtype = ' MONO'
             endif
         endif
         open(nunit, file = 'ph_buffering.opt', status = 'UNKNOWN')
         write(nunit, '(A/)')'Enhanced pH buffering output file'
         write(nunit, '(2A)')'Ammonia buffering: ', ADJUSTL(TRIM(nh4bufc))
         write(nunit, '(2A)')'Phosphate buffering: ', ADJUSTL(TRIM(po4bufc))
         write(nunit, '(2A)')'OM buffering: ', ADJUSTL(TRIM(ombufc))
         if(om_buffering)then
             write(nunit, '(2A)')'POM buffering: ', ADJUSTL(TRIM(pombufc))
             write(nunit, '(2A)')'OM buffer type: ', ADJUSTL(TRIM(omtype))
             write(nunit, '(/A/A)')'Inputs:', 'Group Site density pKa std.dev.'
             do ja = 1, nagi
                 if(omtype==' DIST')then
                     write(nunit, '(1X,I3,5X,F8.4,4X,F6.3,3X,F6.3)')ja,        &
                         & sdeni(ja), pki(ja), pksd(ja)
                 else
                     write(nunit, '(1X,I3,5X,F8.4,4X,F6.3,3X,A)')ja, sdeni(ja),&
                         & pki(ja), ' N/A'
                 endif
             enddo
             write(nunit, '(/A/A)')'Modeled:', 'Group Site density pKa'
             do ja = 1, nag
                 write(nunit, '(1X,I3,5X,F8.4,4X,F6.3)')ja, sden(ja), pk(ja)
             enddo
         endif
         close(nunit)
     else
         ammonia_buffering = .FALSE.
         phosphate_buffering = .FALSE.
         om_buffering = .FALSE.
         pom_buffering = .FALSE.
     endif
     ngctdg = 0
     ngn2 = 0
     ndc = ndc - 1
                ! SW 10/20/15
     do jg = 1, ngc
                ! SW 10/16/2015
         if(cgcs(jg)== - 1.0)then
             ngctdg = jg
                       ! JG COUNTER
             ndc = ndc + 1
                       ! THERE CAN BE ONLY 1 TDG AMONG THE CONSITUENTS
             ngn2 = ngcs + jg - 1
                                ! GLOBAL COUNTER
             exit
         endif
     enddo
     if(ngctdg==0)cdwbc(ndc + 1, :) = '     OFF'  ! SR 7/15/17
     end subroutine INPUT
