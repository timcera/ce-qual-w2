!*==time_varying_data.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                  S U B R O U T I N E   T I M E   V A R Y I N G   D A T A                                      **
!***********************************************************************************************************************************
 
     subroutine TIME_VARYING_DATA
     use f77kinds
     use GLOBAL
     use SURFHE
     use SCREENC
     use TVDC
     use LOGICC
     use SELWC
     use STRUCTURES
     use NAMESC
     use KINETIC, ONLY:EXH2O
     use SHADEC
     use MAIN, ONLY:pumps, segnum
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real :: Nxtvd
     intent (inout) Nxtvd
!
! Local variables
!
     real, allocatable, dimension(:), save :: bgtnx, bpnx, cloudnx, cloudo,    &
       & eoffpu2, eonpu2, epu2, extnx, exto, nxcdh1, nxcdh2, nxcdt1, nxcdt2,   &
       & nxcin1, nxcin2, nxcpr1, nxcpr2, nxctr1, nxctr2, nxcuh1, nxcuh2,       &
       & nxdyns, nxedh1, nxedh2, nxeuh1, nxeuh2, nxext1, nxext2, nxmet1,       &
       & nxmet2, nxpr1, nxpr2, nxpump, nxqdt1, nxqdt2, nxqin1, nxqin2, nxqot1, &
       & nxqot2, nxqtr1, nxqtr2, nxtdh1, nxtdh2, nxtdt1, nxtdt2, nxtin1,       &
       & nxtin2, nxtpr1, nxtpr2, nxttr1, nxttr2, nxtuh1, nxtuh2, phinx, phio,  &
       & qpu2, sronx, sroo, tairnx, tairo, tdewnx, tdewo, windnx, windo, wscnx,&
       & xx
     integer, allocatable, dimension(:), save :: cdhf, cuhf, dhc, dhe, dht,    &
          & dtc, dtq, dtt, edhf, euhf, ext, inc, inft, inq, jjs, met, odyns,   &
          & otq, prc, pre, prt, pumpd, tdhf, trc, trq, trt, tuhf, uhc, uhe, uht
     real, allocatable, dimension(:, :, :), save :: cdhnx, cdho, cuhnx, cuho
     real(R8KIND), allocatable, dimension(:, :), save :: cdtrnx, cdtro, cinnx,     &
           & cino, cprnx, ctrnx, ctro, qoutnx, qouto, qstrnx, qstro, tdhnx,    &
           & tdho, tuhnx, tuho
     real(R8KIND), save :: cratio, hratio, qratio, ratio, tratio
     logical, allocatable, dimension(:), save :: dtcf, dtqf, dtrib_const, dttf,&
          & dynef, dynpumpf, extf, incf, inflow_const, inqf, intf, metf, otqf, &
          & prcf, precip_const, prqf, prtf, trcf, trib_const, trqf, trtf
     real(R8KIND), allocatable, dimension(:), save :: eldhnx, eldho, eluhnx, eluho,&
           & prnx, qdtrnx, qdtro, qinnx, qino, qtrnx, qtro, qwdnx, qwdo,       &
           & tdtrnx, tdtro, tinnx, tino, tprnx, ttrnx, ttro
     logical, save :: gatef, wdqf, wshf
     integer, save :: gtq, iopenpipe, j, jac, jg, js, jt, jwd, k, l, njs, npt, &
                    & piped, shd, wdq, wsh
     character(1), save :: informat
     character(2), save :: informat2
     real, allocatable, dimension(:, :), save :: nxestrt
     real, save :: nxqgt, nxqgt2, nxqpt, nxqwd1, nxqwd2, nxwsc
!
!*** End of declarations rewritten by SPAG
!
 
!    Type declaration
 
                                                               ! SW 9/26/2017
!    Allocation declarations
 
     if(nac<1)then                                                           ! =0 Old format for head BCs, =1 Time series format no vertical variation, =2 csv format vertical variation             SW 2/28/17
         allocate(xx(1))
     else
         allocate(xx(nct))
     endif
 
     allocate(nxqtr1(ntr), nxttr1(ntr), nxctr1(ntr), nxqin1(nbr), nxtin1(nbr), &
            & nxcin1(nbr), nxqdt1(nbr), nxtdt1(nbr), nxcdt1(nbr))
     allocate(nxpr1(nbr), nxtpr1(nbr), nxcpr1(nbr), nxeuh1(nbr), nxtuh1(nbr),  &
            & nxcuh1(nbr), nxedh1(nbr), nxtdh1(nbr), nxcdh1(nbr))
     allocate(nxqot1(nbr), nxmet1(nwb), nxqtr2(ntr), nxttr2(ntr), nxctr2(ntr), &
            & nxqin2(nbr), nxtin2(nbr), nxcin2(nbr), nxqdt2(nbr))
     allocate(nxtdt2(nbr), nxcdt2(nbr), nxpr2(nbr), nxtpr2(nbr), nxcpr2(nbr),  &
            & nxeuh2(nbr), nxtuh2(nbr), nxcuh2(nbr), nxedh2(nbr))
     allocate(nxtdh2(nbr), nxcdh2(nbr), nxqot2(nbr), nxmet2(nwb), dynpumpf(npu)&
            & )
     allocate(wscnx(imx), prcf(nbr), metf(nwb), odyns(nbr), nxdyns(nbr),       &
            & nxestrt(nst, nbr), prtf(nbr), prqf(nbr), extf(nwb))
     allocate(qdtro(nbr), tdtro(nbr), eluho(nbr), eldho(nbr), qwdo(nwd),       &
            & qtro(ntr), ttro(ntr), qino(nbr), tino(nbr))
     allocate(qdtrnx(nbr), tdtrnx(nbr), prnx(nbr), tprnx(nbr), eluhnx(nbr),    &
            & eldhnx(nbr), qwdnx(nwd), qtrnx(ntr), ttrnx(ntr))
     allocate(qinnx(nbr), tinnx(nbr), sroo(nwb), tairo(nwb), tdewo(nwb),       &
            & cloudo(nwb), phio(nwb), windo(nwb), tairnx(nwb))
     allocate(tdewnx(nwb), cloudnx(nwb), phinx(nwb), windnx(nwb), sronx(nwb),  &
            & bgtnx(ngt), bpnx(npi))
     allocate(trq(ntr), trt(ntr), trc(ntr), inq(nbr), dtq(nbr), pre(nbr),      &
            & uhe(nbr), dhe(nbr), inft(nbr))
     allocate(dtt(nbr), prt(nbr), uht(nbr), dht(nbr), inc(nbr), dtc(nbr),      &
            & prc(nbr), uhc(nbr), dhc(nbr))
     allocate(otq(nbr), met(nwb), ext(nwb))
     allocate(nxext1(nwb), nxext2(nwb), extnx(nwb), exto(nwb))
     allocate(ctro(nct, ntr), cino(nct, nbr), qouto(kmx, nbr), cdtro(nct, nbr),&
            & tuho(kmx, nbr), tdho(kmx, nbr), qstro(nst, nbr))
     allocate(ctrnx(nct, ntr), cinnx(nct, nbr), qoutnx(kmx, nbr),              &
            & cdtrnx(nct, nbr), cprnx(nct, nbr), tuhnx(kmx, nbr),              &
            & tdhnx(kmx, nbr))
     allocate(qstrnx(nst, nbr), pumpd(npu), nxpump(npu), epu2(npu), eonpu2(npu)&
            & , eoffpu2(npu), qpu2(npu), dynef(nbr))
     allocate(cuho(kmx, nct, nbr), cdho(kmx, nct, nbr), cuhnx(kmx, nct, nbr),  &
            & cdhnx(kmx, nct, nbr), jjs(nst))
     allocate(inflow_const(nbr), trib_const(ntr), dtrib_const(nbr),            &
            & precip_const(nbr), otqf(nbr), trcf(ntr), dtcf(nbr), incf(nbr),   &
            & trqf(ntr), trtf(ntr), dttf(nbr), dtqf(nbr), inqf(nbr), intf(nbr))
     allocate(euhf(nbr), tuhf(nbr), cuhf(nbr), edhf(nbr), tdhf(nbr), cdhf(nbr))
 
     nxpr1 = 0.0
 
     nxqtr1 = 0.0
 
     nxttr1 = 0.0
 
     nxctr1 = 0.0
 
     nxqin1 = 0.0
 
     nxtin1 = 0.0
 
     nxcin1 = 0.0
 
     nxqdt1 = 0.0
 
     nxtdt1 = 0.0
     nxcdt1 = 0.0
     nxtpr1 = 0.0
     nxcpr1 = 0.0
     nxeuh1 = 0.0
     nxtuh1 = 0.0
     nxcuh1 = 0.0
     nxedh1 = 0.0
     nxtdh1 = 0.0
     nxcdh1 = 0.0
     nxqot1 = 0.0
     nxmet1 = 0.0
     qstrnx = 0.0
     cdtrnx = 0.0
     ctrnx = 0.0
     cinnx = 0.0
     cprnx = 0.0
     cuhnx = 0.0
     cdhnx = 0.0
     qinnx = 0.0
     tinnx = 0.0
     cinnx = 0.0
     nxwsc = 0.0
 
!    Set logical variables
 
     inflow_const = constituents .AND. nacin>0
 
!    Set logical variables
 
     trib_const = constituents .AND. nactr>0
     dtrib_const = constituents .AND. nacdt>0
     precip_const = constituents .AND. nacpr>0
     otqf = .FALSE.
     wshf = .FALSE.
     wdqf = .FALSE.
     trcf = .FALSE.
     trqf = .FALSE.
     trtf = .FALSE.
     dtcf = .FALSE.
     incf = .FALSE.
     prcf = .FALSE.
     prtf = .FALSE.
     prqf = .FALSE.
     metf = .FALSE.
     dynef = .FALSE.
     dttf = .FALSE.
     dtqf = .FALSE.
     inqf = .FALSE.
     intf = .FALSE.
     extf = .FALSE.
     dynpumpf = .FALSE.
     gatef = .FALSE.                   ! SW 9/26/2017
     euhf = 0
     tuhf = 0
     cuhf = 0
     edhf = 0
     tdhf = 0
     cdhf = 0                                                            ! SW 2/28/17
 
!    Open input files
 
     npt = nunit
     shd = npt
     npt = npt + 1
     open(shd, file = shdfn, status = 'OLD')
     read(shd, '(A1)')informat
     if(informat=='$')then
         read(shd, '(/)')
         do i = 1, imx
             read(shd, *)j, shadei(i)
                                     ! SW 3/14/2018 ADDED TO BE COMPATIBLE WITH PREPROCESSOR
             if(shadei(i)<0.0)then
                 backspace(shd)
                 read(shd, *)j, shadei(i), ttlb(i), ttrb(i), cllb(i), clrb(i), &
                           & srlb1(i), srlb2(i), srrb1(i), srrb2(i),           &
                           & (topo(i, j), j = 1, iang), srfjd1(i), srfjd2(i)
             endif
         enddo
     else
         read(shd, '(//(8X,29F8.0))')(shadei(i), ttlb(i), ttrb(i), cllb(i),    &
                                   & clrb(i), srlb1(i), srlb2(i), srrb1(i),    &
                                   & srrb2(i), (topo(i, j), j = 1, iang),      &
                                   & srfjd1(i), srfjd2(i), i = 1, imx)
     endif
     shade = shadei
     wsh = npt
     npt = npt + 1
     open(wsh, file = wscfn, status = 'OLD')
     read(wsh, '(A1)')informat
     if(informat=='$')wshf = .TRUE.
     if(wshf)then
         read(wsh, '(/)')
         read(wsh, *)nxwsc, (wscnx(i), i = 1, imx)
         wsc = wscnx
         read(wsh, *)nxwsc, (wscnx(i), i = 1, imx)
     else
         read(wsh, '(//10F8.0:/(8X,9F8.0))')nxwsc, (wscnx(i), i = 1, imx)
         wsc = wscnx
         read(wsh, '(10F8.0:/(8X,9F8.0))')nxwsc, (wscnx(i), i = 1, imx)
     endif
     do jw = 1, nwb
         met(jw) = npt
         npt = npt + 1
         open(met(jw), file = METFN(jw), status = 'OLD')
         read(met(jw), '(A1)')informat
         if(informat=='$')metf(jw) = .TRUE.
         if(READ_RADIATION(jw))then
             if(metf(jw))then
                 read(met(jw), '(/)')
                 read(met(jw), *)nxmet2(jw), tairnx(jw), tdewnx(jw), windnx(jw)&
                               & , phinx(jw), cloudnx(jw), sronx(jw)
             else
                 read(met(jw), '(//10F8.0)')nxmet2(jw), tairnx(jw), tdewnx(jw),&
                    & windnx(jw), phinx(jw), cloudnx(jw), sronx(jw)
             endif
             sronx(jw) = sronx(jw)*refl
             SRON(jw) = sronx(jw)
             sroo(jw) = SRON(jw)
         elseif(metf(jw))then
             read(met(jw), '(/)')
             read(met(jw), *)nxmet2(jw), tairnx(jw), tdewnx(jw), windnx(jw),   &
                           & phinx(jw), cloudnx(jw)
         else
             read(met(jw), '(//10F8.0)')nxmet2(jw), tairnx(jw), tdewnx(jw),    &
                                      & windnx(jw), phinx(jw), cloudnx(jw)
         endif
         TAIR(jw) = tairnx(jw)
         TDEW(jw) = tdewnx(jw)
         WIND(jw) = windnx(jw)
         PHI(jw) = phinx(jw)
         CLOUD(jw) = cloudnx(jw)
         tairo(jw) = tairnx(jw)
         tdewo(jw) = tdewnx(jw)
         windo(jw) = windnx(jw)
         phio(jw) = phinx(jw)
         cloudo(jw) = cloudnx(jw)
         if(phiset>0)PHI(jw) = phiset
         if(phiset>0)phio(jw) = phiset
         if(READ_RADIATION(jw))then
             if(metf(jw))then
                 read(met(jw), *)nxmet1(jw), tairnx(jw), tdewnx(jw), windnx(jw)&
                               & , phinx(jw), cloudnx(jw), sronx(jw)
             else
                 read(met(jw), '(10F8.0)')nxmet1(jw), tairnx(jw), tdewnx(jw),  &
                    & windnx(jw), phinx(jw), cloudnx(jw), sronx(jw)
             endif
             sronx(jw) = sronx(jw)*refl
         elseif(metf(jw))then
             read(met(jw), *)nxmet1(jw), tairnx(jw), tdewnx(jw), windnx(jw),   &
                           & phinx(jw), cloudnx(jw)
         else
             read(met(jw), '(10F8.0)')nxmet1(jw), tairnx(jw), tdewnx(jw),      &
                                    & windnx(jw), phinx(jw), cloudnx(jw)
         endif
         if(READ_EXTINCTION(jw))then
             ext(jw) = npt
             npt = npt + 1
             open(ext(jw), file = EXTFN(jw), status = 'OLD')
             read(ext(jw), '(A1)')informat
             if(informat=='$')extf(jw) = .TRUE.
             if(extf(jw))then
                 read(ext(jw), '(/)')
                 read(ext(jw), *)nxext2(jw), extnx(jw)
             else
                 read(ext(jw), '(///2F8.0)')nxext2(jw), extnx(jw)
             endif
 
             EXH2O(jw) = extnx(jw)
             exto(jw) = extnx(jw)
 
             if(extf(jw))then
                 read(ext(jw), *)nxext1(jw), extnx(jw)
             else
                 read(ext(jw), '(2F8.0)')nxext1(jw), extnx(jw)
             endif
         endif
    !DO I=CUS(BS(JW)),DS(BE(JW))   ! SW CODE FIX 5-21-15
    !  WIND2(I) = WIND(JW)*WSC(I)*DLOG(2.0D0/Z0(JW))/DLOG(WINDH(JW)/Z0(JW))
    !END DO
     enddo
     if(nwd>0)then
         wdq = npt
         npt = npt + 1
         open(wdq, file = qwdfn, status = 'OLD')
         read(wdq, '(A1)')informat
         if(informat=='$')wdqf = .TRUE.
         if(wdqf)then
             read(wdq, '(/)')
             read(wdq, *)nxqwd2, (qwdnx(jw), jw = 1, nwd)
             do jw = 1, nwd
                 QWD(jw) = qwdnx(jw)
                 qwdo(jw) = qwdnx(jw)
             enddo
             read(wdq, *)nxqwd1, (qwdnx(jw), jw = 1, nwd)
         else
             read(wdq, '(//10F8.0:/(8X,9F8.0))')nxqwd2,                        &
                & (qwdnx(jw), jw = 1, nwd)
             do jw = 1, nwd
                 QWD(jw) = qwdnx(jw)
                 qwdo(jw) = qwdnx(jw)
             enddo
             read(wdq, '(10F8.0:/(8X,9F8.0))')nxqwd1, (qwdnx(jw), jw = 1, nwd)
         endif
     endif
     if(tributaries)then
         do jt = 1, ntr
             trq(jt) = npt
             npt = npt + 1
             trt(jt) = npt
             npt = npt + 1
             open(trq(jt), file = QTRFN(jt), status = 'OLD')
             open(trt(jt), file = TTRFN(jt), status = 'OLD')
 
             read(trq(jt), '(A1)')informat
             if(informat=='$')trqf(jt) = .TRUE.
             if(trqf(jt))then
                 read(trq(jt), '(/)')
                 read(trq(jt), *)nxqtr2(jt), qtrnx(jt)
             else
                 read(trq(jt), '(//2F8.0)')nxqtr2(jt), qtrnx(jt)
             endif
             read(trt(jt), '(A1)')informat
             if(informat=='$')trtf(jt) = .TRUE.
             if(trtf(jt))then
                 read(trt(jt), '(/)')
                 read(trt(jt), *)nxttr2(jt), ttrnx(jt)
             else
                 read(trt(jt), '(//2F8.0)')nxttr2(jt), ttrnx(jt)
             endif
 
 !     READ (TRQ(JT),'(///2F8.0)') NXQTR2(JT),QTRNX(JT)
 !     READ (TRT(JT),'(///2F8.0)') NXTTR2(JT),TTRNX(JT)
             if(trib_const(jt))then
                 trc(jt) = npt
                 npt = npt + 1
                 open(trc(jt), file = CTRFN(jt), status = 'OLD')
                 read(trc(jt), '(A1)')informat
                 if(informat=='$')trcf(jt) = .TRUE.
                 if(trcf(jt))then
                     read(trc(jt), '(/)')
                     read(trc(jt), *)nxctr2(jt),                               &
                                   & (ctrnx(TRCN(jac, jt), jt), jac = 1,       &
                                   & nactr(jt))
                 else
                     read(trc(jt), '(//1000F8.0)')nxctr2(jt),                  &
                        & (ctrnx(TRCN(jac, jt), jt), jac = 1, nactr(jt))
                 endif
             endif
         enddo
         qtr(1:ntr) = qtrnx(1:ntr)
         qtro(1:ntr) = qtrnx(1:ntr)
         ttr(1:ntr) = ttrnx(1:ntr)
         ttro(1:ntr) = ttrnx(1:ntr)
         ctr(:, 1:ntr) = ctrnx(:, 1:ntr)
         ctro(:, 1:ntr) = ctrnx(:, 1:ntr)
         do jt = 1, ntr
             if(trqf(jt))then
                 read(trq(jt), *)nxqtr1(jt), qtrnx(jt)
                                                ! cb 5/22/14
             else
                 read(trq(jt), '(2F8.0)')nxqtr1(jt), qtrnx(jt)
             endif
             if(trtf(jt))then
                 read(trt(jt), *)nxttr1(jt), ttrnx(jt)
                                               ! cb 5/22/14
             else
                 read(trt(jt), '(2F8.0)')nxttr1(jt), ttrnx(jt)
             endif
             if(trib_const(jt))then
                 if(trcf(jt))then
                     read(trc(jt), *)nxctr1(jt),                               &
                                   & (ctrnx(TRCN(jac, jt), jt), jac = 1,       &
                                   & nactr(jt))
                 else
                     read(trc(jt), '(1000F8.0)')nxctr1(jt),                    &
                        & (ctrnx(TRCN(jac, jt), jt), jac = 1, nactr(jt))
                 endif
             endif
         enddo
     endif
     do jw = 1, nwb
         do jb = BS(jw), BE(jw)
             if(UP_FLOW(jb))then
                 if(.NOT.INTERNAL_FLOW(jb) .AND. .NOT.DAM_INFLOW(jb))then                                             !TC 08/03/04 RA 1/13/06
                     inq(jb) = npt
                     npt = npt + 1
                     inft(jb) = npt
                     npt = npt + 1
                     open(inq(jb), file = QINFN(jb), status = 'OLD')
                     open(inft(jb), file = TINFN(jb), status = 'OLD')
                     read(inq(jb), '(A1)')informat
                     if(informat=='$')inqf(jb) = .TRUE.
                     if(inqf(jb))then
                         read(inq(jb), '(/)')
                         read(inq(jb), *)nxqin2(jb), qinnx(jb)
                     else
                         read(inq(jb), '(//2F8.0)')nxqin2(jb), qinnx(jb)
                     endif
                     read(inft(jb), '(A1)')informat
                     if(informat=='$')intf(jb) = .TRUE.
                     if(intf(jb))then
                         read(inft(jb), '(/)')
                         read(inft(jb), *)nxtin2(jb), tinnx(jb)
                     else
                         read(inft(jb), '(//2F8.0)')nxtin2(jb), tinnx(jb)
                     endif
 
     !     READ (INQ(JB), '(///2F8.0)') NXQIN2(JB),QINNX(JB)
     !     READ (INFT(JB),'(///2F8.0)') NXTIN2(JB),TINNX(JB)
                     if(inflow_const(jb))then
                         inc(jb) = npt
                         npt = npt + 1
                         open(inc(jb), file = CINFN(jb), status = 'OLD')
                         read(inc(jb), '(A1)')informat
                         if(informat=='$')incf(jb) = .TRUE.
                         if(incf(jb))then
                             read(inc(jb), '(/)')
                             read(inc(jb), *)nxcin2(jb),                       &
                                & (cinnx(INCN(jc, jb), jb), jc = 1, nacin(jb))
                         else
                             read(inc(jb), '(//1000F8.0)')nxcin2(jb),          &
                                & (cinnx(INCN(jc, jb), jb), jc = 1, nacin(jb))
                         endif
                     endif
                 endif
                 QIN(jb) = qinnx(jb)
                 QIND(jb) = qinnx(jb)
                 qino(jb) = qinnx(jb)
                 TIN(jb) = tinnx(jb)
                 TIND(jb) = tinnx(jb)
                 tino(jb) = tinnx(jb)
                 cin(:, jb) = cinnx(:, jb)
                 cind(:, jb) = cinnx(:, jb)
                 cino(:, jb) = cinnx(:, jb)
                 if(.NOT.INTERNAL_FLOW(jb) .AND. .NOT.DAM_INFLOW(jb))then                                             !TC 08/03/04  RA 1/13/06
                     if(inqf(jb))then
                         read(inq(jb), *)nxqin1(jb), qinnx(jb)
                     else
                         read(inq(jb), '(2F8.0)')nxqin1(jb), qinnx(jb)
                     endif
                     if(intf(jb))then
                         read(inft(jb), *)nxtin1(jb), tinnx(jb)
                     else
                         read(inft(jb), '(2F8.0)')nxtin1(jb), tinnx(jb)
                     endif
 
                     if(inflow_const(jb))then
                         if(incf(jb))then
                             read(inc(jb), *)nxcin1(jb),                       &
                                & (cinnx(INCN(jc, jb), jb), jc = 1, nacin(jb))
                         else
                             read(inc(jb), '(1000F8.0)')nxcin1(jb),            &
                                & (cinnx(INCN(jc, jb), jb), jc = 1, nacin(jb))
                         endif
                     endif
                 endif
             endif
             if(DN_FLOW(jb))then
                 if(NSTR(jb)>0)then
                     otq(jb) = npt
                     npt = npt + 1
                     open(otq(jb), file = QOTFN(jb), status = 'OLD')
                     read(otq(jb), '(A1)')informat
                     if(informat=='$')otqf(jb) = .TRUE.
                     if(otqf(jb))then
                         read(otq(jb), '(/)')
                         read(otq(jb), *)nxqot2(jb),                           &
                            & (qstrnx(js, jb), js = 1, NSTR(jb))
                         qstr(:, jb) = qstrnx(:, jb)
                         qstro(:, jb) = qstrnx(:, jb)
                         read(otq(jb), *)nxqot1(jb),                           &
                            & (qstrnx(js, jb), js = 1, NSTR(jb))
                     else
                         read(otq(jb), '(//10F8.0:/(8X,9F8.0))')nxqot2(jb),    &
                            & (qstrnx(js, jb), js = 1, NSTR(jb))
                         qstr(:, jb) = qstrnx(:, jb)
                         qstro(:, jb) = qstrnx(:, jb)
                         read(otq(jb), '(10F8.0:/(8X,9F8.0))')nxqot1(jb),      &
                            & (qstrnx(js, jb), js = 1, NSTR(jb))
                     endif
 
                     if(DYNSTRUC(jb)=='      ON')then
                         odyns(jb) = npt
                         npt = npt + 1
                         write(segnum, '(I0)')jb
                         segnum = ADJUSTL(segnum)
                         l = LEN_TRIM(segnum)
                         open(odyns(jb), file = 'dynselev' // segnum(1:l)      &
                             & // '.npt', status = 'OLD')
                         read(odyns(jb), '(A1)')informat
                         if(informat=='$')dynef(jb) = .TRUE.
                         read(odyns(jb), *)njs
                         do j = 1, njs
                             read(odyns(jb), *)jjs(j)
                         enddo
                         if(dynef(jb))then
                             read(odyns(jb), *)
                             read(odyns(jb), *)nxdyns(jb),                     &
                                & (estr(jjs(j), jb), j = 1, njs)
                             read(odyns(jb), *)nxdyns(jb),                     &
                                & (nxestrt(jjs(j), jb), j = 1, njs)
                         else
                             read(odyns(jb), '(/10F8.0:/(8X,9F8.0))')nxdyns(jb)&
                                & , (estr(jjs(j), jb), j = 1, njs)
                             read(odyns(jb), '(10F8.0:/(8X,9F8.0))')nxdyns(jb),&
                                & (nxestrt(jjs(j), jb), j = 1, njs)
                         endif
                     endif
                 endif
             endif
             if(PRECIPITATION(jw))then
                 pre(jb) = npt
                 npt = npt + 1
                 prt(jb) = npt
                 npt = npt + 1
                 open(pre(jb), file = PREFN(jb), status = 'OLD')
                 open(prt(jb), file = TPRFN(jb), status = 'OLD')
                 read(pre(jb), '(A1)')informat
                 if(informat=='$')prqf(jb) = .TRUE.
                 if(prqf(jb))then
                     read(pre(jb), '(/)')
                     read(pre(jb), *)nxpr2(jb), prnx(jb)
                 else
                     read(pre(jb), '(//2F8.0)')nxpr2(jb), prnx(jb)
                 endif
                 read(prt(jb), '(A1)')informat
                 if(informat=='$')prtf(jb) = .TRUE.
                 if(prtf(jb))then
                     read(prt(jb), '(/)')
                     read(prt(jb), *)nxtpr2(jb), tprnx(jb)
                 else
                     read(prt(jb), '(//2F8.0)')nxtpr2(jb), tprnx(jb)
                 endif
 
        !READ (PRE(JB),'(///2F8.0)') NXPR2(JB), PRNX(JB)
        !READ (PRT(JB),'(///2F8.0)') NXTPR2(JB),TPRNX(JB)
                 if(precip_const(jb))then
                     prc(jb) = npt
                     npt = npt + 1
                     open(prc(jb), file = CPRFN(jb), status = 'OLD')
                     read(prc(jb), '(A1)')informat
                     if(informat=='$')prcf(jb) = .TRUE.
                     if(prcf(jb))then
                         read(prc(jb), '(/)')
                         read(prc(jb), *)nxcpr2(jb),                           &
                            & (cprnx(PRCN(jac, jb), jb), jac = 1, nacpr(jb))
                     else
                         read(prc(jb), '(//1000F8.0)')nxcpr2(jb),              &
                            & (cprnx(PRCN(jac, jb), jb), jac = 1, nacpr(jb))
                     endif
                 endif
                 PR(jb) = prnx(jb)
                 TPR(jb) = tprnx(jb)
                 cpr(:, jb) = cprnx(:, jb)
                 if(prqf(jb))then
                     read(pre(jb), *)nxpr1(jb), prnx(jb)
                 else
                     read(pre(jb), '(2F8.0)')nxpr1(jb), prnx(jb)
                 endif
 
                 if(prtf(jb))then
                     read(prt(jb), *)nxtpr1(jb), tprnx(jb)
                 else
                     read(prt(jb), '(2F8.0)')nxtpr1(jb), tprnx(jb)
                 endif
        !READ (PRE(JB),'(2F8.0)') NXPR1(JB), PRNX(JB)
        !READ (PRT(JB),'(2F8.0)') NXTPR1(JB),TPRNX(JB)
                 if(precip_const(jb))then
                     if(prcf(jb))then
                         read(prc(jb), *)nxcpr1(jb),                           &
                            & (cprnx(PRCN(jac, jb), jb), jac = 1, nacpr(jb))
                     else
                         read(prc(jb), '(1000F8.0)')nxcpr1(jb),                &
                            & (cprnx(PRCN(jac, jb), jb), jac = 1, nacpr(jb))
                     endif
                 endif
             endif
             if(DIST_TRIBS(jb))then
                 dtq(jb) = npt
                 npt = npt + 1
                 dtt(jb) = npt
                 npt = npt + 1
                 open(dtq(jb), file = QDTFN(jb), status = 'OLD')
                 open(dtt(jb), file = TDTFN(jb), status = 'OLD')
 
                 read(dtq(jb), '(A1)')informat
                 if(informat=='$')dtqf(jb) = .TRUE.
                 if(dtqf(jb))then
                     read(dtq(jb), '(/)')
                     read(dtq(jb), *)nxqdt2(jb), qdtrnx(jb)
                 else
                     read(dtq(jb), '(//2F8.0)')nxqdt2(jb), qdtrnx(jb)
                 endif
                 read(dtt(jb), '(A1)')informat
                 if(informat=='$')dttf(jb) = .TRUE.
                 if(dttf(jb))then
                     read(dtt(jb), '(/)')
                     read(dtt(jb), *)nxtdt2(jb), tdtrnx(jb)
                 else
                     read(dtt(jb), '(//2F8.0)')nxtdt2(jb), tdtrnx(jb)
                 endif
 
 
  !      READ (DTQ(JB),'(///2F8.0)') NXQDT2(JB),QDTRNX(JB)
  !      READ (DTT(JB),'(///2F8.0)') NXTDT2(JB),TDTRNX(JB)
                 if(dtrib_const(jb))then
                     dtc(jb) = npt
                     npt = npt + 1
                     open(dtc(jb), file = CDTFN(jb), status = 'OLD')
                     read(dtc(jb), '(A1)')informat
                     if(informat=='$')dtcf(jb) = .TRUE.
                     if(dtcf(jb))then
                         read(dtc(jb), '(/)')
                         read(dtc(jb), *)nxcdt2(jb),                           &
                            & (cdtrnx(DTCN(jac, jb), jb), jac = 1, nacdt(jb))
                     else
                         read(dtc(jb), '(//1000F8.0)')nxcdt2(jb),              &
                            & (cdtrnx(DTCN(jac, jb), jb), jac = 1, nacdt(jb))
                     endif
                 endif
                 QDTR(jb) = qdtrnx(jb)
                 qdtro(jb) = qdtrnx(jb)
                 TDTR(jb) = tdtrnx(jb)
                 tdtro(jb) = tdtrnx(jb)
                 cdtr(:, jb) = cdtrnx(:, jb)
                 cdtro(:, jb) = cdtrnx(:, jb)
                 if(dtqf(jb))then
                     read(dtq(jb), *)nxqdt1(jb), qdtrnx(jb)
                 else
                     read(dtq(jb), '(2F8.0)')nxqdt1(jb), qdtrnx(jb)
                 endif
 
                 if(dttf(jb))then
                     read(dtt(jb), *)nxtdt1(jb), tdtrnx(jb)
                 else
                     read(dtt(jb), '(2F8.0)')nxtdt1(jb), tdtrnx(jb)
                 endif
 
                 if(dtrib_const(jb))then
                     if(dtcf(jb))then
                         read(dtc(jb), *)nxcdt1(jb),                           &
                            & (cdtrnx(DTCN(jac, jb), jb), jac = 1, nacdt(jb))
                     else
                         read(dtc(jb), '(1000F8.0)')nxcdt1(jb),                &
                            & (cdtrnx(DTCN(jac, jb), jb), jac = 1, nacdt(jb))
                     endif
                 endif
             endif
             if(UH_EXTERNAL(jb))then
                 uhe(jb) = npt
                 npt = npt + 1
                 uht(jb) = npt
                 npt = npt + 1
                 open(uhe(jb), file = EUHFN(jb), status = 'OLD')
                 open(uht(jb), file = TUHFN(jb), status = 'OLD')
 
                 read(uhe(jb), '(A1)')informat
                 if(informat=='$')euhf(jb) = 1
 
                 read(uht(jb), '(A2)')informat2
                 if(informat2=='$T')then
                     tuhf(jb) = 1
                 elseif(informat2(1:1)=='$')then
                     tuhf(jb) = 2
                 endif
 
                 if(euhf(jb)>0)then
                     read(uhe(jb), '(/)')
                     read(uhe(jb), *)nxeuh2(jb), eluhnx(jb)
                 else
                     read(uhe(jb), '(//2F8.0)')nxeuh2(jb), eluhnx(jb)
                 endif
 
                 if(tuhf(jb)==1)then
                     read(uht(jb), '(/)')
                     read(uht(jb), *)nxtuh2(jb), xx(1)
                     tuhnx(2:KB(US(jb)), jb) = xx(1)
                 elseif(tuhf(jb)==2)then
                     read(uht(jb), '(/)')
                     read(uht(jb), *)nxtuh2(jb),                               &
                                   & (tuhnx(k, jb), k = 2, KB(US(jb)))
                 else
                     read(uht(jb), '(//10F8.0:/(8X,9F8.0))')nxtuh2(jb),        &
                        & (tuhnx(k, jb), k = 2, KB(US(jb)))
                 endif
 
       ! READ (UHE(JB),'(///2F8.0)')              NXEUH2(JB), ELUHNX(JB)
       ! READ (UHT(JB),'(///10F8.0:/(8X,9F8.0))') NXTUH2(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
                 if(constituents)then
                     uhc(jb) = npt
                     npt = npt + 1
                     open(uhc(jb), file = CUHFN(jb), status = 'OLD')
 
                     read(uhc(jb), '(A2)')informat2
                     if(informat2=='$T')then
                         cuhf(jb) = 1
                     elseif(informat2(1:1)=='$')then
                         cuhf(jb) = 2
                     endif
 
                     read(uhc(jb), '(/)')
 
        !  READ (UHC(JB),'(//)')
                     if(cuhf(jb)==1)then
                         read(uhc(jb), *)nxcuh2(jb),                           &
                            & (xx(CN(jac)), jac = 1, nac)
                         do jac = 1, nac
                             cuhnx(2:KB(US(jb)), CN(jac), jb) = xx(CN(jac))
                         enddo
 
                     else
 
 
                         do jac = 1, nac
            !IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH2(JB),(CUHNX(K,CN(JAC),JB),     &
            !                                                  K=2,KB(US(JB)))
 
                             if(cuhf(jb)==2)then
                                 read(uhc(jb), *)nxcuh2(jb),                   &
                                    & (cuhnx(k, CN(jac), jb), k = 2, KB(US(jb))&
                                    & )
                             else
                                 read(uhc(jb), '(10F8.0:/(8X,9F8.0))')         &
                                    & nxcuh2(jb),                              &
                                    & (cuhnx(k, CN(jac), jb), k = 2, KB(US(jb))&
                                    & )
                             endif
 
 
                         enddo
                     endif
 
                 endif
                 ELUH(jb) = eluhnx(jb)
                 eluho(jb) = eluhnx(jb)
                 tuh(:, jb) = tuhnx(:, jb)
                 tuho(:, jb) = tuhnx(:, jb)
                 cuh(:, :, jb) = cuhnx(:, :, jb)
                 cuho(:, :, jb) = cuhnx(:, :, jb)
 
                 if(euhf(jb)>0)then
                     read(uhe(jb), *)nxeuh1(jb), eluhnx(jb)
                 else
                     read(uhe(jb), '(2F8.0)')nxeuh1(jb), eluhnx(jb)
                 endif
 
        !READ (UHE(JB),'(2F8.0)')              NXEUH1(JB), ELUHNX(JB)
 
                 if(tuhf(jb)==1)then
                     read(uht(jb), *)nxtuh1(jb), xx(1)
                     tuhnx(2:KB(US(jb)), jb) = xx(1)
                 elseif(tuhf(jb)==2)then
                     read(uht(jb), *)nxtuh1(jb),                               &
                                   & (tuhnx(k, jb), k = 2, KB(US(jb)))
                 else
                     read(uht(jb), '(10F8.0:/(8X,9F8.0))')nxtuh1(jb),          &
                        & (tuhnx(k, jb), k = 2, KB(US(jb)))
                 endif
 
       ! READ (UHT(JB),'(10F8.0:/(8X,9F8.0))') NXTUH1(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
                 if(constituents)then
 
                     if(cuhf(jb)==1)then
                         read(uhc(jb), *)nxcuh1(jb),                           &
                            & (xx(CN(jac)), jac = 1, nac)
                         do jac = 1, nac
                             cuhnx(2:KB(US(jb)), CN(jac), jb) = xx(CN(jac))
                         enddo
 
                     else
 
 
                         do jac = 1, nac
 
                             if(cuhf(jb)==2)then
                                 read(uhc(jb), *)nxcuh1(jb),                   &
                                    & (cuhnx(k, CN(jac), jb), k = 2, KB(US(jb))&
                                    & )
                             else
                                 read(uhc(jb), '(10F8.0:/(8X,9F8.0))')         &
                                    & nxcuh1(jb),                              &
                                    & (cuhnx(k, CN(jac), jb), k = 2, KB(US(jb))&
                                    & )
                             endif
 
                         enddo
                     endif
 
 
          !DO JAC=1,NAC
          !  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH1(JB),(CUHNX(K,CN(JAC),JB),     &
          !                                                    K=2,KB(US(JB)))
          !END DO
                 endif
             endif
             if(DH_EXTERNAL(jb))then
                 dhe(jb) = npt
                 npt = npt + 1
                 dht(jb) = npt
                 npt = npt + 1
                 open(dhe(jb), file = EDHFN(jb), status = 'OLD')
                 open(dht(jb), file = TDHFN(jb), status = 'OLD')
 
                 read(dhe(jb), '(A1)')informat
                 if(informat=='$')edhf(jb) = 1
 
                 read(dht(jb), '(A2)')informat2
                 if(informat2=='$T')then
                     tdhf(jb) = 1
                 elseif(informat2(1:1)=='$')then
                     tdhf(jb) = 2
                 endif
 
                 if(edhf(jb)>0)then
                     read(dhe(jb), '(/)')
                     read(dhe(jb), *)nxedh2(jb), eldhnx(jb)
                 else
                     read(dhe(jb), '(//2F8.0)')nxedh2(jb), eldhnx(jb)
                 endif
 
                 if(tdhf(jb)==1)then
                     read(dht(jb), '(/)')
                     read(dht(jb), *)nxtdh2(jb), xx(1)
                     tdhnx(2:KB(DS(jb)), jb) = xx(1)
                 elseif(tdhf(jb)==2)then
                     read(dht(jb), '(/)')
                     read(dht(jb), *)nxtdh2(jb),                               &
                                   & (tdhnx(k, jb), k = 2, KB(DS(jb)))
                 else
                     read(dht(jb), '(//10F8.0:/(8X,9F8.0))')nxtdh2(jb),        &
                        & (tdhnx(k, jb), k = 2, KB(DS(jb)))
                 endif
 
        !READ (DHE(JB),'(///10F8.0)')             NXEDH2(JB),ELDHNX(JB)
        !READ (DHT(JB),'(///10F8.0:/(8X,9F8.0))') NXTDH2(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
                 if(constituents)then
                     dhc(jb) = npt
                     npt = npt + 1
                     open(dhc(jb), file = CDHFN(jb), status = 'OLD')
 
                     read(dhc(jb), '(A2)')informat2
                     if(informat2=='$T')then
                         cdhf(jb) = 1
                     elseif(informat2(1:1)=='$')then
                         cdhf(jb) = 2
                     endif
 
                     read(dhc(jb), '(/)')
 
                     if(cdhf(jb)==1)then
                         read(dhc(jb), *)nxcdh2(jb),                           &
                            & (xx(CN(jac)), jac = 1, nac)
                         do jac = 1, nac
                             cdhnx(2:KB(DS(jb)), CN(jac), jb) = xx(CN(jac))
                         enddo
                     else
                         do jac = 1, nac
                             if(cdhf(jb)==2)then
                                 read(dhc(jb), *)nxcdh2(jb),                   &
                                    & (cdhnx(k, CN(jac), jb), k = 2, KB(DS(jb))&
                                    & )
                             else
                                 read(dhc(jb), '(10F8.0:/(8X,9F8.0))')         &
                                    & nxcdh2(jb),                              &
                                    & (cdhnx(k, CN(jac), jb), k = 2, KB(DS(jb))&
                                    & )
                             endif
                         enddo
                     endif
 
          !READ (DHC(JB),'(//)')
          !DO JAC=1,NAC
          !  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH2(JB),(CDHNX(K,CN(JAC),JB),     &
          !                                                    K=2,KB(DS(JB)))
          !END DO
                 endif
                 ELDH(jb) = eldhnx(jb)
                 eldho(jb) = eldhnx(jb)
                 tdh(:, jb) = tdhnx(:, jb)
                 tdho(:, jb) = tdhnx(:, jb)
                 cdh(:, :, jb) = cdhnx(:, :, jb)
                 cdho(:, :, jb) = cdhnx(:, :, jb)
 
                 if(edhf(jb)>0)then
                     read(dhe(jb), *)nxedh1(jb), eldhnx(jb)
                 else
                     read(dhe(jb), '(2F8.0)')nxedh1(jb), eldhnx(jb)
                 endif
 
                 if(tdhf(jb)==1)then
                     read(dht(jb), *)nxtdh1(jb), xx(1)
                     tdhnx(2:KB(DS(jb)), jb) = xx(1)
                 elseif(tdhf(jb)==2)then
                     read(dht(jb), *)nxtdh1(jb),                               &
                                   & (tdhnx(k, jb), k = 2, KB(DS(jb)))
                 else
                     read(dht(jb), '(10F8.0:/(8X,9F8.0))')nxtdh1(jb),          &
                        & (tdhnx(k, jb), k = 2, KB(DS(jb)))
                 endif
 
        !READ (DHE(JB),'(10F8.0)')             NXEDH1(JB),ELDHNX(JB)
        !READ (DHT(JB),'(10F8.0:/(8X,9F8.0))') NXTDH1(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
                 if(constituents)then
 
                     if(cdhf(jb)==1)then
                         read(dhc(jb), *)nxcdh1(jb),                           &
                            & (xx(CN(jac)), jac = 1, nac)
                         do jac = 1, nac
                             cdhnx(2:KB(DS(jb)), CN(jac), jb) = xx(CN(jac))
                         enddo
                     else
                         do jac = 1, nac
                             if(cdhf(jb)==2)then
                                 read(dhc(jb), *)nxcdh1(jb),                   &
                                    & (cdhnx(k, CN(jac), jb), k = 2, KB(DS(jb))&
                                    & )
                             else
                                 read(dhc(jb), '(10F8.0:/(8X,9F8.0))')         &
                                    & nxcdh1(jb),                              &
                                    & (cdhnx(k, CN(jac), jb), k = 2, KB(DS(jb))&
                                    & )
                             endif
                         enddo
                     endif
 
          !DO JAC=1,NAC
          !  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH1(JB),(CDHNX(K,CN(JAC),JB),     &
          !                                                    K=2,KB(DS(JB)))
          !END DO
                 endif
             endif
         enddo
     enddo
     if(gates)then
         gtq = npt
         npt = npt + 1
         open(gtq, file = qgtfn, status = 'OLD')
         read(gtq, '(A1)')informat
         if(informat=='$')gatef = .TRUE.
    ! Added code is to allow for computing flow based on BGT elevation - this can be a target water level in the reservoir - but allowing that flow to be removed from a different elevation
    !READ(GTQ,*)                           ! SW 2/25/11
         if(gatef)then
             read(gtq, '(A8)')gt2char
 
             if(gt2char=='EGT2ELEV')then
                 rewind(gtq)
                 read(gtq, *)
                 read(gtq, *)gt2char, (egt2(jg), jg = 1, ngt)
             endif
 
             read(gtq, *)
             read(gtq, *)nxqgt2, (bgtnx(jg), jg = 1, ngt)
             where(dyngtc=='     ZGT')
                 egt = bgtnx
                 egto = bgtnx
                 bgt = 1.0
                 g1gt = 1.0
                 g2gt = 1.0
             elsewhere
                 bgt = bgtnx
                 bgto = bgtnx
             endwhere
             read(gtq, *)nxqgt, (bgtnx(jg), jg = 1, ngt)
 
         else
 
             read(gtq, '(A8)')gt2char
 
             if(gt2char=='EGT2ELEV')then
                 rewind(gtq)
                 read(gtq, *)
                 read(gtq, '(8X,1000F8.0)')(egt2(jg), jg = 1, ngt)
             endif
 
             read(gtq, *)
             read(gtq, '(1000F8.0)')nxqgt2, (bgtnx(jg), jg = 1, ngt)
!            READ (GTQ,'(///1000F8.0)') NXQGT2,(BGTNX(JG),JG=1,NGT)
             where(dyngtc=='     ZGT')
                 egt = bgtnx
                 egto = bgtnx
                 bgt = 1.0
                 g1gt = 1.0
                 g2gt = 1.0
             elsewhere
                 bgt = bgtnx
                 bgto = bgtnx
             endwhere
             read(gtq, '(1000F8.0)')nxqgt, (bgtnx(jg), jg = 1, ngt)
         endif
 
     endif
     if(pipes)then                                        ! SW 5/5/10
         iopenpipe = 0
         do j = 1, npi
             if(DYNPIPE(j)=='      ON')then
                 iopenpipe = 1
                 exit
             endif
         enddo
         if(iopenpipe==1)then
             piped = npt
             npt = npt + 1
             open(piped, file = 'dynpipe.npt', status = 'OLD')
             read(piped, '(///1000F8.0)')nxqpt, (bpnx(j), j = 1, npi)
             bp = bpnx
             read(piped, '(1000F8.0)')nxqpt, (bpnx(j), j = 1, npi)
         endif
     endif
     if(pumps)then
         do j = 1, npu
             if(DYNPUMP(j)=='      ON')then
                 write(segnum, '(I0)')j
                 segnum = ADJUSTL(segnum)
                 l = LEN_TRIM(segnum)
                 pumpd(j) = npt
                 npt = npt + 1
                 open(pumpd(j), file = 'dynpump' // segnum(1:l) // '.npt',     &
                    & status = 'OLD')
                 read(pumpd(j), '(A1)')informat
                 if(informat=='$')dynpumpf(j) = .TRUE.
                 if(dynpumpf(j))then
                     read(pumpd(j), '(/)')
                     read(pumpd(j), *)nxpump(j), epu2(j), eonpu2(j), eoffpu2(j)&
                                    & , qpu2(j)
                     EPU(j) = epu2(j)
                     EONPU(j) = eonpu2(j)
                     EOFFPU(j) = eoffpu2(j)
                     QPU(j) = qpu2(j)
                     read(pumpd(j), *)nxpump(j), epu2(j), eonpu2(j), eoffpu2(j)&
                                    & , qpu2(j)
                 else
                     read(pumpd(j), '(//1000F8.0)')nxpump(j), epu2(j),         &
                        & eonpu2(j), eoffpu2(j), qpu2(j)
                     EPU(j) = epu2(j)
                     EONPU(j) = eonpu2(j)
                     EOFFPU(j) = eoffpu2(j)
                     QPU(j) = qpu2(j)
                     read(pumpd(j), '(1000F8.0)')nxpump(j), epu2(j), eonpu2(j),&
                        & eoffpu2(j), qpu2(j)
                 endif
      !READ (PUMPD(J),'(///1000F8.0)') NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
      !  EPU(J)=EPU2(J)
      !  EONPU(J)=EONPU2(J)
      !  EOFFPU(J)=EOFFPU2(J)
      !  QPU(J)=QPU2(J)
      !READ (PUMPD(J),'(1000F8.0)') NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
             endif
         enddo
     endif
 
     nopen = npt - 1
     dynamic_shade = shadei<0
     nunit = npt
     return
 
!***********************************************************************************************************************************
!**  R E A D  I N P U T  D A T A                                              
!***********************************************************************************************************************************
 
!    **
     entry READ_INPUT_DATA(Nxtvd)
     Nxtvd = 1.0E10
 
!    Meteorological data
 
     do while (jday>=nxwsc)
         wsc = wscnx
         if(wshf)then
             read(wsh, *)nxwsc, (wscnx(i), i = 1, imx)
         else
             read(wsh, '(10F8.0:/(8X,9F8.0))')nxwsc, (wscnx(i), i = 1, imx)
         endif
     enddo
     do jw = 1, nwb
         do while (jday>=nxmet1(jw))
             TDEW(jw) = tdewnx(jw)
             tdewo(jw) = tdewnx(jw)
             WIND(jw) = windnx(jw)
             windo(jw) = windnx(jw)
             PHI(jw) = phinx(jw)
             phio(jw) = phinx(jw)
             if(phiset>0)PHI(jw) = phiset
             if(phiset>0)phio(jw) = phiset
             TAIR(jw) = tairnx(jw)
             tairo(jw) = tairnx(jw)
             CLOUD(jw) = cloudnx(jw)
             cloudo(jw) = cloudnx(jw)
             nxmet2(jw) = nxmet1(jw)
             if(READ_RADIATION(jw))then
                 SRON(jw) = sronx(jw)
                 sroo(jw) = SRON(jw)
                 if(metf(jw))then
                     read(met(jw), *)nxmet1(jw), tairnx(jw), tdewnx(jw),       &
                                   & windnx(jw), phinx(jw), cloudnx(jw),       &
                                   & sronx(jw)
                 else
                     read(met(jw), '(7F8.0)')nxmet1(jw), tairnx(jw), tdewnx(jw)&
                        & , windnx(jw), phinx(jw), cloudnx(jw), sronx(jw)
                 endif
                 sronx(jw) = sronx(jw)*refl
             elseif(metf(jw))then
                 read(met(jw), *)nxmet1(jw), tairnx(jw), tdewnx(jw), windnx(jw)&
                               & , phinx(jw), cloudnx(jw)
             else
                 read(met(jw), '(6F8.0)')nxmet1(jw), tairnx(jw), tdewnx(jw),   &
                    & windnx(jw), phinx(jw), cloudnx(jw)
             endif
         enddo
         Nxtvd = MIN(Nxtvd, nxmet1(jw))
         if(READ_EXTINCTION(jw))then
             do while (jday>=nxext1(jw))
                 EXH2O(jw) = extnx(jw)
                 exto(jw) = extnx(jw)
                 nxext2(jw) = nxext1(jw)
                 if(extf(jw))then
                     read(ext(jw), *)nxext1(jw), extnx(jw)
                 else
                     read(ext(jw), '(2F8.0)')nxext1(jw), extnx(jw)
                 endif
             enddo
         endif
    !DO I=CUS(BS(JW)),DS(BE(JW))   ! SW CODE FIX 5-21-15
    !  WIND2(I) = WIND(JW)*WSC(I)*DLOG(2.0D0/Z0(JW))/DLOG(WINDH(JW)/Z0(JW))    ! old value  z0 == 0.003
    !END DO
     enddo
 
!    Withdrawals
 
     if(nwd>0)then
         do while (jday>=nxqwd1)
             nxqwd2 = nxqwd1
             do jwd = 1, nwd
                 QWD(jwd) = qwdnx(jwd)
                 qwdo(jwd) = qwdnx(jwd)
             enddo
             if(wdqf)then
                 read(wdq, *)nxqwd1, (qwdnx(jwd), jwd = 1, nwd)
             else
                 read(wdq, '(10F8.0:/(8X,9F8.0))')nxqwd1,                      &
                    & (qwdnx(jwd), jwd = 1, nwd)
             endif
         enddo
         Nxtvd = MIN(Nxtvd, nxqwd1)
     endif
 
!    Tributaries
 
     if(tributaries)then
         do jt = 1, ntr
 
!****        Inflow
 
             do while (jday>=nxqtr1(jt))
                 qtr(jt) = qtrnx(jt)
                 qtro(jt) = qtrnx(jt)
                 nxqtr2(jt) = nxqtr1(jt)
 
                 if(trqf(jt))then
                     read(trq(jt), *)nxqtr1(jt), qtrnx(jt)
                 else
                     read(trq(jt), '(2F8.0)')nxqtr1(jt), qtrnx(jt)
                 endif
 
 
  !      READ (TRQ(JT),'(2F8.0)') NXQTR1(JT),QTRNX(JT)
             enddo
             Nxtvd = MIN(Nxtvd, nxqtr1(jt))
 
!****        Inflow temperatures
 
             if(jday>=nxttr1(jt))then
                 do while (jday>=nxttr1(jt))
                     ttr(jt) = ttrnx(jt)
                     ttro(jt) = ttrnx(jt)
                     nxttr2(jt) = nxttr1(jt)
 
                     if(trtf(jt))then
                         read(trt(jt), *)nxttr1(jt), ttrnx(jt)
                     else
                         read(trt(jt), '(2F8.0)')nxttr1(jt), ttrnx(jt)
                     endif
 
 
   !       READ (TRT(JT),'(2F8.0)') NXTTR1(JT),TTRNX(JT)
                 enddo
             endif
             Nxtvd = MIN(Nxtvd, nxttr1(jt))
 
!****        Inflow constituent concentrations
 
             if(trib_const(jt))then
                 do while (jday>=nxctr1(jt))
                     ctr(TRCN(1:nactr(jt), jt), jt)                            &
                       & = ctrnx(TRCN(1:nactr(jt), jt), jt)
                     ctro(TRCN(1:nactr(jt), jt), jt)                           &
                       & = ctrnx(TRCN(1:nactr(jt), jt), jt)
                     nxctr2(jt) = nxctr1(jt)
                     if(trcf(jt))then
                         read(trc(jt), *)nxctr1(jt),                           &
                            & (ctrnx(TRCN(jac, jt), jt), jac = 1, nactr(jt))
                     else
                         read(trc(jt), '(1000F8.0)')nxctr1(jt),                &
                            & (ctrnx(TRCN(jac, jt), jt), jac = 1, nactr(jt))
                     endif
                 enddo
                 Nxtvd = MIN(Nxtvd, nxctr1(jt))
             endif
         enddo
     endif
 
!    Branch related inputs
 
     do jw = 1, nwb
         do jb = BS(jw), BE(jw)
 
!****        Inflow
 
             if(UP_FLOW(jb))then
                 if(.NOT.INTERNAL_FLOW(jb) .AND. .NOT.DAM_INFLOW(jb))then                                             !TC 08/03/04 RA 1/13/06
                     do while (jday>=nxqin1(jb))
                         QIND(jb) = qinnx(jb)
                         qino(jb) = qinnx(jb)
                         nxqin2(jb) = nxqin1(jb)
 
                         if(inqf(jb))then
                             read(inq(jb), *)nxqin1(jb), qinnx(jb)
                         else
                             read(inq(jb), '(2F8.0)')nxqin1(jb), qinnx(jb)
                         endif
 
        !    READ (INQ(JB),'(2F8.0)') NXQIN1(JB),QINNX(JB)
                     enddo
                     Nxtvd = MIN(Nxtvd, nxqin1(jb))
 
!********            Inflow temperature
 
                     do while (jday>=nxtin1(jb))
                         TIND(jb) = tinnx(jb)
                         tino(jb) = tinnx(jb)
                         nxtin2(jb) = nxtin1(jb)
 
                         if(intf(jb))then
                             read(inft(jb), *)nxtin1(jb), tinnx(jb)
                         else
                             read(inft(jb), '(2F8.0)')nxtin1(jb), tinnx(jb)
                         endif
 
     !       READ (INFT(JB),'(2F8.0)') NXTIN1(JB),TINNX(JB)
                     enddo
                     Nxtvd = MIN(Nxtvd, nxtin1(jb))
 
!********            Inflow constituent concentrations
 
                     if(inflow_const(jb))then
                         do while (jday>=nxcin1(jb))
                             cind(INCN(1:nacin(jb), jb), jb)                   &
                               & = cinnx(INCN(1:nacin(jb), jb), jb)
                             cino(INCN(1:nacin(jb), jb), jb)                   &
                               & = cinnx(INCN(1:nacin(jb), jb), jb)
                             nxcin2(jb) = nxcin1(jb)
                             if(incf(jb))then
                                 read(inc(jb), *)nxcin1(jb),                   &
                                    & (cinnx(INCN(jac, jb), jb), jac = 1,      &
                                    & nacin(jb))
                             else
                                 read(inc(jb), '(1000F8.0)')nxcin1(jb),        &
                                    & (cinnx(INCN(jac, jb), jb), jac = 1,      &
                                    & nacin(jb))
                             endif
                         enddo
                         Nxtvd = MIN(Nxtvd, nxcin1(jb))
                     endif
                 endif
             endif
 
!****        Outflow
 
             if(DN_FLOW(jb) .AND. NSTR(jb)>0)then
                 do while (jday>=nxqot1(jb))
                     qstr(1:NSTR(jb), jb) = qstrnx(1:NSTR(jb), jb)
                     qstro(1:NSTR(jb), jb) = qstrnx(1:NSTR(jb), jb)
                     nxqot2(jb) = nxqot1(jb)
                     if(otqf(jb))then
                         read(otq(jb), *)nxqot1(jb),                           &
                            & (qstrnx(js, jb), js = 1, NSTR(jb))
                     else
                         read(otq(jb), '(10F8.0:/(8X,9F8.0))')nxqot1(jb),      &
                            & (qstrnx(js, jb), js = 1, NSTR(jb))
                     endif
                 enddo
                 if(DYNSTRUC(jb)=='      ON')then
                     do while (jday>=nxdyns(jb))
                         do j = 1, njs
                             estr(jjs(j), jb) = nxestrt(jjs(j), jb)
                         enddo
                         if(dynef(jb))then
                             read(odyns(jb), *)nxdyns(jb),                     &
                                & (nxestrt(jjs(j), jb), j = 1, njs)
                         else
                             read(odyns(jb), '(10F8.0:/(8X,9F8.0))')nxdyns(jb),&
                                & (nxestrt(jjs(j), jb), j = 1, njs)
                         endif
                     enddo
                     Nxtvd = MIN(Nxtvd, nxdyns(jb))
                 endif
                 Nxtvd = MIN(Nxtvd, nxqot1(jb))
             endif
 
!****        Distributed tributaries
 
             if(DIST_TRIBS(jb))then
 
!******          Inflow
 
                 do while (jday>=nxqdt1(jb))
                     QDTR(jb) = qdtrnx(jb)
                     qdtro(jb) = qdtrnx(jb)
                     nxqdt2(jb) = nxqdt1(jb)
 
                     if(dtqf(jb))then
                         read(dtq(jb), *)nxqdt1(jb), qdtrnx(jb)
                     else
                         read(dtq(jb), '(2F8.0)')nxqdt1(jb), qdtrnx(jb)
                     endif
 
  !        READ (DTQ(JB),'(2F8.0)') NXQDT1(JB),QDTRNX(JB)
                 enddo
                 Nxtvd = MIN(Nxtvd, nxqdt1(jb))
 
!******          Temperature
 
                 do while (jday>=nxtdt1(jb))
                     TDTR(jb) = tdtrnx(jb)
                     tdtro(jb) = tdtrnx(jb)
                     nxtdt2(jb) = nxtdt1(jb)
 
                     if(dttf(jb))then
                         read(dtt(jb), *)nxtdt1(jb), tdtrnx(jb)
                     else
                         read(dtt(jb), '(2F8.0)')nxtdt1(jb), tdtrnx(jb)
                     endif
 
 !         READ (DTT(JB),'(2F8.0)') NXTDT1(JB),TDTRNX(JB)
                 enddo
                 Nxtvd = MIN(Nxtvd, nxtdt1(jb))
 
!******          Constituent concentrations
 
                 if(dtrib_const(jb))then
                     do while (jday>=nxcdt1(jb))
                         cdtr(DTCN(1:nacdt(jb), jb), jb)                       &
                           & = cdtrnx(DTCN(1:nacdt(jb), jb), jb)
                         cdtro(DTCN(1:nacdt(jb), jb), jb)                      &
                           & = cdtrnx(DTCN(1:nacdt(jb), jb), jb)
                         nxcdt2(jb) = nxcdt1(jb)
                         if(dtcf(jb))then
                             read(dtc(jb), *)nxcdt1(jb),                       &
                                & (cdtrnx(DTCN(jac, jb), jb), jac = 1,         &
                                & nacdt(jb))
                         else
                             read(dtc(jb), '(1000F8.0)')nxcdt1(jb),            &
                                & (cdtrnx(DTCN(jac, jb), jb), jac = 1,         &
                                & nacdt(jb))
                         endif
                     enddo
                     Nxtvd = MIN(Nxtvd, nxcdt1(jb))
                 endif
             endif
 
!****        Precipitation
 
             if(PRECIPITATION(jw))then
                 do while (jday>=nxpr1(jb))
                     PR(jb) = prnx(jb)
                     nxpr2(jb) = nxpr1(jb)
 
                     if(prqf(jb))then
                         read(pre(jb), *)nxpr1(jb), prnx(jb)
                     else
                         read(pre(jb), '(2F8.0)')nxpr1(jb), prnx(jb)
                     endif
 
          !READ (PRE(JB),'(2F8.0)') NXPR1(JB),PRNX(JB)
                 enddo
                 Nxtvd = MIN(Nxtvd, nxpr1(jb))
 
!******          Temperature
 
                 do while (jday>=nxtpr1(jb))
                     TPR(jb) = tprnx(jb)
                     nxtpr2(jb) = nxtpr1(jb)
 
                     if(prtf(jb))then
                         read(prt(jb), *)nxtpr1(jb), tprnx(jb)
                     else
                         read(prt(jb), '(2F8.0)')nxtpr1(jb), tprnx(jb)
                     endif
 
          !READ (PRT(JB),'(2F8.0)') NXTPR1(JB),TPRNX(JB)
                 enddo
                 Nxtvd = MIN(Nxtvd, nxtpr1(jb))
 
!******          Constituent concentrations
 
                 if(precip_const(jb))then
                     do while (jday>=nxcpr1(jb))
                         cpr(PRCN(1:nacpr(jb), jb), jb)                        &
                           & = cprnx(PRCN(1:nacpr(jb), jb), jb)
                         nxcpr2(jb) = nxcpr1(jb)
                         if(prcf(jb))then
                             read(prc(jb), *)nxcpr1(jb),                       &
                                & (cprnx(PRCN(jac, jb), jb), jac = 1, nacpr(jb)&
                                & )
                         else
                             read(prc(jb), '(1000F8.0)')nxcpr1(jb),            &
                                & (cprnx(PRCN(jac, jb), jb), jac = 1, nacpr(jb)&
                                & )
                         endif
                     enddo
                     Nxtvd = MIN(Nxtvd, nxcpr1(jb))
                 endif
             endif
 
!****        Upstream head conditions
 
             if(UH_EXTERNAL(jb))then
 
!******          Elevations
 
                 do while (jday>=nxeuh1(jb))
                     ELUH(jb) = eluhnx(jb)
                     eluho(jb) = eluhnx(jb)
                     nxeuh2(jb) = nxeuh1(jb)
 
                     if(euhf(jb)>0)then
                         read(uhe(jb), *)nxeuh1(jb), eluhnx(jb)
                     else
                         read(uhe(jb), '(2F8.0)')nxeuh1(jb), eluhnx(jb)
                     endif
          !READ (UHE(JB),'(2F8.0)') NXEUH1(JB),ELUHNX(JB)
                 enddo
                 Nxtvd = MIN(Nxtvd, nxeuh1(jb))
 
!******          Temperatures
 
                 do while (jday>=nxtuh1(jb))
                     do k = 2, kmx - 1
                         tuh(k, jb) = tuhnx(k, jb)
                         tuho(k, jb) = tuhnx(k, jb)
                     enddo
                     nxtuh2(jb) = nxtuh1(jb)
 
                     if(tuhf(jb)==1)then
                         read(uht(jb), *)nxtuh1(jb), xx(1)
                         tuhnx(2:KB(US(jb)), jb) = xx(1)
                     elseif(tuhf(jb)==2)then
                         read(uht(jb), *)nxtuh1(jb),                           &
                            & (tuhnx(k, jb), k = 2, KB(US(jb)))
                     else
                         read(uht(jb), '(10F8.0:/(8X,9F8.0))')nxtuh1(jb),      &
                            & (tuhnx(k, jb), k = 2, KB(US(jb)))
                     endif
 
 
          !READ (UHT(JB),'(10F8.0:/(8X,9F8.0))') NXTUH1(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
                 enddo
                 Nxtvd = MIN(Nxtvd, nxtuh1(jb))
 
!******          Constituent concentrations
 
                 if(constituents)then
                     do while (jday>=nxcuh1(jb))
                         do k = 2, kmx - 1
                             cuh(k, CN(1:nac), jb) = cuhnx(k, CN(1:nac), jb)
                             cuho(k, CN(1:nac), jb) = cuhnx(k, CN(1:nac), jb)
                         enddo
                         nxcuh2(jb) = nxcuh1(jb)
 
                         if(cuhf(jb)==1)then
                             read(uhc(jb), *)nxcuh1(jb),                       &
                                & (xx(CN(jac)), jac = 1, nac)
                             do jac = 1, nac
                                 cuhnx(2:KB(US(jb)), CN(jac), jb) = xx(CN(jac))
                             enddo
                         else
 
                             do jac = 1, nac
                                 if(cuhf(jb)==2)then
                                     read(uhc(jb), *)nxcuh1(jb),               &
                                       & (cuhnx(k, CN(jac), jb), k = 2,        &
                                       & KB(US(jb)))
                                 else
                                     read(uhc(jb), '(10F8.0:/(8X,9F8.0))')     &
                                       & nxcuh1(jb),                           &
                                       & (cuhnx(k, CN(jac), jb), k = 2,        &
                                       & KB(US(jb)))
                                 endif
                             enddo
                         endif
            !DO JAC=1,NAC
            !  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH1(JB),(CUHNX(K,CN(JAC),JB),   &
            !                                                    K=2,KB(US(JB)))
            !END DO
                     enddo
                     Nxtvd = MIN(Nxtvd, nxcuh1(jb))
                 endif
             endif
 
!****        Downstream head
 
             if(DH_EXTERNAL(jb))then
 
!******          Elevation
 
                 do while (jday>=nxedh1(jb))
                     ELDH(jb) = eldhnx(jb)
                     eldho(jb) = eldhnx(jb)
                     nxedh2(jb) = nxedh1(jb)
 
                     if(edhf(jb)>0)then
                         read(dhe(jb), *)nxedh1(jb), eldhnx(jb)
                     else
                         read(dhe(jb), '(2F8.0)')nxedh1(jb), eldhnx(jb)
                     endif
          !READ (DHE(JB),'(2F8.0)') NXEDH1(JB),ELDHNX(JB)
                 enddo
                 Nxtvd = MIN(Nxtvd, nxedh1(jb))
 
!******          Temperature
 
                 do while (jday>=nxtdh1(jb))
                     do k = 2, kmx - 1
                         tdh(k, jb) = tdhnx(k, jb)
                         tdho(k, jb) = tdhnx(k, jb)
                     enddo
                     nxtdh2(jb) = nxtdh1(jb)
                     if(tdhf(jb)==1)then
                         read(dht(jb), *)nxtdh1(jb), xx(1)
                         tdhnx(2:KB(DS(jb)), jb) = xx(1)
                     elseif(tdhf(jb)==2)then
                         read(dht(jb), *)nxtdh1(jb),                           &
                            & (tdhnx(k, jb), k = 2, KB(DS(jb)))
                     else
                         read(dht(jb), '(10F8.0:/(8X,9F8.0))')nxtdh1(jb),      &
                            & (tdhnx(k, jb), k = 2, KB(DS(jb)))
                     endif
          !READ (DHT(JB),'(10F8.0:/(8X,9F8.0))') NXTDH1(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
                 enddo
                 Nxtvd = MIN(Nxtvd, nxtdh1(jb))
 
!******          Constituents
 
                 if(constituents)then
                     do while (jday>=nxcdh1(jb))
                         do k = 2, kmx - 1
                             cdh(k, CN(1:nac), jb) = cdhnx(k, CN(1:nac), jb)
                             cdho(k, CN(1:nac), jb) = cdhnx(k, CN(1:nac), jb)
                         enddo
                         nxcdh2(jb) = nxcdh1(jb)
                         if(cdhf(jb)==1)then
                             read(dhc(jb), *)nxcdh1(jb),                       &
                                & (xx(CN(jac)), jac = 1, nac)
                             do jac = 1, nac
                                 cdhnx(2:KB(DS(jb)), CN(jac), jb) = xx(CN(jac))
                             enddo
                         else
                             do jac = 1, nac
                                 if(cdhf(jb)==2)then
                                     read(dhc(jb), *)nxcdh1(jb),               &
                                       & (cdhnx(k, CN(jac), jb), k = 2,        &
                                       & KB(DS(jb)))
                                 else
                                     read(dhc(jb), '(10F8.0:/(8X,9F8.0))')     &
                                       & nxcdh1(jb),                           &
                                       & (cdhnx(k, CN(jac), jb), k = 2,        &
                                       & KB(DS(jb)))
                                 endif
                             enddo
                         endif
            !DO JAC=1,NAC
            !  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH1(JB),(CDHNX(K,CN(JAC),JB),   &
            !                                                    K=2,KB(DS(JB)))
            !END DO
                     enddo
                     Nxtvd = MIN(Nxtvd, nxcdh1(jb))
                 endif
             endif
         enddo
     enddo
 
!    Gate height opening
 
     if(gates)then
         do while (jday>=nxqgt)
             nxqgt2 = nxqgt
             where(dyngtc=='     ZGT')
                 egt = bgtnx
                 egto = bgtnx
                 bgt = 1.0
                 g1gt = 1.0
                 g2gt = 1.0
             elsewhere
                 bgt = bgtnx
                 bgto = bgtnx
             endwhere
             if(gatef)then
                 read(gtq, *)nxqgt, (bgtnx(jg), jg = 1, ngt)
             else
                 read(gtq, '(1000F8.0)')nxqgt, (bgtnx(jg), jg = 1, ngt)
             endif
         enddo
         Nxtvd = MIN(Nxtvd, nxqgt)
     endif
 
!    Pipe reduction factor
 
     if(pipes .AND. iopenpipe==1)then
         do while (jday>=nxqpt)
             bp = bpnx
             read(piped, '(1000F8.0)')nxqpt, (bpnx(j), j = 1, npi)
         enddo
         Nxtvd = MIN(Nxtvd, nxqpt)
     endif
 
!    DYNAMIC PUMPS
 
     if(pumps)then
 
         do j = 1, npu
             if(DYNPUMP(j)=='      ON')then
                 do while (jday>=nxpump(j))
                     EPU(j) = epu2(j)
                     EONPU(j) = eonpu2(j)
                     EOFFPU(j) = eoffpu2(j)
                     QPU(j) = qpu2(j)
                     if(dynpumpf(j))then
                         read(pumpd(j), *)nxpump(j), epu2(j), eonpu2(j),       &
                            & eoffpu2(j), qpu2(j)
                     else
                         read(pumpd(j), '(1000F8.0)')nxpump(j), epu2(j),       &
                            & eonpu2(j), eoffpu2(j), qpu2(j)
                     endif
 !     READ (PUMPD(J),'(1000F8.0)') NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
                 enddo
                 Nxtvd = MIN(Nxtvd, nxpump(j))
             endif
         enddo
     endif
 
!    Dead sea case
 
     do jw = 1, nwb
         if(NO_INFLOW(jw))then
             QIN(BS(jw):BE(jw)) = 0.0
             qino(BS(jw):BE(jw)) = 0.0
             QIND(BS(jw):BE(jw)) = 0.0
             qinnx(BS(jw):BE(jw)) = 0.0
             QDTR(BS(jw):BE(jw)) = 0.0
             qdtro(BS(jw):BE(jw)) = 0.0
             qdtrnx(BS(jw):BE(jw)) = 0.0
             PR(BS(jw):BE(jw)) = 0.0
             prnx(BS(jw):BE(jw)) = 0.0
         endif
         if(NO_OUTFLOW(jw))then
             qstr(:, BS(jw):BE(jw)) = 0.0
             qstro(:, BS(jw):BE(jw)) = 0.0
             qstrnx(:, BS(jw):BE(jw)) = 0.0
         endif
     enddo
     where(no_wind)
         WIND = 0.0
         windo = 0.0
         windnx = 0.0
     endwhere
     where(READ_RADIATION .AND. no_heat)
         SRON = 0.0
         sroo = 0.0
         sronx = 0.0
     endwhere
     if(ANY(NO_INFLOW))then
         qtr = 0.0
         qtro = 0.0
         qtrnx = 0.0
         QWD = 0.0
         qwdo = 0.0
         qwdnx = 0.0
     endif
     return
 
!***********************************************************************************************************************************
!**  I N T E R P O L A T E  I N P U T S                                       
!***********************************************************************************************************************************
 
!    **
     entry INTERPOLATE_INPUTS
 
!    Meteorological/light extinction data
 
     do jw = 1, nwb
         if(INTERP_METEOROLOGY(jw))then
             ratio = (nxmet1(jw) - jday)/(nxmet1(jw) - nxmet2(jw))
             TDEW(jw) = (1.0 - ratio)*tdewnx(jw) + ratio*tdewo(jw)
             WIND(jw) = (1.0 - ratio)*windnx(jw) + ratio*windo(jw)
      ! CONVERT PHIO AND PHINX TO LESS THAN 2*PI     SW 2/13/15
             do while (phio(jw)>2.*pi)
                 phio(jw) = phio(jw) - 2.*pi
             enddo
             do while (phinx(jw)>2.*pi)
                 phinx(jw) = phinx(jw) - 2.*pi
             enddo
             if(phio(jw) - phinx(jw)>pi)then
                 PHI(jw) = (1.0 - ratio)*(phinx(jw) + 2.0*pi) + ratio*phio(jw)
             elseif(phio(jw) - phinx(jw)< - pi)then                ! WX 2/13/15
                 PHI(jw) = (1.0 - ratio)*phinx(jw) + ratio*(phio(jw) + 2.0*pi)
                                                                   ! WX 2/13/15
             else
                 PHI(jw) = (1.0 - ratio)*phinx(jw) + ratio*phio(jw)
             endif
 
      !IF (ABS(PHIO(JW)-PHINX(JW)) > PI) THEN
      !  PHI(JW) = (1.0-RATIO)*(PHINX(JW)+2.0*PI)+RATIO*PHIO(JW)
      !ELSE
      !  PHI(JW) = (1.0-RATIO)*PHINX(JW)+RATIO*PHIO(JW)
      !END IF
             TAIR(jw) = (1.0 - ratio)*tairnx(jw) + ratio*tairo(jw)
             CLOUD(jw) = (1.0 - ratio)*cloudnx(jw) + ratio*cloudo(jw)
             if(READ_RADIATION(jw))SRON(jw) = (1.0 - ratio)*sronx(jw)          &
              & + ratio*sroo(jw)
         endif
         if(READ_EXTINCTION(jw) .AND. INTERP_EXTINCTION(jw))then
                                                               ! 6/30/15 SW
             ratio = (nxext1(jw) - jday)/(nxext1(jw) - nxext2(jw))
             EXH2O(jw) = (1.0 - ratio)*extnx(jw) + ratio*exto(jw)
         endif
     enddo
 
!    Withdrawals
 
     if(nwd>0)then
         qratio = (nxqwd1 - jday)/(nxqwd1 - nxqwd2)
         do jwd = 1, nwd
             if(INTERP_WITHDRAWAL(jwd))QWD(jwd) = (1.0 - qratio)*qwdnx(jwd)    &
              & + qratio*qwdo(jwd)
         enddo
     endif
 
!    Gates  adding interpolation cb 8/13/2010
     if(gates)then
         qratio = (nxqgt - jday)/(nxqgt - nxqgt2)
         do jg = 1, ngt
             if(INTERP_GATE(jg))then
                 if(dyngtc(jg)=='     ZGT')then
                     egt(jg) = (1.0 - qratio)*bgtnx(jg) + qratio*egto(jg)
                 else
                     bgt(jg) = (1.0 - qratio)*bgtnx(jg) + qratio*bgto(jg)
                 endif
             endif
         enddo
     endif
 
!    Tributaries
 
     if(ntr>0)then
         do jt = 1, ntr
             if(INTERP_TRIBS(jt))then
                 qratio = (nxqtr1(jt) - jday)/(nxqtr1(jt) - nxqtr2(jt))
                 tratio = (nxttr1(jt) - jday)/(nxttr1(jt) - nxttr2(jt))
                 if(trib_const(jt))cratio = (nxctr1(jt) - jday)                &
                  & /(nxctr1(jt) - nxctr2(jt))
                 qtr(jt) = (1.0 - qratio)*qtrnx(jt) + qratio*qtro(jt)
                 ttr(jt) = (1.0 - tratio)*ttrnx(jt) + tratio*ttro(jt)
                 ctr(TRCN(1:nactr(jt), jt), jt) = (1.0 - cratio)               &
                   & *ctrnx(TRCN(1:nactr(jt), jt), jt)                         &
                   & + cratio*ctro(TRCN(1:nactr(jt), jt), jt)
             endif
         enddo
     endif
 
!    Branch related inputs
 
     do jb = 1, nbr
 
!**      Inflow
 
         if(UP_FLOW(jb))then
             if(.NOT.INTERNAL_FLOW(jb) .AND. .NOT.DAM_INFLOW(jb))then                                                 !TC 08/03/04 RA 1/13/06
                 if(INTERP_INFLOW(jb))then
                     qratio = (nxqin1(jb) - jday)/(nxqin1(jb) - nxqin2(jb))
                     tratio = (nxtin1(jb) - jday)/(nxtin1(jb) - nxtin2(jb))
                     if(inflow_const(jb))cratio = (nxcin1(jb) - jday)          &
                      & /(nxcin1(jb) - nxcin2(jb))
                     QIND(jb) = (1.0 - qratio)*qinnx(jb) + qratio*qino(jb)
                     TIND(jb) = (1.0 - tratio)*tinnx(jb) + tratio*tino(jb)
                     cind(INCN(1:nacin(jb), jb), jb) = (1.0 - cratio)          &
                       & *cinnx(INCN(1:nacin(jb), jb), jb)                     &
                       & + cratio*cino(INCN(1:nacin(jb), jb), jb)
                 endif
             endif
         endif
 
!**      Outflow
 
         if(DN_FLOW(jb) .AND. NSTR(jb)>0)then
             qratio = (nxqot1(jb) - jday)/(nxqot1(jb) - nxqot2(jb))
             do js = 1, NSTR(jb)
                 if(INTERP_OUTFLOW(js, jb))qstr(js, jb) = (1.0 - qratio)       &
                  & *qstrnx(js, jb) + qratio*qstro(js, jb)
             enddo
         endif
 
!**      Distributed tributaries
 
         if(DIST_TRIBS(jb))then
             if(INTERP_DTRIBS(jb))then
                 qratio = (nxqdt1(jb) - jday)/(nxqdt1(jb) - nxqdt2(jb))
                 tratio = (nxtdt1(jb) - jday)/(nxtdt1(jb) - nxtdt2(jb))
                 if(dtrib_const(jb))cratio = (nxcdt1(jb) - jday)               &
                  & /(nxcdt1(jb) - nxcdt2(jb))
                 QDTR(jb) = (1.0 - qratio)*qdtrnx(jb) + qratio*qdtro(jb)
                 TDTR(jb) = (1.0 - tratio)*tdtrnx(jb) + tratio*tdtro(jb)
                 cdtr(DTCN(1:nacdt(jb), jb), jb) = (1.0 - cratio)              &
                   & *cdtrnx(DTCN(1:nacdt(jb), jb), jb)                        &
                   & + cratio*cdtro(DTCN(1:nacdt(jb), jb), jb)
             endif
         endif
 
!**      Upstream head
 
         if(UH_EXTERNAL(jb))then
             if(INTERP_HEAD(jb))then
                 hratio = (nxeuh1(jb) - jday)/(nxeuh1(jb) - nxeuh2(jb))
                 tratio = (nxtuh1(jb) - jday)/(nxtuh1(jb) - nxtuh2(jb))
                 if(constituents)cratio = (nxcuh1(jb) - jday)                  &
                  & /(nxcuh1(jb) - nxcuh2(jb))
                 ELUH(jb) = (1.0 - hratio)*eluhnx(jb) + hratio*eluho(jb)
                 do k = 2, kmx - 1
                     tuh(k, jb) = (1.0 - tratio)*tuhnx(k, jb)                  &
                                & + tratio*tuho(k, jb)
                     cuh(k, CN(1:nac), jb) = (1.0 - cratio)                    &
                       & *cuhnx(k, CN(1:nac), jb)                              &
                       & + cratio*cuho(k, CN(1:nac), jb)
                 enddo
             endif
         endif
 
!**      Downstream head
 
         if(DH_EXTERNAL(jb))then
             if(INTERP_HEAD(jb))then
                 hratio = (nxedh1(jb) - jday)/(nxedh1(jb) - nxedh2(jb))
                 tratio = (nxtdh1(jb) - jday)/(nxtdh1(jb) - nxtdh2(jb))
                 if(constituents)cratio = (nxcdh1(jb) - jday)                  &
                  & /(nxcdh1(jb) - nxcdh2(jb))
                 ELDH(jb) = (1.0 - hratio)*eldhnx(jb) + hratio*eldho(jb)
                 do k = 2, kmx - 1
                     tdh(k, jb) = (1.0 - tratio)*tdhnx(k, jb)                  &
                                & + tratio*tdho(k, jb)
                     cdh(k, CN(1:nac), jb) = (1.0 - cratio)                    &
                       & *cdhnx(k, CN(1:nac), jb)                              &
                       & + cratio*cdho(k, CN(1:nac), jb)
                 enddo
             endif
         endif
     enddo
     return
     entry DEALLOCATE_TIME_VARYING_DATA
     deallocate(nxqtr1, nxttr1, nxctr1, nxqin1, nxtin1, nxcin1, nxqdt1, nxtdt1,&
              & nxcdt1, nxpr1, nxtpr1, nxcpr1, nxeuh1, nxtuh1)
     deallocate(nxcuh1, nxedh1, nxtdh1, nxcdh1, nxqot1, nxmet1, nxqtr2, nxttr2,&
              & nxctr2, nxqin2, nxtin2, nxcin2, nxqdt2, nxtdt2)
     deallocate(nxcdt2, nxpr2, nxtpr2, nxcpr2, nxeuh2, nxtuh2, nxcuh2, nxedh2, &
              & nxtdh2, nxcdh2, nxqot2, nxmet2, wscnx, dynpumpf)
     deallocate(qdtro, tdtro, eluho, eldho, qwdo, qtro, ttro, qino, tino,      &
              & qdtrnx, tdtrnx, prnx, tprnx, eluhnx)
     deallocate(eldhnx, qwdnx, qtrnx, ttrnx, qinnx, tinnx, sroo, tairo, tdewo, &
              & cloudo, phio, windo, tairnx, bgtnx, bpnx)
     deallocate(tdewnx, cloudnx, phinx, windnx, sronx, trq, trt, trc, inq, dtq,&
              & pre, uhe, dhe, inft)
     deallocate(dtt, prt, uht, dht, inc, dtc, prc, uhc, dhc, otq, met, ext,    &
              & extnx, exto, extf)
     deallocate(nxext1, nxext2, ctro, cino, qouto, cdtro, tuho, tdho, qstro,   &
              & ctrnx, cinnx, qoutnx, cdtrnx, cprnx)
     deallocate(tuhnx, tdhnx, qstrnx, cuho, cdho, cuhnx, cdhnx, incf, dtcf,    &
              & prcf, metf, odyns, nxdyns, nxestrt, dynef, jjs, prtf, prqf)
     deallocate(euhf, tuhf, cuhf, edhf, tdhf, cdhf, xx)
     deallocate(inflow_const, trib_const, dtrib_const, precip_const, pumpd,    &
              & nxpump, epu2, eonpu2, eoffpu2, qpu2, otqf, trcf, trqf, trtf,   &
              & dttf, dtqf, inqf, intf)
     end subroutine TIME_VARYING_DATA
