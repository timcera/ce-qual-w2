!*==kinetic.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module KINETIC
     use PREC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real, allocatable, dimension(:) :: ac, achla, ae, ag, ahsn, ahsp, ahssi,  &
                                      & ak1, ak2, ak3, ak4, am, an, anpr, ap,  &
                                      & apom, ar, as, asat, asi, at1, at2, at3,&
                                      & at4, bodc, bodn, bodp, cadk, caq10,    &
                                      & cas, cbiom, cbods, cg0dk, cg1dk, cgcs, &
                                      & cgklf, cgldk, cgq10, cgs, co2r, dsir,  &
                                      & eb, ec, echla, ee, eg, ehs, ehsn, ehsp,&
                                      & ehssi, ek1, ek2, ek3, ek4, em, en,     &
                                      & enpr, ep, epom, er, esat, esi, et1,    &
                                      & et2, et3, et4, fer, fes, fno3sed,      &
                                      & fract, kbod, ldomdk, lpomdk, lrddk,    &
                                      & lrpdk, nbiom, nfluxin, nh4dk, nh4k1,   &
                                      & nh4k2
     real, pointer, dimension(:, :, :) :: aer, agr, amr, arr, asr, ebr, eer,   &
       & egr, emr, err
     real(r8), pointer, dimension(:, :, :) :: alg, ass, cbod, cbodn, cbodnss,  &
           & cbodp, cbodpss, cbodss, cg, cgss, ss, ssss
     real(r8), pointer, dimension(:, :) :: alk, alkss, cass, col, colss, doss, &
           & dsi, dsiss, fe, fess, ldom, ldomn, ldomnss, ldomp, ldompss,       &
           & ldomss, lpom, lpomn, lpomnss, lpomp, lpompss, lpomss, nh4, nh4ss, &
           & no3, no3ss, o2, po4, po4ss, psi, psiss, rdom, rdomn, rdomnss,     &
           & rdomp, rdompss, rdomss, rpom, rpomn, rpomnss, rpomp, rpompss,     &
           & rpomss, tds, tic, ticss
     logical :: ammonia_buffering, noncon_alkalinity, om_buffering,            &
              & phosphate_buffering, ph_buffering, pom_buffering
     integer, allocatable, dimension(:) :: aneqn, eneqn, naf, neqn
     real, pointer, dimension(:, :) :: apr, atot, cboddk, cbodns, cbodnsn,     &
                                     & cbodnsp, cbodu, ch4d, ch4reaer, chla,   &
                                     & co2, co2reaer, co3, doae, doap, doar,   &
                                     & dobod, doc, doch4, dodom, doep, doer,   &
                                     & dofe2, doh2s, domn2, don, donit, doom,  &
                                     & dop, dopom, dosed, dosedia, dosod,      &
                                     & dsiag, dsid, dsieg, dsis, dsisd, dsisr, &
                                     & fe2d, fens, fesr, h2sd, h2sreaer,       &
                                     & h2ssedd, hco3, ldomap, ldomd, ldomep,   &
                                     & ldomnap, ldomnep, ldompap, ldompep,     &
                                     & lpomap, lpomd, lpomep, lpomepc, lpomepn,&
                                     & lpomepp, lpomnap, lpomnns, lpomns,      &
                                     & lpompap, lpompns, lrdomd, lrpomd, mn2d, &
                                     & nh4ag, nh4ap, nh4ar, nh4d, nh4dom,      &
                                     & nh4eg, nh4ep, nh4er, nh4om, nh4pom,     &
                                     & nh4sd, nh4sr, no3ag
     real(r8), allocatable, dimension(:) :: beta, cz, exa, exh2o, exom, exss,  &
           & qc, qerr, wind10
     character(8), allocatable, dimension(:) :: cac, reaerc
     real, allocatable, dimension(:, :, :) :: cbodd, epc, epd, epm
     real, allocatable, dimension(:, :) :: do1, do2, do3, fpfe, fpss, gamma,   &
       & ldomnmp, ldompmp, lpomnmp, lpompmp, lpzooinn, lpzooinp, lpzoooutn,    &
       & lpzoooutp, orgnld, orgnlp, orgnrd, orgnrp, orgpld, orgplp, orgprd,    &
       & orgprp, rpomnmp, rpompmp, sdkv, sed, sed1, sed1ic, sed2, sed2ic, sedc,&
       & seddktot, sedn, sedninflux, sedp, sedpinflux, sedvpc, sedvpn, sedvpp
     real :: kdo, r8
     integer, allocatable, dimension(:, :) :: kfcn
     character(10), allocatable, dimension(:, :) :: lfpr
     integer :: nag, nagi, nldomn, nldomp, nlpomn, nlpomp, nrdomn, nrdomp,     &
              & nrpomn, nrpomp
     character(8) :: ncalkc, nh4bufc, ombufc, omtype, phbufc, po4bufc, pombufc
     real, allocatable, dimension(:) :: nh4r, nh4t1, nh4t2, no3dk, no3k1,      &
                                      & no3k2, no3s, no3t1, no3t2, o2ag, o2ar, &
                                      & o2eg, o2er, o2nh4, o2om, omk1, omk2,   &
                                      & omt1, omt2, orgc, orgn, orgp, orgsi,   &
                                      & partp, partsi, pbiom, pfluxin, pk, pki,&
                                      & pksd, po4r, poms, psidk, psis, rbod,   &
                                      & rcoef1, rcoef2, rcoef3, rcoef4, rdomdk,&
                                      & reaer, rpomdk, sden, sdeni, sdk, sdk1, &
                                      & sdk2, sedb, seds, sod, sodk1, sodk2,   &
                                      & sodt1, sodt2, sss, taucr, tbod
     real, pointer, dimension(:, :) :: no3d, no3eg, no3sed, o2dg, ph, po4ag,   &
                                     & po4ap, po4ar, po4dom, po4eg, po4ep,     &
                                     & po4er, po4ns, po4om, po4pom, po4sd,     &
                                     & po4sr, poc, pon, pop, psiam, psid,      &
                                     & psins, rdomd, rpomd, rpomnns, rpomns,   &
                                     & rpompns, sdinc, sdinfeooh, sdinmno2,    &
                                     & sdinn, sdinp, sedas, sedasc, sedasn,    &
                                     & sedasp, sedbr, sedbrc, sedbrn, sedbrp,  &
                                     & sedcb, sedcbc, sedcbn, sedcbp, sedd,    &
                                     & sedd1, sedd2, seddc, seddn, seddp,      &
                                     & sedno3, sedns, sednsc, sednsn, sednsp,  &
                                     & sedoms, sedomsc, sedomsn, sedomsp, sodd,&
                                     & sssi, ssso, tdg, ticap, ticep, tiss,    &
                                     & tkn, tn, toc, ton, top, totss, tp
     real :: PREC
     logical, allocatable, dimension(:, :) :: sdfirstadd
     logical, allocatable, dimension(:) :: sediment_resuspension
!
!*** End of declarations rewritten by SPAG
!
                                                                       ! enhanced pH buffering
                                                                   ! SW 10/17/15
                                                                      ! Amaila
                                                                                                            !CB 11/30/06
                                                                              ! CB 6/6/10
                                                                                  ! CB 6/6/10
                                                                                   ! cb 6/6/10
                                                                                                         !LCJ 2/26/15 SW 10/16/15
                                                                     ! amaila
                                                                                                 !  , SSFLOC   cb 11/27/06   SR 04/21/13
  !REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SED1,sed2   ! Amaila
                                                                           ! Amaila, cb 6/8/17
                                                                                  ! cb 6/17/17
                                                                                                     ! SW 4/2016
  ! enhanced pH buffering start
  ! enhanced pH buffering end
     contains                                                                          !, FLOCEQN   ! SR 04/21/13
     function SATO(T, Sal, P, Salt_water)
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(R8) :: P, Sal, T
     logical :: Salt_water
     real :: SATO
     intent (in) P, Sal, Salt_water, T
!
!*** End of declarations rewritten by SPAG
!
     SATO = EXP(7.7117 - 1.31403*(LOG(T + 45.93)))*P
     if(Salt_water)then
         SATO = EXP(LOG(SATO) - Sal*(1.7674E-2 - 1.0754E1/(T + 273.15)         &
              & + 2.1407E3/(T + 273.15)**2))                                                          ! SAL is in ppt
     elseif(Sal>100.)then
         SATO = EXP(LOG(SATO) - (Sal/1000.)                                    &
              & *(1.7674E-2 - 1.0754E1/(T + 273.15) + 2.1407E3/(T + 273.15)**2)&
              & )                                                                                     ! SAL is in mg/l
     endif
     end function SATO
     function FR(Tt, Tt1, Tt2, Sk1, Sk2)
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real :: Sk1, Sk2, Tt1, Tt2
     real(R8) :: Tt
     real :: FR
     intent (in) Sk1, Sk2, Tt, Tt1, Tt2
!
!*** End of declarations rewritten by SPAG
!
     FR = Sk1*EXP(LOG(Sk2*(1.0 - Sk1)/(Sk1*(1.0-Sk2)))/(Tt2 - Tt1)*(Tt - Tt1))
     end function FR
     function FF(Tt, Tt3, Tt4, Sk3, Sk4)
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real :: Sk3, Sk4, Tt3, Tt4
     real(R8) :: Tt
     real :: FF
     intent (in) Sk3, Sk4, Tt, Tt3, Tt4
!
!*** End of declarations rewritten by SPAG
!
     FF = Sk4*EXP(LOG(Sk3*(1.0 - Sk4)/(Sk4*(1.0-Sk3)))/(Tt4 - Tt3)*(Tt4 - Tt))
     end function FF
     end module KINETIC
