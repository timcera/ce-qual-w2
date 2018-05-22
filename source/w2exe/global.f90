!*==global.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module GLOBAL
     use PREC
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
     real(r8), parameter :: DAY = 86400.0D0, NONZERO = 1.0D-20, REFL = 0.94D0, &
                          & FRAZDZ = 0.14D0, DZMIN = 1.4D-7, AZMIN = 1.4D-6,   &
                          & DZMAX = 1.0D3, RHOW = 1000.0D0
!
! Local variables
!
     real(r8), pointer, dimension(:, :) :: admx, admz, az, dltlim, dm, grav,   &
           & hdg, hpg, rho, sb, st, t2, u, vsh, w
     real, allocatable, target, dimension(:, :, :, :) :: af, ef
     real, allocatable, dimension(:, :, :) :: allim, anlim, aplim, aslim,      &
       & ellim, enlim, eplim, eslim, kfs
     real(r8), allocatable, dimension(:) :: azmax, cdmult, cmult, elkt, hmult, &
           & iceth, palt, wind2, z0
     integer, allocatable, dimension(:) :: be, bs, cdhs, cus, dhs, dqb, ds,    &
          & jbdn, kb, kbmin, kti, ktwb, nbodc, nbodn, nbodp, skti, uhs, uqb, us
     real(r8) :: betabr, current, dlt, dltmin, dlttvd, hmax2, start
     logical, allocatable, dimension(:) :: br_inactive, ice, ice_calc,         &
          & layerchange
     real(r8), allocatable, target, dimension(:, :, :) :: c1, c1s, c2, cssb,   &
           & cssk, hyd
     character(10) :: cctime
     real, allocatable, target, dimension(:, :, :) :: cd, kf
     character(12) :: cdate
     integer, pointer, dimension(:) :: cpl, flx, flx2, prf, snp, spr, vpl
     real(r8), allocatable, dimension(:, :), save :: curz1, curz2, curz3, ratz
     real(r8), external :: DENSITY
     real(r8), save :: g, pi
     integer :: i, id, imx, iu, jb, jc, jjb, jw, jz, kmx, kt, nal, nbod, nbr,  &
              & nct, nep, nept, ngce, ngcs, ngt, niktsr, nmc, nmct, nod, npi,  &
              & npu, nsp, nss, nst, ntr, nunit, nwb, nwd, nwdo, nzooe, nzoos,  &
              & nzp, nzpt
     character(180) :: moddir
     integer, save :: ndc, nfl, nhy, w2err, wrn
     integer, allocatable, target, dimension(:, :) :: opt
     real :: PREC
     real(r8), allocatable, dimension(:, :) :: qdh1, qss, quh1, uxbr, uybr,    &
           & vol, voldh2, voluh2
     real :: r8
     character(72) :: rsifn
     real(r8), allocatable, target, dimension(:, :) :: t1, tss
     real, allocatable, dimension(:) :: tndtrib, tnin, tnout, tnpr, tntrib,    &
                                      & tnwd, tn_sedsod_nh4, tpdtrib, tpin,    &
                                      & tpout, tppr, tptrib, tpwd,             &
                                      & tp_sedsod_po4
     real(R4KIND) :: w2ver = 4.1
!
!*** End of declarations rewritten by SPAG
!
                                                                                                                                                              !TP_SEDBURIAL,TN_SEDBURIAL,
                                                                                                           ! number of zooplankton groups, CONSTIUENT NUMBER FOR ZOOPLANKTON, START AND END
                                                                     ! CEMA
                                                                                   ! CB 6/6/10
                                                                                   ! CURRENT WORKING DIRECTORY
                                                                                         ! SW 5/15/06
     data ndc/23/, nhy/15/, nfl/139/                                             ! CEMA and Amaila -altered 'NFL'
     data g/9.81D0/, pi/3.14159265359D0/
     data wrn/32/, w2err/33/
     end module GLOBAL
