!*==rstart.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module RSTART
     use PREC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(r8), allocatable, dimension(:, :) :: cmbrt, savh2, savhr, saz, su,   &
           & sw, tssdh2, tssuh2
     integer, allocatable, dimension(:) :: cpldp, flxdp, nsprf, prfdp, scrdp,  &
          & snpdp, sprdp, vpldp
     real(r8), allocatable, dimension(:, :, :) :: cssdh2, cssuh2
     integer :: cuf, dltdp, rsodp, tsrdp, wdodp
     real(r8) :: curmax, dltff, dltmaxx, dlts, eltm, nxtmrs, nxtmts, nxtmwd,   &
               & nxtmwd_sec
     real(r8), allocatable, dimension(:) :: ebri, eltmf, esbr, etbr, icebank,  &
           & sbkt, sz, tssb, tssdh, tssdt, tssev, tssice, tssin, tssout, tsspr,&
           & tsss, tsstr, tssuh, tsswd, voldh, voldt, volev, volice, volin,    &
           & volout, volpr, volsbr, volsr, voltbr, voltr, voltrb, voluh, volwd
     real, allocatable, dimension(:) :: nxtmcp, nxtmfl, nxtmpr, nxtmsc, nxtmsn,&
                                      & nxtmsp, nxtmvp
     real :: PREC
     real :: r8
     integer :: rso = 31
!
!*** End of declarations rewritten by SPAG
!
  !REAL                                               :: NXTMRS, NXTMWD, NXTMTS
                                                                                          ! SW 7/13/2010
     end module RSTART                                                                         ! cb 4/6/17
