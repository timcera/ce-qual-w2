!*==structures.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module STRUCTURES
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND), allocatable, dimension(:) :: a1gt, a1sp, a2gt, a2sp, b1gt, b1sp,&
           & b2gt, b2sp, bgt, bp, dlxpi, dtp, dtps, edpi, egt, egt2, endpu,    &
           & eoffpu, eonpu, epu, esp, eupi, fminpi, fpi, g1gt, g2gt, gta1,     &
           & gta2, gtb1, gtb2, qgt, qold, qolds, qpi, qpu, qsp, strtpu, vmax,  &
           & wpi
     logical, allocatable, dimension(:) :: begin, lateral_gate, lateral_pipe,  &
          & lateral_pump, lateral_spillway, wlflag
     real, allocatable, dimension(:) :: bgto, egto, ekbwr, ektwr
     real(R8KIND) :: clen, closs, dia, dnie, fman, upie
     character(8), allocatable, dimension(:) :: dyngtc, dynpipe, dynpump, gtic,&
               & latgtc, latpic, latpuc, latspc
     real(R8KIND), save :: eps2, omega, thr
     character(8) :: gt2char
     integer, allocatable, dimension(:) :: idgt, idpi, idpu, idsp, iugt, iupi, &
          & iupu, iusp, iwr, jbdgt, jbdpi, jbdpu, jbdsp, jbugt, jbupi, jbupu,  &
          & jbusp, jwdgt, jwdpi, jwdpu, jwdsp, jwugt, jwupi, jwupu, jwusp,     &
          & kbpu, kbwr, ktpu, ktwr
     integer, save :: nc, nn, nnpipe
     real(R8KIND), allocatable, dimension(:, :) :: vs, vss, vst, vsts, ys, yss,    &
           & yst, ysts
!
!*** End of declarations rewritten by SPAG
!
                                                                          ! SW 3/18/16
                                                                                                                                         ! SW 5/10/10
                                                                        ! CB/8/13/ 2010
     data thr/0.01D0/, omega/0.8D0/, eps2/0.0001D0/                     ! cb/8/13/ 2010
     data nn/19/, nnpipe/19/, nc/7/
     end module STRUCTURES
