!*==geomc.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module GEOMC
     use PREC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(r8), allocatable, dimension(:) :: alpha, bkt, cosa, dlx, dlxr, elws, &
           & sina, sinac, slope, slopec, z
     real(r8), allocatable, dimension(:, :) :: avh1, avh2, avhr, b, bb, bh,    &
           & bh1, bh2, bhr, bhr1, bhr2, bi, bnew, br, depthb, depthm, el,      &
           & fetchd, fetchu, h, h1, h2
     integer, allocatable, dimension(:) :: jbdh, jbuh, jwdh, jwuh
     real :: PREC
     real :: r8
!
!*** End of declarations rewritten by SPAG
!
     end module GEOMC                                                                                                              ! SW 1/23/06
