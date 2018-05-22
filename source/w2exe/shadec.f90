!*==shadec.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module SHADEC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
     integer, parameter :: IANG = 18
     real, parameter :: GAMA = (3.1415926*2.)/REAL(IANG)
!
! Local variables
!
     real, allocatable, dimension(:) :: a00, cllb, clrb, decl, hh, shadei,     &
                                      & srfjd1, srfjd2, srlb1, srlb2, srrb1,   &
                                      & srrb2, ttlb, ttrb
     real, dimension(IANG), save :: ang
     logical, allocatable, dimension(:) :: dynamic_shade
     real, allocatable, dimension(:, :) :: topo
!
!*** End of declarations rewritten by SPAG
!
                                                                                                           ! SW 10/17/05
                                                                                                           ! SW 10/17/05
                                                                                                           ! SW 10/17/05
     data ang/0.00000, 0.34907, 0.69813, 1.04720, 1.39626, 1.74533, 2.09440,   &
        & 2.44346, 2.79253, 3.14159, 3.49066, 3.83972, 4.18879, 4.53786,       &
        & 4.88692, 5.23599, 5.58505, 5.93412/                                                              ! SW 10/17/05
     end module SHADEC
