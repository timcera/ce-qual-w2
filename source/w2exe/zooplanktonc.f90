!*==zooplanktonc.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module ZOOPLANKTONC
     use PREC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real, allocatable, dimension(:, :, :, :) :: agz, zgz
     real, allocatable, dimension(:, :, :) :: agzt, tgraze, zmt, zmu, zoorm,   &
       & zoormf, zoormr, zrt
     real, allocatable, dimension(:, :) :: dozr, lpzooin, lpzooout, nh4zr,     &
       & po4zr, prefa, prefz, ticzr
     real, allocatable, dimension(:) :: exz, o2zr, prefp, zc, zeff, zg, zk1,   &
                                      & zk2, zk3, zk4, zm, zn, zoomin, zp, zr, &
                                      & zs2p, zt1, zt2, zt3, zt4
     real :: PREC
     real :: r8
     real(r8), pointer, dimension(:, :, :) :: zoo, zooss
     logical :: zooplankton_calc
!
!*** End of declarations rewritten by SPAG
!
                                                                     ! OMNIVOROUS ZOOPLANKTON
     end module ZOOPLANKTONC                                     ! OMNIVOROUS ZOOPLANKTON
