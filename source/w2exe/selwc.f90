!*==selwc.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module SELWC
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND), allocatable, dimension(:, :, :) :: cavg, cdavg
     real(R8KIND), allocatable, dimension(:, :) :: cavgw, cdavgw, qstr, qsw
     character(8), allocatable, dimension(:) :: dynstruc
     real, allocatable, dimension(:, :) :: estr, wstr
     real, allocatable, dimension(:) :: ewd
     integer, allocatable, dimension(:, :) :: kbsw, kout, ktsw
     integer, allocatable, dimension(:) :: kbw, kbwd, ktw, ktwd, nout, nstr
     real(R8KIND), allocatable, dimension(:) :: qnew, vnorm
     real(R8KIND), allocatable, dimension(:, :) :: tavg
     real(R8KIND), allocatable, dimension(:) :: tavgw
!
!*** End of declarations rewritten by SPAG
!
  !REAL,                  ALLOCATABLE, DIMENSION(:)    :: EWD, TAVGW                   ! cb 1/16/13
  !REAL,                  ALLOCATABLE, DIMENSION(:,:)  :: ESTR,   WSTR, TAVG            ! SW Selective 7/30/09
                                                                                       ! cb 1/16/13
                                                                               ! cb 1/16/13
     end module SELWC                                                           ! cb 1/16/13
