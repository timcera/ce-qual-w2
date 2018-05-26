!*==tvdc.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module TVDC
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(R8KIND), allocatable, dimension(:, :, :) :: cdh, cuh
     character(72), allocatable, dimension(:) :: cdhfn, cdtfn, cinfn, cprfn,   &
                & ctrfn, cuhfn, edhfn, euhfn, extfn, metfn, prefn, qdtfn,      &
                & qinfn, qotfn, qtrfn, tdhfn, tdtfn, tinfn, tprfn, ttrfn, tuhfn
     real(R8KIND), allocatable, dimension(:, :) :: cdtr, cin, cind, cpr, ctr, qout,&
           & tdh, tuh
     real(R8KIND), allocatable, dimension(:) :: cloud, eldh, eluh, phi, pr, qdtr,  &
           & qin, qind, qsum, qtr, qwd, sron, tair, tdew, tdtr, tin, tind,     &
           & tout, tpr, ttr, twdo
     integer, allocatable, dimension(:) :: cn, nacd, nacdt, nacin, nacpr, nactr
     logical :: constituents
     integer, allocatable, dimension(:, :) :: dtcn, incn, prcn, trcn
     integer :: nac, nopen
     character(72) :: qgtfn, qwdfn, shdfn, wscfn
!
!*** End of declarations rewritten by SPAG
!
     end module TVDC
