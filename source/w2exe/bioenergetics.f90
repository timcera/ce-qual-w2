!*==bioenergetics.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module BIOENERGETICS
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     character(8), dimension(20) :: bhead
     character(8) :: bioc
     real, allocatable, dimension(:) :: biod, biof, volroos
     integer, allocatable, dimension(:) :: biodp, bioexpfn, ibio, weightnum
     logical :: bioexp, fishbio
     character(72) :: biofn, weightfn
     real, allocatable, dimension(:, :) :: c2w
     real, allocatable, dimension(:, :, :) :: c2zoo
     integer, save :: fishbiofn
     real(R8KIND) :: gammab, nxbio, nxtbio
     integer :: klim, nbio, nibio
!
!*** End of declarations rewritten by SPAG
!
     data fishbiofn/9502/
     end module BIOENERGETICS
