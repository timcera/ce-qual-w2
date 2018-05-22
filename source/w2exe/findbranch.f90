!*==findbranch.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine FINDBRANCH(Seg)
 
 
     use FISHY
     use GLOBAL
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real :: Seg
     intent (in) Seg
!
! Local variables
!
     integer :: iseg, j
!
!*** End of declarations rewritten by SPAG
!
 
 
     iseg = INT(Seg)
     do j = 1, nbr
         if(iseg>=US(j) .AND. iseg<=DS(j))exit
     enddo
     fnbp = j
 
     end subroutine FINDBRANCH
