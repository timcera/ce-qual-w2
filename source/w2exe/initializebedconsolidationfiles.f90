!*==initializebedconsolidationfiles.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     subroutine INITIALIZEBEDCONSOLIDATIONFILES(Tempfilnum, Tempfilname)
 
     use MAIN
     use GLOBAL
     use SCREENC
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     character(256) :: Tempfilname
     integer(4) :: Tempfilnum
     intent (in) Tempfilname, Tempfilnum
!
! Local variables
!
     character(256) :: messagetemp
     logical :: skiploop
!
!*** End of declarations rewritten by SPAG
!
 
 
    !Open File
     open(Tempfilnum, file = Tempfilname(1:LEN_TRIM(Tempfilname) - 1))
!
!!Read Header
     skiploop = .FALSE.
     do while (.NOT.skiploop)
         read(Tempfilnum, '(a)')messagetemp
         if(INDEX(messagetemp, "$")==0)skiploop = .TRUE.
     enddo
!
     end subroutine INITIALIZEBEDCONSOLIDATIONFILES
