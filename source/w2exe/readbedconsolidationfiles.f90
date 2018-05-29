!*==readbedconsolidationfiles.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine READBEDCONSOLIDATIONFILES(Tempfilnum, Consolidratetemp)
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
     real(8) :: Consolidratetemp
     integer(4) :: Tempfilnum
     intent (out) Consolidratetemp
!
! Local variables
!
     real(8) :: consolidratetemp1, consolidratetemp2, factorinterp, timejd1,   &
              & timejd2
     real :: EOF
     logical :: skiploop
!
!*** End of declarations rewritten by SPAG
!
 
 
    !Read Data
     skiploop = .FALSE.
     do while (.NOT.skiploop .OR. EOF(Tempfilnum))
!
         read(Tempfilnum, '(2F8.0)')timejd1, consolidratetemp1
         read(Tempfilnum, '(2F8.0)')timejd2, consolidratetemp2
 
         if(jday>=timejd1 .AND. jday<=timejd2)skiploop = .TRUE.
         backspace(Tempfilnum)
 
     enddo
!
     backspace(Tempfilnum)
!
     factorinterp = (jday - timejd1)/(timejd2 - timejd1)
 
     Consolidratetemp = consolidratetemp1*(1 - factorinterp)                   &
                      & + consolidratetemp2*factorinterp
!
     end subroutine READBEDCONSOLIDATIONFILES
