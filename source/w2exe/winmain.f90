!*==winmain.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
!***********************************************************************************************************************************
!**                                                 F U N C T I O N   W I N M A I N                                               **
!***********************************************************************************************************************************
 
     function WINMAIN(Hinstance, Hprevinstance, Lpszcmdline, Ncmdshow)
     use DFLIB
     use DFWIN, RENAMED => DLT
     use DFLOGM
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     integer(4) :: Hinstance, Hprevinstance, Lpszcmdline, Ncmdshow
     integer(I4KIND) :: WINMAIN
!
! Local variables
!
     real :: DFLIB, DFLOGM, DFWIN
!
!*** End of declarations rewritten by SPAG
!
  !DEC$ IF DEFINED(_X86_)
  !DEC$ ATTRIBUTES STDCALL, ALIAS : '_WinMain@16' :: WINMAIN
  !DEC$ ELSE
  !DEC$ ATTRIBUTES STDCALL, ALIAS: 'WinMain':: WinMain
  !DEC$ ENDIF
  !DEC$ IF DEFINED(_X86_)
  !DEC$ ATTRIBUTES STDCALL, ALIAS : '_WinMain@16' :: WINMAIN
  !DEC$ ELSE
  !DEC$ ATTRIBUTES STDCALL, ALIAS: 'WinMain':: WinMain
  !DEC$ ENDIF
  !DEC$ IF DEFINED(_X86_)
  !DEC$ ATTRIBUTES STDCALL, ALIAS : '_WinMain@16' :: WINMAIN
  !DEC$ ELSE
  !DEC$ ATTRIBUTES STDCALL, ALIAS: 'WinMain':: WinMain
  !DEC$ ENDIF
 
     call W2_DIALOG
     WINMAIN = 0
     end function WINMAIN
