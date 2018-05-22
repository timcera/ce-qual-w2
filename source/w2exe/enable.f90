!*==enable.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine ENABLE(Dlg)
     use DFWIN, RENAMED => DLT
     use DFLOGM
     use MSCLIB
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     type(DIALOG) :: Dlg
!
! Local variables
!
     real :: DFLOGM, DFWIN
     real :: DLGSET
     real :: dlg_enable, high, highest
     integer :: idle, low, lowest, normal, result
     integer :: MSCLIB
!
!*** End of declarations rewritten by SPAG
!
                                                                                                   !Rename DLT in DFWIN
 
                                                                                                   ! 9/27/07 SW
     result = DLGSET(Dlg, idle, .TRUE., dlg_enable)                                                !Enable 'Idle'    button
     result = DLGSET(Dlg, highest, .TRUE., dlg_enable)                                             !Enable 'Highest' button
     result = DLGSET(Dlg, high, .TRUE., dlg_enable)                                                !Enable 'High'    button
     result = DLGSET(Dlg, normal, .TRUE., dlg_enable)                                              !Enable 'Normal'  button
     result = DLGSET(Dlg, low, .TRUE., dlg_enable)                                                 !Enable 'Low'     button
     result = DLGSET(Dlg, lowest, .TRUE., dlg_enable)                                              !Enable 'Lowest'  button
     end subroutine ENABLE
