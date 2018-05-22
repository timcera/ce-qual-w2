!*==w2_dialog.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                                P R O G R A M   W 2   D I A L O G                                              **
!***********************************************************************************************************************************
 
     subroutine W2_DIALOG
     use DFLOGM
     use MSCLIB
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real :: DFLOGM
     type(DIALOG) :: dlg
     real :: DLGINIT, DLGMODAL, DLGSET, DLGSETSUB, RUN_W2
     real :: dlg_enable, high, highest, output_dialog, restart, run, status,   &
           & stop_execution
     integer :: idle, idok, low, lowest, normal, result
     integer :: MSCLIB
     logical :: restarted
     logical :: rso_exists = .FALSE.
!
!*** End of declarations rewritten by SPAG
!
 
     restarted = .FALSE.
     open(1, file = 'rso.opt', status = 'OLD', iostat = result)
     close(1)                                                                                     !Does restart file exist?
     rso_exists = result==0                                                                       !Restart file exists
     result = DLGINIT(output_dialog, dlg)                                                         !Initialize dialog box
     result = DLGSET(dlg, status, 'Pending execution')                                            !Display execution status
     result = DLGSET(dlg, restart, rso_exists, dlg_enable)                                        !Enable  'Restart' button
     result = DLGSET(dlg, run, .FALSE., dlg_enable)
     result = DLGSET(dlg, stop_execution, .FALSE., dlg_enable)                                    !Disable 'Stop'    button
     result = DLGSETSUB(dlg, output_dialog, RUN_W2)
     result = DLGSETSUB(dlg, run, RUN_W2)                                                         !Set 'Run'     callback
     result = DLGSETSUB(dlg, stop_execution, RUN_W2)                                              !Set 'Stop'    callback
     result = DLGSETSUB(dlg, restart, RUN_W2)                                                     !Set 'Restart' callback
     result = DLGSETSUB(dlg, highest, RUN_W2)                                                     !Set 'Highest' callback
     result = DLGSETSUB(dlg, high, RUN_W2)                                                        !Set 'High'    callback
     result = DLGSETSUB(dlg, normal, RUN_W2)                                                      !Set 'Normal'  callback
     result = DLGSETSUB(dlg, low, RUN_W2)                                                         !Set 'Low'     callback
     result = DLGSETSUB(dlg, lowest, RUN_W2)                                                      !Set 'Lowest'  callback
     result = DLGSETSUB(dlg, idle, RUN_W2)                                                        !Set 'Idle'    callback
     result = DLGSETSUB(dlg, idok, RUN_W2)
     result = DLGMODAL(dlg)                                                                       !Show dialog box
     call DLGUNINIT(dlg)                                                                          !Close dialog box
     end subroutine W2_DIALOG
