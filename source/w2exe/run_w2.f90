!*==run_w2.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!*                                                  S U B R O U T I N E   R U N   W 2                                             **
!***********************************************************************************************************************************
 
     subroutine RUN_W2(Dlg, Control_name, Action)
     use DFWIN, RENAMED => DLT                                                                                               !Rename DLT in DFWIN
     use DFLOGM
     use MSCLIB
     use GLOBAL, ONLY:cdate, cctime
     use MAIN, ONLY:end_run
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     integer :: Action, Control_name
     type(DIALOG) :: Dlg
     character(72) :: Text
     intent (in) Control_name
!
! Local variables
!
     integer :: i, idthread, result
     type(T_SECURITY_ATTRIBUTES), pointer :: null_sa
     character(1000) :: text1
     character(8) :: time
!
!*** End of declarations rewritten by SPAG
!
     interface
     integer(4) function ce_qual_w2(H)
    !DEC$ATTRIBUTES STDCALL :: CE_QUAL_W2
     integer H
     end function ce_qual_w2
     endinterface
 
     selectcase(Control_name)
     case(RUN)
         call BLANK_DIALOG(Dlg)
         call DATE_AND_TIME(cdate, cctime)
         restart_pushed = .FALSE.
         stop_pushed = .FALSE.
         hthread = CREATETHREAD(null_sa, 0, LOC(ce_qual_w2), LOC(Dlg), 0,      &
                 & LOC(idthread))                                                                  !Start W2 in a new thread
         time = cctime(1:2) // ':' // cctime(3:4) // ':' // cctime(5:6)
         result = DLGSET(Dlg, RUN, .FALSE., dlg_enable)                                            !Disable 'Run'     button
         result = DLGSET(Dlg, close, .FALSE., dlg_enable)                                          !Disable 'Run'     button
         result = DLGSET(Dlg, RESTART, .FALSE., dlg_enable)                                        !Disable 'Restart' button
         result = DLGSET(Dlg, STOP_EXECUTION, .TRUE., dlg_enable)                                  !Enable  'Stop'    button
         result = DLGSET(Dlg, starting_time, time)                                                 !Display starting time
         result = DLGSET(Dlg, status, 'Executing')                                                 !Display execution status
     case(OUTPUT_DIALOG)
         call BLANK_DIALOG(Dlg)
         call DATE_AND_TIME(cdate, cctime)
         restart_pushed = .FALSE.
         stop_pushed = .FALSE.
         hthread = CREATETHREAD(null_sa, 0, LOC(ce_qual_w2), LOC(Dlg), 0,      &
                 & LOC(idthread))                                                                  !Start W2 in a new thread
         time = cctime(1:2) // ':' // cctime(3:4) // ':' // cctime(5:6)
         result = DLGSET(Dlg, RUN, .FALSE., dlg_enable)                                            !Disable 'Run'     button
         result = DLGSET(Dlg, close, .FALSE., dlg_enable)                                          !Disable 'Run'     button
         result = DLGSET(Dlg, RESTART, .FALSE., dlg_enable)                                        !Disable 'Restart' button
         result = DLGSET(Dlg, STOP_EXECUTION, .TRUE., dlg_enable)                                  !Enable  'Stop'    button
         result = DLGSET(Dlg, starting_time, time)                                                 !Display starting time
         result = DLGSET(Dlg, status, 'Executing')                                                 !Display execution status
     case(STOP_EXECUTION)
         stop_pushed = .TRUE.
         end_run = .TRUE.
         result = DLGSET(Dlg, STOP_EXECUTION, .FALSE., dlg_enable)                                 !Disable 'Stop'    button
         result = DLGSET(Dlg, RUN, .TRUE., dlg_enable)                                             !Enable  'Run'     button
         result = CLOSEHANDLE(hthread)                                                             !Close thread handle
     case(RESTART)
         stop_pushed = .FALSE.
         restart_pushed = .TRUE.
         end_run = .FALSE.
         hthread = CREATETHREAD(null_sa, 0, LOC(ce_qual_w2), LOC(Dlg), 0,      &
                 & LOC(idthread))                                                                  !Start W2 in a new thread
         result = DLGSET(Dlg, ending_time, ' ')
         result = DLGSET(Dlg, RUN, .FALSE., dlg_enable)                                            !Disable 'Run'     button
         result = DLGSET(Dlg, close, .FALSE., dlg_enable)                                          !Disable 'Run'     button
         result = DLGSET(Dlg, RESTART, .FALSE., dlg_enable)                                        !Disable 'Restart' button
         result = DLGSET(Dlg, STOP_EXECUTION, .TRUE., dlg_enable)                                  !Enable  'Stop'    button
         result = DLGSET(Dlg, status, 'Executing')                                                 !Display execution status
     case(HIGHEST)
         call ENABLE(Dlg)                                                                          !Enable  priority buttons
         i = SETTHREADPRIORITY(hthread, thread_priority_highest)                                   !Set highest priority
         result = DLGSET(Dlg, HIGHEST, .FALSE., dlg_enable)                                        !Disable 'Highest' button
     case(HIGH)
         call ENABLE(Dlg)                                                                          !Enable  priority buttons
         i = SETTHREADPRIORITY(hthread, thread_priority_above_normal)                              !Set high    priority
         result = DLGSET(Dlg, HIGH, .FALSE., dlg_enable)                                           !Disable 'High'   button
     case(NORMAL)
         call ENABLE(Dlg)                                                                          !Enable  priority buttons
         i = SETTHREADPRIORITY(hthread, thread_priority_normal)                                    !Set normal  priority
         result = DLGSET(Dlg, NORMAL, .FALSE., dlg_enable)                                         !Disable 'Normal' button
     case(LOW)
         call ENABLE(Dlg)                                                                          !Enable  priority buttons
         i = SETTHREADPRIORITY(hthread, thread_priority_below_normal)                              !Set low     priority
         result = DLGSET(Dlg, LOW, .FALSE., dlg_enable)                                            !Disable 'Low'    button
     case(LOWEST)
         call ENABLE(Dlg)                                                                          !Enable  priority buttons
         i = SETTHREADPRIORITY(hthread, thread_priority_lowest)                                    !Set lowest  priority
         result = DLGSET(Dlg, LOWEST, .FALSE., dlg_enable)                                         !Disable 'Lowest' button
     case(IDLE)
         call ENABLE(Dlg)                                                                          !Enable  priority buttons
         i = SETTHREADPRIORITY(hthread, thread_priority_idle)                                      !Set idle priority
         result = DLGSET(Dlg, IDLE, .FALSE., dlg_enable)                                           !Disable 'Idle'   button
     endselect
     return
 
     entry EXITDIALOG(Dlg, Text)
     call DLGEXIT(Dlg)
     return
 
     entry STOP_W2(Dlg, Text)
     call DATE_AND_TIME(cdate, cctime)
     text1 = cctime(1:2) // ':' // cctime(3:4) // ':' // cctime(5:6)
     result = DLGSET(Dlg, ending_time, text1)                                                      !Display ending time
     result = DLGSET(Dlg, status, Text)                                                            !Execution status
     result = DLGSET(Dlg, STOP_EXECUTION, .FALSE., dlg_enable)                                     !Disable 'Stop'    button
     result = DLGSET(Dlg, RUN, .TRUE., dlg_enable)                                                 !Enable  'Run'     button
     result = DLGSET(Dlg, close, .TRUE., dlg_enable)                                               !Enable  'Close'    button
     if(stop_pushed)result = DLGSET(Dlg, RESTART, .TRUE., dlg_enable)                              !Enable  'Restart' button
     end subroutine RUN_W2
