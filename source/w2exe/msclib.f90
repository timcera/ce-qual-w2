!*==msclib.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
!***********************************************************************************************************************************
!**                                                                                                                               **
!**                                                         CE-QUAL-W2                                                            **
!**                                            A Two-dimensional, Laterally Averaged,                                             **
!**                                             Hydrodynamic and Water Quality Model                                              **
!**                                                            for                                                                **
!**                                           Rivers, Lakes, Reservoirs, and Estuaries                                            **
!**                                                                                                                               **
!**                                                       Version 4.1                                                             **
!**                                                                                                                               **
!**                                                  Thomas M. Cole, Retired                                                      **
!**                                                Water Quality Modeling Group                                                   **
!**                                                U.S. Army Corps of Engineers                                                   **
!**                                                Waterways Experiment Station                                                   **
!**                                                Vicksburg, Mississippi 39180                                                   **
!**                                                                                                                               **
!**                                                        Scott A. Wells                                                         **
!**                                       Department of Civil and Environmental Engineering                                       **
!**                                                  Portland State University                                                    **
!**                                                         PO Box 751                                                            **
!**                                                 Portland, Oregon  97207-0751                                                  **
!**                                                 phone number: (503) 725-4276                                                  **
!**                                                 fax   number: (503) 725-5950                                                  **
!**                                                 e-mail: wellss@pdx.edu                                                        **
!**                                                                                                                               **
!***********************************************************************************************************************************
 
!***********************************************************************************************************************************
!**                                                                                                                               **
!**                  The long arm of the lawyers has found its way into the water quality modeling arena, so:                     **
!**                                                                                                                               **
!**  This model was developed by the U.S. Army Engineer Waterways Experiment Station, Vicksburg, MS and is maintained by          **
!**  Portland State University.  Portland State University and the US government and its components are not responsible           **
!**  for any damages,including incidental or consequential damages, arising                                                       **
!**  from use or misuse of this model, or from results achieved or conclusions drawn by others.  Distribution of this model is    **
!**  restricted by the Export Administration Act of 1969,  50 app. USC subsections 2401-2420, as amended, and other applicable    **
!**  laws or regulations.                                                                                                         **
!**                                                                                                                               **
!***********************************************************************************************************************************
 
!***********************************************************************************************************************************
!**                                                      Module Declaration                                                       **
!***********************************************************************************************************************************
     module MSCLIB
     implicit none
     include "RESOURCE.FD"
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: hthread
     logical :: restart_exists, restart_pushed, stopped, stop_pushed
!
!*** End of declarations rewritten by SPAG
!
     interface
     function $BEGINTHREADEX(SECURITY, STACK_SIZE, START_ADDRESS, ARGLIST,     &
                           & INITFLAG, THRDADDR)
     use DFWINTY, RENAMED => DLT
      !DEC$ ATTRIBUTES C,ALIAS : "__BEGINTHREADEX" :: $BEGINTHREADEX
      !DEC$ ATTRIBUTES REFERENCE,ALLOW_NULL        :: SECURITY
      !DEC$ ATTRIBUTES REFERENCE,IGNORE_LOC        :: THRDADDR
     integer(UINT) :: $BEGINTHREADEX
     integer(UINT), intent(in) :: STACK_SIZE, INITFLAG
     integer(PVOID), intent(in) :: START_ADDRESS, ARGLIST
     integer(UINT), intent(out) :: THRDADDR
     type(T_SECURITY_ATTRIBUTES), INTENT(IN)::SECURITY
     end function $BEGINTHREADEX
     endinterface
     interface
     subroutine $ENDTHREADEX(RETVAL)
     use DFWINTY, RENAMED => DLT
      !DEC$ ATTRIBUTES C, ALIAS : "__ENDTHREADEX" :: $ENDTHREADEX
     integer(UINT), intent(in) :: RETVAL
     end subroutine $ENDTHREADEX
     endinterface
     end module MSCLIB
