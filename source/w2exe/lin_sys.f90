!*==lin_sys.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!
!
     subroutine LIN_SYS(A11, A12, A21, A22, B1, B2, X1, X2, Nflog, Nfcle)
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(8) :: A11, A12, A21, A22, B1, B2, X1, X2
     integer(I1KIND) :: Nfcle, Nflog
     intent (in) A11, A12, A21, A22, B1, B2, Nfcle, Nflog
     intent (out) X1, X2
!
! Local variables
!
     real :: w2err
!
!*** End of declarations rewritten by SPAG
!
!!from 03-Nov-2003 version of Q2KMaster
!!This subroutine solves a linear system of 2 equations and 2 unknowns
     if(A11*A22==A12*A21)then
!	!MsgBox "The sediment flux solution matrix is singular: " & a11 & ", " & a12 & ", " & a21 & ", " & a22
         write(Nflog, '(a)')'The sediment flux solution matrix is singular: '
         write(Nflog, *)'a11  == ', A11, 'a12 = ', A12, 'a21 = ', A21,         &
                      & 'a22 = ', A22
         write(Nfcle, '(a)')'The sediment flux solution matrix is singular: '
         write(Nfcle, *)'a11  == ', A11, 'a12 = ', A12, 'a21 = ', A21,         &
                      & 'a22 = ', A22
         write(Nfcle, *)                                                       &
     &'Error in the solution of linear system of equations used in sediment dia&
     &genesis model'
         write(Nflog, *)                                                       &
     &'Error in the solution of linear system of equations used in sediment dia&
     &genesis model'
         write(w2err, *)                                                       &
     &'Error in Sediment Diagenesis - see log files: Please review the sediment&
     & diagensis parameters'
         stop 'Please review the sediment diagensis parameters'
     endif
     X1 = (A22*B1 - A12*B2)/(A11*A22 - A12*A21)
     X2 = (A11*B2 - A21*B1)/(A11*A22 - A12*A21)
 
     end subroutine LIN_SYS
