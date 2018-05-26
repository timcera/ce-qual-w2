!*==tridiag.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!*                                              S U B R O U T I N E    T R I D I A G                                              **
!***********************************************************************************************************************************
 
     subroutine TRIDIAG(A, V, C, D, S, E, N, U)
     use f77kinds
     use TRIDIAG_V
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     integer :: E, N, S
     real(R8KIND), dimension(E) :: A, C, D, V
     real(R8KIND), dimension(N) :: U
     intent (in) A, C, D, E, N, S, V
     intent (out) U
!
! Local variables
!
     integer :: i
!
!*** End of declarations rewritten by SPAG
!
 ! REAL(R8), ALLOCATABLE, DIMENSION(:)              :: BTA, GMA
 ! REAL(R8), DIMENSION(1000)              :: BTA, GMA
!    ALLOCATE (BTA(N),GMA(N))
 
     BTA1(S) = V(S)
     GMA1(S) = D(S)
     do i = S + 1, E
         BTA1(i) = V(i) - A(i)/BTA1(i - 1)*C(i - 1)
         GMA1(i) = D(i) - A(i)/BTA1(i - 1)*GMA1(i - 1)
     enddo
     U(E) = GMA1(E)/BTA1(E)
     do i = E - 1, S, -1
         U(i) = (GMA1(i) - C(i)*U(i + 1))/BTA1(i)
     enddo
!    DEALLOCATE (BTA, GMA)                                                    
!    ! SW 10/17/05
     end subroutine TRIDIAG
