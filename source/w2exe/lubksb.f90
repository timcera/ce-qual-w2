!*==lubksb.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                              S U B R O U T I N E   L U B K S B                                                **
!***********************************************************************************************************************************
 
     subroutine LUBKSB(A, N, Np, Indx, B)
     use F77KINDS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     integer :: N, Np
     real(R8KIND), dimension(Np, Np) :: A
     real(R8KIND), dimension(N) :: B
     integer, dimension(Np) :: Indx
     intent (in) A, Indx, N, Np
     intent (inout) B
!
! Local variables
!
     integer :: i, ii, j, ll
     real(R8KIND) :: sum
!
!*** End of declarations rewritten by SPAG
!
 
     ii = 0
     do i = 1, N
         ll = Indx(i)
         sum = B(ll)
         B(ll) = B(i)
         if(ii/=0)then
             do j = ii, i - 1
                 sum = sum - A(i, j)*B(j)
             enddo
         elseif(sum/=0.0)then
             ii = i
         endif
         B(i) = sum
     enddo
     do i = N, 1, -1
         sum = B(i)
         do j = i + 1, N
             sum = sum - A(i, j)*B(j)
         enddo
         B(i) = sum/A(i, i)
     enddo
     end subroutine LUBKSB
