!*==ludcmp.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                              S U B R O U T I N E   L U D C M P                                                **
!***********************************************************************************************************************************
 
     subroutine LUDCMP(A, N, Np, Indx, D)
     use F77KINDS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
     real, parameter :: TINY = 1.0E-20
!
! Dummy arguments
!
     real(R8KIND) :: D
     integer :: N, Np
     real(R8KIND), dimension(Np, Np) :: A
     integer, dimension(Np) :: Indx
     intent (in) N, Np
     intent (out) Indx
     intent (inout) A, D
!
! Local variables
!
     real(R8KIND) :: aamax, dum, sum
     integer :: i, imax, j, k
     real(R8KIND), dimension(500) :: vv
!
!*** End of declarations rewritten by SPAG
!
 
 
     D = 1.0
     do i = 1, N
         aamax = 0.0
         do j = 1, N
             if(ABS(A(i, j))>aamax)aamax = ABS(A(i, j))
         enddo
         vv(i) = 1.0/aamax
     enddo
     do j = 1, N
         do i = 1, j - 1
             sum = A(i, j)
             do k = 1, i - 1
                 sum = sum - A(i, k)*A(k, j)
             enddo
             A(i, j) = sum
         enddo
         aamax = 0.0
         do i = j, N
             sum = A(i, j)
             do k = 1, j - 1
                 sum = sum - A(i, k)*A(k, j)
             enddo
             A(i, j) = sum
             dum = vv(i)*ABS(sum)
             if(dum>=aamax)then
                 imax = i
                 aamax = dum
             endif
         enddo
         if(j/=imax)then
             do k = 1, N
                 dum = A(imax, k)
                 A(imax, k) = A(j, k)
                 A(j, k) = dum
             enddo
             D = -D
             vv(imax) = vv(j)
         endif
         Indx(j) = imax
         if(A(j, j)==0.0)A(j, j) = TINY
         if(j/=N)then
             dum = 1.0/A(j, j)
             do i = j + 1, N
                 A(i, j) = A(i, j)*dum
             enddo
         endif
     enddo
     end subroutine LUDCMP
