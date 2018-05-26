!*==zbrent1.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                                  F U N C T I O N   Z B R E N T                                                **
!***********************************************************************************************************************************
 
     function ZBRENT1(X1, X2, Tol, Barg)
     use f77kinds
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
     real, parameter :: FACTOR = 0.1, NTRY = 50, ITMAX = 100, EPS = 3.E-8
!
! Dummy arguments
!
     real(R8KIND) :: Barg, Tol, X1, X2
     real :: ZBRENT1
     intent (in) Tol
     intent (inout) X1, X2
!
! Local variables
!
     real(R8KIND) :: b, ba, c, d, e, f1, f2, fa, fb, fc, p, q, r, s, tol1, xm
     real(R8KIND), external :: CDFUNC
     integer :: i, iter, j
!
!*** End of declarations rewritten by SPAG
!
 
     f1 = CDFUNC(X1, Barg)
     f2 = CDFUNC(X2, Barg)
     if(f1<=0.0)then
         do i = 1, 40
             X1 = X1/10.0
             f1 = CDFUNC(X1, Barg)
             if(f1>0.0)exit
         enddo
     endif
     do j = 1, NTRY
         if(f1*f2<0.0)exit
         if(ABS(f1)<ABS(f2))then
             X1 = X1 + FACTOR*(X1 - X2)
             f1 = CDFUNC(X1, Barg)
         else
             X2 = X2 + FACTOR*(X2 - X1)
             f2 = CDFUNC(X2, Barg)
         endif
     enddo
     ba = X1
     b = X2
     fa = CDFUNC(ba, Barg)
     fb = CDFUNC(b, Barg)
     fc = fb
     do iter = 1, ITMAX
         if(fb*fc>0.0)then
             c = ba
             fc = fa
             d = b - ba
             e = d
         endif
         if(ABS(fc)<ABS(fb))then
             ba = b
             b = c
             c = ba
             fa = fb
             fb = fc
             fc = fa
         endif
         tol1 = 2.0*EPS*ABS(b) + 0.5*Tol
         xm = 0.5*(c - b)
         if(ABS(xm)<=tol1 .OR. fb==0.0)then
             ZBRENT1 = b
             exit
         endif
         if(ABS(e)>=tol1 .AND. ABS(fa)>ABS(fb))then
             s = fb/fa
             if(ba==c)then
                 p = 2.0*xm*s
                 q = 1.0 - s
             else
                 q = fa/fc
                 r = fb/fc
                 p = s*(2.*xm*q*(q - r) - (b - ba)*(r - 1.0))
                 q = (q - 1.0)*(r - 1.0)*(s - 1.0)
             endif
             if(p>0.0)q = -q
             p = ABS(p)
             if(2.0*p<MIN(3.0*xm*q - ABS(tol1*q), ABS(e*q)))then
                 e = d
                 d = p/q
             else
                 d = xm
                 e = d
             endif
         else
             d = xm
             e = d
         endif
         ba = b
         fa = fb
         if(ABS(d)>tol1)then
             b = b + d
         else
             b = b + SIGN(tol1, xm)
         endif
         fb = CDFUNC(b, Barg)
     enddo
     ZBRENT1 = b
     end function ZBRENT1
