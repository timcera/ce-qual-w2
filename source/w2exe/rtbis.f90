!*==rtbis.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     function RTBIS(X1, X2, Xacc, U, Y, V)
     use GLOBAL, ONLY:wrn, jw
     use EDDY, ONLY:AZC
     use SCREENC, ONLY:jday
     use PREC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! PARAMETER definitions
!
     integer, parameter :: MAXIT = 70
!
! Dummy arguments
!
     real(r8) :: U, V, X1, X2, Xacc, Y
     real(r8) :: RTBIS
     intent (in) U, V, X1, X2, Xacc, Y
!
! Local variables
!
     real(r8) :: dx, f, fmid, xmid
     integer :: j
     real(r8), external :: SEMILOG
!
!*** End of declarations rewritten by SPAG
!
!    USING BISECTION, FIND THE ROOT OF A FUNCTION FUNC KNOWN TO LIE BETWEEN X1
!    AND X2. THE ROOT, RETURNED AS RTBIS, WILL BE REFINED UNTIL ITS ACCURACY
!    IS ±XACC. PARAMETER: MAXIT IS THE MAXIMUM ALLOWED NUMBER OF BISECTIONS.
     fmid = SEMILOG(X2, U, Y, V)                                                  !80
     f = SEMILOG(X1, U, Y, V)
     if(f*fmid>=0.0)then
         write(wrn, *)'RTBIS: ROOT MUST BE BRACKETED IN TKE1 ROUTINE ON JDAY:',&
                    & jday
         write(wrn, *)'Setting TKE calculation to TKE from TKE1'
         RTBIS = X1
         AZC(jw) = '     TKE'
         return
     endif
     if(f<0.0)then  !ORIENT THE SEARCH SO THAT F>0 LIES AT X+DX.
         RTBIS = X1
         dx = X2 - X1
     else
         RTBIS = X2
         dx = X1 - X2
     endif
     do j = 1, MAXIT
               !BISECTION LOOP.
         dx = dx*0.5
         xmid = RTBIS + dx
         fmid = SEMILOG(xmid, U, Y, V)
         if(fmid<=0.0)RTBIS = xmid
         if(ABS(dx)<Xacc .OR. fmid==0.0)return
     enddo
     write(wrn, *)'RTBIS: TOO MANY BISECTIONS IN TKE1 ROUTINE ON JDAY:', jday
     write(wrn, *)'Setting TKE calculation to TKE from TKE1'
     AZC(jw) = '     TKE'
     end function RTBIS
