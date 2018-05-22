!*==gasdev.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     function GASDEV(Idum3)
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     integer :: Idum3
     real :: GASDEV
     intent (in) Idum3
!
! Local variables
!
     real :: fac, gset, r, v1, v2
     integer, save :: iset
!
!*** End of declarations rewritten by SPAG
!
 
!    From S. Li, PSU, PArticle Transport Algorithm
     data iset/0/
     if(iset==0)then
         do
             v1 = 2.*RAN(Idum3) - 1
             v2 = 2.*RAN(Idum3) - 1
             r = v1**2 + v2**2
             if(r<1)then
                 fac = SQRT( - 2.*LOG(r)/r)
                 gset = v1*fac
                 GASDEV = v2*fac
                 iset = 1
                 exit
             endif
         enddo
     else
         GASDEV = gset
         iset = 0
     endif
 
     end function GASDEV
