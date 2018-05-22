!*==output.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                           S U B R O U T I N E   O U T P U T                                                   **
!***********************************************************************************************************************************
 
     subroutine OUTPUT(Jday, Iupr, Idpr, Kbr, Isnp, Bl, Nbl)
     use GLOBAL
     use GDAYC
     use GEOMC
     use KINETIC
     use TVDC
     use NAMESC
     use LOGICC
     use MACROPHYTEC
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     integer :: Idpr, Iupr, Kbr, Nbl
     real :: Jday
     integer, dimension(imx) :: Bl
     integer, dimension(imx, nwb) :: Isnp
     intent (in) Idpr, Jday, Kbr
     intent (out) Bl
     intent (inout) Iupr, Nbl
!
! Local variables
!
     character(10), save :: blank
     integer :: j, ja, jac, jd, je, jh, jj, jt, k, l, nlines
     character(8) :: lfac
     real :: limit
     logical :: new_page
!
!*** End of declarations rewritten by SPAG
!
 
!    Type declaration
 
 
!    Data declaration
 
     data blank/'          '/
 
!    Variable initialization
 
     new_page = .TRUE.
 
!    Blank inactive cells
 
     Nbl = 1
     jb = 1
     Iupr = 1
     do i = 1, Idpr - 1
         if(CUS(jb)>Isnp(i, jw))then
             Bl(Nbl) = i
             Nbl = Nbl + 1
             if(jb==1)Iupr = i + 1
         endif
         if(Isnp(i + 1, jw)>DS(jb))jb = jb + 1
     enddo
     Nbl = Nbl - 1
 
!    Water surface elevation, water surface deviation, ice cover, and sediment
 
!    oxygen demand
     conv(1, :) = blank
     do i = Iupr, Idpr
         do jjb = 1, nbr
             if(Isnp(i, jw)>=US(jjb) - 1 .AND. Isnp(i, jw)<=DS(jjb) + 1)exit
         enddo
         write(conv(1, i), '(F10.3)')EL(KTWB(jw), Isnp(i, jw)) - Z(Isnp(i, jw))&
                                   & *COSA(jjb)
     enddo
     write(SNP(jw), '(/A//2X,1000I10)')'          Water Surface, m',           &
                                     & (Isnp(i, jw), i = Iupr, Idpr)
     write(SNP(jw), '(2X,1000A10/)')(conv(1, i), i = Iupr, Idpr)
     do i = Iupr, Idpr
         write(conv(1, i), '(F10.4)')SNGL(Z(Isnp(i, jw)))
     enddo
     write(SNP(jw), '(/A//2X,1000I10)')                                        &
                   &'          Water Surface Deviation (positive downwards), m'&
                  & , (Isnp(i, jw), i = Iupr, Idpr)
     write(SNP(jw), '(2X,1000A10/)')(conv(1, i), i = Iupr, Idpr)
     if(ICE_CALC(jw))then
         do i = Iupr, Idpr
             write(conv(1, i), '(F10.3)')ICETH(Isnp(i, jw))
         enddo
         write(SNP(jw), '(/A//3X,1000A10)')'          Ice Thickness, m',       &
             & (conv(1, i), i = Iupr, Idpr)
     endif
     if(constituents)then
         do i = Iupr, Idpr
             write(conv(1, i), '(F10.3)')SOD(Isnp(i, jw))*day
         enddo
         if(oxygen_demand)write(SNP(jw), '(/A//3X,1000A10/)')                  &
                               &'          Sediment Oxygen Demand, g/m^2/day', &
                              & (conv(1, i), i = Iupr, Idpr)
     endif
 
!    Hydrodynamic variables and temperature
 
     conv = blank
     do jh = 1, nhy
         l = LEN_TRIM(FMTH(jh))
         if(PRINT_HYDRO(jh, jw))then
             do i = Iupr, Idpr
                 if(jh==1)then
                     do k = KTWB(jw), KB(Isnp(i, jw))
                         write(conv(k, i), FMTH(jh))                           &
                             & INT(HYD(k, Isnp(i, jw), jh))
                     enddo
                 elseif(jh>6)then
                     do k = KTWB(jw), KB(Isnp(i, jw))
                         write(conv(k, i), FMTH(jh)(1:l))                      &
                             & HYD(k, Isnp(i, jw), jh)*dlt
                     enddo
                 else
                     do k = KTWB(jw), KB(Isnp(i, jw))
                         write(conv(k, i), FMTH(jh)(1:l))                      &
                             & HYD(k, Isnp(i, jw), jh)*HMULT(jh)
                     enddo
                 endif
             enddo
             if(new_page)then
                 write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                 nlines = kmx - KTWB(jw) + 14
             endif
             nlines = nlines + kmx - KTWB(jw) + 11
             new_page = nlines>72
             write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,A/)')month, gday, ',',     &
                 & year, '    Julian day = ', INT(Jday), ' days ',             &
                 & (Jday - INT(Jday))*24.0, ' hours   ' // HNAME(jh)
             write(SNP(jw), '(1X,A,1000I10)')'Layer  Depth',                   &
                 & (Isnp(i, jw), i = Iupr, Idpr)
             do k = KTWB(jw), Kbr
                 write(SNP(jw), '(1X,I4,F8.2,1000A10)')k, DEPTHM(k, DS(BS(jw)))&
                     & , (conv(k, i), i = Iupr, Idpr)
             enddo
         endif
     enddo
 
!    Constituent concentrations
 
     if(constituents)then
         do jac = 1, nac
             jc = CN(jac)
             l = LEN_TRIM(FMTC(jc))
             if(PRINT_CONST(jc, jw))then
                 do i = Iupr, Idpr
                     do k = KTWB(jw), KB(Isnp(i, jw))
                         write(conv(k, i), FMTC(jc)(1:l))C2(k, Isnp(i, jw), jc)&
                             & *CMULT(jc)
                     enddo
                 enddo
                 if(new_page)then
                     write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                     nlines = kmx - KTWB(jw) + 14
                 endif
                 nlines = nlines + kmx - KTWB(jw) + 11
                 new_page = nlines>72
                 write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,A/)')month, gday, ',', &
                     & year, '    Julian Date ', INT(Jday), ' days ',          &
                     & (Jday - INT(Jday))*24.0, ' hours   ' // CNAME(jc)
                 write(SNP(jw), '(1X,A,1000I10)')'Layer  Depth',               &
                     & (Isnp(i, jw), i = Iupr, Idpr)
                 do k = KTWB(jw), Kbr
                     write(SNP(jw), '(1X,I4,F8.2,1000A10)')k,                  &
                         & DEPTHM(k, DS(BS(jw))), (conv(k, i), i = Iupr, Idpr)
                 enddo
             endif
         enddo
 
!**      Derived constituent concentrations
 
         do jd = 1, ndc
             l = LEN_TRIM(FMTCD(jd))
             if(PRINT_DERIVED(jd, jw))then
                 do i = Iupr, Idpr
                     do k = KTWB(jw), KB(Isnp(i, jw))
                         write(conv(k, i), FMTCD(jd)(1:l))                     &
                             & CD(k, Isnp(i, jw), jd)*CDMULT(jd)
                     enddo
                 enddo
                 if(new_page)then
                     write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                     nlines = kmx - KTWB(jw) + 14
                 endif
                 nlines = nlines + kmx - KTWB(jw) + 11
                 new_page = nlines>72
                 write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,A/)')month, gday, ',', &
                     & year, '    Julian Date ', INT(Jday), ' days ',          &
                     & (Jday - INT(Jday))*24.0, ' hours    ' // CDNAME(jd)
                 write(SNP(jw), '(1X,A,1000I10)')'Layer  Depth',               &
                     & (Isnp(i, jw), i = Iupr, Idpr)
                 do k = KTWB(jw), Kbr
                     write(SNP(jw), '(1X,I4,F8.2,1000A10)')k,                  &
                         & DEPTHM(k, DS(BS(jw))), (conv(k, i), i = Iupr, Idpr)
                 enddo
             endif
         enddo
 
!**      Sediment
 
         if(PRINT_SEDIMENT(jw))then
             do i = Iupr, Idpr
                 do k = KTWB(jw), KB(Isnp(i, jw))
                     write(conv(k, i), '(F10.2)')SED(k, Isnp(i, jw))
                 enddo
             enddo
             if(new_page)then
                 write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                 nlines = kmx - KTWB(jw) + 14
             endif
             nlines = nlines + kmx - KTWB(jw) + 11
             new_page = nlines>72
             write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,A/)')month, gday, ',',     &
                 & year, '    Julian Date ', INT(Jday), ' days ',              &
                 & (Jday - INT(Jday))*24.0,                                    &
                  &' hours     Organic sediments, g/m^3'
             write(SNP(jw), '(1X,A,1000I10)')'Layer  Depth',                   &
                 & (Isnp(i, jw), i = Iupr, Idpr)
             do k = KTWB(jw), Kbr
                 write(SNP(jw), '(1X,I4,F8.2,1000A10)')k, DEPTHM(k, DS(BS(jw)))&
                     & , (conv(k, i), i = Iupr, Idpr)
             enddo
         endif
 
         if(PRINT_SEDIMENT(jw))then
             do i = Iupr, Idpr
                 do k = KTWB(jw), KB(Isnp(i, jw))
                     write(conv(k, i), '(F10.2)')SEDP(k, Isnp(i, jw))
                 enddo
             enddo
             if(new_page)then
                 write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                 nlines = kmx - KTWB(jw) + 14
             endif
             nlines = nlines + kmx - KTWB(jw) + 11
             new_page = nlines>72
             write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,A/)')month, gday, ',',     &
                 & year, '    Julian Date ', INT(Jday), ' days ',              &
                 & (Jday - INT(Jday))*24.0,                                    &
                  &' hours     Organic phosphorus sediments, g/m^3'
             write(SNP(jw), '(1X,A,1000I10)')'Layer  Depth',                   &
                 & (Isnp(i, jw), i = Iupr, Idpr)
             do k = KTWB(jw), Kbr
                 write(SNP(jw), '(1X,I4,F8.2,1000A10)')k, DEPTHM(k, DS(BS(jw)))&
                     & , (conv(k, i), i = Iupr, Idpr)
             enddo
         endif
 
         if(PRINT_SEDIMENT(jw))then
             do i = Iupr, Idpr
                 do k = KTWB(jw), KB(Isnp(i, jw))
                     write(conv(k, i), '(F10.2)')SEDN(k, Isnp(i, jw))
                 enddo
             enddo
             if(new_page)then
                 write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                 nlines = kmx - KTWB(jw) + 14
             endif
             nlines = nlines + kmx - KTWB(jw) + 11
             new_page = nlines>72
             write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,A/)')month, gday, ',',     &
                 & year, '    Julian Date ', INT(Jday), ' days ',              &
                 & (Jday - INT(Jday))*24.0,                                    &
                  &' hours     Organic nitrogen sediments, g/m^3'
             write(SNP(jw), '(1X,A,1000I10)')'Layer  Depth',                   &
                 & (Isnp(i, jw), i = Iupr, Idpr)
             do k = KTWB(jw), Kbr
                 write(SNP(jw), '(1X,I4,F8.2,1000A10)')k, DEPTHM(k, DS(BS(jw)))&
                     & , (conv(k, i), i = Iupr, Idpr)
             enddo
         endif
 
         if(PRINT_SEDIMENT(jw))then
             do i = Iupr, Idpr
                 do k = KTWB(jw), KB(Isnp(i, jw))
                     write(conv(k, i), '(F10.2)')SEDC(k, Isnp(i, jw))
                 enddo
             enddo
             if(new_page)then
                 write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                 nlines = kmx - KTWB(jw) + 14
             endif
             nlines = nlines + kmx - KTWB(jw) + 11
             new_page = nlines>72
             write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,A/)')month, gday, ',',     &
                 & year, '    Julian Date ', INT(Jday), ' days ',              &
                 & (Jday - INT(Jday))*24.0,                                    &
                  &' hours     Organic carbon sediments, g/m^3'
             write(SNP(jw), '(1X,A,1000I10)')'Layer  Depth',                   &
                 & (Isnp(i, jw), i = Iupr, Idpr)
             do k = KTWB(jw), Kbr
                 write(SNP(jw), '(1X,I4,F8.2,1000A10)')k, DEPTHM(k, DS(BS(jw)))&
                     & , (conv(k, i), i = Iupr, Idpr)
             enddo
         endif
 
!        Amaila Start
         if(PRINT_SEDIMENT1(jw))then
             do i = Iupr, Idpr
                 do k = KTWB(jw), KB(Isnp(i, jw))
                     write(conv(k, i), '(F10.2)')SED1(k, Isnp(i, jw))
                 enddo
             enddo
             if(new_page)then
                 write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                 nlines = kmx - KTWB(jw) + 14
             endif
             nlines = nlines + kmx - KTWB(jw) + 11
             new_page = nlines>72
             write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,A/)')month, gday, ',',     &
                 & year, '    Julian Date ', INT(Jday), ' days ',              &
                 & (Jday - INT(Jday))*24.0,                                    &
                  &' hours     Tree Organic Matter-Labile, g/m^3'
             write(SNP(jw), '(1X,A,1000I10)')'Layer  Depth',                   &
                 & (Isnp(i, jw), i = Iupr, Idpr)
             do k = KTWB(jw), Kbr
                 write(SNP(jw), '(1X,I4,F8.2,1000A10)')k, DEPTHM(k, DS(BS(jw)))&
                     & , (conv(k, i), i = Iupr, Idpr)
             enddo
         endif
 
         if(PRINT_SEDIMENT2(jw))then
             do i = Iupr, Idpr
                 do k = KTWB(jw), KB(Isnp(i, jw))
                     write(conv(k, i), '(F10.2)')SED2(k, Isnp(i, jw))
                 enddo
             enddo
             if(new_page)then
                 write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                 nlines = kmx - KTWB(jw) + 14
             endif
             nlines = nlines + kmx - KTWB(jw) + 11
             new_page = nlines>72
             write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,A/)')month, gday, ',',     &
                 & year, '    Julian Date ', INT(Jday), ' days ',              &
                 & (Jday - INT(Jday))*24.0,                                    &
                  &' hours     Tree Organic Matter-Refactory, g/m^3'
             write(SNP(jw), '(1X,A,1000I10)')'Layer  Depth',                   &
                 & (Isnp(i, jw), i = Iupr, Idpr)
             do k = KTWB(jw), Kbr
                 write(SNP(jw), '(1X,I4,F8.2,1000A10)')k, DEPTHM(k, DS(BS(jw)))&
                     & , (conv(k, i), i = Iupr, Idpr)
             enddo
         endif
 
!        Amaila End
 
!**      Epiphyton
 
         do je = 1, nep
             if(PRINT_EPIPHYTON(jw, je))then
                 do i = Iupr, Idpr
                     do k = KTWB(jw), KB(Isnp(i, jw))
                         write(conv(k, i), '(F10.2)')EPD(k, Isnp(i, jw), je)
                     enddo
                 enddo
                 if(new_page)then
                     write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                     nlines = kmx - KTWB(jw) + 14
                 endif
                 nlines = nlines + kmx - KTWB(jw) + 11
                 new_page = nlines>72
                 write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,A/)')month, gday, ',', &
                     & year, '    Julian Date ', INT(Jday), ' days ',          &
                     & (Jday - INT(Jday))*24.0, ' hours     Epiphyton, g/m^2'
                 write(SNP(jw), '(1X,A,1000I10)')'Layer  Depth',               &
                     & (Isnp(i, jw), i = Iupr, Idpr)
                 do k = KTWB(jw), Kbr
                     write(SNP(jw), '(1X,I4,F8.2,1000A10)')k,                  &
                         & DEPTHM(k, DS(BS(jw))),                              &
                         & (ADJUSTR(conv(k, i)), i = Iupr, Idpr)
                 enddo
             endif
         enddo
 
!********* macrophytes
         do l = 1, nmc
             conv = blank
             if(PRINT_MACROPHYTE(jw, l))then
                 do i = Iupr, Idpr
                     do k = KTWB(jw), KB(Isnp(i, jw))
                         write(conv(k, i), '(F10.2)')MAC(k, Isnp(i, jw), l)
                     enddo
                 enddo
                 if(new_page)then
                     write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                     nlines = kmx - KTWB(jw) + 14
                 endif
                 nlines = nlines + kmx - KTWB(jw) + 11
                 new_page = nlines>72
                 write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,2A,i0,A/)')month, gday,&
                      &',', year, '    Julian Date ', INT(Jday), ' days ',     &
                     & (Jday - INT(Jday))*24.0, ' hours   ',                   &
                      &'  Macrophyte Group #', l, ' g/m^3'                                                                                !    ! cb 12/20/11
                 write(SNP(jw), '(1X,A,1000I10)')'Layer  Depth',               &
                     & (Isnp(i, jw), i = Iupr, Idpr)
                 do k = KTWB(jw), Kbr
                     write(SNP(jw), '(1X,I4,F8.2,1000A10)')k,                  &
                         & DEPTHM(k, DS(BS(jw))), (conv(k, i), i = Iupr, Idpr)
                 enddo
                 do i = Iupr, Idpr
                     conv2 = blank
                     do k = KTWB(jw), KB(Isnp(i, jw))
                         if(k==KTWB(jw))then
                             jt = KTI(Isnp(i, jw))
                         else
                             jt = k
                         endif
                         je = KB(Isnp(i, jw))
                         do jj = jt, je
                             write(conv2(k, jj), '(F10.2)')                    &
                                 & MACRC(jj, k, Isnp(i, jw), l)
                         enddo
                     enddo
                     if(new_page)then
                         write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                         nlines = kmx - KTWB(jw) + 14
                     endif
                     nlines = nlines + kmx - KTWB(jw) + 11
                     new_page = nlines>72
!                    WRITE (SNP(jw),'(/1X,3(A,1X,I0),A,F0.2,A/)') 
!                    MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days
!                    ',(JDAY-INT(JDAY))    & *24.0,' hours   ','  Macrophyte
!                    Columns g/m^3'
                     write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,2A,i0,A/)')month,  &
                         & gday, ',', year, '    Julian Date ', INT(Jday),     &
                          &' days ', (Jday - INT(Jday))*24.0, ' hours   ',     &
                          &'  Macrophyte Group #', l, ' Columns g/m^3'                                                                      !  ! cb 12/20/11
                     write(SNP(jw), 9001)Isnp(i, jw)
9001                 format(7x, "SEGMENT   ", i8)
                     jt = KTI(Isnp(i, jw))
                     je = KB(Isnp(i, jw))
                     write(SNP(jw), '(A,200I10)')' LAYER  DEPTH',              &
                         & (jj, jj = jt, je)
                     do k = KTWB(jw), KB(Isnp(i, jw))
                         write(SNP(jw), '(1X,I4,F8.2,1000A10)')k,              &
                             & DEPTHM(k, DS(BS(jw))),                          &
                             & (conv2(k, jj), jj = jt, je)
                     enddo
                     write(SNP(jw), 9002)l, Isnp(i, jw)
9002                 format(7x, "MACROPHYTE GROUP #", i0,                      &
                           &" LIMITATION SEGMENT   ", i8)                  ! CB 12/20/11
                     write(SNP(jw), '(A,200I10)')' LAYER  DEPTH',              &
                         & (jj, jj = jt, je)
                     mlfpr = blank
                     do k = KTWB(jw), KB(Isnp(i, jw))
                         if(k==KTWB(jw))then
                             jt = KTI(Isnp(i, jw))
                         else
                             jt = k
                         endif
                         do jj = jt, je
                             limit = MIN(MPLIM(k, Isnp(i, jw), l),             &
                                   & MNLIM(k, Isnp(i, jw), l),                 &
                                   & MCLIM(k, Isnp(i, jw), l),                 &
                                   & MLLIM(jj, k, Isnp(i, jw), l))
                             if(limit==MPLIM(k, Isnp(i, jw), l))then
                                 write(lfac, '(F8.4)')MPLIM(k, Isnp(i, jw), l)
                                 mlfpr(jj, k, Isnp(i, jw), l) = ' P' // lfac
                             elseif(limit==MNLIM(k, i, l))then
                                 write(lfac, '(F8.4)')MNLIM(k, Isnp(i, jw), l)
                                 mlfpr(jj, k, Isnp(i, jw), l) = ' N' // lfac
                             elseif(limit==MCLIM(k, i, l))then
                                 write(lfac, '(F8.4)')MCLIM(k, Isnp(i, jw), l)
                                 mlfpr(jj, k, Isnp(i, jw), l) = ' C' // lfac
                             elseif(limit==MLLIM(jj, k, Isnp(i, jw), l))then
                                 write(lfac, '(F8.4)')                         &
                                     & MLLIM(jj, k, Isnp(i, jw), l)
                                 mlfpr(jj, k, Isnp(i, jw), l) = ' L' // lfac
                             endif
                         enddo
                     enddo
                     jt = KTI(Isnp(i, jw))
                     je = KB(Isnp(i, jw))
                     do k = KTWB(jw), KB(Isnp(i, jw))
                         write(SNP(jw), '(1X,I4,F8.2,1000A10)')k,              &
                             & DEPTHM(k, DS(BS(jw))),                          &
                             & (mlfpr(jj, k, Isnp(i, jw), l), jj = jt, je)
                     enddo
                 enddo
             endif
         enddo
 
!**      Algal nutrient limitations
 
         do ja = 1, nal
             if(LIMITING_FACTOR(ja))then
                 if(new_page)then
                     write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                     nlines = kmx - KTWB(jw) + 14
                 endif
                 nlines = nlines + kmx - KTWB(jw) + 11
                 new_page = nlines>72
                 do i = Iupr, Idpr
                           !mlm 6/30/06
                     do k = KTWB(jw), KB(Isnp(i, jw))
                                         !mlm  6/30/06
                         limit = MIN(APLIM(k, Isnp(i, jw), ja),                &
                               & ANLIM(k, Isnp(i, jw), ja),                    &
                               & ASLIM(k, Isnp(i, jw), ja),                    &
                               & ALLIM(k, Isnp(i, jw), ja))
                         if(limit==APLIM(k, Isnp(i, jw), ja))then
                             write(lfac, '(F8.4)')APLIM(k, Isnp(i, jw), ja)
                             LFPR(k, Isnp(i, jw)) = ' P' // lfac
                         elseif(limit==ANLIM(k, Isnp(i, jw), ja))then
                             write(lfac, '(F8.4)')ANLIM(k, Isnp(i, jw), ja)
                             LFPR(k, Isnp(i, jw)) = ' N' // lfac
                         elseif(limit==ASLIM(k, Isnp(i, jw), ja))then
                             write(lfac, '(F8.4)')ASLIM(k, Isnp(i, jw), ja)
                             LFPR(k, Isnp(i, jw)) = ' S' // lfac
                         elseif(limit==ALLIM(k, Isnp(i, jw), ja))then
                             write(lfac, '(F8.4)')ALLIM(k, Isnp(i, jw), ja)
                             LFPR(k, Isnp(i, jw)) = ' L' // lfac
                         endif
                     enddo
                 enddo
                 write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,A,I0,A/)')month, gday, &
                      &',', year, '    Julian Date', INT(Jday), ' days ',      &
                     & (Jday - INT(Jday))*24.0, ' hours    Algal group ', ja,  &
                      &' limiting factor'
                 write(SNP(jw), '(1X,A,1000I10)')'Layer  Depth',               &
                     & (Isnp(i, jw), i = Iupr, Idpr)
                 do k = KTWB(jw), Kbr
                     write(SNP(jw), '(1X,I4,F8.2,1000A)')k,                    &
                         & DEPTHM(k, DS(BS(jw))),                              &
                         & (LFPR(k, Isnp(i, jw)), i = Iupr, Idpr)
                 enddo
             endif
         enddo
 
!**      Epiphyton nutrient limitations
 
         do je = 1, nep
             if(PRINT_EPIPHYTON(jw, je))then
                 if(new_page)then
                     write(SNP(jw), '("1",11(A/1X))')(TITLE(j), j = 1, 11)
                     nlines = kmx - KTWB(jw) + 14
                 endif
                 nlines = nlines + kmx - KTWB(jw) + 11
                 new_page = nlines>72
                 do i = Iupr, Idpr                                !mlm   6/30/2006
                     do k = KTWB(jw), KB(Isnp(i, jw))             !mlm  6/30/2006
                         limit = MIN(EPLIM(k, Isnp(i, jw), je),                &
                               & ENLIM(k, Isnp(i, jw), je),                    &
                               & ESLIM(k, Isnp(i, jw), je),                    &
                               & ELLIM(k, Isnp(i, jw), je))
                         if(limit==EPLIM(k, Isnp(i, jw), je))then
                             write(lfac, '(F8.4)')EPLIM(k, Isnp(i, jw), je)
                             LFPR(k, Isnp(i, jw)) = ' P' // lfac
                         elseif(limit==ENLIM(k, Isnp(i, jw), je))then
                             write(lfac, '(F8.4)')ENLIM(k, Isnp(i, jw), je)
                             LFPR(k, Isnp(i, jw)) = ' N' // lfac
                         elseif(limit==ESLIM(k, Isnp(i, jw), je))then
                             write(lfac, '(F8.4)')ESLIM(k, Isnp(i, jw), je)
                             LFPR(k, Isnp(i, jw)) = ' S' // lfac
                         elseif(limit==ELLIM(k, Isnp(i, jw), je))then
                             write(lfac, '(F8.4)')ELLIM(k, Isnp(i, jw), je)
                             LFPR(k, Isnp(i, jw)) = ' L' // lfac
                         endif
                     enddo
                 enddo
                 write(SNP(jw), '(/1X,3(A,1X,I0),A,F0.2,A,I0,A/)')month, gday, &
                      &',', year, '    Julian Date', INT(Jday), ' days ',      &
                     & (Jday - INT(Jday))*24.0, ' hours    Epiphyton group ',  &
                     & je, ' limiting factor'
                 write(SNP(jw), '(1X,A,1000I10)')'Layer  Depth',               &
                     & (Isnp(i, jw), i = Iupr, Idpr)
                 do k = KTWB(jw), Kbr
                     write(SNP(jw), '(1X,I4,F8.2,1000A10)')k,                  &
                         & DEPTHM(k, DS(BS(jw))),                              &
                         & (LFPR(k, Isnp(i, jw)), i = Iupr, Idpr)
                 enddo
             endif
         enddo
     endif
     end subroutine OUTPUT
