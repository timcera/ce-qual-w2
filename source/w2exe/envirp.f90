!*==envirp.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
!*******************************************************************
!**           S U B R O U T I N E   E N V I R P
!*******************************************************************
     subroutine ENVIRP
 
 
 
 
 
 
 
     use GLOBAL
     use MAIN
     use NAMESC
     use SCREENC, ONLY:nit, jday
     use TVDC, ONLY:constituents
     use RSTART, ONLY:ELTM
     use GEOMC, ONLY:DEPTHB
     use ENVIRPMOD
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     character(1), save :: i_int
     integer, save :: n
!
!*** End of declarations rewritten by SPAG
!
 
 
!    initializing variables at first call
     if(nit==1 .OR. iopenfish==0)then
         dltt = dlt
         allocate(cc_e(nct), c_int(nct), c_top(nct), cd_e(ndc), cd_int(ndc),   &
                & cd_top(ndc), c_avg(nct), cd_avg(ndc), cn_e(nct), cdn_e(ndc))
         cc_e = '   '
         c_int = 0.0
         c_top = 0.0
         cd_e = '   '
         cd_int = 0.0
         cd_top = 0.0
         c_avg = 0.0
         cd_avg = 0.0
         cn_e = 0.0
         cdn_e = 0.0
         nac_e = 0
         nacd_e = 0
         cone = nunit
         nunit = nunit + 1
         open(cone, file = 'w2_envirprf.npt', status = 'old')
         read(cone, 9001)i_segint, numclass, selectivec, sjday1, sjday2,       &
                       & istart(1), iend(1),                                   &
                       & (istart(i), iend(i), i = 2, i_segint)
9001     format( // i1, 7x, i8, 5x, a3, f8.0, f8.0, 9(i8, i8))
                                                         ! UP TO A MAX OF 9 INTERVALS
         if(i_segint==0)i_segint = 1
         read(cone, 9002)vel_vpr, vel_int, vel_top, temp_vpr, temp_int,        &
                       & temp_top, depth_vpr, d_int, d_top
9002     format( // 8x, 3(5x, a3, f8.3, f8.3))
         read(cone, 9008)(cc_e(jc), c_int(jc), c_top(jc), jc = 1, nct)
         read(cone, 9008)(cd_e(jd), cd_int(jd), cd_top(jd), jd = 1, ndc)
         close(cone)
 
         do jc = 1, nct
             if(cc_e(jc)==' ON')then
                 nac_e = nac_e + 1
                 cn_e(nac_e) = jc
             endif
         enddo
         do jd = 1, ndc
             if(cd_e(jd)==' ON')then
                 nacd_e = nacd_e + 1
                 cdn_e(nacd_e) = jd
             endif
         enddo
 
         allocate(c_cnt(i_segint, nct), cd_cnt(i_segint, ndc),                 &
                & c_class(i_segint, nct, numclass),                            &
                & cd_class(i_segint, ndc, numclass), c_tot(i_segint, nct),     &
                & cd_tot(i_segint, ndc), t_class(i_segint, numclass),          &
                & v_class(i_segint, numclass), c_sum(nct), cd_sum(ndc))
         allocate(conc_c(nct, numclass), conc_cd(ndc, numclass))
         allocate(d_class(i_segint, numclass))
         allocate(d_tot(i_segint), d_cnt(i_segint), t_tot(i_segint),           &
                & t_cnt(i_segint))
         allocate(v_tot(i_segint), v_cnt(i_segint), volgl(i_segint),           &
                & sumvolt(i_segint))
         cd_cnt = 0.0
         c_cnt = 0.0
         c_class = 0.0
         cd_class = 0.0
         c_tot = 0.0
         cd_tot = 0.0
 
         v_cnt = 0.0
         v_class = 0.0
         v_tot = 0.0
 
         t_cnt = 0.0
         t_class = 0.0
         t_tot = 0.0
 
         d_cnt = 0.0
         d_class = 0.0
         d_tot = 0.0
 
         sumvolt = 0.0
     else
         dltt = (jday - timlast)*86400.
     endif
 
     if(iopenfish/=3)then         !iopenfish=3 is end of simulation deallocate arrays
 
         if(selectivec==' ON')then
             if(jday<sjday1 .OR. jday>sjday2)goto 100
         endif
 
         do n = 1, i_segint
                           ! LOOP OVER SEGMENT INTERVALS BASED ON INPUT DATA
 
             volgl(n) = 0.0
!            start loop for succeeding calls to subroutine
 
             do jw = 1, nwb
                 do jb = BS(jw), BE(jw)
                     do i = CUS(jb), DS(jb)
                         if(selectivec==' ON')then
            !if(i < istart(N) .or. i > iend(N))then
                             if(i<istart(n))cycle
                         ! SW 2/16/2017   exit
                             if(i>iend(n))exit
                         endif
 
!                        Depth
                         if(depth_vpr==' ON')then
                             d_tot(n) = d_tot(n) + DEPTHB(KB(i), i)*dltt
                             d_cnt(n) = d_cnt(n) + dltt
                             d_crit = d_top
                             if(DEPTHB(KB(i), i)>=d_top)d_class(n, 1)          &
                              & = d_class(n, 1) + dltt
                             do jj = 2, numclass
                                 if(DEPTHB(KB(i), i)<d_crit .AND.              &
                                  & DEPTHB(KB(i), i)>=d_crit - d_int)then
                                     d_class(n, jj) = d_class(n, jj) + dltt
                                     exit
                                 else
                                     d_crit = d_crit - d_int
                                 endif
                                 if(jj==numclass .AND. DEPTHB(KB(i), i)        &
                                  & <d_crit + d_int)d_class(n, jj)             &
                                  & = d_class(n, jj) + dltt
                             enddo
                         endif
 
 
 
 
 
                         do k = KTWB(jw), KB(i)
                             volgl(n) = volgl(n) + VOL(k, i)
 
!                            Temperature
 
                             if(temp_vpr==' ON')then
 
                                 t_tot(n) = t_tot(n) + T2(k, i)*VOL(k, i)*dltt
                                 t_cnt(n) = t_cnt(n) + VOL(k, i)*dltt
                                 t_crit = temp_top
                                 if(T2(k, i)>=temp_top)t_class(n, 1)           &
                                  & = t_class(n, 1) + dltt*VOL(k, i)
                                 do jj = 2, numclass
                                     if(T2(k, i)<t_crit .AND. T2(k, i)         &
                                      & >=t_crit - temp_int)then
                                         t_class(n, jj) = t_class(n, jj)       &
                                           & + dltt*VOL(k, i)
                                         exit
                                     else
                                         t_crit = t_crit - temp_int
                                     endif
                                     if(jj==numclass .AND. T2(k, i)            &
                                      & <t_crit + temp_int)t_class(n, jj)      &
                                      & = t_class(n, jj) + dltt*VOL(k, i)
                                 enddo
                             endif
 
 
!                            Velocity
 
                             if(vel_vpr==' ON')then
 
                                 v_tot(n) = v_tot(n) + U(k, i)*VOL(k, i)*dltt
                                 v_cnt(n) = v_cnt(n) + VOL(k, i)*dltt
                                 v_crit = vel_top
                                 if(U(k, i)>=vel_top)v_class(n, 1)             &
                                  & = v_class(n, 1) + dltt*VOL(k, i)
                                 do jj = 2, numclass
                                     if(U(k, i)<v_crit .AND. U(k, i)>=v_crit - &
                                      & vel_int)then
                                         v_class(n, jj) = v_class(n, jj)       &
                                           & + dltt*VOL(k, i)
                                         exit
                                     else
                                         v_crit = v_crit - vel_int
                                     endif
                                     if(jj==numclass .AND. U(k, i)             &
                                      & <v_crit + vel_int)v_class(n, jj)       &
                                      & = v_class(n, jj) + dltt*VOL(k, i)
                                 enddo
                             endif
 
 
!                            Constituents
                             if(constituents)then
                                 do jc = 1, nac_e
                                     jac = cn_e(jc)
 
                                     c_tot(n, jc) = c_tot(n, jc)               &
                                       & + C2(k, i, jac)*CMULT(jac)*VOL(k, i)  &
                                       & *dltt
                                     c_cnt(n, jc) = c_cnt(n, jc) + VOL(k, i)   &
                                       & *dltt
                                     c_crit = c_top(jac)
                                     if(C2(k, i, jac)*CMULT(jac)>=c_top(jac))  &
                                      & c_class(n, jc, 1) = c_class(n, jc, 1)  &
                                      & + dltt*VOL(k, i)
                                     do jj = 2, numclass
                                         if(C2(k, i, jac)*CMULT(jac)           &
                                           & <c_crit .AND. C2(k, i, jac)       &
                                           & *CMULT(jac)>=c_crit - c_int(jac)) &
                                           & then
                                         c_class(n, jc, jj)                    &
                                           & = c_class(n, jc, jj)              &
                                           & + dltt*VOL(k, i)
                                         exit
                                         else
                                         c_crit = c_crit - c_int(jac)
                                         endif
                                         if(jj==numclass .AND. C2(k, i, jac)   &
                                           & *CMULT(jac)<c_crit + c_int(jac))  &
                                           & c_class(n, jc, jj)                &
                                           & = c_class(n, jc, jj)              &
                                           & + dltt*VOL(k, i)
                                     enddo
 
 
 
                                 enddo
 
 
!                                Derived Constituents
 
                                 do jc = 1, nacd_e
 
                                     jacd = cdn_e(jc)
 
                                     cd_tot(n, jc) = cd_tot(n, jc)             &
                                       & + CD(k, i, jacd)*CDMULT(jacd)         &
                                       & *VOL(k, i)*dltt
                                     cd_cnt(n, jc) = cd_cnt(n, jc) + VOL(k, i) &
                                       & *dltt
                                     cd_crit = cd_top(jacd)
                                     if(CD(k, i, jacd)*CDMULT(jacd)            &
                                      & >=cd_top(jacd))cd_class(n, jc, 1)      &
                                      & = cd_class(n, jc, 1) + dltt*VOL(k, i)
                                     do jj = 2, numclass
                                         if(CD(k, i, jacd)*CDMULT(jacd)        &
                                           & <cd_crit .AND. CD(k, i, jacd)     &
                                           & *CDMULT(jacd)>=cd_crit -          &
                                           & cd_int(jacd))then
                                         cd_class(n, jc, jj)                   &
                                           & = cd_class(n, jc, jj)             &
                                           & + dltt*VOL(k, i)
                                         exit
                                         else
                                         cd_crit = cd_crit - cd_int(jacd)
                                         endif
                                         if(jj==numclass .AND. CD(k, i, jacd)  &
                                           & *CDMULT(jacd)<cd_crit +           &
                                           & cd_int(jacd))cd_class(n, jc, jj)  &
                                           & = cd_class(n, jc, jj)             &
                                           & + dltt*VOL(k, i)
                                     enddo
 
                                 enddo
                             endif
                         enddo
                     enddo
                 enddo
             enddo
 
!            sum of volgL*dltt for volume fraction calculation
 
             sumvolt(n) = sumvolt(n) + volgl(n)*dltt
         enddo
          ! END OF INTERVAL FOR SEGMENTS
     endif
 
 
100   if(iopenfish==3)then
!        calculating average violation concentration and writing to file
         do n = 1, i_segint
             cd_sum = 0.0
             c_sum = 0.0
             t_sum = 0.0
             v_sum = 0.0
             write(i_int, '(I1)')n
 
             if(temp_vpr==' ON')then
                 if(t_cnt(n)>0.0)then
                     t_avg = t_tot(n)/t_cnt(n)
                 else
                     t_avg = 0.0
                 endif
 
                 open(cone, file = 'envrprf_t_' // i_int // '.csv',            &
                     &status = 'unknown')
                 write(cone, *)'"Temperature interval,","Fraction of volume"'
                 temp_c = temp_top
                 do i = 1, numclass
                     write(cone, 9009)temp_c, t_class(n, i)/sumvolt(n)
                     temp_c = temp_c - temp_int
                     t_sum = t_sum + t_class(n, i)/sumvolt(n)
                 enddo
                 write(cone, '(1x)')
                 write(cone, '(" 0, ",e12.4)')t_sum
                 write(cone, '(1x)')
                 write(cone, '(" 0, ",e12.4)')t_avg
                 close(cone)
             endif
 
             if(vel_vpr==' ON')then
                 if(v_cnt(n)>0.0)then
                     v_avg = v_tot(n)/v_cnt(n)
                 else
                     v_avg = 0.0
                 endif
                 open(cone, file = 'envrprf_v' // i_int // '.csv',             &
                     &status = 'unknown')
                 write(cone, *)'"Velocity interval,","Fraction of volume"'
                 vel_c = vel_top
                 do i = 1, numclass
                     write(cone, 9009)vel_c, v_class(n, i)/sumvolt(n)
                     vel_c = vel_c - vel_int
                     v_sum = v_sum + v_class(n, i)/sumvolt(n)
                 enddo
                 write(cone, '(1x)')
                 write(cone, '(" 0, ",e12.4)')v_sum
                 write(cone, '(1x)')
                 write(cone, '(" 0, ",e12.4)')v_avg
                 close(cone)
             endif
 
 
             if(depth_vpr==' ON')then
                 if(d_cnt(n)>0.0)then
                     d_avg = d_tot(n)/d_cnt(n)
                 else
                     d_avg = 0.0
                 endif
                 open(cone, file = 'envrprf_depth' // i_int // '.csv',         &
                     &status = 'unknown')
                 write(cone, *)'"Depth interval,","Fraction of time"'
                 d_c = d_top
                 do i = 1, numclass
                     write(cone, 9009)d_c, d_class(n, i)/d_cnt(n)
                     d_c = d_c - d_int
                     d_sum = d_sum + d_class(n, i)/d_cnt(n)
                 enddo
                 write(cone, '(1x)')
                 write(cone, '(" 0, ",e12.4)')d_sum
                 write(cone, '(1x)')
                 write(cone, '(" 0, ",e12.4)')d_avg
                 close(cone)
             endif
 
             if(nac_e>0)then
                 open(cone, file = 'envrprf_c' // i_int // '.csv',             &
                     &status = 'unknown')
                 write(cone, 9003)(CNAME2(cn_e(jc)), jc = 1, nac_e)
9003             format((a8, '_interval, Fraction_of_volume, '))
                 do jc = 1, nac_e
                     if(c_cnt(n, jc)>0.0)then
                         c_avg(jc) = c_tot(n, jc)/c_cnt(n, jc)
                     else
                         c_avg(jc) = 0.0
                     endif
                     do i = 1, numclass
                         if(i==1)then
                             conc_c(jc, i) = c_top(cn_e(jc))
                         else
                             conc_c(jc, i) = conc_c(jc, i - 1)                 &
                               & - c_int(cn_e(jc))
                         endif
                         c_sum(jc) = c_sum(jc) + c_class(n, jc, i)/sumvolt(n)
                     enddo
                 enddo
                 do i = 1, numclass
                     write(cone, 9004)(conc_c(jc, i), c_class(n, jc, i)/sumvolt&
                                    & (n), jc = 1, nac_e)
9004                 format((f10.4, ',', e12.4, ','))
                 enddo
                 write(cone, '(1x)')
                 write(cone, 9011)(c_sum(jc), jc = 1, nac_e)
                 write(cone, '(1x)')
                 write(cone, 9011)(c_avg(jc), jc = 1, nac_e)
                 close(cone)
 
                 if(nacd_e>0)then
                     open(cone, file = 'envrprf_cd' // i_int // '.csv',        &
                        & status = 'unknown')
                     write(cone, 9005)(CDNAME2(cdn_e(jc)), jc = 1, nacd_e)
9005                 format((a8, '_interval, Fraction_of_volume,'))
                     do jc = 1, nacd_e
                         if(cd_cnt(n, jc)>0.0)then
                             cd_avg(jc) = cd_tot(n, jc)/cd_cnt(n, jc)
                         else
                             cd_avg(jc) = 0.0
                         endif
 
                         do i = 1, numclass
                             if(i==1)then
                                 conc_cd(jc, i) = cd_top(cdn_e(jc))
                             else
                                 conc_cd(jc, i) = conc_cd(jc, i - 1)           &
                                   & - cd_int(cdn_e(jc))
                             endif
                             cd_sum(jc) = cd_sum(jc) + cd_class(n, jc, i)      &
                               & /sumvolt(n)
                         enddo
                     enddo
                     do i = 1, numclass
                         write(cone, 9006)                                     &
                             & (conc_cd(jc, i), cd_class(n, jc, i)/sumvolt(n), &
                             & jc = 1, nacd_e)
 
9006                     format((f10.4, ",", e12.4, ","))
                     enddo
                     write(cone, '(1x)')
                     write(cone, 9010)(cd_sum(jc), jc = 1, nacd_e)
                     write(cone, '(1x)')
                     write(cone, 9010)(cd_avg(jc), jc = 1, nacd_e)
                     close(cone)
                 endif
             endif
         enddo
         deallocate(c_cnt, cd_cnt, c_class, cd_class, c_tot, cd_tot, t_class,  &
                  & v_class, c_sum, cd_sum)
         deallocate(conc_c, conc_cd)
         deallocate(cc_e, c_int, c_top, cd_e, cd_int, cd_top, c_avg, cd_avg,   &
                  & cn_e, cdn_e)
         deallocate(t_tot, t_cnt, v_tot, v_cnt, d_tot, d_cnt, d_class, volgl,  &
                  & sumvolt)
9007     format((f10.4, ',', e12.4, ','))
 
     endif
 
     timlast = jday
 
9008  format( // (8x, (5x, a3, f8.0, f8.0)))
9009  format((f6.2, ',', e12.4, ','))
9010  format((" 0, ", e12.4, ","))
9011  format((" 0, ", e12.4, ","))
 
     end subroutine ENVIRP
