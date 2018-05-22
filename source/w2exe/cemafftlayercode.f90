!*==cemafftlayercode.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     subroutine CEMAFFTLAYERCODE
 
 
     use MAIN
     use GLOBAL
     use GEOMC
     use SCREENC
     use RSTART
     use PREC
     use EDDY
     use LOGICC
     use TVDC
     use KINETIC
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: itemp
!
!*** End of declarations rewritten by SPAG
!
 
 
     if(firsttimeinfftcode)then
         firsttimeinfftcode = .FALSE.
         fftactive = .TRUE.
         do jw = 1, nwb
             do jb = BS(jw), BE(jw)
                 do i = CUS(jb), DS(jb)
                    !JAC = NSSS
                     jac = nmft     ! cb 2/18/13
                     k = KB(i)
                     C1(k, i, jac) = initfftlayerconc
                     C1S(k, i, jac) = initfftlayerconc
                     C2(k, i, jac) = initfftlayerconc
                 enddo
             enddo
         enddo
         return
     endif
 
 
     do itemp = 1, numfftactiveprds
 
         if(.NOT.fftactive)then
             if(jday>FFTACTPRDST(itemp) .AND. jday<FFTACTPRDEN(itemp))then
 
                 fftactive = .TRUE.
                 do jw = 1, nwb
                     do jb = BS(jw), BE(jw)
                         do i = CUS(jb), DS(jb)
                            !JAC = NSSS
                             jac = nmft ! cb 2/18/13
                             k = KB(i)
                             C1(k, i, jac) = FFTLAYCONC(i)
                             C1S(k, i, jac) = FFTLAYCONC(i)
                             C2(k, i, jac) = FFTLAYCONC(i)
                         enddo
                     enddo
                 enddo
                 fftactprd = itemp
                 exit
             endif
         endif
 
     enddo
 
     if(fftactive)then
         itemp = fftactprd
         if(jday>FFTACTPRDEN(itemp))then
             fftactive = .FALSE.
             do jw = 1, nwb
                 do jb = BS(jw), BE(jw)
                     do i = CUS(jb), DS(jb)
                        !JAC = NSSS
                         jac = nmft ! cb 2/18/13
                         k = KB(i)
                         FFTLAYCONC(i) = C1(k, i, jac)
                         C1(k, i, jac) = 0.D00
                         C1S(k, i, jac) = 0.D00
                         C2(k, i, jac) = 0.D00
                     enddo
                 enddo
             enddo
         endif
     endif
 
     return
 
     entry MOVEFFTLAYERCONSOLID
 
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
                 if(CEMALAYERADDED(segnumi))then
                    !JAC = NSSS
                     jac = nmft
                              ! cb 2/18/13
                     do k = KB(segnumi) - 1, kt, -1
                         C1(k + 1, segnumi, jac) = C1(k, segnumi, jac)
                         C1S(k + 1, segnumi, jac) = C1S(k, segnumi, jac)
                         C2(k + 1, segnumi, jac) = C2(k, segnumi, jac)
                     enddo !K
                     C1(kt, segnumi, jac) = 0.D00
                     C1S(kt, segnumi, jac) = 0.D00
                     C2(kt, segnumi, jac) = 0.D00
                 endif
             enddo !SegNumI
         enddo !JB
     enddo !JW
 
     end subroutine CEMAFFTLAYERCODE
