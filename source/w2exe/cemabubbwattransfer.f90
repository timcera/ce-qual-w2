!*==cemabubbwattransfer.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
 
     subroutine CEMABUBBWATTRANSFER
 
     use MAIN
     use GLOBAL
     use GEOMC
     use SCREENC
     use KINETIC
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real(8) :: bubbdisssrcsnk, bubsedt, eqbdissconcentration
     integer :: bubblnumber, ngas, ngasconst, nrelarr
     real(8), allocatable, dimension(:) :: henry, kvalue, mw
!
!*** End of declarations rewritten by SPAG
!
 
     allocate(henry(numgas), kvalue(numgas), mw(numgas))
 
     henry(1) = henryconst_h2s          !L atm/M  H2S
     henry(2) = henryconst_ch4          !L atm/M  CH4
     henry(3) = henryconst_nh3          !L atm/M  NH3
     henry(4) = henryconst_co2          !L atm/M  CO2
     bubsedt = 273.15 + SD_T(1)          !K
     mw(1) = 36.                         !H2S gm/mol
     mw(2) = 16.                         !CH4 gm/mol
     mw(3) = 17.                         !NH3 gm/mol
     mw(4) = 44.                         !CO2 gm/mol
 
     do jw = 1, nwb
         kt = KTWB(jw)
         do jb = BS(jw), BE(jw)
             iu = CUS(jb)
             id = DS(jb)
             do segnumi = iu, id
                 do nrelarr = 1, numbubrelarr
                     if(BUBBLESSTATUS(segnumi, nrelarr)==1)then
 
                         do ngas = 1, numgas
!                            Do nGas = 1, NumGas-1 ! debug
 
                             if(ngas==1)ngasconst = nh2s  ! cb 2/18/13
                             if(ngas==2)ngasconst = nch4
                             if(ngas==3)ngasconst = nso4
                             if(ngas==4)ngasconst = ntic
 
                             bubblnumber = BUBBLESLNUMBER(segnumi, nrelarr)
                            !KValue(1)        = Henry(1)/GasConst_R/(T1(BubbLNumber,SegNumI) + 273.15)
                            !KValue(2)        = Henry(2)/GasConst_R/(T1(BubbLNumber,SegNumI) + 273.15)
                            !KValue(3)        = Henry(3)/GasConst_R/(T1(BubbLNumber,SegNumI) + 273.15)
                            !KValue(4)        = Henry(4)/GasConst_R/(T1(BubbLNumber,SegNumI) + 273.15)
                             kvalue(ngas) = henry(ngas)                        &
                               & /gasconst_r/(T1(bubblnumber, segnumi)         &
                               & + 273.15)
 
                             eqbdissconcentration = BUBBLESGASCONC(segnumi,    &
                               & nrelarr, ngas)/kvalue(ngas)
                            !BubbDissSrcSnk = BubbWatGasExchRate*(EqbDissConcentration - C1(BubbLNumber,SegNumI,1+nGas))    !g/m³/s
                             bubbdisssrcsnk = bubbwatgasexchrate*              &
                               & (eqbdissconcentration -                       &
                               & C1(bubblnumber, segnumi, ngasconst))             !g/m³/s  cb 2/18/13
                            !CGSS(BubbLNumber,SegNumI,nGasconst) = CGSS(BubbLNumber,SegNumI,nGasconst) + BubbDissSrcSnk    !BubbDissSrcSnk > 0 Bubbles --> Water
                             C1(bubblnumber, segnumi, ngasconst)               &
                               & = C1(bubblnumber, segnumi, ngasconst)         &
                               & + bubbdisssrcsnk                                 !BubbDissSrcSnk > 0 Bubbles --> Water
                             BUBBLESGASCONC(segnumi, nrelarr, ngas)            &
                               & = BUBBLESGASCONC(segnumi, nrelarr, ngas)      &
                               & - bubbdisssrcsnk*dlt                             !BubbDissSrcSnk > 0 Bubbles --> Water
                             if(BUBBLESGASCONC(segnumi, nrelarr, ngas)<0.D0)   &
                              & BUBBLESGASCONC(segnumi, nrelarr, ngas) = 0.D0
 
                         enddo !nGas
 
                     endif
                 enddo
             enddo
         enddo
     enddo
 
     end subroutine CEMABUBBWATTRANSFER
