!*==gasbubblesformation.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     subroutine GASBUBBLESFORMATION(Radius, Deltat, Volume)
 
 
    ! Type declarations
     use GLOBAL
     use SCREENC
     use CEMAVARS
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     real(8) :: Deltat, Radius, Volume
     intent (in) Deltat, Volume
     intent (inout) Radius
!
! Local variables
!
     real(8) :: bubsedt, c0t, c1t, cgt, ctt, diffvolume, dismass, gasmass,     &
              & nbubbles, nbubblesp, nbubblost, netmass, p0, pbubb, pbubbt,    &
              & pcrit, porosity, ro, rsi, source, sourcet, vbub, vbubbles
     real(8), allocatable, dimension(:) :: c0b, c1b, cgb, ctot, henry, k, mw,  &
          & source0
     logical :: foundopenarray
     integer :: ngas, nrelarr, temp
!
!*** End of declarations rewritten by SPAG
!
 
 
     allocate(ctot(numgas), cgb(numgas), c0b(numgas), c1b(numgas))
                                  ! nTry
     allocate(source0(numgas), henry(numgas), k(numgas), mw(numgas))
 
     porosity = BEDPOROSITY(segnumi)
     henry(1) = henryconst_h2s          !L atm/M  H2S
     henry(2) = henryconst_ch4          !L atm/M  CH4
     henry(3) = henryconst_nh3          !L atm/M  NH3
     henry(4) = henryconst_co2          !L atm/M  CO2
     bubsedt = 273.15 + SD_T(1)          !K           ! cb 5/22/15
     k(1) = henry(1)/gasconst_r/bubsedt
     k(2) = henry(2)/gasconst_r/bubsedt
     k(3) = henry(3)/gasconst_r/bubsedt
     k(4) = henry(4)/gasconst_r/bubsedt
     ro = 0.D0              !m
     rsi = 8318.78          !l-N/m²/mol/K
     p0 = 9800.              !N/m²
     mw(1) = 36.             !H2S gm/mol
     mw(2) = 16.             !CH4 gm/mol
     mw(3) = 17.             !NH3 gm/mol
     mw(4) = 44.             !CO2 gm/mol
     if(CRACKOPEN(segnumi))nbubblost = MFTBUBBRELEASED(segnumi)
 
     do ngas = 1, numgas
         ctot(ngas) = TCONC(ngas, segnumi)
         source0(ngas) = SCONC(ngas, segnumi)
     enddo !nGas
 
     if(firsttimeinbubbles)then
 
         cgt = 0.D0
         c1t = 0.D0
         c0t = 0.D0
         ctt = 0.D0
         do ngas = 1, numgas
             c0b(ngas) = ctot(ngas)/(1 + k(ngas))
             cgb(ngas) = c0b(ngas)*k(ngas)
             c1b(ngas) = ctot(ngas)
             cgt = cgt + cgb(ngas)
             c1t = c1t + c1b(ngas)
             c0t = c0t + c0b(ngas)
             ctt = ctt + ctot(ngas)
         enddo !nGas
         Radius = SQRT(2.0*porosity*gasdiff_sed*Deltat*(c1t - c0t)/cgt + ro**2)
         vbub = (4.0/3.0)*3.1415927*Radius**3
         netmass = ctt*Volume*porosity
         dismass = c0t*Volume*porosity
         gasmass = netmass - dismass
         nbubbles = gasmass/(vbub*cgt)
         pcrit = 1.32*(critstressif**6/(youngmodulus*nbubbles*vbub))**0.2 + p0
         pbubbt = 0.D0
         do ngas = 1, numgas
             pbubb = cgb(ngas)*rsi*0.001*bubsedt/mw(ngas)
             pbubbt = pbubbt + pbubb
         enddo !nGas
 
     else
 
         sourcet = 0.D0
         cgt = 0.D0
         c1t = 0.D0
         c0t = 0.D0
         ctt = 0.D0
         do ngas = 1, numgas
 
             source = source0(ngas)
             sourcet = sourcet + source
             ctot(ngas) = ctot(ngas) + source*Deltat
             c0b(ngas) = ctot(ngas)/(1 + k(ngas))
             cgb(ngas) = c0b(ngas)*k(ngas)
             c1b(ngas) = c0b(ngas)
             cgt = cgt + cgb(ngas)
             c1t = c1t + c1b(ngas)
             c0t = c0t + c0b(ngas)
             ctt = ctt + ctot(ngas)
 
         enddo !nGas
 
         Radius = Radius + porosity*gasdiff_sed/(Radius*cgt)                   &
                & *(sourcet*calibparam_r1**2/(6*gasdiff_sed) + (c1t - c0t))    &
                & *Deltat
 
         if(limbubbsize)then
             if(Radius>maxbubbrad/1000.0)Radius = maxbubbrad/1000.0
         endif
 
 
     endif
 
     vbub = (4.0/3.0)*3.1415927*Radius**3
     netmass = 0.D0
     dismass = 0.D0
     cgt = 0.D0
     do ngas = 1, numgas
         netmass = netmass + ctot(ngas)*Volume*porosity
         dismass = dismass + c0b(ngas)*Volume*porosity
         cgt = cgt + cgb(ngas)
     enddo !nGas
     gasmass = netmass - dismass
     nbubbles = gasmass/(vbub*cgt)
 
     pcrit = 1.324*(critstressif**6/(youngmodulus*nbubbles*vbub))**0.2 + p0
     pbubbt = 0.D0
     do ngas = 1, numgas
         pbubb = cgb(ngas)*rsi*0.001*bubsedt/mw(ngas)
         pbubbt = pbubbt + pbubb
     enddo !nGas
 
     PRESBUBBSED(segnumi) = pbubbt
     PRESCRITSED(segnumi) = pcrit
 
     if(pbubbt<pcrit*crackclosefraction)then
         CRACKOPEN(segnumi) = .FALSE.
         nbubblost = 0
     endif
 
     if(LASTDIFFVOLUME(segnumi)<0)LASTDIFFVOLUME(segnumi) = 0.D00
 
     if(pbubbt>pcrit .AND. .NOT.CRACKOPEN(segnumi))then
 
         CRACKOPEN(segnumi) = .TRUE.
         vbubbles = critstressif**6/(youngmodulus*((pbubbt - p0)/1.32)**5)
         diffvolume = vbub*nbubbles - vbubbles
         diffvolume = diffvolume*bubbrelscale
         if(usereleasefraction)then
             vbubbles = critstressif**6/                                       &
                      & (youngmodulus*((pcrit*crackclosefraction)/1.32)**5)
             diffvolume = bubbrelfraction*(vbub*nbubbles - vbubbles)
         endif
         LASTDIFFVOLUME(segnumi) = diffvolume
         nbubblost = (diffvolume)/vbub
         nbubblesp = nbubbles
         nbubbles = nbubbles - nbubblost
         do ngas = 1, numgas
             cgb(ngas) = cgb(ngas)*(vbub*nbubblesp - diffvolume)               &
                       & /(vbub*nbubblesp)
             c0b(ngas) = cgb(ngas)/k(ngas)
             ctot(ngas) = c0b(ngas)*(1 + k(ngas))
         enddo !nGas
 
     endif
 
!!!!!!!!!!!!!!!!! debug
!    CrackOpen(SegNumI)=.false.
!!!!!!!!!!!!!!!!!!!!!!! debug
     if(CRACKOPEN(segnumi))then
 
         nbubblesp = nbubbles
        !Nbubbles = Nbubbles - NbubbLost  ! cb 2/21/13
         do ngas = 1, numgas
             cgb(ngas) = cgb(ngas)*(vbub*nbubblesp - LASTDIFFVOLUME(segnumi))  &
                       & /(vbub*nbubblesp)
             c0b(ngas) = cgb(ngas)/k(ngas)
             ctot(ngas) = c0b(ngas)*(1 + k(ngas))
 
             TCONCP(ngas, segnumi) = ctot(ngas)
             TCONC(ngas, segnumi) = ctot(ngas)
 
         enddo !nGas
 
         temp = INT4(nbubblost)
 
         foundopenarray = .FALSE.
         MFTBUBBRELEASED(segnumi) = KIDINT(nbubblost)
        !nTry = 0
         do nrelarr = 1, numbubrelarr
            !nTry = nTry + 1
             if(BUBBLESSTATUS(segnumi, nrelarr)==0)then
                 foundopenarray = .TRUE.
                 BUBBLESSTATUS(segnumi, nrelarr) = 1
                 BUBBLESCARRIED(segnumi, nrelarr) = MFTBUBBRELEASED(segnumi)
                 BUBBLESRADIUS(segnumi, nrelarr) = Radius
                 BUBBLESLNUMBER(segnumi, nrelarr) = KB(segnumi)
                 do ngas = 1, numgas
                     BUBBLESGASCONC(segnumi, nrelarr, ngas) = cgb(ngas)
                 enddo
                 exit
             endif
         enddo
         if(.NOT.foundopenarray)then
             write(cemalogfiln, *)                                             &
                      &"Insufficient array size for bubbles release at JDAY = "&
                     & , jday
             stop
         endif
     endif
 
     CGSED(segnumi) = cgt
     C0SED(segnumi) = c0t
     CTSED(segnumi) = ctt
 
     end subroutine GASBUBBLESFORMATION
