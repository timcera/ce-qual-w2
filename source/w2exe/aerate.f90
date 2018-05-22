!*==aerate.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
     subroutine AERATE
     use GLOBAL
     use MAIN
     use KINETIC
     use TRANS
     use SCREENC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     real, allocatable, dimension(:), save :: actual_mass, atimoff, atimon,    &
       & cumdomass, dooff, doon, dzmulta, smass
     logical, allocatable, dimension(:), save :: aerateo2
     character(16), save :: conaer
     real, allocatable, dimension(:, :), save :: dzmult
     integer, allocatable, dimension(:), save :: iaseg, iprb, kbota, kprb,     &
          & ktopa
     integer, save :: kbot, ktop, naer, nlayers
!
!*** End of declarations rewritten by SPAG
!
 
     allocate(dzmult(kmx, imx))
     dzmult = 1.0
     open(aeratefn, file = 'W2_AERATE.NPT', status = 'OLD')
     read(aeratefn, '(//I8,A16)')naer, conaer
 
     if(naer==0)naer = 1
     allocate(dzmulta(naer), iaseg(naer), ktopa(naer), kbota(naer), smass(naer)&
            & , atimon(naer), atimoff(naer), dooff(naer), doon(naer),          &
            & iprb(naer), kprb(naer), actual_mass(naer), aerateo2(naer),       &
            & cumdomass(naer))
     dzmulta = 0.0
     actual_mass = 0.0
     cumdomass = 0.0
     read(aeratefn, 9001)
9001  format(/)
     do i = 1, naer
         read(aeratefn, '(I8,I8,I8,F8.0,F8.0,F8.0,3F8.0,2I8)')iaseg(i),        &
            & ktopa(i), kbota(i), smass(i), atimon(i), atimoff(i), dzmulta(i), &
            & dooff(i), doon(i), iprb(i), kprb(i)
     enddo
 
     close(aeratefn)
 
     if(restart_in)then
         open(aeratefn, file = conaer, position = 'APPEND')
         jday1 = 0.0
         rewind(aeratefn)
         read(aeratefn, '(//)')
         do while (jday1<jday)
             read(aeratefn, '(F9.0)', end = 50)jday1
         enddo
         backspace(aeratefn)
50       jday1 = 0.0
     else
         open(aeratefn, file = conaer, status = 'UNKNOWN')
         write(aeratefn, '(A,I3,A)')'OUTPUT FILE FOR AERATION INPUT WITH',     &
                                  & naer, ' INPUT(S).'
         write(aeratefn, '(A114)')                                             &
     &'JDAY,INSTMASSRATE#1(KGO2/D),CUMMASS#1(KGO2),DOPROBE#1(MG/L),INSTMASSRATE&
     &#2(KGO2/D),CUMMASS#2(KGO2),DOPROBE#2(MG/L)'
     endif
     dzmult = 1.0
              ! ALWAYS RESET MIXING COEFFICIENT FOR AERATION SW IPC 2/01/01
     return
 
     entry DZAERATE
                  ! FROM W2_ MAIN CODE
     do i = 1, naer
         dz(ktopa(i):kbota(i), iaseg(i)) = dz(ktopa(i):kbota(i), iaseg(i))     &
           & *dzmult(ktopa(i):kbota(i), iaseg(i))
     enddo
     return
 
     entry AERATEMASS
                    ! FROM WQ_CONSTITUENTS
 
!    SECTION FOR HYPOLIMNETIC AERATION
     dzmult = 1.0
              ! ALWAYS RESET TO 1.0 IN CASE NO AERATION
     do ii = 1, naer
         if(jday>=atimon(ii) .AND. jday<=atimoff(ii))then
 
           ! FIND BRANCH AND WATERBODY FOR ISEG
             do jjb = 1, nbr
                 if(BR_INACTIVE(jjb))cycle     ! SW 6/12/2017
                 if(iaseg(ii)>=US(jjb) .AND. iaseg(ii)<=DS(jjb))exit
             enddo
 
             if(jjb==jb)then          ! IF THIS BRANCH DOESN'T HAVE AERATION SKIP IT
 
                 ktop = MAX(KTWB(jw), ktopa(ii))
                 kbot = MIN(KB(iaseg(ii)), kbota(ii))
 
                 nlayers = kbot - ktop + 1
 
                 aerateo2(ii) = .TRUE.
 
                 if(O2(kprb(ii), iprb(ii))>doon(ii) .AND.                      &
                  & O2(kprb(ii), iprb(ii))>dooff(ii))then
                     aerateo2(ii) = .FALSE.
                     actual_mass(ii) = 0.0
                 endif
 
                 if(aerateo2(ii))then
                     actual_mass(ii) = smass(ii) ! SW 12/28/01
                     cumdomass(ii) = cumdomass(ii) + smass(ii)*dlt/86400.
                     do k = ktop, kbot
 
                         dzmult(k, iaseg(ii)) = dzmulta(ii)
                                                          ! THIS MEANS THE INCREASE IN DZ IS LAGGED ONE TIME STEP
                         CSSB(k, iaseg(ii), ndo) = CSSB(k, iaseg(ii), ndo)     &
                           & + smass(ii)/(86.4*REAL(nlayers))
!                        UNITS OF SMASS ARE IN KG/DAY
!                        TYPICAL UNITS OF CSSB: (MG/L)*(M3/S) TO OONVERT
 
!                        MULTIPLY BY
!                        (1KG/10^6MG)*(1000L/1M3)*(86400S/DAY)==86.4
                     enddo
                 endif
             endif
         else
             aerateo2(ii) = .FALSE.
             actual_mass(ii) = 0.0
         endif
 
     enddo
 
     return
 
     entry AERATEOUTPUT
                       ! FROM W2_MAIN CODE
 
     write(aeratefn, '(F9.3,<NAER>(1X,E12.3,1X,E12.3,1X,F8.3))')jday,          &
         & (actual_mass(ii), cumdomass(ii), O2(kprb(ii), iprb(ii)), ii = 1,    &
         & naer)
 
     return
     entry DEALLOCATE_AERATE
     deallocate(dzmulta, iaseg, ktopa, kbota, smass, atimon, atimoff, dooff,   &
              & doon, iprb, kprb, actual_mass, aerateo2, cumdomass)
     deallocate(dzmult)
     close(aeratefn)
     end subroutine AERATE
