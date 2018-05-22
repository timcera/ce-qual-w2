!*==read_fish_data.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
! ************************ READ FISH INPUT DATA FILES
 
     subroutine READ_FISH_DATA                                             !FISH
 
 
 
     use FISHY
     use SCREENC, ONLY:jday
     use MAIN, ONLY:fish_particle_exist
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     character(3) :: aline
     real, allocatable, dimension(:) :: delaydateinit, fxlocinit, fzlocinit,   &
                                      & sedvelinit
     integer :: htstbot, htstside, i, idebug, ilinear, j, na, nfishlast, ng
!
!*** End of declarations rewritten by SPAG
!
 
 
     open(diagfn, file = 'particle.csv', status = 'old')
     read(diagfn, *)
     read(diagfn, *)
     read(diagfn, *)parton, nfishseg, nfishpcel, aline, dxtheory, outfreqp,    &
                  & ilinear, htstbot, htstside, idebug                                                !,ALPHAX,ALPHAZ     !'(/10x,7x,a3,i10,i10,f10.0,(7x,a3),1f10.0,7X,A3)
     if(parton=='ON')then
         particle = .TRUE.
     else
         particle = .FALSE.
         fish_particle_exist = .FALSE.
                                 ! STOPS ALL FURTHER
         return
     endif
     if(aline=='ON')then
         line = .TRUE.
     else
         line = .FALSE.
     endif
     if(outfreqp==0.0)outfreqp = 1.0
                                    ! error trapping
 
!    LINE: If LINE, distribute NFISHCEL not at a point but linearly along the
!    line defined by IFISH and IFISHT and IFISHB; IF LINE, FZLOC=0 PARTICLE =
!    ON only "dumb" particle transport with random fluid motion +
!    Sedimentation Vel (SEDVEL) NFISHPCEL: # of particles added per cel
!    NFISHSEG: # of segements to add particles to DATE:  Date particles are
!    deposited
     if(nfishseg==0)nfishseg = 1
     allocate(ifish(nfishseg), ifisht(nfishseg), ifishb(nfishseg),             &
            & delaydateinit(nfishseg), fxlocinit(nfishseg), fzlocinit(nfishseg)&
            & , sedvelinit(nfishseg))
     read(diagfn, *)
     do i = 1, nfishseg
         read(diagfn, *)na, ifish(i), ifisht(i), ifishb(i), fxlocinit(i),      &
                      & fzlocinit(i), sedvelinit(i), delaydateinit(i)
 
         if(delaydateinit(i)<jday)delaydateinit(i) = jday
     enddo
 
    !READ(DIAGFN,'(//(10X,9I10))')(IFISH(I),I=1,NFISHSEG)
    !READ (DIAGFN,'(//(10X,9I10))')(IFISHT(I),I=1,NFISHSEG)
    !READ (DIAGFN,'(//(10X,9I10))')(IFISHB(I),I=1,NFISHSEG)
    ! IFISH: SEG#s where fish are added
    ! IFISHT:Uppermost top Layer where Fish are added in Segment IFISH
    ! IFISHB: Bottommost Layer where Fish are added in Segment IFISH
     nfish = 0
     do i = 1, nfishseg
         nfish = nfish + (ifishb(i) - ifisht(i) + 1)*nfishpcel
     enddo
 
     allocate(group(nfish), fxloci(nfish), fzloci(nfish), delaydate(nfish),    &
            & sedvel(nfish))
 
     nfishlast = 1
     do i = 1, nfishseg
         do j = nfishlast, (nfishlast - 1) + (ifishb(i) - ifisht(i) + 1)       &
           & *nfishpcel
             group(j) = i
             fzloci(j) = fzlocinit(i)
             fxloci(j) = fxlocinit(i)
             delaydate(j) = delaydateinit(i)
             sedvel(j) = sedvelinit(i)/86400.
                                          ! convert m/d to m/s
         enddo
         nfishlast = j
     enddo
 
!    READ(DIAGFN,'(//3F10.0,I10)')FXLOC,FZLOC,OUTFREQP,IDEBUG
 
!    ! FXLOC   = Location of fish within segment IMP from upstream side
!    ! FZLOC   = Location of fish within layer KMP from top side
!    ! FNBP    = Branch where fish is released
 
     xrefl = 0.1                       ! Reflect this percentage of segment length
                                       !   when fish encounters horizontal obstacle
     yrefl = 0.1                       ! Reflect this % of cell width when encounter stream bank
     zbotrefl = 2.0                    ! Reflect this percentage of bottom layer height
                                       !   when fish encounters bottom
     zsurrefl = 0.5                    ! Reflect this percentage of surface layer height
                                       !   when fish encounters water surface
 
 
     debug = .FALSE.
     linear = .FALSE.
     hitstickbottom = .FALSE.
     hitstickside = .FALSE.
     if(htstbot==1)hitstickbottom = .TRUE.
     if(htstside==1)hitstickside = .TRUE.
     if(ilinear==1)linear = .TRUE.
     if(idebug==1)debug = .TRUE.
 
     close(diagfn)
 
     end subroutine READ_FISH_DATA
