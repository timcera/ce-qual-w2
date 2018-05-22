!*==gregorian_date.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!**                                     S U B R O U T I N E    G R E G O R I A N   D A T E                                        **
!***********************************************************************************************************************************
 
     subroutine GREGORIAN_DATE
     use GDAYC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     integer :: incr
!
!*** End of declarations rewritten by SPAG
!
 
!    Determine if new year (regular or leap) and increment year
 
     do while (jdayg>=366)
         if(.NOT.leap_year .AND. jdayg>=366)then
             jdayg = jdayg - 365
             year = year + 1
             leap_year = MOD(year, 4)==0
         elseif(jdayg>=367)then
             jdayg = jdayg - 366
             year = year + 1
             leap_year = MOD(year, 4)==0
         else
             exit
         endif
     enddo
     incr = 0
     if(leap_year)incr = 1
 
!    Determine month and day of year
 
     if(jdayg>=1 .AND. jdayg<32)then
         gday = jdayg
         daym = 31.0
         month = '  January'
         imon = 1
     elseif(jdayg>=32 .AND. jdayg<60 + incr)then
         gday = jdayg - 31
         daym = 29.0
         month = ' February'
         imon = 2
     elseif(jdayg>=60 .AND. jdayg<91 + incr)then
         gday = jdayg - 59 - incr
         daym = 31.0
         month = '    March'
         imon = 3
     elseif(jdayg>=91 .AND. jdayg<121 + incr)then
         gday = jdayg - 90 - incr
         daym = 30.0
         month = '    April'
         imon = 4
     elseif(jdayg>=121 .AND. jdayg<152 + incr)then
         gday = jdayg - 120 - incr
         daym = 31.0
         month = '      May'
         imon = 5
     elseif(jdayg>=152 .AND. jdayg<182 + incr)then
         gday = jdayg - 151 - incr
         daym = 30.0
         month = '     June'
         imon = 6
     elseif(jdayg>=182 .AND. jdayg<213 + incr)then
         gday = jdayg - 181 - incr
         daym = 31.0
         month = '     July'
         imon = 7
     elseif(jdayg>=213 .AND. jdayg<244 + incr)then
         gday = jdayg - 212 - incr
         daym = 31.0
         month = '   August'
         imon = 8
     elseif(jdayg>=244 .AND. jdayg<274 + incr)then
         gday = jdayg - 243 - incr
         daym = 30.0
         month = 'September'
         imon = 9
     elseif(jdayg>=274 .AND. jdayg<305 + incr)then
         gday = jdayg - 273 - incr
         daym = 31.0
         month = '  October'
         imon = 10
     elseif(jdayg>=305 .AND. jdayg<335 + incr)then
         gday = jdayg - 304 - incr
         daym = 30.0
         month = ' November'
         imon = 11
     elseif(jdayg>=335 .AND. jdayg<366 + incr)then
         gday = jdayg - 334 - incr
         daym = 31.0
         month = ' December'
         imon = 12
     endif
     end subroutine GREGORIAN_DATE
