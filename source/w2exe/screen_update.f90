!*==screen_update.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
 
!***********************************************************************************************************************************
!*                                                E N T R Y   S C R E E N  U P D A T E                                            **
!***********************************************************************************************************************************
 
     subroutine SCREEN_UPDATE(Dlg)
     use DFLOGM
     use DFLIB
     use MSCLIB
     use GEOMC
     use GLOBAL
     use GDAYC
     use SCREENC
     use SURFHE
     use TVDC
     use LOGICC
     use NAMESC
     use STRUCTURES
     use MAIN, ONLY:tmstrt, tmend
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Dummy arguments
!
     type(DIALOG) :: Dlg
!
! Local variables
!
     integer :: iprog, result
     character(3000) :: text1
     character(8) :: time
!
!*** End of declarations rewritten by SPAG
!
     time = cctime(1:2) // ':' // cctime(3:4) // ':' // cctime(5:6)
     call DATE_AND_TIME(cdate, cctime)
     call CPU_TIME(current)
     write(text1, '(A," ",I0,", ",I0)')month, gday, year
     result = DLGSET(Dlg, gregorian_day, text1)
     write(text1, '(I0)')INT(jday)
     result = DLGSET(Dlg, julian_day, text1)
     write(text1, '(F0.2)')(jday - INT(jday))*24.0
     result = DLGSET(Dlg, julian_hour, text1)
     write(text1, '(I0)')INT(eltmjd)
     result = DLGSET(Dlg, elapsed_day, text1)
     write(text1, '(F0.2)')(eltmjd - INT(eltmjd))*24.0
     result = DLGSET(Dlg, elapsed_hour, text1)
     write(text1, '(I0)')INT(dlts1)
     result = DLGSET(Dlg, timestep, text1)
     write(text1, '("(",I0,",",I0,")")')kloc, iloc
     result = DLGSET(Dlg, timestep_location, text1)
     write(text1, '(I0)')INT(mindlt)
     result = DLGSET(Dlg, min_timestep, text1)
     write(text1, '("(",I0,",",I0,")")')kmin, imin
     result = DLGSET(Dlg, min_timestep_location, text1)
     write(text1, '(I0)')INT(jdmin)
     result = DLGSET(Dlg, min_timestep_day, text1)
     write(text1, '(F0.2)')(jdmin - INT(jdmin))*24.0
     result = DLGSET(Dlg, min_timestep_hour, text1)
     write(text1, '(I0)')INT(dltav)
     result = DLGSET(Dlg, average_timestep, text1)
     write(text1, '(I0)')nit
     result = DLGSET(Dlg, iterations, text1)
     write(text1, '(I0)')nv
     result = DLGSET(Dlg, timestep_violations, text1)
     write(text1, '(F0.2)')FLOAT(nv)/FLOAT(nit)*100.0
     result = DLGSET(Dlg, percent, text1)
     write(text1, '(F0.2)')TAIR(jw)
     result = DLGSET(Dlg, air_temperature, text1)
     write(text1, '(F0.2)')TDEW(jw)
     result = DLGSET(Dlg, dew_point_temperature, text1)
     write(text1, '(F0.2)')WIND(jw)
     result = DLGSET(Dlg, wind_speed, text1)
     write(text1, '(F0.2)')PHI(jw)
     result = DLGSET(Dlg, wind_direction, text1)
     write(text1, '(F0.2)')CLOUD(jw)
     result = DLGSET(Dlg, cloud_cover, text1)
     write(text1, '(F0.1)')ET(DS(JBDN(jw)))
     result = DLGSET(Dlg, equilibrium_temp, text1)
     write(text1, '(F0.1)')CSHE(DS(JBDN(jw)))*rhowcp
     result = DLGSET(Dlg, surface_heat_exchange, text1)
     write(text1, '(F0.1)')SRON(jw)
     result = DLGSET(Dlg, solar_radiation, text1)
     write(text1, '(A)')time
     result = DLGSET(Dlg, current_time, text1)
     write(text1, '(I0)')kt
     result = DLGSET(Dlg, surface_layer, text1)
     write(text1, '(F0.2)')ELKT(jw)
     result = DLGSET(Dlg, surface_elevation, text1)
     write(text1, '(F0.2)')ZMIN(jw)
     result = DLGSET(Dlg, min_deviation, text1)
     write(text1, '(I0)')IZMIN(jw)
     result = DLGSET(Dlg, min_deviation_segment, text1)
     write(text1, '(1000(F0.2,2X))')qin
     result = DLGSET(Dlg, branch_inflow, text1)
     write(text1, '(1000(F0.2,2X))')tin
     result = DLGSET(Dlg, branch_inflow_temperature, text1)
     write(text1, '(1000(F0.2,2X))')qdtr
     result = DLGSET(Dlg, dist_tributary_inflow, text1)
     write(text1, '(1000(F0.2,2X))')tdtr
     result = DLGSET(Dlg, dist_tributary_inflow_temperature, text1)
     write(text1, '(1000(F0.3,2X))')qtr
     result = DLGSET(Dlg, tributary_inflow, text1)
     write(text1, '(1000(F0.2,2X))')ttr
     result = DLGSET(Dlg, tributary_inflow_temperature, text1)
     write(text1, '(1000(F0.2,2X))')qsum
     result = DLGSET(Dlg, outflow, text1)
     write(text1, '(1000(F0.2,2X))')qwd
     result = DLGSET(Dlg, withdrawal, text1)
     write(text1, '(1000(F0.2,2X))')qpu
     result = DLGSET(Dlg, pumpflow, text1)
     write(text1, '(1000(F0.2,2X))')qsp
     result = DLGSET(Dlg, spillwayflow, text1)
     write(text1, '(1000(F0.2,2X))')qpi
     result = DLGSET(Dlg, pipeflow, text1)
     write(text1, '(1000(F0.2,2X))')qgt
     result = DLGSET(Dlg, gateflow, text1)
     write(text1, '(A180)')moddir
     result = DLGSET(Dlg, modeldirectory, text1)
 
     write(text1, '(F0.2)')(current)/60.0
 
     result = DLGSET(Dlg, cpu_times, text1)
     if(mindlt>=1.0)write(text1, '(I0)')INT(mindlt)
     if(mindlt<1.0)write(text1, '(F0.3)')mindlt
     result = DLGSET(Dlg, min_timestep, text1)
 
     iprog = INT(((jday - tmstrt)/(tmend - tmstrt))*100)
                                                  ! range is 0 to 100
     result = DLGSET(Dlg, progressbar, iprog, dlg_position)
     return
 
     entry BLANK_DIALOG(Dlg)
     result = DLGSET(Dlg, gregorian_day, ' ')
     result = DLGSET(Dlg, julian_day, ' ')
     result = DLGSET(Dlg, julian_hour, ' ')
     result = DLGSET(Dlg, elapsed_day, ' ')
     result = DLGSET(Dlg, elapsed_hour, ' ')
     result = DLGSET(Dlg, timestep, ' ')
     result = DLGSET(Dlg, timestep_location, ' ')
     result = DLGSET(Dlg, min_timestep, ' ')
     result = DLGSET(Dlg, min_timestep_location, ' ')
     result = DLGSET(Dlg, min_timestep_day, ' ')
     result = DLGSET(Dlg, min_timestep_hour, ' ')
     result = DLGSET(Dlg, average_timestep, ' ')
     result = DLGSET(Dlg, iterations, ' ')
     result = DLGSET(Dlg, timestep_violations, ' ')
     result = DLGSET(Dlg, percent, ' ')
     result = DLGSET(Dlg, air_temperature, ' ')
     result = DLGSET(Dlg, dew_point_temperature, ' ')
     result = DLGSET(Dlg, wind_speed, ' ')
     result = DLGSET(Dlg, wind_direction, ' ')
     result = DLGSET(Dlg, cloud_cover, ' ')
     result = DLGSET(Dlg, equilibrium_temp, ' ')
     result = DLGSET(Dlg, surface_heat_exchange, ' ')
     result = DLGSET(Dlg, solar_radiation, ' ')
     result = DLGSET(Dlg, current_time, ' ')
     result = DLGSET(Dlg, surface_layer, ' ')
     result = DLGSET(Dlg, surface_elevation, ' ')
     result = DLGSET(Dlg, min_deviation, ' ')
     result = DLGSET(Dlg, min_deviation_segment, ' ')
     result = DLGSET(Dlg, branch_inflow, ' ')
     result = DLGSET(Dlg, branch_inflow_temperature, ' ')
     result = DLGSET(Dlg, dist_tributary_inflow, ' ')
     result = DLGSET(Dlg, dist_tributary_inflow_temperature, ' ')
     result = DLGSET(Dlg, tributary_inflow, ' ')
     result = DLGSET(Dlg, tributary_inflow_temperature, ' ')
     result = DLGSET(Dlg, outflow, ' ')
     result = DLGSET(Dlg, ending_time, ' ')
     result = DLGSET(Dlg, withdrawal, ' ')
     result = DLGSET(Dlg, pipeflow, ' ')
     result = DLGSET(Dlg, spillwayflow, ' ')
     result = DLGSET(Dlg, gateflow, ' ')
     result = DLGSET(Dlg, pumpflow, ' ')
     result = DLGSET(Dlg, modeldirectory, ' ')
     result = DLGSET(Dlg, progressbar, 0, dlg_position)
     end subroutine SCREEN_UPDATE
