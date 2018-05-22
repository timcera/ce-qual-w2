!*==envirpmod.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module ENVIRPMOD
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     character(3), allocatable, dimension(:), save :: cc_e, cd_e
     real, allocatable, dimension(:), save :: cdn_e, cd_avg, cd_int, cd_sum,   &
       & cd_top, cn_e, c_avg, c_int, c_sum, c_top, d_cnt, d_tot, sumvolt,      &
       & t_cnt, t_tot, vel_c, volgl, v_cnt, v_tot
     real, allocatable, dimension(:, :, :), save :: cd_class, c_class
     real, allocatable, dimension(:, :), save :: cd_cnt, cd_tot, conc_c,       &
       & conc_cd, c_cnt, c_tot, d_class, t_class, v_class
     real, save :: cd_crit, c_crit, dltt, d_avg, d_c, d_crit, d_int, d_sum,    &
                 & d_top, sjday1, sjday2, temp_c, temp_int, temp_top, timlast, &
                 & t_avg, t_crit, t_sum, vel_int, vel_top, v_avg, v_crit, v_sum
     integer, save :: cone
     character(3), save :: depth_vpr, selectivec, temp_vpr, vel_vpr
     integer, dimension(9) :: iend, istart
     integer :: iopenfish, i_segint, jacd, jj, nacd_e, nac_e, numclass
!
!*** End of declarations rewritten by SPAG
!
     data cone/1500/
                 !SW 5/26/15
     end module ENVIRPMOD
