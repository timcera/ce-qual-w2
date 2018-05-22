!*==logicc.spg  processed by SPAG 6.70Rc at 14:33 on 22 May 2018
     module LOGICC
     implicit none
!
!*** Start of declarations rewritten by SPAG
!
! Local variables
!
     logical, allocatable, dimension(:) :: celerity_limit, constituent_plot,   &
          & dam_inflow, dam_outflow, derived_plot, dh_external, dh_internal,   &
          & dist_tribs, dn_flow, dq_external, dq_internal, fresh_water,        &
          & hydro_plot, implicit_az, internal_flow, interp_dtribs,             &
          & interp_extinction, interp_gate, interp_head, interp_inflow,        &
          & interp_meteorology, interp_tribs, interp_withdrawal, limiting_dlt, &
          & limiting_factor, mannings_n, no_heat, no_inflow, no_outflow,       &
          & no_wind, one_layer, ph_calc, precipitation, print_sediment,        &
          & print_sediment1, print_sediment2, read_extinction, read_radiation, &
          & salt_water, term_by_term, trapezoidal, uh_external, uh_internal,   &
          & ultimate, upwind, up_flow, uq_external, uq_internal,               &
          & viscosity_limit
     logical :: gates, initialize_graph, oxygen_demand, pipes, susp_solids,    &
              & tributaries, update_graph, withdrawals
     logical, allocatable, dimension(:, :) :: internal_weir, interp_outflow,   &
          & point_sink, print_const, print_derived, print_epiphyton,           &
          & print_hydro
!
!*** End of declarations rewritten by SPAG
!
                                                                                      ! Amaila
                                                                                                                       !TC 08/03/04
                                                                                                                        !SW 07/16/04
     end module LOGICC                                              ! cb 8/13/2010
