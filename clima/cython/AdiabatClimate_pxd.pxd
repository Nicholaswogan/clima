from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

# callback signatures
ctypedef double (*temp_dependent_albedo_fcn)(double T_surf)
ctypedef void (*ocean_solubility_fcn)(double T_surf, int ng, double *P_i, double *m_i, void *args_p)

# allocate and destroy
cdef extern void *allocate_adiabatclimate();
cdef extern void deallocate_adiabatclimate(void *ptr);

# wrappers for functions
cdef extern void adiabatclimate_create_wrapper(void *ptr, char *species_file,
                                                    char *settings_file, char *flux_file, char *data_dir, bool *double_radiative_grid, char *err);

cdef extern void adiabatclimate_make_profile_wrapper(void *ptr, double *T_surf, int *ng, 
                                      double *P_i_surf, char *err)

cdef extern void adiabatclimate_make_column_wrapper(void *ptr, double *T_surf, int *ng, 
                                      double *N_i_surf, char *err)

cdef extern void adiabatclimate_make_profile_bg_gas_wrapper(void *ptr, double *T_surf, int *ng, 
                                      double *P_i_surf, double *P_surf, char *bg_gas, char *err)

cdef extern void adiabatclimate_toa_fluxes_wrapper(void *ptr, double *T_surf, int *ng, 
                                     double *P_i_surf, double *ISR, double *OLR, char *err)

cdef extern void adiabatclimate_toa_fluxes_column_wrapper(void *ptr, double *T_surf, int *ng, 
                                     double *P_i_surf, double *ISR, double *OLR, char *err)

cdef extern void adiabatclimate_toa_fluxes_bg_gas_wrapper(void *ptr, double *T_surf, int *ng, 
                                     double *P_i_surf, double *P_surf, char *bg_gas, double *ISR, 
                                     double *OLR, char *err)

cdef extern void adiabatclimate_surface_temperature_wrapper(void *ptr, int *ng, 
                                      double *P_i_surf, double *T_guess, double *T_surf, char *err)

cdef extern void adiabatclimate_surface_temperature_column_wrapper(void *ptr, int *ng, 
                                      double *N_i_surf, double *T_guess, double *T_surf, char *err)

cdef extern void adiabatclimate_surface_temperature_bg_gas_wrapper(void *ptr, int *ng, 
                                      double *P_i_surf, double *P_surf, char *bg_gas, 
                                      double *T_guess, double *T_surf, char *err)

cdef extern void adiabatclimate_rce_wrapper(void *ptr, int *ng, double *P_i_surf, double *T_surf_guess,
                                      int *dim_T_guess, double *T_guess, bool *convecting_with_below_present, 
                                      int *dim_convecting_with_below, bool *convecting_with_below, bool *converged, char *err)

cdef extern void adiabatclimate_set_ocean_solubility_fcn_wrapper(void *ptr, char *species_c, ocean_solubility_fcn fcn, char *err)

cdef extern void adiabatclimate_to_regular_grid_wrapper(void *ptr, char *err)

cdef extern void adiabatclimate_out2atmosphere_txt_wrapper(void *ptr, char *filename, int *nz, double *eddy, bool *overwrite, bool *clip, char *err)

cdef extern void adiabatclimate_heat_redistribution_parameters_wrapper(void *ptr, double *tau_LW, double *k_term, double *f_term, char *err)

# getters and setters
cdef extern void adiabatclimate_p_top_get(void *ptr, double *val)
cdef extern void adiabatclimate_p_top_set(void *ptr, double *val)

cdef extern void adiabatclimate_t_trop_get(void *ptr, double *val)
cdef extern void adiabatclimate_t_trop_set(void *ptr, double *val)

cdef extern void adiabatclimate_use_make_column_p_guess_get(void *ptr, bool *val)
cdef extern void adiabatclimate_use_make_column_p_guess_set(void *ptr, bool *val)

cdef extern void adiabatclimate_make_column_p_guess_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_make_column_p_guess_get(void *ptr, int *dim1, double *val)
cdef extern void adiabatclimate_make_column_p_guess_set(void *ptr, int *dim1, double *val)

cdef extern void adiabatclimate_solve_for_t_trop_get(void *ptr, bool *val)
cdef extern void adiabatclimate_solve_for_t_trop_set(void *ptr, bool *val)

cdef extern void adiabatclimate_albedo_fcn_set(void *ptr, temp_dependent_albedo_fcn fcn)

cdef extern void adiabatclimate_rh_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_rh_get(void *ptr, int *dim1, double *val)
cdef extern void adiabatclimate_rh_set(void *ptr, int *dim1, double *val)

cdef extern void adiabatclimate_ocean_args_p_set(void *ptr, void *p)

cdef extern void adiabatclimate_tidally_locked_dayside_get(void *ptr, bool *val)
cdef extern void adiabatclimate_tidally_locked_dayside_set(void *ptr, bool *val)

cdef extern void adiabatclimate_l_get(void *ptr, double *val)
cdef extern void adiabatclimate_l_set(void *ptr, double *val)

cdef extern void adiabatclimate_chi_get(void *ptr, double *val)
cdef extern void adiabatclimate_chi_set(void *ptr, double *val)

cdef extern void adiabatclimate_n_lw_get(void *ptr, double *val)
cdef extern void adiabatclimate_n_lw_set(void *ptr, double *val)

cdef extern void adiabatclimate_cd_get(void *ptr, double *val)
cdef extern void adiabatclimate_cd_set(void *ptr, double *val)

cdef extern void adiabatclimate_surface_heat_flow_get(void *ptr, double *val)
cdef extern void adiabatclimate_surface_heat_flow_set(void *ptr, double *val)

cdef extern void adiabatclimate_species_names_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_species_names_get(void *ptr, int *dim1, char* species_names)

cdef extern void adiabatclimate_rad_get(void *ptr, void *ptr1)

cdef extern void adiabatclimate_convecting_with_below_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_convecting_with_below_get(void *ptr, int *dim1, bool *arr)

cdef extern void adiabatclimate_lapse_rate_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_lapse_rate_get(void *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_lapse_rate_intended_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_lapse_rate_intended_get(void *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_convective_newton_step_size_get(void *ptr, double *val)
cdef extern void adiabatclimate_convective_newton_step_size_set(void *ptr, double *val)

cdef extern void adiabatclimate_rtol_get(void *ptr, double *val)
cdef extern void adiabatclimate_rtol_set(void *ptr, double *val)

cdef extern void adiabatclimate_atol_get(void *ptr, double *val)
cdef extern void adiabatclimate_atol_set(void *ptr, double *val)

cdef extern void adiabatclimate_tol_make_column_get(void *ptr, double *val)
cdef extern void adiabatclimate_tol_make_column_set(void *ptr, double *val)

cdef extern void adiabatclimate_epsj_get(void *ptr, double *val)
cdef extern void adiabatclimate_epsj_set(void *ptr, double *val)

cdef extern void adiabatclimate_xtol_rc_get(void *ptr, double *val)
cdef extern void adiabatclimate_xtol_rc_set(void *ptr, double *val)

cdef extern void adiabatclimate_max_rc_iters_get(void *ptr, int *val)
cdef extern void adiabatclimate_max_rc_iters_set(void *ptr, int *val)

cdef extern void adiabatclimate_max_rc_iters_convection_get(void *ptr, int *val)
cdef extern void adiabatclimate_max_rc_iters_convection_set(void *ptr, int *val)

cdef extern void adiabatclimate_radiation_norm_term_get(void *ptr, double *val)
cdef extern void adiabatclimate_radiation_norm_term_set(void *ptr, double *val)

cdef extern void adiabatclimate_verbose_get(void *ptr, bool *val)
cdef extern void adiabatclimate_verbose_set(void *ptr, bool *val)

cdef extern void adiabatclimate_p_surf_get(void *ptr, double *val)

cdef extern void adiabatclimate_p_trop_get(void *ptr, double *val)

cdef extern void adiabatclimate_p_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_p_get(void *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_t_surf_get(void *ptr, double *val)

cdef extern void adiabatclimate_t_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_t_get(void *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_f_i_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void adiabatclimate_f_i_get(void *ptr, int *dim1, int *dim2, double *arr)

cdef extern void adiabatclimate_z_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_z_get(void *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_dz_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_dz_get(void *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_densities_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void adiabatclimate_densities_get(void *ptr, int *dim1, int *dim2, double *arr)

cdef extern void adiabatclimate_n_atmos_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_n_atmos_get(void *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_n_surface_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_n_surface_get(void *ptr, int *dim1, double *arr)

cdef extern void adiabatclimate_n_ocean_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void adiabatclimate_n_ocean_get(void *ptr, int *dim1, int *dim2, double *arr)



