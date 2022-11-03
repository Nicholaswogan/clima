from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

# allocate and destroy
cdef extern void allocate_adiabatclimate(void *ptr);
cdef extern void deallocate_adiabatclimate(void *ptr);

# wrappers for functions
cdef extern void adiabatclimate_create_wrapper(void *ptr, char *data_dir, char *species_file,
                                                    char *settings_file, char *flux_file, char *err);

cdef extern void adiabatclimate_make_profile_wrapper(void *ptr, double *T_surf, int *ng, 
                                      double *P_i_surf, char *err)

cdef extern void adiabatclimate_make_column_wrapper(void *ptr, double *T_surf, int *ng, 
                                      double *N_i_surf, char *err)

cdef extern void adiabatclimate_toa_fluxes_wrapper(void *ptr, double *T_surf, int *ng, 
                                     double *P_i_surf, double *ISR, double *OLR, char *err)

cdef extern void adiabatclimate_toa_fluxes_column_wrapper(void *ptr, double *T_surf, int *ng, 
                                     double *P_i_surf, double *ISR, double *OLR, char *err)

cdef extern void adiabatclimate_surface_temperature_wrapper(void *ptr, int *ng, 
                                      double *P_i_surf, double *T_guess, double *T_surf, char *err)

cdef extern void adiabatclimate_surface_temperature_column_wrapper(void *ptr, int *ng, 
                                      double *N_i_surf, double *T_guess, double *T_surf, char *err)

cdef extern void adiabatclimate_to_regular_grid_wrapper(void *ptr, char *err)

cdef extern void adiabatclimate_out2atmosphere_txt_wrapper(void *ptr, char *filename, int *nz, double *eddy, bool *overwrite, bool *clip, char *err)

# getters and setters
cdef extern void adiabatclimate_p_top_get(void *ptr, double *val)
cdef extern void adiabatclimate_p_top_set(void *ptr, double *val)

cdef extern void adiabatclimate_t_trop_get(void *ptr, double *val)
cdef extern void adiabatclimate_t_trop_set(void *ptr, double *val)

cdef extern void adiabatclimate_rh_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_rh_get(void *ptr, int *dim1, double *val)
cdef extern void adiabatclimate_rh_set(void *ptr, int *dim1, double *val)

cdef extern void adiabatclimate_species_names_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_species_names_get(void *ptr, int *dim1, char* species_names)

cdef extern void adiabatclimate_p_get_size(void *ptr, int *dim1)
cdef extern void adiabatclimate_p_get(void *ptr, int *dim1, double *arr)

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

cdef extern void adiabatclimate_rad_get(void *ptr, void *ptr1)

