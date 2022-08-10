from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

# allocate and destroy
cdef extern void allocate_wateradiabatclimate(void *ptr);
cdef extern void deallocate_wateradiabatclimate(void *ptr);

cdef extern void wateradiabatclimate_create_wrapper(void *ptr, char *data_dir, char *species_file,
                                                    char *settings_file, char *flux_file, char *err);

cdef extern void wateradiabatclimate_make_profile_wrapper(void *ptr, double *T_surf, int *ng, 
                                      double *P_i_surf, char *err)

cdef extern void wateradiabatclimate_make_column_wrapper(void *ptr, double *T_surf, int *ng, 
                                      double *N_i_surf, char *err)

cdef extern void wateradiabatclimate_toa_fluxes_wrapper(void *ptr, double *T_surf, int *ng, 
                                     double *P_i_surf, double *ISR, double *OLR, char *err)

cdef extern void wateradiabatclimate_toa_fluxes_column_wrapper(void *ptr, double *T_surf, int *ng, 
                                     double *P_i_surf, double *ISR, double *OLR, char *err)

cdef extern void wateradiabatclimate_surface_temperature_wrapper(void *ptr, int *ng, 
                                      double *P_i_surf, double *T_guess, double *T_surf, char *err)

cdef extern void wateradiabatclimate_surface_temperature_column_wrapper(void *ptr, int *ng, 
                                      double *N_i_surf, double *T_guess, double *T_surf, char *err)


cdef extern void wateradiabatclimate_p_top_get(void *ptr, double *val)
cdef extern void wateradiabatclimate_p_top_set(void *ptr, double *val)

cdef extern void wateradiabatclimate_t_trop_get(void *ptr, double *val)
cdef extern void wateradiabatclimate_t_trop_set(void *ptr, double *val)

cdef extern void wateradiabatclimate_rh_get(void *ptr, double *val)
cdef extern void wateradiabatclimate_rh_set(void *ptr, double *val)

cdef extern void wateradiabatclimate_p_get_size(void *ptr, int *dim1)
cdef extern void wateradiabatclimate_p_get(void *ptr, int *dim1, double *arr)

cdef extern void wateradiabatclimate_t_get_size(void *ptr, int *dim1)
cdef extern void wateradiabatclimate_t_get(void *ptr, int *dim1, double *arr)

cdef extern void wateradiabatclimate_f_i_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void wateradiabatclimate_f_i_get(void *ptr, int *dim1, int *dim2, double *arr)

cdef extern void wateradiabatclimate_z_get_size(void *ptr, int *dim1)
cdef extern void wateradiabatclimate_z_get(void *ptr, int *dim1, double *arr)

cdef extern void wateradiabatclimate_dz_get_size(void *ptr, int *dim1)
cdef extern void wateradiabatclimate_dz_get(void *ptr, int *dim1, double *arr)

cdef extern void wateradiabatclimate_densities_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void wateradiabatclimate_densities_get(void *ptr, int *dim1, int *dim2, double *arr)


