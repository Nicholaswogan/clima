from libcpp cimport bool
cdef extern from "<stdbool.h>":
  pass

# allocate and destroy
cdef extern void allocate_wateradiabatclimate(void *ptr);
cdef extern void deallocate_wateradiabatclimate(void *ptr);

cdef extern void wateradiabatclimate_create_wrapper(void *ptr, char *data_dir, char *species_file,
                                                    char *settings_file, char *flux_file, char *err);

cdef extern void wateradiabatclimate_olr_wrapper(void *ptr, double *T_surf, int *ng, 
                                      double *P_i_surf, double *OLR, char *err)

cdef extern void wateradiabatclimate_surface_temperature_wrapper(void *ptr, int *ng, 
                                      double *P_i_surf, double *T_guess, double *T_surf, char *err)

cdef extern void wateradiabatclimate_surface_temperature_column_wrapper(void *ptr, int *ng, 
                                      double *N_i_surf, double *T_guess, double *T_surf, char *err)

