
cimport WaterAdiabatClimate_pxd as wa_pxd
  
cdef class WaterAdiabatClimate:
  cdef void *_ptr

  def __init__(self, species_file = None, settings_file = None, 
                     flux_file = None):           
    # Allocate memory
    wa_pxd.allocate_wateradiabatclimate(&self._ptr)
    
    # convert strings to char
    cdef bytes data_dir_b = pystring2cstring(os.path.dirname(os.path.realpath(__file__))+'/data')
    cdef char *data_dir_c = data_dir_b
    cdef bytes species_file_b = pystring2cstring(species_file)
    cdef char *species_file_c = species_file_b
    cdef bytes settings_file_b = pystring2cstring(settings_file)
    cdef char *settings_file_c = settings_file_b
    cdef bytes flux_file_b = pystring2cstring(flux_file)
    cdef char *flux_file_c = flux_file_b
    cdef char err[ERR_LEN+1]
    
    # Initialize
    wa_pxd.wateradiabatclimate_create_wrapper(&self._ptr, data_dir_c, species_file_c,
                                              settings_file_c, flux_file_c,
                                              err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def __dealloc__(self):
    wa_pxd.deallocate_wateradiabatclimate(&self._ptr);

  def OLR(self, double T_surf, ndarray[double, ndim=1] P_i_surf):
    cdef int ng = P_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double OLR
    wa_pxd.wateradiabatclimate_olr_wrapper(&self._ptr, &T_surf,
    &ng, <double *>P_i_surf.data, &OLR, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return OLR
  
  def surface_temperature(self, ndarray[double, ndim=1] P_i_surf, double T_guess = 300):
    cdef int ng = P_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double T_surf;
    wa_pxd.wateradiabatclimate_surface_temperature_wrapper(&self._ptr, 
    &ng, <double *>P_i_surf.data, &T_guess, &T_surf, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return T_surf

  def surface_temperature_column(self, ndarray[double, ndim=1] N_i_surf, double T_guess = 300):
    cdef int ng = N_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double T_surf;
    wa_pxd.wateradiabatclimate_surface_temperature_column_wrapper(&self._ptr, 
    &ng, <double *>N_i_surf.data, &T_guess, &T_surf, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return T_surf





 