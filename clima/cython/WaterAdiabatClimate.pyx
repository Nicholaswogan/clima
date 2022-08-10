
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

  def make_profile(self, double T_surf, ndarray[double, ndim=1] P_i_surf):
    cdef int ng = P_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    wa_pxd.wateradiabatclimate_make_profile_wrapper(&self._ptr, &T_surf,
    &ng, <double *>P_i_surf.data, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    
  def make_column(self, double T_surf, ndarray[double, ndim=1] N_i_surf):
    cdef int ng = N_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    wa_pxd.wateradiabatclimate_make_column_wrapper(&self._ptr, &T_surf,
    &ng, <double *>N_i_surf.data, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def TOA_fluxes(self, double T_surf, ndarray[double, ndim=1] P_i_surf):
    cdef int ng = P_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double ISR, OLR;
    wa_pxd.wateradiabatclimate_toa_fluxes_wrapper(&self._ptr, &T_surf,
    &ng, <double *>P_i_surf.data, &ISR, &OLR, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return ISR, OLR

  def TOA_fluxes_column(self, double T_surf, ndarray[double, ndim=1] N_i_surf):
    cdef int ng = N_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double ISR, OLR;
    wa_pxd.wateradiabatclimate_toa_fluxes_column_wrapper(&self._ptr, &T_surf,
    &ng, <double *>N_i_surf.data, &ISR, &OLR, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return ISR, OLR
  
  def surface_temperature(self, ndarray[double, ndim=1] P_i_surf, double T_guess = 280):
    cdef int ng = P_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double T_surf;
    wa_pxd.wateradiabatclimate_surface_temperature_wrapper(&self._ptr, 
    &ng, <double *>P_i_surf.data, &T_guess, &T_surf, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return T_surf

  def surface_temperature_column(self, ndarray[double, ndim=1] N_i_surf, double T_guess = 280):
    cdef int ng = N_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double T_surf;
    wa_pxd.wateradiabatclimate_surface_temperature_column_wrapper(&self._ptr, 
    &ng, <double *>N_i_surf.data, &T_guess, &T_surf, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return T_surf

  property P_top:
    def __get__(self):
      cdef double val
      wa_pxd.wateradiabatclimate_p_top_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.wateradiabatclimate_p_top_set(&self._ptr, &val)

  property T_trop:
    def __get__(self):
      cdef double val
      wa_pxd.wateradiabatclimate_t_trop_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.wateradiabatclimate_t_trop_set(&self._ptr, &val)

  property RH:
    def __get__(self):
      cdef double val
      wa_pxd.wateradiabatclimate_rh_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.wateradiabatclimate_rh_set(&self._ptr, &val)

  property P:
    def __get__(self):
      cdef int dim1
      wa_pxd.wateradiabatclimate_p_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.wateradiabatclimate_p_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property T:
    def __get__(self):
      cdef int dim1
      wa_pxd.wateradiabatclimate_t_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.wateradiabatclimate_t_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property f_i:
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.wateradiabatclimate_f_i_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wa_pxd.wateradiabatclimate_f_i_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property z:
    def __get__(self):
      cdef int dim1
      wa_pxd.wateradiabatclimate_z_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.wateradiabatclimate_z_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property dz:
    def __get__(self):
      cdef int dim1
      wa_pxd.wateradiabatclimate_dz_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.wateradiabatclimate_dz_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property densities:
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.wateradiabatclimate_densities_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wa_pxd.wateradiabatclimate_densities_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr





 