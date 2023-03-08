
cimport AdiabatClimate_pxd as wa_pxd
  
cdef class AdiabatClimate:
  cdef void *_ptr

  def __init__(self, species_file = None, settings_file = None, 
                     flux_file = None):           
    # Allocate memory
    wa_pxd.allocate_adiabatclimate(&self._ptr)
    
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
    wa_pxd.adiabatclimate_create_wrapper(&self._ptr, data_dir_c, species_file_c,
                                              settings_file_c, flux_file_c,
                                              err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def __dealloc__(self):
    wa_pxd.deallocate_adiabatclimate(&self._ptr);

  def make_profile(self, double T_surf, ndarray[double, ndim=1] P_i_surf):
    cdef int ng = P_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_make_profile_wrapper(&self._ptr, &T_surf,
    &ng, <double *>P_i_surf.data, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    
  def make_column(self, double T_surf, ndarray[double, ndim=1] N_i_surf):
    cdef int ng = N_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_make_column_wrapper(&self._ptr, &T_surf,
    &ng, <double *>N_i_surf.data, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
  
  def make_profile_bg_gas(self, double T_surf, ndarray[double, ndim=1] P_i_surf, double P_surf, str bg_gas):
    cdef int ng = P_i_surf.shape[0]
    cdef bytes bg_gas_b = pystring2cstring(bg_gas)
    cdef char *bg_gas_c = bg_gas_b
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_make_profile_bg_gas_wrapper(&self._ptr, &T_surf,
    &ng, <double *>P_i_surf.data, &P_surf, bg_gas_c, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def TOA_fluxes(self, double T_surf, ndarray[double, ndim=1] P_i_surf):
    cdef int ng = P_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double ISR, OLR;
    wa_pxd.adiabatclimate_toa_fluxes_wrapper(&self._ptr, &T_surf,
    &ng, <double *>P_i_surf.data, &ISR, &OLR, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return ISR, OLR

  def TOA_fluxes_column(self, double T_surf, ndarray[double, ndim=1] N_i_surf):
    cdef int ng = N_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double ISR, OLR;
    wa_pxd.adiabatclimate_toa_fluxes_column_wrapper(&self._ptr, &T_surf,
    &ng, <double *>N_i_surf.data, &ISR, &OLR, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return ISR, OLR

  def TOA_fluxes_bg_gas(self, double T_surf, ndarray[double, ndim=1] P_i_surf, double P_surf, str bg_gas):
    cdef int ng = P_i_surf.shape[0]
    cdef bytes bg_gas_b = pystring2cstring(bg_gas)
    cdef char *bg_gas_c = bg_gas_b
    cdef char err[ERR_LEN+1]
    cdef double ISR, OLR;
    wa_pxd.adiabatclimate_toa_fluxes_bg_gas_wrapper(&self._ptr, &T_surf,
    &ng, <double *>P_i_surf.data, &P_surf, bg_gas_c, &ISR, &OLR, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return ISR, OLR
  
  def surface_temperature(self, ndarray[double, ndim=1] P_i_surf, double T_guess = 280):
    cdef int ng = P_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double T_surf;
    wa_pxd.adiabatclimate_surface_temperature_wrapper(&self._ptr, 
    &ng, <double *>P_i_surf.data, &T_guess, &T_surf, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return T_surf

  def surface_temperature_column(self, ndarray[double, ndim=1] N_i_surf, double T_guess = 280):
    cdef int ng = N_i_surf.shape[0]
    cdef char err[ERR_LEN+1]
    cdef double T_surf;
    wa_pxd.adiabatclimate_surface_temperature_column_wrapper(&self._ptr, 
    &ng, <double *>N_i_surf.data, &T_guess, &T_surf, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return T_surf

  def surface_temperature_bg_gas(self, ndarray[double, ndim=1] P_i_surf, double P_surf, str bg_gas, double T_guess = 280):
    cdef int ng = P_i_surf.shape[0]
    cdef bytes bg_gas_b = pystring2cstring(bg_gas)
    cdef char *bg_gas_c = bg_gas_b
    cdef char err[ERR_LEN+1]
    cdef double T_surf;
    wa_pxd.adiabatclimate_surface_temperature_bg_gas_wrapper(&self._ptr, 
    &ng, <double *>P_i_surf.data, &P_surf, bg_gas_c, &T_guess, &T_surf, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())
    return T_surf

  def to_regular_grid(self):
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_to_regular_grid_wrapper(&self._ptr, err)
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def out2atmosphere_txt(self, filename, ndarray[double, ndim=1] eddy, bool overwrite = False, bool clip = True):
    cdef int nz = eddy.shape[0]
    cdef bytes filename_b = pystring2cstring(filename)
    cdef char *filename_c = filename_b
    cdef char err[ERR_LEN+1]
    wa_pxd.adiabatclimate_out2atmosphere_txt_wrapper(&self._ptr, filename_c, &nz, <double *>eddy.data, &overwrite, &clip, err)  
    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  property P_top:
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_p_top_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_p_top_set(&self._ptr, &val)

  property T_trop:
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_t_trop_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      wa_pxd.adiabatclimate_t_trop_set(&self._ptr, &val)

  property solve_for_T_trop:
    def __get__(self):
      cdef bool val
      wa_pxd.adiabatclimate_solve_for_t_trop_get(&self._ptr, &val)
      return val
    def __set__(self, bool val):
      wa_pxd.adiabatclimate_solve_for_t_trop_set(&self._ptr, &val)

  property albedo_fcn:
    def __set__(self, object fcn):
      cdef bool set_to_null
      cdef unsigned long long int fcn_l
      cdef wa_pxd.temp_dependent_albedo_fcn fcn_c
      if isinstance(fcn, type(None)):
        set_to_null = True
      else:
        set_to_null = False
        argtypes = (ct.c_double,)
        restype = ct.c_double
        if not fcn.ctypes.argtypes == argtypes:
          raise ClimaException("The callback function has the wrong argument types.")
        if not fcn.ctypes.restype == restype:
          raise ClimaException("The callback function has the wrong return type.")
        fcn_l = <unsigned long long int> fcn.address
        fcn_c = <wa_pxd.temp_dependent_albedo_fcn> fcn_l
      wa_pxd.adiabatclimate_albedo_fcn_set(&self._ptr, &set_to_null, fcn_c)

  property RH:
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_rh_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_rh_get(&self._ptr, &dim1, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=1] arr):
      cdef int dim1
      wa_pxd.adiabatclimate_rh_get_size(&self._ptr, &dim1)
      if arr.shape[0] != dim1:
        raise ClimaException('"RH" is the wrong size')
      wa_pxd.adiabatclimate_rh_set(&self._ptr, &dim1, <double *>arr.data)

  property species_names:
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_species_names_get_size(&self._ptr, &dim1)
      cdef ndarray species_names_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      wa_pxd.adiabatclimate_species_names_get(&self._ptr, &dim1, <char *>species_names_c.data)
      return c2stringarr(species_names_c, S_STR_LEN, dim1)

  property rad:
    def __get__(self):
      cdef void *ptr1
      wa_pxd.adiabatclimate_rad_get(&self._ptr, &ptr1)
      var = Radtran()
      var._ptr = ptr1
      return var

  property P_surf:
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_p_surf_get(&self._ptr, &val)
      return val

  property P_trop:
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_p_trop_get(&self._ptr, &val)
      return val

  property P:
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_p_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_p_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property T_surf:
    def __get__(self):
      cdef double val
      wa_pxd.adiabatclimate_t_surf_get(&self._ptr, &val)
      return val

  property T:
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_t_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_t_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property f_i:
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.adiabatclimate_f_i_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wa_pxd.adiabatclimate_f_i_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property z:
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_z_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_z_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property dz:
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_dz_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_dz_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property densities:
    def __get__(self):
      cdef int dim1, dim2
      wa_pxd.adiabatclimate_densities_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      wa_pxd.adiabatclimate_densities_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property N_surface:
    def __get__(self):
      cdef int dim1
      wa_pxd.adiabatclimate_n_surface_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      wa_pxd.adiabatclimate_n_surface_get(&self._ptr, &dim1, <double *>arr.data)
      return arr





 