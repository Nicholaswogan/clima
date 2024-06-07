cimport Radtran_pxd as rad_pxd

cdef class Radtran:
  cdef rad_pxd.Radtran *_ptr

  def __init__(self):
    pass
  
  def __dealloc__(self):
    pass

  def skin_temperature(self, double bond_albedo):
    """The skin temperature

    Parameters
    ----------
    bond_albedo : float
        The bond albedo of a planet

    Returns 
    -------
    float
        The skin temperature
    """
    cdef double T_skin
    rad_pxd.radtran_skin_temperature_wrapper(self._ptr, &bond_albedo, &T_skin)
    return T_skin

  def opacities2yaml(self):
    """Returns a yaml string representing all opacities in the model.

    Returns 
    -------
    str
        string representing all opacities.
    """
    cdef int out_len
    cdef void *out_cp
    rad_pxd.radtran_opacities2yaml_wrapper_1(self._ptr, &out_len, &out_cp)
    cdef ndarray out_c = np.empty(out_len + 1, 'S1')
    rad_pxd.radtran_opacities2yaml_wrapper_2(self._ptr, &out_cp, &out_len, <char *>out_c.data)
    return out_c[:-1].tobytes().decode()
    
  property surface_albedo:
    "ndarray[double,ndim=1]. The surface albedo in each solar wavelength bin"
    def __get__(self):
      cdef int dim1
      rad_pxd.radtran_surface_albedo_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rad_pxd.radtran_surface_albedo_get(self._ptr, &dim1, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=1] arr):
      cdef int dim1
      rad_pxd.radtran_surface_albedo_get_size(self._ptr, &dim1)
      if arr.shape[0] != dim1:
        raise ClimaException('"surface_albedo" is the wrong size')
      rad_pxd.radtran_surface_albedo_set(self._ptr, &dim1, <double *>arr.data)

  property surface_emissivity:
    "ndarray[double,ndim=1]. The surface emissivity in each IR wavelength bin"
    def __get__(self):
      cdef int dim1
      rad_pxd.radtran_surface_emissivity_get_size(self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rad_pxd.radtran_surface_emissivity_get(self._ptr, &dim1, <double *>arr.data)
      return arr
    def __set__(self, ndarray[double, ndim=1] arr):
      cdef int dim1
      rad_pxd.radtran_surface_emissivity_get_size(self._ptr, &dim1)
      if arr.shape[0] != dim1:
        raise ClimaException('"surface_emissivity" is the wrong size')
      rad_pxd.radtran_surface_emissivity_set(self._ptr, &dim1, <double *>arr.data)

  property ir:
    "The OpticalProperties for longwave radiative transfer"
    def __get__(self):
      var = OpticalProperties()
      rad_pxd.radtran_ir_get(self._ptr, &var._ptr)
      return var

  property sol:
    "The OpticalProperties for shortwave radiative transfer"
    def __get__(self):
      var = OpticalProperties()
      rad_pxd.radtran_sol_get(self._ptr, &var._ptr)
      return var

  property wrk_ir:
    "The ClimaRadtranWrk for longwave radiative transfer"
    def __get__(self):
      var = ClimaRadtranWrk()
      rad_pxd.radtran_wrk_ir_get(self._ptr, &var._ptr)
      return var

  property wrk_sol:
    "The ClimaRadtranWrk for shortwave radiative transfer"
    def __get__(self):
      var = ClimaRadtranWrk()
      rad_pxd.radtran_wrk_sol_get(self._ptr, &var._ptr)
      return var
  


  