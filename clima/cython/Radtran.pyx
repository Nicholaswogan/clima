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

  def set_custom_optical_properties(self, ndarray[double, ndim=1] wv, ndarray[double, ndim=1] P, ndarray[double, ndim=2] dtau_dz, ndarray[double, ndim=2] w0, ndarray[double, ndim=2] g0):
    """Sets particle densities and radii.

    Parameters
    ----------
    P_i_surf : ndarray[double,ndim=1]
        Array of pressures in dynes/cm^2
    pdensities : ndarray[double,ndim=2]
        Particle densities in particles/cm^3 at each pressure and for each 
        particle in the model. Shape (nz, np).
    pradii : ndarray[double,ndim=2]
        Particle radii in cm at each pressure and for each particle
        in the model. Shape (nz, np).
    """
    cdef int dim_wv = wv.shape[0]
    cdef int dim_P = P.shape[0]
    cdef int dim1_dtau_dz = dtau_dz.shape[0]
    cdef int dim2_dtau_dz = dtau_dz.shape[1]
    dtau_dz = np.asfortranarray(dtau_dz)
    cdef int dim1_w0 = w0.shape[0]
    cdef int dim2_w0 = w0.shape[1]
    w0 = np.asfortranarray(w0)
    cdef int dim1_g0 = g0.shape[0]
    cdef int dim2_g0 = g0.shape[1]
    g0 = np.asfortranarray(g0)

    cdef char err[ERR_LEN+1]

    rad_pxd.radtran_set_custom_optical_properties(
      self._ptr, &dim_wv, <double *> wv.data, &dim_P, <double *> P.data,  
      &dim1_dtau_dz, &dim2_dtau_dz, <double *> dtau_dz.data, 
      &dim1_w0, &dim2_w0, <double *> w0.data, 
      &dim1_g0, &dim2_g0, <double *> g0.data, 
      err
    )

    if len(err.strip()) > 0:
      raise ClimaException(err.decode("utf-8").strip())

  def unset_custom_optical_properties(self):
    rad_pxd.radtran_unset_custom_optical_properties(self._ptr)
    
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
  


  