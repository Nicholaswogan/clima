cimport Radtran_pxd as rad_pxd

cdef class Radtran:
  cdef void *_ptr

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
    rad_pxd.radtran_skin_temperature_wrapper(&self._ptr, &bond_albedo, &T_skin)
    return T_skin

  property surface_albedo:
    "The surface albedo"
    def __get__(self):
      cdef double val
      rad_pxd.radtran_surface_albedo_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      rad_pxd.radtran_surface_albedo_set(&self._ptr, &val)

  property ir:
    "The OpticalProperties for longwave radiative transfer"
    def __get__(self):
      cdef void *ptr1;
      rad_pxd.radtran_ir_get(&self._ptr, &ptr1)
      var = OpticalProperties()
      var._ptr = ptr1
      return var

  property sol:
    "The OpticalProperties for shortwave radiative transfer"
    def __get__(self):
      cdef void *ptr1;
      rad_pxd.radtran_sol_get(&self._ptr, &ptr1)
      var = OpticalProperties()
      var._ptr = ptr1
      return var

  property wrk_ir:
    "The ClimaRadtranWrk for longwave radiative transfer"
    def __get__(self):
      cdef void *ptr1;
      rad_pxd.radtran_wrk_ir_get(&self._ptr, &ptr1)
      var = ClimaRadtranWrk()
      var._ptr = ptr1
      return var

  property wrk_sol:
    "The ClimaRadtranWrk for shortwave radiative transfer"
    def __get__(self):
      cdef void *ptr1;
      rad_pxd.radtran_wrk_sol_get(&self._ptr, &ptr1)
      var = ClimaRadtranWrk()
      var._ptr = ptr1
      return var
  


  