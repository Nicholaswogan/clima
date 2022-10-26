cimport Radtran_pxd as rad_pxd

cdef class Radtran:
  cdef void *_ptr

  def __init__(self):
    pass
  
  def __dealloc__(self):
    pass

  property ir:
    def __get__(self):
      cdef void *ptr1;
      rad_pxd.radtran_ir_get(&self._ptr, &ptr1)
      var = OpticalProperties()
      var._ptr = ptr1
      return var

  property sol:
    def __get__(self):
      cdef void *ptr1;
      rad_pxd.radtran_sol_get(&self._ptr, &ptr1)
      var = OpticalProperties()
      var._ptr = ptr1
      return var

  property wrk_ir:
    def __get__(self):
      cdef void *ptr1;
      rad_pxd.radtran_wrk_ir_get(&self._ptr, &ptr1)
      var = ClimaRadtranWrk()
      var._ptr = ptr1
      return var

  property wrk_sol:
    def __get__(self):
      cdef void *ptr1;
      rad_pxd.radtran_wrk_sol_get(&self._ptr, &ptr1)
      var = ClimaRadtranWrk()
      var._ptr = ptr1
      return var
  


  