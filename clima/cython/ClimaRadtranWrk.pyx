cimport ClimaRadtranWrk_pxd as rwrk_pxd

cdef class ClimaRadtranWrk:
  cdef void *_ptr

  def __init__(self):
    pass
  
  def __dealloc__(self):
    pass
  
  property fup_a:
    def __get__(self):
      cdef int dim1, dim2
      rwrk_pxd.climaradtranwrk_fup_a_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      rwrk_pxd.climaradtranwrk_fup_a_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property fdn_a:
    def __get__(self):
      cdef int dim1, dim2
      rwrk_pxd.climaradtranwrk_fdn_a_get_size(&self._ptr, &dim1, &dim2)
      cdef ndarray arr = np.empty((dim1, dim2), np.double, order="F")
      rwrk_pxd.climaradtranwrk_fdn_a_get(&self._ptr, &dim1, &dim2, <double *>arr.data)
      return arr

  property fup_n:
    def __get__(self):
      cdef int dim1
      rwrk_pxd.climaradtranwrk_fup_n_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rwrk_pxd.climaradtranwrk_fup_n_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property fdn_n:
    def __get__(self):
      cdef int dim1
      rwrk_pxd.climaradtranwrk_fdn_n_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      rwrk_pxd.climaradtranwrk_fdn_n_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

      