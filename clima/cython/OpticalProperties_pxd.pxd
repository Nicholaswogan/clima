
cdef extern from *:
  struct OpticalProperties:
    pass

# getters and setters

cdef extern void opticalproperties_wavl_get_size(OpticalProperties *ptr, int *dim1)
cdef extern void opticalproperties_wavl_get(OpticalProperties *ptr, int *dim1, double *arr)

cdef extern void opticalproperties_freq_get_size(OpticalProperties *ptr, int *dim1)
cdef extern void opticalproperties_freq_get(OpticalProperties *ptr, int *dim1, double *arr)