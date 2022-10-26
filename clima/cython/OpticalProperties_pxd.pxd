
# getters and setters

cdef extern void opticalproperties_wavl_get_size(void *ptr, int *dim1)
cdef extern void opticalproperties_wavl_get(void *ptr, int *dim1, double *arr)

cdef extern void opticalproperties_freq_get_size(void *ptr, int *dim1)
cdef extern void opticalproperties_freq_get(void *ptr, int *dim1, double *arr)