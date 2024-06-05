
# getters and setters

cdef extern void climaradtranwrk_fup_a_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void climaradtranwrk_fup_a_get(void *ptr, int *dim1, int *dim2, double *arr)

cdef extern void climaradtranwrk_fdn_a_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void climaradtranwrk_fdn_a_get(void *ptr, int *dim1, int *dim2, double *arr)

cdef extern void climaradtranwrk_fup_n_get_size(void *ptr, int *dim1)
cdef extern void climaradtranwrk_fup_n_get(void *ptr, int *dim1, double *arr)

cdef extern void climaradtranwrk_fdn_n_get_size(void *ptr, int *dim1)
cdef extern void climaradtranwrk_fdn_n_get(void *ptr, int *dim1, double *arr)

cdef extern void climaradtranwrk_tau_band_get_size(void *ptr, int *dim1, int *dim2)
cdef extern void climaradtranwrk_tau_band_get(void *ptr, int *dim1, int *dim2, double *arr)