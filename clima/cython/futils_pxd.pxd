cdef extern void futils_rebin_wrapper(int *n_old, double *old_bins, double *old_vals, 
                                      int *n_new, double *new_bins, double *new_vals, int *ierr)