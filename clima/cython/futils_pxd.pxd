cdef extern void futils_rebin_wrapper(
  int *old_bins_len, double *old_bins, 
  int *old_vals_len, double *old_vals, 
  int *new_bins_len, double *new_bins, 
  int *new_vals_len, double *new_vals, 
  int *ierr
)

cdef extern void futils_rebin_with_errors_wrapper(
  int *old_bins_len, double *old_bins, 
  int *old_vals_len, double *old_vals, 
  int *old_errs_len, double *old_errs, 
  int *new_bins_len, double *new_bins, 
  int *new_vals_len, double *new_vals, 
  int *new_errs_len, double *new_errs, 
  int *ierr
)

cdef extern void futils_grid_at_resolution_wrapper1(
  double *wv_min, double *wv_max, double *R, int *wavl_len,
  void *wavl_ptr, int *ierr 
)
cdef extern void futils_grid_at_resolution_wrapper2(int *wavl_len, void *wavl_ptr, double *wavl)

cdef extern void futils_make_bins_wrapper(
  int *wv_len, double *wv, int *wavl_len, double *wavl, int *ierr
)

cdef extern void futils_rebin_error_message_wrapper1(
  int *ierr, int *err_length, void *err_ptr
)
cdef extern void futils_rebin_error_message_wrapper2(
  int *err_length, void *err_ptr, char *err
)