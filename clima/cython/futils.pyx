cimport futils_pxd as f_pxd

cpdef rebin(ndarray[double, ndim=1] old_bins, ndarray[double, ndim=1]  old_vals, ndarray[double, ndim=1]  new_bins):
  cdef int n_old = old_bins.shape[0] - 1
  cdef int n_new = new_bins.shape[0] - 1
  cdef ndarray new_vals = np.empty(n_new, np.double)
  cdef int ierr

  if old_vals.shape[0] != n_old:
    raise ValueError('"old_bins" and "old_vals" arguments to "rebin" have incompatible shapes')

  f_pxd.futils_rebin_wrapper(&n_old, <double *> old_bins.data, <double *> old_vals.data, 
                             &n_new, <double *> new_bins.data, <double *> new_vals.data, &ierr)
  if ierr < 0:
    raise Exception("rebin returned error code: "+str(ierr))

  return new_vals

  
