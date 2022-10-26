! futils

subroutine futils_rebin_wrapper(n_old, old_bins, old_vals, n_new, new_bins, new_vals, ierr) bind(c)
  use futils, only: rebin
  integer(c_int), intent(in) :: n_old
  real(c_double), intent(in) :: old_bins(n_old+1)
  real(c_double), intent(in) :: old_vals(n_old)
  integer(c_int), intent(in) :: n_new
  real(c_double), intent(in) :: new_bins(n_new+1)
  real(c_double), intent(out) :: new_vals(n_new)
  integer(c_int), intent(out) :: ierr
  call rebin(old_bins, old_vals, new_bins, new_vals, ierr)
end subroutine