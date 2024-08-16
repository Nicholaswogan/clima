! futils

subroutine futils_rebin_wrapper(old_bins_len, old_bins, old_vals_len, old_vals, &
                                new_bins_len, new_bins, new_vals_len, new_vals, ierr) bind(c)
  use futils, only: rebin
  integer(c_int), intent(in) :: old_bins_len
  real(c_double), intent(in) :: old_bins(old_bins_len)
  integer(c_int), intent(in) :: old_vals_len
  real(c_double), intent(in) :: old_vals(old_vals_len)
  integer(c_int), intent(in) :: new_bins_len
  real(c_double), intent(in) :: new_bins(new_bins_len)
  integer(c_int), intent(in) :: new_vals_len
  real(c_double), intent(out) :: new_vals(new_vals_len)
  integer(c_int), intent(out) :: ierr
  call rebin(old_bins, old_vals, new_bins, new_vals, ierr)
end subroutine

subroutine futils_rebin_with_errors_wrapper(old_bins_len, old_bins, old_vals_len, old_vals, old_errs_len, old_errs, &
                                            new_bins_len, new_bins, new_vals_len, new_vals, new_errs_len, new_errs, ierr) bind(c)
  use futils, only: rebin_with_errors
  integer(c_int), intent(in) :: old_bins_len
  real(c_double), intent(in) :: old_bins(old_bins_len)
  integer(c_int), intent(in) :: old_vals_len
  real(c_double), intent(in) :: old_vals(old_vals_len)
  integer(c_int), intent(in) :: old_errs_len
  real(c_double), intent(in) :: old_errs(old_errs_len)
  integer(c_int), intent(in) :: new_bins_len
  real(c_double), intent(in) :: new_bins(new_bins_len)
  integer(c_int), intent(in) :: new_vals_len
  real(c_double), intent(out) :: new_vals(new_vals_len)
  integer(c_int), intent(in) :: new_errs_len
  real(c_double), intent(out) :: new_errs(new_errs_len)
  integer(c_int), intent(out) :: ierr
  call rebin_with_errors(old_bins, old_vals, old_errs, new_bins, new_vals, new_errs, ierr)
end subroutine

subroutine futils_grid_at_resolution_wrapper1(wv_min, wv_max, R, wavl_len, wavl_ptr, ierr) bind(c)
  use futils, only: grid_at_resolution, dp
  real(c_double), intent(in) :: wv_min
  real(c_double), intent(in) :: wv_max
  real(c_double), intent(in) :: R
  integer(c_int), intent(out) :: wavl_len
  type(c_ptr), intent(out) :: wavl_ptr
  integer(c_int), intent(out) :: ierr
  
  real(dp), allocatable :: wavl(:)
  real(c_double), pointer :: wavl_p(:)

  call grid_at_resolution(wv_min, wv_max, R, wavl, ierr)

  if (allocated(wavl)) then
    wavl_len = size(wavl)
    allocate(wavl_p(size(wavl)))
    wavl_p = wavl
    wavl_ptr = c_loc(wavl_p)
  else
    wavl_len = 0
    wavl_ptr = c_null_ptr
  endif

end subroutine

subroutine futils_grid_at_resolution_wrapper2(wavl_len, wavl_ptr, wavl) bind(c)
  integer(c_int), intent(in) :: wavl_len
  type(c_ptr), value, intent(in) :: wavl_ptr
  real(c_double), intent(out) :: wavl(wavl_len)
  
  real(c_double), pointer :: wavl_p(:)

  if (.not.c_associated(wavl_ptr)) return

  call c_f_pointer(wavl_ptr, wavl_p, [wavl_len])
  wavl = wavl_p
  deallocate(wavl_p)

end subroutine

subroutine futils_make_bins_wrapper(wv_len, wv, wavl_len, wavl, ierr) bind(c)
  use futils, only: make_bins
  integer(c_int), intent(in) :: wv_len
  real(c_double), intent(in) :: wv(wv_len)
  integer(c_int), intent(in) :: wavl_len
  real(c_double), intent(out) :: wavl(wavl_len)
  integer(c_int), intent(out) :: ierr
  call make_bins(wv, wavl, ierr)
end subroutine

subroutine futils_rebin_error_message_wrapper1(ierr, err_length, err_ptr) bind(c)
  use futils, only: rebin_error_message
  integer(c_int), intent(in) :: ierr
  integer(c_int), intent(out) :: err_length
  type(c_ptr), intent(out) :: err_ptr

  character(:), allocatable :: err
  character(kind=c_char), pointer :: err_p(:)

  err = rebin_error_message(ierr)

  if (allocated(err)) then
    err_length = len(err)
    allocate(err_p(err_length+1))
    call copy_string_ftoc(err, err_p)
    err_ptr = c_loc(err_p)
  else
    err_length = 0
    err_ptr = c_null_ptr
  endif

end subroutine

subroutine futils_rebin_error_message_wrapper2(err_length, err_ptr, err) bind(c)
  integer(c_int), intent(in) :: err_length
  type(c_ptr), value, intent(in) :: err_ptr
  character(c_char), intent(out) :: err(err_length+1)

  character(kind=c_char), pointer :: err_p(:)

  if (.not.c_associated(err_ptr)) return

  call c_f_pointer(err_ptr, err_p, [err_length+1]) 
  err = err_p
  deallocate(err_p)

end subroutine