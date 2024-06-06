! OpticalProperties

subroutine opticalproperties_wavl_get_size(ptr, dim1) bind(c)
  use clima_radtran_types, only: OpticalProperties
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(OpticalProperties), pointer :: op
  call c_f_pointer(ptr, op)
  dim1 = size(op%wavl,1)
end subroutine

subroutine opticalproperties_wavl_get(ptr, dim1, arr) bind(c)
  use clima_radtran_types, only: OpticalProperties
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(OpticalProperties), pointer :: op
  call c_f_pointer(ptr, op)
  arr = op%wavl
end subroutine

subroutine opticalproperties_freq_get_size(ptr, dim1) bind(c)
  use clima_radtran_types, only: OpticalProperties
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(OpticalProperties), pointer :: op
  call c_f_pointer(ptr, op)
  dim1 = size(op%freq,1)
end subroutine

subroutine opticalproperties_freq_get(ptr, dim1, arr) bind(c)
  use clima_radtran_types, only: OpticalProperties
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(OpticalProperties), pointer :: op
  call c_f_pointer(ptr, op)
  arr = op%freq
end subroutine