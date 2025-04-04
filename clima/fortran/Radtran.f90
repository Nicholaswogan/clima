! Radtran

subroutine radtran_skin_temperature_wrapper(ptr, bond_albedo, T_skin) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  real(c_double), intent(in) :: bond_albedo
  real(c_double), intent(out) :: T_skin
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  T_skin = rad%skin_temperature(bond_albedo)
end subroutine

subroutine radtran_opacities2yaml_wrapper_1(ptr, out_len, out_cp) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: out_len
  type(c_ptr), intent(out) :: out_cp
  character(:), allocatable :: out_f
  character(kind=c_char), pointer :: out_cp1(:)
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)

  out_f = rad%opacities2yaml()
  out_len = len(out_f)
  allocate(out_cp1(out_len+1))
  call copy_string_ftoc(out_f, out_cp1)
  out_cp = c_loc(out_cp1)

end subroutine

subroutine radtran_opacities2yaml_wrapper_2(ptr, out_cp, out_len, out_c) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  type(c_ptr), intent(in) :: out_cp
  integer(c_int), intent(in) :: out_len
  character(kind=c_char), intent(out) :: out_c(out_len+1)
  character(kind=c_char), pointer :: out_cp1(:)
  integer :: i

  call c_f_pointer(out_cp, out_cp1, [out_len+1])

  do i = 1,out_len+1
    out_c(i) = out_cp1(i)
  enddo
  deallocate(out_cp1)

end subroutine

subroutine radtran_set_custom_optical_properties(ptr, dim_wv, wv, dim_P, P,  &
                                                 dim1_dtau_dz, dim2_dtau_dz, dtau_dz, &
                                                 dim1_w0, dim2_w0, w0, &
                                                 dim1_g0, dim2_g0, g0, &
                                                 err) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim_wv
  real(c_double), intent(in) :: wv(dim_wv)
  integer(c_int), intent(in) :: dim_P
  real(c_double), intent(in) :: P(dim_P)
  integer(c_int), intent(in) :: dim1_dtau_dz, dim2_dtau_dz
  real(c_double), intent(in) :: dtau_dz(dim1_dtau_dz, dim2_dtau_dz)
  integer(c_int), intent(in) :: dim1_w0, dim2_w0
  real(c_double), intent(in) :: w0(dim1_w0, dim2_w0)
  integer(c_int), intent(in) :: dim1_g0, dim2_g0
  real(c_double), intent(in) :: g0(dim1_g0, dim2_g0)
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(Radtran), pointer :: rad

  call c_f_pointer(ptr, rad)
  call rad%set_custom_optical_properties(wv, P, dtau_dz, w0, g0, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine radtran_unset_custom_optical_properties(ptr) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr

  type(Radtran), pointer :: rad

  call c_f_pointer(ptr, rad)
  call rad%unset_custom_optical_properties()

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! getters and setters !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine radtran_surface_albedo_get_size(ptr, dim1) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  dim1 = size(rad%surface_albedo,1)
end subroutine

subroutine radtran_surface_albedo_get(ptr, dim1, arr) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  arr = rad%surface_albedo
end subroutine

subroutine radtran_surface_albedo_set(ptr, dim1, arr) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(in) :: arr(dim1)
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  rad%surface_albedo = arr
end subroutine

subroutine radtran_surface_emissivity_get_size(ptr, dim1) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  dim1 = size(rad%surface_emissivity,1)
end subroutine

subroutine radtran_surface_emissivity_get(ptr, dim1, arr) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  arr = rad%surface_emissivity
end subroutine

subroutine radtran_surface_emissivity_set(ptr, dim1, arr) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(in) :: arr(dim1)
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  rad%surface_emissivity = arr
end subroutine

subroutine radtran_ir_get(ptr, ptr1) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  type(c_ptr), intent(out) :: ptr1
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  ptr1 = c_loc(rad%ir)
end subroutine

subroutine radtran_sol_get(ptr, ptr1) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  type(c_ptr), intent(out) :: ptr1
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  ptr1 = c_loc(rad%sol)
end subroutine

subroutine radtran_wrk_ir_get(ptr, ptr1) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  type(c_ptr), intent(out) :: ptr1
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  ptr1 = c_loc(rad%wrk_ir)
end subroutine

subroutine radtran_wrk_sol_get(ptr, ptr1) bind(c)
  use clima, only: Radtran
  type(c_ptr), value, intent(in) :: ptr
  type(c_ptr), intent(out) :: ptr1
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  ptr1 = c_loc(rad%wrk_sol)
end subroutine
