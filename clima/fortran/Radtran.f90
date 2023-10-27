! Radtran

subroutine radtran_skin_temperature_wrapper(ptr, bond_albedo, T_skin) bind(c)
  use clima, only: Radtran
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: bond_albedo
  real(c_double), intent(out) :: T_skin
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  T_skin = rad%skin_temperature(bond_albedo)
end subroutine

subroutine radtran_opacities2yaml_wrapper_1(ptr, out_len, out_cp) bind(c)
  use clima, only: Radtran
  type(c_ptr), intent(in) :: ptr
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
  type(c_ptr), intent(in) :: ptr
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! getters and setters !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine radtran_surface_albedo_get(ptr, val) bind(c)
  use clima, only: Radtran
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  val = rad%surface_albedo
end subroutine

subroutine radtran_surface_albedo_set(ptr, val) bind(c)
  use clima, only: Radtran
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  rad%surface_albedo = val
end subroutine

subroutine radtran_ir_get(ptr, ptr1) bind(c)
  use clima, only: Radtran
  type(c_ptr), intent(in) :: ptr
  type(c_ptr), intent(out) :: ptr1
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  ptr1 = c_loc(rad%ir)
end subroutine

subroutine radtran_sol_get(ptr, ptr1) bind(c)
  use clima, only: Radtran
  type(c_ptr), intent(in) :: ptr
  type(c_ptr), intent(out) :: ptr1
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  ptr1 = c_loc(rad%sol)
end subroutine

subroutine radtran_wrk_ir_get(ptr, ptr1) bind(c)
  use clima, only: Radtran
  type(c_ptr), intent(in) :: ptr
  type(c_ptr), intent(out) :: ptr1
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  ptr1 = c_loc(rad%wrk_ir)
end subroutine

subroutine radtran_wrk_sol_get(ptr, ptr1) bind(c)
  use clima, only: Radtran
  type(c_ptr), intent(in) :: ptr
  type(c_ptr), intent(out) :: ptr1
  type(Radtran), pointer :: rad
  call c_f_pointer(ptr, rad)
  ptr1 = c_loc(rad%wrk_sol)
end subroutine
