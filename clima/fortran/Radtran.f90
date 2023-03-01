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
