! Radtran

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! getters and setters !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
