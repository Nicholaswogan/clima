! AdiabatClimate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! allocator and destroyer !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine allocate_adiabatclimate(ptr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(out) :: ptr
  type(AdiabatClimate), pointer :: pc
  allocate(pc)
  ptr = c_loc(pc)
end subroutine

subroutine deallocate_adiabatclimate(ptr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  type(AdiabatClimate), pointer :: pc
  character(:), allocatable :: err_f
  call c_f_pointer(ptr, pc)
  deallocate(pc)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine wrappers  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine adiabatclimate_create_wrapper(ptr, data_dir, species_file, &
                                    settings_file, flux_file, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  character(kind=c_char), intent(in) :: data_dir(*)
  character(kind=c_char), intent(in) :: species_file(*)
  character(kind=c_char), intent(in) :: settings_file(*)
  character(kind=c_char), intent(in) :: flux_file(*)
  character(kind=c_char), intent(out) :: err(err_len+1)
  
  character(len=:), allocatable :: data_dir_f
  character(len=:), allocatable :: species_file_f
  character(len=:), allocatable :: settings_file_f
  character(len=:), allocatable :: flux_file_f
  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c
  
  call c_f_pointer(ptr, c)
  
  allocate(character(len=len_cstring(data_dir))::data_dir_f)
  allocate(character(len=len_cstring(species_file))::species_file_f)
  allocate(character(len=len_cstring(settings_file))::settings_file_f)
  allocate(character(len=len_cstring(flux_file))::flux_file_f)
  
  call copy_string_ctof(data_dir, data_dir_f)
  call copy_string_ctof(species_file, species_file_f)
  call copy_string_ctof(settings_file, settings_file_f)
  call copy_string_ctof(flux_file, flux_file_f)
  
  c = AdiabatClimate(data_dir_f, &
                          species_file_f, &
                          settings_file_f, &
                          flux_file_f, &
                          err_f)
  
  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif
end subroutine

subroutine adiabatclimate_make_profile_wrapper(ptr, T_surf, ng, P_i_surf, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: T_surf
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: P_i_surf(ng)
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  call c%make_profile(T_surf, P_i_surf, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_make_column_wrapper(ptr, T_surf, ng, N_i_surf, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: T_surf
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: N_i_surf(ng)
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  call c%make_column(T_surf, N_i_surf, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_make_profile_bg_gas_wrapper(ptr, T_surf, ng, P_i_surf, P_surf, bg_gas, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: T_surf
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: P_i_surf(ng)
  real(c_double), intent(in) :: P_surf
  character(kind=c_char), intent(in) :: bg_gas(*)
  character(c_char), intent(out) :: err(err_len+1)

  character(len=:), allocatable :: bg_gas_f
  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  allocate(character(len=len_cstring(bg_gas))::bg_gas_f)
  call copy_string_ctof(bg_gas, bg_gas_f)

  call c%make_profile_bg_gas(T_surf, P_i_surf, P_surf, bg_gas_f, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_toa_fluxes_wrapper(ptr, T_surf, ng, P_i_surf, ISR, OLR, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: T_surf
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: P_i_surf(ng)
  real(c_double), intent(out) :: ISR, OLR
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  call c%TOA_fluxes(T_surf, P_i_surf, ISR, OLR, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_toa_fluxes_column_wrapper(ptr, T_surf, ng, N_i_surf, ISR, OLR, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: T_surf
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: N_i_surf(ng)
  real(c_double), intent(out) :: ISR, OLR
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  call c%TOA_fluxes_column(T_surf, N_i_surf, ISR, OLR, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_toa_fluxes_bg_gas_wrapper(ptr, T_surf, ng, P_i_surf, P_surf, bg_gas, ISR, OLR, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: T_surf
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: P_i_surf(ng)
  real(c_double), intent(in) :: P_surf
  character(kind=c_char), intent(in) :: bg_gas(*)
  real(c_double), intent(out) :: ISR, OLR
  character(c_char), intent(out) :: err(err_len+1)

  character(len=:), allocatable :: bg_gas_f
  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  allocate(character(len=len_cstring(bg_gas))::bg_gas_f)
  call copy_string_ctof(bg_gas, bg_gas_f)

  call c%TOA_fluxes_bg_gas(T_surf, P_i_surf, P_surf, bg_gas_f, ISR, OLR, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_surface_temperature_wrapper(ptr, ng, P_i_surf, T_guess, T_surf, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: P_i_surf(ng)
  real(c_double), intent(in) :: T_guess
  real(c_double), intent(out) :: T_surf
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  T_surf = c%surface_temperature(P_i_surf, T_guess, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_surface_temperature_column_wrapper(ptr, ng, N_i_surf, T_guess, T_surf, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: N_i_surf(ng)
  real(c_double), intent(in) :: T_guess
  real(c_double), intent(out) :: T_surf
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  T_surf = c%surface_temperature_column(N_i_surf, T_guess, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_surface_temperature_bg_gas_wrapper(ptr, ng, P_i_surf, P_surf, bg_gas, T_guess, T_surf, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(in) :: ng
  real(c_double), intent(in) :: P_i_surf(ng)
  real(c_double), intent(in) :: P_surf
  character(kind=c_char), intent(in) :: bg_gas(*)
  real(c_double), intent(in) :: T_guess
  real(c_double), intent(out) :: T_surf
  character(c_char), intent(out) :: err(err_len+1)

  character(len=:), allocatable :: bg_gas_f
  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  allocate(character(len=len_cstring(bg_gas))::bg_gas_f)
  call copy_string_ctof(bg_gas, bg_gas_f)

  T_surf = c%surface_temperature_bg_gas(P_i_surf, P_surf, bg_gas_f, T_guess, err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_to_regular_grid_wrapper(ptr, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  character(c_char), intent(out) :: err(err_len+1)

  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)

  call c%to_regular_grid(err_f)

  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif

end subroutine

subroutine adiabatclimate_out2atmosphere_txt_wrapper(ptr, filename, nz, eddy, overwrite, clip, err) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  character(kind=c_char), intent(in) :: filename(*)
  integer(c_int), intent(in) :: nz
  real(c_double), intent(in) :: eddy(nz)
  logical(c_bool), intent(in) :: overwrite, clip
  character(kind=c_char), intent(out) :: err(err_len+1)
  
  character(len=:), allocatable :: filename_f
  logical :: overwrite_f, clip_f
  character(:), allocatable :: err_f
  type(AdiabatClimate), pointer :: c
  
  call c_f_pointer(ptr, c)
  
  allocate(character(len=len_cstring(filename))::filename_f)
  call copy_string_ctof(filename, filename_f)
  overwrite_f = overwrite
  clip_f = clip
  
  call c%out2atmosphere_txt(filename_f, eddy, overwrite_f, clip_f, err_f)
  err(1) = c_null_char
  if (allocated(err_f)) then
    call copy_string_ftoc(err_f, err)
  endif
  
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! getters and setters !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine adiabatclimate_p_top_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%P_top
end subroutine

subroutine adiabatclimate_p_top_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%P_top = val
end subroutine

subroutine adiabatclimate_t_trop_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%T_trop
end subroutine

subroutine adiabatclimate_t_trop_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%T_trop = val
end subroutine

subroutine adiabatclimate_solve_for_t_trop_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  logical(c_bool), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%solve_for_T_trop
end subroutine

subroutine adiabatclimate_solve_for_t_trop_set(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  logical(c_bool), intent(in) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%solve_for_T_trop = val
end subroutine

subroutine adiabatclimate_albedo_fcn_set(ptr, set_to_null, albedo_fcn_c) bind(c)
  use clima, only: AdiabatClimate
  use clima_adiabat, only: temp_dependent_albedo_fcn
  type(c_ptr), intent(in) :: ptr
  logical(c_bool), intent(in) :: set_to_null
  type(c_funptr), value, intent(in) :: albedo_fcn_c

  procedure(temp_dependent_albedo_fcn), pointer :: albedo_fcn_f
  type(AdiabatClimate), pointer :: c

  call c_f_pointer(ptr, c)
  if (set_to_null) then
    c%albedo_fcn => null()
  else
    call c_f_procpointer(albedo_fcn_c, albedo_fcn_f)
    c%albedo_fcn => albedo_fcn_f
  endif
end subroutine

subroutine adiabatclimate_rh_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%RH,1)
end subroutine

subroutine adiabatclimate_rh_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%RH
end subroutine

subroutine adiabatclimate_rh_set(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(in) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  c%RH = arr
end subroutine

subroutine adiabatclimate_species_names_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%species_names)
end subroutine

subroutine adiabatclimate_species_names_get(ptr, dim1, species_names) bind(c)
  use clima, only: AdiabatClimate, s_str_len
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  character(kind=c_char), intent(out) :: species_names(dim1*s_str_len+1)
  type(AdiabatClimate), pointer :: c
  
  integer :: i, j, k
  
  call c_f_pointer(ptr, c)
  do i = 1,dim1
    do j = 1,s_str_len
      k = j + (i - 1) * s_str_len
      species_names(k) = c%species_names(i)(j:j)
    enddo
  enddo
  species_names(dim1*s_str_len+1) = c_null_char
  
end subroutine

subroutine adiabatclimate_p_surf_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%P_surf
end subroutine

subroutine adiabatclimate_p_trop_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%P_trop
end subroutine

subroutine adiabatclimate_p_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%P,1)
end subroutine

subroutine adiabatclimate_p_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%P
end subroutine

subroutine adiabatclimate_t_surf_get(ptr, val) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  real(c_double), intent(out) :: val
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  val = c%T_surf
end subroutine

subroutine adiabatclimate_t_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%T,1)
end subroutine

subroutine adiabatclimate_t_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%T
end subroutine

subroutine adiabatclimate_f_i_get_size(ptr, dim1, dim2) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(out) :: dim1, dim2
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%f_i,1)
  dim2 = size(c%f_i,2)
end subroutine

subroutine adiabatclimate_f_i_get(ptr, dim1, dim2, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(in) :: dim1, dim2
  real(c_double), intent(out) :: arr(dim1,dim2)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%f_i
end subroutine

subroutine adiabatclimate_z_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%z,1)
end subroutine

subroutine adiabatclimate_z_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%z
end subroutine

subroutine adiabatclimate_dz_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%dz,1)
end subroutine

subroutine adiabatclimate_dz_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%dz
end subroutine

subroutine adiabatclimate_densities_get_size(ptr, dim1, dim2) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(out) :: dim1, dim2
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%densities,1)
  dim2 = size(c%densities,2)
end subroutine

subroutine adiabatclimate_n_surface_get_size(ptr, dim1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(out) :: dim1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  dim1 = size(c%N_surface,1)
end subroutine

subroutine adiabatclimate_n_surface_get(ptr, dim1, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(in) :: dim1
  real(c_double), intent(out) :: arr(dim1)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%N_surface
end subroutine

subroutine adiabatclimate_densities_get(ptr, dim1, dim2, arr) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  integer(c_int), intent(in) :: dim1, dim2
  real(c_double), intent(out) :: arr(dim1,dim2)
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  arr = c%densities
end subroutine

subroutine adiabatclimate_rad_get(ptr, ptr1) bind(c)
  use clima, only: AdiabatClimate
  type(c_ptr), intent(in) :: ptr
  type(c_ptr), intent(out) :: ptr1
  type(AdiabatClimate), pointer :: c
  call c_f_pointer(ptr, c)
  ptr1 = c_loc(c%rad)
end subroutine