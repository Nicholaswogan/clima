module WaterAdiabatClimate_wrapper
  use clima, only: WaterAdiabatClimate
  use wrapper_utils, only: copy_string_ftoc, copy_string_ctof, len_cstring, err_len
  use iso_c_binding
  implicit none

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! allocator and destroyer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine allocate_wateradiabatclimate(ptr) bind(c)
    type(c_ptr), intent(out) :: ptr
    type(WaterAdiabatClimate), pointer :: pc
    allocate(pc)
    ptr = c_loc(pc)
  end subroutine

  subroutine deallocate_wateradiabatclimate(ptr) bind(c)
    type(c_ptr), intent(in) :: ptr
    type(WaterAdiabatClimate), pointer :: pc
    character(:), allocatable :: err_f
    call c_f_pointer(ptr, pc)
    deallocate(pc)
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! subroutine wrappers  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine wateradiabatclimate_create_wrapper(ptr, data_dir, species_file, &
                                     settings_file, flux_file, err) bind(c)
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
    type(WaterAdiabatClimate), pointer :: c
    
    call c_f_pointer(ptr, c)
    
    allocate(character(len=len_cstring(data_dir))::data_dir_f)
    allocate(character(len=len_cstring(species_file))::species_file_f)
    allocate(character(len=len_cstring(settings_file))::settings_file_f)
    allocate(character(len=len_cstring(flux_file))::flux_file_f)
    
    call copy_string_ctof(data_dir, data_dir_f)
    call copy_string_ctof(species_file, species_file_f)
    call copy_string_ctof(settings_file, settings_file_f)
    call copy_string_ctof(flux_file, flux_file_f)
    
    c = WaterAdiabatClimate(data_dir_f, &
                            species_file_f, &
                            settings_file_f, &
                            flux_file_f, &
                            err_f)
    
    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif
  end subroutine

  subroutine wateradiabatclimate_make_profile_wrapper(ptr, T_surf, ng, P_i_surf, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: T_surf
    integer(c_int), intent(in) :: ng
    real(c_double), intent(in) :: P_i_surf(ng)
    character(c_char), intent(out) :: err(err_len+1)

    character(:), allocatable :: err_f
    type(WaterAdiabatClimate), pointer :: c

    call c_f_pointer(ptr, c)

    call c%make_profile(T_surf, P_i_surf, err_f)

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif

  end subroutine

  subroutine wateradiabatclimate_make_column_wrapper(ptr, T_surf, ng, N_i_surf, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: T_surf
    integer(c_int), intent(in) :: ng
    real(c_double), intent(in) :: N_i_surf(ng)
    character(c_char), intent(out) :: err(err_len+1)

    character(:), allocatable :: err_f
    type(WaterAdiabatClimate), pointer :: c

    call c_f_pointer(ptr, c)

    call c%make_column(T_surf, N_i_surf, err_f)

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif

  end subroutine

  subroutine wateradiabatclimate_toa_fluxes_wrapper(ptr, T_surf, ng, P_i_surf, ISR, OLR, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: T_surf
    integer(c_int), intent(in) :: ng
    real(c_double), intent(in) :: P_i_surf(ng)
    real(c_double), intent(out) :: ISR, OLR
    character(c_char), intent(out) :: err(err_len+1)

    character(:), allocatable :: err_f
    type(WaterAdiabatClimate), pointer :: c

    call c_f_pointer(ptr, c)

    call c%TOA_fluxes(T_surf, P_i_surf, ISR, OLR, err_f)

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif

  end subroutine

  subroutine wateradiabatclimate_toa_fluxes_column_wrapper(ptr, T_surf, ng, N_i_surf, ISR, OLR, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: T_surf
    integer(c_int), intent(in) :: ng
    real(c_double), intent(in) :: N_i_surf(ng)
    real(c_double), intent(out) :: ISR, OLR
    character(c_char), intent(out) :: err(err_len+1)

    character(:), allocatable :: err_f
    type(WaterAdiabatClimate), pointer :: c

    call c_f_pointer(ptr, c)

    call c%TOA_fluxes_column(T_surf, N_i_surf, ISR, OLR, err_f)

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif

  end subroutine

  subroutine wateradiabatclimate_surface_temperature_wrapper(ptr, ng, P_i_surf, T_guess, T_surf, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: ng
    real(c_double), intent(in) :: P_i_surf(ng)
    real(c_double), intent(in) :: T_guess
    real(c_double), intent(out) :: T_surf
    character(c_char), intent(out) :: err(err_len+1)

    character(:), allocatable :: err_f
    type(WaterAdiabatClimate), pointer :: c

    call c_f_pointer(ptr, c)

    T_surf = c%surface_temperature(P_i_surf, T_guess, err_f)

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif

  end subroutine

  subroutine wateradiabatclimate_surface_temperature_column_wrapper(ptr, ng, N_i_surf, T_guess, T_surf, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: ng
    real(c_double), intent(in) :: N_i_surf(ng)
    real(c_double), intent(in) :: T_guess
    real(c_double), intent(out) :: T_surf
    character(c_char), intent(out) :: err(err_len+1)

    character(:), allocatable :: err_f
    type(WaterAdiabatClimate), pointer :: c

    call c_f_pointer(ptr, c)

    T_surf = c%surface_temperature_column(N_i_surf, T_guess, err_f)

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif

  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! getters and setters !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wateradiabatclimate_p_top_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    val = c%P_top
  end subroutine
  
  subroutine wateradiabatclimate_p_top_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    c%P_top = val
  end subroutine

  subroutine wateradiabatclimate_t_trop_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    val = c%T_trop
  end subroutine
  
  subroutine wateradiabatclimate_t_trop_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    c%T_trop = val
  end subroutine

  subroutine wateradiabatclimate_rh_get(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    val = c%RH
  end subroutine
  
  subroutine wateradiabatclimate_rh_set(ptr, val) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    c%RH = val
  end subroutine

  subroutine wateradiabatclimate_p_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    dim1 = size(c%P,1)
  end subroutine
  
  subroutine wateradiabatclimate_p_get(ptr, dim1, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    arr = c%P
  end subroutine

  subroutine wateradiabatclimate_t_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    dim1 = size(c%T,1)
  end subroutine
  
  subroutine wateradiabatclimate_t_get(ptr, dim1, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    arr = c%T
  end subroutine

  subroutine wateradiabatclimate_f_i_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    dim1 = size(c%f_i,1)
    dim2 = size(c%f_i,2)
  end subroutine
  
  subroutine wateradiabatclimate_f_i_get(ptr, dim1, dim2, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: arr(dim1,dim2)
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    arr = c%f_i
  end subroutine

  subroutine wateradiabatclimate_z_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    dim1 = size(c%z,1)
  end subroutine
  
  subroutine wateradiabatclimate_z_get(ptr, dim1, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    arr = c%z
  end subroutine

  subroutine wateradiabatclimate_dz_get_size(ptr, dim1) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    dim1 = size(c%dz,1)
  end subroutine
  
  subroutine wateradiabatclimate_dz_get(ptr, dim1, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    arr = c%dz
  end subroutine

  subroutine wateradiabatclimate_densities_get_size(ptr, dim1, dim2) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(out) :: dim1, dim2
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    dim1 = size(c%densities,1)
    dim2 = size(c%densities,2)
  end subroutine
  
  subroutine wateradiabatclimate_densities_get(ptr, dim1, dim2, arr) bind(c)
    type(c_ptr), intent(in) :: ptr
    integer(c_int), intent(in) :: dim1, dim2
    real(c_double), intent(out) :: arr(dim1,dim2)
    type(WaterAdiabatClimate), pointer :: c
    call c_f_pointer(ptr, c)
    arr = c%densities
  end subroutine

end module