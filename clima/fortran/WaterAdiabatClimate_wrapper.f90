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

  subroutine wateradiabatclimate_olr_wrapper(ptr, T_surf, ng, P_i_surf, OLR, err) bind(c)
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: T_surf
    integer(c_int), intent(in) :: ng
    real(c_double), intent(in) :: P_i_surf(ng)
    real(c_double), intent(out) :: OLR
    character(c_char), intent(out) :: err(err_len+1)

    character(:), allocatable :: err_f
    type(WaterAdiabatClimate), pointer :: c

    call c_f_pointer(ptr, c)

    OLR = c%OLR(T_surf, P_i_surf, err_f)

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

end module