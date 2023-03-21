module clima_rc
  use clima_types, only: Species
  use clima_radtran, only: Radtran
  use clima_const, only: dp, s_str_len
  implicit none
  private

  public :: RadiativeConvectiveClimate

  type :: InitialConditions
    character(s_str_len), allocatable :: condensible_names(:)
    real(dp), allocatable :: condensible_P(:)
    real(dp), allocatable :: f_i(:,:)
    real(dp), allocatable :: T_init(:)
    real(dp), allocatable :: pdensities(:,:)
    real(dp), allocatable :: radii(:,:)
  end type
  interface InitialConditions
    module procedure :: create_InitialConditions
  end interface

  type :: RadiativeConvectiveClimate

    ! Species
    character(s_str_len), allocatable :: species_names(:) !! copy of gas names
    character(s_str_len), allocatable :: particle_names(:) !! copy of particle names
    type(Species) :: sp

    ! Radiative transfer
    integer :: nz !! number of layers.
    integer :: neq !! number of equations is nz + 1 for ground.
    integer :: nz_r !! (2*nz) number of layers for radiative calculations.
    type(Radtran) :: rad !! radiative transfer object

    ! Planet
    character(:), allocatable :: bg_gas !! background gas
    integer :: bg_gas_ind !! background gas index 
    real(dp) :: planet_mass !! grams
    real(dp) :: planet_radius !! cm
    real(dp) :: P_surf !! dynes/cm^2
    real(dp) :: P_top = 1.0e2_dp !! Pressure at top of the atmosphere (dynes/cm^2)

    type(InitialConditions), allocatable :: init
  
    ! State of the atmosphere
    real(dp), allocatable :: P(:) !! grid-center pressure (dynes/cm^2)
    real(dp), allocatable :: log10_P(:) !! log space 
    real(dp), allocatable :: log10_delta_P(:) !! thickness of each pressure layer (dynes/cm^2)
    real(dp), allocatable :: z(:) !! nz
    real(dp), allocatable :: dz(:) !! nz
    real(dp), allocatable :: f_i(:,:) !! mixing ratios (nz,ng)
    real(dp), allocatable :: densities(:,:) !! molecules/cm^3 (nz,ng)
    real(dp), allocatable :: pdensities(:,:) !! particle densities (particles/cm^3) (nz,np)
    real(dp), allocatable :: radii(:,:) !! particle radii (cm) (nz,np)
    real(dp), allocatable :: T(:) !! (nz)
    real(dp) :: T_surf

  contains
    procedure :: initialize_stepper => RadiativeConvectiveClimate_initialize_stepper

  end type
  interface RadiativeConvectiveClimate
    module procedure :: create_RadiativeConvectiveClimate
  end interface

contains

  function create_RadiativeConvectiveClimate(datadir, species_f, settings_f, star_f, err) result(c)
    use futils, only: linspace
    use clima_types, only: ClimaSettings, AtmosphereFile, unpack_atmospherefile
    use clima_eqns, only: vertical_grid, gravity_z
    character(*), intent(in) :: datadir
    character(*), intent(in) :: species_f
    character(*), intent(in) :: settings_f
    character(*), intent(in) :: star_f
    character(:), allocatable, intent(out) :: err

    type(RadiativeConvectiveClimate) :: c
    
    integer :: i
    character(s_str_len) :: bg_gas_copy
    type(ClimaSettings) :: s

    ! create settings
    s = ClimaSettings(settings_f, err)
    if (allocated(err)) return

    ! Check settings for the proper inputs
    call check_RadiativeConvectiveClimate_settings(s, err)
    if (allocated(err)) return

    ! unpack the settings
    c%nz = s%nz
    c%nz_r  = c%nz*2
    c%neq = c%nz + 1
    c%bg_gas = s%back_gas_name
    c%planet_mass = s%planet_mass
    c%planet_radius = s%planet_radius
    c%P_surf = s%P_surf*1.0e6_dp ! bars to dynes/cm^2 conversion

    c%sp = Species(species_f, err)
    if (allocated(err)) return

    ! make copy of species names
    allocate(c%species_names(c%sp%ng))
    do i = 1,c%sp%ng
      c%species_names(i) = c%sp%g(i)%name
    enddo

    ! make copy of particle names
    allocate(c%particle_names(c%sp%np))
    do i = 1,c%sp%np
      c%particle_names(i) = c%sp%p(i)%name
    enddo

    ! background gas index
    bg_gas_copy = c%bg_gas
    c%bg_gas_ind = findloc(c%species_names, bg_gas_copy, 1)
    if (c%bg_gas_ind == 0) then
      err = 'Background gas "'//c%bg_gas//'" is not in the list of species.'
      return
    endif

    ! create radiative transfer
    c%rad = Radtran(datadir, c%species_names, c%particle_names, s, star_f, &
                    s%number_of_zenith_angles, s%surface_albedo, c%nz_r, err)
    if (allocated(err)) return

    ! Pressure grid
    allocate(c%P(c%nz))
    allocate(c%log10_P(c%nz))
    allocate(c%log10_delta_P(c%nz))
    c%log10_delta_P = (log10(c%P_surf) - log10(c%P_top))/c%nz
    c%log10_P(1) = log10(c%P_surf) - 0.5_dp*c%log10_delta_P(1)
    do i = 2,c%nz
      c%log10_P(i) = c%log10_P(i-1) - c%log10_delta_P(1)
    enddo
    c%P = 10.0_dp**c%log10_P

    ! allocate
    allocate(c%z(c%nz))
    allocate(c%dz(c%nz))
    allocate(c%f_i(c%nz,c%sp%ng))
    allocate(c%densities(c%nz,c%sp%ng))
    allocate(c%pdensities(c%nz,c%sp%np))
    allocate(c%radii(c%nz,c%sp%np))
    allocate(c%T(c%nz))

  end function

  subroutine check_RadiativeConvectiveClimate_settings(s, err)
    use clima_types, only: ClimaSettings
    type(ClimaSettings), intent(in) :: s
    character(:), allocatable, intent(out) :: err

    if (.not. s%atmos_grid_is_present) then
      err = '"'//s%filename//'/atmosphere-grid" does not exist.'
      return
    endif
    
    if (.not. s%planet_is_present) then
      err = '"'//s%filename//'/planet" does not exist.'
      return
    endif

    if (.not. allocated(s%back_gas_name)) then
      err = '"'//s%filename//'/planet/background-gas" does not exist.'
    endif

    if (.not. allocated(s%P_surf)) then
      err = '"'//s%filename//'/planet/surface-pressure" does not exist.'
      return
    endif

    if (.not. allocated(s%surface_albedo)) then
      err = '"'//s%filename//'/planet/surface-albedo" does not exist.'
      return
    endif

    if (.not. allocated(s%number_of_zenith_angles)) then
      err = '"'//s%filename//'/planet/number-of-zenith-angles" does not exist.'
      return
    endif

  end subroutine

  function create_InitialConditions(condensible_names, condensible_P, f_i, pdensities, T_init, radii) result(init)
    character(*), intent(in) :: condensible_names(:)
    real(dp), intent(in) :: condensible_P(:)
    real(dp), intent(in) :: f_i(:,:)
    real(dp), intent(in) :: T_init(:)
    real(dp), optional, intent(in) :: pdensities(:,:)
    real(dp), optional, intent(in) :: radii(:,:)
    type(InitialConditions) :: init

    init%condensible_names = condensible_names
    init%condensible_P = condensible_P
    init%f_i = f_i
    init%T_init = T_init
    if (present(pdensities)) then
      init%pdensities = pdensities
    endif 
    if (present(radii)) then
      init%radii = radii
    endif 

  end function

  subroutine RadiativeConvectiveClimate_initialize_stepper(self, condensible_names, condensible_P, &
                                                           f_i, pdensities, radii, err)
    use clima_rc_adiabat, only: make_profile_rc
    class(RadiativeConvectiveClimate), intent(inout) :: self
    character(*), intent(in) :: condensible_names(:)
    real(dp), intent(in) :: condensible_P(:)
    real(dp), intent(in) :: f_i(:,:)
    real(dp), optional, intent(in) :: pdensities(:,:)
    real(dp), optional, intent(in) :: radii(:,:)
    character(:), allocatable, intent(out) :: err

    integer :: i

    ! check inputs
    call check_inputs(self%nz, self%sp%ng, self%sp%np, f_i, pdensities, radii, err)
    if (allocated(err)) return

    self%f_i = f_i

    ! now draw initial T-profile
    block
    real(dp) :: T_surf, T_trop
    real(dp), allocatable :: RH(:)
    integer, allocatable :: cond_inds(:)
    integer, allocatable :: cond_z_inds(:,:)

    T_surf = 280.0_dp
    T_trop = 1.0_dp
    allocate(RH(size(condensible_names)))
    allocate(cond_inds(size(condensible_names)))
    allocate(cond_z_inds(2,size(condensible_names)))

    RH = 1.0_dp
    do i = 1,size(condensible_names)
      cond_inds(i) = findloc(self%species_names, condensible_names(i), 1)
    enddo

    call make_profile_rc( &
      self%sp, self%planet_mass, self%planet_radius, &
      self%P_surf, T_surf, T_trop, self%nz, self%P, self%log10_delta_P, &
      RH, condensible_P, cond_inds, self%bg_gas_ind, &
      self%f_i, cond_z_inds, self%z, self%dz, self%T, err &
    ) 
    if (allocated(err)) return

    end block

    ! if (allocated(self%init)) deallocate(self%init)
    ! allocate(self%init)
    ! self%init = InitialConditions(condensible_names, condensible_P, f_i, pdensities, T_init, radii)

  end subroutine

  ! function RadiativeConvectiveClimate_step(self) result(tn)
  !   class(RadiativeConvectiveClimate), intent(inout) :: self
  !   real(dp) :: tn
  !   tn = 1.0_dp
  ! end function

  subroutine check_inputs(nz, ng, np, f_i, pdensities, radii, err)
    integer, intent(in) :: nz, ng, np
    real(dp), intent(in) :: f_i(:,:)
    real(dp), optional, intent(in) :: pdensities(:,:), radii(:,:)
    character(:), allocatable, intent(out) :: err

    if (size(f_i,1) /= nz .or. size(f_i,2) /= ng) then
      err = 'Input "f_i" has the wrong dimensions.'
      return
    endif

    if ((present(pdensities) .and. .not. present(radii)) &
      .or.(present(radii) .and. .not. present(pdensities))) then
      err = 'Both pdensities and radii must be arguments.'
      return
    endif
    if (np > 0) then
      if (.not. present(radii)) then
        err = 'The model contains particles but "pdensities" and "radii" are not arguments.'
        return
      endif
    endif
    if (present(radii)) then
      if (size(pdensities,1) /= nz .or. size(pdensities,2) /= np) then
        err = '"pdensities" has the wrong input dimension.'
        return
      endif
      if (size(radii,1) /= nz .or. size(radii,2) /= np) then
        err = '"radii" has the wrong input dimension.'
        return
      endif
    endif

  end subroutine

end module