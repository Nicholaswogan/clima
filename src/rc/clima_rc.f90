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
    real(dp) :: T_trop = 200.0_dp !! in/out
    real(dp) :: P_trop
    integer :: trop_ind

    ! Set
    character(s_str_len), allocatable :: condensible_names(:)
    integer, allocatable :: condensible_inds(:)
    real(dp), allocatable :: condensible_P(:)
    integer, allocatable :: condensible_z_inds(:)

    ! work variables for radiative transfer
    real(dp), allocatable :: T_r(:) !! (nz_r)
    real(dp), allocatable :: P_r(:) !! (nz_r)
    real(dp), allocatable :: densities_r(:,:) !! (nz_r,ng)
    real(dp), allocatable :: dz_r(:) !! nz_r
    real(dp), allocatable :: pdensities_r(:,:) !! (nz_r,np)
    real(dp), allocatable :: radii_r(:,:) !! (nz_r,np)

  contains
    procedure :: make_profile => RadiativeConvectiveClimate_make_profile
    procedure :: TOA_fluxes => RadiativeConvectiveClimate_TOA_fluxes
    procedure :: surface_temperature => RadiativeConvectiveClimate_surface_temperature

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

    allocate(c%T_r(c%nz_r))
    allocate(c%P_r(c%nz_r))
    allocate(c%densities_r(c%nz_r,c%sp%ng))
    allocate(c%dz_r(c%nz_r))
    allocate(c%pdensities_r(c%nz_r,c%sp%np))
    allocate(c%radii_r(c%nz_r,c%sp%np))

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

  subroutine RadiativeConvectiveClimate_make_profile(self, T_surf, condensible_names, condensible_P, condensible_RH, &
                                                     f_i, err)
    use clima_rc_adiabat, only: make_profile_rc
    use clima_const, only: k_boltz
    class(RadiativeConvectiveClimate), intent(inout) :: self
    real(dp), intent(in) :: T_surf
    character(*), intent(in) :: condensible_names(:)
    real(dp), intent(in) :: condensible_P(:)
    real(dp), intent(in) :: condensible_RH(:)
    real(dp), intent(in) :: f_i(:,:)
    character(:), allocatable, intent(out) :: err

    integer :: i
    real(dp), allocatable :: density(:)

    self%T_surf = T_surf
    self%f_i = f_i

    if (allocated(self%condensible_names)) then
      deallocate(self%condensible_names)
      deallocate(self%condensible_inds)
      deallocate(self%condensible_P)
      deallocate(self%condensible_z_inds)
    endif
    allocate(self%condensible_names(size(condensible_names)))
    allocate(self%condensible_inds(size(condensible_names)))
    allocate(self%condensible_P(size(condensible_names)))
    allocate(self%condensible_z_inds(size(condensible_names)))
    self%condensible_names = condensible_names
    self%condensible_P = condensible_P
    do i = 1,size(condensible_names)
      self%condensible_inds(i) = findloc(self%species_names, condensible_names(i), 1)
      if (self%condensible_inds(i) == 0) then
        err = 'Condensible input"'//trim(condensible_names(i))//'" is not in the list of species.'
        return
      endif
    enddo

    call make_profile_rc( &
      self%sp, self%planet_mass, self%planet_radius, &
      self%P_surf, T_surf, self%T_trop, self%nz, self%P, self%log10_delta_P, &
      condensible_RH, condensible_P, self%condensible_inds, self%bg_gas_ind, &
      self%f_i, self%P_trop, self%condensible_z_inds, self%z, self%dz, self%T, err &
    ) 
    if (allocated(err)) return

    self%trop_ind = -1
    if (self%P_trop > 0.0_dp) then
      do i = 1,self%nz
        if (self%P(i) < self%P_trop) then
          self%trop_ind = i - 1
          exit
        endif
      enddo
    endif

    allocate(density(self%nz))
    density = self%P/(k_boltz*self%T)
    do i = 1,self%sp%ng
      self%densities(:,i) = self%f_i(:,i)*density(:)
    enddo

  end subroutine

  subroutine RadiativeConvectiveClimate_TOA_fluxes(self, T_surf, condensible_names, condensible_P, condensible_RH, &
                                                   f_i, pdensities, radii, ISR, OLR, err)
    use clima_rc_adiabat, only: make_profile_rc
    use clima_const, only: k_boltz
    class(RadiativeConvectiveClimate), intent(inout) :: self
    real(dp), intent(in) :: T_surf
    character(*), intent(in) :: condensible_names(:)
    real(dp), intent(in) :: condensible_P(:)
    real(dp), intent(in) :: condensible_RH(:)
    real(dp), intent(in) :: f_i(:,:)
    real(dp), optional, intent(in) :: pdensities(:,:)
    real(dp), optional, intent(in) :: radii(:,:)
    real(dp), intent(out) :: ISR, OLR
    character(:), allocatable, intent(out) :: err

    integer :: i,j

    ! check inputs
    call check_inputs(self%nz, self%sp%ng, self%sp%np, f_i, pdensities, radii, err)
    if (allocated(err)) return

    ! make profile
    call self%make_profile(T_surf, condensible_names, condensible_P, condensible_RH, f_i, err)
    if (allocated(err)) return

    ! copy atmosphere to radiative transfer variables
    do i = 1,self%nz
      j = 2*(i-1)
      self%T_r(j+1) = self%T(i)
      self%T_r(j+2) = self%T(i)

      self%P_r(j+1) = self%P(i)
      self%P_r(j+2) = self%P(i)

      self%densities_r(j+1,:) = self%densities(i,:)
      self%densities_r(j+2,:) = self%densities(i,:)

      self%dz_r(j+1) = 0.5_dp*self%dz(i)
      self%dz_r(j+2) = 0.5_dp*self%dz(i)
    enddo
    if (present(pdensities)) then
      self%pdensities = pdensities
      self%radii = radii
      do i = 1,self%nz
        j = 2*(i-1)
        self%pdensities_r(j+1,:) = self%pdensities(i,:)
        self%pdensities_r(j+2,:) = self%pdensities(i,:)

        self%radii_r(j+1,:) = self%radii(i,:)
        self%radii_r(j+2,:) = self%radii(i,:)
      enddo
    endif

    ! Do radiative transfer
    if (present(pdensities)) then
      call self%rad%TOA_fluxes(T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, &
                               self%pdensities_r, self%radii_r, ISR=ISR, OLR=OLR, err=err)
      if (allocated(err)) return
    else
      call self%rad%TOA_fluxes(T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, &
                               ISR=ISR, OLR=OLR, err=err)
      if (allocated(err)) return
    endif

  end subroutine

  function RadiativeConvectiveClimate_surface_temperature(self, condensible_names, condensible_P, condensible_RH, &
                                                            f_i, pdensities, radii, T_guess, err) result(T_surf)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars
    class(RadiativeConvectiveClimate), intent(inout) :: self
    character(*), intent(in) :: condensible_names(:)
    real(dp), intent(in) :: condensible_P(:)
    real(dp), intent(in) :: condensible_RH(:)
    real(dp), intent(in) :: f_i(:,:)
    real(dp), optional, intent(in) :: pdensities(:,:)
    real(dp), optional, intent(in) :: radii(:,:)
    real(dp), optional, intent(in) :: T_guess
    character(:), allocatable, intent(out) :: err

    real(dp) :: T_surf

    real(dp) :: T_guess_
    type(MinpackHybrd1Vars) :: mv

    if (present(T_guess)) then
      T_guess_ = T_guess
    else
      T_guess_ = 280.0_dp
    endif

    mv = MinpackHybrd1Vars(1)
    mv%x(1) = log10(T_guess_)
    call hybrd1(fcn, mv%n, mv%x, mv%fvec, mv%tol, mv%info, mv%wa, mv%lwa)
    if (mv%info == 0 .or. mv%info > 1) then
      err = 'hybrd1 root solve failed in surface_temperature.'
      return
    elseif (mv%info < 0) then
      err = 'hybrd1 root solve failed in surface_temperature: '//err
      return
    endif

    T_surf = 10.0_dp**mv%x(1)
    call fcn(mv%n, mv%x, mv%fvec, mv%info)

  contains 
    subroutine fcn(n_, x_, fvec_, iflag_)
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(n_)
      integer, intent(inout) :: iflag_
      real(dp) :: T, ISR, OLR
      T = 10.0_dp**x_(1)
      call self%TOA_fluxes(T, condensible_names, condensible_P, condensible_RH, &
                           f_i, pdensities, radii, ISR, OLR, err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif
      fvec_(1) = ISR - OLR
    end subroutine
  end function

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

  subroutine RadiativeConvectiveClimate_initialize_stepper(self, condensible_names, condensible_P, condensible_RH, &
                                                           f_i, pdensities, radii, T_surf_guess, err)
    use clima_rc_adiabat, only: make_profile_rc
    class(RadiativeConvectiveClimate), intent(inout) :: self
    character(*), intent(in) :: condensible_names(:)
    real(dp), intent(in) :: condensible_P(:)
    real(dp), intent(in) :: condensible_RH(:)
    real(dp), intent(in) :: f_i(:,:)
    real(dp), optional, intent(in) :: pdensities(:,:)
    real(dp), optional, intent(in) :: radii(:,:)
    real(dp), optional :: T_surf_guess
    character(:), allocatable, intent(out) :: err

    real(dp) :: T_surf

    ! Initial guess
    T_surf = self%surface_temperature(condensible_names, condensible_P, condensible_RH, &
                                      f_i, pdensities, radii, T_surf_guess, err)
    if (allocated(err)) return
    
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
    elseif (np == 0) then
      if (present(radii)) then
        err = 'The model does not contain particles but "pdensities" and "radii" are arguments.'
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