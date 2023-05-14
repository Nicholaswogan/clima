module clima_adiabat
  use clima_const, only: dp, s_str_len
  use clima_types, only: Species
  use clima_radtran, only: Radtran
  use clima_eqns, only: ocean_solubility_fcn, temp_dependent_albedo_fcn
  use clima_adiabat_general, only: OceanFunction
  use iso_c_binding, only: c_ptr, c_null_ptr
  implicit none
  private

  public :: AdiabatClimate

  type :: AdiabatClimate

    ! settings and free parameters
    integer :: nz
    real(dp) :: P_top = 1.0e-2_dp !! (dynes/cm2)
    real(dp) :: T_trop = 180.0_dp !! (T)
    real(dp), allocatable :: RH(:) !! relative humidity (ng)

    !> If .true., then Tropopause temperature is non-linearly solved for such that
    !> it matches the skin temperature. The initial guess will always be self%T_trop.
    logical :: solve_for_T_trop = .false.
    !> Callback that sets the surface albedo based on the surface temperature.
    !> This can be used to parameterize the ice-albedo feedback.
    procedure(temp_dependent_albedo_fcn), nopass, pointer :: albedo_fcn => null()

    !> Function describing how gases disolve in oceans. This allows for multiple oceans, each
    !> made of different condensed volatiles.
    type(OceanFunction), allocatable :: ocean_fcns(:)
    !> A c-ptr which will be passed to each ocean solubility function.
    type(c_ptr) :: ocean_args_p = c_null_ptr
    
    ! planet properties
    real(dp) :: planet_mass !! (g)
    real(dp) :: planet_radius !! (cm)
    
    ! species in the model
    character(s_str_len), allocatable :: species_names(:) !! copy of species names
    type(Species) :: sp
    
    ! Radiative transfer
    type(Radtran) :: rad
    
    ! State of the atmosphere
    real(dp) :: P_surf !! Surface pressure (dynes/cm^2)
    real(dp) :: P_trop !! Tropopause pressure (dynes/cm^2)
    real(dp), allocatable :: P(:) !! pressure in each grid cell, dynes/cm^2 (nz)
    real(dp) :: T_surf !! Surface Temperature (K)
    real(dp), allocatable :: T(:) !! Temperature in each grid cell, K (nz) 
    real(dp), allocatable :: f_i(:,:) !! mixing ratios of species in each grid cell (nz,ng)
    real(dp), allocatable :: z(:) !! Altitude at the center of the grid cell, cm (nz)
    real(dp), allocatable :: dz(:) !! Thickness of each grid cell, cm (nz)
    real(dp), allocatable :: densities(:,:) !! densities in each grid cell, molecules/cm^3 (nz,ng)
    real(dp), allocatable :: N_atmos(:) !! reservoir of gas in atmosphere mol/cm^2 (ng)
    real(dp), allocatable :: N_surface(:) !! reservoir of gas on surface mol/cm^2 (ng)
    !> reservoir of gas dissolved in oceans in mol/cm^2 (ng, ng). There can be multiple oceans.
    !> The gases dissolved in ocean made of species 1 is given by `N_ocean(:,1)`.
    real(dp), allocatable :: N_ocean(:,:) 
    
  contains
    ! Constructs atmospheres
    procedure :: make_profile => AdiabatClimate_make_profile
    procedure :: make_column => AdiabatClimate_make_column
    procedure :: make_profile_bg_gas => AdiabatClimate_make_profile_bg_gas
    ! Constructs atmosphere and does radiative transfer
    procedure :: TOA_fluxes => AdiabatClimate_TOA_fluxes
    procedure :: TOA_fluxes_column => AdiabatClimate_TOA_fluxes_column
    procedure :: TOA_fluxes_bg_gas => AdiabatClimate_TOA_fluxes_bg_gas
    ! Non-linear solves for equilibrium climate
    procedure :: surface_temperature => AdiabatClimate_surface_temperature
    procedure :: surface_temperature_column => AdiabatClimate_surface_temperature_column
    procedure :: surface_temperature_bg_gas => AdiabatClimate_surface_temperature_bg_gas
    ! Utilities
    procedure :: set_ocean_solubility_fcn => AdiabatClimate_set_ocean_solubility_fcn
    procedure :: to_regular_grid => AdiabatClimate_to_regular_grid
    procedure :: out2atmosphere_txt => AdiabatClimate_out2atmosphere_txt
  end type
  
  interface AdiabatClimate
    module procedure :: create_AdiabatClimate
  end interface
  
contains
  
  function create_AdiabatClimate(species_f, settings_f, star_f, data_dir, err) result(c)
    use clima_types, only: ClimaSettings 
    character(*), intent(in) :: species_f !! Species yaml file
    character(*), intent(in) :: settings_f !! Settings yaml file
    character(*), intent(in) :: star_f !! Star text file
    character(*), intent(in) :: data_dir !! Directory with radiative transfer data
    character(:), allocatable, intent(out) :: err
    
    type(AdiabatClimate) :: c
    
    type(ClimaSettings) :: s
    integer :: i
    character(s_str_len) :: particle_names(0)

    ! species
    c%sp = Species(species_f, err)
    if (allocated(err)) return
    
    ! unpack species
    allocate(c%species_names(c%sp%ng))
    do i = 1,c%sp%ng
      c%species_names(i) = c%sp%g(i)%name
    enddo

    ! default relative humidty is 1
    allocate(c%RH(c%sp%ng))
    c%RH(:) = 1.0_dp
    
    ! There must be more than 1 species
    if (c%sp%ng == 1) then 
      err = 'There must be more than 1 species in "'//species_f//'"'
      return
    endif
    
    ! settings
    s = ClimaSettings(settings_f, err)
    if (allocated(err)) return
    
    ! unpack setting
    if (s%atmos_grid_is_present) then
      c%nz = s%nz
    else
      err = '"atmosphere-grid" is missing from file "'//settings_f//'"'
      return
    endif
    if (.not.s%planet_is_present) then
      err = '"planet" is missing from file "'//settings_f//'"'
      return
    endif
    c%planet_mass = s%planet_mass
    c%planet_radius = s%planet_radius
    if (.not.allocated(s%number_of_zenith_angles)) then
      err = '"number-of-zenith-angles" is missing from file "'//settings_f//'"'
      return
    endif
    if (.not.allocated(s%surface_albedo)) then
      err = '"surface-albedo" is missing from file "'//settings_f//'"'
      return
    endif
    
    ! Make radiative transfer with list of species, from species file
    ! and the optical-properties from the settings file
    c%rad = Radtran(c%species_names, particle_names, s, star_f, s%number_of_zenith_angles, s%surface_albedo, c%nz, data_dir, err)
    if (allocated(err)) return

    ! allocate ocean functions
    allocate(c%ocean_fcns(c%sp%ng))

    ! allocate work variables
    allocate(c%P(c%nz), c%T(c%nz), c%f_i(c%nz,c%sp%ng), c%z(c%nz), c%dz(c%nz))
    allocate(c%densities(c%nz,c%sp%ng))
    allocate(c%N_atmos(c%sp%ng),c%N_surface(c%sp%ng),c%N_ocean(c%sp%ng,c%sp%ng))
    
  end function
  
  !> Constructs an atmosphere using a multispecies pseudoadiabat (Eq. 1 in Graham et al. 2021, PSJ)
  !> troposphere connected to an isothermal stratosphere. Input are surface partial pressures
  !> of each gas.
  subroutine AdiabatClimate_make_profile(self, T_surf, P_i_surf, err)
    use clima_adiabat_general, only: make_profile
    use clima_const, only: k_boltz, N_avo
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: T_surf !! K
    real(dp), intent(in) :: P_i_surf(:) !! dynes/cm^2
    character(:), allocatable, intent(out) :: err
    
    real(dp), allocatable :: P_e(:), z_e(:), T_e(:), f_i_e(:,:)
    real(dp), allocatable :: density(:)
    integer :: i, j

    if (size(P_i_surf) /= self%sp%ng) then
      err = "P_i_surf has the wrong dimension"
      return
    endif
    
    allocate(P_e(2*self%nz+1),  z_e(2*self%nz+1), T_e(2*self%nz+1))
    allocate(f_i_e(2*self%nz+1,self%sp%ng))
    allocate(density(self%nz))

    call make_profile(T_surf, P_i_surf, &
                      self%sp, self%nz, self%planet_mass, &
                      self%planet_radius, self%P_top, self%T_trop, self%RH, &
                      self%ocean_fcns, self%ocean_args_p, &
                      P_e, z_e, T_e, f_i_e, self%P_trop, &
                      self%N_surface, self%N_ocean, &
                      err)
    if (allocated(err)) return

    self%T_surf = T_surf
    self%P_surf = P_e(1)
    do i = 1,self%nz
      self%P(i) = P_e(2*i)
      self%T(i) = T_e(2*i)
      self%z(i) = z_e(2*i)
      self%dz(i) = z_e(2*i+1) - z_e(2*i-1)
      do j =1,self%sp%ng
        self%f_i(i,j) = f_i_e(2*i,j)
      enddo
    enddo
    
    density = self%P/(k_boltz*self%T)
    do j =1,self%sp%ng
      self%densities(:,j) = self%f_i(:,j)*density(:)
    enddo

    do i = 1,self%sp%ng
      ! mol/cm^2 in atmosphere
      self%N_atmos(i) = sum(density*self%f_i(:,i)*self%dz)/N_avo
    enddo
    
  end subroutine

  !> Similar to `make_profile`, but instead the input is column reservoirs 
  !> of each gas (mol/cm^2). 
  subroutine AdiabatClimate_make_column(self, T_surf, N_i_surf, err)
    use clima_adiabat_general, only: make_column
    use clima_const, only: k_boltz, N_avo
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: T_surf !! K
    real(dp), intent(in) :: N_i_surf(:) !! mole/cm^2
    character(:), allocatable, intent(out) :: err
    
    real(dp), allocatable :: P_e(:), z_e(:), T_e(:), f_i_e(:,:)
    real(dp), allocatable :: density(:)
    integer :: i, j

    if (size(N_i_surf) /= self%sp%ng) then
      err = "N_i_surf has the wrong dimension"
      return
    endif
    
    allocate(P_e(2*self%nz+1),  z_e(2*self%nz+1), T_e(2*self%nz+1))
    allocate(f_i_e(2*self%nz+1,self%sp%ng))
    allocate(density(self%nz))
    
    call make_column(T_surf, N_i_surf, &
                     self%sp, self%nz, self%planet_mass, &
                     self%planet_radius, self%P_top, self%T_trop, self%RH, &
                     self%ocean_fcns, self%ocean_args_p, &
                     P_e, z_e, T_e, f_i_e, self%P_trop, &
                     self%N_surface, self%N_ocean, &
                     err)
    if (allocated(err)) return
    
    self%T_surf = T_surf
    self%P_surf = P_e(1)
    do i = 1,self%nz
      self%P(i) = P_e(2*i)
      self%T(i) = T_e(2*i)
      self%z(i) = z_e(2*i)
      self%dz(i) = z_e(2*i+1) - z_e(2*i-1)
      do j =1,self%sp%ng
        self%f_i(i,j) = f_i_e(2*i,j)
      enddo
    enddo
    
    density = self%P/(k_boltz*self%T)
    do j =1,self%sp%ng
      self%densities(:,j) = self%f_i(:,j)*density(:)
    enddo

    do i = 1,self%sp%ng
      ! mol/cm^2 in atmosphere
      self%N_atmos(i) = sum(density*self%f_i(:,i)*self%dz)/N_avo
    enddo
    
  end subroutine

  !> Similar to `make_profile`, but instead imposes a background gas and fixed
  !> surface pressure. We do a non-linear solve for the background gas pressure
  !> so that the desired total surface pressure is satisfied.
  subroutine AdiabatClimate_make_profile_bg_gas(self, T_surf, P_i_surf, P_surf, bg_gas, err)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: T_surf !! K
    real(dp), intent(in) :: P_i_surf(:) !! dynes/cm^2
    real(dp), intent(in) :: P_surf !! dynes/cm^2
    character(*), intent(in) :: bg_gas !! background gas
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: P_i_surf_copy(:)
    integer :: i, ind
    type(MinpackHybrd1Vars) :: mv
    real(dp), parameter :: scale_factors(*) = [1.0_dp, 0.1_dp]

    if (P_surf <= 0.0_dp) then
      err = 'P_surf must be greater than zero.'
      return
    endif

    ind = findloc(self%species_names, bg_gas, 1)
    if (ind == 0) then
      err = 'Gas "'//bg_gas//'" is not in the list of species'
      return
    endif

    allocate(P_i_surf_copy(size(P_i_surf)))
    P_i_surf_copy = P_i_surf

    mv = MinpackHybrd1Vars(1)
  
    do i = 1,size(scale_factors)
      mv%x(1) = log10(P_surf*scale_factors(i))
      call hybrd1(fcn, mv%n, mv%x, mv%fvec, mv%tol, mv%info, mv%wa, mv%lwa)
      if (mv%info == 1) then
        exit
      endif
    enddo
    if (mv%info == 0 .or. mv%info > 1) then
      err = 'hybrd1 root solve failed in make_profile_bg_gas.'
      return
    elseif (mv%info < 0) then
      err = 'hybrd1 root solve failed in make_profile_bg_gas: '//err
      return
    endif

    call fcn(mv%n, mv%x, mv%fvec, mv%info)

  contains
    subroutine fcn(n_, x_, fvec_, iflag_)
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(n_)
      integer, intent(inout) :: iflag_

      P_i_surf_copy(ind) = 10.0_dp**x_(1)

      call self%make_profile(T_surf, P_i_surf_copy, err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif

      fvec_(1) = self%P_surf - P_surf
    end subroutine
  end subroutine
  
  !> Calls `make_profile`, then does radiative transfer on the constructed atmosphere
  subroutine AdiabatClimate_TOA_fluxes(self, T_surf, P_i_surf, ISR, OLR, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), target, intent(in) :: P_i_surf(:) !! dynes/cm^2
    real(dp), intent(out) :: ISR !! Top-of-atmosphere incoming solar radiation (mW/m^2)
    real(dp), intent(out) :: OLR !! Top-of-atmosphere outgoing longwave radiation (mW/m^2)
    character(:), allocatable, intent(out) :: err
    
    ! make atmosphere profile
    call self%make_profile(T_surf, P_i_surf, err)
    if (allocated(err)) return

    ! set albedo
    if (associated(self%albedo_fcn)) then
      self%rad%surface_albedo = self%albedo_fcn(T_surf)
    endif
    
    ! Do radiative transfer
    call self%rad%TOA_fluxes(T_surf, self%T, self%P/1.0e6_dp, self%densities, self%dz, ISR=ISR, OLR=OLR, err=err)
    if (allocated(err)) return
    
  end subroutine

  !> Calls `make_column`, then does radiative transfer on the constructed atmosphere
  subroutine AdiabatClimate_TOA_fluxes_column(self, T_surf, N_i_surf, ISR, OLR, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), target, intent(in) :: N_i_surf(:) !! mole/cm^2
    real(dp), intent(out) :: ISR !! Top-of-atmosphere incoming solar radiation (mW/m^2)
    real(dp), intent(out) :: OLR !! Top-of-atmosphere outgoing longwave radiation (mW/m^2)
    character(:), allocatable, intent(out) :: err    
    
    ! make atmosphere profile
    call self%make_column(T_surf, N_i_surf, err)
    if (allocated(err)) return

    ! set albedo
    if (associated(self%albedo_fcn)) then
      self%rad%surface_albedo = self%albedo_fcn(T_surf)
    endif
    
    ! Do radiative transfer
    call self%rad%TOA_fluxes(T_surf, self%T, self%P/1.0e6_dp, self%densities, self%dz, ISR=ISR, OLR=OLR, err=err)
    if (allocated(err)) return

  end subroutine

  !> Calls `make_profile_bg_gas`, then does radiative transfer on the constructed atmosphere
  subroutine AdiabatClimate_TOA_fluxes_bg_gas(self, T_surf, P_i_surf, P_surf, bg_gas, ISR, OLR, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), target, intent(in) :: P_i_surf(:) !! dynes/cm^2
    real(dp), intent(in) :: P_surf !! dynes/cm^2
    character(*), intent(in) :: bg_gas !! background gas
    real(dp), intent(out) :: ISR !! Top-of-atmosphere incoming solar radiation (mW/m^2)
    real(dp), intent(out) :: OLR !! Top-of-atmosphere outgoing longwave radiation (mW/m^2)
    character(:), allocatable, intent(out) :: err
    
    ! make atmosphere profile
    call self%make_profile_bg_gas(T_surf, P_i_surf, P_surf, bg_gas, err)
    if (allocated(err)) return

    ! set albedo
    if (associated(self%albedo_fcn)) then
      self%rad%surface_albedo = self%albedo_fcn(T_surf)
    endif
    
    ! Do radiative transfer
    call self%rad%TOA_fluxes(T_surf, self%T, self%P/1.0e6_dp, self%densities, self%dz, ISR=ISR, OLR=OLR, err=err)
    if (allocated(err)) return
    
  end subroutine
  
  !> Does a non-linear solve for the surface temperature that balances incoming solar
  !> and outgoing longwave radiation. Uses `make_profile`.
  function AdiabatClimate_surface_temperature(self, P_i_surf, T_guess, err) result(T_surf)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:) !! dynes/cm^2
    real(dp), optional, intent(in) :: T_guess !! K
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: T_surf
    
    real(dp) :: T_guess_
    type(MinpackHybrd1Vars) :: mv
    
    if (present(T_guess)) then
      T_guess_ = T_guess
    else
      T_guess_ = 280.0_dp
    endif
    
    if (self%solve_for_T_trop) then
      mv = MinpackHybrd1Vars(2)
      mv%x(1) = log10(T_guess_)
      mv%x(2) = log10(self%T_trop)
    else
      mv = MinpackHybrd1Vars(1)
      mv%x(1) = log10(T_guess_)
    endif
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
      if (self%solve_for_T_trop) then
        self%T_trop = 10.0_dp**x_(2)
      endif
      call self%TOA_fluxes(T, P_i_surf, ISR, OLR, err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif
      fvec_(1) = ISR - OLR

      if (self%solve_for_T_trop) then; block
        real(dp) :: bond_albedo
        bond_albedo = self%rad%wrk_sol%fup_n(self%nz+1)/self%rad%wrk_sol%fdn_n(self%nz+1)
        fvec_(2) = self%rad%skin_temperature(bond_albedo) - self%T_trop
      endblock; endif
    end subroutine
    
  end function

  function AdiabatClimate_surface_temperature_column(self, N_i_surf, T_guess, err) result(T_surf)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: N_i_surf(:)
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
    
    if (self%solve_for_T_trop) then
      mv = MinpackHybrd1Vars(2)
      mv%x(1) = log10(T_guess_)
      mv%x(2) = log10(self%T_trop)
    else
      mv = MinpackHybrd1Vars(1)
      mv%x(1) = log10(T_guess_)
    endif
    call hybrd1(fcn, mv%n, mv%x, mv%fvec, mv%tol, mv%info, mv%wa, mv%lwa)
    if (mv%info == 0 .or. mv%info > 1) then
      err = 'hybrd1 root solve failed in surface_temperature_column.'
      return
    elseif (mv%info < 0) then
      err = 'hybrd1 root solve failed in surface_temperature_column: '//err
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
      if (self%solve_for_T_trop) then
        self%T_trop = 10.0_dp**x_(2)
      endif
      call self%TOA_fluxes_column(T, N_i_surf, ISR, OLR, err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif
      fvec_(1) = ISR - OLR

      if (self%solve_for_T_trop) then; block
        real(dp) :: bond_albedo
        bond_albedo = self%rad%wrk_sol%fup_n(self%nz+1)/self%rad%wrk_sol%fdn_n(self%nz+1)
        fvec_(2) = self%rad%skin_temperature(bond_albedo) - self%T_trop
      endblock; endif
    end subroutine
    
  end function

  !> Similar to surface_temperature. The difference is that this function imposes
  !> a background gas and fixed surface pressure.
  function AdiabatClimate_surface_temperature_bg_gas(self, P_i_surf, P_surf, bg_gas, T_guess, err) result(T_surf)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:) !! dynes/cm^2
    real(dp), intent(in) :: P_surf !! dynes/cm^2
    character(*), intent(in) :: bg_gas !! background gas
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

    if (self%solve_for_T_trop) then
      mv = MinpackHybrd1Vars(2)
      mv%x(1) = log10(T_guess_)
      mv%x(2) = log10(self%T_trop)
    else
      mv = MinpackHybrd1Vars(1)
      mv%x(1) = log10(T_guess_)
    endif
    call hybrd1(fcn, mv%n, mv%x, mv%fvec, mv%tol, mv%info, mv%wa, mv%lwa)
    if (mv%info == 0 .or. mv%info > 1) then
      err = 'hybrd1 root solve failed in surface_temperature_bg_gas.'
      return
    elseif (mv%info < 0) then
      err = 'hybrd1 root solve failed in surface_temperature_bg_gas: '//err
      return
    endif

    T_surf = 10.0**mv%x(1)
    call fcn(mv%n, mv%x, mv%fvec, mv%info)

  contains
    subroutine fcn(n_, x_, fvec_, iflag_)
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(n_)
      integer, intent(inout) :: iflag_
      real(dp) :: T, ISR, OLR
      T = 10.0_dp**x_(1)
      if (self%solve_for_T_trop) then
        self%T_trop = 10.0_dp**x_(2)
      endif
      call self%TOA_fluxes_bg_gas(T, P_i_surf, P_surf, bg_gas, ISR, OLR, err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif
      fvec_(1) = ISR - OLR

      if (self%solve_for_T_trop) then; block
        real(dp) :: bond_albedo
        bond_albedo = self%rad%wrk_sol%fup_n(self%nz+1)/self%rad%wrk_sol%fdn_n(self%nz+1)
        fvec_(2) = self%rad%skin_temperature(bond_albedo) - self%T_trop
      endblock; endif
    end subroutine
  end function

  !> Sets a function for describing how gases dissolve in a liquid ocean.
  subroutine AdiabatClimate_set_ocean_solubility_fcn(self, species, fcn, err)
    class(AdiabatClimate), intent(inout) :: self
    character(*), intent(in) :: species !! name of species that makes the ocean
    !> Function describing solubility of other gases in ocean
    procedure(ocean_solubility_fcn), pointer, intent(in) :: fcn 
    character(:), allocatable, intent(out) :: err

    integer :: ind

    ind = findloc(self%species_names, species, 1)
    if (ind == 0) then
      err = 'Gas "'//species//'" is not in the list of species'
      return
    endif

    self%ocean_fcns(ind)%fcn => fcn

  end subroutine

  !> Re-grids atmosphere so that each grid cell is equally spaced in altitude.
  subroutine AdiabatClimate_to_regular_grid(self, err)
    use futils, only: rebin, interp
    use clima_eqns, only: vertical_grid
    use clima_const, only: k_boltz
    class(AdiabatClimate), intent(inout) :: self
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: ze(:), ze_new(:)
    real(dp), allocatable :: z_new(:), dz_new(:)
    real(dp), allocatable :: densities_new(:,:)
    real(dp), allocatable :: f_i_new(:,:)
    real(dp), allocatable :: T_new(:)
    real(dp), allocatable :: density_new(:)
    real(dp), allocatable :: P_new(:)

    integer :: i, j, ierr

    allocate(z_new(self%nz), dz_new(self%nz))

    ! compute the new grid
    call vertical_grid(0.0_dp, self%z(self%nz)+0.5_dp*self%dz(self%nz), &
                       self%nz, z_new, dz_new)

    ! rebin of the densities
    allocate(ze(self%nz+1), ze_new(self%nz+1))
    ze_new(1) = z_new(1) - 0.5_dp*dz_new(1)
    do i = 1,self%nz
      ze_new(i+1) = z_new(i) + 0.5_dp*dz_new(i)
    enddo
    ze(1) = self%z(1) - 0.5_dp*self%dz(1)
    do i = 1,self%nz
      ze(i+1) = self%z(i) + 0.5_dp*self%dz(i)
    enddo

    allocate(densities_new(self%nz, self%sp%ng))
    do i = 1,self%sp%ng
      call rebin(ze, self%densities(:,i), ze_new, densities_new(:,i), ierr)
      if (ierr /= 0) then
        err = 'subroutine `rebin` returned an error'
        return
      endif
    enddo

    allocate(T_new(self%nz))
    call interp(self%nz, self%nz, z_new, self%z, self%T, T_new, ierr)
    if (ierr /= 0) then
      err = 'Subroutine `interp` returned an error.'
      return
    endif

    allocate(density_new(self%nz))
    allocate(f_i_new(self%nz,self%sp%ng))
    do i = 1,self%nz
      density_new(i) = sum(densities_new(i,:))
      do j = 1,self%sp%ng
        f_i_new(i,j) = densities_new(i,j)/density_new(i)
      enddo
    enddo
    allocate(P_new(self%nz))
    P_new = density_new(:)*k_boltz*T_new(:)

    self%P(:) = P_new(:)
    self%T(:) = T_new(:)
    self%f_i(:,:) = f_i_new(:,:)
    self%z(:) = z_new(:)
    self%dz(:) = dz_new(:)
    self%densities(:,:) = densities_new(:,:)

  end subroutine

  subroutine AdiabatClimate_out2atmosphere_txt(self, filename, eddy, overwrite, clip, err)
    class(AdiabatClimate), target, intent(inout) :: self
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: eddy(:)
    logical, intent(in) :: overwrite, clip
    character(:), allocatable, intent(out) :: err
    
    character(len=100) :: tmp
    integer :: io, i, j

    call self%to_regular_grid(err)
    if (allocated(err)) return

    if (size(eddy,1) /= self%nz) then
      err = '"eddy" has the wrong size'
      return
    endif
    
    if (overwrite) then
      open(2, file=filename, form='formatted', status='replace', iostat=io)
      if (io /= 0) then
        err = "Unable to overwrite file "//trim(filename)
        return
      endif
    else
      open(2, file=filename, form='formatted', status='new', iostat=io)
      if (io /= 0) then
        err = "Unable to create file "//trim(filename)//" because it already exists"
        return
      endif
    endif
    
    tmp = 'alt'
    write(unit=2,fmt="(3x,a27)",advance='no') tmp
    tmp = 'press'
    write(unit=2,fmt="(a27)",advance='no') tmp
    tmp = 'den'
    write(unit=2,fmt="(a27)",advance='no') tmp
    tmp = 'temp'
    write(unit=2,fmt="(a27)",advance='no') tmp
    tmp = 'eddy'
    write(unit=2,fmt="(a27)",advance='no') tmp
    do j = 1,self%sp%ng
      tmp = self%species_names(j)
      write(unit=2,fmt="(a27)",advance='no') tmp
    enddo
    
    do i = 1,self%nz
      write(2,*)
      write(unit=2,fmt="(es27.17e3)",advance='no') self%z(i)/1.e5_dp
      write(unit=2,fmt="(es27.17e3)",advance='no') self%P(i)/1.e6_dp
      write(unit=2,fmt="(es27.17e3)",advance='no') sum(self%densities(i,:))
      write(unit=2,fmt="(es27.17e3)",advance='no') self%T(i)
      write(unit=2,fmt="(es27.17e3)",advance='no') eddy(i)
      do j = 1,self%sp%ng
        if (clip) then
          write(unit=2,fmt="(es27.17e3)",advance='no') max(self%f_i(i,j),1.0e-40_dp)
        else
          write(unit=2,fmt="(es27.17e3)",advance='no') self%f_i(i,j)
        endif
      enddo
    enddo
    
    close(2)
    
  end subroutine
  
end module