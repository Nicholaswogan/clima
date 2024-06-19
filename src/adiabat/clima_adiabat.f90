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
    real(dp) :: P_top = 1.0_dp !! (dynes/cm2)
    real(dp) :: T_trop = 180.0_dp !! (T)
    real(dp), allocatable :: RH(:) !! relative humidity (ng)
    !> If .true., then any function that calls `make_column` will
    !> use the initial guess in `self%make_column_P_guess`
    logical :: use_make_column_P_guess = .true.
    !> Initial guess for surface pressure of all gases for `make_column`
    !> routine (ng).
    real(dp), allocatable :: make_column_P_guess(:)

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

    ! Heat re-distribution terms for Equation (10) in Koll (2020), ApJ when
    ! considering tidally locked planets. Default values are from the paper.
    !> If true, then will attempt to compute the climate corresponding to the
    !> observed dayside temperature of a tidally locked planet.
    logical :: tidally_locked_dayside = .false.
    real(dp) :: L !! = planet radius. Circulation’s horizontal scale (cm)
    real(dp) :: chi = 0.2_dp !! Heat engine efficiency term (no units)
    real(dp) :: n_LW = 2.0_dp !! = 1 or 2 (no units)
    real(dp) :: Cd = 1.9e-3_dp !! Drag coefficient (no units)

    !> Heat flow from the surface into the atmosphere (mW/m^2)
    real(dp) :: surface_heat_flow = 0.0_dp
    
    ! planet properties
    real(dp) :: planet_mass !! (g)
    real(dp) :: planet_radius !! (cm)
    
    ! species in the model
    character(s_str_len), allocatable :: species_names(:) !! copy of species names
    type(Species) :: sp
    
    ! Radiative transfer
    type(Radtran) :: rad
    logical :: double_radiative_grid
    integer :: nz_r
    ! Work space for radiative transfer
    real(dp), allocatable :: T_r(:)
    real(dp), allocatable :: P_r(:)
    real(dp), allocatable :: densities_r(:,:)
    real(dp), allocatable :: dz_r(:)

    ! Information about convection
    integer :: n_convecting_zones !! number of convecting zones
    !> Describes the lower index of each convecting zone. Note index 1
    !> is the ground layer.
    integer, private, allocatable :: ind_conv_lower(:) 
    !> Describes the upper index of each convecting zone.
    integer, private, allocatable :: ind_conv_upper(:)
    integer, private, allocatable :: ind_conv_lower_x(:)
    !> Indicates if a layer is super-saturated
    logical, private, allocatable :: super_saturated(:)
    !> Another representation of where convection is occuring. If True,
    !> then the layer below is convecting with the current layer. Index 1 determines
    !> if the first atomspheric layer is convecting with the ground.
    logical, allocatable :: convecting_with_below(:)
    integer, allocatable :: inds_Tx(:)
    real(dp), allocatable :: lapse_rate(:) !! The true lapse rate (dlnT/dlnP)
    real(dp), allocatable :: lapse_rate_intended(:) !! The computed lapse rate (dlnT/dlnP)
    !> The size of the newton step.
    real(dp) :: convective_newton_step_size = 1.0e-1_dp

    ! tolerances
    !> Relative tolerance of integration
    real(dp) :: rtol = 1.0e-12_dp
    !> Absolute tolerance of integration
    real(dp) :: atol = 1.0e-12_dp
    !> Tolerance for nonlinear solve in make_column
    real(dp) :: tol_make_column = 1.0e-8_dp
    !> Perturbation for the jacobian
    real(dp) :: epsj = 1.0e-2_dp
    !> xtol for RC equilibrium
    real(dp) :: xtol_rc = 1.0e-5_dp
    !> Max number of iterations in the RCE routine
    integer :: max_rc_iters = 10
    !> Max number of iterations for which convective layers can
    !> be converged to radiative layers in the RCE routine
    integer :: max_rc_iters_convection = 5
    !> A term that weights the importance of maintaining radiative
    !> equilibrium to convection.
    real(dp) :: radiation_norm_term = 1.0e-3_dp
    logical :: verbose = .true. !! verbosity
    
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
    ! Routines for full radiative convective equilibrium
    procedure, private :: make_profile_rc => AdiabatClimate_make_profile_rc
    procedure :: RCE => AdiabatClimate_RCE
    ! Utilities
    procedure :: set_ocean_solubility_fcn => AdiabatClimate_set_ocean_solubility_fcn
    procedure :: to_regular_grid => AdiabatClimate_to_regular_grid
    procedure :: out2atmosphere_txt => AdiabatClimate_out2atmosphere_txt
    ! For tidally locked planets
    procedure :: heat_redistribution_parameters => AdiabatClimate_heat_redistribution_parameters

    ! Private methods
    procedure, private :: copy_atm_to_radiative_grid => AdiabatClimate_copy_atm_to_radiative_grid
  end type
  
  interface AdiabatClimate
    module procedure :: create_AdiabatClimate
  end interface

  interface
    module subroutine AdiabatClimate_make_profile_rc(self, P_i_surf, T_in, err)
      class(AdiabatClimate), intent(inout) :: self
      real(dp), intent(in) :: P_i_surf(:) !! dynes/cm^2
      real(dp), intent(in) :: T_in(:)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Compute full radiative-convective equilibrium.
    module function AdiabatClimate_RCE(self, P_i_surf, T_surf_guess, T_guess, convecting_with_below, err) result(converged)
      class(AdiabatClimate), intent(inout) :: self
      !> Array of surface pressures of each species (dynes/cm^2)
      real(dp), intent(in) :: P_i_surf(:)
      !> A guess for the surface temperature (K)
      real(dp), intent(in) :: T_surf_guess
      !> A guess for the temperature in each atmospheric layer (K)
      real(dp), intent(in) :: T_guess(:)
      !> An array describing a guess for the radiative vs. convective 
      !> regions of the atmosphere
      logical, optional, intent(in) :: convecting_with_below(:)
      character(:), allocatable, intent(out) :: err
      logical :: converged !! Whether the routine converged or not.
    end function
  end interface
  
contains
  
  function create_AdiabatClimate(species_f, settings_f, star_f, data_dir, double_radiative_grid, err) result(c)
    use clima_types, only: ClimaSettings 
    character(*), intent(in) :: species_f !! Species yaml file
    character(*), intent(in) :: settings_f !! Settings yaml file
    character(*), intent(in) :: star_f !! Star text file
    character(*), intent(in) :: data_dir !! Directory with radiative transfer data
    logical, optional, intent(in) :: double_radiative_grid
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

    allocate(c%make_column_P_guess(c%sp%ng))
    c%make_column_P_guess(:) = 1.0_dp
    
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

    if (present(double_radiative_grid)) then
      c%double_radiative_grid = double_radiative_grid
    else
      c%double_radiative_grid = .true.
    endif
    if (c%double_radiative_grid) then
      c%nz_r = c%nz*2
    else
      c%nz_r = c%nz
    endif
    allocate(c%T_r(c%nz_r),c%P_r(c%nz_r),c%densities_r(c%nz_r,c%sp%ng),c%dz_r(c%nz_r))
    ! Make radiative transfer with list of species, from species file
    ! and the optical-properties from the settings file
    c%rad = Radtran(c%species_names, particle_names, s, star_f, s%number_of_zenith_angles, s%surface_albedo, c%nz_r, data_dir, err)
    if (allocated(err)) return

    ! Convection
    allocate(c%super_saturated(c%nz))
    allocate(c%convecting_with_below(c%nz))
    allocate(c%lapse_rate(c%nz),c%lapse_rate_intended(c%nz))

    ! allocate ocean functions
    allocate(c%ocean_fcns(c%sp%ng))

    ! allocate work variables
    allocate(c%P(c%nz), c%T(c%nz), c%f_i(c%nz,c%sp%ng), c%z(c%nz), c%dz(c%nz))
    allocate(c%densities(c%nz,c%sp%ng))
    allocate(c%N_atmos(c%sp%ng),c%N_surface(c%sp%ng),c%N_ocean(c%sp%ng,c%sp%ng))

    ! Heat redistribution parameter
    c%L = c%planet_radius

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
                      self%rtol, self%atol, &
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

    do i = 1,self%nz
      if (self%P(i) > self%P_trop) then
        self%convecting_with_below(i) = .true.
      else
        self%convecting_with_below(i) = .false.
      endif
    enddo

    self%lapse_rate(1) = (log(self%T(1)) - log(self%T_surf))/(log(self%P(1)) - log(self%P_surf))
    do i = 2,self%nz
      self%lapse_rate(i) = (log(self%T(i)) - log(self%T(i-1)))/(log(self%P(i)) - log(self%P(i-1)))
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
                     self%rtol, self%atol, self%tol_make_column, &
                     self%ocean_fcns, self%ocean_args_p, &
                     self%use_make_column_P_guess, self%make_column_P_guess, &
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

    do i = 1,self%nz
      if (self%P(i) > self%P_trop) then
        self%convecting_with_below(i) = .true.
      else
        self%convecting_with_below(i) = .false.
      endif
    enddo

    self%lapse_rate(1) = (log(self%T(1)) - log(self%T_surf))/(log(self%P(1)) - log(self%P_surf))
    do i = 2,self%nz
      self%lapse_rate(i) = (log(self%T(i)) - log(self%T(i-1)))/(log(self%P(i)) - log(self%P(i-1)))
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

  !> Copies the atmosphere to the data used for radiative transfer
  subroutine AdiabatClimate_copy_atm_to_radiative_grid(self)
    class(AdiabatClimate), intent(inout) :: self
    integer :: i

    if (self%double_radiative_grid) then
      do i = 1,self%nz
        self%T_r(2*(i-1)+1) = self%T(i)
        self%T_r(2*(i-1)+2) = self%T(i)

        self%P_r(2*(i-1)+1) = self%P(i)
        self%P_r(2*(i-1)+2) = self%P(i)

        self%densities_r(2*(i-1)+1,:) = self%densities(i,:)
        self%densities_r(2*(i-1)+2,:) = self%densities(i,:)

        self%dz_r(2*(i-1)+1) = 0.5_dp*self%dz(i)
        self%dz_r(2*(i-1)+2) = 0.5_dp*self%dz(i)
      enddo
    else
      self%T_r = self%T
      self%P_r = self%P
      self%densities_r = self%densities
      self%dz_r = self%dz
    endif

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
    call self%copy_atm_to_radiative_grid()
    call self%rad%TOA_fluxes(T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, ISR=ISR, OLR=OLR, err=err)
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
    call self%copy_atm_to_radiative_grid()
    call self%rad%TOA_fluxes(T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, ISR=ISR, OLR=OLR, err=err)
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
    call self%copy_atm_to_radiative_grid()
    call self%rad%TOA_fluxes(T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, ISR=ISR, OLR=OLR, err=err)
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
      real(dp) :: T, ISR, OLR, rad_enhancement
      T = 10.0_dp**x_(1)
      if (self%solve_for_T_trop) then
        self%T_trop = 10.0_dp**x_(2)
      endif
      call self%TOA_fluxes(T, P_i_surf, ISR, OLR, err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif
      rad_enhancement = 1.0_dp
      if (self%tidally_locked_dayside) then; block
        real(dp) :: tau_LW, k_term, f_term
        call self%heat_redistribution_parameters(tau_LW, k_term, f_term, err)
        if (allocated(err)) then
          iflag_ = -1
          return
        endif
        ! Increase the stellar flux, because we are computing the climate of
        ! observed dayside.
        rad_enhancement = 4.0_dp*f_term
        call self%rad%apply_radiation_enhancement(rad_enhancement)
      endblock; endif
      fvec_(1) = ISR*rad_enhancement - OLR + self%surface_heat_flow

      if (self%solve_for_T_trop) then; block
        use clima_eqns, only: skin_temperature
        real(dp) :: bond_albedo, stellar_radiation
        integer :: i
        bond_albedo = self%rad%wrk_sol%fup_n(self%nz_r+1)/self%rad%wrk_sol%fdn_n(self%nz_r+1)
        stellar_radiation = 0.0_dp
        do i = 1,self%rad%sol%nw
          stellar_radiation = stellar_radiation + self%rad%photons_sol(i)*(self%rad%sol%freq(i) - self%rad%sol%freq(i+1))
        enddo
        stellar_radiation = stellar_radiation/1.0e3_dp
        fvec_(2) = skin_temperature(stellar_radiation*rad_enhancement, bond_albedo) - self%T_trop
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
      real(dp) :: T, ISR, OLR, rad_enhancement
      T = 10.0_dp**x_(1)
      if (self%solve_for_T_trop) then
        self%T_trop = 10.0_dp**x_(2)
      endif
      call self%TOA_fluxes_column(T, N_i_surf, ISR, OLR, err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif
      rad_enhancement = 1.0_dp
      if (self%tidally_locked_dayside) then; block
        real(dp) :: tau_LW, k_term, f_term
        call self%heat_redistribution_parameters(tau_LW, k_term, f_term, err)
        if (allocated(err)) then
          iflag_ = -1
          return
        endif
        ! Increase the stellar flux, because we are computing the climate of
        ! observed dayside.
        rad_enhancement = 4.0_dp*f_term
        call self%rad%apply_radiation_enhancement(rad_enhancement)
      endblock; endif
      fvec_(1) = ISR*rad_enhancement - OLR + self%surface_heat_flow

      if (self%solve_for_T_trop) then; block
        use clima_eqns, only: skin_temperature
        real(dp) :: bond_albedo, stellar_radiation
        integer :: i
        bond_albedo = self%rad%wrk_sol%fup_n(self%nz_r+1)/self%rad%wrk_sol%fdn_n(self%nz_r+1)
        stellar_radiation = 0.0_dp
        do i = 1,self%rad%sol%nw
          stellar_radiation = stellar_radiation + self%rad%photons_sol(i)*(self%rad%sol%freq(i) - self%rad%sol%freq(i+1))
        enddo
        stellar_radiation = stellar_radiation/1.0e3_dp
        fvec_(2) = skin_temperature(stellar_radiation*rad_enhancement, bond_albedo) - self%T_trop
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
      real(dp) :: T, ISR, OLR, rad_enhancement
      T = 10.0_dp**x_(1)
      if (self%solve_for_T_trop) then
        self%T_trop = 10.0_dp**x_(2)
      endif
      call self%TOA_fluxes_bg_gas(T, P_i_surf, P_surf, bg_gas, ISR, OLR, err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif
      rad_enhancement = 1.0_dp
      if (self%tidally_locked_dayside) then; block
        real(dp) :: tau_LW, k_term, f_term
        call self%heat_redistribution_parameters(tau_LW, k_term, f_term, err)
        if (allocated(err)) then
          iflag_ = -1
          return
        endif
        ! Increase the stellar flux, because we are computing the climate of
        ! observed dayside.
        rad_enhancement = 4.0_dp*f_term
        call self%rad%apply_radiation_enhancement(rad_enhancement)
      endblock; endif
      fvec_(1) = ISR*rad_enhancement - OLR + self%surface_heat_flow

      if (self%solve_for_T_trop) then; block
        use clima_eqns, only: skin_temperature
        real(dp) :: bond_albedo, stellar_radiation
        integer :: i
        bond_albedo = self%rad%wrk_sol%fup_n(self%nz_r+1)/self%rad%wrk_sol%fdn_n(self%nz_r+1)
        stellar_radiation = 0.0_dp
        do i = 1,self%rad%sol%nw
          stellar_radiation = stellar_radiation + self%rad%photons_sol(i)*(self%rad%sol%freq(i) - self%rad%sol%freq(i+1))
        enddo
        stellar_radiation = stellar_radiation/1.0e3_dp
        fvec_(2) = skin_temperature(stellar_radiation*rad_enhancement, bond_albedo) - self%T_trop
      endblock; endif
    end subroutine
  end function

  !> Sets a function for describing how gases dissolve in a liquid ocean.
  subroutine AdiabatClimate_set_ocean_solubility_fcn(self, sp, fcn, err)
    class(AdiabatClimate), intent(inout) :: self
    character(*), intent(in) :: sp !! name of species that makes the ocean
    !> Function describing solubility of other gases in ocean
    procedure(ocean_solubility_fcn), pointer, intent(in) :: fcn 
    character(:), allocatable, intent(out) :: err

    integer :: ind

    ind = findloc(self%species_names, sp, 1)
    if (ind == 0) then
      err = 'Gas "'//sp//'" is not in the list of species'
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
    use futils, only: FileCloser
    class(AdiabatClimate), target, intent(inout) :: self
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: eddy(:)
    logical, intent(in) :: overwrite, clip
    character(:), allocatable, intent(out) :: err
    type(FileCloser) :: file
    
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
      file%unit = 2
      if (io /= 0) then
        err = "Unable to overwrite file "//trim(filename)
        return
      endif
    else
      open(2, file=filename, form='formatted', status='new', iostat=io)
      file%unit = 2
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
    
  end subroutine

  !> For considering a tidally locked planet. This function computes key parameters for
  !> Equation (10) in Koll (2022), ApJ. The function must be called after calling a function
  !> `self%TOA_fluxes`, because it uses atmosphere properties, and radiative properties.
  subroutine AdiabatClimate_heat_redistribution_parameters(self, tau_LW, k_term, f_term, err)
    use clima_const, only: c_light
    use clima_eqns, only: k_term_heat_redistribution, f_heat_redistribution, &
                          gravity, heat_capacity_eval, planck_fcn
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(out) :: tau_LW !! Long wavelength optical depth
    real(dp), intent(out) :: k_term !! k term in Equation (10)
    real(dp), intent(out) :: f_term !! The heat redistribution parameter, f in Equation (10)
    character(:), allocatable, intent(out) :: err

    integer :: i

    real(dp) :: bond_albedo, Teq
    real(dp) :: grav
    real(dp) :: mubar
    logical :: found
    real(dp) ::  cp_tmp, cp

    ! equilibrium temperature
    bond_albedo = self%rad%wrk_sol%fup_n(self%nz_r+1)/self%rad%wrk_sol%fdn_n(self%nz_r+1)
    Teq = self%rad%equilibrium_temperature(bond_albedo)
    
    ! gravity
    grav = gravity(self%planet_radius, self%planet_mass, 0.0_dp)

    ! mean molecular weight at surface
    mubar = 0.0_dp
    do i = 1,self%sp%ng
      mubar = mubar + self%f_i(1,i)*self%sp%g(i)%mass
    enddo

    ! heat capacity at surface
    cp = tiny(0.0_dp)
    do i = 1,self%sp%ng
      call heat_capacity_eval(self%sp%g(i)%thermo, self%T_surf, found, cp_tmp) ! J/(mol*K)
      if (.not. found) then
        err = "Failed to compute heat capacity"
        return
      endif
      ! J/(mol*K)
      cp = cp + cp_tmp*self%f_i(1,i) ! J/(mol*K)
    enddo
    ! J/(mol*K) * (mol/kg) = J/(kg*K)
    cp = cp*(1.0_dp/(mubar*1.0e-3_dp))
    ! J/(kg*K) * (erg/J) * (kg/g) = erg/(g*K)
    cp = cp*1.0e4_dp

    ! tau_LW. Computed with Equation (13) in Koll (2020), ApJ
    block
      real(dp) :: numerator, denominator, dlambda, tau_lambda, bplank, avg_freq, avg_lam
      numerator = 0.0_dp
      denominator = 0.0_dp
      do i = 1,self%rad%ir%nw
        dlambda = self%rad%ir%wavl(i+1) - self%rad%ir%wavl(i) ! Width of a wavelength bin (nm)
        tau_lambda = sum(self%rad%wrk_ir%tau_band(:,i)) ! IR optical depth

        avg_freq = 0.5_dp*(self%rad%ir%freq(i) + self%rad%ir%freq(i+1))
        avg_lam = (c_light*1.0e9_dp/avg_freq) ! (m/s) * (nm/m) * (s) = nm

        bplank = planck_fcn(avg_freq, self%T_surf) ! (mW sr^−1 m^−2 Hz^-1)
        bplank = bplank * (avg_freq/avg_lam) ! (mW sr^−1 m^−2 Hz^-1) * (Hz/nm) = (mW sr^−1 m^−2 nm^-1)

        ! integrate the numerator
        numerator = numerator + exp(-tau_lambda)*bplank*dlambda
        ! integrate the denominator
        denominator = denominator + bplank*dlambda
      enddo
      tau_LW = - log(numerator/denominator)
    endblock

    k_term = k_term_heat_redistribution(self%L, grav, self%chi, mubar, cp, self%n_LW, self%Cd)
    f_term = f_heat_redistribution(tau_LW, self%P_surf, Teq, k_term)

  end subroutine

end module