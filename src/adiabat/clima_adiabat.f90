module clima_adiabat
  use clima_const, only: dp, s_str_len
  use clima_types, only: Species
  use clima_radtran, only: Radtran
  implicit none
  private

  public :: WaterAdiabatClimate

  type :: WaterAdiabatClimate

    ! settings and free parameters
    integer :: nz
    real(dp) :: P_top = 1.0e-2_dp ! (dynes/cm2)
    real(dp) :: T_trop = 180.0_dp ! (T)
    real(dp) :: RH = 1.0_dp ! relative humidity
    
    ! planet properties
    real(dp) :: planet_mass ! (g)
    real(dp) :: planet_radius ! (cm)
    
    ! species in the model
    character(s_str_len), allocatable :: species_names(:) ! copy of species names
    integer :: LH2O
    type(Species) :: sp
    
    ! Radiative transfer
    type(Radtran) :: rad
    
    ! work variables
    real(dp), allocatable :: P(:), T(:), f_i(:,:), z(:), dz(:)
    real(dp), allocatable :: densities(:,:)
    
  contains
    procedure :: make_profile => WaterAdiabatClimate_make_profile
    procedure :: make_column => WaterAdiabatClimate_make_column
    procedure :: TOA_fluxes => WaterAdiabatClimate_TOA_fluxes
    procedure :: TOA_fluxes_column => WaterAdiabatClimate_TOA_fluxes_column
    procedure :: surface_temperature => WaterAdiabatClimate_surface_temperature
    procedure :: surface_temperature_column => WaterAdiabatClimate_surface_temperature_column
  end type
  
  interface WaterAdiabatClimate
    module procedure :: create_WaterAdiabatClimate
  end interface
  
contains
  
  function create_WaterAdiabatClimate(datadir, species_f, settings_f, star_f, err) result(c)
    use clima_types, only: ClimaSettings 
    character(*), intent(in) :: datadir
    character(*), intent(in) :: species_f
    character(*), intent(in) :: settings_f
    character(*), intent(in) :: star_f
    character(:), allocatable, intent(out) :: err
    
    type(WaterAdiabatClimate) :: c
    
    type(ClimaSettings) :: s
    integer :: i, ind
    character(s_str_len) :: particle_names(0)

    ! species
    c%sp = Species(species_f, err)
    if (allocated(err)) return
    
    ! unpack species
    allocate(c%species_names(c%sp%ng))
    do i = 1,c%sp%ng
      c%species_names(i) = c%sp%g(i)%name
    enddo
    
    ! H2O must be a species
    ind = findloc(c%species_names, 'H2O', 1)
    if (ind == 0) then
      err = '"H2O" must be a species in "'//species_f//'"'
      return
    else
      c%LH2O = ind
    endif
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
    if (s%planet_is_present) then
      c%planet_mass = s%planet_mass
      c%planet_radius = s%planet_radius
      if (.not.allocated(s%solar_zenith)) then
        err = '"solar-zenith-angle" is missing from file "'//settings_f//'"'
        return
      endif
      if (.not.allocated(s%surface_albedo)) then
        err = '"surface-albedo" is missing from file "'//settings_f//'"'
        return
      endif
    else
      err = '"planet" is missing from file "'//settings_f//'"'
      return
    endif
    
    ! make IR radiative transfer with list of species, from species file
    ! and the optical-properties from the settings file
    c%rad = Radtran(datadir, c%species_names, particle_names, s, star_f, s%solar_zenith, s%surface_albedo, c%nz, err)
    if (allocated(err)) return

    ! allocate work variables
    allocate(c%P(c%nz), c%T(c%nz), c%f_i(c%nz,c%sp%ng), c%z(c%nz), c%dz(c%nz))
    allocate(c%densities(c%nz,c%sp%ng))
    
  end function
  
  subroutine WaterAdiabatClimate_make_profile(self, T_surf, P_i_surf, err)
    use clima_adiabat_water, only: make_profile_water
    use clima_const, only: k_boltz
    class(WaterAdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: T_surf !! K
    real(dp), intent(in) :: P_i_surf(:)
    character(:), allocatable, intent(out) :: err
    
    real(dp), allocatable :: P_e(:), z_e(:), T_e(:), f_i_e(:,:)
    real(dp), allocatable :: density(:)
    integer :: i, j

    if (size(P_i_surf) /= self%sp%ng) then
      err = "P_i_surf has the wrong dimension"
      return
    endif
    
    allocate(P_e(self%nz+1),  z_e(self%nz+1), T_e(self%nz+1))
    allocate(f_i_e(self%nz+1,self%sp%ng))
    allocate(density(self%nz))
    
    call make_profile_water(T_surf, P_i_surf, &
                            self%sp, self%nz, self%LH2O, self%planet_mass, &
                            self%planet_radius, self%P_top, self%T_trop, self%RH, &
                            P_e, z_e, T_e, f_i_e, &
                            err)
    if (allocated(err)) return

    
    do i = 1,self%nz
      self%P(i) = P_e(i)
      self%T(i) = T_e(i)
      self%z(i) = 0.5_dp*(z_e(i)+z_e(i+1))
      self%dz(i) = z_e(i+1) - z_e(i)
      
      do j =1,self%sp%ng
        self%f_i(i,j) = f_i_e(i,j)
      enddo
    enddo
    
    density = self%P/(k_boltz*self%T)
    do j =1,self%sp%ng
      self%densities(:,j) = self%f_i(:,j)*density(:)
    enddo
    
  end subroutine

  subroutine WaterAdiabatClimate_make_column(self, T_surf, N_i_surf, err)
    use clima_adiabat_water, only: make_column_water
    use clima_const, only: k_boltz
    class(WaterAdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: T_surf !! K
    real(dp), intent(in) :: N_i_surf(:)
    character(:), allocatable, intent(out) :: err
    
    real(dp), allocatable :: P_e(:), z_e(:), T_e(:), f_i_e(:,:)
    real(dp), allocatable :: density(:)
    integer :: i, j

    if (size(N_i_surf) /= self%sp%ng) then
      err = "N_i_surf has the wrong dimension"
      return
    endif
    
    allocate(P_e(self%nz+1),  z_e(self%nz+1), T_e(self%nz+1))
    allocate(f_i_e(self%nz+1,self%sp%ng))
    allocate(density(self%nz))
    
    call make_column_water(T_surf, N_i_surf, &
                           self%sp, self%nz, self%LH2O, self%planet_mass, &
                           self%planet_radius, self%P_top, self%T_trop, self%RH, &
                           P_e, z_e, T_e, f_i_e, &
                           err)
    if (allocated(err)) return
    
    do i = 1,self%nz
      self%P(i) = P_e(i)
      self%T(i) = T_e(i)
      self%z(i) = 0.5_dp*(z_e(i)+z_e(i+1))
      self%dz(i) = z_e(i+1) - z_e(i)
      
      do j =1,self%sp%ng
        self%f_i(i,j) = f_i_e(i,j)
      enddo
    enddo
    
    density = self%P/(k_boltz*self%T)
    do j =1,self%sp%ng
      self%densities(:,j) = self%f_i(:,j)*density(:)
    enddo
    
  end subroutine
  
  subroutine WaterAdiabatClimate_TOA_fluxes(self, T_surf, P_i_surf, ISR, OLR, err)
    class(WaterAdiabatClimate), intent(inout) :: self
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), target, intent(in) :: P_i_surf(:)
    real(dp), intent(out) :: ISR, OLR
    character(:), allocatable, intent(out) :: err

    if (size(P_i_surf) /= self%sp%ng) then
      err = "P_i_surf has the wrong dimension"
      return
    endif
    
    ! make atmosphere profile
    call self%make_profile(T_surf, P_i_surf, err)
    if (allocated(err)) return
    
    ! Do radiative transfer
    ! MUST CONVERT P TO BARS
    call self%rad%TOA_fluxes(T_surf, self%T, self%P/1.0e6_dp, self%densities, self%dz, ISR=ISR, OLR=OLR, err=err)
    if (allocated(err)) return
    
  end subroutine

  subroutine WaterAdiabatClimate_TOA_fluxes_column(self, T_surf, N_i_surf, ISR, OLR, err)
    class(WaterAdiabatClimate), intent(inout) :: self
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), target, intent(in) :: N_i_surf(:)
    character(:), allocatable, intent(out) :: err    
    real(dp), intent(out) :: ISR, OLR

    if (size(N_i_surf) /= self%sp%ng) then
      err = "N_i_surf has the wrong dimension"
      return
    endif
    
    ! make atmosphere profile
    call self%make_column(T_surf, N_i_surf, err)
    if (allocated(err)) return
    
    ! Do radiative transfer
    ! MUST CONVERT P TO BARS
    call self%rad%TOA_fluxes(T_surf, self%T, self%P/1.0e6_dp, self%densities, self%dz, ISR=ISR, OLR=OLR, err=err)
    if (allocated(err)) return

  end subroutine
  
  function WaterAdiabatClimate_surface_temperature(self, P_i_surf, T_guess, err) result(T_surf)
    use minpack_module, only: hybrd1
    class(WaterAdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), optional, intent(in) :: T_guess
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: T_surf
    
    real(dp) :: T_guess_
    
    integer, parameter :: n = 1
    real(dp) :: x(1)
    real(dp) :: fvec(1)
    real(dp), parameter :: tol = 1.0e-8_dp
    integer :: info
    integer, parameter :: lwa = (n*(3*n+13))/2 + 1
    real(dp) :: wa(lwa)
    
    if (present(T_guess)) then
      T_guess_ = T_guess
    else
      T_guess_ = 280.0_dp
    endif

    if (size(P_i_surf) /= self%sp%ng) then
      err = "P_i_surf has the wrong dimension"
      return
    endif
    
    x(1) = log10(T_guess_)
    call hybrd1(fcn, n, x, fvec, tol, info, wa, lwa)
    if (info == 0 .or. info > 4) then
      err = 'hybrd1 root solve failed'
      return
    elseif (info < 0) then
      ! err already set
      err = 'hybrd1 root solve failed: '//err
      return
    endif
    
    T_surf = 10.0_dp**x(1)
    call fcn(n, x, fvec, info)
    
  contains
    subroutine fcn(n_, x_, fvec_, iflag_)
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(n_)
      integer, intent(inout) :: iflag_
      real(dp) :: T, ISR, OLR
      T = 10.0_dp**x_(1)
      call self%TOA_fluxes(T, P_i_surf, ISR, OLR, err)
      if (allocated(err)) then
        iflag_ = -1
      endif
      fvec_(1) = ISR - OLR
    end subroutine
    
  end function

  function WaterAdiabatClimate_surface_temperature_column(self, N_i_surf, T_guess, err) result(T_surf)
    use minpack_module, only: hybrd1
    class(WaterAdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: N_i_surf(:)
    real(dp), optional, intent(in) :: T_guess
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: T_surf
    
    real(dp) :: T_guess_
    
    integer, parameter :: n = 1
    real(dp) :: x(1)
    real(dp) :: fvec(1)
    real(dp), parameter :: tol = 1.0e-8_dp
    integer :: info
    integer, parameter :: lwa = (n*(3*n+13))/2 + 1
    real(dp) :: wa(lwa)

    if (size(N_i_surf) /= self%sp%ng) then
      err = "N_i_surf has the wrong dimension"
      return
    endif
    
    if (present(T_guess)) then
      T_guess_ = T_guess
    else
      T_guess_ = 280.0_dp
    endif
    
    x(1) = log10(T_guess_)
    call hybrd1(fcn, n, x, fvec, tol, info, wa, lwa)
    if (info == 0 .or. info > 4) then
      err = 'hybrd1 root solve failed'
      return
    elseif (info < 0) then
      ! err already set
      err = 'hybrd1 root solve failed: '//err
      return
    endif
    
    T_surf = 10.0_dp**x(1)
    call fcn(n, x, fvec, info)
    
  contains
    subroutine fcn(n_, x_, fvec_, iflag_)
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(n_)
      integer, intent(inout) :: iflag_
      real(dp) :: T, ISR, OLR
      T = 10.0_dp**x_(1)
      call self%TOA_fluxes_column(T, N_i_surf, ISR, OLR, err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif
      fvec_(1) = ISR - OLR
    end subroutine
    
  end function
  
end module