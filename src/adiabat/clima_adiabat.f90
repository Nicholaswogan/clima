module clima_adiabat
  use clima_const, only: dp, s_str_len
  use clima_types, only: Species
  use clima_radtran, only: Radtran
  implicit none
  private

  public :: AdiabatClimate

  type :: AdiabatClimate

    ! settings and free parameters
    integer :: nz
    real(dp) :: P_top = 1.0e-2_dp ! (dynes/cm2)
    real(dp) :: T_trop = 180.0_dp ! (T)
    real(dp), allocatable :: RH(:) ! relative humidity (ng)
    
    ! planet properties
    real(dp) :: planet_mass ! (g)
    real(dp) :: planet_radius ! (cm)
    
    ! species in the model
    character(s_str_len), allocatable :: species_names(:) ! copy of species names
    integer :: LH2O
    type(Species) :: sp
    
    ! Radiative transfer
    type(Radtran) :: rad
    
    ! State of the atmosphere
    real(dp), allocatable :: P(:) !! pressure in each grid cell, dynes/cm^2 (nz)
    real(dp), allocatable :: T(:) !! Temperature in each grid cell, K (nz) 
    real(dp), allocatable :: f_i(:,:) !! mixing ratios of species in each grid cell (nz,ng)
    real(dp), allocatable :: z(:) !! Altitude at the center of the grid cell, cm (nz)
    real(dp), allocatable :: dz(:) !! Thickness of each grid cell, cm (nz)
    real(dp), allocatable :: densities(:,:) !! densities in each grid cell, molecules/cm^2 (nz,ng)
    real(dp), allocatable :: N_surface(:) !! reservoir of gas on surface mol/cm^2 (ng)
    
  contains
    procedure :: make_profile => AdiabatClimate_make_profile
    procedure :: make_column => AdiabatClimate_make_column
    procedure :: TOA_fluxes => AdiabatClimate_TOA_fluxes
    procedure :: TOA_fluxes_column => AdiabatClimate_TOA_fluxes_column
    procedure :: surface_temperature => AdiabatClimate_surface_temperature
    procedure :: surface_temperature_column => AdiabatClimate_surface_temperature_column
    procedure :: to_regular_grid => AdiabatClimate_to_regular_grid
    procedure :: out2atmosphere_txt => AdiabatClimate_out2atmosphere_txt
  end type
  
  interface AdiabatClimate
    module procedure :: create_AdiabatClimate
  end interface
  
contains
  
  function create_AdiabatClimate(datadir, species_f, settings_f, star_f, err) result(c)
    use clima_types, only: ClimaSettings 
    character(*), intent(in) :: datadir
    character(*), intent(in) :: species_f
    character(*), intent(in) :: settings_f
    character(*), intent(in) :: star_f
    character(:), allocatable, intent(out) :: err
    
    type(AdiabatClimate) :: c
    
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

    ! default relative humidty is 1
    allocate(c%RH(c%sp%ng))
    c%RH(:) = 1.0_dp
    
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
    allocate(c%N_surface(c%sp%ng))
    
  end function
  
  subroutine AdiabatClimate_make_profile(self, T_surf, P_i_surf, err)
    use clima_adiabat_general, only: make_profile
    use clima_const, only: k_boltz
    class(AdiabatClimate), intent(inout) :: self
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
    
    call make_profile(T_surf, P_i_surf, &
                      self%sp, self%nz, self%planet_mass, &
                      self%planet_radius, self%P_top, self%T_trop, self%RH, &
                      P_e, z_e, T_e, f_i_e, self%N_surface, &
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

  subroutine AdiabatClimate_make_column(self, T_surf, N_i_surf, err)
    use clima_adiabat_general, only: make_column
    use clima_const, only: k_boltz
    class(AdiabatClimate), intent(inout) :: self
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
    
    call make_column(T_surf, N_i_surf, &
                     self%sp, self%nz, self%planet_mass, &
                     self%planet_radius, self%P_top, self%T_trop, self%RH, &
                     P_e, z_e, T_e, f_i_e, self%N_surface, &
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
  
  subroutine AdiabatClimate_TOA_fluxes(self, T_surf, P_i_surf, ISR, OLR, err)
    class(AdiabatClimate), intent(inout) :: self
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

  subroutine AdiabatClimate_TOA_fluxes_column(self, T_surf, N_i_surf, ISR, OLR, err)
    class(AdiabatClimate), intent(inout) :: self
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
  
  function AdiabatClimate_surface_temperature(self, P_i_surf, T_guess, err) result(T_surf)
    use minpack_module, only: hybrd1
    class(AdiabatClimate), intent(inout) :: self
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

  function AdiabatClimate_surface_temperature_column(self, N_i_surf, T_guess, err) result(T_surf)
    use minpack_module, only: hybrd1
    class(AdiabatClimate), intent(inout) :: self
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
        err = 'subroutine conserving_rebin returned an error'
        return
      endif
    enddo

    allocate(T_new(self%nz))
    call interp(self%nz, self%nz, z_new, self%z, self%T, T_new, ierr)
    if (ierr /= 0) then
      err = 'Subroutine interp returned an error.'
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