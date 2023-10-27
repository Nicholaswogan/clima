module clima_radtran
  use clima_const, only: dp, s_str_len
  use clima_radtran_types, only: OpticalProperties, RadiateXSWrk, RadiateZWrk
  implicit none
  private

  public :: ClimaRadtranWrk ! Work type that holds results
  public :: RadtranIR ! IR radiative transfer only
  public :: Radtran ! IR and solar radiative transfer
  
  type :: ClimaRadtranWrk
    
    ! work arrays
    type(RadiateXSWrk) :: rx
    type(RadiateZWrk) :: rz
    
    !> (nz+1,nw) mW/m2/Hz in each wavelength bin
    !> at the edges of the vertical grid
    real(dp), allocatable :: fup_a(:,:), fdn_a(:,:)
    !> (nz+1) mW/m2 at the edges of the vertical grid 
    !> (integral of fup_a and fdn_a over wavelength grid)
    real(dp), allocatable :: fup_n(:), fdn_n(:)
    !> Band optical thickness (nz,nw)
    real(dp), allocatable :: tau_band(:,:)
    
  end type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! IR Radiative Transfer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type :: RadtranIR
  
    integer :: ng
    character(s_str_len), allocatable :: species_names(:) ! (ng) copy of species names

    integer :: np
    character(s_str_len), allocatable :: particle_names(:)
  
    ! number of layers
    integer :: nz
  
    ! Optical properties
    type(OpticalProperties) :: ir
    
    ! work
    type(ClimaRadtranWrk) :: wrk_ir
    
  contains
    procedure :: radiate => RadtranIR_radiate
    procedure :: OLR => RadtranIR_OLR
  end type

  interface RadtranIR
    module procedure :: create_RadtranIR_1
    module procedure :: create_RadtranIR_2
  end interface
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! IR and Solar Radiative Transfer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type :: Radtran

    integer :: ng
    character(s_str_len), allocatable :: species_names(:) ! (ng) copy of species names

    integer :: np
    character(s_str_len), allocatable :: particle_names(:)
  
    ! number of layers
    integer :: nz
  
    !!! Optical properties !!!
    type(OpticalProperties) :: ir
    type(OpticalProperties) :: sol
    
    real(dp) :: diurnal_fac = 0.5_dp
    real(dp), allocatable :: zenith_u(:) ! cos(zenith_angle)
    real(dp), allocatable :: zenith_weights(:) ! 
    real(dp) :: surface_albedo
    real(dp), allocatable :: photons_sol(:) ! (nw) mW/m2/Hz in each bin  

    type(ClimaRadtranWrk) :: wrk_ir
    type(ClimaRadtranWrk) :: wrk_sol
    real(dp), allocatable :: f_total(:)
    
  contains
    procedure :: radiate => Radtran_radiate
    procedure :: TOA_fluxes => Radtran_TOA_fluxes
    procedure :: skin_temperature => Radtran_skin_temperature
    procedure :: equilibrium_temperature => Radtran_equilibrium_temperature
    procedure :: opacities2yaml => Radtran_opacities2yaml
    procedure :: print_opacities => Radtran_print_opacities
  end type

  interface Radtran
    module procedure :: create_Radtran_1
    module procedure :: create_Radtran_2
  end interface
  
contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! IR Radiative Transfer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  function create_RadtranIR_1(settings_f, nz, datadir, err) result(rad)
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: settings_f
    integer, intent(in) :: nz
    character(*), intent(in) :: datadir
    character(:), allocatable, intent(out) :: err
    
    type(RadtranIR) :: rad
    
    type(ClimaSettings) :: s
    
    s = ClimaSettings(settings_f, err)
    if (allocated(err)) return
    
    if (.not. allocated(s%gases)) then
      err = '"'//settings_f//'/optical-properties/species/gases" does not exist'
      return
    endif

    if (.not. allocated(s%particles)) allocate(s%particles(0))
    
    rad = create_RadtranIR_2(s%gases, s%particles, s, nz, datadir, err)
    if (allocated(err)) return
    
  end function
  
  function create_RadtranIR_2(species_names, particle_names, s, nz, datadir, err) result(rad)
    use clima_radtran_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: species_names(:)
    character(*), intent(in) :: particle_names(:)
    type(ClimaSettings), intent(in) :: s
    integer, intent(in) :: nz
    character(*), intent(in) :: datadir
    character(:), allocatable, intent(out) :: err
    
    type(RadtranIR) :: rad
    
    if (nz < 1) then
      err = '"nz" can not be less than 1.'
      return
    endif
    
    rad%ng = size(species_names)
    rad%species_names = species_names
    rad%np = size(particle_names)
    if (rad%np > 0) then
      rad%particle_names = particle_names
    endif
    rad%nz = nz
    
    if (.not. allocated(s%ir)) then
      err = '"'//s%filename//'/optical-properties/ir" does not exist.'
      return
    endif
    rad%ir = OpticalProperties(datadir, IROpticalProperties, species_names, particle_names, s%ir, s%wavelength_bins_file, err)
    if (allocated(err)) return
    
    ! work arrays
    rad%wrk_ir%rx = RadiateXSWrk(rad%ir, nz)
    rad%wrk_ir%rz = RadiateZWrk(nz,rad%ir%npart)
    allocate(rad%wrk_ir%fup_a(nz+1, rad%ir%nw))
    allocate(rad%wrk_ir%fdn_a(nz+1, rad%ir%nw))
    allocate(rad%wrk_ir%fup_n(nz+1))
    allocate(rad%wrk_ir%fdn_n(nz+1))
    allocate(rad%wrk_ir%tau_band(nz,rad%ir%nw))

  end function
  
  subroutine RadtranIR_radiate(self, T_surface, T, P, densities, dz, pdensities, radii, err)
    use clima_radtran_radiate, only: radiate
    class(RadtranIR), target, intent(inout) :: self
    real(dp), intent(in) :: T_surface
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    real(dp), optional, target, intent(in) :: pdensities(:,:), radii(:,:) !! (nz,np)
    character(:), allocatable, intent(out) :: err
    
    integer :: ierr
    type(ClimaRadtranWrk), pointer :: wrk 
    
    wrk => self%wrk_ir  
    
    call check_inputs(self%nz, self%ng, self%np, T, P, densities, dz, pdensities, radii, err)
    if (allocated(err)) return
                                         
    ierr = radiate(self%ir, &
                   0.0_dp, 0.0_dp, [0.0_dp], &
                   [0.0_dp], [1.0_dp], &
                   P, T_surface, T, densities, dz, &
                   pdensities, radii, &
                   wrk%rx, wrk%rz, &
                   wrk%fup_a, wrk%fdn_a, wrk%fup_n, wrk%fdn_n, wrk%tau_band)
    if (ierr /= 0) then
      err = 'Input particle radii are outside the data range.'
      return
    endif
    
  end subroutine
  
  function RadtranIR_OLR(self, T_surface, T, P, densities, dz, pdensities, radii, err) result(res)
    class(RadtranIR), target, intent(inout) :: self
    real(dp), intent(in) :: T_surface
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    real(dp), optional, target, intent(in) :: pdensities(:,:), radii(:,:) !! (nz,np)
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: res
    
    call RadtranIR_radiate(self, T_surface, T, P, densities, dz, pdensities, radii, err)
    if (allocated(err)) return
    res = self%wrk_ir%fup_n(self%nz+1)
    
  end function

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! IR and Solar Radiative Transfer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function create_Radtran_1(settings_f, star_f, num_zenith_angles, surface_albedo, nz, datadir, err) result(rad)
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: settings_f
    character(*), intent(in) :: star_f
    integer, intent(in) :: num_zenith_angles
    real(dp), intent(in) :: surface_albedo
    integer, intent(in) :: nz
    character(*), intent(in) :: datadir
    character(:), allocatable, intent(out) :: err
    
    type(Radtran) :: rad
    
    type(ClimaSettings) :: s
    
    s = ClimaSettings(settings_f, err)
    if (allocated(err)) return
    
    if (.not. allocated(s%gases)) then
      err = '"'//settings_f//'/optical-properties/gases" does not exist'
      return
    endif

    if (.not. allocated(s%particles)) allocate(s%particles(0))
    
    rad = create_Radtran_2(s%gases, s%particles, s, star_f, num_zenith_angles, surface_albedo, nz, datadir, err)
    if (allocated(err)) return
    
  end function
  
  function create_Radtran_2(species_names, particle_names, s, star_f, &
                            num_zenith_angles, surface_albedo, nz, datadir, err) result(rad)
    use clima_radtran_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties, &
                                   read_stellar_flux
    use clima_types, only: ClimaSettings
    use clima_eqns, only: zenith_angles_and_weights
    use clima_const, only: pi
    
    character(*), intent(in) :: species_names(:)
    character(*), intent(in) :: particle_names(:)
    type(ClimaSettings), intent(in) :: s
    character(*), intent(in) :: star_f
    integer, intent(in) :: num_zenith_angles
    real(dp), intent(in) :: surface_albedo
    integer, intent(in) :: nz
    character(*), intent(in) :: datadir
    character(:), allocatable, intent(out) :: err
    
    type(Radtran) :: rad
    
    if (nz < 1) then
      err = '"nz" can not be less than 1.'
      return
    endif
    
    rad%ng = size(species_names)
    rad%species_names = species_names
    rad%np = size(particle_names)
    if (rad%np > 0) then
      rad%particle_names = particle_names
    endif
    rad%nz = nz

    allocate(rad%zenith_u(num_zenith_angles))
    allocate(rad%zenith_weights(num_zenith_angles))
    call zenith_angles_and_weights(num_zenith_angles, rad%zenith_u, rad%zenith_weights)
    rad%zenith_u = cos(rad%zenith_u*pi/180.0_dp)

    rad%surface_albedo = surface_albedo
    
    if (.not. allocated(s%ir)) then
      err = '"'//s%filename//'/optical-properties/ir" does not exist.'
      return
    endif
    rad%ir = OpticalProperties(datadir, IROpticalProperties, species_names, particle_names, s%ir, s%wavelength_bins_file, err)
    if (allocated(err)) return

    if (.not. allocated(s%sol)) then
      err = '"'//s%filename//'/optical-properties/solar" does not exist.'
      return
    endif
    rad%sol = OpticalProperties(datadir, SolarOpticalProperties, species_names, particle_names, s%sol, s%wavelength_bins_file, err)
    if (allocated(err)) return

    ! photons hitting the planet
    allocate(rad%photons_sol(rad%sol%nw))
    call read_stellar_flux(star_f, rad%sol%nw, rad%sol%wavl, s%photon_scale_factor, rad%photons_sol, err)
    if (allocated(err)) return

    ! IR work arrays
    rad%wrk_ir%rx = RadiateXSWrk(rad%ir, nz)
    rad%wrk_ir%rz = RadiateZWrk(nz,rad%ir%npart)
    allocate(rad%wrk_ir%fup_a(nz+1, rad%ir%nw))
    allocate(rad%wrk_ir%fdn_a(nz+1, rad%ir%nw))
    allocate(rad%wrk_ir%fup_n(nz+1))
    allocate(rad%wrk_ir%fdn_n(nz+1))
    allocate(rad%wrk_ir%tau_band(nz,rad%ir%nw))

    ! Solar work arrays
    rad%wrk_sol%rx = RadiateXSWrk(rad%sol, nz)
    rad%wrk_sol%rz = RadiateZWrk(nz,rad%sol%npart)
    allocate(rad%wrk_sol%fup_a(nz+1, rad%sol%nw))
    allocate(rad%wrk_sol%fdn_a(nz+1, rad%sol%nw))
    allocate(rad%wrk_sol%fup_n(nz+1))
    allocate(rad%wrk_sol%fdn_n(nz+1))
    allocate(rad%wrk_sol%tau_band(nz,rad%sol%nw))

    ! total flux
    allocate(rad%f_total(nz+1))

  end function

  subroutine Radtran_radiate(self, T_surface, T, P, densities, dz, pdensities, radii , err)
    use clima_radtran_radiate, only: radiate
    use clima_const, only: pi
    class(Radtran), target, intent(inout) :: self
    real(dp), intent(in) :: T_surface
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), optional, target, intent(in) :: pdensities(:,:), radii(:,:)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    character(:), allocatable, intent(out) :: err

    integer :: ierr

    type(ClimaRadtranWrk), pointer :: wrk_ir
    type(ClimaRadtranWrk), pointer :: wrk_sol
    
    wrk_ir => self%wrk_ir
    wrk_sol => self%wrk_sol                                        
    
    call check_inputs(self%nz, self%ng, self%np, T, P, densities, dz, pdensities, radii, err)
    if (allocated(err)) return

    ! IR radiative transfer                                     
    ierr = radiate(self%ir, &
                  0.0_dp, 0.0_dp, [0.0_dp], &
                  [0.0_dp], [1.0_dp], &
                  P, T_surface, T, densities, dz, &
                  pdensities, radii, &
                  wrk_ir%rx, wrk_ir%rz, &
                  wrk_ir%fup_a, wrk_ir%fdn_a, wrk_ir%fup_n, wrk_ir%fdn_n, wrk_ir%tau_band)
    if (ierr /= 0) then
      err = 'Input particle radii are outside the data range.'
      return
    endif
    ! Solar radiative transfer
    ierr = radiate(self%sol, &
                   self%surface_albedo, self%diurnal_fac, self%photons_sol, &
                   self%zenith_u, self%zenith_weights, &
                   P, T_surface, T, densities, dz, &
                   pdensities, radii, &
                   wrk_sol%rx, wrk_sol%rz, &
                   wrk_sol%fup_a, wrk_sol%fdn_a, wrk_sol%fup_n, wrk_sol%fdn_n, wrk_sol%tau_band)
    if (ierr /= 0) then
      err = 'Input particle radii are outside the data range.'
      return
    endif
    ! Total flux at edges of layers (ergs/(cm2 s) which is the same as mW/m2).
    ! Index 1 is bottom. Index nz+1 is top edge of top layer.
    self%f_total = (wrk_sol%fdn_n - wrk_sol%fup_n) + (wrk_ir%fdn_n - wrk_ir%fup_n)

  end subroutine

  subroutine Radtran_TOA_fluxes(self, T_surface, T, P, densities, dz, pdensities, radii, ISR, OLR, err)
    use clima_radtran_radiate, only: radiate
    class(Radtran), target, intent(inout) :: self
    real(dp), intent(in) :: T_surface
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    real(dp), optional, target, intent(in) :: pdensities(:,:), radii(:,:) !! (nz,np)
    real(dp), intent(out) :: ISR, OLR
    character(:), allocatable, intent(out) :: err
      
    call self%radiate(T_surface, T, P, densities, dz, pdensities, radii , err)
    if (allocated(err)) return
    
    ISR = (self%wrk_sol%fdn_n(self%nz+1) - self%wrk_sol%fup_n(self%nz+1)) ! Incoming short wave
    OLR = - (self%wrk_ir%fdn_n(self%nz+1) - self%wrk_ir%fup_n(self%nz+1)) ! Outgoing long wave

  end subroutine

  function Radtran_skin_temperature(self, bond_albedo) result(T_skin)
    use clima_eqns, only: skin_temperature
    class(Radtran), target, intent(inout) :: self
    real(dp), intent(in) :: bond_albedo
    real(dp) :: T_skin

    real(dp) :: stellar_radiation
    integer :: i

    stellar_radiation = 0.0_dp
    do i = 1,self%sol%nw
      stellar_radiation = stellar_radiation + self%photons_sol(i)*(self%sol%freq(i) - self%sol%freq(i+1))
    enddo
    stellar_radiation = stellar_radiation/1.0e3_dp
    T_skin = skin_temperature(stellar_radiation, bond_albedo)
  end function

  function Radtran_equilibrium_temperature(self, bond_albedo) result(T_eq)
    use clima_eqns, only: equilibrium_temperature
    class(Radtran), target, intent(inout) :: self
    real(dp), intent(in) :: bond_albedo
    real(dp) :: T_eq

    real(dp) :: stellar_radiation
    integer :: i

    stellar_radiation = 0.0_dp
    do i = 1,self%sol%nw
      stellar_radiation = stellar_radiation + self%photons_sol(i)*(self%sol%freq(i) - self%sol%freq(i+1))
    enddo
    stellar_radiation = stellar_radiation/1.0e3_dp
    T_eq = equilibrium_temperature(stellar_radiation, bond_albedo)
  end function

  function Radtran_opacities2yaml(self) result(out)
    class(Radtran), target, intent(inout) :: self
    character(:), allocatable :: out

    character(:), allocatable :: line

    out = ''
    
    line = ''
    line = line//'optical-properties:'
    out = out//line

    out = out//'\n'
    line = '  '
    line = line//'ir:'
    out = out//line

    out = out//'\n'
    out = out//self%ir%opacities2yaml()

    out = out//'\n'
    line = '  '
    line = line//'solar:'
    out = out//line

    out = out//'\n'
    out = out//self%sol%opacities2yaml()

  end function

  subroutine Radtran_print_opacities(self)
    use iso_fortran_env, only: output_unit
    class(Radtran), target, intent(inout) :: self
    character(:), allocatable :: out

    integer :: i
    logical :: skip_next

    out = self%opacities2yaml()

    skip_next = .false.
    do i = 1,len(out)-1
      if (out(i:i+1) == '\n') then
        write(output_unit,'(a)')
        skip_next = .true.
      else
        if (.not. skip_next) then
          write(output_unit,'(a)',advance='no') out(i:i)
        else
          skip_next = .false.
        endif
      endif
    enddo
    write(output_unit,'(a)') out(len(out):len(out))

  end subroutine

  !!!!!!!!!!!!!!!!!
  !!! Utilities !!!
  !!!!!!!!!!!!!!!!!

  subroutine check_inputs(nz, ng, np, T, P, densities, dz, pdensities, radii, err)
    integer, intent(in) :: nz, ng, np
    real(dp), intent(in) :: T(:)
    real(dp), intent(in) :: P(:)
    real(dp), intent(in) :: densities(:,:)
    real(dp), intent(in) :: dz(:)
    real(dp), optional, target, intent(in) :: pdensities(:,:), radii(:,:)
    character(:), allocatable, intent(out) :: err

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
    
    call check_dimensions(nz, ng, T, P, densities, dz, err)
    if (allocated(err)) return
    if (present(radii)) then
      call check_dimensions_p(nz, np, pdensities, radii, err)
      if (allocated(err)) return
    endif

  end subroutine

  subroutine check_dimensions_p(nz, np, pdensities, radii, err)
    integer, intent(in) :: nz, np
    real(dp), intent(in) :: pdensities(:,:)
    real(dp), intent(in) :: radii(:,:)
    character(:), allocatable, intent(out) :: err

    if (size(pdensities,1) /= nz .or. size(pdensities,2) /= np) then
      err = '"pdensities" has the wrong input dimension.'
      return
    endif

    if (size(radii,1) /= nz .or. size(radii,2) /= np) then
      err = '"radii" has the wrong input dimension.'
      return
    endif

  end subroutine
  
  subroutine check_dimensions(nz, ng, T, P, densities, dz, err)
    integer, intent(in) :: nz, ng
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    character(:), allocatable, intent(out) :: err
    
    if (size(T) /= nz) then
      err = '"T" has the wrong input dimension.'
      return
    endif
    if (size(P) /= nz) then
      err = '"P" has the wrong input dimension.'
      return
    endif
    if (size(densities,1) /= nz .or. size(densities,2) /= ng) then
      err = '"densities" has the wrong input dimension.'
      return
    endif
    if (size(dz) /= nz) then
      err = '"dz" has the wrong input dimension.'
      return
    endif
    
  end subroutine
  
end module
  
  