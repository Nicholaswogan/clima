module clima_radtran
  use clima_const, only: dp, s_str_len
  use clima_radtran_types, only: OpticalProperties, RadiateXSWrk, RadiateZWrk
  implicit none
  private

  public :: RadtranIR ! IR radiative transfer only
  public :: Radtran ! IR and solar radiative transfer
  
  type :: ClimaRadtranWrk
    
    ! work arrays
    type(RadiateXSWrk) :: rx
    type(RadiateZWrk) :: rz
    
    !! (nz+1,nw) mW/m2/Hz in each wavelength bin
    !! at the edges of the vertical grid
    real(dp), allocatable :: fup_a(:,:)
    real(dp), allocatable :: fdn_a(:,:)
    !! (nz+1) mW/m2 at the edges of the vertical grid 
    !! (integral of fup_a and fdn_a over wavelength grid)
    real(dp), allocatable :: fup_n(:)
    real(dp), allocatable :: fdn_n(:) 
    
  end type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! IR Radiative Transfer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type :: RadtranIR
  
    integer :: ng
    character(s_str_len), allocatable :: species_names(:) ! (ng) copy of species names
  
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
  
    ! number of layers
    integer :: nz
  
    !!! Optical properties !!!
    type(OpticalProperties) :: ir
    type(OpticalProperties) :: sol
    
    real(dp) :: diurnal_fac = 0.5_dp
    real(dp) :: solar_zenith
    real(dp) :: surface_albedo
    real(dp), allocatable :: photons_sol(:) ! (nw) mW/m2/Hz in each bin  

    type(ClimaRadtranWrk) :: wrk_ir
    type(ClimaRadtranWrk) :: wrk_sol
    real(dp), allocatable :: f_total(:)
    
  contains
    procedure :: radiate =>Radtran_radiate
  end type

  interface Radtran
    module procedure :: create_Radtran_1
    module procedure :: create_Radtran_2
  end interface
  
contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! IR Radiative Transfer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  function create_RadtranIR_1(datadir, settings_f, nz, err) result(rad)
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: datadir
    character(*), intent(in) :: settings_f
    integer, intent(in) :: nz
    character(:), allocatable, intent(out) :: err
    
    type(RadtranIR) :: rad
    
    type(ClimaSettings) :: s
    
    s = ClimaSettings(settings_f, err)
    if (allocated(err)) return
    
    if (.not. allocated(s%species)) then
      err = '"'//settings_f//'/optical-properties/species" does not exist'
      return
    endif
    
    rad = create_RadtranIR_2(datadir, s%species, s, nz, err)
    if (allocated(err)) return
    
  end function
  
  function create_RadtranIR_2(datadir, species_names, s, nz, err) result(rad)
    use clima_radtran_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: datadir
    character(*), intent(in) :: species_names(:)
    type(ClimaSettings), intent(in) :: s
    integer, intent(in) :: nz
    character(:), allocatable, intent(out) :: err
    
    type(RadtranIR) :: rad
    
    if (nz < 1) then
      err = '"nz" can not be less than 1.'
      return
    endif
    
    rad%ng = size(species_names)
    rad%species_names = species_names
    rad%nz = nz
    
    if (.not. allocated(s%ir)) then
      err = '"'//s%filename//'/optical-properties/ir" does not exist.'
      return
    endif
    rad%ir = OpticalProperties(datadir, IROpticalProperties, species_names, s%ir, err)
    if (allocated(err)) return
    
    ! work arrays
    rad%wrk_ir%rx = RadiateXSWrk(rad%ir, nz)
    rad%wrk_ir%rz = RadiateZWrk(nz)
    allocate(rad%wrk_ir%fup_a(nz+1, rad%ir%nw))
    allocate(rad%wrk_ir%fdn_a(nz+1, rad%ir%nw))
    allocate(rad%wrk_ir%fup_n(nz+1))
    allocate(rad%wrk_ir%fdn_n(nz+1))

  end function
  
  subroutine RadtranIR_radiate(self, T_surface, T, P, densities, dz, err)
    use clima_radtran_radiate, only: radiate
    class(RadtranIR), target, intent(inout) :: self
    real(dp), intent(in) :: T_surface
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    character(:), allocatable, intent(out) :: err
    
    type(ClimaRadtranWrk), pointer :: wrk 
    
    wrk => self%wrk_ir                                           
    call check_dimensions(self%nz, self%ng, T, P, densities, dz, err)
    if (allocated(err)) return
                                         
    call radiate(self%ir, &
                 0.0_dp, 0.0_dp, 0.0_dp, [0.0_dp], &
                 P, T_surface, T, densities, dz, &
                 wrk%rx, wrk%rz, &
                 wrk%fup_a, wrk%fdn_a, wrk%fup_n, wrk%fdn_n)
    
  end subroutine
  
  function RadtranIR_OLR(self, T_surface, T, P, densities, dz, err) result(res)
    class(RadtranIR), target, intent(inout) :: self
    real(dp), intent(in) :: T_surface
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: res
    
    call RadtranIR_radiate(self, T_surface, T, P, densities, dz, err)
    if (allocated(err)) return
    res = self%wrk_ir%fup_n(self%nz+1)
    
  end function
  
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! IR and Solar Radiative Transfer !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function create_Radtran_1(datadir, settings_f, star_f, solar_zenith, surface_albedo, nz, err) result(rad)
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: datadir
    character(*), intent(in) :: settings_f
    character(*), intent(in) :: star_f
    real(dp), intent(in) :: solar_zenith, surface_albedo
    integer, intent(in) :: nz
    character(:), allocatable, intent(out) :: err
    
    type(Radtran) :: rad
    
    type(ClimaSettings) :: s
    
    s = ClimaSettings(settings_f, err)
    if (allocated(err)) return
    
    if (.not. allocated(s%species)) then
      err = '"'//settings_f//'/optical-properties/species" does not exist'
      return
    endif
    
    rad = create_Radtran_2(datadir, s%species, s, star_f, solar_zenith, surface_albedo, nz, err)
    if (allocated(err)) return
    
  end function
  
  function create_Radtran_2(datadir, species_names, s, star_f, solar_zenith, surface_albedo, nz, err) result(rad)
    use clima_radtran_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties, &
                                   read_stellar_flux
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: datadir
    character(*), intent(in) :: species_names(:)
    type(ClimaSettings), intent(in) :: s
    character(*), intent(in) :: star_f
    real(dp), intent(in) :: solar_zenith, surface_albedo
    integer, intent(in) :: nz
    character(:), allocatable, intent(out) :: err
    
    type(Radtran) :: rad
    
    if (nz < 1) then
      err = '"nz" can not be less than 1.'
      return
    endif
    
    rad%ng = size(species_names)
    rad%species_names = species_names
    rad%nz = nz

    rad%solar_zenith = solar_zenith
    rad%surface_albedo = surface_albedo
    
    if (.not. allocated(s%ir)) then
      err = '"'//s%filename//'/optical-properties/ir" does not exist.'
      return
    endif
    rad%ir = OpticalProperties(datadir, IROpticalProperties, species_names, s%ir, err)
    if (allocated(err)) return

    if (.not. allocated(s%sol)) then
      err = '"'//s%filename//'/optical-properties/solar" does not exist.'
      return
    endif
    rad%sol = OpticalProperties(datadir, SolarOpticalProperties, species_names, s%sol, err)
    if (allocated(err)) return

    ! photons hitting the planet
    allocate(rad%photons_sol(rad%sol%nw))
    call read_stellar_flux(star_f, rad%sol%nw, rad%sol%wavl, rad%photons_sol, err)
    if (allocated(err)) return

    ! IR work arrays
    rad%wrk_ir%rx = RadiateXSWrk(rad%ir, nz)
    rad%wrk_ir%rz = RadiateZWrk(nz)
    allocate(rad%wrk_ir%fup_a(nz+1, rad%ir%nw))
    allocate(rad%wrk_ir%fdn_a(nz+1, rad%ir%nw))
    allocate(rad%wrk_ir%fup_n(nz+1))
    allocate(rad%wrk_ir%fdn_n(nz+1))

    ! Solar work arrays
    rad%wrk_sol%rx = RadiateXSWrk(rad%sol, nz)
    rad%wrk_sol%rz = RadiateZWrk(nz)
    allocate(rad%wrk_sol%fup_a(nz+1, rad%sol%nw))
    allocate(rad%wrk_sol%fdn_a(nz+1, rad%sol%nw))
    allocate(rad%wrk_sol%fup_n(nz+1))
    allocate(rad%wrk_sol%fdn_n(nz+1))

    ! total flux
    allocate(rad%f_total(nz+1))

  end function

  subroutine Radtran_radiate(self, T_surface, T, P, densities, dz, err)
    use clima_radtran_radiate, only: radiate
    use clima_const, only: pi
    class(Radtran), target, intent(inout) :: self
    real(dp), intent(in) :: T_surface
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: u0

    type(ClimaRadtranWrk), pointer :: wrk_ir
    type(ClimaRadtranWrk), pointer :: wrk_sol
    
    wrk_ir => self%wrk_ir
    wrk_sol => self%wrk_sol                                        
    call check_dimensions(self%nz, self%ng, T, P, densities, dz, err)
    if (allocated(err)) return

    ! IR radiative transfer                                     
    call radiate(self%ir, &
                 0.0_dp, 0.0_dp, 0.0_dp, [0.0_dp], &
                 P, T_surface, T, densities, dz, &
                 wrk_ir%rx, wrk_ir%rz, &
                 wrk_ir%fup_a, wrk_ir%fdn_a, wrk_ir%fup_n, wrk_ir%fdn_n)
    
    ! Solar radiative transfer
    u0 = cos(self%solar_zenith*pi/180.0_dp)
    call radiate(self%sol, &
                 self%surface_albedo, u0, self%diurnal_fac, self%photons_sol, &
                 P, T_surface, T, densities, dz, &
                 wrk_sol%rx, wrk_sol%rz, &
                 wrk_sol%fup_a, wrk_sol%fdn_a, wrk_sol%fup_n, wrk_sol%fdn_n)
    
    ! Total flux at edges of layers (ergs/(cm2 s) which is the same as mW/m2).
    ! Index 1 is bottom. Index nz+1 is top edge of top layer.
    self%f_total = (wrk_sol%fup_n - wrk_sol%fdn_n) + (wrk_ir%fup_n - wrk_ir%fdn_n)

  end subroutine
  
end module
  
  