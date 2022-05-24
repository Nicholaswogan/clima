module clima_radtran
  use clima_const, only: dp, s_str_len
  use clima_radtran_types, only: OpticalProperties, RadiateXSWrk, RadiateZWrk
  implicit none
  
  ! initialize with
  ! settings.yaml, sun.txt, and data dir, species names
  ! 
  ! inputs are
  ! surface_albedo, u0, T, densities, dz
  ! ouptuts are
  ! fluxes (fup, fdn)
  
  type :: ClimaRadtranWrk
    
    ! work arrays
    type(RadiateXSWrk) :: rx_ir
    type(RadiateZWrk) :: rz_ir
    
    !! (nz+1,nw) mW/m2/Hz in each wavelength bin
    !! at the edges of the vertical grid
    real(dp), allocatable :: fup_a(:,:)
    real(dp), allocatable :: fdn_a(:,:)
    !! (nz+1) mW/m2 at the edges of the vertical grid 
    !! (integral of fup_a and fdn_a over wavelength grid)
    real(dp), allocatable :: fup_n(:)
    real(dp), allocatable :: fdn_n(:) 
    
  end type
  
  type :: ClimaRadtranIR
  
    integer :: ng
    character(s_str_len), allocatable :: species_names(:) ! (ng) copy of species names
  
    ! number of layers
    integer :: nz
  
    !!! Optical properties !!!
    type(OpticalProperties) :: ir
    
    type(ClimaRadtranWrk) :: wrk_ir
    
  
  end type
  
  type, extends(ClimaRadtranIR) :: ClimaRadtran

    !!! Optical properties !!!
    type(OpticalProperties) :: sol
  
    real(dp) :: diurnal_fac = 0.5_dp
    real(dp), allocatable :: photons_sol(:) ! (nw) mW/m2/Hz in each bin  
    
    type(ClimaRadtranWrk) :: wrk_sol
    
  end type
  
  interface ClimaRadtranIR
    module procedure :: create_ClimaRadtranIR
  end interface
  
contains
  
  function create_ClimaRadtranIR(datadir, settings_f, nz, err) result(rad)
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: datadir
    character(*), intent(in) :: settings_f
    integer, intent(in) :: nz
    character(:), allocatable, intent(out) :: err
    
    type(ClimaRadtranIR) :: rad
    
    type(ClimaSettings) :: s
    
    s = ClimaSettings(settings_f, err)
    if (allocated(err)) return
    
    if (.not. allocated(s%species)) then
      err = '"'//settings_f//'/optical-properties/species" does not exist'
      return
    endif
    
    rad = create_ClimaRadtranIR_(datadir, s%species, s, nz, err)
    if (allocated(err)) return
    
  end function
  
  function create_ClimaRadtranIR_(datadir, species_names, s, nz, err) result(rad)
    use clima_radtran_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: datadir
    character(s_str_len), intent(in) :: species_names(:)
    type(ClimaSettings), intent(in) :: s
    integer, intent(in) :: nz
    character(:), allocatable, intent(out) :: err
    
    type(ClimaRadtranIR) :: rad
    
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
    rad%wrk_ir%rx_ir = RadiateXSWrk(rad%ir, nz)
    rad%wrk_ir%rz_ir = RadiateZWrk(nz)
    allocate(rad%wrk_ir%fup_a(nz+1, rad%ir%nw))
    allocate(rad%wrk_ir%fdn_a(nz+1, rad%ir%nw))
    allocate(rad%wrk_ir%fup_n(nz+1))
    allocate(rad%wrk_ir%fdn_n(nz+1))

  end function
  
  subroutine ClimaRadtranIR_radiate(self, T, P, densities, dz, err)
    use clima_radtran_radiate, only: radiate
    class(ClimaRadtranIR), target, intent(inout) :: self
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    character(:), allocatable, intent(out) :: err
    
    type(ClimaRadtranWrk), pointer :: wrk 
    
    wrk => self%wrk_ir                                           
    call check_dimensions(self, T, P, densities, dz, err)
    if (allocated(err)) return
                                         
    call radiate(self%ir, &
                 0.0_dp, 0.0_dp, 0.0_dp, [0.0_dp], &
                 P, T, densities, dz, &
                 wrk%rx_ir, wrk%rz_ir, &
                 wrk%fup_a, wrk%fdn_a, wrk%fup_n, wrk%fdn_n)
    
  end subroutine
  
  function ClimaRadtranIR_OLR(self, T, P, densities, dz, err) result(res)
    class(ClimaRadtranIR), target, intent(inout) :: self
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: res
    
    call ClimaRadtranIR_radiate(self, T, P, densities, dz, err)
    if (allocated(err)) return
    res = self%wrk_ir%fup_n(self%nz+1)
    
  end function
  
  subroutine check_dimensions(self, T, P, densities, dz, err)
    class(ClimaRadtranIR), intent(inout) :: self
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    character(:), allocatable, intent(out) :: err
    
    if (size(T) /= self%nz) then
      err = '"T" has the wrong input dimension.'
      return
    endif
    if (size(P) /= self%nz) then
      err = '"P" has the wrong input dimension.'
      return
    endif
    if (size(densities,1) /= self%nz .or. size(densities,2) /= self%ng) then
      err = '"densities" has the wrong input dimension.'
      return
    endif
    if (size(dz) /= self%nz) then
      err = '"dz" has the wrong input dimension.'
      return
    endif
    
  end subroutine

  
end module
  
  