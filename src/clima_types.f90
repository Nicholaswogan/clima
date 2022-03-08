
module clima_types
  use clima_const, only: dp, s_str_len
  use linear_interpolation_module, only: linear_interp_1d, linear_interp_2d
  implicit none
  public
  
  !!!!!!!!!!!!!!!!
  !!! Settings !!!
  !!!!!!!!!!!!!!!!
  
  type :: SettingsOpacity
    character(s_str_len), allocatable :: k_distributions(:)
    character(s_str_len), allocatable :: cia(:)
    character(s_str_len), allocatable :: rayleigh(:)
    logical, allocatable :: rayleigh_bool
    character(s_str_len), allocatable :: absorption_xs(:)
    character(s_str_len), allocatable :: photolysis_xs(:)
  end type
  
  type :: ClimaSettings
    
    type(SettingsOpacity) :: op
    
  end type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Optical Properties  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Two main different types: `Ktable` and `Xsection`.
  ! **`Ktable`** is a k-distribution table representative of line opacities.
  ! **`Xsection`** is a continuum opacity source There are 4 types (CIA, Rayleigh,
  ! Photolysis, Absorption), and can be grids, which are interpolated over 
  ! constant (no interpolation), or T or T and log10P.
  
  type :: Ktable
    integer :: sp_ind
    integer :: ngauss 
    integer :: npress
    integer :: ntemp
    integer :: nwav
    real(dp), allocatable :: weights(:) ! (ngauss)
    real(dp), allocatable :: log10P(:) ! (nwav) log10(bars)
    real(dp), allocatable :: temp(:) ! (ntemp) kelvin
    type(linear_interp_2d), allocatable :: log10k(:,:) ! (ngauss, nwav)
  end type
  
  enum, bind(c)
    enumerator :: CIAXsection, RayleighXsection, AbsorptionXsection, PhotolysisXsection
  end enum
  
  type :: Xsection
    integer :: xs_type ! see enum above.
    integer :: dim ! 0, 1, or 2.
    integer, allocatable :: sp_ind(:) ! (1 or 2). species indexes
    integer, allocatable :: rxn ! photolysis reaction number (only for PhotolysisXsection)
    integer, allocatable :: ntemp ! number of temperatures (only for xs_dim = 1 or 2)
    integer, allocatable :: npress ! number of pressure (only for xs_dim = 2)
    real(dp), allocatable :: temp(:) ! (ntemp) Kelvin
    real(dp), allocatable :: log10P(:) ! (npress) log10(bars)
    real(dp), allocatable :: log10_xs_0d(:) ! (nw) 
    type(linear_interp_1d), allocatable :: log10_xs_1d(:) ! (nw) 
    type(linear_interp_2d), allocatable :: log10_xs_2d(:) ! (nw)
  end type
  
  enum, bind(c)
    enumerator :: SolarOpticalProperties, IROpticalProperties
  end enum
  
  type :: OpticalProperties
    integer :: op_type
    integer :: nw
    real(dp), allocatable :: wavenums(:)
    
    ! K-distributions (e.g. H2O)
    integer :: nk
    type(Ktable), allocatable :: k(:)
    ! T-dependent CIA coefficients (e.g. H2-H2)
    integer :: ncia
    type(Xsection), allocatable :: cia(:)
    ! Rayleigh Scattering
    integer :: nray
    type(Xsection), allocatable :: ray(:)
    ! Absorption cross section
    integer :: naxs
    type(Xsection), allocatable :: axs(:)
    ! Photolysis cross sections (e.g. O3 photolysis)
    integer :: npxs
    type(Xsection), allocatable :: pxs(:)
    
  end type
  
  !!!!!!!!!!!!!!!
  !!! Species !!!
  !!!!!!!!!!!!!!!
  
  type :: ThermodynamicData
    integer :: dtype
    integer :: ntemps
    real(dp), allocatable :: temps(:)
    real(dp), allocatable :: data(:,:)
  end type
  
  type :: Species
    character(:), allocatable :: name
    integer, allocatable :: composition(:) ! (natoms)
    real(dp) :: mass
    
    ! thermodynamics
    type(ThermodynamicData) :: thermo
    
  end type

  type :: ClimaData
    
    integer :: natoms
    character(s_str_len), allocatable :: atoms_names(:)
    real(dp), allocatable :: atoms_mass(:)
    
    integer :: ng
    character(s_str_len), allocatable :: species_names(:) ! copy of species name
    type(Species), allocatable :: sp(:)
    
    !!! Optical properties !!!
    type(OpticalProperties) :: sol
    type(OpticalProperties) :: ir
    
  end type
  
  type :: ClimaVars
    
    character(:), allocatable :: data_dir
    integer :: nz
    real(dp), allocatable :: mix(:)
  
  end type
  
  type :: ClimaWrk
    ! work variables
  end type
  
end module