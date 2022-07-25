
module clima_radtran_types
  use iso_c_binding
  use clima_const, only: dp
  use linear_interpolation_module, only: linear_interp_1d, linear_interp_2d
  implicit none
  public

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
    real(dp), allocatable :: weight_e(:) ! (ngauss+1) bin edges
    real(dp), allocatable :: weights(:) ! (ngauss)
    real(dp), allocatable :: log10P(:) ! (nwav) log10(bars)
    real(dp) :: log10P_min
    real(dp) :: log10P_max
    real(dp), allocatable :: temp(:) ! (ntemp) kelvin
    real(dp) :: T_min
    real(dp) :: T_max
    type(linear_interp_2d), allocatable :: log10k(:,:) ! (ngauss, nwav)
  end type
  
  enum, bind(c)
    enumerator :: CIAXsection, RayleighXsection, AbsorptionXsection, PhotolysisXsection
  end enum
  
  type :: Xsection
    integer :: xs_type ! see enum above.
    integer :: dim ! 0, 1
    integer, allocatable :: sp_ind(:) ! (1 or 2). species indexes
    integer, allocatable :: rxn ! photolysis reaction number (only for PhotolysisXsection)
    integer, allocatable :: ntemp ! number of temperatures (only for xs_dim = 1)
    real(dp), allocatable :: temp(:) ! (ntemp) Kelvin
    real(dp), allocatable :: T_min
    real(dp), allocatable :: T_max
    real(dp), allocatable :: xs_0d(:) ! (nw) 
    type(linear_interp_1d), allocatable :: log10_xs_1d(:) ! (nw) 
  end type

  type :: ParticleXsection
    integer :: p_ind
    integer :: nrad
    real(dp), allocatable :: radii(:) ! cm
    real(dp) :: r_min, r_max
    type(linear_interp_1d), allocatable :: w0(:) ! (nw)
    type(linear_interp_1d), allocatable :: qext(:) ! (nw)
    type(linear_interp_1d), allocatable :: gt(:) ! (nw)
  end type
  
  type :: WaterContinuum
    integer :: LH2O ! index of H2O
    integer :: ntemp ! number of temperatures
    real(dp), allocatable :: temp(:) ! (ntemp) Kelvin
    real(dp) :: T_min
    real(dp) :: T_max
    type(linear_interp_1d), allocatable :: log10_xs_H2O(:) ! (nw) 
    type(linear_interp_1d), allocatable :: log10_xs_foreign(:) ! (nw) 
  end type
  
  ! integer :: k_method 
  enum, bind(c)
    enumerator :: K_RandomOverlap, k_RandomOverlapResortRebin
  end enum
  
  type :: Ksettings
    ! approach to combining k-distributions
    integer :: k_method 
    ! if k_method == k_RandomOverlapResortRebin
    ! then these are the weights we are re-binning to.
    integer :: nbin = -1
    real(dp), allocatable :: wbin(:) ! (nbin)
    real(dp), allocatable :: wbin_e(:) ! (nbin+1)
  end type
  
  enum, bind(c)
    enumerator :: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
  end enum
  
  type :: OpticalProperties
    integer :: op_type
    integer :: nw
    real(dp), allocatable :: wavl(:)
    real(dp), allocatable :: freq(:)
    
    ! K-distributions (e.g. H2O)
    integer :: ngauss_max = -1
    integer :: nk
    type(Ksettings) :: kset
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

    ! Particles
    integer :: npart
    type(ParticleXsection), allocatable :: part(:)
    
    ! Water Continuum absorption. Can only be a thing if H2O is a species
    ! and if there are other gases in the atmosphere.
    type(WaterContinuum), allocatable :: cont
    
  end type
  
  interface
    module function create_OpticalProperties(datadir, optype, species_names, particle_names, sop, err) result(op)
      use clima_types, only: SettingsOpacity
      character(*), intent(in) :: datadir
      integer, intent(in) :: optype
      character(*), intent(in) :: species_names(:)
      character(*), intent(in) :: particle_names(:)
      type(SettingsOpacity), intent(in) :: sop
      character(:), allocatable, intent(out) :: err
      type(OpticalProperties) :: op
    end function
  end interface
  interface OpticalProperties
    module procedure :: create_OpticalProperties
  end interface
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Work arrays for optical properties  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type :: Kcoefficients
    real(dp), allocatable :: k(:,:) ! (nz, ngauss)
  end type
  
  type :: RadiateXSWrk
    type(Kcoefficients), allocatable :: ks(:)
    real(dp), allocatable :: cia(:,:)
    real(dp), allocatable :: axs(:,:)
    real(dp), allocatable :: pxs(:,:)
    real(dp), allocatable :: H2O(:), foreign(:)
    
    ! work arrays that are needed only if
    ! k_method == k_RandomOverlapResortRebin
    real(dp), allocatable :: tau_k(:,:) ! (nz,nbin)
    real(dp), allocatable :: tau_xy(:,:) ! (nz,nbin*ngauss_max)
    real(dp), allocatable :: wxy(:) ! (nbin*ngauss_max)
    real(dp), allocatable :: wxy1(:) ! (nbin*ngauss_max)
    real(dp), allocatable :: wxy_e(:) ! (nbin*ngauss_max+1)
    integer, allocatable :: inds(:) ! (nbin*ngauss_max)
  end type
  
  interface
    module function create_RadiateXSWrk(op, nz) result(rw)
      type(OpticalProperties), target, intent(in) :: op
      integer, intent(in) :: nz
      type(RadiateXSWrk) :: rw
    end function
  end interface
  interface RadiateXSWrk
    module procedure :: create_RadiateXSWrk
  end interface
  
  type :: RadiateZWrk
    real(dp), allocatable :: tausg(:)
    real(dp), allocatable :: taua(:)
    real(dp), allocatable :: taua_1(:)
    real(dp), allocatable :: tau(:)
    real(dp), allocatable :: w0(:)
    real(dp), allocatable :: gt(:)
    real(dp), allocatable :: amean(:)
    real(dp), allocatable :: fup1(:)
    real(dp), allocatable :: fdn1(:)
    real(dp), allocatable :: fup(:)
    real(dp), allocatable :: fdn(:)
    real(dp), allocatable :: bplanck(:)
  end type
  
  interface
    module function create_RadiateZWrk(nz) result(rz)
      integer, intent(in) :: nz
      type(RadiateZWrk) :: rz
    end function
  end interface
  interface RadiateZWrk
    module procedure :: create_RadiateZWrk
  end interface
  
  ! for reading the stellar flux
  interface
    module subroutine read_stellar_flux(star_file, nw, wavl, photon_flux, err)
      character(len=*), intent(in) :: star_file
      integer, intent(in) :: nw
      real(dp), intent(in) :: wavl(nw+1)
      real(dp), intent(out) :: photon_flux(nw)
      character(:), allocatable, intent(out) :: err
    end subroutine
  end interface
  
end module