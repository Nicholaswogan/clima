
module clima_radtran_types
  use iso_c_binding
  use clima_const, only: dp, s_str_len
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
    real(dp), allocatable :: temp(:) ! (ntemp) kelvin
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
    real(dp), allocatable :: xs_0d(:) ! (nw) 
    type(linear_interp_1d), allocatable :: log10_xs_1d(:) ! (nw) 
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
    type(Ksettings) :: kset
    integer :: ngauss_max = -1
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
    
    ! work arrays that are needed only if
    ! k_method == k_RandomOverlapResortRebin
    real(dp), allocatable :: tau_k(:,:) ! (nz,nbin)
    real(dp), allocatable :: tau_xy(:,:) ! (nz,nbin*ngauss_max)
    real(dp), allocatable :: wxy(:) ! (nbin*ngauss_max)
    real(dp), allocatable :: wxy1(:) ! (nbin*ngauss_max)
    real(dp), allocatable :: wxy_e(:) ! (nbin*ngauss_max+1)
    integer, allocatable :: inds(:) ! (nbin*ngauss_max)
  end type
  
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
  
contains
  
  function create_RadiateXSWrk(op, kset, nz) result(rw)
    type(OpticalProperties), intent(in) :: op
    type(Ksettings), intent(in) :: kset
    integer, intent(in) :: nz
    
    type(RadiateXSWrk) :: rw
    
    integer :: i
    
    allocate(rw%ks(op%nk))
    do i = 1,op%nk
      allocate(rw%ks(i)%k(nz,op%k(i)%ngauss))
    enddo
    allocate(rw%cia(nz,op%ncia))
    allocate(rw%axs(nz,op%naxs))
    allocate(rw%pxs(nz,op%npxs))
    
    ! if there are k-distributions
    ! then we need to allocate some work arrays
    if (op%nk /= 0) then
      if (kset%k_method == K_RandomOverlap) then
        ! no need to allocate anything
      elseif (kset%k_method == K_RandomOverlapResortRebin) then
        allocate(rw%tau_k(nz,kset%nbin))
        allocate(rw%tau_xy(nz,kset%nbin*op%ngauss_max))
        allocate(rw%wxy(kset%nbin*op%ngauss_max))
        allocate(rw%wxy1(kset%nbin*op%ngauss_max))
        allocate(rw%wxy_e(kset%nbin*op%ngauss_max+1))
        allocate(rw%inds(kset%nbin*op%ngauss_max))
      endif
    endif
    
  end function  
  
  function create_RadiateZWrk(nz) result(rz)
    integer, intent(in) :: nz
    
    type(RadiateZWrk) :: rz
  
    allocate(rz%tausg(nz))
    allocate(rz%taua(nz))
    allocate(rz%taua_1(nz))
    allocate(rz%tau(nz))
    allocate(rz%w0(nz))
    allocate(rz%gt(nz))
    allocate(rz%amean(nz+1))
    allocate(rz%fup1(nz+1))
    allocate(rz%fdn1(nz+1))
    allocate(rz%fup(nz+1))
    allocate(rz%fdn(nz+1))
    allocate(rz%bplanck(nz+1))
  end function
  
end module