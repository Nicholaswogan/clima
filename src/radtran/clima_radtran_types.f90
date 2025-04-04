
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
    character(:), allocatable :: dat_name
    integer :: p_ind
    integer :: nrad
    real(dp), allocatable :: radii(:) ! cm
    real(dp) :: r_min, r_max
    type(linear_interp_1d), allocatable :: w0(:) ! (nw)
    type(linear_interp_1d), allocatable :: qext(:) ! (nw)
    type(linear_interp_1d), allocatable :: gt(:) ! (nw)
  end type
  
  type :: WaterContinuum
    character(:), allocatable :: model
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
    enumerator :: K_RandomOverlap, k_RandomOverlapResortRebin, k_AdaptiveEquivalentExtinction
  end enum
  
  type :: Ksettings
    ! approach to combining k-distributions
    character(:), allocatable :: k_method_name ! name
    integer :: k_method ! enum (see above)
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

    ! Copy of species and particle names
    character(s_str_len), allocatable :: species_names(:)
    character(s_str_len), allocatable :: particle_names(:)
    
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

    ! Custom optical properties
    type(linear_interp_1d), private, allocatable :: dtau_dz_interp(:) ! nw
    type(linear_interp_1d), private, allocatable :: w0_interp(:) ! nw
    type(linear_interp_1d), private, allocatable :: g0_interp(:) ! nw

  contains
    procedure :: opacities2yaml => OpticalProperties_opacities2yaml
    procedure :: set_custom_optical_properties => OpticalProperties_set_custom_optical_properties
    procedure :: unset_custom_optical_properties => OpticalProperties_unset_custom_optical_properties
    procedure :: custom_optical_properties => OpticalProperties_custom_optical_properties
  end type
  
  interface
    module function create_OpticalProperties(datadir, optype, species_names, &
                                             particle_names, sop, wavelength_bins_file, err) result(op)
      use clima_types, only: SettingsOpacity
      character(*), intent(in) :: datadir
      integer, intent(in) :: optype
      character(*), intent(in) :: species_names(:)
      character(*), intent(in) :: particle_names(:)
      type(SettingsOpacity), intent(in) :: sop
      character(:), allocatable, intent(in) :: wavelength_bins_file
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

    ! logical :: interpolated_particles = .false.
    real(dp), allocatable :: w0(:,:)
    real(dp), allocatable :: qext(:,:)
    real(dp), allocatable :: gt(:,:)
    
    ! work arrays that are needed only if
    ! k_method == k_RandomOverlapResortRebin
    real(dp), allocatable :: tau_k(:,:) ! (nz,nbin)
    real(dp), allocatable :: tau_xy(:,:) ! (nz,nbin*ngauss_max)
    real(dp), allocatable :: wxy(:) ! (nbin*ngauss_max)
    real(dp), allocatable :: wxy1(:) ! (nbin*ngauss_max)
    real(dp), allocatable :: wxy_e(:) ! (nbin*ngauss_max+1)
    integer, allocatable :: inds(:) ! (nbin*ngauss_max)
    ! end work arrays that are needed for k_RandomOverlapResortRebin

    ! work arrays that are only needed for if 
    ! k_method == k_AdapativeEquivalentExtinction
    real(dp), allocatable :: tau_grey(:,:) ! (nz,op%nk)
    real(dp), allocatable :: tau_grey_sum(:) ! (nz)
    integer, allocatable :: ind_major(:) ! (nz)
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

    real(dp), allocatable :: tausp(:)
    real(dp), allocatable :: tausp_1(:,:)
    real(dp), allocatable :: taup(:)

    real(dp), allocatable :: tauc(:), tausc(:), w0c(:), g0c(:)

    real(dp), allocatable :: tau(:)
    real(dp), allocatable :: w0(:)
    real(dp), allocatable :: gt(:)

    real(dp), allocatable :: tau_band(:) !! band optical thickness.
  
    real(dp), allocatable :: amean(:)
    real(dp), allocatable :: amean1(:)
    real(dp), allocatable :: amean2(:)
    real(dp), allocatable :: fup1(:)
    real(dp), allocatable :: fdn1(:)
    real(dp), allocatable :: fup2(:)
    real(dp), allocatable :: fdn2(:)
    real(dp), allocatable :: fup(:)
    real(dp), allocatable :: fdn(:)
    real(dp), allocatable :: bplanck(:)
  end type
  
  interface
    module function create_RadiateZWrk(nz, npart) result(rz)
      integer, intent(in) :: nz, npart
      type(RadiateZWrk) :: rz
    end function
  end interface
  interface RadiateZWrk
    module procedure :: create_RadiateZWrk
  end interface
  
  ! for reading the stellar flux
  interface
    module subroutine read_stellar_flux(star_file, nw, wavl, photon_scale_factor, photon_flux, err)
      character(len=*), intent(in) :: star_file
      integer, intent(in) :: nw
      real(dp), intent(in) :: wavl(nw+1)
      real(dp), intent(in) :: photon_scale_factor
      real(dp), intent(out) :: photon_flux(nw)
      character(:), allocatable, intent(out) :: err
    end subroutine
  end interface

contains
  
  function OpticalProperties_opacities2yaml(self) result(out)
    use clima_const, only: 
    class(OpticalProperties), intent(inout) :: self
    character(:), allocatable :: out

    character(:), allocatable :: line
    integer :: i

    out = ''

    line = '    '
    line = line//'k-method: '//self%kset%k_method_name
    out = out//line

    out = out//new_line('(a)')
    line = '    '
    line = line//'opacities:'
    out = out//line

    if (allocated(self%k)) then
      out = out//new_line('(a)')
      line = '      '
      line = line//'k-distributions: ['
      do i = 1,self%nk
        line = line//trim(self%species_names(self%k(i)%sp_ind))
        if (i /= self%nk) then
          line = line//', '
        endif
      enddo
      line = line//']'
      out = out//line
    endif

    if (allocated(self%cia)) then
      out = out//new_line('(a)')
      line = '      '
      line = line//'CIA: ['
      do i = 1,self%ncia
        line = line//trim(self%species_names(self%cia(i)%sp_ind(1)))// &
                '-'//trim(self%species_names(self%cia(i)%sp_ind(2)))
        if (i /= self%ncia) then
          line = line//', '
        endif
      enddo
      line = line//']'
      out = out//line
    endif

    if (allocated(self%ray)) then
      out = out//new_line('(a)')
      line = '      '
      line = line//'rayleigh: ['
      do i = 1,self%nray
        line = line//trim(self%species_names(self%ray(i)%sp_ind(1)))
        if (i /= self%nray) then
          line = line//', '
        endif
      enddo
      line = line//']'
      out = out//line
    endif

    if (allocated(self%pxs)) then
      out = out//new_line('(a)')
      line = '      '
      line = line//'photolysis-xs: ['
      do i = 1,self%npxs
        line = line//trim(self%species_names(self%pxs(i)%sp_ind(1)))
        if (i /= self%npxs) then
          line = line//', '
        endif
      enddo
      line = line//']'
      out = out//line
    endif

    if (allocated(self%cont)) then
      out = out//new_line('(a)')
      line = '      '
      line = line//'water-continuum: '//self%cont%model
      out = out//line
    endif

    if (allocated(self%part)) then
      out = out//new_line('(a)')
      line = '      '
      line = line//'particle-xs: ['
      do i = 1,self%npart
        line = line//'{name: '//trim(self%particle_names(self%part(i)%p_ind))
        line = line//', data: '//trim(self%part(i)%dat_name)//'}'
        if (i /= self%npart) then
          line = line//', '
        endif
      enddo
      line = line//']'
      out = out//line
    endif

  end function

  subroutine OpticalProperties_set_custom_optical_properties(self, wv, P, dtau_dz, w0, g0, err)
    use futils, only: interp
    class(OpticalProperties), intent(inout) :: self
    real(dp), intent(in) :: wv(:) !! Array of of wavelengths in nm
    real(dp), intent(in) :: P(:) !! Array of pressures in dynes/cm^2
    real(dp), intent(in) :: dtau_dz(:,:) !! (size(P),size(wv)), Optical depth per altitude.
    real(dp), intent(in) :: w0(:,:) !! (size(P),size(wv)), Single scattering albedo
    real(dp), intent(in) :: g0(:,:) !! (size(P),size(wv)), Asymetry parameter
    character(:), allocatable, intent(out) :: err
    
    real(dp), allocatable :: wv1(:), log10P(:), dtau_dz_tmp(:,:), w0_tmp(:,:), g0_tmp(:,:)
    integer :: i, j, ierr

    ! Check inputs
    if (any(wv <= 0.0_dp)) then
      err = 'All elements of `wv` must be larger than zero'
      return
    endif
    if (any(P <= 0.0_dp)) then
      err = 'All elements of `P` must be larger than zero'
      return
    endif
    if (size(P) /= size(dtau_dz,1)) then
      err = '`P` and `dtau_dz` have incompatible shapes'
      return
    endif
    if (size(wv) /= size(dtau_dz,2)) then
      err = '`wv` and `dtau_dz` have incompatible shapes'
      return
    endif
    if (size(P) /= size(w0,1)) then
      err = '`P` and `w0` have incompatible shapes'
      return
    endif
    if (size(wv) /= size(w0,2)) then
      err = '`wv` and `w0` have incompatible shapes'
      return
    endif
    if (size(P) /= size(g0,1)) then
      err = '`P` and `g0` have incompatible shapes'
      return
    endif
    if (size(wv) /= size(g0,2)) then
      err = '`wv` and `g0` have incompatible shapes'
      return
    endif

    allocate(wv1(size(self%wavl)-1))
    allocate(dtau_dz_tmp(size(P),size(self%wavl)-1))
    allocate(w0_tmp(size(P),size(self%wavl)-1))
    allocate(g0_tmp(size(P),size(self%wavl)-1))

    ! Median wavelength
    wv1 = 0.5_dp*(self%wavl(2:) + self%wavl(1:size(self%wavl)-1))

    ! At each pressure, interpolate to the wavelength grid
    do i = 1,size(P)
      j = size(P) + 1 - i
      call interp(wv1, wv, dtau_dz(i,:), dtau_dz_tmp(j,:), ierr=ierr)
      if (ierr /= 0) then
        err = 'Interpolation error in `set_custom_optical_properties`'
        return
      endif
      call interp(wv1, wv, w0(i,:), w0_tmp(j,:), ierr=ierr)
      if (ierr /= 0) then
        err = 'Interpolation error in `set_custom_optical_properties`'
        return
      endif
      call interp(wv1, wv, g0(i,:), g0_tmp(j,:), ierr=ierr)
      if (ierr /= 0) then
        err = 'Interpolation error in `set_custom_optical_properties`'
        return
      endif
    enddo

    ! Consider log10 pressure
    log10P = log10(P)
    log10P = log10P(size(log10P):1:-1)

    if (.not.allocated(self%dtau_dz_interp)) then
      allocate(self%dtau_dz_interp(self%nw))
      allocate(self%w0_interp(self%nw))
      allocate(self%g0_interp(self%nw))
    endif
   
    ! At each wavelength bin, create an interpolator for altitude
    do i = 1,size(wv1)
      call self%dtau_dz_interp(i)%initialize(log10P, dtau_dz_tmp(:,i), ierr)
      if (ierr /= 0) then
        err = 'Interpolation initialization error in `set_custom_optical_properties`'
        return
      endif
      call self%w0_interp(i)%initialize(log10P, w0_tmp(:,i), ierr)
      if (ierr /= 0) then
        err = 'Interpolation initialization error in `set_custom_optical_properties`'
        return
      endif
      call self%g0_interp(i)%initialize(log10P, g0_tmp(:,i), ierr)
      if (ierr /= 0) then
        err = 'Interpolation initialization error in `set_custom_optical_properties`'
        return
      endif
    enddo

  end subroutine

  subroutine OpticalProperties_unset_custom_optical_properties(self)
    class(OpticalProperties), intent(inout) :: self
    if (allocated(self%dtau_dz_interp)) then
      deallocate(self%dtau_dz_interp)
      deallocate(self%w0_interp)
      deallocate(self%g0_interp)
    endif
  end subroutine

  subroutine OpticalProperties_custom_optical_properties(self, log10P, dz, l, tau, w0, g0)
    class(OpticalProperties), intent(inout) :: self
    real(dp), intent(in) :: log10P(:) !! Pressure in dynes/cm^2
    real(dp), intent(in) :: dz(:) !! Layer thickness in cm
    integer, intent(in) :: l !! Wavelength index
    real(dp), intent(out) :: tau(:) !! (nz), Optical depth.
    real(dp), intent(out) :: w0(:) !! (nz), Single scattering albedo
    real(dp), intent(out) :: g0(:) !! (nz), Asymetry parameter

    integer :: j

    if (.not.allocated(self%dtau_dz_interp)) then
      tau = tiny(0.0_dp)
      w0 = tiny(0.0_dp)
      g0 = tiny(0.0_dp)
      return
    endif

    do j = 1,size(log10P)
      call self%dtau_dz_interp(l)%evaluate(log10P(j), tau(j))
      tau(j) = tau(j)*dz(j)
      call self%w0_interp(l)%evaluate(log10P(j), w0(j))
      call self%g0_interp(l)%evaluate(log10P(j), g0(j))
    enddo

  end subroutine
  
end module