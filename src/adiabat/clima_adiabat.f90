module clima_adiabat
  use clima_const, only: dp, s_str_len
  use clima_radtran, only: RadtranIR
  implicit none
  
  type :: KastingClimateModel

    ! settings and free parameters
    integer :: nz
    real(dp) :: P_top = 1.0_dp ! (dynes/cm2)
    real(dp) :: T_trop = 180.0_dp ! (T)
    
    ! planet properties
    real(dp) :: planet_mass = 5.972e27_dp ! Earth (g)
    real(dp) :: planet_radius = 6.371e8_dp ! Earth (cm)
    
    ! species in the model
    character(s_str_len) :: species_names(3) = ['H2O','CO2','N2 ']
    
    ! Radiative transfer
    type(RadtranIR) :: rad
    
    ! work variables
    logical :: surface_liquid_H2O
    real(dp), allocatable :: P(:), T(:), f_H2O(:), f_CO2(:), f_N2(:), z(:), dz(:)
    real(dp), allocatable :: densities(:,:)
    
  contains
    
    procedure :: make_profile => KastingClimateModel_make_profile
    procedure :: OLR => KastingClimateModel_OLR
  end type
  
  interface KastingClimateModel
    module procedure :: create_KastingClimateModel
  end interface
  
contains
  
  function create_KastingClimateModel(datadir, nz, err) result(c)
    use clima_types, only: ClimaSettings 
    character(*), intent(in) :: datadir
    integer, intent(in) :: nz
    character(:), allocatable, intent(out) :: err
    
    type(KastingClimateModel) :: c
    
    type(ClimaSettings) :: s
    
    ! create the settings
    allocate(s%ir)
    s%ir%k_method = "RandomOverlapResortRebin"
    s%ir%nbins = 16
    s%ir%k_distributions = ['H2O', 'CO2']  
    s%ir%cia = ['CO2-CO2']
    s%ir%rayleigh = ['H2O', 'CO2', 'N2 ']
    s%ir%water_continuum = "MT_CKD"
    
    c%nz = nz
    ! species names
    c%species_names = ['H2O','CO2','N2 '] 
    ! make the radtranir object
    c%rad = RadtranIR(datadir, c%species_names, s, nz, err)
    if (allocated(err)) return
    
    ! allocate work variables
    allocate(c%P(nz), c%T(nz), c%f_H2O(nz), c%f_CO2(nz), c%f_N2(nz), c%z(nz), c%dz(nz))
    allocate(c%densities(nz,3))
    
  end function
  
  subroutine KastingClimateModel_make_profile(self, T_surf, P_H2O_surf, P_CO2_surf, P_N2_surf, err)
    use clima_const, only: k_boltz
    use clima_adiabat_kasting, only: make_profile
    class(KastingClimateModel), intent(inout) :: self
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), target, intent(in) :: P_H2O_surf, P_CO2_surf, P_N2_surf !! dynes/cm2
    character(:), allocatable, intent(out) :: err
    
    real(dp), allocatable :: P_e(:), T_e(:), f_H2O_e(:), f_CO2_e(:), f_N2_e(:), z_e(:)
    real(dp), allocatable :: density(:)
    integer :: i
    
    allocate(P_e(self%nz+1), T_e(self%nz+1), f_H2O_e(self%nz+1))
    allocate(f_CO2_e(self%nz+1), f_N2_e(self%nz+1), z_e(self%nz+1))
    allocate(density(self%nz))
    
    call make_profile(T_surf, P_H2O_surf, P_CO2_surf, P_N2_surf, self%nz+1, &
                      self%planet_mass, self%planet_radius, self%P_top, self%T_trop, &
                      P_e, z_e, T_e, f_H2O_e, f_CO2_e, f_N2_e, self%surface_liquid_H2O, &
                      err)
    if (allocated(err)) return                  
    
    do i = 1,self%nz
      self%P(i) = sqrt(P_e(i)*P_e(i+1))
      self%T(i) = 0.5_dp*(T_e(i)+T_e(i+1))
      self%f_H2O(i) = sqrt(f_H2O_e(i)*f_H2O_e(i+1))
      self%f_CO2(i) = sqrt(f_CO2_e(i)*f_CO2_e(i+1))
      self%f_N2(i) = sqrt(f_N2_e(i)*f_N2_e(i+1))
      self%z(i) = 0.5_dp*(z_e(i)+z_e(i+1))
      self%dz(i) = z_e(i+1) - z_e(i)
    enddo
    
    density = self%P/(k_boltz*self%T)
    self%densities(:,1) = self%f_H2O*density
    self%densities(:,2) = self%f_CO2*density
    self%densities(:,3) = self%f_N2*density
  
  end subroutine
  
  function KastingClimateModel_OLR(self, T_surf, P_H2O_surf, P_CO2_surf, P_N2_surf, err) result(OLR)
    class(KastingClimateModel), intent(inout) :: self
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), target, intent(in) :: P_H2O_surf, P_CO2_surf, P_N2_surf !! dynes/cm2
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: OLR
    
    ! make atmosphere profile
    call self%make_profile(T_surf, P_H2O_surf, P_CO2_surf, P_N2_surf, err)
    if (allocated(err)) return
    
    ! Do radiative transfer
    ! MUST CONVERT P TO BARS
    OLR = self%rad%OLR(self%T, self%P/1.0e6_dp, self%densities, self%dz, err)
    if (allocated(err)) return
  
  end function

  
end module