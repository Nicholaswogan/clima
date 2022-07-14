submodule(clima_climate) clima_climate_rhs
  implicit none
  
  
contains

  module subroutine right_hand_side(self, T_in, dTdt, err)
    use clima_eqns, only: press_and_den, heat_capacity_eval
    use clima_const, only: k_boltz, N_avo
    class(Climate), intent(inout), target :: self
    real(dp), intent(in) :: T_in(:)
    real(dp), intent(out) :: dTdt(:)
    character(:), allocatable :: err
    
    integer :: i, j
    
    logical :: found
    real(dp) :: cp(self%nz), cp_tmp
    real(dp) :: adiabat_lapse(self%nz)
    real(dp) :: scale_height(self%nz)
    real(dp) :: rho(self%nz)
    real(dp) :: Fc_e(self%nz-1)
    real(dp) :: Latent_heat(self%nz)
    
    real(dp) :: delta_z

    type(ClimateWrk), pointer :: w
    
    w => self%wrk
    
    w%T = T_in
    
    if (self%switch) then
      ! First, we compute P profile
      call press_and_den(self%nz, w%T, self%grav, &
                        self%surface_pressure*1.e6_dp, self%dz, &
                        self%mubar, w%P, w%density)
      ! convert P from cgs to bar
      w%P = w%P/1.0e6_dp
      self%switch = .false.
    endif
    
    ! Compute densities and columns
    do i = 1,self%sp%ng
      w%densities(:,i) = self%mix(:,i)*w%density(:)
    enddo
    ! density (g/cm3)
    rho = w%density*(1.0_dp/N_avo)*self%mubar
    
    ! Radiative transfer
    call self%rad%radiate(w%T(1), w%T, w%P, w%densities, self%dz, err)
    if (allocated(err)) return

    ! Heat capacity in erg/(g*K)
    do j = 1,self%nz
      cp(j) = 0.0_dp
      do i = 1,self%sp%ng
        call heat_capacity_eval(self%sp%g(i)%thermo, w%T(j), found, cp_tmp)
        if (.not. found) then
          err = "not found"
          return
        endif
        ! J/(kg*K) (si units)
        cp(j) = cp(j) + cp_tmp*self%mix(j,i)*(1.0_dp/(self%mubar(j)*1.0e-3_dp))
      enddo
      ! convert to erg/(g*K)
      cp(j) = cp(j)*1.0e4_dp
    enddo
    
    ! Adiabatic lapse rate (K/cm)
    adiabat_lapse(:) = (self%grav(:)/cp(:))
    
    ! Scale height at all altitudes (cm)
    scale_height(:) = (k_boltz*w%T(:)*N_avo)/(self%mubar(:)*self%grav(:))
    
    ! Convection
    call convection_diffusion(w%T, self%grav, self%z, self%dz, cp, &
                              rho, adiabat_lapse, scale_height, Fc_e)
    
    ! Latent heating (ergs/(cm2 s))
    do j = 1,self%nz
      Latent_heat(j) = 0.0_dp
      ! do i = 1,d%ng
      !   Latent_heat(j) = Latent_heat(j) + Latent(i,j)*con_rate(i,j)
      ! enddo 
    enddo

    Fc_e(:) = 0.0_dp
    
    !!! Right hand side (K/s) !!!
    ! center grid points
    do j = 2,self%nz-1
      dTdt(j) = &
        (1.0_dp/(rho(j)*cp(j)))*(self%rad%f_total(j+1) - self%rad%f_total(j))/self%dz(j) + &
        (-1.0_dp/(rho(j)*cp(j)))*(Fc_e(j) - Fc_e(j-1))/self%dz(j) + &
        (1.0_dp/(rho(j)*cp(j)))*Latent_heat(j)
    enddo
    ! lower boundary
    j = 1
    dTdt(j) = (1.0_dp/(rho(j)*cp(j)))*(self%rad%f_total(j+1) - self%rad%f_total(j))/self%dz(j) + &
              (-1.0_dp/(rho(j)*cp(j)))*(Fc_e(1)/self%dz(j) - 0.0_dp) + &
              (1.0_dp/(rho(j)*cp(j)))*Latent_heat(j)
        
    ! upper boundary
    j = self%nz
    dTdt(j) = (1.0_dp/(rho(j)*cp(j)))*(self%rad%f_total(j+1) - self%rad%f_total(j))/self%dz(j) + &
              (-1.0_dp/(rho(j)*cp(j)))*(0.0_dp - Fc_e(j-1)/self%dz(j)) + &
              (1.0_dp/(rho(j)*cp(j)))*Latent_heat(j)
      
      
  end subroutine
  
  subroutine convection_diffusion(T, grav, z, dz, cp, rho, adiabat_lapse, scale_height, Fc_e)
    use clima_const, only: von_karman_const
    use clima_eqns, only: eddy_for_heat
    real(dp), intent(in) :: T(:)
    real(dp), intent(in) :: grav(:)
    real(dp), intent(in) :: z(:)
    real(dp), intent(in) :: dz(:)
    real(dp), intent(in) :: cp(:)
    real(dp), intent(in) :: rho(:)
    real(dp), intent(in) :: adiabat_lapse(:)
    real(dp), intent(in) :: scale_height(:)
    
    ! convection on the edges
    real(dp), intent(out) :: Fc_e(:)
    
    ! heat diffusion coefficient on all inner edges
    real(dp) :: Kh_e(size(T)-1)
    ! mixing length (cm)
    real(dp) :: free_length_prop_const
    real(dp) :: free_mixing_length(size(T))
    real(dp) :: mixing_length(size(T))
    
    real(dp) :: grav_av, T_av, adiabat_av, cp_av, rho_av
    real(dp) :: delta_z
    real(dp) :: mixing_length_av
    real(dp) :: dTdz
    
    integer :: nz
    integer :: i
    
    nz = size(T)
    
    free_length_prop_const = 1.0_dp
    free_mixing_length = free_length_prop_const*scale_height ! cm
    ! cm
    mixing_length = von_karman_const*z/(1.0_dp + von_karman_const*z/free_mixing_length) 
    
    do i = 1,nz-1
      rho_av = sqrt(rho(i)*rho(i+1))
      cp_av = sqrt(cp(i)*cp(i+1))
      grav_av = sqrt(grav(i)*grav(i+1))
      T_av = sqrt(T(i)*T(i+1))
      adiabat_av = sqrt(adiabat_lapse(i)*adiabat_lapse(i+1))
      mixing_length_av = sqrt(mixing_length(i)*mixing_length(i+1))
      delta_z = dz(i)/2.0_dp + dz(i+1)/2.0_dp 
      
      dTdz = (T(i+1)-T(i))/delta_z
      
      Kh_e(i) = eddy_for_heat(mixing_length_av, grav_av, T_av, dTdz, adiabat_av)
      
      Fc_e(i) = - (rho_av*cp_av*Kh_e(i))*(dTdz + adiabat_av)

    enddo
 
  end subroutine
  
end submodule