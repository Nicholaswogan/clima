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
    real(dp) :: Fc_e(self%nz-1), Fc_g
    real(dp) :: Latent_heat(self%nz)
    real(dp) :: dTdt_l(self%nz)
    real(dp) :: dFdz(self%nz)

    real(dp), parameter :: cp_ground = 4.182e7_dp ! H2O, erg/(g*K)
    real(dp), parameter :: rho_ground = 1.0_dp ! H2O, g/cm3
    real(dp), parameter :: dz_ground = 500.0_dp ! cm

    type(ClimateWrk), pointer :: w
    
    w => self%wrk
    
    w%T(:) = T_in(2:)
    w%T_surf = T_in(1)
    
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

    if (self%double_radiative_grid) then
      do i = 1,self%nz
        w%T_r(2*(i-1)+1) = w%T(i)
        w%T_r(2*(i-1)+2) = w%T(i)

        w%P_r(2*(i-1)+1) = w%P(i)
        w%P_r(2*(i-1)+2) = w%P(i)

        w%densities_r(2*(i-1)+1,:) = w%densities(i,:)
        w%densities_r(2*(i-1)+2,:) = w%densities(i,:)
      enddo
    else
      w%T_r = w%T
      w%P_r = w%P
      w%densities_r = w%densities
    endif
    
    ! Radiative transfer
    call self%rad%radiate(w%T_surf, w%T_r(:), w%P_r, w%densities_r, self%dz_r, err=err)
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
    call convection_diffusion(w%T(:), w%T_surf, self%grav, self%z, self%dz, dz_ground, cp, cp_ground, &
                              rho, rho_ground, adiabat_lapse, scale_height, Fc_e, Fc_g)
    
    ! Latent heating (ergs/(cm2 s))
    do j = 1,self%nz
      Latent_heat(j) = 0.0_dp
      ! do i = 1,d%ng
      !   Latent_heat(j) = Latent_heat(j) + Latent(i,j)*con_rate(i,j)
      ! enddo 
    enddo

    ! Fc_e(:) = 0.0_dp

    if (self%double_radiative_grid) then
      do i = 1,self%nz
        dFdz(i) = (self%rad%f_total(2*(i-1)+3) - self%rad%f_total(2*(i-1)+1))/self%dz(i)
      enddo
    else
      do j = 1,self%nz
        dFdz(j) = (self%rad%f_total(j+1) - self%rad%f_total(j))/self%dz(j)
      enddo
    endif
    
    !!! Right hand side (K/s) !!!
    ! center grid points
    do j = 2,self%nz-1
      dTdt_l(j) = &
        (1.0_dp/(rho(j)*cp(j)))*dFdz(j) + &
        (-1.0_dp/(rho(j)*cp(j)))*(Fc_e(j) - Fc_e(j-1))/self%dz(j) + &
        (1.0_dp/(rho(j)*cp(j)))*Latent_heat(j)
    enddo
    ! lower layer
    j = 1
    dTdt_l(j) = (1.0_dp/(rho(j)*cp(j)))*dFdz(j) + &
                (-1.0_dp/(rho(j)*cp(j)))*(Fc_e(1) - Fc_g)/self%dz(j) + &
                (1.0_dp/(rho(j)*cp(j)))*Latent_heat(j)
        
    ! upper boundary
    j = self%nz
    dTdt_l(j) = (1.0_dp/(rho(j)*cp(j)))*dFdz(j) + &
              (-1.0_dp/(rho(j)*cp(j)))*(0.0_dp - Fc_e(j-1)/self%dz(j)) + &
              (1.0_dp/(rho(j)*cp(j)))*Latent_heat(j)
    
    ! ground level
    dTdt(1) = (1.0_dp/(rho_ground*cp_ground))*(self%rad%f_total(1))/dz_ground + &
              (-1.0_dp/(rho_ground*cp_ground))*(Fc_g/dz_ground - 0.0_dp) 

    ! append atmosphere
    dTdt(2:) = dTdt_l
      
  end subroutine
  
  subroutine convection_diffusion(T, T_surf, grav, z, dz, dz_ground, cp, cp_ground, rho, rho_ground, &
                                  adiabat_lapse, scale_height, Fc_e, Fc_g)
    use clima_const, only: von_karman_const
    use clima_eqns, only: eddy_for_heat
    real(dp), intent(in) :: T(:), T_surf
    real(dp), intent(in) :: grav(:)
    real(dp), intent(in) :: z(:)
    real(dp), intent(in) :: dz(:), dz_ground
    real(dp), intent(in) :: cp(:), cp_ground
    real(dp), intent(in) :: rho(:), rho_ground
    real(dp), intent(in) :: adiabat_lapse(:)
    real(dp), intent(in) :: scale_height(:)
    
    ! convection on the edges
    real(dp), intent(out) :: Fc_e(:), Fc_g
    
    ! heat diffusion coefficient on all inner edges
    real(dp) :: Kh_e(size(T)-1), Kh_g
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

    ! surface layer 
    ! (layer between the ground and the first atmospheric layer)
    rho_av = sqrt(rho_ground*rho(1))
    cp_av = sqrt(cp_ground*cp(1))
    grav_av = grav(i)
    T_av = sqrt(T_surf*T(1))
    adiabat_av = adiabat_lapse(1)
    mixing_length_av = mixing_length(1)

    delta_z = dz_ground/2.0_dp + dz(1)/2.0_dp 

    dTdz = (T(1)-T_surf)/delta_z

    Kh_g = eddy_for_heat(mixing_length_av, grav_av, T_av, dTdz, adiabat_av)
    Fc_g = - (rho_av*cp_av*Kh_g)*(dTdz + adiabat_av)

  end subroutine
  
end submodule