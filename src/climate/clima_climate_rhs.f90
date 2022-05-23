submodule(clima_climate) clima_climate_rhs
  implicit none
  
  
contains

  module subroutine right_hand_side(self, T_in, dTdt, err)
    use clima_eqns, only: press_and_den, heat_capacity_eval
    use clima_radtran, only: radiate
    use clima_const, only: k_boltz, N_avo
    class(Climate), intent(inout), target :: self
    real(dp), intent(in) :: T_in(:)
    real(dp), intent(out) :: dTdt(:)
    character(:), allocatable :: err
    
    integer :: i, j
    
    logical :: found
    real(dp) :: cp(self%v%nz), cp_tmp
    real(dp) :: adiabat_lapse(self%v%nz)
    real(dp) :: scale_height(self%v%nz)
    real(dp) :: rho(self%v%nz)
    real(dp) :: Fc_e(self%v%nz-1)
    real(dp) :: Latent_heat(self%v%nz)
    
    real(dp) :: delta_z
    
    type(ClimaData), pointer :: d
    type(ClimaVars), pointer :: v
    type(ClimaWrk), pointer :: w
    
    d => self%d
    v => self%v
    w => self%w
    
    w%rin%T = T_in
    ! print*,minval(T_in),maxval(T_in)
    
    ! First, we compute P profile
    call press_and_den(v%nz, w%rin%T, v%grav, &
                       v%surface_pressure*1.e6_dp, v%dz, &
                       v%mubar, w%rin%P, w%rin%density)
    ! convert P from cgs to bar
    w%rin%P = w%rin%P/1.0e6_dp
    
    ! Compute densities and columns
    do i = 1,d%ng
      w%rin%densities(:,i) = v%mix(:,i)*w%rin%density(:)
      w%rin%cols(:,i) = w%rin%densities(:,i)*v%dz(:)
    enddo
    ! density (g/cm3)
    rho = w%rin%density*(1.0_dp/N_avo)*v%mubar
    
    ! Solar radiative transfer
    call radiate(d%sol, d%kset, &
                 v%surface_albedo, v%u0, v%diurnal_fac, v%photons_sol, &
                 w%rin%P, w%rin%T, w%rin%densities, v%dz, &
                 w%rx_sol, w%rz, &
                 w%fup_a_sol, w%fdn_a_sol, w%fup_sol, w%fdn_sol)
      
    ! IR radiative transfer
    call radiate(d%ir, d%kset, &
                 v%surface_albedo, v%u0, v%diurnal_fac, v%photons_sol, &
                 w%rin%P, w%rin%T, w%rin%densities, v%dz, &
                 w%rx_ir, w%rz, &
                 w%fup_a_ir, w%fdn_a_ir, w%fup_ir, w%fdn_ir)
                 
    ! Total flux at edges of layers (ergs/(cm2 s)).
    ! Index 1 is bottom. Index nz+1 is top edge of top layer.
    ! mW/m2
    w%f_total = w%fdn_sol - w%fup_sol + w%fdn_ir - w%fup_ir
    ! convert to ergs/(cm2 s) (no conversion necessary)
    ! w%f_total = w%f_total*1.0_dp
    
    ! Heat capacity in erg/(g*K)
    do j = 1,v%nz
      cp(j) = 0.0_dp
      do i = 1,d%ng
        call heat_capacity_eval(d%sp(i)%thermo, w%rin%T(j), found, cp_tmp)
        if (.not. found) then
          err = "not found"
          return
        endif
        ! J/(kg*K) (si units)
        cp(j) = cp(j) + cp_tmp*v%mix(j,i)*(1.0_dp/(v%mubar(j)*1.0e-3_dp))
      enddo
      ! convert to erg/(g*K)
      cp(j) = cp(j)*1.0e4_dp
    enddo
    
    ! Adiabatic lapse rate (K/cm)
    adiabat_lapse = (v%grav/cp)
    
    ! Scale height at all altitudes (cm)
    scale_height = (k_boltz*w%rin%T*N_avo)/(v%mubar*v%grav)
    
    ! Convection
    call convection_diffusion(w%rin%T, v%grav, v%z, v%dz, cp, &
                              rho, adiabat_lapse, scale_height, Fc_e)
    
    ! Latent heating (ergs/(cm2 s))
    do j = 1,v%nz
      Latent_heat(j) = 0.0_dp
      ! do i = 1,d%ng
      !   Latent_heat(j) = Latent_heat(j) + Latent(i,j)*con_rate(i,j)
      ! enddo 
    enddo
    
    !!! Right hand side (K/s) !!!
    ! center grid points
    do j = 2,v%nz-1
      
      dTdt(j) = &
        (1.0_dp/(rho(j)*cp(j)))*(w%f_total(j+1) - w%f_total(j))/v%dz(j) + &
        (-1.0_dp/(rho(j)*cp(j)))*(Fc_e(j) - Fc_e(j-1))/v%dz(j) + &
        (1.0_dp/(rho(j)*cp(j)))*Latent_heat(j)
    enddo
    ! lower boundary
    j = 1
    dTdt(j) = (1.0_dp/(rho(j)*cp(j)))*(w%f_total(j+1) - w%f_total(j))/v%dz(j) + &
              (-1.0_dp/(rho(j)*cp(j)))*(Fc_e(1)/v%dz(j) - 0.0_dp) + &
              (1.0_dp/(rho(j)*cp(j)))*Latent_heat(j)
        
        
    ! print*,(-1.0_dp/(rho(j)*cp(j)))*(Fc_e(1)/v%dz(j) - 0.0_dp), &
    !   (-1.0_dp/(rho(j)*cp(j)))*(w%f_total(j+1) - w%f_total(j))/v%dz(j)   
              
    
    ! upper boundary
    j = v%nz
    dTdt(j) = (1.0_dp/(rho(j)*cp(j)))*(w%f_total(j+1) - w%f_total(j))/v%dz(j) + &
              (-1.0_dp/(rho(j)*cp(j)))*(0.0_dp - Fc_e(v%nz-1)/v%dz(j)) + &
              (1.0_dp/(rho(j)*cp(j)))*Latent_heat(j)
      
    ! dTdt(j) = 0.0_dp        
    ! do i = 1,v%nz
    !   print*,T_in(i),dTdt(i)
    ! 
    ! enddo    
    ! print*,''
    
    ! print*,
    
      
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