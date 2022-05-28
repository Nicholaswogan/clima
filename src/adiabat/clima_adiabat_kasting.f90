! submodule (clima_adiabat) clima_adiabat_kasting
module clima_adiabat_kasting
  use clima_const, only: dp
  use dop853_module, only: dop853_class
  implicit none
  ! private
  
  type, extends(dop853_class) :: dop853_custom
    real(dp) :: f_H2O_surf, f_CO2_dry, f_N2_dry
    integer :: j
    real(dp), pointer :: P_out(:), T_out(:)
    real(dp), pointer :: T_strat
    
    real(dp), allocatable :: P_strat
    real(dp), allocatable :: P_bound, T_bound
    integer :: stopping_reason
  end type
  
  !! molar masses (g/mol)
  real(dp), parameter :: MU_H2O = 18.0_dp
  real(dp), parameter :: MU_CO2 = 44.0_dp
  real(dp), parameter :: MU_N2 = 28.0_dp
  
  !! latent heat of H2O vaporization/condensation (erg/g)
  real(dp), parameter :: L_H2O = 2264e7_dp
  
  real(dp), parameter :: T_crit_H2O = 647.0_dp !! K


contains
  
  subroutine make_profile(T_surf, P_H2O_surf, P_CO2_surf, P_N2_surf, nz, P_top, T_strat, &
                          P_out, T_out, err)
    use stdlib_math, only: logspace
    real(dp), intent(in) :: T_surf !! K
    real(dp), intent(in) :: P_H2O_surf, P_CO2_surf, P_N2_surf !! dynes/cm2
    integer, intent(in) :: nz
    real(dp), intent(in) :: P_top, T_strat
    
    real(dp), intent(out) :: P_out(nz), T_out(nz)
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: P_H2O_sat_surf, P_H2O_start
    logical :: dry_start
    real(dp) :: P_surf, P_surf_dry
    real(dp) :: f_H2O_surf, f_CO2_dry, f_N2_dry
    
    ! Check if H2O liquid can exist on the surface
    if (T_surf > T_crit_H2O) then
      ! Above critical point of water.
      ! All H2O must be in the atmosphere.
      dry_start = .true.
      P_H2O_start = P_H2O_surf
    else
      ! Below critical point
      P_H2O_sat_surf = sat_pressure_H2O(T_surf)
      if (P_H2O_sat_surf > P_H2O_surf) then
        ! All water on surface is vaporized
        dry_start = .true.
        P_H2O_start = P_H2O_surf
      else
        ! Water is able to condense at the surface
        dry_start = .false.
        P_H2O_start = P_H2O_sat_surf
      endif
    endif
    
    P_surf = P_H2O_start + P_CO2_surf + P_N2_surf
    if (P_top > P_surf) then
      err = 'make_profile: "P_top" is bigger than the surface pressure'
      return
    endif
    ! mixing ratio of H2O at the surface
    f_H2O_surf = P_H2O_start/P_surf
    ! fraction of the dry atmosphere that is
    ! CO2 and N2
    P_surf_dry = P_CO2_surf + P_N2_surf
    f_CO2_dry = P_CO2_surf/P_surf_dry
    f_N2_dry = P_N2_surf/P_surf_dry
    
    ! Make P profile
    P_out = logspace(log10(P_surf),log10(P_top),nz)
    
    if (dry_start) then
      call make_profile_dry_start(T_surf, P_surf, f_H2O_surf, f_CO2_dry, f_N2_dry, &
                                  T_strat, P_top, P_out, T_out)
    else ! moist start
      call make_profile_moist_start(T_surf, P_surf, f_CO2_dry, f_N2_dry, &
                                    T_strat, P_top, P_out, T_out)
    endif
    
  end subroutine
  
  subroutine make_profile_moist_start(T_surf, P_surf, f_CO2_dry, f_N2_dry, &
                                      T_strat, P_top, P_out, T_out)
    real(dp), intent(in) :: T_surf !! K
    real(dp), intent(in) :: P_surf !! dynes/cm2
    real(dp), intent(in) :: f_CO2_dry, f_N2_dry
    real(dp), target, intent(in) :: T_strat, P_top
    real(dp), target, intent(inout) :: P_out(:)
    
    real(dp), target, intent(out) :: T_out(:)
    
    type(dop853_custom) :: dop_m
    logical :: status_ok
    integer :: idid, ind
    real(dp) :: P
    real(dp) :: T(1)
    
    call dop_m%initialize(fcn=dT_dP_moist_dop, solout=solout_moist, n=1, &
                          iprint=0, icomp=[1], status_ok=status_ok)
    if (.not. status_ok) then
      print*,'init failed'
      stop 1
    endif
    dop_m%f_CO2_dry = f_CO2_dry
    dop_m%f_N2_dry = f_N2_dry
    dop_m%T_strat => T_strat
    dop_m%P_out => P_out
    dop_m%T_out => T_out
    dop_m%T_out(1) = T_surf
    dop_m%j = 2
    dop_m%stopping_reason = 1

    T(1) = T_surf
    P = P_surf
    call dop_m%integrate(P, T, P_top, [1.0e-6_dp], [1.0e-6_dp], iout=2, idid=idid)
    if (idid < 0) then
      print*,'integration failed'
      stop 1
    endif
    if (dop_m%stopping_reason == 1) then
      ! we made it to P_top (finished the integration)
      ! Nothing to do.
    elseif (dop_m%stopping_reason == 2) then
      ! we made it to T_strat
      T_out(dop_m%j:) = T_strat
      
      ind = minloc(abs(P_out - dop_m%P_strat), 1)
      P_out(ind) = dop_m%P_strat
      T_out(ind) = T_strat
    endif
    
  end subroutine
  
  subroutine make_profile_dry_start(T_surf, P_surf, f_H2O_surf, f_CO2_dry, f_N2_dry, &
                                   T_strat, P_top, P_out, T_out)
    real(dp), intent(in) :: T_surf !! K
    real(dp), intent(in) :: P_surf !! dynes/cm2
    real(dp), intent(in) :: f_H2O_surf, f_CO2_dry, f_N2_dry
    real(dp), target, intent(in) :: T_strat, P_top
    real(dp), target, intent(inout) :: P_out(:)
    
    real(dp), target, intent(out) :: T_out(:)
    
    type(dop853_custom) :: dop, dop_m
    logical :: status_ok
    integer :: idid, ind
    
    real(dp) :: P
    real(dp) :: T(1)
    
    ! integrate dry adiabat upward
    call dop%initialize(fcn=dT_dP_dry_dop, solout=solout_dry, n=1, &
                        iprint=0, icomp=[1], status_ok=status_ok)
    if (.not. status_ok) then
      print*,'init failed'
      stop 1
    endif
    dop%f_H2O_surf = f_H2O_surf
    dop%f_CO2_dry = f_CO2_dry
    dop%f_N2_dry = f_N2_dry
    dop%j = 2
    dop%T_strat => T_strat
    dop%P_out => P_out
    dop%T_out => T_out
    dop%T_out(1) = T_surf
    dop%stopping_reason = 1
    
    T(1) = T_surf
    P = P_surf
    call dop%integrate(P, T, P_top, [1.0e-6_dp], [1.0e-6_dp], iout=2, idid=idid)
    if (idid < 0) then
      print*,'integration failed'
      stop 1
    endif
    ! why did we stop?
    if (dop%stopping_reason == 1) then
      ! we made it to P_top (finished the integration)
      ! nothing to do!
    elseif (dop%stopping_reason == 2) then
      ! we made it to T_strat
      T_out(dop_m%j:) = T_strat
      
      ind = minloc(abs(P_out - dop%P_strat), 1)
      P_out(ind) = dop%P_strat
      T_out(ind) = T_strat
      
    elseif (dop%stopping_reason == 3) then
      ! Hit the moist regime
      
      ! do another integration with a moist adiabat
      call dop_m%initialize(fcn=dT_dP_moist_dop, solout=solout_moist, n=1, &
                            iprint=0, icomp=[1], status_ok=status_ok)
      if (.not. status_ok) then
        print*,'init failed'
        stop 1
      endif
      dop_m%f_H2O_surf = f_H2O_surf
      dop_m%f_CO2_dry = f_CO2_dry
      dop_m%f_N2_dry = f_N2_dry
      dop_m%T_strat => T_strat
      dop_m%P_out => P_out
      dop_m%T_out => T_out
      dop_m%j = dop%j
      dop_m%stopping_reason = 1

      T(1) = dop%T_bound
      P = dop%P_bound
      call dop_m%integrate(P, T, P_top, [1.0e-6_dp], [1.0e-6_dp], iout=2, idid=idid)
      if (idid < 0) then
        print*,'integration failed'
        stop 1
      endif
      
      if (dop_m%stopping_reason == 1) then
        ! we made it to P_top (finished the integration)
        ! Nothing to do.
        ind = minloc(abs(P_out - dop%P_bound), 1)
        P_out(ind) = dop%P_bound
        T_out(ind) = dop%T_bound
      elseif (dop_m%stopping_reason == 2) then
        ! we made it to T_strat
        T_out(dop_m%j:) = T_strat
        
        ind = minloc(abs(P_out - dop%P_bound), 1)
        P_out(ind) = dop%P_bound
        T_out(ind) = dop%T_bound
        
        ind = minloc(abs(P_out - dop_m%P_strat), 1)
        P_out(ind) = dop_m%P_strat
        T_out(ind) = T_strat
        
      endif
      
    endif
    
  end subroutine
  
  subroutine find_stratosphere_boundary(dop, P_cur, P_old, P_strat)
    use minpack_module, only: lmdif1
    
    class(dop853_custom),intent(inout) :: dop
    real(dp), intent(in) :: P_cur, P_old
    real(dp), intent(out) :: P_strat
    
    integer, parameter :: n = 1, m = 1
    real(dp) :: x(1)
    real(dp) :: fvec(1)
    real(dp), parameter :: tol = 1.0e-8_dp
    integer :: info
    integer :: iwa(n)
    integer, parameter :: lwa = (n*(3*n+13))/2 + 1
    real(dp) :: wa(lwa)
    
    x(1) = log10(0.5_dp*(P_cur+P_old))
    call lmdif1(fcn, m, n, x, fvec, tol, info, iwa, wa, lwa)
    if (info /= 2) then
      print*,info
      print*,'root solve failed'
      stop
    endif
    P_strat = 10.0_dp**x(1)
    
  contains
    subroutine fcn(m, n, x, fvec, iflag)
      integer, intent(in) :: m
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: fvec(n)
      integer, intent(inout) :: iflag
      real(dp) :: T, P, p_H2O_sat
      P = 10.0_dp**x(1)
      T = dop%contd8(1, P)
      fvec(1) = dop%T_strat - T
    end subroutine
    
  end subroutine
  
  subroutine find_dry_moist_boundary(dop, P_cur, P_old, P_bound)
    use minpack_module, only: lmdif1
    
    class(dop853_custom),intent(inout) :: dop
    real(dp), intent(in) :: P_cur, P_old
    real(dp), intent(out) :: P_bound
      
    integer, parameter :: n = 1, m = 1
    real(dp) :: x(1)
    real(dp) :: fvec(1)
    real(dp), parameter :: tol = 1.0e-8_dp
    integer :: info
    integer :: iwa(n)
    integer, parameter :: lwa = (n*(3*n+13))/2 + 1
    real(dp) :: wa(lwa)
    
    x(1) = log10(0.5_dp*(P_cur+P_old))
    call lmdif1(fcn, m, n, x, fvec, tol, info, iwa, wa, lwa)
    if (info /= 2) then
      print*,info
      print*,'root solve failed'
      stop
    endif
    P_bound = 10.0_dp**x(1)
    
  contains
    subroutine fcn(m, n, x, fvec, iflag)
      integer, intent(in) :: m !! the number of variables.
      integer, intent(in) :: n !! the number of variables.
      real(dp), intent(in) :: x(n) !! independent variable vector
      real(dp), intent(out) :: fvec(n) !! value of function at `x`
      integer, intent(inout) :: iflag !! set to <0 to terminate execution
      real(dp) :: T, P, p_H2O_sat
      P = 10.0_dp**x(1)
      T = dop%contd8(1, P)
      p_H2O_sat = sat_pressure_H2O(T)
      fvec(1) = p_H2O_sat - P*dop%f_H2O_surf
    end subroutine
  end subroutine
  
  !! Dry adiabat
  subroutine solout_dry(self, nr, xold, x, y, irtrn, xout)
    class(dop853_class),intent(inout) :: self
    integer,intent(in)                :: nr
    real(dp),intent(in)               :: xold
    real(dp),intent(in)               :: x
    real(dp),dimension(:),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(dp),intent(out)              :: xout
    
    real(dp) :: p_H2O_sat
    real(dp) :: T_cur, P_cur, P_old
    real(dp) :: PP
    
    P_old = xold
    P_cur = x
    T_cur = y(1)
    
    select type (self)
    class is (dop853_custom)
      
      if (T_cur < self%T_strat) then
        ! entered the stratosphere
        ! do a root solve for the exact P
        ! where we enter the stratosphere
        allocate(self%P_strat)
        call find_stratosphere_boundary(self, P_cur, P_old, self%P_strat)
        self%stopping_reason = 2
        irtrn = -1
        
      elseif (self%T_strat <= T_cur .and. T_cur < T_crit_H2O) then
        p_H2O_sat = sat_pressure_H2O(T_cur)
        
        if (self%f_H2O_surf*P_cur > p_H2O_sat) then
          ! Entered the moist regime. 
          
          ! Do a root solve for the exact P
          ! where we enter the moist regime
          allocate(self%P_bound, self%T_bound)
          call find_dry_moist_boundary(self, P_cur, P_old, self%P_bound)
          ! Temperature at the moist dry boundary
          self%T_bound = self%contd8(1, self%P_bound)
          
          ! Hault integration
          self%stopping_reason = 3
          irtrn = -2
          
        endif
      endif
      
      if (self%P_out(self%j) <= P_old .and. self%P_out(self%j) >= P_cur) then
        if (irtrn == -1) then
          PP = self%P_strat
        elseif (irtrn == -2) then
          PP = self%P_bound
        else
          PP = P_cur
        endif
        
        do while (self%P_out(self%j) >= PP)
          self%T_out(self%j) = self%contd8(1, self%P_out(self%j))
          self%j = self%j + 1
        enddo
      endif
      
    end select
  end subroutine
  
  subroutine dT_dP_dry_dop(self, tn, u, du)
    class(dop853_class), intent(inout) :: self
    real(dp), intent(in) :: tn
    real(dp), intent(in), dimension(:) :: u
    real(dp), intent(out), dimension(:) :: du
    
    select type (self)
    class is (dop853_custom)
      du(1) = dT_dP_dry(tn, u(1), self%f_H2O_surf, self%f_CO2_dry, self%f_N2_dry)
    end select
  
  end subroutine
  
  function dT_dP_dry(P, T, f_H2O, f_CO2_dry, f_N2_dry) result(dTdP)
    use clima_const, only: N_avo, k_boltz
    real(dp), intent(in) :: P !! dynes/cm2
    real(dp), intent(in) :: T !! K
    real(dp), intent(in) :: f_H2O, f_CO2_dry, f_N2_dry
    
    real(dp) :: dTdP
    
    real(dp) :: mubar, cp
    real(dp) :: f_dry, f_CO2, f_N2
    
    ! get CO2 and N2 mixing ratios
    f_dry = 1.0_dp - f_H2O
    f_CO2 = f_CO2_dry*f_dry
    f_N2 = f_N2_dry*f_dry
    
    ! mean molecular weight
    mubar = f_H2O*MU_H2O + f_CO2*MU_CO2 + f_N2*MU_N2
    
    ! heat capacity
    cp = (f_H2O*heat_capacity_H2O(T) + &
          f_CO2*heat_capacity_CO2(T) + &
          f_N2*heat_capacity_N2(T))*(1.0_dp/(mubar*1.0e-3_dp))
    ! convert to erg/(g*K)
    cp = cp*1.0e4_dp
    
    ! adiabat
    dTdP = (N_avo*k_boltz*T)/(P*mubar*cp)
    
  end function
  
  !! Moist adiabat
  
  subroutine solout_moist(self, nr, xold, x, y, irtrn, xout)
    class(dop853_class),intent(inout) :: self
    integer,intent(in)                :: nr
    real(dp),intent(in)               :: xold
    real(dp),intent(in)               :: x
    real(dp),dimension(:),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(dp),intent(out)              :: xout
    
    real(dp) :: T_cur, P_cur, P_old
    real(dp) :: PP
    
    P_old = xold
    P_cur = x
    T_cur = y(1)
    
    select type (self)
    class is (dop853_custom)
      
      if (T_cur < self%T_strat) then
        ! entered the stratosphere
        ! do a root solve for the exact P
        ! where we enter the stratosphere
        allocate(self%P_strat)
        call find_stratosphere_boundary(self, P_cur, P_old, self%P_strat)
        self%stopping_reason = 2
        irtrn = -1
      endif
      
      if (self%P_out(self%j) <= P_old .and. self%P_out(self%j) >= P_cur) then
        if (irtrn == -1) then
          PP = self%P_strat
        else
          PP = P_cur
        endif
        
        do while (self%P_out(self%j) >= PP)
          self%T_out(self%j) = self%contd8(1, self%P_out(self%j))
          self%j = self%j + 1
        enddo
      endif
      
    end select
  end subroutine
  
  subroutine dT_dP_moist_dop(self, tn, u, du)
    class(dop853_class), intent(inout) :: self
    real(dp), intent(in) :: tn
    real(dp), intent(in), dimension(:) :: u
    real(dp), intent(out), dimension(:) :: du
    
    select type (self)
    class is (dop853_custom)
      du(1) = dT_dP_moist(tn, u(1), self%f_CO2_dry, self%f_N2_dry)
    end select
  
  end subroutine
  
  pure function dT_dP_moist(P, T, f_CO2_dry, f_N2_dry) result(dTdP)
    use clima_const, only: N_avo, k_boltz
    real(dp), intent(in) :: P !! dynes/cm2
    real(dp), intent(in) :: T !! K
    real(dp), intent(in) :: f_CO2_dry, f_N2_dry
    
    real(dp) :: dTdP
    
    real(dp) :: mubar, cp
    real(dp) :: f_dry, f_CO2, f_N2
    real(dp) :: f_H2O, mass_frac_H2O
    real(dp) :: P_H2O_sat
    
    ! H2O pressure is determined by saturation
    P_H2O_sat = sat_pressure_H2O(T)
    f_H2O = P_H2O_sat/P
    
    ! get CO2 and N2 mixing ratios
    f_dry = 1.0_dp - f_H2O
    f_CO2 = f_CO2_dry*f_dry
    f_N2 = f_N2_dry*f_dry
    
    ! mean molecular weight
    mubar = f_H2O*MU_H2O + f_CO2*MU_CO2 + f_N2*MU_N2
    
    ! heat capacity
    cp = (f_H2O*heat_capacity_H2O(T) + &
         f_CO2*heat_capacity_CO2(T) + &
         f_N2*heat_capacity_N2(T))*(1.0_dp/(mubar*1.0e-3_dp))
    ! convert to erg/(g*K)
    cp = cp*1.0e4_dp
    
    ! mass fraction of atmosphere that is H2O
    mass_frac_H2O = (MU_H2O/mubar)*f_H2O
    
    ! adiabat
    dTdP = ((N_avo*k_boltz*T)/(P*mubar))* &
            (1.0_dp + (L_H2O*mubar*mass_frac_H2O)/(N_avo*k_boltz*T))/ &
            (cp + (L_H2O**2.0_dp*MU_H2O*mass_frac_H2O)/(N_avo*k_boltz*T**2.0_dp))

  end function
  
  !! Saturation pressure of H2O
  pure function sat_pressure_H2O(T) result(p_H2O_sat)
    real(dp), intent(in) :: T !! Temperature (K)
    real(dp) :: p_H2O_sat !! Saturation pressure (dynes/cm2)
    p_H2O_sat = 1.0e6_dp*exp(5000.0_dp/373.0_dp-5000.0_dp/T)
  end function
  
  !! Heat capacity of H2O, CO2, and N2
  pure function heat_capacity_H2O(T) result(cp)
    use clima_eqns, only: heat_capacity_shomate
    real(dp), intent(in) :: T !! K
    real(dp) :: cp !! J/(mol K)
    
    real(dp), parameter :: coeffs_1(7) = &
          [30.092, 6.832514, 6.793435, -2.53448, 0.082139, -250.881, 223.3967]
    real(dp), parameter :: coeffs_2(7) = &
          [41.96426, 8.622053, -1.49978, 0.098119, -11.15764, -272.1797, 219.7809]
    
    if (T > 0.0_dp .and. T <= 1700.0_dp) then
      cp = heat_capacity_shomate(coeffs_1, T)
    elseif (T > 1700.0_dp .and. T <= 6000.0_dp) then
      cp = heat_capacity_shomate(coeffs_2, T)
    endif
    
  end function
  
  pure function heat_capacity_CO2(T) result(cp)
    use clima_eqns, only: heat_capacity_shomate
    real(dp), intent(in) :: T !! K
    real(dp) :: cp !! J/(mol K)
    
    real(dp), parameter :: coeffs_1(7) = &
          [24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, 228.2431]
    real(dp), parameter :: coeffs_2(7) = &
          [58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186, 263.6125]
    
    if (T > 0.0_dp .and. T <= 1200.0_dp) then
      cp = heat_capacity_shomate(coeffs_1, T)
    elseif (T > 1200.0_dp .and. T <= 6000.0_dp) then
      cp = heat_capacity_shomate(coeffs_2, T)
    endif
    
  end function
  
  pure function heat_capacity_N2(T) result(cp)
    use clima_eqns, only: heat_capacity_shomate
    real(dp), intent(in) :: T !! K
    real(dp) :: cp !! J/(mol K)
    
    real(dp), parameter :: coeffs_1(7) = &
          [26.09, 8.22, -1.98, 0.16, 0.04, -7.99, 221.02]

    if (T > 0.0_dp .and. T <= 6000.0_dp) then
      cp = heat_capacity_shomate(coeffs_1, T)
    endif
    
  end function
  
end module
! end submodule