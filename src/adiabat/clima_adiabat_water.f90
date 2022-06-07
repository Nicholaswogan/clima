submodule(clima_adiabat) clima_adiabat_water
  use dop853_module, only: dop853_class
  implicit none
  
  type :: AdiabatProfileData
    real(dp), pointer :: T_surf
    real(dp) :: P_surf
    real(dp), allocatable :: f_i_surf(:), f_i_dry(:)
    
    ! intent(inout)
    real(dp), pointer :: P(:)
    
    ! intent(out)
    real(dp), pointer :: z(:)
    real(dp), pointer :: T(:)
    real(dp), pointer :: f_i(:,:)
    
    ! work
    integer :: stopping_reason
    real(dp), allocatable :: P_trop, z_trop
    real(dp), allocatable :: P_bound, z_bound, T_bound
    integer :: j
    
    ! error
    character(:), allocatable :: err
    
  end type
  
  ! stopping_reason
  enum, bind(c)
    enumerator :: ReachedPtop, ReachedTropopause, ReachedMoistAdiabat
  end enum
  
  !! latent heat of H2O vaporization/condensation (erg/g)
  real(dp), parameter :: L_H2O = 2307.61404e7_dp
  !! critical point H2O
  real(dp), parameter :: T_crit_H2O = 647.0_dp !! K
  
  type, extends(dop853_class) :: dop853_custom
    type(AdiabatClimateModel), pointer :: c => NULL()
    type(AdiabatProfileData), pointer :: d => NULL()
  end type
  
contains
  
  module subroutine make_profile_a(self, T_surf, P_i_surf, &
                                   P, z, T, f_i, &
                                   err)
    use stdlib_math, only: logspace
    type(AdiabatClimateModel), target, intent(in) :: self
    
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), intent(in) :: P_i_surf(:) !! (ng) dynes/cm2
    
    real(dp), target, intent(out) :: P(:), z(:), T(:) ! (ng)
    real(dp), target, intent(out) :: f_i(:,:) ! (nz,ng)
    character(:), allocatable, intent(out) :: err
    
    type(AdiabatProfileData) :: d
    logical :: moist_start
    real(dp) :: P_H2O_sat_surf, P_H2O_start, P_dry
    integer :: i
    
    real(dp), allocatable :: P_i_surf_(:)
    
    ! check inputs
    if (any(P_i_surf < 0.0_dp)) then
      err = "make_profile: Surface pressures must be positive"
      return
    endif
    if (size(P_i_surf) /= self%sp%ng) then
      err = 'make_profile: Input "P_i_surf" has the wrong shape'
      return
    endif
    if (size(P) /= self%nz+1) then
      err = 'make_profile: Input "P" has the wrong shape'
      return
    endif
    if (size(z) /= self%nz+1) then
      err = 'make_profile: Input "z" has the wrong shape'
      return
    endif
    if (size(T) /= self%nz+1) then
      err = 'make_profile: Input "T" has the wrong shape'
      return
    endif
    if (size(f_i, 1) /= self%nz+1 .or. size(f_i, 2) /= self%sp%ng) then
      err = 'make_profile: Input "f_i" has the wrong shape'
      return
    endif
    
    if (T_surf > T_crit_H2O) then
      ! Above critical point of water.
      ! All H2O must be in the atmosphere.
      moist_start = .false.
      P_H2O_start = P_i_surf(self%LH2O)
    else
      ! Below critical point
      P_H2O_sat_surf = sat_pressure_H2O(T_surf, self%sp%g(self%LH2O)%mass)
      if (P_H2O_sat_surf > P_i_surf(self%LH2O)) then
        ! All water on surface is vaporized
        moist_start = .false.
        P_H2O_start = P_i_surf(self%LH2O)
      else
        ! Water is able to condense at the surface
        moist_start = .true.
        P_H2O_start = P_H2O_sat_surf
      endif
    endif
    
    ! allocate
    allocate(P_i_surf_(self%sp%ng))
    allocate(d%f_i_surf(self%sp%ng))
    allocate(d%f_i_dry(self%sp%ng))
    
    P_i_surf_ = P_i_surf
    P_i_surf_(self%LH2O) = P_H2O_start
    
    d%P_surf = sum(P_i_surf_)
    if (self%P_top > d%P_surf) then
      err = 'make_profile: "P_top" is bigger than the surface pressure'
      return
    endif
    
    ! mixing ratios at the surface
    d%f_i_surf = P_i_surf_/d%P_surf
    P_dry = tiny(1.0_dp)
    do i = 1,self%sp%ng
      if (i /= self%LH2O) then
        P_dry = P_dry + P_i_surf_(i)
      endif
    enddo
    d%f_i_dry = P_i_surf_/P_dry 
    
    ! Make P profile
    P = logspace(log10(d%P_surf),log10(self%P_top),self%nz+1)
    P(1) = d%P_surf
    P(self%nz+1) = self%P_top
    
    ! associate d with inputs
    d%T_surf => T_surf
    d%P => P
    d%z => z
    d%T => T
    d%f_i => f_i
    
    if (moist_start) then
      call integrate_moist_start(self, d, err)
      if (allocated(err)) return
    else
      
      stop 1
      
    endif

  end subroutine
  
  subroutine integrate_moist_start(c, d, err)
    type(AdiabatClimateModel), target, intent(in) :: c
    type(AdiabatProfileData), target, intent(inout) :: d
    
    character(:), allocatable, intent(out) :: err
    
    type(dop853_custom) :: dop
    logical :: status_ok
    integer :: idid, i
    real(dp) :: u(2), Pn
    
    call dop%initialize(fcn=rhs_moist_dop, solout=solout_moist, n=2, &
                        iprint=0, icomp=[1,2], status_ok=status_ok)
    dop%c => c ! associate data
    dop%d => d
    d%j = 2
    call mixing_ratios_moist(c, d, d%P_surf, d%T_surf, d%f_i(1,:))
    d%T(1) = d%T_surf
    d%z(1) = 0.0_dp
    
    d%stopping_reason = ReachedPtop
    if (.not. status_ok) then
      err = 'dop853 initialization failed'
      return
    endif
    
    ! integrate
    u = [d%T_surf, 0.0_dp]
    Pn = d%P_surf
    call dop%integrate(Pn, u, c%P_top, [1.0e-9_dp], [1.0e-6_dp], iout=2, idid=idid)
    
    
    
  end subroutine
  
  subroutine find_tropopause(dop, c, d, P_cur, P_old)
    use minpack_module, only: lmdif1
    
    class(dop853_class),intent(inout) :: dop
    type(AdiabatClimateModel), intent(inout) :: c
    type(AdiabatProfileData), intent(inout) :: d
    real(dp), intent(in) :: P_cur, P_old
    
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
    if (info < 1 .or. info > 4) then
      d%err = 'lmdif1 root solve failed'
    endif
    allocate(d%P_trop)
    d%P_trop = 10.0_dp**x(1)
    
  contains
    subroutine fcn(m, n, x, fvec, iflag)
      integer, intent(in) :: m
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: fvec(n)
      integer, intent(inout) :: iflag
      real(dp) :: T, P
      P = 10.0_dp**x(1)
      T = dop%contd8(1, P)
      fvec(1) = c%T_trop - T
    end subroutine
    
  end subroutine
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Moist adiabat funcitons !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine solout_moist(self, nr, xold, x, y, irtrn, xout)
    class(dop853_class),intent(inout) :: self
    integer,intent(in)                :: nr
    real(dp),intent(in)               :: xold
    real(dp),intent(in)               :: x
    real(dp),dimension(:),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(dp),intent(out)              :: xout
    
    type(AdiabatClimateModel), pointer :: c
    type(AdiabatProfileData), pointer :: d
    
    real(dp) :: T_cur, z_cur, P_cur, P_old
    real(dp) :: PP
    
    P_old = xold
    P_cur = x
    T_cur = y(1)
    z_cur = y(2)
    
    select type (self)
    class is (dop853_custom)
      c => self%c
      d => self%d
    end select
    
    ! Check if the integration needs to be stopped 
    if (T_cur < c%T_trop) then
      ! Hit the tropopause

      ! Solve for the exact pressure of
      ! the tropopause.
      call find_tropopause(self, c, d, P_cur, P_old)
      
      ! altitude of tropopause
      allocate(d%z_trop)
      d%z_trop = self%contd8(2, d%P_trop)
      
      d%stopping_reason = ReachedTropopause
      irtrn = -1
      
    endif
    
    ! save results
    if (d%j <= size(d%P)) then
      if (d%P(d%j) <= P_old .and. d%P(d%j) >= P_cur) then
        if (irtrn == -1) then
          PP = d%P_trop
        else
          PP = P_cur
        endif
        
        do while (d%P(d%j) >= PP)
          d%T(d%j) = self%contd8(1, d%P(d%j))
          d%z(d%j) = self%contd8(2, d%P(d%j))
          call mixing_ratios_moist(c, d, d%P(d%j), d%T(d%j), d%f_i(d%j,:))
          d%j = d%j + 1
          if (d%j > size(d%P)) exit
        enddo
      endif
      
    endif

  end subroutine
  
  subroutine mixing_ratios_moist(c, d, P, T, f_i_layer)
    type(AdiabatClimateModel) :: c
    type(AdiabatProfileData) :: d
    real(dp) :: P, T
    real(dp) :: f_i_layer(:)
    
    real(dp) :: mu_H2O, f_H2O, f_dry
    integer :: i
    
    mu_H2O = c%sp%g(c%LH2O)%mass
    f_H2O = sat_pressure_H2O(T, mu_H2O)/P
    
    f_dry = 1.0_dp - f_H2O
    
    f_i_layer(c%LH2O) = f_H2O
    do i = 1,c%sp%ng
      if (i /= c%LH2O) then
        f_i_layer(i) = f_dry*d%f_i_dry(i)
      endif
    enddo
    print*,f_i_layer
    
  end subroutine
  
  
  subroutine rhs_moist_dop(self, P, u, du)
    use clima_const, only: Rgas
    use clima_eqns, only: heat_capacity_eval, gravity
    class(dop853_class), intent(inout) :: self
    real(dp), intent(in) :: P
    real(dp), intent(in), dimension(:) :: u
    real(dp), intent(out), dimension(:) :: du
    
    type(AdiabatClimateModel), pointer :: c
    type(AdiabatProfileData), pointer :: d
    
    real(dp) :: T, z
    real(dp) :: mubar
    real(dp) :: cp_dry, mu_dry, P_dry, f_dry
    real(dp) :: cp_H2O, mu_H2O, P_H2O, f_H2O
    real(dp) :: dP_dry_dT, dP_H2O_dT
    real(dp) :: dT_dp, dz_dP
    real(dp) :: grav
    
    real(dp) :: cp
    integer :: i
    logical :: found
    
    select type (self)
    class is (dop853_custom)
      c => self%c
      d => self%d
    end select
    
    ! unpack u
    T = u(1)
    z = u(2)
    
    ! Water
    mu_H2O = c%sp%g(c%LH2O)%mass
    P_H2O = sat_pressure_H2O(T, mu_H2O)
    f_H2O = P_H2O/P
    call heat_capacity_eval(c%sp%g(c%LH2O)%thermo, T, found, cp)
    if (.not. found) then
      print*,'not good!'
      stop 1
    endif
    ! convert to erg/(g*K)
    cp_H2O = cp*(1.0_dp/(c%sp%g(c%LH2O)%mass*1.0e-3_dp))*1.0e4_dp
    
    ! Dry atmosphere
    P_dry = P - P_H2O
    f_dry = 1.0_dp - f_H2O
    
    ! mudry and dry heat capacity
    mu_dry = tiny(1.0_dp)**(0.25_dp)
    cp_dry = tiny(1.0_dp)**(0.25_dp)
    do i = 1,c%sp%ng
      if (i /= c%LH2O) then
        mu_dry = mu_dry + d%f_i_dry(i)*c%sp%g(i)%mass
        
        ! J/(mol*K)
        call heat_capacity_eval(c%sp%g(i)%thermo, T, found, cp)
        if (.not. found) then
          print*,'not good!'
          stop 1
        endif
        ! J/(kg*K)
        cp_dry = cp_dry + d%f_i_dry(i)*cp*(1.0_dp/(c%sp%g(i)%mass*1.0e-3_dp))
      endif
    enddo
    ! convert to erg/(g*K)
    cp_dry = cp_dry*1.0e4_dp
           
    dP_dry_dT = (P_dry*mu_dry*cp_dry)/(Rgas*T)*&
                (P_dry + (cp_H2O/cp_dry + (L_H2O*mu_H2O/(Rgas*T) - 1.0_dp)&
                       *(L_H2O/(cp_dry*T)))*((mu_H2O*P_H2O)/(mu_dry)))/ &
                (P_dry + (L_H2O*mu_dry)/(Rgas*T)*(mu_H2O*P_H2O)/(mu_dry))
                
    dP_H2O_dT = (L_H2O*mu_H2O*P_H2O)/(Rgas*T**2.0_dp)
    
    dT_dP = 1.0_dp/(dP_dry_dT + dP_H2O_dT)
    
    ! gravity
    mubar = f_dry*mu_dry + f_H2O*mu_H2O
    grav = gravity(c%planet_radius, c%planet_mass, z)
    dz_dP = -(Rgas*T)/(grav*P*mubar)
    
    ! output
    du = [dT_dP, dz_dP]
    
  end subroutine

  
  !!!!!!!!!!!!!!!!!!!!!!!!
  !!! useful functions !!!
  !!!!!!!!!!!!!!!!!!!!!!!!
  
  pure function altitude_vs_P(P, T, mubar, P0, z0, planet_mass, planet_radius) result(z)
    use clima_const, only: k_boltz, N_avo
    real(dp), intent(in) :: P
    real(dp), intent(in) :: T, mubar, P0, z0, planet_mass, planet_radius
    real(dp) :: z
    real(dp), parameter :: G_grav_cgs = 6.67e-8_dp
    
    z = ((N_avo*k_boltz*T)/(G_grav_cgs*planet_mass*mubar)*log(P/P0) &
        + 1/(planet_radius + z0))**(-1.0_dp) - planet_radius

  end function
  
  function sat_pressure_H2O(T, mu_H2O) result(p_H2O_sat)
    use clima_const, only: Rgas
    real(dp), intent(in) :: T !! Temperature (K)
    real(dp), intent(in) :: mu_H2O !! g/mol
    real(dp) :: p_H2O_sat !! Saturation pressure (dynes/cm2)
    
    p_H2O_sat = 1.0e6_dp*exp((L_H2O*mu_H2O)/(Rgas)*(1.0_dp/373.0_dp - 1.0_dp/T))
    
  end function
  
  
  
end submodule