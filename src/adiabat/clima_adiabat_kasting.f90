module clima_adiabat_kasting
  use clima_const, only: dp
  use dop853_module, only: dop853_class
  implicit none
  private
  
  public :: make_profile
  
  type :: KastingProfileData
    ! intent(in)
    real(dp), pointer :: T_surf
    real(dp) :: P_surf
    real(dp) :: f_H2O_surf, f_CO2_dry, f_N2_dry
    integer, pointer :: nz
    real(dp), pointer :: planet_mass, planet_radius
    real(dp), pointer :: P_top
    real(dp), pointer :: T_trop
    
    ! intent(inout)
    real(dp), pointer :: P(:)
    
    ! intent(out)
    real(dp), pointer :: z(:)
    real(dp), pointer :: T(:)
    real(dp), pointer :: f_H2O(:), f_CO2(:), f_N2(:)
    
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
  
  !! molar masses (g/mol)
  real(dp), parameter :: MU_H2O = 18.0_dp
  real(dp), parameter :: MU_CO2 = 44.0_dp
  real(dp), parameter :: MU_N2 = 28.0_dp
  !! latent heat of H2O vaporization/condensation (erg/g)
  real(dp), parameter :: L_H2O = 2264e7_dp
  !! critical point H2O
  real(dp), parameter :: T_crit_H2O = 647.0_dp !! K
  
  type, extends(dop853_class) :: dop853_custom
    type(KastingProfileData), pointer :: d => NULL()
  end type
  
contains
  
  subroutine make_profile(T_surf, P_H2O_surf, P_CO2_surf, P_N2_surf, nz, &
                          planet_mass, planet_radius, P_top, T_trop, &
                          P, z, T, f_H2O, f_CO2, f_N2, surface_liquid_H2O, &
                          err)
    use stdlib_math, only: logspace
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), target, intent(in) :: P_H2O_surf, P_CO2_surf, P_N2_surf !! dynes/cm2
    integer, target, intent(in) :: nz
    real(dp), target, intent(in) :: planet_mass, planet_radius
    real(dp), target, intent(in) :: P_top, T_trop
    
    real(dp), target, intent(out) :: P(nz), z(nz), T(nz)
    real(dp), target, intent(out) :: f_H2O(nz), f_CO2(nz), f_N2(nz)
    logical, target, intent(out) :: surface_liquid_H2O
    character(:), allocatable, intent(out) :: err
    
    type(KastingProfileData) :: d
    
    real(dp) :: P_H2O_sat_surf, P_H2O_start
    
    ! check inputs
    if (P_H2O_surf < 0.0_dp .or. P_CO2_surf < 0.0_dp .or. P_N2_surf < 0.0_dp) then
      err = "make_profile: Surface pressures must be positive"
      return
    endif
    if (nz < 5) then
      err = 'make_profile: "nz" must be bigger than 4'
      return
    endif
    if (P_top < 0.0_dp .or. T_trop < 0.0_dp) then
      err = "make_profile: P_top and T_strat must be positive"
      return
    endif
    
    if (T_surf > T_crit_H2O) then
      ! Above critical point of water.
      ! All H2O must be in the atmosphere.
      surface_liquid_H2O = .false.
      P_H2O_start = P_H2O_surf
    else
      ! Below critical point
      P_H2O_sat_surf = sat_pressure_H2O(T_surf)
      if (P_H2O_sat_surf > P_H2O_surf) then
        ! All water on surface is vaporized
        surface_liquid_H2O = .false.
        P_H2O_start = P_H2O_surf
      else
        ! Water is able to condense at the surface
        surface_liquid_H2O = .true.
        P_H2O_start = P_H2O_sat_surf
      endif
    endif
    
    d%P_surf = P_H2O_start + P_CO2_surf + P_N2_surf
    if (P_top > d%P_surf) then
      err = 'make_profile: "P_top" is bigger than the surface pressure'
      return
    endif
    
    ! mixing ratios at the surface
    d%f_H2O_surf = P_H2O_start/d%P_surf
    d%f_CO2_dry = P_CO2_surf/max(P_CO2_surf + P_N2_surf, tiny(1.0_dp))
    d%f_N2_dry = P_N2_surf/max(P_CO2_surf + P_N2_surf, tiny(1.0_dp))
    
    ! Make P profile
    P = logspace(log10(d%P_surf),log10(P_top),nz)
    P(1) = d%P_surf
    P(nz) = P_top
    
    ! associate d with inputs
    d%T_surf => T_surf
    d%planet_mass => planet_mass
    d%planet_radius => planet_radius
    d%T_trop => T_trop
    d%P_top => P_top
    d%nz => nz
    d%P => P
    d%z => z
    d%T => T
    d%f_H2O => f_H2O
    d%f_CO2 => f_CO2
    d%f_N2 => f_N2
    
    ! begin a new scope
    if (.not. surface_liquid_H2O) then
      call integrate_dry_start(d, err)
      if (allocated(err)) return
    else
      call integrate_moist_start(d, err)
      if (allocated(err)) return
    endif
              
  end subroutine
  
  subroutine integrate_moist_start(d, err)
    type(KastingProfileData), target, intent(inout) :: d
    character(:), allocatable, intent(out) :: err
    
    type(dop853_custom) :: dop
    logical :: status_ok
    integer :: idid, i, ind
    real(dp) :: u(2), Pn
    real(dp) :: mubar
    
    call dop%initialize(fcn=rhs_moist_dop, solout=solout_moist, n=2, &
                        iprint=0, icomp=[1,2], status_ok=status_ok)
    dop%d => d ! associate data
    d%j = 2
    d%T(1) = d%T_surf
    d%z(1) = 0.0_dp
    call mixing_ratios_moist(d%P_surf, d%T_surf, &
                             d%f_CO2_dry, d%f_N2_dry, &
                             d%f_H2O(1), d%f_CO2(1), d%f_N2(1))
    d%stopping_reason = ReachedPtop
    if (.not. status_ok) then
      err = 'dop853 initialization failed'
      return
    endif
    
    ! integrate
    u = [d%T_surf, 0.0_dp]
    Pn = d%P_surf
    call dop%integrate(Pn, u, d%P_top, [1.0e-6_dp], [1.0e-6_dp], iout=2, idid=idid)
    if (allocated(d%err)) then
      err = d%err
      return
    endif
    if (idid < 0) then
      err = 'dop853 integration failed'
      return
    endif
    
    if (d%stopping_reason == ReachedPtop) then
      ! Nothing to do.
    elseif (d%stopping_reason == ReachedTropopause) then
      ! We must compute results at and above the tropopause
      
      ! Values at tropopause
      d%P(d%j) = d%P_trop
      d%T(d%j) = d%T_trop
      d%z(d%j) = d%z_trop
      call mixing_ratios_moist(d%P(d%j), d%T(d%j), &
                               d%f_CO2_dry, d%f_N2_dry, &
                               d%f_H2O(d%j), d%f_CO2(d%j), d%f_N2(d%j))
      
      ! Values above the tropopause
      if (size(d%T) > d%j) then
        d%T(d%j+1:) = d%T_trop
        d%f_H2O(d%j+1:) = d%f_H2O(d%j)
        d%f_CO2(d%j+1:) = d%f_CO2(d%j)
        d%f_N2(d%j+1:) = d%f_N2(d%j)
        mubar = d%f_H2O(d%j)*MU_H2O + d%f_CO2(d%j)*MU_CO2 + d%f_N2(d%j)*MU_N2
        do i = d%j+1,size(d%z)
          d%z(i) = altitude_vs_P(d%P(i), d%T_trop, mubar, d%P_trop, &
                   d%z_trop, d%planet_mass, d%planet_radius) 
        enddo
      endif
    endif
    
  end subroutine
  
  subroutine integrate_dry_start(d, err)
    type(KastingProfileData), target, intent(inout) :: d
    character(:), allocatable, intent(out) :: err
    
    type(dop853_custom) :: dop
    logical :: status_ok
    integer :: idid, i, ind
    real(dp) :: u(2), Pn
    real(dp) :: mubar
    
    call dop%initialize(fcn=rhs_dry_dop, solout=solout_dry, n=2, &
                        iprint=0, icomp=[1,2], status_ok=status_ok)
    dop%d => d ! associate data
    d%j = 2
    d%T(1) = d%T_surf
    d%z(1) = 0.0_dp
    d%f_H2O(1) = d%f_H2O_surf
    call mixing_ratios_dry(d%f_H2O_surf, d%f_CO2_dry, d%f_N2_dry, &
                           d%f_CO2(1), d%f_N2(1))
    d%stopping_reason = ReachedPtop
    if (.not. status_ok) then
      err = 'dop853 initialization failed'
      return
    endif
    
    ! integrate
    u = [d%T_surf, 0.0_dp]
    Pn = d%P_surf
    call dop%integrate(Pn, u, d%P_top, [1.0e-6_dp], [1.0e-6_dp], iout=2, idid=idid)
    if (allocated(d%err)) then
      err = d%err
      return
    endif
    if (idid < 0) then
      err = 'dop853 integration failed'
      return
    endif
    
    if (d%stopping_reason == ReachedPtop) then
      ! Nothing to do.
    elseif (d%stopping_reason == ReachedTropopause) then
      ! We must compute results at and above the tropopause
      
      ! Values at tropopause
      d%P(d%j) = d%P_trop
      d%T(d%j) = d%T_trop
      d%z(d%j) = d%z_trop
      d%f_H2O(d%j) = d%f_H2O_surf
      call mixing_ratios_dry(d%f_H2O_surf, d%f_CO2_dry, d%f_N2_dry, &
                             d%f_CO2(d%j), d%f_N2(d%j))
      
      ! Values above the tropopause
      if (size(d%T) > d%j) then
        d%T(d%j+1:) = d%T_trop
        d%f_H2O(d%j+1:) = d%f_H2O(d%j)
        d%f_CO2(d%j+1:) = d%f_CO2(d%j)
        d%f_N2(d%j+1:) = d%f_N2(d%j)
        mubar = d%f_H2O(d%j)*MU_H2O + d%f_CO2(d%j)*MU_CO2 + d%f_N2(d%j)*MU_N2
        do i = d%j+1,size(d%z)
          d%z(i) = altitude_vs_P(d%P(i), d%T_trop, mubar, d%P_trop, &
                   d%z_trop, d%planet_mass, d%planet_radius) 
        enddo
      endif

    elseif (d%stopping_reason == ReachedMoistAdiabat) then
      ! Do another integration with a moist adiabat
      call dop%initialize(fcn=rhs_moist_dop, solout=solout_moist, n=2, &
                          iprint=0, icomp=[1,2], status_ok=status_ok)
      dop%d => d ! associate data
      d%stopping_reason = ReachedPtop
      if (.not. status_ok) then
        err = 'dop853 initialization failed'
        return
      endif
      
      u = [d%T_bound, d%z_bound]
      Pn = d%P_bound
      call dop%integrate(Pn, u, d%P_top, [1.0e-6_dp], [1.0e-6_dp], iout=2, idid=idid)
      if (allocated(d%err)) then
        err = d%err
        return
      endif
      if (idid < 0) then
        err = 'dop853 integration failed'
        return
      endif
      
      ! In either case, we must save where the dry-moist
      ! transition occured
      ind = minloc(abs(d%P - d%P_bound), 1)
      d%P(ind) = d%P_bound
      d%T(ind) = d%T_bound
      d%z(ind) = d%z_bound
      d%f_H2O(ind) = d%f_H2O_surf
      call mixing_ratios_dry(d%f_H2O_surf, d%f_CO2_dry, d%f_N2_dry, &
                             d%f_CO2(ind), d%f_N2(ind))
      
      if (d%stopping_reason == ReachedPtop) then
        ! Nothing more to do
      elseif (d%stopping_reason == ReachedTropopause) then
        ! Values at tropopause
        d%P(d%j) = d%P_trop
        d%T(d%j) = d%T_trop
        d%z(d%j) = d%z_trop
        call mixing_ratios_moist(d%P_trop, d%T_trop, &
                                 d%f_CO2_dry, d%f_N2_dry, &
                                 d%f_H2O(d%j), d%f_CO2(d%j), d%f_N2(d%j))
        ! Values above the tropopause
        if (size(d%T) > d%j) then
          d%T(d%j+1:) = d%T_trop
          d%f_H2O(d%j+1:) = d%f_H2O(d%j)
          d%f_CO2(d%j+1:) = d%f_CO2(d%j)
          d%f_N2(d%j+1:) = d%f_N2(d%j)
          mubar = d%f_H2O(d%j)*MU_H2O + d%f_CO2(d%j)*MU_CO2 + d%f_N2(d%j)*MU_N2
          do i = d%j+1,size(d%z)
            d%z(i) = altitude_vs_P(d%P(i), d%T_trop, mubar, d%P_trop, &
                  d%z_trop, d%planet_mass, d%planet_radius) 
          enddo
        endif
      
      endif
      
    endif

  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!
  !!! Boundary finders !!!
  !!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine find_tropopause(dop, d, P_cur, P_old)
    use minpack_module, only: lmdif1
    
    class(dop853_class),intent(inout) :: dop
    type(KastingProfileData), intent(inout) :: d
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
    if (info /= 2) then
      d%err = 'lmdif1 root solve failed'
      return
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
      fvec(1) = d%T_trop - T
    end subroutine
    
  end subroutine
  
  subroutine find_dry_moist_boundary(dop, d, P_cur, P_old)
    use minpack_module, only: lmdif1
    
    class(dop853_class),intent(inout) :: dop
    type(KastingProfileData), intent(inout) :: d
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
    if (info /= 2) then
      d%err = 'lmdif1 root solve failed'
      return
    endif
    allocate(d%P_bound)
    d%P_bound = 10.0_dp**x(1)
    
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
      fvec(1) = p_H2O_sat - P*d%f_H2O_surf
    end subroutine
  end subroutine
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Dry adiabat funcitons !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine solout_dry(self, nr, xold, x, y, irtrn, xout)
    class(dop853_class),intent(inout) :: self
    integer,intent(in)                :: nr
    real(dp),intent(in)               :: xold
    real(dp),intent(in)               :: x
    real(dp),dimension(:),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(dp),intent(out)              :: xout
    
    type(KastingProfileData), pointer :: d
    
    real(dp) :: T_cur, z_cur, P_cur, P_old
    real(dp) :: p_H2O_sat, PP
    
    P_old = xold
    P_cur = x
    T_cur = y(1)
    z_cur = y(2)
    
    select type (self)
    class is (dop853_custom)
      d => self%d
    end select
    
    ! Check if the integration needs to be stopped 
    if (T_cur < d%T_trop) then
      ! Hit the tropopause

      ! Solve for the exact pressure of
      ! the tropopause.
      call find_tropopause(self, d, P_cur, P_old)
      
      ! altitude of tropopause
      allocate(d%z_trop)
      d%z_trop = self%contd8(2, d%P_trop)
      
      d%stopping_reason = ReachedTropopause
      irtrn = -1
      
    elseif (d%T_trop <= T_cur .and. T_cur < T_crit_H2O) then
      p_H2O_sat = sat_pressure_H2O(T_cur)
      if (d%f_H2O_surf*P_cur > p_H2O_sat) then
        ! Entered the moist adiabat regime.
        
        ! Solve for the exact P where we entered the
        ! moist regime.
        call find_dry_moist_boundary(self, d, P_cur, P_old)
        
        allocate(d%T_bound, d%z_bound)
        d%T_bound = self%contd8(1, d%P_bound)
        d%z_bound = self%contd8(2, d%P_bound)
        
        d%stopping_reason = ReachedMoistAdiabat
        irtrn = -2
      endif

    endif
    
    ! save results
    if (d%P(d%j) <= P_old .and. d%P(d%j) >= P_cur) then
      if (irtrn == -1) then
        PP = d%P_trop
      elseif (irtrn == -2) then
        PP = d%P_bound
      else
        PP = P_cur
      endif
      
      do while (d%P(d%j) >= PP .and. d%j <= size(d%P))
        d%T(d%j) = self%contd8(1, d%P(d%j))
        d%z(d%j) = self%contd8(2, d%P(d%j))
        d%f_H2O(d%j) = d%f_H2O_surf
        call mixing_ratios_dry(d%f_H2O_surf, d%f_CO2_dry, d%f_N2_dry, &
                               d%f_CO2(d%j), d%f_N2(d%j))
        d%j = d%j + 1
      enddo
    endif

  end subroutine
  
  subroutine rhs_dry_dop(self, tn, u, du)
    class(dop853_class), intent(inout) :: self
    real(dp), intent(in) :: tn
    real(dp), intent(in), dimension(:) :: u
    real(dp), intent(out), dimension(:) :: du
    
    select type (self)
    class is (dop853_custom)
      du(1:2) = rhs_dry(tn, u, self%d%f_H2O_surf, self%d%f_CO2_dry, &
                        self%d%f_N2_dry, self%d%planet_mass, self%d%planet_radius)
    end select

  end subroutine
  
  pure subroutine mixing_ratios_dry(f_H2O_surf, f_CO2_dry, f_N2_dry, f_CO2, f_N2)
    real(dp), intent(in) :: f_H2O_surf, f_CO2_dry, f_N2_dry
    real(dp), intent(out) :: f_CO2, f_N2
    real(dp) :: f_dry
    f_dry = 1.0_dp - f_H2O_surf
    f_CO2 = f_CO2_dry*f_dry
    f_N2 = f_N2_dry*f_dry
  end subroutine
  
  pure function rhs_dry(P, u, f_H2O_surf, f_CO2_dry, f_N2_dry, planet_mass, planet_radius) result(rhs)
    use clima_const, only: N_avo, k_boltz
    use clima_eqns, only: gravity
    real(dp), intent(in) :: P !! dynes/cm2
    real(dp), intent(in) :: u(2) !! K
    real(dp), intent(in) :: f_H2O_surf, f_CO2_dry, f_N2_dry
    real(dp), intent(in) :: planet_mass, planet_radius
    
    real(dp) :: rhs(2)
    
    real(dp) :: T, z
    real(dp) :: mubar, cp
    real(dp) :: f_CO2, f_N2, grav
    real(dp) :: dTdP, dzdP
    
    ! unpack u
    T = u(1)
    z = u(2)
    
    ! get CO2 and N2 mixing ratios
    call mixing_ratios_dry(f_H2O_surf, f_CO2_dry, f_N2_dry, f_CO2, f_N2)
    
    ! mean molecular weight
    mubar = f_H2O_surf*MU_H2O + f_CO2*MU_CO2 + f_N2*MU_N2
    
    ! heat capacity
    cp = (f_H2O_surf*heat_capacity_H2O(T) + &
          f_CO2*heat_capacity_CO2(T) + &
          f_N2*heat_capacity_N2(T))*(1.0_dp/(mubar*1.0e-3_dp))
    ! convert to erg/(g*K)
    cp = cp*1.0e4_dp
    
    ! How T changes with pressure
    dTdP = (N_avo*k_boltz*T)/(P*mubar*cp)
    
    ! How altitude changes with pressure
    grav = gravity(planet_radius, planet_mass, z)
    dzdP = -(N_avo*k_boltz*T)/(grav*P*mubar)
    
    rhs = [dTdP, dzdP]
    
  end function
  
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
    
    type(KastingProfileData), pointer :: d
    
    real(dp) :: T_cur, z_cur, P_cur, P_old
    real(dp) :: PP
    
    P_old = xold
    P_cur = x
    T_cur = y(1)
    z_cur = y(2)
    
    select type (self)
    class is (dop853_custom)
      d => self%d
    end select
    
    ! Check if the integration needs to be stopped 
    if (T_cur < d%T_trop) then
      ! Hit the tropopause

      ! Solve for the exact pressure of
      ! the tropopause.
      call find_tropopause(self, d, P_cur, P_old)
      
      ! altitude of tropopause
      allocate(d%z_trop)
      d%z_trop = self%contd8(2, d%P_trop)
      
      d%stopping_reason = ReachedTropopause
      irtrn = -1
      
    endif
    
    ! save results
    if (d%P(d%j) <= P_old .and. d%P(d%j) >= P_cur) then
      if (irtrn == -1) then
        PP = d%P_trop
      else
        PP = P_cur
      endif
      
      do while (d%P(d%j) >= PP .and. d%j <= size(d%P))
        d%T(d%j) = self%contd8(1, d%P(d%j))
        d%z(d%j) = self%contd8(2, d%P(d%j))
        call mixing_ratios_moist(d%P(d%j), d%T(d%j), &
                                 d%f_CO2_dry, d%f_N2_dry, &
                                 d%f_H2O(d%j), d%f_CO2(d%j), d%f_N2(d%j))
        d%j = d%j + 1
      enddo
    endif
    
    
  end subroutine
  
  subroutine rhs_moist_dop(self, tn, u, du)
    class(dop853_class), intent(inout) :: self
    real(dp), intent(in) :: tn
    real(dp), intent(in), dimension(:) :: u
    real(dp), intent(out), dimension(:) :: du
    
    select type (self)
    class is (dop853_custom)
      du(1:2) = rhs_moist(tn, u, self%d%f_CO2_dry, self%d%f_N2_dry, &
                          self%d%planet_mass, self%d%planet_radius)
    end select
  
  end subroutine
  
  pure subroutine mixing_ratios_moist(P, T, f_CO2_dry, f_N2_dry, f_H2O, f_CO2, f_N2)
    real(dp), intent(in) :: P, T
    real(dp), intent(in) :: f_CO2_dry, f_N2_dry
    real(dp), intent(out) :: f_H2O, f_CO2, f_N2
    
    real(dp) :: P_H2O_sat, f_dry
    
    ! H2O pressure is determined by saturation
    P_H2O_sat = sat_pressure_H2O(T)
    f_H2O = P_H2O_sat/P
    
    ! get CO2 and N2 mixing ratios
    f_dry = 1.0_dp - f_H2O
    f_CO2 = f_CO2_dry*f_dry
    f_N2 = f_N2_dry*f_dry
    
  end subroutine
  
  pure function rhs_moist(P, u, f_CO2_dry, f_N2_dry, planet_mass, planet_radius) result(rhs)
    use clima_const, only: N_avo, k_boltz
    use clima_eqns, only: gravity
    real(dp), intent(in) :: P !! dynes/cm2
    real(dp), intent(in) :: u(2) !! K
    real(dp), intent(in) :: f_CO2_dry, f_N2_dry
    real(dp), intent(in) :: planet_mass, planet_radius
    
    real(dp) :: rhs(2)
    
    real(dp) :: T, z
    real(dp) :: mubar, cp
    real(dp) :: f_CO2, f_N2
    real(dp) :: f_H2O, mass_frac_H2O, grav
    real(dp) :: dTdP, dzdP
    
    ! unpack u
    T = u(1)
    z = u(2)
  
    call mixing_ratios_moist(P, T, f_CO2_dry, f_N2_dry, f_H2O, f_CO2, f_N2)
    
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
    
    ! How T changes with pressure
    dTdP = ((N_avo*k_boltz*T)/(P*mubar))* &
            (1.0_dp + (L_H2O*mubar*mass_frac_H2O)/(N_avo*k_boltz*T))/ &
            (cp + (L_H2O**2.0_dp*MU_H2O*mass_frac_H2O)/(N_avo*k_boltz*T**2.0_dp))
    
    ! How altitude changes with pressure
    grav = gravity(planet_radius, planet_mass, z)
    dzdP = -(N_avo*k_boltz*T)/(grav*P*mubar)
    
    rhs = [dTdP, dzdP]
    
  end function
  
  !!!!!!!!!!!!!!!!!!!!!!!!
  !!! useful functions !!!
  !!!!!!!!!!!!!!!!!!!!!!!!
  
  !! analytical solution for z(P) for constant mubar and T
  pure function altitude_vs_P(P, T, mubar, P0, z0, planet_mass, planet_radius) result(z)
    use clima_const, only: k_boltz, N_avo
    real(dp), intent(in) :: P
    real(dp), intent(in) :: T, mubar, P0, z0, planet_mass, planet_radius
    real(dp) :: z
    real(dp), parameter :: G_grav_cgs = 6.67e-8_dp
    
    z = ((N_avo*k_boltz*T)/(G_grav_cgs*planet_mass*mubar)*log(P/P0) &
        + 1/(planet_radius + z0))**(-1.0_dp) - planet_radius

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