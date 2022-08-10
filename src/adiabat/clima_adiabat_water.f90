module clima_adiabat_water
  use clima_const, only: dp
  use clima_types, only: Species
  use dop853_module, only: dop853_class
  use futils, only: brent_class
  implicit none

  private
  public :: make_profile_water
  public :: make_column_water
  
  type :: AdiabatProfileData
    ! pointers to input data
    type(Species), pointer :: sp
    integer, pointer :: LH2O
    real(dp), pointer :: planet_mass
    real(dp), pointer :: planet_radius
    real(dp), pointer :: P_top
    real(dp), pointer :: T_trop
    real(dp), pointer :: RH
    ! pointers to input variables
    real(dp), pointer :: T_surf
    
    ! intent(inout)
    real(dp), pointer :: P(:)
    
    ! intent(out)
    real(dp), pointer :: z(:)
    real(dp), pointer :: T(:)
    real(dp), pointer :: f_i(:,:)
    
    ! work
    real(dp) :: P_surf
    real(dp), allocatable :: f_i_surf(:), f_i_dry(:)
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
    type(AdiabatProfileData), pointer :: d => NULL()
  end type
  
contains

  subroutine make_column_water(T_surf, N_i_surf, &
                               sp, nz, LH2O, planet_mass, &
                               planet_radius, P_top, T_trop, RH, &
                               P, z, T, f_i, &
                               err)
    use minpack_module, only: hybrd1
    use clima_eqns, only: gravity
    use clima_const, only: k_boltz, N_avo
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), intent(in) :: N_i_surf(:) !! (ng) dynes/cm2

    type(Species), target, intent(in) :: sp
    integer, intent(in) :: nz
    integer, target, intent(in) :: LH2O
    real(dp), target, intent(in) :: planet_mass, planet_radius
    real(dp), target, intent(in) :: P_top, T_trop, RH
    
    real(dp), target, intent(out) :: P(:), z(:), T(:) ! (ng)
    real(dp), target, intent(out) :: f_i(:,:) ! (nz,ng)
    character(:), allocatable, intent(out) :: err

    integer :: i, j
    real(dp) :: grav
    real(dp), allocatable :: P_av(:), T_av(:), density_av(:), f_i_av(:,:), dz(:)
    real(dp) :: N_i(sp%ng), P_i(sp%ng)
    real(dp) :: P_H2O_ocean, N_H2O_ocean

    integer :: n
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: fvec(:)
    real(dp), parameter :: tol = 1.0e-5_dp
    integer :: info
    integer :: lwa
    real(dp), allocatable :: wa(:)

    n = sp%ng
    allocate(x(n))
    allocate(fvec(n))
    lwa = (n*(3*n+13))/2 + 1
    allocate(wa(lwa))

    allocate(P_av(nz), T_av(nz), density_av(nz), f_i_av(nz,sp%ng), dz(nz))

    ! Initial guess will be crude conversion of moles/cm2 (column) to 
    ! to dynes/cm2 (pressure)
    grav = gravity(planet_radius, planet_mass, 0.0_dp)
    do i = 1,n
      x(i) = log10(max(N_i_surf(i)*sp%g(i)%mass*grav, sqrt(tiny(1.0_dp))))
    enddo

    call hybrd1(fcn, n, x, fvec, tol, info, wa, lwa)
    if (info == 0 .or. info > 1) then
      err = 'hybrd1 root solve failed in make_column_water.'
      return
    elseif (info < 0) then
      err = 'hybrd1 root solve failed in make_column_water: '//err
      return
    endif

    ! call one more time with solution
    call fcn(n, x, fvec, info)

  contains
    subroutine fcn(n_, x_, fvec_, iflag_)
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(n_)
      integer, intent(inout) :: iflag_

      P_i(:) = 10.0_dp**x_(:)
      call make_profile_water(T_surf, P_i, &
                              sp, nz, LH2O, planet_mass, &
                              planet_radius, P_top, T_trop, RH, &
                              P, z, T, f_i, &
                              err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif 

      ! Figure out the "pressure" of the ocean
      P_H2O_ocean = P_i(LH2O) - P(1)*f_i(1,LH2O)
      ! convert to a column (moles/cm2)
      N_H2O_ocean = P_H2O_ocean/(sp%g(LH2O)%mass*grav)

      ! Compute the columns by splitting grid into cells
      do i = 1,nz
        P_av(i) = P(i)
        T_av(i) = T(i)
        dz(i) = z(i+1)-z(i)
        do j = 1,sp%ng
          f_i_av(i,j) = f_i(i,j)
        enddo
      enddo

      density_av = P_av/(k_boltz*T_av)

      do i = 1,sp%ng
        N_i(i) = sum(density_av*f_i_av(:,i)*dz)/N_avo
      enddo
      
      ! Add the ocean to the H2O column in the atmosphere
      N_i(LH2O) = N_i(LH2O) + N_H2O_ocean

      ! residual
      fvec_(:) = N_i - N_i_surf

    end subroutine

  end subroutine
  
  subroutine make_profile_water(T_surf, P_i_surf, &
                                sp, nz, LH2O, planet_mass, &
                                planet_radius, P_top, T_trop, RH, &
                                P, z, T, f_i, &
                                err)
    use futils, only: linspace
    
    real(dp), target, intent(in) :: T_surf !! K
    real(dp), intent(in) :: P_i_surf(:) !! (ng) dynes/cm2

    type(Species), target, intent(in) :: sp
    integer, intent(in) :: nz
    integer, target, intent(in) :: LH2O
    real(dp), target, intent(in) :: planet_mass, planet_radius
    real(dp), target, intent(in) :: P_top, T_trop, RH
    
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
    if (size(P_i_surf) /= sp%ng) then
      err = 'make_profile: Input "P_i_surf" has the wrong shape'
      return
    endif
    if (size(P) /= nz+1) then
      err = 'make_profile: Input "P" has the wrong shape'
      return
    endif
    if (size(z) /= nz+1) then
      err = 'make_profile: Input "z" has the wrong shape'
      return
    endif
    if (size(T) /= nz+1) then
      err = 'make_profile: Input "T" has the wrong shape'
      return
    endif
    if (size(f_i, 1) /= nz+1 .or. size(f_i, 2) /= sp%ng) then
      err = 'make_profile: Input "f_i" has the wrong shape'
      return
    endif
    if (T_surf < T_trop) then
      err = 'T_surf is less than T_trop'
      return
    endif
    
    if (T_surf > T_crit_H2O) then
      ! Above critical point of water.
      ! All H2O must be in the atmosphere.
      moist_start = .false.
      P_H2O_start = P_i_surf(LH2O)
    else
      ! Below critical point
      P_H2O_sat_surf = RH*sat_pressure_H2O(T_surf, sp%g(LH2O)%mass)
      if (P_H2O_sat_surf > P_i_surf(LH2O)) then
        ! All water on surface is vaporized
        moist_start = .false.
        P_H2O_start = P_i_surf(LH2O)
      else
        ! Water is able to condense at the surface
        moist_start = .true.
        P_H2O_start = P_H2O_sat_surf
      endif
    endif
    
    ! allocate
    allocate(P_i_surf_(sp%ng))
    allocate(d%f_i_surf(sp%ng))
    allocate(d%f_i_dry(sp%ng))
    
    P_i_surf_ = P_i_surf
    P_i_surf_(LH2O) = P_H2O_start
    
    d%P_surf = sum(P_i_surf_)
    if (P_top > d%P_surf) then
      err = 'make_profile: "P_top" is bigger than the surface pressure'
      return
    endif
    
    ! mixing ratios at the surface
    d%f_i_surf = P_i_surf_/d%P_surf
    P_dry = tiny(1.0_dp)
    do i = 1,sp%ng
      if (i /= LH2O) then
        P_dry = P_dry + P_i_surf_(i)
      endif
    enddo
    d%f_i_dry = P_i_surf_/P_dry 
    
    ! Make P profile
    call linspace(log10(d%P_surf),log10(P_top),P)
    P(:) = 10.0_dp**P(:)
    P(1) = d%P_surf
    P(nz+1) = P_top
    
    ! associate d with data
    d%sp => sp
    d%LH2O => LH2O
    d%planet_mass => planet_mass
    d%planet_radius => planet_radius
    d%P_top => P_top
    d%T_trop => T_trop
    d%RH => RH
    ! associate d with inputs
    d%T_surf => T_surf
    d%P => P
    d%z => z
    d%T => T
    d%f_i => f_i

    if (moist_start) then
      call integrate_moist_start(d, err)
      if (allocated(err)) return
    else
      call integrate_dry_start(d, err)
      if (allocated(err)) return
    endif
    
  end subroutine
  
  subroutine integrate_moist_start(d, err)
    type(AdiabatProfileData), target, intent(inout) :: d
    
    character(:), allocatable, intent(out) :: err
    
    type(dop853_custom) :: dop
    logical :: status_ok
    integer :: idid, i, j
    real(dp) :: u(2), Pn, mubar
    character(6) :: tmp_char
    
    call dop%initialize(fcn=rhs_moist_dop, solout=solout_moist, n=2, &
                        iprint=0, icomp=[1,2], status_ok=status_ok)
    dop%d => d
    d%j = 2
    call mixing_ratios_moist(d, d%P_surf, d%T_surf, d%f_i(1,:))
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
    call dop%integrate(Pn, u, d%P_top, [1.0e-9_dp], [1.0e-9_dp], iout=2, idid=idid)
    if (allocated(d%err)) then
      err = d%err
      return
    endif
    if (idid < 0) then
      write(tmp_char,'(i6)') idid
      err = 'dop853 integration failed: '//trim(tmp_char)
      return
    endif
    
    if (d%stopping_reason == ReachedPtop) then
      ! Nothing to do.
    elseif (d%stopping_reason == ReachedTropopause) then
      ! Values at tropopause
      d%P(d%j) = d%P_trop
      d%T(d%j) = d%T_trop
      d%z(d%j) = d%z_trop
      
      call mixing_ratios_moist(d, d%P(d%j), d%T(d%j), d%f_i(d%j,:))
      
      ! Values above the tropopause
      do i = d%j+1,size(d%z)
        d%T(i) = d%T_trop
        d%f_i(i,:) = d%f_i(d%j,:)
        mubar = 0.0_dp
        do j = 1,d%sp%ng
          mubar = mubar + d%f_i(d%j,j)*d%sp%g(j)%mass
        enddo
        d%z(i) = altitude_vs_P(d%P(i), d%T_trop, mubar, d%P_trop, &
                 d%z_trop, d%planet_mass, d%planet_radius) 
      enddo
      
    endif
    
  end subroutine
  
  subroutine integrate_dry_start(d, err)
    type(AdiabatProfileData), target, intent(inout) :: d
    
    character(:), allocatable, intent(out) :: err
    
    type(dop853_custom) :: dop
    logical :: status_ok
    integer :: idid, i, j, ind
    real(dp) :: u(2), Pn
    real(dp) :: mubar
    character(6) :: tmp_char
    
    call dop%initialize(fcn=rhs_dry_dop, solout=solout_dry, n=2, &
                        iprint=0, icomp=[1,2], status_ok=status_ok)
    dop%d => d ! associate data
    d%j = 2
    d%T(1) = d%T_surf
    d%z(1) = 0.0_dp
    d%f_i(1,:) = d%f_i_surf(:)
    d%stopping_reason = ReachedPtop
    if (.not. status_ok) then
      err = 'dop853 initialization failed'
      return
    endif
    
    ! integrate
    u = [d%T_surf, 0.0_dp]
    Pn = d%P_surf
    call dop%integrate(Pn, u, d%P_top, [1.0e-9_dp], [1.0e-9_dp], iout=2, idid=idid)
    if (allocated(d%err)) then
      err = d%err
      return
    endif
    if (idid < 0) then
      write(tmp_char,'(i6)') idid
      err = 'dop853 integration failed: '//trim(tmp_char)
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
      d%f_i(d%j,:) = d%f_i_surf(:)
      
      ! Values above the tropopause
      do i = d%j+1,size(d%z)
        d%T(i) = d%T_trop
        d%f_i(i,:) = d%f_i(d%j,:)
        mubar = 0.0_dp
        do j = 1,d%sp%ng
          mubar = mubar + d%f_i(d%j,j)*d%sp%g(j)%mass
        enddo
        d%z(i) = altitude_vs_P(d%P(i), d%T_trop, mubar, d%P_trop, &
                 d%z_trop, d%planet_mass, d%planet_radius) 
      enddo
      

    elseif (d%stopping_reason == ReachedMoistAdiabat) then
      ! ! Do another integration with a moist adiabat
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
      call dop%integrate(Pn, u, d%P_top, [1.0e-9_dp], [1.0e-9_dp], iout=2, idid=idid)
      if (allocated(d%err)) then
        err = d%err
        return
      endif
      if (idid < 0) then
        write(tmp_char,'(i6)') idid
        err = 'dop853 integration failed: '//trim(tmp_char)
        return
      endif
      
      ! In either case, we must save where the dry-moist
      ! transition occured
      ind = minloc(abs(d%P - d%P_bound), 1)
      d%P(ind) = d%P_bound
      d%T(ind) = d%T_bound
      d%z(ind) = d%z_bound
      d%f_i(ind,:) = d%f_i_surf(:)
      
      if (d%stopping_reason == ReachedPtop) then
        ! Nothing more to do
      elseif (d%stopping_reason == ReachedTropopause) then
        ! Values at tropopause
        d%P(d%j) = d%P_trop
        d%T(d%j) = d%T_trop
        d%z(d%j) = d%z_trop
        call mixing_ratios_moist(d, d%P(d%j), d%T(d%j), d%f_i(d%j,:))
        
        ! Values above the tropopause
        do i = d%j+1,size(d%z)
          d%T(i) = d%T_trop
          d%f_i(i,:) = d%f_i(d%j,:)
          mubar = 0.0_dp
          do j = 1,d%sp%ng
            mubar = mubar + d%f_i(d%j,j)*d%sp%g(j)%mass
          enddo
          d%z(i) = altitude_vs_P(d%P(i), d%T_trop, mubar, d%P_trop, &
                   d%z_trop, d%planet_mass, d%planet_radius) 
        enddo
      
      endif
      
    endif
    
  end subroutine
  
  subroutine find_tropopause(dop, d, P_cur, P_old)
    class(dop853_class),intent(inout) :: dop
    type(AdiabatProfileData), intent(inout) :: d
    real(dp), intent(in) :: P_cur, P_old

    real(dp), parameter :: tol = 1.0e-8_dp
    real(dp) :: xzero, fzero
    integer :: info

    type(brent_class) :: brent
    call brent%set_function(fcn)
    call brent%find_zero(P_old, P_cur, tol, xzero, fzero, info)
    if (info /= 0) then
      d%err = 'brent failed in "find_tropopause"'
      return
    endif

    allocate(d%P_trop)
    d%P_trop = xzero
    
  contains
    function fcn(me, x) result(f)
      class(brent_class), intent(inout) :: me
      real(dp), intent(in) :: x
      real(dp) :: f
      f = d%T_trop - dop%contd8(1, x)
    end function
  end subroutine
  
  subroutine find_dry_moist_boundary(dop, d, P_cur, P_old)
    class(dop853_class),intent(inout) :: dop
    type(AdiabatProfileData), intent(inout) :: d
    real(dp), intent(in) :: P_cur, P_old
      
    real(dp), parameter :: tol = 1.0e-8_dp
    real(dp) :: xzero, fzero
    integer :: info
    
    type(brent_class) :: brent
    call brent%set_function(fcn)
    call brent%find_zero(P_old, P_cur, tol, xzero, fzero, info)
    if (info /= 0) then
      d%err = 'brent failed in find_dry_moist_boundary'
      return
    endif

    allocate(d%P_bound)
    d%P_bound = xzero

  contains
    function fcn(me, x) result(f)
      class(brent_class), intent(inout) :: me
      real(dp), intent(in) :: x
      real(dp) :: f
      real(dp) :: T, P, p_H2O_sat
      P = x
      T = dop%contd8(1, P)
      p_H2O_sat = d%RH*sat_pressure_H2O(T, d%sp%g(d%LH2O)%mass)
      f = p_H2O_sat - P*d%f_i_surf(d%LH2O)
    end function
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
    
    type(AdiabatProfileData), pointer :: d
    
    real(dp) :: T_cur, z_cur, P_cur, P_old
    real(dp) :: PP
    real(dp) :: p_H2O_sat
    
    P_old = xold
    P_cur = x
    T_cur = y(1)
    z_cur = y(2)
    
    select type (self)
    class is (dop853_custom)
      d => self%d
    end select

    if (allocated(d%err)) then
      irtrn = -10
      return
    endif
    
    ! Check if the integration needs to be stopped 
    if (T_cur < d%T_trop) then
      ! Hit the tropopause

      ! Solve for the exact pressure of
      ! the tropopause.
      call find_tropopause(self, d, P_cur, P_old)
      if (allocated(d%err)) then
        irtrn = -1
        return
      endif
      
      ! altitude of tropopause
      allocate(d%z_trop)
      d%z_trop = self%contd8(2, d%P_trop)
      
      d%stopping_reason = ReachedTropopause
      irtrn = -1
      
    elseif (d%T_trop <= T_cur .and. T_cur < T_crit_H2O) then
      p_H2O_sat = d%RH*sat_pressure_H2O(T_cur, d%sp%g(d%LH2O)%mass)
      if (d%f_i_surf(d%LH2O)*P_cur > p_H2O_sat) then
        ! Entered the moist adiabat regime.
        
        ! Solve for the exact P where we entered the
        ! moist regime.
        call find_dry_moist_boundary(self, d, P_cur, P_old)
        if (allocated(d%err)) then
          deallocate(d%err)
          d%err = 'Failed to cross the dry-moist boundary.'
          irtrn = -1
          return
        endif
        
        allocate(d%T_bound, d%z_bound)
        d%T_bound = self%contd8(1, d%P_bound)
        d%z_bound = self%contd8(2, d%P_bound)
        
        d%stopping_reason = ReachedMoistAdiabat
        irtrn = -2
  
      endif

    endif
    
    ! save results
    if (d%j <= size(d%P)) then
    
      if (d%P(d%j) <= P_old .and. d%P(d%j) >= P_cur) then
        if (irtrn == -1) then
          PP = d%P_trop
        elseif (irtrn == -2) then
          PP = d%P_bound
        else
          PP = P_cur
        endif
      
        do while (d%P(d%j) >= PP)
          d%T(d%j) = self%contd8(1, d%P(d%j))
          d%z(d%j) = self%contd8(2, d%P(d%j))
          d%f_i(d%j, :) = d%f_i_surf(:)
          d%j = d%j + 1
          if (d%j > size(d%P)) exit
        enddo
      endif
      
    endif
    
  end subroutine
  
  subroutine rhs_dry_dop(self, P, u, du)
    use clima_const, only: Rgas
    use clima_eqns, only: gravity, heat_capacity_eval
    class(dop853_class), intent(inout) :: self
    real(dp), intent(in) :: P
    real(dp), intent(in), dimension(:) :: u
    real(dp), intent(out), dimension(:) :: du
    
    type(AdiabatProfileData), pointer :: d
    
    real(dp) :: T, z
    real(dp) :: mubar, cp, tmp
    logical :: found
    real(dp) :: grav
    real(dp) :: dT_dP, dz_dP
    integer :: i
    
    select type (self)
    class is (dop853_custom)
      d => self%d
    end select
    
    ! unpack u
    T = u(1)
    z = u(2)
    
    ! mean molecular weight
    cp = 0.0_dp
    mubar = 0.0_dp
    do i = 1,d%sp%ng
        mubar = mubar + d%f_i_surf(i)*d%sp%g(i)%mass
        
        ! J/(mol*K)
        call heat_capacity_eval(d%sp%g(i)%thermo, T, found, tmp)
        if (.not. found) then
          d%err = "Failed to compute heat capacity"
          return
        endif
        ! J/(kg*K)
        cp = cp + d%f_i_surf(i)*tmp*(1.0_dp/(d%sp%g(i)%mass*1.0e-3_dp))
    enddo
    ! convert to erg/(g*K)
    cp = cp*1.0e4_dp
    
    ! How T changes with pressure
    dT_dP = (Rgas*T)/(P*mubar*cp)
    
    ! How altitude changes with pressure
    grav = gravity(d%planet_radius, d%planet_mass, z)
    dz_dP = -(Rgas*T)/(grav*P*mubar)
    
    du = [dT_dP, dz_dP]
    
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
    
    type(AdiabatProfileData), pointer :: d
    
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

    if (allocated(d%err)) then
      irtrn = -10
      return
    endif
    
    ! Check if the integration needs to be stopped 
    if (T_cur < d%T_trop) then
      ! Hit the tropopause

      ! Solve for the exact pressure of
      ! the tropopause.
      call find_tropopause(self, d, P_cur, P_old)
      if (allocated(d%err)) then
        irtrn = -1
        return
      endif
      
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
          call mixing_ratios_moist(d, d%P(d%j), d%T(d%j), d%f_i(d%j,:))
          d%j = d%j + 1
          if (d%j > size(d%P)) exit
        enddo
      endif
      
    endif

  end subroutine
  
  subroutine mixing_ratios_moist(d, P, T, f_i_layer)
    type(AdiabatProfileData) :: d
    real(dp) :: P, T
    real(dp) :: f_i_layer(:)
    
    real(dp) :: f_H2O, f_dry
    integer :: i
    
    f_H2O = d%RH*sat_pressure_H2O(T, d%sp%g(d%LH2O)%mass)/P
    f_dry = 1.0_dp - f_H2O
    
    f_i_layer(d%LH2O) = f_H2O
    do i = 1,d%sp%ng
      if (i /= d%LH2O) then
        f_i_layer(i) = f_dry*d%f_i_dry(i)
      endif
    enddo
    
  end subroutine
  
  subroutine rhs_moist_dop(self, P, u, du)
    use clima_const, only: Rgas
    use clima_eqns, only: heat_capacity_eval, gravity
    class(dop853_class), intent(inout) :: self
    real(dp), intent(in) :: P
    real(dp), intent(in), dimension(:) :: u
    real(dp), intent(out), dimension(:) :: du
    
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
      d => self%d
    end select
    
    ! unpack u
    T = u(1)
    z = u(2)
    
    ! Water
    mu_H2O = d%sp%g(d%LH2O)%mass
    P_H2O = d%RH*sat_pressure_H2O(T, mu_H2O)
    f_H2O = P_H2O/P
    call heat_capacity_eval(d%sp%g(d%LH2O)%thermo, T, found, cp)
    if (.not. found) then
      d%err = "Failed to compute heat capacity"
      return
    endif
    ! convert to erg/(g*K)
    cp_H2O = cp*(1.0_dp/(d%sp%g(d%LH2O)%mass*1.0e-3_dp))*1.0e4_dp
    
    ! Dry atmosphere
    P_dry = P - P_H2O
    f_dry = 1.0_dp - f_H2O
    
    ! mudry and dry heat capacity
    mu_dry = tiny(0.0_dp)
    cp_dry = tiny(0.0_dp)
    do i = 1,d%sp%ng
      if (i /= d%LH2O) then
        mu_dry = mu_dry + d%f_i_dry(i)*d%sp%g(i)%mass
        
        ! J/(mol*K)
        call heat_capacity_eval(d%sp%g(i)%thermo, T, found, cp)
        if (.not. found) then
          d%err = "Failed to compute heat capacity"
          return
        endif
        ! J/(kg*K)
        cp_dry = cp_dry + d%f_i_dry(i)*cp*(1.0_dp/(d%sp%g(i)%mass*1.0e-3_dp))
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
    grav = gravity(d%planet_radius, d%planet_mass, z)
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
  
end module