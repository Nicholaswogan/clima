module clima_rc_adiabat
  use clima_const, only: dp
  use clima_types, only: Species
  use linear_interpolation_module, only: linear_interp_1d
  use dop853_module, only: dop853_class
  implicit none
  private

  public :: make_profile_rc

  type :: RCAdiabatData
    ! Input
    type(Species), pointer :: sp
    real(dp) :: planet_mass
    real(dp) :: planet_radius
    real(dp) :: P_surf !! dynes/cm^2
    real(dp) :: T_surf !! K
    real(dp) :: T_trop
    integer :: nz
    real(dp), pointer :: P(:) !! dynes/cm^2
    real(dp), pointer :: log10_delta_P(:) !! dynes/cm^2
    real(dp), pointer :: RH(:) !! RH of condensibles
    real(dp), pointer :: cond_P(:) !! dynes/cm^2
    integer, pointer :: cond_inds(:) !! indexes of condensibles
    integer :: bg_gas_ind !! index of bg gas
    ! Input/Output
    real(dp), pointer :: f_i(:,:)
    ! Output
    integer, pointer :: cond_z_inds(:) !! (size(cond_inds))
    real(dp), pointer :: z(:) !! cm
    real(dp), pointer :: dz(:) !! cm
    real(dp), pointer :: T(:) !! K

    integer :: stopping_reason !! Reason why we stop integration (see enum below)
    !> Index of gas that reached saturation if stopping_reason == ReachedGasSaturation
    integer :: ind_gas_reached_saturation
    logical :: isothermal = .false.
    logical :: search_for_roots = .true.

    ! work space
    real(dp), allocatable :: log10_P(:) !! P grid but in log10 space.
    integer :: ncond
    integer :: ndry !! ng - ncond
    integer, allocatable :: dry_inds(:) !! (ndry)
    type(linear_interp_1d), allocatable :: f_interps(:) !! (ndry)
    integer, allocatable :: sp_type(:) !! (ncond)
    real(dp), allocatable :: L_i_cur(:) !! (ncond)
    real(dp), allocatable :: P_at_start_cond(:) !! (ncond)
    real(dp), allocatable :: f_i_cur(:) !! (ng)
    real(dp), allocatable :: P_i_cur(:) !! (ng)
    real(dp), allocatable :: cp_i_cur(:) !! (ng)

    real(dp), allocatable :: gout(:) !! ncond+1
    real(dp), allocatable :: gout_old(:) !! ncond+1
    real(dp), allocatable :: gout_tmp(:) !! ncond+1
    logical, allocatable :: root_found(:) !! ncond+1
    real(dp), allocatable :: P_roots(:) !! ncond+1
    !> When we find a root, and exit integration, these
    !> guys will give use the root
    real(dp) :: P_root 
    real(dp) :: u_root(2)
    
    real(dp), allocatable :: P_integ(:)
    real(dp), allocatable :: T_integ(:)
    real(dp), allocatable :: z_integ(:)
    real(dp), allocatable :: f_i_integ(:,:)

    !> index for saving T, z and mixing ratios
    integer :: j

    !> This helps us propogate error messages outside of the integration
    character(:), allocatable :: err

  end type

  ! stopping_reason
  enum, bind(c)
    enumerator :: ReachedPtop, ReachedTropopause, ReachedGasSaturation
  end enum

  ! sp_type
  enum, bind(c)
    enumerator :: DrySpeciesType, CondensingSpeciesType
  end enum

  type, extends(dop853_class) :: dop853_custom
    type(RCAdiabatData), pointer :: d => NULL()
  end type

contains

  subroutine make_profile_rc(sp, planet_mass, planet_radius, &
                             P_surf, T_surf, T_trop, nz, P, log10_delta_P, &
                             RH, cond_P, cond_inds, bg_gas_ind, &
                             f_i, cond_z_inds, z, dz, T, err)
    type(Species), target, intent(inout) :: sp
    real(dp), target, intent(in) :: planet_mass, planet_radius
    real(dp), target, intent(in) :: P_surf, T_surf, T_trop
    integer, intent(in) :: nz
    real(dp), target, intent(in) :: P(:), log10_delta_P(:) !! 
    real(dp), target, intent(in) :: RH(:), cond_P(:)
    integer, target, intent(in) :: cond_inds(:), bg_gas_ind
    real(dp), target, intent(inout) :: f_i(:,:)
    integer, target, intent(out) :: cond_z_inds(:)
    real(dp), target, intent(out) :: z(:), dz(:), T(:)
    character(:), allocatable, intent(out) :: err

    type(RCAdiabatData) :: d
    integer :: i, j
    real(dp) :: P_sat

    ! check inputs
    if (P_surf < 0.0_dp) then
      err = 'make_profile_rc: Input "P_surf" is less than 0'
      return
    endif
    if (T_surf < T_trop) then
      err = 'make_profile_rc: Input "T_surf" is less than input "T_trop"'
      return
    endif
    if (T_trop < 0.0_dp) then
      err = 'make_profile_rc: Input "T_trop" is less than 0'
      return
    endif
    if (size(P) /= nz) then
      err = 'make_profile_rc: Input "P" has the wrong shape'
      return
    endif
    if (size(log10_delta_P) /= nz) then
      err = 'make_profile_rc: Input "log10_delta_P" has the wrong shape'
      return
    endif
    if (size(RH) /= size(cond_P)) then
      err = 'make_profile_rc: Input "RH" has the wrong shape'
      return
    endif
    if (size(cond_inds) /= size(cond_P)) then
      err = 'make_profile_rc: Input "cond_inds" has the wrong shape'
      return
    endif
    do i = 1,size(cond_inds)
      if (.not.allocated(sp%g(cond_inds(i))%sat)) then
        err = 'make_profile_rc: Input "cond_inds" indicates species without staturation data.'
        return
      endif
    enddo
    if (bg_gas_ind < 1 .or. bg_gas_ind > sp%ng) then
      err = 'make_profile_rc: Background gas index is outside of the range of species indices.'
      return
    endif
    if (any(bg_gas_ind == cond_inds)) then
      err = 'make_profile_rc: The background gas cannot be a condensing gas.'
      return
    endif
    if (size(f_i,1) /= nz .and. size(f_i,2) /= sp%ng) then
      err = 'make_profile_rc: Input "f_i" has the wrong shape'
      return
    endif
    if (size(cond_z_inds) /= size(cond_P)) then
      err = 'make_profile_rc: Input "cond_z_inds" has the wrong shape'
      return
    endif
    if (size(z) /= nz) then
      err = 'make_profile_rc: Input "z" has the wrong shape'
      return
    endif
    if (size(dz) /= nz) then
      err = 'make_profile_rc: Input "dz" has the wrong shape'
      return
    endif
    if (size(T) /= nz) then
      err = 'make_profile_rc: Input "T" has the wrong shape'
      return
    endif

    ! associate
    ! in
    d%sp => sp
    d%planet_mass = planet_mass
    d%planet_radius = planet_radius
    d%P_surf = P_surf
    d%T_surf = T_surf
    d%T_trop = T_trop
    d%nz = nz
    d%P => P
    d%log10_delta_P => log10_delta_P
    d%RH => RH
    d%cond_P => cond_P
    d%cond_inds => cond_inds
    d%bg_gas_ind = bg_gas_ind
    ! in/out
    d%f_i => f_i
    ! out
    d%cond_z_inds => cond_z_inds
    d%z => z
    d%dz => dz
    d%T => T

    allocate(d%log10_P(nz))
    d%log10_P = log10(d%P)
    d%ncond = size(cond_P)
    d%ndry = sp%ng - d%ncond - 1 ! - 1 for bg gas
    allocate(d%dry_inds(d%ndry))
    allocate(d%f_interps(d%ndry))
    allocate(d%sp_type(d%ncond))
    allocate(d%L_i_cur(d%ncond))
    allocate(d%P_at_start_cond(d%ncond))
    allocate(d%f_i_cur(sp%ng))
    allocate(d%P_i_cur(sp%ng))
    allocate(d%cp_i_cur(sp%ng))
    allocate(d%gout(d%ncond+1))
    allocate(d%gout_old(d%ncond+1))
    allocate(d%gout_tmp(d%ncond+1))
    allocate(d%root_found(d%ncond+1))
    allocate(d%P_roots(d%ncond+1))

    ! Figure out dry indexes. All species that can not condense, and that are not
    ! the background gas are "dry" species.
    j = 1
    do i = 1,sp%ng
      if (.not.any(i == d%cond_inds) .and. .not.(i == d%bg_gas_ind)) then
        d%dry_inds(j) = i
        j = j + 1
      endif
    enddo

    ! Initialize mixing ratio interpolators
    call initialize_mix_interps(d, err)
    if (allocated(err)) return

    ! Dry species mixing ratios and pressures
    do i = 1,d%ndry
      j = d%dry_inds(i)
      d%f_i_cur(j) = d%f_i(1,j)
      d%P_i_cur(j) = d%f_i(1,j)*P_surf
    enddo

    d%P_at_start_cond = -1.0_dp
    ! Check if condensation is happening at surface
    do i = 1,d%ncond
      j = d%cond_inds(i)
      P_sat = huge(1.0_dp)
      if (T_surf < d%sp%g(j)%sat%T_critical) then
        P_sat = d%RH(i)*d%sp%g(j)%sat%sat_pressure(T_surf)
      endif
      if (cond_P(i) > P_sat) then
        d%P_i_cur(j) = P_sat
        d%f_i_cur(j) = P_sat/P_surf
        d%sp_type(i) = CondensingSpeciesType
        d%P_at_start_cond(i) = P_surf
      else
        d%P_i_cur(j) = cond_P(i)
        d%f_i_cur(j) = cond_P(i)/P_surf
        d%sp_type(i) = DrySpeciesType
      endif
    enddo

    ! Background gas. Check if pressure is too big.
    d%f_i_cur(d%bg_gas_ind) = 0.0_dp
    d%f_i_cur(d%bg_gas_ind) = 1.0_dp - sum(d%f_i_cur)
    if (d%f_i_cur(d%bg_gas_ind) < 0.0_dp) then
      err = 'make_profile_rc: The surface pressure is exceeding the input surface pressure'
      return
    endif
    d%P_i_cur(d%bg_gas_ind) = d%f_i_cur(d%bg_gas_ind)*P_surf

    ! Grid of integration. Evaluate at edges and centers of grid cells.
    allocate(d%P_integ(2*d%nz+1))
    allocate(d%z_integ(2*d%nz+1))
    allocate(d%T_integ(2*d%nz+1))
    allocate(d%f_i_integ(2*d%nz+1,d%sp%ng))
    d%P_integ(1) = P_surf
    do i = 1,d%nz
      j = 2*i
      d%P_integ(j) = d%P(i)
      d%P_integ(j+1) = 10.0_dp**(d%log10_P(i) - d%log10_delta_P(i)/2.0_dp)
    enddo

    call integrate(d, err)
    if (allocated(err)) return

  end subroutine

  subroutine initialize_mix_interps(d, err)
    type(RCAdiabatData), intent(inout) :: d
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: tmp_log10_P(:), tmp_log10_f_i(:)
    integer :: i, j, ierr

    allocate(tmp_log10_P(d%nz+2))
    tmp_log10_P(1) = log10(d%P_surf)
    tmp_log10_P(2:d%nz+1) = d%log10_P(:)
    tmp_log10_P(d%nz+2) = d%log10_P(d%nz) - 0.5_dp*d%log10_delta_P(d%nz)
    tmp_log10_P(:) = tmp_log10_P(d%nz+2:1:-1)

    allocate(tmp_log10_f_i(d%nz+2))
    do i = 1,d%ndry
      j = d%dry_inds(i)

      tmp_log10_f_i(1) = log10(d%f_i(1,j))
      tmp_log10_f_i(2:d%nz+1) = log10(d%f_i(:,j))
      tmp_log10_f_i(2:d%nz+2) = log10(d%f_i(d%nz,j))
      tmp_log10_f_i(:) = tmp_log10_f_i(d%nz+2:1:-1)
      call d%f_interps(i)%initialize(tmp_log10_P, tmp_log10_f_i, ierr)
      if (ierr /= 0) then
        err = 'Failed to initialize mixing ratio interpolators'
        return
      endif
    enddo

  end subroutine

  subroutine integrate(d, err)
    type(RCAdiabatData), target, intent(inout) :: d
    character(:), allocatable, intent(out) :: err

    type(dop853_custom) :: dop
    logical :: status_ok
    integer :: idid, i, j
    real(dp) :: Pn, u(2)
    character(6) :: tmp_char

    call dop%initialize(fcn=right_hand_side_dop, solout=solout_dop, n=2, &
                        iprint=0, icomp=[1,2], status_ok=status_ok)
    dop%d => d
    d%j = 2
    d%f_i_integ(1,:) = d%f_i_cur(:)
    d%T_integ(1) = d%T_surf
    d%z_integ(1) = 0.0_dp
    d%stopping_reason = ReachedPtop
    if (.not. status_ok) then
      err = 'dop853 initialization failed'
      return
    endif

    u = [d%T_surf, 0.0_dp]
    Pn = d%P_surf
    do
      call dop%integrate(Pn, u, d%P_integ(size(d%P_integ)), [1.0e-9_dp], [1.0e-9_dp], iout=2, idid=idid)
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
        exit
      elseif (d%stopping_reason == ReachedTropopause) then
        Pn = d%P_root
        u = d%u_root
        d%sp_type(:) = DrySpeciesType
        d%isothermal = .true.
        d%search_for_roots = .false.
        d%stopping_reason = ReachedPtop
      elseif (d%stopping_reason == ReachedGasSaturation) then
        Pn = d%P_root
        u = d%u_root
        d%sp_type(d%ind_gas_reached_saturation) = CondensingSpeciesType
        d%P_at_start_cond(d%ind_gas_reached_saturation) = Pn
        d%stopping_reason = ReachedPtop
      endif

    enddo

    ! Fill in output
    do i = 1,d%nz
      j = 2*i
      d%f_i(i,:) = d%f_i_integ(j,:)
      d%z(i) = d%z_integ(j)
      d%dz(i) = d%z_integ(j+1) - d%z_integ(j-1)
      d%T(i) = d%T_integ(j)
    enddo

    ! Indexes where gases begin saturation
    d%cond_z_inds = -1
    do j = 1,d%ncond
      if (d%P_at_start_cond(j) > 0.0_dp) then
        do i = 1,d%nz
          if (d%P(i) < d%P_at_start_cond(j)) then
            d%cond_z_inds(j) = i
            exit
          endif
        enddo
      endif
    enddo

  end subroutine

  subroutine solout_dop(self, nr, xold, x, y, irtrn, xout)
    class(dop853_class),intent(inout) :: self
    integer,intent(in)                :: nr
    real(dp),intent(in)               :: xold
    real(dp),intent(in)               :: x
    real(dp),dimension(:),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(dp),intent(out)              :: xout

    type(RCAdiabatData), pointer :: d
    real(dp) :: T_cur, z_cur, P_cur, P_old
    real(dp) :: PP, f_dry
    integer :: i

    P_old = xold
    P_cur = x
    T_cur = y(1)
    z_cur = y(2)

    select type (self)
    class is (dop853_custom)
      d => self%d
    end select

    if (allocated(d%err)) then
      ! if there is an error (from RHS function)
      ! we return immediately
      irtrn = -10
      return
    endif

    if (d%search_for_roots) then
    call root_fcn(d, P_cur, T_cur, d%gout)
    if (allocated(d%err)) then
      irtrn = -10
      return
    endif
    if (nr == 1) then
      d%gout_old(:) = d%gout(:)
    endif
    do i = 1,size(d%gout_old)
      if ((d%gout_old(i) < 0 .and. d%gout(i) > 0) .or. &
          (d%gout_old(i) > 0 .and. d%gout(i) < 0)) then
        d%root_found(i) = .true. 
        call find_root(self, d, P_old, P_cur, i, d%P_roots(i))
        if (allocated(d%err)) then
          irtrn = -10
          return
        endif
      else
        d%root_found(i) = .false. 
        d%P_roots(i) = 0.0_dp
      endif
    enddo

    if (any(d%root_found)) then; block
      integer :: ind_root
      ! lets use the highest pressure root.
      ! this makes sure we don't skip an interesting root
      irtrn = -1
      ind_root = maxloc(d%P_roots,1)
      d%P_root = d%P_roots(ind_root)
      d%u_root(1) = self%contd8(1, d%P_root)
      d%u_root(2) = self%contd8(2, d%P_root)

      if (ind_root == 1) then
        ! Tropopause
        d%stopping_reason = ReachedTropopause
      elseif (ind_root >= 2 .and. ind_root <= d%sp%ng+1) then
        ! Some species is supersaturated
        d%stopping_reason = ReachedGasSaturation
        d%ind_gas_reached_saturation = ind_root - 1
      endif

    endblock; endif

    d%gout_old(:) = d%gout(:) ! save for next step
    endif

    ! save the results
    if (d%j <= size(d%P_integ)) then
      if (irtrn == -1) then
        PP = d%P_root
      else
        PP = P_cur
      endif

      if (d%P_integ(d%j) <= P_old .and. d%P_integ(d%j) >= PP) then
        do while (d%P_integ(d%j) >= PP)
          d%T_integ(d%j) = self%contd8(1, d%P_integ(d%j))
          d%z_integ(d%j) = self%contd8(2, d%P_integ(d%j))
          call mixing_ratios(d, d%P_integ(d%j), d%T_integ(d%j), d%f_i_integ(d%j,:), f_dry, d%err)
          if (allocated(d%err)) then
            irtrn = -10
            return
          endif
          d%j = d%j + 1
          if (d%j > size(d%P_integ)) exit
        enddo
      endif

    endif

  end subroutine

  subroutine find_root(dop, d, P_old, P_cur, ind, P_root)
    use futils, only: brent_class
    type(dop853_class), intent(inout) :: dop
    type(RCAdiabatData), intent(inout) :: d
    real(dp), intent(in) :: P_old
    real(dp), intent(in) :: P_cur
    integer, intent(in) :: ind
    real(dp), intent(out) :: P_root

    real(dp), parameter :: tol = 1.0e-8_dp
    real(dp) :: fzero
    integer :: info

    type(brent_class) :: b

    call b%set_function(fcn)
    call b%find_zero(P_old, P_cur, tol, P_root, fzero, info)
    if (info /= 0) then
      d%err = 'brent failed in "find_root"'
      return
    endif

  contains
    function fcn(me_, x_) result(f_)
      class(brent_class), intent(inout) :: me_
      real(dp), intent(in) :: x_
      real(dp) :: f_
      call root_fcn(d, x_, dop%contd8(1, x_), d%gout_tmp)
      f_ = d%gout_tmp(ind)
    end function
  end subroutine

  subroutine root_fcn(d, P, T, gout)
    type(RCAdiabatData), intent(inout) :: d
    real(dp), intent(in) :: P
    real(dp), intent(in) :: T
    real(dp), intent(out) :: gout(:) 

    real(dp) :: f_dry, P_sat
    integer :: i, j

    call mixing_ratios(d, P, T, d%f_i_cur, f_dry, d%err)
    if (allocated(d%err)) return
    d%P_i_cur = d%f_i_cur*P

    do i = 1,d%ncond
      j = d%cond_inds(i)
      P_sat = huge(1.0_dp)
      if (T < d%sp%g(j)%sat%T_critical) then
        P_sat = d%RH(i)*d%sp%g(j)%sat%sat_pressure(T)
      endif

      if (d%sp_type(i) == CondensingSpeciesType) then
        gout(i+1) = 1.0
      elseif (d%sp_type(i) == DrySpeciesType) then
        gout(i+1) = P_sat - d%P_i_cur(i)
      endif
    enddo

    ! also stop at tropopause
    gout(1) = T - d%T_trop

  end subroutine

  subroutine right_hand_side_dop(self, P, u, du)
    class(dop853_class), intent(inout) :: self
    real(dp), intent(in) :: P
    real(dp), intent(in), dimension(:) :: u
    real(dp), intent(out), dimension(:) :: du
    select type (self)
    class is (dop853_custom)
      call right_hand_side(self%d, P, u, du)
    end select
  end subroutine

  subroutine mixing_ratios(d, P, T, f_i_layer, f_dry, err)
    type(RCAdiabatData), intent(inout) :: d
    real(dp), intent(in) :: P, T
    real(dp), intent(out) :: f_i_layer(:), f_dry
    character(:), allocatable, intent(out) :: err

    real(dp) :: log10P, log10_f_i
    integer :: i, j

    log10P = log10(P)

    ! Dry species
    f_dry = 0.0_dp
    do i = 1,d%ndry
      j = d%dry_inds(i)
      call d%f_interps(i)%evaluate(log10P, log10_f_i)
      f_i_layer(j) = 10.0_dp**log10_f_i
      f_dry = f_dry + f_i_layer(j)
    enddo

    ! Species that can condense
    do i = 1,d%ncond
      j = d%cond_inds(i)
      if (d%sp_type(i) == CondensingSpeciesType) then
        f_i_layer(j) = d%RH(i)*d%sp%g(j)%sat%sat_pressure(T)/P
      elseif (d%sp_type(i) == DrySpeciesType) then
        f_i_layer(j) = d%f_i_cur(j) ! Same mixing ratio.
        f_dry = f_dry + f_i_layer(j)
      endif
    enddo

    ! Background gas
    f_i_layer(d%bg_gas_ind) = 0.0_dp
    f_i_layer(d%bg_gas_ind) = 1.0_dp - sum(f_i_layer)
    if (f_i_layer(d%bg_gas_ind) < 0.0_dp) then
      err = 'The background gas is negative.'
      return
    endif
    f_dry = f_dry + f_i_layer(d%bg_gas_ind)

  end subroutine

  subroutine right_hand_side(d, P, u, du)
    use clima_eqns, only: heat_capacity_eval, gravity
    use clima_eqns_water, only: Rgas_cgs => Rgas
    type(RCAdiabatData), intent(inout) :: d
    real(dp), intent(in) :: P
    real(dp), intent(in) :: u(:)
    real(dp), intent(out) :: du(:)

    real(dp) :: T, z
    real(dp) :: f_dry, cp_dry
    real(dp) :: cp, mubar, L
    real(dp) :: first_sumation, second_sumation, Beta_i
    real(dp) :: grav
    real(dp) :: dlnT_dlnP, dT_dP, dz_dP
    logical :: found
    integer :: i, j
    real(dp), parameter :: Rgas_si = Rgas_cgs/1.0e7_dp ! ideal gas constant in SI units (J/(mol*K))

    T = u(1)
    z = u(2)

    call mixing_ratios(d, P, T, d%f_i_cur, f_dry, d%err)
    if (allocated(d%err)) return

    mubar = 0.0_dp
    do i = 1,d%sp%ng
      mubar = mubar + d%f_i_cur(i)*d%sp%g(i)%mass
    enddo

    if (d%isothermal) then
    dT_dP = 0.0_dp
    else
    ! Latent heat
    do i = 1,d%ncond
      j = d%cond_inds(i)
      if (d%sp_type(i) == CondensingSpeciesType) then
        L = d%sp%g(j)%sat%latent_heat(T) ! erg/g
        L = L*d%sp%g(j)%mass*1.0e-7_dp ! convert to J/mol
        d%L_i_cur(i) = L
      endif
    enddo

    ! Heat capacity
    do i = 1,d%sp%ng
      ! heat capacity in J/(mol*K)
      call heat_capacity_eval(d%sp%g(i)%thermo, T, found, cp)
      if (.not. found) then
        d%err = "Failed to compute heat capacity"
        return
      endif
      d%cp_i_cur(i) = cp
    enddo

    ! dry heat capacity
    cp_dry = tiny(0.0_dp)
    do i = 1,d%ndry
      j = d%dry_inds(i)
      cp_dry = cp_dry + d%cp_i_cur(j)*(d%f_i_cur(j)/f_dry)
    enddo
    do i = 1,d%ncond
      j = d%cond_inds(i)
      if (d%sp_type(i) == DrySpeciesType) then
        cp_dry = cp_dry + d%cp_i_cur(j)*(d%f_i_cur(j)/f_dry)
      endif
    enddo
    cp_dry = cp_dry + d%cp_i_cur(d%bg_gas_ind)*(d%f_i_cur(d%bg_gas_ind)/f_dry)

    
    ! Two sums in Equation 1 of Graham et al. (2021)
    first_sumation = 0.0_dp
    second_sumation = 0.0_dp
    do i = 1,d%ncond
      j = d%cond_inds(i)
      if (d%sp_type(i) == CondensingSpeciesType) then
        Beta_i = d%L_i_cur(i)/(Rgas_si*T)
        first_sumation = first_sumation + &
          d%f_i_cur(j)*(d%cp_i_cur(j) - Rgas_si*Beta_i + Rgas_si*Beta_i**2.0_dp)

        second_sumation = second_sumation + &
          Beta_i*d%f_i_cur(j)
      endif
    enddo

    ! Equation 1 in Graham et al. (2021), except simplified to assume no condensate present.
    ! Units are SI
    dlnT_dlnP = 1.0_dp/(f_dry*((cp_dry*f_dry + first_sumation)/(Rgas_si*(f_dry + second_sumation))) + second_sumation)
    ! Convert to dT/dP. Here we introduce CGS units
    ! P = [dynes/cm2]
    ! T = [K]
    dT_dP = dlnT_dlnP*(T/P)
    endif

    ! rate of change of altitude
    grav = gravity(d%planet_radius, d%planet_mass, z)
    dz_dP = -(Rgas_cgs*T)/(grav*P*mubar)

    du(:) = [dT_dP, dz_dP]

  end subroutine
  
end module