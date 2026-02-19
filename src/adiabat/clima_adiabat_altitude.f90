submodule(clima_adiabat) clima_adiabat_altitude
  use dop853_module, only: dop853_class
  implicit none

  type :: AltitudeIntegrationData
    type(linear_interp_1d) :: T_interp
    type(linear_interp_1d) :: mubar_interp
    real(dp), pointer :: P_e(:) => NULL()
    real(dp), pointer :: z_e(:) => NULL()
    real(dp) :: planet_mass
    real(dp) :: planet_radius
    integer :: j
  end type

  type, extends(dop853_class) :: dop853_alt
    type(AltitudeIntegrationData), pointer :: d => NULL()
  end type

contains

  module subroutine AdiabatClimate_compute_altitude(self, err)
    use clima_eqns, only: gravity
    use clima_eqns_water, only: Rgas_cgs => Rgas
    class(AdiabatClimate), intent(inout) :: self
    character(:), allocatable, intent(out) :: err

    type(AltitudeIntegrationData), target :: d
    type(dop853_alt) :: dop
    logical :: status_ok
    integer :: i, j, idid, ierr
    real(dp) :: Pn
    real(dp), allocatable :: T_grid(:), mubar_grid(:), P_grid(:)
    real(dp), allocatable :: P_interp(:), T_interp(:), mubar_interp(:)
    real(dp), allocatable :: u(:)
    integer, allocatable :: icomp(:)
    real(dp), allocatable, target :: P_e(:), z_e(:)
    character(6) :: tmp_char

    if (self%nz < 2) then
      err = 'compute_altitude requires nz >= 2'
      return
    endif
    if (self%P_surf <= 0.0_dp .or. self%P_top <= 0.0_dp) then
      err = 'compute_altitude requires positive P_surf and P_top'
      return
    endif

    allocate(P_e(2*self%nz+1), z_e(2*self%nz+1))
    allocate(T_grid(self%nz+1), mubar_grid(self%nz+1), P_grid(self%nz+1))
    allocate(P_interp(self%nz+1), T_interp(self%nz+1), mubar_interp(self%nz+1))
    allocate(u(1), icomp(1))

    P_e(1) = self%P_surf
    do i = 1,self%nz
      P_e(2*i) = self%P(i)
      if (i < self%nz) P_e(2*i+1) = sqrt(self%P(i)*self%P(i+1))
    enddo
    P_e(2*self%nz+1) = self%P_top

    if (minval(P_e(1:2*self%nz) - P_e(2:2*self%nz+1)) <= 10.0_dp*spacing(P_e(1))) then
      err = 'compute_altitude: pressure grid is not strictly decreasing.'
      return
    endif

    P_grid(1) = self%P_surf
    T_grid(1) = self%T_surf
    mubar_grid(1) = 0.0_dp
    do j = 1,self%sp%ng
      mubar_grid(1) = mubar_grid(1) + self%f_i_surf(j)*self%sp%g(j)%mass
    enddo
    do i = 1,self%nz
      P_grid(i+1) = self%P(i)
      T_grid(i+1) = self%T(i)
      mubar_grid(i+1) = 0.0_dp
      do j = 1,self%sp%ng
        mubar_grid(i+1) = mubar_grid(i+1) + self%f_i(i,j)*self%sp%g(j)%mass
      enddo
    enddo

    do i = 1,self%nz+1
      P_interp(i) = log10(P_grid(self%nz+2-i))
      T_interp(i) = T_grid(self%nz+2-i)
      mubar_interp(i) = mubar_grid(self%nz+2-i)
    enddo

    call d%T_interp%initialize(P_interp, T_interp, ierr)
    if (ierr /= 0) then
      err = 'compute_altitude: failed to initialize temperature interpolator.'
      return
    endif
    call d%mubar_interp%initialize(P_interp, mubar_interp, ierr)
    if (ierr /= 0) then
      err = 'compute_altitude: failed to initialize mubar interpolator.'
      return
    endif

    d%P_e => P_e
    d%z_e => z_e
    d%planet_mass = self%planet_mass
    d%planet_radius = self%planet_radius
    d%j = 2

    z_e(1) = 0.0_dp
    u(1) = 0.0_dp
    icomp = [1]
    Pn = P_e(1)

    call dop%initialize(fcn=right_hand_side_altitude_dop, solout=solout_altitude_dop, n=1, &
                        iprint=0, icomp=icomp, status_ok=status_ok)
    if (.not. status_ok) then
      err = 'dop853 initialization failed in compute_altitude'
      return
    endif
    dop%d => d

    call dop%integrate(Pn, u, P_e(2*self%nz), [self%rtol], [self%atol], iout=2, idid=idid)
    if (idid < 0) then
      write(tmp_char,'(i6)') idid
      err = 'dop853 integration failed in compute_altitude: '//trim(tmp_char)
      return
    endif
    z_e(2*self%nz) = u(1)
    if (d%j <= 2*self%nz) then
      err = 'compute_altitude: dense output stalled.'
      return
    endif

    z_e(2*self%nz+1) = z_e(2*self%nz) + (z_e(2*self%nz) - z_e(2*self%nz-1))

    do i = 1,self%nz
      self%z(i) = z_e(2*i)
      self%dz(i) = z_e(2*i+1) - z_e(2*i-1)
    enddo

  end subroutine

  subroutine right_hand_side_altitude_dop(self_, P, u, du)
    use clima_eqns, only: gravity
    use clima_eqns_water, only: Rgas_cgs => Rgas
    class(dop853_class), intent(inout) :: self_
    real(dp), intent(in) :: P
    real(dp), intent(in) :: u(:)
    real(dp), intent(out) :: du(:)

    real(dp) :: T, mubar, grav

    select type (self_)
    class is (dop853_alt)
      call self_%d%T_interp%evaluate(log10(P), T)
      call self_%d%mubar_interp%evaluate(log10(P), mubar)
      grav = gravity(self_%d%planet_radius, self_%d%planet_mass, u(1))
      du(1) = -(Rgas_cgs*T)/(grav*P*mubar)
    end select

  end subroutine

  subroutine solout_altitude_dop(self_, nr, xold, x, y, irtrn, xout)
    class(dop853_class),intent(inout) :: self_
    integer,intent(in)                :: nr
    real(dp),intent(in)               :: xold
    real(dp),intent(in)               :: x
    real(dp),dimension(:),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(dp),intent(out)              :: xout

    type(AltitudeIntegrationData), pointer :: d
    real(dp) :: P_old, P_cur, P_tol

    P_old = xold
    P_cur = x
    xout = x

    select type (self_)
    class is (dop853_alt)
      d => self_%d
    end select

    if (d%j > size(d%P_e)-1) return

    P_tol = max(1.0e-10_dp*max(abs(P_old), abs(P_cur), abs(d%P_e(d%j))), &
                10.0_dp*spacing(max(P_old, P_cur)))
    if (d%P_e(d%j) <= P_old .and. d%P_e(d%j) >= P_cur - P_tol) then
      do while (d%P_e(d%j) >= P_cur - P_tol)
        d%z_e(d%j) = self_%contd8(1, d%P_e(d%j))
        d%j = d%j + 1
        if (d%j > size(d%P_e)-1) exit
      enddo
    endif

  end subroutine

end submodule
