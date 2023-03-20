module clima_rc_adiabat
  use clima_const, only: dp
  use clima_types, only: Species
  use linear_interpolation_module, only: linear_interp_1d
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
    integer, pointer :: cond_z_inds(:,:) !! (2,size(cond_inds))
    real(dp), pointer :: z(:) !! cm
    real(dp), pointer :: dz(:) !! cm
    real(dp), pointer :: T(:) !! K

    ! work space
    real(dp), allocatable :: log10_P(:) !! P grid but in log10 space.
    integer :: ncond
    integer :: ndry !! ng - ncond
    integer, allocatable :: dry_inds(:) !! (ndry)
    type(linear_interp_1d), allocatable :: f_interps(:) !! (ndry)
    integer, allocatable :: sp_type(:) !! (ncond)
    real(dp), allocatable :: f_i_cur(:) !! (ng)
    real(dp), allocatable :: P_i_cur(:) !! (ng)

    real(dp), allocatable :: P_integ(:)
    real(dp), allocatable :: z_integ(:)
    real(dp), allocatable :: f_i_integ(:,:)

    !> This helps us propogate error messages outside of the integration
    character(:), allocatable :: err

  end type

  ! sp_type
  enum, bind(c)
    enumerator :: DrySpeciesType, CondensingSpeciesType
  end enum

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
    integer, target, intent(out) :: cond_z_inds(:,:)
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
    if (size(cond_z_inds,1) /= 2 .and. size(cond_z_inds,2) /= size(cond_P)) then
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
    allocate(d%f_i_cur(sp%ng))
    allocate(d%P_i_cur(sp%ng))

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
    allocate(d%f_i_integ(2*d%nz+1,d%sp%ng))
    d%P_integ(1) = P_surf
    do i = 1,d%nz
      j = 2*i
      d%P_integ(j) = d%P(i)
      d%P_integ(j+1) = 10.0_dp**(d%log10_P(i) - d%log10_delta_P(i)/2.0_dp)
    enddo
    d%z_integ(1) = 0.0_dp
    d%f_i_integ(1,:) = d%f_i_cur(:)

    call integrate(d, err)
    if (allocated(err)) return

  end subroutine

  subroutine integrate(d, err)
    type(RCAdiabatData), intent(inout) :: d
    character(:), allocatable, intent(out) :: err

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
  
end module