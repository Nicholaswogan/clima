submodule(clima_rc) clima_rc_rhs
  implicit none

contains

  !> Computes dT/dt from radiative effects.
  module subroutine RadiativeConvectiveClimate_rhs(self, T, dT_dt, err)
    use clima_eqns, only: heat_capacity_eval, gravity
    use clima_rc_adiabat, only: make_profile_rc
    use clima_const, only: k_boltz
    class(RadiativeConvectiveClimate), intent(inout) :: self
    real(dp), intent(in) :: T(:)
    real(dp), intent(out) :: dT_dt(:)
    character(:), allocatable :: err

    integer :: i, j
    logical :: tropopause_was_reached
    real(dp) :: cp_tmp
    logical :: found
    real(dp), allocatable :: density(:)
    real(dp), allocatable :: dF_dP(:), grav(:), mubar(:), cp(:)

    ! Unpack input T
    self%T_surf = T(1)
    self%T = T(2:)

    !~~ Construct the atmosphere ~~!
    call make_profile_rc( &
      .true., .true., &
      self%sp, self%planet_mass, self%planet_radius, &
      self%P_surf, self%T_surf, self%nz, self%P, self%log10_delta_P, &
      self%condensible_RH, self%condensible_P, self%condensible_inds, self%bg_gas_ind, &
      self%f_i, self%T, self%T_trop, self%P(self%trop_ind), &
      tropopause_was_reached, self%condensible_z_inds, self%z, self%dz, err &
    ) 
    if (allocated(err)) return

    ! Compute densities
    allocate(density(self%nz))
    density = self%P/(k_boltz*self%T)
    do i = 1,self%sp%ng
      self%densities(:,i) = self%f_i(:,i)*density(:)
    enddo

    ! Copy atmosphere to radiative transfer variables 
    do i = 1,self%nz
      j = 2*(i-1)
      self%T_r(j+1) = self%T(i)
      self%T_r(j+2) = self%T(i)

      self%P_r(j+1) = self%P(i)
      self%P_r(j+2) = self%P(i)

      self%densities_r(j+1,:) = self%densities(i,:)
      self%densities_r(j+2,:) = self%densities(i,:)

      self%dz_r(j+1) = 0.5_dp*self%dz(i)
      self%dz_r(j+2) = 0.5_dp*self%dz(i)
    enddo

    !~~ Radiative transfer ~~!
    if (size(self%pdensities_r,2) > 0) then
      call self%rad%radiate(self%T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, &
                            self%pdensities_r, self%radii_r, err=err)
      if (allocated(err)) return
    else
      call self%rad%radiate(self%T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, &
                            err=err)
      if (allocated(err)) return
    endif
    ! The heating rates
    allocate(dF_dP(self%nz))
    do i = 1,self%nz
      dF_dP(i) = (self%rad%f_total(2*(i-1)+3) - self%rad%f_total(2*(i-1)+1))/(-self%delta_P(i))
    enddo

    !~~ Gravity in cm/s^2 ~~!
    allocate(grav(self%nz))
    do i = 1,self%nz
      grav(i) = gravity(self%planet_radius, self%planet_mass, self%z(i))
    enddo

    !~~ Mean molecular weight ~~!
    allocate(mubar(self%nz))
    mubar = 0.0_dp
    do i = 1,self%sp%ng
      do j = 1,self%nz
        mubar(j) = mubar(j) + self%f_i(j,i)*self%sp%g(i)%mass
      enddo
    enddo

    !~~ Heat capacity in erg/(g*K) ~~!
    allocate(cp(self%nz))
    do j = 1,self%nz
      cp(j) = 0.0_dp
      do i = 1,self%sp%ng
        call heat_capacity_eval(self%sp%g(i)%thermo, self%T(j), found, cp_tmp) ! J/(mol*K)
        if (.not. found) then
          err = "Failed to compute heat capacity"
          return
        endif
        ! J/(mol*K)
        cp(j) = cp(j) + cp_tmp*self%f_i(j,i) ! J/(mol*K)
      enddo
      ! J/(mol*K) * (mol/kg) = J/(kg*K)
      cp(j) = cp(j)*(1.0_dp/(mubar(j)*1.0e-3_dp))
      ! J/(kg*K) * (erg/J) * (kg/g) = erg/(g*K)
      cp(j) = cp(j)*1.0e4_dp
    enddo

    !~~ Rate of change from radiative processes ~~!
    dT_dt(1) = (1.0_dp/(self%rho_ground*self%cp_ground))*(self%rad%f_total(1))/self%dz_ground
    dT_dt(2:) = - (grav/cp)*dF_dP ! Layers of the atmosphere

  end subroutine

  !> Adjusts tropopause lapse rate so that it matches an adiabat.
  !> The adjustment conserves energy in the troposphere.
  module subroutine RadiativeConvectiveClimate_convective_adjustment(self, T, err)
    use clima_rc_adiabat, only: make_profile_rc
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars
    class(RadiativeConvectiveClimate), intent(inout) :: self
    real(dp), intent(inout) :: T(:)
    character(:), allocatable :: err

    real(dp) :: E_trop, E_trop1
    integer :: i
    logical :: tropopause_was_reached
    type(MinpackHybrd1Vars) :: mv

    self%T_surf = T(1)
    self%T = T(2:)

    ! note that self%z, self%f_i are out of date at this point.
    E_trop = RadiativeConvectiveClimate_troposphere_energy(self, T(1), T(2:), self%trop_ind, self%z, self%f_i, err)
    E_trop = E_trop + self%rad%f_total(1)

    mv = MinpackHybrd1Vars(1, tol=1.0e-5_dp)
    mv%x(1) = T(1)
    call hybrd1(fcn, mv%n, mv%x, mv%fvec, mv%tol, mv%info, mv%wa, mv%lwa)
    if (mv%info == 0 .or. mv%info > 1) then
      err = 'hybrd1 root solve failed in make_profile_bg_gas.'
      return
    elseif (mv%info < 0) then
      err = 'hybrd1 root solve failed in make_profile_bg_gas: '//err
      return
    endif

    ! Return a convectively adjusted T profile.
    call fcn(mv%n, mv%x, mv%fvec, mv%info)
    T(1) = self%T_surf
    T(2:) = self%T

  contains
    subroutine fcn(n_, x_, fvec_, iflag_)
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(n_)
      integer, intent(inout) :: iflag_

      

      self%T_surf = x_(1)

      ! Get T above the tropopause
      do i = self%trop_ind+1,self%nz
        self%T(i) = T(i+1)
      enddo

      ! Adjust profile to adiabat
      call make_profile_rc( &
        .true., .true., &
        self%sp, self%planet_mass, self%planet_radius, &
        self%P_surf, self%T_surf, self%nz, self%P, self%log10_delta_P, &
        self%condensible_RH, self%condensible_P, self%condensible_inds, self%bg_gas_ind, &
        self%f_i, self%T, self%T_trop, self%P(self%trop_ind), &
        tropopause_was_reached, self%condensible_z_inds, self%z, self%dz, err &
      ) 
      if (allocated(err)) return

      ! Compute energy in the troposphere
      E_trop1 = RadiativeConvectiveClimate_troposphere_energy(self, self%T_surf, &
                self%T, self%trop_ind, self%z, self%f_i, err)

      ! residual is the difference between energy and what it should be.
      fvec_(1) = E_trop - E_trop1

    end subroutine

  end subroutine

  function RadiativeConvectiveClimate_troposphere_energy(self, T_surf, T, trop_ind, z, f_i, err) result(E_trop)
    use clima_eqns, only: heat_capacity_eval, gravity
    class(RadiativeConvectiveClimate), intent(inout) :: self
    real(dp), intent(in) :: T_surf
    real(dp), intent(in) :: T(:)
    integer, intent(in) :: trop_ind
    real(dp), intent(in) :: z(:)
    real(dp), intent(in) :: f_i(:,:)
    character(:), allocatable :: err
    real(dp) :: E_trop

    integer :: i, j
    real(dp) :: cp_tmp
    logical :: found
    real(dp), allocatable :: grav(:), mubar(:), cp(:)

    !~~ Gravity in cm/s^2 ~~!
    allocate(grav(trop_ind))
    do i = 1,trop_ind
      grav(i) = gravity(self%planet_radius, self%planet_mass, z(i))
    enddo

    !~~ Mean molecular weight ~~!
    allocate(mubar(trop_ind))
    mubar = 0.0_dp
    do i = 1,self%sp%ng
      do j = 1,trop_ind
        mubar(j) = mubar(j) + f_i(j,i)*self%sp%g(i)%mass
      enddo
    enddo

    !~~ Heat capacity in erg/(g*K) ~~!
    allocate(cp(trop_ind))
    do j = 1,trop_ind
      cp(j) = 0.0_dp
      do i = 1,self%sp%ng
        call heat_capacity_eval(self%sp%g(i)%thermo, T(j), found, cp_tmp) ! J/(mol*K)
        if (.not. found) then
          err = "Failed to compute heat capacity"
          return
        endif
        ! J/(mol*K)
        cp(j) = cp(j) + cp_tmp*f_i(j,i) ! J/(mol*K)
      enddo
      ! J/(mol*K) * (mol/kg) = J/(kg*K)
      cp(j) = cp(j)*(1.0_dp/(mubar(j)*1.0e-3_dp))
      ! J/(kg*K) * (erg/J) * (kg/g) = erg/(g*K)
      cp(j) = cp(j)*1.0e4_dp
    enddo

    ! (g/cm^3) * (erg/(g*K)) * (K) * (cm) = erg/cm^2
    E_trop = self%rho_ground*self%cp_ground*T_surf*self%dz_ground ! surface
    do j = 1,trop_ind
      ! (erg/(g*K)) * (s^2/cm) * (K) * (g/(cm*s^2)) = erg/cm^2
      E_trop = E_trop + (cp(j)/grav(j))*T(j)*self%delta_P(j) ! layers
    enddo

  end function
  
end submodule