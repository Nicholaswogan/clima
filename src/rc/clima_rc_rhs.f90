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
    dT_dt(1) = - (grav(1)/cp(1))*(self%rad%f_total(1)/(-self%delta_P(1))) ! Ground
    dT_dt(2:) = - (grav/cp)*dF_dP ! Layers of the atmosphere

  end subroutine

  !> Adjusts tropopause lapse rate so that it matches an adiabat.
  !> The adjustment conserves energy in the troposphere.
  module subroutine RadiativeConvectiveClimate_convective_adjustment()

    ! Start at surface.


  end subroutine
  
end submodule