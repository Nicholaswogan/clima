submodule(clima_rc) clima_rc_rhs
  implicit none

contains

  !> Computes dT/dt from radiative effects.
  module subroutine RadiativeConvectiveClimate_rhs(self, T, dTdt, err)
    use clima_rc_adiabat, only: make_profile_rc
    use clima_const, only: k_boltz
    class(RadiativeConvectiveClimate), intent(inout) :: self
    real(dp), intent(in) :: T(:)
    real(dp), intent(out) :: dTdt(:)
    character(:), allocatable :: err

    integer :: i, j
    logical :: tropopause_was_reached
    real(dp), allocatable :: density(:)
    real(dp), allocatable :: dF_dP(:)

    ! Unpack input T
    self%T_surf = T(1)
    self%T = T(2:)

    ! Construct the atmosphere
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

    ! Radiative transfer
    if (size(self%pdensities_r,2) > 0) then
      call self%rad%radiate(self%T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, &
                            self%pdensities_r, self%radii_r, err=err)
      if (allocated(err)) return
    else
      call self%rad%radiate(self%T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, &
                            err=err)
      if (allocated(err)) return
    endif
    ! allocate(dFdP(self%nz))
    ! do i = 1,self%nz
    !   dF_dP(i) = (self%rad%f_total(2*(i-1)+3) - self%rad%f_total(2*(i-1)+1))/self%dz(i)
    ! enddo

    ! Heat capacity
    
    ! Rate of change from radiative processes
    ! dT_dt = 

  end subroutine

  !> Adjusts tropopause lapse rate so that it matches an adiabat.
  !> The adjustment conserves energy in the troposphere.
  module subroutine RadiativeConvectiveClimate_convective_adjustment()

    ! Start at surface.


  end subroutine
  
end submodule