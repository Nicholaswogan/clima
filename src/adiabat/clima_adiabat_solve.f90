
submodule(clima_adiabat) clima_adiabat_solve
  implicit none

contains

  module subroutine AdiabatClimate_make_profile_rc(self, P_i_surf, T_surf, T, err)
    use clima_const, only: k_boltz, N_avo
    use clima_adiabat_rc, only: make_profile_rc
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:) !! dynes/cm^2
    real(dp), intent(in) :: T_surf, T(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: P_e(:), z_e(:), f_i_e(:,:), lapse_rate_e(:)
    real(dp), allocatable :: density(:)
    integer :: i, j

    allocate(P_e(2*self%nz+1),z_e(2*self%nz+1),f_i_e(2*self%nz+1,self%sp%ng),lapse_rate_e(2*self%nz+1))
    allocate(density(self%nz))

    if (size(P_i_surf) /= self%sp%ng) then
      err = "P_i_surf has the wrong dimension"
      return
    endif
    if (size(T) /= self%nz) then
      err = "T has the wrong dimension"
    endif

    self%T_surf = T_surf
    self%T = T
    call make_profile_rc(self%T_surf, self%T, P_i_surf, &
                         self%sp, self%nz, self%planet_mass, &
                         self%planet_radius, self%P_top, self%RH, &
                         self%rtol, self%atol, &
                         self%ocean_fcns, self%ocean_args_p, &
                         P_e, z_e, f_i_e, lapse_rate_e, &
                         self%N_surface, self%N_ocean, &
                         err)
    if (allocated(err)) return

    self%P_surf = P_e(1)
    do i = 1,self%nz
      self%P(i) = P_e(2*i)
      self%z(i) = z_e(2*i)
      self%dz(i) = z_e(2*i+1) - z_e(2*i-1)
      do j =1,self%sp%ng
        self%f_i(i,j) = f_i_e(2*i,j)
      enddo
    enddo

    density = self%P/(k_boltz*self%T)
    do j =1,self%sp%ng
      self%densities(:,j) = self%f_i(:,j)*density(:)
    enddo

    do i = 1,self%sp%ng
      ! mol/cm^2 in atmosphere
      self%N_atmos(i) = sum(density*self%f_i(:,i)*self%dz)/N_avo
    enddo

    ! lapse rates
    self%rcwrk%lapse_rate(1) = lapse_rate_e(1)
    do i = 2,self%nz
      self%rcwrk%lapse_rate(i) = lapse_rate_e(2*i-1)
    enddo

    self%rcwrk%lapse_rate_obs(1) = (log(self%T(1)) - log(self%T_surf))/(log(self%P(1)) - log(self%P_surf))
    do i = 2,self%nz
      self%rcwrk%lapse_rate_obs(i) = (log(self%T(i)) - log(self%T(i-1)))/(log(self%P(i)) - log(self%P(i-1)))
    enddo

  end subroutine

  module subroutine AdiabatClimate_RCE(self, P_i_surf, T_guess, err)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrj
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_guess
    character(:), allocatable, intent(out) :: err

    real(dp) :: T_surf
    type(MinpackHybrj) :: mv

    ! Create the initial temperature profile
    T_surf = self%surface_temperature(P_i_surf, T_guess, err)
    if (allocated(err)) return

    ! call self%make_profile(T_guess, P_i_surf, err)
    ! if (allocated(err)) return

    ! Establish the convective zone
    if (self%P_trop > 0.0_dp) then
      self%rcwrk%n_convecting_zones = 1
      allocate(self%rcwrk%ind_conv_lower(1))
      allocate(self%rcwrk%ind_conv_upper(1))
      self%rcwrk%ind_conv_lower(1) = 1
      self%rcwrk%ind_conv_upper(1) = minloc(abs(self%P_trop - [self%P_surf, self%P]),dim=1)
    else
      err = 'error!'
      return
    endif

    open(unit=2,file='test.dat',form='unformatted',status='replace')
    write(unit=2) self%nz
    write(unit=2) self%sp%ng
  
    ! Non-linear solve
    mv = MinpackHybrj(fcn,self%nz+1)
    mv%x(1) = self%T_surf
    mv%x(2:) = self%T
    mv%nprint = 1 ! enable printing
    mv%xtol = 1.0e-8_dp
    call mv%hybrj()
    if (mv%info == 0 .or. mv%info > 1) then
      print*,mv%info
      err = 'hybrj root solve failed in surface_temperature.'
      return
    elseif (mv%info < 0) then
      err = 'hybrj root solve failed in surface_temperature: '//err
      return
    endif

    close(2)

  contains

    subroutine fcn(n_, x_, fvec_, fjac_, ldfjac_, iflag_)
      implicit none
      integer, intent(in) :: n_
      real(dp), dimension(n_), intent(in) :: x_
      integer, intent(in) :: ldfjac_
      real(dp), dimension(n_), intent(inout) :: fvec_
      real(dp), dimension(ldfjac_, n_), intent(inout) :: fjac_
      integer, intent(inout) :: iflag_

      if (allocated(err)) return

      if (iflag_ == 1) then
        ! Compute right-hand-side
        call AdiabatClimate_objective(self, P_i_surf, x_, fvec_, err)
        if (allocated(err)) then
          iflag_ = -1
          return
        endif
      elseif (iflag_ == 2) then
        ! Compute jacobian
        call AdiabatClimate_jacobian(self, P_i_surf, x_, fjac_, err)
        if (allocated(err)) then
          iflag_ = -1
          return
        endif
      endif

      if (iflag_ == 0) then
        print*,maxval(x_),minval(x_),sum(fvec_**2.0_dp),mv%nfev, mv%njev
        write(2) [self%P_surf,self%P]
        write(2) x_
        write(2) fvec_
        write(2) sum(fvec_**2.0_dp)
        write(2) self%f_i
      endif

    end subroutine

  end subroutine

  subroutine AdiabatClimate_objective(self, P_i_surf, T_in, res, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_in(:)
    real(dp), intent(out) :: res(:)
    character(:), allocatable, intent(out) :: err

    ! Sets, self%T_surf, self%T, self%N_surface, self%N_ocean, self%P_surf, self%P,
    ! self%z, self%dz, self%f_i, self%densities, self%N_atmos, self%rcwrk%lapse_rate
    ! self%rcwrk%lapse_rate_obs
    self%T_surf = T_in(1)
    self%T = T_in(2:)
    call self%make_profile_rc(P_i_surf, self%T_surf, self%T, err)
    if (allocated(err)) return

    ! resets self%T_surf, self%T, self%densities, self%rcwrk%lapse_rate_obs
    ! also does radiative transfer and computes res
    call AdiabatClimate_objective_(self, P_i_surf, T_in, res, err)
    if (allocated(err)) return

  end subroutine

  subroutine AdiabatClimate_objective_(self, P_i_surf, T_in, res, err)
    use clima_const, only: k_boltz
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_in(:)
    real(dp), intent(out) :: res(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: fluxes(:)
    real(dp), allocatable :: density(:)
    integer :: i, j

    ! work storage
    allocate(fluxes(self%nz+1))
    allocate(density(self%nz))

    ! Stuff set in make_profile_rc that does not change
    ! - self%N_surface, self%N_ocean, self%P_surf, self%P, self%z, 
    ! - self%dz, self%f_i, self%N_atmos, self%rcwrk%lapse_rate
    ! Stuff that does change
    ! - self%T_surf, self%T, self%densities, self%rcwrk%lapse_rate_obs
    self%T_surf = T_in(1)
    self%T = T_in(2:)

    density = self%P/(k_boltz*self%T)
    do j =1,self%sp%ng
      self%densities(:,j) = self%f_i(:,j)*density(:)
    enddo

    self%rcwrk%lapse_rate_obs(1) = (log(self%T(1)) - log(self%T_surf))/(log(self%P(1)) - log(self%P_surf))
    do i = 2,self%nz
      self%rcwrk%lapse_rate_obs(i) = (log(self%T(i)) - log(self%T(i-1)))/(log(self%P(i)) - log(self%P(i-1)))
    enddo

    ! radiative transfer
    call self%copy_atm_to_radiative_grid()
    call self%rad%radiate(self%T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, err=err)
    if (allocated(err)) return

    ! Energy going into each layer (ergs/(cm^2*s))
    fluxes(1) = self%rad%f_total(1)
    do i = 1,self%nz
      fluxes(i+1) = (self%rad%f_total(2*(i-1)+3) - self%rad%f_total(2*(i-1)+1))
    enddo

    ! Radiative layers
    res(j+1:) = fluxes(j+1:)

    ! Convective layers
    j = self%rcwrk%ind_conv_upper(1)
    i = j - 1 ! i is the atmospheric layer
    res(j) = self%rad%f_total(2*(i-1)+3)
    do i = 1,j-1
      res(i) = self%rcwrk%lapse_rate_obs(i) - self%rcwrk%lapse_rate(i)
    enddo

  end subroutine

  subroutine AdiabatClimate_jacobian(self, P_i_surf, T_in, jac, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_in(:)
    real(dp), intent(out) :: jac(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: res(:)
    real(dp), allocatable :: T_perturb(:), res_perturb(:)
    real(dp) :: deltaT

    integer :: i

    ! allocate work
    allocate(res(self%nz+1))
    allocate(T_perturb(self%nz+1),res_perturb(self%nz+1))

    ! First evaluate res at T.
    call AdiabatClimate_objective(self, P_i_surf, T_in, res, err)
    if (allocated(err)) return

    ! Now compute Jacobian
    T_perturb = T_in
    do i = 1,size(T_perturb)

      ! Perturb temperature
      deltaT = self%epsj*abs(T_in(i))
      T_perturb(i) = T_in(i) + deltaT

      ! Compute rhs_perturb
      call AdiabatClimate_objective_(self, P_i_surf, T_perturb, res_perturb, err)
      if (allocated(err)) return

      ! Compute jacobian
      jac(:,i) = (res_perturb(:) - res(:))/deltaT
      
      ! unperturb T
      T_perturb(i) = T_in(i)
    enddo

  end subroutine


end submodule