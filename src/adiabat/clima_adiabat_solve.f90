
submodule(clima_adiabat) clima_adiabat_solve
  implicit none

contains

  module subroutine AdiabatClimate_make_profile_rc(self, P_i_surf, T_in, err)
    use clima_const, only: k_boltz, N_avo
    use clima_adiabat_rc, only: make_profile_rc
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:) !! dynes/cm^2
    real(dp), intent(in) :: T_in(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: P_e(:), z_e(:), f_i_e(:,:), lapse_rate_e(:)
    logical, allocatable :: super_saturated_e(:)
    real(dp), allocatable :: density(:)
    integer :: i, j

    allocate(P_e(2*self%nz+1),z_e(2*self%nz+1),f_i_e(2*self%nz+1,self%sp%ng),lapse_rate_e(2*self%nz+1))
    allocate(super_saturated_e(2*self%nz+1))
    allocate(density(self%nz))

    if (size(P_i_surf) /= self%sp%ng) then
      err = "P_i_surf has the wrong dimension"
      return
    endif
    if (size(T_in) /= self%nz+1) then
      err = "T_in has the wrong dimension"
      return
    endif

    self%T_surf = T_in(1)
    self%T = T_in(2:)
    call make_profile_rc(self%T_surf, self%T, P_i_surf, &
                         self%convecting_with_below, &
                         self%sp, self%nz, self%planet_mass, &
                         self%planet_radius, self%P_top, self%RH, &
                         self%rtol, self%atol, &
                         self%ocean_fcns, self%ocean_args_p, &
                         P_e, z_e, f_i_e, lapse_rate_e, super_saturated_e, &
                         self%N_surface, self%N_ocean, &
                         err)
    if (allocated(err)) return

    z_e(2*self%nz+1) = z_e(2*self%nz) + (z_e(2*self%nz) - z_e(2*self%nz-1))

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
    self%lapse_rate_intended(1) = lapse_rate_e(1)
    do i = 2,self%nz
      self%lapse_rate_intended(i) = lapse_rate_e(2*i-2)
    enddo

    self%lapse_rate(1) = (log(self%T(1)) - log(self%T_surf))/(log(self%P(1)) - log(self%P_surf))
    do i = 2,self%nz
      self%lapse_rate(i) = (log(self%T(i)) - log(self%T(i-1)))/(log(self%P(i)) - log(self%P(i-1)))
    enddo

    do i = 1,self%nz
      self%super_saturated(i) = super_saturated_e(2*i)
    enddo

  end subroutine

  module function AdiabatClimate_RCE(self, P_i_surf, T_surf_guess, T_guess, convecting_with_below, err) result(converged)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrj, linear_solve
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_surf_guess
    real(dp), intent(in) :: T_guess(:)
    logical, optional, intent(in) :: convecting_with_below(:)
    character(:), allocatable, intent(out) :: err
    logical :: converged

    type(MinpackHybrj) :: mv
    logical, allocatable :: convecting_with_below_save(:,:)
    real(dp), allocatable :: difference(:)
    real(dp), allocatable :: T_in(:), x_init(:)
    real(dp) :: perturbation
    integer :: i, j, k

    if (.not.self%double_radiative_grid) then
      err = 'AdiabatClimate must be initialized with "double_radiative_grid" set to True '// &
            'in order to call RCE.'
      return
    endif
    if (size(T_guess) /= self%nz) then
      err = "T_guess has the wrong dimension"
      return
    endif
    
    allocate(convecting_with_below_save(self%nz,0))
    allocate(difference(self%nz))
    allocate(T_in(self%nz+1))

    T_in(1) = T_surf_guess
    T_in(2:) = T_guess(:)
    self%T_surf = T_surf_guess
    self%T(:) = T_guess(:)

    ! setup the convecting zone
    if (present(convecting_with_below)) then
      call AdiabatClimate_set_convecting_zones(self, convecting_with_below, err)
      if (allocated(err)) return
    else
      self%convecting_with_below = .false.
      call AdiabatClimate_update_convecting_zones(self, P_i_surf, T_in, .false., err)
      if (allocated(err)) return
    endif

    converged = .false.
    do i = 1,self%max_rc_iters
      j = i

      if (self%verbose) then
        print"(1x,'Iteration =',i3)", i
      endif

      if (allocated(x_init)) deallocate(x_init)
      allocate(x_init(size(self%inds_Tx)))
      x_init(1) = self%T_surf
      do k = 2,size(self%inds_Tx)
        x_init(k) = self%T(self%inds_Tx(k)-1)
      enddo

      mv = MinpackHybrj(fcn, size(self%inds_Tx))
      mv%xtol = self%xtol_rc
      mv%nprint = 1

      k = 0
      do
        if (mod(k,2) == 0) then
          perturbation = real(k,dp)*1.0_dp
        else
          perturbation = -real(k,dp)*1.0_dp
        endif

        if (self%verbose .and. k > 0) then
          print'(3x,"Perturbation = ",f7.1)',perturbation
        endif

        mv%x = x_init + perturbation
        call mv%hybrj()
        if (mv%info == 1) then
          exit
        else
          if (mv%info < 0) deallocate(err)
        endif

        if (k > 6) then
          err = 'hybrj root solve failed in RCE.'
          return
        endif

        k = k + 1
      enddo

      ! Update all variables to the current root
      call AdiabatClimate_objective(self, P_i_surf, mv%x, .false., mv%fvec, err)
      if (allocated(err)) return

      ! Save the current convective zones
      convecting_with_below_save = reshape(convecting_with_below_save,shape=[self%nz,i],pad=self%convecting_with_below)

      ! Update the convective zones
      if (i < self%max_rc_iters_convection) then
        ! Here, we permit 
        ! - convective layers to become radiative
        ! - radiative layers to become convective
        call AdiabatClimate_update_convecting_zones(self, P_i_surf, [self%T_surf, self%T], .false., err)
        if (allocated(err)) return
      else
        ! Here, we only permit radiative layers to become convective
        call AdiabatClimate_update_convecting_zones(self, P_i_surf, [self%T_surf, self%T], .true., err)
        if (allocated(err)) return
      endif
      
      ! Check for convergence
      if (all(convecting_with_below_save(:,i) .eqv. self%convecting_with_below)) then
        ! Convective zones did not change between iterations, so converged
        converged = .true.
      endif
      if (converged) then
        if (self%verbose) then
          print'(1x,A)','CONVERGED'
        endif
        exit
      endif

    enddo

    ! Return all information to what is was prior to checking for the root
    call AdiabatClimate_set_convecting_zones(self, convecting_with_below_save(:,j), err)
    if (allocated(err)) return
    call AdiabatClimate_objective(self, P_i_surf, mv%x, .false., mv%fvec, err)
    if (allocated(err)) return

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
        call AdiabatClimate_objective(self, P_i_surf, x_, .false., fvec_, err)
        if (allocated(err)) then
          iflag_ = -1
          return
        endif
      elseif (iflag_ == 2) then
        ! Compute jacobian
        call AdiabatClimate_jacobian(self, P_i_surf, x_, .false., fjac_, err)
        if (allocated(err)) then
          iflag_ = -1
          return
        endif
      endif

      if (iflag_ == 0 .and. self%verbose) then
        print"(3x,'step =',i3,3x,'njev =',i3,3x,'||y|| = ',es8.2,3x,'max(T) = ',f7.1,3x,'min(T) = ',f7.1)", &
              mv%nfev, mv%njev, sum(fvec_**2.0_dp), maxval(x_), minval(x_)
      endif

    end subroutine

  end function

  subroutine AdiabatClimate_objective(self, P_i_surf, x, ignoring_convection, res, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: x(:)
    logical, intent(in) :: ignoring_convection
    real(dp), intent(out) :: res(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: T_in(:)
    integer :: i

    allocate(T_in(self%nz+1))
    do i = 1,size(self%inds_Tx)
      T_in(self%inds_Tx(i)) = x(i)
    enddo

    ! Sets, self%T_surf, self%T, self%N_surface, self%N_ocean, self%P_surf, self%P,
    ! self%z, self%dz, self%f_i, self%densities, self%N_atmos, self%lapse_rate_intended
    ! self%lapse_rate
    call self%make_profile_rc(P_i_surf, T_in, err)
    if (allocated(err)) return

    T_in(1) = self%T_surf
    T_in(2:) = self%T

    ! resets self%T_surf, self%T, self%densities, self%lapse_rate
    ! also does radiative transfer and computes res
    call AdiabatClimate_objective_(self, P_i_surf, T_in, ignoring_convection, res, err)
    if (allocated(err)) return

  end subroutine

  subroutine AdiabatClimate_objective_(self, P_i_surf, T_in, ignoring_convection, res, err)
    use clima_const, only: k_boltz
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_in(:)
    logical, intent(in) :: ignoring_convection
    real(dp), intent(out) :: res(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: f_total(:)
    real(dp), allocatable :: density(:)
    integer :: i, j

    ! work storage
    allocate(f_total(self%nz+1))
    allocate(density(self%nz))

    ! Stuff set in make_profile_rc that does not change
    ! - self%N_surface, self%N_ocean, self%P_surf, self%P, self%z, 
    ! - self%dz, self%f_i, self%N_atmos, self%lapse_rate_intended
    ! Stuff that does change
    ! - self%T_surf, self%T, self%densities, self%lapse_rate
    self%T_surf = T_in(1)
    self%T = T_in(2:)

    density = self%P/(k_boltz*self%T)
    do j =1,self%sp%ng
      self%densities(:,j) = self%f_i(:,j)*density(:)
    enddo

    self%lapse_rate(1) = (log(self%T(1)) - log(self%T_surf))/(log(self%P(1)) - log(self%P_surf))
    do i = 2,self%nz
      self%lapse_rate(i) = (log(self%T(i)) - log(self%T(i-1)))/(log(self%P(i)) - log(self%P(i-1)))
    enddo

    ! radiative transfer
    call self%copy_atm_to_radiative_grid()
    call self%rad%radiate(self%T_surf, self%T_r, self%P_r/1.0e6_dp, self%densities_r, self%dz_r, err=err)
    if (allocated(err)) return
    if (self%tidally_locked_dayside) then; block
      real(dp) :: tau_LW, k_term, f_term, rad_enhancement
      call self%heat_redistribution_parameters(tau_LW, k_term, f_term, err)
      if (allocated(err)) return

      rad_enhancement = 4.0_dp*f_term
      call self%rad%apply_radiation_enhancement(rad_enhancement)
    endblock; endif

    do i = 1,self%nz+1
      f_total(i) = self%rad%f_total(2*i-1)
    enddo
    f_total(1) = f_total(1) + self%surface_heat_flow

    if (ignoring_convection) then
      call AdiabatClimate_residuals_with_heat_capacity(self, f_total, res, err)
      if (allocated(err)) return
    else
      call AdiabatClimate_residuals_with_convection(self, f_total, self%lapse_rate, self%lapse_rate_intended, res)
    endif

  end subroutine

  subroutine AdiabatClimate_jacobian(self, P_i_surf, x, ignoring_convection, jac, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: x(:)
    logical, intent(in) :: ignoring_convection
    real(dp), intent(out) :: jac(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: res(:)
    real(dp), allocatable :: res_perturb(:), T_perturb(:), T_in(:)
    real(dp) :: deltaT

    integer :: i, ind

    ! Check inputs
    if (size(jac,1) /= size(self%inds_Tx) .or. size(jac,2) /= size(self%inds_Tx)) then
      err = "jac has the wrong dimension"
      return
    endif

    ! allocate work
    allocate(res(size(self%inds_Tx)))
    allocate(T_in(self%nz+1),res_perturb(size(self%inds_Tx)))

    ! First evaluate res at T.
    call AdiabatClimate_objective(self, P_i_surf, x, ignoring_convection, res, err)
    if (allocated(err)) return

    T_in(1) = self%T_surf
    T_in(2:) = self%T

    T_perturb = T_in
    do i = 1,size(x)

      ! Perturb temperature
      deltaT = self%epsj*abs(x(i))
      T_perturb(self%inds_Tx(i)) = T_in(self%inds_Tx(i)) + deltaT

      ! We perturb the temp also in a convecting region
      ind = findloc(self%ind_conv_lower_x, i, 1)
      if (ind /= 0) then
        T_perturb(self%ind_conv_lower(ind):self%ind_conv_upper(ind)) = &
        T_in(self%ind_conv_lower(ind):self%ind_conv_upper(ind)) + deltaT
      endif

      call AdiabatClimate_objective_(self, P_i_surf, T_perturb, ignoring_convection, res_perturb, err)
      if (allocated(err)) return

      ! Compute jacobian
      jac(:,i) = (res_perturb(:) - res(:))/deltaT

      ! unperturb T
      T_perturb(:) = T_in(:)
    enddo

  end subroutine

  subroutine AdiabatClimate_set_convecting_zones(self, convecting_with_below, err)
    class(AdiabatClimate), intent(inout) :: self
    !> Length nz. States whether the layer below is convecting
    !> with the current layer. `convecting(1)` determines whether
    !> the first atmospheric layer is convecting with the ground.
    logical, intent(in) :: convecting_with_below(:)
    character(:), allocatable, intent(out) :: err

    integer :: i,j,k,ind

    if (size(convecting_with_below) /= self%nz) then
      err = 'Input "convecting_with_below" has the wrong dimension'
      return
    endif

    self%convecting_with_below = convecting_with_below

    self%n_convecting_zones = 0
    if (allocated(self%ind_conv_lower)) deallocate(self%ind_conv_lower)
    if (allocated(self%ind_conv_upper)) deallocate(self%ind_conv_upper)
    allocate(self%ind_conv_lower(0),self%ind_conv_upper(0))
    i = 1
    do while(i <= self%nz)
      if (convecting_with_below(i)) then
        ! identified a convecting zone
        self%n_convecting_zones = self%n_convecting_zones + 1
        ! Save the lower index
        self%ind_conv_lower = [self%ind_conv_lower, i]
        ! Search for the upper bound of convecting zone
        do j = i,self%nz
          if (convecting_with_below(j)) then
            k = j + 1
          else
            ! convecting zone stops
            exit
          endif
        enddo
        self%ind_conv_upper = [self%ind_conv_upper, k]
        i = j
        cycle
      endif
      i = i + 1
    enddo

    if (allocated(self%inds_Tx)) deallocate(self%inds_Tx)
    allocate(self%inds_Tx(1))
    self%inds_Tx(1) = 1
    do i = 1,size(convecting_with_below)
      if (convecting_with_below(i)) then
        ! nothing
      else
        self%inds_Tx = [self%inds_Tx, i+1]
      endif
    enddo

    if (allocated(self%ind_conv_lower_x)) deallocate(self%ind_conv_lower_x)
    allocate(self%ind_conv_lower_x(self%n_convecting_zones))
    do i = 1,size(self%ind_conv_lower)
      ind = findloc(self%inds_Tx, self%ind_conv_lower(i), 1)
      if (ind == 0) then
        err = 'Problem setting a convective zone'
        return
      endif
      self%ind_conv_lower_x(i) = ind
    enddo

  end subroutine

  !> Given P_surf and T profile, this routine determines which layers are convecting
  !> and which are purely radiative. To accomplish this, the routine does a small
  !> Newton step for a purely radiative atmosphere. In a layer, if the step tends towards
  !> a T profile stable to convection, then a layer is radiative. If instead, the step
  !> tends to a T profile unstable to convection, then the layer is convective. The
  !> routine specifically updates self%convecting_with_below, self%n_convecting_zones, 
  !> self%ind_conv_lower and self%ind_conv_upper. It also updates all atmosphere varibles.
  subroutine AdiabatClimate_update_convecting_zones(self, P_i_surf, T_in, no_convection_to_radiation, err)
    use clima_useful, only: linear_solve
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_in(:)
    logical, intent(in) :: no_convection_to_radiation
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: F(:), dFdT(:,:), deltaT(:), T_perturb(:)
    real(dp), allocatable :: lapse_rate_perturb(:), difference(:)
    logical, allocatable :: convecting_with_below_save(:)
    integer :: i, ierr

    ! work storage
    allocate(F(size(T_in)),dFdT(size(T_in),size(T_in)),deltaT(size(T_in)),T_perturb(size(T_in)))
    allocate(lapse_rate_perturb(self%nz),difference(self%nz))

    convecting_with_below_save = self%convecting_with_below
    self%convecting_with_below = .false.
    call AdiabatClimate_set_convecting_zones(self, self%convecting_with_below, err)
    if (allocated(err)) return

    call AdiabatClimate_objective(self, P_i_surf, T_in, .true., F, err)
    if (allocated(err)) return
    call AdiabatClimate_jacobian(self, P_i_surf, T_in, .true., dFdT, err)
    if (allocated(err)) return

    deltaT = -F
    call linear_solve(dFdT, deltaT, ierr)
    if (ierr /= 0) then
      err = 'Linear solved failed in "update_convecting_zones"'
      return
    endif

    ! Newton step
    T_perturb = deltaT*self%convective_newton_step_size + T_in

    ! Compute perturbed lapse rate
    call self%make_profile_rc(P_i_surf, T_perturb, err)
    if (allocated(err)) return
    lapse_rate_perturb = self%lapse_rate

    ! Re-update all variables at T_in, including self%lapse_rate_intended
    call AdiabatClimate_objective(self, P_i_surf, T_in, .true., F, err)
    if (allocated(err)) return

    difference = lapse_rate_perturb - self%lapse_rate_intended

    if (no_convection_to_radiation) then
      self%convecting_with_below = convecting_with_below_save
      do i = 1,self%nz
        if (.not.self%convecting_with_below(i)) then
          difference(i) = self%lapse_rate(i) - self%lapse_rate_intended(i)
        endif
        if (difference(i) >= 0.0_dp) then
          self%convecting_with_below(i) = .true.
        endif
      enddo
    else
      do i = 1,self%nz
        if (difference(i) >= 0.0_dp) then
          self%convecting_with_below(i) = .true.
        else
          self%convecting_with_below(i) = .false.
        endif
      enddo
    endif

    call AdiabatClimate_set_convecting_zones(self, self%convecting_with_below, err)
    if (allocated(err)) return

  end subroutine

  subroutine AdiabatClimate_residuals_with_convection(self, f_total, lapse_rate, lapse_rate_intended, res)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: f_total(:) !! fluxes at the edges of layers
    real(dp), intent(in) :: lapse_rate(:)
    real(dp), intent(in) :: lapse_rate_intended(:)
    real(dp), intent(out) :: res(:)

    real(dp), allocatable :: fluxes(:)
    real(dp) :: f_lower, f_upper
    integer :: i, ind_lower, ind_upper

    ! work storage
    allocate(fluxes(self%nz+1))

    ! Radiative energy going into each layer (ergs/(cm^2*s))
    fluxes(1) = f_total(1)
    do i = 2,self%nz+1
      fluxes(i) = (f_total(i) - f_total(i-1))
    enddo

    ! Radiative equilibrium
    do i = 1,size(self%inds_Tx)
      res(i) = fluxes(self%inds_Tx(i))
    enddo

    ! Change residual to account for convection
    do i = 1,self%n_convecting_zones

      ind_lower = self%ind_conv_lower(i)
      if (ind_lower == 1) then
        f_lower = 0.0_dp
      else
        f_lower = f_total(ind_lower-1)
      endif

      ind_upper = self%ind_conv_upper(i)
      if (ind_lower == 1) then
        f_upper = f_total(ind_upper) + self%surface_heat_flow
      else
        f_upper = f_total(ind_upper)
      endif

      ! Radiative energy going into the convective layer      
      res(self%ind_conv_lower_x(i)) = f_upper - f_lower

    enddo

  end subroutine

  subroutine AdiabatClimate_residuals_with_heat_capacity(self, f_total, res, err)
    use clima_eqns, only: heat_capacity_eval
    use clima_const, only: k_boltz, N_avo
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: f_total(:) !! fluxes at the edges of layers (ergs/(cm^2*s))
    real(dp), intent(out) :: res(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: fluxes(:), mubar(:), cp(:), rho(:), density(:)
    real(dp) :: cp_tmp
    integer :: i, j
    logical :: found

    if (self%n_convecting_zones /= 0) then
      err = 'residuals_with_heat_capacity can not be called when there is convection'
      return
    endif
    if (size(res) /= size(f_total)) then
      err = 'res has the wrong shape in residuals_with_heat_capacity'
      return
    endif

    ! work storage
    allocate(fluxes(self%nz+1), mubar(self%nz), cp(self%nz), rho(self%nz), density(self%nz))

    ! Radiative energy going into each layer (ergs/(cm^2*s))
    fluxes(1) = f_total(1)
    do i = 2,self%nz+1
      fluxes(i) = (f_total(i) - f_total(i-1))
    enddo

    ! mean molecular weight (g/mol)
    do j = 1,self%nz
      mubar(j) = 0.0_dp
      do i = 1,self%sp%ng
        mubar(j) = mubar(j) + self%f_i(j,i)*self%sp%g(i)%mass
      enddo
    enddo

    ! molecules/cm^3
    density = self%P/(k_boltz*self%T)
    ! [molecules/cm^3]*[mol/molecules]*[g/mol] = [g/cm^3]
    rho = density*(1.0_dp/N_avo)*mubar 

    ! Heat capacity in erg/(g*K)
    do j = 1,self%nz
      cp(j) = 0.0_dp
      do i = 1,self%sp%ng
        call heat_capacity_eval(self%sp%g(i)%thermo, self%T(j), found, cp_tmp) ! J/(mol*K)
        if (.not. found) then
          err = "not found"
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

    ! [ergs/(cm^2*s)] * [1/cm] * [cm^3/g] * [g*K/erg] = [K/s]
    res(1) = (fluxes(1)/self%dz(1))*(1.0_dp/(rho(1)*cp(1)))
    do i = 1,self%nz
      res(i+1) = (fluxes(i+1)/self%dz(i))*(1.0_dp/(rho(i)*cp(i)))
    enddo

  end subroutine

end submodule