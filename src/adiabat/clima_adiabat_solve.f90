
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
    real(dp), allocatable :: density(:)
    integer :: i, j

    allocate(P_e(2*self%nz+1),z_e(2*self%nz+1),f_i_e(2*self%nz+1,self%sp%ng),lapse_rate_e(2*self%nz+1))
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
    self%lapse_rate_intended(1) = lapse_rate_e(1)
    do i = 2,self%nz
      self%lapse_rate_intended(i) = lapse_rate_e(2*i-1)
    enddo

    self%lapse_rate(1) = (log(self%T(1)) - log(self%T_surf))/(log(self%P(1)) - log(self%P_surf))
    do i = 2,self%nz
      self%lapse_rate(i) = (log(self%T(i)) - log(self%T(i-1)))/(log(self%P(i)) - log(self%P(i-1)))
    enddo

  end subroutine

  module subroutine AdiabatClimate_RCE(self, P_i_surf, T_guess, err)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrj, linear_solve
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_guess
    character(:), allocatable, intent(out) :: err

    real(dp) :: T_surf
    type(MinpackHybrj) :: mv

    logical, allocatable :: convecting(:)
    integer :: i

    if (.not.self%double_radiative_grid) then
      err = 'AdiabatClimate must be initialized with "double_radiative_grid" set to True '// &
            'in order to call RCE.'
      return
    endif

    ! Create the initial temperature profile
    T_surf = self%surface_temperature(P_i_surf, T_guess, err)
    if (allocated(err)) return

    ! call self%make_profile(T_guess, P_i_surf, err)
    ! if (allocated(err)) return

    if (self%P_trop <= 0.0_dp) then
      err = 'error!'
      return
    endif

    ! is layer i convecting with the layer below it?
    allocate(convecting(self%nz))
    do i = 1,self%nz
      if (self%P(i) > self%P_trop) then
        convecting(i) = .true.
      else
        convecting(i) = .false.
      endif
    enddo

    ! call AdiabatClimate_set_convecting_zones(self, convecting, err)
    ! if (allocated(err)) return
    call AdiabatClimate_update_convecting_zones(self, P_i_surf, [self%T_surf,self%T], err)
    if (allocated(err)) return

    ! self%ind_conv_upper(1) = self%ind_conv_upper(1)
    print*,self%n_convecting_zones
    print*,self%ind_conv_upper(1)

    open(unit=2,file='test.dat',form='unformatted',status='replace')
    write(unit=2) self%nz
    write(unit=2) self%sp%ng
    write(unit=2) self%ind_conv_upper(1)
  
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

      if (iflag_ == 0) then; block
        real(dp), allocatable :: F(:), dFdT(:,:), deltaT(:) 
        print*,maxval(x_),minval(x_),sum(fvec_**2.0_dp),mv%nfev, mv%njev
        write(2) [self%P_surf,self%P]
        write(2) x_
        write(2) fvec_
        write(2) self%lapse_rate
        write(2) self%lapse_rate_intended
        write(2) sum(fvec_**2.0_dp)
        write(2) self%f_i

        allocate(F(self%nz+1),dFdT(self%nz+1,self%nz+1),deltaT(self%nz+1))

        call AdiabatClimate_objective(self, P_i_surf, x_, .true., F, err)
        if (allocated(err)) then
          iflag_ = -1
          return
        endif
        call AdiabatClimate_jacobian(self, P_i_surf, x_, .true., dFdT, err)
        if (allocated(err)) then
          iflag_ = -1
          return
        endif
        deltaT = -F
        call linear_solve(dFdT, deltaT, i)
        if (i /= 0) then
          err = 'linear solve failed'
          iflag_ = -1
          return
        endif
        write(2) deltaT

      endblock; endif

    end subroutine

  end subroutine

  subroutine AdiabatClimate_objective(self, P_i_surf, T_in, ignore_convection, res, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_in(:)
    logical, intent(in) :: ignore_convection
    real(dp), intent(out) :: res(:)
    character(:), allocatable, intent(out) :: err

    ! Sets, self%T_surf, self%T, self%N_surface, self%N_ocean, self%P_surf, self%P,
    ! self%z, self%dz, self%f_i, self%densities, self%N_atmos, self%lapse_rate_intended
    ! self%lapse_rate
    call self%make_profile_rc(P_i_surf, T_in, err)
    if (allocated(err)) return

    ! resets self%T_surf, self%T, self%densities, self%lapse_rate
    ! also does radiative transfer and computes res
    call AdiabatClimate_objective_(self, P_i_surf, T_in, ignore_convection, res, err)
    if (allocated(err)) return

  end subroutine

  subroutine AdiabatClimate_objective_(self, P_i_surf, T_in, ignore_convection, res, err)
    use clima_const, only: k_boltz
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_in(:)
    logical, intent(in) :: ignore_convection
    real(dp), intent(out) :: res(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: f_total(:), fluxes(:)
    real(dp), allocatable :: density(:)
    integer :: i, j

    ! work storage
    allocate(f_total(self%nz+1))
    allocate(fluxes(self%nz+1))
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

    do i = 1,self%nz+1
      f_total(i) = self%rad%f_total(2*i-1)
    enddo

    ! Radiative energy going into each layer (ergs/(cm^2*s))
    fluxes(1) = f_total(1)
    do i = 2,self%nz+1
      fluxes(i) = (f_total(i) - f_total(i-1))
    enddo

    if (ignore_convection) then
      res = fluxes
      return
    endif

    ! j = self%ind_conv_upper(1)
    ! ! Radiative layers
    ! res(j+1:) = fluxes(j+1:)
    ! ! Convective layers
    ! i = j - 1 ! i is the atmospheric layer
    ! res(j) = f_total(i+1)
    ! do i = 1,j-1
    !   res(i) = self%lapse_rate(i) - self%lapse_rate_intended(i)
    ! enddo

    call AdiabatClimate_residuals_with_convection(self, f_total, self%lapse_rate, self%lapse_rate_intended, res)

  end subroutine

  subroutine AdiabatClimate_jacobian(self, P_i_surf, T_in, ignore_convection, jac, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_in(:)
    logical, intent(in) :: ignore_convection
    real(dp), intent(out) :: jac(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: res(:)
    real(dp), allocatable :: T_perturb(:), res_perturb(:)
    real(dp) :: deltaT

    integer :: i

    ! Check inputs
    if (size(jac,1) /= self%nz+1 .or. size(jac,2) /= self%nz+1) then
      err = "jac has the wrong dimension"
      return
    endif

    ! allocate work
    allocate(res(self%nz+1))
    allocate(T_perturb(self%nz+1),res_perturb(self%nz+1))

    ! First evaluate res at T.
    call AdiabatClimate_objective(self, P_i_surf, T_in, ignore_convection, res, err)
    if (allocated(err)) return

    ! Now compute Jacobian
    T_perturb = T_in
    do i = 1,size(T_perturb)

      ! Perturb temperature
      deltaT = self%epsj*abs(T_in(i))
      T_perturb(i) = T_in(i) + deltaT

      ! Compute rhs_perturb
      call AdiabatClimate_objective_(self, P_i_surf, T_perturb, ignore_convection, res_perturb, err)
      if (allocated(err)) return

      ! Compute jacobian
      jac(:,i) = (res_perturb(:) - res(:))/deltaT
      
      ! unperturb T
      T_perturb(i) = T_in(i)
    enddo

  end subroutine

  subroutine AdiabatClimate_set_convecting_zones(self, convecting_with_below, err)
    class(AdiabatClimate), intent(inout) :: self
    !> Length nz. States whether the layer below is convecting
    !> with the current layer. `convecting(1)` determines whether
    !> the first atmospheric layer is convecting with the ground.
    logical, intent(in) :: convecting_with_below(:)
    character(:), allocatable, intent(out) :: err

    integer :: i,j,k

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

  end subroutine

  !> Given P_surf and T profile, this routine determines which layers are convecting
  !> and which are purely radiative. To accomplish this, the routine does a small
  !> Newton step for a purely radiative atmosphere. In a layer, if the step tends towards
  !> a T profile stable to convection, then a layer is radiative. If instead, the step
  !> tends to a T profile unstable to convection, then the layer is convective. The
  !> routine specifically updates self%convecting_with_below, self%n_convecting_zones, 
  !> self%ind_conv_lower and self%ind_conv_upper. It also updates all atmosphere varibles.
  subroutine AdiabatClimate_update_convecting_zones(self, P_i_surf, T_in, err)
    use clima_useful, only: linear_solve
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: T_in(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: F(:), dFdT(:,:), deltaT(:), T_perturb(:)
    real(dp), allocatable :: lapse_rate_perturb(:), difference(:)
    integer :: i, ierr

    ! work storage
    allocate(F(size(T_in)),dFdT(size(T_in),size(T_in)),deltaT(size(T_in)),T_perturb(size(T_in)))
    allocate(lapse_rate_perturb(self%nz),difference(self%nz))

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

    do i = 1,self%nz
      if (difference(i) >= 0.0_dp) then
        self%convecting_with_below(i) = .true.
      else
        self%convecting_with_below(i) = .false.
      endif
    enddo

    call AdiabatClimate_set_convecting_zones(self, self%convecting_with_below, err)
    if (allocated(err)) return

  end subroutine

  !> Given
  subroutine AdiabatClimate_residuals_with_convection(self, f_total, lapse_rate, lapse_rate_intended, res)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), intent(in) :: f_total(:) !! fluxes at the edges of layers
    real(dp), intent(in) :: lapse_rate(:)
    real(dp), intent(in) :: lapse_rate_intended(:)
    real(dp), intent(out) :: res(:)

    real(dp), allocatable :: fluxes(:)
    real(dp) :: f_lower, f_upper
    integer :: i, j, ind_lower, ind_upper

    ! work storage
    allocate(fluxes(self%nz+1))

    ! Radiative energy going into each layer (ergs/(cm^2*s))
    fluxes(1) = f_total(1)
    do i = 2,self%nz+1
      fluxes(i) = (f_total(i) - f_total(i-1))
    enddo

    ! Radiative equilibrium
    res = fluxes

    ! Change residual to account for convection
    do i = 1,self%n_convecting_zones
      ind_lower = self%ind_conv_lower(i)
      f_lower = f_total(ind_lower)
      if (ind_lower == 1) then
        f_lower = 0.0_dp
      else
        f_lower = f_total(ind_lower-1)
      endif

      ind_upper = self%ind_conv_upper(i)
      f_upper = f_total(ind_upper)

      ! Radiative energy going into the convective layer      
      res(ind_upper) = f_upper - f_lower

      ! Enforcing the convective lapse rate
      do j = ind_lower,ind_upper-1
        res(j) = lapse_rate(j) - lapse_rate_intended(j)
      enddo

    enddo

  end subroutine

end submodule