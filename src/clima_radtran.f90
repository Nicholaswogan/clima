module clima_radtran
  use clima_const, only: dp
  use clima_types, only: ClimaVars
  implicit none
  
contains
  
  ! v: nz, dz, u0, surface_albedo, photons_sol, diurnal_fac
  ! rin: T, P, densities
  
  ! We should switch to this interface instead. It will make the radiative
  ! Transfer more generalized.
  
  ! subroutine radiate_explicit(op, kset, &
  !                    dz, surface_albedo, u0, diurnal_fac, photons_sol, &
  !                    T, P, densities, &
  !                    rw, rz, &
  !                    fup_a, fdn_a, fup_n, fdn_n)
  !   use clima_types, only: OpticalProperties
  !   use clima_types, only: RadiateXSWrk, RadiateZWrk, Ksettings
  ! 
  !   type(OpticalProperties), intent(inout) :: op !! Optical properties
  !   type(Ksettings), intent(in) :: kset !! Settings for how to deal with k-coefficients
  ! 
  !   real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
  !   real(dp), intent(in) :: surface_albedo !! Surface albedo
  !   real(dp), optional, intent(in) :: u0 !! Cosine of solar zenith angle
  !   real(dp), optional, intent(in) :: diurnal_fac !! Diurnal averaging factor (0.5)
  !   real(dp), optional, intent(in) :: photons_sol(:) !! (nw) Solar flux (mW/m2)
  ! 
  !   real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
  !   real(dp), intent(in) :: P(:) !! (nz) Pressure (bar) 
  !   real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
  !                                          !! molecule in each layer (molcules/cm3)
  ! 
  !   type(RadiateXSWrk), intent(inout) :: rw !! work arrays
  !   type(RadiateZWrk), intent(inout) :: rz !! work arrays
  ! 
  !   real(dp), intent(out) :: fup_a(:,:), fdn_a(:,:) !! (nz+1,nw)
  !   real(dp), intent(out) :: fup_n(:), fdn_n(:) !! (nz+1)
  ! 
  !   integer :: nz
  ! 
  !   nz = size(dz)
  ! 
  ! end subroutine
             
    
  subroutine radiate(op, kset, v, rin, rw, rz, fup_a, fdn_a, fup_n, fdn_n)
    use clima_types, only: RadiateXSWrk, RadiateZWrk, Ksettings, RadiateInputs
    use clima_types, only: OpticalProperties, Kcoefficients
    use clima_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_types, only: k_RandomOverlap, k_RandomOverlapResortRebin
    use clima_eqns, only: planck_fcn, ten2power
    
    type(OpticalProperties), intent(inout) :: op
    type(Ksettings), intent(in) :: kset
    type(ClimaVars), intent(in) :: v
    type(RadiateInputs), intent(in) :: rin
    type(RadiateXSWrk), intent(inout) :: rw
    type(RadiateZWrk), intent(inout) :: rz
    real(dp), intent(out) :: fup_a(:,:), fdn_a(:,:) ! (nz+1,nw)
    real(dp), intent(out) :: fup_n(:), fdn_n(:) ! (nz+1)
    
    real(dp), parameter :: max_w0 = 0.99999_dp
    
    integer :: i, j, k, l, n, jj
    real(dp) :: avg_freq
    
    ! array of indexes for recursive
    ! correlated-k
    integer :: iks(op%nk)
    
    ! other work
    real(dp) :: dfreq
    
    !$omp parallel private(i, j, k, l, n, jj, &
    !$omp& iks, &
    !$omp& rw, rz)
    
    !$omp do
    do l = 1,op%nw
      
      ! interpolate to T and P grid
      ! k-distributions
      do i = 1,op%nk
        do k = 1,op%k(i)%ngauss
          do j = 1,v%nz
            call op%k(i)%log10k(k,l)%evaluate(log10(rin%P(j)), rin%T(j), rw%ks(i)%k(j,k))
            rw%ks(i)%k(j,k) = ten2power(rw%ks(i)%k(j,k))
          enddo
        enddo
      enddo
    
      ! CIA
      do i = 1,op%ncia
        call interpolate_Xsection(op%cia(i), l, rin%P, rin%T, rw%cia(:,i))
      enddo
      
      ! Absorption xs
      do i = 1,op%naxs
        call interpolate_Xsection(op%axs(i), l, rin%P, rin%T, rw%axs(:,i))
      enddo
      
      ! Photolysis xs
      do i = 1,op%npxs
        call interpolate_Xsection(op%pxs(i), l, rin%P, rin%T, rw%pxs(:,i))
      enddo
    
      ! compute tau
      ! rayleigh scattering
      rz%tausg(:) = 0.0_dp
      do i = 1,op%nray
        j = op%ray(i)%sp_ind(1)
        do k = 1,v%nz
          n = v%nz+1-k
          rz%tausg(n) = rz%tausg(n) + op%ray(i)%xs_0d(l)*rin%cols(k,j)
        enddo
      enddo
      
      ! CIA
      rz%taua(:) = 0.0_dp
      do i = 1,op%ncia
        j = op%cia(i)%sp_ind(1)
        jj = op%cia(i)%sp_ind(2)
        do k = 1,v%nz
          n = v%nz+1-k
          rz%taua(n) = rz%taua(n) + rw%cia(k,i)*rin%densities(k,j)*rin%densities(k,jj)*v%dz(k)
        enddo
      enddo
      
      ! absorption
      do i = 1,op%naxs
        j = op%cia(i)%sp_ind(1)
        do k = 1,v%nz
          n = v%nz+1-k
          rz%taua(n) = rz%taua(n) + rw%axs(k,i)*rin%cols(k,j)
        enddo
      enddo
      
      ! photolysis
      do i = 1,op%npxs
        j = op%cia(i)%sp_ind(1)
        do k = 1,v%nz
          n = v%nz+1-k
          rz%taua(n) = rz%taua(n) + rw%pxs(k,i)*rin%cols(k,j)
        enddo
      enddo
      
      ! plank function, only if in the IR
      ! bplanck has units [W sr^−1 m^−2 Hz^-1]
      if (op%op_type == IROpticalProperties) then
        avg_freq = 0.5_dp*(op%freq(l) + op%freq(l+1))
        rz%bplanck(v%nz+1) = planck_fcn(avg_freq, rin%T(1)) ! ground level
        do j = 1,v%nz
          n = v%nz+1-j
          rz%bplanck(n) = planck_fcn(avg_freq, rin%T(j))
        enddo
      endif
      
      ! asymetry factor
      rz%gt = 0.0_dp
      
      rz%fup1 = 0.0_dp
      rz%fdn1 = 0.0_dp
      
      if (op%nk /= 0) then
        
        ! Deal with k-distributions
        if (kset%k_method == k_RandomOverlap) then
          ! Random Overlap method. Slow for >2 k species.
          ! k_loops uses recursion to make an arbitrary
          ! number of nested loops.
          call k_loops(v, op, rw%ks, rin, rz, iks, 1)
        elseif (kset%k_method == k_RandomOverlapResortRebin) then
          ! Random Overlap with Resorting and Rebinning.
          call k_rorr(v, op, kset, rin, rw, rz)
        endif
        
      else
        ! there are no k-distributions
        print*,"there are no k-distributions"
        stop 1
      endif ! endif k-dist
      
      fup_a(:,l) =  rz%fup1
      fdn_a(:,l) =  rz%fdn1
      
    enddo
    !$omp enddo
    !$omp end parallel
    
    ! Integrate fluxes over wavelength grid. Move values to grid-center,
    ! by averaging over edge values.
    fup_n = 0.0_dp
    fdn_n = 0.0_dp

    ! THIS IS THE RIGHT UNITS FOR fup_n, and fdn_n
    ! mW/m2
    if (op%op_type == FarUVOpticalProperties .or. &
        op%op_type == SolarOpticalProperties) then
        
      do l = 1,op%nw
        do i = 1,v%nz+1
          n = v%nz+2-i
          fup_n(i) = fup_n(i) + fup_a(n,l)*v%photons_sol(l)*v%diurnal_fac
          fdn_n(i) = fdn_n(i) + fdn_a(n,l)*v%photons_sol(l)*v%diurnal_fac
        enddo
      enddo
      
    elseif (op%op_type == IROpticalProperties) then
      
      do l = 1,op%nw
        dfreq = op%freq(l)-op%freq(l+1)
        do i = 1,v%nz+1
          n = v%nz+2-i
          fup_n(i) = fup_n(i) + fup_a(n,l)*dfreq*1.0e3_dp
          fdn_n(i) = fdn_n(i) + fdn_a(n,l)*dfreq*1.0e3_dp
        enddo
      enddo
      
    endif
    
  end subroutine
  
  subroutine k_rorr(v, op, kset, rin, rw, rz)
    use futils, only: rebin
    use mrgrnk_mod, only: mrgrnk
    
    use clima_types, only: OpticalProperties, Ksettings
    use clima_types, only: ClimaVars, RadiateZWrk, RadiateXSWrk, RadiateInputs
    use clima_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_eqns, only: weights_to_bins
    use clima_twostream, only: two_stream_solar, two_stream_ir
    
    type(ClimaVars), intent(in) :: v
    type(OpticalProperties), intent(in) :: op
    type(Ksettings), intent(in) :: kset
    type(RadiateInputs), target, intent(in) :: rin
    type(RadiateXSWrk), target, intent(in) :: rw
    type(RadiateZWrk), intent(inout) :: rz
    
    real(dp), pointer :: tau_k(:,:) ! (nz,nbin)
    real(dp), pointer :: tau_xy(:,:) ! (nz,nbin*ngauss_max)
    real(dp), pointer :: wxy(:) ! (nbin*ngauss_max)
    real(dp), pointer :: wxy1(:) ! (nbin*ngauss_max)
    real(dp), pointer :: wxy_e(:) ! (nbin*ngauss_max+1)
    integer, pointer :: inds(:) ! (nbin*ngauss_max)
      
    integer :: i, j, k, n, jj
    integer :: j1, j2, ngauss, ngauss_mix
    real(dp) :: surf_rad
      
    ! pointers to work arrays
    tau_k => rw%tau_k
    tau_xy => rw%tau_xy
    wxy => rw%wxy
    wxy1 => rw%wxy1
    wxy_e => rw%wxy_e
    inds => rw%inds
    
    ! rebin k-coeffs of first species to new grid
    do i = 1,v%nz
      call rebin(op%k(1)%weight_e, rw%ks(1)%k(i,:), kset%wbin_e, tau_k(i,:))
    enddo
    ! combine k-coefficients with species concentrations (molecules/cm2)
    j1 = op%k(1)%sp_ind
    do i = 1,kset%nbin
      tau_k(:,i) = tau_k(:,i)*rin%cols(:,j1)
    enddo
    
    ! Mix rest of k-coeff species with the first species
    do jj = 2,op%nk
      
      j2 = op%k(jj)%sp_ind
      ngauss = op%k(jj)%ngauss
      ngauss_mix = kset%nbin*ngauss
      
      do i = 1,kset%nbin
        do j = 1,ngauss
          tau_xy(:,j + (i-1)*ngauss) = tau_k(:,i) &
                                       + rw%ks(jj)%k(:,j)*rin%cols(:,j2)
          wxy(j + (i-1)*ngauss) = kset%wbin(i)*op%k(jj)%weights(j) 
        enddo
      enddo
      
      ! Here we only use the relevant space in 
      ! each work array.
      do i = 1,v%nz
        ! sort tau_xy and the weights
        wxy1(1:ngauss_mix) = tau_xy(i,1:ngauss_mix)
        call mrgrnk(wxy1(1:ngauss_mix), inds(1:ngauss_mix))
        tau_xy(i,1:ngauss_mix) = tau_xy(i,inds(1:ngauss_mix))
        wxy1(1:ngauss_mix) = wxy(inds(1:ngauss_mix))
        ! rebin to smaller grid
        call weights_to_bins(wxy1(1:ngauss_mix), wxy_e(1:ngauss_mix+1))
        call rebin(wxy_e(1:ngauss_mix+1), tau_xy(i,1:ngauss_mix), kset%wbin_e, tau_k(i,:))
      enddo
      
    enddo
    
    ! loop over mixed k-coeffs
    do i = 1,kset%nbin

      ! tau_k(:,i) is optical depth of ith mixed k-coeff.
      ! tau_k(1,i) is ground level. Need to reorder so that
      ! first element in array is top of the atmosphere.
      do k = 1,v%nz
        n = v%nz+1-k
        rz%taua_1(n) = tau_k(k,i)
      enddo
      
      ! sum all optical depths
      ! total = gas scattering + continumm/gray opacities + k-coeff
      rz%tau = rz%tausg + rz%taua + rz%taua_1
      do j = 1,v%nz
        rz%w0(j) = min(0.99999_dp,rz%tausg(j)/rz%tau(j))
      enddo
      
      if (op%op_type == FarUVOpticalProperties .or. &
          op%op_type == SolarOpticalProperties) then
        call two_stream_solar(v%nz, rz%tau, rz%w0, rz%gt, v%u0, v%surface_albedo, &
                              rz%amean, surf_rad, rz%fup, rz%fdn)
      elseif (op%op_type == IROpticalProperties) then
        call two_stream_ir(v%nz, rz%tau, rz%w0, rz%gt, v%surface_albedo, rz%bplanck, &
                           rz%fup, rz%fdn)
      endif
      
      ! weight upward and downward fluxes by
      ! k-coeff weights (kset%wbin(i))
      rz%fup1 = rz%fup1 + rz%fup*kset%wbin(i)
      rz%fdn1 = rz%fdn1 + rz%fdn*kset%wbin(i)
      
    enddo
    
  end subroutine
  
  recursive subroutine k_loops(v, op, ks, rin, rz, iks, ik)
    use clima_types, only: OpticalProperties, Kcoefficients
    use clima_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_types, only: ClimaVars, RadiateZWrk, RadiateInputs
    
    use clima_twostream, only: two_stream_solar, two_stream_ir
    
    type(ClimaVars), intent(in) :: v
    type(OpticalProperties), intent(in) :: op
    type(Kcoefficients), intent(in) :: ks(:)
    type(RadiateInputs), intent(in) :: rin
    type(RadiateZWrk), intent(inout) :: rz
    
    integer, intent(inout) :: iks(:)
    integer, intent(in) :: ik
    
    real(dp), parameter :: max_w0 = 0.99999_dp
    
    real(dp) :: gauss_weight, surf_rad
    integer :: i, j, k, n
    
    
    do i = 1,op%k(ik)%ngauss
      iks(ik) = i
      
      ! if ik = number of k coeffients
      ! then we are in the lowest loop
      ! and can actually do radiative transfer
      if (ik == op%nk) then
        ! do radiative transfer
        rz%taua_1(:) = 0.0_dp
        
        do k = 1,v%nz
          n = v%nz+1-k
          do j = 1,op%nk
            rz%taua_1(n) = rz%taua_1(n) + &
                           ks(j)%k(k,iks(j))*rin%cols(k,op%k(j)%sp_ind)
          enddo
        enddo
        
        ! sum
        rz%tau = rz%tausg + rz%taua + rz%taua_1
        do j = 1,v%nz
          rz%w0(j) = min(max_w0,rz%tausg(j)/rz%tau(j))
        enddo
        
        if (op%op_type == FarUVOpticalProperties .or. &
            op%op_type == SolarOpticalProperties) then
          call two_stream_solar(v%nz, rz%tau, rz%w0, rz%gt, v%u0, v%surface_albedo, &
                                rz%amean, surf_rad, rz%fup, rz%fdn)
        elseif (op%op_type == IROpticalProperties) then
          call two_stream_ir(v%nz, rz%tau, rz%w0, rz%gt, v%surface_albedo, rz%bplanck, &
                             rz%fup, rz%fdn)
        endif
        
        gauss_weight = 1.0_dp
        do j = 1,op%nk
          gauss_weight = gauss_weight*op%k(j)%weights(iks(j))
        enddo
        
        rz%fup1 = rz%fup1 + rz%fup*gauss_weight
        rz%fdn1 = rz%fdn1 + rz%fdn*gauss_weight
        
      else
        ! go into a deeper loop
        call k_loops(v, op, ks, rin, rz, iks, ik + 1)
      
      endif
      
    enddo
    
  end subroutine
  
  subroutine interpolate_Xsection(xs, l, P, T, res)
    use clima_types, only: Xsection
    use clima_eqns, only: ten2power
    
    type(Xsection), intent(inout) :: xs
    integer, intent(in) :: l
    real(dp), intent(in) :: P(:), T(:)
    real(dp), intent(out) :: res(size(P))
    
    integer :: j
    real(dp) :: val

    if (xs%dim == 0) then
      do j = 1,size(P)
        res(j) = xs%xs_0d(l)
      enddo
    elseif (xs%dim == 1) then
      do j = 1,size(P)
        call xs%log10_xs_1d(l)%evaluate(T(j), val)
        res(j) = ten2power(val)
      enddo
    endif
    
  end subroutine
  
  
end module
  
  