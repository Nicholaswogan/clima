module clima_radtran_radiate
  use clima_const, only: dp
  implicit none

  real(dp), parameter :: max_w0 = 0.99999_dp
  real(dp), parameter :: max_gt = 0.999999_dp
  
contains
  
  
  !> Does Solar or IR radiative transfer, depending on value of 
  !> `op%op_type`. `u0`, `diurnal_fac` and `photons_sol` are all
  !> only used during Solar radiative transfer.
  function radiate(op, &
                   surface_emissivity, &
                   surface_albedo, diurnal_fac, photons_sol, &
                   zenith_u, zenith_weights, &
                   P, T_surface, T, densities, dz, &
                   pdensities, radii, &
                   rw, rz, &
                   fup_a, fdn_a, fup_n, fdn_n, amean, tau_band) result(ierr)
    use clima_radtran_types, only: RadiateXSWrk, RadiateZWrk, Ksettings
    use clima_radtran_types, only: OpticalProperties, Kcoefficients
    use clima_radtran_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_radtran_types, only: k_RandomOverlap, k_RandomOverlapResortRebin, k_AdaptiveEquivalentExtinction
    use clima_eqns, only: planck_fcn, ten2power
    use clima_const, only: k_boltz, pi, plank, c_light
  
    type(OpticalProperties), target, intent(inout) :: op !! Optical properties
    
    real(dp), intent(in) :: surface_emissivity(:) !! Surface emissivity (Needed only for IR)

    real(dp), intent(in) :: surface_albedo(:) !! Surface albedo (Needed only for solar)
    real(dp), intent(in) :: diurnal_fac !! Diurnal averaging factor (0.5) (Needed only for solar)
    real(dp), intent(in) :: photons_sol(:) !! (nw) Average solar flux in each bin (mW/m2/Hz) (Needed only for solar)

    real(dp), intent(in) :: zenith_u(:) !! cos(zenith angles) to average stellar flux over.
                                        !! SHOULD BE SIZE 1 array for IR with any value.
    real(dp), intent(in) :: zenith_weights(:) !! Weights for integrating over zenith angle. 
                                              !! SHOULD BE SIZE 1 array for IR equal to 1.0.
  
    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: T_surface !! Surface Tempeature (K) (Needed only for solar)
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), optional, intent(in) :: pdensities(:,:) !! (nz,np) particles/cm3
    real(dp), optional, intent(in) :: radii(:,:) !! (nz,np) cm
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    
    type(RadiateXSWrk), intent(inout) :: rw !! work arrays
    type(RadiateZWrk), intent(inout) :: rz !! work arrays
  
    real(dp), intent(out) :: fup_a(:,:), fdn_a(:,:) !! (nz+1,nw) mW/m2/Hz in each wavelength bin
                                                    !! at the edges of the vertical grid
    real(dp), intent(out) :: fup_n(:), fdn_n(:) !! (nz+1) mW/m2 at the edges of the vertical grid 
                                                !! (integral of fup_a and fdn_a over wavelength grid)
    real(dp), intent(out) :: amean(:,:) !! (nz+1,nw) Mean intensity in photons/cm^2/s (Only for Solar)
    real(dp), intent(out) :: tau_band(:,:) !! (nz,nw) The optical depth of each layer
    integer :: ierr !! if ierr /= 0 on return, then there was an error

    type(Ksettings), pointer :: kset
    integer :: nz, ng
    integer :: i, j, k, l, n, jj
    integer :: ierr_0
    real(dp) :: avg_freq, avg_wavl
    
    ! array of indexes for recursive
    ! correlated-k
    integer :: iks(op%nk)
    
    ! other work
    real(dp) :: dfreq, TT, log10PP
    real(dp) :: albedo, emissivity
    real(dp), allocatable :: cols(:,:), log10P(:), log10P_cgs(:)
    real(dp), allocatable :: foreign_col(:)

    ! check to make sure the zenith anges and weights are OK
    ! for IR scenario.
    if (op%op_type == IROpticalProperties) then
      if (size(zenith_u) /= 1 .or. size(zenith_weights) /= 1) then
        ierr = -1
        return
      endif
      if (zenith_weights(1) /= 1.0_dp) then
        ierr = -1
        return
      endif
    endif
    
    kset => op%kset
    nz = size(dz)
    ng = size(densities, 2)
    allocate(log10P(nz))
    log10P = log10(P)
    log10P_cgs = log10(P*1.0e6_dp)
    allocate(cols(nz,ng))
    do i = 1,ng
      cols(:,i) = densities(:,i)*dz(:)
    enddo
    if (allocated(op%cont)) then
      allocate(foreign_col(nz))
      do j = 1,nz
        foreign_col(j) = 0.0_dp
        do i = 1,ng
          if (i /= op%cont%LH2O) then
            foreign_col(j) = foreign_col(j) + cols(j,i)
          endif
        enddo
      enddo
    endif

    ierr = 0
    
    !$omp parallel private(i, j, k, l, n, jj, &
    !$omp& iks, avg_freq, TT, log10PP, albedo, emissivity, &
    !$omp& ierr_0, &
    !$omp& rw, rz)
    
    ierr_0 = 0

    !$omp do
    do l = 1,op%nw
      
      ! interpolate to T and P grid
      ! k-distributions
      do i = 1,op%nk
        do k = 1,op%k(i)%ngauss
          do j = 1,nz
            
            if (T(j) < op%k(i)%T_min) then
              TT = op%k(i)%T_min
            elseif (T(j) > op%k(i)%T_max) then
              TT = op%k(i)%T_max
            else
              TT = T(j)
            endif
            
            if (log10P(j) < op%k(i)%log10P_min) then
              log10PP = op%k(i)%log10P_min
            elseif (log10P(j) > op%k(i)%log10P_max) then
              log10PP =  op%k(i)%log10P_max
            else
              log10PP = log10P(j)
            endif
  
            call op%k(i)%log10k(k,l)%evaluate(log10PP, TT, rw%ks(i)%k(j,k))
            rw%ks(i)%k(j,k) = ten2power(rw%ks(i)%k(j,k))
          enddo
        enddo
      enddo
    
      ! CIA
      do i = 1,op%ncia
        call interpolate_Xsection(op%cia(i), l, P, T, rw%cia(:,i))
      enddo
      
      ! Absorption xs
      do i = 1,op%naxs
        call interpolate_Xsection(op%axs(i), l, P, T, rw%axs(:,i))
      enddo
      
      ! Photolysis xs
      do i = 1,op%npxs
        call interpolate_Xsection(op%pxs(i), l, P, T, rw%pxs(:,i))
      enddo
      
      ! continuum absorption
      if (allocated(op%cont)) then
        call interpolate_WaterContinuum(op%cont, l, T, rw%H2O, rw%foreign)
      endif

      ! particles
      do i = 1,op%npart
        call interpolate_Particle(op%part(i), l, radii, rw%w0(:,i), rw%qext(:,i), rw%gt(:,i), ierr_0)
      enddo
      
      ! compute tau
      ! rayleigh scattering
      rz%tausg(:) = 0.0_dp
      do i = 1,op%nray
        j = op%ray(i)%sp_ind(1)
        do k = 1,nz
          n = nz+1-k
          rz%tausg(n) = rz%tausg(n) + op%ray(i)%xs_0d(l)*cols(k,j)
        enddo
      enddo
      
      ! CIA
      rz%taua(:) = 0.0_dp
      do i = 1,op%ncia
        j = op%cia(i)%sp_ind(1)
        jj = op%cia(i)%sp_ind(2)
        do k = 1,nz
          n = nz+1-k
          rz%taua(n) = rz%taua(n) + rw%cia(k,i)*densities(k,j)*densities(k,jj)*dz(k)
        enddo
      enddo
      
      ! absorption
      do i = 1,op%naxs
        j = op%axs(i)%sp_ind(1)
        do k = 1,nz
          n = nz+1-k
          rz%taua(n) = rz%taua(n) + rw%axs(k,i)*cols(k,j)
        enddo
      enddo
      
      ! photolysis
      do i = 1,op%npxs
        j = op%pxs(i)%sp_ind(1)
        do k = 1,nz
          n = nz+1-k
          rz%taua(n) = rz%taua(n) + rw%pxs(k,i)*cols(k,j)
        enddo
      enddo
      
      ! continuum absorption
      if (allocated(op%cont)) then
        do k = 1,nz
          n = nz+1-k
          rz%taua(n) = rz%taua(n) &
                       + rw%H2O(k)*densities(k,op%cont%LH2O)*cols(k,op%cont%LH2O) &
                       + rw%foreign(k)*densities(k,op%cont%LH2O)*foreign_col(k)
        enddo
      endif

      ! Custom opacity
      call op%custom_optical_properties(log10P_cgs, dz, l, rz%tauc, rz%w0c, rz%g0c)
      rz%tauc = rz%tauc(nz:1:-1)
      rz%w0c = rz%w0c(nz:1:-1)
      rz%g0c = rz%g0c(nz:1:-1)
      rz%tausc = rz%w0c*rz%tauc

      ! particles
      block
      real(dp) :: taup_1

      rz%tausp(:) = 0.0_dp
      rz%taup(:) = 0.0_dp
      do i = 1,op%npart
        j = op%part(i)%p_ind
        do k = 1,nz
          n = nz+1-k
          taup_1 = rw%qext(k,i)*pi*radii(k,j)**2.0_dp*pdensities(k,j)*dz(k)
          rz%taup(n) = rz%taup(n) + taup_1
          rz%tausp_1(n,i) = rw%w0(k,i)*taup_1
          rz%tausp(n) = rz%tausp(n) + rz%tausp_1(n,i)
        enddo
      enddo

      rz%gt(:) = 0.0_dp
      do i = 1,op%npart
        do k = 1,nz
          n = nz+1-k
          rz%gt(n) = rz%gt(n) + rw%gt(k,i)*rz%tausp_1(n,i)/(rz%tausp(n) + rz%tausg(n) + rz%tausc(n))
        enddo
      enddo
      rz%gt = rz%gt + rz%g0c*rz%tausc/(rz%tausp + rz%tausg + rz%tausc)
      do k = 1,nz
        n = nz+1-k
        rz%gt(n) = min(rz%gt(n), max_gt)
      enddo

      end block

      ! plank function, only if in the IR
      ! bplanck has units [mW sr^−1 m^−2 Hz^-1]
      if (op%op_type == IROpticalProperties) then
        avg_freq = 0.5_dp*(op%freq(l) + op%freq(l+1))
        rz%bplanck(nz+1) = planck_fcn(avg_freq, T_surface)
        do j = 1,nz
          n = nz+1-j
          rz%bplanck(n) = planck_fcn(avg_freq, T(j))
        enddo
        emissivity = surface_emissivity(l)
        albedo = 0.0_dp
      elseif (op%op_type == FarUVOpticalProperties .or. &
              op%op_type == SolarOpticalProperties) then
        emissivity = 0.0_dp
        albedo = surface_albedo(l)
      endif
      
      rz%fup2 = 0.0_dp
      rz%fdn2 = 0.0_dp
      rz%amean2 = 0.0_dp
      
      if (op%nk /= 0) then
        
        ! Deal with k-distributions
        if (kset%k_method == k_RandomOverlap) then
          ! Random Overlap method. Slow for >2 k species.
          ! k_loops uses recursion to make an arbitrary
          ! number of nested loops.
          do j = 1,size(zenith_u)
            rz%fup1 = 0.0_dp
            rz%fdn1 = 0.0_dp
            rz%amean1 = 0.0_dp
            rz%tau_band = 0.0_dp
            call k_loops(op, zenith_u(j), albedo, emissivity, cols, rw%ks, rz, iks, 1)
            rz%fup2(:) = rz%fup2(:) + rz%fup1(:)*zenith_weights(j)
            rz%fdn2(:) = rz%fdn2(:) + rz%fdn1(:)*zenith_weights(j)
            if (op%op_type == FarUVOpticalProperties .or. &
                op%op_type == SolarOpticalProperties) then
              rz%amean2 = rz%amean2 + rz%amean1*zenith_weights(j)
            endif
          enddo
        elseif (kset%k_method == k_RandomOverlapResortRebin) then
          ! Random Overlap with Resorting and Rebinning.
          call k_rorr(op, kset, zenith_u, zenith_weights, albedo, emissivity, cols, rw, rz)
        elseif (kset%k_method == k_AdaptiveEquivalentExtinction) then
          ! Adaptive Equivalent Extinction method
          call k_aee(op, kset, zenith_u, zenith_weights, albedo, emissivity, cols, rw, rz)
        endif
        
      else
        ! there are no k-distributions
        print*,"there are no k-distributions"
        stop 1
      endif ! endif k-dist
      
      ! save results. Here I reverse order so that
      ! fup_a(1,l) is ground level.
      do i = 1,nz+1
        n = nz+2-i
        fup_a(i,l) = rz%fup2(n)
        fdn_a(i,l) = rz%fdn2(n)
      enddo
      if (op%op_type == FarUVOpticalProperties .or. &
          op%op_type == SolarOpticalProperties) then
        do i = 1,nz+1
          n = nz+2-i
          amean(i,l) = rz%amean2(n)
        enddo
      endif
      do i = 1,nz
        n = nz+1-i
        tau_band(i,l) = rz%tau_band(n) ! band optical thickness
      enddo
      
    enddo
    !$omp enddo
    !$omp critical
    ierr = ierr + ierr_0
    !$omp end critical
    !$omp end parallel
    
    ! In Solar case, units for fup_a and fdn_a are unit-less.
    ! need to multiply by photons_sol (mW/m2/Hz). Then
    ! fup_a and fdn_a are in units mW/m2/Hz.
    if (op%op_type == FarUVOpticalProperties .or. &
        op%op_type == SolarOpticalProperties) then
      do l = 1,op%nw
        fup_a(:,l) = fup_a(:,l)*photons_sol(l)*diurnal_fac
        fdn_a(:,l) = fdn_a(:,l)*photons_sol(l)*diurnal_fac
        amean(:,l) = amean(:,l)*photons_sol(l)*diurnal_fac

        ! Convert from from mW/m^2/Hz to mW/m^2/nm
        avg_freq = 0.5_dp*(op%freq(l) + op%freq(l+1))
        avg_wavl = 1.0e9_dp*c_light/avg_freq ! nm
        amean(:,l) = amean(:,l)*(avg_freq/avg_wavl)
        ! Convert from mW/m^2/nm to photons/cm^2/s
        amean(:,l) = amean(:,l)*(avg_wavl/(plank*c_light*1.0e16_dp))*(op%wavl(l+1) - op%wavl(l))
      enddo
    endif
    
    ! Integrate fluxes over wavelength or frequency grid.
    ! Units for fup_n and fdn_n are mW/m^2.
    fup_n = 0.0_dp
    fdn_n = 0.0_dp
    do l = 1,op%nw
      dfreq = op%freq(l)-op%freq(l+1)
      do i = 1,nz+1
        fup_n(i) = fup_n(i) + fup_a(i,l)*dfreq
        fdn_n(i) = fdn_n(i) + fdn_a(i,l)*dfreq
      enddo
    enddo
  
  end function

  subroutine k_aee(op, kset, zenith_u, zenith_weights, surface_albedo, surface_emissivity, cols, rw, rz)
    use clima_radtran_types, only: OpticalProperties, Ksettings
    use clima_radtran_types, only: RadiateZWrk, RadiateXSWrk
    use clima_radtran_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_eqns, only: weights_to_bins
    use clima_radtran_twostream, only: two_stream_solar, two_stream_ir
    
    type(OpticalProperties), intent(in) :: op
    type(Ksettings), intent(in) :: kset
    real(dp), intent(in) :: zenith_u(:), zenith_weights(:)
    real(dp), intent(in) :: surface_albedo, surface_emissivity
    real(dp), intent(in) :: cols(:,:)
    type(RadiateXSWrk), target, intent(inout) :: rw
    type(RadiateZWrk), intent(inout) :: rz

    integer :: nz
    integer :: ii, i, j, k, n
    real(dp) :: surf_rad

    nz = size(cols, 1)

    ! we first need to compute the grey opacity of all k-coefficients
    rw%tau_grey(:,:) = 0.0_dp
    do i = 1,op%nk
      ii = op%k(i)%sp_ind
      do k = 1,op%k(i)%ngauss
        do j = 1,nz
          rw%tau_grey(j,i) = rw%tau_grey(j,i) + rw%ks(i)%k(j,k)*op%k(i)%weights(k)*cols(j,ii)
        enddo
      enddo
    enddo

    ! determine which species is the major absorber
    do j = 1,nz
      n = nz+1-j
      rw%ind_major(j) = maxloc(rw%tau_grey(j,:), 1)
      rw%tau_grey_sum(n) = sum(rw%tau_grey(j,:)) - rw%tau_grey(j,rw%ind_major(j))
    enddo

    ! integration over multiple zenith angles
    do ii = 1,size(zenith_u)
    rz%fup1 = 0.0_dp
    rz%fdn1 = 0.0_dp
    rz%amean1 = 0.0_dp
    rz%tau_band = 0.0_dp
    
    do i = 1,op%k(1)%ngauss
      ! This will only work if all ktables have the exact same number of
      ! gausian points, and the same weightings. This is checked while
      ! creating optical properties.

      ! Consider the k-coeffs from major species
      rz%taua_1(:) = 0.0_dp
      do k = 1,nz
        n = nz+1-k
        j = rw%ind_major(k)
        rz%taua_1(n) = rz%taua_1(n) + &
                       rw%ks(j)%k(k,i)*cols(k,op%k(j)%sp_ind)
      enddo

      ! sum
      rz%tau(:) = rz%tausg(:) + rz%taua(:) + rz%taup(:) + rz%taua_1(:) &
                  + rw%tau_grey_sum(:) + rz%tauc(:)
      do j = 1,nz
        rz%w0(j) = min(max_w0,(rz%tausg(j) + rz%tausp(j) + rz%tausc(j))/rz%tau(j))
      enddo

      if (op%op_type == FarUVOpticalProperties .or. &
          op%op_type == SolarOpticalProperties) then
        call two_stream_solar(nz, rz%tau, rz%w0, rz%gt, zenith_u(ii), surface_albedo, &
                              rz%amean, surf_rad, rz%fup, rz%fdn)
      elseif (op%op_type == IROpticalProperties) then
        call two_stream_ir(nz, rz%tau, rz%w0, rz%gt, surface_emissivity, rz%bplanck, &
                           rz%fup, rz%fdn)
      endif
      
      ! weight upward and downward fluxes by
      ! k-coeff weights (kset%wbin(i))
      rz%fup1(:) = rz%fup1(:) + rz%fup(:)*op%k(1)%weights(i)
      rz%fdn1(:) = rz%fdn1(:) + rz%fdn(:)*op%k(1)%weights(i)
      if (op%op_type == FarUVOpticalProperties .or. &
          op%op_type == SolarOpticalProperties) then
        rz%amean1 = rz%amean1 + rz%amean(:)*op%k(1)%weights(i)
      endif

      ! Save the optical thickness in the band.
      rz%tau_band(:) = rz%tau_band(:) + rz%tau(:)*op%k(1)%weights(i)

    enddo

    rz%fup2 = rz%fup2 + rz%fup1*zenith_weights(ii)
    rz%fdn2 = rz%fdn2 + rz%fdn1*zenith_weights(ii)
    if (op%op_type == FarUVOpticalProperties .or. &
          op%op_type == SolarOpticalProperties) then
      rz%amean2 = rz%amean2 + rz%amean1*zenith_weights(ii)
    endif
    enddo

  end subroutine
  
  subroutine k_rorr(op, kset, zenith_u, zenith_weights, surface_albedo, surface_emissivity, cols, rw, rz)
    use futils, only: rebin
    use futils_mrgrnk, only: mrgrnk
    
    use clima_radtran_types, only: OpticalProperties, Ksettings
    use clima_radtran_types, only: RadiateZWrk, RadiateXSWrk
    use clima_radtran_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_eqns, only: weights_to_bins
    use clima_radtran_twostream, only: two_stream_solar, two_stream_ir
    
    type(OpticalProperties), intent(in) :: op
    type(Ksettings), intent(in) :: kset
    real(dp), intent(in) :: zenith_u(:), zenith_weights(:)
    real(dp), intent(in) :: surface_albedo, surface_emissivity
    real(dp), intent(in) :: cols(:,:)
    type(RadiateXSWrk), target, intent(in) :: rw
    type(RadiateZWrk), intent(inout) :: rz
    
    real(dp), pointer :: tau_k(:,:) ! (nz,nbin)
    real(dp), pointer :: tau_xy(:,:) ! (nz,nbin*ngauss_max)
    real(dp), pointer :: wxy(:) ! (nbin*ngauss_max)
    real(dp), pointer :: wxy1(:) ! (nbin*ngauss_max)
    real(dp), pointer :: wxy_e(:) ! (nbin*ngauss_max+1)
    integer, pointer :: inds(:) ! (nbin*ngauss_max)
      
    integer :: nz
    integer :: i, j, k, n, jj, ii
    integer :: j1, j2, ngauss, ngauss_mix
    real(dp) :: surf_rad
      
    nz = size(cols, 1)
    ! pointers to work arrays
    tau_k => rw%tau_k
    tau_xy => rw%tau_xy
    wxy => rw%wxy
    wxy1 => rw%wxy1
    wxy_e => rw%wxy_e
    inds => rw%inds
    
    ! rebin k-coeffs of first species to new grid
    do i = 1,nz
      call rebin(op%k(1)%weight_e, rw%ks(1)%k(i,:), kset%wbin_e, tau_k(i,:))
    enddo
    ! combine k-coefficients with species concentrations (molecules/cm2)
    j1 = op%k(1)%sp_ind
    do i = 1,kset%nbin
      tau_k(:,i) = tau_k(:,i)*cols(:,j1)
    enddo
    
    ! Mix rest of k-coeff species with the first species
    do jj = 2,op%nk
      
      j2 = op%k(jj)%sp_ind
      ngauss = op%k(jj)%ngauss
      ngauss_mix = kset%nbin*ngauss
      
      do i = 1,kset%nbin
        do j = 1,ngauss
          tau_xy(:,j + (i-1)*ngauss) = tau_k(:,i) &
                                       + rw%ks(jj)%k(:,j)*cols(:,j2)
          wxy(j + (i-1)*ngauss) = kset%wbin(i)*op%k(jj)%weights(j) 
        enddo
      enddo
      
      ! Here we only use the relevant space in 
      ! each work array.
      do i = 1,nz
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

    ! integration over multiple zenith angles
    do ii = 1,size(zenith_u)
    rz%fup1 = 0.0_dp
    rz%fdn1 = 0.0_dp
    rz%amean1 = 0.0_dp
    rz%tau_band = 0.0_dp
    
    ! loop over mixed k-coeffs
    do i = 1,kset%nbin

      ! tau_k(:,i) is optical depth of ith mixed k-coeff.
      ! tau_k(1,i) is ground level. Need to reorder so that
      ! first element in array is top of the atmosphere.
      do k = 1,nz
        n = nz+1-k
        rz%taua_1(n) = tau_k(k,i)
      enddo
      
      ! sum all optical depths
      ! total = gas scattering + continumm opacities + particle absorption + k-coeff + custom opacity
      rz%tau(:) = rz%tausg(:) + rz%taua(:) + rz%taup(:) + rz%taua_1(:) + rz%tauc(:)
      do j = 1,nz
        rz%w0(j) = min(max_w0,(rz%tausg(j) + rz%tausp(j) + rz%tausc(j))/rz%tau(j))
      enddo
      
      if (op%op_type == FarUVOpticalProperties .or. &
          op%op_type == SolarOpticalProperties) then
        call two_stream_solar(nz, rz%tau, rz%w0, rz%gt, zenith_u(ii), surface_albedo, &
                              rz%amean, surf_rad, rz%fup, rz%fdn)
      elseif (op%op_type == IROpticalProperties) then
        call two_stream_ir(nz, rz%tau, rz%w0, rz%gt, surface_emissivity, rz%bplanck, &
                           rz%fup, rz%fdn)
      endif
      
      ! weight upward and downward fluxes by
      ! k-coeff weights (kset%wbin(i))
      rz%fup1 = rz%fup1 + rz%fup*kset%wbin(i)
      rz%fdn1 = rz%fdn1 + rz%fdn*kset%wbin(i)
      if (op%op_type == FarUVOpticalProperties .or. &
          op%op_type == SolarOpticalProperties) then
        rz%amean1 = rz%amean1 + rz%amean*kset%wbin(i)
      endif

      ! Save the optical thickness in the band.
      rz%tau_band(:) = rz%tau_band(:) + rz%tau(:)*kset%wbin(i)
      
    enddo

    rz%fup2 = rz%fup2 + rz%fup1*zenith_weights(ii)
    rz%fdn2 = rz%fdn2 + rz%fdn1*zenith_weights(ii)
    if (op%op_type == FarUVOpticalProperties .or. &
        op%op_type == SolarOpticalProperties) then
      rz%amean2 = rz%amean2 + rz%amean1*zenith_weights(ii)
    endif
    enddo
    
  end subroutine
  
  recursive subroutine k_loops(op, u0, surface_albedo, surface_emissivity, cols, ks, rz, iks, ik)
    use clima_radtran_types, only: OpticalProperties, Kcoefficients
    use clima_radtran_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_radtran_types, only: RadiateZWrk
    
    use clima_radtran_twostream, only: two_stream_solar, two_stream_ir
    
    type(OpticalProperties), intent(in) :: op
    real(dp), intent(in) :: u0
    real(dp), intent(in) :: surface_albedo, surface_emissivity
    real(dp), intent(in) :: cols(:,:)
    type(Kcoefficients), intent(in) :: ks(:)
    type(RadiateZWrk), intent(inout) :: rz
    
    integer, intent(inout) :: iks(:)
    integer, intent(in) :: ik
    
    real(dp) :: gauss_weight, surf_rad
    integer :: nz
    integer :: i, j, k, n
    
    nz = size(cols, 1)
    
    do i = 1,op%k(ik)%ngauss
      iks(ik) = i
      
      ! if ik = number of k coeffients
      ! then we are in the lowest loop
      ! and can actually do radiative transfer
      if (ik == op%nk) then
        ! do radiative transfer
        rz%taua_1(:) = 0.0_dp
        
        do k = 1,nz
          n = nz+1-k
          do j = 1,op%nk
            rz%taua_1(n) = rz%taua_1(n) + &
                           ks(j)%k(k,iks(j))*cols(k,op%k(j)%sp_ind)
          enddo
        enddo
        
        ! sum
        rz%tau(:) = rz%tausg(:) + rz%taua(:) + rz%taup(:) + rz%taua_1(:) + rz%tauc(:)
        do j = 1,nz
          rz%w0(j) = min(max_w0,(rz%tausg(j) + rz%tausp(j) + rz%tausc(j))/rz%tau(j))
        enddo
        
        if (op%op_type == FarUVOpticalProperties .or. &
            op%op_type == SolarOpticalProperties) then
          call two_stream_solar(nz, rz%tau, rz%w0, rz%gt, u0, surface_albedo, &
                                rz%amean, surf_rad, rz%fup, rz%fdn)
        elseif (op%op_type == IROpticalProperties) then
          call two_stream_ir(nz, rz%tau, rz%w0, rz%gt, surface_emissivity, rz%bplanck, &
                             rz%fup, rz%fdn)
        endif
        
        gauss_weight = 1.0_dp
        do j = 1,op%nk
          gauss_weight = gauss_weight*op%k(j)%weights(iks(j))
        enddo
        
        rz%fup1 = rz%fup1 + rz%fup*gauss_weight
        rz%fdn1 = rz%fdn1 + rz%fdn*gauss_weight
        if (op%op_type == FarUVOpticalProperties .or. &
          op%op_type == SolarOpticalProperties) then
          rz%amean1 = rz%amean1 + rz%amean*gauss_weight
        endif

        ! Save the optical thickness in the band.
        rz%tau_band(:) = rz%tau_band(:) + rz%tau(:)*gauss_weight
        
      else
        ! go into a deeper loop
        call k_loops(op, u0, surface_albedo, surface_emissivity, cols, ks, rz, iks, ik + 1)
      
      endif
      
    enddo
    
  end subroutine
  
  subroutine interpolate_Xsection(xs, l, P, T, res)
    use clima_radtran_types, only: Xsection
    use clima_eqns, only: ten2power
    
    type(Xsection), intent(inout) :: xs
    integer, intent(in) :: l
    real(dp), intent(in) :: P(:), T(:)
    real(dp), intent(out) :: res(size(P))
    
    integer :: j
    real(dp) :: val, TT

    if (xs%dim == 0) then
      do j = 1,size(T)
        res(j) = xs%xs_0d(l)
      enddo
    elseif (xs%dim == 1) then
      do j = 1,size(T)
        
        if (T(j) < xs%T_min) then
          TT = xs%T_min
        elseif (T(j) > xs%T_max) then
          TT = xs%T_max
        else
          TT = T(j)
        endif
        
        call xs%log10_xs_1d(l)%evaluate(TT, val)
        res(j) = ten2power(val)
      enddo
    endif
    
  end subroutine
  
  subroutine interpolate_WaterContinuum(cont, l, T, H2O, foreign)
    use clima_radtran_types, only: WaterContinuum
    use clima_eqns, only: ten2power
    
    type(WaterContinuum), intent(inout) :: cont
    integer, intent(in) :: l
    real(dp), intent(in) :: T(:)
    real(dp), intent(out) :: H2O(size(T))
    real(dp), intent(out) :: foreign(size(T))
    
    integer :: j
    real(dp) :: val, TT
    
    do j = 1,size(T)
      if (T(j) < cont%T_min) then
        TT = cont%T_min
      elseif (T(j) > cont%T_max) then
        TT = cont%T_max
      else
        TT = T(j)
      endif
      
      call cont%log10_xs_H2O(l)%evaluate(TT, val)
      H2O(j) = ten2power(val)
      
      call cont%log10_xs_foreign(l)%evaluate(TT, val)
      foreign(j) = ten2power(val)
    enddo
    
  end subroutine

  subroutine interpolate_Particle(part, l, radii, w0, qext, gt, ierr)
    use clima_radtran_types, only: ParticleXsection
    type(ParticleXsection), intent(inout) :: part
    integer, intent(in) :: l
    real(dp), intent(in) :: radii(:,:)
    real(dp), intent(out) :: w0(:)
    real(dp), intent(out) :: qext(:)
    real(dp), intent(out) :: gt(:)
    integer, intent(out) :: ierr

    integer :: j
    real(dp) :: rp

    ierr = 0

    do j = 1,size(radii,1)
      rp = radii(j,part%p_ind)

      ! make sure radii is within data
      if (rp < part%r_min) then
        ierr = 1
      elseif (rp > part%r_max) then
        ierr = 1
      endif

      call part%w0(l)%evaluate(rp, w0(j))
      call part%qext(l)%evaluate(rp, qext(j))
      call part%gt(l)%evaluate(rp, gt(j))
    enddo

  end subroutine
  
  
end module
  
  