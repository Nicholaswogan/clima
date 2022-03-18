module clima_climate
  use clima_const, only: dp
  use clima_types, only: ClimaData, ClimaVars, ClimaWrk
  implicit none
  private
  
  public :: Climate, radiative_transfer
  
  type :: Climate
    type(ClimaData) :: d
    type(ClimaVars) :: v
    type(ClimaWrk) :: w
  end type
  
  interface Climate
    module procedure :: create_Climate
  end interface
  
contains
  
  function create_Climate(data_dir, species_f, settings_f, star_f, atmosphere_f, err) result(c)
    use clima_const, only: k_boltz
    use clima_types, only: ClimaSettings
    use clima_input, only: create_ClimaVars, create_ClimaData, create_ClimaSettings
    character(*), intent(in) :: data_dir
    character(*), intent(in) :: species_f
    character(*), intent(in) :: settings_f
    character(*), intent(in) :: star_f
    character(*), intent(in) :: atmosphere_f
    character(:), allocatable, intent(out) :: err
    
    type(Climate) :: c
    
    type(ClimaSettings) :: s
    integer :: i
    
    s = create_ClimaSettings(settings_f, err)
    if (allocated(err)) return
    c%d = create_ClimaData(species_f, data_dir, s, err)
    if (allocated(err)) return
    c%v = create_ClimaVars(atmosphere_f, star_f, c%d, err)
    if (allocated(err)) return
    
    ! after file read-in
    allocate(c%v%density(c%v%nz))
    allocate(c%v%densities(c%v%nz,c%d%ng))
    allocate(c%v%cols(c%v%nz,c%d%ng))
    c%v%density = 1.0e6_dp*c%v%P/(k_boltz*c%v%T)
    do i = 1,c%d%ng
      c%v%densities(:,i) = c%v%mix(:,i)*c%v%density
      c%v%cols(:,i) = c%v%densities(:,i)*c%v%dz(:)
    enddo
    
    ! work
    c%w = ClimaWrk(c%d, c%v%nz)
    
  end function
  
  subroutine radiative_transfer(d, v, w)
    type(ClimaData), intent(inout) :: d
    type(ClimaVars), intent(in) :: v
    type(ClimaWrk), intent(inout) :: w
    
    real(dp), allocatable :: fup_sol(:)
    real(dp), allocatable :: fdn_sol(:)
    
    real(dp), allocatable :: fup_ir(:)
    real(dp), allocatable :: fdn_ir(:)
    
    real(dp), allocatable :: fup_a(:,:)
    real(dp), allocatable :: fdn_a(:,:)
    
    real(dp), allocatable :: ftot(:)
    
    
    allocate(fup_sol(v%nz),fdn_sol(v%nz),fup_ir(v%nz),fdn_ir(v%nz),ftot(v%nz))
    
    allocate(fup_a(v%nz+1,d%ir%nw))
    allocate(fdn_a(v%nz+1,d%ir%nw))
    
    ! Solar radiative transfer
    call radiate(d%ir, v, w%rx_ir, w%rz, fup_a, fdn_a, fup_ir, fdn_ir)
    ! call radiate(d%sol, v, fup_sol, fdn_sol)
    

  end subroutine
  
  subroutine radiate(op, v, rw, rz, fup_a, fdn_a, fup_n, fdn_n)
    use clima_types, only: RadiateXSWrk, RadiateZWrk
    use clima_types, only: OpticalProperties, Kcoefficients
    use clima_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_types, only: k_RandomOverlap, k_RandomOverlapResortRebin
    use clima_eqns, only: planck_fcn
    
    use futils, only: Timer
    type(OpticalProperties), intent(inout) :: op
    type(ClimaVars), intent(in) :: v
    type(RadiateXSWrk), intent(inout) :: rw
    type(RadiateZWrk), intent(inout) :: rz
    real(dp), intent(out) :: fup_a(:,:), fdn_a(:,:) ! (nz+1,nw)
    real(dp), intent(out) :: fup_n(:), fdn_n(:) ! (nz)
    
    real(dp), parameter :: max_w0 = 0.99999_dp
    
    integer :: i, j, k, l, n, jj
    
    ! array of indexes for recursive
    ! correlated-k
    integer :: iks(op%nk)
    
    ! other work
    real(dp) :: dfreq
  
  
    type(Timer) :: tm
    
    call tm%start()

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
            call op%k(i)%log10k(k,l)%evaluate(log10(v%P(j)), v%T(j), rw%ks(i)%k(j,k))
            rw%ks(i)%k(j,k) = 10.0_dp**rw%ks(i)%k(j,k)
          enddo
        enddo
      enddo
    
      ! CIA
      do i = 1,op%ncia
        call interpolate_Xsection(op%cia(i), l, v%P, v%T, rw%cia(:,i))
      enddo
      
      ! Absorption xs
      do i = 1,op%naxs
        call interpolate_Xsection(op%axs(i), l, v%P, v%T, rw%axs(:,i))
      enddo
      
      ! Photolysis xs
      do i = 1,op%npxs
        call interpolate_Xsection(op%pxs(i), l, v%P, v%T, rw%pxs(:,i))
      enddo
    
      ! compute tau
      ! rayleigh scattering
      rz%tausg(:) = 0.0_dp
      do i = 1,op%nray
        j = op%ray(i)%sp_ind(1)
        do k = 1,v%nz
          n = v%nz+1-k
          rz%tausg(n) = rz%tausg(n) + op%ray(i)%xs_0d(l)*v%densities(k,j)*v%dz(k)
        enddo
      enddo
      
      ! CIA
      rz%taua(:) = 0.0_dp
      do i = 1,op%ncia
        j = op%cia(i)%sp_ind(1)
        jj = op%cia(i)%sp_ind(2)
        do k = 1,v%nz
          n = v%nz+1-k
          rz%taua(n) = rz%taua(n) + rw%cia(k,i)*v%densities(k,j)*v%densities(k,jj)*v%dz(k)
        enddo
      enddo
      
      ! absorption
      do i = 1,op%naxs
        j = op%cia(i)%sp_ind(1)
        do k = 1,v%nz
          n = v%nz+1-k
          rz%taua(n) = rz%taua(n) + rw%axs(k,i)*v%densities(k,j)*v%dz(k)
        enddo
      enddo
      
      ! photolysis
      do i = 1,op%npxs
        j = op%cia(i)%sp_ind(1)
        do k = 1,v%nz
          n = v%nz+1-k
          rz%taua(n) = rz%taua(n) + rw%pxs(k,i)*v%densities(k,j)*v%dz(k)
        enddo
      enddo
      
      ! plank function, only if in the IR
      ! bplanck has units [W sr^−1 m^−2 Hz^-1]
      if (op%op_type == IROpticalProperties) then
        rz%bplanck(v%nz+1) = planck_fcn(op%freq(l), v%T(1)) ! ground level
        do j = 1,v%nz
          n = v%nz+1-j
          rz%bplanck(n) = planck_fcn(op%freq(l), v%T(j))
        enddo
      endif
      
      ! asymetry factor
      rz%gt = 0.0_dp
      
      rz%fup1 = 0.0_dp
      rz%fdn1 = 0.0_dp
      
      if (v%k_method == k_RandomOverlap) then
        ! Random Overlap method. Very slow.
        ! k_loops uses recursion to make an arbitrary
        ! number of nested loops.
        call k_loops(v, op, rw%ks, rz, iks, 1)
      elseif (v%k_method == k_RandomOverlapResortRebin) then
        ! Random Overlap with Resorting and Rebinning.
        
        block
          use futils, only: rebin, argsort
          use clima_eqns, only: weights_to_bins
          use clima_twostream, only: two_stream_solar, two_stream_ir
          real(dp), allocatable :: tau_k(:,:) ! (nz,nbin)
          real(dp), allocatable :: tau_xy(:,:) ! (nz,nbin*ngauss_max)
          real(dp), allocatable :: wxy(:)
          real(dp), allocatable :: wxy1(:)
          real(dp), allocatable :: wxy_e(:)
          integer, allocatable :: inds(:)
          
          integer :: j1, j2, ngauss_max, ngauss, ngauss_mix
          real(dp) :: surf_rad
          
          
          allocate(tau_k(v%nz,v%nbin))
          ! find biggest ngauss
          ngauss_max = 0
          do i = 1,op%nk
            ngauss_max = max(ngauss_max,op%k(i)%ngauss)
          enddo
          ! allocate enough space in these arrays for mixing
          ! all k-coeffs
          allocate(tau_xy(v%nz,v%nbin*ngauss_max))
          allocate(wxy(v%nbin*ngauss_max))
          allocate(wxy1(v%nbin*ngauss_max))
          allocate(wxy_e(v%nbin*ngauss_max+1))
          allocate(inds(v%nbin*ngauss_max))
          
          ! rebin k-coeffs of first species to new grid
          do i = 1,v%nz
            call rebin(op%k(1)%weight_e, rw%ks(1)%k(i,:), v%wbin_e, tau_k(i,:))
          enddo
          ! combine k-coefficients with species concentrations (molecules/cm2)
          j1 = op%k(1)%sp_ind
          do i = 1,v%nbin
            tau_k(:,i) = tau_k(:,i)*v%cols(:,j1)
          enddo
          
          ! Mix rest of k-coeff species with the first species
          do jj = 2,op%nk
            
            j2 = op%k(jj)%sp_ind
            ngauss = op%k(jj)%ngauss
            ngauss_mix = v%nbin*ngauss
            
            do i = 1,v%nbin
              do j = 1,ngauss
                tau_xy(:,j + (i-1)*ngauss) = tau_k(:,i) &
                                                      + rw%ks(jj)%k(:,j)*v%cols(:,j2)
                wxy(j + (i-1)*ngauss) = v%wbin(i)*op%k(jj)%weights(j) 
              enddo
            enddo
            
            ! Here we only use the relevant space in 
            ! each work array.
            do i = 1,v%nz
              ! sort tau_xy and the weights
              inds(1:ngauss_mix) = argsort(tau_xy(i,1:ngauss_mix))
              tau_xy(i,1:ngauss_mix) = tau_xy(i,inds(1:ngauss_mix))
              wxy1(1:ngauss_mix) = wxy(inds(1:ngauss_mix))
              ! rebin to smaller grid
              call weights_to_bins(wxy1(1:ngauss_mix), wxy_e(1:ngauss_mix+1))
              call rebin(wxy_e(1:ngauss_mix+1), tau_xy(i,1:ngauss_mix), v%wbin_e, tau_k(i,:))
            enddo
            
          enddo
          
          ! loop over mixed k-coeffs
          do i = 1,v%nbin

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
            ! k-coeff weights (v%wbin(i))
            rz%fup1 = rz%fup1 + rz%fup*v%wbin(i)
            rz%fdn1 = rz%fdn1 + rz%fdn*v%wbin(i)
            
          enddo
                
          
        end block
        
      endif
      
      fup_a(:,l) =  max(rz%fup1,0.0_dp)
      fdn_a(:,l) =  max(rz%fdn1,0.0_dp)
      
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
        do i = 1,v%nz
          n = v%nz+1-i
          fup_n(i) = fup_n(i) + (0.5_dp*(fup_a(n,l)+fup_a(n+1,l)))*v%photons_sol(l)
          fdn_n(i) = fdn_n(i) + (0.5_dp*(fdn_a(n,l)+fdn_a(n+1,l)))*v%photons_sol(l)
        enddo
      enddo
      
    elseif (op%op_type == IROpticalProperties) then
      
      do l = 1,op%nw
        dfreq = op%freq(l)-op%freq(l+1)
        do i = 1,v%nz
          n = v%nz+1-i
          fup_n(i) = fup_n(i) + (0.5_dp*(fup_a(n,l)+fup_a(n+1,l)))*dfreq*1.0e3_dp
          fdn_n(i) = fdn_n(i) + (0.5_dp*(fdn_a(n,l)+fdn_a(n+1,l)))*dfreq*1.0e3_dp
        enddo
      enddo
      
    endif
    
    call tm%finish('')
    
    ! open(unit=1,file='../fup_rorr8.dat',form='formatted',status='replace')
    ! do i = 1,op%nw
    !   write(1,*) op%freq(i),op%wavl(i)*1.0e-3_dp,fup_a(1,i)
    ! enddo
    ! close(1)  
    
  end subroutine
  
  
  recursive subroutine k_loops(v, op, ks, rz, iks, ik)
    use clima_types, only: OpticalProperties, Kcoefficients
    use clima_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_types, only: ClimaVars, RadiateZWrk
    
    use clima_twostream, only: two_stream_solar, two_stream_ir
    
    type(ClimaVars), intent(in) :: v
    type(OpticalProperties), intent(in) :: op
    type(Kcoefficients), intent(in) :: ks(:)
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
                           ks(j)%k(k,iks(j))*v%densities(k,op%k(j)%sp_ind)*v%dz(k)
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
        call k_loops(v, op, ks, rz, iks, ik + 1)
      
      endif
      
    enddo
    
  end subroutine
  
  subroutine interpolate_Xsection(xs, l, P, T, res)
    use clima_types, only: Xsection
    
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
        res(j) = 10.0_dp**val
      enddo
    elseif (xs%dim == 2) then
      print*,'NOOOOO'
      stop 1
    endif
    
  end subroutine
  
  
end module