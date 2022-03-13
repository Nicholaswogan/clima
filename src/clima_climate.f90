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
    c%v%density = 1.0e6_dp*c%v%P/(k_boltz*c%v%T)
    do i = 1,c%d%ng
      c%v%densities(:,i) = c%v%mix(:,i)*c%v%density
    enddo
    
  end function
  
  subroutine radiative_transfer(d, v)
    type(ClimaData), intent(inout) :: d
    type(ClimaVars), intent(in) :: v
    
    real(dp), allocatable :: fup_sol(:)
    real(dp), allocatable :: fdn_sol(:)
    
    real(dp), allocatable :: fup_ir(:)
    real(dp), allocatable :: fdn_ir(:)
    
    real(dp), allocatable :: ftot(:)
    
    integer :: i
    
    allocate(fup_sol(v%nz),fdn_sol(v%nz),fup_ir(v%nz),fdn_ir(v%nz),ftot(v%nz))
    
    ! Solar radiative transfer
    call radiate(d%ir, v, fup_ir, fdn_ir)
    call radiate(d%sol, v, fup_sol, fdn_sol)
    
    ftot = (fup_sol-fdn_sol) + (fup_ir-fdn_ir)

    
    
  end subroutine
  
  subroutine radiate(op, v, fup_n, fdn_n)
    use clima_const, only: pi
    use clima_types, only: OpticalProperties, Kcoefficients
    use clima_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_twostream, only: two_stream_solar, two_stream_ir
    use futils, only: Timer
    type(OpticalProperties), intent(inout) :: op
    type(ClimaVars), intent(in) :: v
    real(dp), intent(inout) :: fup_n(:), fdn_n(:)
    
    real(dp), parameter :: max_w0 = 0.99999e0_dp
    
    integer :: i, j, k, l, n, jj
    integer :: i1, i2, i3
    integer :: j1, j2, j3
    real(dp) :: u0
    real(dp) :: surface_albedo
    real(dp) :: gauss_weight, surf_rad, weights, dfreq
    
    ! work
    real(dp) :: val
    type(Kcoefficients), allocatable :: ks(:)
    real(dp), allocatable :: cia(:,:)
    real(dp), allocatable :: axs(:,:)
    real(dp), allocatable :: pxs(:,:)
    
    real(dp), allocatable :: tausg(:)
    real(dp), allocatable :: taua(:)
    real(dp), allocatable :: taua_1(:)
    real(dp), allocatable :: tau(:)
    real(dp), allocatable :: w0(:)
    real(dp), allocatable :: gt(:)
    real(dp), allocatable :: amean(:)
    real(dp), allocatable :: fup(:)
    real(dp), allocatable :: fdn(:)
    real(dp), allocatable :: bplanck(:)
    
    real(dp), allocatable :: fup1(:)
    real(dp), allocatable :: fdn1(:)
    real(dp), allocatable :: fup_a(:,:)
    real(dp), allocatable :: fdn_a(:,:)
    
    allocate(ks(op%nk))
    do i = 1,op%nk
      allocate(ks(i)%k(v%nz,op%k(i)%ngauss))
    enddo
    
    allocate(cia(v%nz,op%ncia))
    allocate(axs(v%nz,op%naxs))
    allocate(pxs(v%nz,op%npxs))
    allocate(tausg(v%nz))
    allocate(taua(v%nz))
    allocate(taua_1(v%nz))
    allocate(tau(v%nz))
    allocate(w0(v%nz))
    allocate(gt(v%nz))
    allocate(amean(v%nz+1))
    allocate(fup(v%nz+1))
    allocate(fdn(v%nz+1))
    allocate(bplanck(v%nz+1))
    
    allocate(fup1(v%nz+1))
    allocate(fdn1(v%nz+1))
    
    allocate(fup_a(v%nz+1,op%nw))
    allocate(fdn_a(v%nz+1,op%nw))
    
    ! solar zenith angle
    u0 = cos(50.0_dp*pi/180.0_dp)
    surface_albedo = 0.25_dp
    
    !$omp parallel private(i, j, k, l, n, jj, &
    !$omp& i1, i2, i3, j1, j2, j3, surf_rad, val, &
    !$omp& ks, cia, axs, pxs, tausg, taua, taua_1, &
    !$omp& tau, w0, gt, amean, fup, fdn, bplanck, fup1, fdn1)
    
    !$omp do
    do l = 1,op%nw
      
      fup1 = 0.0_dp
      fdn1 = 0.0_dp
    
      ! interpolate to T and P grid
      ! k-distributions
      do i = 1,op%nk
        ks(i)%ngauss = op%k(i)%ngauss
        do k = 1,op%k(i)%ngauss
          do j = 1,v%nz
            call op%k(i)%log10k(k,l)%evaluate(log10(v%P(j)), v%T(j), val)
            ks(i)%k(j,k) = 10.0_dp**val
          enddo
        enddo
      enddo
    
      ! CIA
      do i = 1,op%ncia
        call interpolate_Xsection(op%cia(i), l, v%P, v%T, cia(:,i))
      enddo
      
      ! Absorption xs
      do i = 1,op%naxs
        call interpolate_Xsection(op%axs(i), l, v%P, v%T, axs(:,i))
      enddo
      
      ! Photolysis xs
      do i = 1,op%npxs
        call interpolate_Xsection(op%pxs(i), l, v%P, v%T, pxs(:,i))
      enddo
    
      ! compute tau
      ! rayleigh scattering
      tausg(:) = 0.0_dp
      do i = 1,op%nray
        j = op%ray(i)%sp_ind(1)
        do k = 1,v%nz
          n = v%nz+1-k
          tausg(n) = tausg(n) + op%ray(i)%xs_0d(l)*v%densities(k,j)*v%dz(k)
        enddo
      enddo
      
      ! CIA
      taua(:) = 0.0_dp
      do i = 1,op%ncia
        j = op%cia(i)%sp_ind(1)
        jj = op%cia(i)%sp_ind(2)
        do k = 1,v%nz
          n = v%nz+1-k
          taua(n) = taua(n) + cia(k,i)*v%densities(k,j)*v%densities(k,jj)*v%dz(k)
        enddo
      enddo
      
      ! absorption
      do i = 1,op%naxs
        j = op%cia(i)%sp_ind(1)
        do k = 1,v%nz
          n = v%nz+1-k
          taua(n) = taua(n) + axs(k,i)*v%densities(k,j)*v%dz(k)
        enddo
      enddo
      
      ! photolysis
      do i = 1,op%npxs
        j = op%cia(i)%sp_ind(1)
        do k = 1,v%nz
          n = v%nz+1-k
          taua(n) = taua(n) + pxs(k,i)*v%densities(k,j)*v%dz(k)
        enddo
      enddo
      
      ! plank function, only if in the IR
      ! bplanck has units [W sr^−1 m^−2 Hz^-1]
      if (op%op_type == IROpticalProperties) then
        bplanck(v%nz+1) = planck_fcn(op%freq(l), v%T(1)) ! ground level
        do j = 1,v%nz
          n = v%nz+1-j
          bplanck(n) = planck_fcn(op%freq(l), v%T(j))
        enddo
      endif
      
      gt = 0.0_dp
        
      ! two species with k coefficients  
      if (op%nk == 2) then
        
        j1 =op%k(1)%sp_ind
        j2 =op%k(2)%sp_ind
        
        do i1 = 1,ks(1)%ngauss
          do i2 = 1,ks(2)%ngauss
            
            taua_1(:) = 0.0_dp
            do k = 1,v%nz
              n = v%nz+1-k
              taua_1(n) = taua_1(n) + &
                          ks(1)%k(k,i1)*v%densities(k,j1)*v%dz(k) + &
                          ks(2)%k(k,i2)*v%densities(k,j2)*v%dz(k)
            enddo
            
            ! sum
            tau = tausg + taua + taua_1
            do i = 1,v%nz
              w0(i) = min(max_w0,tausg(i)/tau(i))
            enddo
            
            if (op%op_type == FarUVOpticalProperties .or. &
                op%op_type == SolarOpticalProperties) then
              call two_stream_solar(v%nz, tau, w0, gt, u0, surface_albedo, &
                                    amean, surf_rad, fup, fdn)
            elseif (op%op_type == IROpticalProperties) then
              call two_stream_ir(v%nz, tau, w0, gt, surface_albedo, bplanck, &
                                 fup, fdn)
            endif
            
            gauss_weight = op%k(1)%weights(i1)*op%k(2)%weights(i2)
            fup1 = fup1 + fup*gauss_weight
            fdn1 = fdn1 + fdn*gauss_weight
          
          enddo
        enddo
        
      ! three species with k coefficients
      elseif (op%nk == 3) then
      
        j1 =op%k(1)%sp_ind
        j2 =op%k(2)%sp_ind
        j3 =op%k(3)%sp_ind
      
        do i1 = 1,ks(1)%ngauss
          do i2 = 1,ks(2)%ngauss
            do i3 = 1,ks(3)%ngauss
      
              taua_1(:) = 0.0_dp
              do k = 1,v%nz
                n = v%nz+1-k
                taua_1(n) = taua_1(n) + &
                            ks(1)%k(k,i1)*v%densities(k,j1)*v%dz(k) + &
                            ks(2)%k(k,i2)*v%densities(k,j2)*v%dz(k) + &
                            ks(2)%k(k,i3)*v%densities(k,j3)*v%dz(k)
              enddo
      
              ! sum
              tau = tausg + taua + taua_1
              do i = 1,v%nz
                w0(i) = min(max_w0,tausg(i)/tau(i))
              enddo
      
              if (op%op_type == FarUVOpticalProperties .or. &
                  op%op_type == SolarOpticalProperties) then
                call two_stream_solar(v%nz, tau, w0, gt, u0, surface_albedo, &
                                      amean, surf_rad, fup, fdn)
              elseif (op%op_type == IROpticalProperties) then
                call two_stream_ir(v%nz, tau, w0, gt, surface_albedo, bplanck, &
                                   fup, fdn)
              endif

              gauss_weight = op%k(1)%weights(i1)*op%k(2)%weights(i2)*op%k(3)%weights(i3)
              fup1 = fup1 + fup*gauss_weight
              fdn1 = fdn1 + fdn*gauss_weight
              
            enddo
          enddo
        enddo
        
      else
        print*,'ERR'
        stop 1 
      endif
      
      fup_a(:,l) =  max(fup1,0.0_dp)
      fdn_a(:,l) =  max(fdn1,0.0_dp)
      
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
  
  ! function planck_fcn(lambda_nm, T) result(B)
  !   use clima_const, only: c_light, k_boltz_si, plank
  !   real(dp), intent(in) :: lambda_nm ! (nm) 
  !   real(dp), intent(in) :: T ! (K)
  !   real(dp) :: B ! W sr^−1 m^−2 m^-1
  ! 
  !   real(dp) :: lambda
  ! 
  !   lambda = lambda_nm*1.0e-9_dp ! convert to meters
  ! 
  !   B = ((2.0_dp*plank*c_light**2.0_dp)/(lambda**5.0_dp)) * &
  !       ((1.0_dp)/(exp((plank*c_light)/(lambda*k_boltz_si*T)) - 1.0_dp)) 
  ! end function
  
  function planck_fcn(nu, T) result(B)
    use clima_const, only: c_light, k_boltz => k_boltz_si, plank
    real(dp), intent(in) :: nu ! (1/s) 
    real(dp), intent(in) :: T ! (K)
    real(dp) :: B ! W sr^−1 m^−2 Hz^-1
    
    B = ((2.0_dp*plank*nu**3.0_dp)/(c_light**2.0_dp)) * &
        ((1.0_dp)/(exp((plank*nu)/(k_boltz*T)) - 1.0_dp)) 
    
  end function
  
  
end module