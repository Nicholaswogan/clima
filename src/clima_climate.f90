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
    c%v%density = c%v%P/(k_boltz*c%v%T)
    do i = 1,c%d%ng
      c%v%densities(:,i) = c%v%mix(:,i)*c%v%density
    enddo
    
  end function
  
  subroutine radiative_transfer(d, v)
    type(ClimaData), intent(inout) :: d
    type(ClimaVars), intent(in) :: v
    
    ! Solar radiative transfer
    call radiate(d%ir, v)
    
    
    
  end subroutine
  
  subroutine radiate(op, v)
    use clima_const, only: pi
    use clima_types, only: OpticalProperties, Kcoefficients
    use clima_twostream, only: two_stream
    use futils, only: Timer
    type(OpticalProperties), intent(inout) :: op
    type(ClimaVars), intent(in) :: v
    
    real(dp), parameter :: max_w0 = 0.99999e0_dp
    
    integer :: i, j, k, l, n, jj
    integer :: i1, i2, i3
    integer :: j1, j2, j3
    integer :: ie
    real(dp) :: u0
    real(dp) :: surface_albedo, surf_rad
    type(Timer) :: tm
    
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
    
    
    ! solar zenith angle
    u0 = cos(50.0_dp*pi/180.0_dp)
    surface_albedo = 0.25
    
    call tm%start()
    

    !$omp parallel private(i, j, k, l, n, jj, &
    !$omp& i1, i2, i3, j1, j2, j3, ie, surf_rad, val, &
    !$omp& ks, cia, axs, pxs, tausg, taua, taua_1, &
    !$omp& tau, w0, gt, amean)
    
    !$omp do
    do l = 1,op%nw
    
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
      
      ! gt
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
            
            call two_stream(v%nz, tau, w0, gt, u0, surface_albedo, amean, surf_rad, ie)
            
            
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
      
              call two_stream(v%nz, tau, w0, gt, u0, surface_albedo, amean, surf_rad, ie)

            enddo
          enddo
        enddo
        
      
      else
        print*,'ERR'
        stop 1 
      endif
      
      
    enddo
    !$omp enddo
    !$omp end parallel
    
    call tm%finish('')

    
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
  
  pure function planck_fcn(v, T) result(B)
    use clima_const, only: c_light, k_boltz_si, plank
    real(dp), intent(in) :: v, T
    real(dp) :: B
    B = ((2*plank*v**3.0_dp)/(c_light**2.0_dp)) * &
        ((1.0_dp)/(exp((plank*v)/(k_boltz_si*T)) - 1)) 
  end function
  
  
end module