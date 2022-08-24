
program test
  implicit none

  call test_RadtranIR()
  call test_Radtran()

contains

  subroutine test_RadtranIR()
    use clima, only: RadtranIR
    use clima_types, only: AtmosphereFile, unpack_atmospherefile
    use clima_eqns, only: vertical_grid
    use clima_const, only: k_boltz, dp
    implicit none
  
    type(RadtranIR) :: rad
    type(AtmosphereFile) :: atm
    character(:), allocatable :: err
    
    integer :: nz, i
    real(dp) :: OLR
    real(dp), allocatable :: T(:), P(:), densities(:,:), mix(:,:), dz(:), z(:), density(:)
    real(dp), allocatable :: pdensities(:,:), radii(:,:)
    
    atm = AtmosphereFile("../templates/ModernEarth/atmosphere.txt", err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
  
    nz = size(atm%columns(1,:))
    
    rad = RadtranIR("../clima/data","../templates/ModernEarth/settings.yaml",nz, err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    
    ! step up the atmosphere
    allocate(T(nz), P(nz), densities(nz,rad%ng), mix(nz,rad%ng), dz(nz), z(nz), density(nz))
    allocate(pdensities(nz,rad%np), radii(nz,rad%np))

    call vertical_grid(0.0_dp, 1.0e7_dp, nz, z, dz)
    call unpack_atmospherefile(atm, rad%species_names, z, mix, T, P, err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    
    density = (P*1.0e6_dp)/(k_boltz*T)
    do i = 1,rad%ng
      densities(:,i) = mix(:,i)*density
    enddo

    pdensities(:,:) = 1.0_dp
    radii(:,:) = 1.0e-5_dp
    
    OLR = rad%OLR(289.0_dp, T, P, densities, dz, pdensities, radii, err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    
    print*,OLR*1.0e-3_dp
    
    open(unit=1,file='ModernEarthIR.dat',form='unformatted',status='replace')
    write(1) rad%ir%freq
    write(1) rad%wrk_ir%fup_a(nz+1,:)
    close(1)

  end subroutine

  subroutine test_Radtran()
    use clima, only: Radtran
    use clima_types, only: AtmosphereFile, unpack_atmospherefile
    use clima_eqns, only: vertical_grid
    use clima_const, only: k_boltz, dp
    implicit none
  
    type(Radtran) :: rad
    type(AtmosphereFile) :: atm
    character(:), allocatable :: err
    
    integer :: nz, i
    real(dp) :: solar_zenith, surface_albedo, T_shift
    real(dp), allocatable :: T(:), P(:), densities(:,:), mix(:,:), dz(:), z(:), density(:)
    real(dp), allocatable :: pdensities(:,:), radii(:,:)
    
    atm = AtmosphereFile("../templates/ModernEarth/atmosphere.txt", err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
  
    nz = size(atm%columns(1,:))

    solar_zenith = 60.0_dp
    surface_albedo = 0.15_dp
    rad = Radtran("../clima/data","../templates/ModernEarth/settings.yaml",&
    "../templates/ModernEarth/Sun_now.txt", solar_zenith, surface_albedo, nz, err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    ! rad%diurnal_fac = 1.0_dp
    
    ! step up the atmosphere
    allocate(T(nz), P(nz), densities(nz,rad%ng), mix(nz,rad%ng), dz(nz), z(nz), density(nz))
    allocate(pdensities(nz,rad%np), radii(nz,rad%np))
    
    call vertical_grid(0.0_dp, 1.0e7_dp, nz, z, dz)
    call unpack_atmospherefile(atm, rad%species_names, z, mix, T, P, err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    T_shift = 0.0_dp
    T = T + T_shift
    
    density = (P*1.0e6_dp)/(k_boltz*T)
    do i = 1,rad%ng
      densities(:,i) = mix(:,i)*density
    enddo

    pdensities(:,:) = 1.0_dp
    radii(:,:) = 1.0e-5_dp
    
    call rad%radiate(T(1), T, P, densities, dz, pdensities, radii, err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    print*,rad%wrk_sol%fdn_n(nz+1)*1.0e-3_dp

    open(unit=1,file='ModernEarth.dat',form='unformatted',status='replace')
    write(1) rad%ir%freq
    write(1) rad%wrk_ir%fup_a(nz+1,:)
    write(1) rad%sol%freq
    write(1) rad%wrk_sol%fup_a(nz+1,:)
    write(1) rad%photons_sol
    close(1)

  end subroutine
  
end program