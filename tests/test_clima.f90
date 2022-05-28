
program test_clima
  use clima_radtran, only: RadtranIR, dp, RadtranIR_OLR
  use clima_types, only: AtmosphereFile, unpack_atmospherefile
  use clima_eqns, only: vertical_grid
  use clima_const, only: k_boltz
  implicit none

  type(RadtranIR) :: rad
  type(AtmosphereFile) :: atm
  character(:), allocatable :: err
  
  integer :: nz, i
  real(dp) :: OLR
  real(dp), allocatable :: T(:), P(:), densities(:,:), mix(:,:), dz(:), z(:), density(:)
  
  nz = 200
  
  rad = RadtranIR("../data","../settings.yaml",nz, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  ! step up the atmosphere
  allocate(T(nz), P(nz), densities(nz,rad%ng), mix(nz,rad%ng), dz(nz), z(nz), density(nz))
  
  atm = AtmosphereFile("../atmosphere.txt", err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
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
  
  OLR = RadtranIR_OLR(rad, T, P, densities, dz, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  print*,OLR,rad%wrk_ir%fup_n(1),rad%wrk_ir%fdn_n(1)
  
  open(unit=1,file='../OLR_1.dat',form='unformatted',status='replace')
  write(1) 0.5_dp*(rad%ir%freq(1:rad%ir%nw)+rad%ir%freq(2:rad%ir%nw+1))
  write(1) rad%wrk_ir%fup_a(nz+1,:)
  close(1)
  

end program