program validate
  use clima, only: dp, Climate
  use clima_radtran, only: radiate
  use clima_const, only: k_boltz, N_avo
  implicit none

  type(Climate) :: c
  character(:), allocatable :: err
  
  integer, parameter :: nz = 100
  integer :: i, j
  real(dp) :: ftotal(nz+1)
  real(dp) :: cgas(nz,8)
  real(dp) :: T(nz), T_1(nz), P_1(nz)
  real(dp) :: P(nz), density(nz), dz(nz), dz_1(nz)
  real(dp) :: densities(nz,7)
  
  c = Climate("../data", &
              "../species.yaml", &
              "../settings.yaml", &
              "../Sun_now.txt", &
              "../atmosphere.txt", &
              err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  open(1,file='../radtran.dat',status='old',form='unformatted')
  read(1) i
  read(1) ftotal
  read(1) ftotal
  read(1) ftotal
  read(1) ftotal
  read(1) ftotal
  read(1) cgas
  read(1) T
  read(1) P
  close(1)
  
  density = P/(k_boltz*T)
  dz = cgas(:,1)/density
  ! convert to bar
  P = P/1.0e6_dp
  
  
  do i = 1,nz
    j = nz-i+1
    densities(i,1) = cgas(j,6)/dz(j)
    densities(i,2) = cgas(j,5)/dz(j)
    densities(i,3) = (cgas(j,1)*0.78_dp)/dz(j)
    densities(i,4) = cgas(j,4)/dz(j)
    densities(i,5) = cgas(j,3)/dz(j)
    densities(i,6) = 1.0e-10_dp
    densities(i,7) = cgas(j,2)/dz(j)
    
    dz_1(i) = dz(j)
    P_1(i) = P(j)
    T_1(i) = T(j)
    
  enddo
  
  
  
  call radiate(c%d%sol, c%d%kset, &
               c%v%surface_albedo, c%v%u0, c%v%diurnal_fac, c%v%photons_sol, &
               P_1, T_1, densities, dz_1, &
               c%w%rx_sol, c%w%rz, &
               c%w%fup_a_sol, c%w%fdn_a_sol, c%w%fup_sol, c%w%fdn_sol)
               
  call radiate(c%d%ir, c%d%kset, &
               c%v%surface_albedo, c%v%u0, c%v%diurnal_fac, c%v%photons_sol, &
               P_1, T_1, densities, dz_1, &
               c%w%rx_ir, c%w%rz, &
               c%w%fup_a_ir, c%w%fdn_a_ir, c%w%fup_ir, c%w%fdn_ir)
  
  
  c%w%f_total = (-c%w%fup_sol + c%w%fdn_sol) + (-c%w%fup_ir + c%w%fdn_ir)
  
  open(1,file='../clima_radtran.dat',status='replace',form='unformatted')
  write(1) c%w%f_total
  write(1) c%w%fdn_sol
  write(1) c%w%fup_sol
  write(1) c%w%fdn_ir
  write(1) c%w%fup_ir
  close(1)
  
  
  ! print*,c%w%f_total
  
  

end program