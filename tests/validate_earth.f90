program validate
  use clima, only: dp, Climate
  use clima_types, only: ClimaData, ClimaVars, ClimaWrk
  use clima_radiate, only: radiate
  use clima_eqns, only: press_and_den
  implicit none

  type(Climate), target :: c
  type(ClimaData), pointer :: d
  type(ClimaVars), pointer :: v
  type(ClimaWrk), pointer :: w
  character(:), allocatable :: err
  
  integer :: i
  
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
  
  d => c%d
  v => c%v
  w => c%w
  
  w%rin%T = v%T_init
  
  ! First, we compute P profile
  call press_and_den(v%nz, w%rin%T, v%grav, &
                     v%surface_pressure*1.e6_dp, v%dz, &
                     v%mubar, w%rin%P, w%rin%density)
  ! convert P from cgs to bar
  w%rin%P = w%rin%P/1.0e6_dp
  
  ! Compute densities and columns
  do i = 1,d%ng
    w%rin%densities(:,i) = v%mix(:,i)*w%rin%density(:)
    w%rin%cols(:,i) = w%rin%densities(:,i)*v%dz(:)
  enddo
  
  call radiate(d%sol, d%kset, &
               v%surface_albedo, v%u0, v%diurnal_fac, v%photons_sol, &
               w%rin%P, w%rin%T, w%rin%densities, v%dz, &
               w%rx_sol, w%rz, &
               w%fup_a_sol, w%fdn_a_sol, w%fup_sol, w%fdn_sol)
    
  ! IR radiative transfer
  call radiate(d%ir, d%kset, &
               v%surface_albedo, v%u0, v%diurnal_fac, v%photons_sol, &
               w%rin%P, w%rin%T, w%rin%densities, v%dz, &
               w%rx_ir, w%rz, &
               w%fup_a_ir, w%fdn_a_ir, w%fup_ir, w%fdn_ir)
  
  w%f_total = (-w%fup_sol + w%fdn_sol) + (-w%fup_ir + w%fdn_ir)
  
  open(1,file='../fup_earth.dat',status='replace',form='unformatted')
  ! write(1) w%rin%T
  ! write(1) w%rin%P
  ! write(1) [v%z-(v%dz/2.0_dp), v%z(v%nz)+v%dz(v%nz)/2.0_dp]
  ! write(1) w%f_total
  ! write(1) w%fdn_sol
  ! write(1) w%fup_sol
  ! write(1) w%fdn_ir
  ! write(1) w%fup_ir
  write(1) d%ir%freq
  write(1) w%fup_a_ir(v%nz+1,:)
  close(1)


end program