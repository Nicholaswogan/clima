program test
  use clima_const, only: dp
  use clima_adiabat, only: WaterAdiabatClimateIR
  implicit none
  
  type(WaterAdiabatClimateIR) :: c
  character(:), allocatable :: err
  real(dp) :: T
  integer :: i
  
  c = WaterAdiabatClimateIR('../data', &
                          '../templates/runaway_greenhouse/species.yaml', &
                          '../templates/runaway_greenhouse/settings.yaml', &
                          err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  c%T_trop = 200.0_dp
    
  T = c%surface_temperature(260.0_dp, &
      [270.0e6_dp, 400e-6_dp*1.0e6_dp, 1.0e6_dp], &
      T_guess = 300.0_dp, err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  print*,T

  call c%make_column(280.0_dp, &
                    [15.0e3_dp, 0.0*23.0_dp, 1.0*36.0_dp], &
                    err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

end program