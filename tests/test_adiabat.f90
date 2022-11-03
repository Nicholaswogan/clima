program test
  use clima, only: dp
  use clima, only: AdiabatClimate
  implicit none
  
  type(AdiabatClimate) :: c
  character(:), allocatable :: err
  real(dp) :: T, OLR
  integer :: i
  
  c = AdiabatClimate('../clima/data', &
                     '../templates/runaway_greenhouse/species.yaml', &
                     '../templates/runaway_greenhouse/settings.yaml', &
                     '../templates/ModernEarth/Sun_now.txt', &
                     err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
    
  T = c%surface_temperature( &
      [270.0e6_dp, 400e-6_dp*1.0e6_dp, 1.0e6_dp], &
      T_guess = 280.0_dp, err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  print*,T

  T = c%surface_temperature_column( &
      [15.0e3_dp, 400e-6_dp*23.0_dp, 1.0*36.0_dp], &
      T_guess = 280.0_dp, err=err)
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

  call c%to_regular_grid(err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  ! A case where CO2 will also begin condensing
  c%T_trop = 150.0_dp
  call c%make_profile(280.0_dp, &
      [270.0e6_dp, 10.0e6_dp, 1.0e6_dp], &
      err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

end program