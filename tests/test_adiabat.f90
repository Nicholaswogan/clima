program test_adiabat
  use clima_const, only: dp
  use clima_adiabat, only: AdiabatClimateModel
  implicit none
  
  type(AdiabatClimateModel) :: c
  character(:), allocatable :: err
  
  c = AdiabatClimateModel('../data', &
                          '../templates/runaway_greenhouse/species.yaml', &
                          '../templates/runaway_greenhouse/settings.yaml', &
                          err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  call c%make_profile(288.0_dp, [270.0e6_dp, 0.0e6_dp, 1.0e6_dp], err)
    
end program