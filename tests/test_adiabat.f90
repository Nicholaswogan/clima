program test_adiabat
  use clima_const, only: dp
  use clima_adiabat, only: AdiabatClimateModel
  use stdlib_math, only: linspace
  implicit none
  
  type(AdiabatClimateModel) :: c
  character(:), allocatable :: err
  real(dp) :: T
  integer :: i
  
  c = AdiabatClimateModel('../data', &
                          '../templates/runaway_greenhouse/species.yaml', &
                          '../templates/runaway_greenhouse/settings.yaml', &
                          err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  c%T_trop = 180.0_dp
    
  T = c%surface_temperature(240.0_dp, &
      [270.0e6_dp, 400e-6_dp*1.0e6_dp, 1.0e6_dp], &
      T_guess = 300.0_dp, err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  print*,c%P(1)/1.0e6,T - 273.15_dp, c%T(size(c%T))
  print*,c%f_i(1,:)
  
end program