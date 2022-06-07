program test_adiabat
  use clima_const, only: dp
  use clima_adiabat, only: AdiabatClimateModel
  implicit none
  
  type(AdiabatClimateModel) :: c
  character(:), allocatable :: err
  real(dp) :: OLR
  integer :: i
  
  c = AdiabatClimateModel('../data', &
                          '../templates/runaway_greenhouse/species.yaml', &
                          '../templates/runaway_greenhouse/settings.yaml', &
                          err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  OLR = c%OLR(325.0_dp, [270.0e6_dp, 0.0e0_dp, 1.0e6_dp], err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  print*,OLR
  
  ! do i = 1,c%nz
  !   print*,c%z(i)/1.0e5_dp,c%dz(i)/1.0e5_dp,c%P(i)/1.0e6_dp, c%f_i(i,2), c%densities(i,2)
  ! enddo
  ! 
  ! print*,any(c%densities /= c%densities)
  ! print*,any(c%T /= c%T)
  ! print*,any(c%P /= c%P)
    
end program