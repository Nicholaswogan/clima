program test_adiabat
  use clima_adiabat_kasting
  implicit none
  real(dp) :: T_surf, P_H2O_surf, P_CO2_surf, P_N2_surf
  integer :: nz
  real(dp) :: P_top, T_strat
  real(dp), allocatable :: P(:), T(:)
  character(:), allocatable :: err
  
  integer :: i
  
  T_surf = 1000.0_dp
  P_H2O_surf = 50_dp*1e6_dp
  P_CO2_surf = 1e6_dp
  P_N2_surf = 1e6_dp
  nz = 100
  P_top = 1.0e0_dp
  T_strat = 200.0_dp
  allocate(P(nz),T(nz))
  
  call make_profile(T_surf, P_H2O_surf, P_CO2_surf, P_N2_surf, nz, P_top, T_strat, P, T, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  do i =1,nz
    print*,P(i),T(i)
  enddo
  
end program