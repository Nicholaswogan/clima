program test_adiabat
  use clima_adiabat_kasting, only: make_profile, dp
  implicit none
  real(dp) :: T_surf, P_H2O_surf, P_CO2_surf, P_N2_surf
  integer :: nz
  real(dp) :: P_top, T_strat
  real(dp), allocatable :: P(:), T(:), f_H2O(:), f_CO2(:), f_N2(:)
  logical :: surface_liquid_H2O
  character(:), allocatable :: err
  
  integer :: i
  
  T_surf = 1000.0_dp
  P_H2O_surf = 50_dp*1e6_dp
  P_CO2_surf = 1e6_dp
  P_N2_surf = 1e6_dp
  nz = 100
  P_top = 1.0e0_dp
  T_strat = 200.0_dp
  allocate(P(nz),T(nz),f_H2O(nz),f_CO2(nz),f_N2(nz))
  
  call make_profile(T_surf, P_H2O_surf, P_CO2_surf, P_N2_surf, nz, P_top, T_strat, &
                    P, T, f_H2O, f_CO2, f_N2, surface_liquid_H2O, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  do i =1,nz
    print*,P(i),T(i),f_H2O(i),f_CO2(i),f_N2(i)
  enddo
  print*,surface_liquid_H2O
  
end program