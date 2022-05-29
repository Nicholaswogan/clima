program test_adiabat
  use clima_const, only: dp
  use clima_adiabat_kasting, only: make_profile
  implicit none
  real(dp) :: T_surf, P_H2O_surf, P_CO2_surf, P_N2_surf
  integer :: nz
  real(dp) :: planet_mass, planet_radius
  real(dp) :: P_top, T_trop
  real(dp), allocatable :: P(:), T(:), f_H2O(:), f_CO2(:), f_N2(:), z(:)
  logical :: surface_liquid_H2O
  character(:), allocatable :: err
  
  integer :: i
  
  T_surf = 1000.0_dp
  P_H2O_surf = 200.0_dp*1e6_dp
  P_CO2_surf = 0.001_dp*1e6_dp
  P_N2_surf = 1.0_dp*1e6_dp
  nz = 100
  planet_mass = 5.972e27_dp
  planet_radius = 6.371e8_dp
  P_top = 1.0e0_dp
  T_trop = 200.0_dp
  allocate(P(nz),T(nz),f_H2O(nz),f_CO2(nz),f_N2(nz),z(nz))
  
  call make_profile(T_surf, P_H2O_surf, P_CO2_surf, P_N2_surf, nz, &
                    planet_mass, planet_radius, P_top, T_trop, &
                    P, z, T, f_H2O, f_CO2, f_N2, surface_liquid_H2O, &
                    err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  do i =1,nz
    print*,P(i),z(i),T(i),f_H2O(i),f_CO2(i),f_N2(i)
  enddo
  print*,surface_liquid_H2O
  
end program