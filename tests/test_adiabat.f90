program test
  use clima, only: dp
  use clima, only: AdiabatClimate
  use clima, only: ocean_solubility_fcn
  implicit none
  
  type(AdiabatClimate) :: c
  character(:), allocatable :: err
  real(dp) :: T, OLR
  integer :: i
  procedure(ocean_solubility_fcn), pointer :: ocean_fcn_ptr
  
  c = AdiabatClimate('../templates/runaway_greenhouse/species.yaml', &
                     '../templates/runaway_greenhouse/settings.yaml', &
                     '../templates/ModernEarth/Sun_now.txt', &
                     '../clima/data', &
                     err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  c%solve_for_T_trop = .true.
    
  T = c%surface_temperature( &
      [270.0e6_dp, 400e-6_dp*1.0e6_dp, 1.0e6_dp], &
      T_guess = 280.0_dp, err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  print*,T, c%T_trop

  T = c%surface_temperature_column( &
      [15.0e3_dp, 400e-6_dp*23.0_dp, 1.0*36.0_dp], &
      T_guess = 280.0_dp, err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  print*,T, c%T_trop

  T = c%surface_temperature_bg_gas( &
      [270.0e6_dp, 400e-6_dp*1.0e6_dp, 1.0e6_dp], &
      P_surf = 1.0e6_dp, bg_gas='N2', T_guess = 280.0_dp, err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  print*,T, c%T_trop

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

  ocean_fcn_ptr => ocean_fcn
  call c%set_ocean_solubility_fcn('H2O', ocean_fcn_ptr, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  call c%make_column(280.0_dp, &
      [15.0e4_dp, 400e-6_dp*23.0_dp, 1.0*36.0_dp], &
      err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  print*,c%N_surface+c%N_atmos+c%N_ocean(:,1)
  print*,15.0e3_dp, 400e-6_dp*23.0_dp, 1.0*36.0_dp

contains

  subroutine ocean_fcn(T_surf, ng, P_i, m_i) 
    use iso_c_binding, only: c_double, c_int
    real(c_double), value, intent(in) :: T_surf !! K
    integer(c_int), value, intent(in) :: ng
    real(c_double), intent(in) :: P_i(ng) !! bars
    real(c_double), intent(out) :: m_i(ng) !! mol/kg
    m_i(1) = 0.0_dp
    m_i(2) = P_i(2)*3.4e-2_dp
    m_i(3) = P_i(3)*6.1e-4_dp
  end subroutine

end program