program test
  use clima, only: dp
  use clima, only: AdiabatClimate
  use clima, only: ocean_solubility_fcn
  implicit none
  
  type(AdiabatClimate) :: c
  character(:), allocatable :: err
  real(dp) :: T, OLR, ISR, OLR1, ISR1
  real(dp), allocatable :: P_i_surf(:), N_i_surf(:), f_i(:,:)
  integer :: i
  logical :: converged
  procedure(ocean_solubility_fcn), pointer :: ocean_fcn_ptr
  
  c = AdiabatClimate('../templates/AdiabatClimate/species.yaml', &
                     '../templates/AdiabatClimate/Earth/settings.yaml', &
                     '../templates/ModernEarth/Sun_now.txt', &
                     '../clima/data', &
                     err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  ! Test tidally locked dayside
  c%solve_for_T_trop = .true.
  c%tidally_locked_dayside = .true.
  P_i_surf = [270.0_dp, 400e-6_dp, 1.0_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp]
  P_i_surf = P_i_surf*1.0e6_dp
  T = c%surface_temperature( &
      P_i_surf, &
      T_guess = 280.0_dp, err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  print*,T, c%T_trop

  ! Test surface_temperature_column
  c%tidally_locked_dayside = .false.
  N_i_surf = [15.0e3_dp, 400e-6_dp*23.0_dp, 1.0*36.0_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp]
  T = c%surface_temperature_column( &
      N_i_surf, &
      T_guess = 280.0_dp, err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  print*,T, c%T_trop

  ! Test surface_temperature_bg_gas
  P_i_surf = [270.0_dp, 400e-6_dp, 1.0_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp]
  P_i_surf = P_i_surf*1.0e6_dp
  T = c%surface_temperature_bg_gas( &
      P_i_surf, &
      P_surf = 1.0e6_dp, bg_gas='N2', T_guess = 280.0_dp, err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  print*,T, c%T_trop

  ! Test to_regular_grid
  call c%to_regular_grid(err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  ! Test a case where CO2 begins to condense
  c%T_trop = 150.0_dp
  P_i_surf = [270.0_dp, 10.0_dp, 1.0_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp]
  P_i_surf = P_i_surf*1.0e6_dp
  call c%TOA_fluxes(280.0_dp, &
      P_i_surf, ISR, OLR, &
      err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  ! Test make_profile_dry
  allocate(f_i(c%nz+1,c%sp%ng))
  f_i(1,:) = c%f_i(1,:)
  do i = 1,c%nz
    f_i(i+1,:) = c%f_i(i,:)
  enddo
  call c%TOA_fluxes_dry( &
    [c%P_surf,c%P], &
    [c%T_surf,c%T], &
    f_i, &
    ISR1, OLR1, &
    err &
  )
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  print*,ISR1/ISR, OLR1/OLR

  ! Test ocean solubility functionality
  ocean_fcn_ptr => ocean_fcn
  call c%set_ocean_solubility_fcn('H2O', ocean_fcn_ptr, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  c%use_make_column_P_guess = .true.
  N_i_surf = [15.0e4_dp, 400e-6_dp*23.0_dp, 1.0*36.0_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp]
  call c%make_column(280.0_dp, &
      N_i_surf, &
      err=err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  print*,c%N_surface+c%N_atmos+c%N_ocean(:,1)
  print*,N_i_surf
  
  call test_RCE()

  ! deallocate
  deallocate(P_i_surf, N_i_surf)
contains

  subroutine ocean_fcn(T_surf, ng, P_i, m_i, args_p) 
    use iso_c_binding, only: c_double, c_int, c_ptr
    real(c_double), value, intent(in) :: T_surf !! K
    integer(c_int), value, intent(in) :: ng
    real(c_double), intent(in) :: P_i(ng) !! bars
    real(c_double), intent(out) :: m_i(ng) !! mol/kg
    type(c_ptr), value, intent(in) :: args_p
    m_i(:) = 0.0_dp
    m_i(2) = P_i(2)*3.4e-2_dp
    m_i(3) = P_i(3)*6.1e-4_dp
  end subroutine

  subroutine test_RCE()

    c = AdiabatClimate('../templates/AdiabatClimate/species.yaml', &
                     '../tests/settings_RCE_test.yaml', &
                     '../templates/ModernEarth/Sun_now.txt', &
                     '../clima/data', &
                     err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    c%P_top = 1.0e2_dp
    c%T_trop = 200.0_dp
    P_i_surf = [270.0_dp, 400e-6_dp, 1.0_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp, 1.0e-10_dp]
    P_i_surf = P_i_surf*1.0e6_dp
    T = c%surface_temperature(P_i_surf, T_guess=280.0_dp, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    c%max_rc_iters = 1
    converged = c%RCE( &
        P_i_surf, &
        c%T_surf, &
        c%T, &
        c%convecting_with_below, &
        err=err&
    )
    if (allocated(err)) then
      print*,err
      stop 1
    endif

  end subroutine

end program