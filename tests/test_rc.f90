program main
  use clima, only: dp, s_str_len
  use clima_rc, only: RadiativeConvectiveClimate
  call test()
contains
  subroutine test()
    type(RadiativeConvectiveClimate) :: c
    character(:), allocatable :: err
    character(s_str_len), allocatable :: condensible_names(:)
    real(dp), allocatable :: condensible_P(:)
    real(dp), allocatable :: f_i(:,:), T_init(:)

    c = RadiativeConvectiveClimate('../clima/data', &
                     '../templates/ModernEarth/species.yaml', &
                     '../templates/ModernEarth/settings.yaml', &
                     '../templates/ModernEarth/Sun_now.txt', &
                     err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    allocate(f_i(c%nz,c%sp%ng), T_init(c%nz))

    condensible_names = ['H2O', 'CO2']
    condensible_P = [270.0e6_dp, 400.0_dp]
    f_i(:,1) = 1.0e-2_dp
    f_i(:,2) = 1.0e-2_dp
    f_i(:,3) = 9.0e-1_dp
    f_i(:,4) = 1.0e-10_dp
    f_i(:,5) = 1.0e-10_dp
    f_i(:,6) = 1.0e-10_dp
    call c%initialize_stepper(condensible_names, condensible_P, f_i, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

  end subroutine
end program