program main
  use clima, only: dp
  use clima_rc, only: RadiativeConvectiveClimate
  call test()
contains
  subroutine test()
    type(RadiativeConvectiveClimate) :: c
    character(:), allocatable :: err

    c = RadiativeConvectiveClimate('../clima/data', &
                     '../templates/ModernEarth/species.yaml', &
                     '../templates/ModernEarth/settings.yaml', &
                     '../templates/ModernEarth/Sun_now.txt', &
                     err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
  end subroutine
end program