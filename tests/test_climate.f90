program test
  use clima, only: dp
  implicit none

  call test_Climate_evolve()

contains

  subroutine test_Climate_evolve()
    use clima, only: Climate
    use futils, only: linspace

    type(Climate) :: c
    character(:), allocatable :: err
    logical :: success
    real(dp), allocatable :: t_eval(:)

    c = Climate( &
    "../templates/ModernEarth/species.yaml", &
    "../templates/ModernEarth/settings.yaml", &
    "../templates/ModernEarth/Sun_now.txt", &
    "../templates/ModernEarth/atmosphere.txt", &
    "../clima/data", &
    err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    allocate(t_eval(100))
    call linspace(1.0_dp,10.0_dp,t_eval)
    t_eval = 10.0_dp**t_eval

    success = c%evolve("ModernEarthClimate.dat", 0.0_dp, c%T_init, t_eval, .true., err) 
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    
  end subroutine

end program