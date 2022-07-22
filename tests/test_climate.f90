program test_climate
  use clima_const, only: dp
  implicit none

  call test_Climate_evolve()

contains

  subroutine test_Climate_evolve()
    use clima_climate, only: Climate
    use futils, only: linspace

    type(Climate) :: c
    character(:), allocatable :: err
    logical :: success
    real(dp), allocatable :: t_eval(:)


    c = Climate( &
    "../data", &
    "../species.yaml", &
    "../settings.yaml", &
    "../Sun_now.txt", &
    "../atmosphere.txt", &
    err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    allocate(t_eval(100))
    call linspace(1.0_dp,10.0_dp,t_eval)
    t_eval = 10.0_dp**t_eval

    ! c%T_init = 250.0_dp

    success = c%evolve("../test1.dat", 0.0_dp, c%T_init, t_eval, .true., err) 
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    
  end subroutine

end program