
program test_clima
  use clima, only: dp, Climate
  use stdlib_math, only: linspace
  implicit none

  type(Climate) :: c
  real(dp), allocatable :: dTdt(:), t_eval(:)
  character(:), allocatable :: err
  logical :: success
  
  c = Climate("../data", &
              "../species.yaml", &
              "../settings.yaml", &
              "../Sun_now.txt", &
              "../atmosphere.txt", &
              err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  allocate(dTdt(c%v%nz))
  
  call c%right_hand_side(c%v%T_init,dTdt, err)
  
  allocate(t_eval(100))
  
  t_eval = 10.0_dp**linspace(0.0_dp,9.0_dp,400)
  
  
  success = c%evolve("../test1.dat", 0.0_dp, c%v%T_init, t_eval, .true., err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  

end program