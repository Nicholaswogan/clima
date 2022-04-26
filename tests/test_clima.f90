
program test_clima
  use clima, only: dp, Climate
  implicit none

  type(Climate) :: c
  real(dp), allocatable :: dTdt(:)
  character(:), allocatable :: err
  
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
  
  call c%right_hand_side(c%v%T_init,dTdt)
  

end program