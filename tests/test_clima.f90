
program test_clima
  use clima
  implicit none
  type(Climate) :: c
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
  
  call radiative_transfer(c%d, c%v)

end program