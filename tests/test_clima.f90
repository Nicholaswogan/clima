
program test_clima
  use clima_input, only: create_ClimaSettings
  use clima_radtran, only: ClimaRadtranIR
  implicit none

  type(ClimaRadtranIR) :: rad
  character(:), allocatable :: err
  
  rad = ClimaRadtranIR("../data","../settings.yaml",100, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  
end program