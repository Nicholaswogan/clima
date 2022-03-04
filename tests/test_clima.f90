
program test_clima
  use clima
  implicit none
  type(ClimaSettings) :: s
  type(ClimaData) :: dat
  character(:), allocatable :: err
  
  s = create_ClimaSettings('../settings.yaml', err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  dat = create_ClimaData('../species.yaml', '../data', s, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

end program