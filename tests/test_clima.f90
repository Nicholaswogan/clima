
program test_clima
  use clima
  use clima_const
  implicit none
  type(ktable) :: kt
  character(:), allocatable :: err
  
  kt = create_ktable('../data/kdistributions/HITRAN2016_ir_CO2.dat', ir_wavenums, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  kt = create_ktable('../data/kdistributions/HITRAN2016_ir_H2O.dat', ir_wavenums, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  kt = create_ktable('../data/kdistributions/HITRAN2016_solar_CO2.dat', sol_wavenums, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  kt = create_ktable('../data/kdistributions/HITRAN2016_solar_H2O.dat', sol_wavenums, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

end program