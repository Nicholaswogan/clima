
program test_clima
  use clima
  use clima_const
  implicit none
  type(ktable) :: kt
  type(CIAtable) :: cia
  character(:), allocatable :: err
  
  kt = create_ktable('../data/kdistributions/CO2_ir.dat', 0, ir_wavenums, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  kt = create_ktable('../data/kdistributions/H2O_ir.dat', 0, ir_wavenums, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  kt = create_ktable('../data/kdistributions/CO2_solar.dat', 0, sol_wavenums, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif

  kt = create_ktable('../data/kdistributions/H2O_solar.dat', 0, sol_wavenums, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  cia = create_CIAtable('../data/CIA/CO2_CO2_CIA_ir.dat', [0,0], ir_wavenums, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  cia = create_CIAtable('../data/CIA/H2_H2_CIA_ir.dat', [0,0], ir_wavenums, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  cia = create_CIAtable('../data/CIA/H2_N2_CIA_ir.dat', [0,0], ir_wavenums, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  cia = create_CIAtable('../data/CIA/O2_O2_CIA_ir.dat', [0,0], ir_wavenums, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
end program