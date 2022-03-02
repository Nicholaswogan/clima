module clima_const
  use iso_fortran_env, only: dp => real64
  implicit none
  public
  
  integer, parameter :: s_str_len = 20
  
  ! physical constants
  real(dp), parameter :: c_light = 299792458.0_dp
  
  ! For now, we will hard-code wavelength bins. Doing this
  ! because the correlated-k data is specific to this binning.
  
  ! IR wavenumbers (1/cm)
  integer, parameter :: n_ir = 55
  real(dp), parameter :: ir_wavenums(n_ir+1) =  &
      [0.0_dp,                 40.0_dp,                100.0_dp,               &
       160.0_dp,               220.0_dp,               280.0_dp,               &
       330.0_dp,               380.0_dp,               440.0_dp,               &
       495.0_dp,               545.0_dp,               617.0_dp,               &
       667.0_dp,               720.0_dp,               800.0_dp,               &
       875.0_dp,               940.0_dp,               1000.0_dp,              &
       1065.0_dp,              1108.0_dp,              1200.0_dp,              &
       1275.0_dp,              1350.0_dp,              1450.0_dp,              &
       1550.0_dp,              1650.0_dp,              1750.0_dp,              &
       1850.0_dp,              1950.0_dp,              2050.0_dp,              &
       2200.0_dp,              2397.0_dp,              2494.0_dp,              &
       2796.0_dp,              3087.0_dp,              3425.0_dp,              &
       3760.0_dp,              4030.0_dp,              4540.0_dp,              &
       4950.0_dp,              5370.0_dp,              5925.0_dp,              &
       6390.0_dp,              6990.0_dp,              7650.0_dp,              &
       8315.0_dp,              8850.0_dp,              9350.0_dp,              &
       9650.0_dp,              10400.0_dp,             11220.0_dp,             &
       11870.0_dp,             12790.0_dp,             13300.0_dp,             &
       14470.0_dp,             15000.0_dp]
  ! IR wavelengths (micrometer)
  real(dp), parameter :: ir_wavelengths(n_ir+1) = &
      [huge(1.0_dp),            1.0e4_dp/ir_wavenums(2),  1.0e4_dp/ir_wavenums(3),   &
       1.0e4_dp/ir_wavenums(4),  1.0e4_dp/ir_wavenums(5),  1.0e4_dp/ir_wavenums(6),  &
       1.0e4_dp/ir_wavenums(7),  1.0e4_dp/ir_wavenums(8),  1.0e4_dp/ir_wavenums(9),  &
       1.0e4_dp/ir_wavenums(10), 1.0e4_dp/ir_wavenums(11), 1.0e4_dp/ir_wavenums(12), &
       1.0e4_dp/ir_wavenums(13), 1.0e4_dp/ir_wavenums(14), 1.0e4_dp/ir_wavenums(15), &
       1.0e4_dp/ir_wavenums(16), 1.0e4_dp/ir_wavenums(17), 1.0e4_dp/ir_wavenums(18), &
       1.0e4_dp/ir_wavenums(19), 1.0e4_dp/ir_wavenums(20), 1.0e4_dp/ir_wavenums(21), &
       1.0e4_dp/ir_wavenums(22), 1.0e4_dp/ir_wavenums(23), 1.0e4_dp/ir_wavenums(24), &
       1.0e4_dp/ir_wavenums(25), 1.0e4_dp/ir_wavenums(26), 1.0e4_dp/ir_wavenums(27), &
       1.0e4_dp/ir_wavenums(28), 1.0e4_dp/ir_wavenums(29), 1.0e4_dp/ir_wavenums(30), &
       1.0e4_dp/ir_wavenums(31), 1.0e4_dp/ir_wavenums(32), 1.0e4_dp/ir_wavenums(33), &
       1.0e4_dp/ir_wavenums(34), 1.0e4_dp/ir_wavenums(35), 1.0e4_dp/ir_wavenums(36), &
       1.0e4_dp/ir_wavenums(37), 1.0e4_dp/ir_wavenums(38), 1.0e4_dp/ir_wavenums(39), &
       1.0e4_dp/ir_wavenums(40), 1.0e4_dp/ir_wavenums(41), 1.0e4_dp/ir_wavenums(42), &
       1.0e4_dp/ir_wavenums(43), 1.0e4_dp/ir_wavenums(44), 1.0e4_dp/ir_wavenums(45), &
       1.0e4_dp/ir_wavenums(46), 1.0e4_dp/ir_wavenums(47), 1.0e4_dp/ir_wavenums(48), &
       1.0e4_dp/ir_wavenums(49), 1.0e4_dp/ir_wavenums(50), 1.0e4_dp/ir_wavenums(51), &
       1.0e4_dp/ir_wavenums(52), 1.0e4_dp/ir_wavenums(53), 1.0e4_dp/ir_wavenums(54), &
       1.0e4_dp/ir_wavenums(55), 1.0e4_dp/ir_wavenums(56)]
         
  ! Solar wavenumbers (1/cm)
  ! Note solar wavenumbers go the OPPOSITE direction than IR.
  integer, parameter :: n_sol = 38
  real(dp), parameter :: sol_wavenums(n_sol+1) =  &
      [42087.0_dp,             36363.0_dp,             35087.0_dp,             &
       32562.0_dp,             30376.0_dp,             29308.0_dp,             &
       25641.0_dp,             22222.0_dp,             18518.0_dp,             &
       18198.0_dp,             17649.0_dp,             16528.0_dp,             &
       16000.0_dp,             15000.0_dp,             14470.0_dp,             &
       13300.0_dp,             12790.0_dp,             11870.0_dp,             &
       11220.0_dp,             10400.0_dp,             9650.0_dp,              &
       9350.0_dp,              8850.0_dp,              8315.0_dp,              &
       7650.0_dp,              6990.0_dp,              6390.0_dp,              &
       5925.0_dp,              5370.0_dp,              4950.0_dp,              &
       4540.0_dp,              4030.0_dp,              3760.0_dp,              &
       3425.0_dp,              3087.0_dp,              2796.0_dp,              &
       2494.0_dp,              2397.0_dp,              2200.0_dp]
  ! Solar wavelengths (micrometer)
  real(dp), parameter :: sol_wavelengths(n_sol+1) = (1.0e4_dp/sol_wavenums)
  
end module