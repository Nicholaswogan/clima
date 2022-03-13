module clima_const
  use iso_fortran_env, only: dp => real64
  implicit none
  public
  
  integer, parameter :: s_str_len = 20
  
  ! physical constants
  real(dp), parameter :: Rgas = 8.31446261815324e0_dp ! ideal gas constant (j/(mol*K))
  real(dp), parameter :: k_boltz = 1.380649e-16_dp ! boltzmann's constant cgs units (egs/K)
  real(dp), parameter :: k_boltz_si = 1.380649e-23_dp ! boltzmann's constant si units (J/K)
  real(dp), parameter :: G_grav = 6.67430e-11_dp ! gravitational constant (N * m2 / kg)
  real(dp), parameter :: plank = 6.62607004e-34_dp ! planks constant (m2 kg / s)
  real(dp), parameter :: c_light = 299792458.0_dp ! Speed of light (m / s)
  real(dp), parameter :: N_avo = 6.02214076e23_dp ! avagadros number
  real(dp), parameter :: pi = 3.14159265358979323846e0_dp
  
  ! useful
  real(dp), parameter :: log10tiny = log10(sqrt(tiny(1.0_dp)))
  
  ! For now, we will hard-code wavelength bins. Doing this
  ! because the correlated-k data is specific to this binning.
  
  ! IR wavenumbers (1/cm)
  integer, parameter :: n_ir = 55
  real(dp), parameter :: ir_wavenums(n_ir+1) =  &
      [15000.0_dp,             14470.0_dp,             13300.0_dp,             &
       12790.0_dp,             11870.0_dp,             11220.0_dp,             &
       10400.0_dp,             9650.0_dp,              9350.0_dp,              &
       8850.0_dp,              8315.0_dp,              7650.0_dp,              &
       6990.0_dp,              6390.0_dp,              5925.0_dp,              &
       5370.0_dp,              4950.0_dp,              4540.0_dp,              &
       4030.0_dp,              3760.0_dp,              3425.0_dp,              &
       3087.0_dp,              2796.0_dp,              2494.0_dp,              &
       2397.0_dp,              2200.0_dp,              2050.0_dp,              &
       1950.0_dp,              1850.0_dp,              1750.0_dp,              &
       1650.0_dp,              1550.0_dp,              1450.0_dp,              &
       1350.0_dp,              1275.0_dp,              1200.0_dp,              &
       1108.0_dp,              1065.0_dp,              1000.0_dp,              &
       940.0_dp,               875.0_dp,               800.0_dp,               &
       720.0_dp,               667.0_dp,               617.0_dp,               &
       545.0_dp,               495.0_dp,               440.0_dp,               &
       380.0_dp,               330.0_dp,               280.0_dp,               &
       220.0_dp,               160.0_dp,               100.0_dp,               &
       40.0_dp,                1.0_dp]
  ! IR wavelengths (nanometers)
  real(dp), parameter :: ir_wavl(n_ir+1) = (1.0e7_dp/ir_wavenums)
         
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
  ! Solar wavelengths (nanometers)
  real(dp), parameter :: sol_wavl(n_sol+1) = (1.0e7_dp/sol_wavenums)
  
end module