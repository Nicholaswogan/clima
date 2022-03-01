module clima_const
  use iso_fortran_env, only: dp => real64
  implicit none
  real(dp), parameter :: c_light = 299792458.0_dp
  
  
  
  ! Wavenumbers at the ends of the spectral intervals (1/cm)
  integer, parameter :: n_ir = 55
  real(dp), parameter :: ir_wavenums(n_ir+1) =  &
      [0., 40., 100., 160., 220., 280., 330., 380., 440., 495., &
       545., 617., 667., 720., 800., 875., 940., 1000., 1065., &
       1108., 1200., 1275., 1350., 1450., 1550., 1650., 1750., 1850., &
       1950., 2050., 2200., 2397., 2494., 2796., 3087., 3425., 3760., &
       4030., 4540., 4950., 5370., 5925., 6390., 6990., 7650., 8315., &
       8850., 9350., 9650., 10400., 11220., 11870., 12790., 13300., &
       14470., 15000.]
       
  integer, parameter :: n_sol = 38
  ! real(dp), parameter :: sol_wavenums(n_sol+1) =  &
      ! []
  
end module