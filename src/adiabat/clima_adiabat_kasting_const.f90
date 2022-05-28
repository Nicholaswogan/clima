module clima_adiabat_kasting_const
  use clima_const, only: dp
  implicit none
  
  !! molar masses (g/mol)
  real(dp), parameter :: mu_H2O = 18.0_dp
  real(dp), parameter :: mu_CO2 = 44.0_dp
  real(dp), parameter :: mu_N2 = 28.0_dp
  
  !! latent heat of H2O vaporization/condensation (erg/g)
  real(dp), parameter :: L_H2O = 2264e7_dp
  
end module
