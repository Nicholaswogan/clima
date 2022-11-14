
module clima
  ! version
  use clima_version, only: version

  ! constants that matter
  use clima_const, only: dp, s_str_len 

  ! Radiative transfer classes
  use clima_radtran, only: Radtran, RadtranIR

  ! Adiabatic climate model
  use clima_adiabat, only: AdiabatClimate 

  ! Complicated climate model (experimental)
  use clima_climate, only: Climate 
end module