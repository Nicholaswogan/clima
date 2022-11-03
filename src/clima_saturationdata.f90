module clima_saturationdata
  use clima_const, only: dp, s_str_len
  implicit none
  private
  
  public :: SaturationData

  !! This module specifies the SaturationData class, which is used
  !! to compute latents heats of sublimiation/vaporization, as well as 
  !! saturation vapor pressures above solids or liquids

  type :: SaturationData
    real(dp) :: mu !! molar mass (g/mol)
    real(dp) :: T_ref !! Reference Temperature (K)
    real(dp) :: P_ref !! Reference pressure (dynes/cm^2)
    real(dp) :: T_triple !! Triple point temperature
    real(dp) :: T_critical !! Critical temperature (K)

    !> Latent heat fit parameters for vaporization (no units)
    real(dp) :: A_v, B_v, C_v
    !> Latent heat fit parameters for sublimation (no units)
    real(dp) :: A_s, B_s, C_s

  contains
    procedure :: latent_heat_vap
    procedure :: latent_heat_sub
    procedure :: latent_heat
    procedure :: sat_pressure_vap
    procedure :: sat_pressure_sub
    procedure :: sat_pressure
  end type
  interface SaturationData
    procedure create_SaturationData
  end interface

contains

  !> Latent heat of vaporization
  function latent_heat_vap(self, T) result(L)
    class(SaturationData), intent(inout) :: self
    real(dp), intent(in) :: T !! K
    real(dp) :: L !! erg/g
    L = self%A_v*exp(self%B_v*T)+self%C_v
  end function

  !> Latent heat of sublimation
  function latent_heat_sub(self, T) result(L)
    class(SaturationData), intent(inout) :: self
    real(dp), intent(in) :: T !! K
    real(dp) :: L !! erg/g
    L = self%A_s*exp(self%B_s*T)+self%C_s
  end function

  !> Latent heat of vaporization or sublimation
  function latent_heat(self, T) result(L)
    class(SaturationData), intent(inout) :: self
    real(dp), intent(in) :: T !! K
    real(dp) :: L !! erg/g
    if (T > self%T_triple) then
      L = self%latent_heat_vap(T)
    else ! (T <= T_triple) then
      L = self%latent_heat_sub(T)
    endif
  end function

  !> Saturation pressure over liquid
  function sat_pressure_vap(self, T) result(p_sat)
    use clima_const, only: Rgas
    class(SaturationData), intent(inout) :: self
    real(dp), intent(in) :: T !! K
    real(dp) :: p_sat !! dynes/cm2
    real(dp) :: tmp
    tmp = integral_fcn(self%A_v, self%B_v, self%C_v, T) - integral_fcn(self%A_v, self%B_v, self%C_v, self%T_ref)
    p_sat = self%P_ref*exp((self%mu/Rgas)*(tmp))
  end function

  !> Saturation pressure over solid
  function sat_pressure_sub(self, T) result(p_sat)
    use clima_const, only: Rgas
    class(SaturationData), intent(inout) :: self
    real(dp), intent(in) :: T !! K
    real(dp) :: p_sat !! dynes/cm2
    real(dp) :: tmp
    tmp = (integral_fcn(self%A_v, self%B_v, self%C_v, self%T_triple) - integral_fcn(self%A_v, self%B_v, self%C_v, self%T_ref)) + &
          (integral_fcn(self%A_s, self%B_s, self%C_s, T) - integral_fcn(self%A_s, self%B_s, self%C_s, self%T_triple))
    p_sat = self%P_ref*exp((self%mu/Rgas)*(tmp))
  end function

  !> Saturation pressure over liquid or solid
  function sat_pressure(self, T) result(p_sat)
    class(SaturationData), intent(inout) :: self
    real(dp), intent(in) :: T !! K
    real(dp) :: p_sat !! dynes/cm2
    if (T > self%T_triple) then
      p_sat = self%sat_pressure_vap(T)
    else ! (T <= T_triple) then
      p_sat = self%sat_pressure_sub(T)
    endif
  end function

  !> This is $\int L/T^2 dT$
  function integral_fcn(A, B, C, T) result(out)
    use futils, only: expi
    real(dp), intent(in) :: A, B, C, T !! K
    real(dp) :: out
    out = -(-A*B*T*expi(B*T) + A*exp(B*T) + C)/T
  end function

  function create_SaturationData(s, name, filename, err) result(sat)
    use fortran_yaml_c_types, only: type_dictionary, type_error
    type(type_dictionary), intent(in) :: s
    character(*), intent(in) :: name
    character(*), intent(in) :: filename
    character(:), allocatable, intent(out) :: err
    type(SaturationData) :: sat

    type(type_error), allocatable :: io_err
    class(type_dictionary), pointer :: tmpdict
    character(len=:), allocatable :: model
    
    model = s%get_string("model",error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (model /= "Wogan") then
      err = 'Saturation "model" must be "Wogan" for species "'//name//'" in '//filename
      return
    endif

    sat%mu = s%get_real('mu',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%T_ref = s%get_real('T-ref',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%P_ref = s%get_real('P-ref',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%T_triple = s%get_real('T-triple',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%T_critical = s%get_real('T-critical',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%T_triple = s%get_real('T-triple',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    tmpdict => s%get_dictionary("vaporization",.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%A_v = tmpdict%get_real('A',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%B_v = tmpdict%get_real('B',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%C_v = tmpdict%get_real('C',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    tmpdict => s%get_dictionary("sublimation",.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%A_s = tmpdict%get_real('A',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%B_s = tmpdict%get_real('B',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%C_s = tmpdict%get_real('C',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

  end function

end module