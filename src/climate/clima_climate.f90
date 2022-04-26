module clima_climate
  use clima_const, only: dp
  use clima_types, only: ClimaData, ClimaVars, ClimaWrk
  implicit none
  private
  
  public :: Climate
  
  type :: Climate
    type(ClimaData) :: d
    type(ClimaVars) :: v
    type(ClimaWrk) :: w
    
  contains
    procedure :: right_hand_side
  end type
  
  interface Climate
    module procedure :: create_Climate
  end interface
  
  interface 
    
    module subroutine right_hand_side(self, T_in, dTdt)
      class(Climate), intent(inout), target :: self
      real(dp), intent(in) :: T_in(:)
      real(dp), intent(out) :: dTdt(:)
    end subroutine
  
  
  end interface
  
  
contains
  
  function create_Climate(data_dir, species_f, settings_f, star_f, atmosphere_f, err) result(c)
    use clima_const, only: k_boltz
    use clima_types, only: ClimaSettings
    use clima_input, only: create_ClimaVars, create_ClimaData, create_ClimaSettings
    character(*), intent(in) :: data_dir
    character(*), intent(in) :: species_f
    character(*), intent(in) :: settings_f
    character(*), intent(in) :: star_f
    character(*), intent(in) :: atmosphere_f
    character(:), allocatable, intent(out) :: err
    
    type(Climate) :: c
    
    type(ClimaSettings) :: s
    integer :: i
    
    s = create_ClimaSettings(settings_f, err)
    if (allocated(err)) return
    c%d = create_ClimaData(species_f, data_dir, s, err)
    if (allocated(err)) return
    c%v = create_ClimaVars(atmosphere_f, star_f, s, c%d, err)
    if (allocated(err)) return
    c%w = ClimaWrk(c%d, c%v%nz)
    
  end function
  
end module