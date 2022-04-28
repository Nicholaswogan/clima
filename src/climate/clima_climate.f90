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
    procedure :: evolve
  end type
  
  interface Climate
    module procedure :: create_Climate
  end interface
  
  interface 
    
    module subroutine right_hand_side(self, T_in, dTdt, err)
      class(Climate), intent(inout), target :: self
      real(dp), intent(in) :: T_in(:)
      real(dp), intent(out) :: dTdt(:)
      character(:), allocatable :: err
    end subroutine
    
    module function evolve(self, filename, tstart, T_start, t_eval, overwrite, err) result(success)
      use, intrinsic :: iso_c_binding
      class(Climate), target, intent(inout) :: self
      character(len=*), intent(in) :: filename
      real(c_double), intent(in) :: tstart
      real(dp), intent(in) :: T_start(:)
      real(c_double), intent(in) :: t_eval(:)
      logical, optional, intent(in) :: overwrite
      logical :: success
      character(:), allocatable, intent(out) :: err
    end function
  
  
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