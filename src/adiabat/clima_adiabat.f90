module clima_adiabat
  use clima_const, only: dp, s_str_len
  use clima_radtran, only: RadtranIR
  implicit none
  
  type :: KastingClimateModel
    ! species in the model
    character(s_str_len), allocatable :: species_names(:)
    type(RadtranIR) :: rad
  end type
  
  interface KastingClimateModel
    module procedure :: create_KastingClimateModel
  end interface
  
contains
  
  function create_KastingClimateModel(datadir, nz, err) result(c)
    use clima_types, only: ClimaSettings 
    character(*), intent(in) :: datadir
    integer, intent(in) :: nz
    character(:), allocatable, intent(out) :: err
    
    type(KastingClimateModel) :: c
    
    type(ClimaSettings) :: s
    
    ! create the settings
    allocate(s%ir)
    s%ir%k_method = "RandomOverlapResortRebin"
    s%ir%nbins = 16
    s%ir%k_distributions = ['H2O', 'CO2']  
    s%ir%cia = ['CO2-CO2', 'N2-N2  ', 'H2O-H2O', 'H2O-N2 ']
    s%ir%rayleigh = ['H2O', 'CO2', 'N2 ']
    
    ! species names
    c%species_names = ['H2O','CO2','N2 '] 
    ! make the radtranir object
    c%rad = RadtranIR(datadir, c%species_names, s, nz, err)
    if (allocated(err)) return
    
  end function
  


  
end module