
module clima_types
  use clima_const, only: dp
  implicit none
  
  type :: Ktable
    integer :: ngauss  
    integer :: npress
    integer :: ntemp
    real(dp), allocatable :: kappa(:,:,:,:) ! (ngauss, npress, ntemp, nwav)
    real(dp), allocatable :: weights(:) ! (ngauss)
  end type
  
  type :: ThermodynamicData
    integer :: dtype ! shomate = 1
    integer :: ntemps
    real(dp), allocatable :: temps(:)
    real(dp), allocatable :: data(:,:)
  end type
  
  type :: Atom
    character(:), allocatable :: name
    real(dp) :: mass
  end type
  
  enum, bind(c)
    enumerator :: KTableGas, BackgroundGas
  end enum
  
  type :: Species
    integer :: type
    
    character(:), allocatable :: name
    integer, allocatable :: composition(:) ! (natoms)
    real(dp) :: mass
    
    ! thermodynamics
    type(ThermodynamicData) :: thermo
    
    ! Optical properties
    type(Ktable), allocatable :: k_solar
    type(Ktable), allocatable :: k_ir
    
  end type

  type :: ClimaData
    
    integer :: natoms
    type(Atom), allocatable :: atoms(:)
    
    integer :: ng
    type(Species), allocatable :: sp(:)
        
  end type




end module