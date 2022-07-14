module clima_climate
  use clima_types, only: Species
  use clima_radtran, only: Radtran
  use clima_const, only: dp, s_str_len
  implicit none
  private

  public :: Climate

  type :: ClimateWrk
    real(dp), allocatable :: T(:) ! (nz) K
    real(dp), allocatable :: P(:) ! (nz) bars
    real(dp), allocatable :: density(:)
    real(dp), allocatable :: densities(:,:) ! (nz,ng) molecules/cm3

    integer :: nsteps_previous = -1
  end type

  type :: Climate
    logical :: switch = .true.

    ! species
    character(s_str_len), allocatable :: species_names(:)
    type(Species) :: sp

    ! radiative transfer
    integer :: nz ! number of layers, copy
    type(Radtran) :: rad

    ! planet
    real(dp) :: planet_mass
    real(dp) :: planet_radius
    real(dp) :: surface_pressure

    ! atmospheric composition
    real(dp), allocatable :: z(:) ! (nz) cm
    real(dp), allocatable :: dz(:) ! (nz) cm
    real(dp), allocatable :: grav(:)
    real(dp), allocatable :: mix(:,:) ! (nz,ng) mixing ratios
    real(dp), allocatable :: mubar(:) ! (nz) mean molecular weight

    ! initial T
    real(dp), allocatable :: T_init(:) ! (nz) K

    ! work
    type(ClimateWrk) :: wrk

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
      class(Climate), target, intent(inout) :: self
      character(*), target, intent(in) :: filename
      real(dp), intent(in) :: tstart
      real(dp), intent(in) :: T_start(:)
      real(dp), target, intent(in) :: t_eval(:)
      logical, intent(in) :: overwrite
      character(:), allocatable, intent(out) :: err

      logical :: success
    end function
  end interface

contains

  function create_Climate(datadir, species_f, settings_f, star_f, atmosphere_f, err) result(c)
    use clima_types, only: ClimaSettings, AtmosphereFile, unpack_atmospherefile
    use clima_eqns, only: vertical_grid, gravity_z
    character(*), intent(in) :: datadir
    character(*), intent(in) :: species_f
    character(*), intent(in) :: settings_f
    character(*), intent(in) :: star_f
    character(*), intent(in) :: atmosphere_f
    character(:), allocatable, intent(out) :: err

    type(Climate) :: c

    integer :: i, j
    type(ClimaSettings) :: s
    type(AtmosphereFile) :: atm
    real(dp), allocatable :: P_dum(:)

    ! create settings
    s = ClimaSettings(settings_f, err)
    if (allocated(err)) return

    ! Check settings for the proper inputs
    call check_Climate_settings(s, err)
    if (allocated(err)) return

    ! unpack settings
    c%nz = s%nz
    c%planet_mass = s%planet_mass
    c%planet_radius = s%planet_radius
    c%surface_pressure = s%P_surf

    ! create species
    c%sp = Species(species_f, err)
    if (allocated(err)) return

    ! make copy of species names
    allocate(c%species_names(c%sp%ng))
    do i = 1,c%sp%ng
      c%species_names(i) = c%sp%g(i)%name
    enddo

    ! create radiative transfer
    c%rad = Radtran(datadir, c%species_names, s, star_f, s%solar_zenith, s%surface_albedo, s%nz, err)
    if (allocated(err)) return

    ! allocate memory
    allocate(c%z(c%nz))
    allocate(c%dz(c%nz))
    allocate(c%grav(c%nz))
    allocate(c%mix(c%nz,c%sp%ng))
    allocate(c%mubar(c%nz))
    allocate(c%T_init(c%nz))
    ! set up vertical grid
    call vertical_grid(s%bottom, s%top, &
                       s%nz, c%z, c%dz)
    call gravity_z(c%planet_radius, c%planet_mass, &
                   c%nz, c%z, c%grav)
    ! Initial atmophere
    atm = AtmosphereFile("../atmosphere.txt", err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    allocate(P_dum(c%nz))
    call unpack_atmospherefile(atm, c%species_names, c%z, c%mix, c%T_init, P_dum, err)
    if (allocated(err)) return
    
    ! compute mean molecular weight everywhere
    do i = 1,c%nz
      c%mubar(i) = 0
      do j = 1,c%sp%ng
        c%mubar(i) = c%mubar(i) + c%mix(i,j)*c%sp%g(j)%mass
      enddo
    enddo

    ! allocate work variables
    allocate(c%wrk%T(c%nz))
    allocate(c%wrk%P(c%nz))
    allocate(c%wrk%density(c%nz))
    allocate(c%wrk%densities(c%nz,c%sp%ng))

  end function

  subroutine check_Climate_settings(s, err)
    use clima_types, only: ClimaSettings
    type(ClimaSettings), intent(in) :: s
    character(:), allocatable, intent(out) :: err

    if (.not. s%atmos_grid_is_present) then
      err = '"'//s%filename//'/atmosphere-grid" does not exist.'
      return
    endif

    if (.not. allocated(s%bottom)) then
      err = '"'//s%filename//'/atmosphere-grid/bottom" does not exist.'
      return
    endif
    
    if (.not. allocated(s%top)) then
      err = '"'//s%filename//'/atmosphere-grid/top" does not exist.'
      return
    endif
    
    if (.not. s%planet_is_present) then
      err = '"'//s%filename//'/planet" does not exist.'
      return
    endif

    if (.not. allocated(s%P_surf)) then
      err = '"'//s%filename//'/planet/surface-pressure" does not exist.'
      return
    endif

    if (.not. allocated(s%surface_albedo)) then
      err = '"'//s%filename//'/planet/surface-albedo" does not exist.'
      return
    endif

    if (.not. allocated(s%solar_zenith)) then
      err = '"'//s%filename//'/planet/solar-zenith-angle" does not exist.'
      return
    endif

  end subroutine

end module