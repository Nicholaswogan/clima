
module clima_input
  use clima_const, only: dp, s_str_len
  use clima_types, only: ClimaSettings, ClimaData, ClimaVars
  use yaml_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
  implicit none
  private
  
  public :: create_ClimaSettings, create_ClimaVars, create_ClimaData
  public :: create_AtmosphereFile
  

contains
  
  function create_ClimaVars(atm_file, star_file, s, dat, err) result(v)
    use clima_types, only: k_RandomOverlap, k_RandomOverlapResortRebin, AtmosphereFile
    use clima_const, only: pi
    use clima_eqns, only: vertical_grid, gravity
    character(*), intent(in) :: atm_file
    character(*), intent(in) :: star_file
    type(ClimaSettings), intent(in) :: s
    type(ClimaData), intent(in) :: dat
    character(:), allocatable, intent(out) :: err
    
    type(AtmosphereFile) :: atm
    integer :: i, j
    
    type(ClimaVars) :: v
    
    ! read the atmosphere file
    atm = create_AtmosphereFile(atm_file, err)
    if (allocated(err)) return
    
    allocate(v%photons_sol(dat%sol%nw))
    call read_stellar_flux(star_file, dat%sol%nw, dat%sol%wavl, v%photons_sol, err)
    if (allocated(err)) return
    
    !!! AFTER FILE READ IN !!!
    ! settings
    v%nz = s%nz
    v%surface_pressure = s%P_surf
    v%surface_albedo = s%surface_albedo
    v%diurnal_fac = s%diurnal_fac
    v%solar_zenith = s%solar_zenith
    ! allocate memory
    allocate(v%z(v%nz))
    allocate(v%dz(v%nz))
    allocate(v%grav(v%nz))
    allocate(v%mix(v%nz,dat%ng))
    allocate(v%mubar(v%nz))
    allocate(v%T_init(v%nz))
    ! set up vertical grid
    call vertical_grid(s%bottom, s%top, &
                       v%nz, v%z, v%dz)
    call gravity(dat%planet_radius, dat%planet_mass, &
                 v%nz, v%z, v%grav)      
    call unpack_atmospherefile(atm_file, atm, dat, v, err)
    if (allocated(err)) return
    
    ! solar zenith
    v%u0 = cos(v%solar_zenith*pi/180.0_dp)
    
    ! compute mean molecular weight everywhere
    do i = 1,v%nz
      v%mubar(i) = 0
      do j = 1,dat%ng
        v%mubar(i) = v%mubar(i) + v%mix(i,j)*dat%sp(j)%mass
      enddo
    enddo

  end function
  
  function create_AtmosphereFile(atm_file, err) result(atm)
    use clima_types, only: AtmosphereFile
    character(*), intent(in) :: atm_file
    character(:), allocatable, intent(out) :: err
    
    type(AtmosphereFile) :: atm
    
    character(len=10000) :: line
    character(len=s_str_len) :: arr1(1000)
    character(len=s_str_len) :: arr11(1000)
    integer :: io, n, nn, i, ii
    
    open(4, file=trim(atm_file),status='old',iostat=io)
    if (io /= 0) then
      err = 'Can not open file '//trim(atm_file)
      return
    endif
    read(4,'(A)') line
    
    atm%nz = -1
    io = 0
    do while (io == 0)
      read(4,*,iostat=io)
      atm%nz = atm%nz + 1
    enddo
    
    rewind(4)
    read(4,'(A)') line
    n = 0
    nn = 0
    do i=1,1000
      read(line,*,iostat=io) arr1(1:i)
      if (io==-1) exit
      n = n+1
    enddo
    read(4,'(A)') line
    do i=1,1000
      read(line,*,iostat=io) arr11(1:i)
      if (io==-1) exit
      nn = nn+1
    enddo
    if (n /= nn) then
      err = 'There is a missing column label in the file '//trim(atm_file)
      return
    endif
    
    atm%nlabels = n
    allocate(atm%labels(n))
    allocate(atm%columns(n,atm%nz))
    rewind(4)
    
    ! read labels
    read(4,'(A)') line
    read(line,*) (atm%labels(i),i=1,n)
    
    ! First read in all the data into big array
    do i = 1,atm%nz
      read(4,*,iostat=io) (atm%columns(ii,i),ii=1,n)
      if (io /= 0) then
        err = 'Problem reading in initial atmosphere in "'//trim(atm_file)//'"'
        return
      endif
    enddo
    close(4)
    
    
  end function
  
  subroutine unpack_atmospherefile(atm_file, atm, dat, v, err)
  
    use futils, only: is_close, interp
    use clima_types, only: AtmosphereFile
  
    character(*), intent(in) :: atm_file
    type(AtmosphereFile), intent(in) :: atm
    type(ClimaData), intent(in) :: dat
    type(ClimaVars), intent(inout) :: v
    character(:), allocatable, intent(out) :: err
  
    integer :: i, ind, ind1, ierr
  
    ind1 = findloc(atm%labels,"alt", 1)
    if (ind1 == 0) then
      err = '"alt" was not found in input file "'//trim(atm_file)//'"'
      return
    endif
    
    do i=1,dat%ng
      ind = findloc(atm%labels,dat%species_names(i), 1)
      if (ind /= 0) then

        call interp(v%nz, atm%nz, v%z, atm%columns(ind1,:)*1.0e5_dp, log10(atm%columns(ind,:)), v%mix(:,i), ierr)
        if (ierr /= 0) then
          err = 'Error interpolating "'//trim(atm_file)//'"'
          return
        endif
        
      else
        err = 'Species "'//trim(dat%species_names(i))//'" was not found in "'// &
              trim(atm_file)//'"'
        return
      endif
    enddo
    
    v%mix = 10.0_dp**v%mix
    
    ! check mix sums to 1
    do i = 1,v%nz
      if (.not. is_close(sum(v%mix(i,:)), 1.0_dp, tol=1.0e-1_dp)) then
        err = 'mixing ratios do not sum to close to 1 in "'//trim(atm_file)//'"'
      endif
    enddo

    ind = findloc(atm%labels,'temp',1)
    if (ind /= 0) then
      call interp(v%nz, atm%nz, v%z, atm%columns(ind1,:)*1.0e5_dp, atm%columns(ind,:), v%T_init(:), ierr)
      if (ierr /= 0) then
        err = 'Error interpolating "'//trim(atm_file)//'"'
        return
      endif
    else
      err = '"temp" was not found in input file "'//trim(atm_file)//'"'
      return
    endif
    
  end subroutine
  
  subroutine read_stellar_flux(star_file, nw, wavl, photon_flux, err)
    use futils, only: inter2, addpnt
    use clima_const, only: c_light, plank
    
    character(len=*), intent(in) :: star_file
    integer, intent(in) :: nw
    real(dp), intent(in) :: wavl(nw+1)
    real(dp), intent(out) :: photon_flux(nw)
    character(:), allocatable, intent(out) :: err
    
    real(dp), allocatable :: file_wav(:), file_flux(:)
    real(dp) :: flux(nw)
    real(dp) :: dum1, dum2, wavl_av
    integer :: io, i, n, ierr
    real(dp), parameter :: rdelta = 1.0e-4_dp
    
    open(1,file=star_file,status='old',iostat=io)
    if (io /= 0) then
      err = "The input file "//star_file//' does not exist.'
      return
    endif
    
    ! count lines
    n = -1 
    read(1,*)
    do while (io == 0)
      read(1,*,iostat=io) dum1, dum2
      n = n + 1
    enddo
    
    allocate(file_wav(n+4), file_flux(n+4))
    
    ! read data
    rewind(1)
    read(1,*)
    do i = 1,n
      read(1,*,iostat=io) file_wav(i), file_flux(i)
      if (io /= 0) then
        err = "Problem reading "//star_file
        return
      endif
    enddo
    close(1)
    
    i = n
    ! interpolate 
    call addpnt(file_wav, file_flux, n+4, i, file_wav(1)*(1.0_dp-rdelta), 0.0_dp, ierr)
    call addpnt(file_wav, file_flux, n+4, i, 0.0_dp, 0.0_dp, ierr)
    call addpnt(file_wav, file_flux, n+4, i, file_wav(i)*(1.0_dp+rdelta), 0.0_dp,ierr)
    call addpnt(file_wav, file_flux, n+4, i, huge(rdelta), 0.0_dp,ierr)
    if (ierr /= 0) then
      err = "Problem interpolating "//trim(star_file)
      return
    endif

    call inter2(nw+1, wavl, flux, n+4, file_wav, file_flux, ierr)
    if (ierr /= 0) then
      err = "Problem interpolating "//trim(star_file)
      return
    endif

    ! flux is mW/m2/nm
    ! convert to mW/m2/Hz
    ! I use the wavelength average in each bin.
    do i = 1,nw
      wavl_av = 0.5_dp*(wavl(i) + wavl(i+1))
      photon_flux(i) = flux(i)*(((wavl_av*1.0e-9_dp)*wavl_av)/c_light)
    enddo

  end subroutine
  
  function create_ClimaData(spfile, datadir, s, err) result(dat)
    use clima_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_types, only: k_RandomOverlap, k_RandomOverlapResortRebin
    use clima_types, only: Ksettings
    use clima_eqns, only: weights_to_bins
    use stdlib_quadrature, only: gauss_legendre
    character(*), intent(in) :: spfile
    character(*), intent(in) :: datadir
    type(ClimaSettings), intent(in) :: s
    character(:), allocatable, intent(out) :: err
    
    type(ClimaData), target :: dat
    
    type(Ksettings), pointer :: kset
    real(dp), allocatable :: tmp(:) ! dummy variable.
    
    ! Reads species yaml file
    call read_speciesfile(spfile, s, dat, err)
    if (allocated(err)) return
    
    dat%sol = create_OpticalProperties(datadir, SolarOpticalProperties, &
                                      dat%species_names, s%sol, err)
    if (allocated(err)) return
    
    dat%ir = create_OpticalProperties(datadir, IROpticalProperties, &
                                      dat%species_names, s%ir, err)
    if (allocated(err)) return
  
    ! method for mixing k-distributions.
    kset => dat%kset
    if (s%k_method == "RandomOverlapResortRebin") then
      kset%k_method = k_RandomOverlapResortRebin
      kset%nbin = s%nbins
      allocate(kset%wbin_e(kset%nbin+1))
      allocate(kset%wbin(kset%nbin))
      allocate(tmp(kset%nbin))
      call gauss_legendre(tmp, kset%wbin)
      kset%wbin = kset%wbin/2.0_dp
      call weights_to_bins(kset%wbin, kset%wbin_e)
    elseif (s%k_method == "RandomOverlap") then
      kset%k_method = k_RandomOverlap
    endif
    
    ! planet
    if (s%planet_is_present) then
      dat%back_gas_name = s%back_gas_name
      dat%planet_mass = s%planet_mass
      dat%planet_radius = s%planet_radius
    else
      err = '"planet" is required in the settings file.'
      return
    endif
    
  end function
  
  subroutine read_speciesfile(filename, s, dat, err)
    use fortran_yaml_c, only: parse, error_length
    character(*), intent(in) :: filename
    type(ClimaSettings), intent(in) :: s
    type(ClimaData), intent(inout) :: dat
    character(:), allocatable, intent(out) :: err
    
    character(error_length) :: error
    class(type_node), pointer :: root
    
    root => parse(filename, error=error)
    if (len_trim(error) /= 0) then
      err = trim(error)
      return
    end if
    select type (root)
      class is (type_dictionary)
        call unpack_speciesfile(root, s, filename, dat, err)
      class default
        err = 'yaml file "'//filename//'" must have dictionaries at the root level.'
    end select
    call root%finalize()
    deallocate(root)  
    if (allocated(err)) return
    
  end subroutine
  
  subroutine unpack_speciesfile(root, s, filename, dat, err)
    type(type_dictionary), pointer :: root
    character(*), intent(in) :: filename
    type(ClimaSettings), intent(in) :: s
    type(ClimaData), intent(inout) :: dat
    character(:), allocatable, intent(out) :: err
    
    type(type_list), pointer :: atoms, species
    type(type_dictionary), pointer :: dict
    type(type_key_value_pair), pointer :: key_value_pair
    type(type_list_item), pointer :: item
    type(type_error), allocatable :: io_err

    integer :: i, j, ind

    !!! atoms !!!
    atoms => root%get_list('atoms',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    dat%natoms = atoms%size()
    allocate(dat%atoms_names(dat%natoms))
    allocate(dat%atoms_mass(dat%natoms))
    
    j = 1
    item => atoms%first
    do while(associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        dat%atoms_names(j) = element%get_string("name",error = io_err)
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        dat%atoms_mass(j) = element%get_real("mass",error = io_err)
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      class default
        err = '"atoms" in "'//trim(filename)//'" must made of dictionaries.'
        return 
      end select
      j = j + 1
      item => item%next
    enddo
    !!! done with atoms !!!
    
    !!! species !!!
    species => root%get_list('species',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
    dat%ng = species%size()
    allocate(dat%species_names(dat%ng))
    allocate(dat%sp(dat%ng))
    do j = 1,dat%ng
      allocate(dat%sp(j)%composition(dat%natoms))
      dat%sp(j)%composition = 0
    enddo
    
    j = 1
    item => species%first
    do while(associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        ! name
        dat%sp(j)%name = trim(element%get_string("name",error = io_err))
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        dat%species_names(j) = dat%sp(j)%name
        
        ! composition
        dict => element%get_dictionary("composition",.true.,error = io_err)
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        key_value_pair => dict%first
        do while (associated(key_value_pair))
          ind = findloc(dat%atoms_names,trim(key_value_pair%key), 1)
          if (ind == 0) then
            err = 'The atom "'// trim(key_value_pair%key)//'" in species "'// &
                  dat%sp(j)%name//'" is not in the list of atoms.'
            return
          endif
          key_value_pair =>key_value_pair%next
        enddo
        do i=1,dat%natoms
          dat%sp(j)%composition(i) =  &
              dict%get_integer(dat%atoms_names(i), default = 0, error = io_err)
          if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        enddo
        
        ! mass
        dat%sp(j)%mass = sum(dat%sp(j)%composition*dat%atoms_mass)
        
        ! thermodynamics
        call unpack_thermo(element, dat%sp(j)%name, filename, dat%sp(j)%thermo, err)
        if (allocated(err)) return
        
      class default
        err = '"species" in "'//trim(filename)//'" must made of dictionaries.'
        return 
      end select
      j = j + 1
      item => item%next
    enddo
    !!! done with species !!!
    
    
  end subroutine
  
  subroutine unpack_thermo(molecule, molecule_name, infile, thermo, err)
    use clima_types, only: ThermodynamicData
    class(type_dictionary), intent(in) :: molecule
    character(len=*), intent(in) :: molecule_name
    character(len=*), intent(in) :: infile
    
    type(ThermodynamicData), intent(out) :: thermo
    character(:), allocatable, intent(out) :: err
    
    type (type_error), allocatable :: io_err
    class(type_dictionary), pointer :: tmpdict
    class(type_list), pointer :: tmplist
    class(type_list_item), pointer :: item, item1
    character(len=:), allocatable :: model
    logical :: success
    
    integer :: j, k
    
    tmpdict => molecule%get_dictionary("thermo",.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    ! check thermodynamic model
    model = tmpdict%get_string("model",error = io_err)
    if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (model == "Shomate") then
      thermo%dtype = 1
    elseif (model == "NASA9") then
      thermo%dtype = 2
    else
      err = "Thermodynamic data must be in Shomate or NASA9 format for "//trim(molecule_name)
      return
    endif
    
    ! get temperature ranges
    tmplist =>tmpdict%get_list("temperature-ranges",.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    thermo%ntemps = tmplist%size() - 1
    if (thermo%ntemps /= 1 .and. thermo%ntemps /= 2) then
      err = "Problem reading thermodynamic data for  "//trim(molecule_name)
      return
    endif
    allocate(thermo%temps(thermo%ntemps + 1))
    
    j = 1
    item => tmplist%first
    do while (associated(item))
      select type (listitem => item%node)
      class is (type_scalar)
        thermo%temps(j) = listitem%to_real(-1.0_dp,success)
        if (.not. success) then
          err = "Problem reading thermodynamic data for  "//trim(molecule_name)
          return
        endif
      class default
        err = "Problem reading thermodynamic data for "//trim(molecule_name)
        return
      end select
      item => item%next
      j = j + 1
    enddo
    
    ! get data
    tmplist => tmpdict%get_list("data",.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    if (tmplist%size() /= thermo%ntemps) then
      err = "Problem reading thermodynamic data for "//trim(molecule_name)
      return
    endif
    
    if (thermo%dtype == 1) then
      ! Shomate
      allocate(thermo%data(7,thermo%ntemps))
    elseif (thermo%dtype == 2) then
      ! NASA9
      allocate(thermo%data(9,thermo%ntemps))
    endif
    
    k = 1
    item => tmplist%first
    do while (associated(item))
      select type (listitem => item%node)
      class is (type_list)
        
        if (listitem%size() /= size(thermo%data,1)) then
          err = "Too much or too little thermodynamic data for "//trim(molecule_name)
          return
        endif
        
        j = 1
        item1 => listitem%first
        do while (associated(item1)) 
          select type (listitem1 => item1%node)
          class is (type_scalar)

            thermo%data(j, k) = listitem1%to_real(-1.0_dp,success)
            if (.not.success) then
              err = "Problem reading thermodynamic data for "//trim(molecule_name)
              return
            endif
          class default
            err = "Problem reading thermodynamic data for "//trim(molecule_name)
            return
          end select
        item1 => item1%next
        j = j + 1
        enddo
      class default
        err = "IOError: Problem reading thermodynamic data for "//trim(molecule_name)
        return
      end select
      item => item%next
      k = k + 1
    enddo          
                            
  end subroutine
  
  function create_OpticalProperties(datadir, optype, species_names, sop, err) result(op)
    use fortran_yaml_c, only : parse, error_length
    use clima_const, only: c_light
    use clima_types, only: OpticalProperties, SettingsOpacity
    character(*), intent(in) :: datadir
    integer, intent(in) :: optype
    character(*), intent(in) :: species_names(:)
    type(SettingsOpacity), intent(in) :: sop
    character(:), allocatable, intent(out) :: err
    
    type(OpticalProperties) :: op
    
    character(:), allocatable :: filename
    integer :: i, j, ind1, ind2
    character(error_length) :: error
    class(type_node), pointer :: root
    type(type_dictionary), pointer :: root_dict
    
    op%op_type = optype
    ! get the bins
    filename = datadir//"/kdistributions/bins.h5"
    call read_wavl(filename, op%op_type, op%wavl, err)
    if (allocated(err)) return 
    
    op%nw = size(op%wavl) - 1
    allocate(op%freq(size(op%wavl)))
    op%freq = c_light/(op%wavl*1.0e-9_dp)
    
    !!!!!!!!!!!!!!!!!!!!!!!
    !!! k-distributions !!!
    !!!!!!!!!!!!!!!!!!!!!!!
    if (allocated(sop%k_distributions)) then
      op%nk = size(sop%k_distributions)
      allocate(op%k(op%nk))
      
      do i = 1,op%nk
        filename = datadir//"/kdistributions/"//trim(sop%k_distributions(i))//".h5"
        ind1 = findloc(species_names, trim(sop%k_distributions(i)), 1)
        if (ind1 == 0) then
          err = 'Species "'//trim(sop%k_distributions(i))//'" in optical property '// &
                '"k-distributions" is not in the list of species.'
          return
        endif
        op%k(i) = create_Ktable(filename, ind1, optype, op%wavl, err)
        if (allocated(err)) return
      enddo
      
      ! find maximum number of gauss points
      op%ngauss_max = 0
      do i = 1,op%nk
        op%ngauss_max = max(op%ngauss_max,op%k(i)%ngauss)
      enddo
      
    else
      op%nk = 0
    endif
    
    !!!!!!!!!!!
    !!! CIA !!!
    !!!!!!!!!!!
    if (allocated(sop%cia)) then
      op%ncia = size(sop%cia)
      allocate(op%cia(op%ncia))
      
      do i = 1,op%ncia
        
        j = index(sop%cia(i), "-")
        if (j == 0) then
          err = 'missing "-" in CIA species "'//trim(sop%cia(i))//'"'
          return
        endif
        ind1 = findloc(species_names, trim(sop%cia(i)(1:j-1)), 1)
        ind2 = findloc(species_names, trim(sop%cia(i)(j+1:)), 1)
        if (ind1 == 0) then
          err = 'Species "'//trim(sop%cia(i)(1:j-1))//'" in optical property '// &
                '"CIA" is not in the list of species.'
          return
        endif
        if (ind2 == 0) then
          err = 'Species "'//trim(sop%cia(i)(j+1:))//'" in optical property '// &
                '"CIA" is not in the list of species.'
          return
        endif
        
        filename = datadir//"/CIA/"//trim(sop%cia(i))//".h5"
        op%cia(i) = create_CIAXsection(filename, [ind1, ind2], op%wavl, err)
        if (allocated(err)) return
        
      enddo
    else
      op%ncia = 0
    endif
    
    !!!!!!!!!!!!!!!!
    !!! Rayleigh !!!
    !!!!!!!!!!!!!!!!
    if (allocated(sop%rayleigh)) then
      op%nray = size(sop%rayleigh)
      allocate(op%ray(op%nray)) 
      filename = datadir//"/rayleigh/rayleigh.yaml"
      root => parse(filename, error=error)
      if (len_trim(error) /= 0) then
        err = trim(error)
        return
      end if

      select type(root)
      class is (type_dictionary)
        root_dict => root
      class default
        err = 'There is an issue with formatting in "'//filename//'"'
        call root%finalize()
        deallocate(root)  
        return
      end select
      do i = 1,op%nray
        ind1 = findloc(species_names, trim(sop%rayleigh(i)), 1)
        if (ind1 == 0) then
          err = 'Species "'//trim(sop%rayleigh(i))//'" in optical property '// &
                '"rayleigh" is not in the list of species.'
          exit
        endif
        op%ray(i) = create_RayleighXsection(filename, root_dict, &
                    trim(sop%rayleigh(i)), ind1, op%wavl, err)
        if (allocated(err)) exit
      enddo
      call root%finalize()
      deallocate(root)  
      if (allocated(err)) return
    else
      op%nray = 0
    endif
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Absorption Xsections !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (allocated(sop%absorption_xs)) then
      err = "sop%absorption_xs not implemented"
      return
    else
      op%naxs = 0
    endif
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Photolysis Xsections !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (allocated(sop%photolysis_xs)) then
      err = "sop%photolysis_xs not implemented"
      return
    else
      op%npxs = 0
    endif
    
  end function
  
  subroutine read_wavl(filename, optype, wavl, err)
    use clima_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use h5fortran
    
    character(*), intent(in) :: filename
    integer, intent(in) :: optype
    real(dp), allocatable :: wavl(:)
    character(:), allocatable :: err
    
    type(hdf5_file) :: h
    integer(HSIZE_T), allocatable :: dims(:)
    character(:), allocatable :: optype_str
    
    if (optype == FarUVOpticalProperties) then
      optype_str = "uv_wavl"
    elseif (optype == SolarOpticalProperties) then
      optype_str = "sol_wavl"
    elseif (optype == IROpticalProperties) then
      optype_str = "ir_wavl"
    endif
    
    if (.not. is_hdf5(filename)) then
      err = 'Failed to read "'//filename//'".'
      return
    endif
    
    call h%open(filename,'r')
    
    call check_h5_dataset(h, optype_str, 1, H5T_FLOAT_F, filename//"/"//optype_str, err)
    if (allocated(err)) then
      call h%close()
      return
    endif
    call h%shape(optype_str, dims)
    allocate(wavl(dims(1)))
    call h%read(optype_str, wavl)
    wavl = wavl*1.0e3_dp ! convert to nm

    call h%close()
    
  end subroutine
  
  function create_RayleighXsection(filename, dict, sp, sp_ind, wavl, err) result(xs)
    use clima_types, only: Xsection, RayleighXsection
    character(*), intent(in) :: filename
    type(type_dictionary), intent(in) :: dict
    character(*), intent(in) :: sp
    integer, intent(in) :: sp_ind
    real(dp), intent(in) :: wavl(:)
    character(:), allocatable, intent(out) :: err
    
    type(Xsection) :: xs
    
    integer :: i
    type(type_dictionary), pointer :: tmp1, tmp2
    type (type_error), allocatable :: io_err
    real(dp) :: Delta, A, B
    
    xs%xs_type = RayleighXsection
    allocate(xs%sp_ind(1))
    xs%sp_ind(1) = sp_ind
    xs%dim = 0
    
    tmp1 => dict%get_dictionary(sp, required=.true., error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
    tmp2 => tmp1%get_dictionary("data", required=.true., error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    Delta = tmp2%get_real("Delta",error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    A = tmp2%get_real("A",error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    B = tmp2%get_real("B",error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    ! compute xsection for all lamda
    allocate(xs%xs_0d(size(wavl)-1))
    do i = 1,size(wavl)-1
      xs%xs_0d(i) = rayleigh_vardavas(A, B, Delta, wavl(i))
    enddo
  
  end function
  
  pure function rayleigh_vardavas(A, B, Delta, lambda) result(sigray)
    real(dp), intent(in) :: A, B, Delta, lambda
    real(dp) :: sigray
    sigray = 4.577e-21_dp*((6.0_dp+3.0_dp*Delta)/(6.0_dp-7.0_dp*Delta)) * &
            (A*(1.0_dp+B/(lambda*1.0e-3_dp)**2.0_dp))**2.0_dp * &
            (1.0_dp/(lambda*1.0e-3_dp)**4.0_dp)
  end function
  
  function create_CIAXsection(filename, sp_inds, wavl, err) result(xs)
    use clima_types, only: Xsection, CIAXsection
    character(*), intent(in) :: filename
    integer, intent(in) :: sp_inds(:)
    real(dp), intent(in) :: wavl(:)
    character(:), allocatable, intent(out) :: err
    
    type(Xsection) :: xs
    
    xs%xs_type = CIAXsection
    allocate(xs%sp_ind(2))
    xs%sp_ind = sp_inds
    call read_h5_Xsection(filename, wavl, xs, err)
  
  end function

  subroutine read_h5_Xsection(filename, wavl, xs, err)
    use clima_const, only: log10tiny
    use clima_types, only: Xsection
    use futils, only: addpnt, inter2
    use h5fortran
    character(*), intent(in) :: filename
    real(dp), intent(in) :: wavl(:)
    type(Xsection), intent(inout) :: xs
    character(:), allocatable, intent(out) :: err
    
    type(hdf5_file) :: h
    integer(HSIZE_T), allocatable :: dims(:)
    real(dp), allocatable :: wav_f(:), wav_f_save(:), tmp_xs(:,:)
    real(dp), allocatable :: log10_xs_0d(:)
    real(dp), allocatable :: log10_xs_1d(:,:)
    ! real(dp), allocatable :: log10_xs_2d(:,:,:)
    integer :: i, k, kk, ierr
    real(dp), parameter :: rdelta = 1.0e-4_dp
    
    if (.not. is_hdf5(filename)) then
      err = 'Failed to read "'//filename//'".'
      return
    endif
    
    call h%open(filename,'r')
    
    if (.not. h%exists('log10xs')) then
      call h%close()
      err = filename//': dataset "log10xs" does not exist'
      return
    endif
    
    xs%dim = h%ndims('log10xs') - 1
    
    if (.not. any(xs%dim == [0, 1])) then
      call h%close()
      err = "Issue reading "//filename
      return
    endif
    
    if (xs%dim == 0 .or. xs%dim == 1) then
      call check_h5_dataset(h, "wavelengths", 1, H5T_FLOAT_F, filename, err)
      if (allocated(err)) then
        call h%close()
        return
      endif
      call h%shape("wavelengths", dims)
      allocate(wav_f(dims(1)+4))
      wav_f = 0.0_dp
      call h%read("wavelengths", wav_f(1:dims(1)))
      wav_f = wav_f*1.0e3_dp ! convert from um to nm
      wav_f_save = wav_f
    endif
    
    if (xs%dim == 1) then
      call check_h5_dataset(h, "T", 1, H5T_FLOAT_F, filename, err)
      if (allocated(err)) then
        call h%close()
        return
      endif
      call h%shape("T", dims)
      allocate(xs%ntemp)
      xs%ntemp = dims(1)
      allocate(xs%temp(xs%ntemp))
      call h%read("T", xs%temp)
    endif
    
    ! read in data and interpolate to grid
    if (xs%dim == 0) then
      call check_h5_dataset(h, "log10xs", 1, H5T_FLOAT_F, filename, err)
      if (allocated(err)) then
        call h%close()
        return
      endif
      call h%shape("log10xs", dims)
      allocate(log10_xs_0d(dims(1)+4))
      call h%read("log10xs", log10_xs_0d(1:dims(1)))
      
      kk = dims(1) + 4
      k = dims(1)
      call addpnt(wav_f, log10_xs_0d, kk, k, wav_f(1)*(1.0_dp-rdelta), log10tiny,ierr)
      call addpnt(wav_f, log10_xs_0d, kk, k, 0.0_dp, log10tiny,ierr)
      call addpnt(wav_f, log10_xs_0d, kk, k, wav_f(k)*(1.0_dp+rdelta), log10tiny,ierr)
      call addpnt(wav_f, log10_xs_0d, kk, k, huge(rdelta), log10tiny,ierr)
      if (ierr /= 0) then
        err = 'Problem interpolating data in "'//trim(filename)//'"'
        call h%close()
        return
      endif
      allocate(xs%xs_0d(size(wavl)-1))
      call inter2(size(wavl), wavl, xs%xs_0d, &
                  kk, wav_f, log10_xs_0d, ierr)
      if (ierr /= 0) then
        err = 'Problem interpolating data in "'//trim(filename)//'"'
        call h%close()
        return
      endif
      xs%xs_0d = 10.0_dp**xs%xs_0d
      
    elseif (xs%dim == 1) then
      call check_h5_dataset(h, "log10xs", 2, H5T_FLOAT_F, filename, err)
      if (allocated(err)) then
        call h%close()
        return
      endif
      call h%shape("log10xs", dims)
      if (dims(1) /= xs%ntemp) then
        err = '"log10xs" has a bad dimension in "'//trim(filename)//'"'
        call h%close()
        return
      endif
      allocate(log10_xs_1d(dims(1), dims(2)+4))
      log10_xs_1d = 0.0_dp
      call h%read("log10xs", log10_xs_1d(:,1:dims(2)))
      
      allocate(tmp_xs(xs%ntemp,size(wavl)))
      
      kk = dims(2) + 4
      do i = 1,xs%ntemp
        k = dims(2)
        call addpnt(wav_f, log10_xs_1d(i,:), kk, k, wav_f(1)*(1.0_dp-rdelta), log10tiny,ierr)
        call addpnt(wav_f, log10_xs_1d(i,:), kk, k, 0.0_dp, log10tiny,ierr)
        call addpnt(wav_f, log10_xs_1d(i,:), kk, k, wav_f(k)*(1.0_dp+rdelta), log10tiny,ierr)
        call addpnt(wav_f, log10_xs_1d(i,:), kk, k, huge(rdelta), log10tiny,ierr)
        if (ierr /= 0) then
          err = 'Problem interpolating data in "'//trim(filename)//'"'
          call h%close()
          return
        endif
        call inter2(size(wavl), wavl, tmp_xs(i,:), &
                    kk, wav_f, log10_xs_1d(i,:), ierr)
        if (ierr /= 0) then
          err = 'Problem interpolating data in "'//trim(filename)//'"'
          call h%close()
          return
        endif
        
        wav_f = wav_f_save
      enddo
      
      allocate(xs%log10_xs_1d(size(wavl) - 1))
      do i = 1,size(wavl)-1
        call xs%log10_xs_1d(i)%initialize(xs%temp, tmp_xs(:,i), ierr)
        if (ierr /= 0) then
          err = 'Failed to initialize interpolator for "'//filename//'"'
          call h%close()
          return
        endif
      enddo

    endif

    call h%close()

  end subroutine
  
  function create_Ktable(filename, sp_ind, optype, wavl, err) result(k)
    use clima_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties, &
                           Ktable
    use h5fortran
    use futils, only: is_close
    use clima_eqns, only: weights_to_bins
    
    character(*), intent(in) :: filename
    integer, intent(in) :: sp_ind
    integer, intent(in) :: optype
    real(dp), intent(in) :: wavl(:)
    character(:), allocatable, intent(out) :: err
    
    type(Ktable) :: k
    
    type(hdf5_file) :: h
    integer(HSIZE_T), allocatable :: dims(:)
    real(dp), allocatable :: wavl_f(:)
    real(dp), allocatable :: log10k(:,:,:,:)
    character(:), allocatable :: read_err
    integer :: i, j, iflag, ind1, ind2

    if (.not. is_hdf5(filename)) then
      err = 'Failed to read "'//filename//'".'
      return
    endif
    
    call h%open(filename,'r')
    
    ! weights
    call check_h5_dataset(h, "weights", 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) then
      call h%close()
      return
    endif
    call h%shape("weights", dims)
    k%ngauss = dims(1)
    allocate(k%weights(dims(1)))
    call h%read("weights", k%weights)
    
    ! weight edges
    allocate(k%weight_e(k%ngauss+1))
    call weights_to_bins(k%weights, k%weight_e)
    
    ! Pressure
    call check_h5_dataset(h, "log10P", 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) then
      call h%close()
      return
    endif
    call h%shape("log10P", dims)
    k%npress = dims(1)
    allocate(k%log10P(dims(1)))
    call h%read("log10P", k%log10P)
    
    ! Temperature
    call check_h5_dataset(h, "T", 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) then
      call h%close()
      return
    endif
    call h%shape("T", dims)
    k%ntemp = dims(1)
    allocate(k%temp(dims(1)))
    call h%read("T", k%temp)
    
    ! check to make sure that the wavelengths match 
    ! our used wavlengths
    call check_h5_dataset(h, "wavelengths", 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) then
      call h%close()
      return
    endif
    call h%shape("wavelengths", dims)
    allocate(wavl_f(dims(1)))
    call h%read("wavelengths", wavl_f)
    wavl_f = wavl_f*1.0e3_dp
    ind1 = minloc(abs(wavl(1)-wavl_f), 1)
    ind2 = minloc(abs(wavl(size(wavl))-wavl_f), 1)
    k%nwav = size(wavl) - 1
    
    if (size(wavl_f) < size(wavl)) then
      call h%close()
      err = filename//": there not enough wavlength bins."
      return
    endif

    if (.not. all(is_close(wavl, wavl_f(ind1:ind2), tol=1.0e-7_dp))) then
      call h%close()
      err = filename//": wavelengths does not match input wavelength bins."
      return
    endif

    ! coefficients
    call check_h5_dataset(h, "log10k", 4, H5T_FLOAT_F, filename, err)
    if (allocated(err)) then
      call h%close()
      return
    endif
    call h%shape("log10k", dims)
    allocate(log10k(dims(1),dims(2),dims(3),dims(4)))
    call h%read("log10k", log10k)

    call h%close()
    
    ! initalize interpolators
    allocate(k%log10k(k%ngauss,k%nwav))
    do i = 1,k%nwav
      do j = 1,k%ngauss
        call k%log10k(j,i)%initialize(k%log10P, k%temp, log10k(j,:,:,ind1-1+i), iflag)
        if (iflag /= 0) then
          if (allocated(read_err)) deallocate(read_err)
          allocate(character(3)::read_err)
          write(read_err,'(i3)') iflag
          err = 'Failed to initialize interpolator for "'//filename//'"'// &
                '. Error code: '//read_err
          return
        endif
      enddo
    enddo
    
    k%sp_ind = sp_ind
    
  end function
  
  subroutine check_h5_dataset(h, dataset, ndims, dtype, prefix, err)
    use h5fortran, only: hdf5_file
    
    type(hdf5_file), intent(in) :: h
    character(*), intent(in) :: dataset
    integer, intent(in) :: ndims
    integer, intent(in) :: dtype
    character(*), intent(in) :: prefix
    character(:), allocatable, intent(out) :: err
    
    if (.not. h%exists(dataset)) then
      err = prefix//': dataset "'//dataset//'" does not exist'
      return
    endif
    
    if (h%ndims(dataset) /= ndims) then
      err = prefix//': dataset "'//dataset//'" has wrong number of dimensions'
      return
    endif
    
    if (h%class(dataset) /= dtype) then
      err = prefix//': dataset "'//dataset//'" has the wrong type'
      return
    endif
    
  end subroutine
  
  
  
  function create_ClimaSettings(filename, err) result(s)
    use fortran_yaml_c, only : parse, error_length
    character(*), intent(in) :: filename
    character(:), allocatable, intent(out) :: err
    
    type(ClimaSettings) :: s
    
    character(error_length) :: error
    class(type_node), pointer :: root
    
    root => parse(filename, error=error)
    if (len_trim(error) /= 0) then
      err = trim(error)
      return
    end if
    select type (root)
      class is (type_dictionary)
        call unpack_ClimaSettings(root, filename, s, err)
      class default
        err = 'yaml file "'//filename//'" must have dictionaries at the root level.'
    end select
    call root%finalize()
    deallocate(root)  
    if (allocated(err)) return
    
  end function
  
  subroutine unpack_ClimaSettings(root, filename, s, err)
    type(type_dictionary), pointer :: root
    character(*), intent(in) :: filename
    type(ClimaSettings), intent(out) :: s
    character(:), allocatable, intent(out) :: err
    
    type(type_dictionary), pointer :: op_prop, planet, grid
    type (type_error), allocatable :: io_err
    
    !!!!!!!!!!!!!!!!!!!!!!!
    !!! atmosphere-grid !!!
    !!!!!!!!!!!!!!!!!!!!!!!
    grid => root%get_dictionary('atmosphere-grid',.false.,error = io_err)
    if (associated(grid)) then
      s%atmos_grid_is_present = .true.
      call unpack_settingsgrid(grid, filename, s, err)
      if (allocated(err)) return
    else
      s%atmos_grid_is_present = .false.
    endif
    
    !!!!!!!!!!!!!!
    !!! planet !!!
    !!!!!!!!!!!!!!
    planet => root%get_dictionary("planet", required=.false., error=io_err)
    if (associated(planet)) then
      s%planet_is_present = .true.
      call unpack_settingsplanet(planet, filename, s, err)
      if (allocated(err)) return
    else
      s%planet_is_present = .false.
    endif
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! optical properties !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    op_prop => root%get_dictionary("optical-properties", required=.true., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
    call unpack_settingsopticalproperties(op_prop, filename, s, err)
    if (allocated(err)) return

  end subroutine
  
  subroutine unpack_settingsgrid(grid, filename, s, err)
    type(type_dictionary), pointer :: grid
    character(*), intent(in) :: filename
    type(ClimaSettings), intent(inout) :: s
    character(:), allocatable, intent(out) :: err
    
    type (type_error), allocatable :: io_err
    
    s%bottom = grid%get_real('bottom',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    s%top = grid%get_real('top',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    s%nz = grid%get_integer('number-of-layers',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
  end subroutine
    
  subroutine unpack_settingsplanet(planet, filename, s, err)
    type(type_dictionary), pointer :: planet
    character(*), intent(in) :: filename
    type(ClimaSettings), intent(inout) :: s
    character(:), allocatable, intent(out) :: err
    
    type (type_error), allocatable :: io_err
    
    s%back_gas_name = trim(planet%get_string("background-gas", error=io_err))
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
    s%P_surf = planet%get_real('surface-pressure',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (s%P_surf <= 0.0_dp) then
      err = 'IOError: Planet surface pressure must be greater than zero.'
      return
    endif
    s%planet_mass = planet%get_real('planet-mass',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (s%planet_mass <= 0.0_dp) then
      err = 'IOError: Planet mass must be greater than zero.'
      return
    endif
    s%planet_radius = planet%get_real('planet-radius',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (s%planet_radius <= 0.0_dp) then
      err = 'IOError: Planet radius must be greater than zero.'
      return
    endif
    s%surface_albedo = planet%get_real('surface-albedo',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (s%surface_albedo < 0.0_dp) then
      err = 'IOError: Surface albedo must be greater than zero.'
      return
    endif
    s%diurnal_fac = planet%get_real('diurnal-averaging-factor',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (s%diurnal_fac < 0.0_dp .or. s%diurnal_fac > 1.0_dp) then
      err = 'IOError: diurnal-averaging-factor must be between 0 and 1.'
      return
    endif
    
    s%solar_zenith = planet%get_real('solar-zenith-angle',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (s%solar_zenith < 0.0_dp .or. s%solar_zenith > 90.0_dp) then
      err = 'IOError: solar zenith must be between 0 and 90.'
      return
    endif
    
  end subroutine
  
  subroutine unpack_settingsopticalproperties(op_prop, filename, s, err)
    type(type_dictionary), pointer :: op_prop
    character(*), intent(in) :: filename
    type(ClimaSettings), intent(inout) :: s
    character(:), allocatable, intent(out) :: err
    
    type(type_dictionary), pointer :: opacities, tmp_dict
    type (type_error), allocatable :: io_err
    
    ! optical property settings
    s%k_method = trim(op_prop%get_string("k-method", error=io_err))
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    if (s%k_method == "RandomOverlapResortRebin") then
      s%nbins = op_prop%get_integer("number-of-bins", error=io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    elseif (s%k_method == "RandomOverlap") then
      ! do nothing
    else
      err = 'k-method "'//s%k_method//'" in "'//filename//'" is not an option.'
      return
    endif
    
    !!!!!!!!!!!!!!!!!
    !!! opacities !!!
    !!!!!!!!!!!!!!!!!
    opacities => op_prop%get_dictionary("opacities", required=.true., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    tmp_dict => opacities%get_dictionary("solar", required=.false., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (associated(tmp_dict)) then
      call unpack_settingsopacity(tmp_dict, filename, s%sol, err)
      if (allocated(err)) return
    endif
    
    tmp_dict => opacities%get_dictionary("ir", required=.false., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (associated(tmp_dict)) then
      call unpack_settingsopacity(tmp_dict, filename, s%ir, err)
      if (allocated(err)) return
    endif
    
  end subroutine
  
  subroutine unpack_settingsopacity(opacities, filename, op, err)
    use clima_types, only: SettingsOpacity
    type(type_dictionary), intent(in) :: opacities
    character(*), intent(in) :: filename
    type(SettingsOpacity), intent(out) :: op
    character(:), allocatable, intent(out) :: err
    
    type(type_list), pointer :: tmp
    class(type_node), pointer :: node
    type (type_error), allocatable :: io_err
    integer :: ind
    logical :: success
    
    ! k-distributions
    tmp => opacities%get_list("k-distributions",required=.false., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (associated(tmp)) then
      call unpack_string_list(filename, tmp, op%k_distributions, err)
      if (allocated(err)) return
      ind = check_for_duplicates(op%k_distributions)
      if (ind /= 0) then
        err = '"'//trim(op%k_distributions(ind))//'" is a duplicate in '//trim(tmp%path)
        return
      endif
    endif
    
    ! CIA
    tmp => opacities%get_list("CIA", required=.false., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (associated(tmp)) then
      call unpack_string_list(filename, tmp, op%cia, err)
      if (allocated(err)) return
      ind = check_for_duplicates(op%cia)
      if (ind /= 0) then
        err = '"'//trim(op%cia(ind))//'" is a duplicate in '//trim(tmp%path)
        return
      endif
    endif
    
    ! rayleigh
    node => opacities%get("rayleigh")
    if (associated(node)) then
      select type (node)
      class is (type_list)
        call unpack_string_list(filename, node, op%rayleigh, err)
        if (allocated(err)) return
        ind = check_for_duplicates(op%rayleigh)
        if (ind /= 0) then
          err = '"'//trim(op%rayleigh(ind))//'" is a duplicate in '//trim(node%path)
          return
        endif
      class is (type_scalar)
        allocate(op%rayleigh_bool)
        op%rayleigh_bool = node%to_logical(default=.true.,success=success)
        if (.not. success) then
          err = 'Failed to convert "'//trim(node%path)//'" to logical'
          return
        endif
      class default
        err = '"'//trim(node%path)//'" must be a list or a scalar.'
        return
      end select
    endif
    
    ! absorption-xs
    tmp => opacities%get_list("absorption-xs", required=.false., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (associated(tmp)) then
      call unpack_string_list(filename, tmp, op%absorption_xs, err)
      if (allocated(err)) return
      ind = check_for_duplicates(op%absorption_xs)
      if (ind /= 0) then
        err = '"'//trim(op%absorption_xs(ind))//'" is a duplicate in '//trim(tmp%path)
        return
      endif
    endif
    
    ! photolysis-xs
    tmp => opacities%get_list("photolysis-xs", required=.false., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (associated(tmp)) then
      call unpack_string_list(filename, tmp, op%photolysis_xs, err)
      if (allocated(err)) return
      ind = check_for_duplicates(op%photolysis_xs)
      if (ind /= 0) then
        err = '"'//trim(op%photolysis_xs(ind))//'" is a duplicate in '//trim(tmp%path)
        return
      endif
    endif

  end subroutine
  
  pure function check_for_duplicates(str_list) result(ind)
    character(len=s_str_len), intent(in) :: str_list(:)
    integer :: ind
    integer :: i, j
    ind = 0
    do i = 1,size(str_list)-1
      do j = i+1,size(str_list)
        if (str_list(i) == str_list(j)) then
          ind = i
          exit
        endif
      enddo
    enddo
  end function
  
  subroutine unpack_string_list(filename, list, str_list, err)
    character(*), intent(in) :: filename
    type(type_list), intent(in) :: list
    character(len=s_str_len), allocatable, intent(out) :: str_list(:)
    character(:), allocatable, intent(out) :: err
    
    integer :: i
    type(type_list_item), pointer :: item
    
    allocate(str_list(list%size()))
    i = 1
    item => list%first
    do while (associated(item))
      select type (it => item%node)
      class is (type_scalar)
        str_list(i) = trim(it%string)
      class default
        err = '"'//trim(it%path)//'" must be a scalar.'
        return
      end select
      i = i + 1
      item => item%next
    enddo
    
  end subroutine
  
end module

