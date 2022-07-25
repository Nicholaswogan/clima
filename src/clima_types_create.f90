
submodule(clima_types) clima_types_create
  use yaml_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
  implicit none
  
contains
  
  module function create_Species(filename, err) result(sp)
    use fortran_yaml_c, only: parse, error_length
    character(*), intent(in) :: filename
    character(:), allocatable, intent(out) :: err 
    
    type(Species) :: sp
    
    character(error_length) :: error
    class(type_node), pointer :: root
    
    root => parse(filename, error=error)
    if (len_trim(error) /= 0) then
      err = trim(error)
      return
    end if
    select type (root)
      class is (type_dictionary)
        call unpack_speciesfile(root, filename, sp, err)
      class default
        err = 'yaml file "'//filename//'" must have dictionaries at the root level.'
    end select
    call root%finalize()
    deallocate(root)  
    if (allocated(err)) return
    
  end function
  
  subroutine unpack_speciesfile(root, filename, sp, err)
    type(type_dictionary), pointer :: root
    character(*), intent(in) :: filename
    type(Species), intent(inout) :: sp
    character(:), allocatable, intent(out) :: err
    
    type(type_list), pointer :: tmp_list
    type(type_dictionary), pointer :: dict
    type(type_key_value_pair), pointer :: key_value_pair
    type(type_list_item), pointer :: item
    type(type_error), allocatable :: io_err
    
    character(:), allocatable :: tmp_str
    integer :: i, j, ind

    !!! atoms !!!
    tmp_list => root%get_list('atoms',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    sp%natoms = tmp_list%size()
    allocate(sp%atoms_names(sp%natoms))
    allocate(sp%atoms_mass(sp%natoms))
    
    j = 1
    item => tmp_list%first
    do while(associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        sp%atoms_names(j) = element%get_string("name",error = io_err)
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        sp%atoms_mass(j) = element%get_real("mass",error = io_err)
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
    tmp_list => root%get_list('species',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
    sp%ng = tmp_list%size()
    allocate(sp%g(sp%ng))
    do j = 1,sp%ng
      allocate(sp%g(j)%composition(sp%natoms))
      sp%g(j)%composition = 0
    enddo
    
    j = 1
    item => tmp_list%first
    do while(associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        ! name
        tmp_str = element%get_string("name",error = io_err)
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        sp%g(j)%name = trim(tmp_str)

        ! composition
        dict => element%get_dictionary("composition",.true.,error = io_err)
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        key_value_pair => dict%first
        do while (associated(key_value_pair))
          ind = findloc(sp%atoms_names,trim(key_value_pair%key), 1)
          if (ind == 0) then
            err = 'The atom "'// trim(key_value_pair%key)//'" in species "'// &
                  sp%g(j)%name//'" is not in the list of atoms.'
            return
          endif
          key_value_pair =>key_value_pair%next
        enddo
        do i=1,sp%natoms
          sp%g(j)%composition(i) =  &
              dict%get_integer(sp%atoms_names(i), default = 0, error = io_err)
          if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        enddo
        
        ! mass
        sp%g(j)%mass = sum(sp%g(j)%composition*sp%atoms_mass)
        
        ! thermodynamics
        call unpack_thermo(element, sp%g(j)%name, filename, sp%g(j)%thermo, err)
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
  
  module function create_AtmosphereFile(atm_file, err) result(atm)
    character(*), intent(in) :: atm_file
    character(:), allocatable, intent(out) :: err
    
    type(AtmosphereFile) :: atm
    
    character(len=10000) :: line
    character(len=s_str_len) :: arr1(1000)
    character(len=s_str_len) :: arr11(1000)
    integer :: io, n, nn, i, ii
    
    atm%filename = atm_file
    
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
  
  module subroutine unpack_atmospherefile(atm, species_names, z, mix, T, P, err)
  
    use futils, only: is_close, interp
    use clima_types, only: AtmosphereFile
  
    type(AtmosphereFile), intent(in) :: atm
    character(*), intent(in) :: species_names(:)
    real(dp), intent(in) :: z(:)
    
    real(dp), intent(out) :: mix(:,:)
    real(dp), intent(out) :: T(:)
    real(dp), intent(out) :: P(:)
    character(:), allocatable, intent(out) :: err
  
    integer :: i, ind, ind1, ierr, ng, nz
    
    ng = size(species_names)
    nz = size(z)
    
    if (size(mix,1) /= nz .or. size(mix,2) /= ng) then
      err = '"mix" has the wrong dimensions in subroutine "unpack_atmospherefile"'
      return
    endif
    
    if (size(T) /= nz) then
      err = '"T" has the wrong dimensions in subroutine "unpack_atmospherefile"'
      return
    endif
    
    if (size(P) /= nz) then
      err = '"P" has the wrong dimensions in subroutine "unpack_atmospherefile"'
      return
    endif
    
    ind1 = findloc(atm%labels,"alt", 1)
    if (ind1 == 0) then
      err = '"alt" was not found in input file "'//trim(atm%filename)//'"'
      return
    endif
    
    do i=1,ng
      ind = findloc(atm%labels,species_names(i), 1)
      if (ind /= 0) then

        call interp(nz, atm%nz, z, atm%columns(ind1,:)*1.0e5_dp, log10(atm%columns(ind,:)), mix(:,i), ierr)
        if (ierr /= 0) then
          err = 'Error interpolating "'//trim(atm%filename)//'"'
          return
        endif
        
      else
        err = 'Species "'//trim(species_names(i))//'" was not found in "'// &
              trim(atm%filename)//'"'
        return
      endif
    enddo
    
    mix = 10.0_dp**mix
    
    ! check mix sums to 1
    do i = 1,nz
      if (.not. is_close(sum(mix(i,:)), 1.0_dp, tol=1.0e-2_dp)) then
        err = 'mixing ratios do not sum to close to 1 in "'//trim(atm%filename)//'"'
      endif
    enddo

    ind = findloc(atm%labels,'temp',1)
    if (ind /= 0) then
      call interp(nz, atm%nz, z, atm%columns(ind1,:)*1.0e5_dp, atm%columns(ind,:), T(:), ierr)
      if (ierr /= 0) then
        err = 'Error interpolating "'//trim(atm%filename)//'"'
        return
      endif
    else
      err = '"temp" was not found in input file "'//trim(atm%filename)//'"'
      return
    endif
    
    ind = findloc(atm%labels,'press',1)
    if (ind /= 0) then
      call interp(nz, atm%nz, z, atm%columns(ind1,:)*1.0e5_dp, log10(atm%columns(ind,:)), P(:), ierr)
      if (ierr /= 0) then
        err = 'Error interpolating "'//trim(atm%filename)//'"'
        return
      endif
    else
      err = '"press" was not found in input file "'//trim(atm%filename)//'"'
      return
    endif
    P = 10.0**P
    
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
  
  module function create_ClimaSettings(filename, err) result(s)
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
    
    s%filename = filename
    
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
    real(dp) :: tmp
    
    ! required
    s%nz = grid%get_integer('number-of-layers',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    ! not required
    tmp = grid%get_real('bottom', error = io_err)
    if (.not. allocated(io_err)) then
      allocate(s%bottom)
      s%bottom = tmp
    else
      deallocate(io_err)
    endif
    
    tmp = grid%get_real('top', error = io_err)
    if (.not. allocated(io_err)) then
      allocate(s%top)
      s%top = tmp
    else
      deallocate(io_err)
    endif
    
  end subroutine
    
  subroutine unpack_settingsplanet(planet, filename, s, err)
    type(type_dictionary), pointer :: planet
    character(*), intent(in) :: filename
    type(ClimaSettings), intent(inout) :: s
    character(:), allocatable, intent(out) :: err
    
    type (type_error), allocatable :: io_err
    real(dp) :: tmp
    character(:), allocatable :: tmp_str
    
    ! not required
    tmp_str = planet%get_string("background-gas", error=io_err)
    if (.not. allocated(io_err)) then
      s%back_gas_name = trim(tmp_str)
      deallocate(tmp_str)
    else
      deallocate(io_err)
    endif
    
    tmp = planet%get_real('surface-pressure',error = io_err)
    if (.not. allocated(io_err)) then
      allocate(s%P_surf)
      s%P_surf = tmp
      if (s%P_surf <= 0.0_dp) then
        err = 'IOError: Planet surface pressure must be greater than zero.'
        return
      endif
    else
      deallocate(io_err)
    endif
    
    ! required
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
    
    ! not required
    tmp = planet%get_real('surface-albedo',error = io_err)
    if (.not. allocated(io_err)) then
      allocate(s%surface_albedo)
      s%surface_albedo = tmp
      if (s%surface_albedo < 0.0_dp) then
        err = 'IOError: Surface albedo must be greater than zero.'
        return
      endif
    else
      deallocate(io_err)
    endif
    
    tmp = planet%get_real('diurnal-averaging-factor',error = io_err)
    if (.not. allocated(io_err)) then
      allocate(s%diurnal_fac)
      s%diurnal_fac = tmp
      if (s%diurnal_fac < 0.0_dp .or. s%diurnal_fac > 1.0_dp) then
        err = 'IOError: diurnal-averaging-factor must be between 0 and 1.'
        return
      endif
    else
      deallocate(io_err)
    endif
    
    tmp = planet%get_real('solar-zenith-angle',error = io_err)
    if (.not. allocated(io_err)) then
      allocate(s%solar_zenith)
      s%solar_zenith = tmp
      if (s%solar_zenith < 0.0_dp .or. s%solar_zenith > 90.0_dp) then
        err = 'IOError: solar zenith must be between 0 and 90.'
        return
      endif
    else
      deallocate(io_err)
    endif
    
  end subroutine
  
  subroutine unpack_settingsopticalproperties(op_prop, filename, s, err)
    type(type_dictionary), pointer :: op_prop
    character(*), intent(in) :: filename
    type(ClimaSettings), intent(inout) :: s
    character(:), allocatable, intent(out) :: err
    
    type(type_dictionary), pointer :: tmp_dict
    type(type_list), pointer :: tmp
    type (type_error), allocatable :: io_err
    integer :: ind

    ! species
    tmp_dict => op_prop%get_dictionary("species", required=.false., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (associated(tmp_dict)) then

      ! gases
      tmp => tmp_dict%get_list("gases",required=.false., error=io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      if (associated(tmp)) then
        call unpack_string_list(filename, tmp, s%gases, err)
        if (allocated(err)) return
        ind = check_for_duplicates(s%gases)
        if (ind /= 0) then
          err = '"'//trim(s%gases(ind))//'" is a duplicate in '//trim(tmp%path)
          return
        endif
      endif

      ! particles
      tmp => tmp_dict%get_list("particles",required=.false., error=io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      if (associated(tmp)) then
        call unpack_string_list(filename, tmp, s%particles, err)
        if (allocated(err)) return
        ind = check_for_duplicates(s%particles)
        if (ind /= 0) then
          err = '"'//trim(s%particles(ind))//'" is a duplicate in '//trim(tmp%path)
          return
        endif
      endif

      ! Check for duplicates between gases and particles
      if (allocated(s%gases) .and. allocated(s%particles)) then
        block
        character(s_str_len), allocatable :: gas_particles(:)
        
        allocate(gas_particles(size(s%gases)+size(s%particles)))
        gas_particles(:size(s%gases)) = s%gases
        gas_particles(size(s%gases)+1:) = s%particles
        ind = check_for_duplicates(gas_particles)
        if (ind /= 0) then
          err = '"'//trim(gas_particles(ind))//'" is in both '// &
                '/optical-properties/species/gases  and /optical-properties/species/particles'
          return
        endif
        end block
      endif

    endif
    
    ! solar
    tmp_dict => op_prop%get_dictionary("solar", required=.false., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (associated(tmp_dict)) then
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      call unpack_settingsopacity(tmp_dict, filename, s%sol, err)
      if (allocated(err)) return
    endif
    
    ! ir
    tmp_dict => op_prop%get_dictionary("ir", required=.false., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (associated(tmp_dict)) then
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      call unpack_settingsopacity(tmp_dict, filename, s%ir, err)
      if (allocated(err)) return
    endif
    
  end subroutine
  
  subroutine unpack_settingsopacity(op_dict, filename, op, err)
    type(type_dictionary), intent(in) :: op_dict
    character(*), intent(in) :: filename
    type(SettingsOpacity), allocatable, intent(out) :: op
    character(:), allocatable, intent(out) :: err
    
    type(type_list), pointer :: tmp
    class(type_node), pointer :: node
    type(type_dictionary), pointer :: opacities
    type (type_error), allocatable :: io_err
    integer :: ind
    logical :: success
    
    if (allocated(op)) deallocate(op)
    allocate(op)
    
    opacities => op_dict%get_dictionary("opacities", required=.true., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    ! k-distributions
    tmp => opacities%get_list("k-distributions",required=.false., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (associated(tmp)) then
      
      ! k-distribution settings
      op%k_method = op_dict%get_string("k-method", error=io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      op%k_method = trim(op%k_method)
      
      if (op%k_method == "RandomOverlapResortRebin") then
        op%nbins = op_dict%get_integer("number-of-bins", error=io_err)
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      elseif (op%k_method == "RandomOverlap") then
        ! do nothing
      else
        err = 'k-method "'//op%k_method//'" in "'//filename//'" is not an option.'
        return
      endif
      
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
    node => opacities%get("photolysis-xs")
    if (associated(node)) then
      select type (node)
      class is (type_list)
        call unpack_string_list(filename, node, op%photolysis_xs, err)
        if (allocated(err)) return
        ind = check_for_duplicates(op%photolysis_xs)
        if (ind /= 0) then
          err = '"'//trim(op%photolysis_xs(ind))//'" is a duplicate in '//trim(node%path)
          return
        endif
      class is (type_scalar)
        allocate(op%photolysis_bool)
        op%photolysis_bool = node%to_logical(default=.true.,success=success)
        if (.not. success) then
          err = 'Failed to convert "'//trim(node%path)//'" to logical'
          return
        endif
      class default
        err = '"'//trim(node%path)//'" must be a list or a scalar.'
        return
      end select
    endif
    
    ! continuum
    node => opacities%get("water-continuum")
    if (associated(node)) then
      op%water_continuum = opacities%get_string("water-continuum", error=io_err)
      op%water_continuum = trim(op%water_continuum)
    endif
    
  end subroutine
  
  pure function check_for_duplicates(str_list) result(ind)
    character(*), intent(in) :: str_list(:)
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
    character(*), allocatable, intent(out) :: str_list(:)
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
  
end submodule

