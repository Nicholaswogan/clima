
module clima_input
  use clima_const, only: dp
  use clima_types, only: ClimaData, ClimaVars, Ktable, CIAtable, Xsection, &
                         ClimaSettings, OpticalProperties, SettingsOpacity, Species
  
  use yaml_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
  implicit none

contains
  
  
  function create_ClimaData(spfile, datadir, s, err) result(dat)
    use clima_const, only: sol_wavenums, ir_wavenums
    use clima_types, only: SolarOpticalProperties, IROpticalProperties
    character(*), intent(in) :: spfile
    character(*), intent(in) :: datadir
    type(ClimaSettings), intent(in) :: s
    character(:), allocatable, intent(out) :: err
    
    type(ClimaData) :: dat
    
    ! Reads species yaml file
    call read_speciesfile(spfile, dat, err)
    if (allocated(err)) return
    
    ! Opacities
    dat%sol = create_OpticalProperties(datadir, SolarOpticalProperties, sol_wavenums, dat%species_names, s%op, err)
    if (allocated(err)) return
    
    dat%ir = create_OpticalProperties(datadir, IROpticalProperties, ir_wavenums, dat%species_names, s%op, err)
    if (allocated(err)) return

    ! other stuff
  end function
  
  subroutine read_speciesfile(filename, dat, err)
    use fortran_yaml_c, only: parse, error_length
    character(*), intent(in) :: filename
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
        call unpack_speciesfile(root, filename, dat, err)
      class default
        err = 'yaml file "'//filename//'" must have dictionaries at the root level.'
    end select
    call root%finalize()
    deallocate(root)  
    if (allocated(err)) return
    
  end subroutine
  
  subroutine unpack_speciesfile(root, filename, dat, err)
    type(type_dictionary), pointer :: root
    character(*), intent(in) :: filename
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
  
  
  function create_OpticalProperties(datadir, optype, wavenums, species_names, sop, err) result(op)
    character(*), intent(in) :: datadir
    integer, intent(in) :: optype
    real(dp), intent(in) :: wavenums(:)
    character(*), intent(in) :: species_names(:)
    type(SettingsOpacity), intent(in) :: sop(:)
    character(:), allocatable, intent(out) :: err
    
    type(OpticalProperties) :: op
    
    character(:), allocatable :: filename
    integer :: i, j, ind1, ind2
    
    op%type = optype
    
    ! wavenums
    op%nw = size(wavenums)-1
    op%wavenums = wavenums
    
    ! count different types of optical properties
    op%nk = 0
    op%ncia = 0
    op%nxs = 0
    op%nray = 0
    do i = 1,size(sop)
      select case (sop(i)%type)
      case("kdistribution")
        op%nk = op%nk + count_opfilename(sop(i), optype)
      case("CIA")
        op%ncia = op%ncia + count_opfilename(sop(i), optype)
      case("xsection")
        op%nxs = op%nxs + count_opfilename(sop(i), optype)
      case("rayleigh")
        op%nray = op%nray + 1
      end select
    enddo
    
    ! allocate
    allocate(op%k(op%nk))
    allocate(op%cia(op%ncia))
    allocate(op%xs(op%nxs))
    allocate(op%sigray(op%nray,op%nw))
    allocate(op%ray_sp_inds(op%nray))
    
    ! get optical properties
    op%nk = 0
    op%ncia = 0
    op%nxs = 0
    op%nray = 0
    do i = 1,size(sop)
      select case (sop(i)%type)
      case("kdistribution")
        j = count_opfilename(sop(i), optype)
        op%nk = op%nk + j
        if (j == 1) then ! if there is data to read
          filename = get_opfilename(sop(i), optype)
          filename = datadir//"/kdistributions/"//filename
          ind1 = findloc(species_names, sop(i)%species(1), 1)
          if (ind1 == 0) then
            err = 'Species "'//trim(sop(i)%species(1))//'" in optical property "'// &
                  sop(i)%name//'" is not in the list of species.'
            return
          endif
          op%k(op%nk) = create_ktable(filename, ind1, op%wavenums, err)
          if (allocated(err)) return
        endif
      case("CIA")
        j = count_opfilename(sop(i), optype)
        op%ncia = op%ncia + j
        if (j == 1) then ! if there is data to read
          filename = get_opfilename(sop(i), optype)
          filename = datadir//"/CIA/"//filename
          ind1 = findloc(species_names, sop(i)%species(1), 1)
          if (ind1 == 0) then
            err = 'Species "'//trim(sop(i)%species(1))//'" in optical property "'// &
                  sop(i)%name//'" is not in the list of species.'
            return
          endif
          ind2 = findloc(species_names, sop(i)%species(2), 1)
          if (ind2 == 0) then
            err = 'Species "'//trim(sop(i)%species(2))//'" in optical property "'// &
                  sop(i)%name//'" is not in the list of species.'
            return
          endif
          op%cia(op%ncia) = create_CIAtable(filename, [ind1,ind2], op%wavenums, err)
          if (allocated(err)) return
        endif
      case("xsection")
        j = count_opfilename(sop(i), optype)
        op%nxs = op%nxs + j
        if (j == 1) then ! if there is data to read
          filename = get_opfilename(sop(i), optype)
          filename = datadir//"/xsections/"//filename
          ind1 = findloc(species_names, sop(i)%species(1), 1)
          if (ind1 == 0) then
            err = 'Species "'//trim(sop(i)%species(1))//'" in optical property "'// &
                  sop(i)%name//'" is not in the list of species.'
            return
          endif
          
          op%xs(op%nxs) = create_Xsection(filename, ind1, op%wavenums, err)
          if (allocated(err)) return
        endif
      case("rayleigh")
        op%nray = op%nray + 1
        
        ind1 = findloc(species_names, sop(i)%species(1), 1)
        if (ind1 == 0) then
          err = 'Species "'//trim(sop(i)%species(1))//'" in optical property "'// &
                sop(i)%name//'" is not in the list of species.'
          return
        endif
        
        op%ray_sp_inds(op%nray) = ind1
        
        ! get rayleigh data
      end select
    enddo
    
    
  end function
  
  function count_opfilename(sop, optype) result(i)
    use clima_types, only: SolarOpticalProperties, IROpticalProperties
    type(SettingsOpacity), intent(in) :: sop
    integer, intent(in) :: optype
    integer :: i
    i = 0
    if (optype == SolarOpticalProperties .and. allocated(sop%solar_filename)) then
      i = 1
    elseif (optype == IROpticalProperties .and. allocated(sop%ir_filename)) then
      i = 1
    endif 
  end function
  
  function get_opfilename(sop, optype) result(filename)
    use clima_types, only: SolarOpticalProperties, IROpticalProperties
    type(SettingsOpacity), intent(in) :: sop
    integer, intent(in) :: optype
    
    character(:), allocatable :: filename
    
    if (optype == SolarOpticalProperties .and. allocated(sop%solar_filename)) then
      filename = sop%solar_filename
    elseif (optype == IROpticalProperties .and. allocated(sop%ir_filename)) then
      filename = sop%ir_filename
    endif 
  
  end function
  
  
  function create_ktable(filename, sp_ind, wavenums, err) result(kt)
    
    character(*), intent(in) :: filename
    integer, intent(in) :: sp_ind
    real(dp), intent(in) :: wavenums(:)
    character(:), allocatable, intent(out) :: err
    
    type(Ktable) :: kt
    
    integer :: i, j, io, iflag
    character(:), allocatable :: read_err
    real(dp), allocatable :: file_wavenums(:)
    real(dp), allocatable :: coeffs(:,:,:,:)
    
    read_err = 'Failed to read "'//filename//'".'
    
    ! Ktables are unformatted binary files
    open(unit=1, file=filename, form='unformatted', status='old', iostat=io)
    if (io /= 0) then
      err = 'Problem opening "'//filename//'".'
      return
    endif
    
    read(1,iostat=io) kt%ngauss
    if (io /= 0) then; err = read_err; return; endif
    read(1,iostat=io) kt%ntemp
    if (io /= 0) then; err = read_err; return; endif
    read(1,iostat=io) kt%npress
    if (io /= 0) then; err = read_err; return; endif
    read(1,iostat=io) kt%nwav
    if (io /= 0) then; err = read_err; return; endif
      
    ! nwav must be equal
    if (kt%nwav /= size(wavenums)-1) then
      err = 'ktable "'//filename// &
            '" has the wrong number of wavelength bins.'
      return
    endif

    allocate(kt%weights(kt%ngauss))
    read(1,iostat=io) kt%weights
    if (io /= 0) then; err = read_err; return; endif
    
    allocate(kt%temp(kt%ntemp))
    read(1,iostat=io) kt%temp
    if (io /= 0) then; err = read_err; return; endif
    
    allocate(kt%log10P(kt%npress))
    read(1,iostat=io) kt%log10P
    if (io /= 0) then; err = read_err; return; endif
    kt%log10P = log10(kt%log10P)
    
    allocate(file_wavenums(kt%nwav+1))
    read(1,iostat=io) file_wavenums
    if (io /= 0) then; err = read_err; return; endif
      
    ! wavenumber must be equal
    if (.not. all(file_wavenums == wavenums)) then
      err = 'ktable "'//filename// &
            '" does not have the correct wavenumber bins.'
      return
    endif
      
    allocate(coeffs(kt%ngauss,kt%npress,kt%ntemp,kt%nwav))
    read(1,iostat=io) coeffs
    if (io /= 0) then; err = read_err; return; endif
    coeffs = log10(coeffs)
    
    ! end of file should be reached
    read(1,iostat=io) i
    if (io /= -1) then; err = read_err; return; endif
    
    close(1)
    
    ! build interpolators
    allocate(kt%log10kappa(kt%ngauss,kt%nwav))
    do i = 1,kt%nwav
      do j = 1,kt%ngauss
        call kt%log10kappa(j,i)%initialize(kt%log10P, kt%temp, coeffs(j,:,:,i), iflag)
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
    
    kt%sp_ind = sp_ind
    
  end function
  
  function create_CIAtable(filename, sp_inds, wavenums, err) result(cia)
    
    character(*), intent(in) :: filename
    integer, intent(in) :: sp_inds(2)
    real(dp), intent(in) :: wavenums(:)
    character(:), allocatable, intent(out) :: err
    
    type(CIAtable) :: cia
    
    integer :: i, io, iflag
    character(:), allocatable :: read_err
    real(dp), allocatable :: file_wavenums(:)
    real(dp), allocatable :: coeffs(:,:)
    
    read_err = 'Failed to read "'//filename//'".'
    
    ! Ktables are unformatted binary files
    open(unit=1, file=filename, form='unformatted', status='old', iostat=io)
    if (io /= 0) then
      err = 'Problem opening "'//filename//'".'
      return
    endif
    
    read(1,iostat=io) cia%ntemp
    if (io /= 0) then; err = read_err; return; endif
    read(1,iostat=io) cia%nwav
    if (io /= 0) then; err = read_err; return; endif
      
    ! nwav must be equal
    if (cia%nwav /= size(wavenums)-1) then
      err = 'ktable "'//filename// &
            '" has the wrong number of wavelength bins.'
      return
    endif

    allocate(cia%temp(cia%ntemp))
    read(1,iostat=io) cia%temp
    if (io /= 0) then; err = read_err; return; endif
    
    allocate(file_wavenums(cia%nwav+1))
    read(1,iostat=io) file_wavenums
    if (io /= 0) then; err = read_err; return; endif
      
    ! wavenumber must be equal
    if (.not. all(file_wavenums == wavenums)) then
      err = 'ktable "'//filename// &
            '" does not have the correct wavenumber bins.'
      return
    endif
      
    allocate(coeffs(cia%ntemp,cia%nwav))
    read(1,iostat=io) coeffs
    if (io /= 0) then; err = read_err; return; endif
    coeffs = log10(coeffs)
    
    ! end of file should be reached
    read(1,iostat=io) i
    if (io /= -1) then; err = read_err; return; endif
    
    close(1)
    
    ! build interpolators
    allocate(cia%log10kappa(cia%nwav))
    do i = 1,cia%nwav
      call cia%log10kappa(i)%initialize(cia%temp, coeffs(:,i), iflag)
      if (iflag /= 0) then
        if (allocated(read_err)) deallocate(read_err)
        allocate(character(3)::read_err)
        write(read_err,'(i3)') iflag
        err = 'Failed to initialize interpolator for "'//filename//'"'// &
              '. Error code: '//read_err
        return
      endif
    enddo
    
    cia%sp_inds = sp_inds
    
  end function
  
  function create_Xsection(filename, sp_ind, wavenums, err) result(xs)
    character(*), intent(in) :: filename
    integer, intent(in) :: sp_ind
    real(dp), intent(in) :: wavenums(:)
    character(:), allocatable, intent(out) :: err
    
    type(Xsection) :: xs
    
    integer :: i, io, iflag
    character(:), allocatable :: read_err
    real(dp), allocatable :: file_wavenums(:)
    real(dp), allocatable :: coeffs(:,:)
    
    read_err = 'Failed to read "'//filename//'".'
    
    ! Ktables are unformatted binary files
    open(unit=1, file=filename, form='formatted', status='old', iostat=io)
    if (io /= 0) then
      err = 'Problem opening "'//filename//'".'
      return
    endif
      
      
      
    
    
    close(1)
    
    
    
    xs%sp_ind = sp_ind
    
    
    
  end function

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
    
    type(type_list), pointer :: opacities
    type(type_list), pointer :: tmp
    type(type_list_item), pointer :: item
    type (type_error), allocatable :: io_err
    integer :: i
    logical :: present
    
    !!!!!!!!!!!!!!!!!
    !!! opacities !!!
    !!!!!!!!!!!!!!!!!
    opacities => root%get_list("opacities", .true., error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    allocate(s%op(opacities%size()))
    i = 1
    item => opacities%first
    do while(associated(item))
      select type (e => item%node)
      class is (type_dictionary)
        s%op(i)%name = trim(e%get_string("name",error = io_err))
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        s%op(i)%type = trim(e%get_string("type",error = io_err))
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
          
        ! filename
        if (s%op(i)%type == "kdistribution" .or. &
            s%op(i)%type == "CIA" .or. &
            s%op(i)%type == "xsection") then
          ! look for filename
          present = .false.
          s%op(i)%solar_filename = trim(e%get_string("solar-filename",error = io_err))
          if (allocated(io_err)) then
            deallocate(io_err)
            deallocate(s%op(i)%solar_filename)
          else
            present = .true. 
          endif
            
          s%op(i)%ir_filename = trim(e%get_string("ir-filename",error = io_err))
          if (allocated(io_err)) then
            deallocate(io_err)
            deallocate(s%op(i)%ir_filename)
          else
            present = .true. 
          endif
          
          if (.not. present) then
            err = 'Opacity "'//s%op(i)%name//'" in file "'//filename//'"'// &
                  ' must contain "solar-filename" and/or "ir-filename" key.'
            return
          endif
        elseif (s%op(i)%type == "rayleigh") then
          ! nothing
        else
          err = 'Opacity "'//s%op(i)%name//'" in file "'//filename//'"'// &
                ' has the following invalid type: "'//s%op(i)%type//'".'
          return
        endif
        
        if (s%op(i)%type == "CIA") then
          ! species key is a list
          tmp => e%get_list("species", .true., error=io_err)
          if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
          
          if (tmp%size() /= 2) then
            err = 'Opacity "'//s%op(i)%name//'" in file "'//filename//'"'// &
                  ' must have a list of "species" that is only 2 in length.'
            return
          endif
          allocate(s%op(i)%species(2))
          select type (ee => tmp%first%node)
          class is (type_scalar)
            s%op(i)%species(1) = trim(ee%string)
          class default
            err = 'Problem reading list of "species" for'// &
                  ' opacity "'//s%op(i)%name//'" in file "'//filename//'".'
            return
          end select
          
          select type (ee => tmp%first%next%node)
          class is (type_scalar)
            s%op(i)%species(2) = trim(ee%string)
          class default
            err = 'Problem reading list of "species" for'// &
                  ' opacity "'//s%op(i)%name//'" in file "'//filename//'".'
            return
          end select
            
        else
          ! species key is a scalar
          allocate(s%op(i)%species(1))
          s%op(i)%species(1) = trim(e%get_string("species",error = io_err))
          if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
          
        endif
        
      class default
        err = 'Opacities in "'//trim(filename)//'" must made of dictionaries.'
        return 
      end select
      i = i + 1
      item => item%next
    enddo
    

  end subroutine
  
  
  
end module

