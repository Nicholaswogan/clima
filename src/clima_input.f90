
module clima_input
  use clima_const, only: dp
  use clima_types, only: ClimaData, ClimaVars, Ktable, CIAtable, ClimaSettings
  
  use yaml_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
  implicit none

contains
  
  
  function create_ClimaData(spfile, s, err) result(dat)
    character(*), intent(in) :: spfile
    type(ClimaSettings), intent(in) :: s
    character(:), allocatable, intent(out) :: err
    
    type(ClimaData) :: dat
    
    ! reads species yaml file
    call read_speciesfile(spfile, dat, err)
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
    
    type (type_list), pointer :: atoms, species
    type (type_list_item), pointer :: item
    type (type_error), allocatable :: io_err

    integer :: j

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
        err = 'atoms in "'//trim(filename)//'" must made of dictionaries.'
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
        
    
    !!! done with species !!!
    
    
  end subroutine
  
  subroutine get_thermodata(molecule, molecule_name, infile, &
                            thermo, err)
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
          else
            present = .true. 
          endif
            
          s%op(i)%ir_filename = trim(e%get_string("ir-filename",error = io_err))
          if (allocated(io_err)) then
            deallocate(io_err)
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

