
submodule(clima_radtran_types) clima_radtran_types_create
  use fortran_yaml_c_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
  implicit none

contains
  
  module subroutine read_stellar_flux(star_file, nw, wavl, photon_scale_factor, photon_flux, err)
    use futils, only: inter2, addpnt
    use clima_const, only: c_light, plank
    
    character(len=*), intent(in) :: star_file
    integer, intent(in) :: nw
    real(dp), intent(in) :: wavl(nw+1)
    real(dp), intent(in) :: photon_scale_factor
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
      ! scale factor for the flux, for changing to different planets in the same
      ! solar system
      photon_flux(i) = photon_flux(i)*photon_scale_factor
    enddo

  end subroutine
  
  module function create_RadiateXSWrk(op, nz) result(rw)

    type(OpticalProperties), target, intent(in) :: op
    integer, intent(in) :: nz
    type(RadiateXSWrk) :: rw
    
    type(Ksettings), pointer :: kset
    
    integer :: i
    
    allocate(rw%ks(op%nk))
    do i = 1,op%nk
      allocate(rw%ks(i)%k(nz,op%k(i)%ngauss))
    enddo
    allocate(rw%cia(nz,op%ncia))
    allocate(rw%axs(nz,op%naxs))
    allocate(rw%pxs(nz,op%npxs))
    allocate(rw%H2O(nz))
    allocate(rw%foreign(nz))
    allocate(rw%w0(nz,op%npart))
    allocate(rw%qext(nz,op%npart))
    allocate(rw%gt(nz,op%npart))
    
    ! if there are k-distributions
    ! then we need to allocate some work arrays
    kset => op%kset
    if (op%nk /= 0) then
      if (kset%k_method == K_RandomOverlap) then
        ! no need to allocate anything
      elseif (kset%k_method == K_RandomOverlapResortRebin) then
        allocate(rw%tau_k(nz,kset%nbin))
        allocate(rw%tau_xy(nz,kset%nbin*op%ngauss_max))
        allocate(rw%wxy(kset%nbin*op%ngauss_max))
        allocate(rw%wxy1(kset%nbin*op%ngauss_max))
        allocate(rw%wxy_e(kset%nbin*op%ngauss_max+1))
        allocate(rw%inds(kset%nbin*op%ngauss_max))
      elseif (kset%k_method == k_AdaptiveEquivalentExtinction) then
        allocate(rw%tau_grey(nz,op%nk))
        allocate(rw%tau_grey_sum(nz))
        allocate(rw%ind_major(nz))
      endif
    endif
    
  end function  
  
  module function create_RadiateZWrk(nz, npart) result(rz)
    integer, intent(in) :: nz, npart
    
    type(RadiateZWrk) :: rz
  
    allocate(rz%tausg(nz))
    allocate(rz%taua(nz))
    allocate(rz%taua_1(nz))

    allocate(rz%tausp(nz))
    allocate(rz%tausp_1(nz,npart))
    allocate(rz%taup(nz))

    allocate(rz%tauc(nz),rz%tausc(nz),rz%w0c(nz),rz%g0c(nz))

    allocate(rz%tau(nz))
    allocate(rz%w0(nz))
    allocate(rz%gt(nz))

    allocate(rz%tau_band(nz))

    allocate(rz%amean(nz+1))
    allocate(rz%amean1(nz+1))
    allocate(rz%amean2(nz+1))
    allocate(rz%fup1(nz+1))
    allocate(rz%fdn1(nz+1))
    allocate(rz%fup2(nz+1))
    allocate(rz%fdn2(nz+1))
    allocate(rz%fup(nz+1))
    allocate(rz%fdn(nz+1))
    allocate(rz%bplanck(nz+1))
  end function
  
  function create_Ksettings(sop) result(kset)
    use clima_types, only: SettingsOpacity
    use clima_eqns, only: weights_to_bins
    use futils, only: gauss_legendre
    
    type(SettingsOpacity), intent(in) :: sop
    type(Ksettings) :: kset

    kset%k_method_name = sop%k_method
    
    ! method for mixing k-distributions.
    if (sop%k_method == "RandomOverlapResortRebin") then
      kset%k_method = k_RandomOverlapResortRebin
      ! We will hard code the bins to be something reasonable
      kset%nbin = 8
      allocate(kset%wbin_e(kset%nbin+1))
      allocate(kset%wbin(kset%nbin))
      kset%wbin(:) = [0.16523105_dp, 0.30976894_dp, 0.30976894_dp, &
                      0.16523105_dp, 0.00869637_dp, 0.01630363_dp, &
                      0.01630363_dp, 0.00869637_dp]
      call weights_to_bins(kset%wbin, kset%wbin_e)
    elseif (sop%k_method == "RandomOverlap") then
      kset%k_method = k_RandomOverlap
    elseif (sop%k_method == "AdaptiveEquivalentExtinction") then
      kset%k_method = k_AdaptiveEquivalentExtinction
    endif
    
  end function
  
  module function create_OpticalProperties(datadir, optype, species_names, &
                                           particle_names, sop, wavelength_bins_file, err) result(op)
    use fortran_yaml_c, only: YamlFile
    use clima_const, only: c_light, s_str_len
    use clima_types, only: SettingsOpacity
    character(*), intent(in) :: datadir
    integer, intent(in) :: optype
    character(*), intent(in) :: species_names(:)
    character(*), intent(in) :: particle_names(:)
    type(SettingsOpacity), intent(in) :: sop
    character(:), allocatable, intent(in) :: wavelength_bins_file
    character(:), allocatable, intent(out) :: err
    
    type(OpticalProperties) :: op
    
    character(:), allocatable :: filename
    integer :: i, j, ind1, ind2
    type(type_dictionary), pointer :: root_dict
    type(type_key_value_pair), pointer :: pair
    character(s_str_len), allocatable :: tmp_str_list(:), cia_list(:)
    logical :: tmp_bool, file_exists
    
    op%op_type = optype
    ! get the bins
    if (allocated(wavelength_bins_file)) then
      filename = wavelength_bins_file ! custom bins file, specified in settings
    else
      filename = datadir//"/kdistributions/bins.h5" ! default
    endif
    call read_wavl(filename, op%op_type, op%wavl, err)
    if (allocated(err)) return 
    
    op%nw = size(op%wavl) - 1
    allocate(op%freq(size(op%wavl)))
    op%freq = c_light/(op%wavl*1.0e-9_dp)

    allocate(op%species_names(size(species_names)))
    op%species_names = species_names
    allocate(op%particle_names(size(particle_names)))
    op%particle_names = particle_names

    !!!!!!!!!!!!!!!!!!!!!!!
    !!! k-distributions !!!
    !!!!!!!!!!!!!!!!!!!!!!!
    tmp_bool = .false.
    if (allocated(sop%k_distributions_bool)) tmp_bool = sop%k_distributions_bool

    if (allocated(sop%k_distributions) .or. tmp_bool) then

      if (tmp_bool) then

        ! Look to see if the files / data exist 
        j = 0
        do i = 1,size(species_names)
          filename = datadir//"/kdistributions/"//trim(species_names(i))//".h5"
          inquire(file=filename, exist=file_exists)
          if (file_exists) j = j + 1
        enddo

        if (j == 0) then
          err = 'No k-distribution data was found, but at least one k-distribution '// &
                'is needed.'
          return
        endif

        ! Make a list of avaliable data files
        op%nk = j
        if (allocated(tmp_str_list)) deallocate(tmp_str_list)
        allocate(tmp_str_list(op%nk))
        j = 1
        do i = 1,size(species_names)
          filename = datadir//"/kdistributions/"//trim(species_names(i))//".h5"
          inquire(file=filename, exist=file_exists)
          if (file_exists) then
            tmp_str_list(j) = species_names(i)
            j = j + 1
          endif
        enddo

      else
        op%nk = size(sop%k_distributions)
        if (allocated(tmp_str_list)) deallocate(tmp_str_list)
        allocate(tmp_str_list(op%nk))
        tmp_str_list = sop%k_distributions
      endif
      
      ! k distributions settings
      op%kset = create_Ksettings(sop)
      
      allocate(op%k(op%nk))
      
      do i = 1,op%nk
        filename = datadir//"/kdistributions/"//trim(tmp_str_list(i))//".h5"
        ind1 = findloc(species_names, trim(tmp_str_list(i)), 1)
        if (ind1 == 0) then
          err = 'Species "'//trim(tmp_str_list(i))//'" in optical property '// &
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

      if (op%kset%k_method == k_AdaptiveEquivalentExtinction) then
        ! we must check that all k-coeff have the same
        ! number k-bins and the same weights
        do i = 2,op%nk
          if (op%k(1)%ngauss /= op%k(i)%ngauss) then
            err = 'For k-method "AdaptiveEquivalentExtinction", all k-coeff bin weights must match.'
            return
          endif

          if (.not.all(op%k(1)%weights(:) == op%k(i)%weights(:))) then
            err = 'For k-method "AdaptiveEquivalentExtinction", all k-coeff bin weights must match.'
            return
          endif
        enddo
      endif
      
    else
      op%nk = 0
    endif

    if (op%nk == 0) then
      err = "You must specify at least one k-distribution in the settings file."
      return
    endif
    
    !!!!!!!!!!!
    !!! CIA !!!
    !!!!!!!!!!!
    tmp_bool = .false.
    if (allocated(sop%cia_bool)) tmp_bool = sop%cia_bool

    if (allocated(sop%cia) .or. tmp_bool) then

      if (tmp_bool) then; block
        character(s_str_len), allocatable :: cia_combos(:)
        character(s_str_len), allocatable :: cia_species(:,:)
        ! Make a list of all possible CIA combos
        allocate(cia_combos(size(species_names)*size(species_names)))
        allocate(cia_species(2,size(species_names)*size(species_names)))
        do i = 1,size(species_names)
          do j = 1,size(species_names)
            cia_combos(j + (i-1)*size(species_names)) = trim(species_names(i))//'-'//trim(species_names(j))
            cia_species(1, j + (i-1)*size(species_names)) = trim(species_names(i))
            cia_species(2, j + (i-1)*size(species_names)) = trim(species_names(j))
          enddo
        enddo

        ! Look to see if the files/data exist. Logic is included that will not include
        ! H2O continuum if water continuum is separately specified.
        j = 0
        do i = 1,size(cia_combos)
          filename = datadir//"/CIA/"//trim(cia_combos(i))//".h5"
          inquire(file=filename, exist=file_exists)
          if (file_exists .and. &
            .not.(allocated(sop%water_continuum) .and. any('H2O' == cia_species(:,i))))  then
            j = j + 1
          endif
        enddo

        ! Make a list of avaliable data files
        op%ncia = j
        if (allocated(cia_list)) deallocate(cia_list)
        allocate(cia_list(op%ncia))
        j = 1
        do i = 1,size(cia_combos)
          filename = datadir//"/CIA/"//trim(cia_combos(i))//".h5"
          inquire(file=filename, exist=file_exists)
          if (file_exists .and. &
            .not.(allocated(sop%water_continuum) .and. any('H2O' == cia_species(:,i))))  then
            cia_list(j) = cia_combos(i)
            j = j + 1
          endif
        enddo

      endblock; else
        op%ncia = size(sop%cia)
        if (allocated(cia_list)) deallocate(cia_list)
        allocate(cia_list(op%ncia))
        cia_list = sop%cia
      endif

      allocate(op%cia(op%ncia))

      do i = 1,op%ncia
        
        j = index(cia_list(i), "-")
        if (j == 0) then
          err = 'missing "-" in CIA species "'//trim(sop%cia(i))//'"'
          return
        endif
        ind1 = findloc(species_names, trim(cia_list(i)(1:j-1)), 1)
        ind2 = findloc(species_names, trim(cia_list(i)(j+1:)), 1)
        if (ind1 == 0) then
          err = 'Species "'//trim(cia_list(i)(1:j-1))//'" in optical property '// &
                '"CIA" is not in the list of species.'
          return
        endif
        if (ind2 == 0) then
          err = 'Species "'//trim(cia_list(i)(j+1:))//'" in optical property '// &
                '"CIA" is not in the list of species.'
          return
        endif
        
        filename = datadir//"/CIA/"//trim(cia_list(i))//".h5"
        op%cia(i) = create_CIAXsection(filename, [ind1, ind2], op%wavl, err)
        if (allocated(err)) return
        
      enddo
    else
      op%ncia = 0
    endif
    
    !!!!!!!!!!!!!!!!
    !!! Rayleigh !!!
    !!!!!!!!!!!!!!!!
    tmp_bool = .false.
    if (allocated(sop%rayleigh_bool)) tmp_bool = sop%rayleigh_bool

    if (allocated(sop%rayleigh) .or. tmp_bool) then; block
      type(YamlFile) :: file
      ! parse the yaml file
      filename = datadir//"/rayleigh/rayleigh.yaml"
      call file%parse(filename, err)
      if (allocated(err)) return
      select type(root => file%root)
      class is (type_dictionary)
        root_dict => root
      class default
        err = 'There is an issue with formatting in "'//filename//'"'
        return
      end select

      if (tmp_bool) then
        i = 0
        pair => root_dict%first
        do while (associated(pair))
          ind1 = findloc(species_names, trim(pair%key), 1)
          if (ind1 /= 0) then
            i = i + 1
          endif
          pair => pair%next
        enddo
        if (allocated(tmp_str_list)) deallocate(tmp_str_list)
        allocate(tmp_str_list(i))
        i = 1
        pair => root_dict%first
        do while (associated(pair))
          ind1 = findloc(species_names, trim(pair%key), 1)
          if (ind1 /= 0) then
            tmp_str_list(i) = trim(pair%key)
            i = i + 1
          endif
          pair => pair%next
        enddo
        op%nray = size(tmp_str_list)
      else
        op%nray = size(sop%rayleigh)
        if (allocated(tmp_str_list)) deallocate(tmp_str_list)
        allocate(tmp_str_list(op%nray))
        tmp_str_list = sop%rayleigh
      endif
      allocate(op%ray(op%nray)) 

      do i = 1,op%nray
        ind1 = findloc(species_names, trim(tmp_str_list(i)), 1)
        if (ind1 == 0) then
          err = 'Species "'//trim(sop%rayleigh(i))//'" in optical property '// &
                '"rayleigh" is not in the list of species.'
          return
        endif
        op%ray(i) = create_RayleighXsection(filename, root_dict, &
                    trim(tmp_str_list(i)), ind1, op%wavl, err)
        if (allocated(err)) return
      enddo
      
      call file%finalize()
    endblock; else
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
    tmp_bool = .false.
    if (allocated(sop%photolysis_bool)) tmp_bool = sop%photolysis_bool
    
    if (allocated(sop%photolysis_xs) .or. tmp_bool) then

      if (tmp_bool) then
        ! need to go see what photolysis data is avaliable
        j = 0
        do i = 1,size(species_names)
          filename = datadir//"/xsections/"//trim(species_names(i))//".h5"
          inquire(file=filename, exist=file_exists)
          if (file_exists) j = j + 1
        enddo
        
        if (allocated(tmp_str_list)) deallocate(tmp_str_list)
        op%npxs = j
        allocate(tmp_str_list(op%npxs))
        j = 1
        do i = 1,size(species_names)
          filename = datadir//"/xsections/"//trim(species_names(i))//".h5"
          inquire(file=filename, exist=file_exists)
          if (file_exists) then
            tmp_str_list(j) = species_names(i)
            j = j + 1
          endif
        enddo

      else
        op%npxs = size(sop%photolysis_xs)
        if (allocated(tmp_str_list)) deallocate(tmp_str_list)
        allocate(tmp_str_list(op%npxs))
        tmp_str_list = sop%photolysis_xs
      endif

      allocate(op%pxs(op%npxs)) 

      do i = 1,op%npxs
        ind1 = findloc(species_names, trim(tmp_str_list(i)), 1)
        if (ind1 == 0) then
          err = 'Species "'//trim(sop%rayleigh(i))//'" in optical property '// &
                '"photolysis-xs" is not in the list of species.'
          return
        endif
        filename = datadir//"/xsections/"//trim(tmp_str_list(i))//".h5"
        op%pxs(i) = create_PhotolysisXsection(filename, trim(tmp_str_list(i)), ind1, op%wavl, err)
        if (allocated(err)) return

      enddo
    else
      op%npxs = 0
    endif

    !!!!!!!!!!!!!!!!!
    !!! Particles !!!
    !!!!!!!!!!!!!!!!!
    if (allocated(sop%particle_xs)) then
      op%npart = size(sop%particle_xs)
      allocate(op%part(op%npart))

      do i = 1,op%npart

        ind1 = findloc(particle_names, sop%particle_xs(i)%name, 1)
        if (ind1 == 0) then
          err = 'Species "'//sop%particle_xs(i)%name//'" in optical property '// &
                '"particle-xs" is not in the list of particles.'
          return
        endif
        
        filename = datadir//"/aerosol_xsections/"//sop%particle_xs(i)%dat//"/mie_"//sop%particle_xs(i)%dat//'.h5'
        op%part(i) = create_ParticleXsection(filename, ind1, sop%particle_xs(i)%dat, op%wavl, err)
        if (allocated(err)) return
        
      enddo
    else
      op%npart = 0
    endif
    
    !!!!!!!!!!!!!!!!!
    !!! Continuum !!!
    !!!!!!!!!!!!!!!!!
    if (allocated(sop%water_continuum)) then

      ! Make sure water continuum is not already accounted for with CIA
      if (allocated(cia_list)) then
        do i = 1,size(cia_list)
          j = index(cia_list(i), "-")
          if ('H2O' == cia_list(i)(1:j-1) .or. 'H2O' == cia_list(i)(j+1:)) then
            err = 'Optical property "water-continuum" is set, but CIA "'//trim(cia_list(i))// &
                  '" is also set. This is not allowed because it would double count opacity.'
            return
          endif
        enddo
      endif

      allocate(op%cont)
      filename = datadir//"/water_continuum/"//trim(sop%water_continuum)//".h5"
      op%cont = create_WaterContinuum(sop%water_continuum, filename, species_names, op%wavl, err)
      if (allocated(err)) return
    endif
    
  end function
  
  subroutine read_wavl(filename, optype, wavl, err)
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

  function create_ParticleXsection(filename, p_ind, dat_name, wavl, err) result(part)
    use h5fortran
    use clima_useful, only: hdf5_file_closer
    use futils, only: interp_discrete_to_bins
    character(*), intent(in) :: filename
    integer, intent(in) :: p_ind
    character(*), intent(in) :: dat_name
    real(dp), intent(in) :: wavl(:)
    character(:), allocatable, intent(out) :: err

    type(ParticleXsection) :: part

    type(hdf5_file), target :: h
    integer(HSIZE_T), allocatable :: dims(:)
    real(dp), allocatable :: wv_tmp(:)
    real(dp), allocatable :: w0_tmp(:,:), qext_tmp(:,:), g0_tmp(:,:)
    real(dp), allocatable :: w0_file(:,:), qext_file(:,:), g0_file(:,:)
    integer :: i, ierr

    part%p_ind = p_ind
    part%dat_name = trim(dat_name)
    
    if (.not.is_hdf5(filename)) then
      err = "Was unable to open mie data file "//trim(filename)
      return
    endif
    
    block
    type(hdf5_file_closer) :: h_closer

    ! Open file
    call h%open(filename,'r')
    h_closer%h => h

    ! Wavelengths
    call check_h5_dataset(h, 'wavelengths', 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("wavelengths", dims)
    allocate(wv_tmp(dims(1)))
    call h%read("wavelengths", wv_tmp)

    ! Radii
    call check_h5_dataset(h, 'radii', 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("radii", dims)
    part%nrad = dims(1)
    allocate(part%radii(dims(1)))
    call h%read("radii", part%radii)

    ! convert from micron to cm
    part%radii = part%radii/1.0e4_dp
    part%r_min = minval(part%radii)
    part%r_max = maxval(part%radii)

    ! w0
    call check_h5_dataset(h, 'w0', 2, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("w0", dims)
    allocate(w0_tmp(dims(1),dims(2)))
    call h%read("w0", w0_tmp)

    ! qext
    call check_h5_dataset(h, 'qext', 2, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("qext", dims)
    allocate(qext_tmp(dims(1),dims(2)))
    call h%read("qext", qext_tmp)

    ! g0
    call check_h5_dataset(h, 'g0', 2, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("g0", dims)
    allocate(g0_tmp(dims(1),dims(2)))
    call h%read("g0", g0_tmp)

    if (size(w0_tmp,1) /= part%nrad .or. size(w0_tmp,2) /= size(wv_tmp)) then
      err = '"w0" has the wrong shape in "'//filename//'"'
      return
    endif
    if (size(qext_tmp,1) /= part%nrad .or. size(qext_tmp,2) /= size(wv_tmp)) then
      err = '"qext" has the wrong shape in "'//filename//'"'
      return
    endif
    if (size(g0_tmp,1) /= part%nrad .or. size(g0_tmp,2) /= size(wv_tmp)) then
      err = '"g0" has the wrong shape in "'//filename//'"'
      return
    endif

    end block

    ! Now lets interpolate a bunch
    allocate(w0_file(part%nrad,size(wavl)-1))
    allocate(qext_file(part%nrad,size(wavl)-1))
    allocate(g0_file(part%nrad,size(wavl)-1))

    do i = 1, part%nrad

      ! w0
      call interp_discrete_to_bins(wavl, wv_tmp, w0_tmp(i,:), w0_file(i,:), 'Constant', err=err)
      if (allocated(err)) return

      ! ext
      call interp_discrete_to_bins(wavl, wv_tmp, qext_tmp(i,:), qext_file(i,:), 'Constant', err=err)
      if (allocated(err)) return

      ! g0
      call interp_discrete_to_bins(wavl, wv_tmp, g0_tmp(i,:), g0_file(i,:), 'Constant', err=err)
      if (allocated(err)) return

    enddo

    allocate(part%w0(size(wavl) - 1))
    allocate(part%qext(size(wavl) - 1))
    allocate(part%gt(size(wavl) - 1))
    do i = 1,size(wavl)-1
      call part%w0(i)%initialize(part%radii, w0_file(:,i), ierr)
      if (ierr /= 0) then
        err = 'Failed to initialize interpolator for "'//filename//'"'
        return
      endif
      call part%qext(i)%initialize(part%radii, qext_file(:,i), ierr)
      if (ierr /= 0) then
        err = 'Failed to initialize interpolator for "'//filename//'"'
        return
      endif
      call part%gt(i)%initialize(part%radii, g0_file(:,i), ierr)
      if (ierr /= 0) then
        err = 'Failed to initialize interpolator for "'//filename//'"'
        return
      endif
    enddo

  end function
  
  function create_WaterContinuum(model, filename, species_names, wavl, err) result(cont)
    use clima_const, only: log10tiny
    use h5fortran
    use futils, only: addpnt, inter2
    
    character(*), intent(in) :: model
    character(*), intent(in) :: filename
    character(*), intent(in) :: species_names(:)
    real(dp), intent(in) :: wavl(:)
    character(:), allocatable, intent(out) :: err
    
    type(WaterContinuum) :: cont
    integer(HSIZE_T), allocatable :: dims(:)
    real(dp), allocatable :: wav_f(:), wav_f_save(:), tmp_xs(:,:)
    real(dp), allocatable :: log10_xs_H2O(:,:)
    real(dp), allocatable :: log10_xs_foreign(:,:)
    integer :: i, k, kk, ierr, ind
    real(dp), parameter :: rdelta = 1.0e-4_dp
    
    type(hdf5_file) :: h

    cont%model = trim(model)
    
    ! Look for H2O
    ind = findloc(species_names,"H2O",1)
    if (ind == 0) then
      err = '"H2O" must be a species to include the "continuum" opacity'
      return
    endif
    if (.not. size(species_names) > 1) then
      err = 'There must be more than 1 species in order to use the "continuum"'// &
            ' opacity'
      return
    endif
    
    cont%LH2O = ind

    if (.not. is_hdf5(filename)) then
      err = 'Continuum "'//model//'" is not avaliable.'
      return
    endif
    
    call h%open(filename,'r')
    
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
    
    call check_h5_dataset(h, "T", 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) then
      call h%close()
      return
    endif
    call h%shape("T", dims)
    cont%ntemp = dims(1)
    allocate(cont%temp(cont%ntemp))
    call h%read("T", cont%temp)
    
    !!!!!!!!!!!
    !!! H2O !!!
    !!!!!!!!!!!
    call check_h5_dataset(h, "log10xs_H2O", 2, H5T_FLOAT_F, filename, err)
    if (allocated(err)) then
      call h%close()
      return
    endif
    call h%shape("log10xs_H2O", dims)
    if (dims(1) /= cont%ntemp) then
      err = '"log10xs_H2O" has a bad dimension in "'//trim(filename)//'"'
      call h%close()
      return
    endif
    allocate(log10_xs_H2O(dims(1), dims(2)+4))
    log10_xs_H2O = 0.0_dp
    call h%read("log10xs_H2O", log10_xs_H2O(:,1:dims(2)))
    
    allocate(tmp_xs(cont%ntemp,size(wavl)))
    
    kk = dims(2) + 4
    do i = 1,cont%ntemp
      k = dims(2)
      call addpnt(wav_f, log10_xs_H2O(i,:), kk, k, wav_f(1)*(1.0_dp-rdelta), log10tiny,ierr)
      call addpnt(wav_f, log10_xs_H2O(i,:), kk, k, 0.0_dp, log10tiny,ierr)
      call addpnt(wav_f, log10_xs_H2O(i,:), kk, k, wav_f(k)*(1.0_dp+rdelta), log10tiny,ierr)
      call addpnt(wav_f, log10_xs_H2O(i,:), kk, k, huge(rdelta), log10tiny,ierr)
      if (ierr /= 0) then
        err = 'Problem interpolating data in "'//trim(filename)//'"'
        call h%close()
        return
      endif
      call inter2(size(wavl), wavl, tmp_xs(i,:), &
                  kk, wav_f, log10_xs_H2O(i,:), ierr)
      if (ierr /= 0) then
        err = 'Problem interpolating data in "'//trim(filename)//'"'
        call h%close()
        return
      endif
      
      wav_f = wav_f_save
    enddo
    
    allocate(cont%log10_xs_H2O(size(wavl) - 1))
    do i = 1,size(wavl)-1
      call cont%log10_xs_H2O(i)%initialize(cont%temp, tmp_xs(:,i), ierr)
      if (ierr /= 0) then
        err = 'Failed to initialize interpolator for "'//filename//'"'
        call h%close()
        return
      endif
    enddo
    deallocate(tmp_xs)
    
    !!!!!!!!!!!!!!!
    !!! foreign !!!
    !!!!!!!!!!!!!!!
    call check_h5_dataset(h, "log10xs_foreign", 2, H5T_FLOAT_F, filename, err)
    if (allocated(err)) then
      call h%close()
      return
    endif
    call h%shape("log10xs_foreign", dims)
    if (dims(1) /= cont%ntemp) then
      err = '"log10xs_foreign" has a bad dimension in "'//trim(filename)//'"'
      call h%close()
      return
    endif
    allocate(log10_xs_foreign(dims(1), dims(2)+4))
    log10_xs_foreign = 0.0_dp
    call h%read("log10xs_foreign", log10_xs_foreign(:,1:dims(2)))
    
    allocate(tmp_xs(cont%ntemp,size(wavl)))
    
    kk = dims(2) + 4
    do i = 1,cont%ntemp
      k = dims(2)
      call addpnt(wav_f, log10_xs_foreign(i,:), kk, k, wav_f(1)*(1.0_dp-rdelta), log10tiny,ierr)
      call addpnt(wav_f, log10_xs_foreign(i,:), kk, k, 0.0_dp, log10tiny,ierr)
      call addpnt(wav_f, log10_xs_foreign(i,:), kk, k, wav_f(k)*(1.0_dp+rdelta), log10tiny,ierr)
      call addpnt(wav_f, log10_xs_foreign(i,:), kk, k, huge(rdelta), log10tiny,ierr)
      if (ierr /= 0) then
        err = 'Problem interpolating data in "'//trim(filename)//'"'
        call h%close()
        return
      endif
      call inter2(size(wavl), wavl, tmp_xs(i,:), &
                  kk, wav_f, log10_xs_foreign(i,:), ierr)
      if (ierr /= 0) then
        err = 'Problem interpolating data in "'//trim(filename)//'"'
        call h%close()
        return
      endif
      
      wav_f = wav_f_save
    enddo
    
    allocate(cont%log10_xs_foreign(size(wavl) - 1))
    do i = 1,size(wavl)-1
      call cont%log10_xs_foreign(i)%initialize(cont%temp, tmp_xs(:,i), ierr)
      if (ierr /= 0) then
        err = 'Failed to initialize interpolator for "'//filename//'"'
        call h%close()
        return
      endif
    enddo
    
    cont%T_min = minval(cont%temp)
    cont%T_max = maxval(cont%temp)
    
    call h%close()
    
  end function
  
  function create_RayleighXsection(filename, dict, sp, sp_ind, wavl, err) result(xs)
    use clima_eqns, only: rayleigh_vardavas
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
  
  function create_CIAXsection(filename, sp_inds, wavl, err) result(xs)
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
      
      allocate(xs%T_min)
      allocate(xs%T_max)
      xs%T_min = minval(xs%temp)
      xs%T_max = maxval(xs%temp)
    
    endif

    call h%close()

  end subroutine
  
  function create_Ktable(filename, sp_ind, optype, wavl, err) result(k)
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
    
    k%log10P_min = minval(k%log10P)
    k%log10P_max = maxval(k%log10P)
    
    k%T_min = minval(k%temp)
    k%T_max = maxval(k%temp)
    
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

  function create_PhotolysisXsection(xsfilename, sp, sp_ind, wavl, err) result(xs)
    use h5fortran
    use clima_useful, only: hdf5_file_closer
    use clima_const, only: log10tiny
    use futils, only: interp_discrete_to_bins
    character(*), intent(in) :: xsfilename
    character(*), intent(in) :: sp
    integer, intent(in) :: sp_ind
    real(dp), intent(in) :: wavl(:)
    character(:), allocatable, intent(out) :: err

    type(Xsection) :: xs
    
    type(hdf5_file), target :: h
    integer(HSIZE_T), allocatable :: dims(:)
    real(dp), allocatable :: wv_tmp(:), xs_tmp(:)

    xs%xs_type = PhotolysisXsection
    allocate(xs%sp_ind(1))
    xs%sp_ind(1) = sp_ind
    xs%dim = 0

    if (.not.is_hdf5(xsfilename)) then
      err = 'Species "'//sp//'" does not have photolysis xsection data'
      return
    endif

    ! Block to close hdf5 file when we leave this function
    block
    type(hdf5_file_closer) :: h_closer

    ! Open file
    call h%open(xsfilename,'r')
    h_closer%h => h

    ! Wavelengths
    call check_h5_dataset(h, 'wavelengths', 1, H5T_FLOAT_F, xsfilename, err)
    if (allocated(err)) return
    call h%shape("wavelengths", dims)
    allocate(wv_tmp(dims(1)))
    call h%read("wavelengths", wv_tmp)

    ! Photoabsorption
    call check_h5_dataset(h, 'photoabsorption', 1, H5T_FLOAT_F, xsfilename, err)
    if (allocated(err)) return
    call h%shape("photoabsorption", dims)
    allocate(xs_tmp(dims(1)))
    call h%read("photoabsorption", xs_tmp)

    xs_tmp = max(xs_tmp,tiny(1.0_dp))
    xs_tmp = log10(xs_tmp)

    ! Interpolate
    allocate(xs%xs_0d(size(wavl)-1))
    call interp_discrete_to_bins(wavl, wv_tmp, xs_tmp, xs%xs_0d, 'FillValue', log10tiny, err)
    if (allocated(err)) return

    xs%xs_0d = 10.0_dp**xs%xs_0d

    end block
    
  end function
  
end submodule

