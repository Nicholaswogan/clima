
module clima_radtran_types_create
  use clima_const, only: dp, s_str_len
  use yaml_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
  implicit none
  private
  
  public :: create_OpticalProperties
  public :: read_stellar_flux
  public :: create_RadiateXSWrk, create_RadiateZWrk
  
contains
  
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
  
  function create_RadiateXSWrk(op, nz) result(rw)
    use clima_radtran_types, only: Ksettings, RadiateXSWrk, OpticalProperties
    use clima_radtran_types, only: K_RandomOverlap, K_RandomOverlapResortRebin
    
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
      endif
    endif
    
  end function  
  
  function create_RadiateZWrk(nz) result(rz)
    use clima_radtran_types, only: RadiateZWrk
    integer, intent(in) :: nz
    
    type(RadiateZWrk) :: rz
  
    allocate(rz%tausg(nz))
    allocate(rz%taua(nz))
    allocate(rz%taua_1(nz))
    allocate(rz%tau(nz))
    allocate(rz%w0(nz))
    allocate(rz%gt(nz))
    allocate(rz%amean(nz+1))
    allocate(rz%fup1(nz+1))
    allocate(rz%fdn1(nz+1))
    allocate(rz%fup(nz+1))
    allocate(rz%fdn(nz+1))
    allocate(rz%bplanck(nz+1))
  end function
  
  function create_OpticalProperties(datadir, optype, species_names, sop, err) result(op)
    use fortran_yaml_c, only : parse, error_length
    use clima_const, only: c_light
    use clima_radtran_types, only: OpticalProperties
    use clima_types, only: SettingsOpacity
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
    use clima_radtran_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
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
    use clima_radtran_types, only: Xsection, RayleighXsection
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
    use clima_radtran_types, only: Xsection, CIAXsection
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
    use clima_radtran_types, only: Xsection
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
    use clima_radtran_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties, &
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
  
end module

