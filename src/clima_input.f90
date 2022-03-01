
module clima_input
  use clima_const, only: dp
  use clima_types, only: ClimaData, ClimaVars, Ktable
  implicit none
  
contains
  
  function create_ktable(filename, wavenums, err) result(kt)
    
    character(*) :: filename
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
    
    allocate(kt%press(kt%npress))
    read(1,iostat=io) kt%press
    if (io /= 0) then; err = read_err; return; endif
    kt%press = log10(kt%press)
    
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
    
    ! end of file should be reached
    read(1,iostat=io) i
    if (io /= -1) then; err = read_err; return; endif
    
    close(1)
    
    ! build interpolators
    allocate(kt%kappa(kt%ngauss,kt%nwav))
    do i = 1,kt%nwav
      do j = 1,kt%ngauss
        call kt%kappa(j,i)%initialize(kt%press, kt%temp, coeffs(j,:,:,i), iflag)
        if (iflag /= 0) then
          if (allocated(read_err)) deallocate(read_err)
          allocate(character(3)::read_err)
          write(read_err,'(i3)') iflag
          err = 'Failed to initialize 2d interpolator for "'//filename//'"'// &
                '. Error code: '//read_err
          return
        endif
      enddo
    enddo
    
  end function
  
  
  ! function create_settings(err) result(s)
  ! 
  ! 
  ! 
  ! 
  ! 
  ! end function
  
  
  
  
  
end module

