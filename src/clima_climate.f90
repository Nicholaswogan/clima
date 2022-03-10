module clima_climate
  use clima_const, only: dp
  use clima_types, only: ClimaData, ClimaVars, ClimaWrk
  implicit none
  private
  
  public :: Climate, radiative_transfer
  
  type :: Climate
    type(ClimaData) :: d
    type(ClimaVars) :: v
    type(ClimaWrk) :: w
  end type
  
  interface Climate
    module procedure :: create_Climate
  end interface
  
contains
  
  function create_Climate(data_dir, species_f, settings_f, star_f, atmosphere_f, err) result(c)
    use clima_const, only: k_boltz
    use clima_types, only: ClimaSettings
    use clima_input, only: create_ClimaVars, create_ClimaData, create_ClimaSettings
    character(*), intent(in) :: data_dir
    character(*), intent(in) :: species_f
    character(*), intent(in) :: settings_f
    character(*), intent(in) :: star_f
    character(*), intent(in) :: atmosphere_f
    character(:), allocatable, intent(out) :: err
    
    type(Climate) :: c
    
    type(ClimaSettings) :: s
    integer :: i
    
    s = create_ClimaSettings(settings_f, err)
    if (allocated(err)) return
    c%d = create_ClimaData(species_f, data_dir, s, err)
    if (allocated(err)) return
    c%v = create_ClimaVars(atmosphere_f, star_f, c%d, err)
    if (allocated(err)) return
    
    ! after file read-in
    allocate(c%v%density(c%v%nz))
    allocate(c%v%densities(c%v%nz,c%d%ng))
    c%v%density = c%v%P/(k_boltz*c%v%T)
    do i = 1,c%d%ng
      c%v%densities(:,i) = c%v%mix(:,i)*c%v%density
    enddo
    
  end function
  
  subroutine radiative_transfer(d, v)
    type(ClimaData), intent(inout) :: d
    type(ClimaVars), intent(in) :: v
    
    ! Solar radiative transfer
    call solar_rt(d%sol, v)
    
  end subroutine
  
  subroutine solar_rt(op, v)
    use clima_types, only: OpticalProperties, Kcoefficients
    use futils, only: Timer
    type(OpticalProperties), intent(inout) :: op
    type(ClimaVars), intent(in) :: v
    
    integer :: i, j, k, l, n, jj, ll
    integer :: totgauss
    type(Timer) :: tm
    
    ! work
    real(dp) :: val
    type(Kcoefficients), allocatable :: ks(:)
    real(dp), allocatable :: cia(:,:)
    real(dp), allocatable :: ray(:,:)
    real(dp), allocatable :: axs(:,:)
    real(dp), allocatable :: pxs(:,:)
    
    real(dp), allocatable :: tausg(:)
    real(dp), allocatable :: taua(:)
    
    allocate(ks(op%nk))
    do i = 1,op%nk
      ks(i)%ngauss = op%k(i)%ngauss
      allocate(ks(i)%k(v%nz,op%k(i)%ngauss))
    enddo
    
    allocate(cia(v%nz,op%ncia))
    allocate(ray(v%nz,op%nray))
    allocate(axs(v%nz,op%naxs))
    allocate(pxs(v%nz,op%npxs))
    allocate(tausg(v%nz))
    allocate(taua(v%nz))
    
    do l = 1,op%nw
    
      ! interpolate to T and P grid
      ! k-distributions
      do i = 1,op%nk
        do k = 1,op%k(i)%ngauss
          do j = 1,v%nz
            call op%k(i)%log10k(k,l)%evaluate(log10(v%P(j)), v%T(j), val)
            ks(i)%k(j,k) = 10.0_dp**val
          enddo
        enddo
      enddo
    
      ! CIA
      do i = 1,op%ncia
        call interpolate_Xsection(op%cia(i), l, v%P, v%T, cia(:,i))
      enddo
      
      ! Rayleigh
      do i = 1,op%nray
        call interpolate_Xsection(op%ray(i), l, v%P, v%T, ray(:,i))
      enddo
      
      ! Absoprtion xs
      do i = 1,op%naxs
        call interpolate_Xsection(op%axs(i), l, v%P, v%T, axs(:,i))
      enddo
      
      ! Photolysis xs
      do i = 1,op%npxs
        call interpolate_Xsection(op%pxs(i), l, v%P, v%T, pxs(:,i))
      enddo
    
      ! compute tau
      ! rayleigh scattering
      tausg(:) = 0.0_dp
      do i = 1,op%nray
        j = op%ray(i)%sp_ind(1)
        do k = 1,v%nz
          n = v%nz+1-k
          tausg(n) = tausg(n) + ray(k,i)*v%densities(k,j)*v%dz(k)
        enddo
      enddo
      
      ! CIA
      taua(:) = 0.0_dp
      do i = 1,op%ncia
        j = op%cia(i)%sp_ind(1)
        jj = op%cia(i)%sp_ind(2)
        do k = 1,v%nz
          n = v%nz+1-k
          taua(n) = taua(n) + cia(k,i)*v%densities(k,j)*v%densities(k,jj)*v%dz(k)
        enddo
      enddo
      
      ! absorption
      do i = 1,op%naxs
        j = op%cia(i)%sp_ind(1)
        do k = 1,v%nz
          n = v%nz+1-k
          taua(n) = taua(n) + axs(k,i)*v%densities(k,j)*v%dz(k)
        enddo
      enddo
      
      ! photolysis
      do i = 1,op%npxs
        j = op%cia(i)%sp_ind(1)
        do k = 1,v%nz
          n = v%nz+1-k
          taua(n) = taua(n) + pxs(k,i)*v%densities(k,j)*v%dz(k)
        enddo
      enddo
      
      ! k-distribution loop
      totgauss = 1
      do i = 1,op%nk
        totgauss = totgauss*op%k(i)%ngauss
      enddo
      
      
      
    enddo

    
  end subroutine
  
  subroutine interpolate_Xsection(xs, l, P, T, res)
    use clima_types, only: Xsection
    
    type(Xsection), intent(inout) :: xs
    integer, intent(in) :: l
    real(dp), intent(in) :: P(:), T(:)
    real(dp), intent(out) :: res(size(P))
    
    integer :: j
    real(dp) :: val
    
    
    
    if (xs%dim == 0) then
      do j = 1,size(P)
        res(j) = 10.0_dp**xs%log10_xs_0d(l)
      enddo
    elseif (xs%dim == 1) then
      do j = 1,size(P)
        call xs%log10_xs_1d(l)%evaluate(T(j), val)
        res(j) = 10.0_dp**val
      enddo
    elseif (xs%dim == 2) then
      print*,'NOOOOO'
      stop 1
    endif
    
  end subroutine
  
  
end module