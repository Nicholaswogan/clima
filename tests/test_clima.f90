
program test_clima
  use clima_radtran, only: ClimaRadtranIR, dp
  implicit none

  type(ClimaRadtranIR) :: rad
  character(:), allocatable :: err
  
  integer :: nz
  real(dp), allocatable :: T(:), P(:), densities(:,:), dz(:)
  
  nz = 100
  
  rad = ClimaRadtranIR("../data","../settings.yaml",nz, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  allocate(T(nz), P(nz), densities(nz,rad%ng), dz(nz))
  print*,rad%species_names
  
  
end program