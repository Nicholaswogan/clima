program test_kasting
  use clima_adiabat
  implicit none
  
  type(KastingClimateModel) :: c
  integer :: nz
  character(:), allocatable :: err
  
  c = KastingClimateModel('../data', nz, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  
  
  
  
  
end program