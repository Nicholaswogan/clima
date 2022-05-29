program test_kasting
  use clima_const, only: dp
  use clima_adiabat, only: KastingClimateModel
  use stdlib_math, only: linspace
  implicit none
  
  type(KastingClimateModel) :: c
  integer :: nz, i
  real(dp) :: OLR(300)
  real(dp) :: T_surf(300)
  character(:), allocatable :: err
  
  nz = 50
  
  c = KastingClimateModel('../data', nz, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  do i = 1,size(c%rad%ir%k)
    print*,maxval(c%rad%ir%k(i)%temp)
  enddo
  print*,''
  do i = 1,size(c%rad%ir%cia)
    print*,maxval(c%rad%ir%cia(i)%temp)
  enddo
  print*,''
  

  T_surf = linspace(220.0_dp, 1800.0_dp, 300)


  do i = 1,size(T_surf)
    OLR(i) = c%OLR(T_surf(i), 270.0e6_dp, 0.000000001e6_dp, 1.0e6_dp, err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
  enddo
  
  open(unit=1,file='../runaway1.dat',form='unformatted',status='replace')
  write(1) T_surf
  write(1) OLR
  close(1)

  ! do i =1,nz
  !   print*,c%z(i)*1.0e-5_dp,c%P(i)*1.0e-6_dp,c%T(i),c%f_H2O(i),c%f_CO2(i),c%f_N2(i)
  ! enddo
  ! print*,OLR

end program