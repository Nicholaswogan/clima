module clima_eqns
  use clima_const, only: dp
  implicit none
    
contains
  
  subroutine weights_to_bins(weights, bins)
    real(dp), intent(in) :: weights(:)
    real(dp), intent(out) :: bins(:)
    
    integer :: i
    
    bins(1) = 0.0_dp
    do i = 2,size(bins)
      bins(i) = weights(i-1) + bins(i-1)
    enddo
    
  end subroutine
  
  subroutine bins_to_weights(bins, weights)
    real(dp), intent(in) :: bins(:)
    real(dp), intent(out) :: weights(:)
    
    weights = bins(2:) - bins(:size(bins)-1)
    
  end subroutine
  
  ! function planck_fcn(lambda_nm, T) result(B)
  !   use clima_const, only: c_light, k_boltz_si, plank
  !   real(dp), intent(in) :: lambda_nm ! (nm) 
  !   real(dp), intent(in) :: T ! (K)
  !   real(dp) :: B ! W sr^−1 m^−2 m^-1
  ! 
  !   real(dp) :: lambda
  ! 
  !   lambda = lambda_nm*1.0e-9_dp ! convert to meters
  ! 
  !   B = ((2.0_dp*plank*c_light**2.0_dp)/(lambda**5.0_dp)) * &
  !       ((1.0_dp)/(exp((plank*c_light)/(lambda*k_boltz_si*T)) - 1.0_dp)) 
  ! end function
  
  function planck_fcn(nu, T) result(B)
    use clima_const, only: c_light, k_boltz => k_boltz_si, plank
    real(dp), intent(in) :: nu ! (1/s) 
    real(dp), intent(in) :: T ! (K)
    real(dp) :: B ! W sr^−1 m^−2 Hz^-1
    
    B = ((2.0_dp*plank*nu**3.0_dp)/(c_light**2.0_dp)) * &
        ((1.0_dp)/(exp((plank*nu)/(k_boltz*T)) - 1.0_dp)) 
    
  end function
  
  pure function ten2power(y) result(res)
    real(dp), intent(in) :: y
    real(dp) :: res  
    real(dp), parameter :: c = log(10.0_dp)
    res = exp(y*c)
  end function
  
  ! coppied from Photochem
  pure subroutine vertical_grid(bottom, top, nz, z, dz)
    real(dp), intent(in) :: bottom, top
    integer, intent(in) :: nz
    real(dp), intent(out) :: z(nz), dz(nz)
  
    integer :: i
  
    dz = (top - bottom)/nz
    z(1) = dz(1)/2.0_dp
    do i = 2,nz
      z(i) = z(i-1) + dz(i)
    enddo
  end subroutine
  
  pure subroutine gravity(radius, mass, nz, z, grav)
    use clima_const, only: G_grav
    real(dp), intent(in) :: radius, mass ! radius in cm, mass in grams
    integer, intent(in) :: nz
    real(dp), intent(in) :: z(nz) ! cm
    real(dp), intent(out) :: grav(nz) ! cm/s2

    integer :: i
    
    do i = 1, nz              
      grav(i) = G_grav * (mass/1.e3_dp) / ((radius + z(i))/1.e2_dp)**2.0_dp
      grav(i) = grav(i)*1.e2_dp ! convert to cgs
    enddo 
    
  end subroutine

end module