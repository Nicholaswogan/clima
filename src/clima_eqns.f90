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

end module