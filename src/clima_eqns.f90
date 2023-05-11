module clima_eqns
  use clima_const, only: dp
  implicit none

  interface
    !> A temperature dependent surface albedo
    function temp_dependent_albedo_fcn(T_surf) result(albedo)
      use iso_c_binding, only: c_double
      real(c_double), value, intent(in) :: T_surf !! K
      real(c_double) :: albedo
    end function
    
    subroutine ocean_solubility_fcn(T_surf, P_i, m_i) 
      use iso_c_binding, only: c_double
      real(c_double), value, intent(in) :: T_surf !! K
      real(c_double), intent(in) :: P_i(:) !! bars
      real(c_double), intent(out) :: m_i(:) !! mol/kg
    end subroutine
  end interface
    
contains

  subroutine zenith_angles_and_weights(ngauss, zenith_angles, weights)
    use clima_const, only: pi
    use futils, only: gauss_legendre
    integer, intent(in) :: ngauss
    real(dp), intent(out) :: zenith_angles(ngauss), weights(ngauss)

    real(dp) :: x(ngauss), w(ngauss)
    real(dp) :: mu(ngauss)
    
    call gauss_legendre(x, w)
    ! The below is a change of integration bounds. See
    ! https://en.wikipedia.org/wiki/Gaussian_quadrature#Change_of_interval
    mu(:) = x(:)/2.0_dp + 1.0_dp/2.0_dp
    zenith_angles(:) = acos(mu(:))*180.0_dp/pi
    weights(:) = w(:)/2.0_dp
  end subroutine
  
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
  
  function planck_fcn(nu, T) result(B)
    use clima_const, only: c_light, k_boltz => k_boltz_si, plank
    real(dp), intent(in) :: nu ! (1/s) 
    real(dp), intent(in) :: T ! (K)
    real(dp) :: B ! mW sr^−1 m^−2 Hz^-1
    
    B = 1.0e3_dp*((2.0_dp*plank*nu**3.0_dp)/(c_light**2.0_dp)) * &
        ((1.0_dp)/(exp((plank*nu)/(k_boltz*T)) - 1.0_dp)) 
    
  end function
  
  pure function ten2power(y) result(res)
    real(dp), intent(in) :: y
    real(dp) :: res  
    real(dp), parameter :: c = log(10.0_dp)
    res = exp(y*c)
  end function
  
  pure function heat_capacity_shomate(coeffs, T) result(cp)
    real(dp), intent(in) :: coeffs(7)
    real(dp), intent(in) :: T
    real(dp) :: cp !! J/(mol K)
    
    real(dp) :: TT
    
    TT = T/1000.0_dp
    cp = coeffs(1) + coeffs(2)*TT + coeffs(3)*TT**2 + &
         coeffs(4)*TT**3 + coeffs(5)/TT**2
  end function
  
  pure subroutine heat_capacity_eval(thermo, T, found, cp)
    use clima_types, only: ShomatePolynomial, Nasa9Polynomial
    use clima_types, only: ThermodynamicData
  
    type(ThermodynamicData), intent(in) :: thermo
    real(dp), intent(in) :: T
    logical, intent(out) :: found
    real(dp), intent(out) :: cp
  
    integer :: k
  
    found = .false.
    do k = 1,thermo%ntemps
      if (T >= thermo%temps(k) .and. &
          T <  thermo%temps(k+1)) then
  
        found = .true.
        if (thermo%dtype == ShomatePolynomial) then
          cp = heat_capacity_shomate(thermo%data(1:7,k), T)
        elseif (thermo%dtype == Nasa9Polynomial) then
          ! gibbs_energy = gibbs_energy_nasa9(thermo%data(1:9,k), T)
          found = .false.         
        endif
  
        exit
  
      endif
    enddo
  
  end subroutine
  
  function eddy_for_heat(l, g, T, dTdz, adiabat) result(Kh)
    real(dp), intent(in) :: l, g, T, dTdz, adiabat
    
    real(dp) :: Kh
    
    real(dp) :: eta, a1, a2
    
    eta = 0.1_dp*abs(adiabat)
    
    if (dTdz < -adiabat-eta) then
      Kh = l**2.0_dp*sqrt(-(g/T)*(dTdz + adiabat))
    elseif (-adiabat-eta < dTdz .and. dTdz < -adiabat) then
      a1 = - adiabat - eta
      a2 = - adiabat
      Kh = (l**2.0_dp*sqrt(-(g/T)*(dTdz + adiabat)))* &
            smoother(dTdz, a1, a2, -2.0_dp)
    elseif (-adiabat < dTdz) then
      Kh = 0.0_dp
    endif

  end function
  
  pure function smoother(x, a1, a2, beta) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in) :: a1
    real(dp), intent(in) :: a2
    real(dp), intent(in) :: beta
    
    real(dp) :: res
    real(dp) :: y
    
    y = (1.0_dp/(a2-a1))*(x-a1)
    res = 1.0_dp/(1.0_dp + (y/(1.0_dp-y))**(-beta))
    
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
  
  pure subroutine gravity_z(radius, mass, nz, z, grav)
    use clima_const, only: G_grav
    real(dp), intent(in) :: radius, mass ! radius in cm, mass in grams
    integer, intent(in) :: nz
    real(dp), intent(in) :: z(nz) ! cm
    real(dp), intent(out) :: grav(nz) ! cm/s2
  
    integer :: i
  
    do i = 1, nz              
      grav(i) = gravity(radius, mass, z(i))
    enddo 
  
  end subroutine
  
  pure function gravity(radius, mass, z) result(grav)
    use clima_const, only: G_grav
    real(dp), intent(in) :: radius, mass ! radius in cm, mass in grams
    real(dp), intent(in) :: z ! cm 
    
    real(dp) :: grav ! cm/s2
           
    grav = G_grav * (mass/1.0e3_dp) / ((radius + z)/1.0e2_dp)**2.0_dp
    grav = grav*1.0e2_dp ! convert to cgs
    
  end function
  
  pure subroutine press_and_den(nz, T, grav, Psurf, dz, &
                                mubar, pressure, density)
    use clima_const, only: k_boltz, N_avo
  
    integer, intent(in) :: nz
    real(dp), intent(in) :: T(nz), grav(nz)
    real(dp), intent(in) :: Psurf, dz(nz), mubar(nz)
  
    real(dp), intent(out) :: pressure(nz)
    real(dp), intent(out) :: density(nz)
  
    real(dp) :: T_temp
    integer :: i
  
    ! first layer
    T_temp = T(1)
    pressure(1) = Psurf * exp(-((mubar(1) * grav(1))/(N_avo * k_boltz * T_temp)) * 0.5e0_dp * dz(1))
    density(1) = pressure(1)/(k_boltz * T(1))
    ! other layers
    do i = 2,nz
      T_temp = (T(i) + T(i-1))/2.0_dp
      pressure(i) = pressure(i-1) * exp(-((mubar(i) * grav(i))/(N_avo * k_boltz * T_temp))* dz(i))
      density(i) = pressure(i)/(k_boltz * T(i))
    enddo
  
  end subroutine
  
  pure function rayleigh_vardavas(A, B, Delta, lambda) result(sigray)
    real(dp), intent(in) :: A, B, Delta, lambda
    real(dp) :: sigray
    sigray = 4.577e-21_dp*((6.0_dp+3.0_dp*Delta)/(6.0_dp-7.0_dp*Delta)) * &
            (A*(1.0_dp+B/(lambda*1.0e-3_dp)**2.0_dp))**2.0_dp * &
            (1.0_dp/(lambda*1.0e-3_dp)**4.0_dp)
  end function

  pure function equilibrium_temperature(stellar_radiation, bond_albedo) result(T_eq)
    use clima_const, only: sigma_si
    real(dp), intent(in) :: stellar_radiation !! Total stellar radiation (W/m^2)
    real(dp), intent(in) :: bond_albedo !! Bond albedo
    real(dp) :: T_eq !! Equilibrium Temperature (K)
    T_eq = ((stellar_radiation*(1.0_dp - bond_albedo))/(4.0_dp*sigma_si))**(0.25_dp)
  end function

  pure function skin_temperature(stellar_radiation, bond_albedo) result(T_skin)
    real(dp), intent(in) :: stellar_radiation !! Total stellar radiation (W/m^2)
    real(dp), intent(in) :: bond_albedo !! Bond albedo
    real(dp) :: T_skin !! Skin Temperature (K)
    T_skin = equilibrium_temperature(stellar_radiation, bond_albedo)*(0.5_dp)**(0.25_dp)
  end function

end module