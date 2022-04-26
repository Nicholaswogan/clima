submodule(clima_climate) clima_climate_rhs
  implicit none
  
  
contains
  
  
  
  module subroutine right_hand_side(self, T_in, dTdt)
    use clima_eqns, only: press_and_den
    use clima_radtran, only: radiate
    class(Climate), intent(inout), target :: self
    real(dp), intent(in) :: T_in(:)
    real(dp), intent(out) :: dTdt(:)
    
    integer :: i
    
    type(ClimaData), pointer :: d
    type(ClimaVars), pointer :: v
    type(ClimaWrk), pointer :: w
    
    d => self%d
    v => self%v
    w => self%w
    
    w%rin%T = T_in
    
    ! First, we compute P profile
    call press_and_den(v%nz, w%rin%T, v%grav, &
                       v%surface_pressure*1.e6_dp, v%dz, &
                       v%mubar, w%rin%P, w%rin%density)
    ! convert P from cgs to bar
    w%rin%P = w%rin%P/1.0e6_dp
    ! Compute densities and columns
    do i = 1,d%ng
      w%rin%densities(:,i) = v%mix(:,i)*w%rin%density(:)
      w%rin%cols(:,i) = w%rin%densities(:,i)*v%dz(:)
    enddo
    ! Radiative transfer
    call radiate(d%sol, d%kset, v, w%rin, w%rx_sol, w%rz,  &
                w%fup_a_sol, w%fdn_a_sol, w%fup_sol, w%fdn_sol)
                
    call radiate(d%ir, d%kset, v, w%rin, w%rx_ir, w%rz,  &
                 w%fup_a_ir, w%fdn_a_ir, w%fup_ir, w%fdn_ir)

    ! open(unit=1,file='../../clima_validation/fup_clima_climate.dat',form='formatted',status='replace')
    ! do i = 1,d%ir%nw
    !   write(1,*) (d%ir%freq(i)+d%ir%freq(i+1))*0.5_dp*1.0e-3_dp, w%fup_a_ir(1,i)
    ! enddo
    ! close(1)  
    ! stop
    
    
    
    
  end subroutine
  
end submodule