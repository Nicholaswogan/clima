submodule(clima_rc) clima_rc_rhs
  implicit none

contains
  module subroutine RadiativeConvectiveClimate_rhs(self, T_in, dTdt, err)
    class(RadiativeConvectiveClimate), intent(inout) :: self
    real(dp), intent(in) :: T_in(:)
    real(dp), intent(out) :: dTdt(:)
    character(:), allocatable :: err

    integer :: i

    ! Compute mixing ratios
    ! moist species

    ! Compute thickness of each layer

    ! Radiative transfer

    ! Heat capacity

    ! Rate of change

  end subroutine

  module subroutine RadiativeConvectiveClimate_convective_adjustment()
  end subroutine
  
end submodule