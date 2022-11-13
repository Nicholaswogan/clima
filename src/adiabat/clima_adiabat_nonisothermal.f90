module clima_adiabat_nonisothermal
  use clima_adiabat, only: AdiabatClimate
  use clima_const, only: dp
  use dop853_module, only: dop853_class
  implicit none
  private

  !! This module combines the adiabat model with a traditional
  !! time-stepping radiative equilibrium model in order to
  !! self-consistently compute the tropopause height and 
  !! stratosphere temperature

  type :: RCEData
    real(dp), pointer :: P_i_surf(:)
    integer :: trop_ind

    character(:), allocatable :: err
  end type

  type, extends(dop853_class) :: dop853_custom
    type(AdiabatClimate), pointer :: c => NULL()
    type(RCEData), pointer :: d => NULL()
  end type

contains

  !~~ User facing routine ~~!

  subroutine AdiabatClimate_RCE(self, P_i_surf, P_trop_guess, err)
    class(AdiabatClimate), intent(inout) :: self
    real(dp), target, intent(in) :: P_i_surf(:)
    real(dp), intent(in) :: P_trop_guess
    character(:), allocatable, intent(out) :: err
  end subroutine

  !~~ All of the following is not user-facing ~~!

  subroutine right_hand_side_dop(self, t, u, du)
    class(dop853_class), intent(inout) :: self
    real(dp), intent(in) :: t
    real(dp), intent(in), dimension(:) :: u
    real(dp), intent(out), dimension(:) :: du
    select type (self)
    class is (dop853_custom)
      call right_hand_side(self%c, self%d, t, u, du)
    end select
  end subroutine

  subroutine right_hand_side(c, d, t, u, du)
    class(AdiabatClimate), intent(inout) :: c
    class(RCEData), intent(inout) :: d
    real(dp), intent(in) :: t
    real(dp), intent(in) :: u(:)
    real(dp), intent(out) :: du(:)

    real(dp), allocatable :: T_strat(:)
    real(dp) :: T_surf
    integer :: neq

    ! neq = ground layer + number of cells in stratosphere
    neq = size(u)
    T_surf = u(1)
    allocate(T_strat(neq-1))
    T_strat(:) = u(2:)

    ! Draw a profile from the surface upward
    call c%make_profile(T_surf, d%P_i_surf, d%err)
    if (allocated(d%err)) return

    ! interpolate this profile to the altitude grid

    if (c%P_trop < 0.0_dp) then
      d%err = 'The profile never reached the tropopause'
      return
    endif

    ! now we replace the stratosphere with the input stratosphere profile
    ! c%T(d%trop_ind+1:) = T_strat(:)

    ! Compute pressure and other similar stuff

    ! radiative transfer

    ! rate of change

  end subroutine

  subroutine solout_dop(self, nr, xold, x, y, irtrn, xout)
    class(dop853_class),intent(inout) :: self
    integer,intent(in)                :: nr
    real(dp),intent(in)               :: xold
    real(dp),intent(in)               :: x
    real(dp),dimension(:),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(dp),intent(out)              :: xout

    type(AdiabatClimate), pointer :: c
    type(RCEData), pointer :: d

    select type (self)
    class is (dop853_custom)
      c => self%c
      d => self%d
    end select

    ! We stop if

  end subroutine


end module