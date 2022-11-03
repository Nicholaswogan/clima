module clima_useful
  use clima_const, only: dp
  implicit none
  private

  public :: MinpackHybrd1Vars

  type :: MinpackHybrd1Vars
    integer :: n
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: fvec(:)
    real(dp) :: tol = 1.0e-8_dp
    integer :: info
    real(dp), allocatable :: wa(:)
    integer :: lwa
  end type
  interface MinpackHybrd1Vars
    procedure :: create_MinpackHybrd1Vars
  end interface

contains

  function create_MinpackHybrd1Vars(n, tol) result(v)
    integer, intent(in) :: n
    real(dp), optional, intent(in) :: tol
    type(MinpackHybrd1Vars) :: v
    v%n = n
    allocate(v%x(n))
    allocate(v%fvec(n))
    if (present(tol)) then
      v%tol = tol
    endif
    v%lwa = (n*(3*n+13))/2 + 1
    allocate(v%wa(v%lwa))
  end function

end module