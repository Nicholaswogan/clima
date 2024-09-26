module clima_useful
  use clima_const, only: dp
  use minpack_module, only: fcn_hybrj
  use h5fortran, only: hdf5_file
  implicit none
  private

  public :: hdf5_file_closer
  public :: MinpackHybrd1Vars, MinpackHybrj
  public :: linear_solve

  ! Helps close hdf5 files
  type :: hdf5_file_closer
    type(hdf5_file), pointer :: h => null()
  contains
    final :: hdf5_file_closer_final
  end type

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

  type :: MinpackHybrj
    procedure(fcn_hybrj), nopass, pointer :: fcn => null()
    integer :: n
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: fvec(:)
    real(dp), allocatable :: fjac(:,:)
    integer :: Ldfjac
    real(dp) :: xtol
    integer :: maxfev
    real(dp), allocatable :: diag(:)
    integer :: mode
    real(dp) :: factor
    integer :: nprint
    integer :: info
    integer :: nfev
    integer :: njev
    real(dp), allocatable :: r(:)
    integer :: Lr
    real(dp), allocatable :: qtf(:)
    real(dp), allocatable :: wa1(:), wa2(:), wa3(:), wa4(:)
    integer :: lwa
  contains
    procedure :: hybrj => MinpackHybrj_hybrj
    procedure :: code_to_message => MinpackHybrj_code_to_message
  end type
  interface MinpackHybrj
    procedure :: create_MinpackHybrj
  end interface

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      integer, intent(in) :: n
      integer, intent(in) :: nrhs
      real(8), intent(inout) :: a(lda,n)
      integer, intent(in) :: lda
      integer, intent(in) :: ipiv(n)
      real(8), intent(inout) :: b(ldb,nrhs)
      integer, intent(in) :: ldb
      integer, intent(out) :: info
    end subroutine
  end interface

contains

  subroutine hdf5_file_closer_final(self)
    type(hdf5_file_closer), intent(inout) :: self
    if (associated(self%h)) then
      if (self%h%is_open) call self%h%close()
    endif
  end subroutine

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

  function create_MinpackHybrj(fcn, n) result(v)
    procedure(fcn_hybrj) :: fcn
    integer, intent(in) :: n
    type(MinpackHybrj) :: v

    v%fcn => fcn
    v%n = n
    allocate(v%x(v%n))
    v%x = 1.0_dp
    allocate(v%fvec(v%n))
    allocate(v%fjac(v%n,v%n))
    v%ldfjac = v%n
    v%xtol = 1.0e-8_dp ! not default in hybrdj1
    v%maxfev = 100*(v%n + 1)
    allocate(v%diag(v%n))
    v%diag(:) = 1.0_dp
    v%mode = 2
    v%factor = 100.0_dp
    v%nprint = 0
    v%lr = (v%n*(v%n+1))/2
    allocate(v%r(v%lr))
    allocate(v%qtf(v%n))
    allocate(v%wa1(v%n), v%wa2(v%n), v%wa3(v%n), v%wa4(v%n))

  end function

  subroutine MinpackHybrj_hybrj(self)
    use minpack_module, only: hybrj
    class(MinpackHybrj), intent(inout) :: self

    if (.not.associated(self%fcn)) return

    call hybrj(self%fcn, self%n, self%x, self%Fvec, self%Fjac, &
               self%Ldfjac, self%Xtol, self%Maxfev, self%Diag, self%Mode, &
               self%Factor, self%Nprint, self%Info, self%Nfev, self%Njev, &
               self%r, self%Lr, self%Qtf, self%Wa1, self%Wa2, &
               self%Wa3, self%Wa4)
    
  end subroutine

  function MinpackHybrj_code_to_message(self, info) result(message)
    class(MinpackHybrj), intent(inout) :: self
    integer, intent(in) :: info
    character(:), allocatable :: message

    if (info < 0) then
      message = 'user terminated execution.'
    elseif (info == 0) then
      message = 'improper input parameters.'
    elseif (info == 1) then
      message = 'relative error between two consecutive iterates is '// &
      'at most xtol.'
    elseif (info == 2) then
      message = 'number of calls to fcn with iflag = 1 has reached maxfev.'
    elseif (info == 3) then
      message = 'xtol is too small. no further improvement in the '// &
      'approximate solution x is possible.'
    elseif (info == 4) then
      message = 'iteration is not making good progress, as measured '// &
      'by the improvement from the last five jacobian evaluations.'
    elseif (info == 5) then
      message = 'iteration is not making good progress, as measured '// &
      'by the improvement from the last ten iterations.'
    else
      message = 'unknown message'
    endif

  end function

  subroutine linear_solve(A, b, info)
    real(dp), intent(inout) :: a(:,:)
    real(dp), target, intent(inout) :: b(:)
    integer, intent(out) :: info

    integer :: n, nrhs, lda
    integer :: ipiv(size(b))
    integer :: ldb

    n = size(b)
    if (size(a,1) /= n .or. size(a,1) /= n) then
      info = -1
      return
    endif
    lda = n
    ldb = n
    nrhs = 1

    call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
  end subroutine

end module