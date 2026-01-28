module clima_useful
  use clima_const, only: dp
  use minpack_module, only: fcn_hybrj
  use h5fortran, only: hdf5_file

  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_long, c_ptr, c_null_ptr, c_associated
  use fsundials_nvector_mod, only: N_Vector
  use fsundials_matrix_mod, only: SUNMatrix
  use fsundials_linearsolver_mod, only: SUNLinearSolver
  implicit none
  private

  public :: hdf5_file_closer
  public :: MinpackHybrd1Vars, MinpackHybrj
  public :: linear_solve
  public :: SundialsCVode

  type SundialsCVode
    integer :: neq !! number of ODEs
    procedure(cvode_rhs_fcn), pointer :: f => NULL() !! right-hand-side of ODEs
    procedure(cvode_jac_fcn), pointer :: jac => NULL() !! jacobian of ODEs
    !> cvode memory
    type(c_ptr) :: cvode_mem = c_null_ptr
    ! solution vector
    real(c_double), allocatable :: yvec(:)
    type(N_Vector), pointer :: sunvec_y => NULL()
    ! absolute tolerance
    real(c_double), allocatable :: abstol(:)
    type(N_Vector), pointer :: abstol_nvec => NULL()
    ! matrix and linear solver
    type(SUNMatrix), pointer :: sunmat => NULL()
    type(SUNLinearSolver), pointer :: sunlin => NULL()
  contains
    procedure :: finalize => SundialsCVode_finalize
    procedure :: initialize => SundialsCVode_initialize
    procedure :: integrate => SundialsCVode_integrate
    procedure :: steadystate => SundialsCVode_steadystate
    final :: SundialsCVode_final
  end type

  abstract interface

    !> Interface for the right-hand-side function defining the system
    !> of ODEs.
    function cvode_rhs_fcn(self, t, y, ydot) result(ierr)
      import :: dp, SundialsCVode
      implicit none
      class(SundialsCVode), intent(inout) :: self
      real(dp), intent(in) :: t !! current time
      real(dp), intent(in) :: y(:) !! state vector
      real(dp), intent(out) :: ydot(:) !! derivative vector
      integer :: ierr !! Set to >= 0 if successful. Set to < 0 to terminate the integration 
    end function

    !> Interface for the user-supplied jacobian function
    function cvode_jac_fcn(self, t, y, jac) result(ierr)
      import :: dp, SundialsCVode
      implicit none
      class(SundialsCVode), intent(inout) :: self
      real(dp), intent(in) :: t !! current time
      real(dp), intent(in) :: y(:) !! state vector
      real(dp), intent(out) :: jac(:,:) !! Jacobian matrix
      integer  :: ierr !! Set to >= 0 if successful. Set to < 0 to terminate the integration
    end function

  end interface

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

  subroutine SundialsCVode_initialize(self, f, jac, neq, t0, y0, rtol, atol, err)
    class(SundialsCVode), target, intent(inout) :: self
    procedure(cvode_rhs_fcn) :: f !! right-hand-side function defining the system of ODEs.
    procedure(cvode_jac_fcn) :: jac !! Jacobian of the right-hand-side function.
    integer, intent(in) :: neq !! number of ODEs
    real(dp), intent(in) :: t0 !! initial time
    real(dp), intent(in) :: y0(:) !! Initial array
    real(dp), intent(in), optional :: rtol
    real(dp), intent(in), optional :: atol(:)
    character(:), allocatable, intent(out) :: err

    call SundialsCVode_initialize_impl(self, f, jac, neq, t0, y0, rtol, atol, err)
  end subroutine

  subroutine SundialsCVode_initialize_impl(self, f, jac, neq, t0, y0, rtol, atol, err)
    use, intrinsic :: iso_c_binding
    use fcvode_mod, only: CV_BDF, CV_NORMAL, CV_ONE_STEP, FCVodeInit, FCVodeSVtolerances, &
                          FCVodeSetLinearSolver, FCVode, FCVodeCreate, FCVodeFree, &
                          FCVodeSetMaxNumSteps, FCVodeSetJacFn, FCVodeSetInitStep, &
                          FCVodeGetCurrentStep, FCVodeSetMaxErrTestFails, FCVodeSetMaxOrd, &
                          FCVodeSetUserData, FCVodeSetErrHandlerFn, FCVodeSetMaxStep
    use fsundials_nvector_mod, only: N_Vector, FN_VDestroy
    use fnvector_serial_mod, only: FN_VMake_Serial   
    use fsunmatrix_band_mod, only: FSUNBandMatrix
    use fsunmatrix_dense_mod, only: FSUNDenseMatrix
    use fsundials_matrix_mod, only: SUNMatrix, FSUNMatDestroy
    use fsundials_linearsolver_mod, only: SUNLinearSolver, FSUNLinSolFree
    use fsunlinsol_band_mod, only: FSUNLinSol_Band
    use fsunlinsol_dense_mod, only: FSUNLinSol_Dense

    type(SundialsCVode), target, intent(inout) :: self
    procedure(cvode_rhs_fcn) :: f !! right-hand-side function defining the system of ODEs.
    procedure(cvode_jac_fcn) :: jac !! Jacobian of the right-hand-side function.
    integer, intent(in) :: neq !! number of ODEs
    real(dp), intent(in) :: t0 !! initial time
    real(dp), intent(in) :: y0(:) !! Initial array
    real(dp), intent(in), optional :: rtol
    real(dp), intent(in), optional :: atol(:)
    character(:), allocatable, intent(out) :: err

    integer(c_int) :: ierr
    integer(c_int64_t) :: neq_c
    real(c_double) :: t0_c, rtol_c
    type(c_ptr) :: user_data

    call self%finalize(err)
    if (allocated(err)) return

    if (neq /= size(y0)) then
      err = 'Error'
      return
    endif

    ! Set inputs
    self%f => f
    self%jac => jac
    self%neq = neq
    
    ! Convert inputs to C
    neq_c = neq
    t0_c = t0
    rtol_c = 1.0e-3_dp
    if (present(rtol)) rtol_c = rtol
    
    ! Solution vector
    allocate(self%yvec(neq_c))
    self%yvec = y0
    self%sunvec_y => FN_VMake_Serial(neq_c, self%yvec)
    if (.not. associated(self%sunvec_y)) then
      err = "CVODE setup error."
      return
    end if

    ! Absolute tolerance
    allocate(self%abstol(neq_c))
    self%abstol = 1.0e-6_dp
    if (present(atol)) self%abstol = atol
    self%abstol_nvec => FN_VMake_Serial(neq_c, self%abstol)
    if (.not. associated(self%abstol_nvec)) then
      err = "CVODE setup error."
      return
    end if

    ! Create CVode memory
    self%cvode_mem = FCVodeCreate(CV_BDF)
    if (.not. c_associated(self%cvode_mem)) then
      err = "CVODE setup error."
      return
    end if

    ! Set self pointer
    user_data = c_loc(self)
    ierr = FCVodeSetUserData(self%cvode_mem, user_data)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if

    ierr = FCVodeInit(self%cvode_mem, c_funloc(SundialsCVode_RhsFn), t0_c, self%sunvec_y)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if

    ierr = FCVodeSVtolerances(self%cvode_mem, rtol_c, self%abstol_nvec)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if

    ! If dense
    self%sunmat = FSUNDenseMatrix(neq_c, neq_c)
    self%sunlin = FSUNLinSol_Dense(self%sunvec_y, self%sunmat)
    ! If banded
    ! self%sunmat => FSUNBandMatrix(neqs_long, mu, ml)
    ! self%sunlin => FSUNLinSol_Band(self%sunvec_y, self%sunmat)

    ierr = FCVodeSetLinearSolver(self%cvode_mem, self%sunlin, self%sunmat)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if

    ierr = FCVodeSetJacFn(self%cvode_mem, c_funloc(SundialsCVode_JacFn))
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if

    !~~ Bunch of optional stuff ~~!

    ! ierr = FCVodeSetMaxNumSteps(wrk%sun%cvode_mem, mxsteps_)
    ! if (ierr /= 0) then
    !   err = "CVODE setup error."
    !   return
    ! end if
    
    ! ierr = FCVodeSetInitStep(wrk%sun%cvode_mem, var%initial_dt)
    ! if (ierr /= 0) then
    !   err = "CVODE setup error."
    !   return
    ! end if

    ! ierr = FCVodeSetMaxStep(wrk%sun%cvode_mem, var%max_dt)
    ! if (ierr /= 0) then
    !   err = "CVODE setup error."
    !   return
    ! end if

    ! ierr = FCVodeSetMaxErrTestFails(wrk%sun%cvode_mem, var%max_err_test_failures)
    ! if (ierr /= 0) then
    !   err = "CVODE setup error."
    !   return
    ! end if
    
    ! ierr = FCVodeSetMaxOrd(wrk%sun%cvode_mem, var%max_order)
    ! if (ierr /= 0) then
    !   err = "CVODE setup error."
    !   return
    ! end if

    ! ierr = FCVodeSetErrHandlerFn(wrk%sun%cvode_mem, c_funloc(ErrHandlerFn_evo), c_null_ptr)
    ! if (ierr /= 0) then
    !   err = "CVODE setup error."
    !   return
    ! end if


  end subroutine

  subroutine SundialsCVode_integrate(self, tout, tret, itask, err)
    use fcvode_mod, only: CV_ONE_STEP, FCVode, FCVodeGetNumSteps
    class(SundialsCVode), intent(inout) :: self
    real(dp), intent(in) :: tout
    real(dp), intent(out) :: tret
    integer, intent(in) :: itask
    character(:), allocatable, intent(out) :: err

    real(c_double) :: tout_c, tret_c(1)
    integer(c_int) :: itask_c, ierr_c

    tout_c = tout
    itask_c = itask

    ierr_c = FCVode(self%cvode_mem, tout_c, self%sunvec_y, tret_c, itask_c)
    if (ierr_c /= 0) then
      err = "CVODE step failed"
      return
    endif
    tret = tret_c(1)

  end subroutine

  subroutine SundialsCVode_steadystate(self, tmax, ftol, xtol, max_steps, err)
    use fcvode_mod, only: CV_ONE_STEP, FCVode, FCVodeGetCurrentTime
    class(SundialsCVode), intent(inout) :: self
    real(dp), intent(in) :: tmax
    real(dp), intent(in), optional :: ftol
    real(dp), intent(in), optional :: xtol
    integer, intent(in), optional :: max_steps
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: fvec(:)
    real(dp), allocatable :: yprev(:)
    real(dp) :: ftol_, xtol_, norm_f, norm_y, norm_dy
    real(c_double) :: tcur_c(1), tout_c, tret_c(1)
    integer :: step, max_steps_
    integer(c_int) :: ierr_c

    if (.not. associated(self%f)) then
      err = "CVODE steadystate error: RHS function not set."
      return
    end if
    if (.not. c_associated(self%cvode_mem)) then
      err = "CVODE steadystate error: solver not initialized."
      return
    end if

    ftol_ = 1.0e-8_dp
    if (present(ftol)) ftol_ = ftol
    xtol_ = 1.0e-8_dp
    if (present(xtol)) xtol_ = xtol
    max_steps_ = 500
    if (present(max_steps)) max_steps_ = max_steps

    allocate(fvec(self%neq))
    allocate(yprev(self%neq))

    ierr_c = FCVodeGetCurrentTime(self%cvode_mem, tcur_c)
    if (ierr_c /= 0) then
      err = "CVODE steadystate error: failed to get current time."
      return
    end if

    do step = 1, max_steps_
      if (tcur_c(1) >= tmax) exit

      yprev = self%yvec
      tout_c = min(real(tcur_c(1) + 1.0_dp, c_double), real(tmax, c_double))
      ierr_c = FCVode(self%cvode_mem, tout_c, self%sunvec_y, tret_c, CV_ONE_STEP)
      if (ierr_c /= 0) then
        err = "CVODE steadystate step failed."
        return
      end if
      tcur_c(1) = tret_c(1)

      ierr_c = self%f(tcur_c(1), self%yvec, fvec)
      if (ierr_c < 0) then
        err = "CVODE steadystate error: RHS function failed."
        return
      end if

      norm_f = maxval(abs(fvec))
      norm_y = maxval(abs(self%yvec))
      norm_dy = maxval(abs(self%yvec - yprev))

      if (norm_f <= ftol_) return
      if (norm_dy <= xtol_*(xtol_ + norm_y)) return
    end do

    err = "CVODE steadystate did not converge."

  end subroutine

  subroutine SundialsCVode_finalize(self, err)
    use iso_c_binding, only: c_associated, c_null_ptr, c_int
    use fcvode_mod, only: FCVodeFree
    use fsundials_nvector_mod, only: FN_VDestroy
    use fsundials_matrix_mod, only: FSUNMatDestroy
    use fsundials_linearsolver_mod, only: FSUNLinSolFree
    class(SundialsCVode), intent(inout) :: self
    character(:), allocatable, intent(out) :: err

    integer(c_int) :: ierr

    if (allocated(self%yvec)) then
      deallocate(self%yvec)
    endif
    if (associated(self%sunvec_y)) then
      call FN_VDestroy(self%sunvec_y)
      nullify(self%sunvec_y)
    endif

    if (allocated(self%abstol)) then
      deallocate(self%abstol)
    endif
    if (associated(self%abstol_nvec)) then
      call FN_VDestroy(self%abstol_nvec)
      nullify(self%abstol_nvec)
    endif

    if (c_associated(self%cvode_mem)) then
      call FCVodeFree(self%cvode_mem)
      self%cvode_mem = c_null_ptr
    endif

    if (associated(self%sunlin)) then
      ierr = FSUNLinSolFree(self%sunlin)
      if (ierr /= 0) then
        err = "Sundials failed to deallocated linear solver"
      end if
      nullify(self%sunlin)
    endif

    if (associated(self%sunmat)) then
      call FSUNMatDestroy(self%sunmat)
      nullify(self%sunmat)
    endif

  end subroutine

  subroutine SundialsCVode_final(self)
    type(SundialsCVode), intent(inout) :: self
    character(:), allocatable :: err
    call SundialsCVode_finalize(self, err)
  end subroutine

  function SundialsCVode_RhsFn(tn, sunvec_y, sunvec_f, user_data) result(ierr) bind(c, name='SundialsCVode_RhsFn')
    use, intrinsic :: iso_c_binding
    use fcvode_mod
    use fsundials_nvector_mod
    real(c_double), value :: tn
    type(N_Vector) :: sunvec_y
    type(N_Vector) :: sunvec_f
    type(c_ptr), value :: user_data
    integer(c_int) :: ierr
    
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: fvec(:)
    type(SundialsCVode), pointer :: self

    call c_f_pointer(user_data, self)
    yvec(1:self%neq) => FN_VGetArrayPointer(sunvec_y)
    fvec(1:self%neq) => FN_VGetArrayPointer(sunvec_f)
    ierr = self%f(tn, yvec, fvec)
    
  end function

  function SundialsCVode_JacFn(tn, sunvec_y, sunvec_f, sunmat_J, user_data, tmp1, tmp2, tmp3) &
                               result(ierr) bind(C,name='SundialsCVode_JacFn')
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fnvector_serial_mod
    use fsunmatrix_dense_mod
    use fsundials_matrix_mod
    real(c_double), value :: tn        ! current time
    type(N_Vector) :: sunvec_y  ! solution N_Vector
    type(N_Vector) :: sunvec_f
    type(SUNMatrix) :: sunmat_J  ! rhs N_Vector
    type(c_ptr), value :: user_data ! user-defined data
    type(N_Vector) :: tmp1, tmp2, tmp3
    integer(c_int) :: ierr
  
    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: jac(:,:)
    real(c_double), pointer :: jac_data(:)
    integer(c_long) :: nrows, ncols
    type(SundialsCVode), pointer :: self
    
    call c_f_pointer(user_data, self)
    yvec(1:self%neq) => FN_VGetArrayPointer(sunvec_y)
    jac_data => FSUNDenseMatrix_Data(sunmat_J)
    nrows = FSUNDenseMatrix_Rows(sunmat_J)
    ncols = FSUNDenseMatrix_Columns(sunmat_J)
    call c_f_pointer(c_loc(jac_data(1)), jac, [nrows, ncols])
    ierr = self%jac(tn, yvec, jac)
  end function
  
end module
