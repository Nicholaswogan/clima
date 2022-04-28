submodule(clima_climate) clima_climate_rhs
  implicit none
  
  
contains
  
  module function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
                        result(ierr) bind(c, name='RhsFn')
    use, intrinsic :: iso_c_binding
    use fcvode_mod
    use fsundials_nvector_mod
    ! calling variables
    real(c_double), value :: tn        ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! rhs N_Vector
    type(c_ptr), value    :: user_data ! user-defined dat
    integer(c_int)        :: ierr
    
    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: fvec(:)
    character(:), allocatable :: err
    integer(c_long) :: nsteps(1)
    integer(c_int) :: loc_ierr
    real(c_double) :: hcur(1)
    real(dp) :: tmp, mx
    integer :: k, i, j, ii
    
    type(Climate), pointer :: self
  
    ierr = 0
    
    call c_f_pointer(user_data, self)
    
    ! get data arrays from SUNDIALS vectors
    yvec(1:self%v%nz) => FN_VGetArrayPointer(sunvec_y)
    fvec(1:self%v%nz) => FN_VGetArrayPointer(sunvec_f)
    
    ! fill RHS vector
    call self%right_hand_side(yvec, fvec, err)
    loc_ierr = FCVodeGetNumSteps(self%w%cvode_mem, nsteps)
    
    if (nsteps(1) /= self%w%nsteps_previous) then
      loc_ierr = FCVodeGetCurrentStep(self%w%cvode_mem, hcur)
    
      print"(1x,'N =',i6,3x,'Time = ',es11.5,3x,'dt = ',es11.5,3x,'max(dy/dt) = ',es11.5)", &
         nsteps, tn, hcur(1), maxval(abs(fvec))
         
      
      self%w%nsteps_previous = nsteps(1)
    endif
    
    if (allocated(err)) then
      print*,trim(err)//". CVODE will attempt to correct the error."
      ierr = 1
    endif

  end function
  
  module function evolve(self, filename, tstart, T_start, t_eval, overwrite, err) result(success)
                                   
    use, intrinsic :: iso_c_binding
    use fcvode_mod, only: CV_ADAMS, CV_NORMAL, CV_ONE_STEP, FCVodeInit, FCVodeSStolerances, &
                          FCVodeSetNonLinearSolver, FCVode, FCVodeCreate, FCVodeFree, &
                          FCVodeSetMaxNumSteps, FCVodeSetJacFn, FCVodeSetInitStep, &
                          FCVodeGetCurrentStep, FCVodeSetMaxErrTestFails, FCVodeSetMaxOrd, &
                          FCVodeSetUserData
    use fsundials_nvector_mod, only: N_Vector, FN_VDestroy
    use fnvector_serial_mod, only: FN_VMake_Serial   
    
    use fsundials_nonlinearsolver_mod ! Fortran interface to generic SUNNonlinearSolver
    use fsunnonlinsol_fixedpoint_mod 
    
    
    ! in/out
    class(Climate), target, intent(inout) :: self
    character(len=*), intent(in) :: filename
    real(c_double), intent(in) :: tstart
    real(dp), intent(in) :: T_start(:)
    real(c_double), intent(in) :: t_eval(:)
    logical, optional, intent(in) :: overwrite
    logical :: success
    character(:), allocatable, intent(out) :: err
    
    ! local variables
    real(c_double) :: tcur(1)    ! current time
    integer(c_int) :: ierr       ! error flag from C functions
    type(N_Vector), pointer :: sunvec_y ! sundials vector
    
    ! solution vector
    real(c_double) :: yvec(self%v%nz)
    integer(c_int64_t) :: neqs_long
    integer(c_long) :: mxsteps_
    type(SUNNonlinearSolver), pointer :: sunnls 
    
    integer :: i, j, k, ii, io
    
    type(c_ptr)    :: user_data
    type(ClimaData), pointer :: d
    type(ClimaVars), pointer :: v
    type(ClimaWrk), pointer :: w
    type(Climate), pointer :: self_ptr
  
    d => self%d
    v => self%v
    w => self%w
    
    ! check dimensions
    if (size(T_start) /= v%nz) then
      err = "Input to evolve has the wrong dimension"
      return
    endif
    
    ! file prep
    if (overwrite) then
      open(1, file = filename, status='replace', form="unformatted")
    else
      open(1, file = filename, status='new', form="unformatted",iostat=io)
      if (io /= 0) then
        err = "Unable to create file "//trim(filename)//" because it already exists"
        return
      endif
    endif
    write(1) v%nz
    write(1) v%z
    write(1) size(t_eval)
    close(1)
    
    ! settings
    mxsteps_ = 10000
    neqs_long = v%nz
    tcur   = tstart
    
    self_ptr => self
    user_data = c_loc(self_ptr)
    
    ! initialize solution vector
    yvec = T_start

    ! create SUNDIALS N_Vector
    sunvec_y => FN_VMake_Serial(neqs_long, yvec)
    if (.not. associated(sunvec_y)) then
      err = "CVODE setup error."
      return
    end if
    
    ! create CVode memory
    w%cvode_mem = FCVodeCreate(CV_ADAMS)
    if (.not. c_associated(w%cvode_mem)) then
      err = "CVODE setup error."
      return
    end if
    
    ! set user data
    ierr = FCVodeSetUserData(w%cvode_mem, user_data)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeInit(w%cvode_mem, c_funloc(RhsFn), tstart, sunvec_y)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSStolerances(w%cvode_mem, 1.0e-3_dp, 1.0e-6_dp)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ! create fixed point nonlinear solver object
    sunnls => FSUNNonlinSol_FixedPoint(sunvec_y, 0)
    if (.not. associated(sunnls)) then
       print *,'ERROR: sunnls = NULL'
       stop 1
    end if

    ! attache nonlinear solver object to CVode
    ierr = FCVodeSetNonlinearSolver(w%cvode_mem, sunnls)
    if (ierr /= 0) then
      print *, 'Error in FCVodeSetNonlinearSolver, ierr = ', ierr, '; halting'
      stop 1
    end if
    
    ! sunmat => FSUNDenseMatrix(neqs_long, neqs_long)
    ! sunlin => FSUNLinSol_Dense(sunvec_y, sunmat)
    ! 
    ! ierr = FCVodeSetLinearSolver(w%cvode_mem, sunlin, sunmat)
    ! if (ierr /= 0) then
    !   err = "CVODE setup error."
    !   return
    ! end if
    ! 
    ! ierr = FCVodeSetJacFn(wrk%cvode_mem, c_funloc(JacFn))
    ! if (ierr /= 0) then
    !   err = "CVODE setup error."
    !   return
    ! end if
    ! 
    ierr = FCVodeSetMaxNumSteps(w%cvode_mem, mxsteps_)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetInitStep(w%cvode_mem, 1.0e-6_dp)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ierr = FCVodeSetMaxErrTestFails(w%cvode_mem, 15)
    if (ierr /= 0) then
      err = "CVODE setup error."
      return
    end if
    
    ! ierr = FCVodeSetMaxOrd(wrk%cvode_mem, var%max_order)
    ! if (ierr /= 0) then
    !   err = "CVODE setup error."
    !   return
    ! end if
    
    do ii = 1, size(t_eval)
      ierr = FCVode(w%cvode_mem, t_eval(ii), sunvec_y, tcur, CV_NORMAL)
      if (ierr /= 0) then
        success = .false.
        exit
      else
        success = .true.

        open(1, file = filename, status='old', form="unformatted",position="append")
        write(1) tcur(1)
        write(1) yvec
        close(1)
    
      endif
    enddo
    
    ! free memory
    call FN_VDestroy(sunvec_y)
    call FCVodeFree(w%cvode_mem)
    ierr = FSUNNonLinSolFree(sunnls)
    if (ierr /= 0) then
      err = "CVODE deallocation error"
      return
    end if

  end function
  
  
  
end submodule