submodule(clima_climate) clima_climate_integrate
  use dop853_module, only: dop853_class
  implicit none
  
  type, extends(dop853_class) :: dop853_custom
    type(Climate), pointer :: c
  end type
    
contains
  
  subroutine rhs_dop853(self, tn, u, du)
    class(dop853_class), intent(inout) :: self
    real(dp), intent(in) :: tn
    real(dp), intent(in), dimension(:) :: u
    real(dp), intent(out), dimension(:) :: du
    
    character(:), allocatable :: err
    integer :: nsteps
    real(dp) :: hcur
    type(Climate), pointer :: c
    
    select type (self)
    class is (dop853_custom)
      c => self%c
    end select
    
    call c%right_hand_side(u, du, err)
    
    call self%info(h=hcur)
    call self%info(nstep=nsteps)
    
    if (nsteps /= c%wrk%nsteps_previous) then
      print"(1x,'N =',i6,3x,'Time = ',es11.5,3x,'dt = ',es11.5,3x,'max(dy/dt) = ',es11.5)", &
            nsteps, tn, hcur, maxval(abs(du))
      ! print*,'FTIR ',c%w%fdn_ir(c%v%nz+1)-c%w%fup_ir(c%v%nz+1)
      ! print*,'FTSO ',c%w%fdn_sol(c%v%nz+1)-c%w%fup_sol(c%v%nz+1)
      ! print*,'SEFF',abs((c%w%fdn_ir(c%v%nz+1)-c%w%fup_ir(c%v%nz+1))/(c%w%fdn_sol(c%v%nz+1)-c%w%fup_sol(c%v%nz+1)))
      
      c%wrk%nsteps_previous = nsteps
    endif

  end subroutine
  
  module function evolve(self, filename, tstart, T_start, t_eval, overwrite, err) result(success)
                                   
    ! in/out
    class(Climate), target, intent(inout) :: self
    character(*), intent(in) :: filename
    real(dp), intent(in) :: tstart
    real(dp), intent(in) :: T_start(:)
    real(dp), intent(in) :: t_eval(:)
    logical, intent(in) :: overwrite
    character(:), allocatable, intent(out) :: err

    logical :: success
    
    integer :: i, j, ii, io
    
    integer :: idid
    logical :: status_ok
    real(dp) :: tn, tout, hcur
    real(dp) :: yvec(self%nz)
    real(dp) :: t_eval_1(size(t_eval)+1)
    type(dop853_custom) :: dop
    
    ! check dimensions
    if (size(T_start) /= self%nz) then
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
    write(1) self%nz
    write(1) self%z
    write(1) size(t_eval)
    close(1)
    
    ! initialize
    call dop%initialize(fcn=rhs_dop853, n=self%nz, status_ok=status_ok)
    if (.not. status_ok) then
      err = "failed to initialize dop853"
      return
    endif
    
    dop%c => self
    
    ! first integrate from tstart to t_teval(1)
    yvec = T_start
    t_eval_1(1) = tstart
    t_eval_1(2:) = t_eval

    do ii = 2, size(t_eval_1)
      
      tn = t_eval_1(ii-1)
      tout = t_eval_1(ii)
      call dop%integrate(tn, yvec, tout, [1.0e-3_dp], [1.0e-6_dp], iout=0, idid=idid)
      call dop%info(h=hcur)
      call dop%initialize(fcn=rhs_dop853, n=self%nz, hinitial=hcur, status_ok=status_ok)
      dop%c => self
      
      if (idid /= 1) then
        print*,idid
        success = .false.
        exit
      else
        success = .true.
    
        open(1, file = filename, status='old', form="unformatted",position="append")
        write(1) tn
        write(1) yvec
        write(1) self%rad%f_total
        write(1) self%rad%wrk_ir%fup_n
        write(1) self%rad%wrk_ir%fdn_n
        write(1) self%rad%wrk_sol%fup_n
        write(1) self%rad%wrk_sol%fdn_n
        write(1) self%wrk%P
        close(1)
    
      endif
    enddo
    
  end function
  

end submodule