submodule(clima_climate) clima_climate_integrate
  use dop853_module, only: dop853_class
  implicit none
  
  type, extends(dop853_class) :: dop853_custom
    type(Climate), pointer :: c

    integer :: j
    real(dp), pointer :: t_eval(:)
    character(:), pointer :: filename
    real(dp), allocatable :: u(:), du(:)
    character(:), allocatable :: err

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

  subroutine solout(self, nr, xold, x, y, irtrn, xout)
    class(dop853_class),intent(inout) :: self
    integer,intent(in)                :: nr
    real(dp),intent(in)               :: xold
    real(dp),intent(in)               :: x
    real(dp),dimension(:),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(dp),intent(out)              :: xout

    type(dop853_custom), pointer :: dop

    real(dp) :: time_old, time_cur
    integer :: i

    time_old = xold
    time_cur = x

    select type (self)
    class is (dop853_custom)

    if (nr == 1) then
      xout = self%t_eval(1)
      self%j = 1
    else

      do while (self%t_eval(self%j) <= x .and. self%j <= size(self%t_eval))
        do i = 1,size(y)
          self%u(i) = self%contd8(i, self%t_eval(self%j))
        enddo

        call self%c%right_hand_side(self%u, self%du, self%err)

        open(1, file = self%filename, status='old', form="unformatted",position="append")
        write(1) self%t_eval(self%j)
        write(1) self%u

        write(1) self%c%rad%f_total
        write(1) self%c%rad%wrk_ir%fup_n
        write(1) self%c%rad%wrk_ir%fdn_n
        write(1) self%c%rad%wrk_sol%fup_n
        write(1) self%c%rad%wrk_sol%fdn_n
        write(1) [self%c%surface_pressure, self%c%wrk%P]
        close(1)

        self%j = self%j + 1
      enddo
      if (self%j <= size(self%t_eval)) then
        xout = self%t_eval(self%j)
      endif

    endif

    end select

  end subroutine
  
  module function evolve(self, filename, tstart, T_start, t_eval, overwrite, err) result(success)
                                   
    ! in/out
    class(Climate), target, intent(inout) :: self
    character(*), target, intent(in) :: filename
    real(dp), intent(in) :: tstart
    real(dp), intent(in) :: T_start(:)
    real(dp), target, intent(in) :: t_eval(:)
    logical, intent(in) :: overwrite
    character(:), allocatable, intent(out) :: err

    logical :: success
    
    integer :: i, j, ii, io
    
    integer :: idid
    logical :: status_ok
    real(dp) :: tn, tout, hcur
    real(dp) :: yvec(self%neq)
    real(dp) :: t_eval_1(size(t_eval)+1)
    integer :: icomp(self%neq)
    type(dop853_custom) :: dop
    
    ! check dimensions
    if (size(T_start) /= self%neq) then
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
    write(1) [0.0_dp, self%z]
    write(1) size(t_eval)
    close(1)

    do i = 1,self%neq
      icomp(i) = i
    enddo
    
    ! initialize
    call dop%initialize(fcn=rhs_dop853, solout=solout, n=self%neq, icomp=icomp, status_ok=status_ok)
    if (.not. status_ok) then
      err = "failed to initialize dop853"
      return
    endif
    
    dop%c => self
    dop%t_eval => t_eval
    dop%filename => filename
    allocate(dop%u(self%neq), dop%du(self%neq))
    
    yvec = T_start
    tn = tstart
    tout = t_eval(size(t_eval))
    call dop%integrate(tn, yvec, tout, [self%rtol], [self%atol], iout=3, idid=idid)
    if (idid /= 1) then
      success = .false.
    else
      success = .true.
    endif
    
  end function

end submodule