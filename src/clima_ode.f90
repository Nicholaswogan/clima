module clima_ode
  use clima_const, only: dp
  implicit none
  private
  public :: ODEStepper

  !> A simple ODE stepper class using Heun's method (https://en.wikipedia.org/wiki/Heun%27s_method).
  !> Adaptive stepsize selection follows pages 167 - 169 in 
  !> "Solving Ordinary Differential Equations I" by Hairer, Norsett and Wanner
  type :: ODEStepper

    procedure(rhs_fcn), pointer :: fcn !! RHS function
    integer :: n !! number of equations

    ! Current state of ODE integrator
    real(dp) :: t !! current time
    real(dp), allocatable :: u(:) !! current value of variables
    real(dp) :: h !! current stepsize

    ! settings
    real(dp) :: rtol = 1.0e-4_dp !! relative tolerance
    real(dp) :: atol = 1.0e-6_dp !! absolute tolerance 
    real(dp) :: facmax = 3.0_dp
    real(dp) :: facmin = 0.3_dp
    real(dp) :: safety_factor = 0.7_dp
    integer :: max_error_test_failures = 10

    ! statistics for solver
    integer :: nfevals !! number of function evalutions
    integer :: naccept !! number of accepted steps
    integer :: nreject !! number of steps that were rejected

    ! work space
    real(dp), allocatable :: du(:)
    real(dp), allocatable :: du1(:)
    real(dp), allocatable :: u1(:)
    real(dp), allocatable :: u_new(:)
    real(dp) :: t_new
    real(dp), allocatable :: sc(:)
    real(dp) :: h_new

  contains
    procedure :: initialize => ODEStepper_initialize
    procedure :: initial_stepsize => ODEStepper_initial_stepsize
    procedure :: step => ODEStepper_step
  end type

  abstract interface
    subroutine rhs_fcn(self, t, u, du, ierr)
      import :: dp, ODEStepper
      class(ODEStepper), intent(inout) :: self
      real(dp), intent(in) :: t
      real(dp), intent(in) :: u(:)
      real(dp), intent(out) :: du(:)
      integer, intent(out) :: ierr
    end subroutine
  end interface

contains

  subroutine ODEStepper_initialize(self, fcn, t0, u0, ierr)
    class(ODEStepper), intent(inout) :: self
    procedure(rhs_fcn) :: fcn
    real(dp), intent(in) :: t0
    real(dp), intent(in) :: u0(:)
    integer, intent(out) :: ierr

    self%fcn => fcn
    self%n = size(u0)

    self%t = t0
    if (allocated(self%u)) deallocate(self%u)
    self%u = u0

    self%nfevals = 0
    self%naccept = 0
    self%nreject = 0

    if (allocated(self%du)) then
      deallocate(self%du)
      deallocate(self%du1)
      deallocate(self%u1)
      deallocate(self%u_new)
      deallocate(self%sc)
    endif
    allocate(self%du(self%n))
    allocate(self%du1(self%n))
    allocate(self%u1(self%n))
    allocate(self%u_new(self%n))
    allocate(self%sc(self%n))

    self%h = self%initial_stepsize(t0, u0, ierr)
    if (ierr < 0) return

  end subroutine

  subroutine ODEStepper_step(self, ierr)
    class(ODEStepper), intent(inout) :: self
    integer, intent(out) :: ierr

    logical :: accept_step
    real(dp) :: facmax, err
    integer :: ierr1
    integer :: i, j

    facmax = self%facmax

    ierr = 0
    do j = 1,self%max_error_test_failures

      ! Try to perform a step of Heun's Method
      self%t_new = self%t + self%h
      call self%fcn(self%t, self%u, self%du, ierr1)
      if (ierr1 < 0) then
        ierr = -2
        return
      endif
      self%nfevals = self%nfevals + 1
      self%u1 = self%u + self%h*self%du ! first order approximation to u_new
      call self%fcn(self%t_new, self%u1, self%du1, ierr1)
      if (ierr1 < 0) then
        ierr = -2
        return
      endif
      self%nfevals = self%nfevals + 1
      self%u_new = self%u + 0.5_dp*self%h*(self%du + self%du1) ! second-order approximation to u_new

      ! Error and timestep control. Methods generally follow Page 167 - 168 in 
      ! "Solving Ordinary Differential Equations I" by Hairer, Norsett and Wanner
      do i = 1,self%n
        self%sc(i) = self%atol + self%rtol*max(abs(self%u(i)), abs(self%u_new(i)))
      enddo

      if (all(abs(self%u_new - self%u1) <= self%sc)) then
        accept_step = .true.
        facmax = self%facmax
      else
        accept_step = .false.
        facmax = 1.0_dp
      endif
      
      ! Compute new stepsize
      err = sqrt((1.0_dp/real(self%n,dp))*sum(((self%u_new - self%u1)/self%sc)**2.0_dp))
      self%h_new = self%safety_factor*self%h*min(facmax, max(self%facmin, (1.0_dp/err)**(0.5_dp)))
      self%h = self%h_new
      
      if (accept_step) then
        self%t = self%t_new
        self%u = self%u_new
        self%naccept = self%naccept + 1
        return
      else
        self%nreject = self%nreject + 1
      endif

    enddo

    ! Will only get here if step failed for `self%max_error_test_failures`
    ierr = -1
    return

  end subroutine

  !> Computes the initial stepsize. This follows 169 in 
  !> "Solving Ordinary Differential Equations I" by Hairer, Norsett and Wanner
  function ODEStepper_initial_stepsize(self, t0, u0, ierr) result(h_init)
    class(ODEStepper), intent(inout) :: self
    real(dp), intent(in) :: t0
    real(dp), intent(in) :: u0(:)
    integer, intent(out) :: ierr
    real(dp) :: h_init

    integer :: i, ierr1
    real(dp) :: d0, d1, d2
    real(dp) :: h0, h1

    ierr = 0

    !~~ Step a) ~~!
    call self%fcn(t0, u0, self%du, ierr1)
    if (ierr1 < 0) then
      ierr = -2
      return
    endif

    do i = 1,self%n
      self%sc(i) = self%atol + self%rtol*abs(u0(i))
    enddo
    d0 = sqrt((1.0_dp/real(self%n,dp))*sum((u0/self%sc)**2.0_dp))
    d1 = sqrt((1.0_dp/real(self%n,dp))*sum((self%du/self%sc)**2.0_dp))

    !~~ Step b) ~~!
    h0 = 0.01_dp*(d0/d1)
    if (d0 < 1.0e-5_dp .or. d1 < 1.0e-5_dp) then
      h0 = 1.0e-6_dp
    endif

    !~~ Step c) ~~!
    self%u1 = u0 + h0*self%du
    call self%fcn(t0+h0, self%u1, self%du1, ierr1)
    if (ierr1 < 0) then
      ierr = -2
      return
    endif

    !~~ Step d) ~~!
    d2 = sqrt((1.0_dp/real(self%n,dp))*sum(((self%du1 - self%du)/self%sc)**2.0_dp))/h0

    !~~ Step e) ~~!
    h1 = (0.01_dp/max(d1,d2))**(0.5_dp)
    if (max(d1,d2) < 1.0e-15_dp) then
      h1 = max(1.0e-6_dp, h0*1.0e-3_dp)
    endif

    !~~ Step f) ~~!
    h_init = min(100.0_dp*h0, h1)

  end function

end module