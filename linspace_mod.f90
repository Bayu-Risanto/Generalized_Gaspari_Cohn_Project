module linspace_mod

  implicit none

  integer, parameter :: r8 = selected_real_kind(15, 307)
  real(r8), parameter :: pi = 4.0_r8 * atan(1.0_r8)

  private
  public :: linspace, r8, pi

contains

  ! Builds a linearly spaced array of length num_close between a and b.
  ! num_close is intended to be the scalar returned by get_close_obs()
  ! in assim_tools_mod.f90, so the spacing (increment) adapts on every
  ! call to however many close locations were found for the current
  ! observation:
  !     step = (b - a) / (num_close - 1)   [endpoint = .true.,  default]
  !     step = (b - a) /  num_close        [endpoint = .false.]
  !
  ! x is allocated here (intent(out), allocatable) rather than
  ! pre-sized by the caller, since num_close is only known at call
  ! time (e.g. immediately after get_close_obs() returns).

  subroutine linspace(a, b, num_close, x, endpoint)
    real(r8), intent(in)  :: a, b
    integer,  intent(in)  :: num_close
    real(r8), allocatable, intent(out) :: x(:)
    logical,  intent(in), optional :: endpoint

    real(r8) :: step
    integer  :: i
    logical  :: include_endpoint

    include_endpoint = .true.
    if (present(endpoint)) include_endpoint = endpoint

    if (num_close <= 0) then
      ! no close locations found for this observation -- nothing to build
      allocate(x(0))
      return
    end if

    allocate(x(num_close))

    if (num_close == 1) then
      x(1) = a
      return
    end if

    if (include_endpoint) then
      step = (b - a) / real(num_close - 1, r8)
    else
      step = (b - a) / real(num_close, r8)
    end if

    x = [(a + step * real(i - 1, r8), i = 1, num_close)]

  end subroutine linspace

end module linspace_mod
