! This linspace_mod.f90 is to create the prescribed distance
! between the point of observations to the neighboring points
! during the assimilation using Generelized GC. 
! The unit is in radians
! Created by C. Bayu Risanto, S.J. (24 August 2026)

module linspace_mod

  implicit none

  integer, parameter :: r8 = selected_real_kind(15, 307)
  real(r8), parameter :: pi = 4.0_r8 * atan(1.0_r8)

contains

  subroutine linspace(a, b, x, endpoint)
    real(r8), intent(in)  :: a, b
    real(r8), intent(out) :: x(:)
    logical,  intent(in), optional :: endpoint

    real(r8) :: step
    integer  :: i, n
    logical  :: include_endpoint

    n = size(x)

    include_endpoint = .true.
    if (present(endpoint)) include_endpoint = endpoint

    if (n == 1) then
      x(1) = a
      return
    end if

    if (include_endpoint) then
      step = (b - a) / real(n - 1, r8)
    else
      step = (b - a) / real(n, r8)
    end if

    do i = 1, n
      x(i) = a + step * real(i - 1, r8)
    end do

  end subroutine linspace

end module linspace_mod
