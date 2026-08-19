module density_of_states
   use stdlib_kinds, only: dp
   implicit none
   private

   public dos

contains

   subroutine green_function(h, E, eta, gf)
      use stdlib_linalg, only: diag, invert
      complex(dp), intent(in) :: h(:, :)
      real(dp), intent(in) :: E, eta
      complex(dp), allocatable, intent(out) :: gf(:, :)

      complex(dp) :: z
      integer :: n, i

      n = size(h, dim=1)
      z = cmplx(E, eta, kind=dp)

      gf = diag([(z, i=1, n)])
      gf = gf - h

      call invert(gf)
   end subroutine green_function

   function dos(h, E, eta) result(rho)
      use stdlib_linalg, only: trace
      use stdlib_constants, only: PI => PI_dp
      complex(dp), intent(in) :: h(:, :)
      real(dp), intent(in) :: E, eta
      real(dp) :: rho

      complex(dp), allocatable :: gf(:, :)
      integer :: n

      n = size(h, dim=1)

      call green_function(h, E, eta, gf)

      rho = -aimag(trace(gf))/(PI*real(n, kind=dp))
   end function dos

end module density_of_states
