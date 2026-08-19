module lead_green_function
   use stdlib_kinds, only: dp
   use matrix_operations, only: matmul3
   implicit none
   private

   public :: surface_gf_1d
   public :: surface_self_energy_left, surface_self_energy_right
   public :: broadening

contains

   function f(x) result(res)
      real(dp), intent(in) :: x
      complex(dp) :: arg, res

      arg = cmplx(x*x - 4.0_dp, kind=dp)
      res = cmplx(x, kind=dp) - sqrt(arg)
   end function f

   function surface_gf_1d(E, tlead, mu, eta) result(gs)
      use stdlib_optval, only: optval
      real(dp), intent(in) :: E, tlead, mu
      real(dp), intent(in), optional :: eta
      complex(dp) :: gs

      real(dp) :: x

      x = (E - mu)/tlead
      gs = f(x)/cmplx(2.0_dp*tlead, kind=dp)

      if (abs(E - mu) >= 2.0_dp*tlead) then
         gs = gs - cmplx(0.0_dp, optval(eta, 0.0_dp), kind=dp)
      end if
   end function surface_gf_1d

   subroutine surface_self_energy_left(surf_gf_l, u_left, sigma_left)
      complex(dp), intent(in), contiguous :: surf_gf_l(:, :), u_left(:, :)
      complex(dp), intent(out), contiguous :: sigma_left(:, :)

      call matmul3(u_left, surf_gf_l, u_left, sigma_left, transc="C")
   end subroutine surface_self_energy_left

   subroutine surface_self_energy_right(surf_gf_r, u_right, sigma_right)
      complex(dp), intent(in), contiguous :: surf_gf_r(:, :), u_right(:, :)
      complex(dp), intent(out), contiguous :: sigma_right(:, :)

      call matmul3(u_right, surf_gf_r, u_right, sigma_right, transa="C")
   end subroutine surface_self_energy_right

   subroutine broadening(sigma, gam)
      use stdlib_linalg, only: hermitian
      complex(dp), intent(in), contiguous :: sigma(:, :)
      complex(dp), intent(out), contiguous :: gam(:, :)

      complex(dp), parameter :: CI = (0.0_dp, 1.0_dp)

      ! if (size(sigma, 1) /= size(sigma, 2)) error stop "broadening: sigma must be square"
      if (size(gam, 1) /= size(sigma, 1) .or. size(gam, 2) /= size(sigma, 2)) then
         error stop "broadening: gam has incompatible dimensions with sigma"
      end if

      gam = CI*(sigma - hermitian(sigma))
   end subroutine broadening

end module lead_green_function
