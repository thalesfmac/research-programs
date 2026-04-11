module lead_green_function
    use :: precision, only : dp
    use :: constants, only : CI
    use :: lapack_blas, only : matmul3
    implicit none

    private
    public :: surface_gf_1d
    public :: surface_self_energy_left, surface_self_energy_right
    public :: broadening

    contains

    function f(x) result(res)
        real(dp), intent(in) :: x
        complex(dp) :: res
        ! real(dp), parameter :: ETA = 1.0e-12_dp
        ! real(dp), parameter :: ETA = 0.0_dp

        if (abs(x) >= 2.0_dp) error stop "f: |x| must be < 2"
        ! res = cmplx(x, kind=dp) - CI * sqrt(4.0_dp - x*x)
        res = cmplx(x, -sqrt(4.0_dp - x*x), kind=dp)
    end function f

    function surface_gf_1d(E, tlead, mu) result(gs)
        real(dp), intent(in) :: E, tlead, mu
        complex(dp) :: gs

        real(dp) :: x

        x = (E - mu) / tlead
        gs = f(x) / cmplx(2.0_dp*tlead, kind=dp)
    end function surface_gf_1d

    subroutine surface_self_energy_left(surf_gf_l, u_left, sigma_left)
        complex(dp), intent(in)  :: surf_gf_l(:, :), u_left(:, :)
        complex(dp), intent(out) :: sigma_left(:, :)

        call matmul3(u_left, surf_gf_l, u_left, sigma_left, transc="C")
    end subroutine surface_self_energy_left

    subroutine surface_self_energy_right(surf_gf_r, u_right, sigma_right)
        complex(dp), intent(in)  :: surf_gf_r(:, :), u_right(:, :)
        complex(dp), intent(out) :: sigma_right(:, :)

        call matmul3(u_right, surf_gf_r, u_right, sigma_right, transa="C")
    end subroutine surface_self_energy_right

    subroutine broadening(sigma, gam)
        complex(dp), intent(in) :: sigma(:, :)
        complex(dp), intent(out) :: gam(:, :)

        ! if (size(sigma, 1) /= size(sigma, 2)) error stop "broadening: sigma must be square"
        if (size(gam, 1) /= size(sigma, 1) .or. size(gam, 2) /= size(sigma, 2)) then
            error stop "broadening: gam has incompatible dimensions with sigma"
        end if

        gam = CI * (sigma - conjg( transpose( sigma ) ))
    end subroutine broadening

end module lead_green_function
