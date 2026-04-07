module lead_green_function
    use :: precision, only : dp
    use :: constants, only : CI
    use :: matrix_operations
    implicit none

    private
    public :: surface_gf_1d
    public :: surface_self_energy_left, surface_self_energy_right

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

        complex(dp), allocatable :: u_left_dagger(:, :)

        u_left_dagger = conjg( transpose( u_left ) )

        call assert_square(surf_gf_l, "surf_gf_l")
        ! call assert_same_shape(surf_gf_l, sigma_left, matmul( G_NN, U_NNp1 ) )
        ! call invert(z_2)
        ! G_0Np1 = matmul( G_0N, matmul( U_NNp1, z_2 ) ) "surf_gf_l", "sigma_left")

        sigma_left = matmul3(u_left, surf_gf_l, u_left_dagger)

    end subroutine surface_self_energy_left

    subroutine surface_self_energy_right(surf_gf_r, u_right, sigma_right)
        complex(dp), intent(in)  :: surf_gf_r(:, :), u_right(:, :)
        complex(dp), intent(out) :: sigma_right(:, :)

        complex(dp), allocatable :: u_right_dagger(:, :)

        u_right_dagger = conjg( transpose( u_right ) )

        call assert_square(surf_gf_r, "surf_gf_r")
        call assert_same_shape(surf_gf_r, sigma_right, "surf_gf_r", "sigma_right")

        sigma_right = matmul3(u_right_dagger, surf_gf_r, u_right)

    end subroutine surface_self_energy_right

    subroutine broadening(sigma, gam)
        complex(dp), intent(in) :: sigma(:, :)
        complex(dp), intent(out) :: gam(:, :)

        call assert_same_shape(sigma, gam, "sigma", "gam")

        gam = CI * (sigma - conjg( transpose( sigma ) ))
    end subroutine broadening

end module lead_green_function
