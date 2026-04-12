module transmittance
    use :: precision, only : dp
    use :: matrix_operations, only : trace, invert, matmul3, matmul4
    implicit none

    private
    public :: caroli_transmission
    public :: rgf_first_step, rgf_step, rgf_last_step

    contains

    function caroli_transmission(gf, Gamma_L, Gamma_R) result(T)
        complex(dp), dimension(:,:), intent(in), contiguous :: gf, Gamma_L, Gamma_R
        real(dp) :: T

        complex(dp), allocatable :: tmp(:,:)

        allocate( tmp(size(Gamma_L, 1), size(gf, 1)) )
        call matmul4(Gamma_L, gf, Gamma_R, gf, tmp, transd="C")
        T = real(trace(tmp), kind=dp)
    end function caroli_transmission

    subroutine rgf_first_step(cE, h_1, U_01, g_L, G_11, G_01)
        complex(dp), dimension(:,:), intent(in), contiguous :: cE, h_1, U_01, g_L
        complex(dp), dimension(:,:), intent(out), contiguous :: G_11, G_01

        if (size(cE, 1) /= size(h_1, 1) .or. size(cE, 2) /= size(h_1, 2)) then
            error stop "rgf_first_step: cE has incompatible dimensions with h_i"
        end if

        ! G_11 will be used to store the matrix product
        call matmul3(U_01, g_L, U_01, G_11, transa="C")

        G_11 = cE - h_1 - G_11
        call invert(G_11)

        call matmul3(g_L, U_01, G_11, G_01)
    end subroutine rgf_first_step

    subroutine rgf_step(cE, h_n, U_nm1_n, G_nm1_nm1, G_0_nm1, G_n_n, G_0_n)
        complex(dp), dimension(:,:), intent(in), contiguous :: cE, h_n, U_nm1_n, G_nm1_nm1, G_0_nm1
        complex(dp), dimension(:,:), intent(out), contiguous :: G_n_n, G_0_n

        if (size(cE, 1) /= size(h_n, 1) .or. size(cE, 2) /= size(h_n, 2)) then
            error stop "rgf_step: cE has incompatible dimensions with h_i"
        end if

        ! G_n_n will be used to store the matrix product
        call matmul3(U_nm1_n, G_nm1_nm1, U_nm1_n, G_n_n, transa="C")
        G_n_n = cE - h_n - G_n_n
        call invert(G_n_n)

        call matmul3(G_0_nm1, U_nm1_n, G_n_n, G_0_n)
    end subroutine rgf_step

    subroutine rgf_last_step(g_R, U_N_Np1, G_N_N, G_0_N, G_Np1_Np1, G_0_Np1)
        complex(dp), dimension(:,:), intent(in), contiguous :: g_R, G_N_N, G_0_N, U_N_Np1
        complex(dp), dimension(:,:), intent(out), contiguous :: G_Np1_Np1, G_0_Np1

        complex(dp), dimension(:,:), allocatable :: g_R_inv

        allocate( g_R_inv(size(g_R, 1), size(g_R, 2)) )

        g_R_inv = g_R
        call invert(g_R_inv)

        ! G_Np1_Np1 will be used to store the matrix product
        call matmul3(U_N_Np1, G_N_N, U_N_Np1, G_Np1_Np1, transa="C")
        G_Np1_Np1 = g_R_inv - G_Np1_Np1
        call invert(G_Np1_Np1)

        call matmul3(G_0_N, U_N_Np1, g_R_inv, G_0_Np1)
    end subroutine rgf_last_step

end module transmittance
