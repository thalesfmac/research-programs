module transmittance
    use :: precision, only : dp
    use :: lapack_blas, only : invert, matmul2
    use :: matrix_operations, only : trace
    implicit none

    private
    public :: caroli_transmission
    public :: rgf_step, rgf_last_step

    contains

    function caroli_transmission(gf, Gamma_L, Gamma_R) result(T)
        complex(dp), intent(in) :: gf(:,:), Gamma_L(:,:), Gamma_R(:,:)
        real(dp) :: T

        complex(dp), allocatable :: tmp1(:,:), tmp2(:,:), tmp3(:,:)

        call matmul2(Gamma_L, gf, tmp1)
        call matmul2(Gamma_R, gf, tmp2, transb="C")
        call matmul2(tmp1, tmp2, tmp3)

        T = real(trace(tmp3), kind=dp)
    end function caroli_transmission

    subroutine rgf_step(cE, h_i, U_couple, G_NN, G_0N)
        complex(dp), intent(in) :: cE(:, :), h_i(:,:), U_couple(:,:)
        complex(dp), intent(inout) :: G_NN(:,:), G_0N(:,:)

        complex(dp), allocatable :: tmp1(:,:), tmp2(:,:), z(:,:)

        call matmul2(G_NN, U_couple, tmp1)
        call matmul2(U_couple, tmp1, tmp2, transa="C")

        z = cE - h_i - tmp2
        call invert(z)

        call matmul2(U_couple, G_NN, tmp1)
        call matmul2(G_0N, tmp1, tmp2)

        G_NN = z
        G_0N = tmp2
    end subroutine rgf_step

    subroutine rgf_last_step(g_R, U_NNp1, G_NN, G_0N, G_0Np1)
        complex(dp), intent(in) :: g_R(:,:), G_NN(:,:), G_0N(:,:), U_NNp1(:,:)
        complex(dp), intent(out) :: G_0Np1(:,:)

        complex(dp), allocatable :: gR_inv(:,:), z(:,:)
        complex(dp), allocatable :: tmp1(:,:), tmp2(:,:)

        gR_inv = g_R
        call invert(gR_inv)

        call matmul2(G_NN, U_NNp1, tmp1)
        call matmul2(U_NNp1, tmp1, tmp2, transa="C")

        z = gR_inv - tmp2
        call invert(z)

        call matmul2(U_NNp1, z, tmp1)
        call matmul2(G_0N, tmp1, tmp2)
        G_0Np1 = tmp2
    end subroutine rgf_last_step

end module transmittance
