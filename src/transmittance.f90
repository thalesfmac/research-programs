module transmittance
    use :: precision, only : dp
    use :: lapack_blas, only : invert, gemm
    use :: matrix_operations
    implicit none

    private
    public :: caroli_transmission
    public :: rgf_step, rgf_last_step

    contains

    function caroli_transmission(gf, Gamma_L, Gamma_R) result(T)
        complex(dp), intent(in) :: gf(:,:), Gamma_L(:,:), Gamma_R(:,:)
        real(dp) :: T

        integer :: k, l, m, n, p, q
        complex(dp), allocatable :: tmp1(:,:), tmp2(:,:), tmp3(:,:)

        k = size(gf, dim=1); l = size(gf, dim=2)
        m = size(Gamma_L, dim=1); n = size(Gamma_L, dim=2)
        p = size(Gamma_R, dim=1); q = size(Gamma_R, dim=2)

        allocate(tmp1(m, l), tmp2(p, k), tmp3(m, k))

        call gemm(Gamma_L, gf, tmp1)
        call gemm(Gamma_R, gf, tmp2, transb="C")
        call gemm(tmp1, tmp2, tmp3)

        T = real(trace(tmp3), kind=dp)
    end function caroli_transmission

    subroutine rgf_step(cE, h_i, U_couple, G_NN, G_0N)
        complex(dp), intent(in) :: cE(:, :), h_i(:,:), U_couple(:,:)
        complex(dp), intent(inout) :: G_NN(:,:), G_0N(:,:)

        complex(dp), allocatable :: z(:,:)

        call assert_square(cE, "cE")
        call assert_square(h_i, "h_i")
        call assert_square(G_NN, "G_NN")

        z = cE - h_i - matmul3(conjg(transpose(U_couple)), G_NN, U_couple)
        call invert(z)

        G_NN = z
        G_0N = matmul3(G_0N, U_couple, G_NN)
    end subroutine rgf_step

    subroutine rgf_last_step(g_R, U_NNp1, G_NN, G_0N, G_0Np1)
        complex(dp), intent(in) :: g_R(:,:), G_NN(:,:), G_0N(:,:), U_NNp1(:,:)
        complex(dp), intent(out) :: G_0Np1(:,:)

        complex(dp), allocatable :: gR_inv(:,:), z_2(:,:)

        gR_inv = g_R
        call invert(gR_inv)

        z_2 = gR_inv - matmul3(conjg(transpose(U_NNp1)), G_NN, U_NNp1)
        call invert(z_2)

        G_0Np1 = matmul3(G_0N, U_NNp1, z_2)
    end subroutine rgf_last_step

end module transmittance
