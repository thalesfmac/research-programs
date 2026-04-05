module transmittance
    use :: precision, only : dp
    use :: lapack_blas, only : invert
    use :: matrix_operations
    implicit none

    private
    public :: caroli_transmission

    contains

    function caroli_transmission(G, Gamma_L, Gamma_R) result(T)
        complex(dp), intent(in) :: G(:,:), Gamma_L(:,:), Gamma_R(:,:)
        real(dp) :: T

        complex(dp), allocatable :: G_a(:,:), tmp(:,:)

        G_a = conjg(transpose(G))
        tmp = matmul4(Gamma_L, G, Gamma_R, G_a)

        T = real(trace(tmp), kind=dp)
    end function caroli_transmission

    subroutine rgf_step(E, eta, Id, h_i, U_couple, G_NN, G_0N)
        real(dp), intent(in) :: E, eta
        complex(dp), intent(in) :: Id(:,:), h_i(:,:), U_couple(:,:)
        complex(dp), intent(inout) :: G_NN(:,:), G_0N(:,:)

        complex(dp), allocatable :: z(:,:)

        call assert_square(Id, "Id")
        call assert_square(h_i, "h_i")
        call assert_square(U_couple, "U_couple")
        call assert_square(G_NN, "G_NN")

        z = cmplx(E, eta, kind=dp) * Id - h_i - matmul3(conjg(transpose(U_couple)), G_NN, U_couple)
        call invert(z)

        G_NN = z
        G_0N = matmul3(G_0N, U_couple, G_NN)
    end subroutine rgf_step

    subroutine rgf_last_step(G_NN, G_0N, g_R, U_NNp1, G_0Np1)
        complex(dp), intent(in) :: G_NN(:,:), G_0N(:,:), g_R(:,:), U_NNp1(:,:)
        complex(dp), intent(out) :: G_0Np1(:,:)

        complex(dp), allocatable :: gR_inv(:,:), z_2(:,:)

        gR_inv = g_R
        call invert(gR_inv)

        z_2 = gR_inv - matmul(conjg(transpose(U_NNp1)), matmul(G_NN, U_NNp1))
        call invert(z_2)

        G_0Np1 = matmul(G_0N, matmul(U_NNp1, z_2))
    end subroutine rgf_last_step

end module transmittance
