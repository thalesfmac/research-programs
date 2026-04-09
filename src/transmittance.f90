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

        allocate( tmp1(size(Gamma_L, 1), size(gf, 2)) )
        allocate( tmp2(size(Gamma_R, 1), size(gf, 1)) )
        allocate( tmp3(size(tmp1, 1), size(tmp2, 2)) )

        call matmul2(Gamma_L, gf, tmp1)
        call matmul2(Gamma_R, gf, tmp2, transb="C")
        call matmul2(tmp1, tmp2, tmp3)

        T = real(trace(tmp3), kind=dp)
    end function caroli_transmission

    subroutine rgf_step(cE, h_i, U_couple, G_NN, G_0N)
        complex(dp), intent(in) :: cE(:, :), h_i(:,:), U_couple(:,:)
        complex(dp), intent(inout) :: G_NN(:,:), G_0N(:,:)

        complex(dp), allocatable :: tmp1(:,:), tmp2(:,:), z(:,:)
        complex(dp), allocatable :: tmp3(:,:), tmp4(:,:)

        if (size(cE, 1) /= size(h_i, 1) .or. size(cE, 2) /= size(h_i, 2)) then
            error stop "rgf_step: cE has incompatible dimensions with h_i"
        end if

        allocate( tmp1(size(G_NN, 1), size(U_couple, 2)) )
        allocate( tmp2(size(U_couple, 2), size(tmp1, 2)) )
        allocate( z(size(cE, 1), size(cE, 2)) )

        call matmul2(G_NN, U_couple, tmp1)
        call matmul2(U_couple, tmp1, tmp2, transa="C")

        z = cE - h_i - tmp2
        call invert(z)

        allocate( tmp3(size(U_couple, 1), size(G_NN, 2)) )
        allocate( tmp4(size(G_0N, 1), size(tmp3, 2)) )

        call matmul2(U_couple, G_NN, tmp3)
        call matmul2(G_0N, tmp3, tmp4)

        G_NN = z
        G_0N = tmp4
    end subroutine rgf_step

    subroutine rgf_last_step(g_R, U_NNp1, G_NN, G_0N, G_0Np1)
        complex(dp), intent(in) :: g_R(:,:), G_NN(:,:), G_0N(:,:), U_NNp1(:,:)
        complex(dp), intent(out) :: G_0Np1(:,:)

        complex(dp), allocatable :: gR_inv(:,:), z(:,:)
        complex(dp), allocatable :: tmp1(:,:), tmp2(:,:), tmp3(:,:), tmp4(:,:)

        gR_inv = g_R
        call invert(gR_inv)

        allocate( tmp1(size(G_NN, 1), size(U_NNp1, 2)) )
        allocate( tmp2(size(U_NNp1, 2), size(tmp1, 2)) )
        allocate( z(size(gR_inv, 1), size(gR_inv, 2)) )

        call matmul2(G_NN, U_NNp1, tmp1)
        call matmul2(U_NNp1, tmp1, tmp2, transa="C")

        z = gR_inv - tmp2
        call invert(z)

        allocate( tmp3(size(U_NNp1, 1), size(z, 2)) )
        allocate( tmp4(size(G_0N, 2), size(tmp3, 2)) )

        call matmul2(U_NNp1, z, tmp3)
        call matmul2(G_0N, tmp3, tmp4)
        G_0Np1 = tmp4
    end subroutine rgf_last_step

end module transmittance
