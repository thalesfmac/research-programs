module disordered_systems
    use :: precision, only : dp
    use :: constants, only : PI
    use :: lapack_blas
    use :: matrix_operations
    use :: lead_green_function
    use :: peierls_operator
    implicit none

    private
    public :: aa_random_phases
    public :: aa_onsite_potential, cavaa_slice_hamiltonian
    public :: energy_grid, cavaa_rgf_transmission

    contains

    subroutine energy_grid(Egrid, Emin, Emax)
        real(dp), intent(out) :: Egrid(:)
        real(dp), intent(in)  :: Emin, Emax

        integer :: i, N
        real(dp) :: dE

        N = size(Egrid)
        if (N <= 0) return

        dE = (Emax - Emin) / real(N + 1, dp)

        do i = 1, N
            Egrid(i) = Emin + real(i, dp) * dE
        end do
    end subroutine energy_grid

    subroutine aa_random_phases(phi_vals)
        real(dp), intent(out) :: phi_vals(:)

        call random_number(phi_vals)
        phi_vals = 2.0_dp * PI * phi_vals
    end subroutine aa_random_phases

    function aa_onsite_potential(V, site, beta, phi) result(V_i)
        real(dp), intent(in) :: V, beta, phi
        integer, intent(in) :: site
        real(dp) :: V_i

        V_i = V * cos(2.0_dp * PI * beta * real(site, kind=dp) + phi)
    end function aa_onsite_potential

    subroutine cavaa_slice_hamiltonian(h_i, i, V, beta, phi, omega)
        integer, intent(in) :: i
        real(dp), intent(in) :: V, beta, phi, omega
        complex(dp), intent(out) :: h_i(0:, 0:)

        integer :: Nph, n
        real(dp) :: V_i

        call assert_square(h_i, "h_i")
        Nph = ubound(h_i, dim=1)

        h_i = (0.0_dp, 0.0_dp)
        V_i = aa_onsite_potential(V, i, beta, phi)

        do n = 0, Nph
            h_i(n, n) = cmplx(V_i + real(n, kind=dp) * omega, kind=dp)
        end do
    end subroutine cavaa_slice_hamiltonian

    function cavaa_rgf_transmission(E, eta, Lx, Nph, t, V, beta, phi, g, omega, tcL, tcR, tlead, muL, muR) result(TT)
        integer, intent(in) :: Lx, Nph
        real(dp), intent(in) :: E, eta, t, V, beta, phi, g, omega
        real(dp), intent(in) :: tcL, tcR, tlead, muL, muR
        real(dp) :: TT

        integer :: i
        complex(dp), dimension(0:Nph, 0:Nph) :: Id, h_i, U, z, G_NN
        complex(dp) :: G_0N(0:0, 0:Nph)
        complex(dp) :: U_01(0:0, 0:Nph)
        complex(dp) :: U_NNp1(0:Nph, 0:0)
        complex(dp), dimension(0:0, 0:0) :: g_L, g_R, gR_inv, z_2, G_0Np1

        complex(dp), dimension(0:0, 0:0) :: sigma_L, sigma_R
        real(dp), dimension(0:0, 0:0) :: gamma_L, gamma_R

        call identity_matrix(Id)

        g_L(0, 0) = surface_gf_1d(E, tlead, muL)
        g_R(0, 0) = surface_gf_1d(E, tlead, muR)
        gR_inv = g_R; call invert(gR_inv)

        U_01 = (0.0_dp, 0.0_dp)
        U_01(0,0) = cmplx(-tcL, kind=dp)
        U_NNp1 = (0.0_dp, 0.0_dp)
        U_NNp1(0,0) = cmplx(-tcR, kind=dp)

        call peierls_exp(U, g)
        U = cmplx(-t, kind=dp) * U

        ! Caso especial: cadeia de um único sítio
        if (Lx == 1) then
            error stop "rgf_transmission: Lx = 1"
        end if

        ! Primeiro sítio
        call cavaa_slice_hamiltonian(h_i, 1, V, beta, phi, omega)
        z = cmplx(E, eta, kind=dp) * Id - h_i - matmul( conjg(transpose(U_01)), matmul(g_L, U_01) )
        call invert(z)
        G_NN = z
        G_0N = matmul( g_L, matmul(U_01, G_NN) )

        ! Sítios internos: 2, ..., Lx
        do i = 2, Lx
            call cavaa_slice_hamiltonian(h_i, i, V, beta, phi, omega)
            z = cmplx(E, eta, kind=dp) * Id - h_i - matmul( conjg(transpose(U)), matmul( G_NN, U ) )
            call invert(z)
            G_NN = z
            G_0N = matmul( G_0N, matmul( U, G_NN ) )
        end do

        ! Lead R
        z_2 = gR_inv - matmul( conjg(transpose(U_NNp1)), matmul( G_NN, U_NNp1 ) )
        call invert(z_2)
        G_0Np1 = matmul( G_0N, matmul( U_NNp1, z_2 ) )

        sigma_L = cmplx(tlead * tlead, kind=dp) * g_L
        sigma_R = cmplx(tlead * tlead, kind=dp) * g_R

        gamma_L = - 2.0_dp * aimag(sigma_L)
        gamma_R = - 2.0_dp * aimag(sigma_R)

        ! Tmat = abs( matmul( matmul(gamma_L, G_0Np1), matmul(gamma_R, conjg(transpose(G_0Np1))) ) )
        ! TT = real(Tmat(0,0), kind=dp)
        TT = gamma_L(0,0) * gamma_R(0,0) * abs(G_0Np1(0,0))**2
    end function cavaa_rgf_transmission

end module disordered_systems
