module disordered_systems
    use :: precision, only : dp
    use :: constants, only : PI
    implicit none

    private
    public :: aa_random_phases
    public :: aa_onsite_potential, aa_slice_hamiltonian

    contains

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

    subroutine aa_slice_hamiltonian(h_i, i, V, beta, phi, omega)
        integer, intent(in) :: i
        real(dp), intent(in) :: V, beta, phi, omega
        complex(dp), intent(out) ::  h_i(0:, 0:)

        integer :: Nph, n
        real(dp) :: V_i

        if (ubound(h_i,1) /= ubound(h_i,2)) then
            error stop "aa_slice_hamiltonian: h_i must be square"
        end if
        Nph = ubound(h_i, dim=1)

        h_i = (0.0_dp, 0.0_dp)
        V_i = aa_onsite_potential(V, i, beta, phi)

        do n = 0, Nph
            h_i(n, n) = cmplx(V_i + real(n, kind=dp) * omega, kind=dp)
        end do
    end subroutine aa_slice_hamiltonian


end module disordered_systems
