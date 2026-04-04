module disordered_systems
    use :: precision, only : dp
    use :: constants, only : PI
    implicit none

    public

    contains

    function aa_onsite_potential(V, site, beta, phi) result(V_i)
        real(dp), intent(in) :: V, beta, phi
        integer, intent(in) :: site

        real(dp) :: V_i

        V_i = V * cos(2.0_dp * PI * beta * real(site, kind=dp) + phi)
    end function aa_onsite_potential

    subroutine aa_slice_hamiltonian(h_i, i, V, beta, phi, omega)
        integer, intent(in) :: i
        real(dp), intent(in) :: V, beta, phi, omega
        real(dp), intent(out) ::  h_i(0:, 0:)

        integer :: Nph, n

        if (ubound(h_i,1) /= ubound(h_i,2)) then
            error stop "aa_slice_hamiltonian: h_i must be square"
        end if
        Nph = ubound(h_i, dim=1)

        h_i = 0.0_dp

        do n = 0, Nph
            h_i(n, n) = aa_onsite_potential(V, i, beta, phi) + real(n, kind=dp) * omega
        end do
    end subroutine aa_slice_hamiltonian


end module disordered_systems
