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


end module disordered_systems
