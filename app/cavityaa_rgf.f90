module system_procedures
    use, intrinsic :: iso_fortran_env, only : int32
    use precision, only : dp
    use constants, only : CI
    use lapack_blas, only : invert
    use matrix_operations, only : identity_matrix
    use peierls_operator
    use rng_utils
    use array_io
    implicit none

    private
    public :: rng_initialize, energy_grid, random_phases, rgf_transmission
    public :: save_array_1d, save_array_2d, save_array_bin

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

    subroutine random_phases(phi_vals)
        real(dp), intent(out) :: phi_vals(:)
        real(dp), parameter :: PI = acos(-1.0_dp)

        call random_number(phi_vals)
        phi_vals = 2.0_dp * PI * phi_vals
    end subroutine random_phases

    subroutine slice_hamiltonian(h_i, i, V, beta, phi, omega)
        integer, intent(in) :: i
        real(dp), intent(in) :: V, beta, phi, omega
        complex(dp), intent(out) ::  h_i(0:, 0:)

        integer :: Nph, n
        real(dp), parameter :: PI = acos(-1.0_dp)

        if (lbound(h_i,1) /= 0 .or. lbound(h_i,2) /= 0) then
            error stop "slice_hamiltonian: h_i deve ter limites inferiores 0,0"
        end if

        if (ubound(h_i,1) /= ubound(h_i,2)) then
            error stop "slice_hamiltonian: h_i deve ser quadrada"
        end if

        h_i = (0.0_dp, 0.0_dp)

        Nph = ubound(h_i, dim=1)
        do n = 0, Nph
            h_i(n, n) = V * cos(2.0_dp * PI * beta * real(i, kind=dp) + phi) + real(n, kind=dp) * omega
        end do

    end subroutine slice_hamiltonian

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

    function rgf_transmission(E, eta, Lx, Nph, t, V, beta, phi, g, omega, tcL, tcR, tlead, muL, muR) result(TT)
        integer, intent(in) :: Lx, Nph
        real(dp), intent(in) :: E, eta, t, V, beta, phi, g, omega
        real(dp), intent(in) :: tcL, tcR, tlead, muL, muR
        real(dp) :: TT

        integer :: i
        complex(dp), dimension(0:Nph, 0:Nph) :: Id, U, z, h_i, G_NN
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
        call slice_hamiltonian(h_i, 1, V, beta, phi, omega)
        z = cmplx(E, eta, kind=dp) * Id - h_i - matmul( conjg(transpose(U_01)), matmul(g_L, U_01) )
        call invert(z)
        G_NN = z
        G_0N = matmul( g_L, matmul(U_01, G_NN) )

        ! Sítios internos: 2, ..., Lx
        do i = 2, Lx
            call slice_hamiltonian(h_i, i, V, beta, phi, omega)
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
    end function rgf_transmission

end module system_procedures


program main
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use precision
    use constants, only : INV_PHI
    use system_procedures
    implicit none

    character(len=80) :: outname
    integer :: Lx, Nph, seed, Ndisorder, NEpoints
    real(dp) :: t, V, beta
    real(dp) :: tcS, tcD, tlead, muS, muD
    real(dp) :: omega, g
    logical :: use_golden

    real(dp), allocatable :: energies(:), phis(:), transmissions(:, :)
    real(dp) :: E, Emin, Emax, phi
    real(dp), parameter :: ETA = 1.0e-10_dp

    integer :: i, j

    call readInput()

    if (use_golden) beta = INV_PHI

    call writeInput("parameters_" // trim(outname) // ".txt")

    call rng_initialize(seed)

    allocate(energies(NEpoints))
    allocate(phis(Ndisorder))
    allocate(transmissions(NEpoints, Ndisorder))

    Emin = max(muS - 2.0_dp*tlead, muD - 2.0_dp*tlead)
    Emax = min(muS + 2.0_dp*tlead, muD + 2.0_dp*tlead)

    call energy_grid( &
        Egrid = energies, &
        Emin = max(muS - 2.0_dp*tlead, muD - 2.0_dp*tlead), &
        Emax = min(muS + 2.0_dp*tlead, muD + 2.0_dp*tlead))
    call random_phases(phis)

    do i = 1, NEpoints
        E = energies(i)

        do j = 1, Ndisorder
            phi = phis(j)

            transmissions(i, j) = rgf_transmission( &
                E      = E,      &
                eta    = ETA,    &
                Lx     = Lx,     &
                Nph    = Nph,    &
                t      = t,      &
                V      = V,      &
                beta   = beta,   &
                phi    = phi,    &
                g      = g,      &
                omega  = omega,  &
                tcL    = tcS,    &
                tcR    = tcD,    &
                tlead  = tlead,  &
                muL    = muS,    &
                muR    = muD )
        end do
    end do

    call save_array_1d("energies_" // trim(outname) // ".dat", energies)
    call save_array_2d("transmissions_" // trim(outname) // ".dat", transmissions)
    call save_array_bin("transmissions_" // trim(outname) // ".bin", transmissions)

    deallocate(energies, phis, transmissions)

    contains

    subroutine readInput()
        use, intrinsic :: iso_fortran_env, only : input_unit
        read(input_unit,*) outname
        read(input_unit,*) Lx, NEpoints, Ndisorder
        read(input_unit,*) Nph
        read(input_unit,*) seed
        read(input_unit,*) t, V, beta
        read(input_unit,*) omega, g
        read(input_unit,*) tcS, tcD, tlead, muS, muD
        read(input_unit,*) use_golden
    end subroutine readInput

    subroutine writeInput(filename)
        character(len=*), intent(in) :: filename
        integer :: unit

        open(newunit=unit, file=filename, status="replace", action="write")

        write(unit, *) "Input data"
        write(unit, *) "outname=", trim(outname)
        write(unit, *) "Lx=", Lx
        write(unit, *) "Nph=", Nph
        write(unit, *) "seed=", seed
        write(unit, *) "Energy grid=", NEpoints
        write(unit, *) "Number of disorder conf.=", Ndisorder
        write(unit, *) "t=", t
        write(unit, *) "V=", V
        write(unit, *) "beta=", beta
        write(unit, *) "omega=", omega
        write(unit, *) "g=", g
        write(unit, *) "tcS=", tcS
        write(unit, *) "tcD=", tcD
        write(unit, *) "tlead=", tlead
        write(unit, *) "muS=", muS
        write(unit, *) "muD=", muD
        write(unit, *) "use_golden=", use_golden

        close(unit)
    end subroutine writeInput

end program main
