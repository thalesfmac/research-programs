program main
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use precision
    use constants, only : INV_PHI
    use rng_utils
    use array_io
    use disordered_systems
    implicit none

    character(len=80) :: outname
    integer :: Lx, Nph, seed, Ndisorder, NEpoints
    real(dp) :: t, V
    real(dp) :: tcS, tcD, tlead, muS, muD
    real(dp) :: omega, g

    real(dp), allocatable :: energies(:), phis(:), transmissions(:, :)
    real(dp) :: Emin, Emax
    real(dp), parameter :: ETA = 1.0e-10_dp

    integer :: i, j

    call readInput()

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
    call aa_random_phases(phis)

    do i = 1, NEpoints
        do j = 1, Ndisorder
            transmissions(i, j) = cavaa_rgf_transmission( &
                E      = energies(i),  &
                eta    = ETA,          &
                Lx     = Lx,           &
                Nph    = Nph,          &
                t      = t,            &
                V      = V,            &
                beta   = INV_PHI,      &
                phi    = phis(j),      &
                g      = g,            &
                omega  = omega,        &
                tcL    = tcS,          &
                tcR    = tcD,          &
                tlead  = tlead,        &
                muL    = muS,          &
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
        read(input_unit,*) t, V
        read(input_unit,*) omega, g
        read(input_unit,*) tcS, tcD, tlead, muS, muD
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
        write(unit, *) "beta=", INV_PHI
        write(unit, *) "omega=", omega
        write(unit, *) "g=", g
        write(unit, *) "tcS=", tcS
        write(unit, *) "tcD=", tcD
        write(unit, *) "tlead=", tlead
        write(unit, *) "muS=", muS
        write(unit, *) "muD=", muD

        close(unit)
    end subroutine writeInput

end program main
