program main
    use precision, only : dp
    use constants, only : INV_PHI
    use rng_utils
    use array_io
    use disordered_systems
    implicit none

    character(len=256) :: outname
    integer :: Nph, seed, Ndisorder, NEpoints, NLpoints
    real(dp) :: t, V
    real(dp) :: tcS, tcD, tlead, muS, muD
    real(dp) :: omega, g

    integer, allocatable :: lengths_int(:)
    real(dp), allocatable :: energies(:), phis(:), lengths_real(:)
    real(dp), allocatable :: transmissions(:, :, :)
    real(dp) :: Emin, Emax
    integer :: Lmin, Lmax
    real(dp), parameter :: ETA = 1.0e-10_dp

    integer :: i, j, k
    ! character(len=32) :: jstr

    call readInput()

    lengths_int = geomspace_int(start=Lmin, stop=Lmax, num=NLpoints)
    lengths_real = real(lengths_int, dp)

    call writeInput("parameters_" // trim(outname) // ".txt")

    call rng_initialize(seed)

    allocate(energies(NEpoints))
    allocate(phis(Ndisorder))
    allocate(transmissions(NLpoints, NEpoints, Ndisorder))

    call energy_grid(Egrid = energies, Emin = Emin, Emax = Emax)
    call aa_random_phases(phis)

    do k = 1, NLpoints
        do i = 1, NEpoints
            do j = 1, Ndisorder
                transmissions(k, i, j) = cavaa_rgf_transmission( &
                    E      = energies(i),     &
                    eta    = ETA,             &
                    Lx     = lengths_int(k),  &
                    Nph    = Nph,             &
                    t      = t,               &
                    V      = V,               &
                    beta   = INV_PHI,         &
                    phi    = phis(j),         &
                    g      = g,               &
                    omega  = omega,           &
                    tcL    = tcS,             &
                    tcR    = tcD,             &
                    tlead  = tlead,           &
                    muL    = muS,             &
                    muR    = muD )
            end do

            write(*, *) "Done: L = ", lengths_int(k), "E = ", energies(i)
        end do
    end do

    call save_array_bin("transmissions_" // trim(outname) // ".bin", transmissions)
    call save_array_bin("energies_" // trim(outname) // ".bin", energies)
    call save_array_bin("lengths_" // trim(outname) // ".bin", lengths_real)

    ! call save_array_1d("energies_" // trim(outname) // ".dat", energies)
    ! call save_array_1d("lengths_" // trim(outname) // ".dat", lengths_real)

    ! do j = 1, Ndisorder
    !     write(jstr, '(I4.4)') j

    !     call save_array_2d("transmissions_" // trim(outname) // "_" // trim(jstr) //".dat", transmissions(:, :, j))
    ! end do


    deallocate(lengths_int, lengths_real, energies, phis, transmissions)

    contains

    subroutine readInput()
        use, intrinsic :: iso_fortran_env, only : input_unit
        read(input_unit,*) outname
        read(input_unit,*) NEpoints, Ndisorder
        read(input_unit,*) Lmin, Lmax, NLpoints
        read(input_unit,*) Emin, Emax
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
        write(unit, *) "Size grid=", NLpoints
        write(unit, *) "Lmin=", Lmin
        write(unit, *) "Lmax=", Lmax
        write(unit, *) "Nph=", Nph
        write(unit, *) "seed=", seed
        write(unit, *) "Energy grid=", NEpoints
        write(unit, *) "Emin=", Emin
        write(unit, *) "Emax=", Emax
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
