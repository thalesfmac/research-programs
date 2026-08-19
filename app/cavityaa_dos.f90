program main
   use stdlib_kinds, only: dp
   use stdlib_random, only: random_seed
   use stdlib_stats_distribution_uniform, only: rvs_uniform
   use stdlib_constants, only: PI => PI_dp
   use stdlib_linalg, only: eigvalsh
   use stdlib_io_npy, only: save_npy
   use stdlib_datetime, only: datetime_type, timedelta_type, now, format_timedelta, operator(-)

   use array_utils, only: geomspace_int
   use aubry_andre, only: energy_grid, cavaa_rgf_transmission, cavaa_hamiltonian
   use density_of_states, only: dos
   implicit none

   character(len=256) :: outname
   integer :: seed, actual_seed
   integer :: Nph, L, Ndisorder, NEpoints
   real(dp) :: t, V
   real(dp) :: tcS, tcD, tlead, muS, muD
   real(dp) :: omega, g

   integer, allocatable :: lengths(:)
   real(dp), allocatable :: energies(:), phis(:)
   real(dp), allocatable :: transmissions(:, :)
   real(dp) :: Emin, Emax
   ! integer :: Lmin, Lmax
   real(dp), parameter :: ETA = 1.0e-13_dp

   complex(dp), allocatable :: H(:, :)
   real(dp), allocatable :: eig(:), w(:, :), d(:, :)

   integer :: i, j

   type(datetime_type) :: time_start, time_end
   type(timedelta_type) :: elapsed

   time_start = now()

   call readInput()

   call writeInput("parameters_"//trim(outname)//".txt")

   call random_seed(put=seed, get=actual_seed)

   allocate (energies(NEpoints))
   allocate (transmissions(NEpoints, Ndisorder))
   allocate (w(L*(Nph + 1), Ndisorder))
   allocate (d(NEpoints, Ndisorder))

   call energy_grid(Egrid=energies, Emin=Emin, Emax=Emax)
   phis = rvs_uniform(loc=0.0_dp, scale=2.0_dp*PI, array_size=Ndisorder)

   do j = 1, Ndisorder
      call cavaa_hamiltonian(H=H, L=L, Nph=Nph, t=t, V=V, phi=phis(j), gam=g, omega=omega)

      eig = eigvalsh(H)
      w(:, j) = eig

      do i = 1, NEpoints
         transmissions(i, j) = cavaa_rgf_transmission( &
                               E=energies(i), &
                               eta=ETA, &
                               Lx=L, &
                               Nph=Nph, &
                               t=t, &
                               V=V, &
                               phi=phis(j), &
                               g=g, &
                               omega=omega, &
                               tcL=tcS, &
                               tcR=tcD, &
                               tlead=tlead, &
                               muL=muS, &
                               muR=muD)

         d(i, j) = dos(H, energies(i), ETA)

         write (*, '("Done: E = ",g0)') energies(i)
      end do

   end do

   call save_npy("transmissions_"//trim(outname)//".npy", transmissions)
   call save_npy("energies_"//trim(outname)//".npy", energies)
   call save_npy("lengths_"//trim(outname)//".npy", lengths)
   call save_npy("aa_phases_"//trim(outname)//".npy", phis)
   call save_npy("eigvals_"//trim(outname)//".npy", w)
   call save_npy("dos_"//trim(outname)//".npy", d)

   time_end = now()
   elapsed = time_end - time_start
   write (*, '("Execution time: ",a)') format_timedelta(elapsed)

contains

   subroutine readInput()
      use, intrinsic :: iso_fortran_env, only: input_unit
      read (input_unit, *) outname
      read (input_unit, *) NEpoints, Ndisorder
      read (input_unit, *) L
      read (input_unit, *) Emin, Emax
      read (input_unit, *) Nph
      read (input_unit, *) seed
      read (input_unit, *) t, V
      read (input_unit, *) omega, g
      read (input_unit, *) tcS, tcD, tlead, muS, muD
   end subroutine readInput

   subroutine writeInput(filename)
      character(len=*), intent(in) :: filename
      integer :: unit

      real(dp), parameter :: INV_PHI = (sqrt(5.0_dp) - 1.0_dp)/2.0_dp

      open (newunit=unit, file=filename, status="replace", action="write")

      write (unit, *) "Input data"
      write (unit, *) "outname=", trim(outname)
      write (unit, *) "L=", L
      write (unit, *) "Nph=", Nph
      write (unit, *) "seed=", seed
      write (unit, *) "Energy grid=", NEpoints
      write (unit, *) "Emin=", Emin
      write (unit, *) "Emax=", Emax
      write (unit, *) "Number of disorder conf.=", Ndisorder
      write (unit, *) "t=", t
      write (unit, *) "V=", V
      write (unit, *) "beta=", INV_PHI
      write (unit, *) "omega=", omega
      write (unit, *) "g=", g
      write (unit, *) "tcS=", tcS
      write (unit, *) "tcD=", tcD
      write (unit, *) "tlead=", tlead
      write (unit, *) "muS=", muS
      write (unit, *) "muD=", muD

      close (unit)
   end subroutine writeInput

end program main
