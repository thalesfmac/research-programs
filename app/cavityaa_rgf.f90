program main
   use stdlib_kinds, only: dp
   use rng_utils
   use array_io, only: geomspace_int
   use aubry_andre, only: aa_random_phases, energy_grid, cavaa_rgf_transmission
   use stdlib_io_npy, only: save_npy
   implicit none

   character(len=256) :: outname
   integer :: Nph, seed, Ndisorder, NEpoints, NLpoints
   real(dp) :: t, V
   real(dp) :: tcS, tcD, tlead, muS, muD
   real(dp) :: omega, g

   integer, allocatable :: lengths(:)
   real(dp), allocatable :: energies(:), phis(:)
   real(dp), allocatable :: transmissions(:, :, :)
   real(dp) :: Emin, Emax
   integer :: Lmin, Lmax
   real(dp), parameter :: INV_PHI = (sqrt(5.0_dp) - 1.0_dp)/2.0_dp
   real(dp), parameter :: ETA = 1.0e-10_dp

   integer :: i, j, k
   ! character(len=32) :: jstr

   call readInput()

   call writeInput("parameters_"//trim(outname)//".txt")

   call rng_initialize(seed)

   allocate (energies(NEpoints))
   allocate (phis(Ndisorder))
   allocate (transmissions(NLpoints, NEpoints, Ndisorder))

   lengths = geomspace_int(start=Lmin, stp=Lmax, num=NLpoints)
   call energy_grid(Egrid=energies, Emin=Emin, Emax=Emax)
   call aa_random_phases(phis)

   do k = 1, NLpoints
      do i = 1, NEpoints
         do j = 1, Ndisorder
            transmissions(k, i, j) = cavaa_rgf_transmission( &
                                     E=energies(i), &
                                     eta=ETA, &
                                     Lx=lengths(k), &
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
         end do

         write (*, *) "Done: L = ", lengths(k), "E = ", energies(i)
      end do
   end do

   call save_npy("transmissions_"//trim(outname)//".npy", transmissions)
   call save_npy("energies_"//trim(outname)//".npy", energies)
   call save_npy("lengths_"//trim(outname)//".npy", lengths)
   call save_npy("aa_phases_"//trim(outname)//".npy", phis)

contains

   subroutine readInput()
      use, intrinsic :: iso_fortran_env, only: input_unit
      read (input_unit, *) outname
      read (input_unit, *) NEpoints, Ndisorder
      read (input_unit, *) Lmin, Lmax, NLpoints
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

      open (newunit=unit, file=filename, status="replace", action="write")

      write (unit, *) "Input data"
      write (unit, *) "outname=", trim(outname)
      write (unit, *) "Size grid=", NLpoints
      write (unit, *) "Lmin=", Lmin
      write (unit, *) "Lmax=", Lmax
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
