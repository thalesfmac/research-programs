program main
   use precision, only: dp
   use constants, only: INV_PHI
   use rng_utils
   use array_io
   use aubry_andre
   use matrix_operations, only: diagonalize
   implicit none

   character(len=256) :: outname
   integer :: L, Nph
   real(dp) :: t, V
   real(dp) :: gam, omega
   real(dp), parameter :: PHI = 0.0_dp

   integer :: NN
   complex(dp), allocatable :: H(:, :)
   real(dp), allocatable :: egv(:), Pph(:, :)
   ! character(len=32) :: jstr

   call readInput()

   NN = L*(Nph + 1)
   ! lengths_int = geomspace_int(start=Lmin, stop=Lmax, num=NLpoints)
   ! lengths_real = real(lengths_int, dp)

   call writeInput("parameters_"//trim(outname)//".txt")

   allocate (H(NN, NN), egv(NN))

   call cavaa_hamiltonian(H, L, Nph, t, V, INV_PHI, PHI, gam, omega)
   call diagonalize(H, egv)
   call photon_probability(Pph, H, L, Nph)

   call save_array_bin("energies_"//trim(outname)//".bin", egv)
   call save_array_bin("photon_prob_"//trim(outname)//".bin", Pph)

   ! call save_array_2d("hamiltonian_re.txt", H%re)
   ! call save_array_2d("hamiltonian_im.txt", H%im)
   ! allocate (energies(NEpoints))
   ! allocate (phis(Ndisorder))
   ! allocate (transmissions(NLpoints, NEpoints, Ndisorder))

   ! energy_grid(Egrid=energies, Emin=Emin, Emax=Emax)
   ! call aa_random_phases(phis)

   deallocate (H, egv, Pph)

contains

   subroutine readInput()
      use, intrinsic :: iso_fortran_env, only: input_unit
      read (input_unit, *) outname
      read (input_unit, *) L, Nph
      read (input_unit, *) t, V
      read (input_unit, *) gam, omega
   end subroutine readInput

   subroutine writeInput(filename)
      character(len=*), intent(in) :: filename
      integer :: unit

      open (newunit=unit, file=filename, status="replace", action="write")

      write (unit, *) "Input data"
      write (unit, *) "outname=", trim(outname)
      write (unit, *) "L=", L
      write (unit, *) "Nph=", Nph
      write (unit, *) "NN=", NN
      write (unit, *) "t=", t
      write (unit, *) "V=", V
      write (unit, *) "gamma=", gam
      write (unit, *) "omega=", omega
      write (unit, *) "beta=", INV_PHI
      write (unit, *) "phi=", PHI

      close (unit)
   end subroutine writeInput

end program main
