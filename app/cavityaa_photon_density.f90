program main
   use stdlib_kinds, only: dp
   use stdlib_io_npy, only: save_npy
   use stdlib_linalg, only: eigh

   use aubry_andre, only: cavaa_hamiltonian, photon_probability
   implicit none

   character(len=256) :: outname
   integer :: L, Nph
   real(dp) :: t, V
   real(dp) :: gam, omega
   real(dp), parameter :: PHI = 0.0_dp
   real(dp), parameter :: INV_PHI = (sqrt(5.0_dp) - 1.0_dp)/2.0_dp

   integer :: NN
   complex(dp), allocatable :: H(:, :)
   real(dp), allocatable :: egv(:), Pph(:, :)
   ! character(len=32) :: jstr

   call readInput()

   NN = L*(Nph + 1)

   call writeInput("parameters_"//trim(outname)//".txt")

   allocate (H(NN, NN), egv(NN))

   call cavaa_hamiltonian(H, L, Nph, t, V, PHI, gam, omega)
   call eigh(H, egv, overwrite_a=.true.)
   call photon_probability(Pph, H, L, Nph)

   call save_npy("energies_"//trim(outname)//".npy", egv)
   call save_npy("photon_prob_"//trim(outname)//".npy", Pph)

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
