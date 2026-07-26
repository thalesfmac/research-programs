program main
   use stdlib_kinds, only: dp
   use stdlib_random, only: random_seed
   use stdlib_stats_distribution_uniform, only: rvs_uniform
   use stdlib_constants, only: PI => PI_dp
   use stdlib_io_npy, only: save_npy

   use aubry_andre, only: cavaa_rgf_transmission
   implicit none

   character(len=256) :: outname
   integer :: seed, actual_seed
   integer :: Lx, Lmin, Lmax, NL, Nph, Ndisorder
   real(dp) :: t, V
   real(dp) :: tcS, tcD, tlead, muS, muD
   real(dp) :: omega, g
   real(dp) :: E

   real(dp), allocatable :: lengths(:), phis(:), transmissions(:, :)
   real(dp), parameter :: INV_PHI = (sqrt(5.0_dp) - 1.0_dp)/2.0_dp
   real(dp), parameter :: ETA = 1.0e-10_dp

   integer :: i, j

   call readInput()

   NL = Lmax - Lmin + 1
   if (NL <= 0) then
      error stop "main: Lmax must be greater than or equal to Lmin"
   end if

   ! Optional, but useful: check if E lies inside the common lead band
   if (E < max(muS - 2.0_dp*tlead, muD - 2.0_dp*tlead) .or. &
       E > min(muS + 2.0_dp*tlead, muD + 2.0_dp*tlead)) then
      error stop "main: E is outside the common propagating band of the leads"
   end if

   call writeInput("parameters_"//trim(outname)//".txt")

   call random_seed(put=seed, get=actual_seed)

   allocate (lengths(NL))
   allocate (transmissions(NL, Ndisorder))

   phis = rvs_uniform(loc=0.0_dp, scale=2.0_dp*PI, array_size=Ndisorder)

   do i = 1, NL
      Lx = Lmin + i - 1
      lengths(i) = real(Lx, dp)

      do j = 1, Ndisorder
         transmissions(i, j) = cavaa_rgf_transmission( &
                               E=E, &
                               eta=ETA, &
                               Lx=Lx, &
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

      write (*, *) "Lx =", Lx
   end do

   call save_npy("lengths_"//trim(outname)//".npy", lengths)
   call save_npy("aa_phases_"//trim(outname)//".npy", phis)
   call save_npy("transmissions_"//trim(outname)//".npy", transmissions)

   deallocate (lengths, phis, transmissions)

contains

   subroutine readInput()
      use, intrinsic :: iso_fortran_env, only: input_unit
      read (input_unit, *) outname
      read (input_unit, *) Lmin, Lmax, Ndisorder
      read (input_unit, *) Nph
      read (input_unit, *) seed
      read (input_unit, *) E
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
      write (unit, *) "Lmin=", Lmin
      write (unit, *) "Lmax=", Lmax
      write (unit, *) "Number of L values=", NL
      write (unit, *) "Nph=", Nph
      write (unit, *) "seed=", seed
      write (unit, *) "E=", E
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
