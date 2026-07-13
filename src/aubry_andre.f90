module aubry_andre
   use :: precision, only:dp
   use :: constants, only:PI, CI
   use :: matrix_operations, only:identity_matrix, assert_square
   use :: lead_green_function, only:surface_gf_1d, surface_self_energy_left, surface_self_energy_right, broadening
   use :: peierls_operator, only:peierls_exp
   use :: transmittance, only:caroli_transmission, rgf_first_step, rgf_step, rgf_last_step
   implicit none

   private
   public :: aa_random_phases
   public :: aa_onsite_potential, cavaa_slice_hamiltonian
   public :: energy_grid, cavaa_rgf_transmission

contains

   subroutine energy_grid(Egrid, Emin, Emax)
      real(dp), intent(out) :: Egrid(:)
      real(dp), intent(in)  :: Emin, Emax

      integer :: i, N
      real(dp) :: dE

      N = size(Egrid)
      if (N <= 0) return

      dE = (Emax - Emin)/real(N + 1, dp)

      do i = 1, N
         Egrid(i) = Emin + real(i, dp)*dE
      end do
   end subroutine energy_grid

   subroutine aa_random_phases(phi_vals)
      real(dp), intent(out) :: phi_vals(:)

      call random_number(phi_vals)
      phi_vals = 2.0_dp*PI*phi_vals
   end subroutine aa_random_phases

   function aa_onsite_potential(V, site, beta, phi) result(V_i)
      real(dp), intent(in) :: V, beta, phi
      integer, intent(in) :: site
      real(dp) :: V_i

      V_i = V*cos(2.0_dp*PI*beta*real(site, kind=dp) + phi)
   end function aa_onsite_potential

   subroutine cav_site_photon_to_index(site, Nph, L, idx)
      integer, intent(in)  :: site, Nph, L
      integer, intent(out) :: idx

      idx = site + L*nph
   end subroutine cav_site_photon_to_index

   subroutine cav_index_to_site_photon(idx, L, site, Nph)
      integer, intent(in)  :: idx, L
      integer, intent(out) :: site, Nph

      nph = (idx - 1)/L
      site = mod(idx - 1, L) + 1
   end subroutine cav_index_to_site_photon

   subroutine cavaa_hamiltonian(H, L, Nph, t, V, beta, phi, gam, omega)
      complex(dp), intent(out) :: H(:, :)
      integer :: L, Nph
      real(dp) :: t, V, beta, phi, gam, omega

      complex(dp) :: PE(0:Nph, 0:Nph)
      real(dp) :: g
      integer :: i, j, n, m

      g = gam/t
      call peierls_exp(PE, g)

      H = (0.0_dp, 0.0_dp)

      do n = 0, Nph
      end do

   end subroutine cavaa_hamiltonian

   subroutine cavaa_slice_hamiltonian(h_i, i, V, beta, phi, omega)
      integer, intent(in) :: i
      real(dp), intent(in) :: V, beta, phi, omega
      complex(dp), intent(out) :: h_i(0:, 0:)

      integer :: Nph, n
      real(dp) :: V_i

      call assert_square(h_i, "h_i", caller="cavaa_slice_hamiltonian")
      Nph = ubound(h_i, dim=1)

      h_i = (0.0_dp, 0.0_dp)
      V_i = aa_onsite_potential(V, i, beta, phi)

      do n = 0, Nph
         h_i(n, n) = cmplx(V_i + real(n, kind=dp)*omega, kind=dp)
      end do
   end subroutine cavaa_slice_hamiltonian

   function cavaa_rgf_transmission(E, eta, Lx, Nph, t, V, beta, phi, g, omega, tcL, tcR, tlead, muL, muR) result(tt)
      integer, intent(in) :: Lx, Nph
      real(dp), intent(in) :: E, eta, t, V, beta, phi, g, omega
      real(dp), intent(in) :: tcL, tcR, tlead, muL, muR
      real(dp) :: tt

      integer :: n
      complex(dp), dimension(0:Nph, 0:Nph) :: cE, h_n, U, G_nm1_nm1, G_n_n
      complex(dp), dimension(0:0, 0:Nph) :: U_01, G_0_nm1, G_0_n
      complex(dp), dimension(0:Nph, 0:0) :: U_N_Np1
      complex(dp), dimension(0:0, 0:0) :: g_L, g_R, G_Np1_Np1, G_0_Np1
      complex(dp), dimension(0:0, 0:0) :: u_left, u_right
      complex(dp), dimension(0:0, 0:0) :: sigma_L, sigma_R, gamma_L, gamma_R

      ! Check if the system has more than 1 site
      if (Lx <= 1) then
         error stop "cavaa_rgf_transmission: Lx must be greater than 1"
      end if

      call identity_matrix(cE)
      cE = cmplx(E, eta, kind=dp)*cE

      g_L(0, 0) = surface_gf_1d(E, tlead, muL)
      g_R(0, 0) = surface_gf_1d(E, tlead, muR)

      u_left(0, 0) = cmplx(-tlead, kind=dp)
      u_right(0, 0) = cmplx(-tlead, kind=dp)

      call surface_self_energy_left(g_L, u_left, sigma_L)
      call surface_self_energy_right(g_R, u_right, sigma_R)

      call broadening(sigma_L, gamma_L)
      call broadening(sigma_R, gamma_R)

      U_01 = (0.0_dp, 0.0_dp)
      U_01(0, 0) = cmplx(-tcL, kind=dp)
      U_N_Np1 = (0.0_dp, 0.0_dp)
      U_N_Np1(0, 0) = cmplx(-tcR, kind=dp)

      call peierls_exp(U, g)
      U = cmplx(-t, kind=dp)*U

      ! First site
      call cavaa_slice_hamiltonian(h_n, 1, V, beta, phi, omega)
      call rgf_first_step(cE, h_n, U_01, g_L, G_nm1_nm1, G_0_nm1)

      ! Internal sites: 2, ..., Lx
      do n = 2, Lx
         call cavaa_slice_hamiltonian(h_n, n, V, beta, phi, omega)
         call rgf_step(cE, h_n, U, G_nm1_nm1, G_0_nm1, G_n_n, G_0_n)
         G_nm1_nm1 = G_n_n
         G_0_nm1 = G_0_n
      end do

      ! Connect to the right lead
      call rgf_last_step(g_R, U_N_Np1, G_n_n, G_0_n, G_Np1_Np1, G_0_Np1)

      ! Calculate the transmission probability
      tt = caroli_transmission(G_0_Np1, gamma_L, gamma_R)
   end function cavaa_rgf_transmission

end module aubry_andre
