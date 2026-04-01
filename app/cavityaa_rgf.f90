module lapack_blas_interfaces
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private

    ! Wrappers “seguros”
    public :: diagonalize, gemm, invert

    interface
        subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
            import :: real64
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            integer, intent(out) :: info
            complex(real64), intent(inout) :: a(lda,*)
            real(real64),    intent(out)   :: w(*)
            complex(real64), intent(inout) :: work(*)
            real(real64),    intent(inout) :: rwork(*)
        end subroutine zheev

        subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: real64
            character(len=1), intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            complex(real64), intent(in) :: alpha, beta
            complex(real64), intent(in) :: a(lda,*), b(ldb,*)
            complex(real64), intent(inout) :: c(ldc,*)
        end subroutine zgemm

        subroutine zgetrf(m,n,a,lda,ipiv,info)
            import :: real64
            integer, intent(in) :: m,n,lda
            integer, intent(out) :: ipiv(*), info
            complex(real64), intent(inout) :: a(lda,*)
        end subroutine zgetrf

        subroutine zgetri(n,a,lda,ipiv,work,lwork,info)
            import :: real64
            integer, intent(in) :: n, lda, lwork
            integer, intent(in) :: ipiv(*)
            integer, intent(out) :: info
            complex(real64), intent(inout) :: a(lda,*), work(*)
        end subroutine zgetri
    end interface

    contains

    subroutine diagonalize(A, w, jobz, uplo)
        ! Diagonaliza Hermitiana complexa via ZHEEV.
        complex(real64), intent(inout), contiguous :: A(:,:)
        real(real64),    intent(out)   :: w(:)
        character(len=1), intent(in), optional :: jobz, uplo

        character(len=1) :: jobz_loc, uplo_loc
        integer :: n, lda, info, lwork
        complex(real64), allocatable :: work(:)
        real(real64),    allocatable :: rwork(:)
        complex(real64) :: workq(1)

        jobz_loc = 'V'; if (present(jobz)) jobz_loc = jobz
        uplo_loc = 'U'; if (present(uplo)) uplo_loc = uplo

        n = size(A,1)
        if (size(A,2) /= n) error stop "diagonalize: A deve ser quadrada"
        if (size(w)   /= n) error stop "diagonalize: w deve ter tamanho n"
        lda = n

        allocate(rwork(max(1, 3*n - 2)))

        ! Workspace query
        lwork = -1
        call zheev(jobz_loc, uplo_loc, n, A, lda, w, workq, lwork, rwork, info)
        if (info /= 0) then
            deallocate(rwork)
            error stop "diagonalize: ZHEEV query falhou"
        end if

        lwork = int(real(workq(1), real64))
        if (lwork < 1) lwork = 1
        allocate(work(lwork))

        ! Diagonalização
        call zheev(jobz_loc, uplo_loc, n, A, lda, w, work, lwork, rwork, info)

        deallocate(work, rwork)
        if (info /= 0) error stop "diagonalize: ZHEEV falhou"
    end subroutine diagonalize


    subroutine gemm(A, B, C, transa, transb, alpha, beta)
        ! Wrapper para ZGEMM: C <- alpha op(A) op(B) + beta C
        complex(real64), intent(in), contiguous    :: A(:,:), B(:,:)
        complex(real64), intent(inout), contiguous :: C(:,:)
        character(len=1), intent(in), optional :: transa, transb
        complex(real64), intent(in), optional :: alpha, beta

        character(len=1) :: ta, tb
        complex(real64) :: a_loc, b_loc
        integer :: m, n, k
        integer :: a_rows, a_cols, b_rows, b_cols
        integer :: lda, ldb, ldc

        ta = 'N'; if (present(transa)) ta = transa
        tb = 'N'; if (present(transb)) tb = transb

        if (present(alpha)) then
            a_loc = alpha
        else
            a_loc = (1.0_real64, 0.0_real64)
        end if
        if (present(beta)) then
            b_loc = beta
        else
            b_loc = (0.0_real64, 0.0_real64)
            C = (0.0_real64, 0.0_real64)
        end if

        select case (ta)
          case ('N','n')
            a_rows = size(A,1); a_cols = size(A,2)
          case ('T','t','C','c')
            a_rows = size(A,2); a_cols = size(A,1)
          case default
            error stop "gemm: transa deve ser 'N','T' ou 'C'"
        end select

        select case (tb)
          case ('N','n')
            b_rows = size(B,1); b_cols = size(B,2)
          case ('T','t','C','c')
            b_rows = size(B,2); b_cols = size(B,1)
          case default
            error stop "gemm: transb deve ser 'N','T' ou 'C'"
        end select

        m = a_rows
        k = a_cols
        n = b_cols

        if (b_rows /= k) error stop "gemm: dimensoes incompativeis"
        if (size(C,1) /= m .or. size(C,2) /= n) error stop "gemm: C tem dimensoes erradas"

        lda = size(A,1)
        ldb = size(B,1)
        ldc = size(C,1)

        call zgemm(ta, tb, m, n, k, a_loc, A, lda, B, ldb, b_loc, C, ldc)
    end subroutine gemm

    subroutine invert(A)
        complex(real64), intent(inout), contiguous :: A(:,:)
        integer :: n, info, lwork
        integer, allocatable :: ipiv(:)
        complex(real64), allocatable :: work(:)

        n = size(A,1)
        if (size(A,2) /= n) error stop "invert: A deve ser quadrada"

        allocate(ipiv(n))

        call zgetrf(n,n,A,n,ipiv,info)
        if (info /= 0) error stop "invert: ZGETRF falhou"

        ! workspace query
        lwork = -1
        allocate(work(1))
        call zgetri(n,A,n,ipiv,work,lwork,info)
        if (info /= 0) error stop "invert: ZGETRI query falhou"
        lwork = int(real(work(1), real64))
        deallocate(work)

        allocate(work(max(1,lwork)))
        call zgetri(n,A,n,ipiv,work,lwork,info)
        if (info /= 0) error stop "invert: ZGETRI falhou"

        deallocate(work, ipiv)
    end subroutine invert


end module lapack_blas_interfaces

module hdf5_io
    use, intrinsic :: iso_fortran_env, only : dp => real64
    use hdf5
    implicit none
    private
    public :: save_results_h5

    contains

    subroutine save_results_h5(filename, energies, transmissions)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: energies(:)
        real(dp), intent(in) :: transmissions(:,:)

        integer(hid_t) :: file_id
        integer(hid_t) :: space1_id, dset1_id
        integer(hid_t) :: space2_id, dset2_id
        integer(hsize_t), dimension(1) :: dims1
        integer(hsize_t), dimension(2) :: dims2
        integer :: error

        dims1 = [ int(size(energies,1), kind=hsize_t) ]
        dims2 = [ int(size(transmissions,1), kind=hsize_t), &
            int(size(transmissions,2), kind=hsize_t) ]

        call h5open_f(error)
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)

        call h5screate_simple_f(1, dims1, space1_id, error)
        call h5dcreate_f(file_id, "energies", H5T_NATIVE_DOUBLE, space1_id, dset1_id, error)
        call h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, energies, dims1, error)

        call h5screate_simple_f(2, dims2, space2_id, error)
        call h5dcreate_f(file_id, "transmissions", H5T_NATIVE_DOUBLE, space2_id, dset2_id, error)
        call h5dwrite_f(dset2_id, H5T_NATIVE_DOUBLE, transmissions, dims2, error)

        call h5dclose_f(dset1_id, error)
        call h5sclose_f(space1_id, error)

        call h5dclose_f(dset2_id, error)
        call h5sclose_f(space2_id, error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)
    end subroutine save_results_h5

end module hdf5_io


module system_procedures
    use, intrinsic :: iso_fortran_env, only : int32, dp => real64
    use lapack_blas_interfaces, only : invert
    implicit none

    complex(dp), parameter :: CI = (0.0_dp, 1.0_dp)

    private
    public :: rng_initialize, energy_grid, random_phases, rgf_transmission
    public :: save_array_1d, save_array_2d, save_array_bin

    contains

    subroutine rng_initialize(seed)
        integer, intent(in) :: seed

        integer :: n, j
        integer, allocatable :: seed_vec(:)

        call random_seed(size=n)
        allocate(seed_vec(n))

        ! seed_vec = seed
        seed_vec = seed + 37 * [(j - 1, j=1,n)]

        call random_seed(put=seed_vec)
        deallocate(seed_vec)
    end subroutine rng_initialize

    subroutine identity_matrix(I)
        complex(dp), intent(out) :: I(:,:)
        integer :: i1, i2, j

        if (lbound(I,1) /= lbound(I,2) .or. ubound(I,1) /= ubound(I,2)) then
            error stop "identity_matrix: I deve ser quadrada com mesmos limites nas duas dimensoes"
        end if

        i1 = lbound(I,1)
        i2 = ubound(I,1)

        I = (0.0_dp, 0.0_dp)
        do j = i1, i2
            I(j,j) = (1.0_dp, 0.0_dp)
        end do
    end subroutine identity_matrix

    subroutine save_array_1d(filename, A)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: A(:)

        integer :: u, i

        open(newunit=u, file=filename, status="replace", action="write")

        do i = 1, size(A)
            ! write(u, '(ES24.16)') A(i)
            write(u, '(E26.16E3)') A(i)
        end do

        close(u)
    end subroutine save_array_1d

    subroutine save_array_2d(filename, A)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: A(:,:)

        integer :: u, i, j, ncol

        open(newunit=u, file=filename, status="replace", action="write")

        ncol = size(A, 2)

        do i = 1, size(A, 1)
            do j = 1, ncol
                ! write(u, '(ES24.16)', advance='no') A(i,j)
                write(u, '(E26.16E3)', advance='no') A(i,j)
                if (j < ncol) write(u, '(A)', advance='no') ' '
            end do
            write(u, *)
        end do

        close(u)
    end subroutine save_array_2d

    subroutine save_array_bin(filename, A)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: A(..)

        integer :: u
        integer(int32) :: r
        integer(int32), allocatable :: shp(:)

        open(newunit=u, file=filename, status="replace", action="write", &
            form="unformatted", access="stream")

        select rank(A)

          rank(1)
            r = 1_int32
            allocate(shp(1))
            shp = int(shape(A), kind=int32)

            write(u) r
            write(u) shp
            write(u) A

          rank(2)
            r = 2_int32
            allocate(shp(2))
            shp = int(shape(A), kind=int32)

            write(u) r
            write(u) shp
            write(u) A

          rank(3)
            r = 3_int32
            allocate(shp(3))
            shp = int(shape(A), kind=int32)

            write(u) r
            write(u) shp
            write(u) A

          rank default
            close(u)
            error stop "save_array_bin: rank nao suportado"
        end select

        close(u)

        if (allocated(shp)) deallocate(shp)
    end subroutine save_array_bin

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

    subroutine build_factorials(fact)
        real(dp), intent(out) :: fact(0:)
        integer :: k, Nph

        Nph = ubound(fact, dim=1)

        fact(0) = 1.0_dp
        do k = 1, Nph
            fact(k) = fact(k-1) * real(k, dp)
        end do
    end subroutine build_factorials

    ! Replica o funcoesP antigo:
    ! P(0,M)=1 e P(j,M)=sqrt(M-(j-1))*P(j-1,M)  => P(j,M)=sqrt(M!/(M-j)!)
    subroutine build_P(P)
        real(dp), intent(out) :: P(0:, 0:)
        integer :: M, j, Nph

        if (lbound(P,1) /= 0 .or. lbound(P,2) /= 0) then
            error stop "build_P: P deve ter limites inferiores 0,0"
        end if

        if (ubound(P,1) /= ubound(P,2)) then
            error stop "build_P: P deve ser quadrada"
        end if

        Nph = ubound(P, 1)

        P = 0.0_dp
        P(0, :) = 1.0_dp

        do M = 1, Nph
            do j = 1, M
                P(j, M) = P(j-1, M) * sqrt(real(M - (j-1), kind=dp))
            end do
        end do
    end subroutine build_P

    subroutine peierls_exp(A, g)
        complex(dp), intent(out) :: A(0:, 0:)
        real(dp), intent(in) :: g

        integer :: Nph
        real(dp), allocatable :: fact(:)
        real(dp), allocatable :: P(:, :)
        integer :: M, N, s, j
        complex(dp) :: hNM
        real(dp) :: pref

        if (lbound(A,1) /= 0 .or. lbound(A,2) /= 0) then
            error stop "peierls_exp: A deve ter limites inferiores 0,0"
        end if

        if (ubound(A,1) /= ubound(A,2)) then
            error stop "peierls_exp: A deve ser quadrada"
        end if

        Nph = ubound(A,1)

        allocate( fact(0:Nph) )
        allocate( P(0:Nph, 0:Nph) )

        call build_factorials(fact)
        call build_P(P)

        pref = exp(-0.5_dp * (g*g))
        A = (0.0_dp, 0.0_dp)

        ! índices físicos N,M = 0..Nph  (no array: +1)
        do N = 0, Nph
            do M = 0, Nph
                hNM = (0.0_dp, 0.0_dp)

                do s = 0, N
                    do j = 0, M
                        ! delta(N-s, M-j) == 1  <=>  N - s == M - j
                        if (N - s == M - j) then
                            hNM = hNM + pref * ( (CI*g)**s ) * ( (CI*g)**j ) &
                                * P(s, N) * P(j, M) / ( fact(s) * fact(j) )
                        end if
                    end do
                end do

                A(N, M) = hNM
            end do
        end do
    end subroutine peierls_exp

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
    use, intrinsic :: iso_fortran_env, only : dp => real64, input_unit
    use system_procedures
    implicit none

    character(len=80) :: outname
    integer :: Lx, Nph, seed, Ndisorder, NEpoints
    real(dp) :: t, V, beta
    real(dp) :: tcS, tcD, tlead, muS, muD
    real(dp) :: omega, g
    logical :: use_golden

    ! complex(dp), parameter :: CI = (0.0_dp, 1.0_dp)
    real(dp), parameter :: GOLDEN = (sqrt(5.0_dp) - 1.0_dp) / 2.0_dp

    real(dp), allocatable :: energies(:), phis(:), transmissions(:, :)
    real(dp) :: E, Emin, Emax, phi
    real(dp), parameter :: ETA = 1.0e-10_dp

    integer :: i, j

    call readInput()

    if (use_golden) beta = GOLDEN

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

    call save_array_1d("energies_"      // trim(outname) // ".dat", energies)
    call save_array_2d("transmissions_" // trim(outname) // ".dat", transmissions)
    call save_array_bin("transmissions_" // trim(outname) // ".bin", transmissions)

    deallocate(energies, phis, transmissions)

    contains

    subroutine readInput()
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
