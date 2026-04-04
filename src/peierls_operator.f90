module peierls_operator
    use :: precision, only : dp
    use :: constants, only : CI
    use :: matrix_operations

    implicit none

    private
    public :: peierls_exp

    contains

    subroutine build_factorials(fact)
        real(dp), intent(out) :: fact(0:)
        integer :: k, Nph

        Nph = ubound(fact, dim=1)

        fact(0) = 1.0_dp
        do k = 1, Nph
            fact(k) = fact(k-1) * real(k, dp)
        end do
    end subroutine build_factorials

    ! P(0,M)=1 e P(j,M)=sqrt(M-(j-1))*P(j-1,M)  => P(j,M)=sqrt(M!/(M-j)!)
    subroutine build_P(P)
        real(dp), intent(out) :: P(0:, 0:)
        integer :: M, j, Nph

        call assert_square(P, "P")

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
        real(dp) :: prefactor

        call assert_square(A, "A")

        Nph = ubound(A,1)

        allocate( fact(0:Nph) )
        allocate( P(0:Nph, 0:Nph) )

        call build_factorials(fact)
        call build_P(P)

        prefactor = exp(-0.5_dp * (g*g))
        A = (0.0_dp, 0.0_dp)

        ! índices físicos N,M = 0..Nph  (no array: +1)
        do N = 0, Nph
            do M = 0, Nph
                hNM = (0.0_dp, 0.0_dp)

                do s = 0, N
                    do j = 0, M
                        ! delta(N-s, M-j) == 1  <=>  N - s == M - j
                        if (N - s == M - j) then
                            hNM = hNM + prefactor * ( (CI*g)**s ) * ( (CI*g)**j ) &
                                * P(s, N) * P(j, M) / ( fact(s) * fact(j) )
                        end if
                    end do
                end do

                A(N, M) = hNM
            end do
        end do
    end subroutine peierls_exp

end module peierls_operator
