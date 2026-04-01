module matrix_operations
    use :: precision, only : dp
    implicit none
    private

    public :: identity_matrix

    contains

    subroutine identity_matrix(A)
        complex(dp), intent(out) :: A(:,:)
        integer :: lo, hi, i

        if (ubound(A,1) /= ubound(A,2)) then
            error stop "identity_matrix: I must be square"
        end if

        lo = lbound(A,1)
        hi = ubound(A,1)

        A = (0.0_dp, 0.0_dp)
        do i = lo, hi
            A(i, i) = (1.0_dp, 0.0_dp)
        end do
    end subroutine identity_matrix

end module matrix_operations
