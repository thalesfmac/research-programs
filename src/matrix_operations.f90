module matrix_operations
    use :: precision, only : dp
    use :: lapack_blas, only : gemm
    implicit none

    private
    public :: identity_matrix, trace
    public :: assert_square, assert_matmul_compatibility, assert_same_shape
    ! public :: matmul3, matmul4

    interface assert_square
        module procedure assert_square_cdp
        module procedure assert_square_rdp
    end interface

    interface assert_matmul_compatibility
        module procedure assert_matmul_compatibility_cdp
        module procedure assert_matmul_compatibility_rdp
    end interface

    interface assert_same_shape
        module procedure assert_same_shape_cdp
        module procedure assert_same_shape_rdp
    end interface

    contains

    subroutine assert_square_cdp(A, name)
        complex(dp), intent(in) :: A(:, :)
        character(len=*), intent(in) :: name

        if (size(A,1) /= size(A,2)) then
            error stop trim(name) // " must be square"
        end if
    end subroutine assert_square_cdp

    subroutine assert_square_rdp(A, name)
        real(dp), intent(in) :: A(:, :)
        character(len=*), intent(in) :: name

        if (size(A,1) /= size(A,2)) then
            error stop trim(name) // " must be square"
        end if
    end subroutine assert_square_rdp

    subroutine assert_matmul_compatibility_cdp(A, B, nameA, nameB)
        complex(dp), intent(in) :: A(:, :), B(:, :)
        character(len=*), intent(in) :: nameA, nameB

        character(len=32) :: a1, a2, b1, b2

        if (size(A,2) /= size(B,1)) then
            write(a1, '(I0)') size(A,1)
            write(a2, '(I0)') size(A,2)
            write(b1, '(I0)') size(B,1)
            write(b2, '(I0)') size(B,2)

            error stop trim(nameA) // " (" // trim(a1) // "," // trim(a2) // &
                ") and " // trim(nameB) // " (" // trim(b1) // "," // trim(b2) // &
                ") are not compatible for matmul"
        end if
    end subroutine assert_matmul_compatibility_cdp

    subroutine assert_matmul_compatibility_rdp(A, B, nameA, nameB)
        real(dp), intent(in) :: A(:, :), B(:, :)
        character(len=*), intent(in) :: nameA, nameB

        character(len=32) :: a1, a2, b1, b2

        if (size(A,2) /= size(B,1)) then
            write(a1, '(I0)') size(A,1)
            write(a2, '(I0)') size(A,2)
            write(b1, '(I0)') size(B,1)
            write(b2, '(I0)') size(B,2)

            error stop trim(nameA) // " (" // trim(a1) // "," // trim(a2) // &
                ") and " // trim(nameB) // " (" // trim(b1) // "," // trim(b2) // &
                ") are not compatible for matmul"
        end if
    end subroutine assert_matmul_compatibility_rdp

    subroutine assert_same_shape_cdp(A, B, nameA, nameB)
        complex(dp), intent(in) :: A(:, :), B(:, :)
        character(len=*), intent(in) :: nameA, nameB
        character(len=32) :: a1, a2, b1, b2

        if (size(A,1) /= size(B,1) .or. size(A,2) /= size(B,2)) then
            write(a1, '(I0)') size(A,1)
            write(a2, '(I0)') size(A,2)
            write(b1, '(I0)') size(B,1)
            write(b2, '(I0)') size(B,2)

            error stop trim(nameA) // " (" // trim(a1) // "," // trim(a2) // &
                ") and " // trim(nameB) // " (" // trim(b1) // "," // trim(b2) // &
                ") do not have the same shape"
        end if
    end subroutine assert_same_shape_cdp

    subroutine assert_same_shape_rdp(A, B, nameA, nameB)
        real(dp), intent(in) :: A(:, :), B(:, :)
        character(len=*), intent(in) :: nameA, nameB
        character(len=32) :: a1, a2, b1, b2

        if (size(A,1) /= size(B,1) .or. size(A,2) /= size(B,2)) then
            write(a1, '(I0)') size(A,1)
            write(a2, '(I0)') size(A,2)
            write(b1, '(I0)') size(B,1)
            write(b2, '(I0)') size(B,2)

            error stop trim(nameA) // " (" // trim(a1) // "," // trim(a2) // &
                ") and " // trim(nameB) // " (" // trim(b1) // "," // trim(b2) // &
                ") do not have the same shape"
        end if
    end subroutine assert_same_shape_rdp

    subroutine identity_matrix(A)
        complex(dp), intent(out) :: A(:,:)
        integer :: lo, hi, i

        call assert_square(A, "A")
        lo = lbound(A,1)
        hi = ubound(A,1)

        A = (0.0_dp, 0.0_dp)
        do i = lo, hi
            A(i, i) = (1.0_dp, 0.0_dp)
        end do
    end subroutine identity_matrix

    function trace(A) result(retval)
        complex(dp), intent(in) :: A(:,:)
        complex(dp) :: retval
        integer :: lo, hi, i

        call assert_square(A, "A")
        lo = lbound(A,1)
        hi = ubound(A,1)

        retval = (0.0_dp, 0.0_dp)
        do i = lo, hi
            retval = retval + A(i, i)
        end do
    end function trace

    subroutine matmul3(A, B, C, D, transa, transb, transc)
        complex(dp), intent(in)  :: A(:, :), B(:, :), C(:, :)
        complex(dp), intent(out) :: D(:, :)
        character(len=1), intent(in), optional :: transa, transb, transc

        character(len=1) :: ta, tb, tc
        complex(dp), allocatable :: tmp(:, :)

        ta = 'N'; if (present(transa)) ta = transa
        tb = 'N'; if (present(transb)) tb = transb
        tc = 'N'; if (present(transc)) tb = transc

        ! call assert_matmul_compatibility(A, B, "A", "B")
        ! call assert_matmul_compatibility(B, C, "B", "C")

        if (size(D,1) /= size(A,1) .or. size(D,2) /= size(C,2)) then
            error stop "matmul3: D has incompatible shape"
        end if

        tmp = matmul(B, C)
        D   = matmul(A, tmp)
    end subroutine matmul3

    ! subroutine matmul4(A, B, C, D, E)
    !     complex(dp), intent(in)  :: A(:, :), B(:, :), C(:, :), D(:, :)
    !     complex(dp), intent(out) :: E(:, :)
    !     complex(dp), allocatable :: tmp1(:, :), tmp2(:, :)

    !     call assert_matmul_compatibility(A, B, "A", "B")
    !     call assert_matmul_compatibility(B, C, "B", "C")
    !     call assert_matmul_compatibility(C, D, "C", "D")

    !     if (size(E,1) /= size(A,1) .or. size(E,2) /= size(D,2)) then
    !         error stop "matmul4: E has incompatible shape"
    !     end if

    !     tmp1 = matmul(A, B)
    !     tmp2 = matmul(C, D)
    !     E    = matmul(tmp1, tmp2)
    ! end subroutine matmul4

end module matrix_operations
