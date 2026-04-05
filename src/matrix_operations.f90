module matrix_operations
    use :: precision, only : dp
    implicit none

    private
    public :: identity_matrix, trace
    public :: assert_square, assert_matmul_compatibility, assert_same_shape
    public :: matmul3, matmul4

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

        if (ubound(A,1) /= ubound(A,2)) then
            error stop "identity_matrix: A must be square"
        end if
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

        if (ubound(A,1) /= ubound(A,2)) then
            error stop "identity_matrix: A must be square"
        end if
        lo = lbound(A,1)
        hi = ubound(A,1)

        retval = (0.0_dp, 0.0_dp)
        do i = lo, hi
            retval = retval + A(i, i)
        end do
    end function trace

    function matmul3(A, B, C) result(D)
        complex(dp), intent(in) :: A(:, :), B(:, :), C(:, :)
        complex(dp), allocatable :: D(:, :), tmp(:, :)

        call assert_matmul_compatibility(A, B, "A", "B")
        call assert_matmul_compatibility(B, C, "B", "C")

        tmp = matmul(B, C)
        D   = matmul(A, tmp)
    end function matmul3

    function matmul4(A, B, C, D) result(E)
        complex(dp), intent(in) :: A(:, :), B(:, :), C(:, :), D(:, :)
        complex(dp), allocatable :: E(:, :), tmp1(:, :), tmp2(:, :)

        call assert_matmul_compatibility(A, B, "A", "B")
        call assert_matmul_compatibility(B, C, "B", "C")
        call assert_matmul_compatibility(C, D, "C", "D")

        tmp1 = matmul(A, B)
        tmp2 = matmul(C, D)
        E    = matmul(tmp1, tmp2)
    end function matmul4

end module matrix_operations
