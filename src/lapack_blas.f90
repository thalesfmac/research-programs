module lapack_blas
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    private
    public :: diagonalize, invert
    public :: matmul2, matmul3, matmul4

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
        if (size(A,2) /= n) error stop "diagonalize: A must be square"
        if (size(w)   /= n) error stop "diagonalize: w must have length n of A"
        lda = n

        allocate(rwork(max(1, 3*n - 2)))

        ! Workspace query
        lwork = -1
        call zheev(jobz_loc, uplo_loc, n, A, lda, w, workq, lwork, rwork, info)
        if (info /= 0) then
            deallocate(rwork)
            error stop "diagonalize: ZHEEV workspace query failed"
        end if

        lwork = int(real(workq(1), real64))
        if (lwork < 1) lwork = 1
        allocate(work(lwork))

        ! Diagonalização
        call zheev(jobz_loc, uplo_loc, n, A, lda, w, work, lwork, rwork, info)

        deallocate(work, rwork)
        if (info /= 0) error stop "diagonalize: ZHEEV failed"
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
            error stop "gemm: transa must be 'N','T' or 'C'"
        end select

        select case (tb)
          case ('N','n')
            b_rows = size(B,1); b_cols = size(B,2)
          case ('T','t','C','c')
            b_rows = size(B,2); b_cols = size(B,1)
          case default
            error stop "gemm: transb must be 'N','T' or 'C'"
        end select

        m = a_rows
        k = a_cols
        n = b_cols

        if (b_rows /= k) error stop "gemm: A and B have imcompatible dimensions"
        if (size(C,1) /= m .or. size(C,2) /= n) error stop "gemm: C have imcompatible dimensions with A*B"

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
        if (size(A,2) /= n) error stop "invert: A must be square"

        allocate(ipiv(n))

        call zgetrf(n,n,A,n,ipiv,info)
        if (info /= 0) error stop "invert: ZGETRF failed"

        ! workspace query
        lwork = -1
        allocate(work(1))
        call zgetri(n,A,n,ipiv,work,lwork,info)
        if (info /= 0) error stop "invert: ZGETRI workspace query failed"
        lwork = int(real(work(1), real64))
        deallocate(work)

        allocate(work(max(1,lwork)))
        call zgetri(n,A,n,ipiv,work,lwork,info)
        if (info /= 0) error stop "invert: ZGETRI failed"

        deallocate(work, ipiv)
    end subroutine invert

    subroutine op_shape(X, trans, nrow, ncol)
        complex(real64), intent(in) :: X(:,:)
        character(len=1), intent(in) :: trans
        integer, intent(out) :: nrow, ncol

        select case (trans)
          case ('N','n')
            nrow = size(X,1)
            ncol = size(X,2)
          case ('T','t','C','c')
            nrow = size(X,2)
            ncol = size(X,1)
          case default
            error stop "op_shape: must be 'N', 'T' or 'C'"
        end select
    end subroutine op_shape

    subroutine matmul2(A, B, C, transa, transb)
        complex(real64), intent(in), contiguous :: A(:,:), B(:,:)
        complex(real64), intent(out), contiguous :: C(:,:)
        character(len=1), intent(in), optional :: transa, transb

        integer a_rows, a_cols, b_rows, b_cols
        character(len=1) :: ta, tb
        integer :: lda, ldb, ldc

        ta = 'N'; if (present(transa)) ta = transa
        tb = 'N'; if (present(transb)) tb = transb

        call op_shape(A, ta, a_rows, a_cols)
        call op_shape(B, tb, b_rows, b_cols)

        if (a_cols /= b_rows) error stop "matmul2: op(A) and op(B) have imcompatible dimensions"
        if (size(C, 1) /= a_rows .or. size(C, 2) /= b_cols) then
            error stop "matmul2: C has incompatible dimensions for the result"
        end if

        C = (0.0_real64, 0.0_real64)

        lda = size(A,1)
        ldb = size(B,1)
        ldc = size(C,1)

        call zgemm( &
            ta, &
            tb, &
            a_rows,  &
            b_cols,  &
            a_cols,  &
            (1.0_real64, 0.0_real64), &
            A, &
            lda, &
            B, &
            ldb, &
            (0.0_real64, 0.0_real64), &
            C, &
            ldc &
            )
    end subroutine matmul2

    subroutine matmul3(A, B, C, D, transa, transb, transc)
        complex(real64), intent(in), contiguous :: A(:,:), B(:,:), C(:,:)
        complex(real64), intent(out), contiguous :: D(:,:)
        character(len=1), intent(in), optional :: transa, transb, transc

        character(len=1) :: ta, tb, tc
        integer :: a_rows, a_cols
        integer :: b_rows, b_cols
        integer :: c_rows, c_cols
        complex(real64), allocatable :: T(:,:)

        ta = 'N'; if (present(transa)) ta = transa
        tb = 'N'; if (present(transb)) tb = transb
        tc = 'N'; if (present(transc)) tc = transc

        call op_shape(A, ta, a_rows, a_cols)
        call op_shape(B, tb, b_rows, b_cols)
        call op_shape(C, tc, c_rows, c_cols)

        if (a_cols /= b_rows) error stop "matmul3: op(A) and op(B) have incompatible dimensions"
        if (b_cols /= c_rows) error stop "matmul3: op(B) and op(C) have incompatible dimensions"
        if (size(D,1) /= a_rows .or. size(D,2) /= c_cols) then
            error stop "matmul3: D has incompatible dimensions for the result"
        end if

        allocate(T(a_rows, b_cols))

        call matmul2(A, B, T, ta, tb)
        call matmul2(T, C, D, 'N', tc)
    end subroutine matmul3

    subroutine matmul4(A, B, C, D, E, transa, transb, transc, transd)
        use iso_fortran_env, only: real64
        implicit none

        complex(real64), intent(in),  contiguous :: A(:,:), B(:,:), C(:,:), D(:,:)
        complex(real64), intent(out), contiguous :: E(:,:)
        character(len=1), intent(in), optional   :: transa, transb, transc, transd

        character(len=1) :: ta, tb, tc, td
        integer :: a_rows, a_cols
        integer :: b_rows, b_cols
        integer :: c_rows, c_cols
        integer :: d_rows, d_cols
        complex(real64), allocatable :: T1(:,:), T2(:,:)

        ta = 'N'; if (present(transa)) ta = transa
        tb = 'N'; if (present(transb)) tb = transb
        tc = 'N'; if (present(transc)) tc = transc
        td = 'N'; if (present(transd)) td = transd

        ! op(A) = (a_rows, a_cols)
        call op_shape(A, ta, a_rows, a_cols)

        ! op(B) = (b_rows, b_cols)
        call op_shape(B, tb, b_rows, b_cols)

        ! op(C) = (c_rows, c_cols)
        call op_shape(C, tc, c_rows, c_cols)

        ! op(D) = (d_rows, d_cols)
        call op_shape(D, td, d_rows, d_cols)

        ! Verificações de compatibilidade:
        ! AB
        if (a_cols /= b_rows) then
            error stop "matmul4: op(A) and op(B) have incompatible dimensions"
        end if

        ! CD
        if (c_cols /= d_rows) then
            error stop "matmul4: op(C) and op(D) have incompatible dimensions"
        end if

        ! (AB)(CD)
        if (b_cols /= c_rows) then
            error stop "matmul4: op(A)op(B) and op(C)op(D) have incompatible dimensions"
        end if

        ! Shape final de E
        if (size(E,1) /= a_rows .or. size(E,2) /= d_cols) then
            error stop "matmul4: E has incompatible dimensions for the result"
        end if

        allocate(T1(a_rows, b_cols))
        allocate(T2(c_rows, d_cols))

        ! T1 = op(A) * op(B)
        call matmul2(A, B, T1, ta, tb)

        ! T2 = op(C) * op(D)
        call matmul2(C, D, T2, tc, td)

        ! E = T1 * T2
        call matmul2(T1, T2, E, 'N', 'N')

    end subroutine matmul4

end module lapack_blas
