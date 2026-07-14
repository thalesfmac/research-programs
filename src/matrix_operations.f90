module matrix_operations
   use :: precision, only:dp
   use :: lapack_blas_interface
   implicit none

   private
   public :: identity_matrix, trace
   public :: assert_square, assert_same_shape
   public :: diagonalize, invert
   public :: matmul2, matmul3, matmul4

   interface assert_square
      module procedure assert_square_cdp
      module procedure assert_square_rdp
   end interface

   interface assert_same_shape
      module procedure assert_same_shape_cdp
      module procedure assert_same_shape_rdp
   end interface

contains

   pure function with_caller(msg, caller) result(full_msg)
      character(len=*), intent(in) :: msg
      character(len=*), intent(in), optional :: caller
      character(len=:), allocatable :: full_msg

      if (present(caller)) then
         full_msg = trim(caller)//": "//trim(msg)
      else
         full_msg = trim(msg)
      end if
   end function with_caller

   subroutine shape_to_strings(nrow, ncol, srow, scol)
      integer, intent(in) :: nrow, ncol
      character(len=32), intent(out) :: srow, scol

      write (srow, '(I0)') nrow
      write (scol, '(I0)') ncol
   end subroutine shape_to_strings

   subroutine assert_square_cdp(A, name, caller)
      complex(dp), intent(in) :: A(:, :)
      character(len=*), intent(in) :: name
      character(len=*), intent(in), optional :: caller

      character(len=32) :: n1, n2

      if (size(A, 1) /= size(A, 2)) then
         call shape_to_strings(size(A, 1), size(A, 2), n1, n2)
         error stop with_caller( &
            trim(name)//" must be square, got shape ("// &
            trim(n1)//","//trim(n2)//")", caller &
            )
      end if
   end subroutine assert_square_cdp

   subroutine assert_square_rdp(A, name, caller)
      real(dp), intent(in) :: A(:, :)
      character(len=*), intent(in) :: name
      character(len=*), intent(in), optional :: caller

      character(len=32) :: n1, n2

      if (size(A, 1) /= size(A, 2)) then
         call shape_to_strings(size(A, 1), size(A, 2), n1, n2)
         error stop with_caller( &
            trim(name)//" must be square, got shape ("// &
            trim(n1)//","//trim(n2)//")", caller &
            )
      end if
   end subroutine assert_square_rdp

   subroutine assert_same_shape_cdp(A, B, nameA, nameB, caller)
      complex(dp), intent(in) :: A(:, :), B(:, :)
      character(len=*), intent(in) :: nameA, nameB
      character(len=*), intent(in), optional :: caller

      character(len=32) :: a1, a2, b1, b2

      if (size(A, 1) /= size(B, 1) .or. size(A, 2) /= size(B, 2)) then
         call shape_to_strings(size(A, 1), size(A, 2), a1, a2)
         call shape_to_strings(size(B, 1), size(B, 2), b1, b2)

         error stop with_caller( &
            trim(nameA)//" has shape ("//trim(a1)//","//trim(a2)// &
            "), "//trim(nameB)//" has shape ("//trim(b1)//","//trim(b2)// &
            "), expected same shape", caller)
      end if
   end subroutine assert_same_shape_cdp

   subroutine assert_same_shape_rdp(A, B, nameA, nameB, caller)
      real(dp), intent(in) :: A(:, :), B(:, :)
      character(len=*), intent(in) :: nameA, nameB
      character(len=*), intent(in), optional :: caller

      character(len=32) :: a1, a2, b1, b2

      if (size(A, 1) /= size(B, 1) .or. size(A, 2) /= size(B, 2)) then
         call shape_to_strings(size(A, 1), size(A, 2), a1, a2)
         call shape_to_strings(size(B, 1), size(B, 2), b1, b2)

         error stop with_caller( &
            trim(nameA)//" has shape ("//trim(a1)//","//trim(a2)// &
            "), "//trim(nameB)//" has shape ("//trim(b1)//","//trim(b2)// &
            "), expected same shape", caller)
      end if
   end subroutine assert_same_shape_rdp

   subroutine identity_matrix(A)
      complex(dp), intent(out), contiguous :: A(:, :)
      integer :: i

      call assert_square(A, "A", "identity_matrix")

      A = (0.0_dp, 0.0_dp)
      do i = 1, size(A, dim=1)
         A(i, i) = (1.0_dp, 0.0_dp)
      end do
   end subroutine identity_matrix

   function trace(A) result(retval)
      complex(dp), intent(in) :: A(:, :)
      complex(dp) :: retval
      integer :: i

      call assert_square(A, "A", "trace")

      retval = (0.0_dp, 0.0_dp)
      do i = 1, size(A, dim=1)
         retval = retval + A(i, i)
      end do
   end function trace

   subroutine diagonalize(A, w, jobz, uplo)
      ! Diagonaliza Hermitiana complexa via ZHEEV.
      complex(dp), intent(inout), contiguous :: A(:, :)
      real(dp), intent(out) :: w(:)
      character(len=1), intent(in), optional :: jobz, uplo

      character(len=1) :: jobz_loc, uplo_loc
      integer :: n, lda, info, lwork
      complex(dp), allocatable :: work(:)
      real(dp), allocatable :: rwork(:)
      complex(dp) :: workq(1)

      jobz_loc = 'V'; if (present(jobz)) jobz_loc = jobz
      uplo_loc = 'U'; if (present(uplo)) uplo_loc = uplo

      call assert_square(A, "A", caller="diagonalize")

      n = size(A, 1)
      if (size(w) /= n) error stop "diagonalize: w must have length n of A"
      lda = n

      allocate (rwork(max(1, 3*n - 2)))

      ! Workspace query
      lwork = -1
      call zheev(jobz_loc, uplo_loc, n, A, lda, w, workq, lwork, rwork, info)
      if (info /= 0) then
         deallocate (rwork)
         error stop "diagonalize: ZHEEV workspace query failed"
      end if

      lwork = int(real(workq(1), dp))
      if (lwork < 1) lwork = 1
      allocate (work(lwork))

      ! Diagonalização
      call zheev(jobz_loc, uplo_loc, n, A, lda, w, work, lwork, rwork, info)

      deallocate (work, rwork)
      if (info /= 0) error stop "diagonalize: ZHEEV failed"
   end subroutine diagonalize

   subroutine invert(A)
      complex(dp), intent(inout), contiguous :: A(:, :)
      integer :: n, info, lwork
      integer, allocatable :: ipiv(:)
      complex(dp), allocatable :: work(:)

      call assert_square(A, "A", caller="invert")

      n = size(A, 1)
      if (size(A, 2) /= n) error stop "invert: A must be square"

      allocate (ipiv(n))

      call zgetrf(n, n, A, n, ipiv, info)
      if (info /= 0) error stop "invert: ZGETRF failed"

      ! workspace query
      lwork = -1
      allocate (work(1))
      call zgetri(n, A, n, ipiv, work, lwork, info)
      if (info /= 0) error stop "invert: ZGETRI workspace query failed"
      lwork = int(real(work(1), dp))
      deallocate (work)

      allocate (work(max(1, lwork)))
      call zgetri(n, A, n, ipiv, work, lwork, info)
      if (info /= 0) error stop "invert: ZGETRI failed"

      deallocate (work, ipiv)
   end subroutine invert

   subroutine op_shape(X, trans, nrow, ncol)
      complex(dp), intent(in) :: X(:, :)
      character(len=1), intent(in) :: trans
      integer, intent(out) :: nrow, ncol

      select case (trans)
      case ('N', 'n')
         nrow = size(X, 1)
         ncol = size(X, 2)
      case ('T', 't', 'C', 'c')
         nrow = size(X, 2)
         ncol = size(X, 1)
      case default
         error stop "op_shape: must be 'N', 'T' or 'C'"
      end select
   end subroutine op_shape

   subroutine matmul2(A, B, C, transa, transb)
      complex(dp), intent(in), contiguous :: A(:, :), B(:, :)
      complex(dp), intent(out), contiguous :: C(:, :)
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

      C = (0.0_dp, 0.0_dp)

      lda = size(A, 1)
      ldb = size(B, 1)
      ldc = size(C, 1)

      call zgemm( &
         ta, &
         tb, &
         a_rows, &
         b_cols, &
         a_cols, &
         (1.0_dp, 0.0_dp), &
         A, &
         lda, &
         B, &
         ldb, &
         (0.0_dp, 0.0_dp), &
         C, &
         ldc &
         )
   end subroutine matmul2

   subroutine matmul3(A, B, C, D, transa, transb, transc)
      complex(dp), intent(in), contiguous :: A(:, :), B(:, :), C(:, :)
      complex(dp), intent(out), contiguous :: D(:, :)
      character(len=1), intent(in), optional :: transa, transb, transc

      character(len=1) :: ta, tb, tc
      integer :: a_rows, a_cols
      integer :: b_rows, b_cols
      integer :: c_rows, c_cols
      complex(dp), allocatable :: T(:, :)

      ta = 'N'; if (present(transa)) ta = transa
      tb = 'N'; if (present(transb)) tb = transb
      tc = 'N'; if (present(transc)) tc = transc

      call op_shape(A, ta, a_rows, a_cols)
      call op_shape(B, tb, b_rows, b_cols)
      call op_shape(C, tc, c_rows, c_cols)

      if (a_cols /= b_rows) error stop "matmul3: op(A) and op(B) have incompatible dimensions"
      if (b_cols /= c_rows) error stop "matmul3: op(B) and op(C) have incompatible dimensions"
      if (size(D, 1) /= a_rows .or. size(D, 2) /= c_cols) then
         error stop "matmul3: D has incompatible dimensions for the result"
      end if

      allocate (T(a_rows, b_cols))

      call matmul2(A, B, T, ta, tb)
      call matmul2(T, C, D, 'N', tc)
   end subroutine matmul3

   subroutine matmul4(A, B, C, D, E, transa, transb, transc, transd)
      complex(dp), intent(in), contiguous :: A(:, :), B(:, :), C(:, :), D(:, :)
      complex(dp), intent(out), contiguous :: E(:, :)
      character(len=1), intent(in), optional :: transa, transb, transc, transd

      character(len=1) :: ta, tb, tc, td
      integer :: a_rows, a_cols
      integer :: b_rows, b_cols
      integer :: c_rows, c_cols
      integer :: d_rows, d_cols
      complex(dp), allocatable :: T1(:, :), T2(:, :)

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
      if (size(E, 1) /= a_rows .or. size(E, 2) /= d_cols) then
         error stop "matmul4: E has incompatible dimensions for the result"
      end if

      allocate (T1(a_rows, b_cols))
      allocate (T2(c_rows, d_cols))

      ! T1 = op(A) * op(B)
      call matmul2(A, B, T1, ta, tb)

      ! T2 = op(C) * op(D)
      call matmul2(C, D, T2, tc, td)

      ! E = T1 * T2
      call matmul2(T1, T2, E, 'N', 'N')

   end subroutine matmul4

end module matrix_operations
