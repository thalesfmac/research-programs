module lapack_blas_interface
   use, intrinsic :: iso_fortran_env, only: real64
   use stdlib_kinds, only: dp
   implicit none

   public

   interface
      subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
         import :: real64
         character(len=1), intent(in) :: transa, transb
         integer, intent(in) :: m, n, k, lda, ldb, ldc
         complex(real64), intent(in) :: alpha, beta
         complex(real64), intent(in) :: a(lda, *), b(ldb, *)
         complex(real64), intent(inout) :: c(ldc, *)
      end subroutine zgemm

      subroutine zgetrf(m, n, a, lda, ipiv, info)
         import :: real64
         integer, intent(in) :: m, n, lda
         integer, intent(out) :: ipiv(*), info
         complex(real64), intent(inout) :: a(lda, *)
      end subroutine zgetrf

      subroutine zgetri(n, a, lda, ipiv, work, lwork, info)
         import :: real64
         integer, intent(in) :: n, lda, lwork
         integer, intent(in) :: ipiv(*)
         integer, intent(out) :: info
         complex(real64), intent(inout) :: a(lda, *), work(*)
      end subroutine zgetri
   end interface

contains

   subroutine diagonalize_d(A, w, jobz, uplo)
      use stdlib_error, only: check
      use stdlib_linalg, only: is_square
      use stdlib_linalg_lapack, only: ilp, heevd
      ! Diagonaliza matriz Hermitiana complexa via ZHEEVD.
      complex(dp), intent(inout), contiguous :: A(:, :)
      real(dp), intent(out) :: w(:)
      character(len=1), intent(in), optional :: jobz, uplo

      character(len=1) :: jobz_loc, uplo_loc
      integer(ilp) :: n, lda, info
      integer(ilp) :: lwork, lrwork, liwork

      complex(dp), allocatable :: work(:)
      real(dp), allocatable :: rwork(:)
      integer(ilp), allocatable :: iwork(:)

      complex(dp) :: workq(1)
      real(dp) :: rworkq(1)
      integer(ilp) :: iworkq(1)

      jobz_loc = 'V'
      if (present(jobz)) jobz_loc = jobz

      uplo_loc = 'U'
      if (present(uplo)) uplo_loc = uplo

      call check(is_square(A), msg="diagonalize_d: A must be a square matrix")

      n = size(A, 1, kind=ilp)

      call check(size(w, kind=ilp) == n, msg="diagonalize_d: w must have length n of A")

      lda = max(1_ilp, n)

      ! Consulta dos tamanhos ótimos dos workspaces.
      lwork = -1_ilp
      lrwork = -1_ilp
      liwork = -1_ilp

      call heevd(jobz_loc, uplo_loc, n, A, lda, w, workq, lwork, rworkq, lrwork, iworkq, liwork, info)

      call check(info == 0_ilp, msg="diagonalize_d: ZHEEVD workspace query failed")

      lwork = max(1_ilp, int(real(workq(1), kind=dp), kind=ilp))
      lrwork = max(1_ilp, int(rworkq(1), kind=ilp))
      liwork = max(1_ilp, iworkq(1))

      allocate (work(lwork))
      allocate (rwork(lrwork))
      allocate (iwork(liwork))

      ! Diagonalização.
      call heevd(jobz_loc, uplo_loc, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)

      deallocate (work, rwork, iwork)

      if (info < 0) then
         error stop "diagonalize_d: ZHEEVD received an invalid argument"
      else if (info > 0) then
         error stop "diagonalize_d: ZHEEVD failed to converge"
      end if
   end subroutine diagonalize_d

end module lapack_blas_interface
