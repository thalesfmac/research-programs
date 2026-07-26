module lapack_blas_interface
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none

   public

   interface
      subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
         import :: real64
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, lda, lwork
         integer, intent(out) :: info
         complex(real64), intent(inout) :: a(lda, *)
         real(real64), intent(out) :: w(*)
         complex(real64), intent(inout) :: work(*)
         real(real64), intent(inout) :: rwork(*)
      end subroutine zheev

      subroutine zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
         import :: real64
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, lda
         integer, intent(in) :: lwork, lrwork, liwork
         integer, intent(out) :: info
         complex(real64), intent(inout) :: a(lda, *)
         real(real64), intent(out) :: w(*)
         complex(real64), intent(inout) :: work(*)
         real(real64), intent(inout) :: rwork(*)
         integer, intent(inout) :: iwork(*)
      end subroutine zheevd

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
      ! Diagonaliza matriz Hermitiana complexa via ZHEEVD.
      complex(real64), intent(inout), contiguous :: A(:, :)
      real(real64), intent(out) :: w(:)
      character(len=1), intent(in), optional :: jobz, uplo

      character(len=1) :: jobz_loc, uplo_loc
      integer :: n, lda, info
      integer :: lwork, lrwork, liwork

      complex(real64), allocatable :: work(:)
      real(real64), allocatable :: rwork(:)
      integer, allocatable :: iwork(:)

      complex(real64) :: workq(1)
      real(real64) :: rworkq(1)
      integer :: iworkq(1)

      jobz_loc = 'V'
      if (present(jobz)) jobz_loc = jobz

      uplo_loc = 'U'
      if (present(uplo)) uplo_loc = uplo

      call check(is_square(A), msg="diagonalize_d: A must be a square matrix")

      n = size(A, 1)

      if (size(w) /= n) then
         error stop "diagonalize_d: w must have length n of A"
      end if

      lda = max(1, n)

      ! Consulta dos tamanhos ótimos dos workspaces.
      lwork = -1
      lrwork = -1
      liwork = -1

      call zheevd(jobz_loc, uplo_loc, n, A, lda, w, workq, lwork, rworkq, lrwork, iworkq, liwork, info)

      if (info /= 0) then
         error stop "diagonalize_d: ZHEEVD workspace query failed"
      end if

      lwork = max(1, int(real(workq(1), kind=real64)))
      lrwork = max(1, int(rworkq(1)))
      liwork = max(1, iworkq(1))

      allocate (work(lwork))
      allocate (rwork(lrwork))
      allocate (iwork(liwork))

      ! Diagonalização.
      call zheevd(jobz_loc, uplo_loc, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)

      deallocate (work, rwork, iwork)

      if (info < 0) then
         error stop "diagonalize_d: ZHEEVD received an invalid argument"
      else if (info > 0) then
         error stop "diagonalize_d: ZHEEVD failed to converge"
      end if
   end subroutine diagonalize_d

end module lapack_blas_interface
