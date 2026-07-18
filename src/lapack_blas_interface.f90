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

end module lapack_blas_interface
