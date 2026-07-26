module transmittance
   use stdlib_kinds, only: dp
   use stdlib_linalg, only: invert

   use matrix_operations, only: matmul3, matmul4
   implicit none

   private
   public caroli_transmission

contains

   function caroli_transmission(gf, Gamma_L, Gamma_R) result(T)
      use stdlib_linalg, only: trace
      complex(dp), dimension(:, :), intent(in), contiguous :: gf, Gamma_L, Gamma_R
      real(dp) :: T

      complex(dp), allocatable :: tmp(:, :)

      allocate (tmp(size(Gamma_L, 1), size(gf, 1)))
      call matmul4(Gamma_L, gf, Gamma_R, gf, tmp, transd="C")
      T = real(trace(tmp), kind=dp)
   end function caroli_transmission

end module transmittance
