module rng_utils
   implicit none
   private
   public :: rng_initialize

contains

   subroutine rng_initialize(seed)
      integer, intent(in) :: seed
      integer :: n, j
      integer, allocatable :: seed_vec(:)

      call random_seed(size=n)
      allocate (seed_vec(n))

      seed_vec = seed + 37*[(j - 1, j=1, n)]

      call random_seed(put=seed_vec)
      deallocate (seed_vec)
   end subroutine rng_initialize

end module rng_utils
