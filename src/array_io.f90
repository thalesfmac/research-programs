module array_io
   use stdlib_kinds, only: dp
   implicit none
   private

   public geomspace_int

contains

   function geomspace_int(start, stp, num) result(values)
      use stdlib_error, only: check
      integer, intent(in) :: start
      integer, intent(in) :: stp
      integer, intent(in) :: num

      integer, allocatable :: values(:)

      real(dp) :: log_start, log_stop, delta_log
      integer :: i

      call check(num < 0, msg="geomspace_int: num must be non-negative")

      allocate (values(num))

      if (num == 0) return

      call check(start <= 0 .or. stp <= 0, msg="geomspace_int: start and stop must be positive")

      if (num == 1) then
         values(1) = start
         return
      end if

      log_start = log(real(start, dp))
      log_stop = log(real(stp, dp))
      delta_log = (log_stop - log_start)/real(num - 1, dp)

      do i = 1, num
         values(i) = nint(exp(log_start + real(i - 1, dp)*delta_log))
      end do

      values(1) = start
      values(num) = stp

   end function geomspace_int

end module array_io
