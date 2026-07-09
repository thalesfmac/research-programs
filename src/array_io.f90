module array_io
   use precision, only: dp
   implicit none

   private
   public :: save_array_1d, save_array_2d, save_array_bin
   public :: arange_int, geomspace_int

contains
   subroutine save_array_1d(filename, A)
      character(len=*), intent(in) :: filename
      real(dp), intent(in) :: A(:)

      integer :: u, i

      open (newunit=u, file=filename, status="replace", action="write")

      do i = 1, size(A)
         ! write(u, '(ES24.16)') A(i)
         write (u, '(E26.16E3)') A(i)
      end do

      close (u)
   end subroutine save_array_1d

   subroutine save_array_2d(filename, A)
      character(len=*), intent(in) :: filename
      real(dp), intent(in) :: A(:, :)

      integer :: u, i, j, ncol

      open (newunit=u, file=filename, status="replace", action="write")

      ncol = size(A, 2)

      do i = 1, size(A, 1)
         do j = 1, ncol
            ! write(u, '(ES24.16)', advance='no') A(i,j)
            write (u, '(E26.16E3)', advance='no') A(i, j)
            if (j < ncol) write (u, '(A)', advance='no') ' '
         end do
         write (u, *)
      end do

      close (u)
   end subroutine save_array_2d

   subroutine save_array_bin(filename, A)
      use, intrinsic :: iso_fortran_env, only: int32
      character(len=*), intent(in) :: filename
      real(dp), intent(in) :: A(..)

      integer :: u
      integer(int32) :: r
      integer(int32), allocatable :: shp(:)

      open (newunit=u, file=filename, status="replace", action="write", &
            form="unformatted", access="stream")

      select rank (A)

      rank (1)
         r = 1_int32
         allocate (shp(1))
         shp = int(shape(A), kind=int32)

         write (u) r
         write (u) shp
         write (u) A

      rank (2)
         r = 2_int32
         allocate (shp(2))
         shp = int(shape(A), kind=int32)

         write (u) r
         write (u) shp
         write (u) A

      rank (3)
         r = 3_int32
         allocate (shp(3))
         shp = int(shape(A), kind=int32)

         write (u) r
         write (u) shp
         write (u) A

      rank default
         close (u)
         error stop "save_array_bin: rank nao suportado"
      end select

      close (u)

      if (allocated(shp)) deallocate (shp)
   end subroutine save_array_bin

   function arange_int(start, stop, step) result(values)
      integer, intent(in) :: start, stop, step
      integer, allocatable :: values(:)

      integer :: n, i

      if (step == 0) then
         error stop "arange_int: step cannot be zero"
      end if

      if ((stop - start)*step <= 0) then
         allocate (values(0))
         return
      end if

      n = (stop - start + step - sign(1, step))/step

      allocate (values(n))

      do i = 1, n
         values(i) = start + (i - 1)*step
      end do

   end function arange_int

   function geomspace_int(start, stop, num) result(values)
      integer, intent(in) :: start
      integer, intent(in) :: stop
      integer, intent(in) :: num

      integer, allocatable :: values(:)

      real(dp) :: log_start, log_stop, dlog
      integer :: i

      if (num < 0) then
         error stop "geomspace_int: num must be non-negative"
      end if

      allocate (values(num))

      if (num == 0) return

      if (start <= 0 .or. stop <= 0) then
         error stop "geomspace_int: start and stop must be positive"
      end if

      if (num == 1) then
         values(1) = start
         return
      end if

      log_start = log(real(start, dp))
      log_stop = log(real(stop, dp))

      dlog = (log_stop - log_start)/real(num - 1, dp)

      do i = 1, num
         values(i) = nint(exp(log_start + real(i - 1, dp)*dlog))
      end do

      values(1) = start
      values(num) = stop

   end function geomspace_int

end module array_io
