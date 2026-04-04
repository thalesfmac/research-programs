module array_io
    use precision, only : dp
    implicit none

    private
    public :: save_array_1d, save_array_2d, save_array_bin

    contains
    subroutine save_array_1d(filename, A)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: A(:)

        integer :: u, i

        open(newunit=u, file=filename, status="replace", action="write")

        do i = 1, size(A)
            ! write(u, '(ES24.16)') A(i)
            write(u, '(E26.16E3)') A(i)
        end do

        close(u)
    end subroutine save_array_1d

    subroutine save_array_2d(filename, A)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: A(:,:)

        integer :: u, i, j, ncol

        open(newunit=u, file=filename, status="replace", action="write")

        ncol = size(A, 2)

        do i = 1, size(A, 1)
            do j = 1, ncol
                ! write(u, '(ES24.16)', advance='no') A(i,j)
                write(u, '(E26.16E3)', advance='no') A(i,j)
                if (j < ncol) write(u, '(A)', advance='no') ' '
            end do
            write(u, *)
        end do

        close(u)
    end subroutine save_array_2d

    subroutine save_array_bin(filename, A)
        use, intrinsic :: iso_fortran_env, only : int32
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: A(..)

        integer :: u
        integer(int32) :: r
        integer(int32), allocatable :: shp(:)

        open(newunit=u, file=filename, status="replace", action="write", &
            form="unformatted", access="stream")

        select rank(A)

          rank(1)
            r = 1_int32
            allocate(shp(1))
            shp = int(shape(A), kind=int32)

            write(u) r
            write(u) shp
            write(u) A

          rank(2)
            r = 2_int32
            allocate(shp(2))
            shp = int(shape(A), kind=int32)

            write(u) r
            write(u) shp
            write(u) A

          rank(3)
            r = 3_int32
            allocate(shp(3))
            shp = int(shape(A), kind=int32)

            write(u) r
            write(u) shp
            write(u) A

          rank default
            close(u)
            error stop "save_array_bin: rank nao suportado"
        end select

        close(u)

        if (allocated(shp)) deallocate(shp)
    end subroutine save_array_bin
end module array_io
