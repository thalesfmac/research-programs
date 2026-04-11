program test_lapack_blas
    use precision, only : dp
    use matrix_operations, only : identity_matrix, invert, matmul2
    implicit none

    complex(dp), allocatable :: A(:,:), B(:,:), C(:,:), I_ref(:,:)
    real(dp), allocatable :: w(:)
    real(dp) :: tol = 1.0e-12_dp

    ! --- Test Inversion ---
    allocate(A(2,2), I_ref(2,2), C(2,2))
    call identity_matrix(I_ref)
    A(1,1) = (2.0_dp, 0.0_dp); A(1,2) = (0.0_dp, 1.0_dp)
    A(2,1) = (0.0_dp, -1.0_dp); A(2,2) = (2.0_dp, 0.0_dp)

    B = A ! Keep copy
    call invert(B)
    call matmul2(A, B, C)

    if (maxval(abs(C - I_ref)) < tol) then
        print *, "Invert test: PASSED"
    else
        print *, "Invert test: FAILED"
    end if

    ! Add similar blocks for matmul2 and diagonalize...
end program test_lapack_blas
