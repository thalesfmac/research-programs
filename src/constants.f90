module constants
    use precision, only : dp
    implicit none
    public

    complex(dp), parameter :: CI = (0.0_dp, 1.0_dp)
    real(dp),    parameter :: PI = acos(-1.0_dp)
    real(dp),    parameter :: INV_PHI = (sqrt(5.0_dp) - 1.0_dp) / 2.0_dp

end module constants
