module hdf5_io
    use precision, only : dp
    use hdf5
    implicit none
    private
    public :: save_results_h5

    contains

    subroutine save_results_h5(filename, energies, transmissions)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: energies(:)
        real(dp), intent(in) :: transmissions(:,:)

        integer(hid_t) :: file_id
        integer(hid_t) :: space1_id, dset1_id
        integer(hid_t) :: space2_id, dset2_id
        integer(hsize_t), dimension(1) :: dims1
        integer(hsize_t), dimension(2) :: dims2
        integer :: error

        dims1 = [ int(size(energies,1), kind=hsize_t) ]
        dims2 = [ int(size(transmissions,1), kind=hsize_t), &
            int(size(transmissions,2), kind=hsize_t) ]

        call h5open_f(error)
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)

        call h5screate_simple_f(1, dims1, space1_id, error)
        call h5dcreate_f(file_id, "energies", H5T_NATIVE_DOUBLE, space1_id, dset1_id, error)
        call h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, energies, dims1, error)

        call h5screate_simple_f(2, dims2, space2_id, error)
        call h5dcreate_f(file_id, "transmissions", H5T_NATIVE_DOUBLE, space2_id, dset2_id, error)
        call h5dwrite_f(dset2_id, H5T_NATIVE_DOUBLE, transmissions, dims2, error)

        call h5dclose_f(dset1_id, error)
        call h5sclose_f(space1_id, error)

        call h5dclose_f(dset2_id, error)
        call h5sclose_f(space2_id, error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)
    end subroutine save_results_h5

end module hdf5_io
