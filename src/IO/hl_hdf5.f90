module HL_HDF5 !High level HDF5 interface


    USE HDF5 ! This module contains all necessary modules 
    USE ISO_C_BINDING
    USE GlobalVariables, only : chunk_dim, n_chunks_in_cache
    IMPLICIT NONE
    INTEGER     ::   hdferr ! Error flag

    interface h5_write
        module procedure h5_write_1d
        module procedure H5_write_2d
        module procedure H5_write_3d
    end interface 

    interface h5_read
        module procedure h5_read_1d
    end interface 

    interface h5_read_block
        module procedure h5_read_block_4d
        module procedure h5_read_block_3d
    end interface 

    interface h5_write_at_step
        module procedure h5_write_1d_at_step
        module procedure H5_write_2d_at_step
        module procedure H5_write_3d_at_step
    end interface 

    interface h5_extend
        module procedure h5_extend_1d
        module procedure h5_extend_2d
        module procedure h5_extend_3d
        module procedure h5_extend_4d
    end interface     

    contains

    subroutine h5_check_version()
        integer :: major = 0;
        integer :: minor = 0;
        integer :: patch = 0;
        integer :: error = 0;

        call h5check_version_f(major, minor, patch, error)

        print*, "This version of h5 is: ", major, ".", minor, ".", patch

    end subroutine h5_check_version

    subroutine h5_file_open(file_name, file_id)
       
        CHARACTER(*), intent(IN) :: file_name
        INTEGER(HID_T), intent(INOUT) :: file_id       ! File identifier
        logical :: dummy

        CALL h5open_f(hdferr)
        dummy = check_return_value(hdferr, "h5_file_create", "h5fopen_f")
        CALL h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr) 
        !H5F_ACC_TRUNC_F overwrite existing file
        dummy = check_return_value(hdferr, "h5_file_create", "h5fopen_f")
        CALL h5fclose_f(file_id, hdferr)
        dummy = check_return_value(hdferr, "h5_file_create", "h5fclose_f")

    end subroutine h5_file_open

    subroutine h5_file_create(file_name, file_id)
       
        CHARACTER(*), intent(IN) :: file_name
        INTEGER(HID_T), intent(INOUT) :: file_id       ! File identifier
        logical :: dummy

        CALL h5open_f(hdferr)
        dummy = check_return_value(hdferr, "h5_file_create", "h5fopen_f")
        CALL h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, hdferr) !H5F_ACC_TRUNC_F overwrite existing file
        dummy = check_return_value(hdferr, "h5_file_create", "h5fcreate_f")
        CALL h5fclose_f(file_id, hdferr)
        dummy = check_return_value(hdferr, "h5_file_create", "h5fclose_f")
    end subroutine h5_file_create

    function h5_dataset_exists(file_name, dataset_name)

        character(*) :: file_name, dataset_name
        integer(HID_T) :: file
        logical :: isExisting, h5_dataset_exists, dummy

        call h5fopen_f(file_name, H5F_ACC_RDWR_F, file, hdferr);
        isExisting = .FALSE.

        call H5LExists_f(file, dataset_name, isExisting, hdferr)
     
        call H5Fclose_f(file, hdferr)

        dummy = check_return_value(hdferr, "h5_dataset_exists", "h5fclose_f")

        h5_dataset_exists = isExisting

    end function h5_dataset_exists

    subroutine h5_dataset_dimension_all(file_name, dataset_name, dataset_dimension)
        ! Returns the dataset dimension in a certain direction
        character(*) :: file_name, dataset_name
        integer(HID_T) :: file_id, dataset_id, dataspace_id
        integer :: rank
        logical :: dummy
        integer(HSIZE_T),allocatable :: maxdims(:), dataset_dimension(:)
        integer(HSIZE_T) :: dimId
        
        call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr);
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5fopen")

        call h5dopen_f(file_id, dataset_name, dataset_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5dopen")

        call h5dget_space_f(dataset_id, dataspace_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5dget_space")

        call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sget_simple_extent_ndims")

        allocate(dataset_dimension(rank), maxdims(rank))
        
        call h5sget_simple_extent_dims_f(dataspace_id, dataset_dimension, maxdims, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sget_simple_extent_dims")

        call h5dclose_f(dataset_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5dclose")

        call h5sclose_f(dataspace_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sclose")

        call h5fclose_f(file_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5fclose")


    end subroutine h5_dataset_dimension_all

    subroutine h5_dataset_dimension(file_name, dataset_name, dimId, dataset_dimension)
        ! Returns the dataset dimension in a certain direction
        character(*) :: file_name, dataset_name
        integer(HID_T) :: file_id, dataset_id, dataspace_id
        integer :: rank
        logical :: dummy
        integer(HSIZE_T),allocatable :: dims(:), maxdims(:)
        integer(HSIZE_T) :: dataset_dimension, dimId

        call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr);
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5fopen")

        call h5dopen_f(file_id, dataset_name, dataset_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5dopen")

        call h5dget_space_f(dataset_id, dataspace_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5dget_space")

        call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sget_simple_extent_ndims")

        allocate(dims(rank), maxdims(rank))
        call h5sget_simple_extent_dims_f(dataspace_id, dims, maxdims, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sget_simple_extent_dims")

        call h5dclose_f(dataset_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5dclose")

        call h5sclose_f(dataspace_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sclose")

        call h5fclose_f(file_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5fclose")

        dataset_dimension = dims(dimId)

        deallocate(dims, maxdims)

    end subroutine h5_dataset_dimension

    subroutine h5_dataset_create_chunked(file_name, dataset_name, rank, dims, maxdims, chunk_dims)

        character(*) :: file_name, dataset_name
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id, prop_id, plist_id, filter_id
        logical :: dummy, avail
        integer(HSIZE_T) :: dims(:), maxdims(:), chunk_dims(:), &
            n_bytes_chunk

        integer :: rank, filter_info
        ! real(kind=8) :: data(:,:,:,:)

        call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr);
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5fopen")

        !  Check if gzip compression is available and can be used for both
        !  compression and decompression.  Normally we do not perform error
        !  checking in these examples for the sake of clarity, but in this
        !  case we will make an exception because this filter is an
        !  optional part of the hdf5 library.
        !  
        CALL h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, avail, hdferr)

        IF (.NOT.avail) THEN
            WRITE(*,'("gzip filter not available.",/)')
            STOP
         ENDIF
         CALL h5zget_filter_info_f(H5Z_FILTER_DEFLATE_F, filter_info, hdferr)

        call h5screate_simple_f(rank, dims, dataspace_id, hdferr, maxdims)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5screate_simple_f")

        call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5p_create_f")        
        ! CALL h5pset_deflate_f(prop_id, 5, hdferr)
        call h5pset_chunk_f(prop_id, size(chunk_dims), chunk_dims, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5pset_chunk")

        call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, &
                    hdferr, prop_id, H5P_DEFAULT_F, H5P_DEFAULT_F)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5dcreate_f")

        call h5dget_access_plist_f(dataset_id, plist_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5dget_access_plist_f")

        ! Storage in double precision: 8 bytes * n_bytes_cache chunks in the cache * 
        ! the size of the chunk
        n_bytes_chunk = n_chunks_in_cache*8*product(chunk_dims)
        ! The second parameter should be hundred times the nr of chunks that can fit
        ! in n_butes_chunk
        call h5pset_chunk_cache_f(plist_id, int(n_chunks_in_cache*101,8), n_bytes_chunk, 1.0 , hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5pset_chunk_cache")

        call h5dclose_f(dataset_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5dclose")

        call h5pclose_f(prop_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5pclose")

        call h5sclose_f(dataspace_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sclose")

        call h5fclose_f(file_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5fclose")

        
    end subroutine h5_dataset_create_chunked

    subroutine h5_create_hard_link(file_name, dataset_target, dataset_source)

            character(*) :: file_name, dataset_target, dataset_source
            integer(HID_T) :: file_id, dataset_target_id, &
                dataspace_id, prop_id, plist_id, filter_id
            logical :: dummy, avail
    
            integer :: rank, filter_info
            

            call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr);
            dummy = check_return_value(hdferr, "h5_create_hard_link", "h5fopen")

            ! call h5dopen_f(file_id, dataset_target, dataset_target_id, hdferr)
            ! dummy = check_return_value(hdferr, "h5_create_hard_link", "h5dopen_f")

            call h5lcreate_hard_f(file_id, dataset_target, file_id, dataset_source, hdferr)
            dummy = check_return_value(hdferr, "h5_create_hard_link", "h5lcreate_hard_f")

            call h5fclose_f(file_id, hdferr)
            dummy = check_return_value(hdferr, "h5_create_hard_link", "h5fclose")
    
    end subroutine h5_create_hard_link

    subroutine h5_extend_4d(file_name, dataset_name, extended_dimension_id, &
        dims_ext, data)

    ! Data must have the same bounds of dims_ext
    ! Routine to insert a 3D array at a location

    character(*) :: file_name, dataset_name
    real(kind=8) :: data(:,:,:,:)
    integer(HID_T) :: file_id, dataset_id, &
        dataspace_id, extended_dimension_id, memspace_id, prop_id, &
        plist_id, n_bytes_chunk
    logical :: dummy
    integer(HSIZE_T) :: dims_ext(:)
    integer(HSIZE_T),allocatable :: maxdims(:), &
            dims_old(:), dims_new(:), offset(:), chunk_dims(:)
    integer :: rank

    include "h5_extend.f90"

end subroutine h5_extend_4d

    ! 3D routines
    subroutine h5_write_3d(file_name, dataset_name, data)

        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:,:,:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id
        logical :: dummy

        include "h5_write.f90"

    end subroutine h5_write_3d
       
    subroutine h5_write_3d_at_step(file_name, dataset_name, extended_dimension_id, &
            step, dims, data)

        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:,:,:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id, extended_dimension_id, memspace_id
        integer(HSIZE_T) :: dims(:)
        integer(HSIZE_T),allocatable :: maxdims(:), &
                dims_old(:), dims_new(:), offset(:)        
        logical :: dummy
        integer :: rank, i, step

        include "h5_write_at_step.f90"

    end subroutine h5_write_3d_at_step

    subroutine h5_extend_3d(file_name, dataset_name, extended_dimension_id, &
            dims_ext, data)

        ! Data must have the same bounds of dims_ext
        ! Routine to insert a 3D array at a location

        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:,:,:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id, extended_dimension_id, memspace_id, prop_id, &
            plist_id, n_bytes_chunk
        logical :: dummy
        integer(HSIZE_T) :: dims_ext(:)
        integer(HSIZE_T),allocatable :: maxdims(:), &
                dims_old(:), dims_new(:), offset(:), chunk_dims(:)
        integer :: rank

        include "h5_extend.f90"

    end subroutine h5_extend_3d

    ! 2D routines
    subroutine h5_write_2d(file_name, dataset_name, data)

        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:,:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id
        logical :: dummy

        include "h5_write.f90"

    end subroutine h5_write_2d

    subroutine h5_write_2d_at_step(file_name, dataset_name, extended_dimension_id, &
        step, dims, data)

    character(*) :: file_name, dataset_name
    real(kind=8) :: data(:,:)
    integer(HID_T) :: file_id, dataset_id, &
        dataspace_id, extended_dimension_id, memspace_id
    integer(HSIZE_T) :: dims(:)
    integer(HSIZE_T),allocatable :: maxdims(:), &
            dims_old(:), dims_new(:), offset(:)        
    logical :: dummy
    integer :: rank, i, step

    include "h5_write_at_step.f90"

end subroutine h5_write_2d_at_step

    subroutine h5_extend_2d(file_name, dataset_name, extended_dimension_id, &
            dims_ext, data)

        ! Data must have the same bounds of dims_ext
        ! Routine to insert a 2D array at a location

        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:,:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id, extended_dimension_id, memspace_id, prop_id, &
            plist_id, n_bytes_chunk
        logical :: dummy
        integer(HSIZE_T) :: dims_ext(:)
        integer(HSIZE_T),allocatable :: maxdims(:), &
                dims_old(:), dims_new(:), offset(:), chunk_dims(:)
        integer :: rank

        include "h5_extend.f90"

    end subroutine h5_extend_2d

    ! 1D routines
    subroutine h5_write_1d(file_name, dataset_name, data)
        ! We can also write a scalar as 1D array
        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id
        logical :: dummy

        include "h5_write.f90"

    end subroutine h5_write_1d
       
    ! 1D routines
    subroutine h5_read_1d(file_name, dataset_name, data)
        ! We can also write a scalar as 1D array
        character(*) :: file_name, dataset_name
        real(kind=8),allocatable :: data(:)
        integer(HSIZE_T) :: one, dataset_dimension
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id
        logical :: dummy

        ! Check the dimension
        one = 1
        call h5_dataset_dimension(file_name, dataset_name, one, dataset_dimension)
        allocate(data(dataset_dimension))
        include "h5_read.f90"

    end subroutine h5_read_1d

    subroutine h5_read_block_4d(file_name, dataset_name, data, block_offset)
        ! We can also write a scalar as 1D array
        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:,:,:,:)
        integer(HSIZE_T) :: one
        integer(HSIZE_T) :: block_offset(:) !, block_size(:)
        integer(HSIZE_T), allocatable ::dataset_dimension(:), &
            h5count(:), h5block(:), h5stride(:), h5start(:), dimsr(:), &
            maxdims(:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id, creation_id, memspace
        logical :: dummy
        integer :: layout, i, rank

        include "h5_read_block.f90"

    end subroutine h5_read_block_4d

    subroutine h5_read_block_3d(file_name, dataset_name, data, block_offset)
        ! We can also write a scalar as 1D array
        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:,:,:)
        integer(HSIZE_T) :: one
        integer(HSIZE_T) :: block_offset(:) !, block_size(:)
        integer(HSIZE_T), allocatable ::dataset_dimension(:), &
            h5count(:), h5block(:), h5stride(:), h5start(:), dimsr(:), &
            maxdims(:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id, creation_id, memspace
        logical :: dummy
        integer :: layout, i, rank

        include "h5_read_block.f90"  

    end subroutine h5_read_block_3d

    subroutine h5_write_1d_at_step(file_name, dataset_name, extended_dimension_id, &
        step, dims, data)

    character(*) :: file_name, dataset_name
    real(kind=8) :: data(:)
    integer(HID_T) :: file_id, dataset_id, &
        dataspace_id, extended_dimension_id, memspace_id
    integer(HSIZE_T) :: dims(:)
    integer(HSIZE_T),allocatable :: maxdims(:), &
            dims_old(:), dims_new(:), offset(:)        
    logical :: dummy
    integer :: rank, i, step

    include "h5_write_at_step.f90"

    end subroutine h5_write_1d_at_step


    subroutine h5_extend_1d(file_name, dataset_name, extended_dimension_id, &
            dims_ext, data)

        ! Data must have the same bounds of dims_ext

        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id, extended_dimension_id, memspace_id, prop_id, &
            plist_id, n_bytes_chunk
        logical :: dummy
        integer(HSIZE_T) :: dims_ext(:)
        integer(HSIZE_T),allocatable :: maxdims(:), &
                dims_old(:), dims_new(:), offset(:), chunk_dims(:)
        integer :: rank

        include "h5_extend.f90"

    end subroutine h5_extend_1d    

    function check_return_value(status, calling_function, returning_function)
        character(*) :: calling_function, returning_function
        integer :: status
        logical :: check_return_value

        if (status < 0) then
            print *, "Error using ", calling_function
            print *, "Failed during call to ", returning_function
            print *, "Status is ", status
            check_return_value = .FALSE.
        end if

        check_return_value = .TRUE.

    end function check_return_value

end module HL_HDF5
