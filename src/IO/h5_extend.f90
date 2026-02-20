if (allocated(offset)) deallocate(offset)
if (allocated(dims_new)) deallocate(dims_new)
if (allocated(dims_old)) deallocate(dims_old)
if (allocated(maxdims)) deallocate(maxdims)
if(allocated(chunk_dims)) deallocate(chunk_dims)

call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr);
dummy = check_return_value(hdferr, "h5_extend", "h5fopen")

call h5dopen_f(file_id, dataset_name, dataset_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5dopen")

call h5dget_space_f(dataset_id, dataspace_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5dget_space")   

! find the old dimensions
! allocate(dims_old, mold=dims)
call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5sget_simple_extent_ndims")

allocate(dims_old(rank))
allocate(maxdims(rank))
allocate(chunk_dims(rank))

! get dataset creation plist
call h5dget_create_plist_f(dataset_id, prop_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5pcreate_f")

! Find out the chunking information
call H5Pget_chunk_f(prop_id, rank, chunk_dims, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "H5Pget_chunk_f")

! Found out the chunking info. Now we can close the dataset and reopen it
! with the correct chunking information

! Create the good access property list, with chunking info
call h5dget_access_plist_f(dataset_id, plist_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5dget_access_plist_f")

call h5dclose_f(dataset_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5dclose")

!!! End finding chunks & dataset dimensions

!!! Now open for writing the new dataset
!!! with efficient chunking caching options.

! Storage in double precision: 8 bytes * 1024 chunks in the cache * 
! the size of the chunk
n_bytes_chunk = n_chunks_in_cache*8*product(chunk_dims)
! The second parameter should be hundred times the nr of chunks that can fit
! in n_butes_chunk
call h5pset_chunk_cache_f(plist_id, int(n_chunks_in_cache*101,8), n_bytes_chunk, .0 , hdferr)
dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5pset_chunk_cache")

call h5dopen_f(file_id, dataset_name, dataset_id, hdferr, plist_id)
dummy = check_return_value(hdferr, "h5_extend", "h5dopen")

call h5sget_simple_extent_dims_f(dataspace_id, dims_old, maxdims, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5sget_simple_extent_dims")

allocate(dims_new(rank))
dims_new = dims_old
! change the one dimension that is different
dims_new(extended_dimension_id) = dims_old(extended_dimension_id) + dims_ext(extended_dimension_id)

call h5dset_extent_f(dataset_id, dims_new, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5dset_extent_f")
! print *, 'dims_new', dims_new
!!Select a hyperslab in extended portion of dataset:

!!close the dataspace:
call h5sclose_f(dataspace_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5sclose_f")

! Reopen and take the info about the dataspace. Now it will have the new size
call h5dget_space_f(dataset_id, dataspace_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5dget_space")

allocate(offset(rank))
offset= 0 ; offset(extended_dimension_id) = dims_old(extended_dimension_id)

call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, dims_ext, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5sselect_hyperslab_f")

! print*, 'offset ',offset
! print*, 'dims_ext', dims_ext

! Create a simple memspace where to store the new data
call h5screate_simple_f(rank, dims_ext, memspace_id, hdferr, maxdims)
dummy = check_return_value(hdferr, "h5_extend", "h5screate_simple_f")

! Write the new data to the file
call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dims_ext, hdferr, memspace_id, dataspace_id)
dummy = check_return_value(hdferr, "h5_extend", "h5dwrite")        

! close resources
call h5dclose_f(dataset_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5dclose")
call h5sclose_f(dataspace_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5sclose_dataspace")
call h5sclose_f(memspace_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5sclose_memspace")
call h5fclose_f(file_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5fclose")    

