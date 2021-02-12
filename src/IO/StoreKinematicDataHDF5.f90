!
! ****************************************************************
!
SUBROUTINE StoreKinematicDataHDF5(Nx,Ny,Nz,io,it)
    !
    ! ****************************************************************
    !
    !>
    !! Write 3D kinematics data to an unformatted binary file
    !!
    !!
    !! By Allan P. Engsig-Karup.
    !<
    USE GlobalVariables
    USE hl_hdf5
    use kinematicsArray
    IMPLICIT NONE
    ! Input parameters
    INTEGER :: Nx, Ny, Nz, io, it
    ! Local variables
    INTEGER ::  i, j, k, i0, i1, is, j0, j1, js
    INTEGER :: FOUT, i3, iWriteAt
    LOGICAL :: extend, hdf5_file_exists, write_kinematics
    REAL(KIND=long), DIMENSION(:,:), POINTER :: x, y, h, hx, hy, eta, etax, etay
    REAL(KIND=long), DIMENSION(:), POINTER   :: z, tmpxval
    ! Automatic work space
    REAL(KIND=long) :: U(Nz,Nx,Ny), V(Nz,Nx,Ny), W(Nz,Nx,Ny), d(Nx,Ny)
    REAL(KIND=long) :: tmpx(Nx), tmpy(Ny)
    CHARACTER(len=30) :: form
    REAL(KIND=long) :: hint, etaint, dint
    REAL(KIND=long) :: Uint(Nz), Vint(Nz)
    REAL(KIND=long) :: Ux(Nz,Nx,Ny),Uy(Nz,Nx,Ny),Uz(Nz,Nx,Ny)
    REAL(KIND=long) ::              Vy(Nz,Nx,Ny),Vz(Nz,Nx,Ny)
    REAL(KIND=long) ::                           Wz(Nz,Nx,Ny)
    CHARACTER(LEN=30) :: h5file
    INTEGER(HID_T) :: extended_dimension_id, nKinSteps, extraSteps = 2, &
                      nOutputKinSteps  
    INTEGER(HSIZE_T), ALLOCATABLE :: dims_ext(:)
    INTEGER(HSIZE_T), SAVE :: maxdims1(1), maxdims2(3), maxdims3(4), &
                        extdims1(1), extdims2(3), extdims3(4), lenTime, &
                        chunkdims1(1), chunkdims2(3), chunkdims3(4)
                        
    INTEGER(HSIZE_T):: nx_save, ny_save, nz_save, onei = 1, zeroi = 0
    REAL(KIND=long) :: x3d(Nz, Ny, Nx), y3d(Nz, Ny, Nx), z3d(Nz, Ny, Nx)
    REAL(KIND=long), ALLOCATABLE, SAVE :: timeStoredInHDF5(:), dummy2d(:,:), &
       dummy3d(:,:,:), dummy4d(:,:,:,:), dummy(:)
    ! Assign the local pointers
    !
    x => FineGrid%x; y => FineGrid%y; z => FineGrid%z; h => FineGrid%h; hx => FineGrid%hx
    hy => FineGrid%hy; eta => WaveField%E; etax => WaveField%Ex; etay => WaveField%Ey
    
    ! Save n_buffer_timesteps times as many time steps as in t stride
    ! We need also to store some extra steps, as we compute the kin acc with a central scheme
    ! So at every time step, we compute the kin acc at two timesteps back
    nKinSteps = Output(io)%tstride * n_buffer_timesteps + extraSteps
    nOutputKinSteps = Output(io)%tstride * n_buffer_timesteps
    
    ! Determine if we have to write at this timestep
    IF (it >= Output(io)%tbeg .and. it <= Output(io)%tend) THEN !
        write_kinematics = .TRUE.
    else
        write_kinematics = .FALSE.
    end if 
    
    !
    ! Shift the horizontal grid point position to account for the ghost points.  
    !
    i0=Output(io)%xbeg+GhostGridX; i1=Output(io)%xend+GhostGridX; is=Output(io)%xstride; 
    j0=Output(io)%ybeg+GhostGridY; j1=Output(io)%yend+GhostGridY; js=Output(io)%ystride; 


    ! Some parameters that are necessary for h5 saving
    nx_save = size((/(i, i=i0,i1,is)/))
    ny_save = size((/(j, j=j0,j1,js)/))
    nz_save = Nz-GhostGridZ

    ! Set max dimensions for h5 writing
    ! FP 20190705: According to my tests, this syntax works on HDF5/1.8.17 and HDF5/1.8.5
    maxdims1 = (/-1/)
    maxdims2 = (/-1, -1, -1/)
    maxdims3 = (/-1, -1, -1, -1/)
    ! FP 20190705: According to my tests, this syntax works only on HDF5/1.8.5
    ! maxdims1 = (/H5S_UNLIMITED_F/)
    ! maxdims2 = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
    ! maxdims3 = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
    
    ! Chunk dims
    ! FP 20190705: The chunk dims is quite big in the time-axis.
    ! This is done to facilitate the access of the data, which are usually read
    ! as whole time series: usually, all the time steps of one dataset are read at one 
    ! particular x,y,z location
    
    chunkdims1 = (/chunk_dim*onei/)
    chunkdims2 = (/ny_save, nx_save, chunk_dim*onei/)
    chunkdims3 = (/nz_save, ny_save, nx_save, chunk_dim*onei/)
    
    WRITE(h5file, "(A18,I3.3,A3)") "WaveKinematicsZone",fntH5(io),".h5"
    if ((write_kinematics).and.(it < 100)) then
        WRITE(*,FMT='(A,A)') '  File output of h5 file number = ',h5file
    end if
    
    IF(it==0)THEN

        ! FabioPierella 20190319
        ! Output in HDF5 files. We need some preparation on the dataset before
        ! we can output it to file.
        ! 1. Fortran is column-major. We want to therefore transpose the arrays before output to have a row-major HDF5 file structure.
        ! 2. This HDF5 file output wants  to be consistent with OW3D-GPU version, therefore some variables need to be made 3D before being output (e.g. the position arrays)

        ! The kinematicsArray will hold 5 kinematics timesteps.
        ! They will be used to compute the kinematic acceleration at the location
        ! and the dynamics pressure, so that we can save it in the .h5 file.
        call allocateZoneKin(Zones(io), nKinSteps, nOutputKinSteps)
        call allocatePointers(Zones(io), nx_save, ny_save, nz_save)

        ! If the time is zero, then we consider we are starting the simulation from scratch.
        ! If we restart, we consider we do have these fields already.
        ! Initialize all datasets.
        WRITE(h5file, "(A18,I3.3,A3)") "WaveKinematicsZone",fntH5(io),".h5"

    ! Dimensions of the extended dataset, when appending to an existing dataset.
        extdims1 = (/onei/)
        extdims2 = (/ny_save, nx_save, onei/)
        extdims3 = (/nz_save, ny_save, nx_save, onei/)

        ! Create
        call h5_dataset_create_chunked(h5file, 'time', 1, &
                    & extdims1, maxdims1, chunkdims1) 
        ! Surface elevation variables
        call h5_dataset_create_chunked(h5file, 'surface_elevation', 3, &
                    & extdims2, maxdims2, chunkdims2)
        call h5_dataset_create_chunked(h5file, 'surface_elevation_derivative_etax', 3, &
                    & extdims2, maxdims2, chunkdims2)
        call h5_dataset_create_chunked(h5file, 'surface_elevation_derivative_etay', 3, &
                    & extdims2, maxdims2, chunkdims2)
        ! Position variables 
        call h5_dataset_create_chunked(h5file, 'position_x', 4, &
                    & extdims3, maxdims3, chunkdims3)
        call h5_dataset_create_chunked(h5file, 'position_y', 4, &
                    & extdims3, maxdims3, chunkdims3)
        call h5_dataset_create_chunked(h5file, 'position_z', 4, &
                    & extdims3, maxdims3, chunkdims3)

        ! Velocity variables 
        call h5_dataset_create_chunked(h5file, 'velocity_u', 4, &
                    & extdims3, maxdims3, chunkdims3)
        call h5_dataset_create_chunked(h5file, 'velocity_v', 4, &
                    & extdims3, maxdims3, chunkdims3)
        call h5_dataset_create_chunked(h5file, 'velocity_w', 4, &
                    & extdims3, maxdims3, chunkdims3)
        
                    ! Velocity U
        call h5_dataset_create_chunked(h5file, 'velocity_derivative_ux', 4, &
                    & extdims3, maxdims3, chunkdims3)
        call h5_dataset_create_chunked(h5file, 'velocity_derivative_uy', 4, &
                    & extdims3, maxdims3, chunkdims3)
        call h5_dataset_create_chunked(h5file, 'velocity_derivative_uz', 4, &
                    & extdims3, maxdims3, chunkdims3)
                    
        ! Velocity V
        call h5_create_hard_link(h5file, 'velocity_derivative_uy', 'velocity_derivative_vx')
        call h5_dataset_create_chunked(h5file, 'velocity_derivative_vy', 4, &
                    & extdims3, maxdims3, chunkdims3)      
        call h5_dataset_create_chunked(h5file, 'velocity_derivative_vz', 4, &
                    & extdims3, maxdims3, chunkdims3)                          
        ! Velocity W
        call h5_create_hard_link(h5file, 'velocity_derivative_uz', 'velocity_derivative_wx')
        call h5_create_hard_link(h5file, 'velocity_derivative_vz', 'velocity_derivative_wy')
        call h5_dataset_create_chunked(h5file, 'velocity_derivative_wz', 4, &
                    & extdims3, maxdims3, chunkdims3)      

        ! Velocity kinematic accelerations
        call h5_dataset_create_chunked(h5file, 'velocity_derivative_ut', 4, &
                    & extdims3, maxdims3, chunkdims3)
        call h5_dataset_create_chunked(h5file, 'velocity_derivative_vt', 4, &
                    & extdims3, maxdims3, chunkdims3)      
        call h5_dataset_create_chunked(h5file, 'velocity_derivative_wt', 4, &
                    & extdims3, maxdims3, chunkdims3)      

        ! dynamic pressure
        call h5_dataset_create_chunked(h5file, 'dynamic_pressure_pdyn', 4, &
                    & extdims3, maxdims3, chunkdims3)      
        ! potential
        call h5_dataset_create_chunked(h5file, 'potential_phi', 4, &
                    & extdims3, maxdims3, chunkdims3)      
        ! Write the first timestep
        ! call h5_write(h5file, 'time', (/it*dt/))
        ! allocate the position arrays

        ! Fluid thickness
        DO j=1,Ny
            DO i=1,Nx
                d(i,j)=h(i,j)+eta(i,j);
            END DO
        END DO

        do k=1,Nz
            do j=1,Ny
                do i=1,Nx
                    x3d(k,j,i) = x(i,j)
                    y3d(k,j,i) = y(i,j)
                    z3d(k,j,i)  = z(k)*d(i,j)-h(i,j)
                end do 
            end do 
        end do

        extended_dimension_id = 4
        call h5_write(h5file, 'position_x', x3d(1+GhostGridZ:, j0:j1:js, i0:i1:is))                  
        call h5_write(h5file, 'position_y', y3d(1+GhostGridZ:, j0:j1:js, i0:i1:is))                          
        call h5_write(h5file, 'position_z', z3d(1+GhostGridZ:, j0:j1:js, i0:i1:is))

        ! call increment_timestep_counter(Zones(io))
        ! n_overwrites(io) = 0 ! we do not expect to overwrite any sample here since we only append

        IF(curvilinearOnOff/=0)THEN
        Print *, 'StoreKinematicData:  Saving horizontal fluid velocities is not yet implemented for curvilinear grids.'
        end if ! if (time0 == 0.0)
    
    ELSE !IF(it==0)THEN
    
       !
       IF(curvilinearOnOff == 0)THEN

            !
            ! Dump this solution slice to the output file
            !
            ! First the free surface elevation and gradient at all points in this slice
            ! 

            !
            ! The fluid thickness d=h+eta
            !
            DO j=1,Ny
            DO i=1,Nx
                d(i,j)=h(i,j)+eta(i,j);
            END DO
            END DO

            do k=1,Nz
            do j=1,Ny
                do i=1,Nx
                    x3d(k,j,i) = x(i,j)
                    y3d(k,j,i) = y(i,j)
                    z3d(k,j,i)  = z(k)*d(i,j)-h(i,j)
                end do 
            end do 
            end do
            
            !
            !
            ! Then the velocities at all points in this horizontal slice and at all sigma 
            ! (vertical) locations.  
            !
            !
            ! Compute dphi/dsigma
            !
            CALL DiffZArbitrary(phi,W,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
            FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)

            IF (FineGrid%Nx>1) THEN
            !
            ! Compute dphi/dx 
            !	
                CALL DiffXEven(phi,U,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
                    FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,alpha)
                IF ( LinearOnOff /= 0) THEN
                    !
                    ! Add in the chain rule contribution to get the velocity
                    !
                    Do j=1,Ny
                        Do i=1,Nx
                            Do k=1,Nz
                            U(k,i,j) = U(k,i,j) + ((1-z(k))/d(i,j)*hx(i,j)-z(k)/d(i,j)*etax(i,j))*W(k,i,j)
                            END Do
                        END Do
                    END Do
                END IF
            ELSE
                U=zero
            END IF

            ! 
            IF (FineGrid%Ny>1) THEN
            ! dphi/dy
                CALL DiffYEven(phi,V,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
                FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,beta)
                IF ( LinearOnOff /= 0) THEN
                    Do j=1,Ny
                        Do i=1,Nx
                            Do k=1,Nz
                            V(k,i,j) = V(k,i,j)+((1-z(k))/d(i,j)*hy(i,j)-z(k)/d(i,j)*etay(i,j))*W(k,i,j)
                            END Do
                        END Do
                    END Do
                END IF
            ELSE
                V=zero
            END IF

            !
            !--------------------------------------- d /dsigma
            !
            CALL DiffZArbitrary(U,Uz,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
            FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
            CALL DiffZArbitrary(V,Vz,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
            FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
            CALL DiffZArbitrary(W,Wz,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
            FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
            ! ---- See below: we still need to divide all three by fluid thickness
            !


            !--------------------------------------- d /dy
            IF (FineGrid%Ny>1) THEN
                CALL DiffYEven(U,Uy,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
                    FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,beta)
                CALL DiffYEven(V,Vy,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
                    FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,beta)

                ! ---- See below: we still need to divide dw/dy by fluid thickness   
            
            ELSE
                Uy=zero
                Vy=zero
            END IF

            !--------------------------------------- d /dx
            IF (FineGrid%Nx>1) THEN
                CALL DiffXEven(U,Ux,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
                    FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,alpha)
                ! ---- See below: we still need to divide dw/dx by fluid thickness                        
            ELSE
                Ux=zero
                Uy=zero
            END IF

            ! sigs ends lines
    
            ! Save the velocity information in the array of Kinematics, for later
            ! computation of the kinematic acceleration
            ! TODO: save all variables 
            call cycle(Zones(io)%Kinematics)

            Zones(io)%Kinematics(nKinSteps)%time = time0+it*dt
            
            Zones(io)%Kinematics(nKinSteps)%X = x3d(1+GhostGridZ:, j0:j1:js, i0:i1:is)
            Zones(io)%Kinematics(nKinSteps)%Y = y3d(1+GhostGridZ:, j0:j1:js, i0:i1:is)
            Zones(io)%Kinematics(nKinSteps)%Z = z3d(1+GhostGridZ:, j0:j1:js, i0:i1:is)

            Zones(io)%Kinematics(nKinSteps)%U = U(1+GhostGridZ:, i0:i1:is, j0:j1:js);
            Zones(io)%Kinematics(nKinSteps)%V = V(1+GhostGridZ:, i0:i1:is, j0:j1:js);
            Zones(io)%Kinematics(nKinSteps)%Ux = Ux(1+GhostGridZ:, i0:i1:is, j0:j1:js);
            Zones(io)%Kinematics(nKinSteps)%Uy = Uy(1+GhostGridZ:, i0:i1:is, j0:j1:js);
            Zones(io)%Kinematics(nKinSteps)%Vy = Vy(1+GhostGridZ:, i0:i1:is, j0:j1:js);
            DO K=1+GhostGridZ,NZ
                Zones(io)%Kinematics(nKinSteps)%W(K-GhostGridZ,:,:) = W(K, i0:i1:is, j0:j1:js)/d(i0:i1:is, j0:j1:js);
                Zones(io)%Kinematics(nKinSteps)%Uz(K-GhostGridZ,:,:) = Uz(K, i0:i1:is, j0:j1:js)/d(i0:i1:is, j0:j1:js);
                Zones(io)%Kinematics(nKinSteps)%Vz(K-GhostGridZ,:,:) = Vz(K, i0:i1:is, j0:j1:js)/d(i0:i1:is, j0:j1:js);               
                !FP20190430: We need to divide Wz by twice the fluid thickness. Before today it was mistakenly divided only once.
                Zones(io)%Kinematics(nKinSteps)%Wz(K-GhostGridZ,:,:) = Wz(K, i0:i1:is, j0:j1:js)/d(i0:i1:is, j0:j1:js)**2; 
                ! FP20200518: Adding the potential to compute the dyn pressure
                Zones(io)%Kinematics(nKinSteps)%phi = phi(1+GhostGridZ:, i0:i1:is, j0:j1:js);               
            END DO            
            Zones(io)%Kinematics(nKinSteps)%Eta = Eta(i0:i1:is, j0:j1:js);
            Zones(io)%Kinematics(nKinSteps)%Etax = Etax(i0:i1:is, j0:j1:js);
            IF (FineGrid%Ny>1) THEN
                Zones(io)%Kinematics(nKinSteps)%Etay = Etay(i0:i1:is, j0:j1:js);
            else
                Zones(io)%Kinematics(nKinSteps)%Etay = 0.
            end if 

            call increment_timestep_counter(Zones(io)) ! signal that we have added another timestep

            if (Zones(io)%number_of_saved_timesteps >= 5) then 
                ! if we have saved more than 5 timesteps  then we calculate the centered 5-points derivative
                ! and also the dynamic pressure
                call calculateKinAcceleration(Zones(io), dt, z(1+GhostGridZ:), rho)
            end if 

            if (write_kinematics) then
                if (mod(Zones(io)%number_of_saved_timesteps, nOutputKinSteps) == 0) then 
                    ! the array is full, save it
                    ! save it but remember to use the tstride to skip the unnecessary files

                    extdims1 = (/onei*n_buffer_timesteps/)
                    extdims2 = (/ny_save, nx_save, onei*n_buffer_timesteps/)
                    extdims3 = (/nz_save, ny_save, nx_save, onei*n_buffer_timesteps/)

                    if (.not.(allocated(dummy4d))) allocate(dummy4d(nz_save, nx_save, ny_save, onei*n_buffer_timesteps))
                    if (.not.(allocated(dummy3d))) allocate(dummy3d(nx_save, ny_save, onei*n_buffer_timesteps))

                    WRITE(h5file, "(A18,I3.3,A3)") "WaveKinematicsZone",fntH5(io),".h5"
    
                    extended_dimension_id = 1
                    ! Save previous 100 time steps
                    call h5_extend(h5file, 'time', extended_dimension_id, extdims1, (/(Zones(io)%Kinematics(j)%time, j=1, nOutputKinSteps, Output(io)%tstride)/))

                    extended_dimension_id = 4
                    
                    call arrayFromVariable4d(zones(io),'X', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'position_x', extended_dimension_id, extdims3, &
                        & dummy4d)         

                    
                    call arrayFromVariable4d(zones(io),'Y', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'position_y', extended_dimension_id, extdims3, &
                        & dummy4d) 

                    call arrayFromVariable4d(zones(io),'Z', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'position_z', extended_dimension_id, extdims3, &
                        & dummy4d)

                    extended_dimension_id = 3
                    call arrayFromVariable3d(zones(io),'Eta', output(io)%tstride, dummy3d)
                    call h5_extend(h5file, 'surface_elevation', extended_dimension_id, extdims2, &
                        & dummy3d)
                    
                    ! Write only if there is more than one point in the x direction
                    if (FineGrid%Nx>1) then
                        call arrayFromVariable3d(zones(io),'Etax', output(io)%tstride, dummy3d)
                        call h5_extend(h5file, 'surface_elevation_derivative_etax', extended_dimension_id, extdims2, &
                        & dummy3d)

                    else ! just write 0s
                        dummy3d = 0.
                        call h5_extend(h5file, 'surface_elevation_derivative_etax', extended_dimension_id, extdims2, &
                        & dummy3d)     
                    end if

                    ! Write only if there is more than one point in the y direction
                    if (FineGrid%Ny>1) then
                        call arrayFromVariable3d(zones(io),'Etay', output(io)%tstride, dummy3d)                        
                        call h5_extend(h5file, 'surface_elevation_derivative_etay', extended_dimension_id, extdims2, &
                        & dummy3d)
                    else ! just write 0s
                        dummy3d = 0.                                            
                        call h5_extend(h5file, 'surface_elevation_derivative_etay', extended_dimension_id, extdims2, &
                        & dummy3d)
                    end if
                
                    extended_dimension_id = 4
                    ! velocities
                    call arrayFromVariable4d(zones(io),'U', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'velocity_u', extended_dimension_id, extdims3, &
                    & dummy4d)

                    call arrayFromVariable4d(zones(io),'V', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'velocity_v', extended_dimension_id, extdims3, &
                    & dummy4d)
                    
                    call arrayFromVariable4d(zones(io),'W', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'velocity_w', extended_dimension_id, extdims3, &
                    & dummy4d)
                    
                    ! Potential
                    call arrayFromVariable4d(zones(io),'phi', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'potential_phi', extended_dimension_id, extdims3, &
                    & dummy4d)    

                    ! Kinematics accelerations and dynamics pressure
                    ! They are a little different from the others. Since we compute them with a centered scheme with 5 points
                    ! We need always to first extend the file, writing a "zero" value, so that it has the same nr. of points
                    ! as the other quantities.
                    
                    call arrayFromVariable4d(zones(io),'Ut', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'velocity_derivative_ut', extended_dimension_id, extdims3, &
                    & dummy4d)    

                    call arrayFromVariable4d(zones(io),'Vt', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'velocity_derivative_vt', extended_dimension_id, extdims3, &
                    & dummy4d)    

                    call arrayFromVariable4d(zones(io),'Wt', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'velocity_derivative_wt', extended_dimension_id, extdims3, &
                    & dummy4d)

                    call arrayFromVariable4d(zones(io),'pdyn', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'dynamic_pressure_pdyn', extended_dimension_id, extdims3, &
                    & dummy4d)
                    
                    ! velocity z gradients
                    call arrayFromVariable4d(zones(io),'Uz', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'velocity_derivative_uz', extended_dimension_id, extdims3, &
                    & dummy4d)    
                    call arrayFromVariable4d(zones(io),'Vz', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'velocity_derivative_vz', extended_dimension_id, extdims3, &
                    & dummy4d)    

                    call arrayFromVariable4d(zones(io),'Wz', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'velocity_derivative_wz', extended_dimension_id, extdims3, &
                    & dummy4d)  


                    ! velocity x gradients
                    call arrayFromVariable4d(zones(io),'Ux', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'velocity_derivative_ux', extended_dimension_id, extdims3, &
                    & dummy4d)    
                    
                    ! velocity y gradients
                    call arrayFromVariable4d(zones(io),'Uy', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'velocity_derivative_uy', extended_dimension_id, extdims3, &
                    & dummy4d) 
                    
                    call arrayFromVariable4d(zones(io),'Vy', output(io)%tstride, dummy4d)
                    call h5_extend(h5file, 'velocity_derivative_vy', extended_dimension_id, extdims3, &
                    & dummy4d)                

                end if
            end if 
       END IF
    END IF

    
    contains
    
    subroutine chainRuleContribution(InOutArray, h_gradient, eta_gradient)
    
       real(kind=long),intent(INOUT) :: InOutArray(:,:,:)
       real(kind=long),intent(IN) :: eta_gradient(:,:), h_gradient(:,:)
    
       Do j=1,Ny
          Do i=1,Nx
             Do k=1,Nz
                InOutArray(k,i,j) = InOutArray(k,i,j) + ((1-z(k))/d(i,j)*h_gradient(i,j)-z(k)/d(i,j)*eta_gradient(i,j))*W(k,i,j)
             END Do
          END Do
       END Do
    
    end subroutine chainRuleContribution
    
    
    END SUBROUTINE StoreKinematicDataHDF5
    
    