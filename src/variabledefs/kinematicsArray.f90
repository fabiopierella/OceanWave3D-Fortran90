module kinematicsArray

! (c) Fabio Pierella 20190402
! In this module we create an array of pointers that will contain
! the velocities for the Wave Kinematics zone data storage (i == 20 or i == 30),
! see ReadInputFileParameters.f90
! We will save 5 time instants and then take their derivative in order to calculate 
! the kinematic acceleration and store it into the WaveKinematicsZone**.h5 file

! TODO: move into datatypes

use hdf5
use precision
implicit none

type kinArray
    real(kind=long), allocatable :: U(:,:,:)
    real(kind=long), allocatable :: V(:,:,:)
    real(kind=long), allocatable :: W(:,:,:)

    real(kind=long), allocatable :: Ux(:,:,:)
    real(kind=long), allocatable :: Uy(:,:,:)
    real(kind=long), allocatable :: Uz(:,:,:)

    real(kind=long), allocatable :: Vy(:,:,:)
    real(kind=long), allocatable :: Vz(:,:,:)
    
    real(kind=long), allocatable :: Wz(:,:,:)

    real(kind=long), allocatable :: Ut(:,:,:)
    real(kind=long), allocatable :: Vt(:,:,:)
    real(kind=long), allocatable :: Wt(:,:,:)

    real(kind=long), allocatable :: pdyn(:,:,:)
    real(kind=long), allocatable :: phi(:,:,:)
    real(kind=long), allocatable :: Eta(:,:)
    real(kind=long), allocatable :: Etax(:,:)
    real(kind=long), allocatable :: Etay(:,:)
    real(kind=long), allocatable :: X(:,:,:)
    real(kind=long), allocatable :: Y(:,:,:)
    real(kind=long), allocatable :: Z(:,:,:)
    real(kind=long), allocatable :: time

end type kinArray

type zoneKin
    type(kinArray), dimension(:), allocatable :: Kinematics !5: newest, 1:oldest
    integer(kind=long) :: number_of_saved_timesteps = 0
    integer :: id
    integer(kind=long) :: nSteps
end type zoneKin

contains

subroutine allocateZoneKin(inZone, nSteps)

implicit none
type(zoneKin) :: inZone
INTEGER(kind=long) :: nSteps

allocate(inZone%Kinematics(nSteps))
inZone%nSteps = nSteps

end subroutine allocateZoneKin

subroutine increment_timestep_counter(inZone)
    ! This subroutine increases a counter
    ! that says how many samples have been already saved
    ! The counter cannot be over 5 (max 5 timesteps are saved)
 
    implicit none
    type(zoneKin) :: inZone

    inZone%number_of_saved_timesteps = inZone%number_of_saved_timesteps +1

end subroutine increment_timestep_counter

subroutine allocatePointers(inZone, nx_save, ny_save, nz_save)

    implicit none
    type(zoneKin)  :: inZone
    integer(kind=long):: nx_save, ny_save, nz_save ! total of points to be saved
    integer        :: i

    DO I=1,inZone%nSteps
        ALLOCATE(inZone%Kinematics(i)%U(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%V(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%W(nz_save, nx_save, ny_save))

        ALLOCATE(inZone%Kinematics(i)%Ut(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Vt(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Wt(nz_save, nx_save, ny_save))
        
        ALLOCATE(inZone%Kinematics(i)%Ux(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Uy(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Uz(nz_save, nx_save, ny_save))
        
        ALLOCATE(inZone%Kinematics(i)%Vy(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Vz(nz_save, nx_save, ny_save))

        ALLOCATE(inZone%Kinematics(i)%Wz(nz_save, nx_save, ny_save))

        ALLOCATE(inZone%Kinematics(i)%X(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Y(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Z(nz_save, nx_save, ny_save))
        
        
        ALLOCATE(inZone%Kinematics(i)%phi(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%pdyn(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Eta(nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Etay(nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Etax(nx_save, ny_save))

        inZone%Kinematics(i)%U = 0.;
        inZone%Kinematics(i)%V = 0.;
        inZone%Kinematics(i)%W = 0.;

        inZone%Kinematics(i)%Ut = 0.;
        inZone%Kinematics(i)%Vt = 0.;
        inZone%Kinematics(i)%Wt = 0.;
        
        inZone%Kinematics(i)%Ux = 0.;
        inZone%Kinematics(i)%Uy = 0.;
        inZone%Kinematics(i)%Uz = 0.;

        inZone%Kinematics(i)%Vy = 0.;
        inZone%Kinematics(i)%Vz = 0.;
        
        inZone%Kinematics(i)%Wz = 0.;
        
        inZone%Kinematics(i)%X = 0.;
        inZone%Kinematics(i)%Y = 0.;
        inZone%Kinematics(i)%Z = 0.;


        inZone%Kinematics(i)%phi = 0.; 
        inZone%Kinematics(i)%pdyn = 0.;
        inZone%Kinematics(i)%Eta = 0.;
        inZone%Kinematics(i)%Etax = 0.
        inZone%Kinematics(i)%Etay = 0.

    END DO

end subroutine allocatePointers

subroutine cycle(inKinArray)
    ! This subroutine cycles the pointers so that 
    ! we can replace the oldest timestep with the newest
    implicit none
    type(kinArray) :: inKinArray(:)

    inKinArray = cshift(inKinArray, 1)
end subroutine cycle

subroutine calculateKinAcceleration(inZone, dt, sigma, rho)
    ! This subroutine calculate the derivative
    ! of eta in time for a certain zone

    implicit none
    type(zoneKin)  :: inZone
    real(kind=long)   :: dt, dummy(5), dummy2
    real(kind=long)   :: sigma(:), rho
    real(kind=long)   :: time(5), Eta_t, &
                      Ut_nocorr, Vt_nocorr, Wt_nocorr, &      ! Accelerations without the correction
                      phit_nocorr, phi_t                            ! Uncorrected dyn pressure
    integer        :: i, j, k, ii, it, &
                      time_steps, &     ! size of the stencil we use for time differenciation
                      alpha, &          ! half-width of the stencil
                      nz_save, nx_save, ny_save , &  ! size of the zone
                      nSteps ! number of total kinArray saved
    real(kind=long), allocatable, save   :: FDstencil(:,:)

    nSteps = inZone%nSteps

    time_steps = inZone%number_of_saved_timesteps
    alpha = 2
    
    ! Construct the time array
    Do i=1,5
        time(i)=dt*(i-1)
    END DO

    if (.NOT.(allocated(FDStencil))) then
        allocate(FDStencil(time_steps, time_steps))
        ! this works only if the computation is done at all time steps
        ! as it is the "ideal" situation (and then save every tstride iterations)
        CALL BuildStencilsGridX(alpha,1,time,FDStencil,5,1)
    end if

    nz_save = size(inZone%Kinematics(1)%U, 1) ! How many points in the z-direction we are saving
    nx_save = size(inZone%Kinematics(1)%U, 2) ! How many points in the x-direction we are saving
    ny_save = size(inZone%Kinematics(1)%U, 3) ! How many points in the y-direction we are saving
    

    DO j=1,ny_save
        DO i=1,nx_save
            Eta_t = &
            Dot_Product(FDStencil(3,:), (/(inZone%Kinematics(it)%Eta(i,j),it=nSteps-4,nSteps)/))
            ! Central, 5-points derivative of the surface elevation
            DO k=1,nz_save
                Ut_nocorr = &
                Dot_Product(FDStencil(3,:), (/(inZone%Kinematics(it)%U(k,i,j),it=nSteps-4,nSteps)/))
                Vt_nocorr = &
                Dot_Product(FDStencil(3,:), (/(inZone%Kinematics(it)%V(k,i,j),it=nSteps-4,nSteps)/))                        
                Wt_nocorr = &
                Dot_Product(FDStencil(3,:), (/(inZone%Kinematics(it)%W(k,i,j),it=nSteps-4,nSteps)/))
                phit_nocorr = &
                Dot_Product(FDStencil(3,:), (/(inZone%Kinematics(it)%phi(k,i,j),it=nSteps-4,nSteps)/))

                inZone%Kinematics(nSteps-2)%Ut(k,i,j) = Ut_nocorr - &
                    sigma(k)*inZone%Kinematics(nSteps-2)%Uz(k,i,j)*Eta_t
                inZone%Kinematics(nSteps-2)%Vt(k,i,j) = Vt_nocorr - &
                    sigma(k)*inZone%Kinematics(nSteps-2)%Vz(k,i,j)*Eta_t
                inZone%Kinematics(nSteps-2)%Wt(k,i,j) = Wt_nocorr - &
                    sigma(k)*inZone%Kinematics(nSteps-2)%Wz(k,i,j)*Eta_t
                ! dynamic pressure
                phi_t = phit_nocorr - &
                    sigma(k)*inZone%Kinematics(nSteps-2)%W(k,i,j)*Eta_t
                inZone%Kinematics(nSteps-2)%pdyn(k,i,j) = -rho*(phi_t + &
                    inZone%Kinematics(nSteps-2)%U(k,i,j)**2 + & 
                    inZone%Kinematics(nSteps-2)%V(k,i,j)**2 + &
                    inZone%Kinematics(nSteps-2)%W(k,i,j)**2)

            END DO
        END DO
    END Do

end subroutine calculateKinAcceleration

subroutine arrayFromVariable3d(inZone, variable, stride, res)

    implicit none
    type(zoneKin)  :: inZone
    character(*) :: variable
    real(kind=long), dimension(:,:,:) :: res
    integer :: i, stride
    integer, allocatable :: k(:)

    allocate(k(size(res, 3)))
    k = (/(i, i=1,inZone%nSteps, stride)/)

    if (variable == 'Eta') then
        do i=1,size(k)
            res(:, :, i) = inZone%Kinematics(k(i))%Eta
        end do
    else if (variable == 'Etax') then
        do i=1,size(k)
            res(:, :, i) = inZone%Kinematics(k(i))%Etax
        end do
    else if (variable == 'Etay') then
        do i=1,size(k)
            res(:, :, i) = inZone%Kinematics(k(i))%Etay
        end do
    end if 

    deallocate(k)

end subroutine arrayFromVariable3d


subroutine arrayFromVariable4d(inZone, variable, stride, res)

    implicit none
    type(zoneKin)  :: inZone
    character(*) :: variable
    real(kind=long), dimension(:,:,:,:) :: res
    integer :: i, stride
    integer, allocatable :: k(:)

    allocate(k(size(res, 4)))
    k = (/(i, i=1,inZone%nSteps, stride)/)

    if (variable == 'U') then
        do i=1,size(k)
            res(:, :, :, i) = inZone%Kinematics(k(i))%U
        end do
    else if (variable == 'V') then
        do i=1,size(k)
            res(:, :, :, i) = inZone%Kinematics(k(i))%V
        end do
    else if (variable == 'W') then
        do i=1,size(k)
            res(:, :, :, i) = inZone%Kinematics(k(i))%W
        end do

    else if (variable == 'Ux') then
        do i=1,size(k)
            res(:, :, :, i) = inZone%Kinematics(k(i))%Ux
        end do
    else if (variable == 'Uy') then
        do i=1,size(k)
            res(:, :, :, i) = inZone%Kinematics(k(i))%Uy
        end do
    else if (variable == 'Uz') then
        do i=1,size(k)
            res(:, :, :, i)= inZone%Kinematics(k(i))%Uz
        end do
    else if (variable == 'Vx') then
        do i=1,size(k)
            res(:, :, :, i)= inZone%Kinematics(k(i))%Uy
        end do
    else if (variable == 'Vy') then
        do i=1,size(k)
            res(:, :, :, i)= inZone%Kinematics(k(i))%Vy
        end do
    else if (variable == 'Vz') then
        do i=1,size(k)
            res(:, :, :, i)= inZone%Kinematics(k(i))%Vz
        end do
    else if (variable == 'Wx') then
        do i=1,size(k)
            res(:, :, :, i)= inZone%Kinematics(k(i))%Uz
        end do
    else if (variable == 'Wy') then
        do i=1,size(k)
            res(:, :, :, i)= inZone%Kinematics(k(i))%Vz
        end do
    else if (variable == 'Wz') then
        do i=1,size(k)
            res(:, :, :, i)= inZone%Kinematics(k(i))%Wz
        end do

    else if (variable == 'X') then
        do i=1,size(k)
            res(:, :, :, i)= inZone%Kinematics(k(i))%X
        end do

    else if (variable == 'Y') then
        do i=1,size(k)
            res(:, :, :, i)= inZone%Kinematics(k(i))%X
        end do

    else if (variable == 'Z') then
        do i=1,size(k)
            res(:, :, :, i)= inZone%Kinematics(k(i))%Z
        end do


    else if (variable == 'Ut') then
        do i=1,size(k)
            res(:, :, :, i) = inZone%Kinematics(k(i))%Ut
        end do
    else if (variable == 'Vt') then
        do i=1,size(k)
            res(:, :, :, i) = inZone%Kinematics(k(i))%Vt
        end do
    else if (variable == 'Wt') then
        do i=1,size(k)
            res(:, :, :, i) = inZone%Kinematics(k(i))%Wt
        end do
    else if (variable == 'phi') then
        do i=1,size(k)
            res(:, :, :, i) = inZone%Kinematics(k(i))%phi
        end do
    else if (variable == 'pdyn') then
        do i=1,size(k)
            res(:, :, :, i) = inZone%Kinematics(k(i))%pdyn
        end do        
    end if 

    deallocate(k)

end subroutine

end module kinematicsArray