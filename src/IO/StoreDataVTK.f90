SUBROUTINE StoreDataVTK(nx,ny,E,P,FineGrid,nr,formattype)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
USE Lib_VTK_IO
IMPLICIT NONE
INTEGER(4) :: nx, ny, i, j, nr, formattype
REAL(KIND=long), DIMENSION(nx,ny) :: E, P
TYPE (Level_def) :: FineGrid
CHARACTER(len=30) :: filename
INTEGER(4) :: E_IO, nz, ntotal
REAL(8), ALLOCATABLE :: X(:), Y(:), Z(:)

nz = 1
ntotal = nx * ny

ALLOCATE(X(ntotal), Y(ntotal), Z(ntotal))

! Fill coordinate arrays from the grid
DO j = 1, ny
   DO i = 1, nx
      X((j-1)*nx + i) = FineGrid%x(i,j)
      Y((j-1)*nx + i) = FineGrid%y(i,j)
      Z((j-1)*nx + i) = zero
   END DO
END DO

WRITE(unit=filename, FMT="(A,I5.5,A)") "E_", nr, ".vtk"

! Write VTK Legacy structured grid file
E_IO = VTK_INI('binary', TRIM(filename), 'Free surface Field', 'STRUCTURED_GRID')
E_IO = VTK_GEO(nx, ny, nz, ntotal, X, Y, Z)
E_IO = VTK_DAT(ntotal, 'node')
E_IO = VTK_VAR(ntotal, 'Elevation', RESHAPE(E, (/ntotal/)))
E_IO = VTK_VAR(ntotal, 'Potential', RESHAPE(P, (/ntotal/)))
E_IO = VTK_END()

DEALLOCATE(X, Y, Z)

END SUBROUTINE StoreDataVTK
