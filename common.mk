# For inclusion in makefile for use with GNU make (gmake)
#
# Purpose: Modify to local configuration by user.
#

# Program name
PROGNAME = OceanWave3D-Fortran
LIBNAME  = libOceanWave3D.so

# Installation directory
INSTALLDIR = $(HOME)/bin
LIBINSTALLDIR = $(HOME)/lib

# Build directory where object files are stored
BUILDDIR = $(PWD)/build

FC = gfortran
PLATFORM = wsl

ifeq ($(PLATFORM),gbar)
  INCLUDEDIRS = -I/usr/lib64/gfortran/modules
  DBFLAGS  = -I/usr/lib64/gfortran/modules
  OPTFLAGS = -I/usr/lib64/gfortran/modules
endif

ifeq ($(PLATFORM),wsl)
  # WSL: install HDF5 libraries with
  #   sudo apt install libhdf5-dev libhdf5-fortran-dev libopenblas-dev
  INCLUDEDIRS = -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include
endif

ifeq ($(FC),gfortran)
  # WSL with gfortran
  LIBDIRS  = -L$(HOME)/lib/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial
  LINLIB   += -lskit -lopenblas -lhdf5 -lhdf5_fortran -lhdf5_hl
  DBFLAGS  += -pg -g -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fno-automatic
  OPTFLAGS += -pg -O3 -ffree-line-length-none -fno-automatic -ffpe-trap=invalid,zero,overflow
  SHLIBFLAGS = -shared -fPIC -O3 -ffree-line-length-none -fno-automatic
endif

ifeq ($(FC),ifort)
  # TODO: recheck with modern intel fortran compiler
  LIBDIRS  = -L$(HOME)/lib/ -L$(MODULE_OPENBLAS_LIB_DIR)
  LINLIB   += -lskit -lopenblas
  DBFLAGS  += -g -CB -fpe0 -fpstkchk -traceback
  OPTFLAGS += -O3 -tpp7
  SHLIBFLAGS = -shared -fPIC -O3
endif

print-%  : ; @echo $* = $($*)
