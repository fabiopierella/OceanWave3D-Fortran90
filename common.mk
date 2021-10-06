# For inclusion in makefile for use with GNU make (gmake)
# 
# Purpose: Modify to local configuration by user.
#

# Program name
PROGNAME = OceanWave3D-Fortran-Release
LIBNAME  = libOceanWave3D_botp.so

# Installation directory
INSTALLDIR = $(HOME)/bin
LIBINSTALLDIR = $(HOME)/lib

# Build directory where object files are stored 
BUILDDIR = $(PWD)/buildRelease

FC = gfortran
PLATFORM = sophia
USER = jess
# To profile the code: add -pg in OPTFLAGS / DBFLAGS
# This will output gmon.out
# Then you can run gprof to make the results human-readable

ifeq ($(PLATFORM),gbar)
  LINLIB = 
  INCLUDEDIRS = -I/zhome/08/f/75737/Downloads/hdf5-1.10.5/hdf5/include
  INCLUDEDIRS = -I/zhome/08/f/75737/Downloads/hdf5-1.10.5/hdf5/include  
endif
ifeq ($(PLATFORM),jess)
  LINLIB = -lma41
  INCLUDEDIRS = -I$(EBROOTHDF5)/include -gdwarf-2
  # INCLUDEDIRS = -I$(EBROOTHDF5)/include
endif
ifeq ($(PLATFORM),sophia)
  LINLIB = -lma41
  INCLUDEDIRS = -I$(EBROOTHDF5)/include -gdwarf-2
  # INCLUDEDIRS = -I$(EBROOTHDF5)/include
endif

ifeq ($(FC),gfortran)
  # fabpi machine, gfortran
  LIBDIRS  = -L$(HOME)/lib/ -L$(MODULE_OPENBLAS_LIB_DIR) 
  LIBDIRS  += -L/zhome/08/f/75737/Downloads/hdf5-1.10.5/hdf5/lib/
  LINLIB   += -lskit -lopenblas -lhdf5 -lhdf5_fortran -lhdf5_hl 
  DBFLAGS  += -g -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none  -fno-automatic
  OPTFLAGS += -O3 -ffree-line-length-none -fno-automatic -ffpe-trap=invalid,zero,overflow 
endif

ifeq ($(FC),ifort)
  # hbb work machine with intel compiler
  LIBDIRS  = -L$(HOME)/lib/  -L$(MODULE_OPENBLAS_LIB_DIR)
  #  LINLIB   = -mkl -lharwell_intel -lskit_intel -llapack_intel -lblas_intel
  LINLIB   += -lopenblas -lma27 -lma41 
  DBFLAGS  += -g -CB -fpe0 -fpstkchk  -traceback 
  OPTFLAGS += -O3 -tpp7
endif

ifeq ($(USER),fabpi)
  LIBDIRS  = -L$(HOME)/lib/  -L/usr/lib/x86_64-linux-gnu/hdf5/serial
  LINLIB   = -lopenblas -lskit -lhdf5 -lhdf5_fortran -lhdf5_hl
  DBFLAGS  = -g -I/usr/include/hdf5/serial
  OPTFLAGS = -O3 -I/usr/include/hdf5/serial
endif

print-%  : ; @echo $* = $($*)
