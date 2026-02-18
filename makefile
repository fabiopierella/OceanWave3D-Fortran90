# makefile for use with GNU make (gmake)
#
# Author: Allan P. Engsig-Karup.
#

include common.mk

.SUFFIXES: .f90 .f

# Source files
# subdirectories - complete directory-tree should be included
SOURCES :=
include src/variabledefs/makefile.inc
include src/utilities/makefile.inc
include src/analyticalsolutions/makefile.inc
include src/curvilinear/makefile.inc
include src/fft/makefile.inc
include src/functions/makefile.inc
include src/multigrid/makefile.inc
include src/pressure/makefile.inc
include src/timeintegration/makefile.inc
include src/wrappers/makefile.inc
include src/initialization/makefile.inc
include src/IO/makefile.inc
include src/main/makefile.inc
include src/iterative/makefile.inc
include src/OpenFoam/makefile.inc

# Search paths for source dependencies
VPATH  = src/variabledefs
VPATH += src/analyticalsolutions
VPATH += src/curvilinear
VPATH += src/fft
VPATH += src/functions
VPATH += src/initialization
VPATH += src/IO
VPATH += src/main
VPATH += src/multigrid
VPATH += src/pressure
VPATH += src/timeintegration
VPATH += src/utilities
VPATH += src/wrappers
VPATH += src/iterative
VPATH += src/OpenFoam
VPATH += src/OpenFoam/IO
VPATH += src/OpenFoam/OpenFoam
VPATH += $(BUILDDIR)

# Include directories
INCLUDEDIRS += -I$(BUILDDIR)

# Object files
OBJECTS := $(SOURCES:.f=.o)
OBJECTS := $(OBJECTS:.f90=.o)
OBJECTSNODIR = $(notdir $(OBJECTS))
OBJECTSBUILDDIR = $(addprefix $(BUILDDIR)/,$(OBJECTSNODIR))

# Third-party .mod files to copy into build dir
VTKMOD = thirdpartylibs/LIB_VTK_IO/static/lib_vtk_io.mod

# Ensure build directory exists (order-only prerequisite)
$(BUILDDIR):
	@mkdir -p $(BUILDDIR)
	@if [ -f "$(VTKMOD)" ]; then cp $(VTKMOD) $(BUILDDIR)/; fi

# Ensure install directory exists
$(INSTALLDIR):
	@mkdir -p $(INSTALLDIR)

# Default target
default:
	@echo "To build OceanWave3D type 'make Release' or 'make Debug'"
	@echo "To clean build artifacts type 'make clean'"

# Targets for linking
Release: FFLAGS = $(OPTFLAGS)
Release: $(INSTALLDIR) $(OBJECTSBUILDDIR)
	@echo "*** Linking OceanWave3D (Release)... ***"
	$(FC) $(FFLAGS) -o $(INSTALLDIR)/$(PROGNAME) $(OBJECTSBUILDDIR) $(LIBDIRS) $(LINLIB) $(INCLUDEDIRS)
	@echo "OceanWave3D has been built successfully."

Debug: FFLAGS = $(DBFLAGS)
Debug: $(INSTALLDIR) $(OBJECTSBUILDDIR)
	@echo "*** Linking OceanWave3D (Debug)... ***"
	$(FC) $(FFLAGS) -o $(INSTALLDIR)/$(PROGNAME) $(OBJECTSBUILDDIR) $(LIBDIRS) $(LINLIB) $(INCLUDEDIRS)
	@echo "OceanWave3D has been built successfully (Debug)."

all: Release

sharedLib: FFLAGS = $(SHLIBFLAGS)
sharedLib: $(OBJECTSBUILDDIR)
	@echo "*** Linking OceanWave3D shared library... ***"
	$(FC) $(FFLAGS) -o $(LIBINSTALLDIR)/$(LIBNAME) $(OBJECTSBUILDDIR) $(LIBDIRS) $(LINLIB) $(INCLUDEDIRS)
	@echo "Shared library for OceanWave3D has been built successfully."

# Compile only
compile: $(OBJECTSBUILDDIR)

clean:
	rm -rf $(BUILDDIR)
	rm -f $(INSTALLDIR)/$(PROGNAME)

# Move any .mod files generated in the source root into the build directory
define move_mods
	@if ls *.mod 1>/dev/null 2>&1; then mv *.mod $(BUILDDIR)/; fi
endef

# Generic compilation rules (order-only prerequisite on BUILDDIR)
$(BUILDDIR)/%.o: %.f | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDEDIRS)
	$(move_mods)

$(BUILDDIR)/%.o: %.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDEDIRS)
	$(move_mods)
