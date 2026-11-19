#=====================================================
# Run directory file
#---------------------------------------------------
RUN = ../run/
#
# Precompiler flags
#---------------------------------------------------
ARCH = process-gfortran
PARA = 0
#=====================================================


#############################################################################
ifeq ($(ARCH),macOS)
#---------------------------------------------------
# Default (macOS X)
F90      = gfortran
CPPFLAGS = -x f95-cpp-input -DSYSTEM
FLAGS    = -fdefault-double-8 -fdefault-real-8 -ffree-line-length-0 -O3
#FLAGS    = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0     \
       -fbounds-check -ffpe-trap=invalid\
       -Wall -g -fcheck=all -fbacktrace -Warray-bounds

# Libraries
LIBS     = -llapack -lblas -lcurl -lnetcdff -lnetcdf 
NCDFINC  = -I/opt/local/include/
NCDFLIB  = -L/opt/local/lib/
#---------------------------------------------------
else ifeq ($(ARCH),dirac_mono)
#---------------------------------------------------
# Parallel on dirac
F90      = ifort
CPPFLAGS = -cpp
FLAGS    = -r8 -O3
#FLAGS    = -traceback -debug all -check all -warn all
# Libraries
# LAPACK_PATH = /cm/shared/apps/ics/mkl/lib/intel64/
LAPACK_PATH = $(MKLROOT)/lib/intel64/
LIBS        = $(LAPACK_PATH)libmkl_intel_lp64.so \
	      $(LAPACK_PATH)libmkl_sequential.so \
	      $(LAPACK_PATH)libmkl_core.so       \
		     -lcurl -lnetcdff -lnetcdf 
NCDFINC  = -I/cm/shared/apps/netcdf/intel/4.2.1.1/include/
NCDFLIB  = -L/cm/shared/apps/netcdf/intel/4.2.1.1/lib/
#---------------------------------------------------
else ifeq ($(ARCH),dirac_mono_gfortran)
#---------------------------------------------------
# Parallel on dirac
F90      = gfortran
CPPFLAGS = -cpp
#FLAGS    = -r8 -O3
FLAGS    = -fdefault-double-8 -fdefault-real-8 -ffree-line-length-0 -O3
# Libraries
# LAPACK_PATH = /cm/shared/apps/ics/mkl/lib/intel64/
LAPACK_PATH = $(MKLROOT)/lib/intel64/
LIBS        = $(LAPACK_PATH)libmkl_intel_lp64.so \
              $(LAPACK_PATH)libmkl_sequential.so \
              $(LAPACK_PATH)libmkl_core.so       \
              /panfs/dirac.panfs.cluster/home/dc-trem1/pnetcdf/lib/libpnetcdf.a \
                     -lcurl -lnetcdff -lnetcdf
NCDFINC  = -I/cm/shared/apps/netcdf/intel/4.2.1.1/include/
NCDFLIB  = -L/cm/shared/apps/netcdf/intel/4.2.1.1/lib/
#---------------------------------------------------
else ifeq ($(ARCH),zen_mono)
#---------------------------------------------------
# Parallel on dirac
F90      = ifort
CPPFLAGS = -cpp
FLAGS    = -r8 -O3
#FLAGS    = -traceback -debug all -check all -warn all
# Libraries
LAPACK_PATH = /sw/sdev/composer_xe_2013.1.117/mkl/lib/intel64/
LIBS        = $(LAPACK_PATH)libmkl_intel_lp64.so \
	      $(LAPACK_PATH)libmkl_sequential.so \
	      $(LAPACK_PATH)libmkl_core.so       \
                 -lcurl -lnetcdf 
NCDFINC     = -I/sw/netcdf-3.6.2/include/
NCDFLIB     = -L/sw/netcdf-3.6.2/lib/ 

#---------------------------------------------------
else ifeq ($(ARCH),isca_mono)
#---------------------------------------------------
# Normal on  isca
F90      = gfortran
CPPFLAGS = -x f95-cpp-input -DSYSTEM
FLAGS    = -fdefault-double-8 -fdefault-real-8 -ffree-line-length-0 -O3
#FLAGS    = -traceback -debug all -check all -warn all
# Libraries
LAPACK_PATH = /gpfs/ts0/shared/software/lapack/3.6.1/gnu/4.8.5/lib/

LIBS        = $(LAPACK_PATH)liblapack.a \
              $(LAPACK_PATH)libtmglib.a \
              $(LAPACK_PATH)libblas.a \
                 -lnetcdff -lnetcdf -lhdf5 -lhdf5_hl  -lsz
NCDFINC     = -I/gpfs/ts0/shared/software/netCDF-Fortran/4.4.3-foss-2016a/include/
NCDFLIB     = -L/gpfs/ts0/shared/software/netCDF-Fortran/4.4.3-foss-2016a/lib/ -L/gpfs/ts0/shared/software/netCDF/4.4.0-foss-2016a/lib64 -L/gpfs/ts0/shared/software/HDF5/1.8.17-foss-2016a/lib -L/gpfs/ts0/shared/software/Szip/2.1-foss-2016a/lib

#-------------------------------------------------
else ifeq ($(ARCH),zen_para)
#---------------------------------------------------
# Parallel on dirac
F90      = mpiifort
CPPFLAGS = -cpp
FLAGS    = -r8 -O3
#FLAGS    = -traceback -debug all -check all -warn all
# Libraries
LAPACK_PATH = /sw/sdev/composer_xe_2013.1.117/mkl/lib/intel64/
LIBS        = $(LAPACK_PATH)libmkl_intel_lp64.so \
              $(LAPACK_PATH)libmkl_sequential.so \
              $(LAPACK_PATH)libmkl_core.so       \
                 -lmpi -lcurl -lnetcdf
NCDFINC     = -I/sw/netcdf-3.6.2/include/
NCDFLIB     = -L/sw/netcdf-3.6.2/lib/



#---------------------------------------------------
else ifeq ($(ARCH),dirac_para)
#---------------------------------------------------
# Parallel on dirac
F90      = mpiifort
CPPFLAGS = -cpp
FLAGS    = -r8 -O3 -diag-disable remark
#FLAGS    = -r8 -check all -traceback -debug all
# Libraries
# LAPACK_PATH = /cm/shared/apps/ics/mkl/lib/intel64/
LAPACK_PATH = $(MKLROOT)/lib/intel64/
LIBS        = $(LAPACK_PATH)libmkl_intel_lp64.so \
	      $(LAPACK_PATH)libmkl_sequential.so \
	      $(LAPACK_PATH)libmkl_core.so       \
	      /cm/shared/apps/netcdf_par/intel/4.6.1/lib/libnetcdf.a \
                 -lmpi -lcurl -lnetcdff -lnetcdf 
NCDFINC     = -I/cm/shared/apps/netcdf/intel/4.2.1.1/include/
NCDFLIB     = -L/cm/shared/apps/netcdf/intel/4.2.1.1/lib/ 

#---------------------------------------------------
else ifeq ($(ARCH),post-zen-gfortran)
#---------------------------------------------------
# Non-parallel on post-zen
F90      = gfortran
CPPFLAGS = -x f95-cpp-input -DSYSTEM
FLAGS    = -fdefault-double-8 -fdefault-real-8 -ffree-line-length-0 -O3
#FLAGS    = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0     \
       -fbounds-check -ffpe-trap=invalid\
       -Wall -g -fcheck=all -fbacktrace -Warray-bounds
# Libraries
LAPACK_PATH = /home/dsa206/Library/lapack-3.5.0_gfortran/
BLAS_PATH = /home/dsa206/Library/BLAS_gfortran/
LIBS        = $(LAPACK_PATH)liblapack.a \
              $(LAPACK_PATH)libtmglib.a \
              $(BLAS_PATH)blas_LINUX.a \
                 -lnetcdff -lnetcdf
NCDFINC     = -I/home/dsa206/Library/netCDF_gfortran/include/
NCDFLIB     = -L/home/dsa206/Library/netCDF_gfortran/lib/
#---------------------------------------------------
else ifeq ($(ARCH),post-zen-ifort)
#---------------------------------------------------
# Non-parallel on post-zen
F90      = ifort
CPPFLAGS = -cpp
FLAGS = -r8 -O3
#FLAGS    = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0     \
       -fbounds-check -ffpe-trap=invalid\
       -Wall -g -fcheck=all -fbacktrace -Warray-bounds
# Libraries
LAPACK_PATH = /home/dsa206/Library/lapack-3.5.0_ifort/
BLAS_PATH = /home/dsa206/Library/BLAS_ifort/
LIBS        = $(LAPACK_PATH)liblapack.a \
	      $(LAPACK_PATH)libtmglib.a \
	      $(BLAS_PATH)blas_LINUX.a \
                 -lnetcdff -lnetcdf
NCDFINC     = -I/home/dsa206/Library/netCDF_ifort/include/
NCDFLIB     = -L/home/dsa206/Library/netCDF_ifort/lib/
#---------------------------------------------------
else ifeq ($(ARCH),um-serve-gfortran)
#---------------------------------------------------
# Non-parallel on post-zen
F90      = gfortran
CPPFLAGS = -x f95-cpp-input -DSYSTEM
FLAGS    = -fdefault-double-8 -fdefault-real-8 -ffree-line-length-O3 -fallow-argument-mismatch
#FLAGS    = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0     \
       -fbounds-check -ffpe-trap=invalid\
       -Wall -g -fcheck=all -fbacktrace -Warray-bounds
# Libraries
LIBS        =  -llapack -lblas -lnetcdff -lnetcdf
NCDFINC     = -I/home/bd257/Library/netCDF/include/
NCDFLIB     = -L/home/bd257/Library/netCDF/lib/
#---------------------------------------------------
else ifeq ($(ARCH),process-gfortran)
#---------------------------------------------------
# Non-parallel on post-zen
F90      = gfortran
CPPFLAGS = -x f95-cpp-input -DSYSTEM
FLAGS    = -fdefault-double-8 -fdefault-real-8 -ffree-line-length-0 -O3 -fallow-argument-mismatch
#FLAGS    = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0     \
       -fbounds-check -ffpe-trap=invalid\
       -Wall -g -fcheck=all -fbacktrace -Warray-bounds
# Libraries
LIBS        =  -llapack -lblas -lnetcdff -lnetcdf
NCDFINC     = -I/usr/include/
NCDFLIB     = -L/usr/lib/
#---------------------------------------------------
else ifeq ($(ARCH),gitlab_CI)
#---------------------------------------------------
# Non-parallel on gitlab_CI
F90      = gfortran
CPPFLAGS = -x f95-cpp-input -DSYSTEM
#FLAGS    = -fdefault-integer-8 -fdefault-double-8 -fdefault-real-8 -ffree-line-length-0 -O3
FLAGS    = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0     \
       -fbounds-check -ffpe-trap=invalid\
       -Wall -g -fcheck=all -fbacktrace -Warray-bounds

# Libraries
LIBS     = -llapack -lblas -lnetcdff
NCDFINC  = -I/usr/include/
NCDFLIB  = -L/usr/include/
#---------------------------------------------------
endif
#---------------------------------------------------
# directories and vpath for source files
#---------------------------------------------------
SRCDIR = ../src/
ODEDIR = ../src/odepack/
CASDIR = $(RUN)cas/
VPATH  = $(CASDIR):$(SRCDIR):$(ODEDIR)
#---------------------------------------------------
# Source objects
#---------------------------------------------------
SRCOBJ = param.o constant.o util.o opkda2.o opkda1.o opkdmain.o grid.o \
radtrans.o chemistry.o eos.o opacity.o convection.o solver.o chem_neq.o \
cs2006_relaxation.o tsai2017_relaxation.o chem_relax.o

#############################################################################

atmo.x:		$(SRCOBJ) atmo.o
		$(F90) $(FLAGS) $(SRCOBJ) atmo.o -o \
		$(RUN)atmo.x $(NCDFLIB) $(LIBS) -DPARA=$(PARA)

#############################################################################

%.o:		%.f90
		$(F90) $(NCDFINC) $(CPPFLAGS) $(FLAGS) -DPARA=$(PARA) -c $< -o $@

%.o:		%.f
		$(F90) $(NCDFINC) $(CPPFLAGS) $(FLAGS) -DPARA=$(PARA) -c $< -o $@

atmo.o:	 	$(SRCDIR)atmo.f90 Makefile
		$(F90) $(CPPFLAGS) $(FLAGS) -DPARA=$(PARA) -c $(SRCDIR)atmo.f90

#############################################################################

clean :
		rm *.o *.mod

#############################################################################




