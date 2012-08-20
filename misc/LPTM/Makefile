# MAKEFILE for merge_particles.f90
# Assumes using gcc-compilers
# ---------------------------------

# Production
FCFLAGS = -g -O3 -fdefault-real-8 -funroll-all-loops
# Debug
#FCFLAGS = -g -O0 -fdefault-real-8 -finit-real=nan -Wall -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow

# MPIM-PCs
#INCDIRS = -I/sw/squeeze-x64/netcdf-latest-static-gcc45/include
#LIBDIRS =  -L/sw/squeeze-x64/netcdf-latest-static-gcc45/lib -L/sw/squeeze-x64/szip-2.1-static/lib -L/sw/squeeze-x64/hdf5-1.8.7-static/lib 
#LIBS = -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lm -lcurl  -lsz -lm 

# BLIZZARD 
INCDIRS = -I/sw/aix61/netcdf-4.1.2-hdf5-threadsafe/include
LIBDIRS =  -L/sw/aix61/netcdf-4.1.2-hdf5-threadsafe/lib -L/sw/aix61/zlib-1.2.3/lib -L/sw/aix61/hdf5-1.8.6-threadsafe/lib 
LIBS = -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lm -lcurl  -lsz -lm 

PROGRAMS = merge_particles
all: $(PROGRAMS)

merge_particles: merge_particles.o 
	gfortran $(FCFLAGS) -o merge_particles merge_particles.o $(LIBDIRS) $(LIBS) 

merge_particles.o: merge_particles.f90
	gfortran $(FCFLAGS) -c merge_particles.f90 $(INCDIRS)

# Utility targets
.PHONY: clean

clean:
	rm -f *.o
	rm -f *~ $(PROGRAMS)

