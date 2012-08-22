# MAKEFILE for merge_particles.f90
# Assumes using gcc-compilers
# ---------------------------------

# MPI desktops:
FC = gfortran
FCFLAGS = -g -O3 -fdefault-real-8 -funroll-all-loops
#FCFLAGS = -g -O0 -fdefault-real-8 -finit-real=nan -Wall -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
INCDIRS = -I/sw/squeeze-x64/netcdf-latest-static-gcc45/include
LIBS =  -L/sw/squeeze-x64/netcdf-latest-static-gcc45/lib -L/sw/squeeze-x64/szip-2.1-static/lib -L/sw/squeeze-x64/hdf5-1.8.7-static/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lm -lcurl  -lsz -lm

# Blizzard:
#FC = xlf90
#FCFLAGS = -g -O3 
#INCDIRS = -I/sw/aix61/netcdf-4.1.2-hdf5-threadsafe/include
#LIBS =  -L/sw/aix61/netcdf-4.1.2-hdf5-threadsafe/lib -lnetcdff -lnetcdf -L/sw/aix61/hdf5-1.8.6-threadsafe/lib -lhdf5_hl -lhdf5  -L/sw/aix61/zlib-1.2.3/lib -lz -lm -lhdf5_hl -lhdf5 -lz  -L/sw/aix61/szip-2.1/lib -lsz -lhdf5_hl -lhdf5 -lz 


PROGRAMS = merge_particles
all: $(PROGRAMS)

merge_particles: merge_particles.o 
	$(FC) $(FCFLAGS) -o merge_particles merge_particles.o  $(LIBS) 

merge_particles.o: merge_particles.f90
	$(FC) $(FCFLAGS) -c merge_particles.f90 $(INCDIRS)

# Utility targets
.PHONY: clean

clean:
	rm -f *.o
	rm -f *~ $(PROGRAMS)

