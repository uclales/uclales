# MAKEFILE for merge_particles.f90
# Assumes using gcc-compilers
# ---------------------------------

## MPI desktops:
#FC = gfortran
#FCFLAGS = -g -O3 -fdefault-real-8 -funroll-all-loops
##FCFLAGS = -g -O0 -fdefault-real-8 -finit-real=nan -Wall -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
#INCDIRS = -I/sw/squeeze-x64/netcdf-latest-static-gcc45/include
#LIBS =  -L/sw/squeeze-x64/netcdf-latest-static-gcc45/lib -L/sw/squeeze-x64/szip-2.1-static/lib -L/sw/squeeze-x64/hdf5-1.8.7-static/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lm -lcurl  -lsz -lm

# Blizzard:
#FC = xlf90
#FCFLAGS = -g -O3 
#INCDIRS = -I/sw/aix61/netcdf-4.1.2-hdf5-threadsafe/include
#LIBS =  -L/sw/aix61/netcdf-4.1.2-hdf5-threadsafe/lib -lnetcdff -lnetcdf -L/sw/aix61/hdf5-1.8.6-threadsafe/lib -lhdf5_hl -lhdf5  -L/sw/aix61/zlib-1.2.3/lib -lz -lm -lhdf5_hl -lhdf5 -lz  -L/sw/aix61/szip-2.1/lib -lsz -lhdf5_hl -lhdf5 -lz 

########################################################
# Thunder using Intel:
########################################################
FC = ifort
#FCFLAGS = -O3 -no-prec-div -xHost -fp-model source -traceback -r8 -ftz -extend_source
FCFLAGS = -O3 -no-prec-div -xSSE4.2 -fp-model source -traceback -r8 -ftz -extend_source

#FCFLAGS = -traceback -r8 -ftz -extend_source -fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created
INCDIRS = -I/sw/squeeze-x64/netcdf_fortran-4.2-static-intel13/include
LIBS = -L/sw/squeeze-x64/netcdf_fortran-4.2-static-intel13/lib -lnetcdff -L/sw/squeeze-x64/szip-2.1-static/lib -L/sw/squeeze-x64/hdf5-1.8.8-static/lib -L/sw/squeeze-x64/netcdf-4.2.1.1-static/lib -lnetcdf -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz -lm 


PROGRAMS = merge_drops
all: $(PROGRAMS)

merge_drops: merge_drops.o 
	$(FC) $(FCFLAGS) -o merge_drops merge_drops.o  $(LIBS) 

merge_drops.o: merge_drops.f90
	$(FC) $(FCFLAGS) -c merge_drops.f90 $(INCDIRS)

# Utility targets
.PHONY: clean

clean:
	rm -f *.o
	rm -f *~ $(PROGRAMS)

