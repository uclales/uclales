###############################################################
#
# Type of code to build $CNFG, location of code (in $ROOT) and 
# location where model is to be built $BIN
#
ROOT      = ..
BIN       = $(ROOT)/bin
ARCH      := $(shell uname)
#
# Note that on the Linux systems of any mpipc, the software environment
# is such that it requires the additional netcdf flag -l netcdff, in order
# for all netcdf libraries to be found
#
#
#
SRC     =$(ROOT)/src
ECHO    = /bin/echo
RM      = /bin/rm -f

ifeq ($(ARCH),Linux)
  COMPILER  := 'GFORTRAN'
  ARCHIVE   = ar rs
  RANLIB    =:
  ifeq ($(COMPILER),'NAG')
     F90      = /sw/etch-ia32/NAG/NAGWare-5.1.365/bin/nagf95
     NCDF     = /sw/etch-ia32/netcdf-3.6.2-nag
     DEBUG    = -C=all -gline -nan
     FFLAGS   = $(DEBUG) -O3 -kind=byte $(NCDFINC) -Oassumed=always_contig
     SEQFFLAGS = -I$(SRC)/seq -mdir $(SRC)/seq
     MPIFFLAGS = -I$(SRC)/mpi -mdir $(SRC)/mpi
     F77FLAGS = -fixed -dusty
  endif
  ifeq ($(COMPILER),'G95')
     F90      = g95
     NCDF     = /sw/g95/netcdf
     FFLAGS   = -O3 -r8   $(NCDFINC)
     SEQFFLAGS = -I$(SRC)/seq -fmod=$(SRC)/seq
     MPIFFLAGS = -I$(SRC)/mpi -fmod=$(SRC)/mpi
     F77FLAGS = -ffixed-form
  endif
  ifeq ($(COMPILER),'GFORTRAN')
     F90      = gfortran
     MPIF90   = mpif90
     NCDF     = /usr
     FFLAGS   = -W -Wall -fdefault-real-8 -ffree-line-length-none -O3  $(NCDFINC)
#      SEQFFLAGS = -I$(SRC)/seq -fmod=$(SRC)/seq
#      MPIFFLAGS = 
     F77FLAGS = -ffixed-form
  endif
  ifeq ($(COMPILER),'IFORT')
     F90      = ifort
     MPIF90    = /sw/etch-ia32/mpich2-1.0.8-intel/bin/mpif90
     NCDF     = /sw/etch-ia32/netcdf-3.6.2-intel
     SEQFFLAGS = -I$(SRC)/seq -module $(SRC)/seq
     MPIFFLAGS = -I$(SRC)/mpi -module $(SRC)/mpi
     FFLAGS   = -O3 -r8 -diag-disable remark $(NCDFINC)
     F77FLAGS = -fixed
  endif
  NCDFLIB   = -L$(NCDF)/lib64  -lnetcdff -lnetcdf 
  NCDFINC   = -I$(NCDF)/include 
  LIBS      = $(NCDFLIB)
  MPILIBS   = $(LIBS)
endif


LES_ARC_MPI=$(BIN)/libles_mpi.a
LES_OUT_MPI=$(BIN)/les.mpi

LES_ARC_SEQ=$(BIN)/libles_seq.a
LES_OUT_SEQ=$(BIN)/les.seq

default: mpi seq

all:  mpi seq

seq: $(LES_OUT_SEQ)

mpi: $(LES_OUT_MPI)

$(LES_OUT_SEQ): $(SRC)/LES.f90 $(LES_ARC_SEQ)  
	$(F90) $(FFLAGS) $(SEQFFLAGS) -o $@ $< $(LES_ARC_SEQ) $(LIBS)

$(LES_OUT_MPI): $(SRC)/LES.f90 $(LES_ARC_MPI)  
	$(MPIF90) $(FFLAGS) $(MPIFFLAGS) -o $@ $< $(LES_ARC_MPI) $(MPILIBS)

$(LES_ARC_SEQ): FORCE
	cd $(SRC) ; $(MAKE) LES_ARC=$(LES_ARC_SEQ) \
        FFLAGS='$(FFLAGS) $(SEQFFLAGS)' F90=$(F90) ARCHIVE='$(ARCHIVE)' \
	RANLIB=$(RANLIB) F77FLAGS='$(F77FLAGS)'

$(LES_ARC_MPI): FORCE
	cd $(SRC) ; $(MAKE) LES_ARC=$(LES_ARC_MPI) \
        FFLAGS='$(FFLAGS) $(MPIFFLAGS)' F90='$(MPIF90)' ARCHIVE='$(ARCHIVE)' \
	RANLIB='$(RANLIB)' F77FLAGS='$(F77FLAGS)'
#
# cleaning
# --------------------
#
clean: cleanmpi cleanseq 
	$(RM) $(SRC)/*mod $(SRC)/*.o

cleanmpi:
	$(ECHO) "cleaning mpi model"
	$(RM) core $(LES_OUT_MPI) $(SRC)/mpi/*mod $(LES_ARC_MPI)

cleanseq:
	$(ECHO) "clean sequential model"
	$(RM) core $(LES_OUT_SEQ) $(SRC)/seq/*mod $(LES_ARC_SEQ)

FORCE: 
.PRECIOUS: $(LIBS)
