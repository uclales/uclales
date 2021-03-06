###############################################################
#
# Type of code to build $CNFG, location of code (in $ROOT) and 
# location where model is to be built $BIN
#
ROOT      = ..
BIN       = $(ROOT)/bin
ARCH      := $(shell uname)
#
# This Makefile can be used for the IBM compiler on blizzard
#
#
SRC     =$(ROOT)/src
ECHO    = /bin/echo
RM      = /bin/rm -f

ifeq ($(ARCH),AIX)
  NCDF = /sw/aix61/netcdf-4.1.1-rc1
  NCDFLIB = -L$(NCDF)/lib -lnetcdf
  NCDFINC = -I$(NCDF)/include
  HDF5 = /sw/aix61/hdf5-1.8.4-patch1
  HDF5LIB = -L$(HDF5)/lib -lhdf5 -lhdf5_hl
  HDF5INC = -I$(HDF5)/include
  SZIP = /sw/aix53/szip-2.1
  SZIPLIB = -L$(SZIP)/lib -lsz
  SZIPINC = -I$(SZIP)/include
  ZLIB = /sw/aix53/zlib-1.2.3
  ZLIBLIB = -L$(ZLIB)/lib -lz
  ZLIBINC = -I$(ZLIB)/include
  ARCHIVE   = ar rs
  RANLIB    =:
  F90       = xlf95_r 
  MPIF90    = mpxlf95_r
  FFLAGS   = -qextname -q64 -qrealsize=8 -O4 -qnoipa -qstrict=none:exceptions $(NCDFINC)
  #FFLAGS    = -qextname -q64 -qrealsize=8 -O0 -qfullpath -qflttrap=overflow:nanq:zerodivide:enable -qsigtrap $(NCDFINC)
  #FFLAGS   = -qextname -qrealsize=8  -qwarn64            -qflttrap=en:ov:zero:inv:imp -qflag=w:e  $(NCDFINC)
  #FFLAGS   = -qfullpath -qflttrap=overflow:nanq:zerodivide:enable -qsigtrap
  MPIFFLAGS = -qmoddir=$(SRC)/mpi -I$(SRC)/mpi
  SEQFFLAGS = -qmoddir=$(SRC)/seq -I$(SRC)/seq
  LIBS      = $(NCDFLIB) $(HDF5LIB) $(SZIPLIB) $(ZLIBLIB) -lmassv -lmass
  MPILIBS   = -lblacssmp -lpesslsmp $(LIBS)
  F77FLAGS  = -qextname -qfixed
endif

LES_ARC_MPI=$(BIN)/libles_mpi.a
LES_OUT_MPI=$(BIN)/les.mpi

LES_ARC_SEQ=$(BIN)/libles_seq.a
LES_OUT_SEQ=$(BIN)/les.seq

default: mpi

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
