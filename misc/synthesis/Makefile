FC:=ftn
FFLAGS=`nc-config --fflags`
LIBS=`nc-config --libs`

ALL=$(patsubst %.f90, %, $(wildcard *.f90))

all: $(ALL)

%: %.f90
	$(FC) $(FFLAGS) $(LIBS) $< -o $@


