# Makefile for reshape_hist
F90 = ifort
FFLAGS = -r8
objects = read_hist_1.o read_hist_2.o read_hist_3.o\
          read_hist_srfc.o \
          write_hist_1.o write_hist_2.o write_hist_3.o\
          write_hist_srfc.o reshape_hist.o
reshape_hist : $(objects)
	$(F90) ${FFLAGS} -o reshape_hist $(objects)
%.o : %.f90
	$(F90) ${FFLAGS} -c $<
clean:
	rm -f reshape_hist *.o *.mod

