# Generated automatically from Makefile.in by configure.
.SUFFIXES : .f .o .a .f90
.PRECIOUS : $(LES_ARC)
#
#     Define model archive file name.
#
RM       = /bin/rm -f
ECHO     = /bin/echo
#
#     Define all necessary source code
#
ARCOBJ_base :=       \
$(LES_ARC)(defs.o)   \
$(LES_ARC)(rad_solver.o) \
$(LES_ARC)(rad_corkds.o) \
$(LES_ARC)(rad_cldwtr.o) \
$(LES_ARC)(rad_rndnmb.o) \
$(LES_ARC)(rad_d4strm.o) \
$(LES_ARC)(rad_driver.o) \
$(LES_ARC)(ncio.o)   \
$(LES_ARC)(util.o)   \
$(LES_ARC)(grid.o)   \
$(LES_ARC)(stat.o)   \
$(LES_ARC)(thrm.o)   \
$(LES_ARC)(forc.o)   \
$(LES_ARC)(advf.o)   \
$(LES_ARC)(advl.o)   \
$(LES_ARC)(srfc.o)   \
$(LES_ARC)(prss.o)   \
$(LES_ARC)(sgsm.o)   \
$(LES_ARC)(mcrp.o)   \
$(LES_ARC)(step.o)   \
$(LES_ARC)(init.o)   \
$(LES_ARC)(rfft.o)   


ifneq (,$(findstring mpi,$(LES_ARC)))
	ARCOBJ_f := $(LES_ARC)(mpi_interface.o) \
                    $(ARCOBJ_base)
else
	ARCOBJ_f := $(LES_ARC)(seq_interface.o) \
                    $(ARCOBJ_base)
endif

.f90.a:
	$(F90) -c $(FFLAGS) $(%:.o=.f90)
	$(ARCHIVE)  $@ $%
	$(RM) $%

.f.a:
	$(F90) -c $(FFLAGS) $(F77FLAGS) $<
	$(ARCHIVE)  $@ $%
	$(RM) $%

$(LES_ARC): $(ARCOBJ_f) 
	$(ECHO)
	$(ECHO) "$(LES_ARC) up-to-date"
	$(ECHO)


#
# Dependencies
# --------------------
#
rad_rndnmb.o: rad_rndnmb.f90
rad_solver.o: rad_solver.f90 defs.f90
rad_corkds.o: rad_corkds.f90 defs.f90
rad_cldwtr.o: rad_cldwtr.f90 defs.f90 rad_corkds.f90
rad_d4strm.o: rad_d4strm.f90 defs.f90 rad_solver.f90 rad_rndnmb.f90
	rad_corkds.f90 rad_cldwtr.f90
rad_driver.o: rad_driver.f90 defs.f90 rad_d4strm.f90