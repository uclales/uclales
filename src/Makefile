# Generated automatically from Makefile.in by configure.
.SUFFIXES : .f .o .a .f90 .F90
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
$(LES_ARC)(util.o)   \
$(LES_ARC)(grid.o)   \
$(LES_ARC)(thrm.o)   \
$(LES_ARC)(rad_gcss.o)   \
$(LES_ARC)(modnetcdf.o)   \
$(LES_ARC)(modcross.o)   \
$(LES_ARC)(modnudge.o)   \
$(LES_ARC)(modtimedep.o)   \
$(LES_ARC)(ice_sb.o)   \
$(LES_ARC)(icemcrp.o)   \
$(LES_ARC)(lsmdata.o)   \
$(LES_ARC)(srfc.o)   \
$(LES_ARC)(ncio.o)   \
$(LES_ARC)(stat.o)   \
$(LES_ARC)(sgsm.o)   \
$(LES_ARC)(prss.o)   \
$(LES_ARC)(advf.o)   \
$(LES_ARC)(advl.o)   \
$(LES_ARC)(forc.o)   \
$(LES_ARC)(lsvar.o)   \
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


.F90.a:
	$(F90) -c $(FFLAGS) $(%:.o=.F90)
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

grid.o:		
advf.o:		grid.o stat.o util.o
rad_rndnmb.o:	defs.o
rad_solver.o:	defs.o
rad_corkds.o:	defs.o mpi_interface.o
rad_cldwtr.o:	defs.o rad_corkds.o
rad_d4strm.o:	defs.o rad_solver.o rad_rndnmb.o rad_corkds.o rad_cldwtr.o
rad_driver.o:	defs.o rad_d4strm.o
forc.o: 	defs.o grid.o rad_gcss.o stat.o
icemcrp.o: 	defs.o grid.o ice_sb.o
util.o: 	defs.o grid.o ice_sb.o
step.o:		icemcrp.o mpi_interface.o
init.o:		grid.o ncio.o step.o icemcrp.o
ncio.o:		grid.o mpi_interface.o icemcrp.o
lsvar.o:	forc.o
srfc.o: 	lsmdata.o
