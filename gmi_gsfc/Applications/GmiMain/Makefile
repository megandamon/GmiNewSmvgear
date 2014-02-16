include $(CONFIG)/libraries.mk
include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

gmi.x: $(OBJS)
	$(F90) -o gmi.x $(OBJS) $(LIBS) $(GmiInternalLibs) $(EsmfLibs) $(NETCDFLibs) $(MsgPassLibs)
	cp *.mod $(MODULE_DIR)
	mv gmi.x $(BIN_DIR)

clean:
	rm -f *.o *.x *.mod *.f90 *.i90

##dependencies
GmiMain.o: GmiControlInitialize_mod.o GmiControlAdvance_mod.o GmiControlFinalize_mod.o
