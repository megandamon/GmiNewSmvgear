include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiIOutilities.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiIOutilities.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
GmiCheckRange_mod.o: GmiPrintError_mod.o
GmiReadList_mod.o: GmiPrintError_mod.o
GmiFileUnit_mod.o: GmiPrintError_mod.o
GmiASCIIoperations_mod.o: GmiPrintError_mod.o GmiFileUnit_mod.o
