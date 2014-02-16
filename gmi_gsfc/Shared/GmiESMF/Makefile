include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	cp *.h   $(INCLUDE_DIR)
	ar crs libGmiESMF.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiESMF.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
GmiESMFclock_mod.o: GmiESMF_LatLonGrid_mod.o GmiTimeControl_mod.o
GmiESMF_LatLonGrid_mod.o: GmiESMF_ErrorChecking_mod.o
