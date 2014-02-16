include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiAdvecMethod.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiAdvecMethod.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
GmiAdvectionMethod_mod.o: advec_update.o
