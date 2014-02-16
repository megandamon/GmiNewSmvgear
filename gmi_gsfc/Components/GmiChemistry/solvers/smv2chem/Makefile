include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiChemSmv2.a $(OBJS)
#	cp *.mod $(MODULE_DIR)
	cp libGmiChemSmv2.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
