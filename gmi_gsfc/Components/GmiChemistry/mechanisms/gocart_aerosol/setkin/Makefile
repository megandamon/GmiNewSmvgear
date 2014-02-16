include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiChemSetkin.a $(OBJS)
#	cp *.mod $(MODULE_DIR)
	cp libGmiChemSetkin.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
