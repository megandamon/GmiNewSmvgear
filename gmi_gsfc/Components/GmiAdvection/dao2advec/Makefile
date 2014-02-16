include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiAdvecDao2.a $(OBJS)
#	cp *.mod $(MODULE_DIR)
	cp libGmiAdvecDao2.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
#xyz_ppm_dao2.o: GmiAdvectionFlux_mod.o
