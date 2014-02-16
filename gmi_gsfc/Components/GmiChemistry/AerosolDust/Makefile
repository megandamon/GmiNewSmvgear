include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiAerosol_Dust.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiAerosol_Dust.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
GmiAerDustColumnMass_mod.o: GmiReadAerosolDust_mod.o
