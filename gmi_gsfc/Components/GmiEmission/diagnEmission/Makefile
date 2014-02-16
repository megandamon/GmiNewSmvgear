include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libdiagnEmission.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libdiagnEmission.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
#GmiEmissionMethod_mod.o: GmiEmissionReadData_mod.o emiss_update.o
