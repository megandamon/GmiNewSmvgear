include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiSpecConcentration.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiSpecConcentration.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
GmiTaggedCO_AgeOfAir_mod.o: GmiSpcConcentrationMethod_mod.o \
                            GmiTracerMethod_mod.o
GmiProdLossDiagnostics_mod.o: GmiSpcConcentrationMethod_mod.o
