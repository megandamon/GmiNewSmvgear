include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiGOCART.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiGOCART.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
GmiSeaSaltMethod_mod.o: GocartDerivedVariables_mod.o
GmiDustMethod_mod.o: GocartDerivedVariables_mod.o GmiWaterMethod_mod.o
GmiReadGocartSourceFiles_mod.o: GmiSeaSaltMethod_mod.o GmiDustMethod_mod.o
