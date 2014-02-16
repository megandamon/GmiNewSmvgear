include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiChemMethod.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiChemMethod.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
GmiSolverInterface_mod.o: GmiSurfaceEmissionInChemistry_mod.o
GmiUpdateChemistry_mod.o: GmiShipEmission_mod.o GmiPhotolysisRateConstants_mod.o \
                          GmiThermalRateConstants_mod.o GmiCloudPropertiesGT_mod.o \
                          GmiUpdateSAD_mod.o GmiUpdateForcingBC_mod.o \
                          GmiStratosphericLoss_mod.o GmiSolverInterface_mod.o
GmiChemistryMethod_mod.o: GmiUpdateChemistry_mod.o
