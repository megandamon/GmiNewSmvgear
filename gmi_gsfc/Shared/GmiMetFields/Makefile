include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiMetFields.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiMetFields.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
GmiPressureFixer_mod.o: GmiCalcDivergence_mod.o GmiMassFluxes_mod.o
GmiUtilsMetFields_mod.o: GmiPressureFixer_mod.o GmiCalcDivergence_mod.o 
GmiMetFieldsControl_mod.o: GmiUtilsMetFields_mod.o GmiModelData_mod.o GmiSurfaceAlbedo_mod.o  GmiMetFieldsAttribute_mod.o GmiTotalMass_mod.o
GmiSurfaceTemperature_mod.o: GmiMetFieldsControl_mod.o
