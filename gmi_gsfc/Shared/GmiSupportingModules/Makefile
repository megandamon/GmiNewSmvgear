include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiSupportingModules.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiSupportingModules.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
#GmiSyntheticSpecies_mod.o: GmiChemCaseName_mod.o
#GmiTropWater_mod.o: GmiChemCaseName_mod.o
GmiFixedConstituent_mod.o: GmiTimeControl_mod.o
GmiSpeciesRegistry_mod.o: GmiStringManipulation_mod.o
