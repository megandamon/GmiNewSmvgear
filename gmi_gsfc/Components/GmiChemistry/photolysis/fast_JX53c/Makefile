include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiPhotFastJX53c.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiPhotFastJX53c.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
FastJX53cMetVars_mod.o: FastJX53cCTMparameters_mod.o
FastJX53cJvaluesVars_mod.o: FastJX53cCTMparameters_mod.o
FastJX53cCoreScattering_mod.o: FastJX53cMIEparameters_mod.o 
FastJX53cCoreFastj_mod.o: FastJX53cJvaluesVars_mod.o FastJX53cCTMparameters_mod.o FastJX53cCoreScattering_mod.o
FastJX53cCTM_FastjLinking_mod.o: FastJX53cMetVars_mod.o FastJX53cJvaluesVars_mod.o FastJX53cCTMparameters_mod.o FastJX53cMIEparameters_mod.o FastJX53cCoreFastj_mod.o
FastJX53cCTMspecific_mod.o: FastJX53cMetVars_mod.o FastJX53cCTMparameters_mod.o FastJX53cCoreFastj_mod.o FastJX53cJvaluesVars_mod.o FastJX53cCTMparameters_mod.o
FastJX53cMethod_mod.o: FastJX53cMetVars_mod.o FastJX53cCTM_FastjLinking_mod.o FastJX53cCTMspecific_mod.o FastJX53cCTMparameters_mod.o FastJX53cJvaluesVars_mod.o
