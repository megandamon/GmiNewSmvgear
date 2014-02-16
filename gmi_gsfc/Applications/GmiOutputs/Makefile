include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiOutputFiles.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiOutputFiles.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
GmiControlASCII_mod.o: m_ncGeneralOpsOutput.o
GmiControlConcentration_mod.o: m_ncGeneralOpsOutput.o
GmiControlTendencies_mod.o: m_ncGeneralOpsOutput.o
GmiControlColumn_mod.o: m_ncGeneralOpsOutput.o
GmiControlRestart_mod.o: m_ncGeneralOpsOutput.o
GmiControlFlux_mod.o: m_ncGeneralOpsOutput.o
GmiControlAerosolDust_mod.o: m_ncGeneralOpsOutput.o
GmiControlSAD.o: m_ncGeneralOpsOutput.o
GmiControlQqjk.o: m_ncGeneralOpsOutput.o
GmiControlQj.o: m_ncGeneralOpsOutput.o
GmiControlQk.o: m_ncGeneralOpsOutput.o
GmiControlFreq1.o: m_ncGeneralOpsOutput.o
GmiControlFreq2.o: m_ncGeneralOpsOutput.o
GmiControlFreq3.o: m_ncGeneralOpsOutput.o
GmiControlFreq4.o: m_ncGeneralOpsOutput.o
GmiControlOverpass1.o: m_ncGeneralOpsOutput.o
GmiControlOverpass2.o: m_ncGeneralOpsOutput.o
GmiControlOverpass3.o: m_ncGeneralOpsOutput.o
GmiControlOverpass4.o: m_ncGeneralOpsOutput.o
GmiControlOutput_mod.o: GmiControlASCII_mod.o GmiControlConcentration_mod.o GmiControlTendencies_mod.o GmiControlColumn_mod.o GmiControlRestart_mod.o GmiControlFlux_mod.o GmiControlAerosolDust_mod.o GmiControlSAD_mod.o GmiControlQqjk_mod.o GmiControlQj_mod.o GmiControlQk_mod.o GmiControlOverpass1_mod.o GmiControlOverpass2_mod.o GmiControlOverpass3_mod.o GmiControlOverpass4_mod.o GmiControlFreq1_mod.o GmiControlFreq2_mod.o GmiControlFreq3_mod.o GmiControlFreq4_mod.o
