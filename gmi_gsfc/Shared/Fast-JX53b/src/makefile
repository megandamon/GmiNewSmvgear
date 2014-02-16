# Test problem for sparinv
#

OBJS = m_FastjMIEparameters.o    m_FastjCTMparameters.o \
       m_FastjJvaluesVars.o      m_FastjMetVars.o       \
       m_FastjCoreScattering.o   m_FastjCoreFastj.o     \
       m_FastjCTM_FastjLinking.o m_FastjCTMspecific.o   \
       m_FastjReadNamelistFile.o m_FastjMethod.o

FC = f90

#OPT  =  -static
OPT  =  

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) $(OPT) -c $*.F90

FastjStandalone.ex: $(OBJS) FastjStandalone.o 
	$(FC) -o FastjStandalone.ex  FastjStandalone.o $(OBJS)
	cp FastjStandalone.ex ../bin

clean: 
	'rm' *.ex *.o *.mod
	'rm' ../bin/*.ex
