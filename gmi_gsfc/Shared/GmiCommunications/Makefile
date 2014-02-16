include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiCommunications.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiCommunications.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
GmiGhostZones_mod.o: GmiMessagePassing_mod.o 
GmiBroadcast_mod.o: GmiMessagePassing_mod.o
GmiSubdomainsBC_mod.o: GmiGhostZones_mod.o
GmiSub2Glob_mod.o: GmiSubDomainsBC_mod.o GmiMessagePassing_mod.o
GmiGlob2Sub_mod.o: GmiMessagePassing_mod.o
GmiDomainDecomposition_mod.o: GmiGrid_mod.o GmiMessagePassing_mod.o 
