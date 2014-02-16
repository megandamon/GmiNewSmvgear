include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libioChemistry.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libioChemistry.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
#GmiUpdateChemistry_mod.o: GmiAccumQqjk_mod.o               GmiUpdateQj_mod.o \
#                          GmiGravitationalSettling_mod.o   GmiUpdateQk_mod.o \
#                          GmiUpdateForcingBCs_mod.o        GmiUpdateLoss_mod.o        
