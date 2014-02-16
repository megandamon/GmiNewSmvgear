include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiConvecMethod.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiConvecMethod.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
convec_update.o: convec_ncar_old.o convec_ncar.o convec_dao2.o convec_gmao.o
