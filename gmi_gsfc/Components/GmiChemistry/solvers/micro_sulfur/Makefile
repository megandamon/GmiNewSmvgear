include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiSulfur.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiSulfur.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
modbox_umaer.o: precision.o
umaerosol.o: precision.o modpar_umaer.o modbox_umaer.o
sulfchem.o: modinit_umaer.o
modpar_umaer.o: precision.o
modinit_umaer.o: precision.o modpar_umaer.o modbox_umaer.o
modiag_umaer.o: precision.o
