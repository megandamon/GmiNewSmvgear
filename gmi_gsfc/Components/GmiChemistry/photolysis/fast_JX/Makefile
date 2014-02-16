include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiPhotFast_JX.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiPhotFast_JX.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
fastJX_Bundle_mod.o: fastJX_Parameters_mod.o
fast_JX_control.o: fastJX_Bundle_mod.o fastJX_Parameters_mod.o
