include $(CONFIG)/rules.mk

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)


all: $(OBJS)
	ar crs libGmiPhotFast_JX53b.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libGmiPhotFast_JX53b.a $(LIB_DIR)

clean:
	$(RM) *.o *.mod *.a *.f90

#dependencies
