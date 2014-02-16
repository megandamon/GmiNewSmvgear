include $(CONFIG)/rules.mk
#include $(CONFIG)/libraries.mk

#%.o:%.F90
#	$(F90) -c $(FFLAGS) $(INCLUDES_NETCDF) $<

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all: $(OBJS)
	ar crs libNcUtilsSingle.a $(OBJS)
	cp *.mod $(MODULE_DIR)
	cp libNcUtilsSingle.a $(LIB_DIR)

clean: 
	$(RM) *.o *.mod *.a *i90 *f90


#dependencies
