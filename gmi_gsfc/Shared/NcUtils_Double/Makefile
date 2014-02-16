include $(CONFIG)/rules.mk

NETCDFIncs=-I/usr/ulocal/include
NETCDFLibs=-L/usr/ulocal/lib -lnetcdf

%.o:%.F90
	$(F90) -c $(FFLAGS) $(NETCDFIncs) $<

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

libNcUtilsDouble.a: $(OBJS)
	ar crs libNcUtilsDouble.a $(OBJS)
	cp *.mod $(MODULE_DIR)
#	cp include/*.h   $(INCLUDE_DIR)
	cp libNcUtilsDouble.a $(LIB_DIR)

#.PHONY: libtest_chemistry.a

clean: 
	$(RM) *.o *.mod *.a
#	$(MAKE) -C tests clean


#dependencies
