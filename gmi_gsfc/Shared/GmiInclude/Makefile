include $(CONFIG)/rules.mk

%.o:%.F90
	$(F90) -c $(FFLAGS) $<

SRCS_F = $(wildcard *.F90)
OBJS   = $(SRCS_F:.F90=.o)

all: 
	cp *.h   $(INCLUDE_DIR)
	cp $(GMIHOME)/Shared/trsa/src/MAPL/MAPL_Base/*.h $(INCLUDE_DIR)
	cp $(SetkinInc)/*.h $(INCLUDE_DIR)
	cp $(GMIHOME)/Config/*.h $(INCLUDE_DIR)

clean: 
	$(RM) *.o *.mod *.a


#dependencies
