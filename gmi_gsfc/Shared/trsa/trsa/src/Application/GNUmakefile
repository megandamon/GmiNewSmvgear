#
# Makefile for FVAdvStandalone Tracer Advection Test Application
#
# REVISION HISTORY:
#
# 20060922   Sawyer      Taken from fvsa/src/Application/GNUmakefile
# 20061203   Sawyer      Naming scheme revisions
# 20070608   Sawyer      Upgraded for MAPL implementation
# 20071224   Sawyer      Moved FV_Shared files to separate directory
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

THIS = FV_AdvectionEx
BIN  = $(THIS).x

esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"

THIS = FV_AdvectionEx
BIN  = $(THIS).x
TEXS = $(THIS)

PROTEX_FLAGS = -bM


esma_install install: $(BIN)
	$(MKDIR) $(ESMABIN) $(ESMAETC) 
	$(CP) -p *.x    $(ESMABIN)
	$(CP) -p *.rc   $(ESMAETC)

esma_clean esma_distclean clean distclean:
	-$(RM) *~ *.[aox] *.mod *.x

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"


#                  --------------------
#                  User Defined Targets
#                  --------------------

FREAL =  # Turn off this nasty -r4/-r8 stuff
SRCS = $(THIS).F90
OBJS = $(SRCS:.F90=.o)

USER_FINCS = -I$(ESMAINC)/MAPL
USER_FFLAGS = $(BIG_ENDIAN) $(M). \
              $(M)$(ESMAINC)/MAPL \
              $(M)$(ESMAINC)/GMAO_pilgrim \
              $(M)$(ESMAINC)/FV_Shared \
              $(M)$(ESMAINC)/FVadvcore_GridComp \

LIB_COMP = $(ESMALIB)/libFVadvcore_GridComp.a \
           $(ESMALIB)/libFV_Shared.a \
           $(ESMALIB)/libMAPL.a \
           $(ESMALIB)/libMAPL_cfio_r4.a \
           $(ESMALIB)/libMAPL_cfio_r8.a \
           $(ESMALIB)/libGMAO_pilgrim.a

$(BIN) bin : $(OBJS)
	echo "should include stubs $(LIB_ESMF)"
	$(LD) -o $(BIN) $(LDFLAGS) $(OBJS) $(LIB_COMP) $(LIB_ESMF) \
                                   $(LIB_NETCDF) $(LIB_MPI) $(LIB_SCI) $(LIB_SYS)

#                  --------------------
#                      Dependencies
#                  --------------------

$(THIS).o: $(LIB_COMP)

#.
