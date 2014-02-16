#
# Makefile for FV common files
#
#
# REVISION HISTORY:
#
# 20071009  Sawyer    based on standard ESMA makefile
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
#  include               tracer_advection_arch.mk  # arch dependent flags 

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"

THIS = FV_Shared
LIB  = lib$(THIS).a

esma_install install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS)
	$(CP) -p *.a    $(ESMALIB)
	$(CP) -p *.mod  $(ESMAINC)/$(THIS)
#	$(CP) -p *.rc   $(ESMAETC)

esma_clean esma_distclean clean distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] *.x

esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS = shr_kind_mod.F90 fill_module.F90 tp_core.F90 \
       mapz_module.F90 TracerMod.F90 TransposeMod.F90

OBJS = $(SRCS:.F90=.o) 
DEPS = $(SRCS:.F90=.d) 

INC_DIRS = $(INC_MPI)
MOD_DIRS = . \
           $(ESMAINC)/GMAO_pilgrim \
           $(ESMAINC)/MAPL

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
vpath % $(INC_DIRS) $(MOD_DIRS)

USER_FFLAGS += $(BIG_ENDIAN) 
USER_FDEFS  = 

PPAWK=./pp.awk
FREAL = # avoid the nasty -r4 and -r8 flags at all costs

$(LIB) lib : $(DEPS) $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

#                  --------------------
#                      Dependencies
#                  --------------------

# Make sure dep files are not remade during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.
