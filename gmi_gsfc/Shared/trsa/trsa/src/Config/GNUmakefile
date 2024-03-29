#
# Makefile for ESMA build configuration system
#

# Generic Base/Arch + User Defined Fragments
# ------------------------------------------
  ESMA_BASE = ESMA_base.mk $(wildcard $(HOME)/.esma_base.mk) 
  ESMA_ARCH = ESMA_arch.mk $(wildcard $(HOME)/.esma_arch.mk) 
  ESMA_POST = ESMA_post.mk $(wildcard $(HOME)/.esma_post.mk) 

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMA_BASE)            # Generic stuff
  include $(ESMA_ARCH)            # System dependencies

# Make sure ESMADIR/BASEDIR are defined
# -------------------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../..
endif

  BINS   = # f90aib.x f90split.x

# User can provide their own assert script
# ----------------------------------------
  ASSERT_MK = 
  ifdef ASSERT
        RC = $(shell $(ASSERT))
        ifeq ($(RC), 0)
              $(warning Enabling build environment check)
	      ASSERT_MK = Assert.mk
        else 
           ifeq ($(RC), -1)
             $(warning - Node not know, disabling build environment check)
           else 
             $(error Please correct your build environment and try again)
           endif
        endif
  endif


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

esma_install install : fdp.o $(BINS)
	$(MKDIR) $(ESMABIN) $(ESMALIB) $(ESMAINC) \
		 $(ESMAETC) $(ESMADOC) $(ESMACFG)

	@if [ "$(ESMACFG)" != "$(ESMADIR)/Config" ]; then (\
	   $(RM) -rf $(ESMADIR)/Config;  \
	   $(LN) -s $(ESMACFG) $(ESMADIR)/Config  \
	) fi

	$(CAT) $(ESMA_BASE) > $(ESMACFG)/ESMA_base.mk
	$(CAT) $(ESMA_ARCH) > $(ESMACFG)/ESMA_arch.mk
	$(CAT) $(ESMA_POST) > $(ESMACFG)/ESMA_post.mk
	$(CAT) GMAO_base.mk $(ASSERT_MK) > $(ESMACFG)/GMAO_base.mk

	$(CP) -p protex $(MODINFO) $(BINS) $(ESMABIN)

ifdef ASSERT
	$(CP) $(ASSERT) $(ESMABIN)/Assert
endif
	./install_fdp ./fdp $(FDP)
	./CVSTAG_install.pl $(ESMAETC)

esma_clean clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"

FOPT =

%.x : %.F90
	$(FC) -o $@ $(F90FLAGS) $<

  include $(ESMA_POST)

#.

