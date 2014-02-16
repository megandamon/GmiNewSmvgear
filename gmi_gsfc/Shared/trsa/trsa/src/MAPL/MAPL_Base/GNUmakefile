#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 20oct2004  da Silva  First ESMA version.
# 23Aug2006  da Silva  Driver for making either MAPL or GEOS_Base.
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(shell pwd)/../../..
endif

THIS := $(shell basename `pwd`)

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help :
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
	@echo "         SITE = $(SITE) "
	@echo "         THIS = $(THIS) "


#                  --------------------------------
#                   Recurse Make in Sub-directories
#                  --------------------------------

ifeq ($(THIS),MAPL_Base)
     SUBDIRS = MAPL
     THIS_CFIO = MAPL_cfio_r4
else
     SUBDIRS = MAPL #GEOS_Base
     THIS_CFIO = GMAO_cfio_r4
endif

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          help install test clean

export THIS_CFIO

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( $(MKDIR) $$d            		;\
	      $(MAKE) -C $$d -f ../MAPL_base.mk -e $$t ) \
	  done

local_esma_install local_install: $(LIB)
	@echo No local install in here

#                  --------------------
#                     Documentation
#                  --------------------

GGUG_SRCS = MAPL_Generic.F90 MAPL_GenericCplComp.F90  MAPL_sun_uc.P90 MAPL_SatVapor.F90

esma_doc doc:
	$(PROTEX) -g  MAPL_GenericIntro.tex  > $(ESMADOC)/MAPL_Intro.tex
	$(PROTEX) -g  $(GGUG_SRCS)           > $(ESMADOC)/MAPL_UsersGuide.tex

#                  --------------------
#                  User Defined Targets
#                  --------------------

MAPL_install:
	$(MKDIR) MAPL
	$(MAKE) -C MAPL -f ../MAPL_base.mk -e install 

GEOS_Base_install:
	$(MKDIR) GEOS_Base
	$(MAKE) -C GEOS_Base -f ../MAPL_base.mk -e install 

distclean:
	$(MAKE) esma_distclean
	-$(RM) -r *~ *.[aoxd] *.mod $(SUBDIRS)


  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
