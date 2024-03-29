#
# recursive makefile for ESMA.
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
	@echo "         SITE = $(SITE)"


#                  --------------------------------
#                   Recurse Make in Sub-directories
#                  --------------------------------

ALLDIRS = arpack nagpack \
          GMAO_mpeu GMAO_pilgrim GMAO_mfhdf3 \
          GMAO_gfio GMAO_gfioeos GMAO_perllib \
          GMAO_ods  GMAO_cfio MAPL_cfio GEOS_Base \
          MAPL_Base GEOS_Shared Chem_Base Chem_Shared\
          GMAO_etc  GMAO_hermes GMAO_transf GMAO_psas GMAO_gems \
          GMAO_iret GEOS_Util

SUBDIRS = $(wildcard $(ALLDIRS) )

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          help install doc clean distclean

export ESMADIR BASEDIR ARCH SITE

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( cd $$d				;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done


local_esma_install local_install: $(LIB)
	@echo No local install in here

#              ------------------------------------------
#              Package Dependencies for Parallel Install
#              ------------------------------------------

ifeq ($(wildcard GMAO_mfhdf3), $(null))
   GMAO_MFHDF3 =
else
   GMAO_MFHDF3 = GMAO_mfhdf3_install
endif

Chem_Base_install   : GEOS_Base_install   \
                      GEOS_Shared_install \
                      GMAO_mpeu_install   \
                      GMAO_gfio_install   \
                      $(GMAO_MFHDF3)
Chem_Shared_install : Chem_Base_install   \
                      GEOS_Base_install   \
                      GEOS_Shared_install \
                      GMAO_gfio_install   \
                      GMAO_mpeu_install   \
                      GMAO_pilgrim_install
GEOS_Base_install   : GMAO_cfio_install
MAPL_Base_install   : GMAO_cfio_install
GEOS_Shared_install : GEOS_Base_install   \
                      GMAO_cfio_install
GMAO_cfio_install   : $(GMAO_MFHDF3)
MAPL_cfio_install   : $(GMAO_MFHDF3) 
GMAO_etc_install    : GMAO_mpeu_install
GMAO_gfio_install   : $(GMAO_MFHDF3)
GMAO_hermes_install : GMAO_mpeu_install   \
                      GMAO_cfio_install   \
                      GMAO_gfio_install   \
                      GEOS_Base_install   \
                      GEOS_Shared_install
GMAO_iret_install   : GMAO_hermes_install \
                      GMAO_mpeu_install   \
                      Chem_Base_install
GMAO_ods_install    : GMAO_mpeu_install   \
                      GMAO_pilgrim_install
GMAO_psas_install   : GMAO_mpeu_install
GMAO_transf_install : GMAO_hermes_install \
                      GMAO_mpeu_install   \
                      GMAO_hermes_install

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
