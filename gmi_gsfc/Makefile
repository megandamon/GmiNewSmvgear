include Config/base.mk
include Config/architecture.mk
include Config/compiler.mk
include Config/libraries.mk

VPATH=$(LIB_DIR)
export VPATH

LIBS = -L$(LIB_DIR)

# ---------------
# Include Modules
# ---------------
IncMods  =-I$(MODULE_DIR)
IncMods +=$(INCLUDES_Esmf)
IncMods +=$(advecCoreINC)
FFLAGS   +=$(IncMods)

# -------------
# Include Files
# -------------
IncFiles  =-I$(INCLUDE_DIR)
IncFiles +=$(INCLUDES_NETCDF)
IncFiles +=$(INCLUDES_MsgPass)
IncFiles +=$(MAPL_INCdir)
CPPFLAGS +=$(IncFiles)

export LIBS

echo_config:
	@echo "UNAME:         " $(UNAME)
	@echo "ARCH:          " $(ARCH)
	@echo "VENDOR:        " $(VENDOR)
	@echo "BASEDIR:       " $(BASEDIR)
	@echo "F90_VENDOR:    " $(F90_VENDOR)
	@echo "OPT            " $(OPT)
	@echo "FFLAGS:        " $(FFLAGS)
	@echo "HOSTMACH:      " $(HOSTMACH)

new: packageddir shared components applications
all: packageddir shared components applications

packageddir:
	if [ ! -e $(INSTALL_DIR) ] ; then \
	  mkdir -p $(INSTALL_DIR) ; \
	  mkdir -p $(INSTALL_DIR)/mod ; \
	  mkdir -p $(INSTALL_DIR)/lib ; \
	  mkdir -p $(INSTALL_DIR)/include ; \
	  mkdir -p $(INSTALL_DIR)/bin ; \
	fi

components:
	@$(MAKE) -C $(COMPONENTS) all

shared:
	@$(MAKE) -C $(SHARED) all

applications:
	@$(MAKE) -C $(APPLICATIONS) all

distclean: clean
	$(RM) $(BIN_DIR)/*.x
	$(RM) $(LIB_DIR)/*.a
	$(RM) $(MODULE_DIR)/*.mod
	$(RM) $(INCLUDE_DIR)/*.h

clean:
	$(MAKE) -C $(APPLICATIONS) clean
	$(MAKE) -C $(SHARED) clean
	$(MAKE) -C $(COMPONENTS) clean

