TOP_DIR=$(shell pwd)
ifndef $(BASEDIRe)
  BASEDIRe=$(TOP_DIR)
endif

#BASEDIR=/usr/local/other/baselibs/ESMF220rp2_NetCDF362b6_9.1.052
#BASEDIR=/usr/local/other/baselibs/ESMF310r_NetCDF362b6_10.1.017

CONFIG=$(TOP_DIR)/Config
COMPONENTS=$(TOP_DIR)/Components
APPLICATIONS=$(TOP_DIR)/Applications
CHEMISTRY=$(COMPONENTS)/GmiChemistry

SetkinInc=$(CHEMISTRY)/mechanisms/$(CHEMCASE)/include_setkin
SETKINDIR=$(CHEMISTRY)/mechanisms/$(CHEMCASE)

SHARED=$(TOP_DIR)/Shared

INSTALL_DIR=$(BASEDIRe)/../$(ARCHi)
BIN_DIR    =$(INSTALL_DIR)/bin
LIB_DIR    =$(INSTALL_DIR)/lib
INCLUDE_DIR=$(INSTALL_DIR)/include
MODULE_DIR =$(INSTALL_DIR)/mod

export CONFIG
export INSTALL_DIR
export BIN_DIR
export LIB_DIR
export INCLUDE_DIR
export MODULE_DIR
export APPLICATIONS

export SetkinInc
export SETKINDIR
export BASEDIR
