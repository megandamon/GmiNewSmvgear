# use $(ARCHi) to determine default compiler
# then use $(COMPILER) and $(OPT) to determine the default flags

# -----------------------------------------------------
# This option is used to set some non used species indices
# to a non zero value to allow the code to be compiled with
# array bound checking. This does not affect at all the
# calculations.
# -----------------------------------------------------
deadInd=
deadInd_tracers=

ifeq ($(CHEMCASE), strat_trop)
deadInd=nonZeroInd
endif

ifeq ($(CHEMCASE), troposphere)
deadInd=nonZeroInd
endif

ifeq ($(CHEMCASE), stratosphere)
deadInd=nonZeroInd
endif

ifeq ($(CHEMCASE), tracers)
deadInd_tracers=nonZeroInd_tracers
endif
# -----------------------------------------------------
# This option is to select the Georgia Tech cloud module
# -----------------------------------------------------
GTcloud=

#---------------------------------------------------------
# Comment out the line below if the module is not used.
# The module can only be used if the chemical mechanism is
# aerosol, gocart_aerosol or micro_aerosol.
#---------------------------------------------------------
#GTcloud=GTmodule

# -----------------------------------------------------

ifndef $(OPT)
  OPT=DEBUG
endif

ifeq ($(ARCHi), Darwin)
  VENDOR=Apple
  ifndef $(COMPILER)
    F90_VENDOR=IBM
  endif
endif

ifeq ($(ARCHi), IRIX64)
  VENDOR=SGI
  F90_VENDOR=MIPSpro
endif

ifeq ($(ARCHi), OSF1)
  VENDOR=HP
  F90_VENDOR=Compaq
endif

ifeq ($(ARCHi), Linux)
  VENDOR=SGI
  F90_VENDOR=Intel
endif

ifeq ($(CHEMCASE), micro_aerosol)
  SULFCASE=micro_sulfur
  EXTRAFLAGS=1
else
  SULFCASE=sulfur
  EXTRAFLAGS=0
  AerosolCase=
ifeq ($(CHEMCASE), gocart_aerosol)
      AerosolCase=GOCARTaerosol
endif
endif


ifeq ($(F90_VENDOR), IBM)
ifeq ($(USE_MPI),YES)
  F90        = mpif77
else
  F90        = xlf90
endif

  MODULE_INC = 
#  MODULE_INC = -I
  F90DEF     = -WF,-D
  FFLAGS     = -qsuffix=cpp=F90 -qfree=f90 -qzerosize -qnosave
  ifeq ($(OPT),DEBUG)
    FFLAGS += -g -C
  else
    FFLAGS += -O2
  endif
endif

ifeq ($(F90_VENDOR), Compaq)
  F90        = f90
  F90DEF     = -D
  FFLAGS     = -arch host -fast -fpe -assume accuracy_sensitive 
  CPP        = $(F90) -cpp
  ifeq ($(EXTRAFLAGS), 1)
    CPPFLAGS   = -DMICRO_AEROSOL -DKULNEW -DPARNEW -D${HOSTMACH} -D${GTcloud} -D${deadInd_tracers} -D${deadInd} -P -Wp,-P
  else
    CPPFLAGS   = -D${HOSTMACH} -D${AerosolCase} -D${GTcloud} -D${deadInd} -D${deadInd_tracers} -P -Wp,-P
  endif
  CC         = cc
endif

ifeq ($(F90_VENDOR), MIPSpro)
  F90        = f90
  F90DEF     = -D
  FFLAGS     = -c -n32 -mips4 -O3 -OPT:Olimit=0:IEEE_arithmetic=3 -LNO:prefetch=2
  CPP        = /usr/lib/cpp 
  ifeq ($(EXTRAFLAGS), 1)
    CPPFLAGS   = -DMICRO_AEROSOL -DKULNEW -DPARNEW -D${GTcloud} -D${deadInd} -D${deadInd_tracers} -P
  else
    CPPFLAGS   = -D${AerosolCase} -D${GTcloud} -D${deadInd} -D${deadInd_tracers} -P
  endif
endif

ifeq ($(F90_VENDOR), Intel)
  F90DEF     = -D

  F90        = mpif90 
  DebugFLAGS = -g -O0 -check bounds -fpe-all=0 -fpe0
  GenFLAGS   = -warn nogeneral -traceback
  OptFLAGS   = -O2
  gmaoOPT = -O3 -vec-report0 -ftz -align all -fno-alias -fPIC -fpe0 -fp-model precise
#  OptFLAGS   = -O4 -finline-limit=1000 -fno-alias -opt-mem-bandwidth2 -no-prec-div -xP
#  FFLAGS     = $(DebugFLAGS) $(GenFLAGS)
   FFLAGS     = $(OptFLAGS) $(GenFLAGS)
#  FFLAGS     = $(gmaoOPT) -warn nogeneral
  CPP        = $(F90) -cpp -132
  CC         = icc
  ifeq ($(EXTRAFLAGS), 1)
    CPPFLAGS   = -DMICRO_AEROSOL -DKULNEW -DPARNEW -D${GTcloud} -D${deadInd} -D${deadInd_tracers} -E
  else
    CPPFLAGS   = -D${AerosolCase} -D${GTcloud} -D${deadInd} -D${deadInd_tracers} -E
  endif

endif


export F90_VENDOR
export F90
export MODULE_INC
export F90DEF
export FFLAGS
export CPP
export CPPFLAGS
export CC
export EXTRAFLAGS
export SULFCASE
export AerosolCase
