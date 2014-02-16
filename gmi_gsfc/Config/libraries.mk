#include base.mk
#
# use $(ARCHi) to determine the path to include files and libraries.
#
# ----------------------------------
# netCDF include files and libraries
# ----------------------------------

INCLUDES_NETCDF=
NETCDFLibs     =

ifeq ($(F90_VENDOR), IBM)
  INCLUDES_NETCDF=
  NETCDFLibs     =
endif

ifeq ($(F90_VENDOR), Compaq)
  INCLUDES_NETCDF=-I/usr/ulocal/include
  NETCDFLibs     =-L/usr/ulocal/lib -lnetcdf
endif

ifeq ($(F90_VENDOR), MIPSpro)
  INCLUDES_NETCDF=
  NETCDFLibs     =
endif

ifeq ($(F90_VENDOR), Intel)
      INCLUDES_NETCDF = -I$(BASEDIR)/Linux/include/netcdf
      NETCDFLibs      = -L$(BASEDIR)/Linux/lib -lnetcdf -lhdf5_hl -lhdf5 -lmfhdf -ldf -lcurl -lrt -lz -lsz -ljpeg -lm
endif

# -------------------------------------------
# Message Passing include files and libraries
# -------------------------------------------

INCLUDES_MsgPass=
MsgPassLibs     =

ifeq ($(F90_VENDOR), MIPSpro)
  INCLUDES_MsgPass =
  MsgPassLibs      =
endif

ifeq ($(F90_VENDOR), Intel)
      INCLUDES_MsgPass = 
      MsgPassLibs      = 
endif

# -------------------------------------------
# ESMF include files and libraries
# -------------------------------------------

INCLUDES_Esmf=
EsmfLibs     =

ifeq ($(F90_VENDOR), MIPSpro)
  INCLUDES_Esmf =
  EsmfLibs      =
endif

ifeq ($(F90_VENDOR), Intel)
      INCLUDES_Esmf = -I$(BASEDIR)/Linux/include/esmf
      EsmfLibs      = -L$(BASEDIR)/Linux/lib -lesmf -limf -lm -lrt -ldl -threads
      #EsmfLibs      = -L$(BASEDIR)/Linux/lib -lesmf -lcprts -limf -lm -lcxa -lunwind -lrt -ldl -threads
endif


advecCoreINC = $(INCLUDES_Esmf)

advecINCdir = $(GMIHOME)/Shared/trsa/Linux/include
advecLIBSdir = $(GMIHOME)/Shared/trsa/src/../Linux/lib

advecCoreLibs = 

GmiInternalLibs = \
    -lGmiAppUtils -lGmiOutputFiles \
    -ladvecCore \
    -lGmiAdvecMethod    -lGmiAdvecDao2      -lGmiUtilsDao2  \
    -lGmiConvecMethod   -lGmiDeposMethod  -lGmiDiffuMethod          \
    -lGmiChemMethod     -lGmiChemSmv2  -lGmiSulfur  -lGmiAmmonia  -lSolverShared \
    -lGmiChemSetkin     -lGmiChemSad  \
    -lGmiPhotFastj -lGmiPhotFast_JX -lGmiPhotFast_JX65 \
    -lGmiPhotFast_JX53b -lGmiPhotFastJX53c_ref -lGmiPhotLookup \
    -lGmiPhotUtils  -lioChemistry -lGmiAerosol_Dust   \
    -lGmiEmissMethod   -ldiagnEmission -lGmiEmissHarvard -lEmissionMEGAN \
    -lioEmission -lGmiEmissLlnl -lGmiEmissLightning  -lGmiGOCART  -lGmiGCR \
    -lGmiSpecConcentration -lioSpcConcentration  -lGmiDiagnostics \
    -lNcUtilsSingle   -lGmiMetFields -lGmiDomainDecomp -lGmiCommunications \
    -lGmiSupportingModules -lGmiIOutilities -lGmiESMF


export INCLUDES_NETCDF
export NETCDFLibs
export INCLUDES_MsgPass
export MsgPassLibs
export GmiInternalLibs
export INCLUDES_Esmf
export EsmfLibs
export advecCoreINC
