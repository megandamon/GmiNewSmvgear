
/****************************************************************
 *
 * CODE DEVELOPER
 *   John Tannahill, LLNL (Original code from Bill Bosl, LLNL)
 *   jrt@llnl.gov
 *
 * FILE
 *   gem_config.h
 *
 * DESCRIPTION
 *   This file sets configuration parameters for gem Makefiles.
 *
 * HISTORY
 *   - July 1, 2004 - Jules Kouatchou
 *      o New directories for netCDF library and include files
 *        on the COMPAQ architecture.
 *      o Removed the full path when using the cc and f90
 *        compilers on the COMPAQ architecture.
 *      o Point to the ifort compiler when the single processor
 *        version of the code is used on the INTEL architecture.
 ***************************************************************/

# include "gem_sys_options.h"
# include "gem_options.h"

  PROG_NAME       = gmi.x

/* -------------------------------------------------------
   This option is to select the Georgia Tech cloud module
   ------------------------------------------------------- */
GTcloud= 

/* Comment out the line below if the module is not used. */
/* GTcloud=GTmodule  */

/* The module can only be used if the chemical mechanism is
   aerosol, gocart_aerosol or micro_aerosol. */

#ifdef  strat_trop_aerosol
GTcloud=
#endif

#ifdef  strat_trop
GTcloud=
#endif

#ifdef  troposphere
GTcloud=
#endif

#ifdef  stratosphere
GTcloud=
#endif
/* ------------------------------------------------------- */

#ifdef  micro_aerosol
#define SULFCASE micro_sulfur
#define EXTRAFLAGS 1
#else
#define EXTRAFLAGS 0
#define SULFCASE sulfur
   AerosolCase = 
#ifdef  gocart_aerosol
   AerosolCase = GOCARTaerosol
#endif
#endif

#if (ARCH_OPTION == ARCH_IBM_SP)
  MKMF            = (/bin/rm -f Makefile; make -f $(GMIHOME)/Config/Makefile.init.sp  Makefile)
#elif (ARCH_OPTION == ARCH_INTEL)
  MKMF            = (/bin/rm -f Makefile; make -f $(GMIHOME)/Config/Makefile.init.int  Makefile)
#else
  MKMF            = (/bin/rm -f Makefile; make -f $(GMIHOME)/Config/Makefile.init     Makefile)
#endif
 
  MAKE_INC        = -I$(GMIHOME)/Config

  GMIDIR          = ${GMIHOME}

  GMIBINDIR       = $(GMIDIR)/Applications/GmiBin
  GMIINCLUDE      = $(GMIDIR)/Config

  GEMDIR          = ${GMIHOME}/gem

  ESMDIR          = $(GEMDIR)/esm
  ESMINCLUDE      = $(GEMDIR)/esm/include
  ESMLIBDIR       = $(GEMDIR)/esm/lib

  ESMTOOLSDIR     = $(GEMDIR)/esm_tools
  ESMTOOLSINCLUDE = $(GEMDIR)/esm_tools/include
  ESMTOOLSLIBDIR  = $(GEMDIR)/esm_tools/lib

  ACTMDIR         = $(GEMDIR)/actm
  ACTMLIBDIR      = $(GEMDIR)/actm/lib

  CHEMDIR         = $(GMIDIR)/Components/GmiChemistry/mechanisms/$(CHEMCASE)
  SETKINDIR       = $(CHEMDIR)/setkin
  SETKININCLUDE   = $(CHEMDIR)/include_setkin

#if (ARCH_OPTION == ARCH_COMPAQ)
    ARCHi = SC45
#elif (ARCH_OPTION == ARCH_SGI_ORIG)
    ARCHi = Origin
#elif (ARCH_OPTION == ARCH_INTEL)
    ARCHi = Intel
#endif

  GmiLibsDir       =   $(GMIDIR)/$(ARCHi)/lib
  GmiModsDir       = -I$(GMIDIR)/$(ARCHi)/mod
  GmiIncDir        = -I$(GMIDIR)/$(ARCHi)/include
  Gmi_CLEAN_TARGET =    Gmi_clean_target

/****************************************************
 * Define some target names for use in the Makefile.
 ****************************************************/

  GEM_CLEAN_TARGET            = gem_clean_target  
  GEM_COMPILE_ONLY_TARGET     = gem_compile_only_target
  GEM_CPP_TARGET              = gem_cpp_target
  GEM_LINK_TARGET             = gem_link_target

  ESM_CLEAN_TARGET            = esm_clean_target
  ESM_COMPILE_TARGET          = esm_compile_target
  ESM_COMPILE_ONLY_TARGET     = esm_compile_only_target
  ESM_CPP_TARGET              = esm_cpp_target

  ESMTOOLS_CLEAN_TARGET       = esmtools_clean_target
  ESMTOOLS_COMPILE_TARGET     = esmtools_compile_target

  ACTM_CLEAN_TARGET           = actm_clean_target
  ACTM_COMPILE_TARGET         = actm_compile_target
  ACTM_COMPILE_ONLY_TARGET    = actm_compile_only_target
  ACTM_CPP_TARGET             = actm_cpp_target


/**************
 * GEM macros.
 **************/

  GEM_CLEAN         = $(GEM_CLEAN_TARGET)  
  GEM_COMPILE       = $(GEM_COMPILE_TARGET)
  GEM_COMPILE_ONLY  = $(GEM_COMPILE_ONLY_TARGET)
  GEM_CPP           = $(GEM_CPP_TARGET)
  GEM_LINK_PREREQ   = $(GEM_LINK_TARGET)


/**************
 * ESM macros.
 **************/

  ESM_CLEAN         = $(ESM_CLEAN_TARGET)
  ESM_COMPILE       = $(ESM_COMPILE_TARGET)
  ESM_COMPILE_ONLY  = $(ESM_COMPILE_ONLY_TARGET)
  ESM_CPP           = $(ESM_CPP_TARGET)

  ESM1Libs          = -L$(ESMLIBDIR) -lESMcontrol
  ESM2Libs          = -L$(ESMLIBDIR) \
    -lESMcomm  -lESMin_out  -lESMmem_manage  -lESMutils


/*******************
 * ESMTOOLS macros.
 *******************/

  ESMTOOLS_CLEAN   = $(ESMTOOLS_CLEAN_TARGET)
  ESMTOOLS_COMPILE = $(ESMTOOLS_COMPILE_TARGET) 

#if ((ARCH_OPTION == ARCH_COMPAQ) || \
     (ARCH_OPTION == ARCH_INTEL)  || \
     (ARCH_OPTION == ARCH_SGI_ORIG))
  ESMTOOLSLibs = -L$(ESMTOOLSLIBDIR) -lESMTOOLSfi       -lESMTOOLSbaseline
#else
  ESMTOOLSLibs = -L$(ESMTOOLSLIBDIR) -lESMTOOLSbaseline -lESMTOOLSfi
#endif


/***************
 * ACTM macros.
 ***************/

  ACTM_Package_Dir  = gmimod

  ACTM_CLEAN        = $(ACTM_CLEAN_TARGET)
  ACTM_COMPILE      = $(ACTM_COMPILE_TARGET)
  ACTM_COMPILE_ONLY = $(ACTM_COMPILE_ONLY_TARGET)
  ACTM_CPP          = $(ACTM_CPP_TARGET)

  ACTMLibs = -L$(ACTMLIBDIR)  \
    -lACTMControl            \
    -lGmiAdvecMethod    -lGmiAdvecDao2      -lGmiUtilsDao2          \
    -lACTMIn_Out                        \
    -lACTMComm          -lGMIncOutputFiles

  GmiLibs = -L$(GmiLibsDir) \
    -ladvecCore -lGmiAdvecMethod    -lGmiAdvecDao2      -lGmiUtilsDao2    \
    -lGmiConvecMethod   -lGmiDeposMethod  -lGmiDiffuMethod          \
    -lGmiChemMethod     -lGmiChemSmv2  -lGmiSulfur  -lGmiAmmonia  -lGmiChemSetkin  \
    -lGmiChemSad        -lGmiPhotFastj -lGmiPhotFast_JX             \
    -lGmiPhotFast_JX53b -lGmiPhotFastJX53c_ref -lGmiPhotLookup \
    -lGmiPhotUtils  -lioChemistry            \
    -lGmiAerosol_Dust   \
    -lGmiEmissMethod   -ldiagnEmission -lGmiEmissHarvard -lEmissionMEGAN -lioEmission \
    -lGmiEmissLlnl     -lGmiEmissLightning  -lGmiGOCART  -lGmiGCR    \
    -lGmiSpecConcentration -lioSpcConcentration  -lGmiDiagnostics \
    -lGmiMetFields \
    -lNcUtilsSingle -lGmiDomainDecomp -lGmiCommunications \
    -lGmiSupportingModules -lGmiIOutilities -lGmiESMF

/*******************
 * Combined macros.
 *******************/

/*
  CLEAN_PREREQ        = $(GEM_CLEAN)          $(ESM_CLEAN) \
                        $(ESMTOOLS_CLEAN)     $(ACTM_CLEAN)

  COMPILE_PREREQ      = $(GEM_COMPILE)        $(ESM_COMPILE) \
                        $(ESMTOOLS_COMPILE)   $(ACTM_COMPILE)

  COMPILE_ONLY_PREREQ = $(GEM_COMPILE_ONLY)   $(ESM_COMPILE_ONLY) \
                        $(ACTM_COMPILE_ONLY)

  CPP_PREREQ          = $(GEM_CPP)   $(ESM_CPP) \
                        $(ACTM_CPP)
*/
  CLEAN_PREREQ        = $(GEM_CLEAN) $(ACTM_CLEAN)

  COMPILE_PREREQ      = $(GEM_COMPILE)  $(ACTM_COMPILE)

  COMPILE_ONLY_PREREQ = $(GEM_COMPILE_ONLY) $(ACTM_COMPILE_ONLY)

  CPP_PREREQ          = $(GEM_CPP) $(ACTM_CPP)


/*************************
 * Tools external to gem.
 *************************/ 

/***********************************************************************
 * BENCHLib is a Cray MPP library for fast evaluation of transcendental
 * functions.
 ***********************************************************************/

#if (ARCH_OPTION == ARCH_T3E)
  BENCHLibs = -L/usr/local/benchlib -l_scalar -l_vect

#else
  BENCHLibs =

#endif


#if (ARCH_OPTION == ARCH_IBM_SP)

/*********************************************************************
 * MASSLib is an IBM SP library for fast evaluation of transcendental
 * functions.
 *********************************************************************/

#  if (Debug_Option != 1)

#    if (cheetah_machine || frost_machine)
       MASSLibs = -lmass
#    elif (seaborg_machine)
       MASSLibs = -L/usr/common/usg/MASS/3.0/lib -lmass
#    endif

#  else
     MASSLibs =

#  endif

#else
  MASSLibs =

#endif


/*****************************
 * Message passing libraries.
 *****************************/

INCLUDES_MSG =
MSGLibs      =

#if (MSG_OPTION == MSG_MPI)

# if (ARCH_OPTION == ARCH_COMPAQ)
    INCLUDES_MSG = -I/usr/local/apps/mpi/include
    MSGLIBDIR    = /usr/local/apps/mpi/lib

#   if (east_machine || west_machine || north_machine || south_machine)
      MSGLibs    = -L$(MSGLIBDIR) -L/usr/local/lib -lfmpi -lmpi
#   else
      MSGLibs    = -L$(MSGLIBDIR) -L/usr/local/lib \
                   -lfmpi -lmpi -lelan -lelan3 -lrmscall -lmach
#   endif

# elif (ARCH_OPTION == ARCH_IBM_SP)
    INCLUDES_MSG = -I/usr/lpp/ppe.poe/include/thread

# elif (ARCH_OPTION == ARCH_SGI_ORIG)
    MSGLibs      = -lmpi

# elif (ARCH_OPTION == ARCH_INTEL)

# if (HOST_MACH == DISCOVER)
    INCLUDES_MSG = 
    MSGLibs      =
# else
    INCLUDES_MSG = -I/opt/sgi/mpt/1.11-100/include
    MSGLibs      = -L/opt/sgi/mpt/1.11-100/lib -lmpi  -Vaxlib
#endif

# endif

# endif


/********************
 * netCDF libraries.
 ********************/

INCLUDES_NETCDF = -I/usr/local/include
NETCDFLibs      = -L/usr/local/lib -lnetcdf

#if (ARCH_OPTION == ARCH_CRAY)
  INCLUDES_NETCDF = -I/usr/local/pkg/usg/netcdf/3.5/include
  NETCDFLibs      = -L/usr/local/pkg/usg/netcdf/3.5/lib -lnetcdf

#elif (ARCH_OPTION == ARCH_COMPAQ)

# if (halem_machine)
    INCLUDES_NETCDF = -I/usr/local/unsupported/include
    NETCDFLibs      = -L/usr/local/unsupported/lib -lnetcdf
# endif

#elif (ARCH_OPTION == ARCH_IBM_SP) 

# if (seaborg_machine)
  /*INCLUDES_NETCDF = -I/usr/common/usg/netcdf64/3.5/include
    NETCDFLibs      = -L/usr/common/usg/netcdf64/3.5/lib -lnetcdf -lnetcdf_c++*/
    INCLUDES_NETCDF = -I/usr/common/usg/netcdf/3.5/include
    NETCDFLibs      = -L/usr/common/usg/netcdf/3.5/lib -lnetcdf

# endif

#elif (ARCH_OPTION == ARCH_INTEL)

# if (HOST_MACH == DISCOVER)
      INCLUDES_NETCDF = -I/usr/local/other/baselibs/ESMF220rp2_NetCDF362b6_9.1.052/Linux/include/netcdf
      NETCDFLibs      = -L/usr/local/other/baselibs/ESMF220rp2_NetCDF362b6_9.1.052/Linux/lib -lnetcdf
# else
     INCLUDES_NETCDF = -I/local/LinuxIA64/sivo_baselibs/v2_2rp2_10.0.025/LinuxIA64/include/netcdf
     NETCDFLibs      = -L/local/LinuxIA64/sivo_baselibs/v2_2rp2_10.0.025/LinuxIA64/lib -lnetcdf
/*
  INCLUDES_NETCDF = -I/share/ESMA/baselibs/v1_8r1p/LinuxIA64/include/netcdf
  NETCDFLibs      = -L/share/ESMA/baselibs/v1_8r1p/LinuxIA64/lib -lnetcdf
*/
# endif

#elif (ARCH_OPTION == ARCH_SGI_ORIG)

# if (hopper_machine || steger_machine || turing_machine || \
      jimpf0_machine || jimpf1_machine)
    INCLUDES_NETCDF = -I/u/jrt/netcdf/include
    NETCDFLibs      = -L/u/jrt/netcdf/lib -lnetcdf
# endif

#endif

/********************
 * ESMF libraries.
 ********************/

INCLUDES_Esmf =
EsmfLibs      =

#if (ARCH_OPTION == ARCH_CRAY)

  INCLUDES_Esmf =
  EsmfLibs      =

#elif (ARCH_OPTION == ARCH_COMPAQ)

    INCLUDES_Esmf =
    EsmfLibs      =

#elif (ARCH_OPTION == ARCH_IBM_SP)

    INCLUDES_Esmf =
    EsmfLibs      =

#elif (ARCH_OPTION == ARCH_INTEL)

# if (HOST_MACH == DISCOVER)
  INCLUDES_Esmf = -I/usr/local/other/baselibs/ESMF220rp2_NetCDF362b6_9.1.052/Linux/include/esmf
  EsmfLibs      = -L/usr/local/other/baselibs/ESMF220rp2_NetCDF362b6_9.1.052/Linux/lib -lesmf -lcprts -limf -lm -lcxa -lunwind -lrt -ldl -threads
# else
  INCLUDES_Esmf = -I/local/LinuxIA64/sivo_baselibs/v2_2rp2_10.0.025/LinuxIA64/include/esmf
  EsmfLibs      = -L/local/LinuxIA64/sivo_baselibs/v2_2rp2_10.0.025/LinuxIA64/lib -lesmf -lcprts -limf -lm -lcxa -lunwind -lrt -ldl -threads
# endif

#elif (ARCH_OPTION == ARCH_SGI_ORIG)

    INCLUDES_Esmf =
    EsmfLibs      =

#endif

advecINCdir = $(GMIHOME)/Shared/trsa/Linux/include
advecLIBSdir = $(GMIHOME)/Shared/trsa/src/../Linux/lib

advecCoreINC = -I$(advecINCdir)/FVadvcore_GridComp -I$(advecINCdir)/MAPL
advecCoreLibs = $(advecLIBSdir)/libFVadvcore_GridComp.a \
	$(advecLIBSdir)/libFV_Shared.a \
	$(advecLIBSdir)/libMAPL.a \
	$(advecLIBSdir)/libGMAO_pilgrim.a
/*
	$(advecLIBSdir)/libMAPL_cfio_r4.a $(advecLIBSdir)/libMAPL_cfio_r8.a \
*/

/***************
 * Load macros.
 ***************/
/*
GEMMAIN  = $(ESMDIR)/main/esm_main.o
*/

GMIMAIN  = $(GMIDIR)/gem/actm/gmimod/main/GmiMain.o

LoadLibraries = \
  $ $(ACTMLibs)  $(GmiLibs)  \
  $(BENCHLibs)  $(MASSLibs)  $(MSGLibs) $(advecCoreLibs) $(EsmfLibs) $(NETCDFLibs)
/*
  $(BENCHLibs)  $(MASSLibs)  $(MSGLibs) $(advecCoreLibs) $(EsmfLibs) $(NETCDFLibs)
*/


/***************************************************************
 * MACHINE OR SYSTEM DEPENDENT PARAMETERS NEEDED FOR COMPILING.
 ***************************************************************/

/*
INCLUDES = -I$(ESMTOOLSINCLUDE)  $(INCLUDES_MSG)  $(INCLUDES_NETCDF)
*/
INCLUDES = $(INCLUDES_MSG)  $(INCLUDES_NETCDF) $(advecCoreINC)


/**********************
 * COMPAQ / DEC ALPHA.
 **********************/

#if (ARCH_OPTION == ARCH_COMPAQ)

# define RUN_RANLIB  1

# if (Debug_Option == 1)
    CC_DEBUG      = -g
    DEBUG_FLAG    = -g
    LDR_OPTS      = -g
# endif

# if (Optimization_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =
    CC_OPTIMIZE   = -O2
    OPTIMIZE_FLAG = -arch host -fast -fpe -assume accuracy_sensitive
    LDR_OPTS      =
# endif

# if (Profile_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =
    CC_OPTIMIZE   = -O2
    OPTIMIZE_FLAG = -arch host -fast -fpe -assume accuracy_sensitive
    PROFILE_FLAG  = -g3
    LDR_OPTS      =
# endif

  MAKE       = make

  AR         = ar r

  CC         = cc
  CFLAGS     = -c $(CC_DEBUG) $(CC_OPTIMIZE)

  CC_MAKE    = $(CC) -oldcomment -D${CHEMCASE} -DMAKING_MAKEFILE -E

  FC         = f90
  FFLAGS     = -c $(DEBUG_FLAG) $(OPTIMIZE_FLAG) $(PROFILE_FLAG) $(GmiModsDir) $(advecCoreINC) $(INCLUDES_Esmf)

  CPP        = $(FC) -cpp
#if (EXTRAFLAGS == 1)
  CPP_OPTS   = -D${GTcloud} -DMICRO_AEROSOL -DKULNEW -DPARNEW -D${HOSTMACH} -P -Wp,-P $(INCLUDES) \
               $(GmiIncDir)
#endif
#if (EXTRAFLAGS == 0)
  CPP_OPTS   = -D${AerosolCase} -D${HOSTMACH} -P -Wp,-P $(INCLUDES) \
               $(GmiIncDir)
#endif

  LDR        = $(FC)
  LDFLAGS    = $(LDR_OPTS)

#endif


/**************************************
 * CRAY vector machines (SV1/J90/C90).
 **************************************/

#if (ARCH_OPTION == ARCH_CRAY)

# define RUN_RANLIB  0

# if (Debug_Option == 1)
    CC_DEBUG      = -g
    DEBUG_FLAG    = -G0
    LDR_OPTS      =
# endif

# if (Optimization_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =
    CC_OPTIMIZE   = -03
    OPTIMIZE_FLAG = -Oscalar3,task0,vector3,negmsgs,inline1
    LDR_OPTS      =
# endif

# if (Profile_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =
    CC_OPTIMIZE   = -O3
    OPTIMIZE_FLAG = -Oscalar3,task0,vector3,negmsgs,inline1
    PROFILE_FLAG  = -ef -Oscalar3,task0,vector3,negmsgs,inline1
    LDR_OPTS      = -l perf
# endif

  MAKE       = make

  AR         = bld r

  CC         = /opt/ctl/bin/cc
  CFLAGS     = -c $(CC_DEBUG) $(CC_OPTIMIZE)

  CC_MAKE    = $(CC) -D${CHEMCASE} -DMAKING_MAKEFILE -E -Wp"-N -P"

  FC         = /usr/local/Modules/modules/fortran/f90_551J/f90
  FFLAGS     = -c -dp $(DEBUG_FLAG) $(OPTIMIZE_FLAG) $(PROFILE_FLAG) $(GmiModsDir) $(advecCoreINC) $(INCLUDES_Esmf)

  CPP        = /opt/ctl/bin/cpp
#if (EXTRAFLAGS == 1)
  CPP_OPTS   = -D${GTcloud} -DMICRO_AEROSOL -DKULNEW -DPARNEW -D${HOSTMACH} -N -P $(INCLUDES) \
               $(GmiIncDir)
#endif
#if (EXTRAFLAGS == 0)
  CPP_OPTS   = -D${GTcloud} -D${AerosolCase} -D${HOSTMACH} -N -P $(INCLUDES) \
               $(GmiIncDir)
#endif

  LDR        = $(FC)
  LDFLAGS    = $(LDR_OPTS)

#endif


/************
 * CRAY T3E.
 ************/
 
#if (ARCH_OPTION == ARCH_T3E)

# define RUN_RANLIB  0

# if (Debug_Option == 1)
    CC_DEBUG      = -g
    DEBUG_FLAG    = -G0
    LDR_OPTS      =
# endif

# if (Optimization_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =
    CC_OPTIMIZE   = -O2
    OPTIMIZE_FLAG = -O3 -Oaggress -Oinline1
    LDR_OPTS      =
# endif

# if (Profile_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =
    CC_OPTIMIZE   = -O2
    OPTIMIZE_FLAG = -O3 -Oaggress -Oinline1
    PROFILE_FLAG  = -eA
    LDR_OPTS      = -l app
# endif

  MAKE       = make

  AR         = ar r

  CC         = /opt/ctl/bin/cc
  CFLAGS     = -c $(CC_DEBUG) $(CC_OPTIMIZE)

  CC_MAKE    = $(CC) -D${CHEMCASE} -DMAKING_MAKEFILE -E -Wp"-N -P"

  FC         = /opt/ctl/bin/f90
  FFLAGS     = -c -dp $(DEBUG_FLAG) $(OPTIMIZE_FLAG) $(PROFILE_FLAG) $(GmiModsDir) $(advecCoreINC) $(INCLUDES_Esmf)

  CPP        = /opt/ctl/bin/cpp
#if (EXTRAFLAGS == 1)
  CPP_OPTS   = -D${GTcloud} -DMICRO_AEROSOL -DKULNEW -DPARNEW -D${HOSTMACH} -N -P $(INCLUDES) \
               $(GmiIncDir)
#endif
#if (EXTRAFLAGS == 0)
  CPP_OPTS   = -D${GTcloud} -D${AerosolCase} -D${HOSTMACH} -N -P $(INCLUDES) \
               $(GmiIncDir)
#endif

  LDR        = $(FC)
  LDFLAGS    = $(LDR_OPTS)

#endif


/**********
 * IBM SP.
 **********/

#if (ARCH_OPTION == ARCH_IBM_SP)

# define RUN_RANLIB  1

# if (Debug_Option == 1)
    CC_DEBUG      = -g
/*  DEBUG_FLAG    = -g -qfullpath -qflttrap=overflow:zerodivide:invalid:enable -C*/
    DEBUG_FLAG    = -g -qfullpath
    LDR_OPTS      = -g -qsmp=omp -qsmp=noopt
# endif

# if (Optimization_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =
    CC_OPTIMIZE   =
    OPTIMIZE_FLAG = -O3 -qstrict -qarch=auto -qtune=auto -qmaxmem=-1
    LDR_OPTS      = -qsmp=omp
# endif

# if (Profile_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =
    CC_OPTIMIZE   =
    OPTIMIZE_FLAG = -O3 -qstrict -qarch=auto -qtune=auto -qmaxmem=-1
    PROFILE_FLAG  = -g -qfullpath -pg
    LDR_OPTS      = -g -pg -qsmp=omp
# endif

  MAKE       = make

/*AR         = ar -r -X 64*/
  AR         = ar -r

# if (MSG_OPTION == MSG_MPI)
    CC       = /usr/bin/mpcc_r
# else
    CC       = /usr/bin/xlc_r
# endif

/*CFLAGS     = -c -q64 $(CC_DEBUG) $(CC_OPTIMIZE)*/
  CFLAGS     = -c $(CC_DEBUG) $(CC_OPTIMIZE)

  CC_MAKE    = $(CC) -D${CHEMCASE} -DMAKING_MAKEFILE -E -w

# if (MSG_OPTION == MSG_MPI)
    FC       = /usr/bin/mpxlf90_r
# else
    FC       = /usr/bin/xlf90_r
# endif

/*FFLAGS     = -c -q64 -qfixed=132 \
               $(DEBUG_FLAG) $(PROFILE_FLAG) $(OPTIMIZE_FLAG)*/
  FFLAGS     = -c -qfixed=132 $(DEBUG_FLAG) $(PROFILE_FLAG) $(OPTIMIZE_FLAG) $(GmiModsDir) $(advecCoreINC) $(INCLUDES_Esmf)

  CPP        = /usr/ccs/lib/cpp
#if (EXTRAFLAGS == 1)
  CPP_OPTS   = -D${GTcloud} -DMICRO_AEROSOL -DKULNEW -DPARNEW -D${HOSTMACH} -P $(INCLUDES) \
               $(GmiIncDir)
#endif
#if (EXTRAFLAGS == 0)
  CPP_OPTS   = -D${GTcloud} -D${AerosolCase} -D${HOSTMACH} -P $(INCLUDES) \
               $(GmiIncDir)
#endif

  LDR        = $(FC)

/*LDFLAGS  = -q64 $(LDR_OPTS)*/
  LDFLAGS  = -bmaxdata:0x40000000 -bmaxstack:0x4000000 $(LDR_OPTS)

#endif


/*********
 * INTEL.
 *********/

#if (ARCH_OPTION == ARCH_INTEL)

# define RUN_RANLIB  1

# if (Debug_Option == 1)
    CC_DEBUG      = -g -O0
    DEBUG_FLAG    = -g -O0 -traceback -CB
    LDR_OPTS      = -g
# endif

# if (Optimization_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =

# if (HOST_MACH == DISCOVER)
    CC_OPTIMIZE   = -O2
    OPTIMIZE_FLAG = -O2 
    LDR_OPTS      = -O2
# else
    CC_OPTIMIZE   = -O3 -tpp2 -Zp16 -ip
    OPTIMIZE_FLAG = -O2 -tpp2 -Zp16 -ip
    LDR_OPTS      = -O2 -tpp2 -Zp16 -ip 
# endif

# endif

# if (Profile_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =
    CC_OPTIMIZE   =
    OPTIMIZE_FLAG =
    PROFILE_FLAG  =
    LDR_OPTS      =
# endif

  MAKE       = make

  AR         = ar r

  CC         = icc
  CFLAGS     = -c $(CC_DEBUG) $(CC_OPTIMIZE)

  CC_MAKE    = $(CC) -D${CHEMCASE} -DMAKING_MAKEFILE -E

# if (HOST_MACH == DISCOVER)
  FC         = mpif90
# else
  FC         = ifort
# endif 


  FFLAGS     = -c -cm -w95 -WB -warn nogeneral -traceback  \
               $(DEBUG_FLAG) $(OPTIMIZE_FLAG) $(PROFILE_FLAG) $(GmiModsDir) $(advecCoreINC) $(INCLUDES_Esmf)

  CPP        = $(FC) -cpp -132
#if (EXTRAFLAGS == 1)
  CPP_OPTS   = -D${GTcloud} -DMICRO_AEROSOL -DKULNEW -DPARNEW -D${HOSTMACH} -E $(INCLUDES) \
               $(GmiIncDir)
#endif
#if (EXTRAFLAGS == 0)
  CPP_OPTS   = -D${GTcloud} -D${AerosolCase} -D${HOSTMACH} -E $(INCLUDES) \
               $(GmiIncDir)
#endif

  LDR        = $(FC)
  LDFLAGS    = $(LDR_OPTS)

#endif


/**************
 * SGI Origin.
 **************/

#if (ARCH_OPTION == ARCH_SGI_ORIG)

# define RUN_RANLIB  0

# if (Debug_Option == 1)
    CC_DEBUG      = -g 
    DEBUG_FLAG    = -g -n32 -mips4
    LDR_OPTS      = $(DEBUG_FLAG)
# endif

# if (Optimization_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =
    CC_OPTIMIZE   = -O3
    OPTIMIZE_FLAG = -n32 -mips4 -O3 -OPT:Olimit=0:IEEE_arithmetic=3 \
                    -LNO:prefetch=2
    LDR_OPTS      = $(OPTIMIZE_FLAG)
# endif

# if (Profile_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =
    CC_OPTIMIZE   = -O3
    OPTIMIZE_FLAG = -n32 -mips4 -O3 -OPT:Olimit=5000:IEEE_arithmetic=3 \
                    -LNO:prefetch=2
    PROFILE_FLAG  = 
    LDR_OPTS      = $(OPTIMIZE_FLAG) $(PROFILE_FLAG)
# endif

  MAKE       = make

  AR         = ar r

  CC         = cc
  CFLAGS     = -c $(CC_DEBUG) $(CC_OPTIMIZE)

  CC_MAKE    = $(CC) -D${CHEMCASE} -DMAKING_MAKEFILE -E

  FC         = f90
  FFLAGS     = -c $(DEBUG_FLAG) $(OPTIMIZE_FLAG) $(PROFILE_FLAG) $(GmiModsDir) $(advecCoreINC) $(INCLUDES_Esmf)

  CPP        = /usr/lib/cpp
#if (EXTRAFLAGS == 1)
  CPP_OPTS   = -D${GTcloud} -DMICRO_AEROSOL -DKULNEW -DPARNEW -D${HOSTMACH} -P $(INCLUDES) \
               $(GmiIncDir)
#endif
#if (EXTRAFLAGS == 0)
  CPP_OPTS   = -D${GTcloud} -D${AerosolCase} -D${HOSTMACH} -P $(INCLUDES) \
               $(GmiIncDir)
#endif

  LDR        = $(FC)
  LDFLAGS    = $(LDR_OPTS)

#endif


/*********************
 * SUN 4 Workstation.
 *********************/

#if (ARCH_OPTION == ARCH_SUN4)

# define RUN_RANLIB  1

# if (Debug_Option == 1)
    CC_DEBUG      = -g
    DEBUG_FLAG    = -g
    LDR_OPTS      =
# endif

# if (Optimization_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =
    CC_OPTIMIZE   = -O2
    OPTIMIZE_FLAG = -O3
    LDR_OPTS      =
# endif

# if (Profile_Option == 1)
    CC_DEBUG      =
    DEBUG_FLAG    =
    CC_OPTIMIZE   =
    OPTIMIZE_FLAG =
    PROFILE_FLAG  =
    LDR_OPTS      =
# endif

  MAKE       = make

  AR         = ar r

  CC         = /usr/ucb/cc
  CFLAGS     = -c $(CC_DEBUG) $(CC_OPTIMIZE)

  CC_MAKE    = $(CC) -D${CHEMCASE} -DMAKING_MAKEFILE -E

  FC         = /usr/local/SUNWspro/bin/f90
  FFLAGS     = -c -e $(DEBUG_FLAG) $(OPTIMIZE_FLAG) $(PROFILE_FLAG) $(GmiModsDir) $(advecCoreINC) $(INCLUDES_Esmf)

  CPP        = /usr/ccs/lib/cpp
#if (EXTRAFLAGS == 1)
  CPP_OPTS   = -D${GTcloud} -DMICRO_AEROSOL -DKULNEW -DPARNEW -D${HOSTMACH} -P $(INCLUDES) \
               $(GmiIncDir)
#endif
#if (EXTRAFLAGS == 0)
  CPP_OPTS   = -D${GTcloud} -D${AerosolCase} -D${HOSTMACH} -P $(INCLUDES) \
               $(GmiIncDir)
#endif

  LDR        = $(FC)
  LDFLAGS    = $(LDR_OPTS)

#endif


/********************************
 * Set RANLIB to null if needed.
 ********************************/

#if (RUN_RANLIB == 0)
  RANLIB = 
#endif

