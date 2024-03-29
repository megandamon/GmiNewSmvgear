#
# Earth System Modeling Applications (ESMA) base makefile fragment.
# This fragment costumizes ESMF_base.mk for each each architecture. 
#
# REVISION HISTORY:
#
# 06Jun2003  da Silva  First Crack
# 07aug2003  Zaslavsky Added -64 option for IRIX64
# 29Oct2003  Sawyer    Merged back most material lost in 1.26
# 01Dec2003  da Silva  Mods to IRIX64 for building GMAO legacy code
# 02Dec2003  da Silva  Changed INC_MPI to /usr/include (blank causes -I )
# 28Mar2005  Owens/RT  Added courant.
# 01Apr2005  Todling   Changed location of perl on IRIX64 to /usr/bin
# 04Apr2005  da Silva  Changed back perl on IRIX as /usr/bin/perl breaks fdp
# 26Apr2005  da Silva  Merged in Carlos mods for building FVGSI.
# 27Mar2006  da Silva  On IRIX, removed $(OMPFLAG) from FFLAGS and kept its
#                      "-mp" definition; user is supposed to activated where
#                      desired, not by default. Also, added OMPFLAG to Altix.
# 18Apr2006  Stassi    Added MPFLAG for Intel Fortran Compiler
# 05Jan2007  Todling   Add -lsz from palm current baselibs
# 26Jan2007  da Silva  Merged in Kokron's discover port
# 23May2007  Stassi    Removed lf95 as possible Linux default f90 compiler
# 13Aug2008  da Silva  Removed g77 from LIBSCI under Linux.
# 25Oct2008  da Silva  Improved Mac support, made openmpi default on Linux, 
#                      added detection of intelmpi. Slowly, moving towards
#                      use of mpif90. Added F2PY customization.
#
#--------------------------------------------------------------------------

# -----
# TO DO: Remove default for BASEDIR, FC from here; user must set BASEDIR
# -----  and ESMA_FC to control these parameters.
#
  ifndef BASEDIR
     $(warning BASEDIR is undefined --- this will cause an error in the next release!)
  endif

#                             ---------------------
#                             User defined defaults
#                             ---------------------
  ifdef ESMA_FC
	FC := $(ESMA_FC)
  endif
  ifdef ESMA_F2PY
	F2PY := $(ESMA_F2PY)
  endif

  ifeq ($(ESMA_PROFILE),TAU)
##         TAUROOTDIR      = /usr/local/other/tau/tau-2.16.1b3
#         TAUROOTDIR      = /gpfsm/dhome/dkokron/play/tau/tau-2.16.3p1
#         include $(TAUROOTDIR)/x86_64/lib/Makefile.tau-9.1.042-callpath-icpc-mpi-compensate-pdt
##        OPTS  = -optPdtF90Parser=f95parse -optPreProcess -optVerbose -optCompile=$(TAU_F90_SUFFIX) -optKeepFiles -optNoRevert -optCPPOpts="-P -traditional-cpp" -optTauSelectFile=$(ESMADIR)/src/Config/select.tau
#        OPTS  = -optPdtF90Parser=f95parse -optPreProcess -optVerbose -optCompile=$(TAU_F90_SUFFIX) -optKeepFiles -optNoRevert -optDetectMemoryLeaks -optCPPOpts="-P -traditional-cpp" 
#        FC = $(TAU_COMPILER) $(OPTS) $(TAU_F90)
##       FC = tau_f90.sh
  endif

#                               ----
#                               OSF1
#                               ----

ifeq ($(ARCH),OSF1)

  SITE := $(patsubst halem%,nccs,$(SITE))

  ifndef BASEDIR
    ifeq ($(SITE),nccs)
       BASEDIR = /share/ESMA/baselibs/v1_8r1p
    endif
  endif

  BIG_ENDIAN = -convert big_endian

  CC  = cc
  CXX = cxx -x cxx

  BYTERECLEN = -assume byterecl
  EXTENDED_SOURCE = -extend_source
  FREE_SOURCE = -free
  FIXED_SOURCE = -fixed -col72
  fFLAGS   += $(EXTENDED_SOURCE) -OPT:Olimit=0:roundoff=3:reorg_common=OFF -LNO:prefetch=2 
  f90FLAGS += -free -OPT:Olimit=0:roundoff=3:reorg_common=OFF -LNO:prefetch=2 
  FFLAGS   += $(EXTENDED_SOURCE) -cpp 
  F90FLAGS += -free -cpp -OPT:Olimit=0:roundoff=3:reorg_common=OFF -LNO:prefetch=2 

  LIB_SCI  = -lcxml
  LIB_SYS  = -L/usr/ccs/lib/cmplrs/cxx -lcxx

  FREAL4 = -r4
  FREAL8 = -r8

  FDEFS  += $(D)OSF1
  CDEFS  += $(D)OSF1

  OMPFLAG  = -omp

endif  #    OSF1

#                               ---
#                               AIX
#                               ---

ifeq ($(ARCH),AIX)

  SITE := $(patsubst bs%,ncar,$(SITE))

  ifeq ($(SITE),ncar)
    ifndef BASEDIR
	  BASEDIR = /home/bluesky/cacruz/ESMA/baselibs
    endif
  endif

#  CC  = mpcc_r
  CC  = gcc
  CXX = mpcc_r
  CPP = /lib/cpp
  FC  = mpxlf95_r
  AR += -X64

  D = -WF,-D
  DC = -D

  EXTENDED_SOURCE = -qfixed=132
  FIXED_SOURCE = -qfixed=72
  FREE_SOURCE = -qfree

  FREAL4 = -qrealsize=4
  FREAL8 = -qrealsize=8

#  CFLAGS   += -q64 -DAIX 
  fFLAGS   += -q64 $(EXTENDED_SOURCE) -qstrict -qarch=auto 
  f90FLAGS += -q64 -qsuffix=f=f90 -qstrict -qarch=auto
  FFLAGS   += -q64 $(EXTENDED_SOURCE) -qsuffix=cpp=F -qinit=f90ptr -qarch=auto
  F90FLAGS += -q64 -qsuffix=cpp=F90 -qinit=f90ptr -qarch=auto

  INC_MPI = $(BASEDIR)/$(ARCH)/include/mpi

  LDFLAGS  += -q64 -qsmp=omp,noauto
  LIB_MPI   = -lmpi_r -lvtd_r -lmpci_r
  LIB_SCI   = -L/usr/local/lib64/r4i4 -llapack -lessl
  LIB_SYS   = -lxlf90_r -lC_r -brtl

endif  #  AIX


#                               ------
#                               IRIX64
#                               ------

ifeq ($(ARCH),IRIX64)

  SITE := $(patsubst jimpf,nccs,$(SITE))
  SITE := $(patsubst daley,nccs,$(SITE))
  SITE := $(patsubst mintz,nccs,$(SITE))
  SITE := $(patsubst tropic,nccs,$(SITE))
  SITE := $(patsubst courant,nccs,$(SITE))

  ifndef BASEDIR
    ifeq ($(SITE),nccs)
       BASEDIR = /share/ESMA/baselibs/v1_8r1p
    else
       BASEDIR = /share/ESMA/baselibs/latest
    endif
  endif

  ifeq ($(SITE),nccs)
                                               # /usr/bin/perl breaks fdp
     PERL = /ford1/local/bin/perl
                                               # GEOS-5 macros need GNU cpp
     CPP = /ford1/local/irix6.2/gnu/bin/cpp
  endif

  CC  = cc
  CXX = c++

  FPP = f90 
  FDEFS += $(D)HAVE_SHMEM
  FPPFLAGS = -E $(FDEFS) $(FINCS) 

  OMPFLAG  = -mp

  FOPT += -OPT:Olimit=0

  EXTENDED_SOURCE = -extend_source
  FREE_SOURCE = -free
  FIXED_SOURCE = -fixed
  fFLAGS   += -64 -fixedform $(EXTENDED_SOURCE) -cpp
  f90FLAGS += -64 $(EXTENDED_SOURCE) -cpp
  FFLAGS   += -64 $(EXTENDED_SOURCE) -cpp -macro_expand 
  F90FLAGS += -64 -cpp -macro_expand

  CFLAGS   += -64
  CXXFLAGS += -64
  LDFLAGS  += -64 

  LIB_MPI = -L$(MPT_SGI)/usr/lib64 -lmpi -lmpi++ 

#  LIB_SCI =  -lscs -lscs_blas -lftn -lc -lC -lcomplib.sgimath -lm 
#  LIB_SCI =  -lftn -lc -lC -lCio -lcomplib.sgimath -lm
  LIB_SCI =  -lscs

  LIB_SYS = -L/usr/lib64 -lc -lC -lCio

  INC_MPI = ${MPT_SGI}/usr/include

  FDEFS  += $(D)SGI $(D)IRIX64
  CDEFS  += $(D)SGI $(D)IRIX64

endif  #    IRIX64

#                               -----
#                               Linux
#                               -----

ifeq ($(ARCH),Linux)

# Linux default compilers
# -----------------------
  ifndef ESMA_FC
     FC := ifort
  endif
  CC  = cc
  CXX = c++
  CPP = cpp

# Determine which site we are at
# ------------------------------
  SITE := $(patsubst cerebus.gsfc.nasa.gov,nsipp,$(SITE))
  SITE := $(patsubst nsipp02.gsfc.nasa.gov,nsipp,$(SITE))
  SITE := $(patsubst ping-ge.sci.gsfc.nasa.gov,nsipp,$(SITE))
  SITE := $(patsubst boygo-ge.sci.gsfc.nasa.gov,nsipp,$(SITE))
  SITE := $(patsubst boygo.gsfc.nasa.gov,nsipp,$(SITE))
# SITE := $(patsubst columbia.nas.nasa.gov,nas,$(SITE))

  SITE := $(patsubst columbia,nas,$(SITE))
  SITE := $(patsubst cfe1,nas,$(SITE))
  SITE := $(patsubst cfe2,nas,$(SITE))
  SITE := $(patsubst palm,nccs,$(SITE))

#
#                    Linux Site Specific
#                    -------------------

# NSIPP specific
# --------------
  ifeq ($(SITE),nsipp)
    ifndef BASEDIR
       BASEDIR = /home/trayanov/baselibs/latest
    endif
                # Absoft f90 is the default at NSIPP
    ifndef ESMA_FC
      FC := f90
    endif
  endif
  ifeq ($(SITE),cumulus.gsfc.nasa.gov)
       BASEDIR = /home/bacmj/esmf_latest
       ifndef ESMA_FC
           FC := f90
       endif
  endif
  ifeq ($(SITE),beer.gsfc.nasa.gov)
       BASEDIR = /home/suarez/lib/baselibs/latest
       ifndef ESMA_FC
           FC := f90
       endif
  endif

# NAS specific
# ------------
  ifeq ($(SITE),nas)
    ifndef BASEDIR
       BASEDIR = /u/mirvis/v1_8r1p
    endif
  endif

  ifeq ($(SITE),nccs)
    ifndef BASEDIR
       BASEDIR = /share/ESMA/baselibs/v1_8r1p
    endif
  endif

# Linux MPI default: openmpi
# --------------------------
  INC_MPI := $(dir $(shell which mpif90))../include
  LIB_MPI := -L$(dir $(shell which mpif90))../lib -lmpi -lmpi_cxx

#
#                    Linux Compiler Specific
#                    -----------------------

# Absoft compiler
# ---------------
  ifeq ($(FC), f90) 
      EXTENDED_SOURCE = -W 132
      FREE_SOURCE = -f free
      FIXED_SOURCE = -f fixed
      LIB_ESMF = $(BASELIB)/esmf/libesmf.a #$(BASELIB)/esmf/libnetcdf_stubs.a
      LIB_MPI = -L$(BASELIB)/mpi -lpmpich++ -lfmpich -lmpich -lmpichfsup
      LIB_SYS = -lstdc++ -lpthread -lU77 -lrt
      FREAL4 =
      FREAL8 = -N113
      M = -p
      FINCS += -I$(BASEINC)/mpi
      FDEFS += -DABSOFT -DNO_R16 
      XFLAGS += -YEXT_NAMES=LCS -YEXT_SFX=_ $(EXTENDED_SOURCE)
      FOPT   = -g -trap=INVALID,DIVBYZERO,OVERFLOW
  endif

# Intel Fortran Compiler (ifort or mpiifort)
# ------------------------------------------
  ifeq ($(subst mpi,,$(FC)), ifort) 

    EXTENDED_SOURCE := -extend_source
    FREE_SOURCE := -free
    FIXED_SOURCE := -fixed
    MPFLAG  := -mp
    OMPFLAG  := -openmp
    BIG_ENDIAN := -convert big_endian
    BYTERECLEN := -assume byterecl
    FPE = -fpe0
    ALIGNCOM = -align dcommons
    FREAL4 =
    FREAL8 = -r8
    FOPT2 += 
    ifeq ("$(BOPT)","g")
       FOPT = $(FOPTG)
    else
       FOPT = $(FOPT3)
    endif

    LIB_ESMF = $(BASELIB)/libesmf.a

    CC  = gcc
    CXX = g++

#   Handle MPI; need more robust solution for i686
#   ----------------------------------------------
    ifdef I_MPI_ROOT
        FC := mpiifort
        ifeq ($(MACH), x86_64) 
          INC_MPI := $(I_MPI_ROOT)/include64
          LIB_MPI := -L$(I_MPI_ROOT)/lib64  -lmpi -lmpiif # Intel MPI
        else
          INC_MPI := $(I_MPI_ROOT)/include
          LIB_MPI := -L$(I_MPI_ROOT)/lib  -lmpi -lmpiif # Intel MPI
        endif
    else
    ifdef FPATH
        FPATHS := $(subst :, ,$(FPATH))
        ifeq ($(MACH), x86_64) 
          INC_MPI := $(filter /opt/scali%,$(FPATHS)) 
          LIB_MPI := -L$(subst include,lib,$(INC_MPI)) -lfmpi -lmpi -lmpio 
        endif
        ifeq ($(MACH), ia64)
          INC_MPI := $(filter /opt/sgi/mpt%,$(FPATHS)) \
                     $(filter /nasa/sgi/mpt%,$(FPATHS)) 
          INC_MPI := $(word 1,$(INC_MPI))
          LIB_MPI := -L$(subst include,lib,$(INC_MPI)) -lmpi -lmpi++
         endif
    else 
    endif
    endif

#   Define LIB_SYS
#   --------------
    LIB_SCI := 
    LIB_SYS := -ldl -lc -lpthread -lrt 

#   Determine compiler version
#   --------------------------
    IFORT_VER := $(subst ., ,$(word 3,$(shell ifort --version)))
    IFORT_MAJOR := $(word 1,$(IFORT_VER))
    IFORT_MINOR := $(word 2,$(IFORT_VER))

    ifeq ($(IFORT_MAJOR), 10)
          LIB_SYS := -lirc -lguide $(LIB_SYS)
          ifneq ($(MACH), i686)
              FPE += -fp-model precise
          endif
    else
          LIB_SYS +=  -lunwind -lcprts
    endif

#   MKL math library
#   ----------------
    ifdef MKLPATH
       LIB_SCI += -L$(MKLPATH) -lmkl_lapack -lmkl 
    else
    ifeq ($(MACH), ia64)
       LIB_SCI += -lscs 
    endif
    endif

#   Costumize for each MACH
#   -----------------------
    GCC_DIR = $(shell dirname `gcc --print-libgcc-file-name`)
    ifeq ($(MACH), x86_64) 
       OVERRIDE_LIMITS =
       OMPFLAG =
       LOOP_VECT = 
       FDEFS += $(D)HAVE_SHMEM
    else
    ifeq ($(MACH), ia64)
      OVERRIDE_LIMITS = -override_limits 
      FDEFS += $(D)HAVE_SHMEM
    else
    endif # x86_64
    endif # ia64

    LIB_SYS += -L$(GCC_DIR) -lstdc++

    fFLAGS += $(EXTENDED_SOURCE) $(FPE) $(OVERRIDE_LIMITS) $(ALIGNCOM)
    FFLAGS += $(EXTENDED_SOURCE) $(FPE) $(OVERRIDE_LIMITS) $(ALIGNCOM)
    f90FLAGS += $(FPE) $(OVERRIDE_LIMITS) $(ALIGNCOM)
    F90FLAGS += $(FPE) $(OVERRIDE_LIMITS) $(ALIGNCOM)

#   This should be done in general...
#   ---------------------------------
    ifeq ($(MACH), i686)
	FC := mpif90
    endif

  endif

# Lahey F95 Compiler
# ------------------
  ifeq ($(ESMA_FC), lf95) 

      EXTENDED_SOURCE := --wide
      FREE_SOURCE = -free
      FIXED_SOURCE = -fixed
      FREAL4   := 
      FREAL8   := -CcdRR8  

      OMPFLAG = --openmp

      fFLAGS   += --warn $(EXTENDED_SOURCE)
      f90FLAGS += --warn $(EXTENDED_SOURCE)
      FFLAGS   += --warn $(EXTENDED_SOURCE)
      F90FLAGS += --warn $(EXTENDED_SOURCE)
      LDFLAGS  += --warn $(EXTENDED_SOURCE)

      LIB_SCI  = -llapackmt -lblasmt 
      LIB_SYS = -ldl -lc -lpthread -lrt 

#     AssumeOpenMPI is installed
#     --------------------------
      INC_MPI = /usr/include
      LIB_MPI = -lmpi
      override FC = mpif90
      OMPI_FC = lf95
      export OMPI_FC 

  endif

# GNU Fortran Compiler
# --------------------
  ifeq ($(FC), gfortran) 

      CC = gcc

      LIB_ESMF = $(BASELIB)/libesmf.a

      EXTENDED_SOURCE := -ffixed-line-length-132
      FREE_SOURCE = 
      FIXED_SOURCE = -ffixed-form
      FREAL4   := 
      FREAL8   := -fdefault-real-8

      OMPFLAG = 

      fFLAGS   += $(D)__GFORTRAN__ $(EXTENDED_SOURCE)
      FFLAGS   += $(D)__GFORTRAN__ $(EXTENDED_SOURCE)
      f90FLAGS += $(D)__GFORTRAN__ -ffree-line-length-256
      F90FLAGS += $(D)__GFORTRAN__ -ffree-line-length-256

      INC_MPI = /usr/include
      LIB_MPI = -lmpi

#      LIB_SCI  = -llapackmt -lblasmt 
      LIB_SYS = -ldl -lc -lpthread -lrt -lstdc++

  endif


# Wrapper around PGI compiler for the Cray XT-4 platform
# ------------------------------------------------------
  ifeq ($(FC),ftn)

#   Query: what is this? Why are you duplicating ESMA_base here?
#   Please clean up.
#   ------------------------------------------------------------
    ESMF_ARCH = Linux
    ESMABIN = $(ESMADIR)/bin
    ESMALIB = $(ESMADIR)/lib
    ESMAINC = $(ESMADIR)/include
    ESMAMOD = $(ESMADIR)/include
    ESMAETC = $(ESMADIR)/etc
    ESMADOC = $(ESMADIR)/doc
    BASEBIN = $(BASEDIR)/bin
    BASELIB = $(BASEDIR)/lib
    BASEINC = $(BASEDIR)/include
    BASEMOD = $(BASEDIR)/include
    BASEETC = $(BASEDIR)/etc
    DIR_NETCDF = $(BASEDIR)/
    INC_NETCDF = $(DIR_NETCDF)/include/netcdf
    DIR_HDF = $(BASEDIR)/
    INC_HDF = $(DIR_HDF)/include/hdf
#    INC_SDF = $(INC_NETCDF)
#    LIB_SDF = $(LIB_NETCDF)
    INC_SDF = $(INC_HDF)
    LIB_SDF = $(LIB_HDF)
    INC_ESMF = $(DIR_ESMF)/include/esmf
    MOD_ESMF = $(DIR_ESMF)/include/esmf
    LIB_ESMF = $(BASELIB)/libesmf.a
    INC_MPI = $(MPICH_DIR)/include
    LIB_MPI   = -L$(MPICH_DIR)/lib -lmpichf90

    LIB_SYS = -lrt -lC
    EXTENDED_SOURCE = -Mextend
    FREAL4   = -r4
    FREAL8   = -r8
    ifeq ("$(BOPT)","g")
#       FOPT   = $(FOPTG) -Ktrap=fp -Mbounds -Mchkptr
#       FOPT   = $(FOPTG) -O0 -Ktrap=fp -Mbounds
       FOPT   = $(FOPTG) -O0 -Ktrap=fp
    else
       FOPT   = $(FOPT3)
#       FOPT   = -fastsse $(FOPT3)
#       FOPT   = $(FOPT2)
    endif
    BIG_ENDIAN =
    fFLAGS += $(EXTENDED_SOURCE)
    FFLAGS += $(EXTENDED_SOURCE)
    f90FLAGS +=
    F90FLAGS +=
    CFLAGS   += -DpgiFortran
    CXXFLAGS +=
    CC  := cc
    CXX := CC
#    TAUROOTDIR      = /ccs/home/dkokron/play/tau/tau-2.16.5
#    include $(TAUROOTDIR)/xt3/lib/Makefile.tau-callpath-mpi-compensate-pdt-pgi
#    CC              = $(TAU_CC)
#    OPTS            = -optPdtF90Parser=f95parse -optPreProcess -optVerbose -optCompile=$(TAU_F90_SUFFIX) -optKeepFiles -optNoRevert -optCPPOpts="-P -traditional-cpp" -optTauSelectFile=${ESMADIR}/src/Config/select.tau
#    FC              = $(TAU_COMPILER) $(OPTS) $(TAU_F90)
#    CF              = $(FC)
  endif

endif  #    Linux

#                               -----------------
#                               Darwin (Mac OS X)
#                               -----------------

ifeq ($(ARCH),Darwin)

# Linux default compilers
# -----------------------
  ifndef ESMA_FC
     FC = gfortran
  endif
  CC  = gcc
  CXX = g++
  CPP = /usr/bin/cpp -xc++
  SED = /usr/bin/sed
  MAKEFLAGS = -w
  DLLEXT = dylib
  LIB_HDF = -lmfhdf -ldf -ljpeg -lz -lsz 
  LIB_SYS = 
  OMPFLAG = 

#
#                    Darwin Site Specific
#                    -------------------

# Darwin MPI default: assumes openmpi
# -----------------------------------
  INC_MPI := $(dir $(shell which mpif90))../include
  LIB_MPI := -L$(dir $(shell which mpif90))../lib -lmpi -lmpi_cxx

#                    Darwin Compiler Specific
#                    -----------------------


# GNU Fortran Compiler
# --------------------
  ifeq ($(FC), gfortran) 

      LIB_ESMF = $(BASELIB)/libesmf.a

      EXTENDED_SOURCE := -ffixed-line-length-132
      FREE_SOURCE = 
      FIXED_SOURCE = -ffixed-form
      FREAL4   := 
      FREAL8   := -fdefault-real-8

      fFLAGS   += $(D)__GFORTRAN__ $(EXTENDED_SOURCE)
      FFLAGS   += $(D)__GFORTRAN__ $(EXTENDED_SOURCE)
      f90FLAGS += $(D)__GFORTRAN__ -ffree-line-length-256
      F90FLAGS += $(D)__GFORTRAN__ -ffree-line-length-256

  endif


# G95 Fortran Compiler
# --------------------
  ifeq ($(FC), g95) 

      LIB_ESMF = $(BASELIB)/libesmf.a

      EXTENDED_SOURCE := -ffixed-line-length-132
      FREE_SOURCE = 
      FIXED_SOURCE = -ffixed-form
      FREAL4   := 
      FREAL8   := -r8 -i4

      fFLAGS   +=  $(EXTENDED_SOURCE)
      FFLAGS   +=  $(EXTENDED_SOURCE)
      f90FLAGS +=  -ffree-line-length-huge
      F90FLAGS +=  -ffree-line-length-huge

  endif # g95

# Intel Fortran Compiler
# ----------------------
  ifeq ($(FC), ifort) 

      EXTENDED_SOURCE := -extend-source
      FREE_SOURCE = -free
      FIXED_SOURCE = -fixed
      FREAL4 := 
      FREAL8 := -r8 -i4

      fFLAGS += $(EXTENDED_SOURCE)
      FFLAGS += $(EXTENDED_SOURCE)

      LIB_SYS = -limf -lm -ldl -lirc -lguide -lstdc++ -lgcc_s.1

  endif # ifort

  FC := mpif90

endif  #    Darwin
