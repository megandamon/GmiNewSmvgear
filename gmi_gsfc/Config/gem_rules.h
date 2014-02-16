
/**************************************************************************
 *
 * $Id: gem_rules.h,v 1.2 2011-08-08 20:17:08 mrdamon Exp $
 *
 * CODE DEVELOPER
 *   John Tannahill, LLNL (Original code from Bill Bosl, LLNL)
 *   jrt@llnl.gov
 *
 * FILE
 *   gem_rules.h
 *
 * DESCRIPTION
 *   This file contains suffix dependencies and compilation rules for all
 *   Makefile.cpp files in directories that contain source files.
 *
 *************************************************************************/

/***********************
 * Suffix dependencies.
 ***********************/

.SUFFIXES:
.SUFFIXES: .o .c .f90 .F90 .h


/**********************************************
 * Rules for Fortran source file dependencies.
 **********************************************/

.F90.f90:
	rm -f $*.f90 $*.i90
#if (ARCH_OPTION == ARCH_COMPAQ)
	$(CPP) $(CPP_OPTS) $(DEPEND_INC) $*.F90
	sed '/^ *$$/d' $*.i90 > $*.f90
#elif (ARCH_OPTION == ARCH_INTEL)
	$(CPP) $(CPP_OPTS) $(DEPEND_INC) $*.F90 > $*.i90
	sed '/^ *$$/d' $*.i90 > $*.f90
#else
	$(CPP) $(CPP_OPTS) $(DEPEND_INC) $*.F90 > $*.i90
	sed '/^ *$$/d' $*.i90 > $*.f90
#endif
	rm -f $*.i90
	chmod -w $*.f90


/****************************************************
 * Rules for Fortran and C object file dependencies.
 ****************************************************/

.c.o:
	$(CC) $(CFLAGS) $(DEPEND_INC) $*.c


.F90.o:
	rm -f $*.f90 $*.i90
#if (ARCH_OPTION == ARCH_COMPAQ)
	$(CPP) $(CPP_OPTS) $(DEPEND_INC) $*.F90
	sed '/^ *$$/d' $*.i90 > $*.f90
#elif (ARCH_OPTION == ARCH_INTEL)
	$(CPP) $(CPP_OPTS) $(DEPEND_INC) $*.F90
#else
	$(CPP) $(CPP_OPTS) $(DEPEND_INC) $*.F90 > $*.i90
	sed '/^ *$$/d' $*.i90 > $*.f90
#endif
	rm -f $*.i90
	chmod -w $*.f90

	$(FC) $(FFLAGS) $*.f90


.f90.o:
	$(FC) $(FFLAGS) $*.f90 $(PIPE_DESTINATION)

