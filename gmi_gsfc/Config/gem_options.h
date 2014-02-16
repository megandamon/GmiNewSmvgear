
/****************************************************************************
 *
 * $Id: gem_options.h,v 1.2 2011-08-08 20:17:08 mrdamon Exp $
 *
 * CODE DEVELOPER
 *   John Tannahill, LLNL (Original code from Bill Bosl, LLNL)
 *   jrt@llnl.gov
 *
 * FILE
 *   gem_options.h
 *
 * DESCRIPTION
 *   The options specified in the following #define, etc. clauses are used
 *   to build the desired executable program.
 *
 ****************************************************************************/


/* Choose package options. */

#define ACTM_Package    1       /* Current options:     0 = no package
                                                        1 = Gmimod   */

/* Select Gmimod compilation options using either #define or #undef. */

#undef  UCITRANS
#define SMV2CHEM


/* Select debug, optimization, and profiling options. */

#define Debug_Option         0
#define Optimization_Option  1
#define Profile_Option       0

