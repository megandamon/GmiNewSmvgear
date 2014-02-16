#!/usr/bin/perl
#
# Simple utility to parse gmake output and print error summary;
# loosely inspired on Perl's Harness module.
#
#----------------------------------------------------------------------------

use Env;                 # make env vars readily available
use Getopt::Std;         # command line options

# Command line options
# --------------------
  getopts('hp:s');
  usage() if ( $opt_h || $#ARGV < 0 );

  $NAME = $ARGV[0];
  $GCPROXY = "GEOS_GenericMod" unless $GCPROXY = "$opt_p";
  unless ( $SETSERVICES = $opt_s ) {
      if ( "$GCPROXY" eq "GEOS_GenericMod" ) {
           $SETSERVICES = "GEOS_GenericSetServices" ;
      } else {  
           $SETSERVICES = "SetServices" ; # work for most components
      }
  }

  print <<EOF;
!
! Stub code automatically generated on the fly by geos_stub.pl; 
! do not edit or check in, change GNUmakefile instead.
!
module $NAME
   use $GCPROXY,  only:  SetServices => $SETSERVICES
end module $NAME
EOF

exit 0;


#......................................................................

sub usage {

   print <<"EOF";

NAME
     stub.pl - Creates stub ESMF Grid Component
          
SYNOPSIS

     stub.pl [-hps] module_name
          
DESCRIPTION

     Creates a simple stub ESMF Grid Component which inherits the
     Initialize/Run/Finalize methods from another component, usually
     GEOS Generic.

OPTIONS
     -p proxy_name  The proxy Grid Component name; default is 
                    "GEOS_GenericMod"
     -s             SetServices entry point in proxy; default depends on the
                    proxy name: if proxy name is "GEOS_GenericMod" it 
                    defaults to "GEOS_GenericSetServices", otherwise it 
                    defaults to simply "SetServices".

BUGS
     For consistency with current practice in GEOS-5, the file name
     does not correspond to the actual f90 module name, or even the
     package name, example:
       Package name: GEOSagcs_GridComp
        Module name: GEOS_AgcsGridCompMod
          File name: GEOS_AgcsGridComp
     The actual ESMF convention calls for the following convention:
       Package name: GEOSagcs_GridComp
        Module name: GEOSagcs_GridCompMod
          File name: GEOSagcs_GridCompMod.F90
     or even
       Package name: GEOS_AgcsGridComp
        Module name: GEOS_AgcsGridCompMod
          File name: GEOS_AgcsGridCompMod.F90
     For this reason, this utility echoes the results to stdout and
     the actual file name convention is enforced there.

AUTHOR
     Arlindo.daSilva@nasa.gov 

EOF

  exit(1)

 }

