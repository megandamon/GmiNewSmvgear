#!/usr/bin/perl
#
# Reads a Grid Component Registry with definitions of INTERNAL, IMPORT and
# EXPORT states and generate the necessary code fragments and resource files
# with the appropriate MAPL calls for setting up such states, producing
# output, etc.
#
# REVISION HISTORY:
# 16Aug2006  da Silva  Initial version borrowing code from fdp.
#
#........................................................................

use File::Basename;
use Getopt::Std;         # command line options

$Iam = "mapl_acg";

# Defaults and such
# -----------------
  init();

# Main loop over file
# -------------------
$pFlag = 0;
$iFlag = 0;
$xFlag = 0;

LINE: while (<>) {

    chomp($_);

    next LINE if (/^\#/ );  # skip comment lines

#   Beginning of table?
#   -------------------
    if ( /^\<$gcname\:\:InternalSpec/ ) {
	$pFlag = 1;
	next LINE;
    }
    if ( /^\<$gcname\:\:ImportSpec/ ) {
	$iFlag = 1;
	next LINE;
    }
    if ( /^\<$gcname\:\:ExportSpec/ ) {
	$xFlag = 1;
	next LINE;
    }

#   End of table?
#   ------------
    if ( /^\<\/$gcname\:\:InternalSpec/ ) {
	$pFlag = 0;
	next LINE;
    }
    if ( /^\<\/$gcname\:\:ImportSpec/ ) {
	$iFlag = 0;
	next LINE;
    }
    if ( /^\<\/$gcname\:\:ExportSpec/ ) {
	$xFlag = 0;
	next LINE;
    }

#   Record content of table
#   -----------------------
    if ( $pFlag or $iFlag ) {

	die "Internal/Import tables not implemented yet";

#   Export table record
#   -------------------
    } elsif ( $xFlag ) {

         ($name,$units,$dims,$vloc,$stat,$rfr,$avg,$subt,$long) 
         = split('\|');

         $record = { name  => trim($name),
                     units => trim($units),
                     dims  => trim($dims),
                     vloc  => trim($vloc),
                     stat  => trim($stat),
                     rfr   => trim($rfr),
                     avg   => trim($avg),
                     subt  => trim($subt),
                     long  => trim($long) };

	 push @xTable, $record;

	 next LINE;
       
    } 

}

#   Parsing error
#   -------------
    die "Internal Table not terminated" if ( $pFlag );
    die   "Import Table not terminated" if ( $iFlag );
    die   "Export Table not terminated" if ( $xFlag );

#   Write code fragments with State Specs
#   -------------------------------------
    write_InternalSpec() if ( @pTable );
    write_ImportSpec()   if ( @iTable );
    write_ExportSpec()   if ( @xTable );

#   Write fragment with Declarations and GEOS_GetPointer()
#   ------------------------------------------------------
    write_GetPointer()   if ( @pTable+@iTable+@xTable );

#   Write fragments for history 
#   ---------------------------
    write_History() if ( @pTable+@iTable+@xTable );

exit(0);

#.........................................................................

sub init {

# Command line options
# --------------------
  getopts('CFvmdn:N:i:x:g:h:');
  usage() if ( $opt_h || $#ARGV < 0 );

 $ifilen = $ARGV[0] unless ( $ifilen = $opt_i );
($base,$path,$suffix) = fileparse($ifilen,'\..*');

$gcname  = $base unless ($gcname = $opt_n );

$gcname =~ s/_REGISTRY//g;
$gcname =~ s/_Registry//g;
$gcname =~ s/_registry//g;

$gcName = $gcname unless ( $gcName = $opt_N );
$mangle = $opt_m;

$pFilen = "$gcname"."_InternalSpec___.h" unless ( $pFilen=$opt_p );
$iFilen = "$gcname"."_ImportSpec___.h"   unless ( $pFilen=$opt_i );
$xFilen = "$gcname"."_ExportSpec___.h"   unless ( $pFilen=$opt_x );

$gFilen = "$gcname"."_GetPointer___.h"   unless ( $pFilen=$opt_g );
$hFilen = "$gcname"."_History___.rc"     unless ( $pFilen=$opt_h );

}

#.........................................................................

sub write_InternalSpec {

	die "Internal Spec not implemented yet";
}

#.........................................................................

sub write_ImportSpec {

	die "Import Spec not implemented yet";

}

#.........................................................................

sub write_ExportSpec {

    open(FILE,">$xFilen") or die "cannot open ExportSpec file";

    print "Building EXPORT Spec file $xFilen\n" if ( $opt_v ); 

    Preamble('!');

    $n = @xTable - 1;
    for $i ( 0..$n ) {

        %r = %{$xTable[$i]};

        $name  = $r{name};
        $long  = $r{long};
        $units = $r{units};

        if ( $r{dims} eq "xyz" ) { $dims = "GEOS_DimsHorzVert"; }
	else {                     $dims = "GEOS_DimsHorzOnly"; }

        if ( $r{vloc} eq "C" ) { $vloc = "GEOS_VLocationCenter"; }
     elsif ( $r{vloc} eq "E" ) { $vloc = "GEOS_VLocationEdge"; }
	else                   { $vloc = "GEOS_VLocationNone"; }

        print FILE <<EOF;

     call GEOS_StateAddExportSpec(GC,  &
        SHORT_NAME         = '$name',  &
        LONG_NAME          = '$long',  &
        UNITS              = '$units', &
        DIMS               = $dims,    &
        VLOCATION          = $vloc,    &
EOF

     $stat = $r{stat};

     if ( "$stat" eq "B" ) { $stat = 'GEOS_BundleItem'; } # what else?


     if ( "$stat" ne "" ) {
        print FILE <<EOF;
        STAT               = $stat,    &
EOF
     }

     $rfr  = $r{rfr};
     if ( "$rfr" ne "" ) {
        print FILE <<EOF;
        REFRESH_INTERVAL   = $rfr,     &
EOF
     }

     $avg  = $r{avg};
     if ( "$avg" ne "" ) {
        print FILE <<EOF;
        AVERAGING_INTERVAL = $avg,     &
EOF
     }

     $subt = $r{subt};
     if ( "$subt" ne "" ) {
        print FILE <<EOF;
        NUM_SUBTILES       = $subt,    &
EOF
     }

#    End of item
#    -----------
     print FILE <<EOF
                                                       RC=STATUS  )
     VERIFY_(STATUS)

EOF

    }

    close(FILE);

}



#.........................................................................
sub trim {
    $str = shift;
    @words = split(' ',$str);
    $str = join ' ', @words;
    return $str;
}

#.........................................................................

sub write_GetPointer {

    if ( $opt_F ) { write_GetPointer_FlatArray(); }
    else          { write_GetPointer_ChemArray(); }

}

#.........................................................................

sub write_GetPointer_FlatArray {

    open(FILE,">$gFilen") or die "cannot open ExportSpec file";

    print "Building F90 code fragment file $gFilen (flat arrays)" if ( $opt_v ); 

    Preamble('!');

#   Declarations:
#   -------------
    print FILE <<EOF;

!       Local arrays referencing the Import/Export states
!       -------------------------------------------------
EOF

    $n = @xTable - 1;
    for $i ( 0..$n ) {

        %r = %{$xTable[$i]};

        $name  = $r{name};
        $long  = $r{long};

        if ( $r{dims} eq "xyz" ) { $rank = "(:,:,:)"; }
	else {                     $rank = "(:,:)  "; }

        print FILE <<EOF;
        real, pointer, dimension$rank :: $name ! Export: $long        
EOF
        print "#define $name \n" if ( $opt_d );

       }

        print FILE <<EOF;

!       Get pointers to data in state
!       -----------------------------
EOF

#   Get the pointer - non-binned variables
#   --------------------------------------
    $n = @xTable - 1;
    for $i ( 0..$n ) {

        %r = %{$xTable[$i]};
        $name  = $r{name};

        print FILE <<EOF;
        call GEOS_GetPointer ( EXPORT, $name,  '$name', RC=STATUS )
        VERIFY_(STATUS)
EOF
    }

    close(FILE);

}

#.........................................................................

sub write_GetPointer_ChemArray {

    open(FILE,">$gFilen") or die "cannot open ExportSpec file";

    if ( $opt_C ) {

	print "Building F90 code fragment file $gFilen (Chem Arrays)\n" 
              if ( $opt_v ); }

    else {

	print "Building F90 code fragment file $gFilen\n" 
              if ( $opt_v ); 

    }

    Preamble('!');

#   First pass: look for binned variables
#   -------------------------------------
    $n = @xTable - 1;
    for $i ( 0..$n ) {

        %r = %{$xTable[$i]};

        $name  = $r{name};
        $long  = $r{long};
        $dims  = $r{dims};

#       Binned variables, a little trickier
#       -----------------------------------
        if ( $long =~ / Bin / ) {

	    @tokens = split('Bin',$long);
            $title = $tokens[0];
            $nBin = $tokens[1];
            $nBin =~ s/ //g;           

            $name =~ s/$nBin//g;
            $nBin = $nBin + 0;

#           For this work bins must be in ascending order
#           --------------------------------------------
            if ( "$dims" eq "xyz" ) { $v3Bin{$name} = $nBin; }
            else                    { $v2Bin{$name} = $nBin; }

            $vBin{$name} = $nBin; 
            $tBin{$name} = $title; 

            $xName[$i] = $name;
            $xBin[$i] = $nBin;

        } else {

            $xName[$i] = $name;  # needed for -C option
        }
    }

#   Declarations: bin sizes
#   -----------------------
    print FILE <<EOF;

!     Bin sizes
!     ---------
EOF
    foreach $name ( keys %vBin ) {
	$n = $vBin{$name};
        print FILE <<EOF;
      integer, parameter              :: NBIN_$name = $n ! $tBin{$name}
EOF
      print "#define $name \n" if ( $opt_d );
      }

#   Declarations: binned variables
#   ------------------------------
    print FILE <<EOF;

!     Bin-indexed Chem Arrays
!     -----------------------
EOF

    foreach $name ( keys %vBin ) {
      $n = $vBin{$name};

      print FILE <<EOF;
      type(Chem_Array), target        ::    $name(NBIN_$name) ! Export: $tBin{$name}
      type(Chem_Array), pointer       :: ptr$name(:)  ! Export: $tBin{$name}
EOF
      print "#define ptr$name \n" if ( $opt_d );

  }


#   Declarations: non-binned variables
#   ----------------------------------
    print FILE <<EOF;

!     Local array referencing the Import/Export states
!     ------------------------------------------------
EOF

    $n = @xTable - 1;
    for $i ( 0..$n ) {

        %r = %{$xTable[$i]};

        $name  = $r{name};
        $long  = $r{long};

        if ( $r{dims} eq "xyz" ) { $rank = "(:,:,:)"; }
	else {                     $rank = "(:,:)  "; }

#       Binned variables
#       ----------------
        unless ( $long =~ / Bin / ) {

          if ( $opt_C ) {
	      print FILE <<EOF;
      type(Chem_Array), target        ::    $name ! Export: $long
      type(Chem_Array), pointer       :: ptr$name ! Export: $long
EOF
          } else {         
              print FILE <<EOF;
      real, pointer, dimension$rank :: $name ! Export: $longEOF
EOF
          }
          print "#define ptr$name \n" if ( $opt_d );
          }
    
    }

        print FILE <<EOF;

!     Get pointers to data in state
!     -----------------------------
EOF

#   Get the pointer - non-binned variables
#   --------------------------------------
    $n = @xTable - 1;
    for $i ( 0..$n ) {

        %r = %{$xTable[$i]};

        $name  = $r{name};
        $long  = $r{long};
        $dims  = $r{dims};

        if ( $r{dims} eq "xyz" ) { $rank = "3d"; }
	else {                     $rank = "2d"; }

#       Binned variables are Chem Arrays
#       --------------------------------
        if ( $long =~ / Bin / ) {
	   print FILE "\n      ptr$xName[$i] => $xName[$i]   ! $long\n"
                                                       if ( $xBin[$i] == 1 ); 
           print FILE <<EOF;
      call GEOS_GetPointer ( EXPORT, $xName[$i]($xBin[$i])%data$rank,  '$name', RC=STATUS )
      VERIFY_(STATUS)
EOF

#       If desired, even non-binned variables are Chem Arrays
#       -----------------------------------------------------
        } elsif ( $opt_C ) {
           print FILE <<EOF;

      ptr$xName[$i] => $xName[$i]   ! $long
      call GEOS_GetPointer ( EXPORT, $xName[$i]%data$rank,  '$name', RC=STATUS )
      VERIFY_(STATUS)
EOF

        } else {
            print FILE <<EOF;
      call GEOS_GetPointer ( EXPORT, $name,  '$name', RC=STATUS )
      VERIFY_(STATUS)
EOF
        }
    }

    close(FILE);

}

#.........................................................................

sub write_History {

    open(FILE,">$hFilen") or die "cannot open ExportSpec file";

    print "Building History Spec fragment file $hFilen\n" if ( $opt_v ); 

    $Name = "$gcName"."::"."$name";

    Preamble('#');

#   2D quantities
#   -------------

        print FILE <<EOF;

  list(#)%filename:   '/dev/null/%s.$gcName.sfc.%y4%m2%d2_%h2%n2z',
  list(#)%format:     'CFIO',
  list(#)%mode:       'time-averaged',
  list(#)%frequency:  030000,
  list(#)%duration:   030000,
EOF

     $label =   'list(#)%fields:    ';
     $n = @xTable - 1;
ONE: for $i ( 0..$n ) {

        %r = %{$xTable[$i]};
        $name  = $r{name};
        $Name = "$gcName"."::"."$name";

        if ( "$r{dims}" eq "xyz" ) { 
	    $has3D = 1;
            next ONE;
        }

        if ( $opt_m ) {
	    print FILE "  $label '$Name'     , '$gcName'      ,  '$name'   ,\n";
        } else {
	    print FILE "  $label '$name'     , '$gcName'      ,\n";
        }

    $label =   '                   '

    }

# 3D quantities
# -------------
  if ( $has3D ) {

        print FILE <<EOF;

  list(#)%filename:   '/dev/null/%s.$gcName.eta.%y4%m2%d2_%h2%n2z',
  list(#)%format:     'CFIO',
  list(#)%mode:       'time-averaged',
  list(#)%frequency:  030000,
  list(#)%duration:   030000,
EOF

     $label =   'list(#)%fields:    ';
     $n = @xTable - 1;
VAR: for $i ( 0..$n ) {

        %r = %{$xTable[$i]};
        $name  = $r{name};
        $Name = "$gcName"."::"."$name";

        next VAR if ( "$r{dims}" eq "xy" );

        if ( $opt_m ) {
	    print FILE "  $label '$Name'     , '$gcName'      ,  '$name'   ,\n";
        } else {
	    print FILE "  $label '$name'     , '$gcName'      ,\n";
        }

    $label =   '                   '

    }

 }

    close(FILE);

}

#.........................................................................

sub Preamble {

    my $c = shift;

    print FILE <<EOF;
$c                          -------------------
$c                          W  A  R  N  I  N  G
$c                          -------------------
$c
$c   This code fragment is automatically generated by $Iam.
$c   Please DO NOT edit it. Any modification made in here will be overwritten
$c   next time this file is auto-generated. Instead, enter your additions
$c   or deletions in the $gcname_Registry.rc file. 
$c

EOF

}

#.........................................................................

sub usage {

   print <<"EOF";

NAME
     maple_acg - generates code fragments for setting up MAPL states
          
SYNOPSIS

     maple_acg [OPTIONS]  Registry_filename
          
DESCRIPTION

     Reads a Grid Component Registry with definitions of INTERNAL,
     IMPORT and EXPORT states and generate the necessary code
     fragments and resource files. The following code fragments/
     resorce file are generated:

        Code Fragment Type     Default Output File Name
        -------------------    ------------------------
        Internal Spec Setup    GC_NAME_InternalSpec___.h
        Import   Spec Setup    GC_NAME_ImportSpec___.h
        Export   Spec Setup    GC_NAME_ExportSpec___.h
        GEOS_GetPointer()      GC_NAME_GetPointer___.h
        History streams        GC_NAME_History___.rc

     The Grid Component name "GC_NAME" is either specified with he -n
     option or derived from the input file name. Notice that files
     matching "*___.h" and "*___.rc" are ignored by the CVS ESMA
     repository at NASA/GSFC.

     The F90 code fragments contain the appropriate MAPL calls for
     setting up import/export/external states. In addition, it also
     generates a code fragment (*_GetPointer.h) for automatically 
     declaring and retrienving fields from these states in the form 
     of simple F90 arrays, or "arrays of arrays" - read on.

     ... to be completed ...

OPTIONS
     -C             make all arrays Chem Arrays
     -d             print variables #defines on the screen.
     -F             make all arrays flat FORTRAN arrays
     -g             file name for code fragment with calls to GEOS_GetPointer()
     -h             file name for code fragment with History.rc fragmet
     -i             file name for code fragment with Import   state spec
     -m             Mangle variable name in history stream
     -N name        component name to be used in history stream
     -n name        name of grid component; the default is derived from
                    the imput file name
     -p             file name for code fragment with Internal state spec
     -v             verbose mode
     -x             file name for code fragment with Export   state spec

SEE ALSO
     The MAPL User's Guide.

BUGS
     No parsing of the table schema at this point.
     At this point only handling of the Export state has been implemented.

AUTHOR
     Arlindo da Silva, NASA/GSFC.

EOF

  exit(1)

 }

__DATA__
#
# Sample Grid Component Registry Dust (DU) Grid Component Registry from GOCART.
# This sample file defines Export state for this component.
#
# -----------------------------------------------------------------

  COMP_NAME: DU

# Only change the Registry version when major structural changes
# occurs, not changes in content
# --------------------------------------------------------------
  GEOS_REGISTRY_VERSION: 1.00  

<DU::ExportSpec cols="SHORT UNITS DIMS VLOC STAT REFRESH AVG SUBTILES LONG">
# ---------|------------|-----|---|----|---|---|-----|---------------------------------
#  Short   |            |     | V |Item|Intervl| Sub |          Long
#  Name    |   Units    | Dim |Loc|Type| R | A |Tiles|          Name
# ---------|------------|-----|---|----|---|---|-----|---------------------------------
  DUMASS   |            | xyz | C |    |   |   |     | Dust Mass Mixing Ratio 
  DUMASS25 |            | xyz | C |    |   |   |     | Dust Mass Mixing Ratio - PM 2.5 
# .........|............|.....|...|....|...|...|.....|..................................
  DUEM001  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Emission Bin 001     
  DUEM002  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Emission Bin 002     
  DUEM003  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Emission Bin 003     
  DUEM004  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Emission Bin 004     
  DUEM005  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Emission Bin 005    
  DUSD001  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Sedimentation Bin 001     
  DUSD002  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Sedimentation Bin 002     
  DUSD003  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Sedimentation Bin 003     
  DUSD004  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Sedimentation Bin 004     
  DUSD005  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Sedimentation Bin 005    
  DUDP001  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Dry Deposition Bin 001     
  DUDP002  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Dry Deposition Bin 002     
  DUDP003  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Dry Deposition Bin 003     
  DUDP004  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Dry Deposition Bin 004     
  DUDP005  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Dry Deposition Bin 005     
  DUWT001  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Wet Deposition Bin 001     
  DUWT002  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Wet Deposition Bin 002     
  DUWT003  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Wet Deposition Bin 003     
  DUWT004  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Wet Deposition Bin 004     
  DUWT005  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Wet Deposition Bin 005     
# DUSV001  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Convective Scavenging Bin 001     
# DUSV002  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Convective Scavenging Bin 002     
# DUSV003  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Convective Scavenging Bin 003     
# DUSV004  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Convective Scavenging Bin 004     
# DUSV005  | kg m-2 s-1 | xy  |   |    |   |   |     | Dust Convective Scavenging Bin 005     
  DUSMASS  | kg m-3     | xy  |   |    |   |   |     | Dust Surface Mass Concentration    
  DUCMASS  | kg m-2     | xy  |   |    |   |   |     | Dust Column Mass Density    
  DUEXTTAU |    1       | xy  |   |    |   |   |     | Dust Extinction AOT [550 nm]    
  DUSCATAU |    1       | xy  |   |    |   |   |     | Dust Scattering AOT [550 nm]    
  DUSMASS25| kg m-3     | xy  |   |    |   |   |     | Dust Surface Mass Concentration - PM 2.5    
  DUCMASS25| kg m-2     | xy  |   |    |   |   |     | Dust Column Mass Density - PM 2.5 
  DUEXTT25 |    1       | xy  |   |    |   |   |     | Dust Extinction AOT [550 nm] - PM 2.5    
  DUSCAT25 |    1       | xy  |   |    |   |   |     | Dust Scattering AOT [550 nm] - PM 2.5   
  DUAERIDX |    1       | xy  |   |    |   |   |     | Dust TOMS UV Aerosol Index    
# ---------|------------|-----|---|----|---|---|-----|---------------------------------
</DU::ExportSpec>
