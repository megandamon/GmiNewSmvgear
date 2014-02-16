#
# $Id: Makefile.init.sp,v 1.2 2011-08-08 20:17:08 mrdamon Exp $
#
# Gem -- A climate systems model
# Lawrence Livermore National Laboratory
# For a successful make procedure,
# please set the following UNIX environment variables
# GEMHOME HOSTMACH
# See the file gem/doc/README.environ for further instructions
#
# Note that the difference between this file, Makefile.init.sp, &
# Makefile.init is the "-w" in the "cc" line below.  Without the
# -w, the code was producing a repititious warning message that
# was of no consequence.  The -w disables all warning messages.
#


Makefile:
	cp Makefile.cpp Makefile.c
	cc -w -E -D${HOSTMACH} -DMAKING_MAKEFILE -I${GEMHOME}/Config Makefile.c | cat -s > Makefile
	rm -f Makefile.c
	sleep 1
	touch Makefile.cpp
	make -s Makefile

