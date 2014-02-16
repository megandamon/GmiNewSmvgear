%.o:%.F90
	rm -f $*.f90 $*.i90
  ifeq ($(F90_VENDOR), Compaq) # CPP is a separate stage
	$(CPP) $(CPPFLAGS) $*.F90
	sed '/^ *$$/d' $*.i90 > $*.f90
  else
	$(CPP) $(CPPFLAGS) $*.F90 > $*.i90
	sed '/^ *$$/d' $*.i90 > $*.f90
  endif
	rm -f $*.i90
	chmod -w $*.f90
	$(F90) -c $(FFLAGS) $*.f90 

%.o:%.c
	$(CC) $(CFLAGS) $(DEPEND_INC) $*.c
