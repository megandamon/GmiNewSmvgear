! $Id: FVAdvStandalone.F90,v 1.2 2011-08-09 22:12:59 mrdamon Exp $
!
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: FVAdvStandalone.F90 - Main program source file for tracer advection 
!                                 test
!
! !DESCRIPTION:
!
! ESMF Application Wrapper for the tracer advection test.  This file contains 
! the main program, and creates a top level ESMF Gridded Component to contain
! all other Components.
!
!  
!
!
!EOP
#define ASSERT_(cond_) IF(.NOT. (cond_)) PRINT *, "assert failed in file: ", __FILE__, " at line: ",  __LINE__
    program FVAdvStandalone

    ! ESMF module, defines all ESMF data types and procedures
    use ESMF_Mod
    
    ! Devel Component registration routines
    use FVadvcore_GridCompMod, only : FVadvcore_GridCompRegister
    use FVAdvStateMod, only : FVAdvStateCreate, FVAdvStateDestroy, &
                                  FVAdvStateData

    implicit none
    
    ! Local variables

    ! Components
    type(ESMF_GridComp) :: compGridded

    ! State, Virtual Machine, and DELayout
    type(ESMF_VM) :: vm
    type(ESMF_State) :: stateIn, stateOut
    integer :: pet_id

    ! A common grid
    type(ESMF_Grid)     :: grid
    type(ESMF_DELayout) :: layoutTop

    ! A clock, a calendar, and timesteps
    type(ESMF_Clock) :: clock
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_Time) :: startTime
    type(ESMF_Time) :: stopTime

!
! Constants (should actually be declared in a constants module)
!
    real(ESMF_KIND_R8), parameter   :: PI = 3.141592654
    real(ESMF_KIND_R8)              :: xmin, xmax, ymin, ymax
    real(ESMF_KIND_R8), allocatable :: delta(:)

    ! Variables related to grid and clock
    integer :: counts(2)


    integer :: dt_seconds
    integer :: s_month, s_day, s_hour, s_min
    integer :: e_month, e_day, e_hour, e_min
    integer :: nprxy_x, nprxy_y
    integer :: im, jm, km

    integer :: imxy(2), jmxy(2)

    ! Read in from config file
    namelist /input/ dt_seconds, s_month, s_day, s_hour, s_min,       &
                     e_month, e_day, e_hour, e_min,                   &
                     nprxy_x, nprxy_y, im, jm, km

!BOP
!
! !DESCRIPTION:
! \subsubsection{Namelist Input Parameters for App:}
!     The following variables must be input to the FV Advection 
!     Application to run.  They are located in a file called 
!     "FVAdv\_input."  In addition, the tracer advection
!     component has its own resource file (see documentation in
!     that component).
!
!     The variables are:
!     \begin{description}
!     \item [dt\_seconds]
!           Delta time in seconds (integer).
!     \item [s\_month]
!           Simulation start time month (integer).
!     \item [s\_day]
!           Simulation start time day (integer).
!     \item [s\_hour]
!           Simulation start time hour (integer).
!     \item [s\_min]
!           Simulation start time minute (integer).
!     \item [e\_month]
!           Simulationendt time month (integer).
!     \item [e\_day]
!           Simulation end time day (integer).
!     \item [e\_hour]
!           Simulation end time hour (integer).
!     \item [e\_min]
!           Simulation end time minute (integer).
!     \end{description}
!
!EOP

    ! Return codes for error checks
    integer :: rc

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    ESMF_Initialize
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!

!BOP
!
! !DESCRIPTION:
! \subsubsection{Example of Initializing the Framework:}
!
!     The first call to ESMF must be the initialize method.   As part of
!     initialization the default Calendar can be specified, some options
!     for logging can be set, and the default global VM can be returned.
!     Here we are setting the default Calendar to be Gregorian, and getting
!     back the global VM:
!\begin{verbatim}
    ! Initialize ESMF, get the default Global VM, and set
    ! the default calendar to be Gregorian.
    call ESMF_Initialize(vm=vm, defaultCalendar=ESMF_CAL_GREGORIAN, rc=rc)
!\end{verbatim}
!EOP 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!   Read in configuration data - will be replaced by Config routines
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!
      !
      ! Read in input file
      !
      open(9, status="old", file="FVAdv_input",action="read",iostat=rc)
      if (rc .ne. 0) then
        print *, "Error!  Failed to open namelist file 'FVAdv_input' "
        stop
      endif
      read(9, input, end=20)
   20 continue
      close(9)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Create section
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!

!BOP
!\begin{verbatim}
    ! Create the top level Gridded Component.
    compGridded = ESMF_GridCompCreate(name="Tracer Advection Comp", rc=rc)
!\end{verbatim}
!EOP 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  Register section
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      call ESMF_GridCompSetServices(compGridded, FVadvcore_GridCompRegister, rc)
      ASSERT_(rc==0)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Initialize the grid
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    layoutTop = ESMF_DELayoutCreate( vm, deCountList=(/nprxy_x, nprxy_y/), &
                                     rc=rc)
    ASSERT_(rc==0)

!
! Get minimum and maximum coordinates of grid
!
    xmin = 0.0
    xmax = 2.0*PI
    ymin = -PI/2.0
    ymax = +PI/2.0

    grid = ESMF_GridCreateHorzLatLonUni(counts=(/im,jm/),                   &
                                    minGlobalCoordPerDim=(/xmin,ymin/),     &
                                    maxGlobalCoordPerDim=(/xmax,ymax/),     &
                                    horzStagger=ESMF_GRID_HORZ_STAGGER_A,   &
                                    periodic=(/ ESMF_TRUE, ESMF_FALSE /),   &
                                    name="A-grid", rc=rc)
    ASSERT_(rc==0)
    allocate( delta(km) )
    delta = 1.0
    call ESMF_GridAddVertHeight( grid, delta=delta,                        &
                                 vertStagger=ESMF_GRID_VERT_STAGGER_CENTER,&
                                 rc=rc )
    ASSERT_(rc==0)
    deallocate( delta )

    imxy = (/18,18/)
    jmxy = (/12,12/)
    call ESMF_GridDistribute(grid, delayout=layoutTop,                     &
!!!!                             countsPerDEDim1=imxy, countsPerDEDim2=jmxy,   &
                             decompIds=(/1,2,0/), rc=rc)
    ASSERT_(rc==0)

    call ESMF_GridCompSet(compGridded, grid=grid, rc=rc)
    ASSERT_(rc==0)


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  Create and initialize a clock, and a grid.
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
!
! !DESCRIPTION:
! \subsubsection{Example of Calendar and Clock Creation and Usage:}
!
!     The following piece of code provides an example of Clock creation used in
!     the Demo.  Note that the Gregorian calendar was set as the default in
!     the ESMF\_Initialize() call above.  As shown in this example, we first
!     initialize a time interval (timestep) to 2 seconds:
!\begin{verbatim}
      call ESMF_TimeIntervalSet(timeStep, s=dt_seconds, rc=rc)
!\end{verbatim}
!     And then we set the start time and stop time to input values for the month,
!     day, and hour (assuming the year to be 2003):
!\begin{verbatim}
      call ESMF_TimeSet(startTime, yy=2003, mm=s_month, dd=s_day, &
                        h=s_hour, m=s_min, s=0, rc=rc)

      call ESMF_TimeSet(stopTime, yy=2003, mm=e_month, dd=e_day, &
                        h=e_hour, m=e_min, s=0, rc=rc)
!\end{verbatim}
!     With the time interval, start time, and stop time set above, the Clock can
!     now be created:
!\begin{verbatim}
      clock = ESMF_ClockCreate(timeStep=timeStep, startTime=startTime, &
                               stopTime=stopTime, rc=rc)
!\end{verbatim}
!     Subsequent calls to ESMF\_ClockAdvance with this clock will increment the
!     current time from the start time by the timestep.
!EOP 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  Create and initialize a State to use for both import and export.
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      StateIn  = ESMF_StateCreate("Tracer Adv State In", rc=rc)
      ASSERT_(rc==0)
      StateOut = ESMF_StateCreate("Tracer Adv State Out", rc=rc)
      ASSERT_(rc==0)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  Initialize the Tracer Advection Component
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 
      call ESMF_GridCompInitialize(compGridded, StateIn, StateOut, &
                                   clock, rc=rc)
      ASSERT_(rc==0)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  Allocate the Tracer Advection Component states
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 
      call FVAdvStateCreate( compGridded, stateIn, stateOut )

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  Set the test data (service provided by Component)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      call FVAdvStateData( compGridded, stateIn, stateOut )

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  Run the code until the final time is reached
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      call ESMF_ClockPrint(clock, "starttime string", rc)
      do while (.not. ESMF_ClockIsStopTime(clock, rc))
        call ESMF_GridCompRun(compGridded, stateIn, stateOut, &
                              clock, rc=rc)
        call ESMF_ClockAdvance(clock, rc=rc)
      enddo
      call ESMF_ClockPrint(clock, "finaltime string", rc)

      call ESMF_GridCompFinalize(compGridded, stateIn, stateOut, &
                                 clock, rc=rc)
      ASSERT_(rc==0)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!     Destroy section
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!     Clean up

      call  FVAdvStateDestroy( compGridded, StateIn, StateOut )

      call ESMF_StateDestroy(stateIn, rc)
      ASSERT_(rc==0)
      call ESMF_StateDestroy(stateOut, rc)
      ASSERT_(rc==0)

      call ESMF_GridDestroy(grid, rc)
      ASSERT_(rc==0)

      call ESMF_ClockDestroy(clock, rc)
      ASSERT_(rc==0)

      call ESMF_GridCompDestroy(compGridded, rc)
      ASSERT_(rc==0)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      ! This output goes into the log file (standard output, unit 6)
      print *, "**********************************************************"
      print *, "SUCCESS!  Your ESMF  Application Demo ", &
               "ran to completion!"
      print *, "See the output files in the Demo source directory for ", &
               "the generated data."
      print *, "**********************************************************"

      ! Get our PET number from the VM
      call ESMF_VMGet(vm, localPET=pet_id, rc=rc)

      ! This output goes to the console/screen (standard error) where
      ! hopefully the user will see it without needing to inspect the log file.
      if (pet_id .eq. 0) then
        write(0, *) ""
        write(0, *) "SUCCESS!  Your ESMF  Application Demo ", &
                 "ran to completion!"
        write(0, *) ""
      endif

      ! Finalize must be after the message, since this call shuts down MPI if
      ! it is being used and on some platforms that prevents messages from
      ! reaching their destination files.

      call ESMF_Finalize(rc=rc)

      end program FVAdvStandalone
    
