!  $Id: MAPL_Cap.F90,v 1.2 2011-08-09 22:13:00 mrdamon Exp $

#include "MAPL_Generic.h"

module MAPL_CapMod

  use ESMF_Mod
  use MAPL_BaseMod
  use MAPL_ConstantsMod
  use MAPL_GenericMod
  use MAPL_HistoryGridCompMod, only : Hist_SetServices => SetServices
  implicit none

  private
  public MAPL_Cap

contains
!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine MAPL_CAP(SS, RC)
    external :: SS
    integer, optional, intent(out) :: rc


! Handles to the CAP's Gridded Components GCs
! -------------------------------------------

   integer                      :: ROOT
   integer                      :: HIST
   character(len=ESMF_MAXSTR)   :: ROOT_NAME

! A MAPL object for the cap
!--------------------------

   type(MAPL_MetaComp)          :: META

! The children's GCs and IM/Ex states
!------------------------------------

   type(ESMF_GridComp), pointer :: GCS(:)
   type(ESMF_State),    pointer :: IMPORTS(:)
   type(ESMF_State),    pointer :: EXPORTS(:)

! ESMF stuff
!-----------

   type(ESMF_VM)                :: VM
   type(ESMF_Config)            :: config
   type(ESMF_Clock)             :: clock
   type(ESMF_Grid)              :: grid

! ErrLog variables
!-----------------

   integer                      :: STATUS
   character(len=ESMF_MAXSTR), parameter :: Iam="MAPL_Cap"

! Misc locals
!------------

   character(len=ESMF_MAXSTR)   :: cf_file
   character(len=ESMF_MAXSTR)   :: cf_label
   integer                      :: tick_first
   integer                      :: n
   logical                      :: tick_at_top

! Begin
!------

!  Initialize ESMF
!-----------------

   call ESMF_Initialize (vm=vm, rc=status)
   VERIFY_(STATUS)

!  Open the CAP's configuration from CAP.rc
!------------------------------------------

   config = ESMF_ConfigCreate (                   rc=STATUS )
   VERIFY_(STATUS)
   call ESMF_ConfigLoadFile   ( config, 'CAP.rc', rc=STATUS )
   VERIFY_(STATUS)

!  CAP's MAPL MetaComp
!---------------------
!BOR
   call MAPL_Set (META, name='CAP', cf=CONFIG,    rc=STATUS )
   VERIFY_(STATUS)

!  Create Clock
!--------------

   call MAPL_ClockInit ( META, clock,             rc=STATUS )
   VERIFY_(STATUS)

!  Create grid
! -----------

   GRID = MAPL_LogRectGridCreate( META,            rc=STATUS ) 
   VERIFY_(STATUS)

!  Create Root child
!-------------------

! !RESOURCE_ITEM: none :: Name of ROOT's config file 
   call MAPL_GetResource(META, CF_FILE,"ROOT_RC_FILE:",default="ROOT.rc", &
                                                  rc=STATUS ) 
   VERIFY_(STATUS)

   call MAPL_GetResource(META, ROOT_NAME,"MAPLROOT_COMPNAME:",default="MAPLROOT", &
                                                  rc=STATUS ) 
   VERIFY_(STATUS)

   ROOT = MAPL_AddChild ( META,        &
        name       = ROOT_NAME,        &
        grid       = GRID,             &
        configfile = CF_FILE,          &
        SS         = SS,               &
                             rc=STATUS )  
   VERIFY_(STATUS)

!  Create History child
!----------------------

! !RESOURCE_ITEM: none :: Name of HISTORY's config file 
   call MAPL_GetResource(META, CF_FILE,"HIST_CONFIG:", default='HISTORY.rc',&
                         RC=STATUS ) 
   VERIFY_(STATUS)

   HIST = MAPL_AddChild ( META,        &
        name       = 'HIST',           &
        grid       = grid,             &
        configfile = cf_file,          &
        SS         = HIST_SetServices, &
                             rc=STATUS )  
   VERIFY_(STATUS)

!  Query MAPL for the the children's for GCS, IMPORTS, EXPORTS
!-------------------------------------------------------------

   call MAPL_Get(META, GCS=GCS, GIM=IMPORTS, GEX=EXPORTS, RC=STATUS)
   VERIFY_(STATUS)

!  Initialize the Computational Hierarchy
!----------------------------------------

   call ESMF_GridCompInitialize ( GCS(ROOT), IMPORTS(ROOT), EXPORTS(ROOT), CLOCK, rc=STATUS )
   VERIFY_(STATUS)

! All the EXPORTS of the Hierachy are made IMPORTS of History
!------------------------------------------------------------

   call ESMF_StateAddState(IMPORTS(HIST), EXPORTS(ROOT), rc=status)
   VERIFY_(STATUS)

! this probably is a "bug" in HISTORY but it is required :(
   call ESMF_StateAddState(IMPORTS(HIST), EXPORTS(HIST), rc=status)
   VERIFY_(STATUS)

! The history component takes as IMPORT the "all" EXPORTS (and has "fake" EXPORT)

   call ESMF_GridCompInitialize ( GCS(HIST), IMPORTS(HIST), EXPORTS(HIST), clock, &
                                  rc=status )
   VERIFY_(STATUS)
 
! Advance the clock 1 timestep to set proper "ringing" state for all alarms
!--------------------------------------------------------------------------

   call ESMF_ClockAdvance ( clock, rc=status )
   VERIFY_(STATUS)


! *********************************************************************
! *****                     Main Time Loop                         ****
! *********************************************************************

! !RESOURCE_ITEM: 1 or 0 :: Determines when clock is advanced 
   call MAPL_GetResource(META, TICK_FIRST, "TICK_FIRST:",  &
                         default=0,              RC=STATUS ) 
   VERIFY_(STATUS)

   TICK_AT_TOP = tick_first /= 0

  do

! Advance the Clock
! -----------------

     if (TICK_AT_TOP) then
        call ESMF_ClockAdvance ( clock, rc=status )
        VERIFY_(STATUS)
     end if

! Run the Gridded Component
! --------------------------

     call ESMF_GridCompRun ( GCS(ROOT), IMPORTS(ROOT), EXPORTS(ROOT), clock, rc=status )
     VERIFY_(STATUS)

! Advance the Clock
! -----------------

     if (.not.TICK_AT_TOP) then
        call ESMF_ClockAdvance ( clock, rc=status )
        VERIFY_(STATUS)
     end if


! Call History Run for Output
! ---------------------------

     call ESMF_GridCompRun ( GCS(HIST), IMPORTS(HIST), EXPORTS(HIST), clock, rc=status )
     VERIFY_(STATUS)

! Call Record for intermediate checkpoint (if desired)
! ----------------------------------------------------
      call ESMF_GridCompFinalize( GCS(ROOT),IMPORTS(ROOT),EXPORTS(ROOT),clock, &
                                  phase=MAPL_RecordPhase, rc=status )
      VERIFY_(STATUS)
! Check for Segment Ending Time
! -----------------------------

      if ( ESMF_ClockIsStopTime( clock, rc=status ) ) exit
      VERIFY_(STATUS)

! Synchronize for Next TimeStep
! -----------------------------

      call ESMF_VMBarrier(VM, rc=status)
      VERIFY_(STATUS)

  enddo ! end of time loop

!  Finalize
!  --------

   do n=1,size(GCS)
      call ESMF_GridCompFinalize( GCS(n),IMPORTS(n),EXPORTS(n),clock,rc=status )
      VERIFY_(STATUS)
   enddo

!  Finalize framework
!  ------------------

   call ESMF_Finalize (RC=status)
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)

 END subroutine MAPL_CAP


  subroutine MAPL_ClockInit ( META, Clock, rc)

     type(MAPL_MetaComp), intent(inout) :: META
     type(ESMF_Clock),    intent(  out) :: Clock
     integer, optional,   intent(  out) :: rc

     type(ESMF_Time)          :: StartTime    ! Initial     Begin  Time of Experiment
     type(ESMF_Time)          ::  StopTime    ! Final       Ending Time of Experiment
     type(ESMF_TimeInterval)  ::  timeStep    ! HEARTBEAT
     type(ESMF_Calendar)      ::  cal

     integer        :: STATUS
     character(ESMF_MAXSTR) :: IAM="MAPL_ClockInit"

     integer        :: BEG_YY
     integer        :: BEG_MM
     integer        :: BEG_DD
     integer        :: BEG_H
     integer        :: BEG_M
     integer        :: BEG_S

     integer        :: END_YY
     integer        :: END_MM
     integer        :: END_DD
     integer        :: END_H
     integer        :: END_M
     integer        :: END_S

     integer        :: RUN_DT

! Begin
!------

! Read Times From Config
! ----------------------

! !RESOURCE_ITEM: year :: Beginning year (integer)
     call MAPL_GetResource( META, BEG_YY, label='BEG_YY:', DEFAULT=1, rc=STATUS )
     VERIFY_(STATUS)
! !RESOURCE_ITEM: month :: Beginning month (integer 1-12)
     call MAPL_GetResource( META, BEG_MM, label='BEG_MM:', default=1, rc=STATUS )
     VERIFY_(STATUS)
! !RESOURCE_ITEM: day  :: Beginning day of month (integer 1-31)
     call MAPL_GetResource( META, BEG_DD, label='BEG_DD:', default=1, rc=STATUS )
     VERIFY_(STATUS)
! !RESOURCE_ITEM: hour :: Beginning hour of day (integer 0-23)
     call MAPL_GetResource( META, BEG_H , label='BEG_H:' , default=0, rc=STATUS )
     VERIFY_(STATUS)
! !RESOURCE_ITEM: minute :: Beginning minute (integer 0-59)
     call MAPL_GetResource( META, BEG_M , label='BEG_M:' , default=0, rc=STATUS )
     VERIFY_(STATUS)
! !RESOURCE_ITEM: second :: Beginning second (integer 0-59)
     call MAPL_GetResource( META, BEG_S , label='BEG_S:' , default=0, rc=STATUS )
     VERIFY_(STATUS)

! !RESOURCE_ITEM: year :: Ending year (integer)
     call MAPL_GetResource( META, END_YY, label='END_YY:', DEFAULT=1, rc=STATUS )
     VERIFY_(STATUS)
! !RESOURCE_ITEM: month :: Ending month (integer 1-12)
     call MAPL_GetResource( META, END_MM, label='END_MM:', default=1, rc=STATUS )
     VERIFY_(STATUS)
! !RESOURCE_ITEM: day  :: Ending day of month (integer 1-31)
     call MAPL_GetResource( META, END_DD, label='END_DD:', default=1, rc=STATUS )
     VERIFY_(STATUS)
! !RESOURCE_ITEM: hour :: Ending hour of day (integer 0-23)
     call MAPL_GetResource( META, END_H , label='END_H:' , default=0, rc=STATUS )
     VERIFY_(STATUS)
! !RESOURCE_ITEM: minute :: Ending minute (integer 0-59)
     call MAPL_GetResource( META, END_M , label='END_M:' , default=0, rc=STATUS )
     VERIFY_(STATUS)
! !RESOURCE_ITEM: second :: Ending second (integer 0-59)
     call MAPL_GetResource( META, END_S , label='END_S:' , default=0, rc=STATUS )
     VERIFY_(STATUS)

! !RESOURCE_ITEM: seconds :: Interval of the application clock (the Heartbeat)
     call MAPL_GetResource( META, RUN_DT, label='RUN_DT:',            rc=STATUS )
     VERIFY_(STATUS)

! initialize calendar to be Gregorian type
! ----------------------------------------

     cal = ESMF_CalendarCreate( "GregorianCalendar",ESMF_CAL_GREGORIAN,rc=status )
     VERIFY_(STATUS)

! initialize start time for Alarm frequencies
! -------------------------------------------

     call ESMF_TimeSet( StartTime, YY = BEG_YY, &
                                   MM = BEG_MM, &
                                   DD = BEG_DD, &
                                    H = BEG_H , &
                                    M = BEG_M , &
                                    S = BEG_S , &
                    calendar=cal,  rc = STATUS  )
     VERIFY_(STATUS)

! initialize final stop time
! --------------------------

     call ESMF_TimeSet(  StopTime, YY = END_YY, &
                                   MM = END_MM, &
                                   DD = END_DD, &
                                    H = END_H , &
                                    M = END_M , &
                                    S = END_S , &
                    calendar=cal,  rc = STATUS  )
     VERIFY_(STATUS)


! initialize model time step
!    N.B. CONVERTED TO INTEGER!!
! ------------------------------

     call ESMF_TimeIntervalSet( timeStep, S=RUN_DT, rc=STATUS )
     VERIFY_(STATUS)

! Create Clock and set it to one time step before StartTime.
! After Initialize has created all alarms, we will advance the
! clock to ensure the proper ringing state of all alarms
!-------------------------------------------------------------

     StartTime = StartTime - TimeStep

     clock = ESMF_ClockCreate( "ApplClock",timeStep,StartTime,StopTime, rc=STATUS )
     VERIFY_(STATUS)

     call ESMF_ClockSet ( clock, CurrTime=StartTime, rc=STATUS )
     VERIFY_(STATUS)

     RETURN_(ESMF_SUCCESS)
   end subroutine MAPL_ClockInit



!EOR
!BOP

! !IROUTINE: MAPL_LogRectGridCreate

! !DESCRIPTION:  Creates certain standard logically-rectangular 3-D grids.
!      The nature of the grid can be controlled from the configuration.
!      In all cases the horizontal is represented by the first two
!      dimensions. Uniform Lat-Lon grids that spans the sphere are created
!      directly, relying on the configuration for grid parameters.
!      Other, more complicated grids require that a special grid create
!      function (AppGridCreate) with the same signature as MAPL\_LogRectGridCreate
!      be suppied by the user. A stub of this function lives in MAPL.
!      
!      The grid is created in the current M with a layout obtained from the 
!      configuration. If the number of grid points per PET is not specified
!      in the configuration a default distribution is used.
!
! !INTERFACE:

  function MAPL_LogRectGridCreate (META, RC) result(GRID)

! !ARGUMENTS:

    type(MAPL_MetaComp), intent(INOUT) :: META
    integer, optional,   intent(  OUT) :: rc
    type (ESMF_Grid)                   :: grid

!EOP

! Local vars

    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: IAm='MAPL_LogRectGridCreate'

    integer                         :: IM_WORLD
    integer                         :: JM_WORLD
    integer                         :: LM
    integer                         :: L
    integer                         :: NX, NY
    integer                         :: POLEEDGE
    integer, allocatable            :: IMS(:), JMS(:)
    character(len=ESMF_MAXSTR)      :: gridname
    real(ESMF_KIND_R8)              :: minCoord(3)
    real(ESMF_KIND_R8)              :: deltaX, deltaY, deltaZ
    real                            :: LON0
    type (ESMF_VM)                  :: VM
    type(ESMF_DELayout)             :: layout
    integer                         :: latlon

!BOR

! The default is a lat-lon grid, which MAPL can create. If it is not
! lat-lon, a special funcion must be provided in the executable to create
! the APP grid. (Later this can be extended to exchange grid creation,
! leaving the special function for more exotic cases.)
!------------------------------------------------------------------------

! !RESOURCE_ITEM: 0 or 1 :: 1 -> regular lat-lon; 0 -> custom grid
    call MAPL_GetResource( META, LATLON, 'LATLON:', default=1, rc = status )
    VERIFY_(STATUS)

    if(LATLON /= 0) then

! We need the VM to create the grid ??
!------------------------------------

       call ESMF_VMGetCurrent(vm, rc=STATUS)
       VERIFY_(STATUS)

! Get Decomposition from CF
!--------------------------

! !RESOURCE_ITEM: none :: Processing elements in 1st dimension
       call MAPL_GetResource( META, NX,       label ='NX:', default=1, rc = status )
       VERIFY_(STATUS)
! !RESOURCE_ITEM: none :: Processing elements in 2nd dimension
       call MAPL_GetResource( META, NY,       label ='NY:', default=1, rc = status )
       VERIFY_(STATUS)

! Get World problem size from CF
!-------------------------------

! !RESOURCE_ITEM: none :: Grid size in 1st dimension
       call MAPL_GetResource( META, IM_WORLD, 'IM:',            rc = status )
       VERIFY_(STATUS)
! !RESOURCE_ITEM: none :: Grid size in 2nd dimension
       call MAPL_GetResource( META, JM_WORLD, 'JM:',            rc = status )
       VERIFY_(STATUS)
! !RESOURCE_ITEM: none :: Grid size in 3rd dimension
       call MAPL_GetResource( META, LM,       'LM:', default=1, rc = status )
       VERIFY_(STATUS)

! The grid's name is optional
!----------------------------

! !RESOURCE_ITEM: none :: Optional grid name
       call MAPL_GetResource( META, GRIDNAME, 'GRIDNAME:', default='APPGRID', rc = status )
       VERIFY_(STATUS)

! Give the IMS and JMS the MAPL default distribution
! --------------------------------------------------

       allocate( IMS(0:NX-1) )  
       allocate( JMS(0:NY-1) )  

       call MAPL_DecomposeDim ( IM_WORLD, IMS, NX )
       call MAPL_DecomposeDim ( JM_WORLD, JMS, NY )

! Override them with alues in CF, if any
!---------------------------------------

! !RESOURCE_ITEM: none :: gridpoints in each PE along 1st dimension
       call MAPL_GetResource( META, IMS, 'IMS:', default=IMS, rc = status )
       VERIFY_(STATUS)
! !RESOURCE_ITEM: none :: gridpoints in each PE along 2nd dimension
       call MAPL_GetResource( META, JMS, 'JMS:', default=JMS, rc = status )
       VERIFY_(STATUS)

! Lat-Lon grids cover the sphere, but the first grid box can have the pole on the
!  center or on the edge. The latter is the FV way, in which there is an "extra"
!  grid box in the meridional direction. This is the default.
!--------------------------------------------------------------------------------

! !RESOURCE_ITEM: 0 or 1 :: 1->gridedge at pole; 0->gridpoint at pole
       call MAPL_GetResource( META, POLEEDGE, 'POLEEDGE:', default=0       , rc = status )
       VERIFY_(STATUS)
! !RESOURCE_ITEM: degrees :: Longituce of center of first gridbox
       call MAPL_GetResource( META,     LON0, 'BEGLON:'  , default=-90., rc = status )
       VERIFY_(STATUS)

       LON0 = LON0 * (MAPL_PI/180.)
!EOR

! Lat-Lon Grid definition
!-------------------------
       
       deltaX      = 2.0*MAPL_PI/IM_WORLD
       minCoord(1) = LON0-deltaX/2 

       if(POLEEDGE==0) then
          deltaY = MAPL_PI/(JM_WORLD-1)
          minCoord(2) = -MAPL_PI/2-deltaY/2
       else
          deltaY = MAPL_PI/(JM_WORLD  )
          minCoord(2) = -MAPL_PI/2
       end if

       deltaZ = 1.0D0
       minCoord(3) = deltaZ/2

! We should have a much simpler create with the next ESMF grid design
!--------------------------------------------------------------------

       layout = ESMF_DELayoutCreate(vm, deCountList=(/NX, NY/), rc=status)
       VERIFY_(STATUS)

       grid = ESMF_GridCreateHorzLatLonUni(         &
            counts = (/IM_WORLD, JM_WORLD/),        &
            minGlobalCoordPerDim=minCoord(1:2),     &
            deltaPerDim=(/deltaX, deltaY /),        &
            horzStagger=ESMF_Grid_Horz_Stagger_A,   &
            periodic=(/ESMF_TRUE, ESMF_FALSE/),     &
            name=gridname,                          &
            rc=status)
       VERIFY_(STATUS)

       call ESMF_GridAddVertHeight(grid,            &
            delta=(/(deltaZ, L=1,LM) /),            &
            vertStagger=ESMF_GRID_VERT_STAGGER_TOP, &
            rc=status)
       VERIFY_(STATUS)

       call ESMF_GridDistribute(grid,               &
            deLayout=layout,                        &
            countsPerDEDim1=ims,                    &
            countsPerDEDim2=jms,                    &
            rc=status)
       VERIFY_(STATUS)

    else

! Can we have a stubb for appgridcreate in mapl ?
       RETURN_(ESMF_FAILURE)
!       Grid = AppGridCreate(META, rc=STATUS)
!       VERIFY_(STATUS)

    endif

! All Done
!---------

    deallocate(ims)
    deallocate(jms)

    RETURN_(STATUS)
  end function MAPL_LogRectGridCreate

end module MAPL_CapMod

! stub functions here
 
