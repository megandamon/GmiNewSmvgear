!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiESMFclock_mod
!
! !INTERFACE:
!
      module GmiESMFclock_mod
!
! !USES:
      use GmiTimeControl_mod, only : t_GmiClock, Set_leapYearFlag
      use GmiTimeControl_mod, only : Set_begGmiDate, Set_begGmiTime
      use GmiTimeControl_mod, only : Set_endGmiDate, Set_endGmiTime
      use GmiTimeControl_mod, only : Set_curGmiDate, Set_curGmiTime
      use GmiTimeControl_mod, only : Set_gmiSeconds, Set_gmiTimeStep
      use GmiTimeControl_mod, only : Set_numTimeSteps
      use GmiTimeControl_mod, only : Set_totNumDays
      use GmiPrintError_mod , only : GmiPrintError

      ! ESMF module, defines all ESMF data types and procedures
      use ESMF_Mod

      implicit none

! !PUBLIC DATA MEMBERS:
      private
      public  :: createESMFclock
      public  :: advanceESMFclock
! !DESCRIPTION:
! Routines for creating and advancing the ESMF clock using the GMI clock.
!
! !AUTHOR:
!  Jules Kouatchou, NASA/GSFC, Jules.Kouatchou@nasa.gov
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: createESMFclock
!
! !INTERFACE:
!
      subroutine createESMFclock(config, esmfClock, gmiClock, STATUS)
!
      implicit none
!
! !OUPUT PARAMETERS:
      integer         , intent(out) :: STATUS
      type(ESMF_Clock), intent(out) :: esmfClock
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Config), intent(inOut) :: config
      type(t_gmiClock ), intent(inOut) :: gmiClock
!
!
! !DESCRIPTION:
! Create the ESMF clock using information from the resource file.
!
! !LOCAL VARIABLES:
      type(ESMF_Time)         :: StartTime  ! Initial Begin  Time of Experiment
      type(ESMF_Time)         :: StopTime   ! Final   Ending Time of Experiment
      type(ESMF_TimeInterval) :: timeStep   ! HEARTBEAT
      type(ESMF_Calendar)     :: cal
      character(ESMF_MAXSTR) :: IAM="createClock"
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
      integer        :: leapYearFlag, num_time_steps, ndt

      integer        :: begGmiDate, begGmiTime, endGmiDate, endGmiTime
      real*8         :: gmi_sec, tdt
      real           :: totSec
      character(ESMF_MAXSTR) :: err_msg

      integer            :: totNumDays
      real(ESMF_KIND_R8) :: runTimeStepCount
      INTEGER, PARAMETER :: SEC_PER_DAY = 24*60*60

      ! This variable is used to force all the years to be leap years.
      ! It is used when leapYearFlag is greater than zero.
      integer, parameter :: nDaysPerMonth(12) = &
     &                        (/ 31, 29, 31, 30, 31, 30, &
     &                           31, 31, 30, 31, 30, 31  /)
!EOP
!------------------------------------------------------------------------------
!BOC
! Read Times From Config
! ----------------------

! !RESOURCE_ITEM: year :: Beginning year (integer)
      call ESMF_ConfigGetAttribute(config, BEG_YY, label='BEG_YY:',            &
     &                default=1989, rc=STATUS )
!      VERIFY_(STATUS)

! !RESOURCE_ITEM: month :: Beginning month (integer 1-12)
      call ESMF_ConfigGetAttribute(config, BEG_MM, label='BEG_MM:',            &
     &                default=1, rc=STATUS )
!      VERIFY_(STATUS)

! !RESOURCE_ITEM: day  :: Beginning day of month (integer 1-31)
      call ESMF_ConfigGetAttribute(config, BEG_DD, label='BEG_DD:', default=1, &
     &                rc=STATUS )
!      VERIFY_(STATUS)

! !RESOURCE_ITEM: hour :: Beginning hour of day (integer 0-23)
      call ESMF_ConfigGetAttribute(config, BEG_H, label='BEG_H:', default=0,   &
     &                rc=STATUS )
!      VERIFY_(STATUS)

! !RESOURCE_ITEM: minute :: Beginning minute (integer 0-59)
      call ESMF_ConfigGetAttribute(config, BEG_M, label='BEG_M:', default=0,   &
     &                rc=STATUS )
!      VERIFY_(STATUS)

! !RESOURCE_ITEM: second :: Beginning second (integer 0-59)
      call ESMF_ConfigGetAttribute(config, BEG_S, label='BEG_S:', default=0,   &
     &                rc=STATUS )
!      VERIFY_(STATUS)

! !RESOURCE_ITEM: year :: Ending year (integer)
      call ESMF_ConfigGetAttribute(config, END_YY, label='END_YY:',            &
     &                default=1989, rc=STATUS )
!      VERIFY_(STATUS)

! !RESOURCE_ITEM: month :: Ending month (integer 1-12)
      call ESMF_ConfigGetAttribute(config, END_MM, label='END_MM:', default=1, &
     &                rc=STATUS )
!      VERIFY_(STATUS)

! !RESOURCE_ITEM: day  :: Ending day of month (integer 1-31)
      call ESMF_ConfigGetAttribute(config, END_DD, label='END_DD:', default=1, &
     &                rc=STATUS )
!      VERIFY_(STATUS)

! !RESOURCE_ITEM: hour :: Ending hour of day (integer 0-23)
      call ESMF_ConfigGetAttribute(config, END_H, label='END_H:', default=0,   &
     &                rc=STATUS )
!      VERIFY_(STATUS)

! !RESOURCE_ITEM: minute :: Ending minute (integer 0-59)
      call ESMF_ConfigGetAttribute(config, END_M, label='END_M:', default=0,   &
     &                rc=STATUS )
!      VERIFY_(STATUS)

! !RESOURCE_ITEM: second :: Ending second (integer 0-59)
      call ESMF_ConfigGetAttribute(config, END_S, label='END_S:', default=0,   &
     &                rc=STATUS )
!      VERIFY_(STATUS)

! !RESOURCE_ITEM: seconds :: Interval of the application clock (the Heartbeat)
      call ESMF_ConfigGetAttribute(config, RUN_DT, label='RUN_DT:', &
     &                default = 180, rc=STATUS )
!      VERIFY_(STATUS

! !RESOURCE_ITEM: integer flag to determine the type of calendar to use
      call ESMF_ConfigGetAttribute(config, leapYearFlag, &
                      label='leapYearFlag:', default=0, rc=STATUS )
!      VERIFY_(STATUS)

! initialize calendar
! -------------------

      ! Determine the type of calendar to be used
      ! If leapYearFlag < 0 then no year is a leap year
      !                 = 0 then we have the Gregorian calendar
      !                 > 0 then all years are leap years

      if (leapYearFlag < 0) then
         cal = ESMF_CalendarCreate("noLeapYear", ESMF_CAL_NOLEAP, rc=STATUS)
      elseif (leapYearFlag == 0) then
         cal = ESMF_CalendarCreate("Gregorian", ESMF_CAL_GREGORIAN, rc=STATUS)
      else
         cal = ESMF_CalendarCreate("allLeap", daysPerMonth=nDaysPerMonth,      &
     &                                     rc=STATUS)
      end if
!      VERIFY_(STATUS)

! initialize start time for Alarm frequencies
! -------------------------------------------

      call ESMF_TimeSet(StartTime, YY = BEG_YY, MM = BEG_MM, DD = BEG_DD,      &
     &         H = BEG_H, M = BEG_M, S = BEG_S, calendar=cal, rc = STATUS)
!     VERIFY_(STATUS)

! initialize final stop time
! --------------------------

      call ESMF_TimeSet(StopTime, YY = END_YY, MM = END_MM, DD = END_DD,       &
     &         H = END_H, M = END_M, S = END_S, calendar=cal, rc = STATUS)
!     VERIFY_(STATUS)


! initialize model time step
!    N.B. CONVERTED TO INTEGER!!
! ------------------------------

     call ESMF_TimeIntervalSet( timeStep, S=RUN_DT, rc=STATUS )
!     VERIFY_(STATUS)

     ! Create the ESMF clock

      esmfClock = ESMF_ClockCreate(timeStep=timeStep, startTime=startTime, &
     &                             stopTime=stopTime, rc=STATUS)
!     VERIFY_(STATUS)


      ! Determine the total number of time steps
      call ESMF_ClockGet(esmfClock, runTimeStepCount=runTimeStepCount)

      ! Determine the length of the experiment in days
      totNumDays = runTimeStepCount * RUN_DT / SEC_PER_DAY

      !#####################
      ! Create the GMI Clock
      !#####################

       call ESMF_ConfigGetAttribute(config, totSec, label='gmi_sec:', &
     &                default=0.0, rc=STATUS )
!     VERIFY_(STATUS)

      gmi_sec        = 1.0d0 * totSec
      tdt            = 1.0d0 * RUN_DT
      ndt            = RUN_DT
      begGmiDate     = 10000 * BEG_YY + 100 * BEG_MM + BEG_DD
      begGmiTime     = 10000 * BEG_H  + 100 * BEG_M  + BEG_S
      endGmiDate     = 10000 * END_YY + 100 * END_MM + END_DD
      endGmiTime     = 10000 * END_H  + 100 * END_M  + END_S
      num_time_steps = Nint (gmi_sec) / ndt
      
      
      call Set_leapYearFlag (leapYearFlag)
      call Set_gmiTimeStep (gmiClock, tdt           )
      call Set_curGmiDate  (gmiClock, begGmiDate    )
      call Set_curGmiTime  (gmiClock, begGmiTime    )
      call Set_gmiSeconds  (gmiClock, gmi_sec       )
      call Set_begGmiDate  (gmiClock, begGmiDate    )
      call Set_begGmiTime  (gmiClock, begGmiTime    )
      call Set_endGmiDate  (gmiClock, endGmiDate    )
      call Set_endGmiTime  (gmiClock, endGmiTime    )
      call Set_totNumDays  (gmiClock, totNumDays    )
      call Set_numTimeSteps(gmiClock, num_time_steps)

      if (tdt <= 0.0d0) then     
         err_msg = 'RUN_DT <= 0.0 in the resource file.'
         call GmiPrintError (err_msg, .true., 0, 0, 0, 1, tdt, 0.0d0)
      end if
      
      if (Mod (Nint (gmi_sec), ndt) /= 0) then
         err_msg = 'gmi_sec is not evenly div. by ndt.'
         call GmiPrintError (err_msg, .true., 1, ndt, 0, 1, gmi_sec, 0.0d0)
      end if

      return

      end subroutine createESMFclock
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine advanceESMFclock(esmfClock, gmiClock)
!
! !USES
      use GmiTimeControl_mod, only : updateGmiClock

      implicit none
!
! !INPUT/OUPUT PARAMETERS:
      type (t_GmiClock), intent(inOut) :: gmiClock
      type(ESMF_Clock ), intent(inOut) :: esmfClock
!
! !DESCRIPTION:
! Advances the ESMF clock and updates the GMI clock.
!     
! !LOCAL VARIABLES:
      integer               :: rc
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Advance the ESMF clock
      call ESMF_ClockAdvance(esmfClock, rc=rc)
 
      ! Update the GMI clock
      call updateGmiClock (gmiClock, esmfClock)

      return
      end subroutine advanceESMFclock
!EOC
!------------------------------------------------------------------------------

      end module GmiESMFclock_mod
