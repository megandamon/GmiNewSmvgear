!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiTimeControl_mod
!
! !INTERFACE:
!
  module GmiTimeControl_mod
!
! !USES:
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  ::  Get_leapYearFlag, Set_leapYearFlag, initializeDate
  public  ::  GmiAdvanceClock, updateGmiClock, incrementDate
  public  ::  GetSecondsFromJanuary1
  public  ::  GetDaysFromJanuary1
  public  ::  ConvertTimeToSeconds
  public  ::  ConvertSecondstoTime
  public  ::  ConvertDateToSeconds
  public  ::  GmiSplitDateTime
!
  public  :: Get_begGmiDate, Get_begGmiTime
  public  :: Get_endGmiDate, Get_endGmiTime
  public  :: Get_curGmiDate, Get_curGmiTime
  public  :: Get_totNumDays, Set_totNumDays
  public  :: Get_gmiTimeStep
  public  :: Get_gmiSeconds
  public  :: Get_numTimeSteps, Set_numTimeSteps
  public  :: Get_totNumTimeSteps, Set_totNumTimeSteps
  public  :: Set_begGmiDate, Set_begGmiTime
  public  :: Set_endGmiDate, Set_endGmiTime
  public  :: Set_curGmiDate, Set_curGmiTime
  public  :: Set_gmiTimeStep
  public  :: Set_gmiSeconds
  public  :: isOutTime, isOutFreqTime
!
! !PUBLIC DATA MEMBERS:
!
  public  :: SIZTUSML
  public  :: leapYearFlag
  public  :: t_GmiClock

# include "gmi_time_constants.h"

  type t_GmiClock
     integer  :: begGmiDate    ! beginning date of the experiment year/month/day (YYYYMMDD)
     integer  :: begGmiTime    ! beginning time of the experiment hour/min/sec   (HHMMSS)
     integer  :: endGmiDate    ! end       date of the experiment year/month/day (YYYYMMDD)
     integer  :: endGmiTime    ! end       time of the experiment hour/min/sec   (HHMMSS)
     real*8   :: gmiTimeStep   ! GMI model time step (s)
     integer  :: curGmiDate    ! current   date of the experiment year/month/day (YYYYMMDD)
                               ! This is updated at each iteration
     integer  :: curGmiTime    ! current   time of the experiment hour/min/sec   (HHMMSS)
                               ! This is updated at each iteration
     real*8   :: gmiSeconds    ! current integration time in seconds
     integer  :: numTimeSteps  ! current number of time steps
     integer  :: totNumTimeSteps  ! total number of time steps
     integer  :: totNumDays       ! total number of days
  end type t_GmiClock

  integer, parameter :: SIZTUSML = 16
  integer, save      :: leapYearFlag ! Used to determine how leap years are derived
  integer, parameter :: refYear = 0 ! reference year used to date/time init.
  integer, parameter :: refDay  = 1 ! reference day  used to date/time init.
  real*8 , parameter :: secsPerMin      = 60.
  real*8 , parameter :: secsPerHour     = secsPerMin *60.
  real*8 , parameter :: secsPerDay      = secsPerHour*24.
  real*8 , parameter :: secsPerWeek     = secsPerDay * 7.
  real*8 , parameter :: secsp400yrs     = secsPerDay * (365 * 400 + 24 * 4 + 1)
  real*8 , parameter :: secspnrmlyr     = secsPerDay * 365
  real*8 , parameter :: secspleapyr     = secsPerDay * (365 +1)
  real*8 , parameter :: secspnrmlcentry = secsPerDay * (365 * 100 + 24)
  real*8 , parameter :: secspnrml4yrblk = secsPerDay * (365 * 4 + 1)
  character (len=9), data :: monthsYear(12) = &
                           (/'JANUARY  ', 'FEBRUARY ', 'MARCH    ', &
                             'APRIL    ', 'MAY      ', 'JUNE     ', &
                             'JULY     ', 'AUGUST   ', 'SEPTEMBER', &
                             'OCTOBER  ', 'NOVEMBER ', 'DECEMBER ' /)
!
! !DESCRIPTION:
!
! !AUTHOR:
!   John Tannahill (jrt@llnl.gov) and Jules Kouatchou (kouatchou@gsfc.nasa.gov)
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!
      CONTAINS
!
!-------------------------------------------------------------------------
      subroutine Get_begGmiDate (self, ymd)
          integer         , intent(out) :: ymd
          type(t_GmiClock), intent(in ) :: self
          ymd = self%begGmiDate
          return
      end subroutine Get_begGmiDate
!-------------------------------------------------------------------------
      subroutine Set_begGmiDate (self, ymd)
          integer         , intent(in   ) :: ymd
          type(t_GmiClock), intent(inout) :: self
           self%begGmiDate = ymd
          return
      end subroutine Set_begGmiDate
!-------------------------------------------------------------------------
      subroutine Get_begGmiTime (self, hms)
          integer         , intent(out) :: hms
          type(t_GmiClock), intent(in ) :: self
          hms = self%begGmiTime
          return
      end subroutine Get_begGmiTime
!-------------------------------------------------------------------------
      subroutine Set_begGmiTime (self, hms)
          integer         , intent(in   ) :: hms
          type(t_GmiClock), intent(inout) :: self
           self%begGmiTime = hms
          return
      end subroutine Set_begGmiTime
!-------------------------------------------------------------------------
      subroutine Get_endGmiDate (self, ymd)
          integer         , intent(out) :: ymd
          type(t_GmiClock), intent(in ) :: self
          ymd = self%endGmiDate
          return
      end subroutine Get_endGmiDate
!-------------------------------------------------------------------------
      subroutine Set_endGmiDate (self, ymd)
          integer         , intent(in   ) :: ymd
          type(t_GmiClock), intent(inout) :: self
           self%endGmiDate = ymd
          return
      end subroutine Set_endGmiDate
!-------------------------------------------------------------------------
      subroutine Get_endGmiTime (self, hms)
          integer         , intent(out) :: hms
          type(t_GmiClock), intent(in ) :: self
          hms = self%endGmiTime
          return
      end subroutine Get_endGmiTime
!-------------------------------------------------------------------------
      subroutine Set_endGmiTime (self, hms)
          integer         , intent(in   ) :: hms
          type(t_GmiClock), intent(inout) :: self
           self%endGmiTime = hms
          return
      end subroutine Set_endGmiTime
!-------------------------------------------------------------------------
      subroutine Get_totNumDays (self, totNumDays)
          integer         , intent(out) :: totNumDays
          type(t_GmiClock), intent(in ) :: self
          totNumDays = self%totNumDays
          return
      end subroutine Get_totNumDays
!-------------------------------------------------------------------------
      subroutine Set_totNumDays (self, totNumDays)
          integer         , intent(in   ) :: totNumDays
          type(t_GmiClock), intent(inout) :: self
           self%totNumDays = totNumDays
          return
      end subroutine Set_totNumDays
!-------------------------------------------------------------------------
      subroutine Get_curGmiDate (self, ymd)
          integer         , intent(out) :: ymd
          type(t_GmiClock), intent(in ) :: self
          ymd = self%curGmiDate
          return
      end subroutine Get_curGmiDate
!-------------------------------------------------------------------------
      subroutine Set_curGmiDate (self, ymd)
          integer         , intent(in   ) :: ymd
          type(t_GmiClock), intent(inout) :: self
           self%curGmiDate = ymd
          return
      end subroutine Set_curGmiDate
!-------------------------------------------------------------------------
      subroutine Get_curGmiTime (self, hms)
          integer         , intent(out) :: hms
          type(t_GmiClock), intent(in ) :: self
          hms = self%curGmiTime
          return
      end subroutine Get_curGmiTime
!-------------------------------------------------------------------------
      subroutine Set_curGmiTime (self, hms)
          integer         , intent(in   ) :: hms
          type(t_GmiClock), intent(inout) :: self
           self%curGmiTime = hms
          return
      end subroutine Set_curGmiTime
!-------------------------------------------------------------------------
      subroutine Get_totNumTimeSteps (self, totNumTS)
          integer         , intent(out) :: totNumTS
          type(t_GmiClock), intent(in ) :: self
          totNumTS = self%totNumTimeSteps
          return
      end subroutine Get_totNumTimeSteps
!-------------------------------------------------------------------------
      subroutine Set_totNumTimeSteps (self, totNumTS)
          integer         , intent(in   ) :: totNumTS
          type(t_GmiClock), intent(inout) :: self
           self%totNumTimeSteps = totNumTS
          return
      end subroutine Set_totNumTimeSteps
!-------------------------------------------------------------------------
      subroutine Get_numTimeSteps (self, numTS)
          integer         , intent(out) :: numTS
          type(t_GmiClock), intent(in ) :: self
          numTS = self%numTimeSteps
          return
      end subroutine Get_numTimeSteps
!-------------------------------------------------------------------------
      subroutine Set_numTimeSteps (self, numTS)
          integer         , intent(in   ) :: numTS
          type(t_GmiClock), intent(inout) :: self
           self%numTimeSteps = numTS
          return
      end subroutine Set_numTimeSteps
!-------------------------------------------------------------------------
      subroutine Get_gmiTimeStep (self, tdt)
          real*8          , intent(out) :: tdt
          type(t_GmiClock), intent(in ) :: self
          tdt = self%gmiTimeStep
          return
      end subroutine Get_gmiTimeStep
!-------------------------------------------------------------------------
      subroutine Set_gmiTimeStep (self, tdt)
          real*8          , intent(in   ) :: tdt
          type(t_GmiClock), intent(inout) :: self
           self%gmiTimeStep = tdt
          return
      end subroutine Set_gmiTimeStep
!-------------------------------------------------------------------------
      subroutine Get_gmiSeconds (self, gmi_secs)
          real*8          , intent(out) :: gmi_secs
          type(t_GmiClock), intent(in ) :: self
          gmi_secs = self%gmiSeconds
          return
      end subroutine Get_gmiSeconds
!-------------------------------------------------------------------------
      subroutine Set_gmiSeconds (self, gmi_secs)
          real*8          , intent(in   ) :: gmi_secs
          type(t_GmiClock), intent(inout) :: self
           self%gmiSeconds = gmi_secs
          return
      end subroutine Set_gmiSeconds
!-------------------------------------------------------------------------
      subroutine Set_leapYearFlag (leap_year_flag)
          integer, intent(in) :: leap_year_flag
          leapYearFlag = leap_year_flag
          return
      end subroutine Set_leapYearFlag
!-------------------------------------------------------------------------
      subroutine Get_leapYearFlag (leap_year_flag)
          integer, intent(out) :: leap_year_flag
          leap_year_flag = leapYearFlag
          return
      end subroutine Get_leapYearFlag
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GmiAdvanceClock
!
! !INTERFACE:
!
  subroutine GmiAdvanceClock (self)
!
  implicit none

   type(t_GmiClock), intent(inout) :: self
!
! !LOCAL VARIABLES:
  integer :: nsec
  real*8  :: tdt             ! model time step (s)
  integer :: num_time_steps  ! number of time steps
  integer :: nymd            ! year/month/day (YYYYMMDD)
  integer :: nhms            ! hour/min/sec   (HHMMSS)
  real*8  :: gmi_sec         ! total Gmimod seconds (s)
!
! !DESCRIPTION:
!  This routine steps the Gmimod clock one time step.
!
! !AUTHOR:
!
! !REVISION HISTORY:
!  Initial code.
!EOP
!-------------------------------------------------------------------------
!BOC

  call Get_curGmiDate  (self, nymd          )
  call Get_curGmiTime  (self, nhms          )
  call Get_numTimeSteps(self, num_time_steps)
  call Get_gmiSeconds  (self, gmi_sec       )
  call Get_gmiTimeStep (self, tdt           )
 
  num_time_steps = num_time_steps + 1
  gmi_sec        = num_time_steps * tdt

  nsec = ConvertTimeToSeconds (nhms) + Nint (tdt)

  if (nsec >= SECPDY) then
     nsec = nsec - SECPDY
     nymd = IncrementDate (nymd, 1)
  end if
 
  if (nsec < 00000) then
     nsec = SECPDY + nsec
     nymd = IncrementDate (nymd, -1)
  end if
 
  nhms = ConvertSecondstoTime (nsec)
 
  call Set_curGmiDate  (self, nymd          )
  call Set_curGmiTime  (self, nhms          )
  call Set_numTimeSteps(self, num_time_steps)
  call Set_gmiSeconds  (self, gmi_sec       )

  return
 
  end subroutine GmiAdvanceClock
!
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: updateGmiClock
!
! !INTERFACE:
!
  subroutine updateGmiClock (self, esmfClock)
!
! !USES:
      use ESMF_Mod

      implicit none
!
! !INPUT PARAMETERS:
      type(ESMF_Clock), intent(in) :: esmfClock
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_GmiClock), intent(inOut) :: self
!
! !LOCAL VARIABLES:
      integer(ESMF_KIND_I4) :: cYear, cMonth, cDay
      integer(ESMF_KIND_I4) :: cHour, cMin  , cSec
      integer               :: rc
      type(ESMF_Time)       :: currTime
!
! !DESCRIPTION:
!  Update the model clock information using ESMF clock.
!
! !AUTHOR:
! Jules Kouatchou, NASA/GSFC, Jules.Kouatchou@nasa.gov
!
! !REVISION HISTORY:
!  Initial code.
!EOP
!-------------------------------------------------------------------------
!BOC
      ! Get the ESMF clock's current time
      call ESMF_ClockGet(esmfClock, currTime=currTime, rc=rc)

      call ESMF_TimeGet(currTime, yy=cYear, mm=cMonth, dd=cDay, h=cHour, &
     &                  m=cMin, s=cSec, rc=rc)

      self%curGmiDate   = 10000*cYear + 100*cMonth + cDay
      self%curGmiTime   = 10000*cHour + 100*cMin   + cSec
      self%numTimeSteps = self%numTimeSteps + 1
      self%gmiSeconds   = self%numTimeSteps * self%gmiTimeStep

      return
 
      end subroutine updateGmiClock
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetSecondsFromJanuary1
!
! !INTERFACE:
!
  subroutine GetSecondsFromJanuary1 (nsec_jan1, nymd, nhms)
!
  implicit none
!
! !INPUT PARAMETERS:
  integer, intent(in)   :: nymd      ! year/month/day (YYYYMMDD)
  integer, intent(in)   :: nhms      ! hour/min/sec   (HHMMSS)
!
! !OUTPUT PARAMETERS:
  integer, intent(out)  :: nsec_jan1 ! seconds from Jan. 1st (s)
!
! !DESCRIPTION:
!  This routine returns the offset in seconds from Jan. 1st, 
!  given any nymd/nhms.
!
! !LOCAL VARIABLES:
  integer :: nsec_ymd
  integer :: nsec_hms
!
! !AUTHOR:
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
  nsec_ymd = ConvertDateToSeconds (nymd)

  nsec_hms = ConvertTimeToSeconds (nhms) 

  nsec_jan1 = nsec_ymd + nsec_hms

  return

  end subroutine GetSecondsFromJanuary1
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetDaysFromJanuary1
!
  subroutine GetDaysFromJanuary1 (nday_jan1, nymd)
!
  implicit none

! !INPUT PARAMETERS:
  integer, intent(in)  :: nymd           ! year/month/day (YYYYMMDD)
!
! !OUTPUT PARAMETERS:
  integer, intent(out) :: nday_jan1      ! days from Jan. 1st (days)
!
! !DESCRIPTION:
!   This routine returns the offset in days from Jan. 1st, given any
!   nymd/nhms.
!
! !LOCAL VARIABLES:
  integer :: nsec_ymd
!
! !AUTHOR:
! !HISTORY:
!
!EOP 
!-----------------------------------------------------------------------------
!BOC
      nsec_ymd  = ConvertDateToSeconds (nymd)

      nday_jan1 = (nsec_ymd / SECPDY) + 1

      return

      end subroutine GetDaysFromJanuary1
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IncrementDate


!-----------------------------------------------------------------------------
!
! ROUTINE
!   IncrementDate
!
! DESCRIPTION
!   This routine increments/decrements nymd by one day.
!
! ARGUMENTS
!   nymd  : year/month/day (YYYYMMDD)
!   dyinc : the day increment/decrement (1 or -1)
!
!-----------------------------------------------------------------------------

  function IncrementDate (nymd, dyinc)

  use GmiPrintError_mod   , only : GmiPrintError
 
  implicit none
 
! ----------------------
! Argument declarations.
! ----------------------

  integer :: nymd
  integer :: dyinc


! ----------------------
! Function declarations.
! ----------------------

  integer :: IncrementDate
 
! ----------------------
! Variable declarations.
! ----------------------

  character (len=75) :: err_msg

  logical :: is_ldy
  logical :: is_lyr

  integer :: ndy, nmon, nyr

! ----------------
! Begin execution.
! ----------------
 
  if ((dyinc /= -1) .and. (dyinc /= 1)) then
    err_msg = 'dyinc range error in IncrementDate.'
    call GmiPrintError (err_msg, .true., 1, dyinc, 0, 0, 0.0d0, 0.0d0)
  end if

  nyr  = nymd / 10000
  nmon = Mod (nymd, 10000) / 100
  ndy  = Mod (nymd, 100) + dyinc
 
  if (ndy == 0) then

!   -------------------------
!   Return to previous month.
!   -------------------------

    nmon = nmon - 1
 
    if (nmon == 0) then
      nmon = MONTHS_PER_YEAR
      nyr = nyr - 1
!      if (nyr == 0) then
!        nyr = 99
!      else
!        nyr = nyr - 1
!      end if
    end if
 
    is_lyr = IsLeapYear (nyr)

    if (is_lyr) then
      ndy = DAYS_PER_MONTH_LY(nmon)
    else
      ndy = DAYS_PER_MONTH(nmon)
    end if

  else

    is_lyr = IsLeapYear (nyr)

    if (is_lyr .and. (nmon == 2) .and. (ndy == DAYS_PER_MONTH_LY(nmon))) then
      is_ldy = .true.
    else
      is_ldy = .false.
    end if

    if (.not. is_ldy) then

      if (ndy > DAYS_PER_MONTH(nmon)) then

!       ----------------------
!       Proceed to next month.
!       ----------------------
 
        ndy  = 1
        nmon = nmon + 1
 
        if (nmon > MONTHS_PER_YEAR) then
          nmon = 1
          nyr = nyr + 1
!          if (nyr == 99) then
!            nyr = 0
!          else
!            nyr = nyr + 1
!          end if
        end if
 
      end if

    end if

  end if

  IncrementDate = (nyr * 10000) + (nmon * 100) + ndy
 
  return
 
  end function IncrementDate
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IsLeapYear
 

!-----------------------------------------------------------------------------
!
! ROUTINE
!   IsLeapYear
!
! DESCRIPTION
!   This routine always returns false if leapYearFlag < 0,
!                always returns true  if leapYearFlag > 0,
!                and calculates whether or not it is really a leap year if
!                  leapYearFlag = 0.
!
! ARGUMENTS
!   the_year : the year to check (last 2 digits)
!
!-----------------------------------------------------------------------------

      function IsLeapYear (the_year)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: the_year
 
!     ----------------------
!     Function declarations.
!     ----------------------

      logical :: IsLeapYear

!     ----------------------
!     Variable declarations.
!     ----------------------

!c    integer, save :: the_year00 = 1900
      integer, save :: the_year00 = 2000


!     ----------------
!     Begin execution.
!     ----------------

      IsLeapYear = .false.


      if (leapYearFlag > 0) then

        IsLeapYear = .true.

      else if (leapYearFlag == 0) then

        if (Mod (the_year, 4) == 0) then

          if ((the_year /= 0) .or. (Mod (the_year00, 400) == 0)) then

            IsLeapYear = .true.

          end if

        end if

      end if

      return

      end function IsLeapYear
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConvertTimeToSeconds


!-----------------------------------------------------------------------------
!
! ROUTINE
!   ConvertTimeToSeconds
!
! DESCRIPTION
!   This routine converts nhms to total seconds.
!
! ARGUMENTS
!   nhms : hour/min/sec (HHMMSS)
!
!-----------------------------------------------------------------------------

      function ConvertTimeToSeconds (nhms)

     use GmiPrintError_mod   , only : GmiPrintError
 
      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: nhms

 
!     ----------------------
!     Function declarations.
!     ----------------------

      integer :: ConvertTimeToSeconds


!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      integer :: nsec

!     ----------------
!     Begin execution.
!     ----------------

      nsec = (nhms / 10000) * SECPHR +              &
             (Mod (nhms, 10000) / 100) * SECPMN +   &
             Mod (nhms, 100)

      if (nsec >= SECPDY) then

        err_msg = 'nsec too big in ConvertTimeToSeconds.'
        call GmiPrintError (err_msg, .true., 1, nsec, 0, 0, 0.0d0, 0.0d0)

      end if

      ConvertTimeToSeconds = nsec

      return

      end function ConvertTimeToSeconds
 

!-----------------------------------------------------------------------------
!
! ROUTINE
!   ConvertSecondstoTime
!
! DESCRIPTION
!   This routine converts total seconds to hour/min/sec.
!
! ARGUMENTS
!   nsec : total seconds
!
!-----------------------------------------------------------------------------

  function ConvertSecondstoTime (nsec)

    use GmiPrintError_mod   , only : GmiPrintError
 
      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: nsec

 
!     ----------------------
!     Function declarations.
!     ----------------------

      integer :: ConvertSecondstoTime


!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      integer :: nhms


!     ----------------
!     Begin execution.
!     ----------------

      if (nsec >= SECPDY) then

        err_msg = 'nsec too big in ConvertSecondstoTime.'
        call GmiPrintError (err_msg, .true., 1, nsec, 0, 0, 0.0d0, 0.0d0)

      end if

      nhms = (nsec / SECPHR) * 10000 +                 &
             (Mod (nsec, SECPHR) / SECPMN) * 100 +     &
             Mod (nsec, SECPMN)


      ConvertSecondstoTime = nhms

      return

      end function ConvertSecondstoTime
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConvertDateToSeconds
!
! !INTERFACE:
!
  function ConvertDateToSeconds (nymd)
!
! !USES:
  use GmiPrintError_mod   , only : GmiPrintError
!
  implicit none

!
! !INPUT PARAMETERS:
  integer, intent(in)   :: nymd   ! year/month/day (YYYYMMDD)
!
! !RETURNED VALUE:
  integer :: ConvertDateToSeconds
!
! !DESCRIPTION:
!  This routine converts nymd to total seconds.
!
! !LOCAL VARIABLES:
  character (len=75) :: err_msg
  logical            :: is_lyr
  integer            :: idays, isecs
  integer            :: iyrsec
  integer            :: ndy, nmon, nyr
!
! !AUTHOR:
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

  nyr  = nymd / 10000
  nmon = Mod (nymd, 10000) / 100
  ndy  = Mod (nymd, 100)
 
  is_lyr = IsLeapYear (nyr)

  if (is_lyr) then
     idays  = (START_DAY_OF_MONTH_LY(nmon) - 1) + (ndy - 1)
     iyrsec = SECPYR_LY
  else
     idays  = (START_DAY_OF_MONTH(nmon)    - 1) + (ndy - 1)
     iyrsec = SECPYR
  end if
 
  isecs = idays * SECPDY

  if (isecs >= iyrsec) then
     err_msg = 'isecs too big in ConvertDateToSeconds.'
     call GmiPrintError (err_msg, .true., 1, isecs, 0, 0, 0.0d0, 0.0d0)
  end if

  ConvertDateToSeconds = isecs

  return

  end function ConvertDateToSeconds
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GmiSplitDateTime 
!
! !INTERFACE:
!
  subroutine GmiSplitDateTime (nlmr, left2, middle2, right2)
!
  implicit none
!
! !INPUT PARAMETERS:
  integer, intent(in)   :: nlmr    ! nymd or nhms ("lmr" = "left" "middle" "right")
                                   ! (YYYYMMDD or HHMMSS)
!
! !OUTPUT PARAMETERS:
  integer, intent(out)  :: left2   ! left   field (year  or hours)
  integer, intent(out)  :: middle2 ! middle field (month or minutes)
  integer, intent(out)  :: right2  ! right  field (day   or seconds)
!
! !DESCRIPTION:
!  This routine extracts the three fields from either nymd 
!  (year, month, day) or nhms (hours, minutes, seconds).
!
! !AUTHOR:
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  left2   = nlmr / 10000
  middle2 = Mod (nlmr, 10000) / 100
  right2  = Mod (nlmr, 100)

  return

  end subroutine GmiSplitDateTime
!EOC
!-------------------------------------------------------------------------
!BOP
      subroutine initializeDate &
                  (tbegin_days, tbegin, tfinal_days, tfinal, time_esm, &
                   communicator)

      implicit none

      real*8 , intent(in) :: tbegin_days
      real*8 , intent(in) :: tfinal_days
      real*8 , intent(out) :: tbegin
      real*8 , intent(out) :: tfinal
      real*8   :: time_esm
      integer, intent(in) :: communicator

! !DESCRIPTION:
!   This routine sets the starting time and date of the simulation from
!   user input data.  The finishing time is calculated from user input data
!   and reconciled with restrictions required by the constituent models.
!
! !LOCAL VARIABLES:
      integer :: day_of_week, week_of_year
      integer :: month_of_year

      integer :: days_since0, weeks_since0
      integer :: months_since0, years_since0
      real*8  :: smaltime, gmt_esm

      integer  :: leapYear, year_esm, imonth_esm
      integer  :: dayomnth_esm, dayoyear_esm
      character (len=9)   :: month_esm
      character (len=27)  :: mmddyy_esm
      integer :: daypmnth(14)

!EOP
!-------------------------------------------------------------------------
!BOC

      smaltime     = 0.0d0 ! a small time to prevent exactly even hours

      if (tbegin_days /= 0.0d0) tbegin = tbegin_days * SECPDY
      if (tfinal_days /= 0.0d0) tfinal = tfinal_days * SECPDY

      time_esm = tbegin

      if (tfinal <= tbegin) then

        Write (6,800)

 800    format ('ERROR: The final time is <= the ', //,  &
     &          'beginning time.  Check tbegin_days and tfinal_days.')

        call stopCode (communicator, "Check the namelist file.")

      end if

      if ((tbegin_days < 0.0d0) .or.  &
     &    (tfinal_days < 0.0d0) .or.  &
     &    (tbegin      < 0.0d0) .or.  &
     &    (tfinal      < 0.0d0) .or.  &
     &    (refDay      < 0)) then

        Write (6,810)
        Write (6,820)
        Write (6,830)

 810    format ('ERROR: Negative values of time are not allowed.')
 820    format ('       The default origin of time is day0 = 1.')
 830    format ('       Please reset tbegin_days and tfinal_days.')

        call stopCode (communicator, "See the section of the input deck.")

      end if

!     ==============
      call setClock  &
!     ==============
     &  (time_esm, leapYear,  &
     &   mmddyy_esm, gmt_esm, year_esm, dayomnth_esm,  &
     &   imonth_esm, month_esm,  &
     &   dayoyear_esm, day_of_week, week_of_year, month_of_year,  &
     &   days_since0, weeks_since0, months_since0, years_since0, &
     &   daypmnth, smaltime, communicator)

      return

      end subroutine initializeDate
!EOC
!-------------------------------------------------------------------------
!BOP
!
! DESCRIPTION
!   This routine tells time and reads the calendar.  It does this by
!   calculating the time of day at the Prime Meridian (known as Greenwich
!   Mean Time).
!
!   For perpetual month calculations the day of the month is held constant
!   and equal to the first day of the simulation.
!
!   The algorithms used to calculate the time/date were adopted from
!   Ken Caldeira's calendar function, written in Mathematica.  It's kludgy and
!   served the purpose for a quick addition to the esm_clock function.  This
!   should be replaced eventually.
!
!   Also note that the variables, that were returned before the arguments
!   yearx, day_of_year, day_of_week, week_of_year, month_of_year,
!   days_since_year0, weeks_since_year0, & years_since_year0 were added have
!   to be "manipulated" to return the correct values which are different from
!   those without the calendar functions.  This required some adjustments in
!   routines that use those values.  The same adjustments need to be made when
!   other packages start using the calendar scheme.
!
! ARGUMENTS
!   telapsed : elapsed time in the package this routine was called from (s)
!   leapYear : integer flag to indicate leap years
!   mmddyy   : time, day, month year in character form
!   gmt      : Greenwich Mean Time
!   year     : the years since the beginning of the calculation
!   dayomnth : the day of the month
!   imonth   : the month in integer   form
!   month    : the month in character form
!   day_of_year        : the day of the year
!   day_of_week        : tbd
!   week_of_year       :
!   month_of_year      :
!   days_since_refYear   :
!   weeks_since_refYear  :
!   months_since_refYear :
!   years_since_refYear  :
!
!-----------------------------------------------------------------------------

      subroutine setClock  &
     &  (telapsed, leapYear, mmddyy, gmt, year, dayomnth,  &
     &   imonth, month, day_of_year, day_of_week, week_of_year,  &
     &   month_of_year, days_since_refYear, weeks_since_refYear,  &
     &   months_since_refYear, years_since_refYear, &
     &   daypmnth, smaltime, communicator)

      implicit none

! !INPUT PARAMETERS:
      real*8 , intent(in) :: smaltime
      integer, intent(in) :: communicator
! !OUTPUT PARAMETERS:
      real*8  :: telapsed
      integer, intent(out) :: leapYear
      character (len=*) :: mmddyy
      integer :: daypmnth(14)
      real*8  :: gmt
      integer :: year
      integer :: dayomnth
      integer :: imonth
      character (len=*) :: month
      integer :: day_of_year
      integer :: day_of_week
      integer :: week_of_year
      integer :: month_of_year
      integer :: days_since_refYear
      integer :: weeks_since_refYear
      integer :: months_since_refYear
      integer :: years_since_refYear
! !LOCAL VARIABLES:
      character (len=2) :: day2, hour2, minutes2
      character (len=9) :: month9
      character (len=4) :: year4

      integer :: hour, minutes
      integer :: i, iday

      real*8  :: day, hours
      real*8  :: rday_of_year, ryear
!EOP
!------------------------------------------------------------------------------
!BOC
!     ---------------------------------
!     Calculate the clock and calendar.
!     ---------------------------------

      day_of_year   = 0
      day_of_week   = 0
      week_of_year  = 0
      month_of_year = 0

      days_since_refYear   = 0
      weeks_since_refYear  = 0
      months_since_refYear = 0
      years_since_refYear  = 0

!       -----------------------------------------------------------------
!       Realistic calendar with leap years and leap centuries (preferred).
!       -----------------------------------------------------------------

!       ==================
        call findDate  &
!       ==================
     &    (telapsed, day_of_year, day_of_week, week_of_year,  &
     &     month_of_year, days_since_refYear, weeks_since_refYear,  &
     &     months_since_refYear, years_since_refYear, hour, minutes,  &
     &     dayomnth, leapYear, communicator)


        imonth = month_of_year
        month9 = monthsYear(month_of_year)
        hours  = telapsed / SECPHR
        iday   = hours / HRS_PER_DAY
        gmt    = hours + smaltime - (iday * HRS_PER_DAY)
        year   = years_since_refYear + refYear

        daypmnth(1)  = 0
        daypmnth(2)  = 31

        if (leapYear == 0) daypmnth(3) = daypmnth(2) + 28
        if (leapYear == 1) daypmnth(3) = daypmnth(2) + 29

        daypmnth(4)  = daypmnth(3)  + 31
        daypmnth(5)  = daypmnth(4)  + 30
        daypmnth(6)  = daypmnth(5)  + 31
        daypmnth(7)  = daypmnth(6)  + 30
        daypmnth(8)  = daypmnth(7)  + 31
        daypmnth(9)  = daypmnth(8)  + 31
        daypmnth(10) = daypmnth(9)  + 30
        daypmnth(11) = daypmnth(10) + 31
        daypmnth(12) = daypmnth(11) + 30
        daypmnth(13) = daypmnth(12) + 31
        daypmnth(14) = daypmnth(13) + 31

        Write (day2,     '(i2.2)') dayomnth
        Write (year4,    '(i4.4)') year
        Write (hour2,    '(i2.2)') hour
        Write (minutes2, '(i2.2)') minutes

        mmddyy = month9 // ' ' // day2 // ',' // year4 // ' ' //  &
     &           hour2  // ':' // minutes2 // ' GMT'

      return

      end subroutine setClock
!EOC
!-----------------------------------------------------------------------------
!BOP
! DESCRIPTION
!   This routine finds the date given secs after midnight between year -1 and
!   0.  It was adopted from Ken Caldeira's function calendar, written in
!   Mathematica.
!
!   This routine takes leap years into account.
!
!   Note that the variables hour, minutes, dayomnth, & year are redundant;
!   they should be consolidated.
!
! ARGUMENTS
!   time_since_refYear   : time (s)
!   refYear              : any year you want to use as base
!   day_of_year        : the day   of current year
!   day_of_week        : the day   of current week
!   week_of_year       : the week  of current year
!   month_of_year      : the month of current year
!   days_since_refYear   : days   relative to Jan1, refYear
!   weeks_since_refYear  : weeks  relative to Jan1, refYear
!   months_since_refYear : months relative to Jan1, refYear
!   years_since_refYear  : years  relative to Jan1, refYear
!   hour               : tbd
!   minutes            :
!   dayomnth           : the day of the month
!   leapYear           : integer flag to indicate leap years
!
!-----------------------------------------------------------------------------

      subroutine findDate  &
     &  (time_since_refYear, day_of_year, day_of_week,  &
     &   week_of_year, month_of_year, days_since_refYear,  &
     &   weeks_since_refYear, months_since_refYear, years_since_refYear,  &
     &   hour, minutes, dayomnth, leapYear, communicator)

      implicit none

! !INPUT PARAMETERS:
      real*8 , intent(in) :: time_since_refYear
      integer, intent(in) :: communicator
! !OUTPUT PARAMETERS:
      integer, intent(out) :: day_of_year,  day_of_week
      integer, intent(out) :: week_of_year, month_of_year
      integer, intent(out) :: days_since_refYear,   weeks_since_refYear
      integer, intent(out) :: months_since_refYear, years_since_refYear
      integer, intent(out) :: hour, minutes
      integer, intent(out) :: dayomnth
      integer, intent(out) :: leapYear

! !LOCAL VARIABLES:
      integer :: days_since_0,   weeks_since_0  ! relative to Jan1, 0 AD
      integer :: months_since_0, years_since_0  ! relative to Jan1, 0 AD
      integer :: days_at_refYear,   weeks_at_refYear
      integer :: months_at_refYear, years_at_refYear
      real*8  :: time_of_refYear
      real*8  :: time_since_0
!EOP
!------------------------------------------------------------------------------
!BOC
      days_since_refYear  = 0
      weeks_since_refYear = 0
      years_since_refYear = 0

!     =======================
      call getSecsYear0  &
!     =======================
     &  (time_since_refYear, time_since_0)

      time_of_refYear = time_since_0 - time_since_refYear

!     ----------------------------------------------------------------
!     This call is made to determine the reference at Jan 1, refYear AD.
!     ----------------------------------------------------------------

!     =======================
      call getTimeYear0  &
!     =======================
     &  (time_of_refYear,  &
     &   days_at_refYear, weeks_at_refYear, months_at_refYear, years_at_refYear,  &
     &   day_of_year, day_of_week, week_of_year, month_of_year,  &
     &   hour, minutes, dayomnth, leapYear)


!     --------------------------------------
!     This call is made at the current time.
!     --------------------------------------

!     =======================
      call getTimeYear0  &
!     =======================
     &  (time_since_0,  &
     &   days_since_0, weeks_since_0, months_since_0, years_since_0,  &
     &   day_of_year, day_of_week, week_of_year, month_of_year,  &
     &   hour, minutes, dayomnth, leapYear)


      days_since_refYear   = days_since_0   - days_at_refYear
      weeks_since_refYear  = weeks_since_0  - weeks_at_refYear
      months_since_refYear = months_since_0 - months_at_refYear
      years_since_refYear  = years_since_0  - years_at_refYear


      return

      end subroutine findDate


!-----------------------------------------------------------------------------
!
! ROUTINE
!   getSecsYear0
!
! DESCRIPTION
!   This routine, given the time in seconds and the year passed in as
!   arguments, calculates and returns the number of seconds since year 0 to
!   the current year and seconds in that year.
!
! ARGUMENTS
!   secs       : base time (s)
!   timeInSecs : the time (s)
!
!-----------------------------------------------------------------------------

      subroutine getSecsYear0 (secs, timeInSecs)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      real*8 , intent(in   ) :: secs
      real*8 , intent(out) :: timeInSecs

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: leapYear
      integer :: n4yrs, n100yrs, n400yrs
      integer :: year0remaining

      real*8  :: secs0
      real*8  :: secsremaining


!     ----------------
!     Begin execution.
!     ----------------

      year0remaining = refYear
      n400yrs = int(year0remaining/400.)
      secs0 = secsp400yrs*n400yrs
      year0remaining = year0remaining - 400*n400yrs

!  is current year a leap year?

      if(mod(refYear,400) .eq. 0) then
        leapYear = 1
      elseif(mod(refYear,4).eq.0.and.mod(refYear,100).ne.0) then
        leapYear = 1
      else
        leapYear = 0
      end if

!  account for the first year of the 400 year block having an extra
!  leap day

      if (year0remaining .gt. 0) then
        secs0 = secs0 + secsPerDay
      end if

      n100yrs = int(year0remaining/100.)
      secs0 = secs0  + secspnrmlcentry*n100yrs
      year0remaining = year0remaining - 100*n100yrs

      n4yrs = int(year0remaining/4.)
      secs0 = secs0 + secspnrml4yrblk*n4yrs
      year0remaining = year0remaining - 4*n4yrs

      secs0 = secs0 + year0remaining*secspnrmlyr
      secsremaining = secs + secs0
      timeInSecs = secsremaining

      if (leapYear .eq. 1 .and. refYear .ne. 0) then
        timeInSecs = timeInSecs - secsPerDay
      end if


      return

      end subroutine getSecsYear0


!-----------------------------------------------------------------------------
!
! ROUTINE
!   getTimeYear0
!
! DESCRIPTION
!   This routine tbd
!
! ARGUMENTS
!
!-----------------------------------------------------------------------------

      subroutine getTimeYear0  &
     &  (secsremaining, days_since0, weeks_since0, months_since0,  &
     &   years_since0, day_of_year, day_of_week, week_of_year,  &
     &   month_of_year, hourout, minout, dayout, leapYear)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      real*8 , intent(in ) :: secsremaining
      integer, intent(out) :: leapYear
      integer, intent(out) :: days_since0,   weeks_since0
      integer, intent(out) :: months_since0, years_since0
      integer, intent(out) :: day_of_year,   day_of_week
      integer, intent(out) :: week_of_year,  month_of_year
      integer, intent(out) :: hourout, minout, dayout


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: centuries
      integer :: days, month, years
      integer :: daysinmnth
      integer :: fouryrblks
      integer :: i
      integer :: monthout, yearout
      integer :: n400yrs
      integer :: stopflag

      integer :: dayspmnth(12) =  &
     &             (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

      real*8  :: secsleft

!     ----------------
!     Begin execution.
!     ----------------

      yearout = 0
      monthout = 0
      dayout = 0
      hourout = 0
      minout = 0

      n400yrs = int(secsremaining/secsp400yrs)
      yearout = yearout + 400*n400yrs
      secsleft = secsremaining -n400yrs*secsp400yrs

!  takes care of leapyear every 400 years

      if (secsleft .ge. secspleapyr) then
        yearout = yearout + 1
        secsleft = secsleft - secspleapyr

!  now remaining centuries are normal

        centuries = int(secsleft/secspnrmlcentry)
        yearout = yearout + 100*centuries
        secsleft = secsleft -  centuries*secspnrmlcentry

!  now the remaining 4-yr blocks are normal

        fouryrblks = int(secsleft/secspnrml4yrblk)
        yearout = yearout + 4*fouryrblks
        secsleft = secsleft - fouryrblks*secspnrml4yrblk

! since we took out year zero, if any year is a leap year, it
! will be the 4th year remaining

        years = int(secsleft/secspnrmlyr)

! If years = 4, the current year is a leap year

        if (years.eq.4) years = 3
        yearout = yearout + years

        secsleft = secsleft - years*secspnrmlyr

      end if

!  is current year a leap year?

      if(mod(yearout,400) .eq. 0) then
        leapYear = 1
      elseif(mod(yearout,4).eq.0.and.mod(yearout,100).ne.0) then
        leapYear = 1
      else
        leapYear = 0
      end if

!   find current month

      stopflag = 0
      month = 0

      do i = 1,12

        if (stopflag .eq. 0) then
          month = month +1
          daysinmnth = dayspmnth(month)
          if(month .eq. 2.) daysinmnth = daysinmnth + leapYear
          if(secsleft .ge. daysinmnth*secsPerDay) then
            secsleft = secsleft - daysinmnth*secsPerDay
          else
            stopflag = 1
          end if
        end if

      end do

      monthout = month

!  get rest of date and time

      days = int(secsleft/secsPerDay)
      secsleft = secsleft - days*secsPerDay
      dayout = days + 1

      hourout = int(secsleft/secsPerHour)
      secsleft = secsleft - hourout*secsPerHour

      minout = int(secsleft/secsPerMin)
      secsleft = secsleft - minout*secsPerMin

      days_since0   = secsremaining / secsPerDay
      weeks_since0  = secsremaining / secsPerWeek
      years_since0  = yearout
      months_since0 = years_since0 * 12 + monthout

!  calculate the day of the year

      day_of_year = 0

      if (monthout .eq. 1) then
        day_of_year = dayout
      else

        do i = 2, monthout
          day_of_year =  day_of_year + dayspmnth(i-1)
        end do

        day_of_year = day_of_year + dayout

        if (monthout .ge. 3) then
          day_of_year = day_of_year + leapYear
        end if

      end if

!  calculate day of the week

      day_of_week = mod(days_since0,7)

      if(day_of_week .eq. 0) day_of_week = 7

! adjust for the very first 7 days in refYear and for days after leapyear
! for year 0

      if(days_since0 .le. 7) day_of_week = days_since0 + 1

      if(day_of_week .eq. 8) day_of_week = 1

      week_of_year = day_of_year/7 + 1

      month_of_year = monthout

      return

      end subroutine getTimeYear0
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: isOutTime
!
! !INTERFACE:
!
      subroutine isOutTime(outTime, printed_on_this_day, &
     &             month_save, month, day, nhms, ndt, gmi_sec, pr_period)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: month_save ! the month during the last time step
      integer, intent(in) :: month      ! the current month
      integer, intent(in) :: day        ! the current day
      integer, intent(in) :: nhms       ! the current time HHMMSS
      integer, intent(in) :: ndt        ! model time step (s)
      real*8 , intent(in) :: gmi_sec    ! total GMI seconds (s)
      real*8 , intent(in) :: pr_period  ! periodic printing interval (s)
!
! !OUTPUT PARAMETERS:
      logical, intent(out) :: outTime
!
! !INPUT/OUTPUT PARAMETERS:
      logical, intent(inOut) :: printed_on_this_day
!
! !DESCRIPTION:
! Determines if it is time to do some output based on the
!   following values for pr_period:
! \begin{itemize}
! \item >0.0d0:  output at specified interval (s)
! \item -1.0d0:  output at monthly intervals
! \item -2.0d0:  output on 1st & 15th of each month
! \end{itemize}
!
! !LOCAL VARIABLES:
      logical :: is_time
!      logical, save :: printed_on_this_day = .false.
      integer :: prflg
      real*8  :: rsecpdy
!EOP
!------------------------------------------------------------------------------
!BOC
      outTime = .false.
      rsecpdy = SECPDY

!     ======================
      if (pr_period < 0.0d0) then
!     ======================
         prflg = Nint (pr_period / rsecpdy)
         if (prflg == -1) then

            !Monthly output.
            if (month /= month_save) outTime = .true.
         else if (prflg == -2) then

            ! Bimonthly output (1st & 15th).
            if ((day == 1) .or. (day == 15)) then
               if ((nhms >= 0) .and. (.not. (printed_on_this_day))) then
                  outTime = .true.
                  printed_on_this_day = .true.
               end if
            else
               printed_on_this_day = .false.
            end if
         end if

!     ======================================================
      else if (Mod (Nint (gmi_sec), Nint (pr_period)) < ndt) then
!     ======================================================

         ! Specified time for periodic output.
         outTime = .true.
!     ======
      end if
!     ======

      return

      end subroutine isOutTime
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: isOutFreqTime
!
! !INTERFACE:
!
      subroutine isOutFreqTime (outFreqTime, printed_on_this_day, month_save, &
     &               month, day, nhms, ndt, gmi_sec, pr_period, pr_at_time)
!
      implicit none

!
! !INPUT PARAMETERS:
      integer, intent(in) :: month_save  ! the month during the last time step
      integer, intent(in) :: month       ! the current month
      integer, intent(in) :: day         ! the current day
      integer, intent(in) :: nhms        ! the current hour/minute/second
      integer, intent(in) :: ndt         ! model time step (s)
      real*8 , intent(in) :: gmi_sec     ! total Gmimod seconds (s)
      real*8 , intent(in) :: pr_period   ! periodic printing interval (s)
      real*8 , intent(in) :: pr_at_time
!
! !OUTPUT PARAMETERS:
      logical, intent(out) :: outFreqTime
!
! !INPUT/OUTPUT PARAMETERS:
      logical, intent(inOut) :: printed_on_this_day
!
! !DESCRIPTION:
!   Determines if it is time to do some output based on the
!   following values for pr_period:
! \begin{itemize}
! \item >0.0d0:  output at specified interval (days)
! \item -1.0d0:  output at monthly intervals
! \item -2.0d0:  output on 1st & 15th of each month
! \item -3.0d0:  output on 1st, 6th, 11th, 16th, 21st & 26th of each month
! \end{itemize}
!
! !LOCAL VARIABLES:
      integer :: print_days(6), ipr
      integer :: prflg
      real*8  :: rsecpdy
!EOP
!------------------------------------------------------------------------------
!BOC
      outFreqTime = .false.

      rsecpdy = SECPDY

!     ======================
      if (pr_period < 0.0d0) then
!     ======================

        prflg = Nint (pr_period / rsecpdy)

        if (prflg == -1) then

!         ---------------
!         ---------------
!         Monthly output.
!         ---------------

          if (month /= month_save) outFreqTime = .true.

        else if (prflg == -2) then

!         ------------------------------
!         Bimonthly output (1st & 15th).
!         ------------------------------

          if ((day == 1) .or. (day == 15)) then

            if ((nhms >= 0) .and. (.not. (printed_on_this_day))) then
              outFreqTime = .true.
              printed_on_this_day = .true.
            end if

          else

            printed_on_this_day = .false.

          end if

        else if (prflg == -3) then
!         ------------------------------
!         Instantaneous output (1st, 6th, 11th, 16th, 21st, 26th).
!         (Note: this section added by Bigyani)
!         ------------------------------
          data (print_days(ipr), ipr=1,6) /1, 6, 11, 16, 21, 26/

          if (any(print_days == day)) then

            if ((nhms >= 0) .and. (.not. (printed_on_this_day))) then
              outFreqTime = .true.
              printed_on_this_day = .true.
            end if

          else

            printed_on_this_day = .false.

          end if

        end if

!     ======================================================
      else if (Mod (Nint (gmi_sec+pr_at_time), Nint (pr_period)) < ndt) then
!     ======================================================

!       -----------------------------------
!       Specified time for periodic output.
!       -----------------------------------

        outFreqTime = .true.

!     ======
      end if
!     ======

      return

      end subroutine isOutFreqTime
!EOC
!-----------------------------------------------------------------------------
!BOP
! !IROUTINE: stopCode
!
! !INTERFACE:
!
      subroutine stopCode (commu, msg)

      implicit none
!
! !INPUT PARAMETERS:
      integer          , intent(in) :: commu
      character (len=*), intent(in) :: msg
!
! !DESCRIPTION:
! Stops all the processes associated with the communicator.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: icode, ierr

      Write (6,*) msg

      call MPI_Abort (commu, icode, ierr)

      Stop

      return

      end subroutine stopCode
!EOC
!------------------------------------------------------------------------------


  end module GmiTimeControl_mod
 
