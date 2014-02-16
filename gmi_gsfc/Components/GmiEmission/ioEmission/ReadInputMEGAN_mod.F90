!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: ReadInputMEGAN_mod
!
! !INTERFACE:
!
     module ReadInputMEGAN_mod
!
! !USES:
     use m_netcdf_io_open  , only : Ncop_Rd
     use m_netcdf_io_close , only : Nccl
     use m_netcdf_io_read  , only : Ncrd_3d, Ncrd_4d
     use GmiTimeControl_mod, only : GmiSplitDateTime, GetDaysFromJanuary1
!
     implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
     private
     public  :: setMEGANisoLAI
     public  :: readMEGANannualEmissFactor

#     include "GmiParameters.h"
#     include "gmi_time_constants.h"

! !DESCRIPTION:
! Reads and stores AVHRR LAI for calculating MEGAN biogenic VOC emissions.
! (dsa, tmf, bmy, 10/20/05)
!
! !AUTHOR:
!  Bob Yantosca
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!------------------------------------------------------------------------------
     contains
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setMEGANisoLAI
!
! !INTERFACE:
!
      subroutine setMEGANisoLAI &
     &           (isoLai, isoLaiCurr, isoLaiPrev, isoLaiNext, &
     &            days_btw_m, laiMEGAN_InfileName, &
     &            nymd, i1, i2, ju1, j2, i1_gl, ju1_gl, currMonthMEGAN)

      implicit none

! !INPUT PARAMETERS:
      integer            , intent(in) :: i1, i2, ju1, j2, i1_gl, ju1_gl
      integer            , intent(in) :: nymd
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: laiMEGAN_InfileName
!
! !OUTPUT PARAMETERS:
                              ! days between midmonths in the LAI data
      integer, intent(out) :: days_btw_m 
!
! !INPUT/OUTPUT PARAMETERS:
                                ! Daily LAI data (interpolate)
      real*8 , intent(inOut) :: isoLai    (i1:i2,ju1:j2)
                                ! AVHRR LAI data for the current  month
      real*8 , intent(inOut) :: isoLaiCurr(i1:i2,ju1:j2)
                                ! AVHRR LAI data for the previous month
      real*8 , intent(inOut) :: isoLaiPrev(i1:i2,ju1:j2)
                                ! AVHRR LAI data for the next     month
      real*8 , intent(inOut) :: isoLaiNext(i1:i2,ju1:j2)
      integer, intent(inOut) :: currMonthMEGAN
!
! !DESCRIPTION: 
!  Sets isoLai daily.  The stored monthly LAI are used for
!  the middle day in the month and LAIs are interpolated for other days.
!  (dsa, tmf, bmy, 10/20/05)
!
! !LOCAL VARIABLES:
      integer                :: thisMonth, calMonth, jday
      integer                :: ddummy, ydummy, i, j
      integer                :: imul ! days since midmonth
      real*8                 :: fraction

! !DEFINED PARAMETERS:
      ! specify midmonth day for year 2000
      INTEGER, PARAMETER     :: startDay(MONTHS_PER_YEAR+1) = &
     &                         (/  15,  45,  74, 105, 135, 166, &
     &                            196, 227, 258, 288, 319, 349, 380/)
!
! !AUTHOR:
!  Bob Yantosca and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Extract the month
      call GmiSplitDateTime    (nymd, ydummy, calMonth, ddummy)

      ! Extract the number of days since January 1st
      call GetDaysFromJanuary1 (jday, nymd)

      thisMonth = currLAImonth(jday, calMonth, startDay)

      ! Read new data if it is a new LAI month
      if (currMonthMEGAN /= thisMonth) then
         currMonthMEGAN = thisMonth
         call readIsoLai (isoLaiCurr, isoLaiPrev, isoLaiNext, laiMEGAN_InfileName, &
     &            currMonthMEGAN, i1, i2, ju1, j2, i1_gl, ju1_gl)
      end if

      ! IMUL is days since midmonth
      ! ITD  is days between midmonths
      if ( jday < startDay(1) ) then
         imul = 365 + jday - startDay(12)
         days_btw_m  = 31
      ELSE
         imul = jday                  - startDay(currMonthMEGAN)
         days_btw_m  = startDay(currMonthMEGAN+1) - startDay(currMonthMEGAN)
      END IF

      ! Fraction of the LAI month that we are in
      fraction       = dble( imul ) / dble( days_btw_m )

      ! Interpolate to daily LAI value
      
      isoLai(:,:) = isoLaiCurr(:,:) + fraction * (isoLaiNext(:,:) - isoLaiCurr(:,:)) 
      return

      end subroutine setMEGANisoLAI
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readIsoLai
!
! !INTERFACE:
!
      subroutine readIsoLai &
     &           (isoLaiCurr, isoLaiPrev, isoLaiNext, laiMEGAN_InfileName, &
     &            currMonth, i1, i2, ju1, j2, i1_gl, ju1_gl)

      implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: i1, i2, ju1, j2, i1_gl, ju1_gl
      integer, intent(in) :: currMonth
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: laiMEGAN_InfileName
!
! !OUTPUT PARAMETERS:
                              ! AVHRR LAI data for the current  month
      real*8 , intent(out) :: isoLaiCurr(i1:i2, ju1:j2)
                              ! AVHRR LAI data for the previous month
      real*8 , intent(out) :: isoLaiPrev(i1:i2, ju1:j2)
                              ! AVHRR LAI data for the next     month
      real*8 , intent(out) :: isoLaiNext(i1:i2, ju1:j2)

! !DESCRIPTION: 
!  Reads AVHRR LAI data from bpch file for the current month, 
!  the previous month, and the next month. (dsa, tmf, bmy, 10/18/05).
!
! !LOCAL VARIABLES:
      integer :: prevMonth, nextMonth
      integer :: count4d (4), start4d(4)
      integer :: ncid
      integer :: ilen, jlen, il, ij, inb, jnb
      character (len=MAX_LENGTH_VAR_NAME), parameter :: laiVarName = 'AVHRR__AVHRR'
      real*8 :: tempVar(i2-i1+1,j2-ju1+1)
!
! !AUTHOR:
!  
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      start4d(1) = i1  -  i1_gl + 1
      start4d(2) = ju1 - ju1_gl + 1
      start4d(3) = 1
      start4d(4) = 1

      ilen = i2 - i1  + 1
      jlen = j2 - ju1 + 1

      count4d (:) = (/ ilen, jlen, 1, 1 /)

      call Ncop_Rd (ncid, laiMEGAN_InfileName)

      ! Reading the current month data
      start4d(4) = currMonth
      call Ncrd_4d (tempVar, ncid, laiVarName, start4d, count4d)

      do ij = ju1, j2
         jnb = ij - ju1 + 1
         do il = i1, i2
            inb = il - i1 + 1
            isoLaiCurr(il,ij) = tempVar(inb,jnb)
         end do
      end do

      ! Reading the previous month data
      prevMonth = currMonth - 1
      if (prevMonth == 0) prevMonth = 12

      start4d(4) = prevMonth

      call Ncrd_4d (tempVar, ncid, laiVarName, start4d, count4d)

      do ij = ju1, j2
         jnb = ij - ju1 + 1
         do il = i1, i2
            inb = il - i1 + 1
            isoLaiPrev(il,ij) = tempVar(inb,jnb)
         end do
      end do

      ! Reading the next month data
      nextMonth = currMonth + 1
      if (nextMonth == 13) nextMonth = 1

      start4d(4) = nextMonth
      call Ncrd_4d (tempVar, ncid, laiVarName, start4d, count4d)

      do ij = ju1, j2
         jnb = ij - ju1 + 1
         do il = i1, i2
            inb = il - i1 + 1
            isoLaiNext(il,ij) = tempVar(inb,jnb)
         end do
      end do

      call Nccl (ncid)

      return

      end subroutine readIsoLai
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: currLAImonth
!
! !INTERFACE:
!
      function currLAImonth(jday, inMonth, startDay)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: jday         ! days since January 1st
      integer, intent(in) :: inMonth      ! current calendar month
      integer, intent(in) :: startDay(MONTHS_PER_YEAR+1)
!
! !RETURN VALUE:
      integer :: currLAImonth
!
! !DESCRIPTION:
! Determines the current LAI month (different than the calendar month)
!
!EOP
!------------------------------------------------------------------------------
!BOC

     if ( jday < startDay(1) ) THEN
        currLAImonth = 12
     else if ( jday < startDay(inMonth) ) then
        currLAImonth = inMonth-1
     else
        currLAImonth = inMonth
     endif 

     end function currLAImonth
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readMEGANannualEmissFactor
!
! !INTERFACE:
!
      subroutine readMEGANannualEmissFactor &
     &     (aefMboMEGAN_InfileName, aefIsopMEGAN_InfileName, aefMonotMEGAN_InfileName, &
     &      aefOvocMEGAN_InfileName, aefMbo, aefIsop, aefOvoc, aefMonot, &
     &      tdt, mcor, i1, i2, ju1, j2, i1_gl, ju1_gl)
      
      implicit none
      
! !INPUT PARAMETERS:
      integer, intent(in) :: i1, i2, ju1, j2, i1_gl, ju1_gl
      real*8 , intent(in) :: tdt
      real*8 , intent(in) :: mcor(i1:i2,ju1:j2)
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: aefMboMEGAN_InfileName
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: aefIsopMEGAN_InfileName
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: aefOvocMEGAN_InfileName
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: aefMonotMEGAN_InfileName
!     
! !OUTPUT PARAMETERS:
     real*8 , intent(out) :: aefMbo  (i1:i2,ju1:j2)
     real*8 , intent(out) :: aefIsop (i1:i2,ju1:j2)
     real*8 , intent(out) :: aefOvoc (i1:i2,ju1:j2)
     real*8 , intent(out) :: aefMonot(i1:i2,ju1:j2)
!
! !DESCRIPTION:
! Reads Annual Emission Factor for all biogenic VOC species.
! Reference: Guenther et al, 2004.
!
!LOCAL VARIABLES:
      INTEGER                 :: I, J , il, ij, inb, jnb
      REAL*8                  :: DTSRCE
      real*8                  :: factor(i1:i2, ju1:j2)
      integer                 :: ilen, jlen, ncid
      character (len=MAX_LENGTH_VAR_NAME), parameter :: aefMboVarName   = 'BIOGSRCE__AEF'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: aefIsopVarName  = 'BIOGSRCE__AEF'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: aefOvocVarName  = 'BIOGSRCE__AEF'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: aefMonotVarName = 'BIOGSRCE__AEF'
      integer :: count3d (3), start3d(3)
      real*8  :: tempVar(i2-i1+1, j2-ju1+1)
!
! !AUTHOR:
!  
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Emission timestep [min]
      DTSRCE =  tdt / 60.0d0 ! ... get from GMI timing routines

      ! Conversion factor from [ug C/m2/hr] to [kg C/box]
      factor(i1:i2,ju1:j2) = 1.d-9 / 3600.d0 * mcor(i1:i2,ju1:j2) * DTSRCE * 60.d0

      !---------------------------------------------
      ! Read in ISOPRENE Annual Emission Factors
      !---------------------------------------------

      start3d(1) = i1  -  i1_gl + 1
      start3d(2) = ju1 - ju1_gl + 1
      start3d(3) = 1

      ilen = i2 - i1  + 1
      jlen = j2 - ju1 + 1

      count3d (:) = (/ ilen, jlen, 1 /)

      call Ncop_Rd (ncid, aefIsopMEGAN_InfileName)

      call Ncrd_3d (tempVar, ncid, aefIsopVarName, start3d, count3d)

      call Nccl (ncid)
!
      do ij = ju1, j2
         jnb = ij - ju1 + 1
         do il = i1, i2
            inb = il - i1 + 1
            aefIsop(il,ij) = tempVar(inb,jnb)
         end do
      end do

      aefIsop(i1:i2,ju1:j2) = aefIsop(i1:i2,ju1:j2) * factor(i1:i2,ju1:j2)

      !---------------------------------------------
      ! Read in MONOTERPENE Annual Emission Factors
      !---------------------------------------------

      call Ncop_Rd (ncid, aefMonotMEGAN_InfileName)

      call Ncrd_3d (tempVar, ncid, aefMonotVarName, start3d, count3d)

      call Nccl (ncid)
!
      do ij = ju1, j2
         jnb = ij - ju1 + 1
         do il = i1, i2
            inb = il - i1 + 1
            aefMonot(il,ij) = tempVar(inb,jnb)
         end do
      end do

      aefMonot(i1:i2,ju1:j2) = aefMonot(i1:i2,ju1:j2) * factor(i1:i2,ju1:j2)

      !---------------------------------------------
      ! Read in MBO Annual Emission Factors
      !---------------------------------------------

      call Ncop_Rd (ncid, aefMboMEGAN_InfileName)

      call Ncrd_3d (tempVar, ncid, aefMboVarName, start3d, count3d)

      call Nccl (ncid)
!
      do ij = ju1, j2
         jnb = ij - ju1 + 1
         do il = i1, i2
            inb = il - i1 + 1
            aefMbo(il,ij) = tempVar(inb,jnb)
         end do
      end do

      aefMbo(i1:i2,ju1:j2) = aefMbo(i1:i2,ju1:j2) * factor(i1:i2,ju1:j2)

      !---------------------------------------------
      ! Read in other VOC Annual Emission Factors
      !---------------------------------------------

      call Ncop_Rd (ncid, aefOvocMEGAN_InfileName)

      call Ncrd_3d (tempVar, ncid, aefOvocVarName, start3d, count3d)

      call Nccl (ncid)
!
      do ij = ju1, j2
         jnb = ij - ju1 + 1
         do il = i1, i2
            inb = il - i1 + 1
            aefOvoc(il,ij) = tempVar(inb,jnb)
         end do
      end do

      aefOvoc(i1:i2,ju1:j2) = aefOvoc(i1:i2,ju1:j2) * factor(i1:i2,ju1:j2)

      return

      end subroutine readMEGANannualEmissFactor
!EOC
!------------------------------------------------------------------------------

     end module ReadInputMEGAN_mod
