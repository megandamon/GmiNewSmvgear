module GmiFileUnit_mod

      implicit none

      private
      public  :: GetFileUnitNumber
      public  :: InitializeFileUnitNumbers
      public  :: ReleaseFileUnitNumber 

      integer, parameter :: MAX_UNIT_NUMBER = 300
      integer  availUnitNumbers(MAX_UNIT_NUMBER)

      contains

!------------------------------------------------------------------
!------------------------------------------------------------------

      subroutine GetFileUnitNumber (lun, error)

      implicit none

      integer, intent(out) :: lun
      integer, intent(out) :: error

      logical  open

      integer  i

!---------------
! Begin execution
!----------------

      lun = -1

      do i = 7, MAX_UNIT_NUMBER, 1

         if (availUnitNumbers(i) /= 1) then
            inquire (unit = i, OPENED = open)
!                           =====
            if (.not. open) go to 20
!                           =====
         end if

      end do

      if (lun == -1) then
         error = -1
         PRINT*, "GetFileUnitNumber:  No logical unit numbers available"
!        ======
         return
!        ======
      end if

!     ========
 20   continue
!     ========

      error = 0

      lun          = i
      availUnitNumbers(i) = 1

      if ((lun < 7) .or. (lun > MAX_UNIT_NUMBER)) then
         lun   = -1
         error = -1
         PRINT*, "GetFileUnitNumber:  No logical unit numbers available"
!        ======
         return
!        ======
      end if

      return

      end subroutine GetFileUnitNumber 

!------------------------------------------------------------------
!------------------------------------------------------------------

      subroutine ReleaseFileUnitNumber (lun, error)

      implicit none

      integer, intent(in ) :: lun
      integer, intent(out) :: error

!---------------
! Begin execution
!----------------

      if ((lun >= 7) .and. (lun < MAX_UNIT_NUMBER)) then
         error          = 0
         availUnitNumbers(lun) = 0
         close (lun)
      else
         error = -2
         PRINT*, "ReleaseFileUnitNumber:  Invalid input"
      end if

      return

      end subroutine ReleaseFileUnitNumber

!------------------------------------------------------------------
!------------------------------------------------------------------

      subroutine InitializeFileUnitNumbers ( )

      implicit none

      integer  i

!---------------
! Begin execution
!----------------

      do i = 1, MAX_UNIT_NUMBER
        availUnitNumbers(i) = 0
      end do

      availUnitNumbers(5) = -1
      availUnitNumbers(6) = -1

      return

      end subroutine InitializeFileUnitNumbers

end module GmiFileUnit_mod
