!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiStationFinder_mod
!
! !INTERFACE:
!
      module GmiStationFinder_mod
!
      implicit none
!
#     include "GmiParameters.h"
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: getStationPosition
!
! !DESCRIPTION:
!  Routine to determine if a station exists.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  18April2007 - Initial code.
!
!EOP
!-------------------------------------------------------------------------

   CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getStationPosition
!
! !INTERFACE:

     subroutine getStationPosition(inputSite, lat, lon, fileName)
!
! !USES:
     use GmiPrintError_mod     , only : GmiPrintError
     use GmiASCIIoperations_mod, only : AsciiOpenRead
     use GmiFileUnit_mod       , only : ReleaseFileUnitNumber
!
     implicit none
!
! !INPUT PARAMETERS:
      character(len=*), intent(in )  :: inputSite      ! name of the site
      character(len=*), intent(in )  :: fileName       ! file name containing the list of
                                                       ! sites.
!
! !OUTPUT PARAMETERS:
      real*8          , intent(out)  :: lon            ! Longitude position of the site.
      real*8          , intent(out)  :: lat            ! Latitude  position of the site.

!
! !DESCRIPTION:
!  Given the inputSite, this routine searches the file fileName
!  if the inputSite exists. If the inputSite is in the file,
!  the routine returns its position (lon, lat), otherwise the
!  code will stop.
!
! !LOCAL VARIABLES:
      real*8             :: latLoc, lonLoc
      character(len=MAX_LENGTH_SITE_NAME)  :: siteName
      character(len=80)  :: cline, err_msg
      integer            :: lun, ierr
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

      lon = -999.d0
      lat = -999.d0

      call AsciiOpenRead(lun, fileName)

      read(lun,905) cline
      read(lun,905) cline
      read(lun,905) cline

      readLoop: do
        read(lun, 900, iostat = ierr) siteName, latLoc, lonLoc

        if (ierr < 0) exit readLoop

        if (trim(siteName) == trim(inputSite)) then
           lat = latLoc
           lon = lonLoc
           exit readLoop
        end if
      end do readLoop
 
 905  format(a80)
 900  format(a16, 2f9.2)
      call ReleaseFileUnitNumber(lun, ierr)

      if (lon == -999.d0) then
         err_msg = 'The site does not exist: ' // inputSite
         call GmiPrintError (err_msg, .true., 1, 1, 0, 0, 0.0d0, 0.0d0)
      end if

     return

     end subroutine getStationPosition

!EOC
!-------------------------------------------------------------------------

      end module GmiStationFinder_mod
