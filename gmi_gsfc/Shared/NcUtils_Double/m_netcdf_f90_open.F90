!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_f90_open
!
! !INTERFACE:
!
      module m_netcdf_f90_open
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  NcF90_Open_Read
      public  NcF90_Open_Write

#     include "GmiParameters.h"
!
! !DESCRIPTION:
!  Routines to open a netCDF file.
!
! !AUTHOR: 
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Open_Read
!
! !INTERFACE:
!
      subroutine NcF90_Open_Read (ncid, filname)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
!
! !INPUT PARAMETERS:
!!    filname : name of NetCDF file to open for reading
      character (len=*), intent (in)    :: filname
!
! !OUTPUT PARAMETERS:
!!    ncid    : opened NetCDF file id
      integer          , intent (out)   :: ncid
!
! !DESCRIPTION:
!  Open a NetCDF file for reading and do some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
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
      ierr = NF90_Open (path=filname, mode=NF90_NOWRITE, ncid=ncid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Open_Read, cannot open:  ' // Trim (filname)
        call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Open_Read
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Open_Write
!
! !INTERFACE:
!
      subroutine NcF90_Open_Write (ncid, filname)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
!
! !INPUT PARAMETERS:
!!    filname : name of NetCDF file to open for reading
      character (len=*), intent (in)    :: filname
!
! !OUTPUT PARAMETERS:
!!    ncid    : opened NetCDF file id
      integer          , intent (out)   :: ncid
!
! !DESCRIPTION:
!  Open a NetCDF file for reading/writing and do some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
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
      ierr = NF90_Open (path=filname, mode=NF90_WRITE, ncid=ncid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Open_Write, cannot open:  ' // Trim (filname)
        call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Open_Write
!EOC
!------------------------------------------------------------------------
end module m_netcdf_f90_open

