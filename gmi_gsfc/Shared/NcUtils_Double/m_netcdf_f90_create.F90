!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_f90_create
!
! !INTERFACE:
!
      module m_netcdf_f90_create
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  NcF90_Create_Write
      public  NcF90_Do_Sync

#     include "GmiParameters.h"
!
! !DESCRIPTION:
!  Routines for creating and syncronizing netCDF files.
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
! !IROUTINE: NcF90_Create_Write
!
! !INTERFACE:
!
      subroutine NcF90_Create_Write (ncid, filname)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
!
! !INPUT PARAMETERS:
!     filname : name of NetCDF file to open for writing
      character (len=*), intent(in)   :: filname
!
! !OUTPUT PARAMETERS:
!!    ncid    : opened NetCDF file id
      integer          , intent(out)   ::  ncid
!
! !DESCRIPTION:
!  Create a NetCDF file for writing and do some error checking.
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
!
      ierr = NF90_Create (path=filname, cmode=NF90_64BIT_OFFSET, ncid=ncid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Create_Write, cannot create:  ' // Trim (filname)
        call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Create_Write
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Do_Sync
!
! !INTERFACE:
!
      subroutine NcF90_Do_Sync (ncid)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid : NetCDF file id
      integer, intent(in)   :: ncid
!
! !DESCRIPTION:
!  Synchronize a NetCDF file. 
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
!
      ierr = NF90_Sync (ncid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Do_Sync:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Do_Sync
!EOC
!------------------------------------------------------------------------
end module m_netcdf_f90_create
