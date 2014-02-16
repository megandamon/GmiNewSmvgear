!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_f90_close
!
! !INTERFACE:
!
      module m_netcdf_f90_close
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  NcF90_Close
      public  NcF90_Close_Noerr

#     include "GmiParameters.h"
!
! !DESCRIPTION:
!  Routines to close a netCDF file.
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
! !IROUTINE: NcF90_Close
!
! !INTERFACE:
!
      subroutine NcF90_Close (ncid)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid : netCDF file id
      integer, intent (in)   :: ncid
!
! !DESCRIPTION:
!  Close a NetCDF file with file id ncid.
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
      ierr = NF90_Close (ncid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Close:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Close
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Close_Noerr
!
! !INTERFACE:
!
      subroutine NcF90_Close_Noerr (ncid)
!
! !USES:
!
      use netcdf
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid : netCDF file id
      integer, intent (in)   :: ncid
!
! !DESCRIPTION:
!  Close a NetCDF file (with file id ncid) if it is open and
!  suppresses NcF90_Close error messages/exit if it is not.
!
! !LOCAL VARIABLES:
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
      ierr = NF90_Close (ncid)

      return

      end subroutine NcF90_Close_Noerr
!EOC
!------------------------------------------------------------------------
end module m_netcdf_f90_close

