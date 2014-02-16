!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_create
!
! !INTERFACE:
!
      module m_netcdf_io_create
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  Nccr_Wr
      public  Ncdo_Sync

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
! !IROUTINE: Nccr_Wr
!
! !INTERFACE:
!
      subroutine Nccr_Wr (ncid, filname)
!
! !USES:
!
      use GmiPrintError_mod, only : GmiPrintError
!
      implicit none
!
#     include "netcdf.inc"
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
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      ierr = Nf_Create (filname, NF_CLOBBER, ncid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Nccr_Wr, cannot create:  ' // Trim (filname)
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Nccr_Wr
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncdo_Sync
!
! !INTERFACE:
!
      subroutine Ncdo_Sync (ncid)
!
! !USES:
!
      use GmiPrintError_mod, only : GmiPrintError
!
      implicit none
!
#     include "netcdf.inc"
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
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      ierr = Nf_Sync (ncid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncdo_Sync:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncdo_Sync
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_create
