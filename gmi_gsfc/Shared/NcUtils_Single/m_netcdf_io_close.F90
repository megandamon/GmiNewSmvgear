!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_close
!
! !INTERFACE:
!
      module m_netcdf_io_close
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  Nccl
      public  Nccl_Noerr

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
! !IROUTINE: Nccl
!
! !INTERFACE:
!
      subroutine Nccl (ncid)
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
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      ierr = Nf_Close (ncid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Nccl:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Nccl
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nccl_Noerr
!
! !INTERFACE:
!
      subroutine Nccl_Noerr (ncid)
!
      implicit none
!
#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!!    ncid : netCDF file id
      integer, intent (in)   :: ncid
!
! !DESCRIPTION:
!  Close a NetCDF file (with file id ncid) if it is open and
!  suppresses Ncclos error messages/exit if it is not.
!
! !LOCAL VARIABLES:
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
      ierr = Nf_Close (ncid)

      return

      end subroutine Nccl_Noerr
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_close

