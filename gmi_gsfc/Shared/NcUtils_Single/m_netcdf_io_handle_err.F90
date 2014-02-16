!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_handle_err
!
! !INTERFACE:
!
      module m_netcdf_io_handle_err
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  Nchandle_Err

#     include "GmiParameters.h"
!
! !DESCRIPTION:
!  Provide a routine to handle error messages.
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
! !IROUTINE: Nchandle_Err
!
! !INTERFACE:
!
      subroutine Nchandle_Err (ierr)
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
!     ierr : NetCDF error number
      integer, intent (in)   :: ierr
!
! !DESCRIPTION:
!  Handle NetCDF errors. Print out a message and then exit. 
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
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
      err_msg = 'In Nchandle_Err:  ' // Nf_Strerror (ierr)

      call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

      return

      end subroutine Nchandle_Err
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_handle_err

