!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_f90_handle_err
!
! !INTERFACE:
!
      module m_netcdf_f90_handle_err
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  NcF90_handle_Err

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
! !IROUTINE: NcF90_handle_Err
!
! !INTERFACE:
!
      subroutine NcF90_handle_Err (ierr)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
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
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      err_msg = 'In NcF90_handle_Err:  ' // NF90_Strerror (ierr)

      call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

      return

      end subroutine NcF90_handle_Err
!EOC
!------------------------------------------------------------------------
end module m_netcdf_f90_handle_err

