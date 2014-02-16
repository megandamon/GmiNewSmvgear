!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_open
!
! !INTERFACE:
!
      module m_netcdf_io_open
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  Ncop_Rd
      public  Ncop_Wr

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
! !IROUTINE: Ncop_Rd
!
! !INTERFACE:
!
      subroutine Ncop_Rd (ncid, filname)
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
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      ierr = Nf_Open (filname, NF_NOWRITE, ncid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncop_Rd, cannot open:  ' // Trim (filname)
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncop_Rd
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncop_Wr
!
! !INTERFACE:
!
      subroutine Ncop_Wr (ncid, filname)
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
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      ierr = Nf_Open (filname, NF_WRITE, ncid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncop_Rd, cannot open:  ' // Trim (filname)
        call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncop_Wr
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_open

