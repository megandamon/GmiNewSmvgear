!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_netcdf_io_get_dimlen
!
! !INTERFACE:
!
      module m_netcdf_io_get_dimlen
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  Ncget_Dimlen
      public  Ncget_Unlim_Dimlen

#     include "GmiParameters.h"
!
! !DESCRIPTION:
!  Provide routines to obtain the length of a given dimension.
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
! !IROUTINE: Ncget_Dimlen
!
! !INTERFACE:
!
      subroutine Ncget_Dimlen (ncid, dim_name, dim_len)
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
!!  dim_name : netCDF dimension name
!!  ncid     : netCDF file id
      character (len=*), intent(in) :: dim_name
      integer,           intent(in) :: ncid
!
! !OUTPUT PARAMETERS:
!!  dim_len: NetCDF dimension length
      integer,           intent(out)   :: dim_len
!
! !DESCRIPTION:
!  Return the length of a given NetCDF dimension.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: dimid
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
      ierr = Nf_Inq_Dimid  (ncid, dim_name, dimid)

      if (ierr /= NF_NOERR) then 
        err_msg = 'In Ncget_Dimlen #1:  ' // Trim (dim_name) // &
                   ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Inq_Dimlen (ncid, dimid, dim_len)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncget_Dimlen #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, dimid, 0, 0.0d0, 0.0d0)
      end if

      return
      end subroutine Ncget_Dimlen
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncget_Unlim_Dimlen
!
! !INTERFACE:
!
      subroutine Ncget_Unlim_Dimlen (ncid, udim_len)
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
!!  ncid     : NetCDF file id
      integer,           intent(in) :: ncid
!
! !OUTPUT PARAMETERS:
!!    udim_len : NetCDF unlimited dimension length
      integer,           intent(out) :: udim_len
!
! !DESCRIPTION:
!  Return the length of the unlimited NetCDF dimension.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: udimid
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
      ierr = Nf_Inq_Unlimdim (ncid, udimid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncget_Unlim_Dimlen #1:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Inq_Dimlen (ncid, udimid, udim_len)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncget_Unlim_Dimlen #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, udimid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncget_Unlim_Dimlen
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_get_dimlen
