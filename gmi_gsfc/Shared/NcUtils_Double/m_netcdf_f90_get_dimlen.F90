!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_netcdf_f90_get_dimlen
!
! !INTERFACE:
!
      module m_netcdf_f90_get_dimlen
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  NcF90_Get_Dimlen
      public  NcF90_Get_Unlim_Dimlen

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
! !IROUTINE: NcF90_Get_Dimlen
!
! !INTERFACE:
!
      subroutine NcF90_Get_Dimlen (ncid, dim_name, dim_len)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
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
!  Checks if a netCDF dimension exists and returns its length.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: dimid
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
      ierr = NF90_Inq_Dimid  (ncid, dim_name, dimid)

      if (ierr /= NF90_NOERR) then 
        err_msg = 'In NcF90_Get_Dimlen #1:  ' // Trim (dim_name) // &
                   ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Inquire_Dimension (ncid, dimid, len=dim_len)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Get_Dimlen #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, dimid, 0, 0.0d0, 0.0d0)
      end if

      return
      end subroutine NcF90_Get_Dimlen
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Get_Unlim_Dimlen
!
! !INTERFACE:
!
      subroutine NcF90_Get_Unlim_Dimlen (ncid, dim_len)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
!
! !INPUT PARAMETERS:
!!  dim_name : netCDF dimension name
!!  ncid     : netCDF file id
      integer,           intent(in) :: ncid
!
! !OUTPUT PARAMETERS:
!!  dim_len: NetCDF dimension length
      integer,           intent(out)   :: dim_len
!
! !DESCRIPTION:
!  Checks if a netCDF dimension exists and returns its length.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: dimid
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
      ierr = NF90_Inquire (ncid, unlimitedDimId = dimid)

      if (ierr /= NF90_NOERR) then 
        err_msg = 'In NcF90_Get_Unlim_Dimlen #1: no unlimited dim ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Inquire_Dimension (ncid, dimid, len=dim_len)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Get_Unlim_Dimlen #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, dimid, 0, 0.0d0, 0.0d0)
      end if

      return
      end subroutine NcF90_Get_Unlim_Dimlen
!EOC
!-------------------------------------------------------------------------

end module m_netcdf_f90_get_dimlen
