!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_checks
!
! !INTERFACE:
!
      module m_netcdf_io_checks
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  Ncdoes_Udim_Exist
      public  Ncdoes_Var_Exist
!
! !DESCRIPTION:
!  Routines to check if a netCDF file contains a specified variable.
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
! !FUNCTION: Ncdoes_Udim_Exist
!
! !INTERFACE:
!
      function Ncdoes_Udim_Exist (ncid)
!
      implicit none
!
#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!!    ncid : NetCDF file id to check
      integer, intent (in)   :: ncid
!
! !DESCRIPTION:
!  Checks a given NetCDF file to see if it contains an unlimited dimension.
!
! !RETURN VALUE:
      logical :: Ncdoes_Udim_Exist
!
! !LOCAL VARIABLES:
      integer :: ierr
      integer :: udimid
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

      if (ierr == NF_NOERR) then
         Ncdoes_Udim_Exist = .true.
      else
         Ncdoes_Udim_Exist = .false.
      end if

      return

      end function Ncdoes_Udim_Exist
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !FUNCTION: Ncdoes_Var_Exist
!
! !INTERFACE:
!
      function Ncdoes_Var_Exist (ncid, varname)
!
      implicit none
!
#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!!    ncid    : NetCDF file id       to check
!!    varname : NetCDF variable name to check
      integer,           intent (in)   :: ncid
      character (len=*), intent (in)   :: varname
!
! !DESCRIPTION:
!  Checks a given NetCDF file to see if a given NetCDF variable 
!  exists in it.
!
! !RETURN VALUE:
      logical :: Ncdoes_Var_Exist
!
! !LOCAL VARIABLES:
      integer :: ierr
      integer :: varid
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
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr == NF_NOERR) then
         Ncdoes_Var_Exist = .true.
      else
         Ncdoes_Var_Exist = .false.
      end if

      return

      end function Ncdoes_Var_Exist
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_checks

