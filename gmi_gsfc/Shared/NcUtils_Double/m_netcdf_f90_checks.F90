!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_f90_checks
!
! !INTERFACE:
!
      module m_netcdf_f90_checks
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  NcF90_does_Udim_Exist
      public  NcF90_does_Var_Exist
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
! !FUNCTION: NcF90_does_Udim_Exist
!
! !INTERFACE:
!
      function NcF90_does_Udim_Exist (ncid)
!
! !USES:
      use netcdf
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid : NetCDF file id to check
      integer, intent (in)   :: ncid
!
! !DESCRIPTION:
!  Checks a given NetCDF file to see if it contains an unlimited dimension.
!
! !RETURN VALUE:
      logical :: NcF90_does_Udim_Exist
!
! !LOCAL VARIABLES:
      integer :: ierr
      integer :: udimid
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
      ierr = NF90_Inquire (ncid, unlimitedDimId = udimid)

      if (ierr == NF90_NOERR) then
         NcF90_does_Udim_Exist = .true.
      else
         NcF90_does_Udim_Exist = .false.
      end if

      return

      end function NcF90_does_Udim_Exist
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !FUNCTION: NcF90_does_Var_Exist
!
! !INTERFACE:
!
      function NcF90_does_Var_Exist (ncid, varname)
!
! !USES:
      use netcdf
!
      implicit none
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
      logical :: NcF90_does_Var_Exist
!
! !LOCAL VARIABLES:
      integer :: ierr
      integer :: varid
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
      ierr = NF90_Inq_Varid (ncid, varname, varid)

      if (ierr == NF90_NOERR) then
         NcF90_does_Var_Exist = .true.
      else
         NcF90_does_Var_Exist = .false.
      end if

      return

      end function NcF90_does_Var_Exist
!EOC
!------------------------------------------------------------------------
end module m_netcdf_f90_checks

