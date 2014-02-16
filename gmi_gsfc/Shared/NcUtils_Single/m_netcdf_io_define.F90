!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_netcdf_io_define
!
! !INTERFACE:
!
      module m_netcdf_io_define
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public  NcDef_dimension
      public  NcDef_variable
      public  NcDef_var_attributes
      public  NcDef_glob_attributes
      public  NcSetFill
      public  NcEnd_def

#     include "GmiParameters.h"
!
! !DESCRIPTION:
!  Provide netcdf utilities to define dimension, variable and attributes
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
! !IROUTINE: NcDef_dimension
!
! !INTERFACE:
!
      subroutine NcDef_dimension(ncid,name,len,id)
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
!!    ncid  : netCDF file id
!!    name  : dimension name
!!    len   : dimension number
      character (len=*), intent(in) :: name
      integer,           intent(in) :: ncid, len
!
! !OUTPUT PARAMETERS:
!!    id    : dimension id
      integer,           intent(out) :: id
!
! !DESCRIPTION:
!  Define dimension.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer :: ierr
!
! !AUTHOR: 
!  Jules Kouatchou and Maharaj Bhat
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      ierr = Nf_Def_Dim (ncid, name, len, id)

      if (ierr.ne.NF_NOERR) then
         err_msg = 'Nf_Def_Dim: can not define dimension : '// Trim (name)
         call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return
      end subroutine NcDef_dimension
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_variable
!
! !INTERFACE:
!
      subroutine NcDef_variable(ncid,name,type,ndims,dims,var_id)
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
!
!!    ncid   : netCDF file id
!!    name   : name of the variable
!!    type   : type of the variable 
!!             (NF_FLOAT, NF_CHAR, NF_INT, NF_DOUBLE, NF_BYTE, NF_SHORT)
!!    ndims  : number of dimensions of the variable
!!    dims   : netCDF dimension id of the variable
!!    varid  : netCDF varid id

      character (len=*), intent(in) :: name
      integer,           intent(in) :: ncid, ndims, var_id
      integer,           intent(in) :: dims(ndims)
      integer,           intent(in) :: type
!
! !DESCRIPTION:
!  Define a netCDF variable.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer ::  ierr
!
! !AUTHOR: 
!  Jules Kouatchou and Maharaj Bhat
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      ierr = Nf_Def_Var (ncid, name, type, ndims, dims, var_id)

      if (ierr.ne.NF_NOERR) then
         err_msg = 'Nf_Def_Var: can not define variable : '// Trim (name)
         call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_variable
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_var_attributes
!
! !INTERFACE:
!
      subroutine NcDef_var_attributes(ncid,var_id,att_name,att_val)
!
! !USES:
!
      use GmiPrintError_mod, only : GmiPrintError
!
      implicit none
#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!!    ncid    : netCDF file id
!!    var_id  : netCDF variable id
!!    att_name: attribute name
!!    att_val : attribute value
      character (len=*), intent(in) :: att_name, att_val
      integer,           intent(in) :: ncid, var_id
!
! !DESCRIPTION:
!  Define netCDF attributes.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             ::  mylen, ierr
!
! !AUTHOR: 
!  Jules Kouatchou and Maharaj Bhat
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      mylen = len(att_val)
      ierr = Nf_Put_Att_Text (ncid, var_id, att_name, mylen, att_val)

      if (ierr.ne.NF_NOERR) then
         err_msg = 'Nf_Put_Att_Text: can not define attribute : ' // Trim (att_name)
         call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_var_attributes
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcDef_glob_attributes
!
! !INTERFACE:
!
      subroutine NcDef_glob_attributes(ncid,att_name,att_val)
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
!!    ncid    : netCDF file id
!!    att_name: attribute name
!!    att_val : attribute value
!
      character (len=*), intent(in) :: att_name, att_val
      integer,           intent(in) :: ncid
!
! !DESCRIPTION:
!  Define attributes
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             ::  mylen, ierr
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
      mylen = len(att_val)
      ierr = Nf_Put_Att_Text (ncid, NF_GLOBAL, att_name, mylen, att_val)

      if (ierr.ne.NF_NOERR) then
         err_msg = 'Nf_Put_Att_Text: can not define attribute : ' // Trim (att_name)
         call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcDef_glob_attributes
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcSetFill
!
! !INTERFACE:
!
      subroutine NcSetFill(ncid,ifill,omode)
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
!
      integer,           intent(in) :: ncid, ifill
!
! !OUTPUT PARAMETERS:
!
      integer,           intent(out) :: omode
!
! !DESCRIPTION:
!  Set fill method
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             ::  mylen, ierr
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
      ierr = Nf_Set_Fill (ncid, NF_NOFILL, omode)

      if (ierr.ne.NF_NOERR) then
         err_msg = 'Nf_Put_Att_Text: Error in omode  '
         call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcSetFill
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcEnd_def
!
! !INTERFACE:
!
      subroutine NcEnd_def(ncid)
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
!
      integer,           intent(in) :: ncid
!
! !DESCRIPTION:
!  Ends definitions of variables and their attributes.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             ::  ierr
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
      ierr = Nf_Enddef (ncid)

      if (ierr.ne.NF_NOERR) then
         err_msg = 'Nf_Put_Att_Text: Error in closing global attribute'
         call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcEnd_def
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_define
