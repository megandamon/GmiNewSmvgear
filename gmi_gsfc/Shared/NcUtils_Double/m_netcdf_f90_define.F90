!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_netcdf_f90_define
!
! !INTERFACE:
!
      module m_netcdf_f90_define
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  NcF90_Def_dimension
      public  NcF90_Def_variable
      public  NcF90_Def_var_attributes
      public  NcF90_Def_glob_attributes
      public  NcF90_SetFill
      public  NcF90_End_def

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
! !IROUTINE: NcF90_Def_dimension
!
! !INTERFACE:
!
      subroutine NcF90_Def_dimension(ncid,name,len,id)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
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
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      ierr = NF90_Def_Dim (ncid, name, len, id)

      if (ierr.ne.NF90_NOERR) then
         err_msg = 'NF90_Def_Dim: can not define dimension : '// Trim (name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return
      end subroutine NcF90_Def_dimension
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Def_variable
!
! !INTERFACE:
!
      subroutine NcF90_Def_variable(ncid,name,type,ndims,dims,var_id)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
!
! !INPUT PARAMETERS:
!
!!    ncid   : netCDF file id
!!    name   : name of the variable
!!    type   : type of the variable 
!!             (NF90_FLOAT, NF90_CHAR, NF90_INT, NF90_DOUBLE, NF90_BYTE, NF90_SHORT)
!!    ndims  : number of dimensions of the variable
!!    dims   : netCDF dimension id of the variable

      character (len=*), intent(in) :: name
      integer,           intent(in) :: ncid, ndims
      integer,           intent(in) :: dims(ndims)
      integer,           intent(in) :: type
!
! !OUTPUT PARAMETERS:
!
!!    var_id  : netCDF varid id
      integer,           intent(out) :: var_id
!
! !DESCRIPTION:
!  Define a netCDF variable.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer ::  ierr
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
      ierr = NF90_Def_Var (ncid, name, type, dims, var_id)

      if (ierr.ne.NF90_NOERR) then
         err_msg = 'NF90_Def_Var: can not define variable : '// Trim (name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Def_variable
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Def_var_attributes
!
! !INTERFACE:
!
      subroutine NcF90_Def_var_attributes(ncid,var_id,att_name,att_val)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
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
      ierr = NF90_Put_Att (ncid, var_id, att_name, att_val)

      if (ierr.ne.NF90_NOERR) then
         err_msg = 'NF90_Put_Att_Text: can not define attribute : ' // Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Def_var_attributes
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Def_glob_attributes
!
! !INTERFACE:
!
      subroutine NcF90_Def_glob_attributes(ncid,att_name,att_val)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
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
      ierr = NF90_Put_Att (ncid, NF90_GLOBAL, att_name, att_val)

      if (ierr.ne.NF90_NOERR) then
         err_msg = 'NF90_Put_Att: can not define attribute : ' // Trim (att_name)
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Def_glob_attributes
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_SetFill
!
! !INTERFACE:
!
      subroutine NcF90_SetFill(ncid,ifill,omode)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
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
      ierr = NF90_Set_Fill (ncid, NF90_NOFILL, omode)

      if (ierr.ne.NF90_NOERR) then
         err_msg = 'NF90_Put_Att_Text: Error in omode  '
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_SetFill
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_End_def
!
! !INTERFACE:
!
      subroutine NcF90_End_def(ncid)
!
! !USES:
!
      use m_do_err_out, ONLY: Do_Err_Out
      use netcdf
!
      implicit none
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
      ierr = NF90_Enddef (ncid)

      if (ierr.ne.NF90_NOERR) then
         err_msg = 'NF90_Put_Att_Text: Error in closing global attribute'
         call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_End_def
!EOC
!------------------------------------------------------------------------
end module m_netcdf_f90_define
