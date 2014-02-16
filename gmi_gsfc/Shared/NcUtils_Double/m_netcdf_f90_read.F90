!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_f90_read
!
! !INTERFACE:
!
      module m_netcdf_f90_read
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private  
      public  NcF90_Read

      interface NcF90_Read
         module procedure  NcF90_Read_Scal
         module procedure  NcF90_Read_Scal_Int
         module procedure  NcF90_Read_1d
         module procedure  NcF90_Read_1d_Int
         module procedure  NcF90_Read_2d
         module procedure  NcF90_Read_2d_Int
         module procedure  NcF90_Read_3d
         module procedure  NcF90_Read_3d_Int
         module procedure  NcF90_Read_4d
         module procedure  NcF90_Read_5d
         module procedure  NcF90_Read_1d_Char
         module procedure  NcF90_Read_2d_Char
      end interface

#     include "GmiParameters.h"
!
! !DESCRIPTION:
!  Routines for reading variables in a netCDF file.
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
! !IROUTINE: NcF90_Read_Scal
!
! !INTERFACE:
!
      subroutine NcF90_Read_Scal (varRead_scal, ncid, varname)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid       : NetCDF file id to read variable from
!!    varname    : NetCDF variable name
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
!
! !OUTPUT PARAMETERS:
!!    varRead_scal : variable to fill
      real*8           , intent(out)  :: varRead_scal
!
! !DESCRIPTION:
!  Reads in a NetCDF scalar variable.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
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
      ierr = NF90_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_Scal #1:  ' // Trim (varname) // &
                 ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Get_Var (ncid, varid, varRead_scal)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_Scal #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Read_Scal
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Read_Scal_Int
!
! !INTERFACE:
!
      subroutine NcF90_Read_Scal_Int (varRead_scali, ncid, varname)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid       : NetCDF file id to read variable from
!!    varname    : NetCDF variable name
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
!
! !OUTPUT PARAMETERS:
!!    varRead_scali : integer variable to fill
      integer          , intent(out)  :: varRead_scali
!
! !DESCRIPTION:
!  Reads in a NetCDF integer scalar variable.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
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
      ierr = NF90_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_Scal_Int #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Get_Var (ncid, varid, varRead_scali)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_Scal_Int #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Read_Scal_Int
!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Read_1d
!
! !INTERFACE:
!
      subroutine NcF90_Read_1d (varRead_1d, ncid, varname, strt1d, cnt1d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    strt1d   : vector specifying the index in varRead_1d where 
!!               the first of the data values will be read 
!!    cnt1d    : varRead_1d dimension
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt1d(1)
      integer          , intent(in)   :: cnt1d (1)
!
! !OUTPUT PARAMETERS:
!!    varRead_1d : array to fill
      real*8           , intent(out)  :: varRead_1d(cnt1d(1))
!
! !DESCRIPTION:
!  Reads in a 1D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
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
      ierr = NF90_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_1d #1:  ' // Trim (varname) // &
                   ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr =  NF90_Get_Var (ncid, varid, varRead_1d, start=strt1d, count=cnt1d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_1d #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Read_1d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Read_1d_Int
!
! !INTERFACE:
!
      subroutine NcF90_Read_1d_Int (varRead_1di, ncid, varname, strt1d, cnt1d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    strt1d   : vector specifying the index in varRead_1di where
!!               the first of the data values will be read
!!    cnt1d    : varRead_1di dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt1d(1)
      integer          , intent(in)   :: cnt1d (1)
!
! !OUTPUT PARAMETERS:
!!    varRead_1di : intger array to fill
      integer          , intent(out)  :: varRead_1di(cnt1d(1))
!
! !DESCRIPTION:
!  Reads in a 1D NetCDF integer array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
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
      ierr = NF90_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_1d_Int #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr =  NF90_Get_Var (ncid, varid, varRead_1di, start=strt1d, count=cnt1d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_1d_Int #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Read_1d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Read_2d
!
! !INTERFACE:
!
      subroutine NcF90_Read_2d (varRead_2d, ncid, varname, strt2d, cnt2d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    strt2d   : vector specifying the index in varRead_2d where
!!               the first of the data values will be read
!!    cnt2d    : varRead_2d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
!
! !OUTPUT PARAMETERS:
!!    varRead_2d : array to fill
      real*8           , intent(out)  :: varRead_2d(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION:
!  Reads in a 2D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
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
      ierr = NF90_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_2d #1:  ' // Trim (varname) // & 
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr =  NF90_Get_Var (ncid, varid, varRead_2d, start=strt2d, count=cnt2d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_2d #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Read_2d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Read_2d_Int
!
! !INTERFACE:
!
      subroutine NcF90_Read_2d_Int (varRead_2di, ncid, varname, strt2d, cnt2d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    strt2d   : vector specifying the index in varRead_2di where
!!               the first of the data values will be read
!!    cnt2d    : varRead_2di dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
!
! !OUTPUT PARAMETERS:
!!    varRead_2di : intger array to fill
      integer          , intent(out)  :: varRead_2di(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION:
!  Reads in a 2D NetCDF integer array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
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
      ierr = NF90_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_2d_Int #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr =  NF90_Get_Var (ncid, varid, varRead_2di, start=strt2d, count=cnt2d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_2d_Int #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Read_2d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Read_3d
!
! !INTERFACE:
!
      subroutine NcF90_Read_3d (varRead_3d, ncid, varname, strt3d, cnt3d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    strt3d   : vector specifying the index in varRead_3d where
!!               the first of the data values will be read
!!    cnt3d    : varRead_3d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt3d(3)
      integer          , intent(in)   :: cnt3d (3)
!
! !OUTPUT PARAMETERS:
!!    varRead_3d : array to fill
      real*8           , intent(out)  :: varRead_3d(cnt3d(1), cnt3d(2), &
                                                  cnt3d(3))
!
! !DESCRIPTION:
!  Reads in a 3D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
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
      ierr = NF90_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_3d #1:  ' // Trim (varname) // &
                 ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr =  NF90_Get_Var (ncid, varid, varRead_3d, start=strt3d, count=cnt3d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_3d #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Read_3d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Read_3d_Int
!
! !INTERFACE:
!
      subroutine NcF90_Read_3d_Int (varRead_3di, ncid, varname, strt3d, cnt3d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    strt3d   : vector specifying the index in varRead_3di where
!!               the first of the data values will be read
!!    cnt3d    : varRead_3di dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt3d(3)
      integer          , intent(in)   :: cnt3d (3)
!
! !OUTPUT PARAMETERS:
!!    varRead_3di : intger array to fill
      integer          , intent(out)  :: varRead_3di(cnt3d(1), cnt3d(2), &
                                                   cnt3d(3))
!
! !DESCRIPTION:
!  Reads in a 3D NetCDF integer array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
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
      ierr = NF90_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_3d_Int #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr =  NF90_Get_Var (ncid, varid, varRead_3di, start=strt3d, count=cnt3d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_3d_Int #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Read_3d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Read_4d
!
! !INTERFACE:
!
      subroutine NcF90_Read_4d (varRead_4d, ncid, varname, strt4d, cnt4d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    strt4d   : vector specifying the index in varRead_4d where
!!               the first of the data values will be read
!!    cnt4d    : varRead_4d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt4d(4)
      integer          , intent(in)   :: cnt4d (4)
!
! !OUTPUT PARAMETERS:
!!    varRead_4d : array to fill
      real*8           , intent(out)  :: varRead_4d(cnt4d(1), cnt4d(2), &
                                                  cnt4d(3), cnt4d(4))
!
! !DESCRIPTION:
!  Reads in a 4D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
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
      ierr = NF90_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_4d #1:  ' // Trim (varname) // &
                    ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr =  NF90_Get_Var (ncid, varid, varRead_4d, start=strt4d, count=cnt4d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_4d #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Read_4d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Read_5d
!
! !INTERFACE:
!
      subroutine NcF90_Read_5d (varRead_5d, ncid, varname, strt5d, cnt5d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    strt5d   : vector specifying the index in varRead_5d where
!!               the first of the data values will be read
!!    cnt5d    : varRead_5d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt5d(5)
      integer          , intent(in)   :: cnt5d (5)
!
! !OUTPUT PARAMETERS:
!!    varRead_5d : array to fill
      real*8         , intent(out)  :: varRead_5d(cnt5d(1), cnt5d(2), &
                                                cnt5d(3), cnt5d(4), &
                                                cnt5d(5))
!
! !DESCRIPTION:
!  Reads in a 5D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
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
      ierr = NF90_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_5d #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr =  NF90_Get_Var (ncid, varid, varRead_5d, start=strt5d, count=cnt5d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_5d #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Read_5d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Read_1d_Char
!
! !INTERFACE:
!
      subroutine NcF90_Read_1d_Char (varRead_1dc, ncid, varname, strt1d, cnt1d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    strt1d   : vector specifying the index in varRead_1di where 
!!               the first of the data values will be read 
!!    cnt1d    : varRead_1dc dimension
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt1d(1)
      integer          , intent(in)   :: cnt1d (1)
!
! !OUTPUT PARAMETERS:
!!    varRead_1dc : intger array to fill
      character (len=1), intent(out)  :: varRead_1dc(cnt1d(1))
!
! !DESCRIPTION:
!  Reads in a 1D NetCDF character array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
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
      ierr = NF90_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_1d_Char #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr =  NF90_Get_Var (ncid, varid, varRead_1dc, start=strt1d, count=cnt1d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_1d_Char #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Read_1d_Char
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Read_2d_Char
!
! !INTERFACE:
!
      subroutine NcF90_Read_2d_Char (varRead_2dc, ncid, varname, strt2d, cnt2d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    strt2d   : vector specifying the index in varRead_2dc where
!!               the first of the data values will be read
!!    cnt2d    : varRead_2dc dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
!
! !OUTPUT PARAMETERS:
!!    varRead_2dc : charcter array to fill
      character        , intent(out)  :: varRead_2dc(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION:
!  Reads in a 2D NetCDF character array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
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
      ierr = NF90_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_2d_Char #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr =  NF90_Get_Var (ncid, varid, varRead_2dc, start=strt2d, count=cnt2d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Read_2d_Char #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Read_2d_Char
!EOC
!------------------------------------------------------------------------
end module m_netcdf_f90_read

