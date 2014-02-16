!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_f90_Write
!
! !INTERFACE:
!
      module m_netcdf_f90_Write
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  NcF90_Write

      interface NcF90_Write
         module procedure  NcF90_Write_Scal
         module procedure  NcF90_Write_Scal_Int
         module procedure  NcF90_Write_1d
         module procedure  NcF90_Write_1d_Int
         module procedure  NcF90_Write_2d
         module procedure  NcF90_Write_2d_Int
         module procedure  NcF90_Write_3d
         module procedure  NcF90_Write_3d_Int
         module procedure  NcF90_Write_4d
         module procedure  NcF90_Write_5d
         module procedure  NcF90_Write_6d
         module procedure  NcF90_Write_2d_Char
      end interface

#     include "GmiParameters.h"
!
! !DESCRIPTION:
!  Routines for Writing variables in a netCDF file.
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
! !IROUTINE: NcF90_Write_Scal
!
! !INTERFACE:
!
      subroutine NcF90_Write_Scal (varWrite_scal, ncid, varname)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!     ncid       : NetCDF file id to Write variable to
!!     varname    : NetCDF variable name
!!     varWrite_scal : variable to Write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      real*8           , intent(in)   :: varWrite_scal
!
! !DESCRIPTION:
!  Writes out a NetCDF real scalar variable.
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
        err_msg = 'In NcF90_Write_Scal #1:  ' // Trim (varname) // &
                 ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Put_Var (ncid, varid, varWrite_scal)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Write_Scal #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Write_Scal
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Write_Scal_Int
!
! !INTERFACE:
!
      subroutine NcF90_Write_Scal_Int (varWrite_scali, ncid, varname)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid       : NetCDF file id to Write variable to
!!    varname    : NetCDF variable name
!!    varWrite_scali : integer variable to Write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: varWrite_scali
!
! !DESCRIPTION:
!  Writes out a NetCDF integer scalar variable.
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
        err_msg = 'In NcF90_Write_Scal_Int #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Put_Var (ncid, varid, varWrite_scali)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Write_Scal_Int #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Write_Scal_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Write_1d
!
! !INTERFACE:
!
      subroutine NcF90_Write_1d (varWrite_1d, ncid, varname, strt1d, cnt1d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to Write array output data to
!!    varname  : NetCDF variable name for array
!!    strt1d   : vector specifying the index in varWrite_1d where
!!               the first of the data values will be written
!!    cnt1d    : varWrite_1d dimension
!!    varWrite_1d : array to Write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt1d(1)
      integer          , intent(in)   :: cnt1d (1)
      real*8           , intent(in)   :: varWrite_1d(cnt1d(1))
!
! !DESCRIPTION:
!  Writes out a 1D NetCDF real array and does some error checking.
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
        err_msg = 'In NcF90_Write_1d #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Put_Var (ncid, varid, varWrite_1d, start=strt1d, count=cnt1d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Write_1d #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Write_1d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Write_1d_Int
!
! !INTERFACE:
!
      subroutine NcF90_Write_1d_Int (varWrite_1di, ncid, varname, strt1d, cnt1d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to Write array output data to
!!    varname  : NetCDF variable name for array
!!    strt1d   : vector specifying the index in varWrite_1di where
!!               the first of the data values will be written
!!    cnt1d    : varWrite_1di dimension
!!    varWrite_1di : intger array to Write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt1d(1)
      integer          , intent(in)   :: cnt1d (1)
      integer          , intent(in)   :: varWrite_1di(cnt1d(1))
!
! !DESCRIPTION:
!  Writes out a 1D NetCDF integer array and does some error checking.
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
        err_msg = 'In NcF90_Write_1d_Int #1:  ' // Trim (varname) // &
                 ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Put_Var (ncid, varid, varWrite_1di, start=strt1d, count=cnt1d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Write_1d_Int #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Write_1d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Write_2d
!
! !INTERFACE:
!
      subroutine NcF90_Write_2d (varWrite_2d, ncid, varname, strt2d, cnt2d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to Write array output data to
!!    varname  : NetCDF variable name for array
!!    strt2d   : vector specifying the index in varWrite_2d where
!!               the first of the data values will be written
!!    cnt2d    : varWrite_2d dimensions
!!    varWrite_2d : array to Write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
      real*8           , intent(in)   :: varWrite_2d(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION:
!  Writes out a 2D NetCDF real array and does some error checking.
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
        err_msg = 'In NcF90_Write_2d #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Put_Var (ncid, varid, varWrite_2d, start=strt2d, count=cnt2d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Write_2d #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Write_2d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Write_2d_Int
!
! !INTERFACE:
!
      subroutine NcF90_Write_2d_Int (varWrite_2di, ncid, varname, strt2d, cnt2d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to Write array output data to
!!    varname  : NetCDF variable name for array
!!    strt2d   : vector specifying the index in varWrite_2di where
!!               the first of the data values will be written
!!    cnt2d    : varWrite_2di dimensions
!!    varWrite_2di : intger array to Write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
      integer          , intent(in)   :: varWrite_2di(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION:
!  Writes out a 2D NetCDF integer array and does some error checking.
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
        err_msg = 'In NcF90_Write_2d_Int #1:  ' // Trim (varname) //  &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Put_Var (ncid, varid, varWrite_2di, start=strt2d, count=cnt2d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Write_2d_Int #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Write_2d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Write_3d
!
! !INTERFACE:
!
      subroutine NcF90_Write_3d (varWrite_3d, ncid, varname, strt3d, cnt3d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to Write array output data to
!!    varname  : NetCDF variable name for array
!!    strt3d   : vector specifying the index in varWrite_3d where
!!               the first of the data values will be written
!!    cnt3d    : varWrite_3d dimensions
!!    varWrite_3d : array to Write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt3d(3)
      integer          , intent(in)   :: cnt3d (3)
      real*8           , intent(in)   :: varWrite_3d(cnt3d(1), cnt3d(2), cnt3d(3))
!
! !DESCRIPTION:
!  Writes out a 3D NetCDF real array and does some error checking.
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
        err_msg = 'In NcF90_Write_3d #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Put_Var (ncid, varid, varWrite_3d, start=strt3d, count=cnt3d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Write_3d #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Write_3d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Write_3d_Int
!
! !INTERFACE:
!
      subroutine NcF90_Write_3d_Int (varWrite_3di, ncid, varname, strt3d, cnt3d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to Write array output data to
!!    varname  : NetCDF variable name for array
!!    strt3d   : vector specifying the index in varWrite_3di where
!!               the first of the data values will be written
!!    cnt3d    : varWrite_3di dimensions
!!    varWrite_3di : intger array to Write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt3d(3)
      integer          , intent(in)   :: cnt3d (3)
      integer          , intent(in)   :: varWrite_3di(cnt3d(1), cnt3d(2), cnt3d(3))
!
! !DESCRIPTION:
!  Writes out a 3D NetCDF integer array and does some error checking.
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
        err_msg = 'In NcF90_Write_3d_Int #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Put_Var (ncid, varid, varWrite_3di, start=strt3d, count=cnt3d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Write_3d_Int #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Write_3d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Write_4d
!
! !INTERFACE:
!
      subroutine NcF90_Write_4d (varWrite_4d, ncid, varname, strt4d, cnt4d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to Write array output data to
!!    varname  : NetCDF variable name for array
!!    strt4d   : vector specifying the index in varWrite_4d where
!!               the first of the data values will be written
!!    cnt4d    : varWrite_4d dimensions
!!    varWrite_4d : array to Write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt4d(4)
      integer          , intent(in)   :: cnt4d (4)
      real*8           , intent(in)   :: varWrite_4d(cnt4d(1), cnt4d(2), &
                                                  cnt4d(3), cnt4d(4))
!
! !DESCRIPTION:
!  Writes out a 4D NetCDF real array and does some error checking.
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
        err_msg = 'In NcF90_Write_4d #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Put_Var (ncid, varid, varWrite_4d, start=strt4d, count=cnt4d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Write_4d #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Write_4d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Write_5d
!
! !INTERFACE:
!
      subroutine NcF90_Write_5d (varWrite_5d, ncid, varname, strt5d, cnt5d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to Write array output data to
!!    varname  : NetCDF variable name for array
!!    strt5d   : vector specifying the index in varWrite_5d where
!!               the first of the data values will be written
!!    cnt5d    : varWrite_5d dimensions
!!    varWrite_5d : array to Write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt5d(5)
      integer          , intent(in)   :: cnt5d (5)
      real*8           , intent(in)   :: varWrite_5d(cnt5d(1), cnt5d(2), &
                                                  cnt5d(3), cnt5d(4), &
                                                  cnt5d(5))
!
! !DESCRIPTION:
!  Writes out a 5D NetCDF real array and does some error checking.
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
        err_msg = 'In NcF90_Write_5d #1:  ' // Trim (varname) // &
                 ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Put_Var (ncid, varid, varWrite_5d, start=strt5d, count=cnt5d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Write_5d #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Write_5d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Write_6d
!
! !INTERFACE:
!
      subroutine NcF90_Write_6d (varWrite_6d, ncid, varname, strt6d, cnt6d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to Write array output data to
!!    varname  : NetCDF variable name for array
!!    strt6d   : vector specifying the index in varWrite_6d where
!!               the first of the data values will be written
!!    cnt6d    : varWrite_6d dimensions
!!    varWrite_6d : array to Write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: strt6d(6)
      integer          , intent(in)   :: cnt6d (6)
      real*8           , intent(in)   :: varWrite_6d(cnt6d(1), cnt6d(2), &
                                                  cnt6d(3), cnt6d(4), &
                                                  cnt6d(5), cnt6d(6))
!
! !DESCRIPTION:
!  Writes out a 6D NetCDF real array and does some error checking.
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
        err_msg = 'In NcF90_Write_6d #1:  ' // Trim (varname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Put_Var (ncid, varid, varWrite_6d, start=strt6d, count=cnt6d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Write_6d #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Write_6d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcF90_Write_2d_Char
!
! !INTERFACE:
!
      subroutine NcF90_Write_2d_Char (char_2d, ncid, tvarname, strt2d, cnt2d)
!
! !USES:
!
      use netcdf
      use m_do_err_out, ONLY: Do_Err_Out
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to Write text to
!!    tvarname : NetCDF variable name for text
!!    strt2d   : vector specifying the index in char_2d where
!!               the first of the data values will be written
!!    cnt2d    : char_2d dimensions
!!    char_2d  : text to Write
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: tvarname
      integer          , intent(in)   :: strt2d(2)
      integer          , intent(in)   :: cnt2d (2)
      character (len=1), intent(in)   :: char_2d(cnt2d(1), cnt2d(2))
!
! !DESCRIPTION:
!  Writes out a 2D NetCDF character array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: tvarid
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
      ierr = NF90_Inq_Varid (ncid, tvarname, tvarid)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Write_2d_Char #1:  ' // Trim (tvarname) // &
                  ', ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = NF90_Put_Var (ncid, tvarid, char_2d, start=strt2d, count=cnt2d)

      if (ierr /= NF90_NOERR) then
        err_msg = 'In NcF90_Write_2d_Char #2:  ' // NF90_Strerror (ierr)
        call Do_Err_Out (err_msg, .true., 2, ncid, tvarid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine NcF90_Write_2d_Char
!EOC
!------------------------------------------------------------------------
end module m_netcdf_f90_Write

