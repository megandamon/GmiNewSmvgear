!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_write
!
! !INTERFACE:
!
      module m_netcdf_io_write
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public :: Ncwr_Scal, Ncwr_Scal_Int, Ncwr_2d_Char
      public :: Ncwr_1d, Ncwr_1d_Int
      public :: Ncwr_2d, Ncwr_2d_Int
      public :: Ncwr_3d, Ncwr_3d_Int
      public :: Ncwr_4d, Ncwr_5d, Ncwr_6d
!      public  Ncwr

      interface Ncwr
         module procedure  Ncwr_Scal
         module procedure  Ncwr_Scal_Int
         module procedure  Ncwr_1d
         module procedure  Ncwr_1d_Int
         module procedure  Ncwr_2d
         module procedure  Ncwr_2d_Int
         module procedure  Ncwr_3d
         module procedure  Ncwr_3d_Int
         module procedure  Ncwr_4d
         module procedure  Ncwr_5d
         module procedure  Ncwr_6d
         module procedure  Ncwr_2d_Char
      end interface

#     include "GmiParameters.h"
!
! !DESCRIPTION:
!  Routines for writing variables in a netCDF file.
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
! !IROUTINE: Ncwr_Scal
!
! !INTERFACE:
!
      subroutine Ncwr_Scal (varwr_scal, ncid, varname)
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
!!     ncid       : NetCDF file id to write variable to
!!     varname    : NetCDF variable name
!!     varwr_scal : variable to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      real*8           , intent(in)   :: varwr_scal
!
! !DESCRIPTION:
!  Writes out a NetCDF real scalar variable.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4              :: varwr_scal_tmp
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
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_Scal #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      varwr_scal_tmp = varwr_scal

      ierr = Nf_Put_Var_Real (ncid, varid, varwr_scal_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_Scal #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_Scal
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_Scal_Int
!
! !INTERFACE:
!
      subroutine Ncwr_Scal_Int (varwr_scali, ncid, varname)
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
!!    ncid       : NetCDF file id to write variable to
!!    varname    : NetCDF variable name
!!    varwr_scali : integer variable to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: varwr_scali
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
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_Scal_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Var_Int (ncid, varid, varwr_scali)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_Scal_Int #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_Scal_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_1d
!
! !INTERFACE:
!
      subroutine Ncwr_1d (varwr_1d, ncid, varname, start1d, count1d)
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
!!    ncid     : NetCDF file id to write array output data to
!!    varname  : NetCDF variable name for array
!!    start1d   : vector specifying the index in varwr_1d where
!!               the first of the data values will be written
!!    count1d    : varwr_1d dimension
!!    varwr_1d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start1d(1)
      integer          , intent(in)   :: count1d (1)
      real*8           , intent(in)   :: varwr_1d(count1d(1))
!
! !DESCRIPTION:
!  Writes out a 1D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4              :: varwr_1d_tmp(count1d(1))
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
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_1d #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      varwr_1d_tmp(:) = varwr_1d(:)
      ierr = Nf_Put_Vara_Real (ncid, varid, start1d, count1d, varwr_1d_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_1d #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_1d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_1d_Int
!
! !INTERFACE:
!
      subroutine Ncwr_1d_Int (varwr_1di, ncid, varname, start1d, count1d)
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
!!    ncid     : NetCDF file id to write array output data to
!!    varname  : NetCDF variable name for array
!!    start1d   : vector specifying the index in varwr_1di where
!!               the first of the data values will be written
!!    count1d    : varwr_1di dimension
!!    varwr_1di : intger array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start1d(1)
      integer          , intent(in)   :: count1d (1)
      integer          , intent(in)   :: varwr_1di(count1d(1))
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
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_1d_Int #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Int (ncid, varid, start1d, count1d, varwr_1di)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_1d_Int #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_1d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_2d
!
! !INTERFACE:
!
      subroutine Ncwr_2d (varwr_2d, ncid, varname, start2d, count2d)
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
!!    ncid     : NetCDF file id to write array output data to
!!    varname  : NetCDF variable name for array
!!    start2d   : vector specifying the index in varwr_2d where
!!               the first of the data values will be written
!!    count2d    : varwr_2d dimensions
!!    varwr_2d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start2d(2)
      integer          , intent(in)   :: count2d (2)
      real*8           , intent(in)   :: varwr_2d(count2d(1), count2d(2))
!
! !DESCRIPTION:
!  Writes out a 2D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4              :: varwr_2d_tmp(count2d(1), count2d(2))
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
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      varwr_2d_tmp(:,:) = varwr_2d(:,:)
      ierr = Nf_Put_Vara_Real (ncid, varid, start2d, count2d, varwr_2d_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_2d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_2d_Int
!
! !INTERFACE:
!
      subroutine Ncwr_2d_Int (varwr_2di, ncid, varname, start2d, count2d)
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
!!    ncid     : NetCDF file id to write array output data to
!!    varname  : NetCDF variable name for array
!!    start2d   : vector specifying the index in varwr_2di where
!!               the first of the data values will be written
!!    count2d    : varwr_2di dimensions
!!    varwr_2di : intger array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start2d(2)
      integer          , intent(in)   :: count2d (2)
      integer          , intent(in)   :: varwr_2di(count2d(1), count2d(2))
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
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d_Int #1:  ' // Trim (varname) //  &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Put_Vara_Int (ncid, varid, start2d, count2d, varwr_2di)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d_Int #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_2d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_3d
!
! !INTERFACE:
!
      subroutine Ncwr_3d (varwr_3d, ncid, varname, start3d, count3d)
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
!!    ncid     : NetCDF file id to write array output data to
!!    varname  : NetCDF variable name for array
!!    start3d   : vector specifying the index in varwr_3d where
!!               the first of the data values will be written
!!    count3d    : varwr_3d dimensions
!!    varwr_3d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start3d(3)
      integer          , intent(in)   :: count3d (3)
      real*8           , intent(in)   :: varwr_3d(count3d(1), count3d(2), count3d(3))
!
! !DESCRIPTION:
!  Writes out a 3D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4              :: varwr_3d_tmp(count3d(1), count3d(2), count3d(3))
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
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_3d #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      varwr_3d_tmp(:,:,:) = varwr_3d(:,:,:)
      ierr = Nf_Put_Vara_Real (ncid, varid, start3d, count3d, varwr_3d_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_3d #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_3d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_3d_Int
!
! !INTERFACE:
!
      subroutine Ncwr_3d_Int (varwr_3di, ncid, varname, start3d, count3d)
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
!!    ncid     : NetCDF file id to write array output data to
!!    varname  : NetCDF variable name for array
!!    start3d   : vector specifying the index in varwr_3di where
!!               the first of the data values will be written
!!    count3d    : varwr_3di dimensions
!!    varwr_3di : intger array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start3d(3)
      integer          , intent(in)   :: count3d (3)
      integer          , intent(in)   :: varwr_3di(count3d(1), count3d(2), count3d(3))
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
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_3d_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if


      ierr = Nf_Put_Vara_Int (ncid, varid, start3d, count3d, varwr_3di)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_3d_Int #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_3d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_4d
!
! !INTERFACE:
!
      subroutine Ncwr_4d (varwr_4d, ncid, varname, start4d, count4d)
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
!!    ncid     : NetCDF file id to write array output data to
!!    varname  : NetCDF variable name for array
!!    start4d   : vector specifying the index in varwr_4d where
!!               the first of the data values will be written
!!    count4d    : varwr_4d dimensions
!!    varwr_4d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start4d(4)
      integer          , intent(in)   :: count4d (4)
      real*8           , intent(in)   :: varwr_4d(count4d(1), count4d(2), &
                                                  count4d(3), count4d(4))
!
! !DESCRIPTION:
!  Writes out a 4D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4              :: varwr_4d_tmp(count4d(1), count4d(2), count4d(3), count4d(4))
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
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_4d #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if


      varwr_4d_tmp(:,:,:,:) = varwr_4d(:,:,:,:)
      ierr = Nf_Put_Vara_Real (ncid, varid, start4d, count4d, varwr_4d_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_4d #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_4d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_5d
!
! !INTERFACE:
!
      subroutine Ncwr_5d (varwr_5d, ncid, varname, start5d, count5d)
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
!!    ncid     : NetCDF file id to write array output data to
!!    varname  : NetCDF variable name for array
!!    start5d   : vector specifying the index in varwr_5d where
!!               the first of the data values will be written
!!    count5d    : varwr_5d dimensions
!!    varwr_5d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start5d(5)
      integer          , intent(in)   :: count5d (5)
      real*8           , intent(in)   :: varwr_5d(count5d(1), count5d(2), &
                                                  count5d(3), count5d(4), &
                                                  count5d(5))
!
! !DESCRIPTION:
!  Writes out a 5D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4              :: varwr_5d_tmp(count5d(1), count5d(2), count5d(3), &
                              count5d(4), count5d(5))
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
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_5d #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      varwr_5d_tmp(:,:,:,:,:) = varwr_5d(:,:,:,:,:)
      ierr = Nf_Put_Vara_Real (ncid, varid, start5d, count5d, varwr_5d_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_5d #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_5d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_6d
!
! !INTERFACE:
!
      subroutine Ncwr_6d (varwr_6d, ncid, varname, start6d, count6d)
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
!!    ncid     : NetCDF file id to write array output data to
!!    varname  : NetCDF variable name for array
!!    start6d   : vector specifying the index in varwr_6d where
!!               the first of the data values will be written
!!    count6d    : varwr_6d dimensions
!!    varwr_6d : array to write out
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start6d(6)
      integer          , intent(in)   :: count6d (6)
      real*8           , intent(in)   :: varwr_6d(count6d(1), count6d(2), &
                                                  count6d(3), count6d(4), &
                                                  count6d(5), count6d(6))
!
! !DESCRIPTION:
!  Writes out a 6D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4              :: varwr_6d_tmp(count6d(1), count6d(2), count6d(3), &
                                     count6d(4), count6d(5), count6d(6))
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
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_6d #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      varwr_6d_tmp(:,:,:,:,:,:) = varwr_6d(:,:,:,:,:,:)
      ierr = Nf_Put_Vara_Real (ncid, varid, start6d, count6d, varwr_6d_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_6d #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_6d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncwr_2d_Char
!
! !INTERFACE:
!
      subroutine Ncwr_2d_Char (char_2d, ncid, tvarname, start2d, count2d)
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
!!    ncid     : NetCDF file id to write text to
!!    tvarname : NetCDF variable name for text
!!    start2d   : vector specifying the index in char_2d where
!!               the first of the data values will be written
!!    count2d    : char_2d dimensions
!!    char_2d  : text to write
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: tvarname
      integer          , intent(in)   :: start2d(2)
      integer          , intent(in)   :: count2d (2)
      character (len=1), intent(in)   :: char_2d(count2d(1), count2d(2))
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
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      ierr = Nf_Inq_Varid (ncid, tvarname, tvarid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d_Char #1:  ' // Trim (tvarname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if


      ierr = Nf_Put_Vara_Text (ncid, tvarid, start2d, count2d, char_2d)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncwr_2d_Char #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, tvarid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncwr_2d_Char
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_write

