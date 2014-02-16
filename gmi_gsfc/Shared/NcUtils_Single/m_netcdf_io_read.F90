!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_netcdf_io_read
!
! !INTERFACE:
!
      module m_netcdf_io_read
!
! !USES:
!
      use GmiPrintError_mod, only : GmiPrintError
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public ::  Ncrd_Scal
      public ::  Ncrd_Scal_Int
      public ::  Ncrd_1d
      public ::  Ncrd_1d_Int
      public ::  Ncrd_2d
      public ::  Ncrd_2d_Int
      public ::  Ncrd_3d
      public ::  Ncrd_3d_Int
      public ::  Ncrd_4d
      public ::  Ncrd_5d
      public ::  Ncrd_1d_Char
      public ::  Ncrd_2d_Char
!      public  NetcdfRead
!
!      interface NetcdfRead
!         module procedure  Ncrd_Scal
!         module procedure  Ncrd_Scal_Int
!         module procedure  Ncrd_1d
!         module procedure  Ncrd_1d_Int
!         module procedure  Ncrd_2d
!         module procedure  Ncrd_2d_Int
!         module procedure  Ncrd_3d
!         module procedure  Ncrd_3d_Int
!         module procedure  Ncrd_4d
!         module procedure  Ncrd_5d
!         module procedure  Ncrd_1d_Char
!         module procedure  Ncrd_2d_Char
!      end interface
!
#     include "GmiParameters.h"
#     include "netcdf.inc"
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
! !IROUTINE: Ncrd_Scal
!
! !INTERFACE:
!
      subroutine Ncrd_Scal (varrd_scal, ncid, varname)
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
!!    varrd_scal : variable to fill
      real*8           , intent(out)  :: varrd_scal
!
! !DESCRIPTION:
!  Reads in a NetCDF scalar variable.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4              :: varrd_scal_tmp
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
        err_msg = 'In Ncrd_Scal #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Var_Real   (ncid, varid, varrd_scal_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_Scal #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      varrd_scal = varrd_scal_tmp

      return

      end subroutine Ncrd_Scal
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_Scal_Int
!
! !INTERFACE:
!
      subroutine Ncrd_Scal_Int (varrd_scali, ncid, varname)
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
!!    varrd_scali : integer variable to fill
      integer          , intent(out)  :: varrd_scali
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
        err_msg = 'In Ncrd_Scal_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Var_Int (ncid, varid, varrd_scali)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_Scal_Int #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncrd_Scal_Int
!
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_1d
!
! !INTERFACE:
!
      subroutine Ncrd_1d (varrd_1d, ncid, varname, start1d, count1d)
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    start1d   : vector specifying the index in varrd_1d where 
!!               the first of the data values will be read 
!!    count1d    : varrd_1d dimension
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start1d(1)
      integer          , intent(in)   :: count1d (1)
!
! !OUTPUT PARAMETERS:
!!    varrd_1d : array to fill
      real*8           , intent(out)  :: varrd_1d(count1d(1))
!
! !DESCRIPTION:
!  Reads in a 1D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4  :: varrd_1d_tmp(count1d(1))
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
        err_msg = 'In Ncrd_1d #1:  ' // Trim (varname) // &
                   ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr =  Nf_Get_Vara_Real   (ncid, varid, start1d, count1d, varrd_1d_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_1d #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      varrd_1d(:) = varrd_1d_tmp(:)

      return

      end subroutine Ncrd_1d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_1d_Int
!
! !INTERFACE:
!
      subroutine Ncrd_1d_Int (varrd_1di, ncid, varname, start1d, count1d)
!
      implicit none
!
! !INPUT PARAMETERS:
!
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    start1d   : vector specifying the index in varrd_1di where 
!!               the first of the data values will be read 
!!    count1d    : varrd_1di dimension
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start1d(1)
      integer          , intent(in)   :: count1d (1)
!
! !OUTPUT PARAMETERS:
!!    varrd_1di : intger array to fill
      integer          , intent(out)  :: varrd_1di(count1d(1))
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
        err_msg = 'In Ncrd_1d_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if


      ierr = Nf_Get_Vara_Int (ncid, varid, start1d, count1d, varrd_1di)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_1d_Int #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncrd_1d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_2d
!
! !INTERFACE:
!
      subroutine Ncrd_2d (varrd_2d, ncid, varname, start2d, count2d)
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    start2d   : vector specifying the index in varrd_2d where
!!               the first of the data values will be read
!!    count2d    : varrd_2d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start2d(2)
      integer          , intent(in)   :: count2d (2)
!
! !OUTPUT PARAMETERS:
!!    varrd_2d : array to fill
      real*8           , intent(out)  :: varrd_2d(count2d(1), count2d(2))
!
! !DESCRIPTION:
!  Reads in a 2D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4              :: varrd_2d_tmp(count2d(1), count2d(2))
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
        err_msg = 'In Ncrd_2d #1:  ' // Trim (varname) // & 
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Real   (ncid, varid, start2d, count2d, varrd_2d_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_2d #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      varrd_2d(:,:) = varrd_2d_tmp(:,:)

      return

      end subroutine Ncrd_2d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_2d_Int
!
! !INTERFACE:
!
      subroutine Ncrd_2d_Int (varrd_2di, ncid, varname, start2d, count2d)
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    start2d   : vector specifying the index in varrd_2d where
!!               the first of the data values will be read
!!    count2d    : varrd_2di dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start2d(2)
      integer          , intent(in)   :: count2d (2)
!
! !OUTPUT PARAMETERS:
!!    varrd_2di : intger array to fill
      integer          , intent(out)  :: varrd_2di(count2d(1), count2d(2))
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
        err_msg = 'In Ncrd_2d_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Int (ncid, varid, start2d, count2d, varrd_2di)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_2d_Int #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncrd_2d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_3d
!
! !INTERFACE:
!
      subroutine Ncrd_3d (varrd_3d, ncid, varname, start3d, count3d)
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    start3d   : vector specifying the index in varrd_3d where
!!               the first of the data values will be read
!!    count3d    : varrd_3d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start3d(3)
      integer          , intent(in)   :: count3d (3)
!
! !OUTPUT PARAMETERS:
!!    varrd_3d : array to fill
      real*8           , intent(out)  :: varrd_3d(count3d(1), count3d(2), &
                                                  count3d(3))
!
! !DESCRIPTION:
!  Reads in a 3D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4              :: varrd_3d_tmp(count3d(1), count3d(2), count3d(3))
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
        err_msg = 'In Ncrd_3d #1:  ' // Trim (varname) // &
                 ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Real (ncid, varid, start3d, count3d, varrd_3d_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_3d #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      varrd_3d(:,:,:) = varrd_3d_tmp(:,:,:)

      return

      end subroutine Ncrd_3d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_3d_Int
!
! !INTERFACE:
!
      subroutine Ncrd_3d_Int (varrd_3di, ncid, varname, start3d, count3d)
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    start3d   : vector specifying the index in varrd_3d where
!!               the first of the data values will be read
!!    count3d    : varrd_3di dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start3d(3)
      integer          , intent(in)   :: count3d (3)
!
! !OUTPUT PARAMETERS:
!!    varrd_3di : intger array to fill
      integer          , intent(out)  :: varrd_3di(count3d(1), count3d(2), &
                                                   count3d(3))
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
        err_msg = 'In Ncrd_3d_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Int (ncid, varid, start3d, count3d, varrd_3di)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_3d_Int #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncrd_3d_Int
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_4d
!
! !INTERFACE:
!
      subroutine Ncrd_4d (varrd_4d, ncid, varname, start4d, count4d)
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    start4d   : vector specifying the index in varrd_4d where
!!               the first of the data values will be read
!!    count4d    : varrd_4d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start4d(4)
      integer          , intent(in)   :: count4d (4)
!
! !OUTPUT PARAMETERS:
!!    varrd_4d : array to fill
      real*8           , intent(out)  :: varrd_4d(count4d(1), count4d(2), &
                                                  count4d(3), count4d(4))
!
! !DESCRIPTION:
!  Reads in a 4D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4              :: varrd_4d_tmp(count4d(1), count4d(2), count4d(3), &
                                                              count4d(4))
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
        err_msg = 'In Ncrd_4d #1:  ' // Trim (varname) // &
                    ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if


      ierr =  Nf_Get_Vara_Real   (ncid, varid, start4d, count4d, varrd_4d_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_4d #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      varrd_4d(:,:,:,:) = varrd_4d_tmp(:,:,:,:)

      return

      end subroutine Ncrd_4d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_5d
!
! !INTERFACE:
!
      subroutine Ncrd_5d (varrd_5d, ncid, varname, start5d, count5d)
!
      implicit none
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    start5d   : vector specifying the index in varrd_5d where
!!               the first of the data values will be read
!!    count5d    : varrd_5d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start5d(5)
      integer          , intent(in)   :: count5d (5)
!
! !OUTPUT PARAMETERS:
!!    varrd_5d : array to fill
      real*8         , intent(out)  :: varrd_5d(count5d(1), count5d(2), &
                                                count5d(3), count5d(4), &
                                                count5d(5))
!
! !DESCRIPTION:
!  Reads in a 5D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4              :: varrd_5d_tmp(count5d(1), count5d(2), count5d(3), &
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
        err_msg = 'In Ncrd_5d #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Real   (ncid, varid, start5d, count5d, varrd_5d_tmp)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_5d #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      varrd_5d(:,:,:,:,:) = varrd_5d_tmp(:,:,:,:,:)

      return

      end subroutine Ncrd_5d
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_1d_Char
!
! !INTERFACE:
!
      subroutine Ncrd_1d_Char (varrd_1dc, ncid, varname, start1d, count1d)
!
      implicit none
!
! !INPUT PARAMETERS:
!
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    start1d   : vector specifying the index in varrd_1di where 
!!               the first of the data values will be read 
!!    count1d    : varrd_1dc dimension
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start1d(1)
      integer          , intent(in)   :: count1d (1)
!
! !OUTPUT PARAMETERS:
!!    varrd_1dc : intger array to fill
      character (len=1), intent(out)  :: varrd_1dc(count1d(1))
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
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_1d_Char #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Text (ncid, varid, start1d, count1d, varrd_1dc)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_1d_Char #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncrd_1d_Char
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_2d_Char
!
! !INTERFACE:
!
      subroutine Ncrd_2d_Char (varrd_2dc, ncid, varname, start2d, count2d)
!
      implicit none
!
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    start2d   : vector specifying the index in varrd_2dc where
!!               the first of the data values will be read
!!    count2d    : varrd_2dc dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start2d(2)
      integer          , intent(in)   :: count2d (2)
!
! !OUTPUT PARAMETERS:
!!    varrd_2dc : charcter array to fill
      character        , intent(out)  :: varrd_2dc(count2d(1), count2d(2))
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
      ierr = Nf_Inq_Varid (ncid, varname, varid)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_2d_Char #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Get_Vara_Text (ncid, varid, start2d, count2d, varrd_2dc)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_2d_Char #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine Ncrd_2d_Char
!EOC
!------------------------------------------------------------------------
end module m_netcdf_io_read

