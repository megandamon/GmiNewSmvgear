!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_ncGeneralOpsOutput
!
! !INTERFACE:
!
module m_ncGeneralOpsOutput
!
! !USES
      use m_netcdf_io_close     , only : Nccl, Nccl_Noerr
      use m_netcdf_io_write     , only : Ncwr_2d_Int, Ncwr_1d_Int
      use m_netcdf_io_write     , only : Ncwr_1d, Ncwr_Scal, Ncwr_2d_Char
      use m_netcdf_io_create    , only : Ncdo_Sync, Nccr_Wr
      use m_netcdf_io_define    , only : NcDef_glob_attributes, NcDef_dimension
      use m_netcdf_io_define    , only : NcDef_var_attributes, NcDef_variable
      use m_netcdf_io_define    , only : NcSetFill, NcEnd_def
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: Define_Netcdf_Out_Gen
  public  :: Write_Netcdf_Hdr_Gen , WriteNetcdfHdrGen
  public  :: Is_Out_Freq_Time
!
#     include "GmiParameters.h"
!
! !DESCRIPTION:
!  Provides routines for the definition/writing of common variables in
!  netCDF output files.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define_Netcdf_Out_Gen
!
! !INTERFACE:
!
      subroutine Define_Netcdf_Out_Gen (ncid, nhdf, numLon, numLat, numVert,   &
     &                 hdfd, lond, latd, prsd, hdf_dim_name, lon_dim_name,     &
     &                 lat_dim_name, prs_dim_name)
!
      implicit none
!
#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
      integer,           intent(in) :: ncid         ! netCDF file identifier
      integer,           intent(in) :: nhdf         ! number of header fields
      integer,           intent(in) :: numLon       ! number of longitudes
      integer,           intent(in) :: numLat       ! number of latitudes
      integer,           intent(in) :: numVert      ! number of vertical levels
      character (len=*), intent(in) :: hdf_dim_name ! hdf       dimension name
      character (len=*), intent(in) :: prs_dim_name ! pressure  dimension name
      character (len=*), intent(in) :: lat_dim_name ! latitude  dimension name
      character (len=*), intent(in) :: lon_dim_name ! longitude dimension name
!
! !OUTPUT PARAMETERS:
      integer, intent(out) :: hdfd(1)
      integer, intent(out) :: lond(1), latd(1), prsd(1)
!
! !DESCRIPTION:
! Makes some general definitions for the NetCDF output file.
!
! !LOCAL VARIABLES:
      integer :: ierr
      integer :: varid
!EOP
!------------------------------------------------------------------------------
!BOC
!     ------------------
!     Define dimensions.
!     ------------------

      call NcDef_dimension(ncid, hdf_dim_name, nhdf, hdfd(1))
!                                 ------------
      call NcDef_dimension(ncid, lon_dim_name, numLon, lond(1))
!                                 ------------
      call NcDef_dimension(ncid, lat_dim_name, numLat, latd(1))
!                                 ------------
      call NcDef_dimension(ncid, prs_dim_name, numVert, prsd(1))
!                                 ------------
!     -----------------------------------------
!     Define variables and variable attributes.
!     -----------------------------------------

      call NcDef_variable (ncid, hdf_dim_name, NF_INT, 1, hdfd, varid)
!                                 ------------
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Header index')
      call NcDef_var_attributes (ncid, varid, 'units', 'unitless')

      call NcDef_variable (ncid, lon_dim_name, NF_FLOAT, 1, lond, varid)
!                                 ------------
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Longitude')
      call NcDef_var_attributes (ncid, varid, 'units', 'degrees_east')

      call NcDef_variable (ncid, lat_dim_name, NF_FLOAT, 1, latd, varid)
!                                 ------------
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Latitude')
      call NcDef_var_attributes (ncid, varid, 'units', 'degrees_north')

      call NcDef_variable (ncid, prs_dim_name, NF_FLOAT, 1, prsd, varid)
!                                 ------------
      call NcDef_var_attributes (ncid, varid, 'long_name', 'Pressure')
      call NcDef_var_attributes (ncid, varid, 'units', 'hPa')

      return

      end subroutine Define_Netcdf_Out_Gen
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Hdr_Gen
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Hdr_Gen  (ncid, latdeg, londeg, pr_diag,         &
     &                 procID, iMin, iMax, jMin, jMax, hdf_dim_name,         &
     &                 lat_dim_name, lon_dim_name)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical,           intent(in) :: pr_diag
      integer,           intent(in) :: ncid
      integer,           intent(in) :: procID
      integer,           intent(in) :: iMin, iMax, jMin, jMax
      real*8 ,           intent(in) :: latdeg(jMin:jMax) ! latitude  (deg)
      real*8 ,           intent(in) :: londeg(iMin:iMax) ! longitude (deg)
      character (len=*), intent(in) :: hdf_dim_name ! hdf       dimension name
      character (len=*), intent(in) :: lat_dim_name ! latitude  dimension name
      character (len=*), intent(in) :: lon_dim_name ! longitude dimension name
!
! !DESCRIPTION:
! Creates some general header information for a netCDF output file and writes
! it out.
!
! !LOCAL VARIABLES:
      integer :: ih
      integer :: cnt1d (1)
      integer :: strt1d(1)
      real*8 , allocatable :: locLatDeg(:)
      real*8 , allocatable :: locLonDeg(:)
      integer              :: ilen, jlen
      integer :: hdfdat(NETCDF_HDF)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Hdr_Gen called by ', procID

      strt1d(1) = 1

!     ----
!     Hdf.
!     ----

      do ih = 1, NETCDF_HDF
        hdfdat(ih) = ih
      end do

      cnt1d(1) = NETCDF_HDF

      call Ncwr_1d_Int (hdfdat, ncid, hdf_dim_name, strt1d, cnt1d)

      ilen = iMax - iMin + 1
      jlen = jMax - jMin + 1

      allocate(locLonDeg(1:ilen))
      allocate(locLatDeg(1:jlen))

      locLonDeg(1:ilen) = londeg(iMin:iMax)
      locLatDeg(1:jlen) = latdeg(jMin:jMax)

!     ----------
!     Longitude.
!     ----------

      cnt1d(1) = ilen

      call Ncwr_1d (locLonDeg, ncid, lon_dim_name, strt1d, cnt1d)

!     ---------
!     Latitude.
!     ---------

      cnt1d(1) = jlen

      call Ncwr_1d (locLatDeg, ncid, lat_dim_name, strt1d, cnt1d)

      deallocate(locLatDeg)
      deallocate(locLonDeg)

      return

      end subroutine Write_Netcdf_Hdr_Gen
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: WriteNetcdfHdrGen
!
! !INTERFACE:
!
      subroutine WriteNetcdfHdrGen  (ncid, locLatDeg, locLonDeg, pr_diag,      &
     &                procID, ilen, jlen, hdf_dim_name, lat_dim_name,          &
     &                lon_dim_name)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: ncid
      integer, intent(in) :: procID
      integer, intent(in) :: ilen, jlen
      real*8 , intent(in) :: locLatDeg(1:jlen)
      real*8 , intent(in) :: locLonDeg(1:ilen)
      character (len=*), intent(in) :: hdf_dim_name ! hdf       dimension name
      character (len=*), intent(in) :: lat_dim_name ! latitude  dimension name
      character (len=*), intent(in) :: lon_dim_name ! longitude dimension name
!
! !LOCAL VARIABLES:
      integer :: ih
      integer :: cnt1d (1)
      integer :: strt1d(1)
      integer :: hdfdat(NETCDF_HDF)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_Hdr_Gen called by ', procID

      strt1d(1) = 1

!     ----
!     Hdf.
!     ----

      do ih = 1, NETCDF_HDF
        hdfdat(ih) = ih
      end do

      cnt1d(1) = NETCDF_HDF

      call Ncwr_1d_Int (hdfdat, ncid, hdf_dim_name, strt1d, cnt1d)

!     ----------
!     Longitude.
!     ----------

      cnt1d(1) = ilen

      call Ncwr_1d (locLonDeg, ncid, lon_dim_name, strt1d, cnt1d)

!     ---------
!     Latitude.
!     ---------

      cnt1d(1) = jlen

      call Ncwr_1d (locLatDeg, ncid, lat_dim_name, strt1d, cnt1d)

      return

      end subroutine WriteNetcdfHdrGen
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Is_Out_Freq_Time
!
! !INTERFACE:
!
      function Is_Out_Freq_Time (month_save, month, day, nhms, ndt, gmi_sec,   &
     &               pr_period, pr_at_time)
!
      implicit none

#     include "gmi_time_constants.h"
!
! !INPUT PARAMETERS:
      integer, intent(in) :: month_save  ! the month during the last time step
      integer, intent(in) :: month       ! the current month
      integer, intent(in) :: day         ! the current day
      integer, intent(in) :: nhms        ! the current hour/minute/second
      integer, intent(in) :: ndt         ! model time step (s)
      real*8 , intent(in) :: gmi_sec     ! total Gmimod seconds (s)
      real*8 , intent(in) :: pr_period   ! periodic printing interval (s)
      real*8 , intent(in) :: pr_at_time
!
! !RETURN VALUE:
      logical :: Is_Out_Freq_Time
!
! !DESCRIPTION:
!   This routine determines if it is time to do some output based on the
!   following values for pr_period =>
! 
!     >0.0d0:  output at specified interval (s)
!     -1.0d0:  output at monthly intervals
!     -2.0d0:  output on 1st & 15th of each month
!
! !LOCAL VARIABLES:
      integer :: print_days(6), ipr
      logical :: is_time
      logical, save :: printed_on_this_day = .false.
      integer :: prflg
      real*8  :: rsecpdy
!EOP
!------------------------------------------------------------------------------
!BOC
      is_time = .false.

      rsecpdy = SECPDY

!     ======================
      if (pr_period < 0.0d0) then
!     ======================

        prflg = Nint (pr_period / rsecpdy)

        if (prflg == -1) then

!         ---------------
!         ---------------
!         Monthly output.
!         ---------------

          if (month /= month_save) is_time = .true.

        else if (prflg == -2) then

!         ------------------------------
!         Bimonthly output (1st & 15th).
!         ------------------------------

          if ((day == 1) .or. (day == 15)) then

            if ((nhms >= 0) .and. (.not. (printed_on_this_day))) then
              is_time = .true.
              printed_on_this_day = .true.
            end if

          else

            printed_on_this_day = .false.

          end if

        else if (prflg == -3) then
!         ------------------------------
!         Instantaneous output (1st, 6th, 11th, 16th, 21st, 26th).
!         (Note: this section added by Bigyani)
!         ------------------------------
          data (print_days(ipr), ipr=1,6) /1, 6, 11, 16, 21, 26/

          if (any(print_days == day)) then

            if ((nhms >= 0) .and. (.not. (printed_on_this_day))) then
              is_time = .true.
              printed_on_this_day = .true.
            end if

          else

            printed_on_this_day = .false.

          end if

        end if

!     ======================================================
      else if (Mod (Nint (gmi_sec+pr_at_time), Nint (pr_period)) < ndt) then
!     ======================================================

!       -----------------------------------
!       Specified time for periodic output.
!       -----------------------------------

        is_time = .true.

!     ======
      end if
!     ======

      Is_Out_Freq_Time = is_time

      return

      end function Is_Out_Freq_Time
!EOC
!-----------------------------------------------------------------------------

end module m_ncGeneralOpsOutput
