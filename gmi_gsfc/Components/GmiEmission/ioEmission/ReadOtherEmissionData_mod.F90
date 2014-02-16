    module ReadOtherEmissionData_mod

    use GmiTimeControl_mod, only : GmiSplitDateTime
    use GmiPrintError_mod , only : GmiPrintError
    use GmiASCIIoperations_mod, only : AsciiOpenRead
    use m_set_NLANDHAR, only : NLANDHAR

    implicit none

!=============================================================================
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
!   Original code from:
!     Harvard tropospheric emissions module for 3D applications;
!       by Yuhang Wang, Gerry Gardner, and Prof. Daniel Jacob
!       of Harvard University (Release V1.0)
!
!
    private
    public  :: readLightData, readSoilData, readFertilizerData
    public  :: readIsopreneConvertData, readMonoterpeneConvertData
    public  :: readPrecipitationData
!
#     include "GmiParameters.h"
#     include "gmi_emiss_constants.h"
#     include "gmi_time_constants.h"

!=============================================================================

    contains

!-----------------------------------------------------------------------------
!
! ROUTINE
!   readLightData

      subroutine readLightData  &
     &  (light_infile_name, sopcoeff, pr_diag, loc_proc)

      implicit none
!
! !INPUT PARAMETERS:
      logical          , intent(in) :: pr_diag
      integer          , intent(in) :: loc_proc
      character (len=*), intent(in) :: light_infile_name
!                                      light input file name
!
! !OUTPUT PARAMETERS:
      real*8           , intent(out) :: sopcoeff(NPOLY)
!                                       coefficients used for polynomial fit
!
! !DESCRIPTION:
! This routine reads in the polynomial coefficients for isoprene emissions.
!
! !LOCAL VARIABLES:
      character (len=80) :: cdummy
      integer :: ii, lun
!
!EOP
!-----------------------------------------------------------------------------
!BOC
      if (pr_diag) then
        Write (6,*) 'readLightData called by ', loc_proc
      end if

!     ================
      call AsciiOpenRead  &
!     ================
     &  (lun, light_infile_name)

      Read (lun,900) cdummy

      Read (lun,910) (sopcoeff(ii), ii=1,NPOLY)

 900  format (a80)
 910  format (8(1pe10.2))

      Close (lun)

      return

      end subroutine readLightData
!EOC
!-----------------------------------------------------------------------------
!BOP
! ROUTINE
!   readIsopreneConvertData
!
! DESCRIPTION
!   This routine reads in the isoprene conversion table and constructs the
!   isoprene base emission.
!
! ARGUMENTS
!   isopconv_infile_name : isoprene conversion infile name
!   convert_isop         : conversion table used to construct isoprene base
!                          emission
!
!-----------------------------------------------------------------------------

      subroutine readIsopreneConvertData  &
     &  (isopconv_infile_name, convert_isop, pr_diag, loc_proc)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      character (len=MAX_LENGTH_FILE_NAME) :: isopconv_infile_name
      real*8  :: convert_isop(NVEGTYPE)

!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=80) :: cdummy

      integer :: idummy
      integer :: ii
      integer :: lun

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'readIsopreneConvertData called by ', loc_proc
      end if

!     ================
      call AsciiOpenRead  &
!     ================
     &  (lun, isopconv_infile_name)

      Read (lun, 900) cdummy

 900  format (a80)

      do ii = 1, NVEGTYPE

        Read (lun, *) idummy, convert_isop(ii)

      end do

      Close (lun)

      return

      end subroutine readIsopreneConvertData


!-----------------------------------------------------------------------------
!
! ROUTINE
!   readMonoterpeneConvertData
!
! DESCRIPTION
!   This routine reads in the monoterpene conversion table and constructs the
!   monoterpene base emission.
!
! ARGUMENTS
!   monotconv_infile_name : monoterpene conversion infile name
!   convert_monot         : conversion table used to construct monoterpene
!                           base emission
!
!-----------------------------------------------------------------------------

      subroutine readMonoterpeneConvertData  &
     &  (monotconv_infile_name, convert_monot, pr_diag, loc_proc)


      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      character (len=*) :: monotconv_infile_name
      real*8  :: convert_monot(NVEGTYPE)

!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=80) :: cdummy

      integer :: idummy
      integer :: ii
      integer :: lun

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'readMonoterpeneConvertData called by ', loc_proc
      end if

!     ================
      call AsciiOpenRead  &
!     ================
     &  (lun, monotconv_infile_name)

      Read (lun, 900) cdummy

 900  format (a80)

      do ii = 1, NVEGTYPE
        Read (lun, *) idummy, convert_monot(ii)
      end do

      Close (lun)

      return

      end subroutine readMonoterpeneConvertData

!-----------------------------------------------------------------------------

!   ncon_soil  : Olson -> soil type
!   soil_infile_name     : soil type        input file name

      subroutine readSoilData &
     &  (soil_infile_name, ncon_soil, pr_diag, loc_proc)

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      character (len=*) :: soil_infile_name
      integer :: ncon_soil (NVEGTYPE)

      integer :: lun, ii, mm
      character (len=80) :: cdummy

      if (pr_diag) then
        Write (6,*) 'readSoilData called by ', loc_proc
      end if

      call AsciiOpenRead (lun, soil_infile_name)

      Read (lun, 900) cdummy

 900  format (a80)

      do mm = 1, NVEGTYPE
         Read (lun, *) ii, ncon_soil(ii)
      end do

      Close (lun)

      return

      end subroutine readSoilData

!-----------------------------------------------------------------------------

!   fertscal_infile_name : fertilizer scale input file name
!   index_soil : i,j of the grid
!   soil_fert  : fertilizers    (ng N/m^2/s)

      subroutine readFertilizerData &
     &  (fertscal_infile_name, index_soil, soil_fert, pr_diag, loc_proc)

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      character (len=*) :: fertscal_infile_name
      integer :: index_soil(2, NLANDHAR)
      real*8  :: soil_fert (NLANDHAR)

      integer :: lun, mm, ii

      if (pr_diag) then
        Write (6,*) 'readFertilizerData called by ', loc_proc
      end if

      call AsciiOpenRead (lun, fertscal_infile_name)

      do mm = 1, NLANDHAR
          Read (lun, *) (index_soil(ii,mm), ii=1,2), soil_fert(mm)
      end do

      Close (lun)

      return

      end subroutine readFertilizerData

!-----------------------------------------------------------------------------
!
! ROUTINE
!   readPrecipitationData
!
! DESCRIPTION
!   This routine reads in three types of soil data needed to calculate soil
!   NOx emissions:
!     definition of the grid ii ,jj (once),
!     ferterlizer                   (once),
!     two months of precipitation   (monthly).
!
! ARGUMENTS
!   precip_infile_name   : precipitation    input file name
!   nymd      : year/month/day (YYYYMMDD)
!   soil_prep  : two months of observed precipitation (mm/day/box)
!
!-----------------------------------------------------------------------------

      subroutine readPrecipitationData  &
     &  (precip_infile_name, index_soil, soil_prep, nymd, soil_month, pr_diag, loc_proc)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      character (len=*) :: precip_infile_name
      integer, intent(in) :: nymd
      integer, intent(in) :: index_soil(2, NLANDHAR)
      real*8              :: soil_prep (2, NLANDHAR)
      integer, intent(inOut) :: soil_month ! to determine if a new month

!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      integer :: ddummy, ydummy
      integer :: ii, jj, kk, mm
      integer :: imonth
      integer :: lun
      integer :: prev_soil_month

      real*8  :: sp_tmp(MONTHS_PER_YEAR)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'readPrecipitationData called by ', loc_proc
      end if

!     ---------------------------
!     Read in precipitation data.
!     ---------------------------

!     ====================
      call GmiSplitDateTime  &
!     ====================
     &  (nymd, ydummy, imonth, ddummy)

      if (soil_month /= imonth) then

        soil_month = imonth

        if (soil_month == 1) then
          prev_soil_month = MONTHS_PER_YEAR
        else
          prev_soil_month = soil_month - 1
        end if

!       ================
        call AsciiOpenRead  &
!       ================
     &    (lun, precip_infile_name)

        do mm = 1, NLANDHAR

!          Read (lun, *) ii, jj, (sp_tmp(kk), kk=1,MONTHS_PER_YEAR)
          Read (lun, '(2I3,12f6.2)') ii, jj, (sp_tmp(kk), kk=1,MONTHS_PER_YEAR)

          if ((index_soil(1,mm) /= ii) .or.  &
     &        (index_soil(2,mm) /= jj)) then

            err_msg = 'Problem in Read_Soil_Data.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, ii, jj, 0, 0.0d0, 0.0d0)

          else

            soil_prep(1,mm) = sp_tmp(prev_soil_month)
            soil_prep(2,mm) = sp_tmp(soil_month)

          end if

        end do

        Close (lun)

      end if


      return

      end subroutine readPrecipitationData

!-------------------------------------------------------------------------
    end module ReadOtherEmissionData_mod

