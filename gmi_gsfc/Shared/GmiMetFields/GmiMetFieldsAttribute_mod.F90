!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiMetFieldsAttribute_mod
!
! !INTERFACE:
!
      module GmiMetFieldsAttribute_mod
!
! !USES:
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: getMetFieldsAttribute
!
! !DESCRIPTION:
!
! !AUTHOR:
! Jules Kouatchou, Jules.Kouatchou-1@nasa.gov
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
      contains
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getMetFieldsAttribute
!
! !INTERFACE:
!
      subroutine getMetFieldsAttribute &
     &  (ncid, metdata_name, metdata_name_org,  &
     &   metdata_name_model, metdata_name_dims, pr_diag, procID)
!
! !USES:
      use GmiPrintError_mod, only : GmiPrintError

      implicit none

#     include "netcdf.inc"
#     include "GmiParameters.h"
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: ncid       ! netCDF file identifier
      integer, intent(in) :: procID
!
! !OUTPUT PARAMETERS:
      character (len=*), intent(out) :: metdata_name
      character (len=*), intent(out) :: metdata_name_org
      character (len=*), intent(out) :: metdata_name_model
      character (len=*), intent(out) :: metdata_name_dims
!
! !DESCRIPTION:
! This routine determines the attributes of the metFields file.
! It reads in the "Met_Data_Name" global attribute from the
! netCDF met input file.  It also parses the complete name string into
! three substings.
!
! !LOCAL VARIABLES:
      character (len=50 ) :: tmp_mname
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      character (len=50 ) :: metdata_set(50)
      logical             :: found
      integer             :: ierr, ii, ltrim
      integer             :: pos1, pos2
      integer             :: siz1, siz2, siz3
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'getMetFieldsAttribute called by ', procID

      ii = 0

      ii = ii + 1
      metdata_set(ii) = 'DAO_FVCCM_4x5x28'
!
      ii = ii + 1
      metdata_set(ii) = 'DAO_FVDAS_4x5x28'
!
      ii = ii + 1
      metdata_set(ii) = 'DAO_FVDAS_2x2.5x55'
!
      ii = ii + 1
      metdata_set(ii) = 'DAO_GEOS3_1x1x48'
      ii = ii + 1
      metdata_set(ii) = 'DAO_GEOS3_2x2.5x48'
!
      ii = ii + 1
      metdata_set(ii) = 'DAO_GS_2x2.5x20'
      ii = ii + 1
      metdata_set(ii) = 'DAO_GS_2x2.5x26'
      ii = ii + 1
      metdata_set(ii) = 'DAO_GS_2x2.5x29'
      ii = ii + 1
      metdata_set(ii) = 'DAO_GS_2x2.5x46'
!
      ii = ii + 1
      metdata_set(ii) = 'DAO_GS_4x5x29'
      ii = ii + 1
      metdata_set(ii) = 'DAO_GS_4x5x46'
!
      ii = ii + 1
      metdata_set(ii) = 'DAO_GU_2x2.5x25'
!
      ii = ii + 1
      metdata_set(ii) = 'GISS_2prime_4x5x23'
!
      ii = ii + 1
      metdata_set(ii) = 'NCAR_CCM2_4x5x44'
!
      ii = ii + 1
      metdata_set(ii) = 'NCAR_CCM3_4x5x18'
      ii = ii + 1
      metdata_set(ii) = 'NCAR_CCM3_T42x18'
!
      ii = ii + 1
      metdata_set(ii) = 'NCAR_MATCH_4x5x52'
!
      ii = ii + 1
      metdata_set(ii) = 'DAO_FVCCM_2.5x2x28'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GCM_2x2%5x28'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GCM_2x2%5x33'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GEOS4GCM_2%5x2x42'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GEOS4GCM_5x4x42'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GEOS4DAS_5x4x42'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GEOS5DAS_4x5x72'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GEOS5DAS_2x2%5x72'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GEOS5DAS_1x1%25x72'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GEOS5MERRA_2x2%5x72'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GEOS5CCM_2x2%5x72'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GEOS5MERRA300_2x2%5x72'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GEOS5DAS_2x2%5x55'
!
      ii = ii + 1
      metdata_set(ii) = 'GISS_2prime_4x5x46'
      ii = ii + 1
      metdata_set(ii) = 'EC_OSLO_2x2x37'
      ii = ii + 1
      metdata_set(ii) = 'EC_OSLO_2x2x40'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GEOS5MERRA300_1x1%25x72'
      ii = ii + 1
      metdata_set(ii) = 'GMAO_GEOS5MERRA_1x1%25x72'

      ierr = NF_Get_Att_Text(ncid, NF_GLOBAL, 'Met_Data_Name', tmp_mname)

      if (ierr /= NF_NOERR) then
         err_msg = 'In Get_Metdata_Name_Attr:  ' // Nf_Strerror (ierr)
         call GmiPrintError(err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      found = .false.

      do ii = 1, Size (metdata_set, 1)
         ltrim = Len_Trim (metdata_set(ii))

         if (Trim (metdata_set(ii)) == (tmp_mname(1:ltrim))) then
            metdata_name = metdata_set(ii)
            found = .true.
!           ====
            exit
!           ====
         end if
      end do

      if (.not. found) then
         Write (6,*) "tmp_mname = ", tmp_mname(1:ltrim)
         Write (6,*) "metdata_name = ", metdata_name
         err_msg = "Problem with metdata_name in getMetFieldsAttribute (unused?)"
         call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

      ltrim = Len_Trim (metdata_name)

      pos1 = Scan (metdata_name, "_")
      pos2 = Scan (metdata_name, "_", BACK = .true.)

      siz1 = pos1  - 1
      siz2 = pos2  - pos1 - 1
      siz3 = ltrim - pos2

      metdata_name_org  (1:siz1) = metdata_name(1:pos1-1)
      metdata_name_model(1:siz2) = metdata_name(pos1+1:pos2-1)
      metdata_name_dims (1:siz3) = metdata_name(pos2+1:ltrim)

      return

      end subroutine getMetFieldsAttribute
!EOC
!------------------------------------------------------------------------------
      end module GmiMetFieldsAttribute_mod
