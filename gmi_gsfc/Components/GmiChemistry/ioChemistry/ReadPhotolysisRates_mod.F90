!-------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:
!
      module ReadPhotolysisRates_mod
!
! !USES:
      use m_netcdf_io_open , only : Ncop_Rd
      use m_netcdf_io_close, only : Nccl
      use m_netcdf_io_checks, only : Ncdoes_Var_Exist
      use m_netcdf_io_read , only : Ncrd_1d, Ncrd_2d, Ncrd_3d, Ncrd_4d, Ncrd_5d, &
     &                              Ncrd_Scal_Int
      use m_netcdf_io_get_dimlen, only : Ncget_Dimlen
      use GmiPrintError_mod, only : GmiPrintError
!
      implicit none

      private
      public  :: readPhotolysisRates, readPhotolysisTable, readSolarCycle, readRaaQaa

#     include "netcdf.inc"
#     include "GmiParameters.h"
#     include "phot_lookup_constants.h"
!
! !AUTHOR:
! Jules Kouatchou, Jules.Kouatchou-1@nasa.gov
!
!EOP
!-------------------------------------------------------------------------
   contains
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readPhotolysisTable
!
! !INTERFACE:
!
      subroutine readPhotolysisTable(phot_opt, qj_infile_name, chem_mecha, &
     &               num_qjs, i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl)
!
      implicit none
!
#     include "phot_lookup_arrays.h"
#     include "phot_lookup.h"
!
! !INPUT PARAMETERS:
      integer,           intent(in) :: phot_opt, num_qjs
      integer,           intent(in) :: i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl
      character (len=*), intent(in) :: chem_mecha
      character (len=*), intent(in) :: qj_infile_name
!
! !DESCRIPTION:
! Reads a table of photolysis data including a radiative
!   source function, cross sections for most reactions, and photolysis rates
!   for O2 and NO.  This routine also reads the coordinates of the table
!   values (column ozone, pressure and solar zenith angle).  This information
!   will be used to calculate photolysis rates using it as a table and
!   interpolating linearly in that table (similar to Randy Kawa's photloysis
!   lookup model).
!
! !DEFINED PARAMETERS:
      character (len=MAX_LENGTH_VAR_NAME), parameter :: COL_OZONE_DNAM   = 'column_ozone'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: PRESS_DNAM       = 'pressure_dim'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: QJ_DNAM          = 'qj_dim'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: SZA_DNAM         = 'solar_zenith_angle'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: TEMP_DNAM        = 'temperature_dim'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: WVLEN_DNAM       = 'wavelength_dim'

      character (len=MAX_LENGTH_VAR_NAME), parameter :: O3_CLIM_LAT_DNAM = 'o3_clim_lat_dim'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: O3_CLIM_LON_DNAM = 'o3_clim_lon_dim'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: O3_CLIM_MON_DNAM = 'o3_clim_mon_dim'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: O3_CLIM_PRS_DNAM = 'o3_clim_p_dim'

      character (len=MAX_LENGTH_VAR_NAME), parameter :: COL_OZONE_VNAM   = 'column_ozone'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: CROSS_SECT_VNAM  = 'cross_sections'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: NO_JRATE_VNAM    = 'NO_jrate'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: NUMACET_VNAM     = 'NUM_QJ_ACET'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: NUMCH2O_VNAM     = 'NUM_QJ_CH2O'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: NUMHACN_VNAM     = 'NUM_QJ_HACN'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: NUMMGLY_VNAM     = 'NUM_QJ_MGLY'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: NUMNO_VNAM       = 'NUM_QJ_NO'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: NUMO2_VNAM       = 'NUM_QJ_O2'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: NUMO3_VNAM       = 'NUM_QJ_O3_TO_2OH'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: O2_JRATE_VNAM    = 'O2_jrate'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: PRESS_VNAM       = 'pressure_dim'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: RAD_SOURCE_VNAM  = 'radiative_source_function'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: SZA_VNAM         = 'solar_zenith_angle'

      character (len=MAX_LENGTH_VAR_NAME), parameter :: O3_CLIM_PRS_VNAM = 'o3_clim_p_dim'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: O3_CLIM_VNAM     = 'o3_clim'
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      character (len=MAX_LENGTH_VAR_NAME) :: cdimary4(4)
      character (len=MAX_LENGTH_VAR_NAME) :: cdimary6(6)
      integer :: ii, ii_num, iilam, iio3, iiprs, iisza
      integer :: ncid_qj, nferr ,var_id_o1d
      integer :: cnt1d (1), cnt2d (2), cnt3d (3), cnt4d (4)
      integer :: strt1d(1), strt2d(2), strt3d(3), strt4d(4)
      integer :: idimary4(4), ilong, ilat
      integer :: idimary6(6)
      real*8, allocatable  :: tmp_rad_source(:, :, :, :)
      real,   allocatable  :: tmp_o1d_coef(:, :, :, :, :)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ilong = i2 - i1  + 1
      ilat  = j2 - ju1 + 1
!     ---------------------------------------------------------
!     Do some error checking to make sure the dimensions in the
!     table are consistant with the dimensions in the include
!     file that describes the table.
!     ---------------------------------------------------------

      cdimary6 = (/ WVLEN_DNAM, PRESS_DNAM, TEMP_DNAM,  &
     &              QJ_DNAM,    SZA_DNAM,   COL_OZONE_DNAM /)

      idimary6 = (/ NUMLAM, NUMPRS, NUMTMP, num_qjs, NUMSZA, NUMO3 /)

      call Ncop_Rd (ncid_qj, qj_infile_name)

      do ii = 1, 6

        call Ncget_Dimlen (ncid_qj, cdimary6(ii), ii_num)

        if (ii_num /= idimary6(ii)) then

          err_msg = 'Error #1 in Read_Phot_Table:  ' // cdimary6(ii)

          call GmiPrintError  &
     &      (err_msg, .true., 2, ii_num, idimary6(ii), 0, 0.0d0, 0.0d0)

        end if

      end do

      if (phot_opt == 5) then

        cdimary4 = (/ O3_CLIM_LAT_DNAM, O3_CLIM_LON_DNAM,  &
     &                O3_CLIM_MON_DNAM, O3_CLIM_PRS_DNAM /)

        idimary4 = (/ j2_gl-ju1_gl+1, i2_gl-i1_gl+1,  &
     &                NUM_O3CLIM_MON, NUM_O3CLIM_PRS /)

        do ii = 1, 4

          call Ncget_Dimlen (ncid_qj, cdimary4(ii), ii_num)

          if (ii_num /= idimary4(ii)) then

            err_msg = 'Error #2 in Read_Phot_Table:  ' // cdimary4(ii)

            call GmiPrintError (err_msg, .true., 2, ii_num, idimary4(ii),  &
     &                       0, 0.0d0, 0.0d0)

          end if

        end do

      end if

      num_qj_acet      = 0
      num_qj_ch2o      = 0
      num_qj_hacn      = 0
      num_qj_mgly      = 0
      num_qj_no        = 0
      num_qj_o2        = 0
      num_qj_o3_to_2oh = 0

      if (Ncdoes_Var_Exist (ncid_qj, NUMACET_VNAM))  &
     &   call Ncrd_Scal_Int (num_qj_acet,      ncid_qj, NUMACET_VNAM)

      if (Ncdoes_Var_Exist (ncid_qj, NUMCH2O_VNAM))  &
     &   call Ncrd_Scal_Int (num_qj_ch2o,      ncid_qj, NUMCH2O_VNAM)

      if (Ncdoes_Var_Exist (ncid_qj, NUMHACN_VNAM))  &
     &   call Ncrd_Scal_Int (num_qj_hacn,      ncid_qj, NUMHACN_VNAM)

      if (Ncdoes_Var_Exist (ncid_qj, NUMMGLY_VNAM))  &
     &   call Ncrd_Scal_Int (num_qj_mgly,      ncid_qj, NUMMGLY_VNAM)

      if (Ncdoes_Var_Exist (ncid_qj, NUMNO_VNAM))  &
     &   call Ncrd_Scal_Int (num_qj_no,        ncid_qj, NUMNO_VNAM)

      if (Ncdoes_Var_Exist (ncid_qj, NUMO2_VNAM))  &
     &   call Ncrd_Scal_Int (num_qj_o2,        ncid_qj, NUMO2_VNAM)

      if (Ncdoes_Var_Exist (ncid_qj, NUMO3_VNAM))  &
     &   call Ncrd_Scal_Int (num_qj_o3_to_2oh, ncid_qj, NUMO3_VNAM)

      strt1d(:) = (/ 1 /)
      strt2d(:) = (/ 1, 1 /)
      strt3d(:) = (/ 1, 1, 1 /)
      strt4d(:) = (/ 1, 1, 1, 1 /)

      cnt1d(1) = NUMPRS

      call Ncrd_1d (prs_phot, ncid_qj, PRESS_VNAM, strt1d, cnt1d)

      cnt1d(1) = NUMSZA

      call Ncrd_1d (sza_phot, ncid_qj, SZA_VNAM,   strt1d, cnt1d)

      cnt2d(:) = (/ NUMO3, NUMPRS /)

      call Ncrd_2d (col_o3, ncid_qj, COL_OZONE_VNAM, strt2d, cnt2d)

      cnt3d(:) = (/ NUMLAM, NUMTMP, num_qjs /)

      call Ncrd_3d (cross_section, ncid_qj,  CROSS_SECT_VNAM, strt3d, cnt3d)

      cnt3d(:) = (/ NUMSZA, NUMO3, NUMPRS /)

      if (num_qj_no /= 0)  &
     &   call Ncrd_3d (no_qj, ncid_qj, NO_JRATE_VNAM, strt3d, cnt3d)

      cnt3d(:) = (/ NUMSZA, NUMO3, NUMPRS /)

      if (num_qj_o2 /= 0)  &
     &   call Ncrd_3d (o2_qj, ncid_qj, O2_JRATE_VNAM, strt3d, cnt3d)

      cnt4d(:) = (/ NUMSZA, NUMO3, NUMPRS, NUMLAM /)

      allocate (tmp_rad_source(NUMSZA, NUMO3, NUMPRS, NUMLAM))

      call Ncrd_4d (tmp_rad_source, ncid_qj, RAD_SOURCE_VNAM, strt4d, cnt4d)

!
!.... Added to read j adjustment factors for O3 + hv = O(1D) + O2
!.... PSC 030923
!

      if (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) then

         allocate(tmp_o1d_coef(NUMSZA, NUMO3, NUMPRS, 5, 2))

         nferr = NF_INQ_VARID    (ncid_qj ,'O1D_correlation' ,var_id_o1d)
         nferr = NF_GET_VAR_REAL (ncid_qj ,var_id_o1d ,tmp_o1d_coef)

         do iio3    = 1 ,2
           do iilam = 1 ,5
             o1d_coef(iilam,:,:,:,iio3) = dble(tmp_o1d_coef(:,:,:,iilam,iio3))
           end do
         end do

         deallocate(tmp_o1d_coef)

      end if

!.... End of added segment for o(1D)

      do iiprs = 1, NUMPRS
        do iio3 = 1, NUMO3
          do iisza = 1, NUMSZA
            do iilam = 1, NUMLAM

              rad_source(iilam,iisza,iio3,iiprs) =  &
     &          tmp_rad_source(iisza,iio3,iiprs,iilam)

            end do
          end do
        end do
      end do

      deallocate(tmp_rad_source)

      if (phot_opt == 5) then

        cnt1d(1) = NUM_O3CLIM_PRS

        call Ncrd_1d (o3_clim_prs, ncid_qj, O3_CLIM_PRS_VNAM, strt1d, cnt1d)

        cnt4d(:)  = (/ NUM_O3CLIM_PRS, ilong, ilat, NUM_O3CLIM_MON /)

        strt4d(:) = (/ 1, i1, ju1, 1 /)

        call Ncrd_4d (o3_clim, ncid_qj, O3_CLIM_VNAM, strt4d, cnt4d)

      end if

      call Nccl (ncid_qj)

      return

      end subroutine readPhotolysisTable
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readSolarCycle
!
! !INTERFACE:
!
      subroutine readSolarCycle(do_solar_cycle, sc_infile_name)
!
! !USES: 
      use solar_cycle_mod
!
      implicit none
!
! !INPUT PARAMETERS:
      logical          , intent(in) :: do_solar_cycle
      character (len=*), intent(in) :: sc_infile_name
!
! !DESCRIPTION:
! Reads a table of relative strength of the solar flux in each
!   phot_table bin. The file also contains the dates and wavelength bin values.
!   Data from Charley Jackman via the Goddard CTM
!
! !DEFINED PARAMETERS:
      character (len=MAX_LENGTH_VAR_NAME), parameter :: WVLEN_DNAM      = 'wavelength_dim'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: TIME_DNAM       = 'time_dim'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: SC_VNAM         = 's_cycle'
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      character (len=MAX_LENGTH_VAR_NAME) :: cdimary2(2)
      integer :: idimary2(2)
      integer :: ii
      integer :: ii_num
      integer :: ncid_sc
      integer :: cnt1d (1), cnt2d (2)
      integer :: strt1d(1), strt2d(2)
!EOP
!-------------------------------------------------------------------------
!BOC
      if (do_solar_cycle) then
         ! ---------------------------------------------------------
         ! Do some error checking to make sure the dimensions in the
         ! table are consistant with the dimensions in the include
         ! file that describes the table and with lookup table file.
         ! ---------------------------------------------------------

         cdimary2 =  (/ WVLEN_DNAM, TIME_DNAM /)

         idimary2 =  (/ NUM_SOLAR_CYC_LAM, NUM_SOLAR_CYC_MON /)

!... open netCDF file
         call Ncop_Rd (ncid_sc, sc_infile_name)

         do ii = 1, 2
            call Ncget_Dimlen (ncid_sc, cdimary2(ii), ii_num)

            if (ii == 1 .and. ii_num /= NUMLAM) then
               err_msg = 'Error #1 in Read_Solar_Cycle 1: ' // cdimary2(ii)
               call GmiPrintError(err_msg, .true., 2, ii_num, NUMLAM, 0, 0.0d0, 0.0d0)
            endif

            if (ii_num /= idimary2(ii)) then
               err_msg = 'Error #1 in Read_Solar_Cycle 2: ' // cdimary2(ii)
               call GmiPrintError(err_msg, .true., 2, ii_num, NUMLAM, 0, 0.0d0, 0.0d0)
            endif
         enddo

!         strt1d(:) = (/ 1 /)
!         cnt1d(1) = NUM_SOLAR_CYC_MON
!
!         call Ncrd_1d (s_cycle_date, ncid_sc, SC_DATE_VNAM, strt1d, cnt1d)

         strt2d(:) = (/ 1, 1 /)
         cnt2d(:) = (/ idimary2(1), idimary2(2) /)

         call Ncrd_2d (s_cycle, ncid_sc, SC_VNAM, strt2d, cnt2d)

         call Nccl (ncid_sc)
      else
         ! set entire s_cycle array to 1 if no solar cycle
         s_cycle = 1.d0
      end if

      return

      end subroutine readSolarCycle
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readPhotolysisRates
!
! !INTERFACE:
!
  subroutine readPhotolysisRates(i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, &
                                 num_qjs, qj_timpyr, &
                                 qjmon, qj_infile_name, qj_var_name)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, qj_timpyr, num_qjs
      integer, intent(in) :: i1_gl, ju1_gl
      character (len=*), intent(in) :: qj_infile_name
      character (len=*), intent(in) :: qj_var_name
!
! !INPUT/OUTPUT PARAMETERS:
      real*8, intent(inOut) :: qjmon(i1:i2, ju1:j2, k1:k2, num_qjs, qj_timpyr)
!
! !LOCAL VARIABLES:
      integer :: il, ij, ik, it, iq, ilong, ilat, ivert
      integer :: inb, jnb, ncid_qj
      integer :: cnt5d (5), strt5d(5)
      real*8, allocatable  :: qj_tmp(:,:,:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      strt5d(1) = i1  -  i1_gl + 1
      strt5d(2) = ju1 - ju1_gl + 1
!      strt5d(1) =  1
!      strt5d(2) =  1
      strt5d(3) = k1

     ilong = i2 - i1  + 1
     ilat  = j2 - ju1 + 1
     ivert = k2 - k1  + 1

     allocate (qj_tmp(ilong, ilat, ivert, 1, 1))

     cnt5d (:) = (/ ilong, ilat, ivert, 1, 1 /)

     call Ncop_Rd (ncid_qj, qj_infile_name)

     do it = 1, qj_timpyr

         strt5d(5) = it

         do iq = 1, num_qjs

            strt5d(4) = iq

            call Ncrd_5d (qj_tmp, ncid_qj, qj_var_name, strt5d, cnt5d)

            do ik = k1, k2
              do ij = ju1, j2
                 jnb = ij - ju1 + 1
                do il = i1, i2
                   inb = il - i1 + 1

                   qjmon(il,ij,ik,iq,it) = qj_tmp(inb,jnb,ik,1,1)

                end do
              end do
            end do
         end do

     end do

     deallocate(qj_tmp)

     call Nccl (ncid_qj)

  return

  end subroutine readPhotolysisRates
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readRaaQaa
!
! !INTERFACE:
!
  subroutine readRaaQaa(raa_data_file, RAA, QAA)

      implicit none
!
! !DESCRIPTION:
! Reads a table of RAA and QAA fastj(x) variables for aerosol opt depth calc
!
# include "gmi_AerDust_const.h"
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: raa_data_file
!
! !OUTPUT PARAMETERS:
      real*8, intent(out) :: RAA(4, NP_b), QAA(4, NP_b)
!
! !LOCAL VARIABLES:
      integer :: cnt2d (2), strt2d(2)
      real*8, allocatable  :: xaa_tmp(:,:)
      integer ncid_raa

!... begin
      allocate (xaa_tmp(4, NP_b))

!... open netCDF file
      call Ncop_Rd (ncid_raa, raa_data_file)

      strt2d(:) = (/ 1, 1 /)
      cnt2d(:) = (/ 4, NP_b /)

      call Ncrd_2d (xaa_tmp, ncid_raa, 'RAA', strt2d, cnt2d)
      RAA = xaa_tmp
      call Ncrd_2d (xaa_tmp, ncid_raa, 'QAA', strt2d, cnt2d)
      QAA = xaa_tmp

      call Nccl (ncid_raa)

      deallocate (xaa_tmp)

      return

  end subroutine readRaaQaa
end module ReadPhotolysisRates_mod
