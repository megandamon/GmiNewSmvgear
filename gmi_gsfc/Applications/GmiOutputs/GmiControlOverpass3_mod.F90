!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlOverpass3_mod
!
! !INTERFACE:
!
      module GmiControlOverpass3_mod
!
! !USES:
      use GmiSub2Glob_mod              , only : subDomain2Global
      use GmiMessagePassing_mod        , only : synchronizeGroup
      use GmiFileOperations_mod        , only : makeOutfileName
      use m_netcdf_io_close            , only : Nccl, Nccl_Noerr
      use m_netcdf_io_write, only : Ncwr_2d_Int, Ncwr_5d, Ncwr_3d, Ncwr_2d, &
     &       Ncwr_4d
      use m_netcdf_io_write            , only : Ncwr_1d, Ncwr_Scal, Ncwr_2d_Char
      use m_netcdf_io_create           , only : Ncdo_Sync, Nccr_Wr
      use m_netcdf_io_define           , only : NcDef_glob_attributes, NcDef_dimension
      use m_netcdf_io_define           , only : NcDef_var_attributes, NcDef_variable
      use m_netcdf_io_define           , only : NcSetFill, NcEnd_def
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentrationSurf
      use GmiSpcConcentrationMethod_mod, only : Get_concentrationColTrop
      use GmiSpcConcentrationMethod_mod, only : Get_concentrationColCombo
      use GmiChemistryMethod_mod, only : t_Chemistry, Get_overheadO3col,       &
     &       Get_num_qjs, Get_num_qjo, Get_num_qks, Get_qj_labels,             &
     &       Get_const_labels, Get_qk_labels

      use GmiArrayBundlePointer_mod    , only : t_GmiArrayBundle

      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_rootProc,        &
     &       Get_procID, Get_iAmRootProc,         &
     &       Get_numDomains, Get_communicatorWorld, Get_map1_u
      use GmiDomainDecomposition_mod, only : Get_mcorGlob, Get_latdeg, Get_londeg

      use GmiGrid_mod, only : t_gmiGrid, Get_numSpecies, Get_i1, Get_i2,       &
     &       Get_ju1, Get_j2, Get_k1, Get_k2, Get_i1_gl, Get_i2_gl, Get_ju1_gl,&
     &       Get_j2_gl, Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jhi_gl, &
     &                        Get_ilo, Get_ihi, Get_julo, Get_jhi

      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Set_do_qqjk_reset,   &
     &       Get_hdf_dim_name, Get_lat_dim_name, Get_lon_dim_name, Get_k1r_gl, &
     &       Get_hdr_var_name, Get_rec_dim_name, Get_prs_dim_name, Get_k2r_gl, &
     &       Get_qj_var_name, Get_qj_dim_name, Get_qqj_var_name,               &
     &       Get_qqj_dim_name, Get_qqk_var_name, Get_qqk_dim_name, Get_do_mean,&
     &       Get_pr_diag, Get_problem_name, Get_pr_overpass3,                  &
     &       Get_pr_kel_overpass3, Get_pr_psf_overpass3, Get_pr_qj_overpass3,  &
     &       Get_pr_qqjk_overpass3, Get_pr_const_overpass3,                    &
     &       Get_pr_metwater_overpass3, Get_pr_totalMass_overpass3,            &
     &       Get_pr_relHumidity_overpass3, Get_pr_gridBoxHeight_overpass3,     &
     &       Get_pr_cloudOptDepth_overpass3, Get_pr_overheadO3col_overpass3,   &
     &       Get_pr_cloudFraction_overpass3, Get_pr_tropopausePress_overpass3, &
     &       Get_pr_lightningNO_overpass3, &
     &       Get_numSpecies_overpass3, Get_species_overpass3,                  &
     &       Get_begTime_overpass3, Get_endTime_overpass3, Get_pr_overpass3_period

      use GmiTimeControl_mod, only : t_GmiClock, GmiSplitDateTime, isOutTime, &
     &       Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds, Get_gmiTimeStep, &
     &       Get_numTimeSteps

      use GmiMetFieldsControl_mod, only : t_metFields, Get_pt, Get_mdt, Get_ai,&
     &       Get_bi, Get_am, Get_bm, Get_kel, Get_pctm2, Get_humidity,         &
     &       Get_mass, Get_tropopausePress, Get_potentialVorticity,            &
     &       Get_gridBoxHeight, Get_relativeHumidity, Get_totalCloudFraction,  &
     &       Get_tau_cloud

      use GmiEmissionMethod_mod, only : t_Emission, Get_lightning_no
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: controlOutputOverpass3, initializeOutputOverpass3
      public  :: finalizeOutputOverpass3
!
#     include "gmi_phys_constants.h"
#     include "gmi_diag_constants_llnl.h"
#     include "GmiParameters.h"
!
! !DESCRIPTION:
! Contains the necessary routines for the overpass3 outputs. Three routines
! are visible from the outside:
! %
! \begin{description}
! \item[initializeOutputOverpass3:] called once in the initialization steps to
!      open the netCDF file, define variables in the file, write out header
!      information and allocate variables.
! \item[controlOutputOverpass3:] called at the end of each model time step to
!      update variables and if needed for communications and writing out
!      data in the file.
! \item[finalizeOutputOverpass3:] called at the end of the integration to
!      deallocate all the variables and close the netCDF file.
! \end{description}
! %
! The above routines have as arguments derived types.
! %
! This module is self-contained. All the main routines are called by all
! the processors (master included) but each processor does not execute
! all portions of the code.
! There are two categories of arrays: (1) the ones manipulate by all the
! processors, and (2) those handle by worker processors only.
!
      ! Variables manipulate by all processors
      real*8, pointer, save :: qj_overpass3                (:,:,:,:) => null()
      real*8, pointer, save :: qqj_overpass3               (:,:,:,:) => null()
      real*8, pointer, save :: qqk_overpass3               (:,:,:,:) => null()
      real*8, pointer, save :: const_overpass3             (:,:,:,:) => null()
      real*8, pointer, save :: qj_nc_overpass3               (:,:,:) => null()
      real*8, pointer, save :: qqj_nc_overpass3              (:,:,:) => null()
      real*8, pointer, save :: qqk_nc_overpass3              (:,:,:) => null()
      real*8, pointer, save :: kel_nc_overpass3              (:,:,:) => null()
      real*8, pointer, save :: psf_nc_overpass3                (:,:) => null()
      real*8, pointer, save :: const_nc_overpass3            (:,:,:) => null()
      real*8, pointer, save :: metwater_nc_overpass3         (:,:,:) => null()
      real*8, pointer, save :: totalMass_nc_overpass3        (:,:,:) => null()
      real*8, pointer, save :: relHumidity_nc_overpass3      (:,:,:) => null()
      real*8, pointer, save :: gridBoxHeight_nc_overpass3    (:,:,:) => null()
      real*8, pointer, save :: cloudOptDepth_nc_overpass3    (:,:,:) => null()
      real*8, pointer, save :: overheadO3col_nc_overpass3    (:,:,:) => null()
      real*8, pointer, save :: cloudFraction_nc_overpass3    (:,:,:) => null()
      real*8, pointer, save :: tropopausePress_nc_overpass3    (:,:) => null()
      real*8, pointer, save :: lightningNO_nc_overpass3      (:,:,:) => null()
!
      ! Variables manipulate by worker processors only
      real*8, pointer, save :: kel_mean_overpass3            (:,:,:) => null()
      real*8, pointer, save :: psf_mean_overpass3              (:,:) => null()
      real*8, pointer, save :: qj_mean_overpass3           (:,:,:,:) => null()
      real*8, pointer, save :: qqj_mean_overpass3          (:,:,:,:) => null()
      real*8, pointer, save :: qqk_mean_overpass3          (:,:,:,:) => null()
      real*8, pointer, save :: const_mean_overpass3        (:,:,:,:) => null()
      real*8, pointer, save :: metwater_mean_overpass3       (:,:,:) => null()
      real*8, pointer, save :: totalMass_mean_overpass3      (:,:,:) => null()
      real*8, pointer, save :: relHumidity_mean_overpass3    (:,:,:) => null()
      real*8, pointer, save :: gridBoxHeight_mean_overpass3  (:,:,:) => null()
      real*8, pointer, save :: cloudOptDepth_mean_overpass3  (:,:,:) => null()
      real*8, pointer, save :: overheadO3col_mean_overpass3  (:,:,:) => null()
      real*8, pointer, save :: cloudFraction_mean_overpass3  (:,:,:) => null()
      real*8, pointer, save :: tropopausePress_mean_overpass3 (:,:) => null()
      real*8, pointer, save :: lightningNO_mean_overpass3      (:,:,:) => null()

      ! netCDF file identifier
      integer,          save :: ncid_overpass3

      ! Counter for the number of records
      integer,          save :: rnum_out_overpass3

      real*8 , pointer, save :: counter_overpass3(:)

      ! Grid information of the variables to be written out
      integer,          save :: ix1, ix2     ! longitude
      integer,          save :: i1, i2
      integer,          save :: i1_gl, i2_gl
      integer,          save :: ilo_gl, ihi_gl, ilo, ihi
      integer,          save :: jx1, jx2   ! latitude
      integer,          save :: ju1, j2
      integer,          save :: ju1_gl, j2_gl
      integer,          save :: julo_gl, jhi_gl, julo, jhi
      integer,          save :: kx1, kx2   ! vertical
      integer,          save :: k1, k2
      integer,          save :: numLon     ! number of longitudes
      integer,          save :: numLat     ! number of latitudes
      integer,          save :: numVert    ! number of vertical levels
      integer,          save :: numSpecies ! number of species

      logical,          save :: iAmRootProc
      logical,          save :: pr_diag
      integer,          save :: procID, rootProc
      integer,          save :: numDomains
      integer,          save :: commuWorld
      integer, pointer, save :: map1_u       (:,:,:) => null()
!      integer, pointer, save :: ewflag           (:) => null()
!      integer, pointer, save :: nspole           (:) => null()
      real*8 , pointer, save :: londeg(:) => null()

      integer,          save :: numQjo, numQjs, numQks

      integer,          save :: nhdf
      logical,          save :: pr_overpass3
      logical,          save :: pr_kel_overpass3
      logical,          save :: pr_psf_overpass3
      logical,          save :: pr_qj_overpass3
      logical,          save :: pr_qqjk_overpass3
      logical,          save :: pr_const_overpass3
      logical,          save :: pr_metwater_overpass3
      logical,          save :: pr_totalMass_overpass3
      logical,          save :: pr_relHumidity_overpass3
      logical,          save :: pr_gridBoxHeight_overpass3
      logical,          save :: pr_cloudOptDepth_overpass3
      logical,          save :: pr_overheadO3col_overpass3
      logical,          save :: pr_cloudFraction_overpass3
      logical,          save :: pr_tropopausePress_overpass3
      logical,          save :: pr_lightningNO_overpass3
      integer,          save :: numSpecies_overpass3
      integer,          save :: species_overpass3 (MAX_NUM_CONST_GIO)
      real*8 ,          save :: begTime_overpass3, endTime_overpass3
      real*8 ,          save :: pr_overpass3_period
      character (len=MAX_LENGTH_VAR_NAME), save :: rec_dim_name, prs_dim_name, lon_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdr_var_name, hdf_dim_name, lat_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: qj_dim_name, qqj_dim_name, qqk_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: qj_var_name, qqj_var_name, qqk_var_name
      logical           , save :: do_mean
      real*8            , save :: mdt
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
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initializeOutputOverpass3
!
! !INTERFACE:
!
      subroutine initializeOutputOverpass3(gmiGrid, gmiDomain, Diagnostics, &
     &                     Chemistry, metFields)
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiGrid  ), intent(in) :: gmiGrid
      type(t_gmiDomain), intent(in) :: gmiDomain
      type(t_metFields), intent(in) :: metFields
      type(t_Chemistry), intent(in) :: Chemistry
      type(t_Diagnostics), intent(in) :: Diagnostics
!
! !DESCRIPTION:
! The routines (1) opens a netCDF file, (2) defines variables in the file,
! (3) writes out header data in the file, and (4) allocates variables.
!
! !LOCAL VARIABLES:
      character (len=75)  :: err_msg
      integer :: in1, k1r_gl, k2r_gl
!      integer, allocatable :: ewflag_org(:)
!      integer, allocatable :: nspole_org(:)
      character (len=MAX_LENGTH_FILE_NAME) :: fname
      character (len=MAX_LENGTH_FILE_NAME) :: problem_name
!EOP
!-----------------------------------------------------------------------------
!BOC
      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) 'initializeOutputOverpass3 called by ', procID

      call Get_ilo (gmiGrid, ilo )
      call Get_ihi (gmiGrid, ihi )
      call Get_julo(gmiGrid, julo )
      call Get_jhi (gmiGrid, jhi )
      call Get_i1    (gmiGrid, i1   )
      call Get_i2    (gmiGrid, i2   )
      call Get_ju1   (gmiGrid, ju1  )
      call Get_j2    (gmiGrid, j2   )
      call Get_k1    (gmiGrid, k1   )
      call Get_k2    (gmiGrid, k2   )
      call Get_i1_gl (gmiGrid, i1_gl )
      call Get_i2_gl (gmiGrid, i2_gl )
      call Get_ju1_gl(gmiGrid, ju1_gl)
      call Get_j2_gl (gmiGrid, j2_gl )
      call Get_ilo_gl (gmiGrid, ilo_gl )
      call Get_ihi_gl (gmiGrid, ihi_gl )
      call Get_julo_gl(gmiGrid, julo_gl )
      call Get_jhi_gl (gmiGrid, jhi_gl )
      call Get_numSpecies (gmiGrid, numSpecies )

      call Get_rootProc         (gmiDomain, rootProc  )
      call Get_numDomains       (gmiDomain, numDomains )
      call Get_iAmRootProc      (gmiDomain, iAmRootProc)
      call Get_communicatorWorld(gmiDomain, commuWorld)

      allocate(londeg(i1_gl:i2_gl))
      call Get_londeg(gmiDomain, londeg)

      allocate(map1_u(2, 2, numDomains))
      call Get_map1_u (gmiDomain, map1_u)

      call Get_k1r_gl      (Diagnostics, k1r_gl)
      call Get_k2r_gl      (Diagnostics, k2r_gl)

      ix1     = i1
      ix2     = i2
      jx1     = ju1
      jx2     = j2
      kx1     = k1r_gl
      kx2     = k2r_gl

      numLon  = i2_gl - i1_gl  + 1
      numLat  = j2_gl - ju1_gl + 1
      numVert = kx2   - kx1    + 1

      call Get_num_qjs(Chemistry, numQjs)
      call Get_num_qjo(Chemistry, numQjo)
      call Get_num_qks(Chemistry, numQks)

      call Get_pr_overpass3(Diagnostics, pr_overpass3)
      call Get_pr_kel_overpass3(Diagnostics, pr_kel_overpass3)
      call Get_pr_psf_overpass3(Diagnostics, pr_psf_overpass3)
      call Get_pr_qj_overpass3(Diagnostics, pr_qj_overpass3)
      call Get_pr_qqjk_overpass3(Diagnostics, pr_qqjk_overpass3)
      call Get_pr_const_overpass3(Diagnostics, pr_const_overpass3)
      call Get_pr_metwater_overpass3(Diagnostics, pr_metwater_overpass3)
      call Get_pr_totalMass_overpass3(Diagnostics, pr_totalMass_overpass3)
      call Get_pr_relHumidity_overpass3(Diagnostics, pr_relHumidity_overpass3)
      call Get_pr_gridBoxHeight_overpass3(Diagnostics, pr_gridBoxHeight_overpass3)
      call Get_pr_cloudOptDepth_overpass3(Diagnostics, pr_cloudOptDepth_overpass3)
      call Get_pr_overheadO3col_overpass3(Diagnostics, pr_overheadO3col_overpass3)
      call Get_pr_cloudFraction_overpass3(Diagnostics, pr_cloudFraction_overpass3)
      call Get_pr_tropopausePress_overpass3(Diagnostics, pr_tropopausePress_overpass3)
      call Get_pr_lightningNO_overpass3(Diagnostics, pr_lightningNO_overpass3)
      call Get_numSpecies_overpass3(Diagnostics, numSpecies_overpass3)
      call Get_species_overpass3(Diagnostics, species_overpass3)
      call Get_begTime_overpass3(Diagnostics, begTime_overpass3)
      call Get_endTime_overpass3(Diagnostics, endTime_overpass3)
      call Get_pr_overpass3_period(Diagnostics, pr_overpass3_period)

      call Get_do_mean     (Diagnostics, do_mean)
      call Get_hdf_dim_name (Diagnostics, hdf_dim_name)
      call Get_lat_dim_name (Diagnostics, lat_dim_name)
      call Get_lon_dim_name (Diagnostics, lon_dim_name)
      call Get_hdr_var_name(Diagnostics, hdr_var_name)
      call Get_rec_dim_name(Diagnostics, rec_dim_name)
      call Get_prs_dim_name(Diagnostics, prs_dim_name)
      call Get_qj_var_name (Diagnostics, qj_var_name)
      call Get_qj_dim_name (Diagnostics, qj_dim_name)
      call Get_qqj_var_name (Diagnostics, qqj_var_name)
      call Get_qqj_dim_name (Diagnostics, qqj_dim_name)
      call Get_qqk_var_name (Diagnostics, qqk_var_name)
      call Get_qqk_dim_name (Diagnostics, qqk_dim_name)

      nhdf = NETCDF_HDF

      call Get_mdt(metFields, mdt)

      if (iAmRootProc) then

         !###########################
         ! Initialize the output file
         !###########################

         ! Determine the file name
         call Get_problem_name(Diagnostics, problem_name)
         call makeOutfileName (fname, '.overpass3.nc', problem_name)

         ! Create the file and assign a file identifier
         call Nccr_Wr (ncid_overpass3, fname)

         ! Define the variables in the file
         call ncDefineOverpass3 ()

         ! Write header data

         call ncWriteHeaderOverpass3 (gmiDomain, Chemistry, metFields)

         call Ncdo_Sync (ncid_overpass3)

         ! Initialize the counter for the number of records
         rnum_out_overpass3 = 1
      end if

      !########################
      ! Allocation of variables
      !########################

      call allocateVariablesOverpass3()

      return

      end subroutine initializeOutputOverpass3
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlOutputOverpass3
!
! !INTERFACE:
!
      subroutine controlOutputOverpass3 (last_tstp, Chemistry, &
     &                  SpeciesConcentration, Diagnostics, gmiClock, Emissions, &
     &                  metFields)
!
! !USES:
      use GmiTimeControl_mod, only : GmiSplitDateTime

      implicit none
!
! !INPUT PARAMETERS:
      logical                     , intent(in) :: last_tstp ! last time step?
      type(t_Chemistry)           , intent(in) :: Chemistry
      type(t_metFields)           , intent(in) :: metFields
      type(t_Emission)            , intent(in) :: Emissions
      type(t_gmiClock )           , intent(in) :: gmiClock
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_Diagnostics), intent(inOut) :: Diagnostics
!
! !DESCRIPTION:
! Controls the overpass3 file output. It is called at each time step.
!
! !LOCAL VARIABLES:
      logical :: time_for_overpass3
      integer :: day, month, idumyear, ndt, num_time_steps, nhms, nymd
      real*8  :: gmi_sec, tstep
      integer, save :: month_save = -999
      logical, save :: printed_on_this_day = .false.
      logical :: do_qqjk_reset
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'controlOutputOverpass3 called by ', procID

      call Get_curGmiDate(gmiClock, nymd)
      call Get_curGmitime(gmiClock, nhms)
      call Get_numTimeSteps(gmiClock, num_time_steps)
      call Get_gmiSeconds(gmiClock, gmi_sec)
      call Get_gmiTimeStep(gmiClock, tstep)

      ndt = Nint(tstep)

      call GmiSplitDateTime (nymd, idumyear, month, day)

      if (month_save == -999) month_save = month

      time_for_overpass3 = .false.

      call isOutTime (time_for_overpass3, printed_on_this_day, &
     &       month_save, month, day, nhms, ndt, gmi_sec, pr_overpass3_period)

      month_save = month

      if (last_tstp .and. (.not. time_for_overpass3)) then
        if (Mod (num_time_steps, Nint (mdt) / ndt) == 0) then

!         ------------------------------------------------------------
!         Always update restart file after the last step if you are at
!         the end of a met record.
!         ------------------------------------------------------------

          time_for_overpass3 = .true.

        end if
      end if


      if (time_for_overpass3 .or. do_mean) then
         call ncPrepOverpass3 (time_for_overpass3, Chemistry, &
     &                         SpeciesConcentration, metFields, Emissions, &
     &                         nhms)
      end if

!     ================
      if (time_for_overpass3) then
!     ================

!        ==========================
         call ncWriteOverpass3 (nymd, nhms, gmi_sec)
!        ==========================

!        ===========================
         call bufferOutOverpass3 ()
!        ===========================

         if (pr_qqjk_overpass3) then
            do_qqjk_reset = .true.
            call Set_do_qqjk_reset(Diagnostics, do_qqjk_reset)
         end if

        if (iAmRootProc) rnum_out_overpass3 = rnum_out_overpass3 + 1

!     ======
      end if
!     ======

      return

      end subroutine controlOutputOverpass3
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finalizeOutputOverpass3
! 
! !INTERFACE:
!
      subroutine finalizeOutputOverpass3()
!
      implicit none
!
! !DESCRIPTION:
! Deallocates variables necessary to produce freq1 outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'finalizeOutputOverpass3 called by ', procID

      if (pr_overpass3) then
         deallocate (const_nc_overpass3)

         if (pr_overheadO3col_overpass3) then
             deallocate (overheadO3col_nc_overpass3)
         end if

         if (pr_cloudFraction_overpass3) then
             deallocate (cloudFraction_nc_overpass3)
         end if

         if (pr_psf_overpass3) then
            deallocate (psf_nc_overpass3 )
         end if

         if (pr_tropopausePress_overpass3) then
            deallocate (tropopausePress_nc_overpass3 )
         end if

         if (pr_kel_overpass3) then
            deallocate (kel_nc_overpass3 )
         end if

         if (pr_cloudOptDepth_overpass3) then
            deallocate (cloudOptDepth_nc_overpass3 )
         end if

         if (pr_lightningNO_overpass3) then
            deallocate (lightningNO_nc_overpass3 )
         end if

         if (pr_metwater_overpass3) then
            deallocate (metwater_nc_overpass3 )
         end if

         if (pr_relHumidity_overpass3) then
            deallocate (relHumidity_nc_overpass3 )
         end if

         if (pr_totalMass_overpass3) then
            deallocate (totalMass_nc_overpass3 )
         end if

         if (pr_gridBoxHeight_overpass3) then
            deallocate (gridBoxHeight_nc_overpass3 )
         end if

         if (pr_qj_overpass3) then
            deallocate (qj_nc_overpass3 )
         end if

         if (pr_qqjk_overpass3) then
            deallocate (qqj_nc_overpass3 )
            deallocate (qqk_nc_overpass3 )
         end if

         deallocate(counter_overpass3)

         deallocate (const_overpass3)

         if (pr_qj_overpass3) then
            deallocate (qj_overpass3)
         end if

         if (pr_qqjk_overpass3) then
            deallocate (qqj_overpass3)
            deallocate (qqk_overpass3)
         end if

         if (do_mean) then

            if (pr_const_overpass3) then
               deallocate(const_mean_overpass3)
            end if

            if (pr_overheadO3col_overpass3) then
               deallocate(overheadO3col_mean_overpass3)
            end if

            if (pr_cloudFraction_overpass3) then
               deallocate(cloudFraction_mean_overpass3)
            end if

            if (pr_kel_overpass3) then
               deallocate (kel_mean_overpass3 )
            end if

            if (pr_metwater_overpass3) then
               deallocate (metwater_mean_overpass3 )
            end if

            if (pr_totalMass_overpass3) then
               deallocate (totalMass_mean_overpass3 )
            end if

            if (pr_relHumidity_overpass3) then
               deallocate (relHumidity_mean_overpass3 )
            end if

            if (pr_gridBoxHeight_overpass3) then
               deallocate (gridBoxHeight_mean_overpass3 )
            end if

            if (pr_cloudOptDepth_overpass3) then
               deallocate (cloudOptDepth_mean_overpass3 )
            end if

            if (pr_lightningNO_overpass3) then
               deallocate (lightningNO_mean_overpass3 )
            end if

            if (pr_psf_overpass3) then
               deallocate (psf_mean_overpass3 )
            end if

            if (pr_tropopausePress_overpass3) then
               deallocate (tropopausePress_mean_overpass3 )
            end if

            if (pr_qj_overpass3) then
               deallocate (qj_mean_overpass3 )
            end if

            if (pr_qqjk_overpass3) then
               deallocate (qqj_mean_overpass3 )
               deallocate (qqk_mean_overpass3 )
            end if

         end if

         if (iAmRootProc) call Nccl_Noerr (ncid_overpass3)
      end if

      return

      end  subroutine finalizeOutputOverpass3
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutputOverpass3
!
! !INTERFACE:
!
      subroutine bufferOutOverpass3 ()
!
! !USES:
      use GmiBuffer_mod, only : bufferOutput4d

      implicit none
!
! !DESCRIPTION:
! Buffers the overpass related data to the output file.
! The worker processors send a slice of the data to the root processor, 
! the root processor writes it out, they then send another slice and it 
! is written out, and so on.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_NOON2_NC = 1001
      integer, parameter :: SG_QJ2_NC   = 1002
      integer, parameter :: SG_QQJ2_NC   = 1003
      integer, parameter :: SG_QQK2_NC   = 1003
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_VAR_NAME), parameter :: COverpass_VNAM = 'const_overpass'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_const_overpass3) then
!        ====================
         call bufferOutput4d  &
!        ====================
     &     (COverpass_VNAM, do_mean, kx1, kx2, commuWorld,  &
     &      SG_NOON2_NC, ncid_overpass3, numSpecies_overpass3, rnum_out_overpass3,  &
     &      const_overpass3, const_mean_overpass3, const_nc_overpass3, &
     &      map1_u, numDomains, rootProc, procID, &
     &       i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)
      end if

      if (pr_qj_overpass3) then
         if (numQjo > 0) then

!           ====================
            call bufferOutput4d  &
!           ====================
     &        (qj_var_name, do_mean, kx1, kx2, commuWorld,  &
     &        SG_QJ2_NC, ncid_overpass3, numQjo, rnum_out_overpass3, &
     &        qj_overpass3, qj_mean_overpass3, qj_nc_overpass3, &
     &       map1_u, numDomains, rootProc, procID, &
     &       i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

         end if
      end if

      if (pr_qqjk_overpass3) then
         if (numQjs > 0) then
!           ====================
            call bufferOutput4d  &
!           ====================
     &        (qqj_var_name, do_mean, kx1, kx2, commuWorld,  &
     &         SG_QQJ2_NC, ncid_overpass3, numQjs, rnum_out_overpass3, &
     &         qqj_overpass3, qqj_mean_overpass3, qqj_nc_overpass3, &
     &      map1_u, numDomains, rootProc, procID, &
     &       i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)
          end if

         if (numQks > 0) then
!           ====================
            call bufferOutput4d  &
!           ====================
     &        (qqk_var_name, do_mean, kx1, kx2, commuWorld,  &
     &         SG_QQK2_NC, ncid_overpass3, numQks, rnum_out_overpass3, &
     &         qqk_overpass3, qqk_mean_overpass3, qqk_nc_overpass3, &
     &      map1_u, numDomains, rootProc, procID, &
     &       i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)
          end if
      end if

      return

      end subroutine bufferOutOverpass3
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocateVariablesOverpass3
! 
! !INTERFACE:
!
      subroutine allocateVariablesOverpass3()
!
      implicit none
!
! !DESCRIPTION:
! Allocates variables necessary to produce Overpass3 outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'allocateVariablesOverpass3 called by ', procID

      Allocate (const_nc_overpass3(i1:i2, ju1:j2, kx1:kx2))
      const_nc_overpass3 = 0.0d0

      if (pr_overheadO3col_overpass3) then
          Allocate (overheadO3col_nc_overpass3(i1:i2, ju1:j2, kx1:kx2))
          overheadO3col_nc_overpass3 = 0.0d0
      end if

      if (pr_cloudFraction_overpass3) then
          Allocate (cloudFraction_nc_overpass3(i1:i2, ju1:j2, kx1:kx2))
          cloudFraction_nc_overpass3 = 0.0d0
      end if

      if (pr_psf_overpass3) then
         Allocate (psf_nc_overpass3 (i1:i2, ju1:j2))
         psf_nc_overpass3 = 0.0d0
      end if

      if (pr_tropopausePress_overpass3) then
         Allocate (tropopausePress_nc_overpass3 (i1:i2, ju1:j2))
         tropopausePress_nc_overpass3 = 0.0d0
      end if

      if (pr_kel_overpass3) then
         Allocate (kel_nc_overpass3 (i1:i2, ju1:j2, kx1:kx2))
         kel_nc_overpass3 = 0.0d0
      end if

      if (pr_cloudOptDepth_overpass3) then
         Allocate (cloudOptDepth_nc_overpass3 (i1:i2, ju1:j2, kx1:kx2))
         cloudOptDepth_nc_overpass3 = 0.0d0
      end if

      if (pr_lightningNO_overpass3) then
         Allocate (lightningNO_nc_overpass3 (i1:i2, ju1:j2, kx1:kx2))
         lightningNO_nc_overpass3 = 0.0d0
      end if

      if (pr_metwater_overpass3) then
         Allocate (metwater_nc_overpass3 (i1:i2, ju1:j2, kx1:kx2))
         metwater_nc_overpass3 = 0.0d0
      end if

      if (pr_relHumidity_overpass3) then
         Allocate (relHumidity_nc_overpass3 (i1:i2, ju1:j2, kx1:kx2))
         relHumidity_nc_overpass3 = 0.0d0
      end if

      if (pr_totalMass_overpass3) then
         Allocate (totalMass_nc_overpass3 (i1:i2, ju1:j2, kx1:kx2))
         totalMass_nc_overpass3 = 0.0d0
      end if

      if (pr_gridBoxHeight_overpass3) then
         Allocate (gridBoxHeight_nc_overpass3 (i1:i2, ju1:j2, kx1:kx2))
         gridBoxHeight_nc_overpass3 = 0.0d0
      end if

      if (pr_qj_overpass3) then
         Allocate (qj_nc_overpass3 (i1:i2, ju1:j2, kx1:kx2))
         qj_nc_overpass3 = 0.0d0
      end if

      if (pr_qqjk_overpass3) then
         Allocate (qqj_nc_overpass3 (i1:i2, ju1:j2, kx1:kx2))
         qqj_nc_overpass3 = 0.0d0
         Allocate (qqk_nc_overpass3 (i1:i2, ju1:j2, kx1:kx2))
         qqk_nc_overpass3 = 0.0d0
      end if

      Allocate (const_overpass3(i1:i2, ju1:j2, kx1:kx2, numSpecies_overpass3))
      const_overpass3 = 0.0d0

      if (pr_qqjk_overpass3) then
         Allocate (qqj_overpass3(i1:i2, ju1:j2, kx1:kx2, numQjs))
         qqj_overpass3 = 0.0d0
         Allocate (qqk_overpass3(i1:i2, ju1:j2, kx1:kx2, numQks))
         qqk_overpass3 = 0.0d0
      end if

      if (pr_qj_overpass3) then
         Allocate (qj_overpass3(i1:i2, ju1:j2, kx1:kx2, numQjo))
         qj_overpass3 = 0.0d0
      end if

      allocate(counter_overpass3(i1:i2))
      counter_overpass3 = 0.0d0

      if (do_mean) then

         if (pr_const_overpass3) then
            Allocate(const_mean_overpass3(i1:i2, ju1:j2, kx1:kx2, numSpecies_overpass3))
            const_mean_overpass3 = 0.0d0
         end if

         if (pr_overheadO3col_overpass3) then
            Allocate(overheadO3col_mean_overpass3(i1:i2, ju1:j2, kx1:kx2))
            overheadO3col_mean_overpass3 = 0.0d0
         end if

         if (pr_cloudFraction_overpass3) then
            Allocate(cloudFraction_mean_overpass3(i1:i2, ju1:j2, kx1:kx2))
            cloudFraction_mean_overpass3 = 0.0d0
         end if

         if (pr_kel_overpass3) then
            Allocate (kel_mean_overpass3 (i1:i2, ju1:j2, kx1:kx2))
            kel_mean_overpass3 = 0.0d0
         end if

         if (pr_metwater_overpass3) then
            Allocate (metwater_mean_overpass3 (i1:i2, ju1:j2, kx1:kx2))
            metwater_mean_overpass3 = 0.0d0
         end if

         if (pr_totalMass_overpass3) then
            Allocate (totalMass_mean_overpass3 (i1:i2, ju1:j2, kx1:kx2))
            totalMass_mean_overpass3 = 0.0d0
         end if

         if (pr_relHumidity_overpass3) then
            Allocate (relHumidity_mean_overpass3 (i1:i2, ju1:j2, kx1:kx2))
            relHumidity_mean_overpass3 = 0.0d0
         end if

         if (pr_gridBoxHeight_overpass3) then
            Allocate (gridBoxHeight_mean_overpass3 (i1:i2, ju1:j2, kx1:kx2))
            gridBoxHeight_mean_overpass3 = 0.0d0
         end if

         if (pr_cloudOptDepth_overpass3) then
            Allocate (cloudOptDepth_mean_overpass3 (i1:i2, ju1:j2, kx1:kx2))
            cloudOptDepth_mean_overpass3 = 0.0d0
         end if

         if (pr_lightningNO_overpass3) then
            Allocate (lightningNO_mean_overpass3 (i1:i2, ju1:j2, kx1:kx2))
            lightningNO_mean_overpass3 = 0.0d0
         end if

         if (pr_psf_overpass3) then
            Allocate (psf_mean_overpass3 (i1:i2, ju1:j2))
            psf_mean_overpass3 = 0.0d0
         end if

         if (pr_tropopausePress_overpass3) then
            Allocate (tropopausePress_mean_overpass3 (i1:i2, ju1:j2))
            tropopausePress_mean_overpass3 = 0.0d0
         end if

         if (pr_qj_overpass3) then
            Allocate (qj_mean_overpass3 (i1:i2, ju1:j2, kx1:kx2, numQjo))
            qj_mean_overpass3 = 0.0d0
         end if

         if (pr_qqjk_overpass3) then
            Allocate (qqj_mean_overpass3 (i1:i2, ju1:j2, kx1:kx2, numQjs))
            qqj_mean_overpass3 = 0.0d0
            Allocate (qqk_mean_overpass3 (i1:i2, ju1:j2, kx1:kx2, numQks))
            qqk_mean_overpass3 = 0.0d0
         end if

      end if

      return

      end  subroutine allocateVariablesOverpass3
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncDefineOverpass3
!
! !INTERFACE:
!
      subroutine ncDefineOverpass3 ()

      use m_ncGeneralOpsOutput, only : Define_Netcdf_Out_Gen
 
      implicit none

#     include "netcdf.inc"
!
! !DESCRIPTION:
! Makes the necessary definitions for the overpass related variables
! in the overpass3 netCDF output file.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_VAR_NAME) :: prsp1_dim_name

      integer :: ierr
      integer :: nchr1, nchr2
      integer :: nstr
      integer :: omode
      integer :: pos1
      integer :: varid

      integer :: chrd1(1), chrd2 (1), chrd3(1)
      integer :: hdfd (1)
      integer :: lond (1), latd  (1)
      integer :: overpassd(1), oneSca(1)
      integer :: prsd (1), prsp1d(1)
      integer :: recd (1), spcd  (1)
      integer :: spcd_2(1), spcd_3(1), spcd_4(1)
      integer :: strd (1), qjsd(1)
      integer :: qqjsd(1), qqksd(1)

      integer :: var2(2)
      integer :: var3(3)
      integer :: var4(4)
      integer :: var5(5)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'ncDefineOverpass3 called by ', procID

      nchr1 = MAX_LENGTH_SPECIES_NAME
      nchr2 = 50

      nstr  =  1

!     ==========================
      call Define_Netcdf_Out_Gen  &
!     ==========================
     &  (ncid_overpass3, nhdf, numLon, numLat, numVert, hdfd, lond, latd, prsd,  &
     &   hdf_dim_name, lon_dim_name, lat_dim_name, prs_dim_name)

!     ------------------
!     Define dimensions.
!     ------------------

      prsp1_dim_name = prs_dim_name
      pos1 = Len_Trim (prs_dim_name) + 1
      prsp1_dim_name(pos1:pos1+2) = 'p1'

      call NcDef_dimension(ncid_overpass3, prsp1_dim_name, numVert+1, prsp1d(1))
!                                      --------------
      call NcDef_dimension(ncid_overpass3, 'chr_dim1', nchr1, chrd1(1))
!                                     ----------
      if (numSpecies_overpass3 > 0) then
         call NcDef_dimension  &
     &         (ncid_overpass3, 'species_overpass', numSpecies_overpass3, overpassd(1))
!                           --------------
      end if

      if (pr_qj_overpass3 .and. do_mean) then
         call NcDef_dimension(ncid_overpass3, qj_dim_name, numQjo, qjsd(1))
!                                         -----------
         call NcDef_dimension(ncid_overpass3, 'chr_dim', 80, chrd2(1))
!                                          ---------
      end if

      if (pr_qqjk_overpass3 .and. do_mean) then
         call NcDef_dimension (ncid_overpass3, qqj_dim_name, numQjs, qqjsd(1))
!                                          ------------
         call NcDef_dimension (ncid_overpass3, qqk_dim_name, numQks, qqksd(1))
!                                          ------------
         call NcDef_dimension (ncid_overpass3, 'chr_dim3', 80, chrd3(1))
!                                          ---------
      end if

      call NcDef_dimension(ncid_overpass3, 'scalar_dim', 1, oneSca(1))

      call NcDef_dimension(ncid_overpass3, rec_dim_name, NF_UNLIMITED, recd(1))
!                                      ------------
!     -----------------------------------------
!     Define variables and variable attributes.
!     -----------------------------------------
!
      call NcDef_variable (ncid_overpass3, 'am', NF_FLOAT, 1, prsd, varid)
!                                     ----
      call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', &
     &          'Hybrid pressure term')
      call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')

      call NcDef_variable (ncid_overpass3, 'bm', NF_FLOAT, 1, prsd, varid)
!                                     ----
      call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', &
     &          'Hybrid sigma term')
      call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')

      call NcDef_variable (ncid_overpass3, 'ai', NF_FLOAT, 1, prsp1d, varid)
!                                     ----
      call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', &
     &          'Hybrid pressure edge term')
      call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')

      call NcDef_variable (ncid_overpass3, 'bi', NF_FLOAT, 1, prsp1d, varid)
!                                     ----
      call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', &
     &          'Hybrid sigma edge term')
      call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')

      call NcDef_variable (ncid_overpass3, 'pt', NF_FLOAT, 0, prsp1d, varid)

         call NcDef_variable  &
     &         (ncid_overpass3, 'begTime_overpass3', NF_FLOAT, 1, oneSca, varid)
         call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', &
     &          'Beginning Time for overpass')
!
         call NcDef_variable  &
     &         (ncid_overpass3, 'endTime_overpass3', NF_FLOAT, 1, oneSca, varid)
         call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', &
     &          'End Time for overpass')

      var2(:) = (/ lond(1), latd(1) /)
      call NcDef_variable (ncid_overpass3, 'mcor', NF_FLOAT, 2, var2, varid)
!                                           ------
!ccc Overpass Species
!
      if (numSpecies_overpass3 > 0) then
         call NcDef_variable  &
     &         (ncid_overpass3, 'species_overpass', NF_FLOAT, 1, overpassd, varid)
!                           --------------
         call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', &
     &          'Overpass Species index')
         call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')
         call NcDef_var_attributes (ncid_overpass3, varid, 'coord_labels', &
     &          'overpass_labels')
         call NcDef_var_attributes (ncid_overpass3, varid, 'selection_category', 'NULL')

         var2(:) = (/ chrd1(1), overpassd(1) /)
         call NcDef_variable (ncid_overpass3, 'overpass_labels', NF_CHAR, 2, var2, varid)
!                                         -------------
         call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', &
     &          'Overpass Species name')
         call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')
         call NcDef_var_attributes (ncid_overpass3, varid, 'selection_category', 'NULL')

         var5(:) = (/ lond(1), latd(1), prsd(1), overpassd(1), recd(1) /)
         call NcDef_variable (ncid_overpass3, 'const_overpass', NF_FLOAT, 5, var5, varid)
!                                          ------------
         call NcDef_var_attributes (ncid_overpass3, varid, 'long_name',&
     &          'Constituent at overpass time range')
         call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'volume mixing ratio')
      end if

      var2(:) = (/ hdfd(1), recd(1) /)
      call NcDef_variable (ncid_overpass3, hdr_var_name, NF_INT, 2, var2, varid)
!                                      ------------
      call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'Header')
      call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'gmi_sec, nymd, nhms')
!
!cc metwater
!
      if (pr_metwater_overpass3 .and. do_mean) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_overpass3, 'metwater', NF_FLOAT, 4, var4, varid)
!                                        ----------
        call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', &
     &          'Meteorological Water')
        call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'volume mixing ratio')
      end if
!
!cc Relative humidity
!
      if (pr_relHumidity_overpass3 .and. do_mean) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_overpass3, 'relHumidity_overpass', NF_FLOAT, 4, var4, varid)
!                                        ----------
        call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'Relative Humidity')
        call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')
      end if
!
!cc metwater
!
      if (pr_totalMass_overpass3 .and. do_mean) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_overpass3, 'totalMass_overpass', NF_FLOAT, 4, var4, varid)
!                                        ----------
        call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', &
     &          'Total Mass of Atmosp.')
        call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'kg')
      end if
!
!cc grid box height
!
      if (pr_gridBoxHeight_overpass3 .and. do_mean) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_overpass3, 'gridBoxHeight_overpass', NF_FLOAT, 4, var4, varid)
!                                        ----------
        call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'Grid Box Height')
        call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'm')
      end if
!
!cc qj
!
      if (pr_qj_overpass3 .and. do_mean) then
         call NcDef_variable (ncid_overpass3, qj_dim_name, NF_FLOAT, 1, qjsd, varid)
!                                         -----------
         call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'Qjs')
         call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')
         call NcDef_var_attributes (ncid_overpass3, varid, 'coord_labels', 'qj_labels')
         call NcDef_var_attributes (ncid_overpass3, varid, 'selection_category', 'NULL')

         var2(:) = (/ chrd2(1), qjsd(1) /)
         call NcDef_variable (ncid_overpass3, 'qj_labels', NF_CHAR, 2, var2, varid)
!                                         -----------
         call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'qj name')
         call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')
         call NcDef_var_attributes (ncid_overpass3, varid, 'selection_category', 'NULL')

         var5(:) = (/ lond(1), latd(1), prsd(1), qjsd(1), recd(1) /)
         call NcDef_variable (ncid_overpass3, qj_var_name, NF_FLOAT, 5, var5, varid)
!                                         -----------
      end if
!
!cc qqjk
!
      if (pr_qqjk_overpass3 .and. do_mean) then
         call NcDef_variable (ncid_overpass3, qqj_dim_name, NF_FLOAT, 1, qqjsd, varid)
!                                         ------------
         call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'Qqjs')
         call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')
         call NcDef_var_attributes (ncid_overpass3, varid, 'coord_labels', 'qqj_labels')
         call NcDef_var_attributes (ncid_overpass3, varid, 'selection_category', 'NULL')

         call NcDef_variable (ncid_overpass3, qqk_dim_name, NF_FLOAT, 1, qqksd, varid)
!                                         ------------
         call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'Qqks')
         call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')
         call NcDef_var_attributes (ncid_overpass3, varid, 'coord_labels', 'qqk_labels')
         call NcDef_var_attributes (ncid_overpass3, varid, 'selection_category', 'NULL')

         var2(:) = (/ chrd3(1), qqjsd(1) /)
         call NcDef_variable (ncid_overpass3, 'qqj_labels', NF_CHAR, 2, var2, varid)
!                                         ------------
         call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'qqj name')
         call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')
         call NcDef_var_attributes (ncid_overpass3, varid, 'selection_category', 'NULL')

         var2(:) = (/ chrd3(1), qqksd(1) /)
         call NcDef_variable (ncid_overpass3, 'qqk_labels', NF_CHAR, 2, var2, varid)
!                                         ------------
         call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'qqk name')
         call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')
         call NcDef_var_attributes (ncid_overpass3, varid, 'selection_category', 'NULL')

         var5(:) = (/ lond(1), latd(1), prsd(1), qqjsd(1), recd(1) /)
         call NcDef_variable (ncid_overpass3, qqj_var_name, NF_FLOAT, 5, var5, varid)
!                                         ------------

         var5(:) = (/ lond(1), latd(1), prsd(1), qqksd(1), recd(1) /)
         call NcDef_variable (ncid_overpass3, qqk_var_name, NF_FLOAT, 5, var5, varid)
!                                         ------------
      end if

!     ----------------------
!     psf and kel printing
!     ---------------------
!    -----------------
      if (pr_psf_overpass3) then
        var3(:) = (/ lond(1), latd(1), recd(1) /)
        call NcDef_variable (ncid_overpass3, 'psf_overpass', NF_FLOAT, 3, var3, varid)
!                                       -----
        call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'Surface pressure')
        call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'mb')
      end if

      if (pr_kel_overpass3) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_overpass3, 'kel_overpass', NF_FLOAT, 4, var4, varid)
!                                       -----
        call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'Temperature')
        call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'K')
      end if

      if (pr_tropopausePress_overpass3) then
        var3(:) = (/ lond(1), latd(1), recd(1) /)
        call NcDef_variable (ncid_overpass3, 'tropopausePress_overpass', NF_FLOAT, 3, var3, varid)
!                                       -----
        call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'Tropopause pressure')
        call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'mb')
      end if

      if (pr_overheadO3col_overpass3) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_overpass3, 'overheadO3col_overpass', NF_FLOAT, 4, var4, varid)
!                                       -----
        call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'Photolysis overhead ozone column')
        call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'DU')
      end if

      if (pr_cloudFraction_overpass3) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_overpass3, 'cloudFraction_overpass', NF_FLOAT, 4, var4, varid)
!                                       -----
        call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', 'Total cloud fraction')
        call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')
      end if

      if (pr_cloudOptDepth_overpass3) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_overpass3, 'cloudOptDepth_overpass', NF_FLOAT, 4, var4, varid)
!                                       -----
        call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', ' Cloud optical depths (1000 nm)')
        call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'unitless')
      end if

      if (pr_lightningNO_overpass3) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_overpass3, 'lightningNO_overpass', NF_FLOAT, 4, var4, varid)
!                                       -----
        call NcDef_var_attributes (ncid_overpass3, varid, 'long_name', ' Lightning NO')
        call NcDef_var_attributes (ncid_overpass3, varid, 'units', 'kg/sec')
      end if
! ==================

!     -------------------------
!     Define global attributes.
!     -------------------------

      call NcDef_glob_attributes (ncid_overpass3, 'title',  &
     &     'Gmimod overpass variable file')

      call NcSetFill (ncid_overpass3, NF_NOFILL, omode)

      call NcEnd_def (ncid_overpass3)

      return

      end subroutine ncDefineOverpass3
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncWriteHeaderOverpass3
!
! !INTERFACE:
!
      subroutine ncWriteHeaderOverpass3 (gmiDomain, Chemistry, metFields)
!
! !USES:
      use m_ncGeneralOpsOutput, only : WriteNetcdfHdrGen

      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiDomain), intent(in) :: gmiDomain
      type(t_metFields), intent(in) :: metFields
      type(t_Chemistry), intent(in) :: Chemistry
!
! !DESCRIPTION:
! This routine creates some header information for the cloud netCDF output
! file and writes it out.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_SPECIES_NAME) :: spclab(numSpecies_overpass3)
      integer :: ic, icx
      integer :: count1d (1), count2d (2)
      integer :: start1d(1), start2d(2)
      real*8  :: prsdat(k1:k2)
      real*8  :: spcdat(numSpecies_overpass3)
      real*8 , allocatable :: qqdat(:)
      real*8 , allocatable :: mcor(:,:)
      real*8 , allocatable :: latdeg(:)
      character (len=MAX_LENGTH_LABELS) :: qj_labels(MAX_NUM_QJ)
      character (len=MAX_LENGTH_LABELS) :: qk_labels(MAX_NUM_QK)
      character (len=MAX_LENGTH_SPECIES_NAME) :: const_labels(numSpecies)      
      real*8 , allocatable :: ai(:), bi(:)
      real*8 , allocatable :: am(:), bm(:)
      real*8 :: pt
!
!EOP
!-------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'ncWriteHeaderOverpass3 called by ', procID

      call Get_qj_labels   (Chemistry, qj_labels)
      call Get_qk_labels   (Chemistry, qk_labels)
      call Get_const_labels(Chemistry, const_labels)

      allocate(latdeg(ju1_gl:j2_gl))
      call Get_latdeg(gmiDomain, latdeg)

      allocate(mcor(i1_gl:i2_gl,ju1_gl:j2_gl))
      call Get_mcorGlob(gmiDomain, mcor)

!     =========================
      call WriteNetcdfHdrGen  &
!     =========================
     &  (ncid_overpass3, latdeg, londeg, pr_diag, procID, numLon, numLat, &
     &   hdf_dim_name, lat_dim_name, lon_dim_name)

!     ---------
!     Pressure.
!     ---------

      allocate(ai(k1-1:k2))
      allocate(bi(k1-1:k2))
      allocate(am(k1:k2))
      allocate(bm(k1:k2))

      call Get_pt(metFields, pt)
      call Get_ai(metFields, ai)
      call Get_bi(metFields, bi)
      call Get_am(metFields, am)
      call Get_bm(metFields, bm)

      prsdat(:) = (am(:) * pt)  + (bm(:) * 1000.0d0)

      start1d(1) = 1
      count1d (1) = numVert

      call Ncwr_1d (prsdat(kx1:kx2), ncid_overpass3, prs_dim_name, start1d, count1d)
      call Ncwr_1d (am(kx1:kx2), ncid_overpass3, 'am', start1d, count1d)
      call Ncwr_1d (bm(kx1:kx2), ncid_overpass3, 'bm', start1d, count1d)

      count1d(1) = numVert + 1

      call Ncwr_1d (ai(kx1-1:kx2), ncid_overpass3, 'ai', start1d, count1d)
      call Ncwr_1d (bi(kx1-1:kx2), ncid_overpass3, 'bi', start1d, count1d)

      call Ncwr_Scal (pt, ncid_overpass3, 'pt')

!     --------------
!     Grid box area.
!     --------------

      start2d(:) = (/ 1, 1 /)
      count2d (:) = (/ numLon, numLat /)

      call Ncwr_2d (mcor, ncid_overpass3, 'mcor', start2d, count2d)


!     Beging and End times

      call Ncwr_Scal(begTime_overpass3, ncid_overpass3, 'begTime_overpass3')
      call Ncwr_Scal(endTime_overpass3, ncid_overpass3, 'endTime_overpass3')

!     -------------
!     Overpass species.
!     -------------

      icx = 0

      do ic = 1, numSpecies

        if (species_overpass3(ic) /= 0) then
          icx = icx + 1
          spcdat(icx) = ic
        end if

      end do

      start1d(1) = 1
      count1d (1) = numSpecies_overpass3

      if (numSpecies_overpass3 > 0) then
        call Ncwr_1d (spcdat, ncid_overpass3, 'species_overpass', start1d, count1d)
      end if

!     ------------
!     Overpass labels.
!     ------------

      icx = 0

      do ic = 1, numSpecies
        if (species_overpass3(ic) /= 0) then
          icx = icx + 1
          spclab(icx) = const_labels(ic)
        end if
      end do

      start2d(:) = (/  1, 1 /)
      count2d (:) = (/ MAX_LENGTH_SPECIES_NAME, numSpecies_overpass3 /)

      if (numSpecies_overpass3 > 0) then
         call Ncwr_2d_Char  &
     &    (spclab, ncid_overpass3, 'overpass_labels', start2d, count2d)
      end if
!
!cc qj
!
      if (pr_qj_overpass3 .and. do_mean) then

         allocate(qqdat(numQjo))
         do ic = 1, numQjo
           qqdat(ic) = ic
         end do

         start1d(1) = 1
         count1d (1) = numQjo

         call Ncwr_1d (qqdat, ncid_overpass3, qj_dim_name, start1d, count1d)
         deallocate(qqdat)

!     ----------
!     qj labels.
!     ----------

         start2d(:) = (/  1, 1 /)
         count2d (:) = (/ 80, numQjo /)

         call Ncwr_2d_Char (qj_labels, ncid_overpass3, 'qj_labels', start2d, count2d)
      end if
!
!cc qqjk
!
      if (pr_qqjk_overpass3 .and. do_mean) then

         allocate(qqdat(numQjs))
         do ic = 1, numQjs
           qqdat(ic) = ic
         end do

         start1d(1) = 1
         count1d (1) = numQjs

         call Ncwr_1d (qqdat, ncid_overpass3, qqj_dim_name, start1d, count1d)
         deallocate(qqdat)

         allocate(qqdat(numQks))
         do ic = 1, numQks
           qqdat(ic) = ic
         end do

         start1d(1) = 1
         count1d (1) = numQks

         call Ncwr_1d (qqdat, ncid_overpass3, qqk_dim_name, start1d, count1d)
         deallocate(qqdat)

!     -----------
!     qqj labels.
!     -----------

         start2d(:) = (/  1, 1 /)
         count2d (:) = (/ 80, numQjs /)

         call Ncwr_2d_Char (qj_labels, ncid_overpass3, 'qqj_labels', start2d, count2d)

!     -----------
!     qqk labels.
!     -----------

         start2d(:) = (/  1, 1 /)
         count2d (:) = (/ 80, numQks /)

         call Ncwr_2d_Char (qk_labels, ncid_overpass3, 'qqk_labels', start2d, count2d)
      end if

      return

      end subroutine ncWriteHeaderOverpass3
!EOC
!-----------------------------------------------------------------------------
!!BOP
!
! !IROUTINE: ncPrepOverpass3
!
! !INTERFACE:
!
      subroutine ncPrepOverpass3 (time_for_overpass3, Chemistry, &
     &             SpeciesConcentration, metFields, Emissions, nhms)
!
! !USES:
      use GmiChemistryMethod_mod, only : Get_qqjgmi, Get_qqkgmi, Get_qjgmi
!
      implicit none

! !INPUT PARAMETERS:
      integer                     , intent(in) :: nhms
      logical                     , intent(in) :: time_for_overpass3
      type(t_Chemistry)           , intent(in) :: Chemistry
      type(t_metFields)           , intent(in) :: metFields
      type(t_Emission)            , intent(in) :: Emissions
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !DESCRIPTION:
! This routine prepares the overpass3 variables output.
!
! !LOCAL VARIABLES:
      logical :: overpass3(i1:i2)
      integer :: ic, icx, il
      integer, save :: counter = 0
      real*8  :: local_time
      real*8  :: rhms
      real*8, allocatable :: mass(:,:,:), kel(:,:,:)
      real*8, allocatable :: tropopausePress(:,:), pctm2(:,:)
      real*8, allocatable :: relativeHumidity(:,:,:), humidity(:,:,:)
      real*8, allocatable :: gridBoxHeight(:,:,:), tau_cloud(:,:,:)
      real*8, allocatable :: overheadO3col(:,:,:), totalCloudFraction(:,:,:)
      real*8, allocatable :: lightningNO(:,:,:)
      real*8, allocatable :: qqjgmi(:,:,:,:), qqkgmi(:,:,:,:), qjgmi(:,:,:,:)
      type (t_GmiArrayBundle), pointer :: concentration(:)
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'ncPrepOverpass3 called by ', procID

      !--------------------------------------------
      ! Extract array values from the derived types
      !--------------------------------------------

      call Get_concentration(SpeciesConcentration, concentration)

      if (pr_overheadO3col_overpass3)  then
         allocate(overheadO3col(i1:i2, ju1:j2, k1:k2))
         call Get_overheadO3col(Chemistry, overheadO3col)
      end if

      if (pr_totalMass_overpass3) then
         allocate(mass(i1:i2, ju1:j2, k1:k2))
         call Get_mass(metFields, mass)
      end if

      if (pr_relHumidity_overpass3) then
         allocate(relativeHumidity(i1:i2, ju1:j2, k1:k2))
         call Get_relativeHumidity(metFields, relativeHumidity)
      end if

      if (pr_psf_overpass3) then
         allocate(pctm2(ilo:ihi, julo:jhi))
         call Get_pctm2(metFields, pctm2)
      end if

      if (pr_kel_overpass3) then
         allocate(kel(ilo:ihi, julo:jhi, k1:k2))
         call Get_kel(metFields, kel)
      end if

      if (pr_metwater_overpass3) then
         allocate(humidity(i1:i2, ju1:j2, k1:k2))
         call Get_humidity(metFields, humidity)
      end if

      if (pr_gridBoxHeight_overpass3) then
         allocate(gridBoxHeight(i1:i2, ju1:j2, k1:k2))
         call Get_gridBoxHeight(metFields, gridBoxHeight)
      end if

      if (pr_tropopausePress_overpass3) then
         allocate(tropopausePress(i1:i2, ju1:j2))
         call Get_tropopausePress(metFields, tropopausePress)
      end if

      if (pr_cloudOptDepth_overpass3) then
         allocate(tau_cloud(i1:i2, ju1:j2, k1:k2))
         call Get_tau_cloud(metFields, tau_cloud)
      end if

      if (pr_cloudFraction_overpass3) then
         allocate(totalCloudFraction(i1:i2, ju1:j2, k1:k2))
         call Get_totalCloudFraction(metFields, totalCloudFraction)
      end if

      if (pr_lightningNO_overpass3) then
         allocate(lightningNO(i1:i2, ju1:j2, k1:k2))
         call Get_lightning_no(Emissions, lightningNO)
      end if

!     -----------------------------------------------------------
!     Determine if it is local overpass3 time at each longitude and
!     increment the overpass3 counter.  This is used in several of the
!     "Prep_Netcdf" subroutine calls below.
!     -----------------------------------------------------------

      if (Any (species_overpass3(1:numSpecies) /= 0) ) then

        overpass3(:) = .false.
        do il = i1, i2
          rhms = nhms / 10000
          local_time = rhms + (londeg(il) / 15.0d0)
          if (local_time > 24.0d0) local_time = local_time - 24.0d0
          if ((local_time >= begTime_overpass3-0.1d0) .and.  &
     &        (local_time <= endTime_overpass3+0.1d0)) then
            overpass3(il) = .true.
            counter_overpass3(il) = counter_overpass3(il) + 1.0d0
          end if
        end do
      end if

!     ----------------
!     Reset mean sums.
!     ----------------

      if (counter == 0) then
        if (pr_const_overpass3) const_overpass3(:,:,:,:) = 0.0d0
        if (pr_qj_overpass3   )   qj_overpass3 (:,:,:,:) = 0.0d0
        if (pr_qqjk_overpass3 )  qqj_overpass3 (:,:,:,:) = 0.0d0
        if (pr_qqjk_overpass3 )  qqk_overpass3 (:,:,:,:) = 0.0d0
        if (do_mean) then
           if (pr_const_overpass3   )    const_mean_overpass3(:,:,:,:) = 0.0d0
           if (pr_metwater_overpass3) metwater_mean_overpass3(:,:,:  ) = 0.0d0
           if (pr_psf_overpass3     )      psf_mean_overpass3(:,:)     = 0.0d0
           if (pr_tropopausePress_overpass3) tropopausePress_mean_overpass3(:,:) = 0.0d0
           if (pr_kel_overpass3     )      kel_mean_overpass3(:,:,:)   = 0.0d0
           if (pr_totalMass_overpass3)     totalMass_mean_overpass3    (:,:,:) = 0.0d0
           if (pr_relHumidity_overpass3)   relHumidity_mean_overpass3  (:,:,:) = 0.0d0
           if (pr_gridBoxHeight_overpass3) gridBoxHeight_mean_overpass3(:,:,:) = 0.0d0
           if (pr_cloudOptDepth_overpass3) cloudOptDepth_mean_overpass3(:,:,:) = 0.0d0
           if (pr_overheadO3col_overpass3) overheadO3col_mean_overpass3(:,:,:) = 0.0d0
           if (pr_cloudFraction_overpass3) cloudFraction_mean_overpass3(:,:,:) = 0.0d0
           if (pr_lightningNO_overpass3) lightningNO_mean_overpass3(:,:,:) = 0.0d0
           if (pr_qj_overpass3      )       qj_mean_overpass3(:,:,:,:) = 0.0d0
           if (pr_qqjk_overpass3    )      qqj_mean_overpass3(:,:,:,:) = 0.0d0
           if (pr_qqjk_overpass3    )      qqk_mean_overpass3(:,:,:,:) = 0.0d0
        end if
      end if

!     -----------------
!     Update overpass3 sums.
!     -----------------

!     ----------------------------------------------
!     If local time is close to overpass3 save const_overpass3
!     ----------------------------------------------
      if (numSpecies_overpass3 > 0) then
         ic = 0
         do icx = 1, numSpecies
           if (species_overpass3(icx) /= 0) then
              ic = ic + 1
              do il = i1, i2
                if (overpass3(il)) then
                   const_overpass3(il,:,:,ic) =  const_overpass3(il,:,:,ic) + &
                           concentration(icx)%pArray3D(il,:,kx1:kx2)
                   if (do_mean) &
     &                const_mean_overpass3(il,:,:,ic) = const_mean_overpass3(il,:,:,ic) +  &
                           concentration(icx)%pArray3D(il,:,kx1:kx2)
                end if
              end do
            end if
         end do
      end if

!     ---------------------------------------
!     If local time is close to overpass3 save qj.
!     ---------------------------------------

      if (pr_qj_overpass3) then
          Allocate (qjgmi(i1:i2, ju1:j2, k1:k2, numQjo))
          call Get_qjgmi(Chemistry, qjgmi)

          do il = i1, i2
               if (overpass3(il)) then
                  qj_overpass3 (il,:,:,:) = qj_overpass3 (il,:,:,:) +          &
     &                                      qjgmi (il,:,kx1:kx2,:)
                  if (do_mean) qj_mean_overpass3 (il,:,:,:) = &
     &               qj_mean_overpass3 (il,:,:,:) + qjgmi (il,:,kx1:kx2,:)
               end if
          end do

          deallocate(qjgmi)
      end if

!     -----------------------------------------
!     If local time is close to overpass3 save qqjk.
!     -----------------------------------------

      if (pr_qqjk_overpass3) then
          Allocate (qqjgmi(i1:i2, ju1:j2, k1:k2, numQjo))
          call Get_qqjgmi(Chemistry, qqjgmi)
      
          Allocate (qqkgmi(i1:i2, ju1:j2, k1:k2, numQks))
          call Get_qqkgmi(Chemistry, qqkgmi)

          do il = i1, i2
               if (overpass3(il)) then
                  qqj_overpass3 (il,:,:,:) = qqj_overpass3 (il,:,:,:) + &
     &                                       qqjgmi (il,:,kx1:kx2,:)
                  qqk_overpass3 (il,:,:,:) = qqk_overpass3 (il,:,:,:) + &
     &                                       qqkgmi (il,:,kx1:kx2,:)
                  if (do_mean) then
                     qqj_mean_overpass3 (il,:,:,:) = qqj_mean_overpass3 (il,:,:,:) + &
     &                                       qqjgmi (il,:,kx1:kx2,:)
                     qqk_mean_overpass3 (il,:,:,:) = qqk_mean_overpass3 (il,:,:,:) + &
     &                                       qqkgmi (il,:,kx1:kx2,:)
                  end if
               end if
          end do

          deallocate(qqjgmi)
          deallocate(qqkgmi)
      end if

      if (do_mean) then

!        --------------------------------------------
!        If local time is close to overpass3 save metwater
!        --------------------------------------------

         if (pr_metwater_overpass3) then
            do il = i1, i2
               if (overpass3(il)) then
                   metwater_mean_overpass3(il,ju1:j2,:) =  &
     &               metwater_mean_overpass3(il,ju1:j2,:) +  &
     &               humidity (il,ju1:j2,kx1:kx2) * MWTAIR / (MWTH2O * GPKG)
               end if
            end do
         end if

         if (pr_relHumidity_overpass3) then
            do il = i1, i2
               if (overpass3(il)) then
                   relHumidity_mean_overpass3(il,:,:) =  &
     &               relHumidity_mean_overpass3(il,:,:) + relativeHumidity (il,:,kx1:kx2)
               end if
            end do
         end if

         if (pr_totalMass_overpass3) then
            do il = i1, i2
               if (overpass3(il)) then
                   totalMass_mean_overpass3(il,:,:) =  &
     &               totalMass_mean_overpass3(il,:,:) + mass (il,:,kx1:kx2)
               end if
            end do
         end if

         if (pr_gridBoxHeight_overpass3) then
            do il = i1, i2
               if (overpass3(il)) then
                   gridBoxHeight_mean_overpass3(il,:,:) =  &
     &               gridBoxHeight_mean_overpass3(il,:,:) + gridBoxHeight (il,:,kx1:kx2) 
               end if
            end do
         end if

        if (pr_psf_overpass3)  then
          do il=i1,i2
             if (overpass3(il)) then
                   psf_mean_overpass3 (il,ju1:j2)   =  &
     &                psf_mean_overpass3 (il,ju1:j2)   + pctm2(il,ju1:j2)
             end if
          end do
        end if

        if (pr_tropopausePress_overpass3)  then
          do il=i1,i2
             if (overpass3(il)) then
                   tropopausePress_mean_overpass3 (il,ju1:j2)   =  &
     &                tropopausePress_mean_overpass3 (il,ju1:j2) + tropopausePress(il,ju1:j2)
             end if
          end do
        end if

        if (pr_kel_overpass3)  then
          do il=i1,i2
             if (overpass3(il)) then
                   kel_mean_overpass3 (il,ju1:j2,:) =  &
     &                kel_mean_overpass3 (il,ju1:j2,:) + kel  (il,ju1:j2,kx1:kx2)
             end if
          end do
        end if

        if (pr_cloudOptDepth_overpass3)  then
          do il=i1,i2
             if (overpass3(il)) then
                   cloudOptDepth_mean_overpass3 (il,ju1:j2,:) =  &
     &                cloudOptDepth_mean_overpass3 (il,ju1:j2,:) + tau_cloud (il,ju1:j2,kx1:kx2)
             end if
          end do
        end if

        if (pr_overheadO3col_overpass3)  then
           do il=i1,i2
              if (overpass3(il)) then
                 overheadO3col_mean_overpass3 (il,ju1:j2,:) =  &
     &                overheadO3col_mean_overpass3 (il,ju1:j2,:) + &
     &                overheadO3col(il,ju1:j2,kx1:kx2)
              end if
           end do
        end if

        if (pr_cloudFraction_overpass3)  then
           do il=i1,i2
              if (overpass3(il)) then
                 cloudFraction_mean_overpass3 (il,ju1:j2,:) =  &
     &                cloudFraction_mean_overpass3 (il,ju1:j2,:) + &
     &                totalCloudFraction(il,ju1:j2,kx1:kx2)
              end if
           end do
        end if

        if (pr_lightningNO_overpass3)  then
           do il=i1,i2
              if (overpass3(il)) then
                 lightningNO_mean_overpass3 (il,ju1:j2,:) =  &
     &                lightningNO_mean_overpass3 (il,ju1:j2,:) + &
     &                lightningNO(il,ju1:j2,kx1:kx2)
              end if
           end do
        end if

      end if

      counter = counter + 1

      if (time_for_overpass3) then

        if (do_mean) then
!         ---------------------
!         Calculate overpass means.
!         ---------------------

          if (numSpecies_overpass3 > 0) then
             do ic = 1, numSpecies_overpass3
               do il = i1, i2
                 if (counter_overpass3(il) > 0) then
                   const_mean_overpass3(il,:,:,ic) =  &
     &                const_mean_overpass3(il,:,:,ic) / counter_overpass3(il)
                 end if
               end do
             end do
          end if

          if (pr_qj_overpass3) then
             do il = i1, i2
                if (counter_overpass3(il) /= 0.d0) then
                   qj_mean_overpass3(il,:,:,:) =  qj_mean_overpass3(il,:,:,:) / &
     &                                            counter_overpass3(il)
                end if
            end do
          end if

          if (pr_qqjk_overpass3) then
             do il = i1, i2
                if (counter_overpass3(il) /= 0.d0) then
                   qqj_mean_overpass3(il,:,:,:) = qqj_mean_overpass3(il,:,:,:) / &
     &                                            counter_overpass3(il)
                   qqk_mean_overpass3(il,:,:,:) = qqk_mean_overpass3(il,:,:,:) / &
     &                                            counter_overpass3(il)
                end if
             end do
          end if

          do il = i1, i2
             if (counter_overpass3(il) > 0) then
                if (pr_psf_overpass3)  psf_mean_overpass3 (il,:) &
     &                  = psf_mean_overpass3 (il,:) /counter_overpass3(il) 
                if (pr_kel_overpass3)  kel_mean_overpass3 (il,:,:) &
     &                  = kel_mean_overpass3 (il,:,:) /counter_overpass3(il) 
                if (pr_cloudOptDepth_overpass3)  cloudOptDepth_mean_overpass3 (il,:,:) &
     &                  = cloudOptDepth_mean_overpass3 (il,:,:) /counter_overpass3(il) 
                if (pr_overheadO3col_overpass3) overheadO3col_mean_overpass3 (il,:,:) &
     &                  = overheadO3col_mean_overpass3 (il,:,:) /counter_overpass3(il) 
                if (pr_cloudFraction_overpass3) cloudFraction_mean_overpass3 (il,:,:) &
     &                  = cloudFraction_mean_overpass3 (il,:,:) /counter_overpass3(il) 
                if (pr_tropopausePress_overpass3)  tropopausePress_mean_overpass3 (il,:) &
     &                  = tropopausePress_mean_overpass3 (il,:) /counter_overpass3(il) 
                if (pr_lightningNO_overpass3) lightningNO_mean_overpass3 (il,:,:) &
     &                  = lightningNO_mean_overpass3 (il,:,:) /counter_overpass3(il) 
             end if
          end do
          psf_nc_overpass3(:,:) = psf_mean_overpass3(:,:)
          kel_nc_overpass3(:,:,:) = kel_mean_overpass3(:,:,:)
          cloudOptDepth_nc_overpass3(:,:,:) = cloudOptDepth_mean_overpass3(:,:,:)
          cloudFraction_nc_overpass3(:,:,:) = cloudFraction_mean_overpass3(:,:,:)
          tropopausePress_nc_overpass3(:,:) = tropopausePress_mean_overpass3(:,:)
          lightningNO_nc_overpass3(:,:,:) = lightningNO_mean_overpass3(:,:,:)
          ! Converting from mol/cm^2 to Dobson unit (DU)
          overheadO3col_mean_overpass3(:,:,:) = overheadO3col_mean_overpass3(:,:,:)/2.690d+16
          overheadO3col_nc_overpass3  (:,:,:) = overheadO3col_mean_overpass3(:,:,:)

          if (pr_metwater_overpass3) then
             do il = i1, i2
                if (counter_overpass3(il) > 0) then
                   metwater_mean_overpass3(il,:,:) =  &
     &              metwater_mean_overpass3(il,:,:) / counter_overpass3(il)
                end if

             end do
             metwater_nc_overpass3(:,:,:) = metwater_mean_overpass3(:,:,:)
          end if

          if (pr_relHumidity_overpass3) then
             do il = i1, i2
                if (counter_overpass3(il) > 0) then
                   relHumidity_mean_overpass3(il,:,:) =  &
     &              relHumidity_mean_overpass3(il,:,:) / counter_overpass3(il)
                end if

             end do
             relHumidity_nc_overpass3(:,:,:) = relHumidity_mean_overpass3(:,:,:)
          end if

          if (pr_totalMass_overpass3) then
             do il = i1, i2
                if (counter_overpass3(il) > 0) then
                   totalMass_mean_overpass3(il,:,:) =  &
     &              totalMass_mean_overpass3(il,:,:) / counter_overpass3(il)
                end if

             end do
             totalMass_nc_overpass3(:,:,:) = totalMass_mean_overpass3(:,:,:)
          end if

          if (pr_gridBoxHeight_overpass3) then
             do il = i1, i2
                if (counter_overpass3(il) > 0) then
                   gridBoxHeight_mean_overpass3(il,:,:) =  &
     &              gridBoxHeight_mean_overpass3(il,:,:) / counter_overpass3(il)
                end if

             end do
             gridBoxHeight_nc_overpass3(:,:,:) = gridBoxHeight_mean_overpass3(:,:,:)
          end if

        end if
        counter = 0
        counter_overpass3(:) = 0.0d0
      end if

      if (pr_overheadO3col_overpass3) deallocate(overheadO3col) 
      if (pr_totalMass_overpass3) deallocate(mass) 
      if (pr_relHumidity_overpass3) deallocate(relativeHumidity) 
      if (pr_psf_overpass3) deallocate(pctm2) 
      if (pr_kel_overpass3) deallocate(kel) 
      if (pr_metwater_overpass3) deallocate(humidity) 
      if (pr_gridBoxHeight_overpass3) deallocate(gridBoxHeight) 
      if (pr_tropopausePress_overpass3) deallocate(tropopausePress)

      return

      end subroutine ncPrepOverpass3
!EOC
!-----------------------------------------------------------------------------
!BOP
! 
! !IROUTINE: ncWriteOverpass3
!
! !INTERFACE:
!
      subroutine ncWriteOverpass3 (nymd, nhms, gmi_sec)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: nymd, nhms
      real*8 , intent(in) :: gmi_sec
!
! !DESCRIPTION:
! This routine writes out variables that have at most 3 dimensions.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_METWATER_NC = 1010 
      integer, parameter :: SG_NOON2_NC    = 1011 
      integer, parameter :: SG_QJ1_NC      = 1012
      integer, parameter :: SG_QQJ1_NC     = 1013
      integer, parameter :: SG_QQK1_NC     = 1014
      integer, parameter :: SG_KEL_NOON_NC = 1015
      integer, parameter :: SG_PSF_NOON_NC = 1016
      integer, parameter :: SG_NOON1_NC    = 1017
!
! !LOCAL VARIABLES:
      integer :: count2d(2), count3d(3), count4d(4)
      integer :: start2d(2), start3d(3), start4d(4)
      integer :: hdr(NETCDF_HDF)
      real*8, allocatable :: arrayGlob2D(:,:)
      real*8, allocatable :: arrayGlob3D(:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'ncWriteOverpass3 called by ', procID

      if (iAmRootProc) then
         start2d(:) = (/ 1, rnum_out_overpass3 /)
         start3d(:) = (/ 1, 1, rnum_out_overpass3 /)
         start4d(:) = (/ 1, 1, 1, rnum_out_overpass3 /)

         count2d (:) = (/ NETCDF_HDF, 1 /)
         count3d (:) = (/ numLon, numLat, 1 /)
         count4d (:) = (/ numLon, numLat, numVert,  1 /)

         hdr(1) = Nint (gmi_sec)
         hdr(2) = nymd
         hdr(3) = nhms

         call Ncwr_2d_Int (hdr, ncid_overpass3, hdr_var_name, start2d, count2d)

         allocate(arrayGlob2D(numLon, numLat))
         allocate(arrayGlob3D(numLon, numLat, numVert))
      end if

      if (pr_tropopausePress_overpass3) then
         call subDomain2Global (arrayGlob2D, tropopausePress_nc_overpass3, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2,  &
     &           rootProc, procID, map1_u, numDomains, SG_PSF_NOON_NC, commuWorld)

         if (iAmRootProc)  then
            call Ncwr_3d (arrayGlob2D, ncid_overpass3, &
     &               'tropopausePress_overpass', start3d, count3d)

            arrayGlob2D = 0.0d0
         end if
      end if

      if (pr_psf_overpass3) then
         call subDomain2Global (arrayGlob2D, psf_nc_overpass3, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, &
     &           rootProc, procID, map1_u, numDomains, SG_PSF_NOON_NC, commuWorld)

         if (iAmRootProc)  then
            call Ncwr_3d (arrayGlob2D, ncid_overpass3, 'psf_overpass', &
     &                  start3d, count3d)
         end if
      end if


      if (pr_totalMass_overpass3) then
         call subDomain2Global (arrayGlob3D, totalMass_nc_overpass3, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, SG_METWATER_NC, commuWorld)

         if (iAmRootProc) then
            call Ncwr_4d (arrayGlob3D, ncid_overpass3, &
     &                    'totalMass_overpass', start4d, count4d)
            arrayGlob3D = 0.0d0
         end if
      end if

      if (pr_overheadO3col_overpass3) then
         call subDomain2Global (arrayGlob3D, overheadO3col_nc_overpass3, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &           rootProc, procID, map1_u, numDomains, SG_METWATER_NC, commuWorld)

         if (iAmRootProc) then
            call Ncwr_4d (arrayGlob3D, ncid_overpass3, &
     &                    'overheadO3col_overpass', start4d, count4d)

            arrayGlob3D = 0.0d0
         end if
      end if

      if (pr_cloudFraction_overpass3) then
         call subDomain2Global (arrayGlob3D, cloudFraction_nc_overpass3, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, SG_METWATER_NC, commuWorld)

         if (iAmRootProc) then
            call Ncwr_4d (arrayGlob3D, ncid_overpass3, &
     &                    'cloudFraction_overpass', start4d, count4d)

            arrayGlob3D = 0.0d0
         end if
      end if

      if (pr_lightningNO_overpass3) then
         call subDomain2Global (arrayGlob3D, lightningNO_nc_overpass3, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, SG_METWATER_NC, commuWorld)

         if (iAmRootProc) then
            call Ncwr_4d (arrayGlob3D, ncid_overpass3, &
     &                    'lightningNO_overpass', start4d, count4d)

            arrayGlob3D = 0.0d0
         end if
      end if


      if (pr_relHumidity_overpass3) then
         call subDomain2Global (arrayGlob3D, relHumidity_nc_overpass3, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, SG_METWATER_NC, commuWorld)

         if (iAmRootProc) then
            call Ncwr_4d (arrayGlob3D, ncid_overpass3, &
     &                    'relHumidity_overpass', start4d, count4d)

            arrayGlob3D = 0.0d0
         end if
      end if

      if (pr_gridBoxHeight_overpass3) then
         call subDomain2Global (arrayGlob3D, gridBoxHeight_nc_overpass3, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, SG_METWATER_NC, commuWorld)

         if (iAmRootProc) then
            call Ncwr_4d (arrayGlob3D, ncid_overpass3, &
     &                 'gridBoxHeight_overpass', start4d, count4d)

            arrayGlob3D = 0.0d0
         end if
      end if

      if (pr_metwater_overpass3) then
         call subDomain2Global (arrayGlob3D, metwater_nc_overpass3, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, SG_METWATER_NC, commuWorld)

         if (iAmRootProc) then
            call Ncwr_4d (arrayGlob3D, ncid_overpass3, &
     &                 'metwater', start4d, count4d)

            arrayGlob3D = 0.0d0
         end if
      end if

      if (pr_kel_overpass3) then
         call subDomain2Global (arrayGlob3D, kel_nc_overpass3, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, SG_KEL_NOON_NC, commuWorld)

         if (iAmRootProc) then
            call Ncwr_4d (arrayGlob3D, ncid_overpass3, 'kel_overpass', &
     &               start4d, count4d)
            arrayGlob3D = 0.0d0
         end if
      end if

      if (pr_cloudOptDepth_overpass3) then
         call subDomain2Global (arrayGlob3D, cloudOptDepth_nc_overpass3, &
     &           i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &           rootProc, procID, map1_u, numDomains, SG_KEL_NOON_NC, commuWorld)

         if (iAmRootProc) then
            call Ncwr_4d (arrayGlob3D, ncid_overpass3, &
     &                'cloudOptDepth_overpass', start4d, count4d)
         end if
      end if

      if (iAmRootProc) then 
         call Ncdo_Sync (ncid_overpass3)
      end if

      return

      end subroutine ncWriteOverpass3
!EOC
!-------------------------------------------------------------------------------
end module GmiControlOverpass3_mod
