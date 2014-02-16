!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlCloudModuleGT_mod
!
! !INTERFACE:
!
module GmiControlCloudModuleGT_mod
!
! !USES:
      use  GmiSub2Glob_mod      , only : subDomain2Global
      use  GmiMessagePassing_mod, only : synchronizeGroup
      use m_ncGeneralOpsOutput  , only: Is_Out_Freq_Time
      use GmiFileOperations_mod , only : makeOutfileName
      use m_netcdf_io_close     , only : Nccl, Nccl_Noerr
      use m_netcdf_io_write     , only : Ncwr_2d_Int, Ncwr_5d, Ncwr_3d, Ncwr_2d
      use m_netcdf_io_write     , only : Ncwr_1d, Ncwr_Scal, Ncwr_2d_Char
      use m_netcdf_io_create    , only : Ncdo_Sync, Nccr_Wr
      use m_netcdf_io_define    , only : NcDef_glob_attributes, NcDef_dimension
      use m_netcdf_io_define    , only : NcDef_var_attributes, NcDef_variable
      use m_netcdf_io_define    , only : NcSetFill, NcEnd_def
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration
      use GmiChemistryMethod_mod , only : t_Chemistry, Get_overheadO3col
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      use GmiDomainDecomposition_mod, only : t_gmiDomain
      use GmiDomainDecomposition_mod, only : Get_rootProc, Get_procID, Get_iAmRootProc
      use GmiDomainDecomposition_mod, only : Get_numDomains
      use GmiDomainDecomposition_mod, only : Get_communicatorWorld
      use GmiDomainDecomposition_mod, only : Get_map1_u
!      use GmiDomainDecomposition_mod, only : Get_eastWestFlag , Get_northSouthPole
      use GmiDomainDecomposition_mod, only : Get_mcorGlob, Get_latdeg, Get_londeg

      use GmiGrid_mod              , only : t_gmiGrid, Get_numSpecies
      use GmiGrid_mod              , only : Get_i1, Get_i2, Get_ju1, Get_j2, Get_k1, Get_k2
      use GmiGrid_mod              , only : Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl
      use GmiGrid_mod              , only : Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jhi_gl
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics
      use GmiDiagnosticsMethod_mod, only : Get_cloudOutputFrequency, Get_hdr_var_name, Get_rec_dim_name
      use GmiDiagnosticsMethod_mod, only : Get_lat_dim_name, Get_hdf_dim_name, Get_lon_dim_name
      use GmiDiagnosticsMethod_mod, only : Get_prs_dim_name, Get_pr_grid_height
      use GmiDiagnosticsMethod_mod, only : Get_do_mean, Get_k1r_gl, Get_k2r_gl
      use GmiDiagnosticsMethod_mod, only : Get_pr_cloud, Get_pr_diag, Get_problem_name
      use GmiDiagnosticsMethod_mod, only : Get_pr_psf, Get_pr_kel, Get_pr_mass

      use GmiTimeControl_mod, only : t_GmiClock, GmiSplitDateTime, isOutTime, &
     &       Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds, Get_gmiTimeStep, &
     &       Get_numTimeSteps

      use GmiMetFieldsControl_mod, only : t_metFields, Get_pt,   &
     &       Get_metdata_name, Get_mdt, Get_ai, Get_bi, Get_am, Get_bm, &
     &       Get_mass, Get_gridBoxHeight, Get_kel, Get_pctm2, Get_lwi_flags
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: controlOutputCloudGT, initializeOutputCloudGT
      public  :: finalizeOutputCloudGT
!
#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
!
! !DESCRIPTION:
! Contains the necessary routines for the cloud outputs. Three routines
! are visible from the outside:
! %
! \begin{description}
! \item[initializeOutputCloudGT:] called once in the initialization steps to
!      open the netCDF file, define variables in the file, write out header
!      information and allocate variables.
! \item[controlOutputCloudGT:] called at the end of each model time step to
!      update variables and if needed for communications and writing out
!      data in the file.
! \item[finalizeOutputCloudGT:] called at the end of the integration to
!      deallocate all the variables and close the netCDF file.
! \end{description}
! %
! The above routines have as arguments derived types.
! %
! This module is self-contained. All the main routines are called by all
! the processors (master included) but each processor does not execute
! all portions of the code.
!
! DEFINED PARAMETERS:
      integer, parameter :: num_cloud = 22
!
      real*8 , pointer, save :: psf_nc             (:,:) => null()
      real*8 , pointer, save :: kel_nc           (:,:,:) => null()
      real*8 , pointer, save :: mass_nc          (:,:,:) => null()
      real*8 , pointer, save :: gridbox_height_nc(:,:,:) => null()
      real*8 , pointer, save :: cloud_data     (:,:,:,:) => null()
      real*8 , pointer, save :: cloud_data_nc    (:,:,:) => null()
      real*8 , pointer, save :: cloud_rh_nc      (:,:,:) => null() ! Relative Humidity
      real*8 , pointer, save :: cloud_cf_nc      (:,:,:) => null() ! Cloud Fraction
      real*8 , pointer, save :: cloud_lwc_nc     (:,:,:) => null() ! Liquid Water Path
      real*8 , pointer, save :: cloud_r_nc       (:,:,:) => null() ! Effective Radii
      real*8 , pointer, save :: cloud_cod_nc     (:,:,:) => null() ! Cloud Optical Depth
      real*8 , pointer, save :: cloud_cod_n_nc   (:,:,:) => null() ! Cloud Optical Depth from SO4
      real*8 , pointer, save :: cloud_so4_nc     (:,:,:) => null() ! Total Aerosol Sulfate
      real*8 , pointer, save :: cloud_cdnc_nc    (:,:,:) => null() ! Cloud Droplet # Concentration
      real*8 , pointer, save :: cloud_effr_nc    (:,:,:) => null() ! Effective Radii from SO4
      real*8 , pointer, save :: cloud_alb_nc     (:,:,:) => null() ! Cloud Albedo
      real*8 , pointer, save :: cloud_lwp_nc     (:,:,:) => null() ! Liquid Water Path
      real*8 , pointer, save :: cloud_alb_par_nc (:,:,:) => null() ! CA from Fountoukis & Nenes, 2005
      real*8 , pointer, save :: cloud_effr_par_nc(:,:,:) => null() ! E R from Fountoukis & Nenes, 2005
      real*8 , pointer, save :: cloud_cdnc_par_nc(:,:,:) => null() ! CDNC from Fountoukis & Nenes, 2005
      real*8 , pointer, save :: cloud_cod_par_nc (:,:,:) => null() ! COD from Fountoukis & Nenes, 2005
      real*8 , pointer, save :: cloud_smax_par_nc(:,:,:) => null() ! Max Super Fountoukis & Nenes, 2005
      real*8 , pointer, save :: cloud_w_par_nc   (:,:,:) => null() ! Updraft Fountoukis & Nenes, 2005
      real*8 , pointer, save :: cloud_AIE_nc     (:,:,:) => null() ! Flux gt
      real*8 , pointer, save :: cloud_AIE12_nc   (:,:,:) => null() ! Flux gt 1_2
      real*8 , pointer, save :: Qaut_NS_KK1      (:,:,:) => null() ! Autoconversion rate K+K Eq. 29
      real*8 , pointer, save :: Qaut_NS_KK2      (:,:,:) => null() ! Autoconversion rate K+K Eq. 30
      real*8 , pointer, save :: Qaut_NS_R        (:,:,:) => null() ! Autoconversion rate R Eq. 15

      integer,          save :: ncid_cloud           ! NetCDF       species output file id
      integer,          save :: rnum_out_cloud       ! next Netcdf output record to write

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

      integer,          save :: nhdf

      character (len=MAX_LENGTH_VAR_NAME), save :: cloud_var_name ! cloud variable name
      character (len=80)       :: cloud_labels(num_cloud)
      character (len=MAX_LENGTH_VAR_NAME), save :: rec_dim_name, prs_dim_name, lon_dim_name
      character (len=MAX_LENGTH_VAR_NAME), save :: hdr_var_name, hdf_dim_name, lat_dim_name
      logical           , save :: pr_cloud, do_mean, grid_height, p
      logical           , save :: pr_psf, pr_kel, pr_grid_height, pr_mass
      real*8            , save :: cloudOutputFrequency
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
!-------------------------------------------------------------------------------
!BOP  
!
! !IROUTINE: initializeOutputCloudGT
!     
! !INTERFACE:
!     
      subroutine initializeOutputCloudGT(gmiGrid, gmiDomain, Diagnostics, &
     &                     metFields)
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiGrid    ), intent(in) :: gmiGrid  
      type(t_gmiDomain  ), intent(in) :: gmiDomain
      type(t_metFields  ), intent(in) :: metFields
      type(t_Diagnostics), intent(in) :: Diagnostics
!     
! !DESCRIPTION:
! The routines (1) opens a netCDF file, (2) defines variables in the file,
! (3) writes out header data in the file, and (4) allocates varaibles.
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

      if (pr_diag) Write (6,*) 'initializeOutputCloudGT called by ', procID

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
      call Get_iAmRootProc      (gmiDomain, iAmRootProc )
      call Get_communicatorWorld(gmiDomain, commuWorld)

!      allocate(ewflag (i1:i2))
!      allocate(ewflag_org (ilo_gl:ihi_gl))
!      call Get_eastWestFlag(gmiDomain, ewflag_org)
!      ewflag(i1:i2) = ewflag_org(i1:i2)
!      deallocate(ewflag_org)

!      allocate(nspole(ju1:j2))
!      allocate(nspole_org(julo_gl:jhi_gl))
!      call Get_northSouthPole(gmiDomain, nspole_org)
!      nspole(ju1:j2) = nspole_org(ju1:j2)
!      deallocate(nspole_org)

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

      numLon  = i2_gl - i1_gl + 1
      numLat  = j2_gl - ju1_gl + 1
      numVert = kx2 - kx1 + 1

      nhdf = NETCDF_HDF

      call Get_do_mean     (Diagnostics, do_mean)
      call Get_cloudOutputFrequency(Diagnostics, cloudOutputFrequency)
      call Get_pr_cloud    (Diagnostics, pr_cloud)
      call Get_pr_kel      (Diagnostics, pr_kel  )
      call Get_pr_psf      (Diagnostics, pr_psf  )
      call Get_pr_mass     (Diagnostics, pr_mass )
      call Get_pr_grid_height(Diagnostics, pr_grid_height )
      call Get_hdr_var_name(Diagnostics, hdr_var_name)
      call Get_rec_dim_name(Diagnostics, rec_dim_name)
      call Get_prs_dim_name(Diagnostics, prs_dim_name)
      call Get_hdf_dim_name (Diagnostics, hdf_dim_name)
      call Get_lat_dim_name (Diagnostics, lat_dim_name)
      call Get_lon_dim_name (Diagnostics, lon_dim_name)

      call Get_mdt (metFields, mdt)

      if (pr_cloud) then

         if (iAmRootProc) then

            !###########################
            ! Initialize the output file
            !###########################

            ! Determine the file name
            call Get_problem_name(Diagnostics, problem_name)
            call makeOutfileName (fname, '.cloud.nc', problem_name)

            ! Create the file and assign a file identifier
            call Nccr_Wr (ncid_cloud, fname)

            ! Define the variables in the file
            call Define_Netcdf_Out_CloudGT()

            ! Write header data
            call Write_Netcdf_Hdr_CloudGT (gmiDomain, metFields)

            call Ncdo_Sync (ncid_cloud)

            ! Initialize the counter for the number of records
            rnum_out_cloud = 1

            cloud_labels(:) = ' '
            cloud_labels(1)  = 'Relative Humidity'
            cloud_labels(2)  = 'Cloud Fraction '
            cloud_labels(3)  = 'Liquid Water Content'
            cloud_labels(4)  = 'Effective Radius'
            cloud_labels(5)  = 'Cloud Optical Depth'
            cloud_labels(6)  = 'Total Aerosol SO4'
            cloud_labels(7)  = 'Cloud Droplet Number Concentration'
            cloud_labels(8)  = 'Droplet Effective Radius'
            cloud_labels(9)  = 'Cloud Albedo'
            cloud_labels(10) = 'Cloud Optical Depth from SO4'
            cloud_labels(11) = 'Cloud Albedo from SO4'

         end if

         cloud_var_name = 'cloud_data'

         !########################
         ! Allocation of variables
         !########################

         call allocateVariablesCloudGT()

      end if

      return

      end subroutine initializeOutputCloudGT
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: controlOutCloudGT
!
! !INTERFACE:
!
      subroutine controlOutputCloudGT (last_tstp, Chemistry, gmiClock, metFields)
!
! !USES:
      use GmiTimeControl_mod    , only : GmiSplitDateTime
      use GmiChemistryMethod_mod, only : t_Chemistry
!
      implicit none
!
! !INPUT PARAMETERS:
      logical          , intent(in) :: last_tstp
      type(t_gmiClock ), intent(in) :: gmiClock 
      type(t_Chemistry), intent(in) :: Chemistry
      type(t_metFields), intent(in) :: metFields
!
! !DESCRIPTION:
! This routine controls the cloud file output. It is called at each time step.
!
! !LOCAL VARIABLES:
      logical :: time_for_cloud
      integer :: day, month, idumyear, ndt, num_time_steps, nhms, nymd
      real*8  :: gmi_sec, tstep
      integer, save :: month_save = -999
      logical, save :: printed_on_this_day = .false.
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'controlOutputCloudGT called by ', procID

      call Get_curGmiDate(gmiClock, nymd)
      call Get_curGmitime(gmiClock, nhms)
      call Get_numTimeSteps(gmiClock, num_time_steps)
      call Get_gmiSeconds(gmiClock, gmi_sec)
      call Get_gmiTimeStep(gmiClock, tstep)

      ndt = Nint(tstep)

      call GmiSplitDateTime (nymd, idumyear, month, day)

      if (month_save == -999) month_save = month

      time_for_cloud = .false.

      call isOutTime(time_for_cloud, printed_on_this_day, &
     &        month_save, month, day, nhms, ndt, gmi_sec, cloudOutputFrequency)

      month_save = month

      if (last_tstp .and. (.not. time_for_cloud)) then
         if (Mod (num_time_steps, Nint (mdt) / ndt) == 0) then

!         ------------------------------------------------------------
!         Always update restart file after the last step if you are at
!         the end of a met record.
!         ------------------------------------------------------------

           time_for_cloud = .true.
         end if
      end if

      if (time_for_cloud) then
         !====================
         call Prep_Netcdf_CloudGT (time_for_cloud, Chemistry, metFields)
         !====================
      endif

!     ================
      if (time_for_cloud) then
!     ================

!        =====================
         call Write_Netcdf_CloudGT (nymd, nhms, gmi_sec)
!        =====================

!        ======================
         call bufferOutput_CloudGT ( )
!        ======================

        if (iAmRootProc) rnum_out_cloud = rnum_out_cloud + 1

!     ======
      end if
!     ======

      return

      end subroutine controlOutputCloudGT
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finalizeOutputCloudGT
!
! !INTERFACE:
!
      subroutine finalizeOutputCloudGT( )
!
      implicit none
!
! !DESCRIPTION:
! Deallocates variables necessary to produce cloud outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if (pr_diag) Write(6,*) 'finalizeOutputCloudGT called by ', procID

      !==========================
      ! Deallocation of variables
      !==========================

      deallocate (cloud_data        )
      deallocate (cloud_data_nc     )
      deallocate (cloud_rh_nc       )
      deallocate (cloud_cf_nc       )
      deallocate (cloud_lwc_nc      )
      deallocate (cloud_r_nc        )
      deallocate (cloud_cod_nc      )
      deallocate (cloud_cod_n_nc    )
      deallocate (cloud_so4_nc      )
      deallocate (cloud_cdnc_nc     )
      deallocate (cloud_effr_nc     )
      deallocate (cloud_alb_nc      )
      deallocate (cloud_lwp_nc      )
      deallocate (cloud_effr_par_nc )
      deallocate (cloud_alb_par_nc  )
      deallocate (cloud_cod_par_nc  )
      deallocate (cloud_cdnc_par_nc )
      deallocate (cloud_w_par_nc    )
      deallocate (cloud_smax_par_nc )
      deallocate (Qaut_NS_KK1       )
      deallocate (Qaut_NS_KK2       )
      deallocate (Qaut_NS_R         )
      deallocate (cloud_AIE_nc      )
      deallocate (cloud_AIE12_nc    )

      !======================
      ! Close the netCDF file
      !======================
      if (iAmRootProc) call Nccl_Noerr (ncid_cloud)

      return

      end subroutine finalizeOutputCloudGT
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define_Netcdf_Out_CloudGT
!
! !INTERFACES:
!
      subroutine Define_Netcdf_Out_CloudGT  ()
!
! !USES:
      use m_ncGeneralOpsOutput, only: Define_Netcdf_Out_Gen
!
      implicit none
!
#     include "netcdf.inc"
!
! !DESCRIPTION:
! This routine makes the necessary definitions for the cloud
! netCDF output file.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_VAR_NAME) :: prsp1_dim_name
      integer :: ierr
      integer :: nchr1, nchr2
      integer :: nstr
      integer :: omode
      integer :: pos1
      integer :: varid
      integer :: chrd1(1), chrd(1)
      integer :: hdfd (1)
      integer :: lond (1), latd  (1)
      integer :: prsd (1), prsp1d(1)
      integer :: recd (1), spcd  (1)
      integer :: spcd_2(1), spcd_3(1), spcd_4(1)
      integer :: spcd_6(1), spcd_7(1), spcd_8(1)
      integer :: strd (1)
      integer :: var2(2), var3(3), var4(4), var5(5)

      character(len=200) :: AttrValue
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Define_Netcdf_Out_CloudGT called by ', procID

!     ==========================
      call Define_Netcdf_Out_Gen  &
!     ==========================
     &  (ncid_cloud, nhdf, numLon, numLat, numVert, hdfd, lond, latd, prsd,  &
     &   hdf_dim_name, lon_dim_name, lat_dim_name, prs_dim_name)

!     ------------------
!     Define dimensions.
!     ------------------

      prsp1_dim_name = prs_dim_name
      pos1 = Len_Trim (prs_dim_name) + 1
      prsp1_dim_name(pos1:pos1+2) = 'p1'

      call NcDef_dimension(ncid_cloud, prsp1_dim_name, numVert+1, prsp1d(1))
!                                     --------------
 
       call NcDef_dimension (ncid_cloud,'cloud_dim', num_cloud, spcd(1))
!                                   -----------
       call NcDef_dimension (ncid_cloud, 'chr_dim', 80, chrd(1))
!                                 --------
       call NcDef_dimension (ncid_cloud, rec_dim_name, NF_UNLIMITED, recd(1))
!                                 ------------
!     -----------------------------------------
!     Define variables and variable attributes.
!     -----------------------------------------

      call NcDef_variable (ncid_cloud, 'am', NF_FLOAT, 1, prsd, varid)
!                                     ----
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Hybrid pressure term')
      call NcDef_var_attributes (ncid_cloud, varid, 'units', 'unitless')

      call NcDef_variable (ncid_cloud, 'bm', NF_FLOAT, 1, prsd, varid)
!                                     ----
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Hybrid sigma term')
      call NcDef_var_attributes (ncid_cloud, varid, 'units', 'unitless')

      call NcDef_variable (ncid_cloud, 'ai', NF_FLOAT, 1, prsp1d, varid)
!                                     ----
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Hybrid pressure edge term')
      call NcDef_var_attributes (ncid_cloud, varid, 'units', 'unitless')

      call NcDef_variable (ncid_cloud, 'bi', NF_FLOAT, 1, prsp1d, varid)
!                                     ----
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Hybrid sigma edge term')
      call NcDef_var_attributes (ncid_cloud, varid, 'units', 'unitless')

      call NcDef_variable (ncid_cloud, 'pt', NF_FLOAT, 0, prsp1d, varid)

      call NcDef_variable (ncid_cloud, 'cloud_dim', NF_FLOAT, 1, spcd, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'cloud Species index')
      call NcDef_var_attributes (ncid_cloud, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_cloud, varid, 'coord_labels', 'cloud_labels')
      call NcDef_var_attributes (ncid_cloud, varid, 'selection_category', 'NULL')

      var2(:) = (/ chrd(1), spcd(1) /)
      call NcDef_variable (ncid_cloud, 'cloud_labels', NF_CHAR, 2, var2, varid)
!                                 -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'cloud name')
      call NcDef_var_attributes (ncid_cloud, varid, 'units', 'unitless')
      call NcDef_var_attributes (ncid_cloud, varid, 'selection_category', 'NULL')

      var2(:) = (/ hdfd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, hdr_var_name, NF_INT, 2, var2, varid)
!                                 ------------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Header')
      call NcDef_var_attributes (ncid_cloud, varid, 'units', 'gmi_sec, nymd, nhms')
!
      var2(:) = (/ lond(1), latd(1) /)
      call NcDef_variable (ncid_cloud, 'mcor', NF_FLOAT, 2, var2, varid)
!                                     ------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Grid box area')
      call NcDef_var_attributes (ncid_cloud, varid, 'units', 'm^2')

      if (pr_psf) then
        var3(:) = (/ lond(1), latd(1), recd(1) /)
        call NcDef_variable (ncid_cloud, 'psf', NF_FLOAT, 3, var3, varid)
!                                       -----
        call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Surface pressure')
        call NcDef_var_attributes (ncid_cloud, varid, 'units', 'mb')
      end if

      if (pr_kel) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_cloud, 'kel', NF_FLOAT, 4, var4, varid)
!                                       -----
        call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Temperature')
        call NcDef_var_attributes (ncid_cloud, varid, 'units', 'K')
      end if

      if (pr_mass) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_cloud, 'mass', NF_FLOAT, 4, var4, varid)
!                                       ------
        call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Mass')
        call NcDef_var_attributes (ncid_cloud, varid, 'units', 'kg')
      end if

      if (pr_grid_height) then
        var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
        call NcDef_variable (ncid_cloud, 'gridbox_height', NF_FLOAT, 4, var4, varid)
!                                       ------
        call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Grid Box Height')
        call NcDef_var_attributes (ncid_cloud, varid, 'units', 'm')
      end if


!      ____________For Relative Humidity______________________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud,  'RH', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Relative humidity')
!      _______________For  Cloud Fracton ---------------

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
!
      call NcDef_variable (ncid_cloud, 'CF', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Cloud Fraction')

!      ________________For Liquid Water Content ________________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'LWC', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Liquid Water Content (g/m3)')

!     _________________For  Effective Radius _______________________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'R', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Effective Radius (m)')

!    ____________________For Cloud Optical Depth _____________________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'COD', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Cloud Optical Depth')

!     ____________________For Cloud Optical Depth using SO4 _____________________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'COD_N', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Cloud Optical Depth')

!    _________For COD calculated using Fountoukis & Nenes, 2005 paramet________________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'COD_par', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Cloud Optical Depth')
!     ____________________For Total Aerosol Sulfate  _____________________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'SO4', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Total Aerosol SO4 (ug/m3)')

!     ____________________For Cloud Droplet Number Concentration  _____________________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'CDNC', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Cloud Droplet Number  &
     &                            Concentration (cm-3)')

!     __________For CDNC  calculated using Fountoukis & Nenes, 2005 paramet__________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'CDNC_par', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Cloud Droplet Number  &
     &                            Concentration (cm-3)')

!     ____________________For Effective Radius using SO4 _____________________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'effR', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Effective Radius (m)')

!     ________For Effective Radius  from Fountoukis & Nenes, 2005 paramet _________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'effR_par', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Effective Radius (m)')

!
!     ____________________For Cloud Albedo  _______________________________________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'ALB', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Cloud Albedo')

!     ____________________ Liquid Water Path _____________________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'LWP', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Cloud Albedo')

!     _________ For Cloud Albedo from Fountoukis & Nenes, 2005 paramet ____________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'ALB_par', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Cloud Albedo')


!      _________ For Updraft Velocity  for Fountoukis & Nenes, 2005 paramet ____________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'W_par', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Updraft Velocity (m/s)')

!   _________ For Maximum Supersaturation from Fountoukis & Nenes, 2005 paramet ________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'Smax_par', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Maximum  &
     &                            Supersaturation')
!!!! Start rsot
!      _________ For Autoconversion rate from Khairoutdinov and Kogan, [2000] Eq.29 ____________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'Qaut_KK1', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Qaut K+K, eq with Nd (1/s)')


!      _________ For Autoconversion rate from Khairoutdinov and Kogan, [2000] Eq.30____________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'Qaut_KK2', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Qaut K+K, eq with r_v (1/s)')


!      _________ For Autoconversion rate from Rotstayn, [2000] Eq.15____________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'Qaut_R', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'Qaut, Rotsatyn (1/s)')
!

!      _________ For AIE added on 071206 ____________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'AIE', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'AIE (W/m^2')
!

!      _________ For AIE12 added on 071206 ____________

      var4(:) = (/ lond(1), latd(1), prsd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'AIE12', NF_FLOAT, 4, var4, varid)
!                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'AIE12 (W/m^2')
!
!!!! End rsot

!     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      var5(:) = (/ lond(1), latd(1), prsd(1), spcd(1), recd(1) /)
      call NcDef_variable (ncid_cloud, 'cloud_data', NF_FLOAT, 5, var5, varid)
!c                        -----------
      call NcDef_var_attributes (ncid_cloud, varid, 'long_name', 'cloud_data')
!     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     -------------------------
!     Define global attributes.
!     -------------------------

      call NcDef_glob_attributes (ncid_cloud, 'title',  &
     &       'cloud information diagnostic file')

      call NcSetFill (ncid_cloud, NF_NOFILL, omode)

      call NcEnd_def (ncid_cloud)

      return

      end subroutine Define_Netcdf_Out_CloudGT
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocateVariablesCloudGT
!
! !INTERFACE:
!
      subroutine allocateVariablesCloudGT()
!
      implicit none
!
! !DESCRIPTION:
! Allocates variables necessary to produce cloud outputs.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'allocateVariablesCloudGT called by ', procID

      if (pr_cloud) then

         if (pr_psf) then
            Allocate(psf_nc         (i1:i2, ju1:j2))
            psf_nc = 0.0d0
         end if

         if (pr_kel) then
            Allocate(kel_nc (i1:i2, ju1:j2, kx1:kx2))
            kel_nc = 0.0d0
         end if

         if (pr_mass) then
            Allocate(mass_nc (i1:i2, ju1:j2, kx1:kx2))
            mass_nc = 0.0d0
         end if

         if (pr_grid_height) then
            Allocate(gridbox_height_nc (i1:i2, ju1:j2, kx1:kx2))
            gridbox_height_nc = 0.0d0
         end if

         Allocate (cloud_data        (i1:i2, ju1:j2, kx1:kx2, num_cloud))
         cloud_data = 0.0d0

         Allocate (cloud_data_nc     (i1:i2, ju1:j2, kx1:kx2))
         cloud_data_nc = 0.0d0
 
         Allocate (cloud_rh_nc       (i1:i2, ju1:j2, kx1:kx2))
         cloud_rh_nc   = 0.0d0
 
         Allocate (cloud_cf_nc       (i1:i2, ju1:j2, kx1:kx2))
         cloud_cf_nc   = 0.0d0
 
         Allocate (cloud_lwc_nc      (i1:i2, ju1:j2, kx1:kx2))
         cloud_lwc_nc   = 0.0d0
 
         Allocate (cloud_r_nc        (i1:i2, ju1:j2, kx1:kx2))
         cloud_r_nc   = 0.0d0
 
         Allocate (cloud_cod_nc      (i1:i2, ju1:j2, kx1:kx2))
         cloud_cod_nc   = 0.0d0
 
         Allocate (cloud_cod_n_nc    (i1:i2, ju1:j2, kx1:kx2))
         cloud_cod_n_nc   = 0.0d0
 
         Allocate (cloud_so4_nc      (i1:i2, ju1:j2, kx1:kx2))
         cloud_so4_nc   = 0.0d0
 
         Allocate (cloud_cdnc_nc     (i1:i2, ju1:j2, kx1:kx2))
         cloud_cdnc_nc   = 0.0d0
 
         Allocate (cloud_effr_nc     (i1:i2, ju1:j2, kx1:kx2))
         cloud_effr_nc   = 0.0d0
 
         Allocate (cloud_alb_nc      (i1:i2, ju1:j2, kx1:kx2))
         cloud_alb_nc   = 0.0d0
 
         Allocate (cloud_lwp_nc      (i1:i2, ju1:j2, kx1:kx2))
         cloud_lwp_nc   = 0.0d0
 
         Allocate (cloud_effr_par_nc (i1:i2, ju1:j2, kx1:kx2))
         cloud_effr_par_nc   = 0.0d0
 
         Allocate (cloud_alb_par_nc  (i1:i2, ju1:j2, kx1:kx2))
         cloud_alb_par_nc   = 0.0d0
 
         Allocate (cloud_cod_par_nc  (i1:i2, ju1:j2, kx1:kx2))
         cloud_cod_par_nc   = 0.0d0
 
         Allocate (cloud_cdnc_par_nc (i1:i2, ju1:j2, kx1:kx2))
         cloud_cdnc_par_nc   = 0.0d0

         Allocate (cloud_w_par_nc    (i1:i2, ju1:j2, kx1:kx2))
         cloud_w_par_nc   = 0.0d0

         Allocate (cloud_smax_par_nc (i1:i2, ju1:j2, kx1:kx2))
         cloud_smax_par_nc   = 0.0d0

         Allocate (Qaut_NS_KK1       (i1:i2, ju1:j2, kx1:kx2))
         Qaut_NS_KK1 = 0.0d0

         Allocate (Qaut_NS_KK2       (i1:i2, ju1:j2, kx1:kx2))
         Qaut_NS_KK2 = 0.0d0

         Allocate (Qaut_NS_R         (i1:i2, ju1:j2, kx1:kx2))
         Qaut_NS_R = 0.0d0

         Allocate (cloud_AIE_nc      (i1:i2, ju1:j2, kx1:kx2))
         cloud_AIE_nc = 0.0d0

         Allocate (cloud_AIE12_nc    (i1:i2, ju1:j2, kx1:kx2))
         cloud_AIE12_nc = 0.0d0
      endif

     return

     end  subroutine allocateVariablesCloudGT

!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_Hdr_CloudGT
!
! !INTERFACE:
!
      subroutine Write_Netcdf_Hdr_CloudGT  (gmiDomain, metFields)
!
! !USES:
      use m_ncGeneralOpsOutput, only : WriteNetcdfHdrGen
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_metFields), intent(in) :: metFields
      type(t_gmiDomain), intent(in) :: gmiDomain
!
! !DESCRIPTION:
! This routine creates some header information for the cloud netCDF output
! file and writes it out.
!
! !LOCAL VARIABLES:
      logical, save :: first = .true.
      character(len=80)  :: spclab(num_cloud)
      real*8             :: aerdat(num_cloud)
      integer :: ic, icx
      integer :: count1d (1), count2d (2)
      integer :: start1d(1), start2d(2)
      real*8  :: prsdat(k1:k2)
      real*8 , allocatable :: locLatDeg(:)
      real*8 , allocatable :: locLonDeg(:)
      real*8 , allocatable :: locMCOR(:,:)
      real*8 , allocatable :: mcor(:,:)
      real*8 , allocatable :: londeg(:), latdeg(:)
      real*8 , allocatable :: ai(:), bi(:)
      real*8 , allocatable :: am(:), bm(:)
      real*8 :: pt
!EOP
!------------------------------------------------------------------------------
!BOC

      if (pr_diag) Write (6,*) 'Write_Netcdf_Hdr_CloudGT called by ', procID

      allocate(latdeg(ju1_gl:j2_gl))
      call Get_latdeg(gmiDomain, latdeg)

      allocate(londeg(i1_gl:i2_gl))
      call Get_londeg(gmiDomain, londeg)

      allocate(mcor(i1_gl:i2_gl,ju1_gl:j2_gl))
      call Get_mcorGlob(gmiDomain, mcor)

      allocate(locLonDeg(1:numLon))
      allocate(locLatDeg(1:numLat))

      locLonDeg(1:numLon) = londeg(ix1:ix2)
      locLatDeg(1:numLat) = latdeg(jx1:jx2)

!     =========================
      call WriteNetcdfHdrGen  &
!     =========================
     &  (ncid_cloud, locLatDeg, locLonDeg, pr_diag, procID, numLon, numLat, &
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

      call Ncwr_1d (prsdat(kx1:kx2), ncid_cloud, prs_dim_name, start1d, count1d)
      call Ncwr_1d (am(kx1:kx2), ncid_cloud, 'am', start1d, count1d)
      call Ncwr_1d (bm(kx1:kx2), ncid_cloud, 'bm', start1d, count1d)

      count1d(1) = numVert + 1

      call Ncwr_1d (ai(kx1-1:kx2), ncid_cloud, 'ai', start1d, count1d)
      call Ncwr_1d (bi(kx1-1:kx2), ncid_cloud, 'bi', start1d, count1d)

      call Ncwr_Scal (pt, ncid_cloud, 'pt')


!     -------------------
!     cloud labels.
!     -------------------

      do ic = 1, num_cloud
        aerdat(ic) = ic
      end do

      start1d(1) = 1
      count1d (1) = num_cloud

      call Ncwr_1d (aerdat, ncid_cloud, 'cloud_dim', start1d, count1d)

      do ic = 1, num_cloud
        spclab(ic) = cloud_labels(ic)
      end do

      start2d(:) = (/  1, 1 /)
      count2d (:) = (/ 80, num_cloud /)

      call Ncwr_2d_Char (spclab, ncid_cloud, 'cloud_labels', start2d, count2d)

!     --------------
!     Grid box area.
!     --------------

      start2d(:) = (/ 1, 1 /)
      count2d (:) = (/ numLon, numLat /)

      allocate(locMCOR(ix1:ix2, jx1:jx2))

      locMCOR(ix1:ix2,jx1:jx2) = mcor(ix1:ix2, jx1:jx2)
      call Ncwr_2d (locMCOR, ncid_cloud, 'mcor', start2d, count2d)

      deallocate(locMCOR)
      deallocate(locLonDeg)
      deallocate(locLatDeg)

      return

      end subroutine Write_Netcdf_Hdr_CloudGT
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Prep_Netcdf_CloudGT
!
! !INTERFACE:
!
      subroutine Prep_Netcdf_CloudGT (time_for_cloud, Chemistry, metFields)
!
! !USES:
      use GmiChemistryMethod_mod , only : Get_cloud_param

      implicit none
!
! !INPUT PARAMETERS:
      logical          , intent(in) :: time_for_cloud ! time for netCDF output?
      type(t_metFields), intent(in) :: metFields
      type(t_Chemistry), intent(in) :: Chemistry
!
! !DESCRIPTION:
! This routine prepares the const netCDF output.
!
! !LOCAL VARIABLES:
      integer :: ic, icx, il, ik, ij
      real*8, save :: counter = 0
      real*8, allocatable :: pctm2(:,:)
      real*8, allocatable :: kel(:,:,:)
      real*8, allocatable :: mass(:,:,:)
      real*8, allocatable :: gridBoxHeight(:,:,:)
      real*8, allocatable :: cloud_param(:,:,:,:)
      real*8, allocatable, save :: countern(:,:,:)
      real*8, allocatable, save :: rcounter(:,:,:,:)
      integer, allocatable :: lwi_flags(:,:)
      logical, save :: first = .true.
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Prep_Netcdf_CloudGT called by ', procID

!     ==========
      if (first) then
!     ==========

      first = .false.

      Allocate (countern (i1:i2, ju1:j2, k1:k2))
      countern = 0.0d0
      Allocate (rcounter (i1:i2, ju1:j2, k1:k2,1))
      rcounter = 0.0d0

      end if

      allocate(cloud_param(i1:i2, ju1:j2, k1:k2, numSpecies))
      call Get_cloud_param (Chemistry, cloud_param)

      allocate(lwi_flags(i1:i2,ju1:j2))
      call Get_lwi_flags(metFields, lwi_flags)

      if (pr_psf) then                        
         allocate(pctm2(ilo:ihi,julo:jhi))
         call Get_pctm2(metFields, pctm2)       
      end if      
     
      if (pr_kel) then                        
         allocate(kel(ilo:ihi,julo:jhi,k1:k2))
         call Get_kel(metFields, kel)       
      end if      
     
      if (pr_mass) then                        
         allocate(mass(i1:i2,ju1:j2,k1:k2))
         call Get_mass(metFields, mass)       
      end if      
     
      if (pr_grid_height) then     
         allocate(gridBoxHeight(i1:i2,ju1:j2,k1:k2))
         call Get_gridBoxHeight(metFields, gridBoxHeight)  
      end if

      if (do_mean) then

         if (counter == 0.d0) then

!         ----------------
!            Reset means
!         ----------------

            cloud_data     (:,:,:,:) = 0.0d0

            cloud_rh_nc      (:,:,:) = 0.d0
            cloud_cf_nc      (:,:,:) = 0.d0
            cloud_lwc_nc     (:,:,:) = 0.d0
            cloud_r_nc       (:,:,:) = 0.d0
            cloud_so4_nc     (:,:,:) = 0.d0
            cloud_alb_nc     (:,:,:) = 0.d0
            cloud_lwp_nc     (:,:,:) = 0.d0
            cloud_alb_par_nc (:,:,:) = 0.d0
            cloud_cod_nc     (:,:,:) = 0.d0
            cloud_cod_n_nc   (:,:,:) = 0.d0
            cloud_cdnc_nc    (:,:,:) = 0.d0
            cloud_effr_nc    (:,:,:) = 0.d0
            cloud_cdnc_par_nc(:,:,:) = 0.d0
            cloud_w_par_nc   (:,:,:) = 0.d0
            cloud_smax_par_nc(:,:,:) = 0.d0
            cloud_effr_par_nc(:,:,:) = 0.d0
            cloud_cod_par_nc (:,:,:) = 0.d0
            Qaut_NS_KK1      (:,:,:) = 0.d0
            Qaut_NS_KK2      (:,:,:) = 0.d0
            Qaut_NS_R        (:,:,:) = 0.d0
            cloud_AIE_nc     (:,:,:) = 0.d0
            cloud_AIE12_nc   (:,:,:) = 0.d0

            if (pr_psf) psf_nc = 0.0d0
            if (pr_kel) kel_nc = 0.0d0
            if (pr_mass) mass_nc = 0.0d0
            if (pr_grid_height) gridbox_height_nc = 0.0d0
         end if

!       -----------------
!        Calculate means
!       -----------------
     
         if (pr_psf) psf_nc(:,:) =  psf_nc(:,:) + pctm2(i1:i2,ju1:j2)
         if (pr_kel) kel_nc(:,:,:) =  kel_nc(:,:,:) + kel(i1:i2,ju1:j2,kx1:kx2)

         if (pr_mass) then
             mass_nc(:,:,:) =  mass_nc(:,:,:) + mass(:,:,kx1:kx2)
         end if

         if (pr_grid_height) then
            gridbox_height_nc(:,:,:) = gridbox_height_nc(:,:,:) + &
     &                           gridBoxHeight(:,:,kx1:kx2)
         end if

         cloud_rh_nc      (:,:,:)= cloud_rh_nc      (:,:,:)+ cloud_param(:,:,:,1)
         cloud_cf_nc      (:,:,:)= cloud_cf_nc      (:,:,:)+ cloud_param(:,:,:,2)
         cloud_lwc_nc     (:,:,:)= cloud_lwc_nc     (:,:,:)+ cloud_param(:,:,:,3)
         cloud_r_nc       (:,:,:)= cloud_r_nc       (:,:,:)+ cloud_param(:,:,:,4)
         cloud_so4_nc     (:,:,:)= cloud_so4_nc     (:,:,:)+ cloud_param(:,:,:,7)
!!       cloud_alb_nc     (:,:,:)= cloud_alb_nc   (:,:,:)+ cloud_param(:,:,:,10)
         cloud_lwp_nc     (:,:,:)= cloud_lwp_nc     (:,:,:)+ cloud_param(:,:,:,12)
         cloud_alb_par_nc (:,:,:)= cloud_alb_par_nc (:,:,:)+ cloud_param(:,:,:,18)
         cloud_AIE_nc     (:,:,:)= cloud_AIE_nc     (:,:,:)+ cloud_param(:,:,:,22)
         cloud_AIE12_nc   (:,:,:)= cloud_AIE12_nc   (:,:,:)+ cloud_param(:,:,:,23)
  
         cloud_data(:,:,:,1) = cloud_data(:,:,:,1)  + cloud_param(:,:,:,1)
         cloud_data(:,:,:,2) = cloud_data(:,:,:,2)  + cloud_param(:,:,:,2)
         cloud_data(:,:,:,3) = cloud_data(:,:,:,3)  + cloud_param(:,:,:,3)
         cloud_data(:,:,:,4) = cloud_data(:,:,:,4)  + cloud_param(:,:,:,4)
         cloud_data(:,:,:,6) = cloud_data(:,:,:,6)  + cloud_param(:,:,:,7)
         cloud_data(:,:,:,9) = cloud_data(:,:,:,9)  + cloud_param(:,:,:,10)
         cloud_data(:,:,:,11)= cloud_data(:,:,:,11) + cloud_param(:,:,:,12)
         cloud_data(:,:,:,17)= cloud_data(:,:,:,17) + cloud_param(:,:,:,18)
         cloud_data(:,:,:,21)= cloud_data(:,:,:,21) + cloud_param(:,:,:,22)
         cloud_data(:,:,:,22)= cloud_data(:,:,:,22) + cloud_param(:,:,:,23)


         counter = counter + 1.d0

         do ik = kx1, kx2
            do ij = ju1, j2
               do il = i1, i2

                  !!! For ocean:
                  if ((lwi_flags(il,ij)==1)           .and.  &
                      (kel(il,ij,ik)>=269.15d0)       .and. &
     &                (cloud_param(il,ij,ik,3)>0.d0)  .and. &
                      (cloud_param(il,ij,ik,13)>0.d0)) then

                     countern(il,ij,ik) = countern(il,ij,ik) + 1.d0
                     cloud_cod_nc     (il,ij,ik)= cloud_cod_nc     (il,ij,ik) + &
     &                                            cloud_param(il,ij,ik,5)
                     cloud_cdnc_nc    (il,ij,ik)= cloud_cdnc_nc    (il,ij,ik) + &
     &                                            cloud_param(il,ij,ik,8)
                     cloud_effr_nc    (il,ij,ik)= cloud_effr_nc    (il,ij,ik) + &
     &                                            cloud_param(il,ij,ik,9)
                     cloud_cod_n_nc   (il,ij,ik)= cloud_cod_n_nc   (il,ij,ik) + &
     &                                            cloud_param(il,ij,ik,11)
                     cloud_cdnc_par_nc(il,ij,ik)= cloud_cdnc_par_nc(il,ij,ik) + &
     &                                            cloud_param(il,ij,ik,13)
                     cloud_w_par_nc   (il,ij,ik)= cloud_w_par_nc   (il,ij,ik) + &
     &                                            cloud_param(il,ij,ik,14)
                     cloud_smax_par_nc(il,ij,ik)= cloud_smax_par_nc(il,ij,ik) + &
     &                                            cloud_param(il,ij,ik,15)
                     cloud_effr_par_nc(il,ij,ik)= cloud_effr_par_nc(il,ij,ik) + &
     &                                            cloud_param(il,ij,ik,16)
                     cloud_cod_par_nc (il,ij,ik)= cloud_cod_par_nc (il,ij,ik) + &
     &                                            cloud_param(il,ij,ik,17)
                     Qaut_NS_KK1      (il,ij,ik)= Qaut_NS_KK1      (il,ij,ik) + &
     &                                            cloud_param(il,ij,ik,19)
                     Qaut_NS_KK2      (il,ij,ik)= Qaut_NS_KK2      (il,ij,ik) + &
     &                                            cloud_param(il,ij,ik,20)
                     Qaut_NS_R        (il,ij,ik)= Qaut_NS_R        (il,ij,ik) + &
     &                                            cloud_param(il,ij,ik,21)

                     cloud_data(il,ij,ik,5) = cloud_data(il,ij,ik,5)  + &
     &                                        cloud_param(il,ij,ik,5)
                     cloud_data(il,ij,ik,7) = cloud_data(il,ij,ik,7)  + &
     &                                        cloud_param(il,ij,ik,8)
                     cloud_data(il,ij,ik,8) = cloud_data(il,ij,ik,8)  + &
     &                                        cloud_param(il,ij,ik,9)
                     cloud_data(il,ij,ik,10)= cloud_data(il,ij,ik,10) + &
     &                                        cloud_param(il,ij,ik,11)
                     cloud_data(il,ij,ik,12)= cloud_data(il,ij,ik,12) + &
     &                                        cloud_param(il,ij,ik,13)
                     cloud_data(il,ij,ik,13)= cloud_data(il,ij,ik,13) + &
     &                                        cloud_param(il,ij,ik,14)
                     cloud_data(il,ij,ik,14)= cloud_data(il,ij,ik,14) + &
     &                                        cloud_param(il,ij,ik,15)
                     cloud_data(il,ij,ik,15)= cloud_data(il,ij,ik,15) + &
     &                                        cloud_param(il,ij,ik,16)
                     cloud_data(il,ij,ik,16)= cloud_data(il,ij,ik,16) + &
     &                                        cloud_param(il,ij,ik,17)
                     cloud_data(il,ij,ik,18)= cloud_data(il,ij,ik,18) + &
     &                                        cloud_param(il,ij,ik,19)
                     cloud_data(il,ij,ik,19)= cloud_data(il,ij,ik,19) + &
     &                                        cloud_param(il,ij,ik,20)
                     cloud_data(il,ij,ik,20)= cloud_data(il,ij,ik,20) + &
     &                                        cloud_param(il,ij,ik,21)

                  !!!    For land:
                  else if ((lwi_flags(il,ij)==2)           .and. &
     &                     (kel(il,ij,ik)>=263.15d0)       .and. &
     &                     (cloud_param(il,ij,ik,3)>0d0)   .and. &
     &                     (cloud_param(il,ij,ik,13)>0.d0)) then

                     countern(il,ij,ik) = countern(il,ij,ik) + 1.d0
                     cloud_cod_nc     (il,ij,ik)= cloud_cod_nc     (il,ij,ik)+ &
     &                                            cloud_param(il,ij,ik,5)
                     cloud_cdnc_nc    (il,ij,ik)= cloud_cdnc_nc    (il,ij,ik)+ &
     &                                            cloud_param(il,ij,ik,8)
                     cloud_effr_nc    (il,ij,ik)= cloud_effr_nc    (il,ij,ik)+ &
     &                                            cloud_param(il,ij,ik,9)
                     cloud_cod_n_nc   (il,ij,ik)= cloud_cod_n_nc   (il,ij,ik)+ &
     &                                            cloud_param(il,ij,ik,11)
                     cloud_cdnc_par_nc(il,ij,ik)= cloud_cdnc_par_nc(il,ij,ik)+ &
     &                                            cloud_param(il,ij,ik,13)
                     cloud_w_par_nc   (il,ij,ik)= cloud_w_par_nc   (il,ij,ik)+ &
     &                                            cloud_param(il,ij,ik,14)
                     cloud_smax_par_nc(il,ij,ik)= cloud_smax_par_nc(il,ij,ik)+ &
     &                                            cloud_param(il,ij,ik,15)
                     cloud_effr_par_nc(il,ij,ik)= cloud_effr_par_nc(il,ij,ik)+ &
     &                                            cloud_param(il,ij,ik,16)
                     cloud_cod_par_nc (il,ij,ik)= cloud_cod_par_nc (il,ij,ik)+ &
     &                                            cloud_param(il,ij,ik,17)
                     Qaut_NS_KK1      (il,ij,ik)= Qaut_NS_KK1      (il,ij,ik)+ &
     &                                            cloud_param(il,ij,ik,19)
                     Qaut_NS_KK2      (il,ij,ik)= Qaut_NS_KK2      (il,ij,ik)+ &
     &                                            cloud_param(il,ij,ik,20)
                     Qaut_NS_R        (il,ij,ik)= Qaut_NS_R        (il,ij,ik)+ &
     &                                            cloud_param(il,ij,ik,21)
              

                     cloud_data(il,ij,ik,5) = cloud_data(il,ij,ik,5)  + &
     &                                            cloud_param(il,ij,ik,5)
                     cloud_data(il,ij,ik,7) = cloud_data(il,ij,ik,7)  + &
     &                                            cloud_param(il,ij,ik,8)
                     cloud_data(il,ij,ik,8) = cloud_data(il,ij,ik,8)  + &
     &                                            cloud_param(il,ij,ik,9)
                     cloud_data(il,ij,ik,10)= cloud_data(il,ij,ik,10) + &
     &                                            cloud_param(il,ij,ik,11)
                     cloud_data(il,ij,ik,12)= cloud_data(il,ij,ik,12) + &
     &                                            cloud_param(il,ij,ik,13)
                     cloud_data(il,ij,ik,13)= cloud_data(il,ij,ik,13) + &
     &                                            cloud_param(il,ij,ik,14)
                     cloud_data(il,ij,ik,14)= cloud_data(il,ij,ik,14) + &
     &                                            cloud_param(il,ij,ik,15)
                     cloud_data(il,ij,ik,15)= cloud_data(il,ij,ik,15) + &
     &                                            cloud_param(il,ij,ik,16)
                     cloud_data(il,ij,ik,16)= cloud_data(il,ij,ik,16) + &
     &                                            cloud_param(il,ij,ik,17)
                     cloud_data(il,ij,ik,18)= cloud_data(il,ij,ik,18) + &
     &                                            cloud_param(il,ij,ik,19)
                     cloud_data(il,ij,ik,19)= cloud_data(il,ij,ik,19) + &
     &                                            cloud_param(il,ij,ik,20)
                     cloud_data(il,ij,ik,20)= cloud_data(il,ij,ik,20) + &
     &                                            cloud_param(il,ij,ik,21)
              
                  end if
               end do
            end do
         end do

         cloud_alb_nc     (:,:,:)= countern  (:,:,:)


         if (time_for_cloud) then
!           ----------------
!           Calculate means.
!           ----------------
            if (pr_psf) psf_nc(:,:) =  psf_nc(:,:) / counter
            if (pr_kel) kel_nc(:,:,:) =  kel_nc(:,:,:) / counter
            if (pr_mass) mass_nc(:,:,:) =  mass_nc(:,:,:) / counter
            if (pr_grid_height) gridbox_height_nc(:,:,:) = gridbox_height_nc(:,:,:) / &
     &                           counter

            cloud_rh_nc      (:,:,:) = cloud_rh_nc      (:,:,:)/counter
            cloud_cf_nc      (:,:,:) = cloud_cf_nc      (:,:,:)/counter
            cloud_lwc_nc     (:,:,:) = cloud_lwc_nc     (:,:,:)/counter
            cloud_r_nc       (:,:,:) = cloud_r_nc       (:,:,:)/counter
            cloud_so4_nc     (:,:,:) = cloud_so4_nc     (:,:,:)/counter
            cloud_lwp_nc     (:,:,:) = cloud_lwp_nc     (:,:,:)/counter
            cloud_alb_par_nc (:,:,:) = cloud_alb_par_nc (:,:,:)/counter
            cloud_AIE_nc     (:,:,:) = cloud_AIE_nc     (:,:,:)/counter  !added on 071206
            cloud_AIE12_nc   (:,:,:) = cloud_AIE12_nc   (:,:,:)/counter  !added on 071206
            cloud_alb_nc     (:,:,:) = cloud_alb_nc     (:,:,:)/1.d0
  
            where (countern(:,:,:) .ne. 0.d0)

               cloud_cod_nc     (:,:,:) = cloud_cod_nc     (:,:,:)/countern(:,:,:)
               cloud_cod_n_nc   (:,:,:) = cloud_cod_n_nc   (:,:,:)/countern(:,:,:)
               cloud_cdnc_nc    (:,:,:) = cloud_cdnc_nc    (:,:,:)/countern(:,:,:)
               cloud_effr_nc    (:,:,:) = cloud_effr_nc    (:,:,:)/countern(:,:,:)
               cloud_cdnc_par_nc(:,:,:) = cloud_cdnc_par_nc(:,:,:)/countern(:,:,:)
               cloud_w_par_nc   (:,:,:) = cloud_w_par_nc   (:,:,:)/countern(:,:,:)
               cloud_smax_par_nc(:,:,:) = cloud_smax_par_nc(:,:,:)/countern(:,:,:)
               cloud_effr_par_nc(:,:,:) = cloud_effr_par_nc(:,:,:)/countern(:,:,:)
               cloud_cod_par_nc (:,:,:) = cloud_cod_par_nc (:,:,:)/countern(:,:,:)
               Qaut_NS_KK1      (:,:,:) = Qaut_NS_KK1      (:,:,:)/countern(:,:,:)
               Qaut_NS_KK2      (:,:,:) = Qaut_NS_KK2      (:,:,:)/countern(:,:,:)
               Qaut_NS_R        (:,:,:) = Qaut_NS_R        (:,:,:)/countern(:,:,:)

            elsewhere

               cloud_cod_nc     (:,:,:) = 0.d0
               cloud_cod_n_nc   (:,:,:) = 0.d0
               cloud_cdnc_nc    (:,:,:) = 0.d0
               cloud_effr_nc    (:,:,:) = 0.d0
               cloud_cdnc_par_nc(:,:,:) = 0.d0
               cloud_w_par_nc   (:,:,:) = 0.d0
               cloud_smax_par_nc(:,:,:) = 0.d0
               cloud_effr_par_nc(:,:,:) = 0.d0
               cloud_cod_par_nc (:,:,:) = 0.d0
               Qaut_NS_KK1      (:,:,:) = 0.d0
               Qaut_NS_KK2      (:,:,:) = 0.d0
               Qaut_NS_R        (:,:,:) = 0.d0

            end where

            cloud_data(:,:,:,1) = cloud_data(:,:,:,1) /counter
            cloud_data(:,:,:,2) = cloud_data(:,:,:,2) /counter
            cloud_data(:,:,:,3) = cloud_data(:,:,:,3) /counter
            cloud_data(:,:,:,4) = cloud_data(:,:,:,4) /counter
            cloud_data(:,:,:,6) = cloud_data(:,:,:,6) /counter
            cloud_data(:,:,:,9) = cloud_data(:,:,:,9) /counter
            cloud_data(:,:,:,11)= cloud_data(:,:,:,11)/counter
            cloud_data(:,:,:,17)= cloud_data(:,:,:,17)/counter
            cloud_data(:,:,:,21)= cloud_data(:,:,:,21)/counter
            cloud_data(:,:,:,22)= cloud_data(:,:,:,22)/counter

            rcounter (:,:,:,1) = countern(:,:,:)

            do ik = kx1, kx2
               do ij = ju1, j2
                  do il = i1, i2
                     if (rcounter(il,ij,ik,1) .ne. 0.d0) then
                        cloud_data(il,ij,ik,5) = cloud_data(il,ij,ik,5) / &
     &                                           rcounter(il,ij,ik,1)
                        cloud_data(il,ij,ik,7) = cloud_data(il,ij,ik,7) / &
     &                                           rcounter(il,ij,ik,1)
                        cloud_data(il,ij,ik,8) = cloud_data(il,ij,ik,8) / &
     &                                           rcounter(il,ij,ik,1)
                        cloud_data(il,ij,ik,10)= cloud_data(il,ij,ik,10)/ &
     &                                           rcounter(il,ij,ik,1)
                        cloud_data(il,ij,ik,12)= cloud_data(il,ij,ik,12)/ &
     &                                           rcounter(il,ij,ik,1)
                        cloud_data(il,ij,ik,13)= cloud_data(il,ij,ik,13)/ &
     &                                           rcounter(il,ij,ik,1)
                        cloud_data(il,ij,ik,14)= cloud_data(il,ij,ik,14)/ &
     &                                           rcounter(il,ij,ik,1)
                        cloud_data(il,ij,ik,15)= cloud_data(il,ij,ik,15)/ &
     &                                           rcounter(il,ij,ik,1)
                        cloud_data(il,ij,ik,16)= cloud_data(il,ij,ik,16)/ &
     &                                           rcounter(il,ij,ik,1)
                        cloud_data(il,ij,ik,18)= cloud_data(il,ij,ik,18)/ &
     &                                           rcounter(il,ij,ik,1)
                        cloud_data(il,ij,ik,19)= cloud_data(il,ij,ik,19)/ &
     &                                           rcounter(il,ij,ik,1)
                        cloud_data(il,ij,ik,20)= cloud_data(il,ij,ik,20)/ &
     &                                           rcounter(il,ij,ik,1)
                     end if
                  end do
               end do
            end do

            countern(:,:,:) = 0.d0
            rcounter(:,:,:,1) = 0.d0
            counter = 0.d0
         end if

!     ====
      else
!     ====

         if (time_for_cloud) then

!           --------------------------------------------------
!           Fill output arrays with current "non-mean" values.
!           --------------------------------------------------
!
           if (pr_cloud) then
              if (pr_psf) psf_nc(:,:) =  psf_nc(:,:) + pctm2(i1:i2,ju1:j2)
              if (pr_kel) kel_nc(:,:,:) =  kel_nc(:,:,:) + kel(i1:i2,ju1:j2,kx1:kx2)
              if (pr_mass) mass_nc(:,:,:) =  mass_nc(:,:,:) + mass(:,:,kx1:kx2)
              if (pr_grid_height) gridbox_height_nc(:,:,:) = gridbox_height_nc(:,:,:) + &
     &                           gridBoxHeight(:,:,kx1:kx2)

              cloud_rh_nc    (:,:,:)  = cloud_rh_nc    (:,:,:)  + cloud_param(:,:,:,1)
              cloud_cf_nc    (:,:,:)  = cloud_cf_nc    (:,:,:)  + cloud_param(:,:,:,2)
              cloud_lwc_nc   (:,:,:)  = cloud_lwc_nc   (:,:,:)  + cloud_param(:,:,:,3)
              cloud_r_nc     (:,:,:)  = cloud_r_nc     (:,:,:)  + cloud_param(:,:,:,4)
              cloud_cod_nc   (:,:,:)  = cloud_cod_nc   (:,:,:)  + cloud_param(:,:,:,5)
              cloud_so4_nc   (:,:,:)  = cloud_so4_nc   (:,:,:)  + cloud_param(:,:,:,7)
              cloud_cdnc_nc  (:,:,:)  = cloud_cdnc_nc  (:,:,:)  + cloud_param(:,:,:,8)
              cloud_effr_nc  (:,:,:)  = cloud_effr_nc  (:,:,:)  + cloud_param(:,:,:,9)
              cloud_cod_n_nc (:,:,:)  = cloud_cod_n_nc (:,:,:)  + cloud_param(:,:,:,11)
              cloud_alb_nc   (:,:,:)  = cloud_alb_nc   (:,:,:)  + cloud_param(:,:,:,10)
              cloud_lwp_nc   (:,:,:)  = cloud_lwp_nc   (:,:,:)  + cloud_param(:,:,:,12)

              cloud_cdnc_par_nc  (:,:,:)  = cloud_cdnc_par_nc  (:,:,:)  +  &
     &                                               cloud_param(:,:,:,13)
              cloud_w_par_nc     (:,:,:)  = cloud_w_par_nc     (:,:,:)  +  &
     &                                               cloud_param(:,:,:,14)
              cloud_smax_par_nc  (:,:,:)  = cloud_smax_par_nc  (:,:,:)  +  &
     &                                               cloud_param(:,:,:,15)
              cloud_effr_par_nc  (:,:,:)  = cloud_effr_par_nc  (:,:,:)  +  &
     &                                               cloud_param(:,:,:,16)
              cloud_cod_par_nc   (:,:,:)  = cloud_cod_par_nc   (:,:,:)  +  &
     &                                               cloud_param(:,:,:,17)
              cloud_alb_par_nc   (:,:,:)  = cloud_alb_par_nc   (:,:,:)  +  &
     &                                               cloud_param(:,:,:,18)
              Qaut_NS_KK1  (:,:,:)  = Qaut_NS_KK1  (:,:,:)  + &
     &                                               cloud_param(:,:,:,19)
              Qaut_NS_KK2  (:,:,:)  = Qaut_NS_KK2  (:,:,:)  + &
     &                                               cloud_param(:,:,:,20)
              Qaut_NS_R    (:,:,:)  = Qaut_NS_R    (:,:,:)  + &
     &                                               cloud_param(:,:,:,21)
              cloud_AIE_nc (:,:,:)  = cloud_AIE_nc    (:,:,:)  + &
     &                                               cloud_param(:,:,:,22)  ! added on 071206
              cloud_AIE12_nc  (:,:,:)  = cloud_AIE12_nc  (:,:,:)  + &
     &                                               cloud_param(:,:,:,23)  ! added on 071206
!
              cloud_data(:,:,:,1) = cloud_data(:,:,:,1)  + cloud_param(:,:,:,1)
              cloud_data(:,:,:,2) = cloud_data(:,:,:,2)  + cloud_param(:,:,:,2)
              cloud_data(:,:,:,3) = cloud_data(:,:,:,3)  + cloud_param(:,:,:,3)
              cloud_data(:,:,:,4) = cloud_data(:,:,:,4)  + cloud_param(:,:,:,4)
              cloud_data(:,:,:,5) = cloud_data(:,:,:,5)  + cloud_param(:,:,:,5)
              cloud_data(:,:,:,6) = cloud_data(:,:,:,6)  + cloud_param(:,:,:,7)
              cloud_data(:,:,:,7) = cloud_data(:,:,:,7)  + cloud_param(:,:,:,8)
              cloud_data(:,:,:,8) = cloud_data(:,:,:,8)  + cloud_param(:,:,:,9)
              cloud_data(:,:,:,9) = cloud_data(:,:,:,9)  + cloud_param(:,:,:,10)
              cloud_data(:,:,:,10)= cloud_data(:,:,:,10) + cloud_param(:,:,:,11)
              cloud_data(:,:,:,11)= cloud_data(:,:,:,11) + cloud_param(:,:,:,12)
              cloud_data(:,:,:,12)= cloud_data(:,:,:,12) + cloud_param(:,:,:,13)
              cloud_data(:,:,:,13)= cloud_data(:,:,:,13) + cloud_param(:,:,:,14)
              cloud_data(:,:,:,14)= cloud_data(:,:,:,14) + cloud_param(:,:,:,15)
              cloud_data(:,:,:,15)= cloud_data(:,:,:,15) + cloud_param(:,:,:,16)
              cloud_data(:,:,:,16)= cloud_data(:,:,:,16) + cloud_param(:,:,:,17)
              cloud_data(:,:,:,17)= cloud_data(:,:,:,17) + cloud_param(:,:,:,18)
              cloud_data(:,:,:,18)= cloud_data(:,:,:,18) + cloud_param(:,:,:,19)
              cloud_data(:,:,:,19)= cloud_data(:,:,:,19) + cloud_param(:,:,:,20)
              cloud_data(:,:,:,20)= cloud_data(:,:,:,20) + cloud_param(:,:,:,21)
              cloud_data(:,:,:,21)= cloud_data(:,:,:,21) + cloud_param(:,:,:,22)
              cloud_data(:,:,:,22)= cloud_data(:,:,:,22) + cloud_param(:,:,:,23)
           end if
        end if
      end if  ! end of if (do_mean)

      deallocate(cloud_param)

      if (pr_psf) deallocate(pctm2)

      if (pr_kel) deallocate(kel)

      if (pr_mass) deallocate(mass)

      if (pr_grid_height) deallocate(gridBoxHeight)

      return

      end subroutine Prep_Netcdf_CloudGT
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Netcdf_CloudGT
!
! !INTERFACE:
!
      subroutine Write_Netcdf_CloudGT (nymd, nhms, gmi_sec)

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
      integer, parameter :: SG_cloud_cod_par_nc  = 7001
      integer, parameter :: SG_cloud_cdnc_par_nc = 7002
      integer, parameter :: SG_cloud_alb_par_nc  = 7003
      integer, parameter :: SG_cloud_smax_par_nc = 7004
      integer, parameter :: SG_cloud_effr_par_nc = 7005
      integer, parameter :: SG_cloud_w_par_nc    = 7006
      integer, parameter :: SG_Cloud_NC          = 7007
      integer, parameter :: SG_cloud_rh_nc       = 7008
      integer, parameter :: SG_cloud_cf_nc       = 7009
      integer, parameter :: SG_cloud_lwc_nc      = 7010
      integer, parameter :: SG_cloud_r_nc        = 7011
      integer, parameter :: SG_cloud_cod_nc      = 7012
      integer, parameter :: SG_cloud_so4_nc      = 7013
      integer, parameter :: SG_cloud_cdnc_nc     = 7014
      integer, parameter :: SG_cloud_effr_nc     = 7015
      integer, parameter :: SG_cloud_cod_n_nc    = 7016
      integer, parameter :: SG_cloud_alb_nc      = 7017
      integer, parameter :: SG_cloud_lwp_nc      = 7018
      integer, parameter :: SG_Qaut_NS_KK1  = 7100
      integer, parameter :: SG_Qaut_NS_KK2  = 7101
      integer, parameter :: SG_Qaut_NS_R    = 7102
      integer, parameter :: SG_AIE_nc       = 7103
      integer, parameter :: SG_AIE12_nc     = 7104
      integer, parameter :: SG_p3           = 7105
      integer, parameter :: SG_p4           = 7106
      integer, parameter :: SG_p5           = 7107
      integer, parameter :: SG_p6           = 7108
      integer, parameter :: SG_p7           = 7109
      integer, parameter :: SG_p8           = 7110
      integer, parameter :: SG_p9           = 7111
      integer, parameter :: SG_p10          = 7112
      integer, parameter :: SG_p11          = 7113
      integer, parameter :: SG_p12          = 7114
      integer, parameter :: SG_p13          = 7115
      integer, parameter :: SG_p14          = 7116
      integer, parameter :: SG_p15          = 7117
      integer, parameter :: SG_p16          = 7118
      integer, parameter :: SG_p17          = 7119
!
! !LOCAL VARIABLES:
      integer :: ic
      integer :: count2d (2), count3d (3), count4d (4)
      integer :: start2d(2), start3d(3), start4d(4)
      integer :: hdr(NETCDF_HDF)
      real*8, allocatable :: arrayGlob2D(:,:)
      real*8, allocatable :: arrayGlob3D(:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Write_Netcdf_CloudGT called by ', procID

      if (iAmRootProc) then
         start2d(:) = (/ 1, rnum_out_cloud /)
         start3d(:) = (/ 1, 1, rnum_out_cloud /)

         count2d (:) = (/ NETCDF_HDF, 1 /)
         count3d (:) = (/ numLon, numLat, 1 /)

         hdr(1) = Nint (gmi_sec)
         hdr(2) = nymd
         hdr(3) = nhms

         call Ncwr_2d_Int (hdr, ncid_cloud, hdr_var_name, start2d, count2d)

         allocate(arrayGlob2D(i1_gl:i2_gl, ju1_gl:j2_gl))
         allocate(arrayGlob3D(i1_gl:i2_gl, ju1_gl:j2_gl, kx1:kx2))

         start4d(:) = (/ 1, 1, 1, rnum_out_cloud /)
         count4d (:) = (/ numLon, numLat, numVert, 1 /)
      end if

      call subDomain2Global (arrayGlob2D, psf_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2,        &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_rh_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_3d (arrayGlob2D, ncid_cloud, 'psf', start3d, count3d)

! temperature

      call subDomain2Global (arrayGlob3D, kel_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,        &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_rh_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'kel', start4d, count4d)

! total mass

      call subDomain2Global (arrayGlob3D, mass_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,        &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_rh_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'mass', start4d, count4d)

! grid box height

      call subDomain2Global (arrayGlob3D, gridbox_height_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,        &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_rh_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'gridbox_height', start4d, count4d)

! relative humidity

      call subDomain2Global (arrayGlob3D, cloud_rh_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,        &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_rh_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'RH', start4d, count4d)

! Cloud Fraction

      call subDomain2Global (arrayGlob3D, cloud_cf_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,     &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_cf_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'CF', start4d, count4d)

! Liquid Water Path

      call subDomain2Global (arrayGlob3D, cloud_lwc_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_lwc_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'LWC', start4d, count4d)

! Effective Radii

      call subDomain2Global (arrayGlob3D, cloud_r_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,   &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_r_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'R', start4d, count4d)

! Cloud Optical Depth

      call subDomain2Global (arrayGlob3D, cloud_cod_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &        rootProc, procID, map1_u, numDomains,  &
     &        SG_cloud_cod_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'COD', start4d, count4d)

! Cloud Optical Depth calculated with SO4

      call subDomain2Global (arrayGlob3D, cloud_cod_n_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &        rootProc, procID, map1_u, numDomains,  &
     &        SG_cloud_cod_n_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'COD_N', start4d, count4d)

! COD calculated with Fountoukis & Nenes, 2005 paramet

      call subDomain2Global (arrayGlob3D, cloud_cod_par_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,   &
     &        rootProc, procID, map1_u, numDomains,  &
     &        SG_cloud_cod_par_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'COD_par', start4d, count4d)

! Total Aerosol SO4

      call subDomain2Global (arrayGlob3D, cloud_so4_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,    &
     &        rootProc, procID, map1_u, numDomains,  &
     &        SG_cloud_so4_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'SO4', start4d, count4d)

! Cloud Droplet Number

      call subDomain2Global (arrayGlob3D, cloud_cdnc_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,   &
     &        rootProc, procID, map1_u, numDomains,  &
     &        SG_cloud_cdnc_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'CDNC', start4d, count4d)

! Cloud Droplet Number calculated by Fountoukis & Nenes, 2005 paramet

      call subDomain2Global (arrayGlob3D, cloud_cdnc_par_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &        rootProc, procID, map1_u, numDomains,  &
     &        SG_cloud_cdnc_par_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'CNDC_par', start4d, count4d)

! Effective Radii from SO4

      call subDomain2Global (arrayGlob3D, cloud_effr_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,   &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_effr_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'effR', start4d, count4d)

! Effective Radii from Fountoukis & Nenes, 2005 paramet

      call subDomain2Global (arrayGlob3D, cloud_effr_par_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_effr_par_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'effR_par', start4d, count4d)

! Cloud Albedo

      call subDomain2Global (arrayGlob3D, cloud_alb_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,   &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_alb_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'ALB', start4d, count4d)

! Cloud Albedo from SO4
      call subDomain2Global (arrayGlob3D, cloud_lwp_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,   &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_lwp_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'LWP', start4d, count4d)

! Cloud Albedo from  Fountoukis & Nenes, 2005 parame

      call subDomain2Global (arrayGlob3D, cloud_alb_par_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,   &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_alb_par_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'ALB_par', start4d, count4d)

! Updraft Velocity from  Fountoukis & Nenes, 2005 paramet

      call subDomain2Global (arrayGlob3D, cloud_w_par_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &     rootProc, procID, map1_u, numDomains,  SG_cloud_w_par_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'W_par', start4d, count4d)

! Maximum Supersaturation  from  Fountoukis & Nenes, 2005 paramet

      call subDomain2Global (arrayGlob3D, cloud_smax_par_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &        rootProc, procID, map1_u, numDomains, SG_cloud_smax_par_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'Smax_par', start4d, count4d)

! Autoconversion Rate from Khairoutdinov and Kogan, 2000 (Equation 29)

      call subDomain2Global (arrayGlob3D, Qaut_NS_KK1, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &        rootProc, procID, map1_u, numDomains, SG_Qaut_NS_KK1, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'Qaut_KK1', start4d, count4d)

! Autoconversion Rate from Khairoutdinov and Kogan, 2000 (Equation 30)

      call subDomain2Global (arrayGlob3D, Qaut_NS_KK2, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,  &
     &        rootProc, procID, map1_u, numDomains, SG_Qaut_NS_KK2, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'Qaut_KK2', start4d, count4d)

! Autoconversion Rate from Rotstayn, 1999 (Equation 15)

      call subDomain2Global (arrayGlob3D, Qaut_NS_R, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2, &
     &        rootProc, procID, map1_u, numDomains, SG_Qaut_NS_R, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'Qaut_R', start4d, count4d)

! AIE  !added on 071206

      call subDomain2Global (arrayGlob3D, cloud_AIE_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,    &
     &        rootProc, procID, map1_u, numDomains, SG_AIE_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'AIE', start4d, count4d)

! AIE12 !added on 071206

      call subDomain2Global (arrayGlob3D, cloud_AIE12_nc, &
     &        i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, kx1, kx2,   &
     &        rootProc, procID, map1_u, numDomains, SG_AIE12_nc, commuWorld)

      if (iAmRootProc) &
         call Ncwr_4d (arrayGlob3D, ncid_cloud, 'AIE12', start4d, count4d)

      if (iAmRootProc) call Ncdo_Sync (ncid_cloud)

      return

      end subroutine Write_Netcdf_CloudGT
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bufferOutput_CoudGT
!
! !INTERFACE:
!
      subroutine bufferOutput_CloudGT ( )
!
! !USES:
      use GmiBuffer_mod, only : bufferOutput4d_Nomean
!
      implicit none
!
! !DESCRIPTION:
! This routine buffers the cloud data to the output file. The worker processors
! send a slice of the data to the root processor, the root processor writes it 
! out, they then send another slice and it is written out, and so on.
!
! !DEFINED PARAMETERS:
      integer, parameter :: SG_Cloud_NC = 7007
!
!EOP
!------------------------------------------------------------------------------
!BOC

      call bufferOutput4d_Nomean  &
     &    (Cloud_var_name, kx1, kx2, commuWorld, SG_Cloud_NC,  &
     &     ncid_cloud, num_cloud, rnum_out_cloud, cloud_data,  &
     &     cloud_data_nc, map1_u, numDomains, rootProc, procID, &
     &     i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2)

      return

      end subroutine bufferOutput_CloudGT
!-------------------------------------------------------------------------------

end module GmiControlCloudModuleGT_mod
