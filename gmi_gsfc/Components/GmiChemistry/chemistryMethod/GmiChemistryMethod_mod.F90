!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiChemistryMethod_mod
!
! !INTERFACE:
!
#include "GmiESMF_ErrLog.h"
!
  module GmiChemistryMethod_mod
!
! !USES:
      use ESMF_Mod
      use GmiESMF_ErrorChecking_mod
!      use GmiSolar_mod,             ONLY : computeSolarZenithAngle_Photolysis
      use GmiESMFrcFileReading_mod, only : rcEsmfReadTable, rcEsmfReadLogical
      use m_netcdf_io_get_dimlen, only : Ncget_Dimlen
      use m_netcdf_io_read      , only : Ncrd_1d_Int
      use m_netcdf_io_close     , only : Nccl
      use m_netcdf_io_open      , only : Ncop_rd
      use GmiEmissionMethod_mod, only : t_Emission
      use GmiEmissionMethod_mod, only : Get_doReadDailyEmiss, Get_emiss_opt,   &
     &       Get_emiss_timpyr, Get_begDailyEmissRec, Get_endDailyEmissRec,     &
     &       Get_curEmissionFileRecord, Get_do_semiss_inchem
      use GmiEmissionMethod_mod, only : Get_do_ShipEmission
      use GmiGrid_mod, only : t_gmiGrid, Get_numSpecies, Get_i1, Get_i2,       &
     &       Get_ju1, Get_j2, Get_k1, Get_k2, Get_i1_gl, Get_i2_gl, Get_ju1_gl,&
     &       Get_j2_gl, Get_ilo, Get_ihi, Get_julo, Get_jhi, Get_ilong,        &
     &       Get_ilat, Get_ivert, Get_itloop
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID,          &
     &       Get_numDomains, Get_mcor, Get_latdeg, Get_londeg, Get_dlatr,      &
             Get_communicatorWorld
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_do_qqjk_reset,   &
     &       Get_constOutputFrequency, Set_do_qqjk_reset, Get_qj_var_name,     &
     &       Get_restart_inrec, Get_restart_infile_name, Get_do_mean,          &
     &       Get_pr_qj, Get_pr_diag, Get_pr_qqjk, Get_pr_smv2, Get_pr_qk,      &
     &       Get_pr_sad, Get_rd_restart, Get_pr_qj_opt_depth, Get_pr_qj_o3_o1d,&
     &       Get_do_qqjk_inchem, Get_do_ftiming, Get_pr_sulf_src, Get_pr_const,&
     &       Get_pr_emiss_3d, Get_pr_surf_emiss, Get_do_aerocom,               &
     &       Get_pr_ascii5, Get_pr_AerDust
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration,        &
     &       Get_concentration, Get_const_init_val, Set_concentration,         &
     &       Get_prod, Get_loss, Set_prod, Set_loss
      use GmiTimeControl_mod   , only : t_GmiClock, Get_gmiTimeStep,           &
     &       Get_begGmiDate, Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds,   &
     &       Get_numTimeSteps, GmiSplitDateTime, GetSecondsFromJanuary1
      use GmiSpeciesRegistry_mod, only : getSpeciesIndex, UNKNOWN_SPECIES
      use GmiPrintError_mod, only : GmiPrintError
      use GmiMetFieldsControl_mod, only : t_metFields, Get_met_opt, Get_mass, &
     &       Get_metdata_name_org, Get_metdata_name_model, Get_gridBoxHeight, &
     &       Get_tropopausePress, Get_relativeHumidity, Get_press3c,          &
     &       Get_press3e, Get_humidity, Get_pctm2, Get_moistq, Get_lwi_flags, &
     &       Get_totalCloudFraction, Get_fracCloudCover, Get_tau_cloud,       &
     &       Get_kel, Get_surf_air_temp, Get_surf_alb_uv, Get_cmf, Get_clwc,  &
     &       Get_radswg, Get_met_opt, Get_taucli, Get_tauclw
      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved, t_SulfSaved
      use GmiSolver_SavedVariables_mod, only : t_CloudParametersGT
      USE FastJX_Bundle_mod,      ONLY : t_fastJXbundle
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: InitializeChemistry, RunChemistry     , FinalizeChemistry
  public  :: initReadChemistry, runReadChemistry, readChemistryResourceFile
  public  :: runCalcAerosolDust, CalcAerosolDust
  public  :: Get_saldif_data, Get_saldir_data
  public  :: Get_sasdif_data, Get_sasdir_data, Get_uvalbedo_data
  public  :: Get_chem_opt       , Get_sad_opt      , Get_forc_bc_opt
  public  :: Get_lbssad_opt     , Get_h2oclim_opt  , Get_phot_opt
  public  :: Get_fastj_opt      , Get_sfalbedo_opt , Get_uvalbedo_opt
  public  :: Get_loss_opt       , Get_oz_eq_synoz_opt, Get_dehyd_opt
  public  :: Get_h2oclim_timpyr , Get_lbssad_timpyr, Get_sadgmi
  public  :: Get_fbc_j1         , Get_fbc_j2       , Get_forc_bc_years
  public  :: Get_loss_freq_opt  , Get_do_chem_grp  , Get_do_rxnr_adjust
  public  :: Get_do_clear_sky   , Get_t_cloud_ice  , Get_synoz_threshold
  public  :: Get_do_solar_cycle , Get_dust         , Get_wAersl
  public  :: Get_dAersl         , Get_optDepth     , Get_odAer
  public  :: Get_odMdust        , Get_eRadius      , Get_tArea
  public  :: Get_hno3cond       , Get_hno3gas      , Get_h2oback
  public  :: Get_h2ocond        , Get_reffice      , Get_reffsts
  public  :: Get_vfall          , Get_dms_no3      , Get_dms_oh
  public  :: Get_so2_h2o2       , Get_so2_o3       , Get_so2_oh
  public  :: Get_qkgmi          , Get_qqjda        , Get_qqkda 
  public  :: Get_qjgmi          , Get_qqjgmi       , Get_qqkgmi
  public  :: Get_rxnr_adjust    , Get_yda          , Get_overheadO3col
  public  :: Get_decay_3d_out  
  public  :: Set_optDepth       , Set_h2ocond
  public  :: Set_dms_no3        , Set_dms_oh
  public  :: Set_so2_h2o2       , Set_so2_o3       , Set_so2_oh
  public  :: Set_qkgmi          , Set_qqjda        , Set_qqkda 
  public  :: Set_qjgmi          , Set_qqjgmi       , Set_qqkgmi
  public  :: Set_rxnr_adjust    , Set_yda          , Set_overheadO3col
  public  :: Set_decay_3d_out
  public  :: Allocate_flux_gt    , Get_flux_gt    , Set_flux_gt    
  public  :: Allocate_cloud_param, Get_cloud_param, Set_cloud_param
  public  :: Allocate_cloud_tau  , Get_cloud_tau  , Set_cloud_tau  
  public  :: Get_cloudDroplet_opt  , Set_cloudDroplet_opt  
  public  :: Get_do_full_chem, Get_num_active, Get_num_chem, Get_num_molefrac
  public  :: Get_num_qjs, Get_num_qjo, Get_num_qks, Get_num_ks_sbc, Get_num_sad
  public  :: Get_num_spc_sbc, Get_surf_bc_map, Get_mw
  public  :: Get_const_labels, Get_qj_labels, Get_qk_labels
  public  :: Get_ibrono2_num, Get_ich4_num, Get_in2o_num, Get_idehyd_num
  public  :: Get_ih2_num, Get_ih2o_num, Get_ih2o2_num, Get_ih2oaircr_num
  public  :: Get_ihcl_num, Get_ihno3_num, Get_imgas_num, Get_initrogen_num
  public  :: Get_ioxygen_num, Get_ihno3cond_num, Get_ih2ocond_num, Get_io3_num
  public  :: Get_ino_num, Get_iacetone_num, Get_iisoprene_num, Get_ipropene_num
  public  :: Get_isynoz_num, Get_do_synoz, Get_do_nodoz, Get_dehydmin, Set_dehydmin
  public  :: Get_num_nox, Get_num_noy, Get_nox_map, Get_noy_map, Get_ico_num
  public  :: Get_do_AerDust_Calc

! !PUBLIC DATA MEMBERS:

  public  :: t_Chemistry

# include "GmiParameters.h"
!# include "gmi_time_constants.h"
# include "gmi_AerDust_const.h"
# include "setkin_par.h"
# include "setkin_lchem.h"
# include "setkin_mw.h"
# include "setkin_surf_bc.h"
 
! Parameter for loss frequency
  integer, parameter :: JDIM = 18
  integer, parameter :: KDIM = 20
  integer, parameter :: MDIM = 12

  type t_Chemistry
    private
    TYPE(t_CloudParametersGT) :: cloudGT
    TYPE(t_fastJXbundle)   :: JXbundle
    !---------------------------------
    type(t_Smv2Saved)   :: smv2SavedVars ! COMMON BLOCK, SAVE and DATA variables
    type(t_SulfSaved)   :: sulfSavedVars ! COMMON BLOCK, SAVE and DATA variables
    real*8              :: chemintv
    real*8              :: rateintv
    integer             :: jno_num  ! for adjusting J(NO)
    integer             :: jnoxnum
    !---------------------------------
    integer             :: prevEmissionFileRecord
    real*8, pointer     :: emiss_ozone(:,:) => null()      ! for Ship Emissions
    real*8              :: dehydmin ! minimum dehyd value (mixing ratio)
    integer             :: ibrono2_num, ich4_num, in2o_num, idehyd_num, ih2_num
    integer             :: ih2o_num, ih2o2_num, ih2oaircr_num, ihcl_num, ico_num
    integer             :: ihno3_num, imgas_num, initrogen_num, ioxygen_num
    integer             :: ihno3cond_num, ih2ocond_num, io3_num, isynoz_num
    integer             :: ino_num, iacetone_num, iisoprene_num, ipropene_num
    integer             :: num_nox ! number of Nodoz NOx species
    integer             :: num_noy ! number of Nodoz NOy species
    integer             :: nox_map(MAX_NUM_SMARRAY) ! mapping of NOx spec. # to const spc. #
    integer             :: noy_map(MAX_NUM_SMARRAY) ! mapping of NOy spec. # to const spc. #
    logical             :: do_synoz  ! do Synoz?
    logical             :: do_nodoz  ! do Nodoz?
    integer             :: cloudDroplet_opt
    real*8, pointer     :: cloud_param(:,:,:,:) => null()
    real*8, pointer     :: cloud_tau  (:,:,:) => null()
    real*8, pointer     :: flux_gt    (:,:,:) => null()
    real*8, pointer     :: tArea    (:,:,:,:) => null()
    real*8, pointer     :: odAer    (:,:,:,:) => null()
    real*8, pointer     :: odMdust  (:,:,:,:) => null()
    real*8, pointer     :: eRadius  (:,:,:,:) => null()
    real*8, pointer     :: optDepth (:,:,:,:) => null()
    real*8, pointer     :: dust     (:,:,:,:) => null()
    real*8, pointer     :: wAersl   (:,:,:,:) => null()
    real*8, pointer     :: dAersl   (:,:,:,:) => null()

    real*8, pointer     :: lossData (:,:,:,:) => null()
    real*8, pointer     :: sadgmi   (:,:,:,:) => null()
    real*8, pointer     :: lbssad   (:,:,:,:) => null()
    real*8, pointer     :: h2oclim  (:,:,:,:) => null()
    real*8, pointer     :: ch4clim  (:,:,:,:) => null()
    real*8, pointer     :: loss_freq(:,:,:,:) => null()
    real*8, pointer     :: hno3cond (:,:,:)   => null()
    real*8, pointer     :: hno3gas  (:,:,:)   => null()
    real*8, pointer     :: h2oback  (:,:,:)   => null()
    real*8, pointer     :: h2ocond  (:,:,:)   => null()
    real*8, pointer     :: reffice  (:,:,:)   => null()
    real*8, pointer     :: reffsts  (:,:,:)   => null()
    real*8, pointer     :: vfall    (:,:,:)   => null()

    real*8, pointer     :: dms_no3  (:,:,:)   => null()
    real*8, pointer     :: dms_oh   (:,:,:)   => null()
    real*8, pointer     :: so2_h2o2 (:,:,:)   => null()
    real*8, pointer     :: so2_o3   (:,:,:)   => null()
    real*8, pointer     :: so2_oh   (:,:,:)   => null()

    real*8, pointer     :: qkgmi   (:,:,:,:) => null()
    real*8, pointer     :: qqjda   (:,:,:,:) => null()
    real*8, pointer     :: qqkda   (:,:,:,:) => null()
    real*8, pointer     :: yda     (:,:,:,:) => null()
    real*8, pointer     :: qjgmi   (:,:,:,:) => null()
    real*8, pointer     :: qqjgmi  (:,:,:,:) => null()
    real*8, pointer     :: qqkgmi  (:,:,:,:) => null()
    real*8, pointer     :: rxnr_adjust(:,:,:,:,:)   => null()
    real*8, pointer     :: overheadO3col(:,:,:)   => null()
    real*8, pointer     :: decay_3d_out(:,:,:,:)  => null()
!
    integer             :: chem_opt
    real*8              :: chem_cycle
    integer             :: chem_mask_klo
    integer             :: chem_mask_khi
    real*8              :: synoz_threshold
    real*8              :: t_cloud_ice
!   
    integer             :: loss_opt
    integer             :: oz_eq_synoz_opt
!
    logical             :: do_chem_grp
    logical             :: do_smv_reord
    logical             :: do_wetchem
!
    integer             :: be_opt
    real*8              :: t_half_be7
    real*8              :: t_half_be10
    real*8              :: yield_be7
    real*8              :: yield_be10
!
    integer             :: forc_bc_opt
    integer             :: fbc_j1
    integer             :: fbc_j2
    integer             :: forc_bc_years
    integer             :: forc_bc_start_num
    integer             :: forc_bc_kmin
    integer             :: forc_bc_kmax
    integer             :: forc_bc_num
    integer             :: forc_bc_map (MAX_NUM_CONST)
    real*8              :: forc_bc_init_val
    real*8              :: forc_bc_incrpyr
    real*8              :: forc_bc_lz_val
    integer, pointer    :: jlatmd(:)
    integer             :: lastYearFBC
    ! forcing boundary condition data (ppmv)
    real*8, pointer     :: forc_bc_data(:,:,:,:)
    character (len=MAX_LENGTH_FILE_NAME) :: forc_bc_infile_name
    integer             :: loss_freq_opt
    integer             :: kmin_loss
    integer             :: kmax_loss
    real*8              :: loss_init_val
    character (len=MAX_LENGTH_FILE_NAME) :: loss_data_infile_name
!
    character (len=MAX_LENGTH_FILE_NAME) :: chem_infile_name
    character (len=MAX_LENGTH_FILE_NAME) :: linoz_infile_name
    character (len=MAX_LENGTH_FILE_NAME) :: sf6_infile_name
    character (len=MAX_LENGTH_FILE_NAME) :: SO3daily_infile_name
    character (len=MAX_LENGTH_FILE_NAME) :: SO3monthly_infile_name
    character (len=MAX_LENGTH_FILE_NAME) :: aqua_infile_name
    logical             :: do_aqu_chem
!
    integer             :: dehyd_opt
    integer             :: sad_opt
    integer             :: h2oclim_opt
    integer             :: h2oclim_timpyr
    real*8              :: ch4clim_init_val
    real*8              :: h2oclim_init_val
    character (len=MAX_LENGTH_FILE_NAME) :: h2oclim_infile_name
    integer             :: lbssad_opt
    integer             :: lbssad_timpyr
    real*8              :: lbssad_init_val
    character (len=MAX_LENGTH_FILE_NAME) :: lbssad_infile_name
!
    integer             :: num_rxnr_adjust
    integer             :: rxnr_adjust_timpyr
    integer, pointer    :: rxnr_adjust_map(:) => null()
    logical             :: do_rxnr_adjust
    character (len=MAX_LENGTH_FILE_NAME) :: rxnr_adjust_infile_name
    character (len=MAX_LENGTH_VAR_NAME)  :: rxnr_adjust_var_name
    integer             :: phot_opt
    integer             :: fastj_opt
    logical             :: do_clear_sky
    real*8              :: fastj_offset_sec
!
    real*8              :: qj_init_val
    integer             :: qj_timpyr
    character (len=MAX_LENGTH_FILE_NAME) :: qj_infile_name
!
    logical             :: do_full_chem
    real*8 , pointer    :: mw (:) => null()
    integer             :: num_sad
    integer             :: num_qjs
    integer             :: num_qjo
    integer             :: num_qks
    integer             :: num_chem
    integer             :: num_active
    integer             :: num_ks_sbc
    integer             :: num_spc_sbc
    integer             :: num_molefrac
    integer, pointer    :: surf_bc_map(:) => null()
    character(len=MAX_LENGTH_SPECIES_NAME), pointer :: const_labels(:) => null() ! constituent string labels
    character(len=MAX_LENGTH_LABELS) :: qj_labels(MAX_NUM_QJ) ! qj (photolysis) string labels
    character(len=MAX_LENGTH_LABELS) :: qk_labels(MAX_NUM_QK) ! qk (thermal)    string labels
!
    integer             :: sfalbedo_opt
    real*8              :: saldif_init_val
    real*8              :: saldir_init_val
    real*8              :: sasdif_init_val
    real*8              :: sasdir_init_val
            ! surface albedo data for diffuse nearIR (fraction 0-1)
    real*8, pointer     :: saldif_data(:,:,:) => null()
            ! surface albedo data for direct  nearIR (fraction 0-1)
    real*8, pointer     :: saldir_data(:,:,:) => null()
            ! surface albedo data for diffuse uv/vis (fraction 0-1)
    real*8, pointer     :: sasdif_data(:,:,:) => null()
            ! surface albedo data for direct  uv/vis (fraction 0-1)
    real*8, pointer     :: sasdir_data(:,:,:) => null()

    character (len=MAX_LENGTH_FILE_NAME) :: sfalbedo_infile_name
    integer             :: uvalbedo_opt
    real*8              :: uvalbedo_init_val
    real*8              :: jno_adjust
    real*8, pointer     :: uvalbedo_data(:,:,:) => null()
    character (len=MAX_LENGTH_FILE_NAME) :: uvalbedo_infile_name
    character (len=MAX_LENGTH_FILE_NAME) :: cross_section_file
    character (len=MAX_LENGTH_FILE_NAME) :: rate_file
    character (len=MAX_LENGTH_FILE_NAME) :: T_O3_climatology_file 
    character (len=MAX_LENGTH_FILE_NAME) :: scattering_data_file
    character (len=MAX_LENGTH_FILE_NAME) :: raa_qaa_data_file
    logical             :: do_solar_cycle
    character (len=MAX_LENGTH_FILE_NAME) :: sc_infile_name
    logical             :: do_ozone_inFastJX
!
    integer             :: AerDust_Effect_opt
    logical             :: do_AerDust_Calc
    character (len=MAX_LENGTH_FILE_NAME) :: AerDust_infile_name
  end type t_Chemistry

  real*8 :: RAA_b(4, NP_b), QAA_b(4, NP_b)
!
! !DESCRIPTION:
!
! !AUTHOR:
!
! !REVISION HISTORY:
!
!EOP
!-------------------------------------------------------------------------
  CONTAINS
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readChemistryResourceFile
!
! !INTERFACE:
!
      subroutine readChemistryResourceFile (self, Diagnostics, gmiGrid, gmiDomain, &
     &                config)
!
! !USES:
      use GmiCheckNamelistFile_mod, only : CheckNamelistOptionRange
      use GmiStringManipulation_mod, only : constructListNames
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiGrid    ), intent(in) :: gmiGrid    
      type(t_gmiDomain  ), intent(in) :: gmiDomain  
      type(t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Chemistry), intent(inOut) :: self
      type (ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in Chemistry related variables from the resource file.
!
! !LOCAL VARIABLES:
      integer :: numSpecies, procID, STATUS, RC, ic, emiss_opt
      integer :: ju1_gl, j2_gl, i1, i2, ju1, j2, k1, k2
      logical :: pr_diag, pr_qj_o3_o1d, pr_qj_opt_depth
      real    :: tempR4
      real*8  :: hugeReal
      character (len=MAX_LENGTH_SPECIES_NAME), pointer :: tempListNames(:)
      character (len=MAX_STRING_LENGTH      ) :: forcedBcSpeciesNames
      character(len=ESMF_MAXSTR) :: IAm, err_msg
!
!EOP
!------------------------------------------------------------------------------
!BOC
      IAm = "readChemistryResourceFile"

      call Get_procID (gmiDomain, procID)
      call Get_pr_diag (Diagnostics, pr_diag)

      if (pr_diag) Write(6,*) IAm, 'called by ', procID

      call Get_i1(gmiGrid, i1)
      call Get_i2(gmiGrid, i2)
      call Get_j2(gmiGrid, j2)
      call Get_ju1(gmiGrid, ju1)
      call Get_ju1_gl(gmiGrid, ju1_gl)
      call Get_j2_gl (gmiGrid, j2_gl)
      call Get_k1 (gmiGrid, k1)
      call Get_k2 (gmiGrid, k2)
      call Get_numSpecies (gmiGrid, numSpecies )

      call Get_pr_qj_o3_o1d(Diagnostics, pr_qj_o3_o1d)
      call Get_pr_qj_opt_depth(Diagnostics, pr_qj_opt_depth)

      allocate(tempListNames(numSpecies))

      !################################
      ! Begin reading the resource file
      !################################

!     -------------------------------------
!     chem_opt
!       0:  no chemistry
!       1:  call Radon/Lead chemistry
!       2:  call Smvgear2
!       3:  call simple loss
!       4:  call forcing boundary condition
!       5:  call Synoz tracer
!       6:  call Beryllium chemistry
!       7:  call Quadchem
!       8:  call Sulfur chemistry
!       9:  call Tracer Package (Radon/Lead and Beryllium) chemistry
!     -------------------------------------

      call ESMF_ConfigGetAttribute(config, self%chem_opt, &
     &                label   = "chem_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

!     -------------------------------------------------------------
!     chem_cycle:  number of time steps to cycle chemistry calls on
!       < 1.0:  chemistry will subcycle
!         1.0:  chemistry called each time step
!     -------------------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%chem_cycle, &
     &                label   = "chem_cycle:", &
     &                default = 1.0d0, rc=STATUS )
      VERIFY_(STATUS)

!     -----------------------------------------------------
!     chem_mask_klo, chem_mask_khi:
!       chemistry turned off where k is outside of range of
!       [chem_mask_klo, chem_mask_khi]
!     -----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%chem_mask_klo, &
     &                label   = "chem_mask_klo:", &
     &                default = k1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%chem_mask_khi, &
     &                label   = "chem_mask_khi:", &
     &                default = k2, rc=STATUS )
      VERIFY_(STATUS)

!     -------------------------------------------------------------------
!     synoz_threshold:  chemistry turned off where synoz > this threshold
!     -------------------------------------------------------------------

      hugeReal = Huge (hugeReal)

      call ESMF_ConfigGetAttribute(config, self%synoz_threshold, &
     &                label   = "synoz_threshold:", &
     &                default = hugeReal, rc=STATUS )
      VERIFY_(STATUS)

!     -----------------------------------------------------
!     loss_opt
!       0:  do not use stratospheric loss in gmi_step.F
!       1:  use stratospheric loss in gmi_step.F
!     -----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%loss_opt, &
     &                label   = "loss_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

!    cloudDroplet_opt = 1: Boucher and LohMan    Correlation (default)
!                     = 2: Nenes and Seinfeld    Parameterization
!                     = 3: Abdul-Razzak and Ghan Parameterization
!                     = 4: Segal amd Khain       Correllation
!   The variable is only use when the GT cloud module is employed.

      call ESMF_ConfigGetAttribute(config, self%cloudDroplet_opt, &
     &                label   = "cloudDroplet_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

!     -----------------------------------------------------
!     oz_eq_synoz_opt
!       0:  do not use stratospheric loss in gmi_step.F
!       1:  use stratospheric loss in gmi_step.F
!     -----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%oz_eq_synoz_opt, &
     &                label   = "oz_eq_synoz_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

!     --------------------------------------
!     Aerosols and Sulfur from Penner et al.
!     --------------------------------------

      call rcEsmfReadLogical(config, self%do_aqu_chem, "do_aqu_chem:", &
     &                       default=.false., rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%aqua_infile_name, &
     &                label   = "aqua_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%chem_infile_name, &
     &                label   = "chem_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%linoz_infile_name, &
     &                label   = "linoz_infile_name:", &
     &                default = 'linoz_input_prather.nc', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%sf6_infile_name, &
     &                label   = "sf6_infile_name:", &
     &                default = 'SF6_CFC115_loss.nc', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%SO3daily_infile_name, &
     &                label   = "SO3daily_infile_name:", &
     &                default = 'gmic_MERRA_2004_daily_tracer_input.nc', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%SO3monthly_infile_name, &
     &                label   = "SO3monthly_infile_name:", &
     &                default = 'gmic_MERRA_2004_monthly_tracer_input.nc', rc=STATUS )
      VERIFY_(STATUS)

!     ------------------------------------------------------------
!     t_cloud_ice:  temperature for ice formation in clouds (degK)
!     ------------------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%t_cloud_ice, &
     &                label   = "t_cloud_ice:", &
     &                default = 263.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%do_chem_grp, "do_chem_grp:", &
     &                       default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%do_smv_reord, "do_smv_reord:", &
     &                       default=.false., rc=STATUS)

      call rcEsmfReadLogical(config, self%do_wetchem, "do_wetchem:", &
     &                       default=.false., rc=STATUS)

!     -----------------------------------------
!     be_opt
!       1:  use Koch  table  for Be-7 and Be-10
!       2:  use Nagai tables for Be-7 and Be-10
!     -----------------------------------------

      call ESMF_ConfigGetAttribute(config, self%be_opt, &
     &                label   = "be_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%t_half_be7, &
     &                label   = "t_half_be7:", &
     &                default = 53.3d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%t_half_be10, &
     &                label   = "t_half_be10:", &
     &                default = 5.84d8, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%yield_be7, &
     &                label   = "yield_be7:", &
     &                default = 4.5d-7, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%yield_be10, &
     &                label   = "yield_be10:", &
     &                default = 2.5d-7, rc=STATUS )
      VERIFY_(STATUS)

!     ---------------------------
!     Forcing boundary condition:
!     ---------------------------

!     --------------------------------------------------
!     forc_bc_opt
!       1:  set all   forc_bc values to forc_bc_init_val
!       2:  read in   forc_bc
!       3:  calculate forc_bc
!     --------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%forc_bc_opt, &
     &                label   = "forc_bc_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%fbc_j1, &
     &                label   = "fbc_j1:", &
     &                default = ju1_gl, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%fbc_j2, &
     &                label   = "fbc_j2:", &
     &                default = j2_gl, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%forc_bc_years, &
     &                label   = "forc_bc_years:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%forc_bc_start_num, &
     &                label   = "forc_bc_start_num:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%forc_bc_kmin, &
     &                label   = "forc_bc_kmin:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%forc_bc_kmax, &
     &                label   = "forc_bc_kmax:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      self%forc_bc_map(:)      = 0

      call rcEsmfReadTable(config, forcedBcSpeciesNames, &
     &                     "forcedBcSpeciesNames::", rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%forc_bc_init_val, &
     &                label   = "forc_bc_init_val:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%forc_bc_incrpyr, &
     &                label   = "forc_bc_incrpyr:", &
     &                default = 0.3d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%forc_bc_lz_val, &
     &                label   = "forc_bc_lz_val:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%forc_bc_infile_name, &
     &                label   = "forc_bc_infile_name:", &
     &                default = 'forc_bc_co2.asc', rc=STATUS )
      VERIFY_(STATUS)

!     ------------
!     Simple loss:
!     ------------

!     -----------------------------------------------
!     loss_freq_opt
!       1:  set all loss_freq values to loss_init_val
!       2:  read in loss data
!       3:  use NCAR loss
!     -----------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%loss_freq_opt, &
     &                label   = "loss_freq_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%kmin_loss, &
     &                label   = "kmin_loss:", &
     &                default = k1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%kmax_loss, &
     &                label   = "kmax_loss:", &
     &                default = k2, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%loss_init_val, &
     &                label   = "loss_init_val:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)


      call ESMF_ConfigGetAttribute(config, self%loss_data_infile_name, &
     &                label   = "loss_data_infile_name:", &
     &                default = 'loss_n2o.asc', rc=STATUS )
      VERIFY_(STATUS)

      
!     ---------------------------
!     Surface Area Density (SAD):
!     ---------------------------

!     ----------------------------------------------------
!     sad_opt
!       0:  do not allocate or process SAD array
!       1:  allocate, but zero out SAD array
!       2:  call Considine code (i.e., Condense)
!       3:  read SAD array from a file of monthly averages
!     ----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%sad_opt, &
     &                label   = "sad_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

!     ------------------------------------------------
!     h2oclim_opt 
!       1:  set all h2oclim values to h2oclim_init_val
!       2:  read in h2oclim
!       3:  h2oclim, ch4clim not used.  Instead, transported
!           H2O and CH4 are provided by the host AGCM.
!     ------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%h2oclim_opt, &
     &                label   = "h2oclim_opt:", &
     &                default = 2, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%h2oclim_timpyr, &
     &                label   = "h2oclim_timpyr:", &
     &                default = MONTHS_PER_YEAR, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%ch4clim_init_val, &
     &                label   = "ch4clim_init_val:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%h2oclim_init_val, &
     &                label   = "h2oclim_init_val:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%h2oclim_infile_name, &
     &                label   = "h2oclim_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

!     ------------------------------------------------
!     dehyd_opt
!       0:  dehyd is not available and is set to zero.
!           In general this is the AGCM case where H20 and
!           condensed H2O are transported species.
!       1:  dehyd is a transported specie and is available.
!           In general, this is the CTM case.
!     ------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%dehyd_opt, &
     &                label   = "dehyd_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

!     ----------------------------------------------
!     lbssad_opt
!       1:  set all lbssad values to lbssad_init_val
!       2:  read in lbssad 3d fields
!       3:  read in lbssad zonal average fields
!       4:  lbssad provided by AGCM, current month only.
!     ----------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%lbssad_opt, &
     &                label   = "lbssad_opt:", &
     &                default = 2, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%lbssad_timpyr, &
     &                label   = "lbssad_timpyr:", &
     &                default = MONTHS_PER_YEAR, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%lbssad_init_val, &
     &                label   = "lbssad_init_val:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%lbssad_infile_name, &
     &                label   = "lbssad_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

!     -------------------------
!     Reaction rate adjustment:
!     -------------------------

      call rcEsmfReadLogical(config, self%do_rxnr_adjust, "do_rxnr_adjust:", &
     &                       default=.false., rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%rxnr_adjust_infile_name, &
     &                label   = "rxnr_adjust_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%rxnr_adjust_var_name, &
     &                label   = "rxnr_adjust_var_name:", &
     &                default = 'reac_rate_adj', rc=STATUS )
      VERIFY_(STATUS)


!     =========
!     nlGmiPhotolysis
!     =========

!     -----------------------------------------------------
!     phot_opt
!       0:  no photolysis
!       1:  set all qj values to qj_init_val
!       2:  read in qj values
!       3:  use fastj routine (for fastJ, fastJx, fastJx53b)
!           This option should be combined with fastj_opt.
!       4:  lookup table for qj (Kawa style)
!       5:  lookup table for qj (Kawa style) +
!           use ozone climatology for column ozone calc.
!       6:  calculate from table and Gmimod data (Quadchem)
!       7:  read in qj values (2-D, 12 months)
!       8:  use fast-JX routine (troposphere/stratosphere)
!     -----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%phot_opt, &
     &                label   = "phot_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

!     -----------------------------------------------------
!     fastj_opt: set when phot_opt=3
!     0: for fastJ
!     1: for fastJx
!     2: for fastJx53b
!     3: for fastJx53c
!    -----------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%fastj_opt, &
     &                label   = "fastj_opt:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%do_clear_sky, "do_clear_sky:", &
     &                       default=.true., rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%fastj_offset_sec, &
     &                label   = "fastj_offset_sec:", &
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

!... do solar cycle in incoming solar flux?

     call rcEsmfReadLogical(config, self%do_solar_cycle, "do_solar_cycle:", &
     &                       default=.false., rc=STATUS)

      call ESMF_ConfigGetAttribute(config, self%sc_infile_name, &
     &                label   = "sc_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

!     ---------
!     qj / qqj:
!     ---------
      
      call ESMF_ConfigGetAttribute(config, self%qj_init_val, &
     &                label   = "qj_init_val:", &
     &                default = 1.0d-30, rc=STATUS )
      VERIFY_(STATUS)

      ! sets of photolysis per year (1 => yearly, 12 => monthly)

      call ESMF_ConfigGetAttribute(config, self%qj_timpyr, &
     &                label   = "qj_timpyr:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%qj_infile_name, &
     &                label   = "qj_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      
!     -------
!     albedo:
!     -------
      
!     ----------------------------------------------------------------------
!     sfalbedo_opt
!       0:  no sfalbedo
!       1:  set each type of sfalbedo to an intial value
!       2:  read in monthly sfalbedo values from a NetCDF file
!       3:  read in values of four types of surface albedo from the met data
!     ----------------------------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%sfalbedo_opt, &
     &                label   = "sfalbedo_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%saldif_init_val, &
     &                label   = "saldif_init_val:", &
     &                default = 0.1d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%saldir_init_val, &
     &                label   = "saldir_init_val:", &
     &                default = 0.1d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%sasdif_init_val, &
     &                label   = "sasdif_init_val:", &
     &                default = 0.1d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%sasdir_init_val, &
     &                label   = "sasdir_init_val:", &
     &                default = 0.1d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%sfalbedo_infile_name, &
     &                label   = "sfalbedo_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

!     --------------------------------------------------------
!     uvalbedo_opt
!       0:  no uvalbedo
!       1:  set all uvalbedo values to uvalbedo_init_val
!       2:  read in monthly uvalbedo values from an ASCII file
!       3:  read in surface albedo values from the met data
!     --------------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%uvalbedo_opt, &
     &                label   = "uvalbedo_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%uvalbedo_init_val, &
     &                label   = "uvalbedo_init_val:", &
     &                default = 0.1d0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%uvalbedo_infile_name, &
     &                label   = "uvalbedo_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%cross_section_file, &
     &                label   = "cross_section_file:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%rate_file, &
     &                label   = "rate_file:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%T_O3_climatology_file, &
     &                label   = "T_O3_climatology_file:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%scattering_data_file, &
     &                label   = "scattering_data_file:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)
!... file with fastj(x) parameters to calc aerosol opt depth
      call ESMF_ConfigGetAttribute(config, self%raa_qaa_data_file, &
     &                label   = "raa_qaa_data_file:", &
     &                default = 'raa_qaa_data_file', rc=STATUS )
      VERIFY_(STATUS)

!... parameter to adjust scaling of J(NO) in older fastJ/fastJX
      call ESMF_ConfigGetAttribute(config, self%jno_adjust, &
     &                label   = "jno_adjust:", &
     &                default=1.0d0, rc=STATUS)

      call rcEsmfReadLogical(config, self%do_ozone_inFastJX, &
     &              "do_ozone_inFastJX:", default=.false., rc=STATUS)

      !=================================================================
      ! do_AerDust_Calc is used to detrmine if aerosol/dust calculations
      !                 are done in the code. If set to FALSE, the code
      !                 will not read global aerosol/dust concentrations
      !                 and not do any aerosol/dust calculations.
      !=================================================================

      call rcEsmfReadLogical(config, self%do_AerDust_Calc, "do_AerDust_Calc:", &
     &                       default=.false., rc=STATUS)

      !=================================================================
      ! AerDust_Effect_opt is used to select if the radiative effects
      !                    or/and heterogeneous chemistry on different
      !                    aerosols/dust are turned on/off.
      !     0: radiative effects on  and heterogeneous chemistry on
      !     1: radiative effects off and heterogeneous chemistry on
      !     2: radiative effects on  and heterogeneous chemistry off
      !     3: radiative effects off and heterogeneous chemistry off
      !=================================================================

      call ESMF_ConfigGetAttribute(config, self%AerDust_Effect_opt, &
     &                label   = "AerDust_Effect_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(config, self%AerDust_infile_name, &
     &                label   = "AerDust_infile_name:", &
     &                default = '', rc=STATUS )
      VERIFY_(STATUS)

      allocate(self%mw(numSpecies))
      self%mw(1:numSpecies) = mw_data(1:numSpecies)

      call rcEsmfReadTable(config, self%mw, "mw::", rc=STATUS)

      allocate(self%const_labels(numSpecies))

      self%const_labels(1:numSpecies) = lchemvar(1:numSpecies)

      call rcEsmfReadTable(config, self%const_labels, &
     &                     "const_labels::", rc=STATUS)

    ! ---------------------------------------------------------------
    ! Check option ranges.  Note that as new options are added, these
    ! range checks will have to be modified.
    ! ---------------------------------------------------------------
   
      call CheckNamelistOptionRange ('dehyd_opt', self%dehyd_opt, 0, 1)
      call CheckNamelistOptionRange ('chem_opt', self%chem_opt, 0, 9)
      call CheckNamelistOptionRange ('be_opt', self%be_opt, 1, 2)
      call CheckNamelistOptionRange ('forc_bc_opt', self%forc_bc_opt, 1, 3)
      call CheckNamelistOptionRange ('h2oclim_opt', self%h2oclim_opt, 1, 3)
      call CheckNamelistOptionRange ('lbssad_opt', self%lbssad_opt, 1, 4)
      call CheckNamelistOptionRange ('phot_opt', self%phot_opt, 0, 7)
      call CheckNamelistOptionRange ('fastj_opt', self%fastj_opt, 0, 4)
      call CheckNamelistOptionRange ('sad_opt', self%sad_opt, 0, 3)
      call CheckNamelistOptionRange ('uvalbedo_opt', self%uvalbedo_opt, 0, 3)
      call CheckNamelistOptionRange ('sfalbedo_opt', self%sfalbedo_opt, 0, 3)
      call CheckNamelistOptionRange ('loss_opt', self%loss_opt, 0, 1)
      call CheckNamelistOptionRange ('loss_freq_opt', self%loss_freq_opt, 1, 3)
      call CheckNamelistOptionRange ('oz_eq_synoz_opt', self%oz_eq_synoz_opt, 0, 1)
      call CheckNamelistOptionRange ('cloudDroplet_opt', self%cloudDroplet_opt, 1, 4)
      call CheckNamelistOptionRange ('AerDust_Effect_opt', self%AerDust_Effect_opt, 0, 3)

!
      if ((self%chem_opt == 2) .or. (self%chem_opt == 7) .or. &
     &    (self%chem_opt == 8)) then
          self%do_full_chem = .true.

          self%mw          (1:numSpecies) = mw_data (1:numSpecies)
          self%const_labels(1:numSpecies) = lchemvar(1:numSpecies)
      else
          self%do_full_chem = .false.

         if ((self%mw(1) == 0.0d0) .and. (numSpecies == NSP)) then
            self%mw(1:numSpecies) = mw_data(1:numSpecies)
         end if
      end if

      ! Set the initial value of the list
      tempListNames(:) = ''

      ! Construct the list of names using the long string
      call constructListNames(tempListNames, forcedBcSpeciesNames)

      self%forc_bc_num = Count (tempListNames(:) /= '')
      if (self%forc_bc_num > 0) then
         do ic = 1, self%forc_bc_num
            self%forc_bc_map(ic) = getSpeciesIndex(tempListNames(ic))
         end do
      end if

      IF (self%lbssad_opt == 4) self%lbssad_timpyr = 1
!
      self%qj_labels(:) = ' '
      self%qk_labels(:) = ' '

      if (self%do_full_chem) then
         self%num_active    = NACT
         self%num_chem      = NCHEM
         self%num_molefrac  = NMF
         self%num_qjs       = NUM_J
         self%num_qjo       = NUM_J
         self%num_qks       = NUM_K

         self%io3_num       = IO3
         self%ibrono2_num   = IBRONO2
         self%ich4_num      = ICH4
         self%in2o_num      = IN2O
         self%idehyd_num    = IDEHYD
         self%ih2_num       = IH2
         self%ih2o_num      = IH2O
         self%ih2o2_num     = IH2O2
         self%ih2oaircr_num = IH2OAIR
         self%ihcl_num      = IHCL
         self%ihno3_num     = IHNO3
         self%imgas_num     = IMGAS
         self%initrogen_num = INITROGEN
         self%ioxygen_num   = IOXYGEN
         self%ico_num       = ICO
         self%ino_num       = INO
         self%ih2ocond_num  = 0 ! getSpeciesIndex('')
         self%ihno3cond_num = 0 ! getSpeciesIndex('')

         call ESMF_ConfigGetAttribute(config, emiss_opt, &
     &                label   = "emiss_opt:", &
     &                default = 0, rc=STATUS )
         VERIFY_(STATUS)

         if (btest(emiss_opt,1)) then
            self%iacetone_num  = IC3H6O
            self%ico_num       = ICO
            self%iisoprene_num = IC5H8
            self%ipropene_num  = IC3H6
            self%ino_num       = INO
         end if

         if (pr_qj_o3_o1d)    self%num_qjo = self%num_qjo + 1
         if (pr_qj_opt_depth) self%num_qjo = self%num_qjo + 1

         self%num_ks_sbc    = K_SBC

         if (NUM_SBC > 0) then
            allocate(self%surf_bc_map(1:NUM_SBC))
            self%surf_bc_map(:)         = 0
            self%surf_bc_map(1:NUM_SBC) = sbc_map(1:NUM_SBC)

            if (self%surf_bc_map(1) /= 0) then
              self%num_spc_sbc = NUM_SBC
            else
              self%num_spc_sbc = 0
            end if
         end if

         if (self%num_qjs > 0) self%qj_labels(1:self%num_qjs) = lqjchem(1:self%num_qjs)
         if (self%num_qks > 0) self%qk_labels(1:self%num_qks) = lqkchem(1:self%num_qks)

         if (pr_qj_o3_o1d)    self%qj_labels(self%num_qjs+1) = 'O3 + hv = O1D + O2'
         if (pr_qj_opt_depth) self%qj_labels(self%num_qjo)   = 'optical depth'

         if (self%do_chem_grp) then
            call setupChemicalGroup (i1, i2, ju1, j2, k1, k2)
         end if

!     ====
      else
!     ====

         self%num_chem      = numSpecies
         self%num_molefrac  = numSpecies
         self%num_qjs       = 0
         self%num_qks       = 0

         self%num_ks_sbc    = 0
         self%num_spc_sbc   = 0

         self%ih2o2_num     = 0
         self%ih2_num       = 0
         self%ih2o_num      = 0
         self%ihno3_num     = 0
         self%idehyd_num    = 0
         self%imgas_num     = 0
         self%ih2oaircr_num = 0
         self%ich4_num      = 0
         self%ih2ocond_num  = 0
         self%ihno3cond_num = 0
         self%ico_num       = 0
         self%ino_num       = 0
         self%initrogen_num = 0
         self%ioxygen_num   = 0
         self%ihcl_num      = 0
         self%ibrono2_num   = 0
         self%in2o_num      = 0

         call ESMF_ConfigGetAttribute(config, self%io3_num, &
     &                label   = "io3_num:", &
     &                default = 0, rc=STATUS )
         VERIFY_(STATUS)

         call ESMF_ConfigGetAttribute(config, self%iacetone_num, &
     &                label   = "iacetone_num:", &
     &                default = 0, rc=STATUS )
         VERIFY_(STATUS)

         call ESMF_ConfigGetAttribute(config, self%ipropene_num, &
     &                label   = "ipropene_num:", &
     &                default = 0, rc=STATUS )
         VERIFY_(STATUS)

         call ESMF_ConfigGetAttribute(config, self%iisoprene_num, &
     &                label   = "iisoprene_num:", &
     &                default = 0, rc=STATUS )
         VERIFY_(STATUS)
      end if

      self%do_synoz   = .false.
      self%do_nodoz   = .false.

      self%isynoz_num = 0

      self%num_nox    = 0
      self%num_noy    = 0

      self%nox_map(:) = 0
      self%noy_map(:) = 0

      if (self%do_full_chem .or. (self%chem_opt == 5) .or. (self%chem_opt == 9)) then
         call Setup_Syn_Spc (self, numSpecies)
      end if

      if (self%sad_opt /= 0) then
         self%num_sad = NSAD
      else
         self%num_sad = 0
      end if

! When NOT using climatological H2O and CH4, the dehydrated specie is absent.
! --------------------------------------------------------------------------
      IF (self%dehyd_opt == 0 .AND. self%h2oclim_opt /= 3) THEN
         PRINT *,"readChemistryResourceFile: dehyd_opt and h2oclim_opt incompatible."
         STOP
      END IF

      return

      end subroutine readChemistryResourceFile
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InitializeChemistry
!
! !INTERFACE:
!
      subroutine InitializeChemistry (self, Emission, gmiGrid, gmiDomain, &
     &                     gmiClock, Diagnostics, metFields, chem_mecha)
!
! !USES:
      use fastj                     , only : InitializeFastj
      use fast_JX                   , only : InitializeFastJX
      use fastjx65_mod              , only : InitializeFastJX65 
      use Fast_JX53b                , only : InitializeFastJX53b
      use Fast_JX53c                , only : InitializeFastJX53c
      use ReadForcedBC_mod          , only : readForcedBcData
      use ReadUValbedoData_mod      , only : readUValbedoData
      use ReadLossFrequency_mod     , only : readLossFrequency
      use ReadPhotolysisRates_mod   , only : readPhotolysisRates, readRaaQaa
      use ReadLiqBinarySulfate_mod  , only : readLiqBinarySulfate
      use ReadWaterClimatology_mod  , only : readWaterClimatology
      use ReadReactionRateAdjFac_mod, only : readReactionRateAdjustmentFactors
!
      implicit none
!
#     include "phot_lookup.h"
#     include "phot_monthly.h"
#     include "gmi_forc_bc.h"
#     include "phot_lookup_constants.h"
#     include "phot_lookup_arrays.h"
#     include "gmi_phys_constants.h"
!
! !INPUT PARAMETERS:
      character (len=*)  , intent(in) :: chem_mecha
      type(t_gmiGrid    ), intent(in) :: gmiGrid  
      type(t_Emission   ), intent(in) :: Emission 
      type(t_metFields  ), intent(in) :: metFields
      type(t_gmiDomain  ), intent(in) :: gmiDomain
      type(t_gmiClock   ), intent(in) :: gmiClock 
      type(t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type  (t_Chemistry), intent(inOut) :: self
!
! !LOCAL VARIABLES:
      integer :: ncid_rra, cnt1d(1), strt1d(1)
      integer :: i1, i2, i1_gl, i2_gl, ju1, j2, ju1_gl, j2_gl
      integer :: k1, k2, itloop, ic, ilong, ilat, ivert, ij, jj, jx
      integer :: numSpecies, procID, emiss_opt, met_opt
      logical :: pr_const, pr_diag, pr_sulf_src, pr_qqjk, pr_qj, pr_qk, pr_AerDust
      logical :: do_qqjk_inchem, do_full_chem, pr_sad, do_mean
      logical :: do_semiss_inchem, rootProc, pr_smv2, do_ShipEmission
      real*8  :: tdt, rjx
      character (len=MAX_LENGTH_FILE_NAME) :: smv_filnam
      character (len=75) :: err_msg
      character(len=MAX_LENGTH_VAR_NAME) :: qj_var_name
      integer :: pr_nc_period
      real*8 , allocatable :: dlatr(:)
      integer :: idumday, month, year, start_ymd
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID  (gmiDomain, procID )
      call Get_pr_diag (Diagnostics, pr_diag )

      if (pr_diag) Write(6,*) 'InitializeChemistry called by ', procID

      call Get_pr_qj          (Diagnostics, pr_qj)
      call Get_pr_qk          (Diagnostics, pr_qk)
      call Get_pr_sad         (Diagnostics, pr_sad)
      call Get_pr_qqjk        (Diagnostics, pr_qqjk)
      call Get_do_mean        (Diagnostics, do_mean)
      call Get_pr_const       (Diagnostics, pr_const)
      call Get_pr_sulf_src    (Diagnostics, pr_sulf_src)
      call Get_do_qqjk_inchem (Diagnostics, do_qqjk_inchem)

      call Get_emiss_opt       (Emission, emiss_opt)
      call Get_do_semiss_inchem(Emission, do_semiss_inchem)

      call Get_i1    (gmiGrid, i1)
      call Get_i2    (gmiGrid, i2)
      call Get_ju1   (gmiGrid, ju1)
      call Get_j2    (gmiGrid, j2)
      call Get_k1    (gmiGrid, k1)
      call Get_k2    (gmiGrid, k2)
      call Get_i1_gl (gmiGrid, i1_gl)
      call Get_i2_gl (gmiGrid, i2_gl)
      call Get_ju1_gl(gmiGrid, ju1_gl)
      call Get_j2_gl (gmiGrid, j2_gl )
      call Get_numSpecies (gmiGrid, numSpecies )

      call Get_met_opt (metFields, met_opt )

      !###############
      ! Error Checking
      !###############

      call rcCheckChemistrySetting(self, pr_sad, pr_qj, pr_qk, pr_qqjk, &
     &       do_qqjk_inchem, do_mean, emiss_opt, do_semiss_inchem, numSpecies, &
     &       met_opt)

#ifdef GTmodule
       call Allocate_flux_gt    (self, i1,i2, ju1,j2, k1,k2)
       call Allocate_cloud_tau  (self, i1,i2, ju1,j2, k1,k2)
       call Allocate_cloud_param(self, i1,i2, ju1,j2, k1,k2, numSpecies)
#endif

      if (self%chem_opt /= 0) then
         if (self%do_full_chem) then
            call Allocate_qkgmi(self, i1,i2, ju1,j2, k1,k2)
            if (pr_qqjk) then
               if (do_qqjk_inchem) then
                  call Allocate_qqjda(self, i1,i2, ju1,j2, k1,k2)
                  call Allocate_qqkda(self, i1,i2, ju1,j2, k1,k2)
                  call Allocate_yda  (self, i1,i2, ju1,j2, k1,k2)
               else
                  call Allocate_qqjgmi(self, i1,i2, ju1,j2, k1,k2)
                  call Allocate_qqkgmi(self, i1,i2, ju1,j2, k1,k2)
               end if
            end if
         end if
      end if

      Allocate (self%decay_3d_out(i1:i2, ju1:j2, k1:k2, 1:numSpecies))
      self%decay_3d_out = 0.0d0

      if (self%phot_opt /= 0) then
         call Allocate_qjgmi(self, i1,i2, ju1,j2, k1,k2)
         call Allocate_overheadO3col(self, i1,i2, ju1,j2, k1,k2)

         if ((self%phot_opt == 4) .or. (self%phot_opt == 5)) then
            Allocate (rad_source(NUMLAM, NUMSZA, NUMO3, NUMPRS))
            rad_source = 0.0d0

            if (TRIM(chem_mecha) == 'strat_trop' .or. &
                TRIM(chem_mecha) == 'strat_trop_aerosol' ) then
               Allocate (o1d_coef(5, NUMSZA, NUMO3, NUMPRS, 2))
               o1d_coef   = 0.0d0
            end if

            if (self%phot_opt == 5) then
               Allocate (o3_clim(NUM_O3CLIM_PRS, i1:i2, ju1:j2, NUM_O3CLIM_MON))
               o3_clim = 0.0d0
            end if
         end if
      endif


      if ((self%uvalbedo_opt == 1) .or. (self%uvalbedo_opt == 2)) then
         Allocate (self%uvalbedo_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         self%uvalbedo_data = 0.0d0

         if (self%uvalbedo_opt == 1) then
            self%uvalbedo_data(:,:,:) = self%uvalbedo_init_val
         else if (self%uvalbedo_opt == 2) then
            call readUValbedoData (self%uvalbedo_data,                         &
     &               self%uvalbedo_infile_name, i1, i2, ju1, j2, i1_gl, i2_gl, &
     &               ju1_gl, j2_gl, pr_diag, procID)
         end if
      end if

      if (self%phot_opt == 3) then
                rootProc = .FALSE.
                if (procID==0) rootProc = .TRUE.

         !===========================
         select case (self%fastj_opt)
         !===========================

            !=======
            case (0)
            !=======
                call InitializeFastJ (self%cross_section_file, self%rate_file, &
     &                   self%T_O3_climatology_file, num_qj_o3_to_2oh, self%num_qjs,&
     &                   self%chem_mask_khi, k2, k1)

            !=======
            case (1)
            !=======
                call InitializeFastJX (self%JXbundle, self%cross_section_file, &
                         self%rate_file, self%T_O3_climatology_file, &
                         num_qj_o3_to_2oh, self%num_qjs, &
     &                   self%chem_mask_khi, k2, k1, rootProc)

            !=======
            case (2)
            !=======
                call InitializeFastJX53b (k1, k2, self%chem_mask_khi, self%num_qjs, &
     &                   self%cross_section_file, self%scattering_data_file, &
     &                   self%rate_file, self%T_O3_climatology_file)

            !=======
            case (3)
            !=======
                call InitializeFastJX53c (k1, k2, self%chem_mask_khi, self%num_qjs, &
     &                   self%cross_section_file, self%scattering_data_file, &
     &                   self%rate_file, self%T_O3_climatology_file)
            !=======
            case (4)
            !=======
                call InitializeFastJX65 (k1, k2, self%chem_mask_khi, self%num_qjs, &
     &                   self%cross_section_file, &
     &                   self%rate_file, self%T_O3_climatology_file, rootProc)
         !=========
         end select
         !=========
      end if

      if ((TRIM(chem_mecha) == 'troposphere') .or. &
          (TRIM(chem_mecha) == 'strat_trop')  .or. &
           TRIM(chem_mecha) == 'strat_trop_aerosol' ) then
         if ((self%phot_opt == 3) .and. self%do_AerDust_Calc) then
            call Allocate_dust  (self, i1,i2, ju1, j2, k1, k2, nSADdust)
            call Allocate_wAersl(self, i1,i2, ju1, j2, k1, k2, nSADaer)
            call Allocate_dAersl(self, i1,i2, ju1, j2, k1, k2)

            call Allocate_odAer   (self, i1,i2, ju1, j2, k1, k2, nSADaer, NRH_b)
            call Allocate_odMdust (self, i1,i2, ju1, j2, k1, k2, nSADdust)
            call Allocate_tArea   (self, i1,i2, ju1, j2, k1, k2, nSADdust, nSADaer)
            call Allocate_eRadius (self, i1,i2, ju1, j2, k1, k2, nSADdust, nSADaer)
            call Allocate_optDepth(self, i1,i2, ju1, j2, k1, k2, num_AerDust)
         end if
      end if

!... set up for AerDust output with gocart
      call Get_pr_AerDust(Diagnostics, pr_AerDust)
      if ( TRIM(chem_mecha) == 'gocart_aerosol' .and. pr_AerDust ) then
         call readRaaQaa(self%raa_qaa_data_file, RAA_b, QAA_b)

         call Allocate_dust  (self, i1,i2, ju1, j2, k1, k2, nSADdust)              
         call Allocate_wAersl(self, i1,i2, ju1, j2, k1, k2, nSADaer)               
         call Allocate_dAersl(self, i1,i2, ju1, j2, k1, k2)                        

         call Allocate_odAer   (self, i1,i2, ju1, j2, k1, k2, nSADaer, NRH_b)      
         call Allocate_odMdust (self, i1,i2, ju1, j2, k1, k2, nSADdust)            
         call Allocate_tArea   (self, i1,i2, ju1, j2, k1, k2, nSADdust, nSADaer)   
         call Allocate_eRadius (self, i1,i2, ju1, j2, k1, k2, nSADdust, nSADaer)   
         call Allocate_optDepth(self, i1,i2, ju1, j2, k1, k2, num_AerDust)         
      end if

      if (self%chem_opt == 3) then

         call Allocate_loss_freq(self, i1, i2, ju1, j2, k1, k2, numSpecies)
         if (self%loss_freq_opt == 2) call Allocate_lossData(self)

         call readLossFrequency (self%loss_freq, self%lossData,                &
     &            self%loss_data_infile_name, self%loss_init_val,              &
     &            self%loss_freq_opt, procID, pr_diag, i1, i2, ju1, j2, k1, k2,&
     &            numSpecies, KDIM, JDIM, MDIM)
      end if

      if (pr_const .and. pr_sulf_src) then
         call Allocate_dms_no3  (self, i1, i2, ju1, j2, k1, k2)
         call Allocate_dms_oh   (self, i1, i2, ju1, j2, k1, k2)
         call Allocate_so2_h2o2 (self, i1, i2, ju1, j2, k1, k2)
         call Allocate_so2_o3   (self, i1, i2, ju1, j2, k1, k2)
         call Allocate_so2_oh   (self, i1, i2, ju1, j2, k1, k2)
         self%dms_oh  (:,:,:) = 0.0d0
         self%dms_no3 (:,:,:) = 0.0d0

         self%so2_oh  (:,:,:) = 0.0d0
         self%so2_h2o2(:,:,:) = 0.0d0
         self%so2_o3  (:,:,:) = 0.0d0
      endif

      ! Boundary forcing Data
      !----------------------
      if (self%forc_bc_num > 0) then
         Allocate (self%forc_bc_data(FBC_LATDIM, FBC_MONDIM,  &
     &            self%forc_bc_years, self%forc_bc_num))
         self%forc_bc_data = 0.0d0
         call  readForcedBcData (pr_diag, procID, self%forc_bc_opt, &
     &             self%forc_bc_years, self%forc_bc_num, &
     &             self%forc_bc_init_val, self%forc_bc_infile_name, &
                   self%forc_bc_data)

         if ((self%forc_bc_opt == 1) .or. (self%forc_bc_opt == 2)) then
            allocate(dlatr (ju1_gl:j2_gl))
            call Get_dlatr (gmiDomain, dlatr )

            Allocate (self%jlatmd(j2_gl-ju1_gl+1))
            self%jlatmd = 0

            call Get_begGmiDate  (gmiClock, start_ymd     )
            call GmiSplitDateTime (start_ymd, year, month, idumday)

            self%lastYearFBC = year

            !-------------------------------------------------------------
            ! The forc_bc values are known only at certain locations.  The
            ! following loop is finding out which of these locations is
            ! closest to the simulation latitude zone.  No interpolation is
            ! done, just pick closest.
            !-------------------------------------------------------------

            do ij = ju1_gl, j2_gl
               rjx = (dlatr(ij) * DEGPRAD / 10.d0) + 10.5d0

              !-----------------------------------------------------
              ! Add a small delta to take care of precision problems.
              !-----------------------------------------------------
              rjx = rjx + 1.0d-06

              jx = rjx
              jj = ij + 1 - ju1_gl

              self%jlatmd(jj) = Min (FBC_LATDIM, Max (1, jx))
            end do
            deallocate(dlatr)
         end if

      end if

      if ((self%sfalbedo_opt == 1) .or. (self%sfalbedo_opt == 2)) then
         allocate(self%sasdir_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         self%sasdir_data = 0.0d0

         allocate(self%sasdif_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         self%sasdif_data = 0.0d0

         allocate(self%saldir_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         self%saldir_data = 0.0d0

         allocate(self%saldif_data(i1:i2, ju1:j2, MONTHS_PER_YEAR))
         self%saldif_data = 0.0d0

         if (self%sfalbedo_opt  == 1) then
            self%sasdir_data(:,:,:) = self%sasdir_init_val
            self%sasdif_data(:,:,:) = self%sasdif_init_val
            self%saldir_data(:,:,:) = self%saldir_init_val
            self%saldif_data(:,:,:) = self%saldif_init_val
         end if
      end if

      call Get_itloop         (gmiGrid, itloop)
      call Get_pr_smv2        (Diagnostics, pr_smv2)
      call Get_do_ShipEmission(Emission, do_ShipEmission)
      call Get_gmiTimeStep    (gmiClock, tdt           )
      call Get_qj_var_name (Diagnostics, qj_var_name)

      self%chemintv = -1.0d0
      self%rateintv = -1.0d0

      self%chemintv = tdt * self%chem_cycle

      if (self%chem_cycle < 1.0d0) then  ! subcycle chemistry
         self%rateintv = tdt
      else if (self%chem_cycle == 1.0d0) then
         self%rateintv = self%chemintv
      end if

      ! find reaction number of NO photolysis
      self%jno_num = 999
      self%jnoxnum = 999

      do ic = 1, self%num_qjs
         if (TRIM(chem_mecha) == 'strat_trop' .or. &
             TRIM(chem_mecha) == 'strat_trop_aerosol' ) then
            if (Trim(self%qj_labels(ic)) =='NO + hv = N + O') self%jno_num = ic
         endif
      end do

      if (do_ShipEmission) then
         do ic=1, self%num_qjs
            if (TRIM(chem_mecha) == 'strat_trop' .or. &
                TRIM(chem_mecha) == 'strat_trop_aerosol' ) then
               if (Trim(self%qj_labels(ic)) =='NO2 + hv = NO + O') self%jnoxnum = ic
            else if (TRIM(chem_mecha) == 'troposphere') then
               if (Trim(self%qj_labels(ic)) =='NO2 + hv = NO + O3') self%jnoxnum = ic
            endif
         end do
         if (self%jnoxnum.eq.999) then
            err_msg = 'jnoxnum not found in updateChemistry'
            call GmiPrintError(err_msg, .true., 1, self%jnoxnum, 0, 0, 0.0d0, 0.0d0)
         endif

         self%prevEmissionFileRecord = -1
         allocate(self%emiss_ozone (i1:i2,ju1:j2))
      endif

      ! Initialize the Solver
      !----------------------

      IF (self%chem_opt == 2) THEN
         call Get_ilong         (gmiGrid, ilong)
         call Get_ivert         (gmiGrid, ivert)
         call Get_ilat          (gmiGrid, ilat )
        if (pr_smv2) then
           smv_filnam = 'smv2chem.asc'
!          call makeOutfileName (smv_filnam, '.smv2.asc', gmi_problem_name)
        end if

        pr_nc_period = -1  ! need to be properly initialized
        call Do_Smv2_Init (self%smv2SavedVars, smv_filnam, self%do_smv_reord, &
                           pr_diag, pr_smv2, procID, self%chemintv, &
                           tdt, i1, i2, ju1, j2, k1, k2, ilong, ilat, ivert, &
                           itloop, self%num_qjs, self%num_qks, self%num_active, pr_nc_period, &
                           pr_qqjk, do_qqjk_inchem)

        if ( (self%phot_opt == 7) .and. (INDMS .ne. 0) ) then
           call Do_Sulf_Init (self%sulfSavedVars, pr_diag, procID, j2_gl, &
                              self%qj_infile_name, qj_var_name)
        end if

      ENDIF

      return

      end subroutine InitializeChemistry
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initReadChemistry
!
! !INTERFACE:
!
  subroutine initReadChemistry (self, gmiGrid, gmiDomain, Diagnostics, &
                          chem_mecha)
!
! !USES:
  use ReadPhotolysisRates_mod   , only : readPhotolysisRates, readSolarCycle,  &
          readPhotolysisTable
  use ReadLiqBinarySulfate_mod  , only : readLiqBinarySulfate
  use ReadWaterClimatology_mod  , only : readWaterClimatology, readH2ocond
  use ReadReactionRateAdjFac_mod, only : readReactionRateAdjustmentFactors
  use ReadForcedBC_mod          , only : readForcedBcData
!
#     include "phot_lookup.h"
#     include "phot_monthly.h"
!
! !INPUT PARAMETERS:
      type (t_gmiGrid    ), intent(in) :: gmiGrid
      type (t_gmiDomain  ), intent(in) :: gmiDomain
      type (t_Diagnostics), intent(in) :: Diagnostics
      character (len=*)   , intent(in) :: chem_mecha
!
! !INPUT/OUTPUT PARAMETERS:
  type  (t_Chemistry), intent(inOut) :: self
!
! !DESCRIPTION:
! This routine reads in chemistry related files that contain geographical
! related data. The read is done once during initialization procedures
! for the Chemistry component.
!
! !LOCAL VARIABLES:
      integer :: ncid_rra, cnt1d(1), strt1d(1)
      integer :: i1, i2, i1_gl, i2_gl, ju1, j2, ju1_gl, j2_gl
      integer :: k1, k2
      integer :: procID, restart_inrec
      logical :: pr_diag, rd_restart
      character(len=MAX_LENGTH_VAR_NAME) :: qj_var_name
      character (len=MAX_LENGTH_FILE_NAME) :: restart_infile_name
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      call Get_procID (gmiDomain, procID)
      call Get_pr_diag (Diagnostics, pr_diag)

      if (pr_diag) Write(6,*) "initReadChemistry called by ", procID

      call Get_qj_var_name (Diagnostics, qj_var_name)

      call Get_i1    (gmiGrid, i1)
      call Get_i2    (gmiGrid, i2)
      call Get_ju1   (gmiGrid, ju1)
      call Get_j2    (gmiGrid, j2)
      call Get_k1    (gmiGrid, k1)
      call Get_k2    (gmiGrid, k2)
      call Get_i1_gl (gmiGrid, i1_gl)
      call Get_i2_gl (gmiGrid, i2_gl)
      call Get_ju1_gl(gmiGrid, ju1_gl)
      call Get_j2_gl (gmiGrid, j2_gl )

      if (self%do_rxnr_adjust) then

         ! Hard-code for now.
         self%rxnr_adjust_timpyr = MONTHS_PER_YEAR

         call Ncop_Rd (ncid_rra, self%rxnr_adjust_infile_name)
         call Ncget_Dimlen (ncid_rra, 'rra_dim', self%num_rxnr_adjust)
         strt1d(1) = 1
         cnt1d (1) = self%num_rxnr_adjust

         call Allocate_rxnr_adjust_map(self)

         call Ncrd_1d_Int (self%rxnr_adjust_map, ncid_rra, 'rra_dim', strt1d, cnt1d)
         call Nccl (ncid_rra)

         call readReactionRateAdjustmentFactors  &
                 (self%rxnr_adjust, self%rxnr_adjust_infile_name, &
                  self%rxnr_adjust_var_name, self%num_rxnr_adjust, &
                  self%rxnr_adjust_timpyr, pr_diag, procID, &
                  i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl)
      else
         self%num_rxnr_adjust = 0
      end if

      if (self%phot_opt == 1) then
         self%qjgmi(:,:,:,:) = self%qj_init_val
      else if (self%phot_opt == 2) then
         Allocate (qjmon (i1:i2, ju1:j2, k1:k2, self%num_qjs, self%qj_timpyr))
         qjmon    = 0.0d0

         call  readPhotolysisRates(i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, &
                                self%num_qjs, self%qj_timpyr, &
                                qjmon, self%qj_infile_name, qj_var_name)
      else if ((self%phot_opt == 4) .or. (self%phot_opt == 5)) then
         call readPhotolysisTable(self%phot_opt, self%qj_infile_name, TRIM(chem_mecha), &
     &               self%num_qjs, i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl) 

         ! solar cycle relative strength data
         call readSolarCycle(self%do_solar_cycle, self%sc_infile_name)
      end if

      if (self%sad_opt /= 0) then
         call Allocate_sadgmi(self, i1, i2, ju1, j2, k1, k2)
         if ((self%sad_opt == 1) .or. (self%sad_opt == 2)) then
    
            call Allocate_hno3gas (self, i1, i2, ju1, j2, k1, k2)
            call Allocate_hno3cond(self, i1, i2, ju1, j2, k1, k2)

            if (self%sad_opt == 2) then
               call Allocate_h2oback(self, i1, i2, ju1, j2, k1, k2)
               call Allocate_h2ocond(self, i1, i2, ju1, j2, k1, k2)

               call Get_rd_restart(Diagnostics, rd_restart)

               if (rd_restart) then
                  call Get_restart_inrec(Diagnostics, restart_inrec)
                  call Get_restart_infile_name(Diagnostics, restart_infile_name)

                  call readH2ocond (self%h2ocond, restart_inrec,               &
     &                     restart_infile_name, i1, i2, ju1, j2, k1, k2,       &
     &                     i1_gl, i2_gl, ju1_gl, j2_gl)
               end if

               call Allocate_reffice(self, i1, i2, ju1, j2, k1, k2)
               call Allocate_reffsts(self, i1, i2, ju1, j2, k1, k2)
               call Allocate_vfall  (self, i1, i2, ju1, j2, k1, k2)

               call Allocate_h2oclim(self, i1, i2, ju1, j2, k1, k2)
               call Allocate_ch4clim(self, i1, i2, ju1, j2, k1, k2)
               call Allocate_lbssad (self, i1, i2, ju1, j2, k1, k2)

               call ReadWaterClimatology (self%h2oclim, self%ch4clim, &
     &                  self%h2oclim_init_val, self%ch4clim_init_val,          &
     &                  self%h2oclim_opt,  self%h2oclim_infile_name, procID,   &
     &                  pr_diag, i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl,       &
     &                  self%h2oclim_timpyr)

               call ReadLiqBinarySulfate (self%lbssad, self%lbssad_init_val,   &
     &                  self%lbssad_opt, self%lbssad_infile_name, procID,      &
     &                  pr_diag, i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl,       &
     &                  self%lbssad_timpyr)
            end if
         elseif (self%sad_opt == 3) then
            call Allocate_lbssad (self, i1, i2, ju1, j2, k1, k2)

            call ReadLiqBinarySulfate (self%lbssad, self%lbssad_init_val,      &
     &               self%lbssad_opt, self%lbssad_infile_name, procID, pr_diag,&
     &               i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, self%lbssad_timpyr)
         end if
      end if

      return

      end subroutine initReadChemistry
!
!EOC
!-------------------------------------------------------------------------
!
!BOP
!
! !IROUTINE: runReadChemistry
!
! !INTERFACE:
!
      subroutine runReadChemistry (self, Emission, gmiClock, gmiGrid, &
     &              gmiDomain, Diagnostics, chem_mecha)
!
! !USES:
      use ReadAerosolDust_mod, only : ReadAerosolDust
!
! !INPUT PARAMETERS:
      character (len=*)   , intent(in) :: chem_mecha
      type (t_gmiGrid    ), intent(in) :: gmiGrid
      type (t_Emission   ), intent(in) :: Emission
      type (t_GmiClock   ), intent(in) :: gmiClock
      type (t_gmiDomain  ), intent(in) :: gmiDomain
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Chemistry), intent(inOut) :: self
!
! !DESCRIPTION:
! This routine reads in chemistry related files that contain geographical
! related data.
! Read in a new record at the beginning of each day (daily file)
! or the beginning of each month (monthly file)
!
! !LOCAL VARIABLES:
      integer       :: i1, i2, i1_gl, ju1, j2, ju1_gl, k1, k2
      integer       :: ydummy, thisDay, thisMonth, thisDate, ddummy
      integer, save :: curDate   = -1
      integer, save :: dailyRec  = -1
      integer, save :: thisRecord  = -1
      logical       :: doReadDailyEmiss, pr_diag
      integer       :: begDailyEmissRec, endDailyEmissRec
      integer       :: nymd, emiss_timpyr, procID
      integer       :: curEmissionFileRecord
!
!EOP
!-------------------------------------------------------------------------
!BOC
      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write(6,*)'runReadChemistry called by ', procID
!
      call Get_i1    (gmiGrid, i1)
      call Get_i2    (gmiGrid, i2)
      call Get_ju1   (gmiGrid, ju1)
      call Get_j2    (gmiGrid, j2)
      call Get_k1    (gmiGrid, k1)
      call Get_k2    (gmiGrid, k2)
      call Get_i1_gl (gmiGrid, i1_gl)
      call Get_ju1_gl(gmiGrid, ju1_gl)

      ! Obtain model time information

      call Get_curGmiDate  (gmiClock, nymd          )

       if ((TRIM(chem_mecha) == 'troposphere').or.(TRIM(chem_mecha) == 'strat_trop')  &
     &         .or. TRIM(chem_mecha) == 'strat_trop_aerosol' ) then
          if ((self%phot_opt == 3) .and. self%do_AerDust_Calc) then

             call Get_curEmissionFileRecord(Emission, curEmissionFileRecord)
             if (thisRecord /= curEmissionFileRecord) then
                thisRecord = curEmissionFileRecord
                call ReadAerosolDust (self%dust, self%wAersl, self%dAersl, &
     &                                self%AerDust_infile_name, thisRecord, &
     &                                i1,i2, ju1, j2, k1, k2, i1_gl, ju1_gl)
            end if
         end if
      end if

      return

      end subroutine runReadChemistry
!
!EOC
!-------------------------------------------------------------------------
!
!BOP
!
! !IROUTINE: runCalcAerosolDust
!
! !INTERFACE:
!
      subroutine runCalcAerosolDust (self,  &
     &              mass, mcor, gridBoxHeight, SpeciesConcentration,  &
     &              gmiGrid, gmiDomain, Diagnostics, chem_mecha)
! !USES:
!    Calculate the fastj aerosol parameters from the gocart aerosol fields
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration, Get_concentration
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      implicit none
!
!
! !INPUT PARAMETERS:
      type(t_SpeciesConcentration), intent(in) :: SpeciesConcentration
      type(t_GmiArrayBundle), pointer :: concentration(:)
      type (t_gmiGrid    ), intent(in) :: gmiGrid
      type (t_gmiDomain  ), intent(in) :: gmiDomain
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Chemistry), intent(inOut) :: self
!
! !INPUT PARAMETERS:
      real*8 , intent(in ) :: mcor (:,:)
      real*8 , intent(in ) :: mass (:,:,:)
      real*8 , intent(in ) :: gridBoxHeight (:,:,:)
      character (len=*)  , intent(in) :: chem_mecha
!
! !DESCRIPTION:
! This routine reads in chemistry related files that contain geographical
! related data.
! Read in a new record at the beginning of each day (daily file)
! or the beginning of each month (monthly file)
!
! !LOCAL VARIABLES:
      integer       :: i1, i2, ju1, j2, k1, k2
      integer       :: procID
      logical       :: pr_diag
!
!EOP
!-------------------------------------------------------------------------
!BOC
      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write(6,*)'runCalcAerosolDust called by ', procID
!
      call Get_i1    (gmiGrid, i1)
      call Get_i2    (gmiGrid, i2)
      call Get_ju1   (gmiGrid, ju1)
      call Get_j2    (gmiGrid, j2)
      call Get_k1    (gmiGrid, k1)
      call Get_k2    (gmiGrid, k2)

      call Get_concentration(SpeciesConcentration, concentration)

!... calculate the aerosol parameters from the gocart fields to be used in Fastj, etc.
      call CalcAerosolDust (self%dust, self%wAersl, self%dAersl, &
     &     mass, mcor, gridBoxHeight, concentration,  &
     &     gmiDomain, Diagnostics, TRIM(chem_mecha), i1, i2, ju1, j2, k1, k2)

      return

      end subroutine runCalcAerosolDust
!
!EOC
!-------------------------------------------------------------------------
!
!BOP
!
! !IROUTINE: CalcAerosolDust
!
! !INTERFACE:
!
      subroutine CalcAerosolDust (dust, wAersl, dAersl, mass, mcor, gridBoxHeight, &
     &              concentration, gmiDomain, Diagnostics, chem_mecha, i1, i2, ju1, j2, k1, k2)
!
! !USES:
!    Calculate the fastj aerosol parameters from the gocart aerosol fields
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      
      implicit none

#     include "setkin_par.h"
!
!
!
! !INPUT PARAMETERS:
      type(t_GmiArrayBundle), pointer :: concentration(:)
!      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)
      type(t_gmiDomain  ), intent(in) :: gmiDomain
      type(t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT PARAMETERS:
      real*8 , intent(in ) :: mcor (i1:i2,ju1:j2)
      real*8 , intent(in ) :: mass (i1:i2,ju1:j2,k1:k2)
      real*8 , intent(in ) :: gridBoxHeight (i1:i2,ju1:j2,k1:k2)
      character (len=*)   , intent(in) :: chem_mecha
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: dust  (i1:i2, ju1:j2, k1:k2, nSADdust)
      real*8 , intent(out) :: wAersl(i1:i2, ju1:j2, k1:k2, nSADaer )
      real*8 , intent(out) :: dAersl(i1:i2, ju1:j2, k1:k2, 2       )
!
! !DESCRIPTION:
! This routine calculstes the aerosol/dust parameters needed for the calculation
!     of the Aerosol Optical Thickness needed for Fastj from the aerosol fields 
!     carried in the mechanism
!
! !LOCAL VARIABLES:
      integer      :: procID
      logical      :: pr_diag
      integer      :: i1, i2, ju1, j2, k1, k2
      integer      :: im, ik
      real*8, allocatable :: conv(:,:,:)
!.sds... set up weighting to split GOCART dust bin 1 into fields for Fastj
      real*8, save :: Dust1Wgt(4)
      data Dust1Wgt / 0.01053d0, 0.08421d0, 0.25263d0, 0.65263d0 /

!.sds   temporary fix for [f,n]SO4n[1-3] being in mass mixing ratio
      real, save :: mwso4, mwm
      data mwso4 / 98.074 /
      data mwm   / 28.960 /
      real*8, save :: vmr2mmr

!.sds.. gocart already in mass mixing ratio
      if(TRIM(chem_mecha) ==  'gocart_aerosol') then
        vmr2mmr = 1.0
      else
        vmr2mmr = ( mwso4/mwm )
      endif
!
!EOP
!-------------------------------------------------------------------------
!BOC
      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write(6,*)'CalcAerosolDust called by ', procID

#ifndef nonZeroInd_tracers
!.sds... There are 7 dust bins, but GMI-GOCARTsimulates only 5 bins. 
!.sds...  separate the first bin into 4 sub-bins withfactors: 0.01053, 0.08421, 0.25263, 0.65263. 
!.sds...  Then use GOCART bins 2, 3 and 4+5 to fill the rest of the 7 bins

      dust(i1:i2, ju1:j2, k1:k2, 1) = Dust1Wgt(1) *  &
     &    concentration(IDUST1)%pArray3D(i1:i2, ju1:j2, k1:k2)

      dust(i1:i2, ju1:j2, k1:k2, 2) = Dust1Wgt(2) *  &
     &    concentration(IDUST1)%pArray3D(i1:i2, ju1:j2, k1:k2)

      dust(i1:i2, ju1:j2, k1:k2, 3) = Dust1Wgt(3) *  &
     &    concentration(IDUST1)%pArray3D(i1:i2, ju1:j2, k1:k2)

      dust(i1:i2, ju1:j2, k1:k2, 4) = Dust1Wgt(4) *  &
     &    concentration(IDUST1)%pArray3D(i1:i2, ju1:j2, k1:k2)

      dust(i1:i2, ju1:j2, k1:k2, 5) =   &
     &    concentration(IDUST2)%pArray3D(i1:i2, ju1:j2, k1:k2)

      dust(i1:i2, ju1:j2, k1:k2, 6) =   &
     &    concentration(IDUST3)%pArray3D(i1:i2, ju1:j2, k1:k2)

      dust(i1:i2, ju1:j2, k1:k2, 7) =   &
     &    concentration(IDUST4)%pArray3D(i1:i2, ju1:j2, k1:k2) &
     &  + concentration(IDUST5)%pArray3D(i1:i2, ju1:j2, k1:k2)

!.sds... set up conversion factor to go from mass mixing ratio to kg/m^3
      allocate(conv(i1:i2, ju1:j2, k1:k2))
      do ik = k1,k2
        conv(i1:i2, ju1:j2, ik) = mass(i1:i2, ju1:j2, ik) /  &
     &     ( mcor(i1:i2, ju1:j2) * gridBoxHeight(i1:i2, ju1:j2, ik) )
      enddo

!...sds aerosol species are carried in the model as mass mixing ratio, convert to kg/m^3
      do im = 1,nSADdust
         dust(i1:i2, ju1:j2, k1:k2, im) = &
     &      dust(i1:i2, ju1:j2, k1:k2, im) * conv(i1:i2, ju1:j2, k1:k2)
      enddo

!.sds... There are 2 options, we are currently only set up for option 1.
!.sds... 
!.sds... 1. Put all GMI BC into hydrophobic mode (ie. all BC into dAersl) and all 
!.sds... GMI OC into hydrophilic mode (ie. all OC into wAersl).
!.sds... 2. Check the GOCART output results to find a global flat ratio of 
!.sds... BC(hydrophilic) / BC(hydrophobic) and use this ratio to scale GMI BC. 
!.sds... Same for OC.
!.sds...

! Sulfate SO4 - add all GOCART SO4 species together
      wAersl(i1:i2, ju1:j2, k1:k2,1) =   &
     &       concentration(IFSO4A )%pArray3D(i1:i2, ju1:j2, k1:k2)  &
     &     + concentration(INSO4A )%pArray3D(i1:i2, ju1:j2, k1:k2)  &
     &     +(concentration(INSO4N1)%pArray3D(i1:i2, ju1:j2, k1:k2)  &
     &     + concentration(INSO4N2)%pArray3D(i1:i2, ju1:j2, k1:k2)  &
     &     + concentration(INSO4N3)%pArray3D(i1:i2, ju1:j2, k1:k2)  &
     &     + concentration(IFSO4N1)%pArray3D(i1:i2, ju1:j2, k1:k2)  &
     &     + concentration(IFSO4N2)%pArray3D(i1:i2, ju1:j2, k1:k2)  &
     &     + concentration(IFSO4N3)%pArray3D(i1:i2, ju1:j2, k1:k2) ) * vmr2mmr

! Hydrophilic BC
      wAersl(i1:i2, ju1:j2, k1:k2,2) = 0.80d0 * (  &
     &      concentration(IBBC)%pArray3D(i1:i2, ju1:j2, k1:k2)  &
     &    + concentration(IFBC)%pArray3D(i1:i2, ju1:j2, k1:k2) )
         
! Hydrophilic OC
      wAersl(i1:i2, ju1:j2, k1:k2,3) =  &
     &      concentration(IFOC)%pArray3D(i1:i2, ju1:j2, k1:k2)  &
     &    + concentration(IBOC)%pArray3D(i1:i2, ju1:j2, k1:k2)  &
     &    + concentration(INOC)%pArray3D(i1:i2, ju1:j2, k1:k2)

! Sea Salt (accum)
      wAersl(i1:i2, ju1:j2, k1:k2,4) =  &
     &      concentration(ISSLT1)%pArray3D(i1:i2, ju1:j2, k1:k2)

! Sea Salt (coarse)

      wAersl(i1:i2, ju1:j2, k1:k2,5) =  &
     &      concentration(ISSLT2)%pArray3D(i1:i2, ju1:j2, k1:k2)
      wAersl(i1:i2, ju1:j2, k1:k2,6) =  &
     &      concentration(ISSLT3)%pArray3D(i1:i2, ju1:j2, k1:k2)
      wAersl(i1:i2, ju1:j2, k1:k2,7) =  &
     &      concentration(ISSLT4)%pArray3D(i1:i2, ju1:j2, k1:k2)

!...sds aerosol species are carried in the model as mass mixing ratio, convert to kg/m^3      
      do im = 1,nSADaer
         wAersl(i1:i2, ju1:j2, k1:k2, im) = &
     &      wAersl(i1:i2, ju1:j2, k1:k2, im) * conv(i1:i2, ju1:j2, k1:k2)
      enddo

! Hydrophobic BC
      dAersl(i1:i2, ju1:j2, k1:k2,1) = 0.20d0 * (  &
     &      concentration(IBBC)%pArray3D(i1:i2, ju1:j2, k1:k2)  &
     &    + concentration(IFBC)%pArray3D(i1:i2, ju1:j2, k1:k2) )

! Hydrophobic OC
      dAersl(i1:i2, ju1:j2, k1:k2,2) = 0.0d0

!...sds aerosol species are carried in the model as mass mixing ratio, convert to kg/m^3
      do im = 1,2
         dAersl(i1:i2, ju1:j2, k1:k2, im) = &
     &      dAersl(i1:i2, ju1:j2, k1:k2, im) * conv(i1:i2, ju1:j2, k1:k2)
      enddo

 
      deallocate(conv)


#endif
      return

      end subroutine CalcAerosolDust
!
!EOC
!-------------------------------------------------------------------------
!BOP
! 
! !IROUTINE: RunChemistry
!
! !INTERFACE:
!
       subroutine RunChemistry (self, Emission, SpeciesConcentration, gmiClock, &
     &               gmiGrid, gmiDomain, Diagnostics, metFields,                &
     &               chem_mecha, averagePressEdge, num_AerDust)
!
! !USES:
      use GmiUpdateChemistry_mod    , only : updateChemistry
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle, CleanArrayPointer
      use GmiEmissionMethod_mod, only : Get_emissionArray
      use GmiEmissionMethod_mod, only : Get_emiss_map       , Get_emissDust 
      use GmiEmissionMethod_mod, only : Get_emiss_opt 
      use GmiEmissionMethod_mod, only : Get_emiss_dust_opt  , Get_ndust 
      use GmiEmissionMethod_mod, only : Get_do_ShipEmission , Get_num_emiss
      use GmiEmissionMethod_mod, only : Get_do_semiss_inchem, Get_emiss_map_dust
      use GmiEmissionMethod_mod, only : Get_doReadDailyEmiss, Get_o3_index
      use GmiEmissionMethod_mod, only : Get_begDailyEmissRec, Get_endDailyEmissRec
      use GmiEmissionMethod_mod, only : Get_emiss_nox, Get_emiss_isop, Get_emiss_monot
      use GmiEmissionMethod_mod, only : Get_emiss_o3 , Get_emiss_hno3, Set_emiss_o3
      use GmiEmissionMethod_mod, only : Get_surf_emiss_out2, Set_surf_emiss_out2
      use GmiEmissionMethod_mod, only : Get_surf_emiss_out , Set_surf_emiss_out
      use GmiEmissionMethod_mod, only : Get_emiss_3d_out   , Set_emiss_3d_out  
      use ReadAerosolDust_mod  , only : ReadAerosolDust
      use GmiAerDustODSA_mod      , only : Aero_OptDep_SurfArea
      use GmiAerDustODSA_mod      , only : Dust_OptDep_SurfArea
!
      implicit none

!.sds not needed...#     include "gmi_sad_constants.h"
#     include "gmi_diag_constants_llnl.h"

!
! !INPUT PARAMETERS:
      integer            , intent(in) :: num_AerDust
      real*8             , intent(in) :: averagePressEdge (:)
      character (len=*)  , intent(in) :: chem_mecha
      type(t_GmiClock   ), intent(in) :: gmiClock
      type(t_gmiGrid    ), intent(in) :: gmiGrid 
      type(t_gmiDomain  ), intent(in) :: gmiDomain
      type(t_metFields  ), intent(in) :: metFields
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_Diagnostics)         , intent(inOut) :: Diagnostics
      type(t_Emission)            , intent(inOut) :: Emission
      type(t_Chemistry)           , intent(inOut) :: self
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !LOCAL VARIABLES:
      real*8 , allocatable :: latdeg(:), londeg(:), dlatr(:), mcor(:,:)
      logical :: do_semiss_inchem, do_ShipEmission
      integer :: num_emiss, emiss_timpyr, ndust
      integer :: emiss_opt, emiss_dust_opt
      integer, allocatable :: emiss_map     (:)
      integer, allocatable :: emiss_map_dust(:)
      real*8 , allocatable :: emiss_dust(:,:,:)
      real*8 , allocatable :: emiss_o3   (:,:)
      real*8 , allocatable :: emiss_nox  (:,:)
      real*8 , allocatable :: emiss_hno3 (:,:)
      real*8 , allocatable :: emiss_isop (:,:)
      real*8 , allocatable :: emiss_monot(:,:)
      type (t_GmiArrayBundle), pointer :: emissionArray(:)
      integer :: ic, o3_index

      real*8 , allocatable :: const_init_val(:)

      real*8, allocatable :: surf_emiss_out (:,:,:), surf_emiss_out2(:,:,:)
      real*8, allocatable :: emiss_3d_out(:,:,:,:)

      integer       :: nymd, nhms, num_time_steps, start_ymd
      real*8        :: tdt, gmi_sec

      integer       :: numSpecies, procID
      integer       :: commuWorld, curEmissionFileRecord
      integer       :: ilo, ihi, julo, jhi
      integer       :: ilong, ilat, ivert, itloop, met_opt
      integer       :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer       :: i1, i2, ju1, j2, k1, k2, numDomains
      type (t_GmiArrayBundle), pointer :: concentration(:)

      character(len=MAX_LENGTH_VAR_NAME) :: qj_var_name
      real*8  :: constOutputFrequency
      logical :: pr_diag, pr_qqjk, pr_smv2, pr_ascii5, pr_AerDust
      logical :: pr_qj_opt_depth, pr_qj_o3_o1d
      logical :: do_qqjk_inchem, do_ftiming
      logical :: pr_const, pr_sulf_src, rd_restart
      logical :: pr_emiss_3d, pr_surf_emiss
      logical :: do_aerocom, do_qqjk_reset

      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model
      real*8 , allocatable :: relativeHumidity (:,:,:)
      real*8 , allocatable :: mass(:,:,:)
      real*8 , allocatable :: press3c(:,:,:), press3e(:,:,:)
      real*8 , allocatable :: tropopausePress(:,:)
      real*8 , allocatable :: gridBoxHeight(:,:,:)

      real*8 , allocatable :: humidity(:,:,:), pctm2(:,:), moistq(:,:,:)
      real*8 , allocatable :: totalCloudFraction(:,:,:), fracCloudCover(:,:)
      real*8 , allocatable :: tau_cloud(:,:,:), kel(:,:,:)
      real*8 , allocatable :: taucli(:,:,:), tauclw(:,:,:)
      integer, allocatable :: lwi_flags(:,:)
      real*8 , allocatable :: surf_air_temp(:,:), surf_alb_uv(:,:)
      real*8 , allocatable :: cmf   (:,:,:), clwc(:,:,:)
      real*8 , allocatable :: radswg (:,:)
      real*8 , allocatable :: prod(:,:,:)
      real*8 , allocatable :: loss(:,:,:)
      logical :: rootProc
!EOP
!-------------------------------------------------------------------------
!EOC
      call Get_procID    (gmiDomain, procID    )
      call Get_pr_diag        (Diagnostics, pr_diag)

      if (pr_diag) Write(6,*) 'runChemistry called by ', procID

      call Get_numSpecies(gmiGrid, numSpecies)

      call Get_pr_ascii5      (Diagnostics, pr_ascii5)
      call Get_pr_qqjk        (Diagnostics, pr_qqjk)
      call Get_pr_smv2        (Diagnostics, pr_smv2)
      call Get_pr_const       (Diagnostics, pr_const)
      call Get_rd_restart     (Diagnostics, rd_restart)
      call Get_do_aerocom     (Diagnostics, do_aerocom)
      call Get_do_ftiming     (Diagnostics, do_ftiming)
      call Get_qj_var_name    (Diagnostics, qj_var_name)
      call Get_pr_emiss_3d    (Diagnostics, pr_emiss_3d)
      call Get_pr_sulf_src    (Diagnostics, pr_sulf_src)
      call Get_constOutputFrequency   (Diagnostics, constOutputFrequency)
      call Get_do_qqjk_reset   (Diagnostics, do_qqjk_reset)
      call Get_pr_qj_o3_o1d   (Diagnostics, pr_qj_o3_o1d)
      call Get_pr_surf_emiss  (Diagnostics, pr_surf_emiss)
      call Get_do_qqjk_inchem (Diagnostics, do_qqjk_inchem)
      call Get_pr_qj_opt_depth(Diagnostics, pr_qj_opt_depth)

      ! Get the communicator for slave processors
      call Get_communicatorWorld(gmiDomain, commuWorld)

      ! Get the GMI grid information
      call Get_i1    (gmiGrid, i1   )
      call Get_i2    (gmiGrid, i2   )
      call Get_ju1   (gmiGrid, ju1  )
      call Get_j2    (gmiGrid, j2   )
      call Get_k1    (gmiGrid, k1   )
      call Get_k2    (gmiGrid, k2   )
      call Get_i1_gl (gmiGrid, i1_gl)
      call Get_i2_gl (gmiGrid, i2_gl)
      call Get_ju1_gl(gmiGrid, ju1_gl)
      call Get_j2_gl (gmiGrid, j2_gl )
      call Get_ilo   (gmiGrid, ilo  )
      call Get_ihi   (gmiGrid, ihi  )
      call Get_julo  (gmiGrid, julo )
      call Get_jhi   (gmiGrid, jhi  )
      call Get_ilong (gmiGrid, ilong)
      call Get_ilat  (gmiGrid, ilat )
      call Get_ivert (gmiGrid, ivert)
      call Get_itloop(gmiGrid, itloop)
      call Get_numSpecies(gmiGrid, numSpecies)

      call Get_numDomains(gmiDomain, numDomains)

      allocate(latdeg(ju1_gl:j2_gl))
      allocate(londeg( i1_gl:i2_gl))
      allocate(dlatr (ju1_gl:j2_gl))
      allocate(mcor    (i1:i2,ju1:j2))

      call Get_mcor  (gmiDomain, mcor  )
      call Get_dlatr (gmiDomain, dlatr )
      call Get_latdeg(gmiDomain, latdeg)
      call Get_londeg(gmiDomain, londeg)

      ! Obtain model time information

      call Get_begGmiDate  (gmiClock, start_ymd     )
      call Get_curGmiDate  (gmiClock, nymd          )
      call Get_curGmiTime  (gmiClock, nhms          )
      call Get_numTimeSteps(gmiClock, num_time_steps)
      call Get_gmiSeconds  (gmiClock, gmi_sec       )
      call Get_gmiTimeStep (gmiClock, tdt           )

      call Get_concentration(SpeciesConcentration, concentration)

      allocate(const_init_val(1:numSpecies))
      call Get_const_init_val(SpeciesConcentration, const_init_val, numSpecies)

      call Get_met_opt(metFields, met_opt)
      call Get_metdata_name_org  (metFields, metdata_name_org)
      call Get_metdata_name_model(metFields, metdata_name_model)

      allocate(kel(ilo:ihi, julo:jhi, k1:k2))
      allocate(humidity(i1:i2, ju1:j2, k1:k2))
      allocate(pctm2(ilo:ihi, julo:jhi))

      call Get_kel(metFields, kel)
      call Get_pctm2(metFields, pctm2)
      call Get_humidity(metFields, humidity)

      allocate(lwi_flags(i1:i2, ju1:j2))
      call Get_lwi_flags(metFields, lwi_flags)

      if (met_opt /= 1) then
         allocate(surf_alb_uv(i1:i2, ju1:j2))
         call Get_surf_alb_uv(metFields, surf_alb_uv)
      end if

      if (met_opt == 3) then
         allocate(moistq(i1:i2, ju1:j2, k1:k2))
         call Get_moistq(metFields, moistq)

         if ((self%chem_opt ==2) .or. (self%chem_opt ==7) .or. (self%chem_opt == 8)) then
!.sds         if ((self%chem_opt ==7) .or. (self%chem_opt == 8)) then
            allocate(totalCloudFraction(i1:i2, ju1:j2, k1:k2))
            call Get_totalCloudFraction(metFields, totalCloudFraction)
         end if

         if (metdata_name_org(1:4) == 'GMAO') then
            allocate(taucli(i1:i2, ju1:j2, k1:k2))
            call Get_taucli(metFields, taucli)
            allocate(tauclw(i1:i2, ju1:j2, k1:k2))
            call Get_tauclw(metFields, tauclw)
         end if

         allocate(fracCloudCover(i1:i2, ju1:j2))
         call Get_fracCloudCover(metFields, fracCloudCover)

         allocate(tau_cloud(i1:i2, ju1:j2, k1:k2))
         call Get_tau_cloud(metFields, tau_cloud)

         allocate(surf_air_temp(i1:i2, ju1:j2))
         call Get_surf_air_temp(metFields, surf_air_temp)

         allocate(cmf   (i1:i2, ju1:j2, k1:k2))
         call Get_cmf(metFields, cmf)

         if ((metdata_name_org(1:4) == 'GISS') .or.  &
     &      ((metdata_name_org  (1:4) == 'NCAR') .and.  &
     &       (metdata_name_model(1:4) == 'CCM3'))) then
            allocate(clwc(i1:i2, ju1:j2, k1:k2))
            call Get_clwc(metFields, clwc)
         end if

         allocate(radswg (i1:i2, ju1:j2))
         call Get_radswg(metFields, radswg)
      end if

      allocate(press3c(ilo:ihi, julo:jhi, k1:k2))
      allocate(press3e(ilo:ihi, julo:jhi, k1-1:k2))
      allocate(mass            (i1:i2, ju1:j2, k1:k2))
      allocate(gridBoxHeight   (i1:i2, ju1:j2, k1:k2))
      allocate(tropopausePress (i1:i2, ju1:j2))
      allocate(relativeHumidity(i1:i2, ju1:j2, k1:k2))

      call Get_mass              (metFields, mass)
      call Get_press3c           (metFields, press3c)
      call Get_press3e           (metFields, press3e)
      call Get_gridBoxHeight     (metFields, gridBoxHeight)
      call Get_tropopausePress   (metFields, tropopausePress)
      call Get_relativeHumidity  (metFields, relativeHumidity)

      ! Get the Emission related variables
      call Get_curEmissionFileRecord(Emission, curEmissionFileRecord)
      call Get_emiss_dust_opt (Emission, emiss_dust_opt)
      call Get_ndust (Emission, ndust)
      if (emiss_dust_opt /= 0) then
         allocate(emiss_map_dust(ndust))
         call Get_emiss_map_dust (Emission, emiss_map_dust, ndust)
         allocate(emiss_dust(i1:i2, ju1:j2, ndust))
         call Get_emissDust (Emission, emiss_dust)
      end if

      call Get_do_semiss_inchem (Emission, do_semiss_inchem)
      call Get_do_ShipEmission  (Emission, do_ShipEmission)

      call Get_num_emiss        (Emission, num_emiss)
      allocate(emiss_map(num_emiss))
      call Get_emiss_map        (Emission, emiss_map, num_emiss)
      call Get_emiss_opt        (Emission, emiss_opt)
      call Get_emiss_timpyr     (Emission, emiss_timpyr)
      call Get_emissionArray (Emission, emissionArray)

!.sds      if (emiss_opt == 2) then
      if (btest(emiss_opt,1)) then
         allocate(emiss_nox  (i1:i2,ju1:j2))
         allocate(emiss_isop (i1:i2,ju1:j2))
         allocate(emiss_monot(i1:i2,ju1:j2))
         call Get_emiss_nox  (Emission, emiss_nox  )
         call Get_emiss_isop (Emission, emiss_isop )
         call Get_emiss_monot(Emission, emiss_monot)
         if (do_ShipEmission) then
            allocate(emiss_o3   (i1:i2,ju1:j2))
            allocate(emiss_hno3 (i1:i2,ju1:j2))
            call Get_emiss_o3   (Emission, emiss_o3   )
            call Get_emiss_hno3 (Emission, emiss_hno3 )
            call Get_o3_index   (Emission, o3_index   )
         end if
      end if

      ! Emisssion diagnostics
      if (pr_surf_emiss) then
         allocate(surf_emiss_out (i1:i2, ju1:j2, numSpecies))
         call Get_surf_emiss_out (Emission, surf_emiss_out )

         allocate(surf_emiss_out2(i1:i2, ju1:j2, 6))
         call Get_surf_emiss_out2(Emission, surf_emiss_out2)
      end if

      if (pr_emiss_3d) then
         allocate(emiss_3d_out(i1:i2, ju1:j2, k1:k2, numSpecies))
         call Get_emiss_3d_out(Emission, emiss_3d_out)
      end if

      if (pr_ascii5) then
         Allocate (loss(2, NUM_OPERATORS, 1:numSpecies))
         Allocate (prod(2, NUM_OPERATORS, 1:numSpecies))

         call Get_prod(SpeciesConcentration, prod)
         call Get_loss(SpeciesConcentration, loss)
      endif

!.sds.. getAerDust diagnostics for gocart mechanism
      call Get_pr_AerDust(Diagnostics, pr_AerDust)
      if ( pr_AerDust .and. TRIM(chem_mecha) == 'gocart_aerosol' ) then

         call runCalcAerosolDust (self, &
     &           mass, mcor, gridBoxHeight, SpeciesConcentration,  &
     &           gmiGrid, gmiDomain, Diagnostics, TRIM(chem_mecha))

         call Aero_OptDep_SurfArea ( gridBoxHeight, concentration, &
     &            self%OptDepth, self%Eradius, self%Tarea, self%Odaer, relativeHumidity, self%Daersl, self%Waersl, &
     &            RAA_b, QAA_b, numSpecies, 1.0d50, self%AerDust_Effect_opt, &
     &            i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, numSpecies, num_AerDust)

         call Dust_OptDep_SurfArea ( gridBoxHeight, &
     &            concentration, self%OptDepth, self%Eradius, self%Tarea, self%Odmdust, self%Dust, &
     &            RAA_b, QAA_b, numSpecies, 1.0d50, self%AerDust_Effect_opt, &
     &            i1, i2, ju1, j2, k1, k2, numSpecies, num_AerDust)

      endif
!
! Call the Chemistry control routine
                rootProc = .FALSE.
                if (procID==0) rootProc = .TRUE.
      call updateChemistry (self%smv2SavedVars, self%sulfSavedVars, &
                self%JXbundle, self%cloudGT, &
                self%chemintv, self%rateintv, self%jno_num, self%jnoxnum,&
                TRIM(chem_mecha), do_ftiming, metdata_name_org, metdata_name_model, &
     &          self%lossData, KDIM, JDIM, MDIM, pr_qj_o3_o1d, pr_qj_opt_depth,&
     &           do_qqjk_inchem, do_qqjk_reset, pr_qqjk, emiss_isop,           &
     &           emiss_monot, emiss_nox, emiss_hno3, emiss_o3,                 &
     &           self%emiss_ozone, self%prevEmissionFileRecord, &
     &           curEmissionFileRecord,  tropopausePress, press3c, press3e,    &
     &           averagePressEdge, gridBoxHeight, taucli, tauclw, self%qj_labels, &
     &           self%const_labels, pr_smv2, pr_sulf_src, do_aerocom,&
     &           pr_surf_emiss, pr_emiss_3d, constOutputFrequency, num_AerDust,&
     &           self%aqua_infile_name, self%qj_timpyr, self%qj_infile_name,   &
     &           qj_var_name, const_init_val, dlatr, latdeg, londeg, mcor,     &
     &           mass, self%loss_freq, prod, loss, concentration,              &
     &           self%overheadO3col, self%decay_3d_out,                        &
     &           self%qjgmi, self%qkgmi, emissionArray, self%rxnr_adjust, kel, &
     &           humidity, pctm2, totalCloudFraction, fracCloudCover,          &
     &           tau_cloud, lwi_flags, moistq, clwc, cmf, surf_air_temp,       &
     &           surf_alb_uv, self%dms_oh, self%dms_no3, self%so2_oh,          &
     &           self%so2_h2o2, self%so2_o3, emiss_dust, self%qqjgmi,          &
     &           self%qqkgmi, self%yda, self%qqkda, self%qqjda, surf_emiss_out,&
     &           surf_emiss_out2, emiss_3d_out, self%ch4clim, self%h2oclim,    &
     &           self%hno3gas, self%lbssad, self%sadgmi, self%hno3cond,        &
     &           self%h2oback, self%h2ocond, self%reffice, self%reffsts,       &
     &           self%vfall, self%optDepth, self%eRadius, self%tArea,          &
     &           self%odAer, relativeHumidity, radswg, self%odMdust, self%dust,&
     &           self%wAersl, self%dAersl, self%cloud_tau, self%cloud_param,   &
     &           self%flux_gt, self%cloudDroplet_opt, self%be_opt,             &
     &           self%chem_opt, self%sad_opt, self%phot_opt, self%fastj_opt,   &
     &           self%fastj_offset_sec, self%do_smv_reord, self%do_synoz,      &
     &           do_semiss_inchem, self%do_wetchem, self%do_aqu_chem,          &
     &         self%do_clear_sky, self%do_AerDust_Calc, self%do_ozone_inFastJX,&
     &           do_ShipEmission, o3_index, self%loss_freq_opt, self%kmin_loss,&
     &           self%kmax_loss, self%t_half_be7, self%t_half_be10,            &
     &           self%yield_be7, self%yield_be10, self%dehydmin, self%fbc_j1,  &
     &           self%fbc_j2, self%forc_bc_num, self%forc_bc_kmax,             &
     &           self%forc_bc_kmin, self%forc_bc_opt, self%forc_bc_map,        &
     &           self%forc_bc_incrpyr, self%forc_bc_start_num,                 &
     &           self%forc_bc_years, self%forc_bc_data, &
                 self%lastYearFBC, self%jlatmd, &
                 nymd, nhms, start_ymd, gmi_sec, tdt,      &
     &          num_time_steps, self%mw, pr_diag, procID, self%synoz_threshold,&
     &           self%chem_cycle, self%AerDust_Effect_opt, self%lbssad_timpyr, &
     &           self%h2oclim_timpyr, self%chem_mask_klo, self%chem_mask_khi,  &
     &           emiss_timpyr, emiss_opt, emiss_map, self%dehyd_opt,           &
     &           self%h2oclim_opt, self%lbssad_opt, emiss_dust_opt,            &
     &           emiss_map_dust, ndust, self%rxnr_adjust_map, self%ih2_num,    &
     &          self%io3_num, self%ih2o_num, self%ih2ocond_num, self%ihno3_num,&
     &           self%ihno3cond_num, self%idehyd_num, self%ih2oaircr_num,      &
     &           self%ich4_num, self%imgas_num, self%ico_num, self%ino_num,    &
     &           self%ipropene_num, self%iisoprene_num, self%initrogen_num,    &
     &           self%ioxygen_num, self%isynoz_num, num_emiss, numSpecies,     &
     &           self%num_qks, self%num_qjs, self%num_qjo, self%num_sad,       &
     &           self%num_molefrac, self%num_chem, self%num_active,            &
     &           self%num_rxnr_adjust, self%rxnr_adjust_timpyr, ilong, ilat,   &
     &           ivert, itloop, i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2, ilo, ihi, &
     &           julo, jhi, i1, i2, ju1, j2, k1, k2, numDomains, commuWorld,   &
     &           self%ih2o2_num, self%jno_adjust, rootProc,                    &
     &           self%linoz_infile_name, self%sf6_infile_name,                 &
     &           self%SO3daily_infile_name, self%SO3monthly_infile_name)

      if (pr_ascii5) then
         call Set_prod(SpeciesConcentration, prod)
         call Set_loss(SpeciesConcentration, loss)

         deallocate (loss)
         deallocate (prod)
      endif

      ! Call the Chemistry control routine
      call Set_do_qqjk_reset   (Diagnostics, do_qqjk_reset)

      if (do_ShipEmission) then
         call Set_emiss_o3   (Emission, emiss_o3   )
         deallocate(emiss_o3 )
      end if

      ! Emission diagnostics
      if (pr_surf_emiss) then
         call Set_surf_emiss_out (Emission, surf_emiss_out )
         call Set_surf_emiss_out2(Emission, surf_emiss_out2)
         deallocate(surf_emiss_out )
         deallocate(surf_emiss_out2)
      end if

      if (pr_emiss_3d) then
         call Set_emiss_3d_out(Emission, emiss_3d_out)
         deallocate(emiss_3d_out)
      end if

      if (metdata_name_org(1:4) == 'GMAO') then
         deallocate(taucli)
         deallocate(tauclw)
      end if

      deallocate(gridBoxHeight)
      deallocate(tropopausePress)
      deallocate(relativeHumidity)
      deallocate(mass)
      deallocate(press3e)
      deallocate(press3c)
      deallocate(mcor)
      deallocate(dlatr)
      deallocate(latdeg)
      deallocate(londeg)

  return

  end subroutine RunChemistry

!-------------------------------------------------------------------------

  subroutine FinalizeChemistry (self)

  implicit none

  type (t_Chemistry), intent(inout) :: self

  PRINT*,'  Finalize Chemistry'

  return

  end subroutine FinalizeChemistry

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  subroutine Allocate_rxnr_adjust_map (self)
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%rxnr_adjust_map(self%num_rxnr_adjust))
    return
  end subroutine Allocate_rxnr_adjust_map
!-------------------------------------------------------------------------
  subroutine Set_flux_gt (self, flux_gt)
    real*8            , intent(in)  :: flux_gt(:,:,:)
    type (t_Chemistry), intent(inOut)  :: self
    self%flux_gt(:,:,:) = flux_gt(:,:,:)
    return
  end subroutine Set_flux_gt
!-------------------------------------------------------------------------
  subroutine Get_flux_gt (self, flux_gt)
    real*8            , intent(out)  :: flux_gt(:,:,:)
    type (t_Chemistry), intent(in)  :: self
    flux_gt(:,:,:) = self%flux_gt(:,:,:)
    return
  end subroutine Get_flux_gt
!-------------------------------------------------------------------------
  subroutine Allocate_flux_gt (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%flux_gt(i1:i2, ju1:j2, k1:k2))
    self%flux_gt = 0.0d0
    return
  end subroutine Allocate_flux_gt
!-------------------------------------------------------------------------
  subroutine Set_cloud_tau (self, cloud_tau)
    real*8            , intent(in)  :: cloud_tau(:,:,:)
    type (t_Chemistry), intent(inOut)  :: self
    self%cloud_tau(:,:,:) = cloud_tau(:,:,:)
    return
  end subroutine Set_cloud_tau
!-------------------------------------------------------------------------
  subroutine Get_cloud_tau (self, cloud_tau)
    real*8            , intent(out)  :: cloud_tau(:,:,:)
    type (t_Chemistry), intent(in)  :: self
    cloud_tau(:,:,:) = self%cloud_tau(:,:,:)
    return
  end subroutine Get_cloud_tau
!-------------------------------------------------------------------------
  subroutine Allocate_cloud_tau (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%cloud_tau(i1:i2, ju1:j2, k1:k2))
    self%cloud_tau = 0.0d0
    return
  end subroutine Allocate_cloud_tau
!-------------------------------------------------------------------------
  subroutine Set_cloudDroplet_opt (self, cloudDroplet_opt)
    integer           , intent(in)  :: cloudDroplet_opt
    type (t_Chemistry), intent(inOut)  :: self
    self%cloudDroplet_opt = cloudDroplet_opt
    return
  end subroutine Set_cloudDroplet_opt
!-------------------------------------------------------------------------
  subroutine Get_cloudDroplet_opt (self, cloudDroplet_opt)
    integer           , intent(out)  :: cloudDroplet_opt
    type (t_Chemistry), intent(in)  :: self
    cloudDroplet_opt = self%cloudDroplet_opt
    return
  end subroutine Get_cloudDroplet_opt
!-------------------------------------------------------------------------
  subroutine Set_cloud_param (self, cloud_param)
    real*8            , intent(in)  :: cloud_param(:,:,:,:)
    type (t_Chemistry), intent(inOut)  :: self
    self%cloud_param(:,:,:,:) = cloud_param(:,:,:,:)
    return
  end subroutine Set_cloud_param
!-------------------------------------------------------------------------
  subroutine Get_cloud_param (self, cloud_param)
    real*8            , intent(out)  :: cloud_param(:,:,:,:)
    type (t_Chemistry), intent(in)  :: self
    cloud_param(:,:,:,:) = self%cloud_param(:,:,:,:)
    return
  end subroutine Get_cloud_param
!-------------------------------------------------------------------------
  subroutine Allocate_cloud_param (self, i1, i2, ju1, j2, k1, k2, numSpecies)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2, numSpecies
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%cloud_param(i1:i2, ju1:j2, k1:k2, numSpecies))
    self%cloud_param = 0.0d0
    return
  end subroutine Allocate_cloud_param
!-------------------------------------------------------------------------
  subroutine Allocate_qkgmi (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%qkgmi(i1:i2, ju1:j2, k1:k2, self%num_qks))
    self%qkgmi = 0.0d0
    return
  end subroutine Allocate_qkgmi
!-------------------------------------------------------------------------
  subroutine Allocate_qqjda (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%qqjda(i1:i2, ju1:j2, k1:k2, self%num_qjs))
    self%qqjda = 0.0d0
    return
  end subroutine Allocate_qqjda
!-------------------------------------------------------------------------
  subroutine Allocate_qqkda (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%qqkda(i1:i2, ju1:j2, k1:k2, self%num_qks))
    self%qqkda = 0.0d0
    return
  end subroutine Allocate_qqkda
!-------------------------------------------------------------------------
  subroutine Allocate_qqjgmi (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%qqjgmi(i1:i2, ju1:j2, k1:k2, self%num_qjo))
    self%qqjgmi = 0.0d0
    return
  end subroutine Allocate_qqjgmi
!-------------------------------------------------------------------------
  subroutine Allocate_qqkgmi (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%qqkgmi(i1:i2, ju1:j2, k1:k2, self%num_qks))
    self%qqkgmi = 0.0d0
    return
  end subroutine Allocate_qqkgmi
!-------------------------------------------------------------------------
  subroutine Allocate_yda (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%yda(i1:i2, ju1:j2, k1:k2, self%num_active))
    self%yda = 0.0d0
    return
  end subroutine Allocate_yda
!-------------------------------------------------------------------------
  subroutine Allocate_qjgmi (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%qjgmi(i1:i2, ju1:j2, k1:k2, self%num_qjo))
    self%qjgmi = 0.0d0
    return
  end subroutine Allocate_qjgmi
!-------------------------------------------------------------------------
  subroutine Allocate_rxnr_adjust (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%rxnr_adjust(i1:i2, ju1:j2, k1:k2, self%num_rxnr_adjust, self%rxnr_adjust_timpyr))
    self%rxnr_adjust = 0.0d0
    return
  end subroutine Allocate_rxnr_adjust
!-------------------------------------------------------------------------
  subroutine Allocate_overheadO3col (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%overheadO3col(i1:i2, ju1:j2, k1:k2))
    self%overheadO3col = 0.0d0
    return
  end subroutine Allocate_overheadO3col
!-------------------------------------------------------------------------
  subroutine Allocate_dms_no3 (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%dms_no3(i1:i2, ju1:j2, k1:k2))
    self%dms_no3 = 0.0d0
    return
  end subroutine Allocate_dms_no3
!-------------------------------------------------------------------------
  subroutine Allocate_dms_oh (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%dms_oh(i1:i2, ju1:j2, k1:k2))
    self%dms_oh = 0.0d0
    return
  end subroutine Allocate_dms_oh
!-------------------------------------------------------------------------
  subroutine Allocate_so2_h2o2 (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%so2_h2o2(i1:i2, ju1:j2, k1:k2))
    self%so2_h2o2 = 0.0d0
    return
  end subroutine Allocate_so2_h2o2
!-------------------------------------------------------------------------
  subroutine Allocate_so2_o3 (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%so2_o3(i1:i2, ju1:j2, k1:k2))
    self%so2_o3 = 0.0d0
    return
  end subroutine Allocate_so2_o3
!-------------------------------------------------------------------------
  subroutine Allocate_so2_oh (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%so2_oh(i1:i2, ju1:j2, k1:k2))
    self%so2_oh = 0.0d0
    return
  end subroutine Allocate_so2_oh
!-------------------------------------------------------------------------
  subroutine Allocate_hno3cond (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%hno3cond(i1:i2, ju1:j2, k1:k2))
    self%hno3cond = 0.0d0
    return
  end subroutine Allocate_hno3cond
!-------------------------------------------------------------------------
  subroutine Allocate_hno3gas (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%hno3gas(i1:i2, ju1:j2, k1:k2))
    self%hno3gas = 0.0d0
    return
  end subroutine Allocate_hno3gas
!-------------------------------------------------------------------------
  subroutine Allocate_h2oback (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%h2oback(i1:i2, ju1:j2, k1:k2))
    self%h2oback = 0.0d0
    return
  end subroutine Allocate_h2oback
!-------------------------------------------------------------------------
  subroutine Allocate_h2ocond (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%h2ocond(i1:i2, ju1:j2, k1:k2))
    self%h2ocond = 0.0d0
    return
  end subroutine Allocate_h2ocond
!-------------------------------------------------------------------------
  subroutine Allocate_reffice (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%reffice(i1:i2, ju1:j2, k1:k2))
    self%reffice = 0.0d0
    return
  end subroutine Allocate_reffice
!-------------------------------------------------------------------------
  subroutine Allocate_reffsts (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%reffsts(i1:i2, ju1:j2, k1:k2))
    self%reffsts = 0.0d0
    return
  end subroutine Allocate_reffsts
!-------------------------------------------------------------------------
  subroutine Allocate_vfall (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%vfall(i1:i2, ju1:j2, k1:k2))
    self%vfall = 0.0d0
    return
  end subroutine Allocate_vfall
!-------------------------------------------------------------------------
  subroutine Allocate_ch4clim (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%ch4clim(i1:i2, ju1:j2, k1:k2, self%h2oclim_timpyr))
    self%ch4clim = 0.0d0
    return
  end subroutine Allocate_ch4clim
!-------------------------------------------------------------------------
  subroutine Allocate_h2oclim (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%h2oclim(i1:i2, ju1:j2, k1:k2, self%h2oclim_timpyr))
    self%h2oclim = 0.0d0
    return
  end subroutine Allocate_h2oclim
!-------------------------------------------------------------------------
  subroutine Allocate_lbssad (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%lbssad(i1:i2, ju1:j2, k1:k2, self%lbssad_timpyr))
    self%lbssad = 0.0d0
    return
  end subroutine Allocate_lbssad
!-------------------------------------------------------------------------
  subroutine Allocate_sadgmi (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%sadgmi(i1:i2, ju1:j2, k1:k2, self%num_sad))
    self%sadgmi = 0.0d0
    return
  end subroutine Allocate_sadgmi
!-------------------------------------------------------------------------
  subroutine Allocate_loss_freq (self, i1, i2, ju1, j2, k1, k2, numSpecies)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2, numSpecies
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%loss_freq(i1:i2, ju1:j2, k1:k2, numSpecies))
    self%loss_freq = 0.0d0
    return
  end subroutine Allocate_loss_freq
!-------------------------------------------------------------------------
  subroutine Allocate_lossData (self)
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%lossData(KDIM, JDIM, MDIM, 1))
    self%lossData = 0.0d0
    return
  end subroutine Allocate_lossData
!-------------------------------------------------------------------------
  subroutine Allocate_odAer (self, i1, i2, ju1, j2, k1, k2, nSADaer, NRH_b)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2, nSADaer, NRH_b
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%odAer(i1:i2, ju1:j2, k1:k2, nSADaer*NRH_b))
    self%odAer = 0.0d0
    return
  end subroutine Allocate_odAer
!-------------------------------------------------------------------------
  subroutine Allocate_odMdust (self, i1, i2, ju1, j2, k1, k2, nSADdust)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2, nSADdust
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%odMdust(i1:i2, ju1:j2, k1:k2, nSADdust))
    self%odMdust = 0.0d0
    return
  end subroutine Allocate_odMdust
!-------------------------------------------------------------------------
  subroutine Allocate_eRadius (self, i1, i2, ju1, j2, k1, k2, nSADdust, nSADaer)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2, nSADdust, nSADaer
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%eRadius(i1:i2, ju1:j2, k1:k2, nSADdust+nSADaer))
    self%eRadius = 0.0d0
    return
  end subroutine Allocate_eRadius
!-------------------------------------------------------------------------
  subroutine Allocate_tArea (self, i1, i2, ju1, j2, k1, k2, nSADdust, nSADaer)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2, nSADdust, nSADaer
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%tArea(i1:i2, ju1:j2, k1:k2, nSADdust+nSADaer))
    self%tArea = 0.0d0
    return
  end subroutine Allocate_tArea
!-------------------------------------------------------------------------
  subroutine Allocate_optDepth (self, i1, i2, ju1, j2, k1, k2, num_AerDust)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2, num_AerDust
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%optDepth(i1:i2, ju1:j2, k1:k2, num_AerDust))
    self%optDepth = 0.0d0
    return
  end subroutine Allocate_optDepth
!-------------------------------------------------------------------------
  subroutine Allocate_dust (self, i1, i2, ju1, j2, k1, k2, nSADdust)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2, nSADdust
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%dust(i1:i2, ju1:j2, k1:k2, nSADdust))
    self%dust = 0.0d0
    return
  end subroutine Allocate_dust
!-------------------------------------------------------------------------
  subroutine Allocate_wAersl (self, i1, i2, ju1, j2, k1, k2, nSADaer)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2, nSADaer
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%wAersl(i1:i2, ju1:j2, k1:k2, nSADaer))
    self%wAersl = 0.0d0
    return
  end subroutine Allocate_wAersl
!-------------------------------------------------------------------------
  subroutine Allocate_dAersl (self, i1, i2, ju1, j2, k1, k2)
    integer           , intent(in   )  :: i1, i2, ju1, j2, k1, k2
    type (t_Chemistry), intent(inout)  :: self
    Allocate(self%dAersl(i1:i2, ju1:j2, k1:k2, 2))
    self%dAersl = 0.0d0
    return
  end subroutine Allocate_dAersl
!-------------------------------------------------------------------------
  subroutine Set_rxnr_adjust (self, rxnr_adjust)
    real*8          , intent(in)  :: rxnr_adjust(:,:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%rxnr_adjust(:,:,:,:,:) = rxnr_adjust(:,:,:,:,:)
    return
  end subroutine Set_rxnr_adjust
!-----------------------------------------------------------------------
  subroutine Set_qkgmi (self, qkgmi)
    real*8          , intent(in)  :: qkgmi(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%qkgmi(:,:,:,:) = qkgmi(:,:,:,:)
    return
  end subroutine Set_qkgmi
!-----------------------------------------------------------------------
  subroutine Set_qqjda (self, qqjda)
    real*8          , intent(in)  :: qqjda(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%qqjda(:,:,:,:) = qqjda(:,:,:,:)
    return
  end subroutine Set_qqjda
!-----------------------------------------------------------------------
  subroutine Set_qqkda (self, qqkda)
    real*8          , intent(in)  :: qqkda(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%qqkda(:,:,:,:) = qqkda(:,:,:,:)
    return
  end subroutine Set_qqkda
!-----------------------------------------------------------------------
  subroutine Set_yda (self, yda)
    real*8          , intent(in)  :: yda(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%yda(:,:,:,:) = yda(:,:,:,:)
    return
  end subroutine Set_yda
!-----------------------------------------------------------------------
  subroutine Set_qqjgmi (self, qqjgmi)
    real*8          , intent(in)  :: qqjgmi(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%qqjgmi(:,:,:,:) = qqjgmi(:,:,:,:)
    return
  end subroutine Set_qqjgmi
!-----------------------------------------------------------------------
  subroutine Set_qqkgmi (self, qqkgmi)
    real*8          , intent(in)  :: qqkgmi(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%qqkgmi(:,:,:,:) = qqkgmi(:,:,:,:)
    return
  end subroutine Set_qqkgmi
!-----------------------------------------------------------------------
  subroutine Set_qjgmi (self, qjgmi)
    real*8          , intent(in)  :: qjgmi(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%qjgmi(:,:,:,:) = qjgmi(:,:,:,:)
    return
  end subroutine Set_qjgmi
!-----------------------------------------------------------------------
  subroutine Set_overheadO3col (self, overheadO3col)
    real*8          , intent(in)  :: overheadO3col(:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%overheadO3col(:,:,:) = overheadO3col(:,:,:)
    return
  end subroutine Set_overheadO3col
!-----------------------------------------------------------------------
  subroutine Set_decay_3d_out (self, decay_3d_out)
    real*8          , intent(in)  :: decay_3d_out(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%decay_3d_out(:,:,:,:) = decay_3d_out(:,:,:,:)
    return
  end subroutine Set_decay_3d_out
!-----------------------------------------------------------------------
  subroutine Set_dms_no3 (self, dms_no3)
    real*8          , intent(in)  :: dms_no3(:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%dms_no3(:,:,:) = dms_no3(:,:,:)
    return
  end subroutine Set_dms_no3
!-----------------------------------------------------------------------
  subroutine Set_dms_oh (self, dms_oh)
    real*8          , intent(in)  :: dms_oh(:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%dms_oh(:,:,:) = dms_oh(:,:,:)
    return
  end subroutine Set_dms_oh
!-----------------------------------------------------------------------
  subroutine Set_so2_h2o2 (self, so2_h2o2)
    real*8          , intent(in)  :: so2_h2o2(:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%so2_h2o2(:,:,:) = so2_h2o2(:,:,:)
    return
  end subroutine Set_so2_h2o2
!-----------------------------------------------------------------------
  subroutine Set_so2_o3 (self, so2_o3)
    real*8          , intent(in)  :: so2_o3(:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%so2_o3(:,:,:) = so2_o3(:,:,:)
    return
  end subroutine Set_so2_o3
!-----------------------------------------------------------------------
  subroutine Set_so2_oh (self, so2_oh)
    real*8          , intent(in)  :: so2_oh(:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%so2_oh(:,:,:) = so2_oh(:,:,:)
    return
  end subroutine Set_so2_oh
!-----------------------------------------------------------------------
  subroutine Set_h2ocond (self, h2ocond)
    real*8          , intent(in)  :: h2ocond(:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%h2ocond(:,:,:) = h2ocond(:,:,:)
    return
  end subroutine Set_h2ocond
!-----------------------------------------------------------------------
  subroutine Set_optDepth (self, optDepth)
    real*8          , intent(in)  :: optDepth(:,:,:,:)
    type (t_Chemistry), intent(inout)   :: self
    self%optDepth(:,:,:,:) = optDepth(:,:,:,:)
    return
  end subroutine Set_optDepth
!-----------------------------------------------------------------------
  subroutine Get_rxnr_adjust (self, rxnr_adjust)
    real*8          , intent(out)  :: rxnr_adjust(:,:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    rxnr_adjust(:,:,:,:,:) = self%rxnr_adjust(:,:,:,:,:)
    return
  end subroutine Get_rxnr_adjust
!-----------------------------------------------------------------------
  subroutine Get_qkgmi (self, qkgmi)
    real*8          , intent(out)  :: qkgmi(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    qkgmi(:,:,:,:) = self%qkgmi(:,:,:,:)
    return
  end subroutine Get_qkgmi
!-----------------------------------------------------------------------
  subroutine Get_qqjda (self, qqjda)
    real*8          , intent(out)  :: qqjda(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    qqjda(:,:,:,:) = self%qqjda(:,:,:,:)
    return
  end subroutine Get_qqjda
!-----------------------------------------------------------------------
  subroutine Get_qqkda (self, qqkda)
    real*8          , intent(out)  :: qqkda(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    qqkda(:,:,:,:) = self%qqkda(:,:,:,:)
    return
  end subroutine Get_qqkda
!-----------------------------------------------------------------------
  subroutine Get_yda (self, yda)
    real*8          , intent(out)  :: yda(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    yda(:,:,:,:) = self%yda(:,:,:,:)
    return
  end subroutine Get_yda
!-----------------------------------------------------------------------
  subroutine Get_qqjgmi (self, qqjgmi)
    real*8          , intent(out)  :: qqjgmi(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    qqjgmi(:,:,:,:) = self%qqjgmi(:,:,:,:)
    return
  end subroutine Get_qqjgmi
!-----------------------------------------------------------------------
  subroutine Get_qqkgmi (self, qqkgmi)
    real*8          , intent(out)  :: qqkgmi(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    qqkgmi(:,:,:,:) = self%qqkgmi(:,:,:,:)
    return
  end subroutine Get_qqkgmi
!-----------------------------------------------------------------------
  subroutine Get_qjgmi (self, qjgmi)
    real*8          , intent(out)  :: qjgmi(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    qjgmi(:,:,:,:) = self%qjgmi(:,:,:,:)
    return
  end subroutine Get_qjgmi
!-----------------------------------------------------------------------
  subroutine Get_overheadO3col (self, overheadO3col)
    real*8          , intent(out)  :: overheadO3col(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    overheadO3col(:,:,:) = self%overheadO3col(:,:,:)
    return
  end subroutine Get_overheadO3col
!-----------------------------------------------------------------------
  subroutine Get_decay_3d_out (self, decay_3d_out)
    real*8          , intent(out)  :: decay_3d_out(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    decay_3d_out(:,:,:,:) = self%decay_3d_out(:,:,:,:)
    return
  end subroutine Get_decay_3d_out
!-----------------------------------------------------------------------
  subroutine Get_dms_no3 (self, dms_no3)
    real*8          , intent(out)  :: dms_no3(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    dms_no3(:,:,:) = self%dms_no3(:,:,:)
    return
  end subroutine Get_dms_no3
!-------------------------------------------------------------------------
  subroutine Get_dms_oh (self, dms_oh)
    real*8          , intent(out)  :: dms_oh(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    dms_oh(:,:,:) = self%dms_oh(:,:,:)
    return
  end subroutine Get_dms_oh
!-------------------------------------------------------------------------
  subroutine Get_so2_h2o2 (self, so2_h2o2)
    real*8          , intent(out)  :: so2_h2o2(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    so2_h2o2(:,:,:) = self%so2_h2o2(:,:,:)
    return
  end subroutine Get_so2_h2o2
!-------------------------------------------------------------------------
  subroutine Get_so2_o3 (self, so2_o3)
    real*8          , intent(out)  :: so2_o3(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    so2_o3(:,:,:) = self%so2_o3(:,:,:)
    return
  end subroutine Get_so2_o3
!-------------------------------------------------------------------------
  subroutine Get_so2_oh (self, so2_oh)
    real*8          , intent(out)  :: so2_oh(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    so2_oh(:,:,:) = self%so2_oh(:,:,:)
    return
  end subroutine Get_so2_oh
!-------------------------------------------------------------------------
  subroutine Get_hno3cond (self, hno3cond)
    real*8          , intent(out)  :: hno3cond(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    hno3cond(:,:,:) = self%hno3cond(:,:,:)
    return
  end subroutine Get_hno3cond
!-------------------------------------------------------------------------
  subroutine Get_hno3gas (self, hno3gas)
    real*8          , intent(out)  :: hno3gas(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    hno3gas(:,:,:) = self%hno3gas(:,:,:)
    return
  end subroutine Get_hno3gas
!-------------------------------------------------------------------------
  subroutine Get_h2ocond (self, h2ocond)
    real*8          , intent(out)  :: h2ocond(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    h2ocond(:,:,:) = self%h2ocond(:,:,:)
    return
  end subroutine Get_h2ocond
!-------------------------------------------------------------------------
  subroutine Get_h2oback (self, h2oback)
    real*8          , intent(out)  :: h2oback(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    h2oback(:,:,:) = self%h2oback(:,:,:)
    return
  end subroutine Get_h2oback
!-------------------------------------------------------------------------
  subroutine Get_reffice (self, reffice)
    real*8          , intent(out)  :: reffice(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    reffice(:,:,:) = self%reffice(:,:,:)
    return
  end subroutine Get_reffice
!-------------------------------------------------------------------------
  subroutine Get_reffsts (self, reffsts)
    real*8          , intent(out)  :: reffsts(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    reffsts(:,:,:) = self%reffsts(:,:,:)
    return
  end subroutine Get_reffsts
!-------------------------------------------------------------------------
  subroutine Get_vfall (self, vfall)
    real*8          , intent(out)  :: vfall(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    vfall(:,:,:) = self%vfall(:,:,:)
    return
  end subroutine Get_vfall
!-------------------------------------------------------------------------
  subroutine Get_sadgmi (self, sadgmi)
    real*8          , intent(out)  :: sadgmi(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    sadgmi(:,:,:,:) = self%sadgmi(:,:,:,:)
    return
  end subroutine Get_sadgmi
!-------------------------------------------------------------------------
  subroutine Get_optDepth (self, optDepth)
    real*8          , intent(out)  :: optDepth(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    optDepth(:,:,:,:) = self%optDepth(:,:,:,:)
    return
  end subroutine Get_optDepth
!-------------------------------------------------------------------------
  subroutine Get_odAer (self, odAer)
    real*8          , intent(out)  :: odAer(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    odAer(:,:,:,:) = self%odAer(:,:,:,:)
    return
  end subroutine Get_odAer
!-------------------------------------------------------------------------
  subroutine Get_odMdust (self, odMdust)
    real*8          , intent(out)  :: odMdust(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    odMdust(:,:,:,:) = self%odMdust(:,:,:,:)
    return
  end subroutine Get_odMdust
!-------------------------------------------------------------------------
  subroutine Get_eRadius (self, eRadius)
    real*8          , intent(out)  :: eRadius(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    eRadius(:,:,:,:) = self%eRadius(:,:,:,:)
    return
  end subroutine Get_eRadius
!-------------------------------------------------------------------------
  subroutine Get_tArea (self, tArea)
    real*8          , intent(out)  :: tArea(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    tArea(:,:,:,:) = self%tArea(:,:,:,:)
    return
  end subroutine Get_tArea
!-------------------------------------------------------------------------
  subroutine Get_dust (self, dust)
    real*8          , intent(out)  :: dust(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    dust(:,:,:,:) = self%dust(:,:,:,:)
    return
  end subroutine Get_dust
!-------------------------------------------------------------------------
  subroutine Get_wAersl (self, wAersl)
    real*8          , intent(out)  :: wAersl(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    wAersl(:,:,:,:) = self%wAersl(:,:,:,:)
    return
  end subroutine Get_wAersl
!-------------------------------------------------------------------------
  subroutine Get_dAersl (self, dAersl)
    real*8          , intent(out)  :: dAersl(:,:,:,:)
    type (t_Chemistry), intent(in)   :: self
    dAersl(:,:,:,:) = self%dAersl(:,:,:,:)
    return
  end subroutine Get_dAersl
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  subroutine Get_chem_opt (self, chem_opt)
    integer         , intent(out)  :: chem_opt
    type (t_Chemistry), intent(in)   :: self
    chem_opt = self%chem_opt
    return
  end subroutine Get_chem_opt
!-------------------------------------------------------------------------
  subroutine Get_loss_opt (self, loss_opt)
    integer         , intent(out)  :: loss_opt
    type (t_Chemistry), intent(in)   :: self
    loss_opt = self%loss_opt
    return
  end subroutine Get_loss_opt
!-------------------------------------------------------------------------
  subroutine Get_oz_eq_synoz_opt (self, oz_eq_synoz_opt)
    integer         , intent(out)  :: oz_eq_synoz_opt
    type (t_Chemistry), intent(in)   :: self
    oz_eq_synoz_opt = self%oz_eq_synoz_opt
    return
  end subroutine Get_oz_eq_synoz_opt
!-------------------------------------------------------------------------
  subroutine Get_sad_opt (self, sad_opt)
    integer         , intent(out)  :: sad_opt
    type (t_Chemistry), intent(in)   :: self
    sad_opt = self%sad_opt
    return
  end subroutine Get_sad_opt
!-------------------------------------------------------------------------
  subroutine Get_forc_bc_opt (self, forc_bc_opt)
    integer         , intent(out)  :: forc_bc_opt
    type (t_Chemistry), intent(in)   :: self
    forc_bc_opt = self%forc_bc_opt
    return
  end subroutine Get_forc_bc_opt
!-------------------------------------------------------------------------
  subroutine Get_lbssad_opt (self, lbssad_opt)
    integer         , intent(out)  :: lbssad_opt
    type (t_Chemistry), intent(in)   :: self
    lbssad_opt = self%lbssad_opt
    return
  end subroutine Get_lbssad_opt
!-------------------------------------------------------------------------
  subroutine Get_dehyd_opt (self, dehyd_opt)
    integer         , intent(out)  :: dehyd_opt
    type (t_Chemistry), intent(in)   :: self
    dehyd_opt = self%dehyd_opt
    return
  end subroutine Get_dehyd_opt
!-------------------------------------------------------------------------
  subroutine Get_h2oclim_opt (self, h2oclim_opt)
    integer         , intent(out)  :: h2oclim_opt
    type (t_Chemistry), intent(in)   :: self
    h2oclim_opt = self%h2oclim_opt
    return
  end subroutine Get_h2oclim_opt
!-------------------------------------------------------------------------
  subroutine Get_phot_opt (self, phot_opt)
    integer         , intent(out)  :: phot_opt
    type (t_Chemistry), intent(in)   :: self
    phot_opt = self%phot_opt
    return
  end subroutine Get_phot_opt
!-------------------------------------------------------------------------
  subroutine Get_fastj_opt (self, fastj_opt)
    integer         , intent(out)  :: fastj_opt
    type (t_Chemistry), intent(in)   :: self
    fastj_opt = self%fastj_opt
    return
  end subroutine Get_fastj_opt
!-------------------------------------------------------------------------
  subroutine Get_uvalbedo_data (self, uvalbedo_data)
    real*8         , intent(out)  :: uvalbedo_data(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    uvalbedo_data(:,:,:) = self%uvalbedo_data(:,:,:)
    return
  end subroutine Get_uvalbedo_data
!-------------------------------------------------------------------------
  subroutine Get_sasdir_data (self, sasdir_data)
    real*8         , intent(out)  :: sasdir_data(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    sasdir_data(:,:,:) = self%sasdir_data(:,:,:)
    return
  end subroutine Get_sasdir_data
!-------------------------------------------------------------------------
  subroutine Get_sasdif_data (self, sasdif_data)
    real*8         , intent(out)  :: sasdif_data(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    sasdif_data(:,:,:) = self%sasdif_data(:,:,:)
    return
  end subroutine Get_sasdif_data
!-------------------------------------------------------------------------
  subroutine Get_saldir_data (self, saldir_data)
    real*8         , intent(out)  :: saldir_data(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    saldir_data(:,:,:) = self%saldir_data(:,:,:)
    return
  end subroutine Get_saldir_data
!-------------------------------------------------------------------------
  subroutine Get_saldif_data (self, saldif_data)
    real*8         , intent(out)  :: saldif_data(:,:,:)
    type (t_Chemistry), intent(in)   :: self
    saldif_data(:,:,:) = self%saldif_data(:,:,:)
    return
  end subroutine Get_saldif_data
!-------------------------------------------------------------------------
  subroutine Get_sfalbedo_opt (self, sfalbedo_opt)
    integer         , intent(out)  :: sfalbedo_opt
    type (t_Chemistry), intent(in)   :: self
    sfalbedo_opt = self%sfalbedo_opt
    return
  end subroutine Get_sfalbedo_opt
!-------------------------------------------------------------------------
  subroutine Get_uvalbedo_opt (self, uvalbedo_opt)
    integer         , intent(out)  :: uvalbedo_opt
    type (t_Chemistry), intent(in)   :: self
    uvalbedo_opt = self%uvalbedo_opt
    return
  end subroutine Get_uvalbedo_opt
!-------------------------------------------------------------------------
  subroutine Get_h2oclim_timpyr (self, h2oclim_timpyr)
    integer         , intent(out)  :: h2oclim_timpyr
    type (t_Chemistry), intent(in)   :: self
    h2oclim_timpyr = self%h2oclim_timpyr
    return
  end subroutine Get_h2oclim_timpyr
!-------------------------------------------------------------------------
  subroutine Get_lbssad_timpyr (self, lbssad_timpyr)
    integer         , intent(out)  :: lbssad_timpyr
    type (t_Chemistry), intent(in)   :: self
    lbssad_timpyr = self%lbssad_timpyr
    return
  end subroutine Get_lbssad_timpyr
!-------------------------------------------------------------------------
  subroutine Get_fbc_j1 (self, fbc_j1)
    integer         , intent(out)  :: fbc_j1
    type (t_Chemistry), intent(in)   :: self
    fbc_j1 = self%fbc_j1
    return
  end subroutine Get_fbc_j1
!-------------------------------------------------------------------------
  subroutine Get_fbc_j2 (self, fbc_j2)
    integer         , intent(out)  :: fbc_j2
    type (t_Chemistry), intent(in)   :: self
    fbc_j2 = self%fbc_j2
    return
  end subroutine Get_fbc_j2
!-------------------------------------------------------------------------
  subroutine Get_forc_bc_years (self, forc_bc_years)
    integer         , intent(out)  :: forc_bc_years
    type (t_Chemistry), intent(in)   :: self
    forc_bc_years = self%forc_bc_years
    return
  end subroutine Get_forc_bc_years
!-------------------------------------------------------------------------
  subroutine Get_loss_freq_opt (self, loss_freq_opt)
    integer         , intent(out)  :: loss_freq_opt
    type (t_Chemistry), intent(in)   :: self
    loss_freq_opt = self%loss_freq_opt
    return
  end subroutine Get_loss_freq_opt
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  subroutine Get_do_AerDust_Calc (self, do_AerDust_Calc)
    logical          , intent(out)  :: do_AerDust_Calc
    type (t_Chemistry), intent(in )  :: self
    do_AerDust_Calc = self%do_AerDust_Calc
    return
  end subroutine Get_do_AerDust_Calc
!-------------------------------------------------------------------------
  subroutine Get_do_chem_grp (self, do_chem_grp)
    logical          , intent(out)  :: do_chem_grp
    type (t_Chemistry), intent(in )  :: self
    do_chem_grp = self%do_chem_grp
    return
  end subroutine Get_do_chem_grp
!-------------------------------------------------------------------------
  subroutine Get_do_rxnr_adjust (self, do_rxnr_adjust)
    logical          , intent(out)  :: do_rxnr_adjust
    type (t_Chemistry), intent(in )  :: self
    do_rxnr_adjust = self%do_rxnr_adjust
    return
  end subroutine Get_do_rxnr_adjust
!-------------------------------------------------------------------------
  subroutine Get_do_clear_sky (self, do_clear_sky)
    logical          , intent(out)  :: do_clear_sky
    type (t_Chemistry), intent(in )  :: self
    do_clear_sky = self%do_clear_sky
    return
  end subroutine Get_do_clear_sky
!-------------------------------------------------------------------------
  subroutine Get_do_solar_cycle (self, do_solar_cycle)
    logical          , intent(out)  :: do_solar_cycle
    type (t_Chemistry), intent(in )  :: self
    do_solar_cycle = self%do_solar_cycle
    return
  end subroutine Get_do_solar_cycle
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  subroutine Get_synoz_threshold (self, synoz_threshold)
    real*8           , intent(out)  :: synoz_threshold
    type (t_Chemistry), intent(in )  :: self
    synoz_threshold = self%synoz_threshold
    return
  end subroutine Get_synoz_threshold
!-------------------------------------------------------------------------
  subroutine Get_t_cloud_ice (self, t_cloud_ice)
    real*8           , intent(out)  :: t_cloud_ice
    type (t_Chemistry), intent(in )  :: self
    t_cloud_ice = self%t_cloud_ice
    return
  end subroutine Get_t_cloud_ice
!-------------------------------------------------------------------------
  subroutine Get_do_full_chem (self, do_full_chem)
    implicit none
    logical          , intent(out)  :: do_full_chem
    type (t_Chemistry), intent(in )  :: self
    do_full_chem = self%do_full_chem
    return
  end subroutine Get_do_full_chem
!-------------------------------------------------------------------------
  subroutine Get_num_chem (self, num_chem)
    implicit none
    integer           , intent(out)  :: num_chem
    type (t_Chemistry), intent(in )  :: self
    num_chem = self%num_chem
    return
  end subroutine Get_num_chem
!-------------------------------------------------------------------------
  subroutine Get_num_molefrac (self, num_molefrac)
    implicit none
    integer           , intent(out)  :: num_molefrac
    type (t_Chemistry), intent(in )  :: self
    num_molefrac = self%num_molefrac
    return
  end subroutine Get_num_molefrac
!-------------------------------------------------------------------------
  subroutine Get_num_active (self, num_active)
    implicit none
    integer           , intent(out)  :: num_active
    type (t_Chemistry), intent(in )  :: self
    num_active = self%num_active
    return
  end subroutine Get_num_active
!-------------------------------------------------------------------------
  subroutine Get_num_qjs (self, num_qjs)
    implicit none
    integer           , intent(out)  :: num_qjs
    type (t_Chemistry), intent(in )  :: self
    num_qjs = self%num_qjs
    return
  end subroutine Get_num_qjs
!-------------------------------------------------------------------------
  subroutine Get_num_qjo (self, num_qjo)
    implicit none
    integer           , intent(out)  :: num_qjo
    type (t_Chemistry), intent(in )  :: self
    num_qjo = self%num_qjo
    return
  end subroutine Get_num_qjo
!-------------------------------------------------------------------------
  subroutine Get_num_qks (self, num_qks)
    implicit none
    integer           , intent(out)  :: num_qks
    type (t_Chemistry), intent(in )  :: self
    num_qks = self%num_qks
    return
  end subroutine Get_num_qks
!-------------------------------------------------------------------------
  subroutine Get_num_spc_sbc (self, num_spc_sbc)
    implicit none
    integer           , intent(out)  :: num_spc_sbc
    type (t_Chemistry), intent(in )  :: self
    num_spc_sbc = self%num_spc_sbc
    return
  end subroutine Get_num_spc_sbc
!-------------------------------------------------------------------------
  subroutine Get_num_sad (self, num_sad)
    implicit none
    integer           , intent(out)  :: num_sad
    type (t_Chemistry), intent(in )  :: self
    num_sad = self%num_sad
    return
  end subroutine Get_num_sad
!-------------------------------------------------------------------------
  subroutine Get_num_ks_sbc (self, num_ks_sbc)
    implicit none
    integer           , intent(out)  :: num_ks_sbc
    type (t_Chemistry), intent(in )  :: self
    num_ks_sbc = self%num_ks_sbc
    return
  end subroutine Get_num_ks_sbc
!-------------------------------------------------------------------------
  subroutine Get_surf_bc_map (self, surf_bc_map)
    implicit none
    integer           , intent(out)  :: surf_bc_map(:)
    type (t_Chemistry), intent(in )  :: self
    surf_bc_map(:) = self%surf_bc_map(:)
    return
  end subroutine Get_surf_bc_map
!-------------------------------------------------------------------------
  subroutine Get_const_labels (self, const_labels)
    implicit none
    character (len=MAX_LENGTH_SPECIES_NAME), intent(out)  :: const_labels(:)
    type (t_Chemistry), intent(in )  :: self
    const_labels(:) = self%const_labels(:)
    return
  end subroutine Get_const_labels
!-------------------------------------------------------------------------
  subroutine Get_qk_labels (self, qk_labels)
    implicit none
    character (len=MAX_LENGTH_LABELS), intent(out)  :: qk_labels(:)
    type (t_Chemistry), intent(in )  :: self
    qk_labels(:) = self%qk_labels(:)
    return
  end subroutine Get_qk_labels
!-------------------------------------------------------------------------
  subroutine Get_qj_labels (self, qj_labels)
    implicit none
    character (len=MAX_LENGTH_LABELS), intent(out)  :: qj_labels(:)
    type (t_Chemistry), intent(in )  :: self
    qj_labels(:) = self%qj_labels(:)
    return
  end subroutine Get_qj_labels
!-------------------------------------------------------------------------
  subroutine Get_mw (self, mw)
    implicit none
    real*8, intent(out)  :: mw(:)
    type (t_Chemistry), intent(in )  :: self
    mw(:) = self%mw(:)
    return
  end subroutine Get_mw
!-------------------------------------------------------------------------
  subroutine Get_io3_num (self, io3_num)
    implicit none
    integer           , intent(out)  :: io3_num
    type (t_Chemistry), intent(in )  :: self
    io3_num = self%io3_num
    return
  end subroutine Get_io3_num
!-------------------------------------------------------------------------
  subroutine Get_ih2oaircr_num (self, ih2oaircr_num)
    implicit none
    integer           , intent(out)  :: ih2oaircr_num
    type (t_Chemistry), intent(in )  :: self
    ih2oaircr_num = self%ih2oaircr_num
    return
  end subroutine Get_ih2oaircr_num
!-------------------------------------------------------------------------
  subroutine Get_ih2o2_num (self, ih2o2_num)
    implicit none
    integer           , intent(out)  :: ih2o2_num
    type (t_Chemistry), intent(in )  :: self
    ih2o2_num = self%ih2o2_num
    return
  end subroutine Get_ih2o2_num
!-------------------------------------------------------------------------
  subroutine Get_ih2o_num (self, ih2o_num)
    implicit none
    integer           , intent(out)  :: ih2o_num
    type (t_Chemistry), intent(in )  :: self
    ih2o_num = self%ih2o_num
    return
  end subroutine Get_ih2o_num
!-------------------------------------------------------------------------
  subroutine Get_ih2_num (self, ih2_num)
    implicit none
    integer           , intent(out)  :: ih2_num
    type (t_Chemistry), intent(in )  :: self
    ih2_num = self%ih2_num
    return
  end subroutine Get_ih2_num
!-------------------------------------------------------------------------
  subroutine Get_initrogen_num (self, initrogen_num)
    implicit none
    integer           , intent(out)  :: initrogen_num
    type (t_Chemistry), intent(in )  :: self
    initrogen_num = self%initrogen_num
    return
  end subroutine Get_initrogen_num
!-------------------------------------------------------------------------
  subroutine Get_ino_num (self, ino_num)
    implicit none
    integer           , intent(out)  :: ino_num
    type (t_Chemistry), intent(in )  :: self
    ino_num = self%ino_num
    return
  end subroutine Get_ino_num
!-------------------------------------------------------------------------
  subroutine Get_iacetone_num (self, iacetone_num)
    implicit none
    integer           , intent(out)  :: iacetone_num
    type (t_Chemistry), intent(in )  :: self
    iacetone_num = self%iacetone_num
    return
  end subroutine Get_iacetone_num
!-------------------------------------------------------------------------
  subroutine Get_ipropene_num (self, ipropene_num)
    implicit none
    integer           , intent(out)  :: ipropene_num
    type (t_Chemistry), intent(in )  :: self
    ipropene_num = self%ipropene_num
    return
  end subroutine Get_ipropene_num
!-------------------------------------------------------------------------
  subroutine Get_iisoprene_num (self, iisoprene_num)
    implicit none
    integer           , intent(out)  :: iisoprene_num
    type (t_Chemistry), intent(in )  :: self
    iisoprene_num = self%iisoprene_num
    return
  end subroutine Get_iisoprene_num
!-------------------------------------------------------------------------
  subroutine Get_isynoz_num (self, isynoz_num)
    implicit none
    integer           , intent(out)  :: isynoz_num
    type (t_Chemistry), intent(in )  :: self
    isynoz_num = self%isynoz_num
    return
  end subroutine Get_isynoz_num
!-------------------------------------------------------------------------
  subroutine Get_noy_map (self, noy_map)
    implicit none
    integer           , intent(out)  :: noy_map(:)
    type (t_Chemistry), intent(in )  :: self
    noy_map(:) = self%noy_map(:)
    return
  end subroutine Get_noy_map
!-------------------------------------------------------------------------
  subroutine Get_nox_map (self, nox_map)
    implicit none
    integer           , intent(out)  :: nox_map(:)
    type (t_Chemistry), intent(in )  :: self
    nox_map(:) = self%nox_map(:)
    return
  end subroutine Get_nox_map
!-------------------------------------------------------------------------
  subroutine Get_do_nodoz (self, do_nodoz)
    implicit none
    logical           , intent(out)  :: do_nodoz
    type (t_Chemistry), intent(in )  :: self
    do_nodoz = self%do_nodoz
    return
  end subroutine Get_do_nodoz
!-------------------------------------------------------------------------
  subroutine Get_do_synoz (self, do_synoz)
    implicit none
    logical           , intent(out)  :: do_synoz
    type (t_Chemistry), intent(in )  :: self
    do_synoz = self%do_synoz
    return
  end subroutine Get_do_synoz
!-------------------------------------------------------------------------
  subroutine Get_num_noy (self, num_noy)
    implicit none
    integer           , intent(out)  :: num_noy
    type (t_Chemistry), intent(in )  :: self
    num_noy = self%num_noy
    return
  end subroutine Get_num_noy
!-------------------------------------------------------------------------
  subroutine Get_num_nox (self, num_nox)
    implicit none
    integer           , intent(out)  :: num_nox
    type (t_Chemistry), intent(in )  :: self
    num_nox = self%num_nox
    return
  end subroutine Get_num_nox
!-------------------------------------------------------------------------
  subroutine Get_ico_num (self, ico_num)
    implicit none
    integer           , intent(out)  :: ico_num
    type (t_Chemistry), intent(in )  :: self
    ico_num = self%ico_num
    return
  end subroutine Get_ico_num
!-------------------------------------------------------------------------
  subroutine Get_imgas_num (self, imgas_num)
    implicit none
    integer           , intent(out)  :: imgas_num
    type (t_Chemistry), intent(in )  :: self
    imgas_num = self%imgas_num
    return
  end subroutine Get_imgas_num
!-------------------------------------------------------------------------
  subroutine Get_ihno3_num (self, ihno3_num)
    implicit none
    integer           , intent(out)  :: ihno3_num
    type (t_Chemistry), intent(in )  :: self
    ihno3_num = self%ihno3_num
    return
  end subroutine Get_ihno3_num
!-------------------------------------------------------------------------
  subroutine Get_ihcl_num (self, ihcl_num)
    implicit none
    integer           , intent(out)  :: ihcl_num
    type (t_Chemistry), intent(in )  :: self
    ihcl_num = self%ihcl_num
    return
  end subroutine Get_ihcl_num
!-------------------------------------------------------------------------
  subroutine Get_ih2ocond_num (self, ih2ocond_num)
    implicit none
    integer           , intent(out)  :: ih2ocond_num
    type (t_Chemistry), intent(in )  :: self
    ih2ocond_num = self%ih2ocond_num
    return
  end subroutine Get_ih2ocond_num
!-------------------------------------------------------------------------
  subroutine Get_ihno3cond_num (self, ihno3cond_num)
    implicit none
    integer           , intent(out)  :: ihno3cond_num
    type (t_Chemistry), intent(in )  :: self
    ihno3cond_num = self%ihno3cond_num
    return
  end subroutine Get_ihno3cond_num
!-------------------------------------------------------------------------
  subroutine Get_ioxygen_num (self, ioxygen_num)
    implicit none
    integer           , intent(out)  :: ioxygen_num
    type (t_Chemistry), intent(in )  :: self
    ioxygen_num = self%ioxygen_num
    return
  end subroutine Get_ioxygen_num
!-------------------------------------------------------------------------
  subroutine Set_dehydmin (self, dehydmin)
    implicit none
    real*8            , intent(in)  :: dehydmin
    type (t_Chemistry), intent(inOut)  :: self
    self%dehydmin = dehydmin
    return
  end subroutine Set_dehydmin
!-------------------------------------------------------------------------
  subroutine Get_dehydmin (self, dehydmin)
    implicit none
    real*8            , intent(out)  :: dehydmin
    type (t_Chemistry), intent(in )  :: self
    dehydmin = self%dehydmin
    return
  end subroutine Get_dehydmin
!-------------------------------------------------------------------------
  subroutine Get_idehyd_num (self, idehyd_num)
    implicit none
    integer           , intent(out)  :: idehyd_num
    type (t_Chemistry), intent(in )  :: self
    idehyd_num = self%idehyd_num
    return
  end subroutine Get_idehyd_num
!-------------------------------------------------------------------------
  subroutine Get_in2o_num (self, in2o_num)
    implicit none
    integer           , intent(out)  :: in2o_num
    type (t_Chemistry), intent(in )  :: self
    in2o_num = self%in2o_num
    return
  end subroutine Get_in2o_num
!-------------------------------------------------------------------------
  subroutine Get_ich4_num (self, ich4_num)
    implicit none
    integer           , intent(out)  :: ich4_num
    type (t_Chemistry), intent(in )  :: self
    ich4_num = self%ich4_num
    return
  end subroutine Get_ich4_num
!-------------------------------------------------------------------------
  subroutine Get_ibrono2_num (self, ibrono2_num)
    implicit none
    integer           , intent(out)  :: ibrono2_num
    type (t_Chemistry), intent(in )  :: self
    ibrono2_num = self%ibrono2_num
    return
  end subroutine Get_ibrono2_num
!-------------------------------------------------------------------------
      subroutine setupChemicalGroup (i1, i2, ju1, j2, k1, k2)

      implicit none

#     include "setkin_group.h"

      integer :: num_cgrp
      integer :: num_cgrp_elem  (MAX_NUM_CGRP)
      integer :: cgrp_elem_map  (MAX_NUM_SMARRAY, MAX_NUM_CGRP)
      real*8  :: cgrp_fac  (MAX_NUM_SMARRAY, MAX_NUM_CGRP)

      integer, intent(in) :: i1, i2, ju1, j2, k1, k2

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ii, jj


!     ----------------
!     Begin Execution.
!     ----------------

      num_cgrp = NUMGRP

      num_cgrp_elem(:) = 0


      do jj = 1, num_cgrp

        IILOOP1: do ii = 1, MAXGRP_ELEM

          if (sgrp_elem_map(ii,jj) /= 0) then

            num_cgrp_elem(jj)    = num_cgrp_elem(jj) + 1

            cgrp_elem_map(ii,jj) = sgrp_elem_map(ii,jj)

            cgrp_fac(ii,jj) = sgrp_fac(ii,jj)

          else

            exit IILOOP1

          end if

        end do IILOOP1

      end do

      return

      end subroutine setupChemicalGroup

!-----------------------------------------------------------------------------
!   This routine sets up the synthetic species info.
!-----------------------------------------------------------------------------

      subroutine Setup_Syn_Spc (self, numSpecies)

      implicit none

#     include "setkin_synspc.h"
   
      integer, intent(in) :: numSpecies
      type (t_Chemistry), intent(inOut) :: self

      integer :: ii

!      self%do_synoz   = .false.
!      self%do_nodoz   = .false.
!
!      self%isynoz_num = 0
!
!      self%num_nox    = 0
!      self%num_noy    = 0
!
!      self%nox_map(:) = 0
!      self%noy_map(:) = 0


!     ==================
      if (self%chem_opt == 2) then
!     ==================

        if (USE_SYNOZ) then

          self%do_synoz   = .true.

          self%isynoz_num = ISYNOZ

          if (USE_NODOZ) then

            self%do_nodoz = .true.

            NOXLOOP: do ii = 1, MAXNODOZ_ELEM

              if (NOX_ELEM_MAP(ii) /= 0) then

                self%num_nox = self%num_nox + 1

                self%nox_map(ii) = NOX_ELEM_MAP(ii)

              else

                exit NOXLOOP

              end if

            end do NOXLOOP

            NOYLOOP: do ii = 1, MAXNODOZ_ELEM

              if (NOY_ELEM_MAP(ii) /= 0) then

                self%num_noy = self%num_noy + 1

                self%noy_map(ii) = NOY_ELEM_MAP(ii)

              else

                exit NOYLOOP

              end if

            end do NOYLOOP

          end if

        end if

!     =======================
      else if (self%chem_opt == 5) then
!     =======================

        self%do_synoz   = .true.

        self%isynoz_num = 1

        if (numSpecies == 2) then

          self%do_nodoz = .true.

          self%ihno3_num = 2

        end if

!     =======================
      else if (self%chem_opt == 9) then
!     =======================

        self%do_synoz   = .true.

        self%isynoz_num = ISYNOZ

!     ======
      end if
!     ======

      return

      end subroutine Setup_Syn_Spc



!-----------------------------------------------------------------------------
! Does error checking to make sure that resource file variables are properly
! set.
!-----------------------------------------------------------------------------

      subroutine rcCheckChemistrySetting(self, pr_sad, pr_qj, pr_qk, pr_qqjk, &
     &       do_qqjk_inchem, do_mean, emiss_opt, do_semiss_inchem, numSpecies,&
     &        met_opt)
!
      implicit none

#     include "smv2chem_par.h"
!
      logical, intent(in) :: pr_sad, pr_qj, pr_qk, pr_qqjk
      logical, intent(in) :: do_qqjk_inchem, do_mean
      logical          , intent(in) :: do_semiss_inchem
      integer          , intent(in) :: emiss_opt, met_opt
      integer          , intent(in) :: numSpecies
      type(t_Chemistry), intent(in) :: self
!
      character(len=90) :: err_msg
      integer           :: rc, icycles, ic
      real*8            :: rcycles
!
      !-------------------------------------------------------------------
      ! molecular weights must always be supplied for all selected species
      !-------------------------------------------------------------------

      do ic = 1, numSpecies
         if (self%mw(ic) == 0.0d0) then
            err_msg = 'mw problem in rc File.'
            call GmiPrintError (err_msg, .true., 2, 1, ic, 1, self%mw(ic), 0.0d0)
         end if
      end do

      if (do_semiss_inchem .and.  &
     &    ((emiss_opt == 0) .or. (self%chem_opt /= 2))) then
        err_msg = 'do_semiss_inchem problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, emiss_opt, self%chem_opt, 0, 0.0d0, 0.0d0)
      end if
!
      if (pr_sad .and. ((self%sad_opt == 0) .or. (self%num_sad == 0))) then
        err_msg = 'pr_sad/sad_opt/num_sad problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, self%sad_opt, self%num_sad, 0, 0.0d0, 0.0d0)
      end if

      if (pr_qj .and. ((self%phot_opt == 0) .or. (self%num_qjs == 0))) then
        err_msg = 'pr_qj/phot_opt/num_qj problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, self%phot_opt, self%num_qjs, 0, 0.0d0, 0.0d0)
      end if

      if (pr_qk .and. (self%num_qks == 0)) then
        err_msg = 'pr_qk/num_qk problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 1, self%num_qks, 0, 0, 0.0d0, 0.0d0)
      end if

      if (pr_qqjk) then
        if ((self%num_qjs == 0) .or. (self%num_qks == 0)) then
          err_msg = 'pr_qqjk/num_qj/k problem.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%num_qjs, self%num_qks, 0, 0.0d0, 0.0d0)
        end if

        if (self%phot_opt == 0) then
          err_msg = 'pr_qqjk/phot_opt problem.'
          call GmiPrintError  &
     &      (err_msg, .true., 1, self%phot_opt, 0, 0, 0.0d0, 0.0d0)
        end if

        if (.not. self%do_full_chem) then
          err_msg = 'pr_qqjk/chem_opt problem.'
          call GmiPrintError  &
     &      (err_msg, .true., 1, self%chem_opt, 0, 0, 0.0d0, 0.0d0)
        end if

        if (do_qqjk_inchem .and. do_mean) then
          err_msg = 'do_qqjk_inchem/do_mean problem.'
          call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
        end if
      end if

!     ==================
      if (self%chem_opt == 1) then
!     ==================
         if (numSpecies /= 2) then
            err_msg = 'The Radon/Lead chemistry requires two species.'
            call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
         end if
      end if

!     ==================
      if (self%chem_opt == 2) then
!     ==================

        if (numSpecies < IGAS) then
          err_msg = 'num_species/IGAS problem.'
          call GmiPrintError (err_msg, .true., 2, numSpecies, IGAS, 0, 0.0d0, 0.0d0)
        end if

        if (self%num_qjs /= IPHOT) then
          err_msg = 'num_qjs/IPHOT problem.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%num_qjs, IPHOT, 0, 0.0d0, 0.0d0)
        end if

        if (self%num_qks /= ITHERM) then
          err_msg = 'num_qks/ITHERM problem in.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%num_qks, ITHERM, 0, 0.0d0, 0.0d0)
        end if

      end if

!     ==================
      if (self%chem_opt == 4) then
!     ==================

        if (self%forc_bc_lz_val < 0.0d0) then
          err_msg = 'Impossible forc_bc_lz_val.'
          call GmiPrintError  &
     &      (err_msg, .true., 0, 0, 0, 1, self%forc_bc_lz_val, 0.0d0)
        end if

        if ((self%forc_bc_lz_val /= 0.0d0) .and. (numSpecies > 1)) then
          err_msg = 'Set_3dLoss code needs generalization for this case.'
          call GmiPrintError (err_msg, .true., 1, numSpecies, 0,  &
     &       1, self%forc_bc_lz_val, 0.0d0)
        end if

      end if

!     ==================
      if (self%chem_opt == 6) then
!     ==================

        if (numSpecies /= 2) then
          err_msg = 'The Beryllium chemistry requires two species.'
          call GmiPrintError  &
     &      (err_msg, .true., 1, numSpecies, 0, 0, 0.0d0, 0.0d0)
        end if

      end if

!     ==================
      if (self%chem_opt == 9) then
!     ==================
         if (numSpecies .ne. 18) then
            err_msg = 'The Tracer package CHEMCASE requires eighteen species.'
            call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
         end if
      end if

!     =================
      if (self%do_full_chem) then
!     =================

        if (self%phot_opt == 0) then
          err_msg = 'chem_opt/phot_opt problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%chem_opt, self%phot_opt, 0, 0.0d0, 0.0d0)
        end if

        if ((self%sad_opt == 0) .and. (met_opt /= 3)) then
          err_msg = 'sad_opt/met_opt problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%sad_opt, met_opt, 0, 0.0d0, 0.0d0)
        end if

        if ((self%sad_opt == 1) .or. (self%sad_opt == 2)) then

          if ((self%ich4_num == 0) .or. (self%ich4_num /= ICH4)) then
            err_msg = 'chem_opt/ich4_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, self%chem_opt, self%ich4_num, 0, 0.0d0, 0.0d0)
          end if

          if ((self%idehyd_num == 0) .or. (self%idehyd_num /= IDEHYD)) then
            err_msg = 'chem_opt/idehyd_num problem in the rc File.'
            call GmiPrintError  &
     &       (err_msg, .true., 2, self%chem_opt, self%idehyd_num, 0, 0.0d0, 0.0d0)
          end if

          if ((self%ih2o_num == 0) .or. (self%ih2o_num /= IH2O)) then
            err_msg = 'chem_opt/ih2o_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, self%chem_opt, self%ih2o_num, 0, 0.0d0, 0.0d0)
          end if

        end if

!.sds - the below if was in code twice in a row
        if ((self%imgas_num == 0) .or. (self%imgas_num /= IMGAS)) then
          err_msg = 'chem_opt/imgas_num problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%chem_opt, self%imgas_num, 0, 0.0d0, 0.0d0)
        end if

        if (self%do_synoz) then
          if (self%isynoz_num == 0) then
            err_msg = 'do_synoz/isynoz_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 1, self%isynoz_num, 0, 0, 0.0d0, 0.0d0)
          end if
        end if

        if (self%do_nodoz) then
          if (.not. self%do_synoz) then
             err_msg = 'do_nodoz/do_synoz problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
          end if

          if ((self%num_nox == 0) .or. (self%num_noy == 0)) then
            err_msg = 'num_nox/num_noy problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, self%num_nox, self%num_noy, 0, 0.0d0, 0.0d0)
          end if

          if (self%ihno3_num == 0) then
            err_msg = 'do_nodoz/ihno3_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 1, self%ihno3_num, 0, 0, 0.0d0, 0.0d0)
          end if
        end if

        if (((self%io3_num  == 0) .or. (self%io3_num  /= IO3)) .and.  &
     &      (self%phot_opt == 4)) then
          err_msg = 'io3_num/IO3 problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%io3_num, IO3, 0, 0.0d0, 0.0d0)
        end if

        if ((self%num_chem == 0) .or. (self%num_chem /= NCHEM)) then
          err_msg = 'chem_opt/num_chem problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%chem_opt, self%num_chem, 0, 0.0d0, 0.0d0)
        end if

        if ((self%num_molefrac == 0) .or. (self%num_molefrac /= NMF)) then
          err_msg = 'chem_opt/num_molefrac problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%chem_opt, self%num_molefrac,  &
     &       0, 0.0d0, 0.0d0)
        end if

        if ((self%num_qjs == 0) .or. (self%num_qjs /= NUM_J)) then
          err_msg = 'chem_opt/num_qjs problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%chem_opt, self%num_qjs, 0, 0.0d0, 0.0d0)
        end if

        if ((self%num_qks == 0) .or. (self%num_qks /= NUM_K)) then
          err_msg = 'chem_opt/num_qks problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%chem_opt, self%num_qks, 0, 0.0d0, 0.0d0)
        end if

        if (numSpecies /= NSP) then
          err_msg = 'chem_opt/num_species problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%chem_opt, numSpecies, 0, 0.0d0, 0.0d0)
        end if

        if (met_opt /= 3) then
          if ((self%num_sad == 0) .or. (self%num_sad /= NSAD)) then
            err_msg = 'met_opt/num_sad problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, met_opt, self%num_sad, 0, 0.0d0, 0.0d0)
          end if
        end if

        if ((self%chem_cycle <= 0.0d0) .or. (self%chem_cycle > 1.0d0)) then

          err_msg = 'chem_cycle must be > 0 .and. <= 1.'
          call GmiPrintError  &
     &      (err_msg, .true., 0, 0, 0, 1, self%chem_cycle, 0.0d0)

        else if (self%chem_cycle < 1.0d0) then  ! subcycle chemistry

          rcycles = 1.0d0 / self%chem_cycle
          icycles = Nint (rcycles)

          if ((Abs (rcycles - icycles)) > 0.000001d0) then
            err_msg =  &
     &    'if chem_cycle < 1, 1/chem_cycle must be a whole real number.'
            call GmiPrintError  &
     &        (err_msg, .true., 0, 0, 0, 1, self%chem_cycle, 0.0d0)
          end if

        end if

!        if (self%do_chem_grp) then
!
!          if (NUMGRP > MAX_NUM_CGRP) then
!            err_msg = 'NUMGRP/MAX_NUM_CGRP problem in the rc File.'
!            call GmiPrintError  &
!     &        (err_msg, .true., 2, NUMGRP, MAX_NUM_CGRP,  &
!     &         0, 0.0d0, 0.0d0)
!          end if
!
!          if (MAXGRP_ELEM > MAX_NUM_SMARRAY) then
!            err_msg =  &
!     &        'MAXGRP_ELEM/MAX_NUM_SMARRAY prob. in the rc File.'
!            call GmiPrintError  &
!     &        (err_msg, .true., 2, MAXGRP_ELEM, MAX_NUM_SMARRAY,  &
!     &         0, 0.0d0, 0.0d0)
!          end if
!
!          if (num_cgrp <= 0) then
!            err_msg = 'do_chem_grp/num_cgrp problem in the rc File.'
!            call GmiPrintError  &
!     &        (err_msg, .true., 1, num_cgrp, 0, 0, 0.0d0, 0.0d0)
!          end if
!
!          do ic = 1, num_cgrp
!            if (num_cgrp_elem(ic) <= 0) then
!              err_msg = 'num_cgrp/num_cgrp_elem prob.in the rc File.'
!              call GmiPrintError  &
!     &          (err_msg, .true., 2, ic, num_cgrp_elem(ic),  &
!     &           0, 0.0d0, 0.0d0)
!            end if
!          end do
!
!        end if

        if ((self%phot_opt  == 3) .and. (self%uvalbedo_opt == 0)) then
          err_msg = 'uvalbedo_opt problem in the rc File.'
          call GmiPrintError  &
     &      (err_msg, .true., 2, self%phot_opt, self%uvalbedo_opt, 0, 0.0d0,0.0d0)
        end if

!       ===================
        if (btest(emiss_opt,1)) then
!       ===================

          if ((self%iisoprene_num == 0) .and. (self%ino_num == 0)) then
            err_msg = 'isopropene/no problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
          end if

          if ((self%iisoprene_num == 0) .and. (self%iacetone_num /= 0)) then
            err_msg = 'isoprene/acetone problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 1, self%iacetone_num, 0, 0, 0.0d0, 0.0d0)
          end if

          if ((self%iisoprene_num == 0) .and. (self%ico_num /= 0)) then
            err_msg = 'isoprene/co problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 1, self%ico_num, 0, 0, 0.0d0, 0.0d0)
          end if

          if ((self%iisoprene_num == 0) .and. (self%ipropene_num /= 0)) then
            err_msg = 'isopropene/propene problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 1, self%ipropene_num, 0, 0, 0.0d0, 0.0d0)
          end if

!.sds.. mv'd do_full_chem test from here up ~26 lines
          if (self%iacetone_num  /= IC3H6O) then
            err_msg = 'iacetone_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, self%iacetone_num, IC3H6O, 0, 0.0d0,0.0d0)
          end if

          if (self%ico_num       /= ICO) then
            err_msg = 'ico_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, self%ico_num, ICO, 0, 0.0d0,0.0d0)
          end if

          if (self%iisoprene_num /= IC5H8) then
            err_msg = 'iisoprene_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, self%iisoprene_num, IC5H8, 0, 0.0d0,0.0d0)
          end if

          if (self%ipropene_num  /= IC3H6) then
            err_msg = 'ipropene_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, self%ipropene_num, IC3H6, 0, 0.0d0, 0.0d0)
          end if

          if (self%ino_num       /= INO) then
            err_msg = 'ino_num problem in the rc File.'
            call GmiPrintError  &
     &        (err_msg, .true., 2, self%ino_num, INO, 0, 0.0d0, 0.0d0)
          end if

        end if

      end if

      return

      end subroutine rcCheckChemistrySetting
!-----------------------------------------------------------------------------

  end module GmiChemistryMethod_mod
