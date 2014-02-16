!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiEmissionMethod_mod
!
! !INTERFACE:
!
#include "GmiESMF_ErrLog.h"
!
  module GmiEmissionMethod_mod
!
! !USES:
      use ESMF_Mod, only : ESMF_Config, ESMF_MAXSTR, ESMF_ConfigGetAttribute
      use GmiESMF_ErrorChecking_mod
      use GmiSpeciesRegistry_mod   , only : getSpeciesPosition
      use GmiESMFrcFileReading_mod, only : rcEsmfReadTable, rcEsmfReadLogical
      use GmiTimeControl_mod, only : GmiSplitDateTime, GetDaysFromJanuary1,    &
     &       t_GmiClock, Get_gmiTimeStep, Get_curGmiDate, Get_numTimeSteps,    &
     &       Get_curGmiTime
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration,        &
     &       Get_concentration, Set_concentration, Get_tracer_opt,             &
     &       Get_tr_source_land, Get_tr_source_ocean
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle, CleanArrayPointer
      use m_set_NLANDHAR            , only : NLANDHAR_expected, NLANDHAR
      use GmiPrintError_mod, only : GmiPrintError
      use GmiCheckRange_mod, only : CheckRange3d
      use GmiGrid_mod, only : t_gmiGrid, Get_i1, Get_i2, Get_ju1, Get_j2,      &
     &       Get_k1, Get_k2, Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl,      &
     &       Get_ilo, Get_ihi, Get_julo, Get_jhi, Get_ilo_gl, Get_ihi_gl,      &
     &       Get_julo_gl, Get_jhi_gl, Get_ilong, Get_ilat, Get_ivert,          &
     &       Get_numSpecies
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID, Get_mcor,&
     &       Get_latdeg, Get_londeg
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_pr_diag,         &
     &       Get_do_aerocom, Get_pr_const, Get_pr_emiss_3d, Get_pr_surf_emiss
      use ReadLightningData_mod, only : readLightningCoeffData
      use GmiMetFieldsControl_mod, only : t_metFields, Get_am, Get_bm, Get_pt, &
     &       Get_cmi_flags, Get_metdata_name_org, Get_metdata_name_model,      &
     &       Get_met_opt, Get_lwi_flags, Get_surfTemp15DayAvg, Get_mass,       &
     &       Get_gridBoxHeight, Get_press3e, Get_press3c, Get_met_opt, Get_kel,&
     &       Get_pctm2, Get_radswg,                                            &
     &       Get_surf_air_temp, Get_surf_rough, Get_con_precip, Get_tot_precip,&
     &       Get_ustar, Get_u10m, Get_v10m, Get_gwet, Get_pbl, Get_pardif,     &
     &       Get_pardir, Get_fracCloudCover, Get_dtrn, Get_zmmu, Get_cmf,      &
     &       Set_cmi_flags
      use GmiLandWaterIce_mod, only : computeGlobCMIflags1
!
!
  implicit none
!
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: InitializeEmission, initReadEmission
  public  :: RunEmission       , runReadEmission
  public  :: FinalizeEmission
!
  public  :: Set_emissionArray
  public  :: Set_emiss_isop  , Set_emiss_monot, Set_emiss_nox
  public  :: Set_emiss_o3    , Set_emiss_hno3
  public  :: Set_lightning_no, Set_flashrate, Set_cmi_flags1, Set_cldmas0
  public  :: Set_emissDust   , Set_emissAero
  public  :: Set_emissDust_t , Set_emissAero_t
  public  :: Set_isop_scale  , Set_aerosolSurfEmiss, Set_aerosolEmiss3D
  public  :: Set_emiss_3d_out, Set_surf_emiss_out  , Set_surf_emiss_out2
!
  public  :: Get_emissionArray       , Get_curEmissionFileRecord
  public  :: Get_emiss_isop          , Get_emiss_monot , Get_emiss_nox
  public  :: Get_emiss_o3            , Get_emiss_hno3
  public  :: Get_emiss_opt           , Get_emiss_in_opt
  public  :: Get_emiss_map           , Get_emiss_timpyr
  public  :: Get_lightning_no        , Get_flashrate, Get_cmi_flags1, Get_cldmas0
  public  :: Get_GCR_NOx
  public  :: Get_emissDust           , Get_emissAero
  public  :: Get_emissDust_t         , Get_emissAero_t
  public  :: Get_fertscal_infile_name, Get_light_infile_name
  public  :: Get_lai_infile_name     , Get_precip_infile_name
  public  :: Get_soil_infile_name    , Get_veg_infile_name
  public  :: Get_isopconv_infile_name, Get_monotconv_infile_name
  public  :: Get_lightCoeff_infile_name
  public  :: Get_lightning_opt       , Get_emiss_aero_opt
  public  :: Get_emiss_dust_opt      , Get_num_emiss,       Get_num_emiss3d
  public  :: Get_ndust               , Get_naero
  public  :: Get_nst_dust            , Get_nt_dust
  public  :: Get_emiss_map_dust      , Get_emiss_map_aero
  public  :: Get_isop_scale          , Get_do_gcr
  public  :: Get_do_semiss_inchem    , Get_do_ShipEmission, Get_doMEGANemission
  public  :: Get_o3_index            , Get_hno3_index
  public  :: Get_doReadDailyEmiss    , Get_desired_g_N_prod_rate
  public  :: Get_begDailyEmissRec    , Get_endDailyEmissRec
  public  :: Get_aerosolSurfEmiss    , Get_aerosolSurfEmissMap, Get_aerosolEmiss3D
  public  :: Get_emiss_3d_out        , Get_surf_emiss_out     , Get_surf_emiss_out2
  public  :: Get_ireg        , Get_iuse
  public  :: Get_iland       , Get_index_soil, Get_ncon_soil
  public  :: Get_soil_fert   , Get_soil_precip
  public  :: Get_xlai        , Get_xlai2
  public  :: Get_base_isop   , Get_base_monot
  public  :: Get_coeff_isop  , Get_convert_isop, Get_convert_monot
!
  public  :: readEmissionResourceFile
!
! !PUBLIC MEMBER DATA:
!
  public  :: t_Emission
!
# include "GmiParameters.h"
# include "gmi_emiss_constants.h"
# include "gmi_phys_constants.h"
# include "gmi_time_constants.h"
!
  type t_Emission
    integer             :: curEmissionFileRecord  ! current record number for the
                                                  ! main emission file.
    integer             :: semiss_inchem_flag
    logical             :: doScaleNOffEmiss    ! for fossil fuel scaling
    logical             :: doScaleNObbEmiss    ! for biomass burning scaling
    real*8 , pointer    :: scFacNOff(:,:,:) => null()
    real*8 , pointer    :: scFacNObb(:,:,:) => null()
    character (len=MAX_LENGTH_FILE_NAME) :: scFactorNOff_infile_name ! scale factor NO fossil fuel
                                                    ! emission input file name
    character (len=MAX_LENGTH_FILE_NAME) :: scFactorNObb_infile_name ! scale factor NO biomass burning
                                                    ! emission input file name
    character (len=MAX_LENGTH_FILE_NAME) :: emiss_infile_name        ! emission input file name
    character (len=10 ) :: emittedSpeciesNames (MAX_NUM_CONST) ! emission names of species read in
    character (len=MAX_LENGTH_VAR_NAME)  :: emiss_var_name                ! NetCDF emiss variable name
    integer             :: emiss_conv_flag               ! emission conversion flag
    real*8              :: emiss_conv_fac                ! emission conversion factor (s^-1)
    real*8              :: emiss_init_val                ! value to initialize emiss array to (kg/s)
    integer             :: emiss_aero_opt                ! sulfur aerosol emiss option
    integer             :: naero                         ! number of aerosols
    integer             :: emiss_map_aero(MAX_NUM_CONST) ! mapping of aerosol emiss number to const species #
    character (len=MAX_LENGTH_FILE_NAME) :: emiss_aero_infile_name        ! aerosol (sulf. code) input file name
    character (len=MAX_LENGTH_FILE_NAME) :: emiss_dust_infile_name        ! dust (sulf. code) input file name
    !... variables for calculation of the Galactic Cosmic Ray production of NOx (read in)
    logical         :: do_gcr                                ! do Galactic Cosmic Ray
    character (len=MAX_LENGTH_FILE_NAME) :: gcr_infile_name  ! Galactic Cosmic Ray input filename
    real*8          :: gcr_sunspot                           ! Galactic Cosmic Ray parameter
    real*8, pointer :: gcr_slope (:,:) => null()             ! Galactic Cosmic Ray parameter
    real*8, pointer :: gcr_aintcp (:,:) => null()            ! Galactic Cosmic Ray parameter
    real*8, pointer :: GCR_NOx (:,:,:) => null()               ! Galactic Cosmic Ray parameter
!
    character (len=MAX_LENGTH_FILE_NAME) :: fertscal_infile_name          ! fertilizer scale     input file name
    character (len=MAX_LENGTH_FILE_NAME) :: lai_infile_name               ! leaf area index      input file name
    character (len=MAX_LENGTH_FILE_NAME) :: light_infile_name             ! light                input file name
    character (len=MAX_LENGTH_FILE_NAME) :: precip_infile_name            ! precipitation        input file name
    character (len=MAX_LENGTH_FILE_NAME) :: soil_infile_name              ! soil type            input file name
    character (len=MAX_LENGTH_FILE_NAME) :: veg_infile_name               ! vegetation type      input file name
    character (len=MAX_LENGTH_FILE_NAME) :: isopconv_infile_name          ! isoprene convert     input file name
    character (len=MAX_LENGTH_FILE_NAME) :: monotconv_infile_name         ! monoterpene convert  input file name
    character (len=MAX_LENGTH_FILE_NAME) :: GOCARTerod_infile_name
    character (len=MAX_LENGTH_FILE_NAME) :: GOCARTocean_infile_name
    character (len=MAX_LENGTH_FILE_NAME) :: GOCARTerod_mod_infile_name
    real*8              :: isop_scale (12)               ! array of monthly isoprene scaling coefficients
!
    integer             :: lightning_opt                 ! lightning option
    integer             :: i_no_lgt                      ! lightning option
    integer             :: ik0                           ! level for 3-D convective
    real*8              :: desired_g_N_prod_rate         ! lightning option
    real*8 , pointer    :: emissAero   (:,:,:)   => null() ! used in sulfur chemistry
    real*8 , pointer    :: emissDust_t (:,:,:,:) => null() ! used in sulfur chemistry
    real*8 , pointer    :: emissAero_t (:,:,:,:) => null() ! used in sulfur chemistry
!
    integer          :: emiss_timpyr             ! emission times per year
                                                 ! (1 => yearly, 12 => monthly)
    integer          :: emiss_opt                ! emission option
    integer          :: emiss_in_opt             ! emission input option
    integer          :: emiss_map(MAX_NUM_CONST) ! mapping of emission number to const species #
    integer          :: emiss_map_dust(MAX_NUM_CONST) ! mapping of dust emiss number to const species #
    integer          :: ndust                     ! number of dust
    integer          :: nst_dust                  ! starting index for dust
    integer          :: nt_dust                   ! ending   index for dust
    logical          :: do_ShipEmission           ! do ship emissions?
    logical          :: doMEGANemission           ! do MEGAN emissions?
    character (len=MAX_LENGTH_FILE_NAME) :: laiMEGAN_InfileName    ! Input file name for AVHRR 
                                                  ! leaf-area-indices
    character (len=MAX_LENGTH_FILE_NAME) :: aefMboMEGAN_InfileName ! Annual emission factor for
                                                  ! methyl butenol input file name
    character (len=MAX_LENGTH_FILE_NAME) :: aefIsopMEGAN_InfileName ! Annual emission factor for
                                                   ! isoprene input file name
    character (len=MAX_LENGTH_FILE_NAME) :: aefMonotMEGAN_InfileName ! Annual emission factor for
                                                    ! monoterpenes input file name
    character (len=MAX_LENGTH_FILE_NAME) :: aefOvocMEGAN_InfileName ! Annual emission factor for other
                                                   ! biogenic VOCs input file name
    integer          :: days_btw_m     ! days between midmonths in the LAI data
    real*8 , pointer :: isoLai    (:,:) => null()  ! AVHRR LAI data for the current day
    real*8 , pointer :: isoLaiPrev(:,:) => null()  ! AVHRR LAI data for the previous month
    real*8 , pointer :: isoLaiCurr(:,:) => null()  ! AVHRR LAI data for the current  month
    real*8 , pointer :: isoLaiNext(:,:) => null()  ! AVHRR LAI data for the next     month
    real*8 , pointer :: aefMbo   (:,:) => null()   ! Annual emission factor for methyl butenol
    real*8 , pointer :: aefIsop  (:,:) => null()   ! Annual emission factor for isoprene
    real*8 , pointer :: aefOvoc  (:,:) => null()   ! Annual emission factor for other biogenic VOCs
    real*8 , pointer :: aefMonot (:,:) => null()   ! Annual emission factor for monoterpenes
    logical          :: do_semiss_inchem           ! do surface emissions inside chemistry solver?
    integer          :: emiss_dust_opt                ! sulfur dust emiss option
    integer          :: num_emiss                     ! total number of emission species
    integer          :: num_emiss3d                   ! for mixed 2d and 3d emissions, number of 3D emissions
    integer          :: o3_index   ! index of O3   in emiss_map
    integer          :: hno3_index ! index of HNO3 in emiss_map
    integer, pointer    :: cmi_flags1  (:,:) => null () ! used in lightning module
    real*8,  pointer    :: cldmas0     (:,:) => null ()
    real*8 , pointer :: flashrate      (:,:)   => null() ! used in the lightning module
    real*8 , pointer :: lightning_no   (:,:,:) => null() ! used in the lightning module
    real*8 , pointer :: localCoeff     (:,:) => null() ! local coefficient array for lightning
    real*8 , pointer :: midLatAdj (:,:) => null()            ! mid latitude adjustment array for lightning
    real*8           :: globalCoeff                    ! global coefficent scalar for lightning
    real*8           :: lightThreshold                 ! threshold for lightning
    character (len=MAX_LENGTH_FILE_NAME) :: lightCoeff_infile_name      ! file with local & global lightning coefficents
    integer          :: lightYearDim                   ! specified which yearly record to read lightning ratios from
    real*8 , pointer :: emiss          (:,:,:,:)   => null() ! array of emissions (kg/s)
    real*8 , pointer :: emissDust     (:,:,:)     => null() ! used in sulfur chemistry
    real*8 , pointer :: emiss_isop     (:,:)       => null() ! isoprene    emissions (kg/s)
    real*8 , pointer :: emiss_monot    (:,:)       => null() ! monoterpene emissions (kg/s)
    real*8 , pointer :: emiss_nox      (:,:)       => null() ! NOx         emissions (kg/s)
    real*8 , pointer :: emiss_o3    (:,:)       => null() ! ozone       emissions
    real*8 , pointer :: emiss_hno3     (:,:)       => null() ! hno3        emissions
    integer          :: ncon_soil (NVEGTYPE)           ! Olson -> soil type
                                                       !    1 => water/desert/ice
                                                       !    2 => tropical rain forest
                                                       !    3 => conifers
                                                       !    4 => dry deciduous
                                                       !    5 => other deciduous
                                                       !    6 => woodland
                                                       !    7 => grassland
                                                       !    8 => agriculture (other than rice)
                                                       !    9 => rice paddies
                                                       !   10 => wetland/tundra
    integer, pointer :: index_soil  (:,:)   => null()
    integer, pointer :: ireg        (:,:)   => null()  ! number of land types in a grid square
    integer, pointer :: iuse        (:,:,:) => null()  ! fraction of grid box area occupied by land type
    integer, pointer :: iland       (:,:,:) => null()  ! land type id in grid square for ireg land types
    real*8           :: coeff_isop   (NPOLY)           ! coefficients used for polynomial fit
    real*8           :: convert_isop (NVEGTYPE)        ! isoprene emissions by landtype (atomsC/cm^2leaf/s)
    real*8           :: convert_monot(NVEGTYPE)        ! monoterpene emissions by landtype (atomsC/cm^2leaf/s)
    real*8 , pointer :: soil_fert   (:)     => null()  ! fertilizers (ng N/m^2/s)
    real*8 , pointer :: soil_precip (:,:)   => null()  ! two months of observed precip (mm/day/box)
    real*8 , pointer :: soil_pulse  (:,:)   => null()  ! tracking of wet/dry & three types of pulsing (Y&L, 94)
    real*8 , pointer :: base_isop   (:,:,:) => null()  ! baseline emissions for isoprene     (kgC/box/step?)
    real*8 , pointer :: base_monot  (:,:,:) => null()  ! baseline emissions for monoterpenes (kgC/box/step?)
    real*8 , pointer :: xlai        (:,:,:) => null()  ! leaf area index of land type for month #1
    real*8 , pointer :: xlai2       (:,:,:) => null()  ! leaf area index of land type for month #2
    type(t_GmiArrayBundle), pointer :: emissionArray(:)
!
    logical       :: doReadDailyEmiss
    integer       :: begDailyEmissRec, endDailyEmissRec
!
    real*8 , pointer    :: surf_emiss_out     (:,:,:)     => null()
    real*8 , pointer    :: surf_emiss_out2    (:,:,:)     => null()
    real*8 , pointer    :: emiss_3d_out       (:,:,:,:)   => null()
    real*8 , pointer    :: aerosolEmiss3D     (:,:,:,:)   => null()
    real*8 , pointer    :: aerosolSurfEmiss   (:,:,:)     => null()
                   ! Array storing the surface emission diagnostics for
                   ! aerosols. FSO2, NSO2, DMS, NSO4A, FSO4A, aerosol
                   ! (OC, BC, sea salt) and dust.
    integer, pointer    :: aerosolSurfEmissMap(:)     => null()
                   ! mapping of SurfEmiss number for aerosols to const species #
!
    integer :: soil_month         ! Determine if a new month
    logical :: firstReadEmiss
    logical :: firstReadLaiData
    integer :: lai_day_save
    integer :: soil_day_save
    integer :: currMonthMEGAN
    integer :: currDateEmiss
    integer :: dailyRecEmiss
    integer :: currHourEmiss
    integer :: newRecord

    logical :: firstScaleEmissions
      integer :: iNO_ff   
      integer :: iNO_bb   
      integer :: iCO_bb   
      integer :: iMEK_bb  
      integer :: iPRPE_bb 
      integer :: iC2H6_bb 
      integer :: iC3H8_bb 
      integer :: iALK4_bb 
      integer :: iALD2_bb 
      integer :: iCH2O_bb 
      real*8 , allocatable :: ffNO  (:,:,:)
      real*8 , allocatable :: bbNO  (:,:,:)
      real*8 , allocatable :: bbCO  (:,:,:)
      real*8 , allocatable :: bbMEK (:,:,:)
      real*8 , allocatable :: bbPRPE(:,:,:)
      real*8 , allocatable :: bbC2H6(:,:,:)
      real*8 , allocatable :: bbC3H8(:,:,:)
      real*8 , allocatable :: bbALK4(:,:,:)
      real*8 , allocatable :: bbALD2(:,:,:)
      real*8 , allocatable :: bbCH2O(:,:,:)
  end type t_Emission
!
! !DESCRIPTION:
!
! !AUTHOR:
!
! !REVISION HISTORY:
!
!EOP
!-------------------------------------------------------------------------
!
  CONTAINS
!
!-------------------------------------------------------------------------
!BOP
      subroutine readEmissionResourceFile(self, gmiGrid, gmiDomain, &
     &               Diagnostics, config)
!
      use GmiCheckNamelistFile_mod, only : CheckNamelistOptionRange
      use GmiSpeciesRegistry_mod, only : getSpeciesIndex, UNKNOWN_SPECIES
      use GmiStringManipulation_mod, only : constructListNames
!
      implicit none
!
!
! !INPUT PARAMETERS:
      type (t_gmiGrid    ), intent(in) :: gmiGrid
      type (t_gmiDomain  ), intent(in) :: gmiDomain
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Emission), intent(inOut) :: self
      type (ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in Emission related variables from the resource file.
!
! !LOCAL VARIABLES:
      integer :: numSpecies, procID, STATUS, RC, ic, chem_opt
      integer :: ju1_gl, j2_gl, i1, i2, ju1, j2, k1, k2
      logical :: pr_diag
      character (len=MAX_LENGTH_SPECIES_NAME) :: tempListNames(MAX_NUM_CONST)
      character (len=MAX_STRING_LENGTH      ) :: emissionSpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: emission3DSpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: emissionDustSpeciesNames
      character (len=MAX_STRING_LENGTH      ) :: emissionAeroSpeciesNames
      character(len=ESMF_MAXSTR) :: IAm, err_msg
!
!EOP
!--------------------------------------------------------------------
!BOC
      IAm = "readEmissionResourceFile"
!
      call Get_procID (gmiDomain, procID)
      call Get_pr_diag (Diagnostics, pr_diag)
!
      if (pr_diag) Write(6,*) IAm, 'called by ', procID
!
      call Get_numSpecies(gmiGrid, numSpecies)
!
!      allocate(tempListNames(numSpecies))
!
      !################################
      ! Begin reading the resource file
      !################################
!
      call ESMF_ConfigGetAttribute(config, chem_opt, &
     &                label   = "chem_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
    ! ---------------------------------
    ! emiss_opt
    !   (set up with bit switches, can turn on more than one option by adding
    !    the numbers togeether, ie 3 will turn on llnl and harvard emissions)
    !   0:  no emissions
    !   1:  do LLNL emissions
    !   2:  do Harvard emissions
    !   3:  do LLNL and Harvard emissions
    !   4:  do GSFC emissions (Galactic Cosmic Rays only right now)
    !   5:  do LLNL and GSFC emissions
    !   6:  do Harvard and GSFC emissions
    !   7:  do LLNL and Harvard and GSFC emissions
    ! ---------------------------------
!
      call ESMF_ConfigGetAttribute(config, self%emiss_opt, &
     &                label   = "emiss_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
    ! --------------------------------------------
    ! emiss_in_opt
    !   0:  no emissions data
    !   1:  set all emiss values to emiss_init_val
    !   2:  read in emiss values
    !   3:  set in emiss values using constant source terms
    ! --------------------------------------------
!
      call ESMF_ConfigGetAttribute(config, self%emiss_in_opt, &
     &                label   = "emiss_in_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
    ! ------------------------------------
    ! emiss_conv_flag
    !   0:  no conversion performed
    !   1:  use emiss_conv_fac
    !   2:  convert from kg/km2-hr to kg/s
    ! ------------------------------------
!
      call ESMF_ConfigGetAttribute(config, self%emiss_conv_flag, &
     &                label   = "emiss_conv_flag:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
    ! ------------------------------------------------------------------
    ! semiss_inchem_flag
    !   <0:  if emissions are on, surface emissions will be done in
    !        Smvgear chemistry if it is on; outside of chemistry if
    !        Smvgear chemistry is off (i.e., "auto set")
    !    0:  if emissions are on, surface emissions will be done outside
    !        of chemistry
    !   >0:  if emissions are on, surface emissions will be done in
    !        Smvgear chemistry
    ! ------------------------------------------------------------------
!
      call ESMF_ConfigGetAttribute(config, self%semiss_inchem_flag, &
     &                label   = "semiss_inchem_flag:", &
     &                default = -1, rc=STATUS )
      VERIFY_(STATUS)
!
      ! sets of emissons per year (1 => yearly, 12 => monthly)
!
      call ESMF_ConfigGetAttribute(config, self%emiss_timpyr, &
     &                label   = "emiss_timpyr:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)
!
      self%emiss_map(:) =  0
!
      call rcEsmfReadTable(config, emissionSpeciesNames, &
     &                     "emissionSpeciesNames::", rc=STATUS)
!
      call rcEsmfReadTable(config, emission3DSpeciesNames, &
     &                     "emission3DSpeciesNames::", rc=STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%emiss_conv_fac, &
     &                label   = "emiss_conv_fac:", &
     &                default = 1.0d0, rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%emiss_init_val, &
     &                label   = "emiss_init_val:", &
     &                default = 1.0d0, rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%emiss_infile_name, &
     &                label   = "emiss_infile_name:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%emiss_var_name, &
     &                label   = "emiss_var_name:", &
     &                default = 'emiss', rc=STATUS )
      VERIFY_(STATUS)
!
    !---------------------------------------
    ! Reading of daily emission file options
    !---------------------------------------
!
      call rcEsmfReadLogical(config, self%doReadDailyEmiss, &
     &           "doReadDailyEmiss:", default=.false., rc=STATUS)
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%begDailyEmissRec, &
     &                label   = "begDailyEmissRec:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%endDailyEmissRec, &
     &                label   = "endDailyEmissRec:", &
     &                default = 366, rc=STATUS )
      VERIFY_(STATUS)
!
    ! --------------------------------------
    ! Aerosols and Sulfur from Penner et al.
    ! --------------------------------------
!
     ! emiss_aero_opt   0: for no     aerosol emissions
     !                  1: for GMI    aerosol emissions
     !                  2: for GOCART aerosol emissions
!
      call ESMF_ConfigGetAttribute(config, self%emiss_aero_opt, &
     &                label   = "emiss_aero_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
      ! number of aerosol emissions
!
      call ESMF_ConfigGetAttribute(config, self%naero, &
     &                label   = "naero:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)
!
      self%emiss_map_aero(:) = 0    ! map emissions to species #
!
      call rcEsmfReadTable(config, emissionAeroSpeciesNames, &
     &                     "emissionAeroSpeciesNames::", rc=STATUS)
!
      ! aerosol emission input file
!
      call ESMF_ConfigGetAttribute(config, self%emiss_aero_infile_name, &
     &                label   = "emiss_aero_infile_name:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
     ! emiss_dust_opt   0: for no     dust emissions
     !                  1: for GMI    dust emissions
     !                  2: for GOCART dust emissions
!
      call ESMF_ConfigGetAttribute(config, self%emiss_dust_opt, &
     &                label   = "emiss_dust_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
      ! number of dust bins
!
      call ESMF_ConfigGetAttribute(config, self%ndust, &
     &                label   = "ndust:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)
!
     ! number of starting point for new dust emiss.
!
      call ESMF_ConfigGetAttribute(config, self%nst_dust, &
     &                label   = "nst_dust:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)
!
      ! number of dust emissions per emiss. file
!
      call ESMF_ConfigGetAttribute(config, self%nt_dust, &
     &                label   = "nt_dust:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)
!
      self%emiss_map_dust(:) = 0    ! map emissions to species #
!
      emissionDustSpeciesNames = ""
!
      call rcEsmfReadTable(config, emissionDustSpeciesNames, &
     &                     "emissionDustSpeciesNames::", rc=STATUS)
!
      ! dust emission input file
!
      call ESMF_ConfigGetAttribute(config, self%emiss_dust_infile_name, &
     &                label   = "emiss_dust_infile_name:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
     ! -----------------------------------------------------
     ! Information on scale factors for various NO emissions
     ! -----------------------------------------------------
!
     ! Fossil fuel
      call rcEsmfReadLogical(config, self%doScaleNOffEmiss, &
     &           "doScaleNOffEmiss:", default=.false., rc=STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%scFactorNOff_infile_name, &
     &                label   = "scFactorNOff_infile_name:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
     ! Biomass burning
      call rcEsmfReadLogical(config, self%doScaleNObbEmiss, &
     &           "doScaleNObbEmiss:", default=.false., rc=STATUS)
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%scFactorNObb_infile_name, &
     &                label   = "scFactorNObb_infile_name:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
    ! ---------------------------------------------------
    ! Harvard emissions:  acetone, isoprene, propene, NO.
    ! ---------------------------------------------------
!
      self%isop_scale(:) = 1.0d0
!
      call rcEsmfReadTable(config, self%isop_scale, &
     &                     "isop_scale::", rc=STATUS)
!
    !     --------------------------------------
    !     Aerosols and Sulfur from GOCART
    !     --------------------------------------
!
      call ESMF_ConfigGetAttribute(config, self%GOCARTerod_infile_name, &
     &                label   = "GOCARTerod_infile_name:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%GOCARTocean_infile_name, &
     &                label   = "GOCARTocean_infile_name:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%GOCARTerod_mod_infile_name, &
     &                label   = "GOCARTerod_mod_infile_name:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
    !... do Galactic Cosmic Rays source of NOx?
!
      ! turn on Galactic Cosmic ray emission of N and NO
!
      call rcEsmfReadLogical(config, self%do_gcr, &
     &           "do_gcr:", default=.false., rc=STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%gcr_infile_name, &
     &                label   = "gcr_infile_name:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
    !... do Ship Emission Calculations
!
      call rcEsmfReadLogical(config, self%do_ShipEmission, &
     &           "do_ShipEmission:", default=.false., rc=STATUS)
!
    !... do MEGAN emission calculations
!
      call rcEsmfReadLogical(config, self%doMEGANemission, &
     &           "doMEGANemission:", default=.false., rc=STATUS)
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%laiMEGAN_InfileName, &
     &                label   = "laiMEGAN_InfileName:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%aefMboMEGAN_InfileName, &
     &                label   = "aefMboMEGAN_InfileName:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%aefIsopMEGAN_InfileName, &
     &                label   = "aefIsopMEGAN_InfileName:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%aefOvocMEGAN_InfileName, &
     &                label   = "aefOvocMEGAN_InfileName:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%aefMonotMEGAN_InfileName, &
     &                label   = "aefMonotMEGAN_InfileName:", &
     &                default = ' ', rc=STATUS )
      VERIFY_(STATUS)
!
!     ------------------------------------------------------------
!     Note that if ((emiss_opt == 2) && do_full_chem), the indices
!     below will be automatically set by the setkin files.
!     ------------------------------------------------------------
!
!      call ESMF_ConfigGetAttribute(config, self%iacetone_num, &
!     &                label   = "iacetone_num:", &
!     &                default = 0, rc=STATUS )
!      VERIFY_(STATUS)
!
!      call ESMF_ConfigGetAttribute(config, self%ico_num, &
!     &                label   = "ico_num:", &
!     &                default = 0, rc=STATUS )
!      VERIFY_(STATUS)
!
!      call ESMF_ConfigGetAttribute(config, self%iisoprene_num, &
!     &                label   = "iisoprene_num:", &
!     &                default = 0, rc=STATUS )
!      VERIFY_(STATUS)
!
!      call ESMF_ConfigGetAttribute(config, self%ipropene_num, &
!     &                label   = "ipropene_num:", &
!     &                default = 0, rc=STATUS )
!      VERIFY_(STATUS)
!
!      call ESMF_ConfigGetAttribute(config, self%ino_num, &
!     &                label   = "ino_num:", &
!     &                default = 0, rc=STATUS )
!      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%fertscal_infile_name, &
     &                label   = "fertscal_infile_name:", &
     &                default = 'fertscale_4x5_dao.asc', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%lai_infile_name, &
     &                label   = "lai_infile_name:", &
     &                default = 'lai_4x5_dao.asc', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%light_infile_name, &
     &                label   = "light_infile_name:", &
     &                default = 'lighttable.asc', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%precip_infile_name, &
     &                label   = "precip_infile_name:", &
     &                default = 'precip_4x5_dao.asc', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%soil_infile_name, &
     &                label   = "soil_infile_name:", &
     &                default = 'soiltype.asc', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%veg_infile_name, &
     &                label   = "veg_infile_name:", &
     &                default = 'vegtype_4x5_dao.asc', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%isopconv_infile_name, &
     &                label   = "isopconv_infile_name:", &
     &                default = 'isopconvtable.asc', rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%monotconv_infile_name, &
     &                label   = "monotconv_infile_name:", &
     &                default = 'monotconvtable.asc', rc=STATUS )
      VERIFY_(STATUS)
!
    ! -----------------------------------------------------
    ! lightning_opt = 0 --> default lightning
    !               = 1 --> variable_lightning (Dale Allen's algorithm)
    !               = 2 --> no lightning
    !               = in options 1 and 2, NO_lgt (in gmi_emiss.F) is made zero
    !
    ! i_no_lgt = 0  --> set this index to location of NO_lgt in emiss infile
    !               --> 0 indicates not used, > 0 indicates actual location
    !                                             of NO_lgt in emiss infile
    !
    ! desired_g_N_prod_rate = 5 --> desired global Nitrogen production Rate.
    !                           --> default value ( 5 Tg. of N)
    ! -----------------------------------------------------
!
      call ESMF_ConfigGetAttribute(config, self%lightning_opt, &
     &                label   = "lightning_opt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%i_no_lgt, &
     &                label   = "i_no_lgt:", &
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%desired_g_N_prod_rate, &
     &                label   = "desired_g_N_prod_rate:", &
     &                default = 5.0d0, rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%lightCoeff_infile_name, &
     &                label   = "lightCoeff_infile_name:", &
     &                default = "FlashRatio.nc", rc=STATUS )
      VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%lightYearDim, &
     &                label   = "lightYearDim:", &
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)
!
      ! ---------------------------------------------------------------
      ! Check option ranges.  Note that as new options are added, these
      ! range checks will have to be modified.
      ! ---------------------------------------------------------------
!
      call CheckNamelistOptionRange ('emiss_opt'     , self%emiss_opt     , 0, 7)
      call CheckNamelistOptionRange ('emiss_aero_opt', self%emiss_aero_opt, 0, 2)
      call CheckNamelistOptionRange ('emiss_dust_opt', self%emiss_dust_opt, 0, 2)
      call CheckNamelistOptionRange ('emiss_in_opt'  , self%emiss_in_opt  , 0, 3)
      call CheckNamelistOptionRange ('lightning_opt' , self%lightning_opt , 0, 2)
!
      if (self%emiss_opt == 0) then
         self%emiss_in_opt    = 0
         self%emiss_conv_flag = 0
         self%emiss_map(:)    = 0
         self%num_emiss = 0
      else
!
         ! Set the initial value of the list
         tempListNames(:) = ''
!
         ! Construct the list of names using the long string
         call constructListNames(tempListNames, emissionSpeciesNames)
!
         self%num_emiss = Count (tempListNames(:) /= '')
!
         if (self%num_emiss > 0) then
            do ic = 1, self%num_emiss
               self%emiss_map(ic) = getSpeciesIndex(tempListNames(ic))
            end do
         end if
!
!.. fold in 3D names if they exist
         tempListNames(:) = ''
         call constructListNames(tempListNames, emission3DSpeciesNames)
!
         self%num_emiss3d = Count (tempListNames(:) /= '')
!
         if (self%num_emiss3d > 0) then
            do ic = 1, self%num_emiss3d
               self%emiss_map(self%num_emiss+ic) = getSpeciesIndex(tempListNames(ic))
            end do
         end if
         self%num_emiss = self%num_emiss + self%num_emiss3d
!
         if (self%emiss_dust_opt > 0) then
            ! Set the initial value of the list
            tempListNames(:) = ''
!
            ! Construct the list of names using the long string
            call constructListNames(tempListNames, emissionDustSpeciesNames)
!
            do ic = 1, self%ndust
               self%emiss_map_dust(ic) = getSpeciesIndex(tempListNames(ic))
            end do
         end if
!
         if (self%emiss_aero_opt > 0) then
            ! Set the initial value of the list
            tempListNames(:) = ''
!
            ! Construct the list of names using the long string
            call constructListNames(tempListNames, emissionAeroSpeciesNames)
!
            do ic = 1, self%naero
               self%emiss_map_aero(ic) = getSpeciesIndex(tempListNames(ic))
            end do
         end if
!
      endif
!
      if (self%semiss_inchem_flag < 0) then  ! "auto set"
         if ((self%emiss_opt /= 0) .and. (chem_opt == 2)) then
            self%do_semiss_inchem = .true.
         else
            self%do_semiss_inchem = .false.
         end if
      else if (self%semiss_inchem_flag == 0) then
         self%do_semiss_inchem   = .false.
      else if (self%semiss_inchem_flag  > 0) then
         self%do_semiss_inchem   = .true.
      end if
!
      !###############
      ! Error Checking
      !###############
!
      if ((numSpecies > MAX_NUM_CONST) .or.  &
     &    (self%num_emiss   > MAX_NUM_CONST)) then
         err_msg = 'MAX_NUM_CONST problem in rc File.'
         call GmiPrintError (err_msg, .true., 2, numSpecies, &
     &           self%num_emiss, 0, 0.0d0, 0.0d0)
      end if
!
      if ((self%emiss_opt /= 0) .and. (self%num_emiss == 0)) then
        err_msg = 'emiss_opt/num_emiss problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, self%emiss_opt, self%num_emiss, 0, 0.0d0, 0.0d0)
      end if
!
      if ((self%emiss_in_opt == 0) .and. (self%emiss_opt /= 0)) then
        err_msg = 'emiss_in_opt/emiss_opt problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 2, self%emiss_in_opt, self%emiss_opt, 0, 0.0d0, 0.0d0)
      end if
!
      if ((self%emiss_timpyr /= 1) .and.  &
     &    (self%emiss_timpyr /= MONTHS_PER_YEAR)) then
        err_msg = 'emiss_timpyr range problem.'
        call GmiPrintError  &
     &    (err_msg, .true., 1, self%emiss_timpyr, 0, 0, 0.0d0, 0.0d0)
      end if
!
      if ((self%emiss_in_opt /= 2) .and.  &
     &    (self%emiss_timpyr == MONTHS_PER_YEAR)) then
        err_msg = 'emiss_in_opt/emiss_timpyr problem.'
        call GmiPrintError (err_msg, .true., 2, self%emiss_in_opt,             &
     &                      self%emiss_timpyr, 0, 0.0d0, 0.0d0)
      end if
!
  return
!
  end subroutine readEmissionResourceFile
!
!-------------------------------------------------------------------------
!BOP
!
  subroutine runReadEmission (self, gmiClock, gmiGrid, gmiDomain, &
                        Diagnostics, do_drydep)
!
  use ReadEmissionFiles_mod    , only : ReadEmiss
  use ReadVegLaiData_mod       , only : readLeafAreaIndexData
  use ReadOtherEmissionData_mod, only : readPrecipitationData
  use ReadInputMEGAN_mod       , only : setMEGANisoLAI
!
  implicit none
!
! !INPUT PARAMETERS:
  type (t_gmiGrid  ), intent(in) :: gmiGrid
  type (t_GmiClock ), intent(in) :: gmiClock
  type (t_GmiDomain), intent(in) :: gmiDomain
  type (t_Diagnostics), intent(in) :: Diagnostics
  logical           , intent(in) :: do_drydep
!
! !INPUT/OUTPUT PARAMETERS:
  type (t_Emission)  , intent(inOut) :: self
!
! !DESCRIPTION:
! This routines reads in (daily or monthly) emission related files.
!
! !LOCAL VARIABLES:
  character (len=75) :: err_msg
  integer            :: il, ij, num_emiss, ic
  real*8             :: units_fac
  integer       :: nhms, nymd, ydummy, thisDay, thisMonth, thisDate, ddummy
  integer       :: i1, i2, ju1, j2, k1, k2, numSpecies, procID
  integer       :: i1_gl, i2_gl, ju1_gl, j2_gl, ilong, ilat, ivert
  integer       :: thisHour, ik
  logical       :: pr_diag
  real*8 , allocatable :: mcor(:,:)
!
!EOP
!----------------------------------------------------------------------
!BOC
  call Get_procID(gmiDomain, procID)
  call Get_pr_diag(Diagnostics, pr_diag)
!
  if (pr_diag) then
     Write (6,*) 'runReadEmission called by ', procID
  end if
!
  call Get_i1    (gmiGrid, i1)
  call Get_i2    (gmiGrid, i2)
  call Get_ju1   (gmiGrid, ju1)
  call Get_j2    (gmiGrid, j2)
  call Get_k1    (gmiGrid, k1)
  call Get_k2    (gmiGrid, k2)
  call Get_i1_gl (gmiGrid, i1_gl)
  call Get_i2_gl (gmiGrid, i2_gl)
  call Get_ju1_gl(gmiGrid, ju1_gl)
  call Get_j2_gl (gmiGrid, j2_gl)
  call Get_ilong (gmiGrid, ilong )
  call Get_ilat  (gmiGrid, ilat  )
  call Get_ivert (gmiGrid, ivert )
  call Get_numSpecies (gmiGrid, numSpecies )
!
      allocate(mcor(i1:i2,ju1:j2))
      call Get_mcor(gmiDomain, mcor)
!
  ! Obtain model time information
!
  call Get_curGmiDate  (gmiClock, nymd)
  call Get_curGmiTime  (gmiClock, nhms)
!
  num_emiss = self%num_emiss
!
  ! Read in the a new emission record at the beginning of each day (daily file)
  ! or the beginning of each month (monthly file)
!
  if (self%emiss_in_opt == 2) then
    ! Determine the current day and month
     call GmiSplitDateTime (nymd, ydummy, thisMonth, thisDay)
!
     if ((self%emiss_timpyr == MONTHS_PER_YEAR) .and. (.not. self%doReadDailyEmiss)) then
         thisDate = thisMonth
     elseif (self%doReadDailyEmiss)  then
         thisDate = thisDay
     end if
!
     ! Read in a new record if at the beginning of a day or a month
!
     if (self%currDateEmiss /= thisDate) then
        self%currDateEmiss = thisDate
        if (self%doReadDailyEmiss) then
           self%dailyRecEmiss  = self%dailyRecEmiss + 1
           self%curEmissionFileRecord = self%begDailyEmissRec + self%dailyRecEmiss
           if (self%curEmissionFileRecord > self%endDailyEmissRec) then
              self%dailyRecEmiss  = 0
              self%curEmissionFileRecord = self%begDailyEmissRec
           end if
        else
           self%curEmissionFileRecord = thisDate
        end if
!
        call ReadEmiss (self%emissionArray, self%emiss_infile_name, &
                    self%emiss_var_name, self%emiss_map, self%emittedSpeciesNames, &
                    self%num_emiss, self%num_emiss3d, self%lightning_opt, &
                    self%i_no_lgt, &
                    self%curEmissionFileRecord, i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, &
                    ilong, ilat, ivert, numSpecies, self%firstReadEmiss)
!
        if (self%lightning_opt == 1) then
           call readLightningCoeffData (self%lightcoeff_infile_name, &
                self%globalCoeff, self%localCoeff, self%midLatAdj, &
		self%lightThreshold, &
                self%ik0, &
                pr_diag, procID, thisMonth, self%lightYearDim, &
                i1, i2, ju1, j2, i1_gl, ju1_gl)
        end if
!
        ! Convert emiss, if necessary.
!
        if (self%emiss_conv_flag /= 0) then
           if (self%emiss_conv_flag == 1) then
              do ic = 1, num_emiss
                 self%emissionArray(ic)%pArray3D(:,:,:) = &
                      self%emissionArray(ic)%pArray3D(:,:,:) * self%emiss_conv_fac
              end do
           else if (self%emiss_conv_flag == 2) then
              units_fac = (1.0d0 / SECPHR) * (KMPM * KMPM)
!
              do ic = 1, num_emiss
                 do ij = ju1, j2
                    do il = i1, i2
                       self%emissionArray(ic)%pArray3D(il,ij,:) =  &
           &                self%emissionArray(ic)%pArray3D(il,ij,:) * mcor(il,ij) * units_fac
                    end do
                 end do
              end do
           end if
        end if
!
        do ic = 1, num_emiss
           call CheckRange3d ('emissionArray', procID, i1, i2, ju1, j2, k1, k2, &
                  self%emissionArray(ic)%pArray3D(:,:,:), -1.0d20, 1.0d20)
        end do
!
     end if
!
     ! ---------------------------------------------------
     ! Begin - Scale the NOx emissions on an hourly basis.
     ! ---------------------------------------------------
!
      if (self%doScaleNOffEmiss .or. self%doScaleNObbEmiss) then
         if (self.firstScaleEmissions) then
            self.firstScaleEmissions = .FALSE.

            self%iNO_ff   = -999
            self%iNO_bb   = -999
            self%iCO_bb   = -999
            self%iMEK_bb  = -999
            self%iPRPE_bb = -999
            self%iC2H6_bb = -999
            self%iC3H8_bb = -999
            self%iALK4_bb = -999
            self%iALD2_bb = -999
            self%iCH2O_bb = -999

            if (self%doScaleNOffEmiss) then
               allocate(self%ffNO(i1:i2,ju1:j2,k1:k2))
               self%iNO_ff = getSpeciesPosition('NO_ff', self%emittedSpeciesNames, self%num_emiss)
            end if
!
            if (self%doScaleNObbEmiss) then
               allocate(self%bbNO(i1:i2,ju1:j2,k1:k2))
               self%iNO_bb = getSpeciesPosition('NO_bb', self%emittedSpeciesNames, self%num_emiss)
!
               allocate(self%bbCO(i1:i2,ju1:j2,k1:k2))
               self%iCO_bb = getSpeciesPosition('CO_bb', self%emittedSpeciesNames, self%num_emiss)
!
               allocate(self%bbMEK(i1:i2,ju1:j2,k1:k2))
               self%iMEK_bb = getSpeciesPosition('MEK_bb', self%emittedSpeciesNames, self%num_emiss)
!
               allocate(self%bbPRPE(i1:i2,ju1:j2,k1:k2))
               self%iPRPE_bb = getSpeciesPosition('PRPE_bb', self%emittedSpeciesNames, self%num_emiss)
!
               allocate(self%bbC2H6(i1:i2,ju1:j2,k1:k2))
               self%iC2H6_bb = getSpeciesPosition('C2H6_bb', self%emittedSpeciesNames, self%num_emiss)
!
               allocate(self%bbC3H8(i1:i2,ju1:j2,k1:k2))
               self%iC3H8_bb = getSpeciesPosition('C3H8_bb', self%emittedSpeciesNames, self%num_emiss)
!
               allocate(self%bbALK4(i1:i2,ju1:j2,k1:k2))
               self%iALK4_bb = getSpeciesPosition('ALK4_bb', self%emittedSpeciesNames, self%num_emiss)
!
               allocate(self%bbALD2(i1:i2,ju1:j2,k1:k2))
               self%iALD2_bb = getSpeciesPosition('ALD2_bb', self%emittedSpeciesNames, self%num_emiss)
!
               allocate(self%bbCH2O(i1:i2,ju1:j2,k1:k2))
               self%iCH2O_bb = getSpeciesPosition('CH2O_bb', self%emittedSpeciesNames, self%num_emiss)
            end if

         end if
         call scaleEmissions(self, nhms, i1, i2, ju1, j2, k1, k2)
      endif
!
     ! -------------------------------------------------
     ! End - Scale the NOx emissions on an hourly basis.
     ! -------------------------------------------------
!
  end if
!
  if (do_drydep .or. (btest(self%emiss_opt,1))) then
     call readLeafAreaIndexData &
              (self%lai_infile_name, self%ireg, self%xlai, self%xlai2, &
               pr_diag, procID, nymd, i1, i2, ju1, j2, &
               self%lai_day_save, self%firstReadLaiData)
!
     if (btest(self%emiss_opt,1)) then
        call readPrecipitationData (self%precip_infile_name, &
                 self%index_soil, self%soil_precip, nymd,    &
                 self%soil_month, pr_diag, procID)
     end if
  end if
!
  if (self%doMEGANemission) then
     call setMEGANisoLAI (self%isoLai, self%isoLaiCurr, self%isoLaiPrev, &
     &          self%isoLaiNext, self%days_btw_m, self%laiMEGAN_InfileName, &
     &          nymd, i1, i2, ju1, j2, i1_gl, ju1_gl, self%currMonthMEGAN)
  end if
!
  return
!
  end subroutine runReadEmission
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: scaleEmissions
!
! !INTERFACE:
!
      subroutine scaleEmissions(self, nhms, i1, i2, ju1, j2, k1, k2)
!
! !USES:
      use GmiSpeciesRegistry_mod   , only : getSpeciesPosition
!
     implicit none
!
! !INPUT PARAMETERS:
     integer, intent(in) :: nhms
     integer, intent(in) :: i1, i2, ju1, j2, k1, k2
!
! !INPUT/OUTPUT PARAMETERS:
     type (t_Emission)  , intent(inOut) :: self
!
! !DESCRIPTION:
! Hourly performs the scaling of some emitted species.
! Two scaling factors are used: scFacNOff (for fossil fuel) and scFacNObb
! (for biomass burning). They are utilized to scale the following sets of species:
! \begin{verbatim}
!    Set 1: NO_ff
!    Set 2: NO_bb, CO_bb, MEK_bb, PRPE_bb, C2H6_bb, C3H8_bb, ALK4_bb
!           ALD2_bb, CH2O_bb
! \end{verbatim}
!
! !LOCAL VARIABLES:
      integer                    :: ik
      integer                    :: thisHour

      real*8 , allocatable       :: bbFac (:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
      ! Save out the array read in from the emission file
      if (self%newRecord /= self%curEmissionFileRecord) then
          self%newRecord  = self%curEmissionFileRecord
!
         if (self%doScaleNOffEmiss) &
     &      self%ffNO(:,:,:)     = self%emissionArray(self%iNO_ff)%pArray3D(:,:,:)
!
         if (self%doScaleNObbEmiss) then
            self%bbNO  (:,:,:) = self%emissionArray(self%iNO_bb  )%pArray3D(:,:,:)
            self%bbCO  (:,:,:) = self%emissionArray(self%iCO_bb  )%pArray3D(:,:,:)
            self%bbMEK (:,:,:) = self%emissionArray(self%iMEK_bb )%pArray3D(:,:,:)
            self%bbPRPE(:,:,:) = self%emissionArray(self%iPRPE_bb)%pArray3D(:,:,:)
            self%bbC2H6(:,:,:) = self%emissionArray(self%iC2H6_bb)%pArray3D(:,:,:)
            self%bbC3H8(:,:,:) = self%emissionArray(self%iC3H8_bb)%pArray3D(:,:,:)
            self%bbALK4(:,:,:) = self%emissionArray(self%iALK4_bb)%pArray3D(:,:,:)
            self%bbALD2(:,:,:) = self%emissionArray(self%iALD2_bb)%pArray3D(:,:,:)
            self%bbCH2O(:,:,:) = self%emissionArray(self%iCH2O_bb)%pArray3D(:,:,:)
         end if
      endif
!
      thisHour = nhms/10000
      if (thisHour == 0) thisHour = 24
!
      if (thisHour /= self%currHourEmiss) then
         self%currHourEmiss = thisHour
      end if
!
      do ik = k1, k2
         if (self%doScaleNOffEmiss) then
            self%emissionArray(self%iNO_ff)%pArray3D(:,:,ik) = &
                         self%scFacNOff(:,:,self%currHourEmiss)*self%ffNO(:,:,ik)
         end if
!
         if (self%doScaleNObbEmiss) then
            allocate(bbFac(i1:i2,ju1:j2))
            bbFac(:,:) = self%scFacNObb(:,:,self%currHourEmiss)
!
            self%emissionArray(self%iNO_bb  )%pArray3D(:,:,ik) = bbFac*self%bbNO  (:,:,ik)
            self%emissionArray(self%iCO_bb  )%pArray3D(:,:,ik) = bbFac*self%bbCO  (:,:,ik)
            self%emissionArray(self%iMEK_bb )%pArray3D(:,:,ik) = bbFac*self%bbMEK (:,:,ik)
            self%emissionArray(self%iPRPE_bb)%pArray3D(:,:,ik) = bbFac*self%bbPRPE(:,:,ik)
            self%emissionArray(self%iC2H6_bb)%pArray3D(:,:,ik) = bbFac*self%bbC2H6(:,:,ik)
            self%emissionArray(self%iC3H8_bb)%pArray3D(:,:,ik) = bbFac*self%bbC3H8(:,:,ik)
            self%emissionArray(self%iALK4_bb)%pArray3D(:,:,ik) = bbFac*self%bbALK4(:,:,ik)
            self%emissionArray(self%iALD2_bb)%pArray3D(:,:,ik) = bbFac*self%bbALD2(:,:,ik)
            self%emissionArray(self%iCH2O_bb)%pArray3D(:,:,ik) = bbFac*self%bbCH2O(:,:,ik)
!
            deallocate(bbFac)
         end if
      end do
!
      end subroutine scaleEmissions
!EOC
!-------------------------------------------------------------------------
!BOP
  subroutine initReadEmission (self, gmiClock, gmiGrid, gmiDomain, &
                        Diagnostics, metFields, do_drydep)
!
  use GmiReadGocartSourceFiles_mod, only : GmiReadGocartSourceFiles
  use ReadEmissionFiles_mod, only : ReadEmiss
  use ReadVegLaiData_mod   , only : readVegetationData
  use ReadOtherEmissionData_mod, only : readFertilizerData
  use ReadInputMEGAN_mod       , only : readMEGANannualEmissFactor
  use ReadEmissionFiles_mod    , only : readGalacticCosmisRayPar
!
  implicit none
!
! !INPUT PARAMETERS:
  type (t_gmiGrid  ), intent(in   ) :: gmiGrid
  type (t_GmiClock ), intent(in   ) :: gmiClock
  type (t_gmiDomain ), intent(in   ) :: gmiDomain
  type (t_Diagnostics), intent(in) :: Diagnostics
  logical, intent(in) :: do_drydep
!
! !INPUT/OUTPUT PARAMETERS:
  type (t_metFields ), intent(inOut) :: metFields
  type (t_Emission)  , intent(inOut) :: self
!
! !DESCRIPTION:
! This routines reads in (daily or monthly) emission related files.
!
! !LOCAL VARIABLES:
  character (len=75) :: err_msg
  integer       :: il, ij, num_emiss, ic
  real*8        :: units_fac
  integer       :: ydummy, thisDay, thisMonth, thisDate, ddummy
  integer       :: i1, i2, ju1, j2, k1, k2, procID
  integer       :: i1_gl, i2_gl, ju1_gl, j2_gl, ilong, ilat, ivert
  integer       :: nymd
  real*8        :: tdt
  logical       :: pr_diag, do_aerocom
      real*8 , allocatable :: mcor(:,:)
      real*8 , allocatable :: londeg(:), latdeg(:)
      real*8 , allocatable :: am(:), bm(:)
      integer, allocatable :: cmi_flags(:,:)
      real*8  :: pt
!
!EOP
!----------------------------------------------------------------------
!BOC
      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)
!
      if (pr_diag)  Write (6,*) 'initReadEmission called by ', procID
!
      call Get_i1    (gmiGrid, i1)
      call Get_i2    (gmiGrid, i2)
      call Get_ju1   (gmiGrid, ju1)
      call Get_j2    (gmiGrid, j2)
      call Get_k1    (gmiGrid, k1)
      call Get_k2    (gmiGrid, k2)
      call Get_i1_gl (gmiGrid, i1_gl)
      call Get_i2_gl (gmiGrid, i2_gl)
      call Get_j2_gl (gmiGrid, j2_gl)
      call Get_ju1_gl(gmiGrid, ju1_gl)
!
      ilong = i2 - i1  + 1
      ilat  = j2 - ju1 + 1
!
      allocate(latdeg(ju1_gl:j2_gl))
      allocate(londeg( i1_gl:i2_gl))
      allocate(mcor  (i1:i2,ju1:j2))
      call Get_mcor  (gmiDomain, mcor  )
      call Get_latdeg(gmiDomain, latdeg)
      call Get_londeg(gmiDomain, londeg)
!
      call Get_gmiTimeStep (gmiClock, tdt)
!
      if (do_drydep .or. btest(self%emiss_opt,1)) then
         call readVegetationData  &
              (self%veg_infile_name, self%ireg, self%iland, self%iuse, &
               pr_diag, procID, i1, i2, ju1, j2)
!
         if (btest(self%emiss_opt,1)) then
             call readFertilizerData (self%fertscal_infile_name, &
     &                self%index_soil, self%soil_fert, pr_diag, procID)
         end if
      end if
!
      if (self%emiss_dust_opt == 1 ) then  ! GMI dust emissions
         call Get_do_aerocom(Diagnostics, do_aerocom)
!
         call InitEmissDust (self, mcor, do_aerocom, pr_diag, procID, &
                         i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat)
      end if
!
      if (self%emiss_aero_opt /= 0 ) then
         call InitEmissAero (self, mcor, pr_diag, procID, &
                        i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat)
      end if
!
      if ((self%emiss_dust_opt == 2) .or. (self%emiss_aero_opt == 2))  then
         call GmiReadGocartSourceFiles           &
                (self%GOCARTerod_infile_name, self%GOCARTerod_mod_infile_name, &
                 self%GOCARTocean_infile_name, mcor, latdeg, londeg, &
                 i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl)
      end if
!
      if (self%do_gcr) then
         allocate(am(k1:k2))
         allocate(bm(k1:k2))
         call Get_am(metFields, am)
         call Get_bm(metFields, bm)
         call Get_pt(metFields, pt)
         call Get_curGmiDate  (gmiClock, nymd)
!
         call Allocate_gcr_slope (self, ju1, j2, k1, k2)
         call Allocate_gcr_aintcp (self, ju1, j2, k1, k2)
         call Allocate_GCR_NOx (self, i1, i2, ju1, j2, k1, k2)
!
         call readGalacticCosmisRayPar(self%gcr_infile_name, self%gcr_sunspot, self%gcr_slope, &
           self%gcr_aintcp, nymd, am, bm, pt, latdeg, ju1, j2, k1, k2, ju1_gl, pr_diag, procID)
!
         deallocate(am)
         deallocate(bm)
      end if
!
      ! MEGAN Emissions
!
      if (self%doMEGANemission) then
         call  readMEGANannualEmissFactor &
     &     (self%aefMboMEGAN_InfileName, self%aefIsopMEGAN_InfileName, &
     &     self%aefMonotMEGAN_InfileName, self%aefOvocMEGAN_InfileName, &
     &     self%aefMbo, self%aefIsop, self%aefOvoc, self%aefMonot, &
     &     tdt, mcor, i1, i2, ju1, j2, i1_gl, ju1_gl)
      else
        call Biogenic_Base (self%ireg, self%iland, tdt, &
                    self%convert_isop, self%convert_monot, mcor, &
                    self%base_isop, self%base_monot, &
                    pr_diag, procID, i1, i2, ju1, j2)
      end if
!
      ! Lightning parameterization
!
      if (self%lightning_opt == 1) then
         allocate(cmi_flags(i1:i2,ju1:j2))
         call Get_cmi_flags(metFields, cmi_flags)
!
         call computeGlobCMIflags1(gmiDomain, pr_diag, cmi_flags, &
     &               self%cmi_flags1, i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, &
     &               j2_gl)
!
         call Set_cmi_flags(metFields, cmi_flags)
         deallocate(cmi_flags)
      end if
!
      return
!
      end subroutine initReadEmission
!EOC
!-------------------------------------------------------------------------
!
  subroutine InitEmiss (self, mcor, tr_source_ocean, tr_source_land, &
                        pr_diag, procID, &
                        i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl)
!
  use GmiLandWaterIce_mod, only : setLWIflags
!
  implicit none
!
  type (t_Emission)  , intent(inOut) :: self
  integer, intent(in) :: i1, i2, ju1, j2, k1, k2
  integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
  integer, intent(in) :: procID
  real*8 , intent(in) :: mcor(i1:i2, ju1:j2)
  real*8 , intent(in) :: tr_source_ocean, tr_source_land
  logical, intent(in) :: pr_diag
!
! ----------------------
! Variable declarations.
! ----------------------
!
  character (len=75) :: err_msg
  integer :: il, ij, num_emiss, ic
  real*8  :: tr_source
  real*8  :: units_fac
  integer :: lwi_flags(i1:i2, ju1:j2)
!
! ----------------
! Begin execution.
! ----------------
!
  if (pr_diag) then
     Write (6,*) 'InitEmiss called by ', procID
  end if
!
  num_emiss = self%num_emiss
!
  if (self%emiss_in_opt == 1) then
!
     do ic = 1, num_emiss
        self%emissionArray(ic)%pArray3D(:,:,:) = self%emiss_init_val
     end do
!
! ===========================
  else if (self%emiss_in_opt == 3) then
! ======================
!
    do ic = 1, num_emiss
       self%emissionArray(ic)%pArray3D(:,:,:) = 0.0d0
    end do
!
    call setLWIflags (lwi_flags, i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl)
!
    do ij = ju1, j2
       do il = i1, i2
          if (lwi_flags(il,ij) == 1) then  ! ocean
             tr_source = tr_source_ocean
          else if (lwi_flags(il,ij) == 3) then  ! ice
             tr_source = 0.0d0
          else if (lwi_flags(il,ij) == 2) then  ! land
             tr_source = tr_source_land
          end if
          self%emissionArray(1)%pArray3D(il,ij,k1) = tr_source
       end do
    end do
 end if
!
! ----------------------------
! Convert emiss, if necessary.
! ----------------------------
!
  if (self%emiss_conv_flag /= 0) then
!
     if (self%emiss_conv_flag == 1) then
!
        do ic = 1, num_emiss
           self%emissionArray(ic)%pArray3D(:,:,:) = &
              self%emissionArray(ic)%pArray3D(:,:,:) * self%emiss_conv_fac
        end do
!
     else if (self%emiss_conv_flag == 2) then
!
        units_fac = (1.0d0 / SECPHR) * (KMPM * KMPM)
!
        do ic = 1, num_emiss
           do ij = ju1, j2
              do il = i1, i2
                 self%emissionArray(ic)%pArray3D(il,ij,:) =  &
     &                self%emissionArray(ic)%pArray3D(il,ij,:) * mcor(il,ij) * units_fac
              end do
          end do
         end do
!
        end if
!
      end if
!
  do ic = 1, num_emiss
     call CheckRange3d  &
           ('emissionArray', procID, i1, i2, ju1, j2, k1, k2, &
            self%emissionArray(ic)%pArray3D(:,:,:), -1.0d20, 1.0d20)
  end do
!
  return
!
  end subroutine InitEmiss
!
!-------------------------------------------------------------------------
! This routine sets the aerosol emissions.
!-------------------------------------------------------------------------
!
      subroutine InitEmissAero (self, mcor, pr_diag, procID, &
                                 i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat)
!
      use ReadEmissionFiles_mod, only : ReadEmissAero
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      type (t_Emission), intent(inOut) :: self
      integer          , intent(in   ) :: i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat
      integer          , intent(in   ) :: procID
      real*8           , intent(in   ) :: mcor(i1:i2, ju1:j2)
      logical          , intent(in) :: pr_diag
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'InitEmissAero called by ', procID
      end if
!
      if (self%emiss_aero_opt /= 0) then
!
        call ReadEmissAero(self%emissAero_t, self%emiss_aero_infile_name, &
                           self%naero, self%emiss_timpyr,&
                           pr_diag, procID, i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat)
!
!       -----------------------------------------------------------------
!       Other aero (carbon & sslt) emissions in kg/m^2/s, change to kg/s.
!       -----------------------------------------------------------------
!
        do ij = ju1, j2
           do il = i1, i2
              self%emissAero_t(il,ij,:,:) = self%emissAero_t(il,ij,:,:) * mcor(il,ij)
           end do
        end do
!
      end if
!
      return
!
      end subroutine InitEmissAero
!
!-------------------------------------------------------------------------
! This routine sets the dust emissions.
!-------------------------------------------------------------------------
!
      subroutine InitEmissDust  (self, mcor, do_aerocom, pr_diag, procID, &
                                i1, i2, ju1, j2, i1_gl, ju1_gl, ilong, ilat)
!
  use ReadEmissionFiles_mod, only : ReadEmissDust
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      type (t_Emission)  , intent(inOut) :: self
      integer, intent(in) :: i1, i2, ju1, j2
      integer, intent(in) :: i1_gl, ju1_gl, ilong, ilat
      integer, intent(in) :: procID
      logical, intent(in) :: do_aerocom, pr_diag
      real*8 , intent(in) :: mcor(i1:i2, ju1:j2)
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'InitEmissDust called by ', procID
      end if
!
      if (self%emiss_dust_opt == 1) then  ! GMI dust emissions
!
        call ReadEmissDust(self%emissDust_t, self%emiss_dust_infile_name, &
                      self%ndust, self%nt_dust, self%nst_dust, &
                      do_aerocom, pr_diag, procID, i1, i2, ju1, j2, i1_gl, ju1_gl, &
                      ilong, ilat)
!
!       -------------------------------------------
!       Dust emissions in kg/m^2/s, change to kg/s.
!       -------------------------------------------
!
        do ij = ju1, j2
          do il = i1, i2
            self%emissDust_t(il,ij,:,:) = self%emissDust_t(il,ij,:,:) * mcor(il,ij)
          end do
        end do
!
      end if
!
      return
!
      end subroutine InitEmissDust
!
!-------------------------------------------------------------------------
!
  subroutine InitializeEmission (self, SpeciesConcentration, gmiGrid, gmiDomain, &
                                 Diagnostics, &
                                 do_drydep, ihno3_num, io3_num)
!
  use ReadNOxScalingFactor_mod , only : readNOffScalingFactor
  use ReadNOxScalingFactor_mod , only : readNObbScalingFactor
  use ReadOtherEmissionData_mod, only : readLightData
  use ReadOtherEmissionData_mod, only : readIsopreneConvertData, readSoilData
  use ReadOtherEmissionData_mod, only : readMonoterpeneConvertData
  use GocartDerivedVariables_mod, only : AllocateGocartDerivedVars
  use GmiSeaSaltMethod_mod      , only : InitializationSeaSalt
  use GmiDustMethod_mod         , only : InitializationDust
  use CalcAerosolEmissDiagn_mod , only : setAerosolSurfEmissMap
!
  implicit none
!
  type (t_SpeciesConcentration), intent(in   ) :: SpeciesConcentration
  type (t_gmiGrid )  , intent(in   ) :: gmiGrid
  type (t_gmiDomain )  , intent(in   ) :: gmiDomain
  type (t_Diagnostics), intent(in) :: Diagnostics
  type (t_Emission)  , intent(inOut) :: self
  integer            , intent(in   ) :: io3_num, ihno3_num
  logical            , intent(in   ) :: do_drydep
!
  integer :: ix, ic, id
  integer            i1, i2, ju1, j2, k1, k2
  integer            i1_gl, i2_gl, ju1_gl, j2_gl
  integer   ilong, ilat, ivert, numSpecies, num_emiss, procID
  integer :: tracer_opt
  real*8  :: tr_source_land, tr_source_ocean
  logical :: pr_diag, pr_const, pr_surf_emiss, pr_emiss_3d
  real*8, allocatable :: mcor(:, :)
  logical :: do_aerocom
!EOP
!------------------------------------------------------------------------------
!BOC
  call Get_pr_diag(Diagnostics, pr_diag)
  call Get_pr_const(Diagnostics, pr_const)
  call Get_pr_emiss_3d(Diagnostics, pr_emiss_3d)
  call Get_pr_surf_emiss(Diagnostics, pr_surf_emiss)
  call Get_do_aerocom(Diagnostics, do_aerocom)
!
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
  call Get_ilong (gmiGrid, ilong )
  call Get_ilat  (gmiGrid, ilat  )
  call Get_ivert (gmiGrid, ivert )
  call Get_numSpecies (gmiGrid, numSpecies )
  call Get_num_emiss (self, num_emiss )
!
      allocate(mcor(i1:i2,ju1:j2))
      call Get_mcor(gmiDomain, mcor)
!
  call Get_procID(gmiDomain, procID)
!
  call Get_tracer_opt     (SpeciesConcentration, tracer_opt   )
  call Get_tr_source_land (SpeciesConcentration, tr_source_land )
  call Get_tr_source_ocean(SpeciesConcentration, tr_source_ocean)
!
  self%newRecord        = -1
  self%currHourEmiss    = -1
  self%currDateEmiss    = -1
  self%dailyRecEmiss    = -1
  self%soil_month       = -1
  self%currMonthMEGAN   = -1
  self%lai_day_save     = -1
  self%soil_day_save    = -1
  self%firstReadEmiss   = .TRUE.
  self%firstReadLaiData = .TRUE.
  self%firstScaleEmissions = .TRUE.

  self%curEmissionFileRecord = -1
!

!  if ((tracer_opt == 0) .and. (do_drydep .or. btest(self%emiss_opt,1))) then
  if ((tracer_opt /= 1) .and. &
         (do_drydep .or. btest(self%emiss_opt,1))) then
     if (self%doScaleNOffEmiss) then
        call Allocate_scFacNOff   (self, i1, i2, ju1, j2)
!
        call readNOffScalingFactor(self%scFacNOff, self%scFactorNOff_infile_name, &
                         i1, i2, ju1, j2, i1_gl, ju1_gl)
     end if
!
     if (self%doScaleNObbEmiss) then
        call Allocate_scFacNObb   (self, i1, i2, ju1, j2)
!
        call readNObbScalingFactor(self%scFacNObb, self%scFactorNObb_infile_name, &
                         i1, i2, ju1, j2, i1_gl, ju1_gl)
     end if
!
     call Allocate_ireg       (self, i1, i2, ju1, j2)
     call Allocate_iuse       (self, i1, i2, ju1, j2)
     call Allocate_iland      (self, i1, i2, ju1, j2)
     call Allocate_xlai       (self, i1, i2, ju1, j2)
     call Allocate_xlai2      (self, i1, i2, ju1, j2)
!
     if (btest(self%emiss_opt,1)) then
        NLANDHAR = NLANDHAR_expected(i2_gl)        ! set the value of NLANDHAR
        call Allocate_index_soil (self, NLANDHAR)
        call Allocate_soil_fert  (self, NLANDHAR)
        call Allocate_soil_precip(self, NLANDHAR)
        call Allocate_soil_pulse (self, NLANDHAR)
!
        call readSoilData (self%soil_infile_name, &
     &           self%ncon_soil, pr_diag, procID)
!
        if (self%doMEGANemission) then
           call Allocate_isoLai  (self, i1, i2, ju1, j2)
           call Allocate_isoLaiCurr  (self, i1, i2, ju1, j2)
           call Allocate_isoLaiPrev  (self, i1, i2, ju1, j2)
           call Allocate_isoLaiNext  (self, i1, i2, ju1, j2)
           call Allocate_aefMbo  (self, i1, i2, ju1, j2)
           call Allocate_aefIsop (self, i1, i2, ju1, j2)
           call Allocate_aefOvoc (self, i1, i2, ju1, j2)
           call Allocate_aefMonot(self, i1, i2, ju1, j2)
        else
           call Allocate_base_isop  (self, i1, i2, ju1, j2)
           call Allocate_base_monot (self, i1, i2, ju1, j2)
!
           call readLightData (self%light_infile_name, &
     &              self%coeff_isop, pr_diag, procID)
!
           call readIsopreneConvertData (self%isopconv_infile_name, &
     &              self%convert_isop, pr_diag, procID)
!
           call readMonoterpeneConvertData (self%monotconv_infile_name, &
     &              self%convert_monot, pr_diag, procID)
        end if
!
     end if
  end if
!
  if (self%emiss_in_opt /= 0) then
     call Allocate_emissionArray (self, i1, i2, ju1, j2, k1, k2)
  endif
!
!
  ! Diagnostic variables: surface and 3D emissions
!
  if (pr_const) then
     if (pr_surf_emiss) then
        call Allocate_surf_emiss_out (self, i1, i2, ju1, j2, numSpecies)
        call Allocate_surf_emiss_out2(self, i1, i2, ju1, j2)
!
        if ((self%emiss_aero_opt > 0) .or. (self%emiss_dust_opt > 0)) then
           call Allocate_aerosolSurfEmissMap(self)
           call setAerosolSurfEmissMap(self%aerosolSurfEmissMap, &
                                       self%emiss_map_aero, &
                                       self%emiss_map_dust, &
                                       self%naero, &
                                       self%ndust)
           call Allocate_aerosolSurfEmiss(self, i1, i2, ju1, j2)
        endif
     endif
!
     if (pr_emiss_3d) then
        call Allocate_emiss_3d_out (self, i1, i2, ju1, j2, k1, k2, numSpecies)
        if ((self%emiss_aero_opt > 0) .or. (self%emiss_dust_opt > 0)) then
           call Allocate_aerosolEmiss3D(self, i1, i2, ju1, j2, k1, k2)
        end if
     end if
  end if
!
!.sds  if (self%emiss_opt == 2) then
  if (btest(self%emiss_opt,1)) then
     call Allocate_emiss_monot (self, i1, i2, ju1, j2)
     call Allocate_emiss_isop  (self, i1, i2, ju1, j2)
     call Allocate_emiss_nox   (self, i1, i2, ju1, j2)
  endif
!
  if (self%do_ShipEmission) then
     call Allocate_emiss_o3(self, i1, i2, ju1, j2)
     call Allocate_emiss_hno3 (self, i1, i2, ju1, j2)
!
!     id              = 0
     self%o3_index   = 999
     self%hno3_index = 999
!
     if (io3_num   == 0) stop "the index of   O3 should be non zero"
     if (ihno3_num == 0) stop "the index of HNO3 should be non zero"
!
     do ix = 1, num_emiss
        ic = self%emiss_map(ix)
!        if (ic > 0) id = id + 1
        if (ic == io3_num  ) self%o3_index   = ix
        if (ic == ihno3_num) self%hno3_index = ix
     end do
!
     if (self%o3_index   == 999) stop "O3   is not defined in emiss_map"
     if (self%hno3_index == 999) stop "HNO3 is not defined in emiss_map"
  endif
!
  if (self%emiss_aero_opt /= 0) then
     call Allocate_emissAero   (self, i1, i2, ju1, j2)
     call Allocate_emissAero_t (self, i1, i2, ju1, j2)
  endif
!
  if (self%emiss_dust_opt /= 0) then
     call Allocate_emissDust   (self, i1, i2, ju1, j2)
     call Allocate_emissDust_t (self, i1, i2, ju1, j2)
  endif
!
  if (self%lightning_opt == 1) then
     call Allocate_lightning_no (self, i1, i2, ju1, j2, k1, k2)
     call Allocate_flashrate    (self, i1, i2, ju1, j2)
     call Allocate_cmi_flags1   (self, i1, i2, ju1, j2)
     call Allocate_cldmas0   (self, i1, i2, ju1, j2)
     call Allocate_localCoeff (self, i1, i2, ju1, j2)
     call Allocate_midLatAdj (self, i1, i2, ju1, j2)
     self%globalCoeff = 0.0d0
     self%lightThreshold = 0.0d0
     self%ik0 = -1
  endif
!
  if ((self%emiss_dust_opt == 2) .or. (self%emiss_aero_opt == 2))  then
     call AllocateGocartDerivedVars(i1, i2, ju1, j2, k1, k2)
     if (self%emiss_aero_opt == 2) call InitializationSeaSalt()
     if (self%emiss_dust_opt == 2) call InitializationDust   ()
  end if
!
  ! Set the intial emission array
  call InitEmiss (self, mcor, tr_source_ocean, tr_source_land, &
                        pr_diag, procID, &
                        i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl)
!
  return
!
  end subroutine InitializeEmission
!
!-------------------------------------------------------------------------
!
      subroutine RunEmission &
     &  (self, SpeciesConcentration, gmiClock, gmiGrid, gmiDomain, &
     &   Diagnostics, metFields, &
     &   cosSolarZenithAngle, mw, chem_opt, do_drydep, convec_opt, &
     &   IBOC, IBBC, INOC, IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, ISSLT4, &
     &   IFSO2, INSO2, INDMS, IDUST1, IDUST2, IDUST3, IDUST4, IDUST5, IAN, IMGAS, INO, &
     &   iisoprene_num, ino_num, ico_num, ipropene_num, ihno3_num, io3_num)
!
      use GocartDerivedVariables_mod, only : SetGocartDerivedVars
      use GmiEmissionLightning_mod  , only : emiss_lightning
      use GmiEmissionGCR_mod        , only : Add_Emiss_GCR
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      type (t_Emission)  , intent(inOut) :: self
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
      type(t_GmiClock )  , intent(in   ) :: gmiClock
      type (t_gmiGrid   ), intent(in   ) :: gmiGrid
      type (t_gmiDomain ), intent(in   ) :: gmiDomain
      type (t_metFields ), intent(in   ) :: metFields
      type (t_Diagnostics), intent(in  ) :: Diagnostics
!
      integer, intent(in) :: IBOC, IBBC, INOC, IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3
      integer, intent(in) :: ISSLT4, IFSO2, INSO2, INDMS, IDUST1, IDUST2, IDUST3, IDUST4, IDUST5, IAN, IMGAS, INO
      integer, intent(in) :: iisoprene_num, ino_num, ico_num, ipropene_num
      integer, intent(in) :: ihno3_num, io3_num
      real*8 , intent(in   ) :: mw(:)
      integer, intent(in   ) :: chem_opt, convec_opt
      logical, intent(in   ) :: do_drydep
      real*8 , intent(in   ) :: cosSolarZenithAngle(:, :)
!
!
! Local variables
!
      real*8, allocatable :: mcor(:,:), latdeg(:)
!
      integer :: met_opt
      character (len=MAX_LENGTH_MET_NAME)  :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME)  :: metdata_name_model
!
      real*8 , allocatable :: productionNOx(:,:,:)
      integer, allocatable :: cmi_flags1(:,:)
      real*8 , allocatable :: cldmas0 (:,:)     ! Normalized CLDMAS at desired level
      type(t_GmiArrayBundle), pointer :: concentration(:)
!
      character (len=MAX_LENGTH_ERROR_MSG) err_msg
      integer       :: nhms, nymd, num_time_steps, ndt
      real*8        :: tdt
      integer       :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer       :: i1_gl, i2_gl, ju1_gl, j2_gl, ilong, ilat, ivert
      integer       :: numSpecies, procID
      logical       :: pr_diag, pr_const, pr_emiss_3d, pr_surf_emiss
      logical       :: do_aerocom
!
!
      real*8 , allocatable :: gridBoxHeight (:,:,:)
      real*8 , allocatable :: mass (:,:,:)
      real*8 , allocatable :: press3c (:,:,:)
      real*8 , allocatable :: press3e (:,:,:)
      integer, allocatable :: cmi_flags(:, :)
      integer, allocatable :: lwi_flags(:, :)
      real*8 , allocatable :: pctm2(:, :)
      real*8 , allocatable :: radswg(:, :), surf_air_temp(:, :), surf_rough(:, :)
      real*8 , allocatable :: con_precip(:, :), tot_precip(:, :), ustar(:, :)
      real*8 , allocatable :: u10m(:,:), v10m(:,:), gwet(:,:)
      real*8 , allocatable :: kel(:, :, :), pbl(:, :)
      real*8 , allocatable :: pardif(:,:), pardir(:,:)
      real*8 , allocatable :: fracCloudCover (:,:) , dtrn(:,:,:)
      real*8 , allocatable :: zmmu(:,:,:), cmf(:,:,:)
!
      ! For MEGAN emissions
      real*8 , allocatable :: surfTemp15DayAvg(:,:)
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID (gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)
!
      if (pr_diag) Write(6,*) 'runEmission called by ', procID
!
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
      call Get_numSpecies (gmiGrid, numSpecies)
!
      call Get_pr_const(Diagnostics, pr_const)
      call Get_pr_emiss_3d(Diagnostics, pr_emiss_3d)
      call Get_pr_surf_emiss(Diagnostics, pr_surf_emiss)
      call Get_do_aerocom(Diagnostics, do_aerocom)
!
      allocate(mcor  (i1:i2,ju1:j2))
      allocate(latdeg(ju1_gl:j2_gl))
!
      call Get_mcor  (gmiDomain, mcor  )
      call Get_latdeg(gmiDomain, latdeg)
!
      ! Get metFields variables
!
      call Get_met_opt(metFields, met_opt)
      call Get_metdata_name_org  (metFields, metdata_name_org  )
      call Get_metdata_name_model(metFields, metdata_name_model)
!
      allocate(press3c(ilo:ihi, julo:jhi, k1:k2))
      allocate(press3e(ilo:ihi, julo:jhi, k1-1:k2))
      allocate(mass(i1:i2,ju1:j2,k1:k2))
      allocate(gridBoxHeight(i1:i2,ju1:j2,k1:k2))
!
      call Get_mass(metFields, mass)
      call Get_press3e(metFields, press3e)
      call Get_press3c(metFields, press3c)
      call Get_gridBoxHeight(metFields, gridBoxHeight)
!
      allocate(pctm2(ilo:ihi,julo:jhi))
      allocate(kel(ilo:ihi,julo:jhi,k1:k2))
      call Get_pctm2(metFields, pctm2)
      call Get_kel(metFields, kel)
!
      if (self%lightning_opt == 1) then
         allocate(cmi_flags(i1:i2,ju1:j2))
         call Get_cmi_flags(metFields, cmi_flags)
      end if
!
      if (self%emiss_in_opt /= 0) then
         allocate(lwi_flags(i1:i2,ju1:j2))
         call Get_lwi_flags(metFields, lwi_flags)
      end if
!
      if (met_opt == 3) then
         allocate(con_precip(i1:i2,ju1:j2))
         call Get_con_precip(metFields, con_precip)
!
         allocate(tot_precip(i1:i2,ju1:j2))
         call Get_tot_precip(metFields, tot_precip)
!
         allocate(radswg(i1:i2,ju1:j2))
         call Get_radswg(metFields, radswg)
!
         allocate(surf_air_temp(i1:i2,ju1:j2))
         call Get_surf_air_temp(metFields, surf_air_temp)
!
         allocate(surf_rough(i1:i2,ju1:j2))
         call Get_surf_rough(metFields, surf_rough)
!
         allocate(ustar(i1:i2,ju1:j2))
         call Get_ustar(metFields, ustar)
!
         if ((metdata_name_model(1:5) == 'GEOS4') .or. &
     &       (metdata_name_model(1:5) == 'GEOS5')) then
            allocate(u10m(i1:i2,ju1:j2))
            call Get_u10m(metFields, u10m)
!
            allocate(v10m(i1:i2,ju1:j2))
            call Get_v10m(metFields, v10m)
!
            allocate(gwet(i1:i2,ju1:j2))
            call Get_gwet(metFields, gwet)
!
            allocate(pbl(i1:i2,ju1:j2))
            call Get_pbl(metFields, pbl)
!
            allocate(pardif(i1:i2,ju1:j2))
            call Get_pardif(metFields, pardif)
!
            allocate(pardir(i1:i2,ju1:j2))
            call Get_pardir(metFields, pardir)
!
            if (metdata_name_model(1:5) == 'GEOS4' .and. convec_opt == 3) then
!
               allocate(zmmu(i1:i2,ju1:j2,k1:k2))
               call Get_zmmu(metFields, zmmu)
            end if
         end if
!
         allocate(fracCloudCover(i1:i2,ju1:j2))
         call Get_fracCloudCover(metFields, fracCloudCover)
!
         allocate(dtrn(i1:i2,ju1:j2,k1:k2))
         call Get_dtrn(metFields, dtrn)
!
         allocate(cmf(i1:i2,ju1:j2,k1:k2))
         call Get_cmf(metFields, cmf)
      end if
!
      ! Obtain model time information
!
      call Get_curGmiDate  (gmiClock, nymd          )
      call Get_curGmiTime  (gmiClock, nhms          )
      call Get_numTimeSteps(gmiClock, num_time_steps)
      call Get_gmiTimeStep (gmiClock, tdt           )
      ndt = Nint (tdt)
!
      call Get_concentration  (SpeciesConcentration, concentration  )
!
!.sds      if (self%emiss_opt == 2) then
      if (btest(self%emiss_opt,1)) then
         self%emiss_monot = 0.0d0
         self%emiss_isop  = 0.0d0
         self%emiss_nox   = 0.0d0
      end if
!
      if (self%do_ShipEmission) then
         self%emiss_hno3(:,:) = self%emissionArray(self%hno3_index)%pArray3D(:,:,1)
      end if
!
! For GOCART emission
!
      if ((self%emiss_dust_opt == 2) .or. (self%emiss_aero_opt == 2))  then
         call SetGocartDerivedVars (u10m, v10m, gwet, press3c, &
     &                    kel, mass, mcor, gridBoxHeight, lwi_flags,  &
     &                    i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
      end if
!
      if (btest(self%emiss_opt,1)) then
         if (self%doMEGANemission) then
            allocate(surfTemp15DayAvg(i1:i2,ju1:j2))
            call Get_surfTemp15DayAvg(metFields, surfTemp15DayAvg)
         end if
      end if
!
! Emission control routines
      call Update_Emiss (lwi_flags, cosSolarZenithAngle, latdeg, mcor,        &
     &            self%emiss_isop, self%emiss_monot, self%emiss_nox,           &
     &            self%do_ShipEmission, self%emiss_hno3, self%emiss_o3,        &
     &            ihno3_num, io3_num, self%doMEGANemission, self%days_btw_m,   &
     &            surfTemp15DayAvg, self%aefIsop, self%aefMbo, self%aefMonot,  &
     &            self%isoLai, self%isoLaiCurr, self%isoLaiPrev, radswg,       &
     &            surf_air_temp, pardif, pardir, surf_rough, con_precip,       &
     &            tot_precip, ustar, mass, fracCloudCover, kel,                &
     &            self%surf_emiss_out, self%surf_emiss_out2, self%emiss_3d_out,&
     &            self%aerosolEmiss3D, self%aerosolSurfEmiss,                  &
     &            self%aerosolSurfEmissMap, concentration, self%emissionArray, &
     &            self%emissDust_t, self%emissDust, self%emissAero_t,          &
     &            self%emissAero, pbl, gridBoxHeight, self%index_soil,         &
     &            self%ncon_soil, self%soil_fert, self%soil_precip,            &
     &            self%soil_pulse, self%ireg, self%iland, self%iuse,           &
     &            self%convert_isop, self%convert_monot, self%coeff_isop,      &
     &            self%base_isop, self%base_monot, self%xlai, IBOC, IBBC, INOC,&
     &            IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, ISSLT4, IFSO2, INSO2,    &
     &            INDMS, IDUST1, IDUST2, IDUST3, IDUST4, IDUST5, IAN, IMGAS, INO, iisoprene_num, ino_num, ico_num,     &
     &            ipropene_num, pr_surf_emiss, pr_emiss_3d, pr_diag, procID,   &
     &            met_opt, self%emiss_opt, chem_opt, self%emiss_aero_opt,      &
     &            self%emiss_dust_opt, do_aerocom, self%do_semiss_inchem,      &
     &            do_drydep, self%emiss_map, self%emiss_map_dust,              &
     &            self%emiss_map_aero, self%ndust, self%nst_dust, self%nt_dust,&
     &            self%naero, nymd, num_time_steps, mw, tdt, ndt,              &
     &            self%emiss_timpyr, self%num_emiss, self%isop_scale, i1, i2,  &
     &            ju1, j2, k1, k2, ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl,  &
     &            j2_gl, ilong, numSpecies, self%soil_day_save)
!
!=================================================
!   Do Galactic Cosmic Ray Emissions of N and NO
!=================================================
!
      if (self%do_gcr) then
         self%GCR_NOx(:,:,:) = 0.d0
         Allocate (productionNOx(i1:i2, ju1:j2, k1:k2))
         productionNOx = 0.d0
!
!        ===================
         call Add_Emiss_GCR  &
!        ===================
     &     (self%gcr_sunspot, self%gcr_slope, self%gcr_aintcp, concentration, IMGAS, tdt,  &
     &      productionNOx, self%GCR_NOx, mass, pr_diag, procID, i1, i2, ju1, j2, k1, k2, numSpecies)
!
!... 55% to atomic N
         if(IAN.gt.0) concentration(IAN)%pArray3D(:,:,:) = &
     &          concentration(IAN)%pArray3D(:,:,:) + ( 0.55d0 * productionNOx(:,:,:) * tdt )
!
!... 45% to NO
         if(INO.gt.0) concentration(INO)%pArray3D(:,:,:) = &
     &          concentration(INO)%pArray3D(:,:,:) + ( 0.45d0 * productionNOx(:,:,:) * tdt )
!
         deallocate(productionNOx)
!
      end if
!
!=================================================
!
      if (self%lightning_opt == 1) then
         Allocate (cldmas0(i1:i2, ju1:j2))
         Allocate (productionNOx(i1:i2, ju1:j2, k1:k2))
!
         self%flashrate     = 0.d0
         self%lightning_no  = 0.d0
         productionNOx = 0.d0
!
         if (metdata_name_org(1:3) == 'DAO') then
            cldmas0 = cmf(:,:,12)
         else if(metdata_name_org(1:4) == 'NCAR') then
            cldmas0 = cmf(:,:,10)
         else if(metdata_name_org(1:4) == 'GISS') then
            cldmas0 = cmf(:,:,8)
         else if(metdata_name_org(1:4) == 'GMAO') then
            if (metdata_name_model(1:5) == 'GEOS4') cldmas0 = zmmu(:,:,self%ik0-1)
            if (metdata_name_model(1:5) == 'GEOS5') cldmas0 = cmf(:,:,self%ik0-1)
            !if (num_time_steps .lt. 2) then
            !   write (6,*) 'Assigning cmf to cldmas0 layer : ', self%ik0-1
            !endif
         else
            write (6,*) ' Lightning parameterization does',  &
     &             '  not exist for this model'
         endif
!
!
         call emiss_lightning(metdata_name_org, metdata_name_model, nymd, &
     &              i1,i2 ,ju1,j2, k1, k2, numSpecies, ilo, ihi, julo, jhi, &
     &              tdt, procID, pr_diag, self%desired_g_N_prod_rate, cldmas0, &
     &              press3c, press3e, cmi_flags,  self%cmi_flags1, i2_gl, j2_gl, &
     &              mass, dtrn, productionNOx, self%flashrate, self%lightning_no, &
     &              self%lightThreshold, self%globalCoeff, self%localCoeff, &
     &		    self%midLatAdj)
!
         concentration(ino_num)%pArray3D(:,:,:) = &
     &          concentration(ino_num)%pArray3D(:,:,:) + productionNOx(:,:,:) * tdt
!
         deallocate(cldmas0)
         deallocate(productionNOx)
!
      endif
!
      ! Passing the updated variables to the derived types before exiting the routine.
      ! The variables are the arguments of Update_Emiss that are INOUT or OUT.
!
      if (btest(self%emiss_opt,1)) then
         if (self%doMEGANemission) then
             deallocate(surfTemp15DayAvg)
         end if
      end if
!
      deallocate(mcor  )
      deallocate(latdeg)
      deallocate(mass  )
      deallocate(press3c)
      deallocate(press3e)
      deallocate(gridBoxHeight)
!
      return
!
      end subroutine RunEmission
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finalizeEmission
!
! !INTERFACE:
!
      subroutine finalizeEmission (self, Diagnostics, tracer_opt, do_drydep)
!
  implicit none
!
! !INPUT PARAMETERS:
      integer            , intent(in) :: tracer_opt
      logical            , intent(in) :: do_drydep
      type(t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
  type (t_Emission)   , intent(inOut) :: self
!
! !DESCRIPTION:
! Finalize the Emission component by deallocating member variables of the
! Emission derived type.
!
! !LOCAL VARIABLES:
      logical :: pr_const, pr_surf_emiss, pr_emiss_3d
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (self%emiss_in_opt /= 0) then
         deallocate(self%emiss)
      endif
!
!      if ((tracer_opt == 0) .and. (do_drydep .or. btest(self%emiss_opt,1))) then
      if ((tracer_opt /= 1) .and. &
     &     (do_drydep .or. btest(self%emiss_opt,1))) then
         if (self%doScaleNOffEmiss) then
            deallocate(self%scFacNOff   )
         end if
!
         if (self%doScaleNObbEmiss) then
            deallocate(self%scFacNObb   )
         end if
!
         deallocate(self%ireg       )
         deallocate(self%iuse       )
         deallocate(self%iland      )
         deallocate(self%xlai       )
         deallocate(self%xlai2      )
!
         if (btest(self%emiss_opt,1)) then
            deallocate(self%index_soil )
            deallocate(self%soil_fert  )
            deallocate(self%soil_precip)
            deallocate(self%soil_pulse )
!
            if (self%doMEGANemission) then
               deallocate(self%isoLai  )
               deallocate(self%isoLaiCurr  )
               deallocate(self%isoLaiPrev  )
               deallocate(self%isoLaiNext  )
               deallocate(self%aefMbo  )
               deallocate(self%aefIsop )
               deallocate(self%aefOvoc )
               deallocate(self%aefMonot)
            else
               deallocate(self%base_isop  )
               deallocate(self%base_monot )
            end if
!
         end if
      end if
!
      if (self%emiss_in_opt /= 0) then
         deallocate(self%emissionArray)
      endif
!
      ! Diagnostic variables: surface and 3D emissions
!
      call Get_pr_const     (Diagnostics, pr_const     )
      call Get_pr_emiss_3d  (Diagnostics, pr_emiss_3d  )
      call Get_pr_surf_emiss(Diagnostics, pr_surf_emiss)
!
      if (pr_const) then
         if (pr_surf_emiss) then
            deallocate(self%surf_emiss_out)
            deallocate(self%surf_emiss_out2)
!
            if ((self%emiss_aero_opt > 0) .or. (self%emiss_dust_opt > 0)) then
               deallocate(self%aerosolSurfEmissMap)
               deallocate(self%aerosolSurfEmiss)
            endif
         endif
!
         if (pr_emiss_3d) then
            deallocate(self%emiss_3d_out)
            if ((self%emiss_aero_opt > 0) .or. (self%emiss_dust_opt > 0)) then
               deallocate(self%aerosolEmiss3D)
            end if
         end if
      end if
!
!.sds      if (self%emiss_opt == 2) then
      if (btest(self%emiss_opt,1)) then
         deallocate(self%emiss_monot)
         deallocate(self%emiss_isop)
         deallocate(self%emiss_nox)
      endif
!
      if (self%do_ShipEmission) then
         deallocate(self%emiss_o3)
         deallocate(self%emiss_hno3)
      endif
!
      if (self%emiss_aero_opt /= 0) then
         deallocate(self%emissAero)
         deallocate(self%emissAero_t)
      endif
!
      if (self%emiss_dust_opt /= 0) then
         deallocate(self%emissDust)
         deallocate(self%emissDust_t)
      endif
!
      if (self%do_gcr) then
         deallocate(self%gcr_slope)
         deallocate(self%gcr_aintcp)
         deallocate(self%GCR_NOx)
      endif
!
      if (self%lightning_opt == 1) then
         deallocate(self%lightning_no)
         deallocate(self%flashrate)
         deallocate(self%cmi_flags1)
         deallocate(self%cldmas0)
         deallocate(self%localCoeff)
	 deallocate(self%midLatAdj)
      endif
!
      return
!
      end subroutine finalizeEmission
!EOC
!-------------------------------------------------------------------------
subroutine Allocate_cmi_flags1 (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inout) :: self
    Allocate(self%cmi_flags1(i1:i2, ju1:j2))
    self%cmi_flags1 = 0
    return
  end subroutine Allocate_cmi_flags1
!-------------------------------------------------------------------------
  subroutine Get_cmi_flags1 (self, cmi_flags1)
    implicit none
    integer         , intent(out)  :: cmi_flags1 (:,:)
    type (t_Emission), intent(in)   :: self
    cmi_flags1(:,:) = self%cmi_flags1(:,:)
    return
  end subroutine Get_cmi_flags1
!-------------------------------------------------------------------------
  subroutine Set_cmi_flags1 (self, cmi_flags1)
    implicit none
    integer          , intent(in)  :: cmi_flags1 (:,:)
    type (t_Emission), intent(inout) :: self
    self%cmi_flags1(:,:) = cmi_flags1(:,:)
    return
  end subroutine Set_cmi_flags1
!-------------------------------------------------------------------------
  subroutine Get_cldmas0 (self, cldmas0)
    implicit none
    real*8          , intent(out)  :: cldmas0 (:,:)
    type (t_Emission), intent(in)   :: self
    cldmas0(:,:) = self%cldmas0 (:,:)
    return
  end subroutine Get_cldmas0
!-------------------------------------------------------------------------
  subroutine Allocate_cldmas0 (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inout) :: self
    Allocate(self%cldmas0(i1:i2, ju1:j2))
    self%cldmas0 = 0.0d0
    return
  end subroutine Allocate_cldmas0
!-------------------------------------------------------------------------
  subroutine Set_cldmas0 (self, cldmas0)
    implicit none
    real*8          , intent(in)  :: cldmas0 (:,:)
    type (t_Emission), intent(inout) :: self
    self%cldmas0(:,:) = cldmas0(:,:)
    return
  end subroutine Set_cldmas0
!-------------------------------------------------------------------------
  subroutine Allocate_emissionArray (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_Emission), intent(inOut) :: self
    integer                          :: i, num_emiss
    num_emiss = self%num_emiss
    Allocate(self%emissionArray(num_emiss))
    do i = 1, num_emiss
       Allocate(self%emissionArray(i)%pArray3D(i1:i2, ju1:j2, k1:k2))
       self%emissionArray(i)%pArray3D(i1:i2, ju1:j2, k1:k2) = 0.0d0
    end do
    return
  end subroutine Allocate_emissionArray
!-------------------------------------------------------------------------
  subroutine Get_emissionArray (self, emissionArray)
    implicit none
    type (t_GmiArrayBundle), pointer :: emissionArray (:)
    type (t_Emission), intent(in)   :: self
    emissionArray => self%emissionArray
    return
  end subroutine Get_emissionArray
!-------------------------------------------------------------------------
  subroutine Set_emissionArray (self, emissionArray)
    implicit none
    type (t_GmiArrayBundle), pointer :: emissionArray (:)
    type (t_Emission), intent(inOut) :: self
    self%emissionArray => emissionArray
    return
  end subroutine Set_emissionArray
!-------------------------------------------------------------------------
  subroutine Allocate_base_isop (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%base_isop(i1:i2, ju1:j2, NTYPE))
    self%base_isop = 0.0d0
    return
  end subroutine Allocate_base_isop
!-------------------------------------------------------------------------
  subroutine Get_base_isop (self, base_isop)
    implicit none
    real*8          , intent(out)  :: base_isop (:,:,:)
    type (t_Emission), intent(in)   :: self
    base_isop(:,:,:) = self%base_isop(:,:,:)
    return
  end subroutine Get_base_isop
!-------------------------------------------------------------------------
  subroutine Set_base_isop (self, base_isop)
    implicit none
    real*8          , intent(in)  :: base_isop (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%base_isop(:,:,:) = base_isop(:,:,:)
    return
  end subroutine Set_base_isop
!-------------------------------------------------------------------------
  subroutine Allocate_base_monot (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%base_monot(i1:i2, ju1:j2, NTYPE))
    self%base_monot = 0.0d0
    return
  end subroutine Allocate_base_monot
!-------------------------------------------------------------------------
  subroutine Get_base_monot (self, base_monot)
    implicit none
    real*8          , intent(out)  :: base_monot (:,:,:)
    type (t_Emission), intent(in)   :: self
    base_monot(:,:,:) = self%base_monot(:,:,:)
    return
  end subroutine Get_base_monot
!-------------------------------------------------------------------------
  subroutine Set_base_monot (self, base_monot)
    implicit none
    real*8          , intent(in)  :: base_monot (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%base_monot(:,:,:) = base_monot(:,:,:)
    return
  end subroutine Set_base_monot
!-------------------------------------------------------------------------
  subroutine Allocate_xlai (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%xlai(i1:i2, ju1:j2, NTYPE))
    self%xlai = 0.0d0
    return
  end subroutine Allocate_xlai
!-------------------------------------------------------------------------
  subroutine Get_xlai (self, xlai)
    implicit none
    real*8          , intent(out)  :: xlai (:,:,:)
    type (t_Emission), intent(in)   :: self
    xlai(:,:,:) = self%xlai(:,:,:)
    return
  end subroutine Get_xlai
!-------------------------------------------------------------------------
  subroutine Set_xlai (self, xlai)
    implicit none
    real*8          , intent(in)  :: xlai (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%xlai(:,:,:) = xlai(:,:,:)
    return
  end subroutine Set_xlai
!-------------------------------------------------------------------------
  subroutine Allocate_xlai2 (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%xlai2(i1:i2, ju1:j2, NTYPE))
    self%xlai2 = 0.0d0
    return
  end subroutine Allocate_xlai2
!-------------------------------------------------------------------------
  subroutine Get_xlai2 (self, xlai2)
    implicit none
    real*8          , intent(out)  :: xlai2 (:,:,:)
    type (t_Emission), intent(in)   :: self
    xlai2(:,:,:) = self%xlai2(:,:,:)
    return
  end subroutine Get_xlai2
!-------------------------------------------------------------------------
  subroutine Set_xlai2 (self, xlai2)
    implicit none
    real*8          , intent(in)  :: xlai2 (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%xlai2(:,:,:) = xlai2(:,:,:)
    return
  end subroutine Set_xlai2
!-------------------------------------------------------------------------
  subroutine Get_coeff_isop (self, coeff_isop)
    implicit none
    real*8          , intent(out)  :: coeff_isop (:)
    type (t_Emission), intent(in)   :: self
    coeff_isop(:) = self%coeff_isop(:)
    return
  end subroutine Get_coeff_isop
!-------------------------------------------------------------------------
  subroutine Set_coeff_isop (self, coeff_isop)
    implicit none
    real*8          , intent(in)  :: coeff_isop (:)
    type (t_Emission), intent(inOut) :: self
    self%coeff_isop(:) = coeff_isop(:)
    return
  end subroutine Set_coeff_isop
!-------------------------------------------------------------------------
  subroutine Get_convert_isop (self, convert_isop)
    implicit none
    real*8          , intent(out)  :: convert_isop (:)
    type (t_Emission), intent(in)   :: self
    convert_isop(:) = self%convert_isop(:)
    return
  end subroutine Get_convert_isop
!-------------------------------------------------------------------------
  subroutine Set_convert_isop (self, convert_isop)
    implicit none
    real*8           , intent(in)  :: convert_isop (:)
    type (t_Emission), intent(inOut) :: self
    self%convert_isop(:) = convert_isop(:)
    return
  end subroutine Set_convert_isop
!-------------------------------------------------------------------------
  subroutine Get_convert_monot (self, convert_monot)
    implicit none
    real*8           , intent(out)  :: convert_monot (:)
    type (t_Emission), intent(in)   :: self
    convert_monot(:) = self%convert_monot(:)
    return
  end subroutine Get_convert_monot
!-------------------------------------------------------------------------
  subroutine Set_convert_monot (self, convert_monot)
    implicit none
    real*8           , intent(in)  :: convert_monot (:)
    type (t_Emission), intent(inOut) :: self
    self%convert_monot(:) = convert_monot(:)
    return
  end subroutine Set_convert_monot
!-------------------------------------------------------------------------
  subroutine Allocate_soil_fert (self, NLANDHAR)
    implicit none
    integer          , intent(in   ) :: NLANDHAR
    type (t_Emission), intent(inOut) :: self
    Allocate(self%soil_fert(NLANDHAR))
    !self%soil_fert = 0.0d0
    return
  end subroutine Allocate_soil_fert
!-------------------------------------------------------------------------
  subroutine Get_soil_fert (self, soil_fert)
    implicit none
    real*8           , intent(out)  :: soil_fert (:)
    type (t_Emission), intent(in)   :: self
    soil_fert(:) = self%soil_fert(:)
    return
  end subroutine Get_soil_fert
!-------------------------------------------------------------------------
  subroutine Set_soil_fert (self, soil_fert)
    implicit none
    real*8           , intent(in)  :: soil_fert (:)
    type (t_Emission), intent(inOut) :: self
    self%soil_fert(:) = soil_fert(:)
    return
  end subroutine Set_soil_fert
!-------------------------------------------------------------------------
  subroutine Allocate_soil_precip (self, NLANDHAR)
    implicit none
    integer          , intent(in   ) :: NLANDHAR
    type (t_Emission), intent(inOut) :: self
    Allocate(self%soil_precip(2,NLANDHAR))
    !self%soil_precip = 0.0d0
    return
  end subroutine Allocate_soil_precip
!-------------------------------------------------------------------------
  subroutine Get_soil_precip (self, soil_precip)
    implicit none
    real*8           , intent(out)  :: soil_precip (:,:)
    type (t_Emission), intent(in)   :: self
    soil_precip(:,:) = self%soil_precip(:,:)
    return
  end subroutine Get_soil_precip
!-------------------------------------------------------------------------
  subroutine Set_soil_precip (self, soil_precip)
    implicit none
    real*8           , intent(in)  :: soil_precip (:,:)
    type (t_Emission), intent(inOut) :: self
    self%soil_precip(:,:) = soil_precip(:,:)
    return
  end subroutine Set_soil_precip
!-------------------------------------------------------------------------
  subroutine Allocate_soil_pulse (self, NLANDHAR)
    implicit none
    integer          , intent(in   ) :: NLANDHAR
    type (t_Emission), intent(inOut) :: self
    Allocate(self%soil_pulse(NPULSE+1,NLANDHAR))
    self%soil_pulse = 0.0d0
    return
  end subroutine Allocate_soil_pulse
!-------------------------------------------------------------------------
  subroutine Get_soil_pulse (self, soil_pulse)
    implicit none
    real*8           , intent(out)  :: soil_pulse (:,:)
    type (t_Emission), intent(in)   :: self
    soil_pulse(:,:) = self%soil_pulse(:,:)
    return
  end subroutine Get_soil_pulse
!-------------------------------------------------------------------------
  subroutine Set_soil_pulse (self, soil_pulse)
    implicit none
    real*8           , intent(in)  :: soil_pulse (:,:)
    type (t_Emission), intent(inOut) :: self
    self%soil_pulse(:,:) = soil_pulse(:,:)
    return
  end subroutine Set_soil_pulse
!-------------------------------------------------------------------------
  subroutine Allocate_iland (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%iland(i1:i2, ju1:j2, NTYPE))
    self%iland = 0
    return
  end subroutine Allocate_iland
!-------------------------------------------------------------------------
  subroutine Get_iland (self, iland)
    implicit none
    integer          , intent(out)  :: iland (:,:,:)
    type (t_Emission), intent(in)   :: self
    iland(:,:,:) = self%iland(:,:,:)
    return
  end subroutine Get_iland
!-------------------------------------------------------------------------
  subroutine Set_iland (self, iland)
    implicit none
    integer          , intent(in)  :: iland (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%iland(:,:,:) = iland(:,:,:)
    return
  end subroutine Set_iland
!-------------------------------------------------------------------------
  subroutine Allocate_iuse (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%iuse(i1:i2, ju1:j2, NTYPE))
    self%iuse = 0
    return
  end subroutine Allocate_iuse
!-------------------------------------------------------------------------
  subroutine Get_iuse (self, iuse)
    implicit none
    integer          , intent(out)  :: iuse (:,:,:)
    type (t_Emission), intent(in)   :: self
    iuse(:,:,:) = self%iuse(:,:,:)
    return
  end subroutine Get_iuse
!-------------------------------------------------------------------------
  subroutine Set_iuse (self, iuse)
    implicit none
    integer          , intent(in)  :: iuse (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%iuse(:,:,:) = iuse(:,:,:)
    return
  end subroutine Set_iuse
!-------------------------------------------------------------------------
  subroutine Allocate_ireg (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    Allocate(self%ireg(i1:i2, ju1:j2))
    self%ireg = 0
    return
  end subroutine Allocate_ireg
!-------------------------------------------------------------------------
  subroutine Get_ireg (self, ireg)
    implicit none
    integer          , intent(out)  :: ireg (:,:)
    type (t_Emission), intent(in)   :: self
    ireg(:,:) = self%ireg(:,:)
    return
  end subroutine Get_ireg
!-------------------------------------------------------------------------
  subroutine Set_ireg (self, ireg)
    implicit none
    integer          , intent(in)  :: ireg (:,:)
    type (t_Emission), intent(inOut) :: self
    self%ireg(:,:) = ireg(:,:)
    return
  end subroutine Set_ireg
!-------------------------------------------------------------------------
  subroutine Allocate_index_soil (self, NLANDHAR)
    implicit none
    integer          , intent(in   ) :: NLANDHAR
    type (t_Emission), intent(inOut) :: self
    Allocate(self%index_soil(2, NLANDHAR))
    return
  end subroutine Allocate_index_soil
!-------------------------------------------------------------------------
  subroutine Get_index_soil (self, index_soil)
    implicit none
    integer          , intent(out)  :: index_soil (:,:)
    type (t_Emission), intent(in)   :: self
    index_soil(:,:) = self%index_soil(:,:)
    return
  end subroutine Get_index_soil
!-------------------------------------------------------------------------
  subroutine Set_index_soil (self, index_soil)
    implicit none
    integer          , intent(in)  :: index_soil (:,:)
    type (t_Emission), intent(inOut) :: self
    self%index_soil(:,:) = index_soil(:,:)
    return
  end subroutine Set_index_soil
!-------------------------------------------------------------------------
  subroutine Get_ncon_soil (self, ncon_soil)
    implicit none
    integer          , intent(out)  :: ncon_soil (:)
    type (t_Emission), intent(in)   :: self
    ncon_soil(:) = self%ncon_soil(:)
    return
  end subroutine Get_ncon_soil
!-------------------------------------------------------------------------
  subroutine Set_ncon_soil (self, ncon_soil)
    implicit none
    integer          , intent(in)  :: ncon_soil (:)
    type (t_Emission), intent(inOut) :: self
    self%ncon_soil(:) = ncon_soil(:)
    return
  end subroutine Set_ncon_soil
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  subroutine Allocate_emissAero (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)       , intent(inOut) :: self
    Allocate(self%emissAero(i1:i2, ju1:j2, self%naero))
    self%emissAero = 0.0d0
    return
  end subroutine Allocate_emissAero
!-------------------------------------------------------------------------
  subroutine Get_emissAero (self, emiss_aero)
    implicit none
    real*8          , intent(out)  :: emiss_aero (:,:,:)
    type (t_Emission), intent(in)   :: self
    emiss_aero(:,:,:) = self%emissAero(:,:,:)
    return
  end subroutine Get_emissAero
!-------------------------------------------------------------------------
  subroutine Set_emissAero (self, emiss_aero)
    implicit none
    real*8          , intent(in)  :: emiss_aero (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%emissAero(:,:,:) = emiss_aero(:,:,:)
    return
  end subroutine Set_emissAero
!-------------------------------------------------------------------------
  subroutine Allocate_emissDust (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)       , intent(inOut) :: self
    Allocate(self%emissDust(i1:i2, ju1:j2, 1:self%ndust))
    self%emissDust = 0.0d0
    return
  end subroutine Allocate_emissDust
!-------------------------------------------------------------------------
  subroutine Get_emissDust (self, emiss_dust)
    implicit none
    real*8           , intent(out)  :: emiss_dust (:,:,:)
    type (t_Emission), intent(in )  :: self
    emiss_dust(:,:,:) = self%emissDust(:,:,:)
    return
  end subroutine Get_emissDust
!-------------------------------------------------------------------------
  subroutine Set_emissDust (self, emiss_dust)
    implicit none
    real*8           , intent(in   ) :: emiss_dust (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%emissDust(:,:,:) = emiss_dust(:,:,:)
    return
  end subroutine Set_emissDust
!-------------------------------------------------------------------------
  subroutine Allocate_lightning_no (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%lightning_no(i1:i2, ju1:j2, k1:k2))
    self%lightning_no = 0.0d0
    return
  end subroutine Allocate_lightning_no
!-------------------------------------------------------------------------
  subroutine Get_lightning_no (self, lightning_no)
    implicit none
    real*8          , intent(out)  :: lightning_no (:,:,:)
    type (t_Emission), intent(in)   :: self
    lightning_no(:,:,:) = self%lightning_no(:,:,:)
    return
  end subroutine Get_lightning_no
!-------------------------------------------------------------------------
  subroutine Set_lightning_no (self, lightning_no)
    implicit none
    real*8          , intent(in)  :: lightning_no (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%lightning_no(:,:,:) = lightning_no(:,:,:)
    return
  end subroutine Set_lightning_no
!-------------------------------------------------------------------------
  subroutine Get_lightCoeff_infile_name (self, lightCoeff_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: lightCoeff_infile_name
    type (t_Emission)    , intent(in)   :: self
    lightCoeff_infile_name = self%lightCoeff_infile_name
    return
  end subroutine Get_lightCoeff_infile_name
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_isop (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_isop(i1:i2, ju1:j2))
    self%emiss_isop = 0.0d0
    return
  end subroutine Allocate_emiss_isop
!-------------------------------------------------------------------------
  subroutine Get_emiss_isop (self, emiss_isop)
    implicit none
    real*8          , intent(out)  :: emiss_isop (:,:)
    type (t_Emission), intent(in)   :: self
    emiss_isop(:,:) = self%emiss_isop(:,:)
    return
  end subroutine Get_emiss_isop
!-------------------------------------------------------------------------
  subroutine Set_emiss_isop (self, emiss_isop)
    implicit none
    real*8          , intent(in)  :: emiss_isop (:,:)
    type (t_Emission), intent(inOut) :: self
    self%emiss_isop(:,:) = emiss_isop(:,:)
    return
  end subroutine Set_emiss_isop
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_monot (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_monot(i1:i2, ju1:j2))
    self%emiss_monot = 0.0d0
    return
  end subroutine Allocate_emiss_monot
!-------------------------------------------------------------------------
  subroutine Get_emiss_monot (self, emiss_monot)
    implicit none
    real*8          , intent(out)  :: emiss_monot (:,:)
    type (t_Emission), intent(in)   :: self
    emiss_monot(:,:) = self%emiss_monot(:,:)
    return
  end subroutine Get_emiss_monot
!-------------------------------------------------------------------------
  subroutine Set_emiss_monot (self, emiss_monot)
    implicit none
    real*8          , intent(in)  :: emiss_monot (:,:)
    type (t_Emission), intent(inOut) :: self
    self%emiss_monot(:,:) = emiss_monot(:,:)
    return
  end subroutine Set_emiss_monot
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_nox (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_nox(i1:i2, ju1:j2))
    self%emiss_nox = 0.0d0
    return
  end subroutine Allocate_emiss_nox
!-------------------------------------------------------------------------
  subroutine Get_emiss_nox (self, emiss_nox)
    implicit none
    real*8          , intent(out)  :: emiss_nox (:,:)
    type (t_Emission), intent(in)   :: self
    emiss_nox(:,:) = self%emiss_nox(:,:)
    return
  end subroutine Get_emiss_nox
!-------------------------------------------------------------------------
  subroutine Set_emiss_nox (self, emiss_nox)
    implicit none
    real*8          , intent(in)  :: emiss_nox (:,:)
    type (t_Emission), intent(inOut) :: self
    self%emiss_nox(:,:) = emiss_nox(:,:)
    return
  end subroutine Set_emiss_nox
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_hno3 (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_hno3(i1:i2, ju1:j2))
    self%emiss_hno3 = 0.0d0
    return
  end subroutine Allocate_emiss_hno3
!-------------------------------------------------------------------------
  subroutine Get_emiss_hno3 (self, emiss_hno3)
    implicit none
    real*8          , intent(out)  :: emiss_hno3 (:,:)
    type (t_Emission), intent(in)   :: self
    emiss_hno3(:,:) = self%emiss_hno3(:,:)
    return
  end subroutine Get_emiss_hno3
!-------------------------------------------------------------------------
  subroutine Set_emiss_hno3 (self, emiss_hno3)
    implicit none
    real*8          , intent(in)  :: emiss_hno3 (:,:)
    type (t_Emission), intent(inOut) :: self
    self%emiss_hno3(:,:) = emiss_hno3(:,:)
    return
  end subroutine Set_emiss_hno3
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_o3 (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%emiss_o3(i1:i2, ju1:j2))
    self%emiss_o3 = 0.0d0
    return
  end subroutine Allocate_emiss_o3
!-------------------------------------------------------------------------
  subroutine Get_emiss_o3 (self, emiss_o3)
    implicit none
    real*8          , intent(out)  :: emiss_o3 (:,:)
    type (t_Emission), intent(in)   :: self
    emiss_o3(:,:) = self%emiss_o3(:,:)
    return
  end subroutine Get_emiss_o3
!-------------------------------------------------------------------------
  subroutine Set_emiss_o3 (self, emiss_o3)
    implicit none
    real*8          , intent(in)  :: emiss_o3 (:,:)
    type (t_Emission), intent(inOut) :: self
    self%emiss_o3(:,:) = emiss_o3(:,:)
    return
  end subroutine Set_emiss_o3
!-------------------------------------------------------------------------
  subroutine Allocate_surf_emiss_out (self, i1, i2, ju1, j2, numSpecies)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2, numSpecies
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%surf_emiss_out(i1:i2, ju1:j2, numSpecies))
    self%surf_emiss_out = 0.0d0
    return
  end subroutine Allocate_surf_emiss_out
!-------------------------------------------------------------------------
  subroutine Get_surf_emiss_out (self, surf_emiss_out)
    implicit none
    real*8          , intent(out)  :: surf_emiss_out (:,:,:)
    type (t_Emission), intent(in)   :: self
    surf_emiss_out(:,:,:) = self%surf_emiss_out(:,:,:)
    return
  end subroutine Get_surf_emiss_out
!-------------------------------------------------------------------------
  subroutine Set_surf_emiss_out (self, surf_emiss_out)
    implicit none
    real*8          , intent(in)  :: surf_emiss_out (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%surf_emiss_out(:,:,:) = surf_emiss_out(:,:,:)
    return
  end subroutine Set_surf_emiss_out
!-------------------------------------------------------------------------
  subroutine Allocate_surf_emiss_out2 (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%surf_emiss_out2(i1:i2, ju1:j2, 6))
    self%surf_emiss_out2 = 0.0d0
    return
  end subroutine Allocate_surf_emiss_out2
!-------------------------------------------------------------------------
  subroutine Get_surf_emiss_out2 (self, surf_emiss_out2)
    implicit none
    real*8          , intent(out)  :: surf_emiss_out2 (:,:,:)
    type (t_Emission), intent(in)   :: self
    surf_emiss_out2(:,:,:) = self%surf_emiss_out2(:,:,:)
    return
  end subroutine Get_surf_emiss_out2
!-------------------------------------------------------------------------
  subroutine Set_surf_emiss_out2 (self, surf_emiss_out2)
    implicit none
    real*8          , intent(in)  :: surf_emiss_out2 (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%surf_emiss_out2(:,:,:) = surf_emiss_out2(:,:,:)
    return
  end subroutine Set_surf_emiss_out2
!-------------------------------------------------------------------------
  subroutine Allocate_emiss_3d_out (self, i1, i2, ju1, j2, k1, k2, numSpecies)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2, numSpecies
    type (t_Emission), intent(inOut) :: self
    Allocate(self%emiss_3d_out(i1:i2, ju1:j2, k1:k2, numSpecies))
    self%emiss_3d_out = 0.0d0
    return
  end subroutine Allocate_emiss_3d_out
!-------------------------------------------------------------------------
  subroutine Get_emiss_3d_out (self, emiss_3d_out)
    implicit none
    real*8          , intent(out)  :: emiss_3d_out (:,:,:,:)
    type (t_Emission), intent(in)   :: self
    emiss_3d_out(:,:,:,:) = self%emiss_3d_out(:,:,:,:)
    return
  end subroutine Get_emiss_3d_out
!-------------------------------------------------------------------------
  subroutine Set_emiss_3d_out (self, emiss_3d_out)
    implicit none
    real*8          , intent(in)  :: emiss_3d_out (:,:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%emiss_3d_out(:,:,:,:) = emiss_3d_out(:,:,:,:)
    return
  end subroutine Set_emiss_3d_out
!-------------------------------------------------------------------------
  subroutine Allocate_aerosolEmiss3D (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%aerosolEmiss3D(i1:i2, ju1:j2, k1:k2, 5))
    self%aerosolEmiss3D = 0.0d0
    return
  end subroutine Allocate_aerosolEmiss3D
!-------------------------------------------------------------------------
  subroutine Get_aerosolEmiss3D (self, aerosolEmiss3D)
    implicit none
    real*8          , intent(out)  :: aerosolEmiss3D (:,:,:,:)
    type (t_Emission), intent(in)   :: self
    aerosolEmiss3D(:,:,:,:) = self%aerosolEmiss3D(:,:,:,:)
    return
  end subroutine Get_aerosolEmiss3D
!-------------------------------------------------------------------------
  subroutine Set_aerosolEmiss3D (self, aerosolEmiss3D)
    implicit none
    real*8          , intent(in)  :: aerosolEmiss3D (:,:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%aerosolEmiss3D(:,:,:,:) = aerosolEmiss3D(:,:,:,:)
    return
  end subroutine Set_aerosolEmiss3D
!-------------------------------------------------------------------------
  subroutine Allocate_aerosolSurfEmiss (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%aerosolSurfEmiss(i1:i2, ju1:j2, &
        self%ndust + self%naero + 5))
    self%aerosolSurfEmiss = 0.0d0
    return
  end subroutine Allocate_aerosolSurfEmiss
!-------------------------------------------------------------------------
  subroutine Get_aerosolSurfEmiss (self, aerosolSurfEmiss)
    implicit none
    real*8          , intent(out)  :: aerosolSurfEmiss (:,:,:)
    type (t_Emission), intent(in)   :: self
    aerosolSurfEmiss(:,:,:) = self%aerosolSurfEmiss(:,:,:)
    return
  end subroutine Get_aerosolSurfEmiss
!-------------------------------------------------------------------------
  subroutine Set_aerosolSurfEmiss (self, aerosolSurfEmiss)
    implicit none
    real*8          , intent(in)  :: aerosolSurfEmiss (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%aerosolSurfEmiss(:,:,:) = aerosolSurfEmiss(:,:,:)
    return
  end subroutine Set_aerosolSurfEmiss
!-------------------------------------------------------------------------
  subroutine Allocate_aerosolSurfEmissMap (self)
    implicit none
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%aerosolSurfEmissMap &
        (self%ndust + self%naero + 5))
    return
  end subroutine Allocate_aerosolSurfEmissMap
!-------------------------------------------------------------------------
  subroutine Get_aerosolSurfEmissMap (self, aerosolSurfEmissMap)
    implicit none
    integer          , intent(out)  :: aerosolSurfEmissMap (:)
    type (t_Emission), intent(in)   :: self
    aerosolSurfEmissMap(:) = self%aerosolSurfEmissMap(:)
    return
  end subroutine Get_aerosolSurfEmissMap
!-------------------------------------------------------------------------
  subroutine Allocate_flashrate (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%flashrate(i1:i2, ju1:j2))
    self%flashrate = 0.0d0
    return
  end subroutine Allocate_flashrate
!-------------------------------------------------------------------------
  subroutine Get_flashrate (self, flashrate)
    implicit none
    real*8          , intent(out)  :: flashrate (:,:)
    type (t_Emission), intent(in)   :: self
    flashrate(:,:) = self%flashrate(:,:)
    return
  end subroutine Get_flashrate
!-------------------------------------------------------------------------
  subroutine Set_flashrate (self, flashrate)
    implicit none
    real*8          , intent(in)  :: flashrate (:,:)
    type (t_Emission), intent(inOut) :: self
    self%flashrate(:,:) = flashrate(:,:)
    return
  end subroutine Set_flashrate
!-------------------------------------------------------------------------
  subroutine Allocate_emissDust_t (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    integer                          :: ndust, nt_dust, nst_dust
    ndust    = self%ndust
    nst_dust = self%nst_dust
    nt_dust  = self%nt_dust
    Allocate(self%emissDust_t(i1:i2, ju1:j2, 1:ndust, nst_dust:nst_dust+nt_dust-1))
    self%emissDust_t = 0.0d0
    return
  end subroutine Allocate_emissDust_t
!-------------------------------------------------------------------------
  subroutine Get_emissDust_t (self, emiss_dust_t)
    implicit none
    real*8          , intent(out)  :: emiss_dust_t (:,:,:,:)
    type (t_Emission), intent(in)   :: self
    emiss_dust_t(:,:,:,:) = self%emissDust_t(:,:,:,:)
    return
  end subroutine Get_emissDust_t
!-------------------------------------------------------------------------
  subroutine Set_emissDust_t (self, emiss_dust_t)
    implicit none
    real*8          , intent(in)  :: emiss_dust_t (:,:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%emissDust_t(:,:,:,:) = emiss_dust_t(:,:,:,:)
    return
  end subroutine Set_emissDust_t
!-------------------------------------------------------------------------
  subroutine Allocate_emissAero_t (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    integer                          :: naero, emiss_timpyr
    naero        = self%naero
    emiss_timpyr = self%emiss_timpyr
    Allocate(self%emissAero_t(i1:i2, ju1:j2, 1:naero, emiss_timpyr))
    self%emissAero_t = 0.0d0
    return
  end subroutine Allocate_emissAero_t
!-------------------------------------------------------------------------
  subroutine Get_emissAero_t (self, emiss_aero_t)
    implicit none
    real*8          , intent(out)  :: emiss_aero_t (:,:,:,:)
    type (t_Emission), intent(in)   :: self
    emiss_aero_t(:,:,:,:) = self%emissAero_t(:,:,:,:)
    return
  end subroutine Get_emissAero_t
!-------------------------------------------------------------------------
  subroutine Set_emissAero_t (self, emiss_aero_t)
    implicit none
    real*8          , intent(in)  :: emiss_aero_t (:,:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%emissAero_t(:,:,:,:) = emiss_aero_t(:,:,:,:)
    return
  end subroutine Set_emissAero_t
!-------------------------------------------------------------------------
  subroutine Get_emiss_map (self, emiss_map, num_emiss)
    implicit none
    integer        , intent(in )  :: num_emiss
    integer        , intent(out)  :: emiss_map (:)
    type (t_Emission), intent(in)   :: self
    emiss_map(1:num_emiss) = self%emiss_map(1:num_emiss)
    return
  end subroutine Get_emiss_map
!-------------------------------------------------------------------------
  subroutine Get_emiss_timpyr (self, emiss_timpyr)
    implicit none
    integer        , intent(out)  :: emiss_timpyr
    type (t_Emission), intent(in)   :: self
    emiss_timpyr = self%emiss_timpyr
    return
  end subroutine Get_emiss_timpyr
!-------------------------------------------------------------------------
  subroutine Set_emiss_timpyr (self, emiss_timpyr)
    implicit none
    integer        , intent(in)  :: emiss_timpyr
    type (t_Emission), intent(inOut) :: self
    self%emiss_timpyr = emiss_timpyr
    return
  end subroutine Set_emiss_timpyr
!-------------------------------------------------------------------------
  subroutine Get_do_gcr (self, do_gcr)
    implicit none
    logical        , intent(out)  :: do_gcr
    type (t_Emission), intent(in)   :: self
    do_gcr = self%do_gcr
    return
  end subroutine Get_do_gcr
!-------------------------------------------------------------------------
  subroutine Get_do_semiss_inchem (self, do_semiss_inchem)
    implicit none
    logical        , intent(out)  :: do_semiss_inchem
    type (t_Emission), intent(in)   :: self
    do_semiss_inchem = self%do_semiss_inchem
    return
  end subroutine Get_do_semiss_inchem
!-------------------------------------------------------------------------
  subroutine Get_do_ShipEmission (self, do_ShipEmission)
    implicit none
    logical        , intent(out)  :: do_ShipEmission
    type (t_Emission), intent(in)   :: self
    do_ShipEmission = self%do_ShipEmission
    return
  end subroutine Get_do_ShipEmission
!-------------------------------------------------------------------------
  subroutine Get_doMEGANemission (self, doMEGANemission)
    implicit none
    logical        , intent(out)  :: doMEGANemission
    type (t_Emission), intent(in)   :: self
    doMEGANemission = self%doMEGANemission
    return
  end subroutine Get_doMEGANemission
!-------------------------------------------------------------------------
  subroutine Allocate_localCoeff(self, i1, i2, ju1, j2)
    implicit none
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%localCoeff(i1:i2,ju1:j2))
    self%localCoeff = 0.0d0
    return
  end subroutine Allocate_localCoeff
!-------------------------------------------------------------------------
  subroutine Allocate_midLatAdj(self, i1, i2, ju1, j2)
    implicit none
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%midLatAdj(i1:i2,ju1:j2))
    self%midLatAdj = 0.0d0
    return
  end subroutine Allocate_midLatAdj
!-------------------------------------------------------------------------
  subroutine Allocate_isoLaiNext(self, i1, i2, ju1, j2)
    implicit none
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%isoLaiNext(i1:i2,ju1:j2))
    self%isoLaiNext = 0.0d0
    return
  end subroutine Allocate_isoLaiNext
!-------------------------------------------------------------------------
  subroutine Allocate_isoLaiCurr(self, i1, i2, ju1, j2)
    implicit none
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%isoLaiCurr(i1:i2,ju1:j2))
    self%isoLaiCurr = 0.0d0
    return
  end subroutine Allocate_isoLaiCurr
!-------------------------------------------------------------------------
  subroutine Allocate_isoLaiPrev(self, i1, i2, ju1, j2)
    implicit none
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%isoLaiPrev(i1:i2,ju1:j2))
    self%isoLaiPrev = 0.0d0
    return
  end subroutine Allocate_isoLaiPrev
!-------------------------------------------------------------------------
  subroutine Allocate_isoLai(self, i1, i2, ju1, j2)
    implicit none
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%isoLai(i1:i2,ju1:j2))
    self%isoLai = 0.0d0
    return
  end subroutine Allocate_isoLai
!-------------------------------------------------------------------------
  subroutine Get_isoLai(self, isoLai)
    implicit none
    real*8 , intent(out) :: isoLai(:,:)
    type (t_Emission), intent(in)   :: self
    isoLai(:,:) = self%isoLai(:,:)
    return
  end subroutine Get_isoLai
!-------------------------------------------------------------------------
  subroutine Allocate_aefMbo(self, i1, i2, ju1, j2)
    implicit none
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%aefMbo(i1:i2,ju1:j2))
    self%aefMbo = 0.0d0
    return
  end subroutine Allocate_aefMbo
!-------------------------------------------------------------------------
  subroutine Get_aefMbo(self, aefMbo)
    implicit none
    real*8 , intent(out) :: aefMbo(:,:)
    type (t_Emission), intent(in)   :: self
    aefMbo(:,:) = self%aefMbo(:,:)
    return
  end subroutine Get_aefMbo
!-------------------------------------------------------------------------
  subroutine Allocate_aefIsop(self, i1, i2, ju1, j2)
    implicit none
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%aefIsop(i1:i2,ju1:j2))
    self%aefIsop = 0.0d0
    return
  end subroutine Allocate_aefIsop
!-------------------------------------------------------------------------
  subroutine Get_aefIsop(self, aefIsop)
    implicit none
    real*8 , intent(out) :: aefIsop(:,:)
    type (t_Emission), intent(in)   :: self
    aefIsop(:,:) = self%aefIsop(:,:)
    return
  end subroutine Get_aefIsop
!-------------------------------------------------------------------------
  subroutine Allocate_aefOvoc(self, i1, i2, ju1, j2)
    implicit none
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%aefOvoc(i1:i2,ju1:j2))
    self%aefOvoc = 0.0d0
    return
  end subroutine Allocate_aefOvoc
!-------------------------------------------------------------------------
  subroutine Get_aefOvoc(self, aefOvoc)
    implicit none
    real*8 , intent(out) :: aefOvoc(:,:)
    type (t_Emission), intent(in)   :: self
    aefOvoc(:,:) = self%aefOvoc(:,:)
    return
  end subroutine Get_aefOvoc
!-------------------------------------------------------------------------
  subroutine Allocate_aefMonot(self, i1, i2, ju1, j2)
    implicit none
    integer, intent(in) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut)   :: self
    allocate(self%aefMonot(i1:i2,ju1:j2))
    self%aefMonot = 0.0d0
    return
  end subroutine Allocate_aefMonot
!-------------------------------------------------------------------------
  subroutine Get_aefMonot(self, aefMonot)
    implicit none
    real*8 , intent(out) :: aefMonot(:,:)
    type (t_Emission), intent(in)   :: self
    aefMonot(:,:) = self%aefMonot(:,:)
    return
  end subroutine Get_aefMonot
!-------------------------------------------------------------------------
  subroutine Get_emiss_in_opt (self, emiss_in_opt)
    implicit none
    integer        , intent(out)  :: emiss_in_opt
    type (t_Emission), intent(in)   :: self
    emiss_in_opt = self%emiss_in_opt
    return
  end subroutine Get_emiss_in_opt
!-------------------------------------------------------------------------
  subroutine Get_emiss_opt (self, emiss_opt)
    implicit none
    integer        , intent(out)  :: emiss_opt
    type (t_Emission), intent(in)   :: self
    emiss_opt = self%emiss_opt
    return
  end subroutine Get_emiss_opt
!-------------------------------------------------------------------------
  subroutine Get_curEmissionFileRecord (self, curEmissionFileRecord)
    implicit none
    integer        , intent(out)  :: curEmissionFileRecord
    type (t_Emission), intent(in)   :: self
    curEmissionFileRecord = self%curEmissionFileRecord
    return
  end subroutine Get_curEmissionFileRecord
!-------------------------------------------------------------------------
  subroutine Get_lightning_opt (self, lightning_opt)
    implicit none
    integer        , intent(out)  :: lightning_opt
    type (t_Emission), intent(in)   :: self
    lightning_opt = self%lightning_opt
    return
  end subroutine Get_lightning_opt
!-------------------------------------------------------------------------
  subroutine Get_emiss_aero_opt (self, emiss_aero_opt)
    implicit none
    integer        , intent(out)  :: emiss_aero_opt
    type (t_Emission), intent(in)   :: self
    emiss_aero_opt = self%emiss_aero_opt
    return
  end subroutine Get_emiss_aero_opt
!-------------------------------------------------------------------------
  subroutine Get_emiss_dust_opt (self, emiss_dust_opt)
    implicit none
    integer        , intent(out)  :: emiss_dust_opt
    type (t_Emission), intent(in)   :: self
    emiss_dust_opt = self%emiss_dust_opt
    return
  end subroutine Get_emiss_dust_opt
!-------------------------------------------------------------------------
  subroutine Get_num_emiss (self, num_emiss)
    implicit none
    integer        , intent(out)  :: num_emiss
    type (t_Emission), intent(in)   :: self
    num_emiss = self%num_emiss
    return
  end subroutine Get_num_emiss
!-------------------------------------------------------------------------
  subroutine Get_num_emiss3d (self, num_emiss3d)
    implicit none
    integer        , intent(out)  :: num_emiss3d
    type (t_Emission), intent(in)   :: self
    num_emiss3d = self%num_emiss3d
    return
  end subroutine Get_num_emiss3d
!-------------------------------------------------------------------------
  subroutine Get_o3_index (self, o3_index)
    implicit none
    integer        , intent(out)  :: o3_index
    type (t_Emission), intent(in)   :: self
    o3_index = self%o3_index
    return
  end subroutine Get_o3_index
!-------------------------------------------------------------------------
  subroutine Get_hno3_index (self, hno3_index)
    implicit none
    integer        , intent(out)  :: hno3_index
    type (t_Emission), intent(in)   :: self
    hno3_index = self%hno3_index
    return
  end subroutine Get_hno3_index
!-------------------------------------------------------------------------
  subroutine Get_ndust (self, ndust)
    implicit none
    integer        , intent(out)  :: ndust
    type (t_Emission), intent(in)   :: self
    ndust = self%ndust
    return
  end subroutine Get_ndust
!-------------------------------------------------------------------------
  subroutine Get_naero (self, naero)
    implicit none
    integer        , intent(out)  :: naero
    type (t_Emission), intent(in)   :: self
    naero = self%naero
    return
  end subroutine Get_naero
!-------------------------------------------------------------------------
  subroutine Get_nst_dust (self, nst_dust)
    implicit none
    integer        , intent(out)  :: nst_dust
    type (t_Emission), intent(in)   :: self
    nst_dust = self%nst_dust
    return
  end subroutine Get_nst_dust
!-------------------------------------------------------------------------
  subroutine Get_nt_dust (self, nt_dust)
    implicit none
    integer        , intent(out)  :: nt_dust
    type (t_Emission), intent(in)   :: self
    nt_dust = self%nt_dust
    return
  end subroutine Get_nt_dust
!-------------------------------------------------------------------------
  subroutine Get_emiss_map_dust (self, emiss_map_dust, ndust)
    implicit none
    integer        , intent(in )  :: ndust
    integer        , intent(out)  :: emiss_map_dust(:)
    type (t_Emission), intent(in)   :: self
    emiss_map_dust(1:ndust) = self%emiss_map_dust(1:ndust)
    return
  end subroutine Get_emiss_map_dust
!-------------------------------------------------------------------------
  subroutine Get_emiss_map_aero (self, emiss_map_aero, naero)
    implicit none
    integer        , intent(in )  :: naero
    integer        , intent(out)  :: emiss_map_aero(:)
    type (t_Emission), intent(in)   :: self
    emiss_map_aero(1:naero) = self%emiss_map_aero(1:naero)
    return
  end subroutine Get_emiss_map_aero
!-------------------------------------------------------------------------
  subroutine Allocate_scFacNOff (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    allocate(self%scFacNOff(i1:i2,ju1:j2,1:24))
    self%scFacNOff(:,:,:) = 1.0d0
    return
  end subroutine Allocate_scFacNOff
!-------------------------------------------------------------------------
  subroutine Get_scFacNOff (self, scFacNOff)
    implicit none
    real*8, intent(out)  :: scFacNOff(:,:,:)
    type (t_Emission)    , intent(in)   :: self
    scFacNOff(:,:,:) = self%scFacNOff(:,:,:)
    return
  end subroutine Get_scFacNOff
!-------------------------------------------------------------------------
  subroutine Get_scFactorNOff_infile_name (self, scFactorNOff_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: scFactorNOff_infile_name
    type (t_Emission)    , intent(in)   :: self
    scFactorNOff_infile_name = self%scFactorNOff_infile_name
    return
  end subroutine Get_scFactorNOff_infile_name
!-------------------------------------------------------------------------
  subroutine Allocate_scFacNObb (self, i1, i2, ju1, j2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2
    type (t_Emission), intent(inOut) :: self
    allocate(self%scFacNObb(i1:i2,ju1:j2,1:24))
    self%scFacNObb(:,:,:) = 1.0d0
    return
  end subroutine Allocate_scFacNObb
!-------------------------------------------------------------------------
  subroutine Get_scFacNObb (self, scFacNObb)
    implicit none
    real*8, intent(out)  :: scFacNObb(:,:,:)
    type (t_Emission)    , intent(in)   :: self
    scFacNObb(:,:,:) = self%scFacNObb(:,:,:)
    return
  end subroutine Get_scFacNObb
!-------------------------------------------------------------------------
  subroutine Get_scFactorNObb_infile_name (self, scFactorNObb_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: scFactorNObb_infile_name
    type (t_Emission)    , intent(in)   :: self
    scFactorNObb_infile_name = self%scFactorNObb_infile_name
    return
  end subroutine Get_scFactorNObb_infile_name
!-------------------------------------------------------------------------
  subroutine Get_fertscal_infile_name (self, fertscal_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: fertscal_infile_name
    type (t_Emission)    , intent(in)   :: self
    fertscal_infile_name = self%fertscal_infile_name
    return
  end subroutine Get_fertscal_infile_name
!-------------------------------------------------------------------------
  subroutine Get_lai_infile_name (self, lai_infile_name)
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: lai_infile_name
    type (t_Emission)    , intent(in)   :: self
    lai_infile_name = self%lai_infile_name
    return
  end subroutine Get_lai_infile_name
!-------------------------------------------------------------------------
  subroutine Get_light_infile_name (self, light_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: light_infile_name
    type (t_Emission)    , intent(in)   :: self
    light_infile_name = self%light_infile_name
    return
  end subroutine Get_light_infile_name
!-------------------------------------------------------------------------
  subroutine Get_precip_infile_name (self, precip_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: precip_infile_name
    type (t_Emission)    , intent(in)   :: self
    precip_infile_name = self%precip_infile_name
    return
  end subroutine Get_precip_infile_name
!-------------------------------------------------------------------------
  subroutine Get_soil_infile_name (self, soil_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: soil_infile_name
    type (t_Emission)    , intent(in)   :: self
    soil_infile_name = self%soil_infile_name
    return
  end subroutine Get_soil_infile_name
!-------------------------------------------------------------------------
  subroutine Get_veg_infile_name (self, veg_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: veg_infile_name
    type (t_Emission)    , intent(in)   :: self
    veg_infile_name = self%veg_infile_name
    return
  end subroutine Get_veg_infile_name
!-------------------------------------------------------------------------
  subroutine Get_isopconv_infile_name (self, isopconv_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: isopconv_infile_name
    type (t_Emission)    , intent(in)   :: self
    isopconv_infile_name = self%isopconv_infile_name
    return
  end subroutine Get_isopconv_infile_name
!-------------------------------------------------------------------------
  subroutine Get_monotconv_infile_name (self, monotconv_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: monotconv_infile_name
    type (t_Emission)    , intent(in)   :: self
    monotconv_infile_name = self%monotconv_infile_name
    return
  end subroutine Get_monotconv_infile_name
!-------------------------------------------------------------------------
  subroutine Get_GOCARTerod_infile_name (self, GOCARTerod_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: GOCARTerod_infile_name
    type (t_Emission)    , intent(in)   :: self
    GOCARTerod_infile_name = self%GOCARTerod_infile_name
    return
  end subroutine Get_GOCARTerod_infile_name
!-------------------------------------------------------------------------
  subroutine Get_GOCARTocean_infile_name (self, GOCARTocean_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: GOCARTocean_infile_name
    type (t_Emission)    , intent(in)   :: self
    GOCARTocean_infile_name = self%GOCARTocean_infile_name
    return
  end subroutine Get_GOCARTocean_infile_name
!-------------------------------------------------------------------------
  subroutine Get_GOCARTerod_mod_infile_name (self, GOCARTerod_mod_infile_name)
    implicit none
    character (len=MAX_LENGTH_FILE_NAME), intent(out)  :: GOCARTerod_mod_infile_name
    type (t_Emission)    , intent(in)   :: self
    GOCARTerod_mod_infile_name = self%GOCARTerod_mod_infile_name
    return
  end subroutine Get_GOCARTerod_mod_infile_name
!-------------------------------------------------------------------------
  subroutine Get_desired_g_N_prod_rate (self, desired_g_N_prod_rate)
    implicit none
    real*8         , intent(out)  :: desired_g_N_prod_rate
    type (t_Emission), intent(in)   :: self
    desired_g_N_prod_rate = self%desired_g_N_prod_rate
    return
  end subroutine Get_desired_g_N_prod_rate
!-------------------------------------------------------------------------
  subroutine Get_isop_scale (self, isop_scale)
    implicit none
    real*8         , intent(out)  :: isop_scale(:)
    type (t_Emission), intent(in)   :: self
    isop_scale(:) = self%isop_scale(:)
    return
  end subroutine Get_isop_scale
!-------------------------------------------------------------------------
  subroutine Set_isop_scale (self, isop_scale)
    implicit none
    real*8         , intent(in   )  :: isop_scale(:)
    type (t_Emission), intent(inOut)  :: self
    self%isop_scale(:) = isop_scale(:)
    return
  end subroutine Set_isop_scale
!-------------------------------------------------------------------------
  subroutine Get_begDailyEmissRec (self, begDailyEmissRec)
    implicit none
    integer          , intent(out)  :: begDailyEmissRec
    type (t_Emission), intent(in)   :: self
    begDailyEmissRec = self%begDailyEmissRec
    return
  end subroutine Get_begDailyEmissRec
!-------------------------------------------------------------------------
  subroutine Get_endDailyEmissRec (self, endDailyEmissRec)
    implicit none
    integer          , intent(out)  :: endDailyEmissRec
    type (t_Emission), intent(in)   :: self
    endDailyEmissRec = self%endDailyEmissRec
    return
  end subroutine Get_endDailyEmissRec
!-------------------------------------------------------------------------
  subroutine Get_doReadDailyEmiss (self, doReadDailyEmiss)
    implicit none
    logical          , intent(out)  :: doReadDailyEmiss
    type (t_Emission), intent(in)   :: self
    doReadDailyEmiss = self%doReadDailyEmiss
    return
  end subroutine Get_doReadDailyEmiss
!-------------------------------------------------------------------------
   subroutine Allocate_gcr_slope (self, ju1, j2, k1, k2)
    implicit none
    integer     , intent(in)    :: ju1, j2, k1, k2
    type (t_Emission), intent(inOut) :: self
    allocate(self%gcr_slope(ju1:j2, k1:k2))
    self%gcr_slope(:,:) = 0.0d0
    return
  end subroutine Allocate_gcr_slope
!-------------------------------------------------------------------------
   subroutine Allocate_gcr_aintcp (self, ju1, j2, k1, k2)
    implicit none
    integer     , intent(in)    :: ju1, j2, k1, k2
    type (t_Emission), intent(inOut) :: self
    allocate(self%gcr_aintcp(ju1:j2, k1:k2))
    self%gcr_aintcp(:,:) = 0.0d0
    return
  end subroutine Allocate_gcr_aintcp
!-------------------------------------------------------------------------
  subroutine Allocate_GCR_NOx (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer          , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_Emission)      , intent(inOut) :: self
    Allocate(self%GCR_NOx(i1:i2, ju1:j2, k1:k2))
    self%GCR_NOx = 0.0d0
    return
  end subroutine Allocate_GCR_NOx
!-------------------------------------------------------------------------
  subroutine Get_GCR_NOx (self, GCR_NOx)
    implicit none
    real*8          , intent(out)  :: GCR_NOx (:,:,:)
    type (t_Emission), intent(in)   :: self
    GCR_NOx(:,:,:) = self%GCR_NOx(:,:,:)
    return
  end subroutine Get_GCR_NOx
!-------------------------------------------------------------------------
  subroutine Set_GCR_NOx (self, GCR_NOx)
    implicit none
    real*8          , intent(in)  :: GCR_NOx (:,:,:)
    type (t_Emission), intent(inOut) :: self
    self%GCR_NOx(:,:,:) = GCR_NOx(:,:,:) 
    return
  end subroutine Set_GCR_NOx
!-------------------------------------------------------------------------
  end module GmiEmissionMethod_mod
!
