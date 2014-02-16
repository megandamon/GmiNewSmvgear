!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MUDULE: GmiTimeStepping_mod
!
#include "GmiESMF_ErrLog.h"
!
      module GmiTimeStepping_mod
!
! !USES:
      use Ftiming_Dao
      use ESMF_Mod
      use GmiESMF_ErrorChecking_mod
      use GmiSurfaceTemperature_mod, only : updateSurfaceTempMEGAN
      use GmiAdvecCoreWrapper_mod, only : fromGmiToAdvecCore, fromAdvecCoreToGmi
      use GmiESMFderivedType_mod, only : t_gmiESMF
      use GmiReduce_mod,             only : Gmi_Min1_Reduce
      use GmiDiagnosticsMethod_mod, only : t_Diagnostics, Get_do_ftiming,      &
     &       Get_pr_diag, Get_num_const_outrecs, Get_num_tend_outrecs,         &
     &       Get_do_qqjk_reset
      use GmiAdvectionMethod_mod, only : runAdvection, t_Advection,            &
     &       Get_advec_opt, Get_do_grav_set
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_numDomains,      &
     &       Get_communicatorWorld, Get_mcor, Get_mcorGlob, Get_latdeg, Get_londeg, &
     &       Get_procID, Get_dlatr, Get_cosp, Get_coscen
      use GmiMessagePassing_mod        , only : synchronizeGroup
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration,        &
     &       resetFixedConcentration, Get_concentration, Set_concentration,    &
     &       calcConcentrationSurf, calcConcentrationColTrop, Get_tracer_opt,  &
     &       calcConcentrationColCombo
      use GmiTaggedCO_AgeOfAir_mod, only : calcTaggedCO_AgeOfAir
      use GmiSolar_mod, only : CalcCosSolarZenithAngle
      use GmiPressureFixer_mod, only : Calc_Delpm
      use GmiUpdateSyntheticSpecies_mod, only : updateSyntheticSpecies
      use GmiChemistryMethod_mod , only : t_Chemistry, runChemistry,           &
     &       runReadChemistry, runCalcAerosolDust, Get_optDepth, Set_optDepth,    &
     &       Get_loss_opt,                                                     &
     &       Get_oz_eq_synoz_opt, Get_num_qjs, Get_num_qjo, Get_num_qks,       &
     &       Get_num_chem, Get_num_active, Get_num_molefrac, Get_num_sad,      &
     &       Get_num_ks_sbc, Get_num_spc_sbc, Get_surf_bc_map, Get_sad_opt,    &
     &       Get_chem_opt, Get_phot_opt, Get_do_chem_grp, Get_mw,              &
     &       Get_iisoprene_num, Get_ino_num, Get_ico_num, Get_ipropene_num,    &
     &       Get_ihno3_num, Get_io3_num, Get_ih2o2_num, Get_ih2o_num,          &
     &       Get_isynoz_num, Get_do_synoz, Get_do_nodoz, Get_num_nox,          &
     &       Get_num_noy, Get_noy_map, Get_nox_map, Get_ihcl_num,              &
     &       Get_idehyd_num, Get_dehydmin, Set_dehydmin, Get_ibrono2_num,      &
     &       Get_t_cloud_ice, Get_synoz_threshold, Get_do_AerDust_Calc,        &
     &       Get_imgas_num, Get_initrogen_num, Get_ioxygen_num, Get_do_clear_sky
      use GmiDepositionMethod_mod, only : t_Deposition, runDryDeposition,      &
     &       runWetDeposition, runSimpleDeposition, Get_do_simpledep,          &
     &       Get_do_wetdep, Get_do_drydep
!
      use computeMaxMinArray_mod
      use GmiConvectionMethod_mod, only : t_Convection, runConvection, &
     &       Get_convec_opt
      use GmiDiffusionMethod_mod , only : t_Diffusion, runDiffusion, &
     &       Get_diffu_opt
!
      use GmiESMFclock_mod       , only : advanceESMFclock
      use GmiTimeControl_mod, only : t_GmiClock, GetSecondsFromJanuary1,       &
     &       Get_curGmiDate, Get_curGmiTime, Get_gmiSeconds, Get_numTimeSteps, &
     &       Get_gmiTimeStep
      use GmiTotalMass_mod       , only : calcTotalMass
      use GmiRelativeHumidity_mod, only : CalcRelativeHumidity
      use GmiHeightSigmaLevel_mod, only : CalcConvection_bmass
      use GmiHeightSigmaLevel_mod, only : CalcDryDepBoxHeight
      use GmiGridBoxHeight_mod   , only : CalcGridBoxHeight
      use GmiCloudVariables_mod, only : CalcFractionalCloudCover,              &
     &       CalcCloudOpticalDepth, CalcTotalCloudFraction,                    &
     &       DiagCloudOpticalDepth
      use GmiCalcMoisture_mod    , only : CalcMoistureChanges
      use GmiPressure_mod, only : CalcPress3dCenter, CalcPress3dEdge,          &
     &       CalcAveragePressEdge, calcTropopausePress_Stobie,                 &
     &       calcTropopausePress_Ertel, calcTropopausePress_WMO
!
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiTroposphericWater_mod , only : addTroposphericWater
      use GmiTroposphericWater_mod , only : removeTroposphericWater
!
      use GmiEmissionMethod_mod  , only : t_Emission, runEmission, &
     &       runReadEmission, Get_emiss_opt, Get_do_gcr, Get_doMEGANemission
      use GmiGrid_mod, only : t_gmiGrid, Get_numSpecies, Get_i1, Get_i2,       &
     &       Get_ju1, Get_j2, Get_k1, Get_k2, Get_ilo, Get_ihi, Get_julo,      &
     &       Get_jhi, Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl, Get_jvlo,   &
     &       Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jhi_gl
!
      use GmiProdLossDiagnostics_mod, only : doDiagnosticsBefore,              &
     &       doDiagnosticsAfter
!
      use GmiMassFluxes_mod, only : convertMassFlux
      use GmiGravitationalSettling_mod, only : updateGravitationalSettling
!
      use GmiMetFieldsControl_mod, only : t_metFields, Get_pt, Get_met_opt,    &
     &       Get_ai, Get_bi, Get_am, Get_bm, Get_dap, Get_dbk, finishMet,      &
     &       Get_metdata_name_org, Get_metdata_name_model,                     &
     &       Set_mass, Set_tropopausePress, Set_gridBoxHeight,                 &
     &       Set_potentialVorticity, Set_potentialTemp, Set_relativeHumidity,  &
     &       Set_press3e, Set_press3c, Set_tau_cloud, Set_fracCloudCover,      &
     &       Set_totalCloudFraction, Get_humidity, Get_pctm2, Get_pctm1,       &
     &       Get_pctm1Glob, Get_kelGlob, &
     &       Get_max_cloud, Get_ran_cloud, Get_xmass, Get_ymass, Get_zmass,    &
     &       Get_kel, Get_uux, Get_vvx, Get_rain, Get_moistq, Set_moistq,      &
     &       Get_tau_cloud
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: gmiTimeStepping
!
! !DESCRIPTION:
! Time stepping routine.
!
! !AUTHOR:
! JUles Kouatchou, NASA/GSFC, Jules.Kouatchou@nasa.gov
!
!EOP
!------------------------------------------------------------------------------
      contains
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gmiTimeStepping
!
! !INTERFACE:
!
      subroutine gmiTimeStepping (advCoreESMF, esmfClock, &
     &              gmiClock, gmiGrid, gmiDomain, Chemistry, Emission, &
     &              Deposition, Convection, Diffusion, Advection, &
     &              SpeciesConcentration, Diagnostics, metFields, new_met_rec, &
     &              chem_mecha)
!
      implicit none
!
#     include "gmi_diag_constants_llnl.h"
#     include "GmiParameters.h"
#     include "gmi_time_constants.h"
#     include "gmi_AerDust_const.h"
#     include "setkin_par.h"
#     include "gmi_subdomains.h"
#     include "gmi_phys_constants.h"
!
!
! !INPUT PARAMETERS:
      character(len=*)   , intent(in) :: chem_mecha
      logical            , intent(in) :: new_met_rec ! new met record?
      type (t_gmiGrid   ), intent(in) :: gmiGrid
      type (t_gmiDomain ), intent(in) :: gmiDomain
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_Diagnostics), intent(inOut) :: Diagnostics
      type (t_gmiESMF)             , intent(inOut) :: advCoreESMF
      type (ESMF_Clock            ), intent(inOut) :: esmfClock
      type (t_GmiClock            ), intent(inOut) :: gmiClock
      type (t_Chemistry           ), intent(inOut) :: Chemistry
      type (t_Advection           ), intent(inOut) :: Advection
      type (t_Emission            ), intent(inOut) :: Emission
      type (t_Deposition          ), intent(inOut) :: Deposition
      type (t_Convection          ), intent(inOut) :: Convection
      type (t_Diffusion           ), intent(inOut) :: Diffusion
      type (t_metFields           ), intent(inOut) :: metFields
      type (t_SpeciesConcentration), intent(inout) :: SpeciesConcentration
!
! !DESCRIPTION:
!  Steps the simulation one time step.
!
! !DEFINED PARAMETERS:
      real(ESMF_KIND_R4), parameter :: mbToPascal = 100.0
!
! !LOCAL VARIABLES:
      logical, save :: first = .true.
!
      integer :: ic, ik
      integer, save :: curRecord = -1
!
      ! -----------------------------------------------------
      ! diffaer    : diffusivity of aerosol at surface layer.
      ! s_radius   : radius of aerosol at surface layer.
      ! s_velocity : settling velocity at surface layer.
      ! -----------------------------------------------------
!
      real*8, allocatable  :: diffaer   (:,:,:)
      real*8, allocatable  :: s_radius  (:,:,:)
      real*8, allocatable  :: s_velocity(:,:,:)
!
      ! -----------------------------------------------------------------
      ! pchem_water : stratospheric and tropospheric water vapor before
      !               chemistry
      ! strat_water : stratospheric water vapor before chemistry operator
      ! -----------------------------------------------------------------
!
      real*8, allocatable :: pchem_water(:,:,:)
      real*8, allocatable :: strat_water(:,:,:)
!
      real*8, allocatable, save :: averagePressEdge(:)
!
      ! ---------------------------------------------
      ! surf_bc : surface boundary condition (units?)
      ! ---------------------------------------------
!
      real*8, allocatable, save :: surf_bc(:,:,:,:)
!
      real*8, allocatable :: bmass(:,:,:) ! for DAO convection
      real*8, allocatable :: BoxHeightCenter(:,:) ! for dry deposition
      real*8, allocatable :: BoxHeightEdge  (:,:) ! for dry deposition
!
      real*8 loss_level, strat_life, strat_life1, strat_life100
      integer ibr, jbr, kbr
      integer llife(17), lexclude(4)
!
      ! ------------------------------------------------------
      ! REL_SCAV_EFF_new : aerosol scavenging efficiency (0-1)
      ! ------------------------------------------------------
!
      real*8, allocatable :: REL_SCAV_EFF_new(:,:,:,:)
!
      real*8  :: pmin, pmax
      real*8  :: days, time
      integer :: nsec_jan1
      real*8, allocatable :: cosSolarZenithAngle(:,:)
      real*8, allocatable :: OptDepth(:,:,:,:)
!
      integer  :: il, jl, ierr, jc
      integer  :: communicatorWorld, numDomains
      integer  :: i1, i2, ju1, j2, k1, k2, ivert
      integer  :: ilo, ihi, julo, jvlo, jhi
      integer  :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer  :: ilo_gl, ihi_gl, julo_gl, jhi_gl
      integer  :: procID, numSpecies, STATUS, rc
      logical  :: do_ftiming, pr_diag
      integer  :: num_const_outrecs, num_tend_outrecs
      type (t_GmiArrayBundle), pointer :: concentration(:)
      real*8 , allocatable :: mcor(:,:)
      real*8 , allocatable :: mcorGlob(:,:)
      real*8 , allocatable :: londeg(:), latdeg(:)
      real*8 , allocatable :: dlatr(:), coscen(:), cosp(:)
      integer :: num_qks, num_qjs, num_qjo
      integer :: num_sad, num_molefrac, num_chem, num_active
!
      real*8 , allocatable :: ai(:), bi(:)
      real*8 , allocatable :: am(:), bm(:)
      real*8 , allocatable :: dap(:), dbk(:)
!
      real*8            , pointer :: delpm(:,:,:)
      real(ESMF_KIND_R4), pointer :: UC(:,:,:), VC(:,:,:)
      real(ESMF_KIND_R4), pointer :: DPEDT(:,:,:), DP(:,:,:)
      real(ESMF_KIND_R4), pointer :: MX_UR(:,:,:), MY_UR(:,:,:)
      real(ESMF_KIND_R4), pointer :: MX(:,:,:), MY(:,:,:), MZ(:,:,:)
!
      integer :: diffu_opt
      logical :: do_qqjk_reset
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
!
      integer :: emiss_opt, chem_opt, sad_opt, phot_opt, advec_opt
      integer :: convec_opt, tracer_opt, loss_opt, oz_eq_synoz_opt, met_opt
      logical :: do_gcr, doMEGANemission, do_chem_grp, do_clear_sky
      logical :: do_grav_set, do_drydep, do_wetdep, do_simpledep
      logical :: do_synoz, do_nodoz, do_AerDust_Calc
!
      integer :: nox_map(MAX_NUM_SMARRAY), noy_map(MAX_NUM_SMARRAY)
      integer, allocatable :: surf_bc_map(:)
!
      real*8 , allocatable :: mw(:)
      real*8  :: dehydmin, gmi_sec, tdt, pt
      integer :: nhms, nymd, num_time_steps
      real*8  :: t_cloud_ice, synoz_threshold
      real*8  :: total_mass, layer1_mass
!
      integer :: ino_num, iisoprene_num, ipropene_num, ih2o_num , ihcl_num
      integer :: io3_num, ioxygen_num  , isynoz_num  , ih2o2_num, ibrono2_num
      integer :: ico_num, initrogen_num, idehyd_num  , imgas_num, ihno3_num
      integer :: num_nox, num_spc_sbc  , num_ks_sbc  , num_noy
!
      real*8, allocatable :: press3c(:,:,:),             press3e(:,:,:)
      real*8, allocatable :: mass (:,:,:),               tropopausePress(:,:)
      real*8, allocatable :: potentialVorticity(:,:,:),  potentialTemp(:,:,:)
      real*8, allocatable :: relativeHumidity(:,:,:),    gridBoxHeight(:,:,:)
!
      real*8, allocatable :: pctm1(:,:),                pctm2(:,:)
      real*8, allocatable :: pctm1Glob(:,:)
      real*8, allocatable :: kelGlob(:,:,:)
      real*8, allocatable :: uux(:,:,:),                vvx(:,:,:)    
      real*8, allocatable :: kel(:,:,:),                humidity(:,:,:)
      real*8, allocatable :: max_cloud(:,:,:),          ran_cloud(:,:,:)
      real*8, allocatable :: tau_cloud(:,:,:),          fracCloudCover(:,:)
      real*8, allocatable :: moistq(:,:,:),             rain(:,:,:)
      real*8, allocatable :: totalCloudFraction(:,:,:), xmass(:,:,:)
      real*8, allocatable :: ymass(:,:,:),              zmass(:,:,:)
!
      type(ESMF_Grid) :: esmfGrid
      character(len=ESMF_MAXSTR), parameter :: IAm = "gmiTimeStepping"
!EOP
!--------------------------------------------------------------------------------
!BOC
      call Get_procID(gmiDomain, procID)
      call Get_pr_diag   (Diagnostics, pr_diag   )
!
      if (pr_diag) write(6,*) trim(Iam), " called by ", procID
!
      call Get_do_ftiming(Diagnostics, do_ftiming)
!
      if (do_ftiming) call Ftiming_On  ('gmiTimeStepping')
!
      call Get_communicatorWorld (gmiDomain, communicatorWorld)
!
      if (do_ftiming) then
        call Ftiming_On  ('procSyncBegStepping')
        call synchronizeGroup(communicatorWorld)
        call Ftiming_Off ('procSyncBegStepping')
      end if
!
      ! Get grid/domain related variables
!
      call Get_i1         (gmiGrid, i1)
      call Get_i2         (gmiGrid, i2)
      call Get_ju1        (gmiGrid, ju1)
      call Get_j2         (gmiGrid, j2)
      call Get_k1         (gmiGrid, k1)
      call Get_k2         (gmiGrid, k2)
      call Get_ilo        (gmiGrid, ilo)
      call Get_ihi        (gmiGrid, ihi)
      call Get_julo       (gmiGrid, julo)
      call Get_jvlo       (gmiGrid, jvlo)
      call Get_jhi        (gmiGrid, jhi)
      call Get_i1_gl      (gmiGrid, i1_gl)
      call Get_i2_gl      (gmiGrid, i2_gl)
      call Get_ju1_gl     (gmiGrid, ju1_gl)
      call Get_j2_gl      (gmiGrid, j2_gl)
      call Get_ilo_gl     (gmiGrid, ilo_gl)
      call Get_ihi_gl     (gmiGrid, ihi_gl)
      call Get_julo_gl    (gmiGrid, julo_gl)
      call Get_jhi_gl     (gmiGrid, jhi_gl)
      call Get_numSpecies (gmiGrid, numSpecies)
      ivert = k2 - k1 + 1
!
      allocate(dlatr   (ju1_gl:j2_gl))
      allocate(coscen  (ju1_gl:j2_gl))
      allocate(cosp    (ju1_gl:j2_gl))
      allocate(latdeg  (ju1_gl:j2_gl))
      allocate(londeg  (i1_gl :i2_gl))
      allocate(mcor    (i1:i2,ju1:j2))
      allocate(mcorGlob(i1_gl:i2_gl,ju1_gl:j2_gl))
!
      call Get_dlatr      (gmiDomain, dlatr)
      call Get_coscen     (gmiDomain, coscen)
      call Get_cosp       (gmiDomain, cosp)
      call Get_latdeg     (gmiDomain, latdeg)
      call Get_londeg     (gmiDomain, londeg)
      call Get_mcor       (gmiDomain, mcor)
      call Get_mcorGlob   (gmiDomain, mcorGlob)
      call Get_numDomains (gmiDomain, numDomains)
!
      ! Get the species concentration
!      
      call Get_concentration(SpeciesConcentration, concentration)
!
      call Get_tracer_opt   (SpeciesConcentration, tracer_opt)
!
      ! Get metFields related variables
!
      allocate(ai(k1-1:k2))
      allocate(bi(k1-1:k2))
      allocate(am(k1:k2))
      allocate(bm(k1:k2))
      allocate(dap(k1:k2))
      allocate(dbk(k1:k2))
!
      call Get_pt                (metFields, pt)
      call Get_ai                (metFields, ai)
      call Get_bi                (metFields, bi)
      call Get_am                (metFields, am)
      call Get_bm                (metFields, bm)
      call Get_dap               (metFields, dap)
      call Get_dbk               (metFields, dbk)
      call Get_met_opt           (metFields, met_opt)
      call Get_metdata_name_org  (metFields, metdata_name_org)
      call Get_metdata_name_model(metFields, metdata_name_model)
!
      call Get_chem_opt        (Chemistry, chem_opt)
      call Get_idehyd_num      (Chemistry, idehyd_num)
      call Get_num_ks_sbc      (Chemistry, num_ks_sbc )
      call Get_num_spc_sbc     (Chemistry, num_spc_sbc)
!
      if ((chem_opt ==2) .and. (num_ks_sbc > 0) .and. (num_spc_sbc > 0)) then
         allocate(surf_bc_map(1:num_spc_sbc))
         call Get_surf_bc_map(Chemistry, surf_bc_map)
      end if
!
!     ==========
      if (first) then
!     ==========
!
         first = .false.
!
         if (chem_opt == 3 .or. chem_opt == 9) then
             allocate(averagePressEdge(k1-1:k2))
             averagePressEdge = 0.0d0
             call CalcAveragePressEdge (ai, bi, pt, averagePressEdge, k1, k2)
         end if
!
         if ((chem_opt ==2) .and. (num_ks_sbc > 0) .and. (num_spc_sbc > 0)) then
            Allocate (surf_bc(i1:i2, ju1:j2, 1:num_ks_sbc, 1:num_spc_sbc))
            surf_bc = 0.0d0
!
            do ic = 1, num_spc_sbc
               do ik = 1, num_ks_sbc
                  surf_bc(:,:,ik,ic) = &
     &                 concentration(surf_bc_map(ic))%pArray3D(:,:,ik)
               end do
            end do
         end if
!
         if (idehyd_num /= 0) then
!
            dehydmin = Minval (concentration(idehyd_num)%pArray3D(:,:,:))
            dehydmin = Min    (0.0d0, dehydmin)
!
            ! =======================
            call Gmi_Min1_Reduce (dehydmin, numDomains, communicatorWorld)
            ! =======================
!
            call synchronizeGroup(communicatorWorld)
!
            call Set_dehydmin (Chemistry, dehydmin)
!
         end if
!
!     ======
      end if
!     ======
!
      call Get_do_gcr          (Emission, do_gcr   )
      call Get_emiss_opt       (Emission, emiss_opt)
      call Get_doMEGANemission (Emission, doMEGANemission)
!
      ! Get Chemistry related variables
!
      call Get_num_qjs         (Chemistry, num_qjs     )
      call Get_num_qjo         (Chemistry, num_qjo     )
      call Get_num_qks         (Chemistry, num_qks     )
      call Get_num_sad         (Chemistry, num_sad     )
      call Get_num_chem        (Chemistry, num_chem    )
      call Get_num_active      (Chemistry, num_active  )
      call Get_num_molefrac    (Chemistry, num_molefrac)
!
      call Get_do_AerDust_Calc (Chemistry, do_AerDust_Calc)
      call Get_synoz_threshold (Chemistry, synoz_threshold)
      call Get_t_cloud_ice     (Chemistry, t_cloud_ice)
      call Get_num_nox         (Chemistry, num_nox)
      call Get_num_noy         (Chemistry, num_noy)
      call Get_nox_map         (Chemistry, nox_map)
      call Get_noy_map         (Chemistry, noy_map)
      call Get_do_synoz        (Chemistry, do_synoz)
      call Get_do_nodoz        (Chemistry, do_nodoz)
!
      call Get_io3_num         (Chemistry, io3_num)
      call Get_ino_num         (Chemistry, ino_num)
      call Get_ico_num         (Chemistry, ico_num)
      call Get_ih2o_num        (Chemistry, ih2o_num)
      call Get_ihcl_num        (Chemistry, ihcl_num)
      call Get_ihno3_num       (Chemistry, ihno3_num)
      call Get_imgas_num       (Chemistry, imgas_num)
      call Get_ih2o2_num       (Chemistry, ih2o2_num)
      call Get_isynoz_num      (Chemistry, isynoz_num)
      call Get_ibrono2_num     (Chemistry, ibrono2_num)
      call Get_ioxygen_num     (Chemistry, ioxygen_num)
      call Get_ipropene_num    (Chemistry, ipropene_num)
      call Get_initrogen_num   (Chemistry, initrogen_num)
      call Get_iisoprene_num   (Chemistry, iisoprene_num)
!
      call Get_loss_opt        (Chemistry, loss_opt)
      call Get_sad_opt         (Chemistry, sad_opt)
      call Get_phot_opt        (Chemistry, phot_opt)
      call Get_do_chem_grp     (Chemistry, do_chem_grp)
      call Get_do_clear_sky    (Chemistry, do_clear_sky)
      call Get_oz_eq_synoz_opt (Chemistry, oz_eq_synoz_opt)
!
      call Get_dehydmin        (Chemistry, dehydmin)
!
      allocate (mw(numSpecies))
      call Get_mw (Chemistry, mw)
!
      call Get_do_grav_set (Advection, do_grav_set)
      call Get_advec_opt   (Advection, advec_opt)
!
      call Get_diffu_opt (Diffusion, diffu_opt)
!
      call Get_convec_opt (Convection, convec_opt)
!
      ! Get deposition related variables
!
      call Get_do_drydep   (Deposition, do_drydep)
      call Get_do_wetdep   (Deposition, do_wetdep)
      call Get_do_simpledep(Deposition, do_simpledep)
!
      ! Get diagnostics related variables
!
      call Get_do_qqjk_reset    (Diagnostics, do_qqjk_reset)
      call Get_num_tend_outrecs (Diagnostics, num_tend_outrecs)
      call Get_num_const_outrecs(Diagnostics, num_const_outrecs)
!
      ! Get clock related variables
!
      call Get_curGmiDate  (gmiClock, nymd          )
      call Get_curGmiTime  (gmiClock, nhms          )
      call Get_gmiSeconds  (gmiClock, gmi_sec       )
      call Get_gmiTimeStep (gmiClock, tdt           )
      call Get_numTimeSteps(gmiClock, num_time_steps)
!
      allocate(cosSolarZenithAngle(i1:i2,ju1:j2))
!
      allocate(diffaer   (i1:i2, ju1:j2, 1:numSpecies))
      allocate(s_radius  (i1:i2, ju1:j2, 1:numSpecies))
      allocate(s_velocity(i1:i2, ju1:j2, 1:numSpecies))
!
      allocate(pchem_water(i1:i2, ju1:j2, k1:k2))
      allocate(strat_water(i1:i2, ju1:j2, k1:k2))
      pchem_water(:,:,:) = 0.0d0
      strat_water(:,:,:) = 0.0d0
!
      Allocate (press3c(ilo:ihi, julo:jhi, k1:k2))
      Allocate (press3e(ilo:ihi, julo:jhi, k1-1:k2))
      press3c = 0.0d0
      press3e = 0.0d0
!
      Allocate (tropopausePress(i1:i2, ju1:j2))
      tropopausePress = 0.0d0
!
      Allocate (potentialVorticity(i1:i2, ju1:j2, k1:k2))
      potentialVorticity = 0.0d0
!
      Allocate (potentialTemp(i1:i2, ju1:j2, k1:k2))
      potentialTemp = 0.0d0
!
      Allocate (mass(i1:i2, ju1:j2, k1:k2))
      mass = 0.0d0
!
      allocate (pctm1Glob(ilo_gl:ihi_gl, julo_gl:jhi_gl)) 
      Allocate (kelGlob  (ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2))
      allocate (pctm1(ilo:ihi, julo:jhi)) 
      allocate (pctm2(ilo:ihi, julo:jhi))
      allocate (kel    (ilo:ihi, julo:jhi, k1:k2))
      allocate (uux    (ilo:ihi, julo:jhi, k1:k2))
      allocate (vvx    (ilo:ihi, jvlo:jhi, k1:k2))
!
      call Get_pctm1Glob(metFields, pctm1Glob)
      call Get_kelGlob(metFields, kelGlob)
      call Get_pctm1(metFields, pctm1)
      call Get_pctm2(metFields, pctm2)
      call Get_kel  (metFields, kel  )
      call Get_uux  (metFields, uux  )
      call Get_vvx  (metFields, vvx  )
      total_mass = 0.e0
      layer1_mass = 0.e0
!
      if (met_opt == 3) then
         allocate (humidity(i1:i2, ju1:j2,k1:k2))
         call Get_humidity(metFields, humidity)
         Allocate (relativeHumidity(i1:i2, ju1:j2, k1:k2))
         relativeHumidity = 0.0d0
!
         allocate (gridBoxHeight(i1:i2, ju1:j2, k1:k2))
         gridBoxHeight = 0.0d0
!
         allocate(fracCloudCover(i1:i2, ju1:j2))
         allocate(max_cloud(i1:i2, ju1:j2, k1:k2))
         allocate(ran_cloud(i1:i2, ju1:j2, k1:k2))
         allocate(tau_cloud(i1:i2, ju1:j2, k1:k2))
      end if
!
      !==================================================
      ! Compute 3D pressure at the center and at the edge
      !==================================================
!
      call CalcPress3dCenter (pt, am, bm, pctm1, press3c, pr_diag, procID, &
     &         i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
      call Set_press3c(metFields, press3c)
!
      call CalcPress3dEdge (pt, ai, bi, pctm1, press3e, pr_diag, procID, &
     &         i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
      call Set_press3e(metFields, press3e)
!
      if (met_opt == 3) then
         call finishMet  (metFields, new_met_rec, t_cloud_ice)
!
         call Get_humidity(metFields, humidity)
!
         call CalcGridBoxHeight (press3e, pctm1, humidity, kel, gridBoxHeight, &
     &            i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
         call Set_gridBoxHeight(metFields, gridBoxHeight)
!
      end if
!
!==========================================
!  Calculations of the tropopause pressure.
!  Two formulations are available.
!==========================================
!
!  Jim Stobie's formulation
!
!        call calcTropopausePress_Stobie  &
!     &    (press3c, kel, tropopausePress, pr_diag, procID, &
!     &       i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
!
!... Ertel's potential vorticity and potential temperature formulation
      call calcTropopausePress_Ertel  &
     &         (press3c, uux, vvx, kel, tropopausePress, potentialVorticity, &
     &          potentialTemp, dlatr, coscen, cosp, londeg, latdeg, &
     &          pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, &
     &          i1_gl, i2_gl, ju1_gl, j2_gl )
!
!... hyl wants WMO tropopause calc for tracers Rn, Pb, Be
      if(chem_opt.eq.9) then
         call calcTropopausePress_WMO        &
     &     (metdata_name_org, metdata_name_model, &
     &      kel, press3c, gridBoxHeight, latdeg, tropopausePress, &
     &      pr_diag, procID, &
     &      i1, i2, ju1, j2, k1, k2, &
     &      ilo, ihi, julo, jhi, ju1_gl, j2_gl)
      endif
!
!
      call Set_potentialTemp(metFields, potentialTemp)
      call Set_tropopausePress(metFields, tropopausePress)
      call Set_potentialVorticity(metFields, potentialVorticity)
!
      !=======================
      ! Compute the total mass
      !=======================
      call calcTotalMass (pt, ai, bi, pctm1, mcor, mass, pr_diag, procID,      &
     &         i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
      call Set_mass(metFields, mass)
!
      !===================================
      ! Reset fixed species concentrations
      !===================================
!
      call resetFixedConcentration (SpeciesConcentration, gmiClock, gmiGrid,   &
     &          numSpecies)
!
      !=====================================================
      ! Compute relative humidity and fractional cloud cover
      !=====================================================
!
      if (met_opt == 3) then
         call CalcRelativeHumidity (kel, press3c, humidity, relativeHumidity,  &
     &            pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
         call Set_relativeHumidity(metFields, relativeHumidity)
!
         call Get_max_cloud(metFields, max_cloud)
         call Get_ran_cloud(metFields, ran_cloud)
!
         call CalcFractionalCloudCover (max_cloud, ran_cloud, fracCloudCover,  &
     &            i1, i2, ju1, j2, k1, k2)
!
         call Set_fracCloudCover(metFields, fracCloudCover)
      end if
!
      !==============================================
      ! Compute the cosines of the solar zenith angle
      !==============================================
!
      call GetSecondsFromJanuary1(nsec_jan1, nymd, nhms)
      time = nsec_jan1
      days = time / SECPDY
!
      call CalcCosSolarZenithAngle(days, latdeg, londeg, &
               cosSolarZenithAngle, i1, i2, ju1, j2, &
               i1_gl, i2_gl, ju1_gl, j2_gl)
!
      !###########################
      ! Run the Emission Component
      !###########################
!
      if ((emiss_opt /= 0) .or. do_drydep .or. do_gcr) then
         if (do_ftiming) call Ftiming_On  ('gmiEmission')
!
         call doDiagnosticsBefore(Diagnostics, SpeciesConcentration, mass,  &
     &          EMISS_OP, i1, i2, ju1, j2, k1, k2, numSpecies)
!
         !---------------------------------------------------
         ! Update the surface temperature for MEGAN emissions
         !---------------------------------------------------
!
         if (doMEGANemission) then
            call updateSurfaceTempMEGAN (metFields, gmiGrid, gmiDomain, &
     &                 gmiClock, Diagnostics)
         end if
!
         !------------------------
         ! Run the Emission Driver
         !------------------------
!
         call runReadEmission (Emission, gmiClock, gmiGrid, gmiDomain,         &
     &           Diagnostics, do_drydep)
!
         call runEmission (Emission, SpeciesConcentration, gmiClock, gmiGrid,  &
     &           gmiDomain, Diagnostics, metFields, cosSolarZenithAngle, mw,   &
     &           chem_opt, do_drydep, convec_opt, &
     &           IBOC, IBBC, INOC, IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, ISSLT4, &
     &           IFSO2, INSO2, INDMS, IDUST1, IDUST2, IDUST3, IDUST4, IDUST5, IAN, IMGAS, INO, &
     &           iisoprene_num, ino_num, ico_num, ipropene_num, ihno3_num, io3_num)
!
         call doDiagnosticsAfter(SpeciesConcentration, mass, EMISS_OP, &
     &               i1, i2, ju1, j2, k1, k2, numSpecies)
!
         if (do_ftiming) call Ftiming_Off ('gmiEmission')
      end if
!
!
      if(tracer_opt /= 0) then
!
!... calculate total mass of atmosphere (kg)
         total_mass  = 100.d0 * sum((ai(k1-1)*pt+bi(k1-1)*pctm1Glob(i1_gl:i2_gl,ju1_gl:j2_gl))  &
            *mcorGlob(i1_gl:i2_gl, ju1_gl:j2_gl)) / GMI_G
!... calculate total mass of lowest layer of model (kg)
         layer1_mass = 100.d0 * sum((ai(k1  )*pt+bi(k1  )*pctm1Glob(i1_gl:i2_gl,ju1_gl:j2_gl))  &
            *mcorGlob(i1_gl:i2_gl, ju1_gl:j2_gl)) / GMI_G
         layer1_mass = total_mass - layer1_mass
!
         call calcTaggedCO_AgeOfAir (SpeciesConcentration, total_mass, layer1_mass,  &
                gmiClock, gmiDomain, Diagnostics, metFields, mw, londeg, latdeg, coscen,  &
                numSpecies, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi,  &
                i1_gl, i2_gl, ju1_gl, j2_gl)
      endif
!
      !===============================================================
      ! Subtract dehydmin from dehyd; dehyd will then have no negative
      ! values in it; dehydmin will be added back in later.
      !===============================================================
!
      if ((sad_opt == 2) .and. (dehydmin < 0.0d0)) then
         concentration(idehyd_num)%pArray3D(:,:,:) = &
     &               concentration(idehyd_num)%pArray3D(:,:,:) - dehydmin
      end if
!
!     call calcMaxMinArray ('DEHYD', &
!              concentration(idehyd_num)%pArray3D(:,:,:), pmin, pmax, &
!              (i2-i1+1)*(j2-ju1+1), k2, 1.0d0, procID, communicatorWorld)
!
      !############################
      ! Run the Diffusion Component
      !############################
!
      if (diffu_opt /= 0) then
         if (do_ftiming) call Ftiming_On  ('gmiDiffusion')
!
         call doDiagnosticsBefore(Diagnostics, SpeciesConcentration, mass,  &
     &          DIFFU_OP, i1, i2, ju1, j2, k1, k2, numSpecies)
!
!
         call runDiffusion (Diffusion, SpeciesConcentration, gmiClock, gmiGrid, &
     &            gmiDomain, Diagnostics, metFields)
!
         call doDiagnosticsAfter(SpeciesConcentration, mass, DIFFU_OP, &
     &               i1, i2, ju1, j2, k1, k2, numSpecies)
!
         if (do_ftiming) call Ftiming_Off ('gmiDiffusion')
      end if
!
      if (do_ftiming) then
         call Ftiming_On  ('procSyncBeforeAdvection')
         call synchronizeGroup(communicatorWorld)
         call Ftiming_Off ('procSyncBeforeAdvection')
      end if
!
      !############################
      ! Run the Advection Component
      !############################
!
      if (advec_opt /= 0) then
         if (do_ftiming) call Ftiming_On  ('gmiAdvection')
!
         call doDiagnosticsBefore(Diagnostics, SpeciesConcentration, mass,  &
     &          ADVEC_OP, i1, i2, ju1, j2, k1, k2, numSpecies)
!
         if (advec_opt == 1) then
            !-------------------------------
            ! For the old TpCore formulation
            !-------------------------------
!
            call runAdvection (Advection, gmiGrid, gmiDomain, Diagnostics, &
     &              metFields, do_chem_grp, concentration, mw)
!
         elseif (advec_opt == 2) then
            !------------------------
            ! For advCore formulation
            !------------------------
!
            allocate (xmass(ilo:ihi, julo:jhi, k1:k2))
            allocate (ymass(ilo:ihi, julo:jhi, k1:k2))
            allocate (zmass(i1:i2, ju1:j2, k1:k2))
            xmass = 0.0d0
            ymass = 0.0d0 
            zmass = 0.0d0 
!
            call Get_xmass(metFields, xmass)
            call Get_ymass(metFields, ymass)
            call Get_zmass(metFields, zmass)
!
            allocate( MX_UR(i1:i2, ju1:j2, k1:k2) )
            MX_UR = 0.0
!
            allocate( MY_UR(i1:i2, ju1:j2, k1:k2) )
            MY_UR = 0.0
!
            allocate( MX(i1:i2, ju1:j2, k1:k2) )
            MX = 0.0
!
            allocate( MY(i1:i2, ju1:j2, k1:k2) )
            MY = 0.0
!
            allocate( MZ(i1:i2, ju1:j2, k1:k2) )
            MZ = 0.0
!
            allocate( UC(i1:i2, ju1:j2, k1:k2) )
            UC = 0.0
            allocate( VC(i1:i2, ju1:j2, k1:k2) )
            VC = 0.0
!
            do ic = i1, i2
               UC(ic,:,:) = 0.5*(uux(ic-1,:,k2:k1:-1) + uux(ic,:,k2:k1:-1))
            end do
!
            do ic = ju1, j2
               VC(:,ic,:) = 0.5*( vvx(:,ic-1,k2:k1:-1) + vvx(:,ic,k2:k1:-1))
            end do
!
            allocate(delpm(ilo:ihi, julo:jhi, k1:k2))
            call Calc_Delpm(dap, dbk, pctm1, pctm1, delpm, pr_diag, procID, &
     &                k1, k2, ilo, ihi, julo, jhi)
!
            allocate( DP(i1:i2, ju1:j2, k1:k2) )
            DP(i1:i2, ju1:j2, k1:k2) = mbToPascal*delpm(i1:i2, ju1:j2, k2:k1:-1)
!
            allocate( DPEDT(i1:i2, ju1:j2, k1:k2) )
            DPEDT = 0.0
!
            call convertMassFlux (xmass, ymass, zmass, MX, MY, MZ, tdt, &
     &                            i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
            MX_UR = MX
            MY_UR = MY
!
            call ESMF_GridCompGet(advCoreESMF%compGridded, grid = esmfGrid, &
     &                            rc=STATUS)
            VERIFY_(STATUS)
!
            call fromGmiToAdvecCore(advCoreESMF%stateImp, esmfGrid,            &
     &               SpeciesConcentration, Advection, UC, VC, MX_UR, MY_UR,    &
     &               MX, MY, MZ, DPEDT, DP, i1, i2, ju1, j2, k1, k2)
!
            call ESMF_GridCompRun(advCoreESMF%compGridded,    &
     &                    importState = advCoreESMF%stateImp, &
     &                    exportState = advCoreESMF%stateExp, &
     &                    clock       = esmfClock, rc=STATUS)
            VERIFY_(STATUS)
!
            call fromAdvecCoreToGmi (advCoreESMF%stateImp, esmfGrid, &
     &               SpeciesConcentration, Advection)
!
            deallocate(DP)
            deallocate(UC)
            deallocate(VC)
            deallocate(MX)
            deallocate(MY)
            deallocate(MZ)
            deallocate(DPEDT)
            deallocate(MX_UR)
            deallocate(MY_UR)
            deallocate(xmass)
            deallocate(ymass)
            deallocate(zmass)
!
            call synchronizeGroup(communicatorWorld)
         end if
!
         !------------------------------------------------------------
         ! Prior to advection, pctm1 is used for the surface pressure;
         ! after advection, pctm2 is used.
         !------------------------------------------------------------
!
         call CalcPress3dCenter (pt, am, bm, pctm2, press3c, pr_diag, procID, &
     &            i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
         call Set_press3c(metFields, press3c)
!
         call CalcPress3dEdge (pt, ai, bi, pctm2, press3e, pr_diag, procID, &
     &            i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
         call Set_press3e(metFields, press3e)
!
!==========================================
!  Calculations of the tropopause pressure.
!  Two formulations are available.
!==========================================
!
!  Jim Stobie's formulation
!
!        call calcTropopausePress_Stobie  &
!     &    (press3c, kel, tropopausePress, pr_diag, procID, &
!     &       i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
!
!  Ertel's potential vorticity and potential temperature formulation
!
         call calcTropopausePress_Ertel (press3c, uux, vvx, kel,               &
     &            tropopausePress, potentialVorticity, potentialTemp, dlatr,   &
     &            coscen, cosp, londeg, latdeg, pr_diag, procID, i1, i2, ju1,  &
     &            j2, k1, k2, ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl )
!
         call Set_potentialTemp(metFields, potentialTemp)
!
!hyl,12/13/11
         if(chem_opt.eq.9) then
            call calcTropopausePress_WMO        &
     &        (metdata_name_org, metdata_name_model, &
     &         kel, press3c, gridBoxHeight, latdeg, tropopausePress, &
     &         pr_diag, procID, &
     &         i1, i2, ju1, j2, k1, k2, &
     &         ilo, ihi, julo, jhi, ju1_gl, j2_gl)
         endif
!
         call Set_tropopausePress(metFields, tropopausePress)
         call Set_potentialVorticity(metFields, potentialVorticity)
!
         !---------------------
         ! Calculate Total Mass
         !---------------------
         call calcTotalMass (pt, ai, bi, pctm2, mcor, mass, pr_diag, procID, &
     &            i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
         call Set_mass(metFields, mass)
!
         !--------------------------
         ! Compute Relative Humidity
         !--------------------------
!
         if (met_opt == 3) then
            call CalcRelativeHumidity (kel, press3c, humidity,                 &
     &               relativeHumidity, pr_diag, procID, i1, i2, ju1, j2, k1,   &
     &               k2, ilo, ihi, julo, jhi)
!
            call Set_relativeHumidity(metFields, relativeHumidity)
         end if
!
         call doDiagnosticsAfter(SpeciesConcentration, mass, ADVEC_OP, &
     &               i1, i2, ju1, j2, k1, k2, numSpecies)
!
         if (do_ftiming) call Ftiming_Off ('gmiAdvection')
      end if
!
      !========================
      ! Advance the model clock
      !========================
!
      call advanceESMFclock(esmfClock, gmiClock)
!
      call Get_curGmiDate  (gmiClock, nymd          )
      call Get_curGmiTime  (gmiClock, nhms          )
      call Get_numTimeSteps(gmiClock, num_time_steps)
      call Get_gmiSeconds  (gmiClock, gmi_sec       )
!
      !==========================================
      ! Update the aerosol gravitational settling
      !==========================================
!
      if (do_grav_set .or. do_drydep) then
         call doDiagnosticsBefore(Diagnostics, SpeciesConcentration, mass,  &
     &          SETTLING_OP, i1, i2, ju1, j2, k1, k2, numSpecies)
!
         call updateGravitationalSettling (SpeciesConcentration, gmiClock,     &
     &              gridBoxHeight, humidity, mass, press3e, kel, diffaer,      &
     &              s_radius, s_velocity, pr_diag, procID, chem_opt, i1, i2,   &
     &              ju1, j2, k1, k2, ilo, ihi, julo, jhi, numSpecies)
!
         call doDiagnosticsAfter(SpeciesConcentration, mass, SETTLING_OP, &
     &          i1, i2, ju1, j2, k1, k2, numSpecies)
!
!      
      end if
!
      !#############################
      ! Run the Convection Component
      !#############################
!
      if (convec_opt /= 0) then
         if (do_ftiming) call Ftiming_On  ('gmiConvection')
!
         call doDiagnosticsBefore(Diagnostics, SpeciesConcentration, mass,  &
     &          CONVEC_OP, i1, i2, ju1, j2, k1, k2, numSpecies)
!
!
         if (convec_opt == 1) then  ! DAO convection
            allocate (bmass(i1:i2, ju1:j2, k1:k2))
            call CalcConvection_bmass (bmass, pctm2, ai, bi, pt, &
     &               i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
         end if
!
#ifdef MICRO_AEROSOL
       allocate(REL_SCAV_EFF_new(i1:i2, ju1:j2, k1:k2, 1:numSpecies))
#endif
!
         call runConvection (Convection, Deposition, SpeciesConcentration,     &
     &           gmiClock, gmiGrid, gmiDomain, Diagnostics, metFields,         &
     &           chem_opt, ih2o2_num, ihno3_num,  mw, bmass, REL_SCAV_EFF_new, &
     &           chem_mecha, i1, ju1, k1)
!
         if (convec_opt == 1) then  ! DAO convection
             deallocate (bmass)
         endif
!
         call doDiagnosticsAfter(SpeciesConcentration, mass, CONVEC_OP, &
     &          i1, i2, ju1, j2, k1, k2, numSpecies)
!
         if (do_ftiming) call Ftiming_Off ('gmiConvection')
!
      end if
!
      !=============================================
      ! Compute the cosines of the solar zeith angle
      !=============================================
!
       call GetSecondsFromJanuary1(nsec_jan1, nymd, nhms)
       time = nsec_jan1
       days = time / SECPDY
       call CalcCosSolarZenithAngle(days, latDeg, lonDeg, cosSolarZenithAngle, &
     &          i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl)
!
      !#############################
      ! Run the Deposition Component
      !#############################
!
      if (do_drydep .or. do_wetdep) then
!
         if (do_drydep) then
            !----------------------
            ! Dry Deposition Driver
            !----------------------
!
            if (do_ftiming) call Ftiming_On  ('gmiDryDeposition')
!
            call doDiagnosticsBefore(Diagnostics, SpeciesConcentration, mass,  &
     &             DRYDEP_OP, i1, i2, ju1, j2, k1, k2, numSpecies)
!
!
            allocate (BoxHeightCenter(i1:i2,ju1:j2))
            allocate (BoxHeightEdge  (i1:i2,ju1:j2))
!
            call CalcDryDepBoxHeight (BoxHeightCenter, BoxHeightEdge, pctm2,   &
     &               kel, humidity, press3c, press3e, i1, i2, ju1, j2, k1, k2, &
     &               ilo, ihi, julo, jhi)
!
            call runDryDeposition (Deposition, Emission, SpeciesConcentration, &
     &              gmiClock, gmiGrid, gmiDomain, Diagnostics, metFields,      &
     &              cosSolarZenithAngle, diffaer, s_radius, s_velocity,        &
     &              BoxHeightCenter, BoxHeightEdge, chem_opt, mw, i1, i2, ju1, &
     &              j2, numSpecies)
!
            deallocate (BoxHeightCenter)
            deallocate (BoxHeightEdge  )
!
            call doDiagnosticsAfter(SpeciesConcentration, mass, DRYDEP_OP, &
     &             i1, i2, ju1, j2, k1, k2, numSpecies)
!
            if (do_ftiming) call Ftiming_Off ('gmiDryDeposition')
         end if
!
         if (do_wetdep) then
            !----------------------
            ! Wet Deposition Driver
            !----------------------
!
            if (do_ftiming) call Ftiming_On  ('gmiWetDeposition')
!
            call doDiagnosticsBefore(Diagnostics, SpeciesConcentration, mass,  &
     &             WETDEP_OP, i1, i2, ju1, j2, k1, k2, numSpecies)
!
!
            ! Compute the moisture changes
            if ((metdata_name_org(1:4) == 'GISS') .or.  &
     &         ((metdata_name_org  (1:4) == 'NCAR' ) .and.  &
     &          (metdata_name_model(1:5) == 'MATCH'))) then
               allocate(moistq(i1:i2,ju1:j2,k1:k2))
               allocate(rain(i1:i2,ju1:j2,k1:k2))
               call Get_rain(metFields, rain)
!
               call CalcMoistureChanges (moistq, rain, mcor, mass, &
     &                                  i1, i2, ju1, j2, k1, k2, ivert)
!
               call Set_moistq(metFields, moistq)
               deallocate(rain)
               deallocate(moistq)
            end if
!
            call RunWetDeposition (Deposition, SpeciesConcentration, gmiClock, &
     &              gmiGrid, gmiDomain, Diagnostics, metFields, chem_opt,      &
     &              ih2o2_num, ihno3_num, mw, REL_SCAV_EFF_new, i1, ju1, k1)
!
            call doDiagnosticsAfter(SpeciesConcentration, mass, WETDEP_OP, &
     &               i1, i2, ju1, j2, k1, k2, numSpecies)
!
            if (do_ftiming) call Ftiming_Off ('gmiWetDeposition')
!
         end if
!
      else if (do_simpledep) then
         !-------------------------
         ! Simple Deposition Driver
         !-------------------------
!
         if (do_ftiming) call Ftiming_On  ('gmiSimpleDeposition')
!
         call doDiagnosticsBefore(Diagnostics, SpeciesConcentration, mass,  &
     &             SIMPDEP_OP, i1, i2, ju1, j2, k1, k2, numSpecies)
!
         call runSimpleDeposition (Deposition, SpeciesConcentration, gmiClock, &
     &           gmiGrid, gmiDomain, Diagnostics, metFields, ibrono2_num,      &
     &           ih2o2_num, ihcl_num, ihno3_num, imgas_num, initrogen_num,     &
     &           ioxygen_num)
!
         call doDiagnosticsAfter(SpeciesConcentration, mass, SIMPDEP_OP, &
     &          i1, i2, ju1, j2, k1, k2, numSpecies)
!
         if (do_ftiming) call Ftiming_Off ('gmiSimpleDeposition')
!
      end if
!
      !============================================================
      ! Add dehydmin back into dehyd; dehyd will then have negative
      !  values in it again.
      !============================================================
!
      if ((sad_opt == 2) .and. (dehydmin < 0.0d0)) then
         concentration(idehyd_num)%pArray3D(:,:,:) = &
     &              concentration(idehyd_num)%pArray3D(:,:,:) + dehydmin
      end if
!... make sure strat_o3 is > 0.0
      if(tracer_opt == 8) then
         where (concentration(2)%pArray3d(:,:,:) < 1.0d-30)
                concentration(2)%pArray3d(:,:,:) = 1.0d-30
         end where
      endif
!
!
      !############################
      ! Run the Chemistry Component
      !############################
!
      if ((chem_opt /= 0) .or. (phot_opt /= 0) .or. (sad_opt  /= 0)) then
         if (do_ftiming) call Ftiming_On  ('gmiChemistry')
!
         if ((met_opt == 3) .and. (ih2o_num /= 0)) then
            call doDiagnosticsBefore(Diagnostics, SpeciesConcentration, mass,  &
     &             ADDWAT_OP, i1, i2, ju1, j2, k1, k2, numSpecies)
!
            !--------------------------------
            ! Tropospheric water calculations
            !--------------------------------
!
            call addTroposphericWater (chem_mecha, ih2o_num, press3c,          &
     &              concentration, tropopausePress, humidity, pchem_water,     &
     &              strat_water, pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo,&
     &              ihi, julo, jhi, numSpecies)
!
            call doDiagnosticsAfter(SpeciesConcentration, mass, ADDWAT_OP, &
     &             i1, i2, ju1, j2, k1, k2, numSpecies)
         end if
!
         do ic = 1, num_chem
            where (concentration(ic)%pArray3d(:,:,:) < 1.0d-30)
                concentration(ic)%pArray3d(:,:,:) = 1.0d-30
            end where
         end do
!
         call doDiagnosticsBefore(Diagnostics, SpeciesConcentration, mass,  &
     &          CHEM_OP, i1, i2, ju1, j2, k1, k2, numSpecies)
!
         !-----------------------------------------
         ! Calculations for cloud related variables
         !-----------------------------------------
!         
         if ((chem_opt == 2) .or. (chem_opt == 7) .or. (chem_opt == 8)) then
!.sds         if ((chem_opt == 7) .or. (chem_opt == 8)) then
            allocate(totalCloudFraction(i1:i2,ju1:j2,k1:k2))
!
            call CalcTotalCloudFraction (max_cloud, ran_cloud,      &
     &               totalCloudFraction, i1, i2, ju1, j2, k1, k2)
!
            call Set_totalCloudFraction(metFields, totalCloudFraction)
            deallocate(totalCloudFraction)
         end if
!
         if ( .not. do_clear_sky) then
            if (metdata_name_org(1:3) == 'DAO') then
                call CalcCloudOpticalDepth (kel, press3e, pctm2, max_cloud,    &
     &                   ran_cloud, tau_cloud, i1, i2, ju1, j2, k1, k2, ilo,   &
     &                   ihi, julo, jhi)
!
               call Set_tau_cloud(metFields, tau_cloud)
            end if
         else
             tau_cloud = 0.0d0
             call Set_tau_cloud(metFields, tau_cloud)
         end if
!
         if ((chem_mecha == 'troposphere').or. (chem_mecha == 'strat_trop')  &
     &         .or. chem_mecha == 'strat_trop_aerosol' ) then
            if ((phot_opt == 3) .and. do_AerDust_Calc) then
               call Get_tau_cloud(metFields, tau_cloud)
!
               Allocate (optDepth(i1:i2, ju1:j2, k1:k2, num_AerDust))
               call Get_optDepth(Chemistry, optDepth)
!
               call DiagCloudOpticalDepth (max_cloud, ran_cloud, tau_cloud, &
                                 optDepth, i1, i2, ju1, j2, k1, k2, num_AerDust)
!
               call Set_optDepth(Chemistry, optDepth)
               deallocate(OptDepth)
            end if
!
         end if
!
         !---------------------------------------------------
         ! Read in chemistry related files (monthly or daily)
         !---------------------------------------------------
!
!.sds  add in using calculated dust/BC/OC/SeaSalt/SO4 for AOT calculation for photolysis calc
         if ( phot_opt == 3 .and. do_AerDust_Calc ) then
            if ( (chem_mecha == 'strat_trop_aerosol') .or.  &
     &           (chem_mecha == 'strat_trop' .and. INSO4N1 /= 0) ) then
!
               call runCalcAerosolDust (Chemistry, &
     &              mass, mcor, gridBoxHeight, SpeciesConcentration,  &
     &              gmiGrid, gmiDomain, Diagnostics, chem_mecha)
!
!...  reads in aerosol/dust burdens for AOT calc for phot calc
            else
               call runReadChemistry(Chemistry, Emission, gmiClock, gmiGrid, &
     &           gmiDomain, Diagnostics, chem_mecha)
            endif
         endif
!
         !--------------------------
         ! Call the Chemistry driver
         !--------------------------
!
         call runChemistry (Chemistry, Emission, SpeciesConcentration,  &
     &           gmiClock, gmiGrid, gmiDomain, Diagnostics, metFields,  &
     &           chem_mecha, averagePressEdge, num_AerDust)
!
#ifndef MICRO_AEROSOL
!micro-----------begin----------------------------------------------------
!specific call for troposphere
! Option for Stratospheric Loss
!
         if ((loss_opt == 1) .and. (chem_mecha == 'troposphere')) then
            data (lexclude(ibr), ibr=1,4) /IHNO3, INO, INO2, IO3/
!
            data (llife(ibr), ibr=1,17) /  &
     &            IALD2, IALK4, IC2H6, IC3H8, ICH2O, ICO, IH2O2,  &
     &            IHNO4, IMACR,  IMCO3, IMP, IMVK, IPMN, IC3H6, IR4N2,  &
     &            IRCHO,IC3H6O/
!
            loss_level = 150.e-9
            strat_life100 = 100*86400.
            strat_life1 = 1*86400.
!
            do ibr = 1,NACT
!  exclude 100 day lifetime species
               if (.not. any(lexclude == ibr)) then
                  if (any(llife == ibr)) then
                      strat_life=strat_life100
                  else
                      strat_life=strat_life1
                  endif
                  where (concentration(isynoz_num)%pArray3D(:,:,:) > loss_level)
                      concentration(ibr)%pArray3D(:,:,:) = &
     &                    concentration(ibr)%pArray3D(:,:,:) *exp(-tdt/strat_life)
                  end where
               endif
            enddo
         endif
! Stratospheric Loss Option Complete
!micro-----------end------------------------------------------------------
#endif
!
         call doDiagnosticsAfter(SpeciesConcentration, mass, CHEM_OP, &
     &          i1, i2, ju1, j2, k1, k2, numSpecies)
!
         do ic = 1, num_chem
            where (concentration(ic)%pArray3d(:,:,:) < 1.0d-30)
                concentration(ic)%pArray3d(:,:,:) = 1.0d-30
            end where
         end do
!
         if ((met_opt == 3) .and. (ih2o_num /= 0)) then
            call doDiagnosticsBefore(Diagnostics, SpeciesConcentration, mass,  &
     &             REMWAT_OP, i1, i2, ju1, j2, k1, k2, numSpecies)
!
            call removeTroposphericWater (chem_mecha, ih2o_num, concentration, &
     &                 pchem_water, strat_water, pr_diag, procID, i1, i2, ju1, &
     &                 j2, k1, k2, numSpecies)
!
            call doDiagnosticsAfter(SpeciesConcentration, mass, REMWAT_OP, &
     &             i1, i2, ju1, j2, k1, k2, numSpecies)
         end if
!
         if (do_ftiming) call Ftiming_Off ('gmiChemistry')
!
      end if
!
      !================================
      ! Synthetic scpecies calculations
      !================================
!
      if (do_synoz) then
         call doDiagnosticsBefore(Diagnostics, SpeciesConcentration, mass,  &
     &          SYNSPC_OP, i1, i2, ju1, j2, k1, k2, numSpecies)
!
         call updateSyntheticSpecies (do_nodoz, latdeg, press3e, concentration,&
     &              synoz_threshold, tdt, mw, io3_num, isynoz_num, ihno3_num,  &
     &              num_nox, num_noy, nox_map, noy_map, oz_eq_synoz_opt, i1,   &
     &              i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ju1_gl, j2_gl,   &
     &              numSpecies, chem_mecha)
!
         call doDiagnosticsAfter(SpeciesConcentration, mass, SYNSPC_OP, &
     &               i1, i2, ju1, j2, k1, k2, numSpecies)
!
      end if
!
      if ((chem_opt ==2) .and. (num_ks_sbc > 0) .and. (num_spc_sbc > 0)) then
         do ic = 1, num_spc_sbc
            do ik = 1, num_ks_sbc
               concentration(surf_bc_map(ic))%pArray3D(:,:,ik) = &
     &                                                 surf_bc(:,:,ik,ic)
            end do
         end do
      end if
!
      !=============================================================
      ! The following statements alleviate some numeric difficulties
      ! with dehyd.
      !=============================================================
!
      if ((advec_opt == 1) .and. (sad_opt == 2)) then
         where (Abs (concentration(idehyd_num)%pArray3D(:,:,:)) < 1.0d-30)
            concentration(idehyd_num)%pArray3D(:,:,:) = 1.0d-30
         end where
      end if
!
      !============================================
      ! Calculate the surface species concentration
      !============================================
!
      call calcConcentrationSurf(SpeciesConcentration, pr_diag, procID, numSpecies)
!
      !===========================================
      ! Calculate the column species concentration
      !===========================================
!
      if (tracer_opt == 0 .and. chem_opt >= 2) then
         if ((chem_mecha == 'troposphere').or. (chem_mecha == 'strat_trop')  &
     &         .or. chem_mecha == 'strat_trop_aerosol' ) then
            call calcConcentrationColTrop(SpeciesConcentration, mcor, mass, &
     &               synoz_threshold, isynoz_num, pr_diag, procID, &
     &               i1, i2, ju1, j2, k1, k2, numSpecies)
            if (chem_mecha == 'strat_trop'  &
     &         .or. chem_mecha == 'strat_trop_aerosol' ) then
               call calcConcentrationColCombo(SpeciesConcentration, mcor, mass, &
     &                  pr_diag, procID, i1, i2, ju1, j2, k1, k2, numSpecies)
            end if
         end if
      end if
!
      !=====================
      ! Deallocate variables
      !=====================
!
      deallocate(cosSolarZenithAngle)
      deallocate(diffaer)
      deallocate(s_radius)
      deallocate(s_velocity)
      deallocate(pchem_water)
      deallocate(strat_water)
!
      deallocate(mw)
      deallocate(ai)
      deallocate(bi)
      deallocate(am)
      deallocate(bm)
      deallocate(dap)
      deallocate(dbk)
      deallocate(dlatr )
      deallocate(coscen)
      deallocate(cosp  )
      deallocate(latdeg)
      deallocate(londeg)
      deallocate(mcor  )
      deallocate(mcorGlob)
!
      deallocate(kel)
      deallocate(uux)
      deallocate(vvx)
      deallocate(mass)
      deallocate(pctm1Glob)
      deallocate(kelGlob)
      deallocate(pctm1)
      deallocate(pctm2)
      deallocate(press3c)
      deallocate(press3e)
      deallocate(potentialTemp)
      deallocate(tropopausePress)
      deallocate(potentialVorticity)
!
      if (met_opt == 3) then
         deallocate(humidity)
         deallocate(relativeHumidity)
         deallocate(gridBoxHeight)
         deallocate(fracCloudCover)
         deallocate(max_cloud)
         deallocate(ran_cloud)
         deallocate(tau_cloud)
      end if
!
      if (do_ftiming) then
         call Ftiming_On  ('procSyncEndStepping')
         call synchronizeGroup(communicatorWorld)
         call Ftiming_Off ('procSyncEndStepping')
         call Ftiming_Off ('gmiTimeStepping')
      end if
!
      return
!
      end subroutine gmiTimeStepping
!EOC
!------------------------------------------------------------------------------
      end module GmiTimeStepping_mod
!
