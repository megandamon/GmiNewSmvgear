!------------------------------------------------------------------------------
!BOP
! !MODULE: GmiUpdateChemistry_mod
!
      module GmiUpdateChemistry_mod
!
! !USES:
      use GmiFileOperations_mod, only : makeOutfileName
      use GmiUpdateForcingBC_mod, only : updateForcingBC
      use GmiStratosphericLoss_mod, only : updateStratosphericLoss
      use Ftiming_Dao              , only : Ftiming_On, Ftiming_Off
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiPrintError_mod, only : GmiPrintError
      use GmiShipEmission_mod, only : calcShipEmission
      use GmiPhotolysisRateConstants_mod, only : calcPhotolysisRateConstants
      use GmiThermalRateConstants_mod, only : calcThermalRateConstants, Accum_Qqjk
      use GmiCloudPropertiesGT_mod   , only : calcCloudPropertiesGT
      use GmiUpdateSAD_mod           , only : updateSurfaceAreaDensities
      use GmiProdLossDiagnostics_mod, only : doDiagnosticsBefore_BC, doDiagnosticsAfter_BC
      use GmiSolverInterface_mod, only : Update_Smv2chem, Update_Sulfchem
      use GmiTracerMethod_mod, only : Update_Beryl, Update_CH3I, Update_Uniform,  &
           Update_Radon_Lead, Update_LINOZ, Update_SF6, Update_strato3, Update_Radon_Lead_Beryl
      use GmiSolar_mod,        ONLY : computeSolarZenithAngle_Photolysis
      use GmiTimeControl_mod      , only : GmiSplitDateTime
      use GmiTimeControl_mod      , only : GetDaysFromJanuary1
      use GmiTimeControl_mod      , only : ConvertTimeToSeconds
      use GmiTimeControl_mod      , only : GetSecondsFromJanuary1
      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved, t_SulfSaved
      use GmiSolver_SavedVariables_mod, only : t_CloudParametersGT
      use fastJX_Bundle_mod,   ONLY : t_fastJXbundle
!
      implicit none

      private
      public  :: updateChemistry

#     include "GmiParameters.h"
!EOP
!------------------------------------------------------------------------------
      contains

!=============================================================================
!
! $Id: GmiUpdateChemistry_mod.F90,v 1.41 2013-09-16 19:59:22 jkouatch Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   GmiUpdateChemistry_mod.F90
!
! ROUTINES
!   updateChemistry
!   Update_Quadchem
!
! HISTORY
!   - September 23, 2004 * Jules Kouatchou
!     Modified the routine Update_Chem to accomodate for the combined
!     stratosphere/troposphere chemical mechanism.
!     The added variables are: idumday, idumyear, month, h2ocombo, sadcombo.
!     The computations for "sad_opt = 2" are rewritten for the new
!     mechanism.
!
!     The routine Update_Smv2chem was also modified to remove the condition
!     on "do_synoz" for the combined mechanism.
!   - January 24, 2005 * Jules Kouatchou
!     Added the argument "gridBoxHeight" to the subroutine "Update_Qj" for use
!     in the aerosol/dust surface area and optical depth calculations.
!   - September 15, 2005 * Jules Kouatchou
!     Modified the routine Update_Chem to add "phot_opt" as argument of the
!     routine Update_Qk.
!   - December 8, 2005 * Bigyani Das
!     Added 3-D carbon emission through emiss field in Update_Sulfchem by
!     adding an array emiss_dust for DMS emission. This is done by
!     adding the argument "emiss_dust" to the subroutine call 
!     for Update_Sulfchem (DMS emission is in the emiss_dust file).
!   - March 13, 2006 * Jules Kouatchou
!     Added ship emission calculations.
!     They are only done if the user set the namelist variable 
!     "do_ShipEmission" to true.
!   - 2012 * Stephen Steenrod
!     Pulled out tracer routines to their own module and
!     introduced the "tracer" mechanism to run them as a package
!=============================================================================

!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: updateChemistry
!
! !INTERFACE:
!
       subroutine updateChemistry (smv2SavedVars, sulfSavedVars,      &
          JXbundle, cloudGT, &
          chemintv, rateintv, jno_num, jnoxnum, chem_mecha, &
          do_ftiming, metdata_name_org, metdata_name_model, lossData, &
      &   KDIM, JDIM, MDIM, pr_qj_o3_o1d, pr_qj_opt_depth,            &
          do_qqjk_inchem, do_qqjk_reset, pr_qqjk, emiss_isop,         &
          emiss_monot, emiss_nox, emiss_hno3, emiss_o3, emiss_ozone,  &
          prevRecord, curRecord,    &
          tropp, press3c, press3e, averagePressEdge, gridBoxHeight,   &
          taucli, tauclw, qj_labels, const_labels, pr_smv2,           &
          pr_sulf_src, do_aerocom, pr_surf_emiss, pr_emiss_3d,        &
          pr_nc_period, num_AerDust, aqua_infile_name, qj_timpyr,     &
          qj_infile_name, qj_var_name, const_init_val, dlatr, latdeg, &
          londeg, mcor, mass, loss_freq, prod, loss, concentration,   &
          overheadO3col, decay_3d_out, qjgmi, qkgmi, emissionArray, rxnr_adjust,    &
      &   kel, humidity, pctm2,  totalCloudFraction, fracCloudCover,  &
      &   tau_cloud, lwi_flags, moistq, clwc, cmf, surf_air_temp,     &
          surf_alb_uv, dms_oh, dms_no3, so2_oh, so2_h2o2, so2_o3,     &
          emiss_dust, qqjgmi, qqkgmi,  yda, qqkda, qqjda,             &
          surf_emiss_out, surf_emiss_out2, emiss_3d_out, ch4clim,     &
          h2oclim, hno3gas, lbssad, sadgmi, hno3cond, h2oback,        &
          h2ocond, reffice, reffsts, vfall, OptDepth, Eradius, Tarea, &
          Odaer, relativeHumidity, radswg, Odmdust, Dust, Waersl,     &
          Daersl, cloud_tau, cloud_param, flux_gt, cloudDroplet_opt,  &
          be_opt, chem_opt, sad_opt, phot_opt, fastj_opt,             &
          fastj_offset_sec, do_smv_reord, do_synoz, do_semiss_inchem, &
          do_wetchem, do_aqu_chem, do_clear_sky, do_AerDust_Calc,     &
          do_ozone_inFastJX, do_ShipEmission, o3_index, loss_freq_opt,&
          kmin_loss, kmax_loss, t_half_be7, t_half_be10, yield_be7,   &
          yield_be10, dehydmin, fbc_j1, fbc_j2, forc_bc_num,          &
          forc_bc_kmax, forc_bc_kmin, forc_bc_opt, forc_bc_map,       &
          forc_bc_incrpyr, forc_bc_start_num, forc_bc_years,          &
          forc_bc_data, lastYearFBC, jlatmd, nymd,    &
          nhms, start_ymd, gmi_sec, tdt, num_time_steps, mw, pr_diag, &
          loc_proc, synoz_threshold, chem_cycle, AerDust_Effect_opt,  &
      &   lbssad_timpyr, h2oclim_timpyr, chem_mask_klo, chem_mask_khi,&
      &   emiss_timpyr, emiss_opt, emiss_map, dehyd_opt, h2oclim_opt, &
          lbssad_opt, emiss_dust_opt, emiss_map_dust, ndust,          &
          rxnr_adjust_map, ih2_num, io3_num, ih2o_num, ih2ocond_num,  &
          ihno3_num, ihno3cond_num, idehyd_num, ih2oaircr_num,        &
          ich4_num, imgas_num, ico_num, ino_num, ipropene_num,        &
          iisoprene_num, initrogen_num, ioxygen_num, isynoz_num,      &
      &   num_emiss, num_species, num_qks, num_qjs, num_qjo, num_sad, &
          num_molefrac, num_chem, num_active, num_rxnr_adjust,        &
          rxnr_adjust_timpyr, ilong, ilat, ivert, itloop, i1_gl,      &
          i2_gl, ju1_gl, j2_gl, k1_gl, k2_gl, ilo, ihi, julo, jhi, i1,&
          i2, ju1, j2, k1, k2, numDomains, commuWorld, ih2o2_num,     &
          jno_adjust,rootProc, linoz_infile_name, sf6_infile_name,    &
          SO3daily_infile_name, SO3monthly_infile_name)
!
! !USES:
!
      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_sad_constants.h"
#     include "setkin_par.h"
#     include "setkin_depos.h"
#     include "gmi_AerDust_const.h"
#     include "gmi_diag_constants_llnl.h"

!
! !INPUT PARAMETERS:
!!   cross_section_file    : X-Section quantum yield
!!   rate_file             : Master rate file
!!   T_O3_climatology_file : T & O3 climatology
!!   metdata_name_org      : first  part of metdata_name, e.g., "NCAR"
!!   metdata_name_model    : second part of metdata_name, e.g., "MATCH"
!!   pr_qj_o3_o1d          : should special reaction o3->o1d be saved?
!!   pr_qj_opt_depth       : should optical depth be saved at end of
!!                           qj file?
!!   do_qqjk_inchem        : if pr_qqjk is on, should qqj's & qqk's be
!!                           determined inside chemistry solver?
!!   do_qqjk_reset         : reset qqjk accumulators (i.e., qqjk output
!!                           occurred on previous time step)?
!!   pr_qqjk               : should the periodic qqjk output file be
!!                           written?
!!   emiss_isop            : isoprene    emissions (kg/s)
!!   emiss_monot           : monoterpene emissions (kg/s)
!!   emiss_nox             : NOx         emissions (kg/s)
!!   tropp                 : tropopause pressure   (mb)
!!   press3c               : atmospheric pressure at the center of
!!                                                   each grid box (mb)
!!   press3e               : atmospheric pressure at the edge   of
!!                                                   each grid box   (mb)
!!   gridBoxHeight           : height of each grid box (m)
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: ilong, ilat, ivert, itloop
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, k1_gl, k2_gl
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: num_emiss, num_species, num_qks, num_qjs, num_qjo, num_sad
      integer, intent(in) :: num_chem, num_active, num_molefrac
      integer, intent(in) :: num_rxnr_adjust, rxnr_adjust_timpyr
      integer, intent(in) :: rxnr_adjust_map(num_rxnr_adjust)
      integer, intent(in) :: emiss_timpyr, curRecord
      integer, intent(in) :: ndust
      integer, intent(in) :: dehyd_opt, h2oclim_opt, lbssad_opt
      integer, intent(in) :: num_AerDust
      integer, intent(in) :: numDomains, commuWorld
      integer, intent(in) :: qj_timpyr, KDIM, JDIM, MDIM
      character (len=*) ,intent(in) :: chem_mecha
      character (len=MAX_LENGTH_MET_NAME) ,intent(in) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME) ,intent(in) :: metdata_name_model
      character (len=*) ,intent(in) :: qj_labels(num_qjs)
      character (len=MAX_LENGTH_SPECIES_NAME) ,intent(in) :: const_labels(num_species)
      character (len=MAX_LENGTH_FILE_NAME),intent(in) :: qj_infile_name
      character (len=MAX_LENGTH_FILE_NAME),intent(in) :: aqua_infile_name
      character (len=MAX_LENGTH_FILE_NAME),intent(in) :: linoz_infile_name
      character (len=MAX_LENGTH_FILE_NAME),intent(in) :: sf6_infile_name
      character (len=MAX_LENGTH_FILE_NAME),intent(in) :: SO3daily_infile_name
      character (len=MAX_LENGTH_FILE_NAME),intent(in) :: SO3monthly_infile_name
      character (len=MAX_LENGTH_VAR_NAME), intent(in) :: qj_var_name
      real*8 , intent(in   ) :: latdeg (ju1_gl:j2_gl)
      real*8  :: londeg  (i1_gl:i2_gl)
      real*8 , intent(in   ) :: dlatr(ju1_gl:j2_gl)
      real*8 , intent(  out) :: loss_freq(i1:i2, ju1:j2, k1:k2, num_species)
      real*8 , intent(in   ) :: humidity(i1:i2,   ju1:j2,   k1:k2) 
      real*8 , intent(in   ) :: taucli(i1:i2,   ju1:j2,   k1:k2) 
      real*8 , intent(in   ) :: tauclw(i1:i2,   ju1:j2,   k1:k2) 
      real*8 , intent(in   ) :: pctm2   (ilo:ihi, julo:jhi)  
      real*8 , intent(in   ) :: mcor    (i1:i2,   ju1:j2)    ! area of grid box  (m^2)
      real*8 , intent(in   ) :: mass    (i1:i2,   ju1:j2, k1:k2)
      real*8 , intent(inOut) :: qkgmi  (i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inOut) :: qjgmi  (i1:i2, ju1:j2, k1:k2, num_qjo)
      real*8 , intent(inOut) :: overheadO3col (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(inOut) :: decay_3d_out (i1:i2, ju1:j2, k1:k2, num_species)
      real*8 , intent(in   ) :: moistq     (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in   ) :: totalCloudFraction(i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in   ) :: fracCloudCover(i1:i2, ju1:j2)
      real*8 , intent(in   ) :: tau_cloud  (i1:i2, ju1:j2, k1:k2)
      integer, intent(in   ) :: lwi_flags  (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: surf_air_temp (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: surf_alb_uv   (i1:i2, ju1:j2)
      real*8 , intent(inOut) :: kel    (ilo:ihi, julo:jhi, k1:k2) 
      real*8 , intent(inOut) :: qqjgmi(i1:i2, ju1:j2, k1:k2, num_qjo)
      real*8 , intent(inOut) :: qqkgmi(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inOut) :: qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)
      real*8 , intent(inOut) :: qqkda(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inOut) :: yda  (i1:i2, ju1:j2, k1:k2, num_active)
      real*8 , intent(inOut) :: cloud_tau    (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(inOut) :: cloud_param  (i1:i2, ju1:j2, k1:k2,num_species)
      real*8 , intent(inOut) :: flux_gt      (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(inOut) :: prod(2, NUM_OPERATORS, 1:num_species)
      real*8 , intent(inOut) :: loss(2, NUM_OPERATORS, 1:num_species)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
      type (t_GmiArrayBundle), intent(inOut) :: emissionArray(num_emiss  )
      real*8  :: surf_emiss_out(i1:i2, ju1:j2, num_species)
      real*8  :: surf_emiss_out2(i1:i2, ju1:j2, 6)
      real*8  :: emiss_3d_out(i1:i2, ju1:j2, k1:k2, num_species)
      real*8  :: cmf     (i1:i2,  ju1:j2, k1:k2)
      real*8  :: clwc  (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: rxnr_adjust(i1:i2, ju1:j2, k1:k2, num_rxnr_adjust, rxnr_adjust_timpyr)
      real*8  :: dms_oh  (i1:i2, ju1:j2, k1:k2)
      real*8  :: dms_no3 (i1:i2, ju1:j2, k1:k2)
      real*8  :: so2_oh  (i1:i2, ju1:j2, k1:k2)
      real*8  :: so2_h2o2(i1:i2, ju1:j2, k1:k2)
      real*8  :: so2_o3  (i1:i2, ju1:j2, k1:k2)
      real*8  :: emiss_dust(i1:i2, ju1:j2, ndust)
      REAL*8 , intent(in) :: relativeHumidity (i1:i2, ju1:j2, k1:k2)
      REAL*8 , intent(in) :: radswg (i1:i2, ju1:j2)
      real*8 :: OptDepth(i1:i2, ju1:j2, k1:k2, num_AerDust)
      real*8 :: ERADIUS (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
      real*8 :: TAREA   (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
      REAL*8 :: ODAER   (i1:i2, ju1:j2, k1:k2, NSADaer*NRH_b)
      real*8 :: ODmdust (i1:i2, ju1:j2, k1:k2, NSADdust)
      real*8 , intent(in) :: Dust  (i1:i2,ju1:j2,k1:k2,NSADdust)
      real*8 , intent(in) :: Waersl(i1:i2,ju1:j2,k1:k2,NSADaer)
      real*8 , intent(in) :: Daersl(i1:i2,ju1:j2,k1:k2,2      )
      real*8 , intent(in)  :: lossData(KDIM, JDIM, MDIM, 1)
      logical, intent(in) :: pr_diag
      logical, intent(in) :: do_smv_reord, do_synoz, do_semiss_inchem
      logical, intent(in) :: do_AerDust_calc, do_clear_sky, do_ozone_inFastJX
      logical, intent(in) :: do_aqu_chem, do_wetchem
      integer, intent(in) :: loc_proc, AerDust_Effect_opt
      integer, intent(in) :: chem_opt, sad_opt, phot_opt, fastj_opt, be_opt, cloudDroplet_opt
      integer, intent(in) :: nymd, nhms, num_time_steps
      integer, intent(in) :: lbssad_timpyr, h2oclim_timpyr
      real*8 , intent(in) :: tdt
      real*8 , intent(out) :: dehydmin
      integer, intent(in) :: ih2_num, io3_num, ih2o_num, ih2ocond_num, ihno3_num, ihno3cond_num
      integer, intent(in) :: idehyd_num, imgas_num, ih2oaircr_num, ich4_num, ih2o2_num
      integer, intent(in) :: ico_num, ino_num, ipropene_num, iisoprene_num
      integer, intent(in) :: initrogen_num, ioxygen_num, isynoz_num
      integer, intent(in) :: chem_mask_klo, chem_mask_khi
      real*8 , intent(in) :: fastj_offset_sec
      integer, intent(in   ) :: forc_bc_kmax, forc_bc_kmin
      integer, intent(in   ) :: forc_bc_num
      integer, intent(in   ) :: forc_bc_opt
      integer, intent(in   ) :: forc_bc_years
      integer, intent(in   ) :: fbc_j1, fbc_j2
      integer, intent(in   ) :: forc_bc_map (forc_bc_num)
      real*8 , intent(in   ) :: forc_bc_incrpyr
      integer, intent(in   ) :: start_ymd, o3_index
      real*8 , intent(in   ) :: gmi_sec
      real*8 , intent(in   ) :: yield_be7, yield_be10
      real*8 , intent(in   ) :: t_half_be7, t_half_be10
      real*8 , intent(in   ) :: mw(num_species)
      logical, intent(in) :: do_ShipEmission
      integer, intent(in) :: loss_freq_opt, kmin_loss, kmax_loss
      integer, intent(in) :: emiss_opt, emiss_dust_opt
      integer, intent(in) :: emiss_map(num_emiss)
      integer, intent(in) :: emiss_map_dust(ndust)
      real*8,  intent(in) :: synoz_threshold
      real*8,  intent(in) :: chem_cycle
      real*8,  intent(in) :: jno_adjust
      integer, intent(in) :: jlatmd(j2_gl-ju1_gl+1)
      integer, intent(inOut) :: lastYearFBC
      real*8  :: forc_bc_data (:,:,:,:)
      real*8  :: ch4clim (i1:i2,   ju1:j2,   k1:k2, h2oclim_timpyr)
      real*8  :: h2oclim (i1:i2,   ju1:j2,   k1:k2, h2oclim_timpyr)
      real*8  :: hno3cond(i1:i2,   ju1:j2,   k1:k2)
      real*8  :: hno3gas (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: lbssad  (i1:i2,   ju1:j2,   k1:k2, lbssad_timpyr)
      real*8  :: sadgmi  (i1:i2,   ju1:j2,   k1:k2, num_sad)
      real*8  :: h2oback (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: h2ocond (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffice (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: reffsts (i1:i2,   ju1:j2,   k1:k2)
      real*8  :: vfall   (i1:i2,   ju1:j2,   k1:k2)
!
      logical, intent(in) :: rootProc
      logical, intent(in) :: pr_qj_o3_o1d, do_ftiming
      logical, intent(in) :: pr_qj_opt_depth
      logical, intent(in) :: do_qqjk_inchem
      logical, intent(in) :: do_qqjk_reset
      logical, intent(in) :: pr_qqjk
      logical, intent(in) :: pr_smv2
      logical, intent(in) :: pr_sulf_src
      logical, intent(in) :: do_aerocom
      logical, intent(in) :: pr_surf_emiss
      logical, intent(in) :: pr_emiss_3d
      real*8 , intent(in) :: pr_nc_period
      real*8 , intent(in) :: const_init_val(num_species)
      real*8 , intent(in) :: emiss_isop (i1:i2, ju1:j2)
      real*8 , intent(in) :: emiss_monot(i1:i2, ju1:j2)
      real*8 , intent(in) :: emiss_nox  (i1:i2, ju1:j2)
      real*8 , intent(in) :: emiss_hno3  (i1:i2, ju1:j2)
      real*8 , intent(in) :: tropp      (i1:i2, ju1:j2)
      real*8 , intent(in) :: press3e(ilo:ihi, julo:jhi, k1-1:k2)
      real*8 , intent(in) :: averagePressEdge (k1-1:k2)
      real*8 , intent(in) :: gridBoxHeight(i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in) :: chemintv
      real*8 , intent(in) :: rateintv 
      integer, intent(in) :: jno_num
      integer, intent(in) :: jnoxnum

      integer, intent(inOut) :: prevRecord
      integer, intent(inOut) :: forc_bc_start_num
      real*8 , intent(inOut) :: emiss_o3  (i1:i2, ju1:j2)
      real*8 , intent(inOut) :: emiss_ozone  (i1:i2, ju1:j2)
      real*8 , intent(inOut) :: press3c(ilo:ihi, julo:jhi, k1:k2)
      type(t_Smv2Saved), intent(inOut) :: smv2SavedVars
      type(t_SulfSaved), intent(inOut) :: sulfSavedVars
      type (t_fastJXbundle)  , intent(inOut) :: JXbundle
      type (t_CloudParametersGT)  , intent(inOut) :: cloudGT
  
!
! !DESCRIPTION: This routine updates the chemistry.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg

!.sds      integer :: IBe7S, IBe10S, IBe7T, IBe10T, Ie90, IRn, IPb, ICH3I
      integer :: ix
      integer :: nsec_jan1
      integer :: num_loops
!
! bnd for ShipEmission
      real*8        :: jno2val(i1:i2, ju1:j2)
      integer       :: ic
! bnd
!
      real*8  :: rsec_jan1

      integer :: idumday
      integer :: idumyear
      integer :: month
      real*8, allocatable  :: h2ocombo(:, :, :, :)
      real*8, allocatable  :: sadcombo(:, :, :, :)
      real*8, allocatable  :: solarZenithAngle(:,:)
      real*8  :: time_sec
      integer :: jday

!.sds temp fix to stop aerosol chem from updating h2o2 (and oh)
!      real*8, allocatable  :: h2o2gas(:, :, :)
!.sds end

!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'updateChemistry called by ', loc_proc

!     ==========
!      if (first) then
!!     ==========
!        first = .false.
!        chemintv = tdt * chem_cycle
!
!        if (chem_cycle < 1.0d0) then  ! subcycle chemistry
!           rateintv = tdt
!        else if (chem_cycle == 1.0d0) then
!           rateintv = chemintv
!!        end if
!
!!,... find reaction number of NO photolysis
!        do ic=1, num_qjs
!           if (TRIM(chem_mecha) == 'strat_trop' .or. &
!               TRIM(chem_mecha) == 'strat_trop_aerosol' ) then
!              if (Trim(qj_labels(ic)) =='NO + hv = N + O') jno_num = ic
!           endif
!        end do
!
!        if (do_ShipEmission) then
!           do ic=1, num_qjs
!              if (TRIM(chem_mecha) == 'strat_trop' .or. &
!                  TRIM(chem_mecha) == 'strat_trop_aerosol' ) then
!                 if (Trim(qj_labels(ic)) =='NO2 + hv = NO + O') jnoxnum = ic
!              else if (TRIM(chem_mecha) == 'troposphere') then
!                 if (Trim(qj_labels(ic)) =='NO2 + hv = NO + O3') jnoxnum = ic
!              endif
!           end do
!           if (jnoxnum.eq.999) then
!             err_msg = 'jnoxnum not found in updateChemistry'
!             call GmiPrintError  &
!     &         (err_msg, .true., 1, jnoxnum, 0, 0, 0.0d0, 0.0d0)
!           endif
!        endif
!!
!      end if

      call GetSecondsFromJanuary1 (nsec_jan1, nymd, nhms)

      rsec_jan1 = nsec_jan1

      !------------------------------
      ! Update Surface Area Densities
      !------------------------------

      call updateSurfaceAreaDensities (rateintv, tropp, press3c, press3e, kel,   &
     &                 concentration, ch4clim, h2oclim, hno3cond, hno3gas,       &
     &                 lbssad, sadgmi, h2oback, h2ocond, reffice, reffsts, vfall,&
     &                 dehydmin, dehyd_opt, h2oclim_opt, lbssad_opt, sad_opt,    &
     &                 ihno3cond_num, idehyd_num, ih2oaircr_num, ich4_num,       &
     &                 ihno3_num, ih2o_num, ih2ocond_num, nymd, pr_diag, loc_proc, &
     &                 commuWorld, numDomains, num_species, num_sad,             &
     &                 lbssad_timpyr, h2oclim_timpyr, ju1_gl, j2_gl, ilo, ihi,   &
     &                 julo, jhi, i1, i2, ju1, j2, k1, k2, metdata_name_model,   &
     &                 TRIM(chem_mecha))

      if ((phot_opt == 2) .or. (phot_opt == 3) .or. &
          (phot_opt == 4) .or. (phot_opt == 5)) then

        if (do_ftiming) call Ftiming_On  ('gmiPhotolysis')
 
        call GetDaysFromJanuary1 (jday, nymd)
        time_sec = ConvertTimeToSeconds (nhms) 

        ALLOCATE(   solarZenithAngle(i1:i2,ju1:j2))

        solarZenithAngle(i1:i2,ju1:j2) = &
             computeSolarZenithAngle_Photolysis (jday, time_sec, &
                    fastj_offset_sec, latDeg(ju1:j2), lonDeg(i1:i2), &
                    i1, i2, ju1, j2)

        !PRINT*,loc_proc,"CH2O before: ", &
        !       minval(concentration(1)%pArray3D(:,:,1)), &
        !       maxval(concentration(1)%pArray3D(:,:,1))

!       ================================
        call calcPhotolysisRateConstants  &
!       ================================
     &   ( JXbundle, pr_qj_o3_o1d, pr_qj_opt_depth, rateintv, rsec_jan1, pctm2,  &
     &     mass, press3e, press3c, kel, concentration, latdeg, londeg, &
     &     mcor, surf_alb_uv, fracCloudCover, tau_cloud, taucli, tauclw,&
     &     overheadO3col, qjgmi, gridBoxHeight, OptDepth, Eradius,     &
     &     Tarea, Odaer, relativeHumidity, Odmdust, Dust, Waersl,      &
     &     Daersl, humidity, lwi_flags, cloud_tau,     &
     &     cloud_param, flux_gt, num_AerDust, phot_opt, fastj_opt,            &
     &     fastj_offset_sec, do_clear_sky, do_AerDust_Calc,                   &
     &     do_ozone_inFastJX, qj_timpyr, io3_num, ih2o_num, isynoz_num,       &
     &     chem_mask_khi, nymd, nhms, pr_diag, loc_proc, synoz_threshold,     &
     &     AerDust_Effect_opt, num_species, num_qjs, num_qjo, ilo, ihi, julo, &
     &     jhi, i1_gl, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, j2, k1, k2, TRIM(chem_mecha), &
     &     jno_num, jno_adjust, solarZenithAngle, rootProc)

        DEALLOCATE(solarZenithAngle)

        if (do_ftiming) call Ftiming_Off  ('gmiPhotolysis')

      end if

!     ======================
      select case (chem_opt)
!     ======================

!       ========
        case (1)
!       ========
          IRn = 1
          IPb = 2

!         ======================
          call Update_Radon_Lead  &
!         ======================
     &      (lwi_flags, mcor, mass, surf_air_temp, concentration, IRn, IPb, &
     &       pr_diag, loc_proc, &
     &       tdt, i1, i2, ju1, j2, k1, k2, num_species)

!       ========
!       case (2)
!       ========

!         See below.

!       ========
        case (3)
!       ========

!         ================
          call updateStratosphericLoss  &
!         ================
     &      (chemintv, dlatr, loss_freq, lossData, concentration, averagePressEdge, &
     &         loss_freq_opt, nymd, kmin_loss, kmax_loss,   &
     &         const_init_val, num_species, pr_diag, loc_proc, &
     &         i1, i2, ju1, j2, k1, k2, ivert, ju1_gl, j2_gl, KDIM, JDIM, MDIM)

!       ========
        case (4)
!       ========

!         ===================
          call updateForcingBC  &
!         ===================
     &      (concentration, nymd, gmi_sec, fbc_j1, fbc_j2, &
     &       forc_bc_num, forc_bc_kmax, forc_bc_kmin, forc_bc_opt, &
     &       forc_bc_map, forc_bc_incrpyr, forc_bc_start_num, forc_bc_years, &
             forc_bc_data, lastYearFBC, jlatmd, ju1_gl, j2_gl, i1, i2, &
             ju1, j2, k1, k2, num_species)

!       ========
        case (5)
!       ========

!         See Update_Synspc.

!       ========
        case (6)
!       ========
          IBe7  = 1
          IBe10 = 2
!         =================
          call Update_Beryl  &
!         =================
     &      (latdeg, press3c, concentration, IBe7, IBe10, &
     &       pr_diag, loc_proc, tdt, be_opt, t_half_be7, t_half_be10, &
     &       yield_be7, yield_be10, i1, i2, ju1, j2, k1, k2, &
     &       ilo, ihi, julo, jhi, ju1_gl, j2_gl, num_species)

!       ========
!       case (7 and 8)
!       ========

!         See below.

!       ========
        case (9)
!       ========


!... Age, e90, tm25, and clock done with tracer_opt=2 in calcTaggedCO_AgeOfAir


!... Synoz done with do_synoz


!!... Radon/Lead
!!         ======================
!          call Update_Radon_Lead  &
!!         ======================
!     &      (lwi_flags, mcor, mass, surf_air_temp, concentration, IRn, IPb, &
!     &       pr_diag, loc_proc, &
!     &       tdt, i1, i2, ju1, j2, k1, k2, num_species)
!
!!... Beryllium
!... Beryllium
!          Ie90   = 2
!          IBe7S  = 7
!          IBe10S = 8
!          IBe7T  = 9
!          IBe10T = 10
!!         =================
!          call Update_Beryl_2  &
!!         =================
!     &      (latdeg, press3c, concentration, IBe7S, IBe10S, IBe7T, IBe10T, &
!     &       pr_diag, loc_proc, tdt, be_opt, t_half_be7, t_half_be10, &
!     &       yield_be7, yield_be10, Ie90, i1, i2, ju1, j2, k1, k2, &
!     &       ilo, ihi, julo, jhi, ju1_gl, j2_gl, num_species)

!... merged Beryl and Radon/Lead codes (hyl, 08/30/06)
!         ======================
           call Update_Radon_Lead_Beryl  &
!         ======================
     &       (latdeg, tropp, press3c, lwi_flags, mcor, mass, mw, &
     &        surf_air_temp, concentration, IRn, IPb, IPbS, IBe7, IBe7S, IBe10, IBe10S, &
     &        decay_3d_out, &   !jules,hyl
     &        pr_surf_emiss, pr_emiss_3d, surf_emiss_out, emiss_3d_out, &      !hyl
     &        pr_diag, loc_proc, &
     &        tdt, be_opt, t_half_be7, t_half_be10, &
     &        yield_be7, yield_be10, &
     &        i1, i2, ju1, j2, k1, k2, &
     &        ilo, ihi, julo, jhi, ju1_gl, j2_gl, num_species)


!... CH3I
!          ICH3I = 11
!         =================
          call Update_CH3I  &
!         =================
     &      (lwi_flags, mcor, mass, concentration, ICH3I, &
     &       pr_diag, loc_proc, tdt, i1, i2, ju1, j2, ilo, ihi, julo, jhi, k1, k2, num_species)


!... fCO2 passive tracer with emissions
           

!... linoz should be done with do_linoz?
!... 
!         =================
          call Update_LINOZ  &
!         =================
           (linoz_infile_name, ILINOZ, concentration, nymd, nhms, tdt, londeg, latdeg,  &
            mcor, pr_diag, loc_proc, kel, mass, gridBoxHeight, press3c,  &
            averagePressEdge, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, &
            i1_gl, i2_gl, ju1_gl, j2_gl, num_species)


!... SF6 tracer with emissions and parameterized chemistry
!         =================
          call Update_SF6  &
!         =================
           (SF6_infile_name, ISF6, concentration, nymd, tdt, londeg, latdeg,  &
            averagePressEdge, pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2,  &
            ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl, num_species)


!... Uniform tracer with emissions and parameterized chemistry
!         =================
          call Update_Uniform  &
!         =================
           (mcor, mass, concentration, IUniform, tdt,  &
            pr_diag, loc_proc, i1, i2, ju1, j2, ilo, ihi, julo, jhi, k1, k2, num_species)


!... StratO3 tracer with emissions and parameterized chemistry
!         =================
          call Update_stratO3  &
!         =================
           (SO3daily_infile_name, SO3monthly_infile_name, Istrat_O3, Ie90,  &
            concentration, humidity, kel, press3c, nymd, tdt, londeg, latdeg,  &
            pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2,  &
            ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl, num_species)

!       ==============
        case (2, 7, 8)  ! do_full_chem
!       ==============

!         -----------------------------------------------------
!         Impose the forcing boundary conditions if they exist.
!         -----------------------------------------------------

          if (forc_bc_num > 0) then

            call doDiagnosticsBefore_BC(concentration, prod, loss, mass,  &
     &                 FORCED_BC_OP, i1, i2, ju1, j2, k1, k2, num_species)

!           ===================
            call updateForcingBC  &
!           ===================
     &        (concentration, nymd, gmi_sec, fbc_j1, fbc_j2, &
     &         forc_bc_num, forc_bc_kmax, forc_bc_kmin, forc_bc_opt, &
     &         forc_bc_map, forc_bc_incrpyr, forc_bc_start_num, forc_bc_years, &
     &         forc_bc_data, lastYearFBC, jlatmd, ju1_gl, j2_gl, i1, i2,&
               ju1, j2, k1, k2, num_species)

            call doDiagnosticsAfter_BC(concentration, prod, loss, mass,  &
     &                 FORCED_BC_OP, i1, i2, ju1, j2, k1, k2, num_species)
          end if

!         --------------------------------------------------------
!         Calculate the air density at the center of each grid box
!         (molecules/cm^3).
!         --------------------------------------------------------

          concentration(imgas_num)%pArray3D(:,:,:) =  &
     &      press3c(i1:i2,ju1:j2,:) * MB2CGS /  &
     &      (kel (i1:i2,ju1:j2,:) * BOLTZMN_E)

          if (ioxygen_num /= 0) then
            concentration(ioxygen_num)%pArray3D(:,:,:) =  &
     &        concentration(imgas_num)%pArray3D(:,:,:) * MXRO2
          end if

          if (initrogen_num /= 0) then
            concentration(initrogen_num)%pArray3D(:,:,:) =  &
     &        concentration(imgas_num)%pArray3D(:,:,:) * MXRN2
          end if
         
          !!!!! Fixes for H2 - Provided by David Considine
          if (ih2_num /= 0) then
             concentration(ih2_num)%pArray3D(:,:,:) = MXRH2
          end if 
          !!!!! end fixes for H2 
!
!         For Ship Emission computations
          if ((do_ShipEmission) .and. (do_semiss_inchem)) then
             jno2val(:,:) = qjgmi(:,:,1,jnoxnum)
             call calcShipEmission (emiss_o3, emiss_ozone, prevRecord, &
                      curRecord, latdeg, jno2val, emissionArray, &
                      o3_index, i1, i2, ju1, j2, k1, k2, ju1_gl, j2_gl, &
                      num_emiss)
          end if

!         ==================
          if (chem_opt == 2) then
!         ==================
              
!.sds... if aerosols are calculated in SMVGEAR, then need to convert from mmr to vmr              
            if ( TRIM(chem_mecha) == 'strat_trop_aerosol' ) then
               do ic = 1, nact
                  if(aerosol(ic) .ne. 0) then
                     concentration(ic)%pArray3D(:,:,:) = ( mw(imgas_num)/mw(ic) ) &
                                                         * concentration(ic)%pArray3D(:,:,:)
                  endif
               enddo
            endif

!           ==============
!            call Update_Qk  &
            call calcThermalRateConstants  &
!           ==============
     &        (do_wetchem, metdata_name_org, metdata_name_model,  &
     &         num_time_steps, ih2o_num, imgas_num, nymd,  &
     &         rxnr_adjust_map, latdeg, press3c, tropp, kel, clwc,  &
     &         cmf, sadgmi, qkgmi, concentration, rxnr_adjust, &
     &         Eradius, Tarea, relativeHumidity, do_AerDust_Calc,  &
     &         phot_opt, pr_diag, loc_proc, num_rxnr_adjust, rxnr_adjust_timpyr, &
     &         ivert, num_sad, num_qks, num_molefrac, num_species, &
     &         ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &         TRIM(chem_mecha))

!           ------------------------------------------------------
!           Just hno3gas (i.e., no hno3cond) is used by chemistry.
!           ------------------------------------------------------

            if ((sad_opt == 1) .or. (sad_opt == 2)) then
              concentration(ihno3_num)%pArray3D(:,:,:) = hno3gas(:,:,:)
            end if

            if (chem_opt == 2) then

              if (chem_cycle < 1.0d0) then
                num_loops = Nint (1.0d0 / chem_cycle)
              else
                num_loops = 1
              end if

              do ix = 1, num_loops
!               ====================
                call Update_Smv2chem  &
!               ====================
     &            (smv2SavedVars, chemintv, emiss_isop, emiss_monot, emiss_nox,  &
     &             do_ShipEmission, emiss_hno3, emiss_o3, ihno3_num, io3_num, &
     &             mcor, surf_emiss_out, surf_emiss_out2, emiss_3d_out,  &
     &             humidity, qjgmi, qkgmi, emissionArray,  &
     &             press3e, pctm2, kel, concentration, &
     &             pr_diag, pr_qqjk, pr_smv2, pr_surf_emiss, pr_emiss_3d, &
     &             do_smv_reord, do_synoz, do_qqjk_inchem, do_semiss_inchem, &
     &             ico_num, ino_num, ipropene_num, iisoprene_num, &
     &             imgas_num, initrogen_num, ioxygen_num, isynoz_num, &
     &             yda, qqkda, qqjda, mw, pr_nc_period, &
     &             emiss_timpyr, emiss_opt, emiss_map, tdt, nymd, &
     &             chem_mask_klo, chem_mask_khi, &
     &             loc_proc, synoz_threshold, &
     &             ilong, ilat, ivert, itloop,  &
     &             i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, &
     &             num_molefrac, num_emiss, num_qjo, num_qks, num_qjs, &
     &             num_active, num_species, commuWorld, TRIM(chem_mecha))
              end do
              
        !PRINT*,loc_proc,"CH2O after: ", &
        !       minval(concentration(1)%pArray3D(:,:,1)), &
        !       maxval(concentration(1)%pArray3D(:,:,1))

!!.sds... if aerosols are calculated in SMVGEAR, then need to convert from vmr to mmr              
              if ( TRIM(chem_mecha) == 'strat_trop_aerosol' ) then
                 do ic = 1, nact
                    if(aerosol(ic) .ne. 0) then
                       concentration(ic)%pArray3D(:,:,:) = ( mw(ic)/mw(imgas_num) ) &
                                                           * concentration(ic)%pArray3D(:,:,:)
                    endif
                 enddo
              endif

!.sds... do sulfur chem if species nDMS exists
              if(INDMS .ne. 0) then
              
!              print *,'sds-Calling Update_Sulfchem'

!               ====================
                call Update_Sulfchem  &
!               ====================
     &        (cloud_param, chemintv, londeg, latdeg, kel, concentration, humidity, emissionArray,  &
     &         mcor, qjgmi, qkgmi, moistq, cmf, press3c, gridBoxHeight,  &
     &         mass, dms_oh, dms_no3, so2_oh, so2_h2o2, so2_o3, emiss_dust, &
     &         aqua_infile_name, qj_infile_name, qj_var_name, &
     &         pr_diag, pr_sulf_src, do_aerocom, do_clear_sky, &
     &         imgas_num, emiss_dust_opt, emiss_map, emiss_map_dust, mw, &
     &         phot_opt, &
     &         nymd, nhms, num_time_steps, loc_proc, ndust, ilong, itloop, &
     &         i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &         emiss_timpyr, num_molefrac, num_emiss, num_qjs, num_qks, num_species, chem_opt)

             endif

!.sds end

!.sds... do ammonia chem if species index NH3 exists
              if(INH3 .ne. 0) then
              
!              print *,'sds-Calling Do_NH3_Solver'

!               ==================
                call Do_NH3_Solver  &
!               ==================
     &            (loc_proc, pr_diag, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi,  &
     &             kel, relativeHumidity, concentration)

              endif
!.sds end


              if (pr_qqjk .and. (.not. do_qqjk_inchem)) then
!               ===============
                call Accum_Qqjk  &
!               ===============
     &          (do_qqjk_reset, imgas_num, concentration, qjgmi, qkgmi,  &
     &           qqjgmi, qqkgmi, &
     &           num_molefrac, num_species, num_qks, num_qjs, num_qjo, &
     &           pr_diag, loc_proc, ilong, i1, i2, ju1, j2, k1, k2)
              end if

            end if

!           -----------------------------------------
!           Add hno3cond to hno3gas for transporting
!             when not using predicted H2O and CH4.
!           -----------------------------------------

            IF((sad_opt == 1 .OR. sad_opt == 2) .AND. h2oclim_opt /= 3) THEN
              concentration(ihno3_num)%pArray3D(:,:,:) =  &
     &          concentration(ihno3_num)%pArray3D(:,:,:) + hno3cond(:,:,:)
            END IF

!         =======================
          else if (chem_opt == 7) then
!         =======================

!           ====================
            call Update_Quadchem  &
!           ====================
     &        (chemintv, rsec_jan1, latdeg, kel, concentration,  &
     &         humidity, pctm2, mcor, mass, qjgmi, qkgmi, moistq,  &
     &         gridBoxHeight, totalCloudFraction, press3c, lwi_flags, &
     &         do_aqu_chem, do_clear_sky, pr_diag, loc_proc, &
     &         num_qjs, num_qks, num_chem, num_active, num_species, num_molefrac, &
     &         imgas_num, ih2o_num, io3_num, i1, i2, ju1, j2, k1, k2, &
     &         ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl)

!         =======================
          else if (chem_opt == 8) then
!         =======================

#ifdef GTmodule
            call calcCloudPropertiesGT  &
     &       ( cloudGT, londeg, latdeg, kel, concentration, humidity,  &
     &         moistq, cmf, press3c, gridBoxHeight, mass, lwi_flags, &
     &         totalCloudFraction,  tau_cloud, radswg, &
     &         cloud_param, cloud_tau, cloudDroplet_opt, &
     &         pr_diag, loc_proc, num_species, imgas_num, &
     &         i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl)

!           ================================
            call calcPhotolysisRateConstants  &
!           ================================
     &       ( JXbundle, pr_qj_o3_o1d, pr_qj_opt_depth, rateintv, rsec_jan1, pctm2,  &
     &         mass, press3e, press3c, kel, concentration, latdeg, londeg, mcor, surf_alb_uv,  &
     &         fracCloudCover, tau_cloud, overheadO3col, qjgmi, gridBoxHeight, &
     &         OptDepth, Eradius, Tarea, Odaer, relativeHumidity, Odmdust, Dust, &
     &         Waersl, Daersl, humidity, lwi_flags, cloud_tau, cloud_param, flux_gt, &
     &         num_AerDust, phot_opt, fastj_opt, fastj_offset_sec, &
     &         do_clear_sky, do_AerDust_Calc, do_ozone_inFastJX, &
     &         qj_timpyr, io3_num, ih2o_num, isynoz_num, chem_mask_khi, &
     &         nymd, nhms, pr_diag, loc_proc, synoz_threshold, AerDust_Effect_opt, &
     &         num_species, num_qjs, num_qjo, ilo, ihi, julo, jhi, i1_gl, i2_gl, &
     &         ju1_gl, j2_gl, i1, i2, ju1, j2, k1, k2, TRIM(chem_mecha), &
     &         jno_num, jno_adjust)
#endif

!           ====================
            call Update_Sulfchem  &
!           ====================
     &       ( cloud_param, chemintv, londeg, latdeg, kel, concentration, humidity, emissionArray,  &
     &         mcor, qjgmi, qkgmi, moistq, cmf, press3c, gridBoxHeight,  &
     &         mass, dms_oh, dms_no3, so2_oh, so2_h2o2, so2_o3, emiss_dust, &
     &         aqua_infile_name, qj_infile_name, qj_var_name, &
     &         pr_diag, pr_sulf_src, do_aerocom, do_clear_sky, &
     &         imgas_num, emiss_dust_opt, emiss_map, emiss_map_dust, mw, &
     &         phot_opt, &
     &         nymd, nhms, num_time_steps, loc_proc, ndust, ilong, itloop, &
     &         i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &         emiss_timpyr, num_molefrac, num_emiss, num_qjs, num_qks, num_species, chem_opt)

          end if

!       ============
        case DEFAULT
!       ============

          err_msg = 'chem_opt problem in updateChemistry.'
          call GmiPrintError  &
     &      (err_msg, .true., 1, chem_opt, 0, 0, 0.0d0, 0.0d0)

!     ==========
      end select
!     ==========

      return

      end subroutine updateChemistry
!EOC
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update_Quadchem
!
! !INTERFACE:
!
      subroutine Update_Quadchem  &
     &  (chemintv, rsec_jan1, latdeg, kel1, concentration, humidity,  &
     &   pctm2, mcor, mass, hv_gmi, ratek_gmi, moistq, gridBoxHeight,  &
     &   totalCloudFraction, press3c, lwi_flags, &
     &   do_aqu_chem, do_clear_sky, pr_diag, loc_proc, &
     &   num_qjs, num_qks, num_chem, num_active, num_species, num_molefrac, &
     &   imgas_num, ih2o_num, io3_num, i1, i2, ju1, j2, k1, k2, &
     &   ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl)
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!
      implicit none

#     include "gmi_time_constants.h"
#     include "gmi_phys_constants.h"
!
! !INPUT PARAMETERS:
      logical, intent(in   ) :: pr_diag
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: num_qjs, num_qks, num_chem
      integer, intent(in   ) :: num_active, num_species, num_molefrac
      integer, intent(in   ) :: imgas_num, ih2o_num, io3_num
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in   ) :: ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl
      logical, intent(in   ) :: do_aqu_chem, do_clear_sky
      real*8 , intent(in   ) :: chemintv                   ! chemistry time step   (s)
      real*8 , intent(in   ) :: rsec_jan1                  ! seconds from Jan. 1st (s)
      real*8 , intent(in   ) :: latdeg  (ju1_gl:j2_gl)     ! latitude    (deg)
      real*8 , intent(in   ) :: humidity(i1:i2,   ju1:j2,   k1:k2)    ! specific humidity (g/kg)
      real*8 , intent(in   ) :: pctm2   (ilo:ihi, julo:jhi)  ! CTM surface pressure at t1+tdt (mb)
      real*8 , intent(in   ) :: mcor    (i1:i2,   ju1:j2)    ! area of grid box  (m^2)
      real*8 , intent(in   ) :: mass    (i1:i2,   ju1:j2, k1:k2)
!                             ! total mass of the atmosphere within each grid box (kg)
      real*8 , intent(in   ) :: ratek_gmi  (i1:i2, ju1:j2, k1:k2, num_qks)
!                             ! reaction rate constants     (units vary)
      real*8 , intent(in   ) :: moistq     (i1:i2, ju1:j2, k1:k2)
!                             ! moisture changes due to wet processes (g/kg/day)
      real*8 , intent(in   ) :: gridBoxHeight(i1:i2, ju1:j2, k1:k2)
!                             ! height of each grid box (m)
      real*8 , intent(in   ) :: totalCloudFraction  (i1:i2, ju1:j2, k1:k2)
      integer, intent(in   ) :: lwi_flags  (i1:i2, ju1:j2)
!                             ! array of flags that indicate land, water, or ice
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inOut) :: kel1    (ilo:ihi, julo:jhi, k1:k2)    ! temperature (degK)
      real*8 , intent(inOut) :: hv_gmi  (i1:i2, ju1:j2, k1:k2, num_qjs)
!                             ! variables used in chemistry (units vary)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!                             ! species concentration (mixing ratio)
      real*8 , intent(inOut) :: press3c (ilo:ihi, julo:jhi, k1:k2)
!                             ! atmospheric pressure at the center of each grid box (mb)
!
! !DESCRIPTION:
!  This routine updates the quadchem chemistry.
!
! !REVISION HISTORY:
!   Initial code.
!EOP
!-------------------------------------------------------------------------

!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter ::  &
     &  R_61  = (1.0d0 - (MWTH2O / MWTAIR)) / (MWTH2O / MWTAIR)

      real*8, parameter ::  &
     &  R29_3 = (GAS_CONST_J / MWTAIR * GPKG) / GMI_G

      real*8, parameter ::  &
     &  DUCONST = AVOGAD / (MWTAIR * 2.69d+16 * 10.0d0),  &
     &  EXT_EFF = 2.0d0,     & ! extinction efficiency of a water droplet
     &  LWD     = 1.0d+06  ! liquid water densitry (g/m^3)


!     ----------------------
!     Variable declarations.
!     ----------------------

      logical, save :: cloud_top = .true.

      integer :: il, ij, ik, ic

!     ---------------------------------------
!     height_asl : height above sea level (m)
!     ---------------------------------------

      real*8, allocatable :: height_asl(:)

!     ----------------------------------------------
!     cloud_opt_dep : cloud optical depth (unitless)
!     o3du          : ozone column (Dobson units)
!     ----------------------------------------------

      real*8, allocatable :: cloud_opt_dep(:,:)
      real*8, allocatable :: o3du(:,:)

!     ----------------------------------------------------
!     eff_rad     : cloud drop effective radius (m)
!     cloud_depth : depth of cloud (m)
!     lwccol      : liquid water in each grid box (gm/m^3)
!     zq          : height at full sigma levels (m)
!     ----------------------------------------------------

      real*8, allocatable :: cloud_depth(:,:,:)
      real*8, allocatable :: eff_rad    (:,:,:)
      real*8, allocatable :: lwccol     (:,:,:)
      real*8, allocatable :: zq         (:,:,:)

!     ===============================================================
!cc!! Non-standard F90 below!?
!     Only way I could find to share memory between const & const_nb;
!     hopefully it will work on all the machines we run on.
!     Also see below.
!     ===============================================================

      real*8  :: const_nb
      real*8  :: tempc_nb
      real*8  :: presc_nb

      pointer (ptr_const_nb,  &
     &         const_nb(i1:i2, ju1:j2, k1:k2, num_species))
      pointer (ptr_tempc_nb,  &
     &         tempc_nb(i1:i2, ju1:j2, k1:k2))
      pointer (ptr_presc_nb,  &
     &         presc_nb(i1:i2, ju1:j2, k1:k2))


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_Quadchem called by ', loc_proc
      end if


      Allocate (height_asl(k1:k2))
      height_asl = 0.0d0

      Allocate (cloud_opt_dep(i1:i2, ju1:j2))
      Allocate (o3du         (i1:i2, ju1:j2))
      cloud_opt_dep = 0.0d0; o3du = 0.0d0

      Allocate (cloud_depth(i1:i2, ju1:j2, k1:k2))
      Allocate (eff_rad    (i1:i2, ju1:j2, k1:k2))
      cloud_depth = 0.0d0; eff_rad = 0.0d0

      Allocate (lwccol(i1:i2, ju1:j2, k1:k2))
      Allocate (zq    (i1:i2, ju1:j2, k1-1:k2))
      lwccol = 0.0d0; zq = 0.0d0


!     ---------------------------------------
!     Estimate ozone column, in Dobson units.
!     ---------------------------------------

      do il = i1, i2
        do ij = ju1, j2
          do ik = k2, k1, -1
            o3du(il, ij) =  &
     &        o3du(il, ij) +  &
     &        concentration(io3_num)%pArray3D(il, ij, ik) * mass(il, ij, ik) *  &
     &        DUCONST / mcor(il, ij)
          end do
        end do
      end do


!     =================
      if (do_clear_sky) then
!     =================

        hv_gmi(:,:,:,2) = 1.0d0
        hv_gmi(:,:,:,3) = 0.0d0

!     ====
      else
!     ====

        hv_gmi(i1:i2, ju1:j2, k1:k2, 5) = totalCloudFraction(i1:i2, ju1:j2, k1:k2)

!       ----------------------------------------------------------
!       Compute cloud drop effective radius from Kiehl, JGR, 1994.
!       ----------------------------------------------------------

!       -----------------------------------------
!       Ocean, sea ice, and land (kel1 < 243.15).
!       -----------------------------------------

        eff_rad(:,:,:) = 1.0d-5

        do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2

              if (lwi_flags(il,ij) == 2) then  ! land

                if ((kel1(il,ij,ik) >= 243.15d0) .and.  &
     &              (kel1(il,ij,ik) <= 263.15d0)) then

                  eff_rad(il,ij,ik) =  &
     &              5.0d-6 -  &
     &              5.0d-6 *  &
     &              ((kel1(il,ij,ik) + ABS_ZERO + 10.d0) / 20.d0)

                else if (kel1(il,ij,ik) > 263.15d0) then

                  eff_rad(il,ij,ik) = 5.0d-6

                end if

              end if

            end do
          end do
        end do


!       ----------------------------------------------------
!       Estimate cloud above/below index and depth of cloud.
!       ----------------------------------------------------

        do il = i1, i2
          do ij = ju1, j2

            cloud_top = .false.

            do ik = k2, k1, -1

              if (moistq(il,ij,ik) < -1.0d-5) then
                cloud_top = .true.
                cloud_depth(il,ij,ik) = gridBoxHeight(il,ij,ik)
              end if

              if (cloud_top) then

                if (moistq(il,ij,ik) < -1.0d-5) then

!                 --------------------------------------
!                 Cloud above/below index in the middle.
!                 --------------------------------------

                  hv_gmi(il,ij,ik,2) = 0.5d0

                else

!                 ----------------------------------------
!                 Cloud above/below index above the cloud.
!                 ----------------------------------------

                  hv_gmi(il,ij,ik,2) = 1.0d0

                end if

              end if

            end do
          end do
        end do


!       -------------------------------------------------------------------
!       Calculate cloud water from a parameterization of Kiehl found in
!       "Sensitivity of the Simulated Climate to a Diagnostic Formulation
!       for Cloud Liquid Water" by James Hack, Journal of Climate, Vol. 11,
!       July 1998, p 1499.
!       -------------------------------------------------------------------

        do ij = ju1, j2
          do il = i1, i2

            height_asl(:) =  &
     &        -29.3d0 * kel1(il,ij,:) *  &
     &        Log (press3c(il,ij,:) / 1013.25d0)

            height_asl(:) = Max (height_asl(:), 0.0d0)

            lwccol(il, ij, :) =  &
     &        0.18d0 *  &
     &        Exp (-height_asl(:) /  &
     &             (1080.0d0 +  &
     &              (2000.0d0 * Cos (latdeg(ij) * RADPDEG)** 2.0d0)))

          end do
        end do


!       -----------------------------------------
!       In-cloud liquid water mixing ratio (g/g).
!       -----------------------------------------

        if (do_aqu_chem) then

          do ik = k1, k2
            do ij = ju1, j2
              do il = i1, i2

                hv_gmi(il,ij,ik,6) =  &
     &            lwccol(il,ij,ik) * 1.0d-6 /  &
     &            (concentration(imgas_num)%pArray3D(il,ij,ik) * MWTAIR / AVOGAD)

              end do
            end do
          end do

        end if


!       -----------------------------------------
!       Grid-averaged cloud liquid water (g/m^3).
!       -----------------------------------------

        do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2

              lwccol(il,ij,ik) =  &
     &          lwccol(il,ij,ik) * hv_gmi(il,ij,ik,5)

            end do
          end do
        end do


!       ---------------------------------------------------------------
!       Estimate cloud optical depth (Seinfeld & Pandis, p.1173, 1998).
!       ---------------------------------------------------------------

        do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2

              cloud_opt_dep(il,ij) =  &
     &          cloud_opt_dep(il,ij) +  &
     &          (3.0d0 * cloud_depth(il,ij,ik) *  &
     &           lwccol(il,ij,ik) * EXT_EFF) /  &
     &          (4.0d0 * LWD * eff_rad(il,ij,ik))

            end do
          end do
        end do


        do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2
              hv_gmi(il,ij,ik,3) = cloud_opt_dep(il,ij)
            end do
          end do
        end do

!     ======
      end if
!     ======


!     ===============================================================
!cc!! Non-standard F90 below!?
!     Only way I could find to share memory between const & const_nb;
!     hopefully it will work on all the machines we run on.
!     ===============================================================

      ptr_const_nb = Loc (concentration(1)%pArray3D(ilo,julo,k1))
      ptr_tempc_nb = Loc (kel1   (ilo,julo,k1))
      ptr_presc_nb = Loc (press3c(ilo,julo,k1))


!     ------------------------------------------------------
!     Estimate altitude at the center of each grid box (km).
!     ------------------------------------------------------

!     -------------------------------------------------
!     Compute the height of full and half-sigma levels.
!     -------------------------------------------------

      do ik = k1, k2
        zq(i1:i2,ju1:j2,ik) =  &
     &    zq(i1:i2,ju1:j2,ik-1) +  gridBoxHeight(i1:i2,ju1:j2,ik)
      end do


      do ik = k1, k2
        do ij = ju1, j2
          do il = i1, i2

            hv_gmi  (il,ij,ik,1) = o3du(il,ij)
            tempc_nb(il,ij,ik)   = kel1(il,ij,ik)
            presc_nb(il,ij,ik)   = press3c(il,ij,ik)

            hv_gmi(il,ij,ik,4) =  &
     &          0.5d0 * (zq(il,ij,ik) + zq(il,ij,ik-1))

          end do
        end do
      end do


!     ================
      if (do_aqu_chem) then
!     ================

!     -----------------------------------------------------------
!     Calculate relative humidity from Seinfeld (1986) p. 181.
!     The first rh is the temperature dependent parameter a.
!     The second rh is the saturation vapor pressure of water.
!     The third rh is the actual relative humidity as a fraction.
!     Then make sure rh is between 0 and 1.
!     -----------------------------------------------------------

        hv_gmi(:,:,:,7) =  &
     &    1.0d0 - (373.15d0 / tempc_nb(:,:,:))

        hv_gmi(:,:,:,7) =  &
     &    1013.25d0 * Exp (13.3185d0 * hv_gmi(:,:,:,7)    -  &
     &                      1.9760d0 * hv_gmi(:,:,:,7)**2 -  &
     &                      0.6445d0 * hv_gmi(:,:,:,7)**3 -  &
     &                      0.1299d0 * hv_gmi(:,:,:,7)**4)

        hv_gmi(:,:,:,7) =  &
     &    concentration(ih2o_num)%pArray3D(i1:i2,ju1:j2,k1:k2) *  &
     &    presc_nb(:,:,:) / hv_gmi(:,:,:,7)

        hv_gmi(:,:,:,7) =  &
     &    Max (Min (hv_gmi(:,:,:,7), 1.0d0), 0.0d0)

!     ======
      end if
!     ======


!     ------------------------------------------------
!     Change units from mixing ratio to concentration.
!     ------------------------------------------------

      do ic = 1, num_molefrac + 1  ! +1 for CO2
        concentration(ic)%pArray3D(i1:i2,ju1:j2,:) =  &
     &    concentration(ic)%pArray3D(i1:i2,ju1:j2,:) *  &
     &    concentration(imgas_num)%pArray3D(i1:i2,ju1:j2,:)
      end do


      do ic = 1, num_chem
        do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2
              const_nb(il,ij,ik,ic) = concentration(ic)%pArray3D(il,ij,ik)
            end do
          end do
        end do
      end do


      do ic = num_active, 1, -1
        do ik = k2, k1, -1
          do ij = j2, ju1, -1
            do il = i2, i1, -1
              concentration(ic)%pArray3D(il,ij,ik) = const_nb(il,ij,ik,ic)
            end do
          end do
        end do
      end do


      do ik = k2, k1, -1
        do ij = j2, ju1, -1
          do il = i2, i1, -1
            kel1(il,ij,ik) = tempc_nb(il,ij,ik)
            press3c(il,ij,ik) = presc_nb(il,ij,ik)
          end do
        end do
      end do


!     -----------------------------------------------------
!     Change units from concentration back to mixing ratio.
!     -----------------------------------------------------

      do ic = 1, num_molefrac + 1   ! +1 for CO2
        concentration(ic)%pArray3D(i1:i2,ju1:j2,:) =  &
     &    concentration(ic)%pArray3D(i1:i2,ju1:j2,:) /  &
     &    concentration(imgas_num)%pArray3D(i1:i2,ju1:j2,:)
      end do


      Deallocate (cloud_depth)
      Deallocate (cloud_opt_dep)
      Deallocate (eff_rad)
      Deallocate (height_asl)
      Deallocate (lwccol)
      Deallocate (o3du)
      Deallocate (zq)


      return

      end subroutine Update_Quadchem

!-------------------------------------------------------------------------
      end module GmiUpdateChemistry_mod
