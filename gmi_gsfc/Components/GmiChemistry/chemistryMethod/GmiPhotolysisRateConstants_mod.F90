      module GmiPhotolysisRateConstants_mod

      use fastj              , only : Control_Fastj     , GetQAA_RAAinFastj
      use fast_JX            , only : Control_Fast_JX   , GetQAA_RAAinFastJX
      use fastJX_Bundle_mod,   ONLY : t_fastJXbundle
      use fastJX65_mod       , only : controlFastJX65,  getQAA_RAAinFastJX65
      use fast_JX53b         , only : Control_Fast_JX53b, GetQAA_RAAinFastJX53b
      use fast_JX53c         , only : Control_Fast_JX53c, GetQAA_RAAinFastJX53c
!      use FastJX53cMethod_mod, only : RunFastJX53c

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiAerDustODSA_mod      , only : Aero_OptDep_SurfArea
      use GmiAerDustODSA_mod      , only : Dust_OptDep_SurfArea
      use GmiTimeControl_mod      , only : GmiSplitDateTime
      use GmiTimeControl_mod      , only : GetDaysFromJanuary1
      use GmiTimeControl_mod      , only : ConvertTimeToSeconds
!      use GmiSolar_mod, only : computeSolarZenithAngle_Photolysis_2

      implicit none

      private
      public  :: calcPhotolysisRateConstants

      CONTAINS

!=============================================================================
!
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! HISTORY
!   - January 24, 2005 * Jules Kouatchou
!     Modified the routine Update_Qj to compute aerosol/dust
!     surface area and optical depth. The routines used for calculations
!     are Aero_OptDep_SurfArea and Dust_OptDep_SurfArea.
!     The changes are valid only for the troposphere chemical mechanism.
!     The argument "gridBoxHeight" was also added in the subroutine
!     Update_Qj for the above calculations.
!     Pass on the values of optical depht (ODAER_ij, ODMDUST_ij)
!     aerosols/dust to the fastj control routine.
!   - March 9, 2005 * Jules Kouatchou
!     The call of Control_Fastj can now be done for any chemical mechanism.
!     Two arguments (ODAER_ij, ODMDUST_ij) were added to Control_FastJX and
!     the variable "ozone_ij" is now an optional argument of Control_FastJX.
!     For this reason Control_FastJX is declared here as an interface.
!     "ozone_ij" is not used for the tropospheric chemical mechanism because
!     the code uses climatology provide by FastJX instead of the model
!     climatology. For the combo mechanism, there is a namelist variable
!     "do_ozone_inFastJX" for selecting a particular climatology.
!   - March 14, 2005 * Jules Kouatchou
!     Computations of optical depth and surface area for different
!     aerosols/dust are done only if "do_AerDust_Calc" is set to TRUE.
!=============================================================================


!-----------------------------------------------------------------------------
!
! DESCRIPTION
!   This routine updates the photolysis rate constants (i.e., qjs).
!
! ARGUMENTS
!   cross_section_file : X-Section quantum yield
!   rate_file          : Master rate file
!   T_O3_climatology_file : T & O3 climatology
!   metdata_name_org : first  part of metdata_name, e.g., "NCAR"
!   pr_qj_o3_o1d     : should special reaction O3->O1d be saved?
!   pr_qj_opt_depth  : should optical depth be saved at the end of qj file?
!   photintv         : photolysis time step  (s)
!   rsec_jan1        : seconds from Jan. 1st (s)
!   pctm2            : CTM surface pressure at t1+tdt (mb)
!   mass3            : total mass of the atmosphere within   each grid box (kg)
!   pres3c           : atmospheric pressure at the center of each grid box (mb)
!   temp3            : temperature (degK)
!   const            : species concentration, at zone centers (mixing ratio)
!   latdeg           : latitude  (deg)
!   londeg           : longitude (deg)
!   mcor             : area of grid box (m^2)
!   surf_alb_uv      : bulk surface albedo (fraction 0-1)
!   tau_cloud        : optical depth (dimensionless)
!   qjgmi            : photolysis rate constants (s^-1)
!   gridBoxHeight      : height of each grid box (m)
!
!-----------------------------------------------------------------------------

      subroutine calcPhotolysisRateConstants(JXbundle, &
     &   pr_qj_o3_o1d,  pr_qj_opt_depth, photintv, rsec_jan1, pctm2,  &
     &   mass3, pres3e, pres3c, temp3, concentration, latdeg, londeg, mcor,  &
     &   surf_alb_uv, fracCloudCover, &
     &   tau_cloud, taucli, tauclw, overheadO3col, qjgmi, gridBoxHeight, &
     &   OptDepth, Eradius, Tarea, Odaer, relativeHumidity, Odmdust, Dust, &
     &   Waersl, Daersl, humidity, lwi_flags, cloud_tau, cloud_param, flux_gt,  &
     &   num_AerDust, phot_opt, fastj_opt, fastj_offset_sec, &
     &   do_clear_sky, do_AerDust_Calc, do_ozone_inFastJX, &
     &   qj_timpyr, io3_num, ih2o_num, isynoz_num, chem_mask_khi, &
     &   nymd, nhms, pr_diag, loc_proc, synoz_threshold, AerDust_Effect_opt, &
     &   num_species, num_qjs, num_qjo, ilo, ihi, julo, jhi, i1_gl, i2_gl,   &
     &   ju1_gl, j2_gl, i1, i2, ju1, j2, k1, k2, chem_mecha, jno_num,  &
     &   jno_adjust, solarZenithAngle, rootProc)

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
#     include "phot_lookup.h"
#     include "phot_monthly.h"
#     include "setkin_par.h"
#     include "gmi_AerDust_const.h"


!     ----------------------
!     Argument declarations.
!     ----------------------
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: num_species, num_qjs, num_qjo
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: isynoz_num, io3_num, ih2o_num
      integer, intent(in) :: nymd, nhms
      real*8 , intent(in) :: synoz_threshold
      integer, intent(in) :: chem_mask_khi, qj_timpyr
      integer, intent(in) :: AerDust_Effect_opt
      logical, intent(in) :: do_AerDust_calc, do_clear_sky, do_ozone_inFastJX
      integer, intent(in) :: phot_opt, fastj_opt
      real*8 , intent(in) :: fastj_offset_sec
      integer, intent(in) :: jno_num
      real*8 , intent(in) :: jno_adjust
      character(len=*), intent(in) :: chem_mecha

      logical :: pr_qj_o3_o1d
      logical :: pr_qj_opt_depth
      real*8  :: photintv
      real*8  :: rsec_jan1
      integer :: lwi_flags(i1:i2, ju1:j2)
      real*8  :: pctm2 (ilo:ihi, julo:jhi)
      real*8  :: mass3 (i1:i2, ju1:j2, k1:k2)
      real*8  :: pres3c(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: pres3e(ilo:ihi, julo:jhi, k1-1:k2)
      real*8  :: temp3 (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: latdeg(ju1_gl:j2_gl)
      real*8  :: londeg(i1_gl:i2_gl)
      real*8  :: mcor  (i1:i2, ju1:j2)
      real*8  :: solarZenithAngle  (i1:i2, ju1:j2)
      real*8  :: surf_alb_uv(i1:i2, ju1:j2)
      real*8 , intent(in   ) :: fracCloudCover(i1:i2, ju1:j2)
      real*8  :: tau_cloud  (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in) :: taucli  (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in) :: tauclw  (i1:i2, ju1:j2, k1:k2)
      real*8  :: overheadO3col (i1:i2, ju1:j2, k1:k2)
      real*8  :: qjgmi (i1:i2, ju1:j2, k1:k2, num_qjo)
      real*8  :: gridBoxHeight (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in) :: humidity (i1:i2, ju1:j2, k1:k2)
      integer, intent(in) :: num_AerDust

      real*8 cloud_tau  (i1:i2, ju1:j2, k1:k2)
      real*8 cloud_param(i1:i2, ju1:j2, k1:k2,num_species)
      real*8 flux_gt    (i1:i2, ju1:j2, k1:k2)
      REAL*8 , intent(in) :: relativeHumidity (i1:i2, ju1:j2, k1:k2)
      real*8 :: OptDepth(i1:i2, ju1:j2, k1:k2, num_AerDust)
      real*8 :: ERADIUS (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
      real*8 :: TAREA   (i1:i2, ju1:j2, k1:k2, NSADdust+NSADaer)
      REAL*8 :: ODAER   (i1:i2, ju1:j2, k1:k2, NSADaer*NRH_b)
      real*8 :: ODmdust (i1:i2, ju1:j2, k1:k2, NSADdust)
      real*8 , intent(in) :: Dust  (i1:i2,ju1:j2,k1:k2,NSADdust)
      real*8 , intent(in) :: Waersl(i1:i2,ju1:j2,k1:k2,NSADaer)
      real*8 , intent(in) :: Daersl(i1:i2,ju1:j2,k1:k2,2      )
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)
      type (t_fastJXbundle)  , intent(inOut) :: JXbundle

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: idumday, idumyear
      integer :: il, ij, it
      integer :: jday
      integer :: month_gmi
      integer :: photOption

      real*8  :: time_sec

      integer :: latyp_ij
      real*8  :: overheadO3col_ij     (k1:k2)
      real*8  :: kel_ij     (k1:k2)
      real*8  :: optdepth_ij(k1:k2)
      real*8  :: ozone_ij   (k1:k2)
      real*8  :: relativeHumidity_ij   (k1:k2)

      real*8  ::  ODAER_ij  (k1:k2,NSADaer*nrh_b) ! Column optical depth for aerosol
      real*8  ::  ODMDUST_ij(k1:k2,NSADdust)      ! Column optical depth for mineral dust

      real*8  :: qjgmi_ij   (k1:chem_mask_khi, num_qjs)

      real*8  :: sflux_ij_gt_1(k1:chem_mask_khi)
      real*8  :: sflux_ij_gt_2(k1:chem_mask_khi)
      real*8  :: sflux_ij_gt_3(k1:chem_mask_khi)

      real*8  :: n2adj(i1:i2, ju1:j2, k1:k2)
      real*8  :: o2adj(i1:i2, ju1:j2, k1:k2)

      real*8 :: RAA_b(4, NP_b), QAA_b(4, NP_b)

      real*8, allocatable :: r_flux_1    (:,:,:)
      real*8, allocatable :: r_flux_2    (:,:,:)
      real*8, allocatable :: r_flux_3    (:,:,:)
      real*8, allocatable :: tau_cloud_n (:,:,:)
 
      logical :: rootProc

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag)  &
         Write (6,*) 'calcPhotolysisRateConstants called by ', loc_proc

#ifdef GTmodule
      photOption = 3
      Allocate (r_flux_1 (i1:i2, ju1:j2,k1:k2))
      r_flux_2=0.0d0
      Allocate (r_flux_2 (i1:i2, ju1:j2,k1:k2))
      r_flux_2=0.0d0
      Allocate (r_flux_3 (i1:i2, ju1:j2,k1:k2))
      r_flux_3=0.0d0
      Allocate (tau_cloud_n (i1:i2, ju1:j2,k1:k2))
      tau_cloud_n=0.0d0
#else
      photOption = phot_opt
#endif
!      if (first) then
!         first = .false.
         if (photOption == 3) then
            if ((TRIM(chem_mecha) == 'troposphere') .or. &
                (TRIM(chem_mecha) == 'strat_trop')  .or. &
     &           TRIM(chem_mecha) == 'strat_trop_aerosol' ) then
               if (do_AerDust_Calc) then
                  if (fastj_opt == 0) then
                     call  GetQAA_RAAinFastj (RAA_b, QAA_b)
                  end if
                  if (fastj_opt == 1) then
                     call  GetQAA_RAAinFastJX (JXbundle, RAA_b, QAA_b)
                  end if
                  if (fastj_opt == 2) then
                     call  GetQAA_RAAinFastJX53b (RAA_b, QAA_b)
                  end if
                  if (fastj_opt == 3) then
                     call  GetQAA_RAAinFastJX53c (RAA_b, QAA_b)
                  end if
                  if (fastj_opt == 4) then
                     call  GetQAA_RAAinFastJX65 (RAA_b, QAA_b)
                  end if
               end if
            end if
         end if

!      end if

!     ==================
      if (photOption == 2) then
!     ==================

        if (qj_timpyr == MONTHS_PER_YEAR) then
          call GmiSplitDateTime (nymd, idumyear, month_gmi, idumday)
          it = month_gmi
        else
          it = 1
        endif 

          qjgmi(i1:i2,ju1:j2,k1:k2,1:num_qjs) =  & 
     &                          qjmon(i1:i2,ju1:j2,k1:k2,1:num_qjs,it)

!     ==================
      elseif (photOption == 3) then
!     ==================

!       --------------------------------------------------------------
!       First form some non-grid dependent arguments needed for Fastj.
!       --------------------------------------------------------------

        call GmiSplitDateTime (nymd, idumyear, month_gmi, idumday)

        call GetDaysFromJanuary1 (jday, nymd)

        time_sec = ConvertTimeToSeconds (nhms)

#ifdef GTmodule
        tau_cloud_n(:,:,:)= cloud_tau(:,:,:)

        ! For aerosol/dust optical depth/surface area calculations
        if ((TRIM(chem_mecha) == 'troposphere') .or. &
            (TRIM(chem_mecha) == 'strat_trop' ) .or. &
            (TRIM(chem_mecha) == 'strat_trop_aerosol' ) .or. &
            (TRIM(chem_mecha) == 'aerosol'    )) then
#else
        ! For aerosol/dust optical depth/surface area calculations
        if ((TRIM(chem_mecha) == 'troposphere') .or. &
            (TRIM(chem_mecha) == 'strat_trop' ) .or. &
            (TRIM(chem_mecha) == 'strat_trop_aerosol')) then
#endif
           if (do_AerDust_Calc) then
              call Aero_OptDep_SurfArea ( gridBoxHeight, concentration, &
     &                 OptDepth, Eradius, Tarea, Odaer, relativeHumidity, Daersl, Waersl, &
     &                 RAA_b, QAA_b, isynoz_num, synoz_threshold, AerDust_Effect_opt, &
     &                 i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species, num_AerDust)
              call Dust_OptDep_SurfArea ( gridBoxHeight, &
     &                 concentration, OptDepth, Eradius, Tarea, Odmdust, Dust, &
     &                 RAA_b, QAA_b, isynoz_num, synoz_threshold, AerDust_Effect_opt, &
     &                 i1, i2, ju1, j2, k1, k2, num_species, num_AerDust)
           end if
        end if

#ifdef GTmodule
        if (pr_qj_opt_depth) qjgmi(:,:,:,num_qjo) = tau_cloud_n(:,:,:)
#else
        if (pr_qj_opt_depth) qjgmi(:,:,:,num_qjo) = tau_cloud(:,:,:)
#endif

!       ------------------------------------------------------------------
!       Now loop over all latitudes and longitudes for this processor
!       because Fastj is set up to be a column calculation.
!       When doing OpenMP this would be a natural place to split the work.
!       ------------------------------------------------------------------

        do ij = ju1, j2
          do il = i1, i2

            kel_ij(:) = temp3(il,ij,:)
            relativeHumidity_ij(:) = relativeHumidity(il,ij,:)

            optdepth_ij(:) = tau_cloud(il,ij,:)

            ODAER_ij   = 0.0d0
            ODMDUST_ij = 0.0d0

#ifdef GTmodule
            if ((TRIM(chem_mecha) == 'troposphere') .or. &
                (TRIM(chem_mecha) == 'strat_trop' ) .or. &
                (TRIM(chem_mecha) == 'strat_trop_aerosol' ) .or. &
                (TRIM(chem_mecha) == 'aerosol'    )) then

               latyp_ij= lwi_flags(il,ij)
#else
            if ((TRIM(chem_mecha) == 'troposphere') .or. &
                (TRIM(chem_mecha) == 'strat_trop' ) .or. &
                (TRIM(chem_mecha) == 'strat_trop_aerosol')) then
#endif
               if (do_AerDust_calc) then
                  ODAER_ij  (:,:) = ODAER  (il,ij,:,:)
                  ODMDUST_ij(:,:) = ODMDUST(il,ij,:,:)
               end if
            end if

            ozone_ij(:) = concentration(io3_num)%pArray3D(il,ij,:)

            if (fastj_opt == 0) then

!                 ==================
                  call Control_Fastj  &
!                 ==================
     &              (k1, k2, chem_mask_khi,  &
     &               num_qjs, month_gmi, jday, time_sec, fastj_offset_sec,  &
     &               londeg(il), latdeg(ij), pres3e(il,ij,k1-1:k2),  &
     &               kel_ij, optdepth_ij, surf_alb_uv(il,ij), qjgmi_ij, overheadO3col_ij, &
     &               ODAER_ij, ODMDUST_ij)

            elseif (fastj_opt == 1) then

               if (.not. do_ozone_inFastJX) then
!                 ==================
                  call Control_Fast_JX  &
!                 ==================
     &              (JXbundle, k1, k2, chem_mask_khi, num_qjs, month_gmi, jday, &
                     time_sec, fastj_offset_sec, londeg(il), latdeg(ij),&
     &               solarZenithAngle(il,ij),  pres3c(il,ij,k1:k2),     &
     &               pctm2(il,ij), kel_ij, optdepth_ij, surf_alb_uv(il,ij), &
     &               qjgmi_ij, overheadO3col_ij, ODAER_ij, ODMDUST_ij, &
!     &               sflux_ij_gt_1, sflux_ij_gt_2, sflux_ij_gt_3, rootProc,latyp_ij,&
     &               sflux_ij_gt_1, sflux_ij_gt_2, sflux_ij_gt_3, loc_proc,latyp_ij,&
     &               ozone_ij)
               else
!                 ==================
                  call Control_Fast_JX  &
!                 ==================
     &              (JXbundle, k1, k2, chem_mask_khi, num_qjs, month_gmi, jday, &
                     time_sec, fastj_offset_sec, londeg(il), latdeg(ij),&
     &               solarZenithAngle(il,ij),  pres3c(il,ij,k1:k2),     &
     &               pctm2(il,ij), kel_ij, optdepth_ij, surf_alb_uv(il,ij), &
     &               qjgmi_ij, overheadO3col_ij, ODAER_ij, ODMDUST_ij, &
     &               sflux_ij_gt_1, sflux_ij_gt_2, sflux_ij_gt_3, loc_proc,latyp_ij)
!     &               sflux_ij_gt_1, sflux_ij_gt_2, sflux_ij_gt_3, rootProc,latyp_ij)
               endif
            elseif (fastj_opt == 2) then
!                    ==================
                     call Control_Fast_JX53b  &
!                    ==================
     &                 (k1, k2, chem_mask_khi,  &
     &                  num_qjs, month_gmi, jday, time_sec,  &
     &                  londeg(il), latdeg(ij), pres3c(il,ij,k1:k2), pctm2(il,ij),  &
     &                  kel_ij, optdepth_ij, surf_alb_uv(il,ij), qjgmi_ij, overheadO3col_ij, &
     &                  ozone_ij)
            elseif (fastj_opt == 3) then
!                    ==================
                     call Control_Fast_JX53c  &
!                    ==================
     &                 (k1, k2, chem_mask_khi, num_qjs, month_gmi, jday, time_sec,  &
     &                  londeg(il), latdeg(ij), pres3c(il,ij,k1:k2), pctm2(il,ij),  &
     &                  kel_ij, optdepth_ij, surf_alb_uv(il,ij), qjgmi_ij, overheadO3col_ij, &
     &                   ODAER_ij, ODMDUST_ij, ozone_ij)
!!                    ==================
!                     call RunFastJX53c  &
!!                    ==================
!     &                 (k1, k2, jday, time_sec, month_gmi,  &
!     &                  londeg(il), latdeg(ij), pres3c(il,ij,k1:k2), pctm2(il,ij),  &
!     &                  kel_ij, optdepth_ij, surf_alb_uv(il,ij), qjgmi_ij,  &
!     &                  ozone_ij)
            elseif (fastj_opt == 4) then
               if (.not. do_ozone_inFastJX) then
                  call controlFastJX65 (k1, k2, chem_mask_khi, num_qjs, month_gmi, &
     &                        jday, time_sec, londeg(il), latdeg(ij), &
     &                        solarZenithAngle(il,ij),  &
     &                        tauclw(il,ij,k1:k2), taucli(il,ij,k1:k2),    &
     &                        pres3c(il,ij,k1:k2), pctm2(il,ij), kel_ij, optdepth_ij, &
     &                        surf_alb_uv(il,ij), qjgmi_ij, relativeHumidity(il,ij,k1:k2),  &
!     &                        rootProc, overheadO3col_ij, ODAER_ij, ODMDUST_ij, &
     &                        loc_proc, overheadO3col_ij, ODAER_ij, ODMDUST_ij, &
     &                         ozone_ij)
               else
                  call controlFastJX65 (k1, k2, chem_mask_khi, num_qjs, month_gmi, &
     &                        jday, time_sec, londeg(il), latdeg(ij), &
     &                        solarZenithAngle(il,ij),  &
     &                        tauclw(il,ij,k1:k2), taucli(il,ij,k1:k2),    &
     &                        pres3c(il,ij,k1:k2), pctm2(il,ij), kel_ij, optdepth_ij, &
     &                        surf_alb_uv(il,ij), qjgmi_ij, relativeHumidity(il,ij,k1:k2),  &
     &                        loc_proc, overheadO3col_ij, ODAER_ij, ODMDUST_ij)
!     &                        rootProc,overheadO3col_ij, ODAER_ij, ODMDUST_ij)
               endif
            endif

            qjgmi(il,ij,k1:chem_mask_khi,1:num_qjs) =  &
     &        qjgmi_ij(k1:chem_mask_khi,1:num_qjs)

!... adjust the J(NO) rate with value read in from namelist variable
            if(jno_num.ne.999) then           
               qjgmi(il,ij,k1:chem_mask_khi,jno_num) =  &
     &            qjgmi(il,ij,k1:chem_mask_khi,jno_num) * jno_adjust
            endif

            overheadO3col(il,ij,:) = overheadO3col_ij(:)

#ifdef GTmodule
            r_flux_1(il,ij,k1:chem_mask_khi)=sflux_ij_gt_1(k1:chem_mask_khi) !Clear sky solar flux W/m2)
            r_flux_2(il,ij,k1:chem_mask_khi)=sflux_ij_gt_2(k1:chem_mask_khi) !nika All-sky shortwave W/m2)
#endif
          end do
        end do

#ifdef GTmodule
          flux_gt(:,:,:)=  r_flux_2(:,:,:) - r_flux_1(:,:,:)   ! All_sky - Clear_sky

!          cloud_param(:,:,:,1)=   flux_gt(:,:,:)
!          cloud_param(:,:,:,3) = r_flux_2(:,:,:)    ! clear sky
!          cloud_param(:,:,1,3) = r_flux_1(:,:,1)    ! All-sky
! start rsot
          cloud_param(:,:,:,22)=   flux_gt(:,:,:)
          cloud_param(:,:,:,23) = r_flux_2(:,:,:)    ! clear sky
          cloud_param(:,:,1,23) = r_flux_1(:,:,1)    ! All-sky
! end rsot
#endif

!     ==============================================
      else if ((photOption == 4) .or. (photOption == 5)) then
!     ==============================================

!       ==============
        call Lookup_Qj  &
!       ==============
     &    (do_clear_sky, photOption, io3_num, nymd, photintv,  &
     &     rsec_jan1, pres3c, temp3, concentration, latdeg, mcor, mass3,  &
     &     fracCloudCover, qjgmi, pr_diag, loc_proc, num_species, num_qjs, &
     &     num_qjo, ilo, ihi, julo, jhi, i2_gl, ju1_gl, j2_gl, i1, i2, ju1, &
     &     j2, k1, k2, TRIM(chem_mecha))

      end if

      if (TRIM(chem_mecha) == 'troposphere') then
!        ----------------------------------------------------------------------
!        Check to see if the mechanism has a photolysis reaction O3 + hv = 2OH.
!        If so, save the rate in the last entry of qjgmi and then the real rate
!        needs to be adjusted.  Updated to JPL 06-2 (Bryan Duncan 10/2006).
!        ----------------------------------------------------------------------
                   
         if (num_qj_o3_to_2oh > 0) then
                   
            n2adj(:,:,:) = 2.15d-11 * Exp (110.0d0 / temp3(i1:i2,ju1:j2,:)) * MXRN2
            o2adj(:,:,:) = 3.30d-11 * Exp ( 55.0d0 / temp3(i1:i2,ju1:j2,:)) * MXRO2
                   
            if (pr_qj_o3_o1d) then
               qjgmi(:,:,:,num_qjs+1) = qjgmi(:,:,:,num_qj_o3_to_2oh)
            end if
                   
            qjgmi(:,:,:,num_qj_o3_to_2oh) = qjgmi(:,:,:,num_qj_o3_to_2oh) /  &
     &                         (1.0d0 + ((n2adj(:,:,:) + o2adj(:,:,:)) /     &
     &                         (1.63d-10 * Exp(60.0d0/temp3(i1:i2,ju1:j2,:)) &
     &                        * concentration(ih2o_num)%pArray3D(:,:,:))))
                   
         end if
      end if

#ifdef GTmodule
     Deallocate (r_flux_1)
     Deallocate (r_flux_2)
     Deallocate (r_flux_3)
#endif

      return

      end subroutine calcPhotolysisRateConstants

      end module GmiPhotolysisRateConstants_mod
