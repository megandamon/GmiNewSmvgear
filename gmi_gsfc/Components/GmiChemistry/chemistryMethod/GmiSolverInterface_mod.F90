!------------------------------------------------------------------------------
!BOP
! !MODULE: GmiSolverInterface_mod
!
      module GmiSolverInterface_mod
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiFileOperations_mod    , only : makeOutfileName
      use Ftiming_Dao              , only : Ftiming_On, Ftiming_Off
      use GmiTimeControl_mod       , only : GmiSplitDateTime
      use GmiTimeControl_mod       , only : GetSecondsFromJanuary1
      use GmiPrintError_mod        , only : GmiPrintError
      use GmiSurfaceEmissionInChemistry_mod, only : updateSurfEmissionInChemistry
      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved
!
      implicit none
!
      private
!
! !PUBLIC MEMBER FUNCTIONS:
      public  :: Update_Smv2chem, Update_Sulfchem

#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
!
!EOP
!------------------------------------------------------------------------------
      contains
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update_Smv2chem
!
! !INTERFACE:
!
      subroutine Update_Smv2chem (savedVars, chemintv, emiss_isop, emiss_monot, emiss_nox,&
     &                  do_ShipEmission, emiss_hno3, emiss_o3, ihno3_num,      &
     &                  io3_num, mcor, surf_emiss_out, surf_emiss_out2,        &
     &                  emiss_3d_out, humidity, qjgmi, qkgmi, emissionArray,   &
     &                  press3e, pctm2, kel, concentration, pr_diag, pr_qqjk,  &
     &                  pr_smv2, pr_surf_emiss, pr_emiss_3d, do_smv_reord,     &
     &                  do_synoz, do_qqjk_inchem, do_semiss_inchem, ico_num,   &
     &                  ino_num, ipropene_num, iisoprene_num, imgas_num,       &
     &                  initrogen_num, ioxygen_num, isynoz_num, yda, qqkda,    &
     &                  qqjda, mw, pr_nc_period, emiss_timpyr, emiss_opt,      &
     &                  emiss_map, tdt, nymd, chem_mask_klo, chem_mask_khi,    &
     &                  loc_proc, synoz_threshold, ilong, ilat, ivert, itloop, &
     &                  i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi,          &
     &                  num_molefrac, num_emiss, num_qjo, num_qks, num_qjs,    &
     &                  num_active, num_species, commuWorld, chem_mecha)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: ico_num, ino_num, ipropene_num, iisoprene_num
      integer, intent(in) :: imgas_num, initrogen_num, ioxygen_num, isynoz_num
      integer, intent(in) :: ilong, ilat, ivert, itloop, ihno3_num, io3_num
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in) :: chem_mask_klo, chem_mask_khi
      integer, intent(in) :: num_molefrac
      integer, intent(in) :: num_emiss, num_qjo, num_qks, num_qjs
      integer, intent(in) :: num_active, num_species
      integer, intent(in) :: commuWorld
      integer, intent(in) :: loc_proc, nymd, emiss_opt, emiss_timpyr
      integer, intent(in) :: emiss_map(num_emiss)
      real*8,  intent(in) :: mw(num_species)
      real*8,  intent(in) :: pr_nc_period
      real*8,  intent(in) :: synoz_threshold
      real*8,  intent(in) :: tdt
      logical, intent(in) :: pr_diag, do_ShipEmission
      logical, intent(in) :: pr_qqjk
      logical, intent(in) :: pr_smv2
      logical, intent(in) :: do_smv_reord
      logical, intent(in) :: do_synoz
      logical, intent(in) :: do_qqjk_inchem
      logical, intent(in) :: do_semiss_inchem
      logical, intent(in) :: pr_surf_emiss, pr_emiss_3d
      character(len=*), intent(in) :: chem_mecha
      ! chemistry time step (s)
      real*8  :: chemintv
      ! isoprene    emissions (kg/s)
      real*8  :: emiss_isop (i1:i2, ju1:j2)
      ! monoterpene emissions (kg/s)
      real*8  :: emiss_monot(i1:i2, ju1:j2)
      ! NOx         emissions (kg/s)
      real*8  :: emiss_nox  (i1:i2, ju1:j2)
      real*8  :: emiss_hno3 (i1:i2, ju1:j2)
      real*8  :: emiss_o3(i1:i2, ju1:j2)
      ! area of grid box  (m^2)
      real*8  :: mcor       (i1:i2, ju1:j2)
      ! accumulation of surface emission for output (kg/m^2/time)
      real*8  :: surf_emiss_out(i1:i2, ju1:j2, num_species)
      real*8  :: surf_emiss_out2(i1:i2, ju1:j2, 6)
      real*8  :: emiss_3d_out(i1:i2, ju1:j2, k1:k2, num_species)
      ! specific humidity (g/kg)
      real*8  :: humidity(i1:i2, ju1:j2, k1:k2)
      ! photolysis rate constants (s^-1)
      real*8  :: qjgmi(i1:i2, ju1:j2, k1:k2, num_qjo)
      ! thermal rate constants (units vary)
      real*8  :: qkgmi(i1:i2, ju1:j2, k1:k2, num_qks)
      ! CTM surface pressure at t1+tdt (mb)
      real*8  :: pctm2(ilo:ihi, julo:jhi)
      real*8  :: press3e(ilo:ihi, julo:jhi, k1-1:k2)
      ! temperature (degK)
      real*8  :: kel  (ilo:ihi, julo:jhi, k1:k2)
!
! !INPUT/OUPUT PARAMETERS:
      real*8 , intent(inOut) :: qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)
      real*8 , intent(inOut) :: qqkda(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inOut) :: yda  (i1:i2, ju1:j2, k1:k2, num_active)
      ! Species concentration, known at zone centers (mixing ratio)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
      ! Array of emissions (kg/s)
      type (t_GmiArrayBundle), intent(in) :: emissionArray(num_emiss  )
      type (t_Smv2Saved), intent(inOut) :: savedVars
!
! !DESCRIPTION:
!   Updates the Smvgear2 chemistry.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_FILE_NAME) :: smv_filnam
      logical :: do_cell_chem(i1:i2, ju1:j2, k1:k2)
      integer :: ic
      real*8  :: surf_emiss(i1:i2, ju1:j2, num_species)
      real*8  :: surf_emiss2(i1:i2, ju1:j2, 6)
      real*8  :: emiss_3d(i1:i2, ju1:j2, k1:k2, num_species)
      real*8, allocatable  :: speciesConc(:,:,:,:)
!
! !REVISION HISTORY:
!   Initial code.
!EOP
!-------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Update_Smv2chem called by ', loc_proc

!     ------------------------------------------------
!     Change units from mixing ratio to concentration.
!     ------------------------------------------------

      do ic = 1, num_molefrac

        concentration(ic)%pArray3D(:,:,:) =  &
     &    concentration(ic)%pArray3D(:,:,:) * concentration(imgas_num)%pArray3D(:,:,:)

      end do


      if (initrogen_num /= 0) then
        concentration(initrogen_num)%pArray3D(:,:,:) =  &
     &    MXRN2 * concentration(imgas_num)%pArray3D(:,:,:)
      end if

      if (ioxygen_num   /= 0) then
        concentration(ioxygen_num)%pArray3D(:,:,:)   =  &
     &    MXRO2 * concentration(imgas_num)%pArray3D(:,:,:)
      end if


      if (do_semiss_inchem) then

!       =========================
        call updateSurfEmissionInChemistry  &
!       =========================
     &    (pr_surf_emiss, pr_emiss_3d, emiss_isop, emiss_monot, emiss_nox,  &
     &     do_ShipEmission, emiss_hno3, emiss_o3, ihno3_num, io3_num, &
     &     mcor, surf_emiss_out, surf_emiss_out2, emiss_3d_out, humidity,  &
     &     press3e, pctm2, kel, emissionArray, surf_emiss, surf_emiss2, emiss_3d, &
     &     emiss_timpyr, num_emiss, emiss_opt, emiss_map, tdt, nymd, &
     &   ico_num, ino_num, ipropene_num, iisoprene_num, mw, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species)

      end if

      do_cell_chem(:,:,:) = .false.

      do_cell_chem(:,:,chem_mask_klo:chem_mask_khi) = .true.

      IF (chem_mecha /= 'strat_trop'         .and. &
          chem_mecha /= 'strat_trop_aerosol' .and. &
          chem_mecha /= 'tracers' ) THEN
         if (do_synoz) then
           where (concentration(isynoz_num)%pArray3D(:,:,:) > synoz_threshold)
             do_cell_chem(:,:,:) = .false.
           end where

         end if
      ENDIF

       allocate(speciesConc(i1:i2, ju1:j2, k1:k2, num_species))
       do ic = 1, num_species
          speciesConc(:,:,:,ic) = concentration(ic)%pArray3D(:,:,:)
       end do

!c?   OK with chem sub-cycling?  tdt?
!     ===================
      call Do_Smv2_Solver  &
!     ===================
     &  (savedVars, do_qqjk_inchem, do_semiss_inchem, pr_diag, pr_qqjk, pr_smv2,  &
     &   loc_proc, ilat, ilong, ivert, itloop, pr_nc_period, tdt,  &
     &   do_cell_chem, qkgmi, qjgmi, surf_emiss, speciesConc, &
     &   yda, qqkda, qqjda, &
     &   i1, i2, ju1, j2, k1, k2, &
     &   num_qks, num_qjs, num_active, commuWorld)

       do ic = 1, num_species
          concentration(ic)%pArray3D(:,:,:) = speciesConc(:,:,:,ic)
       end do
       deallocate(speciesConc)

!     -----------------------------------------------------
!     Change units from concentration back to mixing ratio.
!     -----------------------------------------------------

      do ic = 1, num_molefrac

        concentration(ic)%pArray3D(:,:,:) =  &
     &    concentration(ic)%pArray3D(:,:,:) / concentration(imgas_num)%pArray3D(:,:,:)

      end do


      return

      end subroutine Update_Smv2chem
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update_Sulfchem
!
! !INTERFACE:
!
      subroutine Update_Sulfchem (cloud_param, chemintv, londeg, latdeg, kel, &
     &                  concentration, humidity, emissionArray, mcor, qjgmi,  &
     &                  qkgmi, moistq, cmf1, press3c, gridBoxHeight, mass,    &
     &                  dms_oh, dms_no3, so2_oh, so2_h2o2, so2_o3, emiss_dust,&
     &                  aqua_infile_name, qj_infile_name, qj_var_name,        &
     &                  pr_diag, pr_sulf_src, do_aerocom,                     &
     &                  do_clear_sky, imgas_num, emiss_dust_opt, emiss_map,   &
     &                  emiss_map_dust, mw, phot_opt, nymd, nhms,             &
     &                  num_time_steps, loc_proc, ndust, ilong, itloop,       &
     &                  i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, i1_gl,  &
     &                  i2_gl, ju1_gl, j2_gl, emiss_timpyr, num_molefrac,     &
     &                  num_emiss, num_qjs, num_qks, num_species, chem_opt)
!

      implicit none

#     include "setkin_par.h"
#     include "setkin_depos.h"
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      logical, intent(in) :: pr_sulf_src
      logical, intent(in) :: do_aerocom
      logical, intent(in) :: do_clear_sky
      integer, intent(in) :: loc_proc
      integer, intent(in) :: ilong, itloop
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: num_molefrac, num_emiss, num_qjs, num_qks, num_species
      integer, intent(in) :: ndust, emiss_timpyr, phot_opt
      integer, intent(in) :: nymd, nhms
      integer, intent(in) :: num_time_steps
      integer, intent(in) :: imgas_num
      integer, intent(in) :: emiss_dust_opt
      integer, intent(in) :: emiss_map(num_emiss)
      integer, intent(in) :: emiss_map_dust(ndust)
      real*8 , intent(in) :: mw(num_species)
      integer, intent(in) :: chem_opt
      ! chemistry time step (s)
      real*8  :: chemintv
      ! longitude (deg)
      real*8 , intent(in) :: londeg  (i1_gl:i2_gl)
      ! latitude (deg)
      real*8 , intent(in) :: latdeg  (ju1_gl:j2_gl)
      ! temperature (degK)
      real*8  :: kel     (ilo:ihi, julo:jhi, k1:k2)
      ! specific humidity  (g/kg)
      real*8 , intent(in) :: humidity(i1:i2, ju1:j2, k1:k2)
      ! area of grid box   (m^2)
      real*8 , intent(in) :: mcor    (i1:i2, ju1:j2)
      ! photolysis rate constants (s^-1)
      real*8  :: qjgmi   (i1:i2, ju1:j2, k1:k2, num_qjs)
      ! thermal    rate constants (units vary)
      real*8  :: qkgmi   (i1:i2, ju1:j2, k1:k2, num_qks)
      ! moisture changes due to wet processes (g/kg/day)
      real*8  :: moistq  (i1:i2, ju1:j2, k1:k2)
      ! convective mass flux used for cloud fraction calculation (kg/m^2/s)
      real*8  :: cmf1    (i1:i2,  ju1:j2, k1:k2)
      ! atmospheric pressure at the center of each grid box (mb)
      real*8 , intent(in) :: press3c (ilo:ihi, julo:jhi, k1:k2)
      ! height of each grid box (m)
      real*8 , intent(in) :: gridBoxHeight(i1:i2, ju1:j2, k1:k2)
      ! total mass of the atmosphere within each grid box (kg)
      real*8 , intent(in) :: mass    (i1:i2, ju1:j2, k1:k2)
      ! SO2 production from DMS + OH  accumulated since last output used for
      ! sulfur budget (kg)
      real*8  :: dms_oh  (i1:i2, ju1:j2, k1:k2)
      ! SO2 production from DMS + NO3 accumulated since last output used for
      ! sulfur budget (kg)
      real*8  :: dms_no3 (i1:i2, ju1:j2, k1:k2)
      ! SO4 production from SO2 + OH  accumulated since last output used for
      ! sulfur budget (kg)
      real*8  :: so2_oh  (i1:i2, ju1:j2, k1:k2)
      ! SO4 production from SO2 + H2O2 (aqueous phase) accumulated since last
      ! utput used for sulfur budget (kg)
      real*8  :: so2_h2o2(i1:i2, ju1:j2, k1:k2)
      ! SO4 production from SO2 + O3   (aqueous phase) accumulated since last
      ! output used for sulfur budget (kg)
      real*8  :: so2_o3  (i1:i2, ju1:j2, k1:k2)
      real*8  :: emiss_dust(i1:i2, ju1:j2, ndust)
      character (len=MAX_LENGTH_VAR_NAME), intent(in) :: qj_var_name
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: qj_infile_name
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: aqua_infile_name
!
! !INPUT/OUTPUT PARAMETERS:
      ! species concentration, known at zone centers (mixing ratio)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
      ! array of emissions (kg/s)
      type (t_GmiArrayBundle), intent(in) :: emissionArray(num_emiss  )
      real*8, intent(inOut) :: cloud_param(i1:i2, ju1:j2, k1:k2,num_species)
!
! !DESCRIPTION:
!   Updates the sulfur chemistry.
!
! !DEFINED PARAMETERS:
      real*8, parameter :: CMPM3 = CMPM * CMPM * CMPM
      real*8, parameter :: M2 = 1.0d12
      real*8, parameter :: WTAIR = 28.9d0
!
! !LOCAL VARIABLES:
      integer :: il, ij, ik, ic, icx
      integer :: idumday, idumyear, imon, month
      real*8  :: fac
      real*8  :: conv_emiss(i1:i2, ju1:j2, k1:k2)
      real*8  :: sulf_emiss(i1:i2, ju1:j2, k1:k2, num_species)
      real*8, allocatable :: c0 (:)
      real*8, allocatable :: c1 (:)
      real*8, allocatable :: c2 (:)
      real*8, allocatable :: cfc(:)
      real*8, allocatable :: cfs(:)
      real*8, allocatable :: rhc(:)
      ! convective mass flux   (mb/hr)
      real*8, allocatable :: cmfcol    (:)
      ! height above sea level (m)
      real*8, allocatable :: height_asl(:)
      ! total cloud fraction [0 - 1]
      real*8, allocatable :: clfrac (:,:,:)
      ! in-cloud liquid water in each grid box (g/g)
      real*8, allocatable :: lwccol (:,:,:)
      ! relative humidity [0 - 1]
      real*8, allocatable :: relhum (:,:,:)

      real*8, allocatable :: settle_so4n (:,:,:)
      real*8  :: tau, lambda

      ! total sulfate (g/g)
      real*8, allocatable :: totso4 (:,:,:)
      real*8, allocatable  :: speciesConc(:,:,:,:)
!
! !REVISION HISTORY:
!   Initial code.
!EOP
!-------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'Update_Sulfchem called by ', loc_proc

      Allocate (c0 (k1:k2))
      Allocate (c1 (k1:k2))
      Allocate (c2 (k1:k2))
      c0 = 0.0d0; c1 = 0.0d0; c2 = 0.0d0

      Allocate (cfc(k1:k2))
      Allocate (cfs(k1:k2))
      cfc = 0.0d0; cfs = 0.0d0

      Allocate (rhc(k1:k2))
      rhc = 0.0d0

      Allocate (cmfcol    (k1:k2))
      Allocate (height_asl(k1:k2))
      cmfcol = 0.0d0; height_asl = 0.0d0

      Allocate (clfrac (i1:i2, ju1:j2, k1:k2))
      clfrac = 0.0d0

      Allocate (lwccol (i1:i2, ju1:j2, k1:k2))
      Allocate (relhum (i1:i2, ju1:j2, k1:k2))
      lwccol = 0.0d0; relhum = 0.0d0

#ifdef GTmodule
      Allocate (totso4 (i1:i2, ju1:j2, k1:k2))
      totso4 = 0.0d0
#endif


!     -----------------------------------------------------------
!     Calculate relative humidity from Seinfeld (1986) p. 181.
!     The first rh is the temperature dependent parameter a.
!     The second rh is the saturation vapor pressure of water.
!     The third rh is the actual relative humidity as a fraction.
!     Then make sure rh is between 0 and 0.95.
!     -----------------------------------------------------------

      do ij = ju1, j2
        do il = i1, i2

          relhum(il,ij,:) =  &
     &      1.0d0 - (373.15d0 / kel(il,ij,:))

          relhum(il,ij,:) =  &
     &      1013.25d0 * Exp (13.3185d0 * relhum(il,ij,:)    -  &
     &                        1.9760d0 * relhum(il,ij,:)**2 -  &
     &                        0.6445d0 * relhum(il,ij,:)**3 -  &
     &                        0.1299d0 * relhum(il,ij,:)**4)

          relhum(il,ij,:) =  &
     &      humidity(il,ij,:) * MWTAIR / 18.0d0 / GPKG *  &
     &      press3c(il,ij,:) / relhum(il,ij,:)

          relhum(il,ij,:) =  &
     &      Max (Min (relhum(il,ij,:), 0.95d0), 0.0d0)

        end do
      end do


!     ======================
      if (do_clear_sky) then
!     ======================

        clfrac(:,:,:) = 0.0d0
        lwccol(:,:,:) = 0.0d0

!     ====
      else
!     ====

!       -------------------------------------------------------------------
!       Calculate cloud water from a parameterization of Kiehl found in
!       "Sensitivity of the Simulated Climate to a Diagnostic Formulation
!       for Cloud Liquid Water" by James Hack, Journal of Climate, Vol. 11,
!       July 1998, p 1499.
!       -------------------------------------------------------------------

        do ij = ju1, j2
          do il = i1, i2

            height_asl(:) =  &
     &        -29.3d0 * kel(il,ij,:) *  &
     &        Log (press3c(il,ij,:) / 1013.25d0)

            height_asl(:) = Max (height_asl(:), 0.0d0)


            lwccol(il,ij,:) =  &
     &        0.18d0 *  &
     &        Exp (-height_asl(:) /  &
     &             (1080.0d0 +  &
     &              (2000.0d0 * Cos (latdeg(ij) * RADPDEG) ** 2.0d0)))


!           ---------------------------------------------------------------
!           Calculate stratiform cloud fractions from a parameterization of
!           Sundqvist et al. (1988).
!           ---------------------------------------------------------------

!           -----------------------------------------------------
!           Now set up the constants to be used in the stratiform
!           parameterization; from Xu and Krueger.
!           -----------------------------------------------------

            where ((height_asl(:) < 15000.0d0) .and.  &
     &             (height_asl(:) >  5900.0d0))
              rhc(:) =  0.445d0
            end where

            where ((height_asl(:) <= 5900.0d0) .and.  &
     &             (height_asl(:) >  2350.0d0))
              rhc(:) =  0.715d0
            end where

            where  (height_asl(:) <= 2350.0d0)
              rhc(:) =  0.780d0
            end where

            where (relhum(il,ij,:) > rhc(:) .and.  &
     &             height_asl(:)   < 15000.0d0)

              cfs(:) =  &
     &          1.0d0 -  &
     &          Sqrt (1.0d0 - (relhum(il,ij,:) - rhc(:)) /  &
     &                         (1.0d0 - rhc(:)))
              cfs(:) = Max (cfs(:), 0.0d0)

            elsewhere

              cfs(:) = 0.0d0

            end where


!           ------------------------------------------------------------
!           Calculate convective cloud fractions from a parameterization
!           of Xu and Krueger, Monthly Weather Review, Vol. 119 p 342.
!           ------------------------------------------------------------

!           --------------------------------------------------
!           Now set up the constants to be used in the cumulus
!           parameterization.
!           --------------------------------------------------

            where ((height_asl(:) < 15000.0d0) .and.  &
     &             (height_asl(:) >  5900.0d0))
              c0(:) = 0.0131d0
              c1(:) = 0.0123d0
              c2(:) = 0.0030d0
            end where

            where ((height_asl(:) <= 5900.0d0) .and.  &
     &             (height_asl(:) >  2350.0d0))
              c0(:) = 0.0722d0
              c1(:) = 0.0760d0
              c2(:) = 0.0254d0
            end where

            where  (height_asl(:) <= 2350.0d0)
              c0(:) = 0.0337d0
              c1(:) = 0.0453d0
              c2(:) = 0.0237d0
            end where


!           -------------------------------------
!           Calculate the cumulus cloud fraction.
!           -------------------------------------

            cmfcol(:) = cmf1(il,ij,:) * 3600.0d0 * 9.8d0 / 100.0d0

            where (cmfcol(:) > 0.01d0)

              cfc(:) = c0(:) +  &
     &                 c1(:) *  Log10 (cmfcol(:)) +  &
     &                 c2(:) * (Log10 (cmfcol(:)))**2
              cfc(:) = Max (cfc(:), 0.0d0)

            elsewhere

              cfc(:) = 0.0d0

            end where


!           ---------------------------------------------------------
!           Compute total cloud fraction for combined large-scale and
!           convective clouds (GEOS-1 GCM descriptive documentation).
!           ---------------------------------------------------------

            clfrac(il,ij,:) =  &
     &        1.0d0 - (1.0d0 - cfc(:)) * (1.0d0 - cfs(:))

          end do
        end do


!       -----------------------------------------
!       In-cloud liquid water mixing ratio (g/g).
!       -----------------------------------------

        do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2

              if (moistq(il,ij,ik) < -1.0d-5) then

                if (clfrac(il,ij,ik) > 1.0d-2) then

                  lwccol(il,ij,ik) =  &
     &              lwccol(il,ij,ik) * 1.0d-6 /  &
     &              (concentration(imgas_num)%pArray3D(il,ij,ik) * MWTAIR / AVOGAD)

!c                -----------------------------------------------------------
!c                Do not allow cloud in the lowest level according to NCAR
!c                CCM3 because LWC is prescribed as the maximum in the lowest
!c                level, Xhliu, 01/17/03.
!c                -----------------------------------------------------------

                  if (ik == k1) lwccol(il,ij,ik) = 0.0d0

!c                ------------------------------------------------------------
!c                Do not allow cloud liquid water when temperature below -25C.
!c                ------------------------------------------------------------

                  if (kel(il,ij,ik) < 248.0d0) then
                    lwccol(il,ij,ik) = 0.0d0
                  end if

                else

                  lwccol(il,ij,ik) = 0.0d0

                end if

              else

                lwccol(il,ij,ik) = 0.0d0

              end if

            end do
          end do
        end do


!c      -------------------------
!c      Added by Xhliu, 02/05/03.
!c      -------------------------

        where ((clfrac(:,:,:) <= 1.0d-2) .or.  &
     &         (lwccol(:,:,:) <= 1.0d-9))
          clfrac(:,:,:) = 0.0d0
          lwccol(:,:,:) = 0.0d0
        end where

!     ======
      end if
!     ======


!     ------------------------------------------------
!     Change units from mixing ratio to concentration.
!     (only for gas-species, not for aerosols).
!     ------------------------------------------------

      do ic = 1, num_molefrac

        if (aerosol(ic) <= 0) then
          concentration(ic)%pArray3D(i1:i2,ju1:j2,:) =  &
     &      concentration(ic)%pArray3D(i1:i2,ju1:j2,:) *  &
     &      concentration(imgas_num)%pArray3D(i1:i2,ju1:j2,:)
        end if

      end do


      if (emiss_timpyr == MONTHS_PER_YEAR) then
        call GmiSplitDateTime (nymd, idumyear, month, idumday)
        imon = month
      else
        imon = 1
      end if


!     --------------------------------------------------------------------
!     Calculate conversion factor to go from kg/box/s to molecules/cm^3/s,
!     but leave out the molecular weight term for each species.
!     --------------------------------------------------------------------

      fac = AVOGAD * GPKG

      do ik = k1, k2
        conv_emiss(:,:,ik) =  &
     &    fac / mcor(:,:) / CMPM3 / gridBoxHeight(:,:,ik)
      end do


!     --------------------------------------------------------------
!     Now do the conversion for each sulfur species that is emitted.
!     --------------------------------------------------------------

      sulf_emiss(:,:,:,:) = 0.0d0

      do icx = 1, num_emiss
        ic = emiss_map(icx)
        if (ic > 0) then

           if ((ic == IFSO2) .or. (ic == INSO2) .or. (ic == INDMS)) then
              sulf_emiss(:,:,:,ic) =  &
     &          emissionArray(icx)%pArray3D(:,:,:) * conv_emiss(:,:,:) / mw(ic)
           endif
        endif
      end do

      if (do_aerocom) then
       if (emiss_dust_opt /= 0) then
        do icx = 1, ndust

          ic = emiss_map_dust(icx)

          if (ic == INDMS) then
            sulf_emiss(:,:,1,ic) =  &
     &        emiss_dust(:,:,icx) * conv_emiss(:,:,1) / mw(ic)
          end if

        end do
       end if
      end if

       allocate(speciesConc(i1:i2, ju1:j2, k1:k2, num_species))
       do ic = 1, num_species
          speciesConc(:,:,:,ic) = concentration(ic)%pArray3D(:,:,:)
       end do

!     ===================
      call Do_Sulf_Solver  &
!     ===================
     &  (chemintv, nymd, nhms, kel(i1:i2,ju1:j2,k1:k2),  &
     &   press3c(i1:i2,ju1:j2,k1:k2), speciesConc, humidity,  &
     &   sulf_emiss, qjgmi, qkgmi, phot_opt, clfrac, lwccol, relhum,  &
     &   aqua_infile_name, num_time_steps, londeg, latdeg, mass,  &
     &   pr_sulf_src, do_aerocom,  &
     &   dms_oh, dms_no3, so2_oh, so2_h2o2, so2_o3, &
     &   pr_diag, loc_proc, ilong, itloop, &
     &   i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   num_qjs, num_qks, num_species, chem_opt)


!.sds... artificial settling for SO4 aerosols - move down one level with a 30-day efold time
       Allocate (settle_so4n (i1:i2, ju1:j2, k1:k2))
       tau = (30.)
       lambda = (1.0/(tau * 86400.))
!... NSO4A
       settle_so4n(:,:,:) = speciesConc(:,:,:,INSO4A) * (1.0d0 - exp(-chemintv * lambda))
       speciesConc(:,:,:,INSO4A) = speciesConc(:,:,:,INSO4A) - settle_so4n(:,:,:)
       do ik = k1+1, k2
         speciesConc(:,:,ik-1,INSO4A) = speciesConc(:,:,ik-1,INSO4A) + settle_so4n(:,:,ik)
       enddo
       speciesConc(:,:,k1,INSO4A) = speciesConc(:,:,k1,INSO4A) + settle_so4n(:,:,k1)
!... NSO4N1
       settle_so4n(:,:,:) = speciesConc(:,:,:,INSO4N1) * (1.0d0 - exp(-chemintv * lambda))
       speciesConc(:,:,:,INSO4N1) = speciesConc(:,:,:,INSO4N1) - settle_so4n(:,:,:)
       do ik = k1+1, k2
         speciesConc(:,:,ik-1,INSO4N1) = speciesConc(:,:,ik-1,INSO4N1) + settle_so4n(:,:,ik)
       enddo
       speciesConc(:,:,k1,INSO4N1) = speciesConc(:,:,k1,INSO4N1) + settle_so4n(:,:,k1)
!... NSO4N2
       settle_so4n(:,:,:) = speciesConc(:,:,:,INSO4N2) * (1.0d0 - exp(-chemintv * lambda))
       speciesConc(:,:,:,INSO4N2) = speciesConc(:,:,:,INSO4N2) - settle_so4n(:,:,:)
       do ik = k1+1, k2
         speciesConc(:,:,ik-1,INSO4N2) = speciesConc(:,:,ik-1,INSO4N2) + settle_so4n(:,:,ik)
       enddo
       speciesConc(:,:,k1,INSO4N2) = speciesConc(:,:,k1,INSO4N2) + settle_so4n(:,:,k1)
!... NSO4N3
       settle_so4n(:,:,:) = speciesConc(:,:,:,INSO4N3) * (1.0d0 - exp(-chemintv * lambda))
       speciesConc(:,:,:,INSO4N3) = speciesConc(:,:,:,INSO4N3) - settle_so4n(:,:,:)
       do ik = k1+1, k2
         speciesConc(:,:,ik-1,INSO4N3) = speciesConc(:,:,ik-1,INSO4N3) + settle_so4n(:,:,ik)
       enddo
       speciesConc(:,:,k1,INSO4N3) = speciesConc(:,:,k1,INSO4N3) + settle_so4n(:,:,k1)
!... FSO4A
       settle_so4n(:,:,:) = speciesConc(:,:,:,IFSO4A) * (1.0d0 - exp(-chemintv * lambda))
       speciesConc(:,:,:,IFSO4A) = speciesConc(:,:,:,IFSO4A) - settle_so4n(:,:,:)
       do ik = k1+1, k2
         speciesConc(:,:,ik-1,IFSO4A) = speciesConc(:,:,ik-1,IFSO4A) + settle_so4n(:,:,ik)
       enddo
       speciesConc(:,:,k1,IFSO4A) = speciesConc(:,:,k1,IFSO4A) + settle_so4n(:,:,k1)
!... FSO4N1
       settle_so4n(:,:,:) = speciesConc(:,:,:,IFSO4N1) * (1.0d0 - exp(-chemintv * lambda))
       speciesConc(:,:,:,IFSO4N1) = speciesConc(:,:,:,IFSO4N1) - settle_so4n(:,:,:)
       do ik = k1+1, k2
         speciesConc(:,:,ik-1,IFSO4N1) = speciesConc(:,:,ik-1,IFSO4N1) + settle_so4n(:,:,ik)
       enddo
       speciesConc(:,:,k1,IFSO4N1) = speciesConc(:,:,k1,IFSO4N1) + settle_so4n(:,:,k1)
!... FSO4N2
       settle_so4n(:,:,:) = speciesConc(:,:,:,IFSO4N2) * (1.0d0 - exp(-chemintv * lambda))
       speciesConc(:,:,:,IFSO4N2) = speciesConc(:,:,:,IFSO4N2) - settle_so4n(:,:,:)
       do ik = k1+1, k2
         speciesConc(:,:,ik-1,IFSO4N2) = speciesConc(:,:,ik-1,IFSO4N2) + settle_so4n(:,:,ik)
       enddo
       speciesConc(:,:,k1,IFSO4N2) = speciesConc(:,:,k1,IFSO4N2) + settle_so4n(:,:,k1)
!... FSO4N3
       settle_so4n(:,:,:) = speciesConc(:,:,:,IFSO4N3) * (1.0d0 - exp(-chemintv * lambda))
       speciesConc(:,:,:,IFSO4N3) = speciesConc(:,:,:,IFSO4N3) - settle_so4n(:,:,:)
       do ik = k1+1, k2
         speciesConc(:,:,ik-1,IFSO4N3) = speciesConc(:,:,ik-1,IFSO4N3) + settle_so4n(:,:,ik)
       enddo
       speciesConc(:,:,k1,IFSO4N3) = speciesConc(:,:,k1,IFSO4N3) + settle_so4n(:,:,k1)
!... nOC
       settle_so4n(:,:,:) = speciesConc(:,:,:,INOC) * (1.0d0 - exp(-chemintv * lambda))
       speciesConc(:,:,:,INOC) = speciesConc(:,:,:,INOC) - settle_so4n(:,:,:)
       do ik = k1+1, k2
         speciesConc(:,:,ik-1,INOC) = speciesConc(:,:,ik-1,INOC) + settle_so4n(:,:,ik)
       enddo
       speciesConc(:,:,k1,INOC) = speciesConc(:,:,k1,INOC) + settle_so4n(:,:,k1)
!... fOC
       settle_so4n(:,:,:) = speciesConc(:,:,:,IFOC) * (1.0d0 - exp(-chemintv * lambda))
       speciesConc(:,:,:,IFOC) = speciesConc(:,:,:,IFOC) - settle_so4n(:,:,:)
       do ik = k1+1, k2
         speciesConc(:,:,ik-1,IFOC) = speciesConc(:,:,ik-1,IFOC) + settle_so4n(:,:,ik)
       enddo
       speciesConc(:,:,k1,IFOC) = speciesConc(:,:,k1,IFOC) + settle_so4n(:,:,k1)
!... fBC
       settle_so4n(:,:,:) = speciesConc(:,:,:,IFBC) * (1.0d0 - exp(-chemintv * lambda))
       speciesConc(:,:,:,IFBC) = speciesConc(:,:,:,IFBC) - settle_so4n(:,:,:)
       do ik = k1+1, k2
         speciesConc(:,:,ik-1,IFBC) = speciesConc(:,:,ik-1,IFBC) + settle_so4n(:,:,ik)
       enddo
       speciesConc(:,:,k1,IFBC) = speciesConc(:,:,k1,IFBC) + settle_so4n(:,:,k1)
!... bOC
       settle_so4n(:,:,:) = speciesConc(:,:,:,IBOC) * (1.0d0 - exp(-chemintv * lambda))
       speciesConc(:,:,:,IBOC) = speciesConc(:,:,:,IBOC) - settle_so4n(:,:,:)
       do ik = k1+1, k2
         speciesConc(:,:,ik-1,IBOC) = speciesConc(:,:,ik-1,IBOC) + settle_so4n(:,:,ik)
       enddo
       speciesConc(:,:,k1,IBOC) = speciesConc(:,:,k1,IBOC) + settle_so4n(:,:,k1)
!... bBC
       settle_so4n(:,:,:) = speciesConc(:,:,:,IBBC) * (1.0d0 - exp(-chemintv * lambda))
       speciesConc(:,:,:,IBBC) = speciesConc(:,:,:,IBBC) - settle_so4n(:,:,:)
       do ik = k1+1, k2
         speciesConc(:,:,ik-1,IBBC) = speciesConc(:,:,ik-1,IBBC) + settle_so4n(:,:,ik)
       enddo
       speciesConc(:,:,k1,IBBC) = speciesConc(:,:,k1,IBBC) + settle_so4n(:,:,k1)
!.sds... end artificial settling for SO4 aerosols

       do ic = 1, num_species
          concentration(ic)%pArray3D(:,:,:) = speciesConc(:,:,:,ic)
       end do
       deallocate(speciesConc)

!     ----------------------------------------------------
!     Change units from concentration back to mixing ratio
!     (only for gas-species).
!     ----------------------------------------------------

      do ic = 1, num_molefrac
        if (aerosol(ic) <= 0) then
          concentration(ic)%pArray3D(i1:i2,ju1:j2,:) =  &
     &      concentration(ic)%pArray3D(i1:i2,ju1:j2,:) /  &
     &      concentration(imgas_num)%pArray3D(i1:i2,ju1:j2,:)
        end if
      end do

#ifdef GTmodule
      totso4(:,:,:) = 0.0d0

     do ic= 5, 12
       totso4(:,:,:) =totso4(:,:,:)+ concentration(ic)%pArray3D(:,:,:)
     end do

!rsot      totso4(:,:,:) =concentration(5)%pArray3D(:,:,:)+concentration(9)%pArray3D(:,:,:)        ! add aqueous fossil and natural SO4

!    Convert aerosol sulfate from g/g to ug/m3
!
      ! M2 converts to ug and m^3
      cloud_param(:,:,:,7)=totso4(:,:,:)*M2  & 
     &          * concentration(imgas_num)%pArray3D(:,:,:)*WTAIR/AVOGAD
#endif

      Deallocate (c0)
      Deallocate (c1)
      Deallocate (c2)
      Deallocate (cfc)
      Deallocate (cfs)
      Deallocate (clfrac)
      Deallocate (cmfcol)
      Deallocate (height_asl)
      Deallocate (lwccol)
      Deallocate (relhum)
      Deallocate (rhc)

      return

      end subroutine Update_Sulfchem

!-------------------------------------------------------------------------

      end module GmiSolverInterface_mod
