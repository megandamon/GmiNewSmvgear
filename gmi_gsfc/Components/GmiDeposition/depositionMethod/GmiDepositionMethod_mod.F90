!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiDepositionMethod_mod
!
! !INTERFACE:
!
#include "GmiESMF_ErrLog.h"
!
  module GmiDepositionMethod_mod
!
! !USES:
      use GmiGrid_mod, only : t_gmiGrid
      use GmiGrid_mod, only : Get_i1, Get_i2, Get_ju1, Get_j2, Get_k1, Get_k2
      use GmiGrid_mod, only : Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl
      use GmiGrid_mod, only : Get_ilo, Get_ihi, Get_julo, Get_jhi
      use GmiGrid_mod, only : Get_ilong, Get_ivert, Get_numSpecies
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID, Get_mcor, &
     &       Get_coscen
      use GmiDiagnosticsMethod_mod  , only : t_Diagnostics, Get_pr_diag
      use GmiDiagnosticsMethod_mod  , only : Get_pr_dry_depos, Get_pr_wet_depos
      use Esmf_Mod
      use GmiESMF_ErrorChecking_mod
      use GmiPrintError_mod, only : GmiPrintError
      use GmiESMFrcFileReading_mod, only : rcEsmfReadLogical, rcEsmfReadTable
      use GmiMetFieldsControl_mod, only : t_metFields, Get_mass,               &
     &       Get_metdata_name_org, Get_metdata_name_model, Get_gridBoxHeight,  &
     &       Get_press3c, Get_press3e, Get_met_opt, Get_lwi_flags, Get_ustar, &
     &       Get_fracCloudCover, Get_radswg, Get_surf_air_temp, Get_surf_rough,&
     &       Get_humidity, Get_kel, Get_tot_precip, Get_con_precip, Get_moistq,&
     &       Get_rain_zm, Get_rain_hk, Get_rain_ls
!
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: InitializeDeposition
  public  :: RunDryDeposition
  public  :: RunWetDeposition
  public  :: RunSimpleDeposition
  public  :: FinalizeDeposition
!
  public  :: Get_num_ks_sdep
  public  :: Get_do_drydep
  public  :: Get_do_wetdep
  public  :: Get_do_simpledep
  public  :: Get_wetdep_eff
  public  :: Get_dry_depos, Get_wet_depos
  public  :: Set_dry_depos, Set_wet_depos
!
! !PUBLIC DATA MEMBERS:
!
  public  :: t_Deposition
!
# include "GmiParameters.h"
!
  type t_Deposition
    private
    real*8  :: pressDryDep
    integer :: num_ks_sdep                        ! number of vertical layers
                                                  ! to apply 2 day loss factor
                                                  ! to in simple deposition
    logical :: do_drydep                          ! do dry    deposition?
    logical :: do_wetdep                          ! do wet    deposition?
    logical :: do_simpledep                       ! do simple deposition?
    real*8  :: wetdep_eff(MAX_NUM_CONST)          ! wet deposition (scavenging)
                                                  ! efficiencies
    real*8, pointer :: dry_depos(:,:,:) => null() ! dry deposition accumulated
                                                  ! since last output (kg/m^2)
    real*8, pointer :: wet_depos(:,:,:) => null() ! wet deposition accumulated
                                                  ! since last output (kg/m^2)
  end type t_Deposition
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
!
! !IROUTINE: readDepositionResourceFile
!
! !INTERFACE:
!
      subroutine readDepositionResourceFile (self, procID, pr_diag, config)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: procID
      logical, intent(in) :: pr_diag
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config), intent(inOut) :: config
      type (t_Deposition), intent(inOut) :: self
!
! !DESCRIPTION:
! Reads the Deposition section of the resource file.
!
! !LOCAL VARIABLES:
      integer :: STATUS, RC
      character(len=ESMF_MAXSTR) :: err_msg, IAm
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = "readDepositionResourceFile"
!
      if (pr_diag) Write(6,*) IAm, 'called by ', procID
!
      !################################
      ! Begin reading the resource file
      !################################
!
      call ESMF_ConfigGetAttribute(config, self%num_ks_sdep, &
     &                label   = "num_ks_sdep:",&
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)
!
      call rcEsmfReadLogical(config, self%do_drydep, &
     &               "do_drydep:", default=.false., rc=STATUS )
      VERIFY_(STATUS)
!
      call rcEsmfReadLogical(config, self%do_wetdep, &
     &               "do_wetdep:", default=.false., rc=STATUS )
      VERIFY_(STATUS)
!
      call rcEsmfReadLogical(config, self%do_simpledep, &
     &               "do_simpledep:", default=.false., rc=STATUS )
      VERIFY_(STATUS)
!
     ! ----------------------------------------------------------------
     ! wetdep_eff : wet deposition (scavenging) efficiencies; should be
     !              set to values between 0.0 and 1.0 for each species
     ! ----------------------------------------------------------------
!
      self%wetdep_eff(:) = 0.0d0
!
      ! 3.12.2012 Megan Rose Damon
      ! The below call was causing problems when running the code.
      ! Specifically, a subsequent call to ESMF_ConfigFindLabel
      ! was returning inconsistent results across processors and causing
      ! the code to read a part of the resource file that doesn't exist.
      !call rcEsmfReadTable(config, self%wetdep_eff, "wetdep_eff::", rc=STATUS)
      !VERIFY_(STATUS)
!
      !##############################
      ! End reading the resource file
      !##############################
!
      if (self%do_simpledep .and. self%do_drydep) then
         err_msg = 'do_simpledep/do_drydep problem in Check_Nlvalue.'
         call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if
!
      if (self%do_simpledep .and. self%do_wetdep) then
         err_msg = 'do_simpledep/do_wetdep problem in Check_Nlvalue.'
         call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if
!
      return
!
      end subroutine readDepositionResourceFile
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InitializeDeposition
!
! !INTERFACE:
!
      subroutine initializeDeposition (self, gmiGrid, gmiDomain, Diagnostics, &
     &                     config)
!
      implicit none
!
! !INPUT PARAMETERS:
      type (t_gmiGrid   ) , intent(in) :: gmiGrid
      type (t_gmiDomain)  , intent(in) :: gmiDomain
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config), intent(inOut) :: config
      type (t_Deposition), intent(inOut) :: self
!
! !DESCRIPTION:
! Initializes the Deposition component.
!
! !LOCAL VARIABLES:
      integer :: i1, i2, ju1, j2, numSpecies
      logical :: pr_diag, pr_dry_depos, pr_wet_depos
      integer :: STATUS, RC, procID
      character(len=ESMF_MAXSTR) :: err_msg, IAm
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = "initializeDeposition"
!
      call Get_procID (gmiDomain, procID)
      call Get_pr_diag (Diagnostics, pr_diag)
!
      if (pr_diag) Write(6,*) IAm, 'called by ', procID
!
      call readDepositionResourceFile (self, procID, pr_diag, config)
!
      call Get_pr_dry_depos(Diagnostics, pr_dry_depos)
      call Get_pr_wet_depos(Diagnostics, pr_wet_depos)
!
      call Get_i1 (gmiGrid, i1 )
      call Get_i2 (gmiGrid, i2 )
      call Get_ju1(gmiGrid, ju1)
      call Get_j2 (gmiGrid, j2 )
      call Get_numSpecies(gmiGrid, numSpecies)
!
      self%pressDryDep = 1.5d5

      if (pr_dry_depos) then
         allocate(self%dry_depos(i1:i2, ju1:j2, 1:numSpecies))
         self%dry_depos(:,:,:) = 0.0d0
      endif
!
      if (pr_wet_depos) then
         allocate(self%wet_depos(i1:i2, ju1:j2, 1:numSpecies))
         self%wet_depos(:,:,:) = 0.0d0
      endif
!
      return
!
      end subroutine InitializeDeposition
!EOC
!-------------------------------------------------------------------------
!BOP
      subroutine RunDryDeposition &
     &  (self, Emission, SpeciesConcentration, gmiClock, gmiGrid, gmiDomain, &
     &   Diagnostics, metFields, cosSolarZenithAngle, &
     &   diffaer, s_radius, s_velocity, BoxHeightCenter, BoxHeightEdge, &
     &   chem_opt, mw, i1, i2, ju1, j2, numSpecies)
!
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiEmissionMethod_mod    , only : t_Emission
      use GmiEmissionMethod_mod    , only : Get_ireg, Get_iland, Get_iuse, Get_xlai
      use GmiTimeControl_mod       , only : t_GmiClock, Get_gmiTimeStep
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration
!
      implicit none
!
#     include "gmi_emiss_constants.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      type(t_Deposition), intent(inout) :: self
      type(t_Emission  ), intent(inout) :: Emission
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
      type(t_gmiGrid   ), intent(in   ) :: gmiGrid
      type(t_gmiDomain ), intent(in   ) :: gmiDomain
      type (t_Diagnostics), intent(in) :: Diagnostics
      integer, intent(in) :: chem_opt, i1, i2, ju1, j2, numSpecies
      real*8 , intent(in) :: mw(1:)
      real*8 , intent(in) :: cosSolarZenithAngle (i1:i2,ju1:j2)
      real*8 , intent(in) :: BoxHeightCenter(i1:i2,ju1:j2)
      real*8 , intent(in) :: BoxHeightEdge  (i1:i2,ju1:j2)
      real*8 , intent(in) :: diffaer   (i1:i2,ju1:j2,1:numSpecies)
      real*8 , intent(in) :: s_radius  (i1:i2,ju1:j2,1:numSpecies)
      real*8 , intent(in) :: s_velocity(i1:i2,ju1:j2,1:numSpecies)
      type (t_GmiClock ), intent(in) :: gmiClock
      type (t_metFields), intent(in) :: metFields
!
! Local variable
!
      real*8 , allocatable :: mcor(:,:)
      real*8 , allocatable :: mass(:,:,:)
      integer, allocatable :: iland(:,:,:), iuse(:,:,:), ireg(:,:)
      real*8 , allocatable :: dry_depos (:,:,:), xlai(:,:,:)
      real*8  :: tdt
      integer :: k1, k2, ilong, procID
      integer :: ilo, ihi, julo, jhi, ju1_gl, j2_gl, i1_gl, i2_gl, met_opt
      type (t_GmiArrayBundle), pointer :: concentration(:)
      logical :: pr_diag, pr_dry_depos
      integer, allocatable :: lwi_flags (:,:)
      real*8 , allocatable :: fracCloudCover (:,:)
      real*8 , allocatable :: radswg(:,:), surf_air_temp(:,:)
      real*8 , allocatable :: surf_rough(:,:), ustar(:,:)
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID (gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)
!
      if (pr_diag) Write(6,*)'runDryDeposition called by ', procID
!
      call Get_pr_dry_depos(Diagnostics, pr_dry_depos)
!
      ! Get the GMI grid information
      !call Get_i1    (gmiGrid, i1   )
      !call Get_i2    (gmiGrid, i2   )
      !call Get_ju1   (gmiGrid, ju1  )
      !call Get_j2    (gmiGrid, j2   )
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
      !call Get_numSpecies (gmiGrid, numSpecies)
!
      ! Obtain the model time step
      call Get_gmiTimeStep(gmiClock, tdt)
!
      call Get_concentration(SpeciesConcentration, concentration)
!
      call Get_met_opt(metFields, met_opt)
!
      allocate(lwi_flags (i1:i2, ju1:j2))
      call Get_lwi_flags(metFields, lwi_flags)
!
      if (met_opt == 3) then
         allocate(fracCloudCover (i1:i2, ju1:j2))
         call Get_fracCloudCover(metFields, fracCloudCover)
!
         allocate(radswg(i1:i2, ju1:j2))
         call Get_radswg(metFields, radswg)
!
         allocate(surf_air_temp(i1:i2, ju1:j2))
         call Get_surf_air_temp(metFields, surf_air_temp)
!
         allocate(surf_rough(i1:i2, ju1:j2))
         call Get_surf_rough(metFields, surf_rough)
!
         allocate(ustar(i1:i2, ju1:j2))
         call Get_ustar(metFields, ustar)
      end if
!
      allocate(mass (i1:i2, ju1:j2, k1:k2))
      call Get_mass(metFields, mass)
!
      allocate(ireg (i1:i2, ju1:j2))
      allocate(iland(i1:i2, ju1:j2, NTYPE))
      allocate(iuse (i1:i2, ju1:j2, NTYPE))
      allocate(xlai (i1:i2, ju1:j2, NTYPE))
!
      call Get_ireg  (Emission, ireg )
      call Get_iland (Emission, iland)
      call Get_iuse  (Emission, iuse )
      call Get_xlai  (Emission, xlai )
!
      allocate(mcor(i1:i2,ju1:j2))
      call Get_mcor(gmiDomain, mcor)
!
!           ==================
      call Update_Drydep  &
!           ==================
     &        (pr_dry_depos, lwi_flags, mcor, cosSolarZenithAngle, &
     &         fracCloudCover, radswg, surf_air_temp, &
     &         surf_rough, ustar, mass, concentration, self%dry_depos, &
     &         diffaer, s_radius, s_velocity, BoxHeightCenter, BoxHeightEdge, &
     &         ireg, iland, iuse, xlai, self%pressDryDep, &
     &         pr_diag, procID, chem_opt, tdt,    &
     &         i1, i2, ju1, j2, k1, k2, ilong, ilo, ihi, julo, jhi, ju1_gl,   &
     &         j2_gl, i1_gl, i2_gl, mw, numSpecies)
!
!
      if (pr_dry_depos) then
         where (self%dry_depos(:,:,:) < 1.0d-30)
              self%dry_depos(:,:,:) = 1.0d-30
         end where
      end if
!
      deallocate(mass)
      deallocate(mcor)
      deallocate(ireg )
      deallocate(iland)
      deallocate(iuse )
      deallocate(xlai )
!
  return
!
  end subroutine RunDryDeposition
!
!-------------------------------------------------------------------------
!
  subroutine RunWetDeposition &
     &        (self, SpeciesConcentration, gmiClock, gmiGrid, gmiDomain, &
     &         Diagnostics, metFields, chem_opt, ih2o2_num, ihno3_num, mw, &
     &         REL_SCAV_EFF_new, i1, ju1, k1)
!
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiTimeControl_mod       , only : t_GmiClock, Get_gmiTimeStep
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration, Set_concentration
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      type(t_Deposition), intent(inout) :: self
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
      type(t_gmiGrid   ), intent(in   ) :: gmiGrid
      type(t_gmiDomain ), intent(in   ) :: gmiDomain
      type(t_metFields ), intent(in   ) :: metFields
      type (t_Diagnostics), intent(in) :: Diagnostics
!
      integer           , intent(in   ) :: chem_opt, i1, ju1, k1
      integer           , intent(in   ) :: ih2o2_num
      integer           , intent(in   ) :: ihno3_num
      real*8            , intent(in   ) :: mw(1:)
      type (t_GmiClock), intent(in   ) :: gmiClock
      real*8           , intent(in   ) :: REL_SCAV_EFF_new(i1:,ju1:,k1:,1:)
!
      real*8, allocatable :: mcor(:,:), coscen(:)
      real*8  :: tdt
      integer :: i2, j2, k2, ilong, ivert, procID, met_opt
      integer :: ilo, ihi, julo, jhi, ju1_gl, j2_gl, numSpecies
      type (t_GmiArrayBundle), pointer :: concentration(:)
      logical :: pr_diag, pr_wet_depos
!
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model
      real*8 , allocatable :: mass(:,:,:)
      real*8 , allocatable :: press3c(:,:,:), press3e(:,:,:)
      real*8 , allocatable :: gridBoxHeight(:,:,:)
      real*8 , allocatable :: con_precip (:,:), tot_precip(:,:)
      real*8 , allocatable :: moistq(:,:,:), rain_zm(:,:,:)
      real*8 , allocatable :: rain_hk(:,:,:), rain_ls(:,:,:)
      real*8 , allocatable :: kel(:,:,:)
      real*8 , allocatable :: humidity   (:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)
!
      if (pr_diag) Write(6,*)'runWetDeposition called by ', procID
!
      call Get_pr_wet_depos(Diagnostics, pr_wet_depos)
!
      ! Get the GMI grid information
      !call Get_i1    (gmiGrid, i1   )
      call Get_i2    (gmiGrid, i2   )
      !call Get_ju1   (gmiGrid, ju1  )
      call Get_j2    (gmiGrid, j2   )
      !call Get_k1    (gmiGrid, k1   )
      call Get_k2    (gmiGrid, k2   )
      call Get_ju1_gl(gmiGrid, ju1_gl)
      call Get_j2_gl (gmiGrid, j2_gl )
      call Get_ilo   (gmiGrid, ilo  )
      call Get_ihi   (gmiGrid, ihi  )
      call Get_julo  (gmiGrid, julo )
      call Get_jhi   (gmiGrid, jhi  )
      call Get_ilong (gmiGrid, ilong)
      call Get_ivert (gmiGrid, ivert)
      call Get_numSpecies (gmiGrid, numSpecies)
!
      ! Obtain the model time step
      call Get_gmiTimeStep(gmiClock, tdt)
!
      call Get_met_opt(metFields, met_opt)
      call Get_metdata_name_org(metFields, metdata_name_org)
      call Get_metdata_name_model(metFields, metdata_name_model)
!
      allocate(kel(ilo:ihi, julo:jhi, k1:k2))
      call Get_kel(metFields, kel)
!
      allocate(humidity   (i1:i2, ju1:j2, k1:k2))
      call Get_humidity(metFields, humidity)
!
      if (met_opt == 3) then
         allocate(con_precip (i1:i2, ju1:j2))
         call Get_con_precip(metFields, con_precip)
!
         allocate(tot_precip(i1:i2, ju1:j2))
         call Get_tot_precip(metFields, tot_precip)
!
         allocate(moistq(i1:i2, ju1:j2, k1:k2))
         call Get_moistq(metFields, moistq)
!
         if (metdata_name_org(1:4) == 'GMAO' .and.  &
     &      (metdata_name_model(1:5) == 'GEOS4' .or. &
     &       metdata_name_model(1:5) == 'GEOS5')) then
            allocate(rain_zm(i1:i2, ju1:j2, k1:k2))
            call Get_rain_zm(metFields, rain_zm)
!
            allocate(rain_hk(i1:i2, ju1:j2, k1:k2))
            call Get_rain_hk(metFields, rain_hk)
!
            allocate(rain_ls(i1:i2, ju1:j2, k1:k2))
            call Get_rain_ls(metFields, rain_ls)
         end if
      end if
!
      allocate(mass(i1:i2, ju1:j2, k1:k2))
      allocate(press3c(ilo:ihi, julo:jhi, k1:k2))
      allocate(press3e(ilo:ihi, julo:jhi, k1-1:k2))
      allocate(gridBoxHeight(i1:i2, ju1:j2, k1:k2))
!
      call Get_mass(metFields, mass)
      call Get_press3c(metFields, press3c)
      call Get_press3e(metFields, press3e)
      call Get_gridBoxHeight(metFields, gridBoxHeight)
!
      call Get_concentration(SpeciesConcentration, concentration)
!
      allocate(mcor(i1:i2,ju1:j2))
      call Get_mcor(gmiDomain, mcor)
!
      allocate(coscen(ju1_gl:j2_gl))
      call Get_coscen(gmiDomain, coscen)
!
!           ==================
      call Update_Wetdep  &
!           ==================
     &        (pr_wet_depos, metdata_name_org, metdata_name_model,  &
     &         chem_opt, ih2o2_num, ihno3_num, tdt, mw, coscen,  &
     &         con_precip, tot_precip, mcor, gridBoxHeight, mass, moistq,  &
     &         rain_zm, rain_hk, rain_ls, &
     &         kel, press3c, press3e, concentration, self%wet_depos,  &
#ifdef MICRO_AEROSOL
     &         humidity, REL_SCAV_EFF_new, &
     &         pr_diag, procID, i1, i2, ju1, j2, k1, k2, &
     &         ilo, ihi, julo, jhi, ilong, ivert, ju1_gl, j2_gl, numSpecies)
#else
     &         pr_diag, procID, i1, i2, ju1, j2, k1, k2, &
     &         ilo, ihi, julo, jhi, ilong, ivert, ju1_gl, j2_gl, numSpecies)
#endif
!
      if (pr_wet_depos) then
         where (self%wet_depos(:,:,:) < 1.0d-30)
              self%wet_depos(:,:,:) = 1.0d-30
         end where
      end if
!
      deallocate(mcor)
      deallocate(coscen)
      deallocate(mass)
      deallocate(press3c)
      deallocate(press3e)
      deallocate(gridBoxHeight)
!
  return
!
  end subroutine RunWetDeposition
!
!-------------------------------------------------------------------------
!
      subroutine RunSimpleDeposition &
     &      (self, SpeciesConcentration, gmiClock, gmiGrid, gmiDomain, &
     &       Diagnostics, metFields, &
     &       ibrono2_num, ih2o2_num, ihcl_num, ihno3_num, imgas_num,  &
     &       initrogen_num, ioxygen_num)
!
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiTimeControl_mod       , only : t_GmiClock, Get_gmiTimeStep
      use GmiGrid_mod              , only : t_gmiGrid
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration, Set_concentration
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      type (t_Deposition), intent(inout) :: self
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
      type (t_GmiClock  ), intent(in   ) :: gmiClock
      type (t_gmiGrid   ), intent(in   ) :: gmiGrid
      type (t_gmiDomain ), intent(in   ) :: gmiDomain
      type (t_metFields), intent(in   ) :: metFields
      type (t_Diagnostics), intent(in) :: Diagnostics
      integer, intent(in   ) :: ibrono2_num, ih2o2_num, ihcl_num, ihno3_num, imgas_num
      integer, intent(in   ) :: initrogen_num, ioxygen_num
!
      integer :: num_ks_sdep, procID
      real*8  :: tdt
      integer :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, numSpecies
      type (t_GmiArrayBundle), pointer :: concentration(:)
      logical :: pr_diag
      real*8 , allocatable :: press3c(:,:,:)
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID (gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)
!
      if (pr_diag) Write(6,*)'runSimpleDeposition called by ', procID
!
      ! Get the GMI grid information
      call Get_i1    (gmiGrid, i1   )
      call Get_i2    (gmiGrid, i2   )
      call Get_ju1   (gmiGrid, ju1  )
      call Get_j2    (gmiGrid, j2   )
      call Get_k1    (gmiGrid, k1   )
      call Get_k2    (gmiGrid, k2   )
      call Get_ilo   (gmiGrid, ilo  )
      call Get_ihi   (gmiGrid, ihi  )
      call Get_julo  (gmiGrid, julo )
      call Get_jhi   (gmiGrid, jhi  )
      call Get_numSpecies (gmiGrid, numSpecies)
!
      ! Obtain the model time step
      call Get_gmiTimeStep(gmiClock, tdt)
!
      call Get_num_ks_sdep(self, num_ks_sdep)
!
      call Get_concentration(SpeciesConcentration, concentration)
!
      allocate(press3c(ilo:ihi, julo:jhi, k1:k2))
      call Get_press3c(metFields, press3c)
!
!         =====================
          call Update_Simpledep  &
!         =====================
     &      (ibrono2_num, ih2o2_num, ihcl_num, ihno3_num, imgas_num,  &
     &       initrogen_num, ioxygen_num, num_ks_sdep, tdt, press3c,  &
     &       concentration, pr_diag, procID, &
     &       i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, numSpecies)
!
!
      return
!
      end subroutine RunSimpleDeposition
!EOC
!-------------------------------------------------------------------------
!
  subroutine FinalizeDeposition (self)
!
  type (t_Deposition)   , intent(inout) :: self
!
  PRINT*,'  Finalize Deposition'
!
  return
!
  end subroutine FinalizeDeposition
!
!-------------------------------------------------------------------------
  subroutine Get_dry_depos (self, dry_depos)
    real*8             , intent(out)  :: dry_depos(:,:,:)
    type (t_Deposition), intent(in)   :: self
    dry_depos(:,:,:) = self%dry_depos(:,:,:)
    return
  end subroutine Get_dry_depos
!-------------------------------------------------------------------------
  subroutine Get_wet_depos (self, wet_depos)
    real*8             , intent(out)  :: wet_depos(:,:,:)
    type (t_Deposition), intent(in)   :: self
    wet_depos(:,:,:) = self%wet_depos(:,:,:)
    return
  end subroutine Get_wet_depos
!-------------------------------------------------------------------------
  subroutine Set_dry_depos (self, dry_depos)
    real*8             , intent(in)    :: dry_depos(:,:,:)
    type (t_Deposition), intent(inout) :: self
    self%dry_depos(:,:,:) = dry_depos(:,:,:)
    return
  end subroutine Set_dry_depos
!-------------------------------------------------------------------------
  subroutine Set_wet_depos (self, wet_depos)
    real*8             , intent(in)    :: wet_depos(:,:,:)
    type (t_Deposition), intent(inout) :: self
    self%wet_depos(:,:,:) = wet_depos(:,:,:)
    return
  end subroutine Set_wet_depos
!-------------------------------------------------------------------------
  subroutine Get_num_ks_sdep (self, num_ks_sdep)
    integer         , intent(out)  :: num_ks_sdep
    type (t_Deposition), intent(in)   :: self
    num_ks_sdep = self%num_ks_sdep
    return
  end subroutine Get_num_ks_sdep
!-------------------------------------------------------------------------
  subroutine Get_do_drydep (self, do_drydep)
    logical          , intent(out)  :: do_drydep
    type (t_Deposition), intent(in )  :: self
    do_drydep = self%do_drydep
    return
  end subroutine Get_do_drydep
!-------------------------------------------------------------------------
  subroutine Get_do_wetdep (self, do_wetdep)
    logical          , intent(out)  :: do_wetdep
    type (t_Deposition), intent(in )  :: self
    do_wetdep = self%do_wetdep
    return
  end subroutine Get_do_wetdep
!-------------------------------------------------------------------------
  subroutine Get_do_simpledep (self, do_simpledep)
    logical          , intent(out)  :: do_simpledep
    type (t_Deposition), intent(in )  :: self
    do_simpledep = self%do_simpledep
    return
  end subroutine Get_do_simpledep
!-------------------------------------------------------------------------
  subroutine Get_wetdep_eff (self, wetdep_eff)
    real*8           , intent(out)  :: wetdep_eff(:)
    type (t_Deposition), intent(in )  :: self
    wetdep_eff(:) = self%wetdep_eff(:)
    return
  end subroutine Get_wetdep_eff
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  end module GmiDepositionMethod_mod
!
