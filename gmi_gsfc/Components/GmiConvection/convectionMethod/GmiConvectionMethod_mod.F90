!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiConvectionMethod_mod 
!
! !INTERFACE:
!
#include "GmiESMF_ErrLog.h"
!
  module GmiConvectionMethod_mod
!
! !USES:
      use Esmf_Mod
      use GmiESMF_ErrorChecking_mod
      use GmiPrintError_mod, only : GmiPrintError
      use GmiESMFrcFileReading_mod, only : rcEsmfReadLogical
      use GmiDiagnosticsMethod_mod  , only : t_Diagnostics
      use GmiDiagnosticsMethod_mod  , only : Get_pr_wet_depos, Get_pr_diag
      use GmiCheckNamelistFile_mod, only : CheckNamelistOptionRange
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID
      use GmiMetFieldsControl_mod, only : t_metFields, Get_mass, &
     &       Get_metdata_name_org, Get_metdata_name_model, Get_gridBoxHeight, &
     &       Get_press3e, Get_press3c, Get_lwi_flags, Get_pbl, Get_cmf,       &
     &       Get_dtrn, Get_eu, Get_ed, Get_md, Get_kel, Get_zmdu, Get_zmeu,   &
     &       Get_zmed, Get_zmmd, Get_zmmu, Get_hkdu, Get_hkeu, Get_hkmu,      &
     &       Get_humidity, Get_met_opt
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: InitializeConvection, RunConvection, FinalizeConvection
      public  :: Get_det_ent, Get_do_downdraft, Get_do_old_ncar, Get_convec_opt
!
! !PUBLIC DATA MEMBERS:
!
      public  :: t_Convection

# include "GmiParameters.h"

      type t_Convection
        private
        integer :: convec_opt   ! Convection option
        logical :: det_ent      ! flag for doing detrainment then entrainment
        logical :: do_downdraft ! flag for doing downdrafts
        logical :: do_old_ncar  ! flag for using old ncar routine
      end type t_Convection

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
! !IROUTINE: readConvectionResourceFile
!     
! !INTERFACE:
!
      subroutine readConvectionResourceFile (self, procID, pr_diag, config)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: procID
      logical, intent(in) :: pr_diag
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config), intent(inOut) :: config
      type (t_Convection), intent(inOut) :: self
!
! !DESCRIPTION:
! Reads the Convection section of the resource file.
!
! !LOCAL VARIABLES: 
      integer :: STATUS, RC
      character(len=ESMF_MAXSTR) :: err_msg, IAm
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = "readConvectionResourceFile"
 
      if (pr_diag) Write(6,*) IAm, 'called by ', procID

      !################################
      ! Begin reading the resource file
      !################################

      ! ------------------------------
      ! convec_opt
      !   0:  no convection
      !   1:  do DAO2 convection
      !   2:  do NCAR convection
      !   3:  do GMAO GEOS4 convection
      !------ ------------------------

      call ESMF_ConfigGetAttribute(config, self%convec_opt, &
     &                label   = "convec_opt:",&
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%det_ent, &
     &               "det_ent:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%do_downdraft, &
     &               "do_downdraft:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(config, self%do_old_ncar, &
     &               "do_old_ncar:", default=.false., rc=STATUS )
      VERIFY_(STATUS)

      !##############################
      ! End reading the resource file
      !##############################

     ! ---------------------------------------------------------------
     ! Check option ranges.  Note that as new options are added, these
     ! range checks will have to be modified.
     ! ---------------------------------------------------------------

      call CheckNamelistOptionRange ('convec_opt' , self%convec_opt, 0, 3)

      return

      end subroutine readConvectionResourceFile
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InitializeConvection
!
! !INTERFACE:
!
      subroutine InitializeConvection (self, gmiDomain, Diagnostics, config)
!
      implicit none
!
! !INPUT PARAMETERS:
      type (t_gmiDomain), intent(in) :: gmiDomain
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config), intent(inOut) :: config
      type (t_Convection), intent(inOut) :: self
!
! !DESCRIPTION:
! Initializes the convection componet.
!
! !LOCAL VARIABLES:
      integer :: STATUS, RC, procID
      logical :: pr_diag
      character(len=ESMF_MAXSTR) :: err_msg, IAm
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = "initializeConvection"

      call Get_procID (gmiDomain, procID)
      call Get_pr_diag (Diagnostics, pr_diag)

      if (pr_diag) Write(6,*) IAm, 'called by ', procID

      call readConvectionResourceFile (self, procID, pr_diag, config)

      return

      end subroutine InitializeConvection
!EOC
!-------------------------------------------------------------------------
!BOP
      subroutine RunConvection (self, Deposition, SpeciesConcentration,        &
     &              gmiClock, gmiGrid, gmiDomain, Diagnostics, metFields,      &
     &              chem_opt, ih2o2_num, ihno3_num, mw, bmass,                 &
     &              REL_SCAV_EFF_new, chem_mecha, i1, ju1, k1)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiDepositionMethod_mod, only : t_Deposition
      use GmiDepositionMethod_mod, only : Get_wet_depos, Set_wet_depos,        &
     &       Get_do_wetdep, Get_do_drydep
      use GmiTimeControl_mod, only : t_GmiClock, Get_gmiTimeStep
      use GmiGrid_mod, only : t_gmiGrid, Get_i1, Get_i2, Get_ju1, Get_j2,      &
     &       Get_k1, Get_k2, Get_ilo, Get_ihi, Get_julo, Get_jhi, Get_ilong,   &
     &       Get_ivert, Get_numSpecies
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration, Set_concentration
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID, Get_mcor

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      type (t_Convection), intent(in) :: self
      type (t_Deposition), intent(inout) :: Deposition
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
      type (t_GmiClock  ), intent(in   ) :: gmiClock
      type (t_gmiGrid   ), intent(in   ) :: gmiGrid 
      type (t_gmiDomain ), intent(in   ) :: gmiDomain 
      type (t_metFields ), intent(in   ) :: metFields
      type (t_Diagnostics), intent(in) :: Diagnostics

      character (len=* ), intent(in   ) :: chem_mecha
      integer           , intent(in   ) :: chem_opt, i1, ju1, k1
      integer           , intent(in   ) :: ih2o2_num
      integer           , intent(in   ) :: ihno3_num
      real*8            , intent(in   ) :: mw(1:)
      real*8            , intent(in   ) :: bmass(i1:,ju1:,k1:)
      real*8            , intent(out  ) :: REL_SCAV_EFF_new(i1:,ju1:,k1:,1:)

      real*8    :: tdt
      logical             :: do_wetdep, do_drydep, pr_diag, pr_wet_depos
      real*8, allocatable :: wet_depos(:,:,:)
      real*8, allocatable :: mcor(:,:)
      integer             :: i2, j2, k2, ilo, ihi, julo, jhi
      integer             :: ilong, ivert, procID, numSpecies, met_opt
      type (t_GmiArrayBundle), pointer :: concentration(:)
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model
      real*8, allocatable :: mass(:,:,:)
      real*8, allocatable :: press3e(:,:,:)
      real*8, allocatable :: press3c(:,:,:)
      real*8, allocatable :: gridBoxHeight(:,:,:)
      integer, allocatable :: lwi_flags(:, :)
      real*8 , allocatable :: pbl        (:, :)
      real*8 , allocatable :: cmf        (:, :, :)
      real*8 , allocatable :: dtrn       (:, :, :)
      real*8 , allocatable :: eu         (:, :, :)
      real*8 , allocatable :: ed         (:, :, :)
      real*8 , allocatable :: md         (:, :, :)
      real*8 , allocatable :: kel        (:, :, :)
      real*8 , allocatable :: zmdu       (:, :, :)
      real*8 , allocatable :: zmeu       (:, :, :)
      real*8 , allocatable :: zmed       (:, :, :)
      real*8 , allocatable :: zmmd       (:, :, :)
      real*8 , allocatable :: zmmu       (:, :, :)
      real*8 , allocatable :: hkdu       (:, :, :)
      real*8 , allocatable :: hkeu       (:, :, :)
      real*8 , allocatable :: hkmu       (:, :, :)
      real*8 , allocatable :: humidity   (:, :, :)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write(6,*) 'runConvection called by ', procID

      call Get_pr_wet_depos(Diagnostics, pr_wet_depos)

      ! Get the GMI grid information
      !call Get_i1   (gmiGrid, i1   )
      call Get_i2   (gmiGrid, i2   )
      !call Get_ju1  (gmiGrid, ju1  )
      call Get_j2   (gmiGrid, j2   )
      !call Get_k1   (gmiGrid, k1   )
      call Get_k2   (gmiGrid, k2   )
      call Get_ilo  (gmiGrid, ilo  )
      call Get_ihi  (gmiGrid, ihi  )
      call Get_julo (gmiGrid, julo )
      call Get_jhi  (gmiGrid, jhi  )
      call Get_ilong(gmiGrid, ilong)
      call Get_ivert(gmiGrid, ivert)
      call Get_numSpecies(gmiGrid, numSpecies)

      allocate(mcor(i1:i2,ju1:j2))
      call Get_mcor(gmiDomain, mcor)

      ! Obtain the model time step
      call Get_gmiTimeStep(gmiClock, tdt)

      call Get_concentration(SpeciesConcentration, concentration)

      call Get_do_wetdep(Deposition, do_wetdep)

      call Get_do_drydep(Deposition, do_drydep)

      if (pr_wet_depos) then
         allocate(wet_depos (i1:i2, ju1:j2, 1:numSpecies))
         call Get_wet_depos(Deposition, wet_depos)
      end if

      call Get_met_opt(metFields, met_opt)
      call Get_metdata_name_org(metFields, metdata_name_org)
      call Get_metdata_name_model(metFields, metdata_name_model)

      allocate(humidity  (i1:i2, ju1:j2, k1:k2))
      call Get_humidity(metFields, humidity)

      allocate(kel(ilo:ihi, julo:jhi, k1:k2))
      call Get_kel(metFields, kel)

      if ((chem_opt /= 0) .or. (self%convec_opt /= 0)) then
         allocate(lwi_flags(i1:i2, ju1:j2))
         call Get_lwi_flags(metFields, lwi_flags)
      end if

      if (met_opt == 3) then
         allocate(pbl(i1:i2, ju1:j2))
         call Get_pbl(metFields, pbl)

         if (self%convec_opt /= 0 .or. do_wetdep) then
            allocate(cmf(i1:i2, ju1:j2, k1:k2))
            call Get_cmf(metFields, cmf)
         end if

         if (self%convec_opt /= 0 .or. do_drydep) then
            allocate(dtrn(i1:i2, ju1:j2, k1:k2))
            call Get_dtrn(metFields, dtrn)
         end if

!
         if (self%convec_opt == 2) then
            allocate(eu (i1:i2, ju1:j2, k1:k2))
            allocate(ed (i1:i2, ju1:j2, k1:k2))
            allocate(md (i1:i2, ju1:j2, k1:k2))
            call Get_eu(metFields, eu)
            call Get_ed(metFields, ed)
            call Get_md(metFields, md)
         end if

         if (self%convec_opt == 3 .and.  &
     &        metdata_name_org(1:4) == 'GMAO' .and.  &
     &        metdata_name_model(1:5) == 'GEOS4') then
            allocate(zmdu      (i1:i2, ju1:j2, k1:k2))
            allocate(zmeu      (i1:i2, ju1:j2, k1:k2))
            allocate(zmed      (i1:i2, ju1:j2, k1:k2))
            allocate(zmmd      (i1:i2, ju1:j2, k1:k2))
            allocate(zmmu      (i1:i2, ju1:j2, k1:k2))
            allocate(hkdu      (i1:i2, ju1:j2, k1:k2))
            allocate(hkeu      (i1:i2, ju1:j2, k1:k2))
            allocate(hkmu      (i1:i2, ju1:j2, k1:k2))

            call Get_zmdu(metFields, zmdu)
            call Get_zmeu(metFields, zmeu)
            call Get_zmed(metFields, zmed)
            call Get_zmmd(metFields, zmmd)
            call Get_zmmu(metFields, zmmu)
            call Get_hkdu(metFields, hkdu)
            call Get_hkeu(metFields, hkeu)
            call Get_hkmu(metFields, hkmu)
         end if
      end if
!
      allocate(mass(i1:i2, ju1:j2, k1:k2))
      allocate(press3c(ilo:ihi, julo:jhi, k1:k2))
      allocate(press3e(ilo:ihi, julo:jhi, k1-1:k2))
      allocate(gridBoxHeight(i1:i2, ju1:j2, k1:k2))

      call Get_mass(metFields, mass)
      call Get_press3e(metFields, press3e)
      call Get_press3c(metFields, press3c)
      call Get_gridBoxHeight(metFields, gridBoxHeight)

!         ==================
          call Update_Convec  &
!         ==================
     &      (chem_mecha, metdata_name_org, metdata_name_model, self%det_ent, &
     &       self%do_downdraft, self%do_old_ncar, do_wetdep,  &
     &       pr_wet_depos, chem_opt, self%convec_opt, ih2o2_num, ihno3_num,  &
     &       lwi_flags, tdt, mw, mcor, pbl, cmf, dtrn, eu,  &
     &       ed, md, zmdu, zmeu, zmed, zmmd, zmmu, hkdu, hkeu, hkmu,    &
     &       gridBoxHeight, mass, wet_depos, kel, press3e, concentration, bmass, &
#ifdef MICRO_AEROSOL
     &       humidity, press3c, REL_SCAV_EFF_new, &
     &       pr_diag, procID, &
     &       i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, ivert, numSpecies)
#else
     &       pr_diag, procID, &
     &       i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, ilong, ivert, numSpecies)
#endif

      if (pr_wet_depos) then
         call Set_wet_depos(Deposition, wet_depos)
         deallocate(wet_depos)
      end if

      deallocate(mcor)
      deallocate(mass)
      deallocate(press3c)
      deallocate(press3e)
      deallocate(gridBoxHeight)

  return

  end subroutine RunConvection

!-------------------------------------------------------------------------

  subroutine FinalizeConvection (self)

  type (t_Convection), intent(inout) :: self

  PRINT*,'  Finalize Convection'

  return

  end subroutine FinalizeConvection

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
  subroutine Get_convec_opt (self, convec_opt)
    integer         , intent(out)  :: convec_opt
    type (t_Convection), intent(in)   :: self
    convec_opt = self%convec_opt
    return
  end subroutine Get_convec_opt
!-------------------------------------------------------------------------
  subroutine Get_det_ent (self, det_ent)
    logical          , intent(out)  :: det_ent
    type (t_Convection), intent(in )  :: self
    det_ent = self%det_ent
    return
  end subroutine Get_det_ent
!-------------------------------------------------------------------------
  subroutine Get_do_downdraft (self, do_downdraft)
    logical          , intent(out)  :: do_downdraft
    type (t_Convection), intent(in )  :: self
    do_downdraft = self%do_downdraft
    return
  end subroutine Get_do_downdraft
!-------------------------------------------------------------------------
  subroutine Get_do_old_ncar (self, do_old_ncar)
    logical          , intent(out)  :: do_old_ncar
    type (t_Convection), intent(in )  :: self
    do_old_ncar = self%do_old_ncar
    return
  end subroutine Get_do_old_ncar
!-------------------------------------------------------------------------
  end module GmiConvectionMethod_mod
