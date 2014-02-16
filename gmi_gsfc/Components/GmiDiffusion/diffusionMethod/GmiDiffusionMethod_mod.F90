!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiDiffusionMethod_mod 
!
! !INTERFACE:
!
#include "GmiESMF_ErrLog.h"
!
      module GmiDiffusionMethod_mod
!
! !USES:
      use ESMF_Mod
      use GmiESMF_ErrorChecking_mod
      use GmiDiagnosticsMethod_mod  , only : t_Diagnostics, Get_pr_diag
      use GmiDomainDecomposition_mod, only : t_gmiDomain, Get_procID
      use GmiMetFieldsControl_mod, only : t_metFields, Get_gridBoxHeight, &
     &       Get_mass, Get_tropopausePress, Get_press3e, Get_press3c, Get_kzz, &
     &       Get_pctm1, Get_pbl, Set_kzz
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private 
      public  :: InitializeDiffusion
      public  :: RunDiffusion
      public  :: FinalizeDiffusion
      public  :: Get_diffu_opt
      public  :: Get_vert_diffu_coef
      public  :: Get_pbl_mixing_tau
!
! !PUBLIC DATA MEMBERS:
      public  :: t_Diffusion
!
# include "GmiParameters.h"

      type t_Diffusion
        private
        integer :: diffu_opt        ! diffusion option
        real*8  :: vert_diffu_coef  ! Scalar vertical diffusion coefficient  (m^2/s)
        real*8  :: pbl_mixing_tau   ! PBL mixing tau
      end type t_Diffusion

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
! !IROUTINE: readDiffusionResourceFile
!
! !INTERFACE:
!
      subroutine readDiffusionResourceFile (self, gmiDomain, Diagnostics, config)
!
! !USES:
      use GmiCheckNamelistFile_mod, only : CheckNamelistOptionRange
!
      implicit none
!
! !INPUT PARAMETERS:
      type (t_gmiDomain  ), intent(in) :: gmiDomain
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config), intent(inOut) :: config
      type (t_Diffusion), intent(inOut) :: self
!
! !DESCRIPTION:
! Reads in Diffusion related variables from the resource file.
!
! !LOCAL VARIABLES:
      integer                    :: procID, STATUS, RC
      logical                    :: pr_diag
      character(len=ESMF_MAXSTR) :: IAm 
!EOP
!-------------------------------------------------------------------------
!BOC
      IAm = "readDiffusionResourceFile"

      call Get_procID (gmiDomain, procID)
      call Get_pr_diag (Diagnostics, pr_diag)

      if (pr_diag) Write(6,*) IAm, 'called by ', procID

      !################################
      ! Begin reading the resource file
      !################################

      ! --------------------------------
      ! diffu_opt
      !   0:  no diffusion
      !   1:  do DAO2 vertical diffusion
      ! --------------------------------

      call ESMF_ConfigGetAttribute(config, self%diffu_opt, &
     &                label   = "diffu_opt:",&
     &                default = 0, rc=STATUS )
      VERIFY_(STATUS)

      ! pbl_mixing_tau : scalar vertical diffusion coefficient (m^2/s)

      call ESMF_ConfigGetAttribute(config, self%pbl_mixing_tau, &
     &                label   = "pbl_mixing_tau:",&
     &                default = 3600.0d0, rc=STATUS )
      VERIFY_(STATUS)

      ! ---------------------------------------------------------------
      ! vert_diffu_coef : scalar vertical diffusion coefficient (m^2/s)
      ! ---------------------------------------------------------------

      call ESMF_ConfigGetAttribute(config, self%vert_diffu_coef, &
     &                label   = "vert_diffu_coef:",&
     &                default = 0.0d0, rc=STATUS )
      VERIFY_(STATUS)

    ! ---------------------------------------------------------------
    ! Check option ranges.  Note that as new options are added, these
    ! range checks will have to be modified.
    ! ---------------------------------------------------------------

      call CheckNamelistOptionRange ('diffu_opt       ', self%diffu_opt, 0, 2)

      return

      end subroutine readDiffusionResourceFile
!EOC
!-------------------------------------------------------------------------
!BOP
! 
! !IROUTINE: InitializeDiffusion
!
! !INTERFACE:
!
      subroutine initializeDiffusion (self, gmiDomain, Diagnostics, config)
!
      implicit none
!
! !INPUT PARAMETERS:
      type (t_gmiDomain  ), intent(in) :: gmiDomain
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config), intent(inOut) :: config
      type (t_Diffusion), intent(inOut) :: self
!
! !DESCRIPTION:
! Initialize the Diffusion component.
!
! !LOCAL VARIABLES:
      integer                    :: procID
      logical                    :: pr_diag
      character(len=ESMF_MAXSTR) :: IAm 
!EOP
!-------------------------------------------------------------------------
!BOC
      IAm = "initializeDiffusion"

      call Get_procID (gmiDomain, procID)
      call Get_pr_diag (Diagnostics, pr_diag)

      if (pr_diag) Write(6,*) IAm, 'called by ', procID

      call readDiffusionResourceFile (self, gmiDomain, Diagnostics, config)

      return

      end subroutine initializeDiffusion
!EOC
!-------------------------------------------------------------------------
!BOP
!
! IROUTINE: runDiffusion
!
! !INTERFACE:
!
      subroutine runDiffusion (self, SpeciesConcentration, gmiClock, gmiGrid,  &
     &              gmiDomain, Diagnostics, metFields)
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiTimeControl_mod       , only : t_GmiClock, Get_gmiTimeStep
      use GmiGrid_mod              , only : t_gmiGrid, Get_numSpecies
      use GmiGrid_mod              , only : Get_i1, Get_i2, Get_ju1, Get_j2,   &
             Get_k1, Get_k2, Get_ilo, Get_ihi, Get_julo, Get_jhi
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiSpcConcentrationMethod_mod, only : Get_concentration
!
      implicit none
!
! !INPUT PARAMETERS:
      type(t_gmiGrid     ), intent(in) :: gmiGrid  
      type(t_GmiClock    ), intent(in) :: gmiClock
      type(t_GmiDomain   ), intent(in) :: gmiDomain
      type(t_Diffusion   ), intent(in) :: self
      type (t_Diagnostics), intent(in) :: Diagnostics
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_metFields),            intent(inOut) :: metFields
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
! Runs the diffusion component.
!
! !LOCAL VARIABLES:
      integer :: numSpecies, procID
      real*8  :: tdt
      integer :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      logical :: pr_diag
      real*8 , allocatable :: tropopausePress(:, :), pbl(:,:)
      real*8 , allocatable :: gridBoxHeight(:,:,:) , pctm1(:,:)
      real*8 , allocatable :: mass   (:,:,:)       , kzz(:,:,:)
      real*8 , allocatable :: press3c(:,:,:)       , press3e(:, :, :)
      type (t_GmiArrayBundle), pointer :: concentration(:)
      character(len=ESMF_MAXSTR) :: IAm 
!EOP
!-------------------------------------------------------------------------
!BOC
      IAm = "runDiffusion"

      call Get_procID(gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write(6,*) Iam, 'called by ', procID

      ! Get the GMI grid information
      call Get_i1   (gmiGrid, i1   )
      call Get_i2   (gmiGrid, i2   )
      call Get_ju1  (gmiGrid, ju1  )
      call Get_j2   (gmiGrid, j2   )
      call Get_k1   (gmiGrid, k1   )
      call Get_k2   (gmiGrid, k2   )
      call Get_ilo  (gmiGrid, ilo  )
      call Get_ihi  (gmiGrid, ihi  )
      call Get_julo (gmiGrid, julo )
      call Get_jhi  (gmiGrid, jhi  )
      call Get_numSpecies  (gmiGrid, numSpecies  )

      ! Obtain model time step
      call Get_gmiTimeStep (gmiClock, tdt)

      call Get_concentration(SpeciesConcentration, concentration)

      if (self%diffu_opt == 1) then
         allocate(pbl(i1:i2, ju1:j2))
         allocate(kzz(i1:i2, ju1:j2, k1:k2))
         allocate(pctm1(ilo:ihi, julo:jhi))
         allocate(press3c(ilo:ihi, julo:jhi, k1:k2))
         allocate(press3e(ilo:ihi, julo:jhi, k1-1:k2))
         allocate(tropopausePress(i1:i2, ju1:j2))

         call Get_pbl(metFields, pbl)
         call Get_kzz(metFields, kzz)
         call Get_pctm1  (metFields, pctm1  )
         call Get_press3c(metFields, press3c)
         call Get_press3e(metFields, press3e)
         call Get_tropopausePress(metFields, tropopausePress)

         call Update_Diffu (tdt, self%vert_diffu_coef, pbl, tropopausePress, &
     &               kzz, press3c, press3e, pctm1, concentration, &
     &               pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, &
     &               jhi, numSpecies)

         call Set_kzz(metFields, kzz)

         deallocate(pbl)
         deallocate(kzz)
         deallocate(pctm1)
         deallocate(press3e)
         deallocate(press3c)
         deallocate(tropopausePress)

      else if (self%diffu_opt == 2) then
         allocate(pbl(i1:i2, ju1:j2))
         allocate(mass(i1:i2, ju1:j2, k1:k2))
         allocate(gridBoxHeight(i1:i2, ju1:j2, k1:k2))

         call Get_pbl(metFields, pbl)
         call Get_mass(metFields, mass)
         call Get_gridBoxHeight(metFields, gridBoxHeight)

         call Update_PBL_Mixing  &
     &        (tdt, self%pbl_mixing_tau, pbl, gridBoxHeight, mass, concentration, &
     &         pr_diag, procID, i1, i2, ju1, j2, k1, k2, numSpecies)

         deallocate(pbl)
         deallocate(mass)
         deallocate(gridBoxHeight)
      end if

      return

      end subroutine RunDiffusion
!EOC
!-------------------------------------------------------------------------
      subroutine FinalizeDiffusion (self)

      type (t_Diffusion), intent(inout) :: self

      PRINT*,'  Finalize Diffusion'

      return

      end subroutine FinalizeDiffusion
!-------------------------------------------------------------------------
  subroutine Get_diffu_opt (self, diffu_opt)
    integer         , intent(out)  :: diffu_opt
    type (t_Diffusion), intent(in)   :: self
    diffu_opt = self%diffu_opt
    return
  end subroutine Get_diffu_opt
!-------------------------------------------------------------------------
  subroutine Get_vert_diffu_coef (self, vert_diffu_coef)
    real*8          , intent(out)  :: vert_diffu_coef
    type (t_Diffusion), intent(in)   :: self
    vert_diffu_coef = self%vert_diffu_coef
    return
  end subroutine Get_vert_diffu_coef
!-------------------------------------------------------------------------
  subroutine Get_pbl_mixing_tau (self, pbl_mixing_tau)
    real*8          , intent(out)  :: pbl_mixing_tau
    type (t_Diffusion), intent(in)   :: self
    pbl_mixing_tau = self%pbl_mixing_tau
    return
  end subroutine Get_pbl_mixing_tau
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  end module GmiDiffusionMethod_mod
