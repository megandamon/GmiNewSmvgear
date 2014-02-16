!-------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!-------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiControlAdvance_mod
!
#include "GmiESMF_ErrLog.h"
!
      module GmiControlAdvance_mod
!
! !USES:
      use ESMF_Mod 
      use GmiESMFderivedType_mod       , only : t_gmiESMF
      use GmiTimeStepping_mod          , only : gmiTimeStepping
      use GmiControlOutput_mod         , only : controlOutputFiles
      use GmiControlMetInput_mod, only : Control_Met1_Input
      use GmiControlMetInput_mod, only : Control_Met2_Input
      use GmiESMFclock_mod             , only : advanceESMFclock
      use GmiTimeControl_mod           , only : t_GmiClock, Get_gmiTimeStep
      use GmiTimeControl_mod           , only : Get_curGmiDate, Get_curGmiTime
      use GmiTimeControl_mod           , only : Get_gmiSeconds, Get_numTimeSteps
      use GmiEmissionMethod_mod        , only : t_Emission
      use GmiDepositionMethod_mod      , only : t_Deposition
      use GmiConvectionMethod_mod      , only : t_Convection
      use GmiDiffusionMethod_mod       , only : t_Diffusion
      use GmiChemistryMethod_mod       , only : t_Chemistry
      use GmiAdvectionMethod_mod       , only : t_Advection
      use GmiPrintError_mod            , only : GmiPrintError
      use GmiDomainDecomposition_mod   , only : t_gmiDomain, Get_iAmRootProc
      use GmiDomainDecomposition_mod   , only : Get_procID
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiGrid_mod, only : t_gmiGrid
      use GmiGrid_mod, only : Get_i1, Get_i2, Get_ju1, Get_j2
      use GmiGrid_mod, only : Get_i1_gl, Get_ju1_gl, Get_i2_gl, Get_j2_gl

      use GmiDiagnosticsMethod_mod     , only : t_Diagnostics, Get_pr_diag
      use GmiMetFieldsControl_mod, only : t_metFields, Get_met_opt, Get_mdt, &
     &       met1Glob2Sub

      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: gmiControlAdvance

#     include "GmiParameters.h"

!-------------------------------------------------------------------------------
      contains
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gmiControlAdvance
!
! !INTERFACE:
!
      subroutine gmiControlAdvance(advCoreESMF, esmfClock, gmiConfigFile,  &
     &              SpeciesConcentration, Emission, Chemistry, Deposition, &
     &              Convection, Diffusion, Advection, Diagnostics, gmiClock, &
     &              gmiGrid, gmiDomain, metFields, chem_mecha)

      implicit none
!
! !INPUT PARAMETERS:
      character(len=*)  , intent(in) :: chem_mecha
      type (t_gmiGrid  ), intent(in) :: gmiGrid   
      type (t_gmiDomain), intent(in) :: gmiDomain           
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config)           , intent(inOut) :: gmiConfigFile
      type (t_gmiESMF)             , intent(inOut) :: advCoreESMF
      type (t_Diagnostics         ), intent(inOut) :: Diagnostics
      type (ESMF_Clock            ), intent(inOut) :: esmfClock
      type (t_GmiClock            ), intent(inOut) :: gmiClock
      type (t_Emission            ), intent(inOut) :: Emission
      type (t_Chemistry           ), intent(inOut) :: Chemistry 
      type (t_Advection           ), intent(inOut) :: Advection 
      type (t_metFields           ), intent(inOut) :: metFields 
      type (t_Deposition          ), intent(inOut) :: Deposition
      type (t_Convection          ), intent(inOut) :: Convection
      type (t_Diffusion           ), intent(inOut) :: Diffusion 
      type (t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
! Carries out the model integration (advance in time at the requested number of
! of time steps).
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg

      logical :: first_tstp
      logical :: last_tstp
      logical :: new_met_rec

      integer :: istep, rc
      integer :: totNumTimeSteps, num_time_steps
      integer :: procID, met_opt, ndt
      real*8  :: mdt, tdt
      logical :: iAmRootProc, pr_diag
      real(ESMF_KIND_R8) :: runTimeStepCount
      character(len=ESMF_MAXSTR), parameter :: IAm = "gmiControlAdvance"
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID     (gmiDomain, procID     )

      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) write(6,*) trim(Iam), " called by ", procID

      ! Determine the total number of time steps
      call ESMF_ClockGet(esmfClock, runTimeStepCount=runTimeStepCount)

      totNumTimeSteps = runTimeStepCount

      call Get_gmiTimeStep(gmiClock, tdt)
      ndt = Nint(tdt)

      call Get_mdt    (MetFields, mdt)
      call Get_met_opt(MetFields, met_opt)

      ! Loop to carry out all the time steps
      
      istep = 0
      do while (.not. ESMF_ClockIsStopTime(esmfClock, rc))
         istep = istep + 1
         first_tstp = (istep == 1)
         last_tstp  = (istep == totNumTimeSteps)

         if (met_opt /= 1) then
            call Get_numTimeSteps(gmiClock, num_time_steps)

            if (first_tstp) then
               new_met_rec = .true.
            else if (Mod (num_time_steps, Nint (mdt) / ndt) == 0) then
               new_met_rec = .true.
            else
               new_met_rec = .false.
            end if
         end if

         call advanceOneTimeStep(advCoreESMF, esmfClock, gmiConfigFile, &
     &               gmiGrid, gmiClock, gmiDomain, Advection, Chemistry, &
     &               Emission, Deposition, Convection, Diffusion, &
     &               SpeciesConcentration, Diagnostics, metFields,     &
     &               first_tstp, last_tstp, new_met_rec, TRIM(chem_mecha))
      end do

      return

      end subroutine gmiControlAdvance
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: advanceOneTimeStep
!
! !INTERFACE:
!
      subroutine advanceOneTimeStep(advCoreESMF, esmfClock, gmiConfigFile, &
     &               gmiGrid, gmiClock, gmiDomain, Advection, Chemistry, &
     &               Emission, Deposition, Convection, Diffusion,              &
     &               SpeciesConcentration, Diagnostics, metFields,             &
     &               first_tstp, last_tstp, new_met_rec, chem_mecha)
!
      implicit none
!
! !INPUT PARAMETERS:
      character(len=*)  , intent(in) :: chem_mecha
      logical           , intent(in) :: first_tstp, last_tstp, new_met_rec
      type (t_gmiDomain), intent(in) :: gmiDomain 
      type (t_gmiGrid  ), intent(in) :: gmiGrid 
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_Config)           , intent(inOut) :: gmiConfigFile
      type (t_gmiESMF)             , intent(inOut) :: advCoreESMF
      type (t_Diagnostics         ), intent(inOut) :: Diagnostics 
      type (ESMF_Clock            ), intent(inOut) :: esmfClock
      type (t_GmiClock            ), intent(inOut) :: gmiClock
      type (t_Emission            ), intent(inOut) :: Emission
      type (t_Chemistry           ), intent(inOut) :: Chemistry
      type (t_Advection           ), intent(inOut) :: Advection
      type (t_metFields           ), intent(inOut) :: metFields
      type (t_Deposition          ), intent(inOut) :: Deposition
      type (t_Convection          ), intent(inOut) :: Convection
      type (t_Diffusion           ), intent(inOut) :: Diffusion 
      type (t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration 
!
! !DESCRIPTION:
!  Prepares the data to carry out one time step.
!
! !LOCAL VARIABLES:
      logical              :: loc_last_tstp, iAmRootProc
      logical              :: pr_diag
      integer              :: procID
      integer              :: i1, i2, ju1, j2, i1_gl, ju1_gl, j2_gl, i2_gl
      integer, save        :: oldDay = -999
      integer              :: curDay, rc, met_opt
      real*8 , allocatable :: surfTemp(:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID      (gmiDomain, procID)
      call Get_pr_diag(Diagnostics, pr_diag)

      if (pr_diag) Write (6,*) 'advanceOneTimeStep called by ', procID

!      print *,'sds-advanceOneTimeStep 00: ',procID


      call Get_iAmRootProc (gmiDomain, iAmRootProc)

      call Get_i1         (gmiGrid, i1         )
      call Get_i2         (gmiGrid, i2         )
      call Get_ju1        (gmiGrid, ju1        )
      call Get_j2         (gmiGrid, j2         )
      call Get_i1_gl      (gmiGrid, i1_gl      )
      call Get_i2_gl      (gmiGrid, i2_gl      )
      call Get_ju1_gl     (gmiGrid, ju1_gl     )
      call Get_j2_gl      (gmiGrid, j2_gl      )

      call Get_met_opt(MetFields, met_opt)

!      print *,'sds-advanceOneTimeStep 01: ',procID


      if (met_opt /= 1) then

!        =======================
         call Control_Met1_Input    & ! pressure fixer, etc.
!        =======================
     &          (metFields, Advection, gmiDomain, &
     &           Diagnostics, first_tstp, new_met_rec)

      end if

!      print *,'sds-advanceOneTimeStep 02: ',procID


      !--------------------------------------------------------------
      ! So that the rootProc can spend its idle time doing
      ! Control_Met1_Input while waiting for the workerProcs to complete
      ! their time step, call controlOutputFiles here, except for the
      ! first time through the time step loop when no output info is
      ! available.  Otherwise if controlOutputFiles were always called at
      ! the end of the time step loop, the rootProc would potentially
      ! have to wait for output messages from the workerProcs and then the
      ! workerProcs would have to wait for the rootProc to complete
      ! Control_Met1_Input.
      !
      ! NOTE ALSO THAT loc_last_tstp MUST ALWAYS BE SET TO FALSE HERE
      ! AS IT IS REALLY INFO FROM THE PREVIOUS TIME STEP THAT IS BEING
      ! WRITTEN OUT.
      ! --------------------------------------------------------------

      if (.not. first_tstp) then
         loc_last_tstp = .false.

!        ===================
         call controlOutputFiles  &
!        ===================
     &       (loc_last_tstp, Chemistry, Deposition, Emission,            &
     &        SpeciesConcentration, Diagnostics, gmiDomain, gmiGrid,     &
     &        gmiClock, metFields, Advection, TRIM(chem_mecha))
      end if

      if (met_opt /= 1) then
!        =====================
         call met1Glob2Sub (metFields)
!        =====================
      end if

      !-------------------------------
      ! Call the time stepping routine
      !-------------------------------

      call Control_Met2_Input(metFields, Chemistry, gmiDomain, gmiGrid, &
     &             Diagnostics, gmiClock, first_tstp, new_met_rec)

      call gmiTimeStepping (advCoreESMF, esmfClock, gmiClock, &
     &              gmiGrid, gmiDomain, Chemistry, Emission, Deposition,    &
     &              Convection, Diffusion, Advection, SpeciesConcentration, &
     &              Diagnostics, metFields, new_met_rec, TRIM(chem_mecha))

      !---------------------------------------------------
      ! Do controlOutputFiles while running on one processor
      ! only. Otherwise do it on the last time step only.
      !---------------------------------------------------

      if (last_tstp) then
!        ===================
         call controlOutputFiles  &
!        ===================
     &          (last_tstp, Chemistry, Deposition, Emission, &
     &           SpeciesConcentration, Diagnostics, gmiDomain, gmiGrid, &
     &           gmiClock, metFields, Advection, TRIM(chem_mecha))
      end if

      return

      end subroutine advanceOneTimeStep
!EOC
!------------------------------------------------------------------------------

      end module GmiControlAdvance_mod
