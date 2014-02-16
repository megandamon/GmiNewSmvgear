#include "GmiESMF_ErrLog.h"
!------------------------------------------------------------------------------
! NASA/GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULES: GmiControlInitialize_mod
!
      module GmiControlInitialize_mod
!
! !USES:
      use ESMF_Mod
      use GmiESMF_ErrorChecking_mod
      use GmiESMFrcFileReading_mod, only : rcEsmfReadLogical
      use GmiESMFderivedType_mod       , only : t_gmiESMF
      ! Devel Component registration routines
      use GmiAdvecCoreStubs_mod, only : advCoreSetServices => SetServices
!!!!!      use FVadvcore_GridCompMod, only : advCoreSetServices => SetServices
      use  GmiAdvecCoreWrapper_mod     , only : initializeAdvecCoreTracers

      use GmiPrintError_mod            , only : GmiPrintError
      use GmiMessagePassing_mod        , only : synchronizeGroup
      use GmiDiagnosticsMethod_mod     , only : t_Diagnostics, initializeDiagnostics
      use GmiDiagnosticsMethod_mod     , only : Get_do_ftiming, Get_pr_diag
      use GmiDomainDecomposition_mod   , only : setSubDomainData
      use GmiDomainDecomposition_mod    , only : t_gmiDomain
      use GmiDomainDecomposition_mod    , only : Get_communicatorWorld
      use GmiDomainDecomposition_mod    , only : Get_map1_u, Get_map1_v
      use GmiDomainDecomposition_mod    , only : Get_iAmRootProc, Get_procID, Get_rootProc
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration,        &
     &       InitializeSpcConcentration, constructConcentrationBundle,         &
     &       Set_num_const_inrecs, ReadSpcConcentrationResourceFile
      use GmiGrid_mod, only : t_gmiGrid, Get_numSpecies, Get_i1, Get_i2,       &
     &       Get_ju1, Get_jv1, Get_j2, Get_k1, Get_k2, Get_ilo, Get_ihi,       &
     &       Get_julo, Get_jvlo, Get_jhi, Get_i1_gl, Get_i2_gl,    &
     &       Get_j2_gl, Get_ilo_gl, Get_ihi_gl, Get_julo_gl,       &
     &       Get_jvlo_gl, Get_jhi_gl, Get_gmi_nborder
      use GmiTimeControl_mod           , only : t_GmiClock
      use GmiEmissionMethod_mod, only : t_Emission, InitializeEmission,        &
     &       readEmissionResourceFile, initReadEmission
      use GmiDiffusionMethod_mod       , only : t_Diffusion, initializeDiffusion
      use GmiChemistryMethod_mod, only : t_Chemistry, initializeChemistry,     &
     &       Get_chem_opt, initReadChemistry, ReadChemistryResourceFile,       &
     &       Get_num_chem, Get_num_molefrac, Get_do_full_chem, Get_phot_opt,   &
     &       Get_const_labels, Get_ihno3_num, Get_io3_num
      use GmiDepositionMethod_mod      , only : t_Deposition, InitializeDeposition
      use GmiConvectionMethod_mod      , only : t_Convection, InitializeConvection
      use GmiAdvectionMethod_mod, only : t_Advection, initializeAdvection, &
     &       Get_advec_opt
      use GmiControlOutput_mod, only : initializeOutputFiles
      use GmiMetFieldsControl_mod, only : t_metFields, initializeMetFields, &
     &       Get_metdata_name_model, Get_metdata_name_org, Get_met_grid_type,  &
     &       Set_Met1, Set_Met2
      use GmiControlMetInput_mod, only : Set_Restart_Partial
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: gmiControlInitialize

#     include "GmiParameters.h"

!EOP
!=============================================================================
      contains
!=============================================================================
!BOP
!
! !IROUTINE: gmiControlInitialize
!
! !INTERFACE:
!
      subroutine gmiControlInitialize (advCoreESMF, esmfClock, grid, gmiConfigFile, &
     &                SpeciesConcentration, Emission, Chemistry, Deposition, &
     &                Convection, Diffusion, Advection, Diagnostics,   &
     &                gmiClock, gmiGrid, gmiDomain, metFields, chem_mecha)
!
      implicit none

! !INPUT PARAMETERS:
      character(len=*) , intent(in) :: chem_mecha
!
! !INPUT/OUTPUT PARAMETERS:    
      type(ESMF_Grid)  , intent(inOut) :: grid
      type(ESMF_Clock) , intent(inOut) :: esmfClock
      type (t_gmiESMF)             , intent(inOut) :: advCoreESMF
      type (ESMF_Config)           , intent(inOut) :: gmiConfigFile
      type (t_Diagnostics)         , intent(inOut) :: Diagnostics
      type (t_GmiClock            ), intent(inOut) :: gmiClock
      type (t_Emission            ), intent(inOut) :: Emission
      type (t_Chemistry           ), intent(inOut) :: Chemistry
      type (t_Advection           ), intent(inOut) :: Advection
      type (t_metFields           ), intent(inOut) :: metFields
      type (t_Deposition          ), intent(inOut) :: Deposition
      type (t_Convection          ), intent(inOut) :: Convection
      type (t_Diffusion           ), intent(inOut) :: Diffusion
      type (t_gmiGrid             ), intent(inOut) :: gmiGrid
      type (t_gmiDomain           ), intent(inOut) :: gmiDomain
      type (t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
!   Initializes all the GMI components.
!
! !LOCAL VARIABLES:
      integer :: i1, i2, ju1, jv1, j2, k1, k2, ilo, ihi, julo, jvlo, jhi
      integer :: i1_gl, i2_gl, j2_gl, gmi_nborder
      integer :: ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl
      integer :: STATUS, rc
      logical :: pr_diag, do_wetdep
      integer :: chem_opt, advec_opt
      integer :: commuWorld, procID, sad_opt
      type(ESMF_VM) :: vm
      character (len=1 ) :: met_grid_type
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model
      character(len=ESMF_MAXSTR) :: IAm = "gmiControlInitialize"
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID  (gmiDomain, procID)
      call Get_pr_diag (Diagnostics, pr_diag)

      if (pr_diag) write(6,*) trim(Iam), " called by ", procID

      !-------------------------------------------
      !  Initialize the Tracer Advection Component
      !-------------------------------------------

      call ESMF_VMGetGlobal( vm, rc=STATUS  )
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute(gmiConfigFile, advec_opt, &
                label="advec_opt:", default=1, rc=STATUS )
      VERIFY_(STATUS)

      IF (advec_opt == 2) THEN
         ! Create advecCore gridded component

         advCoreESMF%compGridded = ESMF_GridCompCreate(    &
                 name   = "Advection Comp",             &
                 grid   = grid,                         &
                 config = gmiConfigFile,                &
                 rc     = STATUS)
         VERIFY_(STATUS)

         call ESMF_GridCompValidate(advCoreESMF%compGridded, rc=STATUS)
         !VERIFY_(STATUS)

         ! Create advecCore import and export states

         advCoreESMF%stateImp = ESMF_StateCreate("advecCoreImpState", &
                                  ESMF_STATE_IMPORT, rc=STATUS)
         VERIFY_(STATUS)

         advCoreESMF%stateExp = ESMF_StateCreate("advecCoreExpState", &
                                  ESMF_STATE_EXPORT, rc=STATUS)
         VERIFY_(STATUS)

         call ESMF_GridCompSetServices(advCoreESMF%compGridded, &
                         advCoreSetServices, STATUS)
         VERIFY_(STATUS)

         call ESMF_GridCompInitialize(advCoreESMF%compGridded, &
                 importState = advCoreESMF%stateImp, &
                 exportState = advCoreESMF%stateExp, &
                 clock       = esmfClock,            &
                 rc          = STATUS)
         VERIFY_(STATUS)
      END IF

      !--------------------------
      ! Initialize the components
      !--------------------------

      call initializeDiagnostics (Diagnostics, gmiConfigFile, gmiGrid, gmiDomain)

      !-----------------------------------
      ! Initialize the metFields component
      !-----------------------------------

      call initializeMetFields(metFields, gmiGrid, gmiDomain, gmiClock, &
     &                     Diagnostics, gmiConfigFile)

      call Get_met_grid_type(metFields, met_grid_type)
      call Get_metdata_name_org(metFields, metdata_name_org)
    
      call setSubDomainData(gmiDomain, gmiGrid, pr_diag, metdata_name_org, &
     &               met_grid_type)

      call Set_Met1 (metFields, gmiDomain)
      
      call Set_Restart_Partial (metFields, Diagnostics, gmiGrid, gmiDomain)
      
      call Set_Met2 (metFields, gmiDomain)

      !------------------------------
      ! Initialize all the components
      !------------------------------

      call initializeGmiComponents(gmiConfigFile, &
     &               gmiGrid, gmiDomain, gmiClock, Diagnostics, metFields,     &
     &               SpeciesConcentration, Chemistry, Emission, Deposition,    &
     &               Diffusion, Convection, Advection, TRIM(chem_mecha))

      !-------------------------------------------------
      ! Initialize the bundle in advecCore import state
      !-------------------------------------------------

      call Get_advec_opt(Advection, advec_opt)
      if (advec_opt == 2) then
         call initializeAdvecCoreTracers (advCoreESMF%stateImp, grid, Advection)
      end if

      !-----------------------------------------
      ! Initialize ASCII and netCDF output files
      !-----------------------------------------

      call initializeOutputFiles(gmiGrid, gmiDomain, Diagnostics,            &
     &               SpeciesConcentration, Chemistry, Emission, &
     &               metFields, TRIM(chem_mecha))

      return

      end subroutine gmiControlInitialize
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initializeGmiComponents
!
! !INTERFACE:
!
      subroutine initializeGmiComponents(gmiConfigFile, &
     &                gmiGrid, gmiDomain, gmiClock, Diagnostics, metFields,   &
     &                SpeciesConcentration, Chemistry, Emission, Deposition,  &
     &                Diffusion, Convection, Advection, chem_mecha)
!
      implicit none
!
! !INPUT PARAMETERS:
      character(len=*)   , intent(in) :: chem_mecha
      type(t_gmiGrid    ), intent(in) :: gmiGrid    
      type(t_gmiDomain  ), intent(in) :: gmiDomain  
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_metFields)           , intent(inOut) :: metFields  
      type (ESMF_Config)          , intent(inOut) :: gmiConfigFile
      type(t_Diagnostics)         , intent(inOut) :: Diagnostics
      type(t_GmiClock            ), intent(inOut) :: gmiClock
      type(t_Chemistry           ), intent(inOut) :: Chemistry           
      type(t_Advection           ), intent(inOut) :: Advection           
      type(t_Deposition          ), intent(inOut) :: Deposition          
      type(t_Convection          ), intent(inOut) :: Convection          
      type(t_Diffusion           ), intent(inOut) :: Diffusion           
      type(t_Emission            ), intent(inOut) :: Emission            
      type(t_SpeciesConcentration), intent(inOut) :: SpeciesConcentration
!
! !DESCRIPTION:
! Initializes all the GMI main components.
!
! !LOCAL VARIABLES:
      integer :: i1, i2, ju1, j2, i1_gl, rc, STATUS
      integer :: num_chem, num_molefrac, procID
      integer :: chem_opt, ihno3_num, io3_num, numSpecies
      logical :: pr_diag, do_full_chem, do_drydep
      character(len=MAX_LENGTH_SPECIES_NAME), pointer :: speciesNames(:)
      character(len=ESMF_MAXSTR) :: IAm = "initializeGmiComponents"
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_procID (gmiDomain, procID)
      call Get_pr_diag (Diagnostics, pr_diag)

      if (pr_diag) Write(6,*) IAm, ' called by ', procID

      call readChemistryResourceFile (Chemistry, Diagnostics, gmiGrid, &
     &         gmiDomain, gmiConfigFile)

      call ReadSpcConcentrationResourceFile(SpeciesConcentration, gmiGrid,     &
     &         gmiDomain, Diagnostics, gmiConfigFile)

      call readEmissionResourceFile(Emission, gmiGrid, gmiDomain, Diagnostics, &
     &         gmiConfigFile)


      call Get_numSpecies(gmiGrid, numSpecies)

      allocate(speciesNames(numSpecies))
      call Get_const_labels(Chemistry, speciesNames)

      call initializeAdvection(Advection, gmiGrid, gmiDomain, Diagnostics, &
     &                      gmiConfigFile, speciesNames)

      deallocate(speciesNames)

      call Get_io3_num     (Chemistry, io3_num)
      call Get_ihno3_num   (Chemistry, ihno3_num)
      call Get_chem_opt    (Chemistry, chem_opt    )
      call Get_num_chem    (Chemistry, num_chem    )
      call Get_num_molefrac(Chemistry, num_molefrac)
      call Get_do_full_chem(Chemistry, do_full_chem)

      call InitializeSpcConcentration &
     &            (SpeciesConcentration, gmiGrid, gmiDomain, Diagnostics, & 
     &             do_full_chem, num_molefrac, num_chem, TRIM(chem_mecha))

      call InitializeDeposition (Deposition, gmiGrid, gmiDomain, Diagnostics, &
     &            gmiConfigFile)

      call InitializeConvection (Convection, gmiDomain, Diagnostics, &
     &             gmiConfigFile)

      call initializeDiffusion (Diffusion, gmiDomain, Diagnostics, gmiConfigFile)

      call rcEsmfReadLogical(gmiConfigFile, do_drydep, "do_drydep:", &
     &           default=.false., rc=STATUS )
      VERIFY_(STATUS)

      call InitializeEmission (Emission, SpeciesConcentration, gmiGrid,     &
     &               gmiDomain, Diagnostics, &
     &               do_drydep, ihno3_num, io3_num)

      call initReadEmission (Emission, gmiClock, gmiGrid, gmiDomain, &
     &            Diagnostics, metFields, do_drydep)

      call InitializeChemistry (Chemistry, Emission, gmiGrid, gmiDomain, &
     &            gmiClock, Diagnostics, metFields, TRIM(chem_mecha))

      call initReadChemistry (Chemistry, gmiGrid, gmiDomain, Diagnostics,   &
     &         TRIM(chem_mecha))

      return

      end subroutine initializeGmiComponents
!EOC
!------------------------------------------------------------------------------
      end module GmiControlInitialize_mod
