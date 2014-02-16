#define I_AM_MAIN

#include "GmiESMF_ErrLog.h"
!------------------------------------------------------------------------------
! NASA/GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !PROGRAM: GmiMain
!
      program GmiMain
!
! !USES:
      use Ftiming_Dao
      use GmiGrid_mod                  , only : t_gmiGrid
      use GmiFileUnit_mod              , only : InitializeFileUnitNumbers
      use GmiTimeControl_mod           , only : t_GmiClock
      use GmiEmissionMethod_mod        , only : t_Emission  
      use GmiControlAdvance_mod        , only : gmiControlAdvance
      use GmiDiffusionMethod_mod       , only : t_Diffusion
      use GmiChemistryMethod_mod       , only : t_Chemistry
      use GmiAdvectionMethod_mod       , only : t_Advection
      use GmiControlFinalize_mod       , only : gmiControlFinalize
      use GmiDepositionMethod_mod      , only : t_Deposition
      use GmiConvectionMethod_mod      , only : t_Convection
      use GmiControlInitialize_mod     , only : gmiControlInitialize
      use GmiDiagnosticsMethod_mod     , only : t_Diagnostics
      use GmiDomainDecomposition_mod   , only : t_gmiDomain, domainDecomposition
      use GmiDomainDecomposition_mod   , only : t_gmiDomain, Set_communicatorWorld
      use GmiDomainDecomposition_mod   , only : Set_iAmRootProc, Set_procID, Set_rootProc
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration
      use GmiMetFieldsControl_mod      , only : t_metFields
      use GmiESMFrcFileReading_mod, only : rcEsmfReadLogical
      use GmiSpeciesRegistry_mod  , only : set_labelsSpecies, set_molWeightSpecies

      use GmiESMFderivedType_mod       , only : t_gmiESMF
      use GmiESMFclock_mod             , only : createESMFclock
      use GmiESMFgrid_mod              , only : createESMFgrid

      ! ESMF module, defines all ESMF data types and procedures
      use ESMF_Mod
      use GmiESMF_ErrorChecking_mod

      implicit none

#     include "GmiParameters.h"

      type (t_gmiGrid             ) :: gmiGrid
      type (t_GmiClock            ) :: gmiClock
      type (t_Emission            ) :: Emission
      type (t_Chemistry           ) :: Chemistry
      type (t_Advection           ) :: Advection
      type (t_Diffusion           ) :: Diffusion
      type (t_gmiDomain           ) :: gmiDomain
      type (t_metFields           ) :: metFields 
      type (t_Deposition          ) :: Deposition
      type (t_Convection          ) :: Convection
      type (t_Diagnostics         ) :: Diagnostics
      type (t_SpeciesConcentration) :: SpeciesConcentration

      type(ESMF_VM)     :: vm
      type(ESMF_Config) :: gmiConfigFile
      type(ESMF_Grid)   :: grid          ! A common grid
      type(ESMF_Clock)  :: esmfClock     ! A clock
      type(t_gmiESMF)   :: advCoreESMF   ! for advCore component

      integer             :: STATUS, rc
      logical             :: do_ftiming
      logical             :: iAmRootProc
      integer             :: commuWorld, procID, numProcessors, ierr, rootProc
      CHARACTER (len=GMI_MAXSTR) :: gmiNamelistFile
      CHARACTER (len=GMI_MAXSTR) :: gmiResourceFile
      character(len=ESMF_MAXSTR), parameter :: IAm = "gmiMainProgram"

      CHARACTER (len=GMI_MAXSTR) :: chem_mecha
      EXTERNAL GETENV
!
! !AUTHOR:
!  Jules Kouatchou, NASA/GSFC, Jules.Kouatchou@nasa.gov
!EOP
!------------------------------------------------------------------------------

      !--------------------------------------------
      ! Initialize ESMF, get the default Global VM.
      !--------------------------------------------

      !call ESMF_Initialize(vm=vm, rc=STATUS)
      call ESMF_Initialize(vm=vm, defaultLogType=ESMF_LOG_NONE, rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_VMGetGlobal( vm, rc=STATUS  )
      VERIFY_(STATUS)

      !------------------------------------------------------------
      ! Get the Virtual machine, processor ID, number of processors
      ! and the message passing communicator
      !------------------------------------------------------------

      call ESMF_VMGet(vm, localpet=procID, PetCount=numProcessors, &
     &                MPICommunicator=commuWorld, rc=STATUS )
      VERIFY_(STATUS)

      iAmRootProc = (procID == 0)

      rootProc = 0
    
      ! Populate the gmiDomain derived type
    
      call Set_procID           (gmiDomain, procID     )
      call Set_rootProc         (gmiDomain, rootProc   )
      call Set_iAmRootProc      (gmiDomain, iamRootProc)
      call Set_communicatorWorld(gmiDomain, commuWorld )

      !----------------------------
      ! Print out some Announcement
      !----------------------------

      if (iAmRootProc) then
         Write(6,*) "---------------------------------------------------------------"
         Write(6,*) "           Global Modeling Initiative (GMI) Code               "
         Write(6,*) "---------------------------------------------------------------"
         Write(6,*) "Initially developed at Lawrence Livermore National Laboratory."
         Write(6,*) "Expended and maintained at NASA Goddard Space Flight Center."
         Write(6,*) "Additional information can be obtained at:"
         Write(6,*) "     - https://gmi.gsfc.nasa.gov"
         Write(6,*) "     - https://modelingguru.nasa.gov"
         Write(6,*) "---------------------------------------------------------------"
      end if

      call ESMF_VMBarrier(vm, rc=STATUS)
      VERIFY_(STATUS)

      !---------------------------
      ! Create the ESMF configFile
      !---------------------------

      gmiResourceFile = 'gmiResourceFile.rc'

      gmiConfigFile = ESMF_ConfigCreate(rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_ConfigLoadFile(gmiConfigFile, gmiResourceFile, rc=STATUS )
      VERIFY_(STATUS)

      call rcEsmfReadLogical(gmiConfigFile, do_ftiming, "do_ftiming:", &
     &                       default=.false., rc=STATUS)
      VERIFY_(STATUS)

      if (do_ftiming) then
         call Ftiming_Init ( )
     
         call Ftiming_On ('whole_GMI')
      end if

      !----------------------------------------------------
      ! Create and initialize a ESMF clock and the gmiClock
      !----------------------------------------------------

      call createESMFclock(gmiConfigFile, esmfClock, gmiClock, STATUS)

      !---------------------------------
      ! Perform the Domain decomposition
      !---------------------------------

      call domainDecomposition(gmiDomain, gmiGrid, numProcessors, gmiConfigFile)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! This section of the code allows us to get species information
      ! (number of species, species labels, species molecular weights, etc.)
      ! regardless of the type of experiments (tracer runs or not) we do.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call set_labelsSpecies   (gmiConfigFile)  ! set the species names
      call set_molWeightSpecies(gmiConfigFile)  ! set the species molecular weights

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !--------------------------------
      ! Get the chemical mechanism name
      !--------------------------------

      CALL GETENV('CHEMCASE', chem_mecha)

      call InitializeFileUnitNumbers ()

      !-------------------
      ! Create a ESMF grid
      !-------------------

      grid = createESMFgrid(gmiConfigFile, vm, rc=STATUS)
      VERIFY_(STATUS)

      !--------------------------
      ! Initialize the components
      !--------------------------

      call gmiControlInitialize(advCoreESMF, esmfClock, grid, gmiConfigFile, &
     &               SpeciesConcentration, Emission, Chemistry, Deposition, &
     &               Convection, Diffusion, Advection, Diagnostics, gmiClock, &
     &               gmiGrid, gmiDomain, metFields, TRIM(chem_mecha))

      call ESMF_VMBarrier(vm, rc=STATUS)
      VERIFY_(STATUS)

      !--------------------------
      ! Advance the model in time
      !--------------------------
     
      call gmiControlAdvance (advCoreESMF, esmfClock, gmiConfigFile, &
     &               SpeciesConcentration, Emission, Chemistry, Deposition, &
     &               Convection, Diffusion, Advection, Diagnostics, gmiClock, &
     &               gmiGrid, gmiDomain, metFields, TRIM(chem_mecha))

      !------------------------
      ! Finalize the components
      !------------------------

      call gmiControlFinalize (advCoreESMF, esmfClock, gmiDomain, Diagnostics)

      call ESMF_VMBarrier(vm, rc=STATUS)

      if (iAmRootProc) then
         Write(6,*) "---------------------------------------------------------------"
         Write(6,*) " --------  Successful completion of the run  -----------       "
         Write(6,*) "---------------------------------------------------------------"
      end if

      ! Finalize ESMF
 
      call ESMF_Finalize(rc=STATUS)

  end program GmiMain
