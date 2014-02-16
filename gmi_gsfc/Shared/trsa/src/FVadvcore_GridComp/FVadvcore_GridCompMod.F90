! $Id: FVadvcore_GridCompMod.F90,v 1.3 2011-08-09 22:12:59 mrdamon Exp $
#include "MAPL_Generic.h"


!
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: FVadvcore_GridCompMod
!
! !DESCRIPTION: 
!    This a MAPL component that can be used in
!    either with offline or online applications to advect an arbitrary set
!    of constituents.
!
! \paragraph{Scientific Description:}
!
!   The advection scheme used is that frome the FVdycore grid-point
!   dynamical core.  It runs on a sphere and uses finite-volume
!   discretization techniques. The advection is time split into a
!   horizontal phase that is assumed to be vertically Lagrangian and a
!   vertical remap phase. A complete description of the core from
!   which this component is taken may be found in:
!
!   \begin{quote}
!   Lin, S.-J. 2004, A vertically Lagrangian Finite-Volume Dynamical 
!   Core for Global Models. {\em Mon. Wea. Rev.}, {\bf 132}, 2293-2307.
!   \end{quote}
!
!  \paragraph{Code Implementation:}
!
!    It code uses the MAPL (http://MAPLCode.org/maplwiki/) to
!    encapsulate the FV advection scheme as an ESMF gridded component
!    using the ESMF paradigm of initialize, run and finalize methods,
!    and their SetServices routine. As in all ESMF codes, only
!    SetServices is public and the interface consists of of a Clock
!    and Import and Export states.  The import state includes a
!    specialized description of the motion field in terms of C-grid
!    winds and mass fluxes. These are assumed to have been accumulated
!    over the time interval specified in the resource file. The
!    default of this interval is 1800 seconds. The layer pressure
!    thicknesses in the import state are assumed to be the
!    instantaneous values valid at the beginning of this interval.  If
!    these thicknesses are friendly they will be updated to values
!    valid at the end of the interval, consistent with the given
!    motion field.  Mixing ratios of the constituents to be advected
!    are placed ESMF Fields within an ESMF Bundle in the Import
!    state. Each Field in the Bundle is tested for ``Friendliness'' to
!    advection; if friendly it is advected and its values updated.
!
!    Currently no Export capability is implemented. 
!
! !INTERFACE:

module FVadvcore_GridCompMod

! !USES:

  use ESMF_Mod
  use MAPL_Mod

  use TracerMod,    only : T_TRACERS, Tracer_Get, Tracer_Set
  use FVAdvMod,     only : T_FVAdv_STATE, FVAdvReal, &
                           FVAdv_Init,  FVAdv_Run, FVAdv_Final, &
                           FVAdv_Remap, FVAdv_PrepVars
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

!------------------------------------------------------------------------------

  type T_FVAdv_STATE_WRAP
     type (T_FVAdv_STATE), pointer :: state
  end type T_FVAdv_STATE_WRAP

  integer, parameter :: r4        = ESMF_KIND_R4
  integer, parameter :: r8        = ESMF_KIND_R8

contains
        
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!BOPI
! !IROUTINE: SetServices - Externally visible registration routine

! !INTERFACE:

  subroutine SetServices(GC, rc)
!
! !ARGUMENTS:
    type(ESMF_GridComp), intent(inout) :: GC
    integer, optional,   intent(  out) :: RC
!
! !DESCRIPTION:
!
!     User-supplied setservices routine.
!     The register routine sets the subroutines to be called
!     as the init, run, and finalize routines.  Note that those are
!     private to the module.
!
!EOPI

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------
    Iam = 'SetServices'

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

!BOP

! !PARAMETERS:
!
!  (IMPORT STATE)

    call MAPL_AddImportSpec(GC,                                  &
       LONG_NAME          = 'eastward_wind_component_on_Cgrid',  &
       UNITS              = 'm s-1',                             &
       SHORT_NAME         = 'UC',                                &
       DIMS               = MAPL_DIMSHORZVERT,                   &
       VLOCATION          = MAPL_VLocationCenter,                &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                  &
       LONG_NAME          = 'northward_wind_component_on_Cgrid', &
       UNITS              = 'm s-1',                             &
       SHORT_NAME         = 'VC',                                &
       DIMS               = MAPL_DIMSHORZVERT,                   &
       VLOCATION          = MAPL_VLocationCenter,                &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                  &
       LONG_NAME          = 'unremapped_eastward_pressure_flux_on_Cgrid',   &
       UNITS              = 'm+2 s-1 Pa',                        &
       SHORT_NAME         = 'MX_UR',                             &
       DIMS               = MAPL_DIMSHORZVERT,                   &
       VLOCATION          = MAPL_VLocationCenter,                &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                  &
       LONG_NAME          = 'unremapped_northward_pressure_flux_on_Cgrid',  &
       UNITS              = 'm+2 s-1 Pa',                        &
       SHORT_NAME         = 'MY_UR',                             &
       DIMS               = MAPL_DIMSHORZVERT,                   &
       VLOCATION          = MAPL_VLocationCenter,                &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                  &
       LONG_NAME          = 'eastward_pressure_flux_on_Cgrid',   &
       UNITS              = 'm+2 s-1 Pa',                        &
       SHORT_NAME         = 'MX',                                &
       DIMS               = MAPL_DIMSHORZVERT,                   &
       VLOCATION          = MAPL_VLocationCenter,                &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                  &
       LONG_NAME          = 'northward_pressure_flux_on_Cgrid',  &
       UNITS              = 'm+2 s-1 Pa',                        &
       SHORT_NAME         = 'MY',                                &
       DIMS               = MAPL_DIMSHORZVERT,                   &
       VLOCATION          = MAPL_VLocationCenter,                &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                  &
       LONG_NAME          = 'upward_mass_flux_due_to_dynamics',  &
       UNITS              = 'kg m-2 s-1',                        &
       SHORT_NAME         = 'MZ',                                &
       DIMS               = MAPL_DIMSHORZVERT,                   &
       VLOCATION          = MAPL_VLocationEdge,                  &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                  &
       LONG_NAME          = 'initial_pressure_thickness',        &
       UNITS              = 'Pa',                                &
       SHORT_NAME         = 'DP',                                &
       DIMS               = MAPL_DIMSHORZVERT,                   &
       VLOCATION          = MAPL_VLocationCenter,                &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

!
! DPNEW currently not used!   PE2 calculated from MX, MY, MZ, DP
!
    call MAPL_AddImportSpec(GC,                                  &
       LONG_NAME          = 'final_pressure_thickness',          &
       UNITS              = 'Pa',                                &
       SHORT_NAME         = 'DPNEW',                             &
       DIMS               = MAPL_DIMSHORZVERT,                   &
       VLOCATION          = MAPL_VLocationCenter,                &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

!
! DPEDT currently not used!
!
    call MAPL_AddImportSpec(GC,                                  &
       LONG_NAME          = 'air_pressure_tendency_due_to_nondynamical_processes', &
       UNITS              = 'Pa s-1',                            &
       SHORT_NAME         = 'DPEDT',                             &
       DIMS               = MAPL_DIMSHORZVERT,                   &
       VLOCATION          = MAPL_VLocationEdge,                  &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                  &
       SHORT_NAME         = 'TRACERS',                           &
       LONG_NAME          = 'advected_quantities',               &
       units              = 'X',                                 &
       DIMS               = MAPL_DIMSHORZVERT,                   &
       VLOCATION          = MAPL_VLocationCenter,                &
       DATATYPE           = MAPL_BundleItem,                     &
                                                      RC=STATUS  )
     VERIFY_(STATUS)

   call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'UC',                                        &
         LONG_NAME  = 'eastward_wind',                            &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,          &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'VC',                                        &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,          &
         RC=STATUS  )
     VERIFY_(STATUS)

!EOP

! Register methods with MAPL
! --------------------------

    call MAPL_GridCompSetEntryPoint(GC, ESMF_SETINIT, Initialize, &
                                    STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint(GC, ESMF_SETRUN, Run,         &
                                    STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint(GC, ESMF_SETFINAL, Finalize,  &
                                    STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint(GC, "ESMF_ReadRestart", TestCase, &
                                    STATUS )
    VERIFY_(STATUS)

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="RUN"    ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-HORZ"  ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-VERT"  ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="INIT"   ,RC=STATUS); VERIFY_(STATUS)

! Ending with a Generic SetServices call is a MAPL requirement 
!-------------------------------------------------------------
    call MAPL_GenericSetServices    ( GC, rc=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine SetServices

!------------------------------------------------------------------------------

!BOPI
! !IROUTINE: Initialize - initialization routine

! !INTERFACE:

  subroutine Initialize(GC, IMPORT, EXPORT, CLOCK, RC)

! !ARGUMENTS:

     type(ESMF_GridComp), intent(inout) :: GC
     type(ESMF_State),    intent(inout) :: IMPORT, EXPORT
     type(ESMF_Clock),    intent(inout) :: CLOCK
     integer, optional,   intent(  out) :: RC
!
! !DESCRIPTION:
!     This initialization routine creates the import and export states,
!     as well as the internal state, which is attached to the component.
!     It also determines the distribution (and therefore the grid) 
!     and performs allocations of persistent data, 
!
!EOPI

    type (T_FVAdv_STATE), pointer      :: state
    type (T_FVAdv_STATE_WRAP)          :: wrap
    type(ESMF_Grid)                    :: grid
    type(ESMF_VM)                      :: vm
    type(MAPL_MetaComp), pointer       :: MAPL
 
    integer                            :: globalCounts(3)
    integer, allocatable               :: cellCounts(:,:)
    integer                            :: lpet
    integer                            :: im, jm, km
    integer                            :: nx, ny
    integer, allocatable               :: ims(:), jms(:)

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME
    integer                            :: comm_esmf

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'Initialize'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Call Generic Initialize 
! -----------------------

    call MAPL_GenericInitialize(GC, IMPORT, EXPORT, CLOCK, RC=STATUS)
    VERIFY_(STATUS)

! Get my instance of the MAPL object
!-----------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,'TOTAL')
    call MAPL_TimerOn(MAPL, 'INIT')

! Allocate the private state, wrap it, and put it in the GC
!----------------------------------------------------------

    allocate( STATE, stat=STATUS )
    VERIFY_(STATUS)

    wrap%state => state

    call ESMF_UserCompSetInternalState ( GC, 'FVAdvState', wrap, STATUS )
    VERIFY_(STATUS)

! Get info from ESMF Gridded Component
!-------------------------------------

    call MAPL_Get(MAPL, nx=nx, ny=ny, rc=status)
    VERIFY_(STATUS)

    call ESMF_GridCompGet(GC, vm=vm, grid=grid, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_VMGet      (vm, localPet=lpet,    rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_GridGet    ( grid,                                      &
                           horzRelloc            = ESMF_CELL_CENTER,  &      
                           Vertrelloc            = ESMF_CELL_CELL,    &      
                           globalCellCountPerDim = globalCounts,      &      
                                                            rc=STATUS )
    VERIFY_(STATUS)

    im  = globalCounts(1)
    jm  = globalCounts(2)
    km  = globalCounts(3)

    allocate(cellCounts(nx*ny,3), stat=STATUS); VERIFY_(STATUS)

! BEWARE:  cellCounts does not translate immediately to IMXY, JMXY
    call ESMF_GridGet    ( grid,                                      &
                           horzRelloc            = ESMF_CELL_CENTER,  &      
                           Vertrelloc            = ESMF_CELL_CELL,    &      
                           cellCountPerDEPerDim  = cellCounts,        &
                                                            rc=STATUS )
    VERIFY_(STATUS)

    allocate(ims(nx), jms(ny))
    ims = cellCounts(1:NX   ,1)
    jms = cellCounts(1:  :NX,2)

    ASSERT_(sum(ims)==im)
    ASSERT_(sum(jms)==jm)

! Initialize State. This also initializes pilgrim
!------------------------------------------------

    call ESMF_VMGet(vm,  MPICommunicator=comm_esmf, rc=STATUS)
    VERIFY_(STATUS)

    call FVAdv_Init     (STATE,   &
         IM      = IM,            &
         JM      = JM,            &
         KM      = KM,            &
         imxy    = ims,           &
         jmxy    = jms,           &
         lpet    = lpet,          &
         comm    = comm_esmf,     &
                        RC=STATUS )
    VERIFY_(STATUS)

! All Done
!---------

    deallocate(cellCounts)
    deallocate(ims, jms)

    call MAPL_TimerOff(MAPL, 'INIT')
    call MAPL_TimerOff(MAPL,'TOTAL')

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOPI
! !IROUTINE: Run - run routine

! !INTERFACE:

  subroutine Run(GC, IMPORT, EXPORT, CLOCK, RC)

!
! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC
    type(ESMF_State),    intent(inout) :: IMPORT, EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC
!
! !DESCRIPTION:
! 
! The Run method advanced the advection one long time step, as
! specified in the configuration.  This may be broken down int a
! number of internal, small steps, also configurable.
!
! Each time step is split into a horizontal advection calculation,
! which is vertically Lagrangian, and a vertical remapping.
!
! The API for remap requires the edge pressures before and after the
! horizontal transport (Lagrangian) step. These are computed at the
! FVAdv precision.
!
! The calculation of PE1 (new pressures on Lagrangian frame) is
! straightforward. The desired final pressures in the Eulerian frame,
! PE2, comes from the continuity equation for Eulerian coordinates
! written for the layer thicknesses $\delta_k dP$:
!
! %
!   \begin{equation*}
!     \frac{\partial (\delta_k P)}{\partial t} = 
!       - \frac{1}{\Delta x \Delta y} (\delta_x (M_x) + \delta_y (M_y)) - g\;\delta_z (M_z)
!   \end{equation*}
! %
! This is time split between Lagrangian and vertical remapping processes:
!   \begin{eqnarray*}
!          \left(\frac{\partial (\delta_k P)}{\partial t}\right) & = & 
!          \left(\frac{\partial (\delta_k P)}{\partial t}\right)_{L} + 
!          \left(\frac{\partial (\delta_k P)}{\partial t}\right)_{R}, \\ 
!  & & \\
!     \left(\frac{\partial (\delta_k P)}{\partial t}\right)_{L} & = &
!       - \frac{1}{\Delta x \Delta y}(\delta_x (M_x) + \delta_y (M_y)),  \\
!     \left(\frac{\partial (\delta_k P)}{\partial t}\right)_{R} & = &
!       - g\;\delta_z (M_z).
!   \end{eqnarray*}
! %
!  In finite difference form, the two time steps are:
!   \begin{eqnarray*}
!      (\delta_k P)^* & = & (\delta_k P)^{(old)}
!            -\frac{\Delta t}{ \Delta x\, \Delta y}(\delta_x (M_x)+\delta_y (M_y)), \\
!      (\delta_k P)^{new} & = & (\delta_k P)^* - \Delta t g\;\delta_z (M_z).
!   \end{eqnarray*}
! %
! The first update, resulting in $(\delta_k P)^*$ is done in the
! internal FVadvRun routine, which implements the Lagrangian
! processes. The PE1 are the edge pressures corresponding to these
! thicknesses. The second update is done below, while computing the
! final thicknesses, $(\delta_k P)^{new}$. The PE2 are the edge
! pressures corresponding to these thicknesses. The vertical remap of
! the tracers is between PE1 and PE2.
!
!EOPI

! Local variables

    type(T_FVAdv_STATE), pointer :: state
    type(T_FVAdv_STATE_WRAP)     :: wrap
    type(T_TRACERS),     pointer :: tr(:)

    type(ESMF_Bundle)      :: TRACERS
    type(ESMF_Field)       :: Field
    type(ESMF_Array)       :: Array
    type(ESMF_Logical)     :: FRIENDLY
    type(ESMF_DataKind)    :: Kind

    type(MAPL_MetaComp), pointer       :: MAPL


    integer                :: i, k, n, im, jm, lm, ntracers
    integer                :: iord, jord, kord, ns
    real                   :: ptop
    real                   :: dt
    real(FVAdvReal)        :: dt_long
    real(FVAdvReal)        :: dt_short

! Pointers to imports

    real,            pointer     :: uc (:,:,:)
    real,            pointer     :: vc (:,:,:)
    real,            pointer     :: mx_ur (:,:,:)
    real,            pointer     :: my_ur (:,:,:)
    real,            pointer     :: mx (:,:,:)
    real,            pointer     :: my (:,:,:)
    real,            pointer     :: mz (:,:,:)
    real,            pointer     :: dp (:,:,:)
    real,            pointer     :: dpnew (:,:,:)
    real,            pointer     :: dpedt (:,:,:)

! Temporary pointer to tracer
    real(r4),        pointer     :: tracer_r4 (:,:,:)
    real(r8),        pointer     :: tracer_r8 (:,:,:)

! Scratch Space at FVAdv precision

    real(FVAdvReal), allocatable :: pe1(:,:,:)
    real(FVAdvReal), allocatable :: pe2(:,:,:)
    real(FVAdvReal), allocatable :: cx (:,:,:)
    real(FVAdvReal), allocatable :: cy (:,:,:)
    real(FVAdvReal), allocatable :: mfx(:,:,:)
    real(FVAdvReal), allocatable :: mfy(:,:,:)
    real(FVAdvReal), allocatable :: dp1(:,:,:)

    character(len=ESMF_MAXSTR)   :: IAm
    integer                      :: STATUS
    character(len=ESMF_MAXSTR)   :: COMP_NAME

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------


    Iam = 'Run'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get MAPL object
!----------------

    call MAPL_GetObjectFromGC(GC,MAPL,rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,'TOTAL')
    call MAPL_TimerOn(MAPL,'RUN'  )

    call MAPL_Get(MAPL, im=im, jm=jm, lm=lm, rc=status)
    VERIFY_(STATUS)

! Query configuration for run parameters
!---------------------------------------


!BOR
! !RESOURCE_ITEM: none :: Number of internal steps per call
    call MAPL_GetResource(MAPL, NS  , 'nsplit:', default=1,     rc=STATUS)
    VERIFY_(STATUS)
! !RESOURCE_ITEM: seconds :: Total time interval for each call
    call MAPL_GetResource(MAPL, DT  , 'dt:'    , default=1800., rc=STATUS)
    VERIFY_(STATUS)
! !RESOURCE_ITEM: none :: Order of zonal advection scheme
    call MAPL_GetResource(MAPL, IORD, 'iord:'  , default=3,     rc=STATUS)
    VERIFY_(STATUS)
! !RESOURCE_ITEM: none :: Order of meridional advection scheme
    call MAPL_GetResource(MAPL, JORD, 'jord:'  , default=3,     rc=STATUS)
    VERIFY_(STATUS)
! !RESOURCE_ITEM: none :: Order of vertical advection scheme
    call MAPL_GetResource(MAPL, KORD, 'kord:'  , default=4,     rc=STATUS)
    VERIFY_(STATUS)
! !RESOURCE_ITEM: Pa :: Constant pressure at the model top
    call MAPL_GetResource(MAPL, PTOP, 'ptop:'  , default=1.0,   rc=STATUS)
    VERIFY_(STATUS)
!EOR


! Get FVAdv's private internal state from the GC
!-----------------------------------------------

    call ESMF_UserCompGetInternalState( GC, 'FVAdvState', wrap, STATUS )
    VERIFY_(STATUS)

    state => wrap%state

! Get the imports. These are all native precision
!------------------------------------------------

    call MAPL_GetPointer(IMPORT,   UC,    "UC",    rc=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   VC,    "VC",    rc=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   MX_UR, "MX_UR", rc=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   MY_UR, "MY_UR", rc=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   MX,    "MX",    rc=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   MY,    "MY",    rc=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   MZ,    "MZ",    rc=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,DPEDT,"DPEDT", rc=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   DP,   "DP", rc=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   DPNEW,"DPNEW", rc=STATUS); VERIFY_(STATUS)

! The quantities to be advected come as friendlies in a bundle
!  in the import state.
!--------------------------------------------------------------

    call ESMF_StateGetBundle(IMPORT, "TRACERS", TRACERS, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_BundleGet(TRACERS, fieldCount=N       ,    rc=STATUS)
    VERIFY_(STATUS)

! We allocate a list of tracers big enough to hold all items in the bundle
!-------------------------------------------------------------------------

    allocate( TR(N),stat=STATUS )
    VERIFY_(STATUS)

! Go through the bundle copying the friendlies into the tracer list.
!  These can be of either precision
!-------------------------------------------------------------------------

    NTRACERS = 0

    do i=1, N
       call ESMF_BundleGetField (TRACERS, I, FIELD, RC=STATUS)
       VERIFY_(STATUS)

! WS 2007-08-24: FRIENDLY code not working: consult Max Suarez
!!!    FRIENDLY = MAPL_VerifyFriendly(FIELD,"ADVECTION",RC=STATUS)
!!!    VERIFY_(STATUS)

!!!    if(FRIENDLY==ESMF_true) then
          NTRACERS = NTRACERS + 1
          call ESMF_FieldGetArray  (FIELD, ARRAY,      RC=STATUS)
          VERIFY_(STATUS)
          call ESMF_ArrayGet       (ARRAY, KIND=Kind,  RC=status)
          VERIFY_(STATUS)
          if ( Kind == ESMF_R4 ) then
             call ESMF_ArrayGetData(array,tracer_r4, rc=status )
             VERIFY_(STATUS)
             call Tracer_Set(tr(NTRACERS),tracer_r4)
          elseif ( Kind == ESMF_R8 ) then
             call ESMF_ArrayGetData(array,tracer_r8, rc=status )
             VERIFY_(STATUS)
             call Tracer_Set(tr(NTRACERS),tracer_r8)
          else
             ASSERT_(.false.)
          endif
!!!    else
!!!      print *, "Tracer not FRIENDLY"
!!!    endif
    enddo

! Temporary space for mass and motion field description
!-----------------------------------------------------

    allocate(CX (IM,JM,1:LM), STAT=STATUS); VERIFY_(STATUS)
    allocate(CY (IM,JM,1:LM), STAT=STATUS); VERIFY_(STATUS)
    allocate(MFX(IM,JM,1:LM), STAT=STATUS); VERIFY_(STATUS)
    allocate(MFY(IM,JM,1:LM), STAT=STATUS); VERIFY_(STATUS)
    allocate(PE1(IM,JM,0:LM), STAT=STATUS); VERIFY_(STATUS)
    allocate(PE2(IM,JM,0:LM), STAT=STATUS); VERIFY_(STATUS)
    allocate(DP1(IM,JM,1:LM), STAT=STATUS); VERIFY_(STATUS)

! The import winds and mass fluxes need to be converted to FVAdv usable forms.
!  FVAdv provides a method to do this, hiding the nature of these forms. They
!  are however, at the FVAdv precision. They may be haloed and they need to be
!  given the proper length and area weights.
!-----------------------------------------------------------------------------

    call MAPL_TimerOn(MAPL,'-HORZ')

    dt_long = dt

    call FVAdv_PrepVars(state,  dt_long, ptop,    MAPL_GRAV,    &
                        uc,     vc,      mx_ur,   my_ur,        &
                        dp,     dpnew,   mx,      my,      mz,  &
                        cx,     cy,      mfx,     mfy,     pe2 )


! Set initial pressure thicknesses from the import state
!-------------------------------------------------------

    DP1 = DP

! We want the capability to split the step over which the mass fluxes
!   are valid into smaller steps for stability. Note only the friendly
!   tracers are passed here. In addition to doing a lagrangian update
!   to the tracers in TR, this also updates the layer thicknesses,
!   producing new values in DP1.
!----------------------------------------------------------------------

    ASSERT_(NS>0)

    DT_SHORT = DT / float(NS)

    do N=1,NS
       call FVAdv_Run(State, TR(1:NTRACERS), DP1, IORD, JORD,  &
                      CX, CY, MFX, MFY )
    enddo

    call MAPL_TimerOff(MAPL,'-HORZ')

! Begin vertical advection through Remap
!---------------------------------------

    call MAPL_TimerOn(MAPL,'-VERT')

    PE1(:,:,0) = PTOP
!!!    PE2(:,:,0) = PTOP

    do K=1,LM
       PE1(:,:,K) = PE1(:,:,K-1) + DP1(:,:,K)
!!!       PE2(:,:,K) = PE2(:,:,K-1) + DP(:,:,K)
    enddo

!
! Make sure PE1 and PE2 are equal at surface to avoid mapn_ppm_tracer bug
!
    PE1(:,:,LM) = PE2(:,:,LM)

! Replace import thicknesses if friendly
!---------------------------------------

    FRIENDLY = MAPL_VerifyFriendly(IMPORT, "DP", "ADVECTION", RC=STATUS)
    VERIFY_(STATUS)

!
! WS:  This is not working!   How do we specify variables as friendly?
!
!!!    if(FRIENDLY==ESMF_true) then
       DP = DP1
!!!    end if

! Remap the tracers. 
!-------------------

    call FVAdv_Remap(STATE, TR(1:NTRACERS), PE1, PE2, KORD)
    VERIFY_(STATUS)

    call MAPL_TimerOff(MAPL,'-VERT')

! Clean-up
!---------

    deallocate( TR, CX, CY, MFX, MFY, PE1, PE2, DP1, STAT=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOff(MAPL,'RUN'  )
    call MAPL_TimerOff(MAPL,'TOTAL')

! All done
!---------
    RETURN_(ESMF_SUCCESS)
  end subroutine Run


!------------------------------------------------------------------------------
!BOPI
! !IROUTINE:  Finalize - user supplied finalize routine

! !INTERFACE:

  subroutine Finalize(GC, IMPORT, EXPORT, CLOCK, RC)
!
! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC
    type(ESMF_State),    intent(inout) :: IMPORT, EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC
!
! !DESCRIPTION:
!    Finalize merely destroys the FVadv object that was created in Initialize
!    and releases the space for the persistent data .
!EOPI

    type (T_FVAdv_STATE), pointer :: state
    type (T_FVAdv_STATE_WRAP)     :: wrap

    character(len=ESMF_MAXSTR)    :: IAm
    integer                       :: STATUS
    character(len=ESMF_MAXSTR)    :: COMP_NAME

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'Finalize'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    call ESMF_UserCompGetInternalState( GC, 'FVAdvState', wrap, STATUS )
    VERIFY_(STATUS)

    state => wrap%state

    call FVAdv_Final    ( state, rc=status )
    VERIFY_(STATUS)

    deallocate(state,stat=STATUS)
    VERIFY_(STATUS)

    call MAPL_GenericFinalize(GC, IMPORT, EXPORT, CLOCK, RC)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine Finalize

!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!BOP
! !IROUTINE: TestCase
!
! !INTERFACE:

subroutine TestCase(gc, import, export, clock, rc)

! !USES:
    use FVAdvUtilsMod, only :  FVAdvUtilsTestCase
    implicit none

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc
    type(ESMF_State),    intent(inout) :: import
    type(ESMF_State),    intent(inout) :: export
    type (ESMF_Clock),   intent(in)    :: clock
    integer, intent(out), optional     :: rc

! !DESCRIPTION:
!   This routine provides a test scenario; usually this routine will
!   not be exercised since another component will provide the import 
!   state.  
!
!   If test data are desired, this routine can be activated
!   by setting the {\bf TEST:} variable in the configuration file 
!   to {\bf 1} (the default is {\bf 0}, or not active).
!
! !REVISION HISTORY:
!   2007-07-09    Sawyer	Creation from FVHS example

!EOP

    character(len=ESMF_MAXSTR)        :: IAm
    character(len=ESMF_MAXSTR)        :: COMP_NAME
    integer                           :: status

    type (T_FVAdv_STATE), pointer :: state
    type (T_FVAdv_STATE_WRAP)     :: wrap

    type (ESMF_State)     :: internal

    type(ESMF_Grid)        :: grid
    type(ESMF_Bundle)      :: tracers
    type(ESMF_Field)       :: field
    type(ESMF_Array)       :: array
    type(ESMF_FieldDataMap):: datamap
    type(ESMF_ArraySpec)   :: arrayspec
    type(ESMF_VM)          :: vm
    type(ESMF_Config)      :: cf


    type(MAPL_MetaComp), pointer    :: MAPL

    integer               :: im, jm, lm, lpet, nx, ny
    integer               :: ntracers
    integer               :: test
    real(FVAdvReal)       :: velocityScale
    real(FVAdvReal)       :: pi, ae, deltaT
    real                  :: dt
    integer               :: globalCounts(3)
    integer, allocatable  :: cellCounts(:,:), ims(:), jms(:)

! Pointers to imports

    real,            pointer     :: uc (:,:,:)
    real,            pointer     :: vc (:,:,:)
    real,            pointer     :: mx_ur (:,:,:)
    real,            pointer     :: my_ur (:,:,:)
    real,            pointer     :: mx (:,:,:)
    real,            pointer     :: my (:,:,:)
    real,            pointer     :: mz (:,:,:)
    real,            pointer     :: dp (:,:,:)
    real,            pointer     :: dpedt (:,:,:)
    real,            pointer     :: content(:,:,:)

    Iam = "Testcase"
    call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Get MAPL object
!----------------
    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"INIT"  )

    call MAPL_Get(MAPL, nx=nx, ny=ny, rc=status)
    VERIFY_(STATUS)

! BEWARE: im, jm are the local extents, not the global!!!
    call MAPL_Get(MAPL, lm=lm, rc=status)
    VERIFY_(STATUS)

! !RESOURCE_ITEM: seconds :: Total time interval for each call
    call MAPL_Get ( MAPL, INTERNAL_ESMF_STATE=INTERNAL, rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, Test, 'TEST:', default=0, rc=STATUS)
    VERIFY_(STATUS)
    if ( Test /= 0 ) then   ! Fill with test data
      call MAPL_GetResource(MAPL, DT  , 'RUN_DT:'    , default=1800., rc=STATUS)
      VERIFY_(STATUS)

! Convert these variables to default real precision
      velocityScale = 2 * MAPL_PI * MAPL_RADIUS / (12 * 86400)
      pi = MAPL_PI
      ae = MAPL_RADIUS
      deltaT = dt 

! WS 2007-07-19  The following statement causes the program to crash,
!                because this routine is performed *before* Initialize,
!                and thus the internal state has not yet been been set.
!                So we have to scrounge for the information in the GC Grid
!!!    call ESMF_UserCompGetInternalState( GC, 'FVAdvState', wrap, STATUS )
!!!    VERIFY_(STATUS)
!!!    state => wrap%state

! Get virtual machine and grid
!-----------------------------
      call ESMF_GridCompGet(GC, vm=vm, grid=grid, rc=STATUS)
      VERIFY_(STATUS)
      call ESMF_VMGet      (vm, localPet=lpet,    rc=STATUS)
      VERIFY_(STATUS)

      allocate(cellCounts(nx*ny,3), stat=STATUS); VERIFY_(STATUS)

      call ESMF_GridGet    ( grid,                                      &
                             horzRelloc            = ESMF_CELL_CENTER,  &      
                             Vertrelloc            = ESMF_CELL_CELL,    &      
                             cellCountPerDEPerDim  = cellCounts,        &
                             rc=STATUS )
      VERIFY_(STATUS)
      allocate(ims(nx), jms(ny))
      ims = cellCounts(1:NX   ,1)
      jms = cellCounts(1:  :NX,2)
      deallocate( cellCounts )
      im = sum(ims)
      jm = sum(jms)

! Get the import/export states
!-----------------------------
      call MAPL_GetPointer(IMPORT,   UC,    "UC",    rc=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT,   VC,    "VC",    rc=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT,   MX_UR, "MX_UR", rc=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT,   MY_UR, "MY_UR", rc=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT,   MX,    "MX",    rc=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT,   MY,    "MY",    rc=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT,   MZ,    "MZ",    rc=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT,DPEDT,"DPEDT", rc=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT,   DP,   "DP", rc=STATUS); VERIFY_(STATUS)


! Add test contents to the tracer bundle
      call ESMF_StateGetBundle(IMPORT, "TRACERS", Tracers, rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArraySpecSet(arrayspec, rank=3, type=ESMF_DATA_REAL, &
                             kind=ESMF_R4, rc=STATUS )
      VERIFY_(STATUS)

      call ESMF_FieldDataMapSetDefault(datamap,3,(/1,2,0/),          &
                                       counts=(/lm/), rc=STATUS)
      VERIFY_(STATUS)

      field = ESMF_FieldCreate(grid, arrayspec, datamap=datamap,&
                              horzRelloc=ESMF_CELL_CENTER,     &
                              haloWidth=0, name="constituent", &
                              rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_BundleAddField(Tracers, Field, rc=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldGetDataPointer(field, content, ESMF_DATA_REF,&
                                    rc=STATUS)
      VERIFY_(STATUS)
      call FVAdvUtilsTestCase( nx, ny, im, jm, lm, ims, jms,     &
                               deltaT, pi, ae, velocityScale,    &
                               uc, vc, mx_ur, my_ur, mx, my, mz, &
                               dp, dpedt, content)
      deallocate(ims, jms)

      call ESMF_BundleGet(Tracers, fieldCount=Ntracers   , rc=STATUS)
      VERIFY_(STATUS)

      print *, "Test Case now set up with Ntracers =", Ntracers
    endif

    call MAPL_TimerOff(MAPL,"INIT")
    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)

end subroutine TestCase

end module FVadvcore_GridCompMod
