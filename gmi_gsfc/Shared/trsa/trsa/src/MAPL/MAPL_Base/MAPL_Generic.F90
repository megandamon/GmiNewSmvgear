!  $Id: MAPL_Generic.F90,v 1.2 2011-08-09 22:13:00 mrdamon Exp $

#include "MAPL_ErrLog.h"
#define GET_POINTER ESMFL_StateGetPointerToData


!=============================================================================


module MAPL_GenericMod

!BOP
! \setlength{\parskip}{12pt};
!
! \part{\Large Module and Function Prologues}
!
! !MODULE: MAPL_GenericMod
!
! !DESCRIPTION:  MAPL\_Generic allows the user to easily build ESMF gridded
!  components.  It has its own SetServices, Initialize, Run, and Finalize
!  (IRF) methods, and thus is itself a valid gridded component, although somewhat
!  non-standard since it makes its IRF methods public. An instance of 
!  MAPL\_Generic does no useful work, but can be used as a null MAPL\_Generic component.
!
!  The standard way to use MAPL\_Generic is as an aid in building ESMF gridded
!  components. A MAPL/ESMF gridded component built in this way will always have
!  its own SetServices, which will call the subroutine MAPL\_GenericSetServices.
!  When MAPL\_GenericSetServices is called it sets the
!  component's IRF methods to the generic versions, MAPL\_GenericInitialize, MAPL\_GenericFinalize, and
!  MAPL\_GenericRun.  Any (or all) of
!  these may be used as default methods by a gridded component. (As we will see below, 
!  using all three default IRF methods in this way need not be equivalent to instanciating
!  a null component.) If for any of the
!  three IRF methods the default version is inadequate, it can simply be overrided
!  by having the component register its own method after the call to MAPL\_GenericSetServices.

!  The generic IRF methods perform a number of useful functions, including 
!  creating, allocating, and initializing the components Import, Export, 
!  and Internal states. It would be a shame to waste this capability when a component
!  needs to write its own version of an IRF method. A common situation is that the component wants support
!  in performing these functions, but needs to do some (usually small) additional specialized
!  work; for example, it may need to do some special initializations. In this case,
!  one would write a light version of the IRF method that does the specialized work
!  and {\it calls directly} the corresponding MAPL\_Generic method to do the boilerplate.
!  This is why MAPL\_Generic, unlike a standard ESMF gridded component, makes its 
!  IRF methods public and why we added the ``Generic'' modifier (i.e., MAPL\_GenericInitialize,
!  rather than MAPL\_Initialize), to emphasize that they are directly callable IRF methods. 
!  
!  MAPL\_Generic may also be viewed as a fairly standard Fortran 90 ``class,'' which
!  defines and makes public an opaque object that we refer to as a ``MAPL\_Generic State.'' 
!  This object can be created only in association with a standard ESMF Gridded Component (GC),
!  by making a MAPL\_GenericSetServices call.  This object can be obtained through an ESMF GC method
!  which is currently provided with MAPL. The MAPL\_Generic State is, therefore, just another thing that
!  lives in the ESMF GC, like the grid and the configuration. The MAPL\_Generic State
!  is private, but user components can access its contents through public
!  MAPL\_Generic methods (Get, Set, etc). The bulk of MAPL_Generic consists of methods that act
!  on this object.
!
!  MAPL\_GenericSetServices and MAPL\_Generic IRF methods cannot create their own ESMF grid.
!  The grid must be inherited from the parent or created by the component
!  either in its own SetServices or in its Initialize, if it is writing one. 
!  In any case, an important assumption of MAPL is that the grid must  already be {\it present
!  in the component and initialized} when MAPL\_GenericSetServices is invoked.
!  The same is true of the configuration.
!
!  In MAPL\_Generic, we distinguish between {\em simple (leaf)}
!  gridded compnents and {\em composite} gridded components, which contain other
!  ({\em child}) gridded components.  We also define three types of services,
!  which can be registered by the component's SetServices routine.

!  \begin{itemize}
!    \item {\bf Functional services}: 
!                   These are the standard EMSF callable IRF methods for
!         the component.
!
!    \item {\bf Data services:}
!                   These are descriptions of the component's import, export,
!         and internal states, which can be manipulated by MAPL\_Generic.
!
!    \item {\bf Child services:} 
!                   These are the services of the component's children and 
!         their connectivity.
!
!    \item {\bf Profiling Services:}
!                   These are profiling counters (clocks) that can be used
!         by the component and are automatically reported by generic finalize.
!  \end{itemize}

!   MAPL\_GenericSetServices provides generic versions of all these, as described below.


! !USES:

use ESMF_Mod
use ESMFL_Mod
use MAPL_BaseMod
use MAPL_IOMod
use MAPL_CFIOMod
use MAPL_ProfMod
use MAPL_CommsMod
use MAPL_ConstantsMod
use MAPL_SunMod
use MAPL_VarSpecMod
use MAPL_GenericCplCompMod
use MAPL_LocStreamMod

! !PUBLIC MEMBER FUNCTIONS:

implicit none
private

public MAPL_GenericSetServices
public MAPL_GenericInitialize
public MAPL_GenericRun
public MAPL_GenericFinalize

public MAPL_AddInternalSpec
public MAPL_AddImportSpec
public MAPL_AddExportSpec

public MAPL_GridCompSetEntryPoint
public MAPL_GetObjectFromGC
public MAPL_Get
public MAPL_Set
!public MAPL_GenericStateGet
!public MAPL_InternalStateGet
!public MAPL_GenericStateSet
!public MAPL_GenericRunCouplers

!public MAPL_StateAddInternalSpec
!public MAPL_StateAddImportSpec
!public MAPL_StateAddExportSpec
!public MAPL_StateAddSpec
!public MAPL_StateGetSpecAttrib
!public MAPL_StateSetSpecAttrib
!public MAPL_StateGetVarSpecs
!public MAPL_StatePrintSpec

! MAPL_Connect
public MAPL_AddChild
public MAPL_AddConnectivity
public MAPL_TerminateImport

! MAPL_Util
!public MAPL_GenericStateClockOn
!public MAPL_GenericStateClockOff
!public MAPL_GenericStateClockAdd
public MAPL_TimerOn
public MAPL_TimerOff
public MAPL_TimerAdd
public MAPL_ReadForcing
public MAPL_GetResource

! Internal public

public MAPL_StateAlarmAdd
public MAPL_StateAlarmGet
public MAPL_StateCreateFromSpec
public MAPL_FriendlyGet
public MAPL_GridCompGetFriendlies
public MAPL_GenericGridSet
public MAPL_SetVarSpecForCC
public MAPL_DateStampGet
public MAPL_ExchangeGridGet
public MAPL_ExchangeGridSet
public MAPL_ImportStateGet
public MAPL_ExportStateGet
public MAPL_GetChildLocstream
public MAPL_CopyFriendliness
public MAPL_VerifyFriendly
public MAPL_GenericMakeXchgNatural
public MAPL_GridCreate

! !INTERFACE:

interface MAPL_AddChild
   module procedure MAPL_AddChild
   module procedure MAPL_AddChildFromMeta
end interface

interface MAPL_StateAddImportSpec
   module procedure MAPL_StateAddImportSpec
   module procedure MAPL_StateAddImportSpecFrmChld
end interface

interface MAPL_StateAddExportSpec
   module procedure MAPL_StateAddExportSpec
   module procedure MAPL_StateAddExportSpecFrmChld
   module procedure MAPL_StateAddExportSpecFrmAll
end interface



interface MAPL_AddExportSpec
   module procedure MAPL_StateAddExportSpec
   module procedure MAPL_StateAddExportSpecFrmChld
   module procedure MAPL_StateAddExportSpecFrmAll
end interface

interface MAPL_AddInternalSpec
   module procedure MAPL_StateAddInternalSpec
end interface

interface MAPL_AddImportSpec
   module procedure MAPL_StateAddImportSpec
   module procedure MAPL_StateAddImportSpecFrmChld
end interface

interface MAPL_Get
   module procedure MAPL_GenericStateGet
end interface

interface MAPL_Set
   module procedure MAPL_GenericStateSet
end interface

interface MAPL_GetObjectFromGC
   module procedure MAPL_InternalStateGet
end interface

interface MAPL_TimerOn
   module procedure MAPL_GenericStateClockOn
end interface

interface MAPL_TimerOff
   module procedure MAPL_GenericStateClockOff
end interface

interface MAPL_TimerAdd
   module procedure MAPL_GenericStateClockAdd
end interface


interface MAPL_TerminateImport
   module procedure MAPL_DoNotConnect
   module procedure MAPL_DoNotConnectMany
   module procedure MAPL_DoNotConnectAnyImport
   module procedure MAPL_TerminateImportAll
end interface

interface MAPL_AddConnectivity
!   module procedure MAPL_AddConnectivityE2E
   module procedure MAPL_AddConnectivityRename
   module procedure MAPL_AddConnectivityRenameMany
   module procedure MAPL_AddConnectivityMany
!   module procedure MAPL_AddConnectivityOld
end interface

interface  MAPL_GridCompGetFriendlies
   module procedure MAPL_GridCompGetFriendlies0
   module procedure MAPL_GridCompGetFriendlies1
end interface

interface  MAPL_GetResource
   module procedure MAPL_GetResourceR4
   module procedure MAPL_GetResourceI4
   module procedure MAPL_GetResourceI41
   module procedure MAPL_GetResourceR8
   module procedure MAPL_GetResourceI8
   module procedure MAPL_GetResourceC
end interface

interface MAPL_CopyFriendliness
   module procedure MAPL_CopyFriendlinessInField
   module procedure MAPL_CopyFriendlinessInState
end interface

interface MAPL_VerifyFriendly
   module procedure MAPL_VerifyFriendlyInField
   module procedure MAPL_VerifyFriendlyInState
end interface

interface MAPL_ReadForcing
   module procedure MAPL_ReadForcing1
   module procedure MAPL_ReadForcing2
end interface

interface MAPL_GenericStateSet
   module procedure MAPL_GenericStateSet
   module procedure MAPL_GenericStateSetFromGC
end interface


! !PUBLIC TYPES:

public MAPL_MetaComp

!  \ev
!  This defines the MAPL\_Generic class. It is an opaque object that
!  can be queried using {\tt MAPL\_GenericStateGet}. An instance of
!  this type is placed in the default internal state location of the
!  ESMF gridded component by 
!  {\tt MAPL\_GenericSetServices}. This instance can be retreived by using  
!  {\tt MAPL\_InternalStateGet}.\bv

!EOP

! =======================================================================


integer, parameter :: LAST_ALARM = 99


type MAPL_GenericWrap
   type(MAPL_MetaComp       ), pointer :: INTERNAL_STATE
end type MAPL_GenericWrap

type MAPL_GenericGrid
   type(ESMF_Grid)                          :: ESMFGRID
   type(ESMF_Grid)                          :: HORZGRID
   type(ESMF_DELayout)                      :: LAYOUT
   real, pointer, dimension(:,:)            :: LATS
   real, pointer, dimension(:,:)            :: LONS
   integer                                  :: VERTDIM
   integer                                  :: IM, JM, LM ! Local counts
   integer                                  :: MYID, NX, NY, NX0, NY0
end type  MAPL_GenericGrid

type MAPL_GenericRecordType
   type(ESMF_Alarm), pointer                :: ALARM(:)
   character (len=ESMF_MAXSTR)              :: IMP_FNAME
   integer                                  :: IMP_LEN
   character (len=ESMF_MAXSTR)              :: INT_FNAME
   integer                                  :: INT_LEN
end type  MAPL_GenericRecordType

type  MAPL_MetaComp
   private
   type (ESMF_GridComp           ), pointer :: GCS(:)           => null()
   type (ESMF_State              ), pointer :: GIM(:)           => null()
   type (ESMF_State              ), pointer :: GEX(:)           => null()
   type (ESMF_CplComp            ), pointer :: CCS(:,:)         => null()
   type (ESMF_State              ), pointer :: CIM(:,:)         => null()
   type (ESMF_State              ), pointer :: CEX(:,:)         => null()
   logical,                         pointer :: CCcreated(:,:)   => null()
   type (MAPL_GenericGrid        )          :: GRID
   type (ESMF_Alarm              )          :: ALARM(0:LAST_ALARM)
   integer                                  :: ALARMLAST=0
   type (ESMF_Clock              )          :: CLOCK
   type (ESMF_Config             )          :: CF
   type (ESMF_State              )          :: INTERNAL
   type (MAPL_SunOrbit           )          :: ORBIT
   type (MAPL_VarSpec            ), pointer :: IMPORT_SPEC(:)   => null()
   type (MAPL_VarSpec            ), pointer :: EXPORT_SPEC(:)   => null()
   type (MAPL_VarSpec            ), pointer :: INTERNAL_SPEC(:) => null()
   type (MAPL_VarSpec            ), pointer :: FRCSPEC(:)  => null()
   type (MAPL_Prof               ), pointer :: TIMES(:)         => null()
   character(len=ESMF_MAXSTR)     , pointer :: GCNameList(:)    => null()
   type(ESMF_GridComp)                      :: RootGC
   type(ESMF_GridComp)            , pointer :: parentGC         => null()
   logical                                  :: ChildInit = .true.
   type (MAPL_Link)               , pointer :: LINK(:)          => null()
   type (MAPL_LocStream)                    :: ExchangeGrid
   type (MAPL_LocStream)                    :: LOCSTREAM
   character(len=ESMF_MAXSTR)               :: COMPNAME
   type (MAPL_GenericRecordType)  , pointer :: RECORD           => null()
   type (ESMF_State)                        :: FORCING
   integer                        , pointer :: phase_init (:)    => null()
   integer                        , pointer :: phase_run  (:)    => null()
   integer                        , pointer :: phase_final(:)    => null()
   integer                        , pointer :: phase_record(:)   => null()
   integer                        , pointer :: phase_coldstart(:)=> null()
end type MAPL_MetaComp

type MAPL_Connectivity
   type (MAPL_VarConn), pointer :: CONNECT(:)       => null()
   type (MAPL_VarConn), pointer :: DONOTCONN(:)     => null()
end type MAPL_Connectivity

type MAPL_ConnectivityWrap
   type(MAPL_Connectivity), pointer :: PTR
end type MAPL_ConnectivityWrap

type MAPL_LinkType
   type (ESMF_GridComp) :: GC
   integer              :: StateType
   integer              :: SpecId
end type MAPL_LinkType

type MAPL_LinkForm
   type (MAPL_LinkType) :: FROM
   type (MAPL_LinkType) :: TO
end type MAPL_LinkForm

type MAPL_Link
   type (MAPL_LinkForm), pointer :: PTR
end type MAPL_Link

type MAPL_MetaPtr
   type(MAPL_MetaComp), pointer  :: PTR
end type MAPL_MetaPtr

contains

#define LOWEST_(c) m=0; do while (m /= c) ;\
 m = c; c=label(c);\
enddo

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================



!BOPI

! !IROUTINE: MAPL_GenericSetServices

! !DESCRIPTION: {\tt MAPL\_GenericSetServices} performs the following tasks:

!\begin{itemize}
!\item
!  Allocate an instance of MAPL\_GenericState, wrap it, and set it as the
!  GC's internal state.
!\item
!  Exract the grid and configuration from the GC and save them in the 
!  generic state.
!\item
!  Set GC's IRF methods to the generic versions
!\item
!  If there are children
!\begin{itemize}
!\item
!   Allocate a gridded comoponent and an import and export state for each child
!\item
!   Create each child's GC using the natural grid and the inherited configuration.
!\item
!  Create each child's Import and Export states. These are named
!  {\tt GCNames(I)//"\_IMPORT"} and {\tt GCNames(I)//"\_EXPORT"}
!\item
!   Invoke each child's set services.
!\item
!   Add each item in each child's export state to GC's export state.
!\item
!   Add each item in each child's import state to GC's import, 
!   eliminating duplicates.  
!\end{itemize}
!\end{itemize}
! Since {\tt MAPL\_GenericSetServices} calls SetServices for the children,
! which may be generic themselves, the routine must be recursive.
!
! The optional arguments describe the component's children. There can be any 
!  number of children but they must be of one of the types specified by the
!  five SetServices entry points passed. If SSptr is not specified there can
!  only be five children, one for each {\tt SSn}, and the names must be in
!  {\tt SSn} order.
!  \newline

! !INTERFACE:

recursive subroutine MAPL_GenericSetServices ( GC,                  &
                                                                 RC )

! !ARGUMENTS:
!  

type(ESMF_GridComp),                  intent(INOUT) :: GC         ! Gridded component
integer,                    optional, intent(  OUT) :: RC         ! Return code
!EOPI


! ErrLog Variables
!-----------------

character(len=ESMF_MAXSTR)        :: IAm
character(len=ESMF_MAXSTR)        :: COMP_NAME
integer                           :: STATUS

! Local variables
! ---------------

character(len=ESMF_MAXSTR)        :: NAME
type (MAPL_MetaComp), pointer     :: INTERNAL_STATE
integer                           :: I
integer                           :: NC
integer                           :: TS
integer                           :: LBL, K, M
integer                           :: fLBL, tLBL
integer                           :: good_label, bad_label
integer, pointer                  :: LABEL(:)
integer, pointer                  :: SSnum(:)
type (MAPL_VarSpec),  pointer     :: SPECS(:)
type(ESMF_Field), pointer         :: FIELD
type(ESMF_VM)                     :: VM

type (MAPL_ConnectivityWrap)      :: connwrap
type (MAPL_VarConn),      pointer :: CONNECT(:)
type (MAPL_VarSpec),      pointer :: IM_SPECS(:)
type (MAPL_VarSpec),      pointer :: EX_SPECS(:)
type (MAPL_VarSpecPtr),   pointer :: ImSpecPtr(:)
type (MAPL_VarSpecPtr),   pointer :: ExSpecPtr(:)
type(ESMF_GridComp)               :: rootGC

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

   Iam = "MAPL_GenericSetServices"
   call ESMF_GridCompGet( GC, name=COMP_NAME, VM=VM, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // trim(Iam)

! Create the generic state, intializing its configuration and grid.
!  The parents grid is probably not defined at this point
!----------------------------------------------------------

   call MAPL_InternalStateRetrieve( GC, INTERNAL_STATE, RC=STATUS)
   VERIFY_(STATUS)

   INTERNAL_STATE%COMPNAME = COMP_NAME


! Set the Component's Total timer
! -------------------------------

    call MAPL_GenericStateClockAdd(GC,    name="TOTAL"  ,RC=STATUS)
    VERIFY_(STATUS)

! Fill in the rootGC
!-------------------
    call MAPL_GetRootGC(GC, rootGC, rc=status)
    VERIFY_(STATUS)
    INTERNAL_STATE%rootGC = rootGC

! If this is a composite component, Setup the children
! ----------------------------------------------------

   CHILDREN: if(associated(INTERNAL_STATE%GCNameList)) then

      NC = size(INTERNAL_STATE%GCNameList)

  
    do I=1,NC
       call MAPL_GenericStateClockAdd(GC, name=trim(INTERNAL_STATE%GCNameList(I))  ,RC=STATUS)
       VERIFY_(STATUS)
    end do


! The child should've been already created by MAPL_AddChild
! and set his services should've been called.
! -------------------------------------

! Create internal couplers and composite
! component's Im/Ex specs.
!---------------------------------------

      call MAPL_WireComponent(GC, RC=STATUS)
      VERIFY_(STATUS)

! Relax connectivity for non-existing imports
      if (NC > 0) then
         call ESMF_UserCompGetInternalState(gc, 'M0APL_Connectivity', &
              connwrap, status)
         if (STATUS == ESMF_FAILURE) then
            NULLIFY(CONNECT)
         else
            CONNECT => connwrap%ptr%CONNECT
         end if

         allocate (ImSpecPtr(NC), ExSpecPtr(NC), stat=status)
         VERIFY_(STATUS)

         DO I = 1, NC
            call MAPL_GridCompGetVarSpecs(INTERNAL_STATE%GCS(I), &
                 IMPORT=IM_SPECS, EXPORT=EX_SPECS, RC=STATUS)
            VERIFY_(STATUS)
            ImSpecPtr(I)%Spec => IM_SPECS
            ExSpecPtr(I)%Spec => EX_SPECS
         END DO

         call MAPL_ConnCheckReq(CONNECT, ImSpecPtr, ExSpecPtr, rc=status)
         VERIFY_(STATUS)

         deallocate (ImSpecPtr, ExSpecPtr)

      end if

! If I am root call Label from here; everybody else
!  will be called recursively from Label
!--------------------------------------------------
   ROOT: if (.not. associated(INTERNAL_STATE%parentGC)) then

      call MAPL_GenericConnCheck(GC, RC=status)
      VERIFY_(STATUS)

! Collect all IMPORT and EXPORT specs in the entire tree in one list
!-------------------------------------------------------------------

      nullify(SPECS)
      call MAPL_GenericSpecEnum(GC, SPECS, RC=STATUS)
      VERIFY_(STATUS)

! Label each spec by its place on the list--sort of.
!--------------------------------------------------

      TS = size(SPECS)
      allocate(LABEL(TS), STAT=STATUS)
      VERIFY_(STATUS)

      do I = 1, TS
         LABEL(I)=I
      end do

! For each spec...
!-----------------

      do I = 1, TS

!  Get the LABEL attribute on the spec
!-------------------------------------

         call MAPL_VarSpecGet(SPECS(I), LABEL = LBL, RC = STATUS)
         VERIFY_(STATUS)
         if (LBL <= 0) then
            RETURN_(ESMF_FAILURE)
         endif

! Do something to sort labels???
!-------------------------------

         LOWEST_(LBL)

         if (LBL < I) then
            good_label = LBL
            bad_label  = I
         else
            good_label = I
            bad_label  = LBL
         end if
         label(bad_label) = good_label
      end do

      if (associated(INTERNAL_STATE%LINK)) then
        do I = 1, size(INTERNAL_STATE%LINK)
           fLBL = MAPL_LabelGet(INTERNAL_STATE%LINK(I)%ptr%FROM, RC=STATUS)
           VERIFY_(STATUS)
           tLBL = MAPL_LabelGet(INTERNAL_STATE%LINK(I)%ptr%TO,   RC=STATUS)
           VERIFY_(STATUS)
           LOWEST_(fLBL)
           LOWEST_(tLBL)

           if (fLBL < tLBL) then
              good_label = fLBL
              bad_label  = tLBL
           else
              good_label = tLBL
              bad_label  = fLBL
           end if
           label(bad_label) = good_label
        end do
      end if

      K=0
      do I = 1, TS
         LBL = LABEL(I)
         LOWEST_(LBL)
         call MAPL_VarSpecGet(SPECS(LBL), FIELDPTR = FIELD, RC=STATUS  )
         VERIFY_(STATUS)

         if (LBL == I) then
            K = K+1
         else
            call MAPL_VarSpecSet(SPECS(I), FIELDPTR = FIELD, RC=STATUS  )
            VERIFY_(STATUS)
         end if

         call MAPL_VarSpecSet(SPECS(I), LABEL=LBL, RC=STATUS)
         VERIFY_(STATUS)
      end do

      DEALLOCATE(LABEL, STAT=status)
      VERIFY_(STATUS)

   end if ROOT

end if CHILDREN  !  Setup children

! Register default services for this component if needed
! --------------------------------------------

   if (.not. associated(INTERNAL_STATE%phase_init)) then
      call MAPL_GridCompSetEntrypoint(GC, ESMF_SETINIT, MAPL_GenericInitialize,  RC=STATUS)
      VERIFY_(STATUS)
   endif

   if (.not. associated(INTERNAL_STATE%phase_run)) then
      call MAPL_GridCompSetEntrypoint(GC, ESMF_SETRUN, MAPL_GenericRun,  RC=STATUS)
      VERIFY_(STATUS)
   endif


   if (.not. associated(INTERNAL_STATE%phase_final)) then
      call MAPL_GridCompSetEntrypoint(GC, ESMF_SETFINAL, MAPL_GenericFinalize,  RC=STATUS)
      VERIFY_(STATUS)
   endif

!ALT check record!!!
   if (.not. associated(INTERNAL_STATE%phase_record)) then
      call MAPL_GridCompSetEntryPoint ( GC, "ESMF_WriteRestart", MAPL_GenericRecord, RC=STATUS)
      VERIFY_(STATUS)
   end if
   ASSERT_(size(INTERNAL_STATE%phase_record)==1)  !ALT: currently we support only 1 record (99)
   

   if (associated(INTERNAL_STATE%phase_coldstart)) then
      ASSERT_(size(INTERNAL_STATE%phase_coldstart)==1) !ALT: currently we support only 1 coldstart (99)
   endif

!ALT ATTENTION HERE!!!!

! Timers for generic initialize and finalize
!-------------------------------------------

   call MAPL_GenericStateClockAdd(GC, name="GenInitTot"     ,RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GenericStateClockAdd(GC, name="--GenInitMine"  ,RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GenericStateClockAdd(GC, name="GenRunTot"      ,RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GenericStateClockAdd(GC, name="--GenRunMine"   ,RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GenericStateClockAdd(GC, name="GenFinalTot"    ,RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GenericStateClockAdd(GC, name="--GenFinalMine" ,RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GenericStateClockAdd(GC, name="GenRecordTot"   ,RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GenericStateClockAdd(GC, name="--GenRecordMine",RC=STATUS)
   VERIFY_(STATUS)


! All done
!---------

   RETURN_(ESMF_SUCCESS)

end subroutine MAPL_GenericSetServices

#undef LOWEST_
!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

!BOPI

! !IROUTINE: MAPL_GenericInitialize -- Initializes the component and its children

! !INTERFACE:

recursive subroutine MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component 
type(ESMF_State),    intent(INOUT) :: IMPORT ! Import state
type(ESMF_State),    intent(INOUT) :: EXPORT ! Export state
type(ESMF_Clock),    intent(INOUT) :: CLOCK  ! The clock
integer, optional,   intent(  OUT) :: RC     ! Error code:

!EOPI

! ErrLog Variables

  character(len=ESMF_MAXSTR)    :: IAm
  character(len=ESMF_MAXSTR)    :: COMP_NAME
  integer                       :: STATUS

! Local derived type aliases

  type (MAPL_MetaComp),pointer     :: STATE 
  type (MAPL_GenericGrid ),pointer :: MYGRID
  
! Local variables

  type (ESMF_Time  )            :: ORIGIN    
  type (ESMF_Time  )            :: ringTime
  type (ESMF_Time  )            :: FILETIME
  type (ESMF_TimeInterval)      :: TIMEINT
  type (ESMF_Calendar)          :: cal
  character(len=ESMF_MAXSTR)    :: FILENAME
  character(len=ESMF_MAXSTR)    :: FILETYPE
  real                          :: DT
  real                          :: OFFSET
  real                          :: DEFDT
  integer                       :: COUNTS(3)
  integer                       :: L
  integer                       :: I, J, N
  integer                       :: I1, IN, J1, JN
  integer                       :: NC
  integer                       :: NSUBTILES
  integer                       :: DIMCOUNT
  integer                       :: DECOUNT(2)
  type (ESMF_Grid)              :: TILEGRID
  type (ESMF_Alarm)             :: recordAlarm
  integer, dimension(:), allocatable :: ref_date, ref_time, ref_freq
  integer                       :: NRA, sec
  character(len=ESMF_MAXSTR)    :: AlarmName
  type(ESMF_Time)               :: CurrTime    ! Current time of the ESMF clock
  type(ESMF_Time)               :: RefTime
  type(ESMF_TimeInterval)       :: Frequency
  character(len=ESMF_MAXSTR)    :: CHILD_NAME
  type(ESMF_IOSpec)             :: iospec
  type(ESMF_Grid)               :: CHLGRID

  integer                          :: nymd,nhms  ! Current Time date and hour/minute
  integer                          :: PHASE
  integer                          :: NUMPHASES
  integer                          :: MAXPHASES
  integer, allocatable             :: CHLDPHASES(:)
  type (MAPL_MetaPtr), allocatable :: CHLDMAPL(:)

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

  Iam = "MAPL_GenericInitialize"
  call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
  VERIFY_(STATUS)
  Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the internal state.
! -------------------------------------------

  call MAPL_InternalStateGet ( GC, STATE, RC=STATUS)
  VERIFY_(STATUS)

! Start my timer
!---------------

  call MAPL_GenericStateClockOn(STATE,"TOTAL")
  call MAPL_GenericStateClockOn(STATE,"GenInitTot")
  call MAPL_GenericStateClockOn(STATE,"--GenInitMine")

! Put the inherited grid in the generic state
!--------------------------------------------

  MYGRID    =>  STATE%GRID
  call ESMF_GridCompGet( GC, GRID = MYGRID%ESMFGRID, RC=STATUS )
  VERIFY_(STATUS)

! At this point, this component must have a valid grid!
!------------------------------------------------------
  call ESMF_GridValidate(MYGRID%ESMFGRID, RC=STATUS)
  VERIFY_(STATUS)

! Check children's grid. If they don't have a valid grid yet, put this one in their GC
! ------------------------------------------------------------------------------------
  if(associated(STATE%GCS)) then
      do I=1, size(STATE%GCS)
         call ESMF_GridCompGet(STATE%GCS(I), grid=ChlGrid, rc=status)
         VERIFY_(STATUS)
         call ESMF_GridValidate(ChlGrid, RC=STATUS)
         if (STATUS /= ESMF_SUCCESS) then
! This child does not have a valid grid
            call ESMF_GridCompset( STATE%GCS(I), GRID = MYGRID%ESMFGRID, RC=STATUS )
            VERIFY_(STATUS)
         end if
      end do
   end if

! Put the horizontal cross-section of the grid in the component's grid
! For now we support only horizontal ESMF grids.
! --------------------------------------------------------------------

  MYGRID%HORZGRID = MYGRID%ESMFGRID
  !call ESMFL_GridGetHorizontal(MYGRID%ESMFGRID, MYGRID%HORZGRID, RC=STATUS)
  !VERIFY_(STATUS)

! We keep these in the component's grid  for convenience
!-------------------------------------------------------

  call ESMF_GridGet(MYGRID%ESMFGRID, deLayout=MYGRID%LAYOUT, RC=STATUS)
  VERIFY_(STATUS)

  call ESMF_DELayoutGet(MYGRID%LAYOUT, dimCount = DIMCOUNT ,RC=STATUS)
  VERIFY_(STATUS)

! For now we support only 2-D layouts
! ---------------------------------

  ASSERT_(DIMCOUNT == 2)

! Processors in each direction
!-----------------------------

  call ESMF_DELayoutGet(MYGRID%LAYOUT, deCountPerDim = DECOUNT, &
                        localDE=MYGRID%MYID, RC=STATUS)
  VERIFY_(STATUS)

  MYGRID%NX = DECOUNT(1)
  MYGRID%NY = DECOUNT(2)

! My processor coordinates
!-------------------------

  call ESMF_DELayoutGetDELocalInfo(delayout=MYGRID%LAYOUT, de=MYGRID%MYID, coord=DECOUNT, rc=status)
  VERIFY_(STATUS)

  MYGRID%NX0 = DECOUNT(1)
  MYGRID%NY0 = DECOUNT(2)

! Vertical coordinate must exist and be THE THIRD DIMENSION
! ---------------------------------------------------------

  MYGRID%VERTDIM = 3

  call ESMF_GridGetDELocalInfo(MYGRID%ESMFGRID, &
       horzRelLoc=ESMF_CELL_CENTER, &
       vertRelLoc=ESMF_CELL_CELL, &
       localCellCountPerDim=COUNTS,RC=STATUS)
  VERIFY_(STATUS)

! Local sizes of three dimensions
!--------------------------------

  MYGRID%IM = COUNTS(1)
  MYGRID%JM = COUNTS(2)
  MYGRID%LM = COUNTS(3)

! Create and initialize factors saved as ESMF arrays in MYGRID
!-------------------------------------------------------------

  call ESMFL_GridCoordGet(   MYGRID%HORZGRID, MYGRID%LATS       , &
                             Name     = "Latitude"              , &
                             Location = ESMF_CELL_CENTER        , &
                             Units    = ESMFL_UnitsRadians      , &
                             RC       = STATUS                    )
  VERIFY_(STATUS)

  call ESMFL_GridCoordGet(   MYGRID%HORZGRID, MYGRID%LONS       , &
                             Name     = "Longitude"             , &
                             Location = ESMF_CELL_CENTER        , &
                             Units    = ESMFL_UnitsRadians      , &
                             RC       = STATUS                    )
  VERIFY_(STATUS)

! Put the clock passed down in the generic state
!-----------------------------------------------

   STATE%CLOCK = CLOCK

! We get our calling interval from the configuration,
! set the alarm, and attach it to the callers clock.
! ---------------------------------------------------
 
   call ESMF_ConfigGetAttribute( state%CF, DEFDT, Label="RUN_DT:", RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetResource( STATE   , DT, Label="DT:", default=DEFDT, RC=STATUS)
   VERIFY_(STATUS)
   call ESMF_ClockGet(clock, calendar = cal, currTime=currTime, rc=status)
   VERIFY_(STATUS)

   call MAPL_GetResource( STATE   , OFFSET, Label="ALARM_OFFSET:", default=0., RC=STATUS)
   VERIFY_(STATUS)
   call ESMF_TimeIntervalSet(TIMEINT,  S=nint(OFFSET) , calendar=cal, RC=STATUS)
   VERIFY_(STATUS)

   ringTime = currTime+TIMEINT

   call ESMF_TimeIntervalSet(TIMEINT,  S=nint(DT) , calendar=cal, RC=STATUS)
   VERIFY_(STATUS)

   STATE%ALARM(0) = ESMF_AlarmCreate(trim(COMP_NAME) // "_Alarm" , &
        CLOCK = CLOCK, &
        RingInterval = TIMEINT  ,  &
        RingTime     = ringTime,  & 
!        Enabled      = .true.   ,  &
!        sticky       = .false.  ,  &
        RC           = STATUS      )
   VERIFY_(STATUS)
!   call ESMF_AlarmRingerOn(STATE%ALARM, rc=status); VERIFY_(STATUS)

! Create tiling for all gridded components with associated LocationStream
! -----------------------------------------------------------------------

   if (MAPL_LocStreamIsAssociated(STATE%LOCSTREAM, RC=STATUS)) then
      NSUBTILES = MAPL_GetNumSubtiles(STATE, RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_LocStreamAttachGrid(STATE%LocStream, MYGRID%ESMFGRID, NSUBTILES, RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_LocStreamGet(STATE%LocStream, TILEGRID=TILEGRID, RC=STATUS)
      VERIFY_(STATUS)

   endif

! Create import and initialize state variables
! --------------------------------------------

   if (associated(STATE%IMPORT_SPEC)) then
      if (MAPL_LocStreamIsAssociated(STATE%LOCSTREAM, RC=STATUS)) then
         call MAPL_StateCreateFromSpec(IMPORT,STATE%IMPORT_SPEC,     &
                                       MYGRID%ESMFGRID,              &
                                       TILEGRID=TILEGRID,            &
                                       RC=STATUS       )
      else
         call MAPL_StateCreateFromSpec(IMPORT,STATE%IMPORT_SPEC,     &
                                       MYGRID%ESMFGRID,              &
                                       RC=STATUS       )
      endif
      VERIFY_(STATUS)

      call MAPL_GetResource( STATE   , FILENAME,         &
                                 LABEL="IMPORT_RESTART_FILE:", &
                                 RC=STATUS)
      if(STATUS==ESMF_SUCCESS) then
         call MAPL_GetResource( STATE   , FILETYPE,         &
                                 default='binary', &
                                 LABEL="IMPORT_RESTART_TYPE:", &
                                 RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_ESMFStateReadFromFile(IMPORT, CLOCK, FILENAME, &
                                         FILETYPE, RC=STATUS)
         if (STATUS /= ESMF_SUCCESS) then
            if (MAPL_AM_I_Root()) then
               call ESMF_StatePrint(Import)
            end if
         end if
         VERIFY_(STATUS)
      endif
   end if

! Create internal and initialize state variables
! -----------------------------------------------

   STATE%INTERNAL = ESMF_StateCreate(statename = trim(COMP_NAME) // "_INTERNAL", &
                              RC=STATUS)
   VERIFY_(STATUS)

   if (associated(STATE%INTERNAL_SPEC)) then
      if (MAPL_LocStreamIsAssociated(STATE%LOCSTREAM, RC=STATUS)) then
         call MAPL_StateCreateFromSpec(STATE%INTERNAL,STATE%INTERNAL_SPEC, &
                                        MYGRID%ESMFGRID,             &
                                        TILEGRID=TILEGRID,           &
                                        RC=STATUS       )
      else
         call MAPL_StateCreateFromSpec(STATE%INTERNAL,STATE%INTERNAL_SPEC, &
                                        MYGRID%ESMFGRID,             &
                                        RC=STATUS       )
      end if
      VERIFY_(STATUS)

      call MAPL_GetResource( STATE   , FILENAME,         &
                                 LABEL="INTERNAL_RESTART_FILE:", &
                                 RC=STATUS)
      if(STATUS==ESMF_SUCCESS) then
         call MAPL_GetResource( STATE   , FILETYPE,         &
                                 default='binary', &
                                 LABEL="INTERNAL_RESTART_TYPE:", &
                                 RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_ESMFStateReadFromFile(STATE%INTERNAL, CLOCK, FILENAME, &
                                         FILETYPE, RC=STATUS)
         if (STATUS /= ESMF_SUCCESS) then
            if (MAPL_AM_I_Root()) then
               call ESMF_StatePrint(STATE%INTERNAL)
            end if
            RETURN_(ESMF_FAILURE)
         end if
      else
! try to coldstart the internal state 
! -------------------------------
         if (associated(STATE%phase_coldstart)) then
            call ESMF_GridCompInitialize(GC, import, export, CLOCK, &
                 PHASE=MAPL_ColdstartPhase, BlockingFlag=ESMF_NONBLOCKING, RC=STATUS)
            VERIFY_(STATUS)
         endif

      endif
   end if

! Create export state variables
!------------------------------

   if (associated(STATE%EXPORT_SPEC)) then
      if (MAPL_LocStreamIsAssociated(STATE%LOCSTREAM, RC=STATUS)) then
         call MAPL_StateCreateFromSpec(EXPORT,STATE%EXPORT_SPEC,     &
                                       MYGRID%ESMFGRID,              &
                                       TILEGRID=TILEGRID,            &
                                       DEFER=.true., RC=STATUS       )
      else
         call MAPL_StateCreateFromSpec(EXPORT,STATE%EXPORT_SPEC,     &
                                       MYGRID%ESMFGRID,              &
                                       DEFER=.true., RC=STATUS       )
      end if
      VERIFY_(STATUS)
   end if

! Create forcing state
   STATE%FORCING = ESMF_StateCreate(statename = trim(COMP_NAME) // "_FORCING", &
                              RC=STATUS)
   VERIFY_(STATUS)
! Check if user wants RECORD feature

   call ESMF_ConfigFindLabel( STATE%CF, LABEL="RECORD_FREQUENCY:", RC=STATUS)
   if (STATUS==ESMF_SUCCESS) then
      nra = ESMF_ConfigGetLen( STATE%CF, RC = STATUS)
      ASSERT_( NRA > 0 .and. NRA < 10)

      allocate (ref_date(NRA), ref_time(NRA), ref_freq(NRA), stat=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigFindLabel( STATE%CF, LABEL="RECORD_FREQUENCY:", RC=STATUS)
      call ESMF_ConfigGetattribute( STATE%CF, valueList=ref_freq, count=NRA, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigFindLabel( STATE%CF, LABEL="RECORD_REF_DATE:", RC=STATUS)
      VERIFY_(STATUS)
!      ASSERT_(NRA == ESMF_ConfigGetLen(STATE%CF))
      call ESMF_ConfigGetattribute( STATE%CF, valueList=ref_date, count=NRA, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigFindLabel( STATE%CF, LABEL="RECORD_REF_TIME:", RC=STATUS)
      VERIFY_(STATUS)
!      ASSERT_(NRA == ESMF_ConfigGetLen(STATE%CF))
      call ESMF_ConfigGetattribute( STATE%CF, valueList=ref_time, count=NRA, RC=STATUS)
      VERIFY_(STATUS)

      allocate(STATE%RECORD, stat=STATUS)
      VERIFY_(STATUS)

      call MAPL_GetResource( STATE, FILENAME,         &
                             LABEL="IMPORT_CHECKPOINT_FILE:", &
                             RC=STATUS)
      if(STATUS==ESMF_SUCCESS) then
         STATE%RECORD%IMP_FNAME = FILENAME
         STATE%RECORD%IMP_LEN = LEN_TRIM(FILENAME)
      else
         STATE%RECORD%IMP_LEN = 0
      end if
         
      call MAPL_GetResource( STATE, FILENAME,         &
                             LABEL="INTERNAL_CHECKPOINT_FILE:", &
                             RC=STATUS)
      if(STATUS==ESMF_SUCCESS) then
         STATE%RECORD%INT_FNAME = FILENAME
         STATE%RECORD%INT_LEN = LEN_TRIM(FILENAME)
      else
         STATE%RECORD%INT_LEN = 0
      end if

      allocate (STATE%RECORD%ALARM(NRA), STAT=STATUS)
      VERIFY_(STATUS)

      DO  I = 1, NRA
         AlarmName = "RecordAlarm" // CHAR(I+ICHAR('0')) 
         call ESMF_ClockGetAlarm(clock, trim(AlarmName), recordAlarm, rc=status)
         if (STATUS/=ESMF_SUCCESS) then
            ! create alarm
            call ESMF_TimeSet( RefTime, YY = ref_date(I)/10000, &
                               MM = mod(ref_date(I),10000)/100, &
                               DD = mod(ref_date(I),100), &
                               H = ref_time(I)/10000, &
                               M = mod(ref_time(I),10000)/100, &
                               S = mod(ref_time(I),100), calendar=cal, rc=status )
            VERIFY_(STATUS)

            nhms = ref_freq(I)
            sec = nhms/10000*3600 + mod(nhms,10000)/100*60 + mod(nhms,100)
            call ESMF_TimeIntervalSet( frequency, S=sec, rc=status )
            VERIFY_(STATUS)
            RingTime = RefTime
            if (RingTime < currTime .and. sec /= 0) then
               RingTime = RingTime + (INT((currTime - RingTime)/frequency)+1)*frequency
            endif

            RecordAlarm = ESMF_AlarmCreate( trim(AlarmName), clock=clock, RingInterval=Frequency, &
                                             RingTime=RingTime, sticky=.false.,rc=status )
            VERIFY_(STATUS)

         end if
         STATE%RECORD%alarm(I) = recordAlarm
      END DO
      deallocate (ref_freq, ref_time, ref_date)

   endif


         
   
! Put the Export state of each child into my export
! -------------------------------------------------

!ALT: export might have to be declared ESMF_STATELIST
   if(associated(STATE%GCS)) then
      do I=1, size(STATE%GCS)
         call ESMF_StateAddState(EXPORT, STATE%GEX(I), RC=STATUS)
         VERIFY_(STATUS)
      end do
   end if

  call MAPL_GenericStateClockOff(STATE,"--GenInitMine")

! Initialize the children
! -----------------------

   if(associated(STATE%GCS)) then
      NC = size(STATE%GCS)
      if (STATE%ChildInit) then
         allocate(CHLDMAPL(NC), stat=status)
         MAXPHASES = 0
         do I=1,NC
            call MAPL_GetObjectFromGC(STATE%GCS(I), CHLDMAPL(I)%PTR, RC=STATUS)
            VERIFY_(STATUS)
            MAXPHASES = MAX(MAXPHASES, SIZE(CHLDMAPL(I)%PTR%PHASE_INIT))
         end do

         do PHASE = 1, MAXPHASES
            do I=1,NC
               NUMPHASES = SIZE(CHLDMAPL(I)%PTR%PHASE_INIT)
               if (PHASE .le. NUMPHASES) then
                  call ESMF_GridCompGet( STATE%GCS(I), NAME=CHILD_NAME, RC=STATUS )
                  VERIFY_(STATUS)
      
                  call MAPL_GenericStateClockOn (STATE,trim(CHILD_NAME))
                  call ESMF_GridCompInitialize (STATE%GCS(I), STATE%GIM(I), STATE%GEX(I), &
                       CLOCK, PHASE=CHLDMAPL(I)%PTR%PHASE_INIT(PHASE), RC=STATUS )
                  VERIFY_(STATUS)
                  call MAPL_GenericStateClockOff(STATE,trim(CHILD_NAME))
               end if
            end do
         end do
         deallocate(CHLDMAPL)
      end if

! Initialize all needed couplers
! ---------------------------------------------------

      do I=1,NC
         do J=1,NC
            if(STATE%CCcreated(J,I)) then
!               call WRITE_PARALLEL( "DEBUG: initilaizing CPL in " // &
!                    trim(comp_name) // " for " // &
!                    trim(STATE%GCNameList(J)) // " and " // &
!                    trim(STATE%GCNameList(I)))
               call ESMF_CplCompInitialize (STATE%CCS(J,I), STATE%GEX(J), STATE%GIM(I), CLOCK, RC=STATUS )
               VERIFY_(STATUS)
            endif
         enddo
! ---------------------------------------------------
      enddo
   endif

   if (.not. associated(STATE%parentGC)) then
      call MAPL_AdjustIsNeeded(GC, EXPORT, RC=STATUS)
      VERIFY_(STATUS)
   end if

  call MAPL_GenericStateClockOff(STATE,"GenInitTot")
  call MAPL_GenericStateClockOff(STATE,"TOTAL")

  RETURN_(ESMF_SUCCESS)

end subroutine MAPL_GenericInitialize

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

!BOPI

! !IROUTINE: MAPL_GenericRun

! !INTERFACE:

recursive subroutine MAPL_GenericRun ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component 
type(ESMF_State),    intent(INOUT) :: IMPORT ! Import state
type(ESMF_State),    intent(INOUT) :: EXPORT ! Export state
type(ESMF_Clock),    intent(INOUT) :: CLOCK  ! The clock
integer, optional,   intent(  OUT) :: RC     ! Error code:

!EOPI

! ErrLog Variables


character(len=ESMF_MAXSTR)    :: IAm
character(len=ESMF_MAXSTR)    :: COMP_NAME
integer                       :: STATUS

! Local derived type aliases

type (MAPL_MetaComp),pointer :: STATE 

character(len=ESMF_MAXSTR)       :: CHILD_NAME
integer                          :: I, J
integer                          :: NC
integer                          :: PHASE
integer                          :: NUMPHASES
integer                          :: MAXPHASES
integer, allocatable             :: CHLDPHASES(:)
type (MAPL_MetaPtr), allocatable :: CHLDMAPL(:)
!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

Iam = "MAPL_GenericRun"
call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
VERIFY_(STATUS)
Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the internal state. It comes in a wrapper.
! ------------------------------------------------------------------

call MAPL_InternalStateGet ( GC, STATE, RC=STATUS)
VERIFY_(STATUS)

call MAPL_GenericStateClockOn (STATE,"TOTAL")
call MAPL_GenericStateClockOn (STATE,"GenRunTot")
call MAPL_GenericStateClockOn (STATE,"--GenRunMine")
call MAPL_GenericStateClockOff(STATE,"--GenRunMine")

! Run the children
! ----------------

if(associated(STATE%GCS)) then
   NC = size(STATE%GCS)
   allocate(CHLDMAPL(NC), stat=status)
   MAXPHASES = 0
   do I=1,NC
      call MAPL_GetObjectFromGC(STATE%GCS(I), CHLDMAPL(I)%PTR, RC=STATUS)
      VERIFY_(STATUS)
      MAXPHASES = MAX(MAXPHASES, SIZE(CHLDMAPL(I)%PTR%PHASE_RUN))
   end do

   do PHASE = 1, MAXPHASES
      do I=1,NC
         NUMPHASES = SIZE(CHLDMAPL(I)%PTR%PHASE_RUN)
         if (PHASE .le. NUMPHASES) then
            call ESMF_GridCompGet( STATE%GCS(I), NAME=CHILD_NAME, RC=STATUS )
            VERIFY_(STATUS)
      
            call MAPL_GenericStateClockOn (STATE,trim(CHILD_NAME))
            call ESMF_GridCompRun (STATE%GCS(I), STATE%GIM(I), STATE%GEX(I), &
                CLOCK, PHASE=CHLDMAPL(I)%PTR%PHASE_RUN(PHASE), RC=STATUS )
            VERIFY_(STATUS)
            call MAPL_GenericStateClockOff(STATE,trim(CHILD_NAME))
         end if

         if (PHASE == NUMPHASES) then
            do J=1,NC
               if(STATE%CCcreated(I,J)) then
                  call ESMF_CplCompRun (STATE%CCS(I,J), STATE%GEX(I), STATE%GIM(J), &
                       CLOCK, RC=STATUS )
                  VERIFY_(STATUS)
               endif
            enddo
         end if
      enddo
   enddo
   deallocate(CHLDMAPL)
endif

call MAPL_GenericStateClockOff(STATE,"GenRunTot")
call MAPL_GenericStateClockOff(STATE,"TOTAL")

RETURN_(ESMF_SUCCESS)

end subroutine MAPL_GenericRun


!BOPI

! !IROUTINE: MAPL_GenericFinalize -- Finalizes the component and its children

! !INTERFACE:

  recursive subroutine MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! composite gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! import state
    type(ESMF_State),    intent(inout) :: EXPORT ! export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! the clock
    integer, optional,   intent(  out) :: RC     ! Error code:
                                                 ! = 0 all is well
                                                 ! otherwise, error
!EOPI

! LOCAL VARIABLES

  character(len=ESMF_MAXSTR)                  :: IAm
  character(len=ESMF_MAXSTR)                  :: COMP_NAME
  integer                                     :: STATUS

  character(len=ESMF_MAXSTR)                  :: FILENAME
  character(len=ESMF_MAXSTR)                  :: FILETYPE
  character(len=ESMF_MAXSTR)                  :: CHILD_NAME
  character(len=ESMF_MAXSTR)                  :: RECFIN
  type (MAPL_MetaComp), pointer               :: STATE
  integer                                     :: I
  logical                                     :: final_checkpoint
  integer                                     :: NC
  integer                                     :: PHASE
  integer                                     :: NUMPHASES
  integer                                     :: MAXPHASES
  integer, allocatable                        :: CHLDPHASES(:)
  type (MAPL_MetaPtr), allocatable            :: CHLDMAPL(:)

!=============================================================================

!  Begin...

  Iam = "MAPL_GenericFinalize"
  call ESMF_GridCompGet(GC, name=COMP_NAME, RC=STATUS )
  VERIFY_(STATUS)
  Iam = trim(COMP_NAME) // Iam


! Retrieve the pointer to the state
!----------------------------------

  call MAPL_InternalStateRetrieve(GC, STATE, RC=STATUS)
  VERIFY_(STATUS)

! Finalize the children
! ---------------------

  call MAPL_GenericStateClockOn(STATE,"TOTAL")
  call MAPL_GenericStateClockOn(STATE,"GenFinalTot")
   if(associated(STATE%GCS)) then
      NC = size(STATE%GCS)
      allocate(CHLDMAPL(NC), stat=status)
      MAXPHASES = 0
      do I=1,NC
         call MAPL_GetObjectFromGC(STATE%GCS(I), CHLDMAPL(I)%PTR, RC=STATUS)
         VERIFY_(STATUS)
         MAXPHASES = MAX(MAXPHASES, SIZE(CHLDMAPL(I)%PTR%PHASE_FINAL))
      end do

      do PHASE = 1, MAXPHASES
         do I=1,NC
            NUMPHASES = SIZE(CHLDMAPL(I)%PTR%PHASE_RUN)
            if (PHASE .le. NUMPHASES) then
               call ESMF_GridCompGet( STATE%GCS(I), NAME=CHILD_NAME, RC=STATUS )
               VERIFY_(STATUS)
      
               call MAPL_GenericStateClockOn (STATE,trim(CHILD_NAME))
               call ESMF_GridCompFinalize (STATE%GCS(I), STATE%GIM(I), STATE%GEX(I), &
                    CLOCK, PHASE=CHLDMAPL(I)%PTR%PHASE_FINAL(PHASE), RC=STATUS )
               VERIFY_(STATUS)
               call MAPL_GenericStateClockOff(STATE,trim(CHILD_NAME))
            end if
         enddo
      end do
      deallocate(CHLDMAPL)
   endif

  call MAPL_GenericStateClockOn(STATE,"--GenFinalMine")

  call MAPL_GetResource( STATE, RECFIN, LABEL="RECORD_FINAL:", &
       RC=STATUS )
  final_checkpoint = .true.
  IF (STATUS == ESMF_SUCCESS) then
     IF (RECFIN == "NO")  final_checkpoint = .false.
  END IF

  if (final_checkpoint) then
! Checkpoint the internal state if required.
!------------------------------------------

     call MAPL_GetResource( STATE   , FILENAME,       &
          LABEL="INTERNAL_CHECKPOINT_FILE:", &
          RC=STATUS )
     if(STATUS==ESMF_SUCCESS) then
        call MAPL_GetResource( STATE   , FILETYPE,         &
             default='binary', &
             LABEL="INTERNAL_CHECKPOINT_TYPE:", &
             RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_ESMFStateWriteToFile(STATE%INTERNAL,CLOCK,FILENAME, &
             FILETYPE, RC=STATUS)
        VERIFY_(STATUS)
     endif

! Checkpoint the import state if required.
!----------------------------------------

     call ESMF_ConfigGetattribute( STATE%CF, FILENAME,         &
          LABEL=trim(COMP_NAME) // "_IMPORT_CHECKPOINT_FILE:", &
          RC=STATUS)
     if(STATUS==ESMF_SUCCESS) then
        call ESMF_ConfigGetattribute( STATE%CF, FILETYPE,         &
             default='binary', &
             LABEL=trim(COMP_NAME) // "_IMPORT_CHECKPOINT_TYPE:", &
             RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_ESMFStateWriteToFile(IMPORT,CLOCK,FILENAME, &
             FILETYPE, RC=STATUS)
        VERIFY_(STATUS)
     endif
  end if

  call MAPL_GenericStateClockOff(STATE,"--GenFinalMine")
  call MAPL_GenericStateClockOff(STATE,"GenFinalTot")
  call MAPL_GenericStateClockOff(STATE,"TOTAL")

! Write summary of profiled times
!--------------------------------

  if (.not. MAPL_ProfIsDisabled()) then
     call WRITE_PARALLEL(" ")
     call WRITE_PARALLEL(" Times for "//trim(COMP_NAME))

     call MAPL_ProfWrite(STATE%TIMES,RC=STATUS)
     VERIFY_(STATUS)

     call WRITE_PARALLEL(" ")
  end if

! Clean-up
!---------
!ALT
 call MAPL_GenericStateDestroy (STATE,  RC=STATUS)
 VERIFY_(STATUS)
! call ESMF_StateDestroy        (IMPORT, RC=STATUS)
! VERIFY_(STATUS)
! call ESMF_StateDestroy        (EXPORT, RC=STATUS)
! VERIFY_(STATUS)

  RETURN_(ESMF_SUCCESS)
end subroutine MAPL_GenericFinalize


! !IROUTINE: MAPL_GenericRecord -- Record the component and its children

! !INTERFACE:

  recursive subroutine MAPL_GenericRecord ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! composite gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! import state
    type(ESMF_State),    intent(inout) :: EXPORT ! export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! the clock
    integer, optional,   intent(  out) :: RC     ! Error code:
                                                 ! = 0 all is well
                                                 ! otherwise, error
!EOPI

! LOCAL VARIABLES

  character(len=ESMF_MAXSTR)                  :: IAm
  character(len=ESMF_MAXSTR)                  :: COMP_NAME
  character(len=ESMF_MAXSTR)                  :: CHILD_NAME
  character(len=18)                           :: datestamp ! YYYYMMDD_HHMMz.bin
  integer                                     :: STATUS
  integer                                     :: I
  type (MAPL_MetaComp), pointer               :: STATE
  logical                                     :: record

!=============================================================================

!  Begin...

  Iam = "MAPL_GenericRecord"
  call ESMF_GridCompGet(GC, name=COMP_NAME, RC=STATUS )
  VERIFY_(STATUS)
  Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
!----------------------------------

  call MAPL_InternalStateRetrieve(GC, STATE, RC=STATUS)
  VERIFY_(STATUS)

  call MAPL_GenericStateClockOn(STATE,"TOTAL")
  call MAPL_GenericStateClockOn(STATE,"GenRecordTot")
! Record the children
! ---------------------
  if(associated(STATE%GCS)) then
     do I=1,size(STATE%GCS)
        call ESMF_GridCompGet( STATE%GCS(I), NAME=CHILD_NAME, RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GenericStateClockOn (STATE,trim(CHILD_NAME))
        call ESMF_GridCompFinalize (STATE%GCS(I), STATE%GIM(I), STATE%GEX(I), CLOCK, &
             PHASE=MAPL_RecordPhase, BlockingFlag=ESMF_NONBLOCKING, RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GenericStateClockOff(STATE,trim(CHILD_NAME))
     enddo
  endif

! Do my "own" record
! ------------------
  call MAPL_GenericStateClockOn(STATE,"--GenRecordMine")
  if (associated(STATE%RECORD)) then
     RECORD = .false.
     DO I = 1, size(STATE%RECORD%ALARM)
        if ( ESMF_AlarmIsRinging(STATE%RECORD%ALARM(I), RC=STATUS) ) then
           VERIFY_(STATUS)
           RECORD = .true.
           cycle
        end if
     END DO
     if (record) then
! add timestamp to filename
        call MAPL_DateStampGet(clock, datestamp, status)
        VERIFY_(STATUS)

        I=STATE%RECORD%IMP_LEN
        if (I > 0) then
           STATE%RECORD%IMP_FNAME(I+1:) = '.' // DATESTAMP
        end if

        I=STATE%RECORD%INT_LEN
        if (I > 0) then
           STATE%RECORD%INT_FNAME(I+1:) = '.' // DATESTAMP
        end if

! call the actual record method
        call MAPL_StateRecord (GC, IMPORT, EXPORT, CLOCK, RC=STATUS )
        VERIFY_(STATUS)
     endif
  endif
  call MAPL_GenericStateClockOff(STATE,"--GenRecordMine")

  call MAPL_GenericStateClockOff(STATE,"GenRecordTot")
  call MAPL_GenericStateClockOff(STATE,"TOTAL")


  RETURN_(ESMF_SUCCESS)
end subroutine MAPL_GenericRecord

subroutine MAPL_StateRecord( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! composite gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! import state
    type(ESMF_State),    intent(inout) :: EXPORT ! export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! the clock
    integer, optional,   intent(  out) :: RC     ! Error code:
                                                 ! = 0 all is well
                                                 ! otherwise, error
!EOPI

! LOCAL VARIABLES

  character(len=ESMF_MAXSTR)                  :: IAm
  character(len=ESMF_MAXSTR)                  :: COMP_NAME
  integer                                     :: STATUS

  type (MAPL_MetaComp), pointer               :: STATE

!=============================================================================

!  Begin...

  Iam = "MAPL_StateRecord"
  call ESMF_GridCompGet(GC, name=COMP_NAME, RC=STATUS )
  VERIFY_(STATUS)
  Iam = trim(COMP_NAME) // Iam


! Retrieve the pointer to the state
!----------------------------------

  call MAPL_InternalStateRetrieve(GC, STATE, RC=STATUS)
  VERIFY_(STATUS)

  if (.not.associated(STATE%RECORD)) then
     RETURN_(ESMF_SUCCESS)
  end if

  if (STATE%RECORD%IMP_LEN > 0) then
     call MAPL_ESMFStateWriteToFile(IMPORT, CLOCK, &
                                    STATE%RECORD%IMP_FNAME, &
                                    'BINARY', RC=STATUS)
     VERIFY_(STATUS)
  end if
  
  if (STATE%RECORD%INT_LEN > 0) then
     call MAPL_ESMFStateWriteToFile(STATE%INTERNAL, CLOCK, &
                                    STATE%RECORD%INT_FNAME, &
                                    'BINARY', RC=STATUS)
     VERIFY_(STATUS)
  end if
  

  RETURN_(ESMF_SUCCESS)
end subroutine MAPL_StateRecord

recursive integer function MAPL_AddChildFromMeta(META, NAME, GRID, CONFIGFILE, SS, PARENTGC, RC)
  type(MAPL_MetaComp),           intent(INOUT) :: META
  character(len=*),              intent(IN   ) :: NAME
  type(ESMF_Grid),  optional,    intent(IN   ) :: GRID
  character(len=*), optional,    intent(IN   ) :: CONFIGFILE
  external                                     :: SS
  type(ESMF_GridComp), optional, intent(IN   ) :: parentGC
  integer,           optional  , intent(  OUT) :: rc
  
  character(len=ESMF_MAXSTR)                  :: IAm='MAPL_AddChild'
  integer                                     :: STATUS

  integer                                     :: I
  type(MAPL_MetaComp), pointer                :: CHILD_META
  type(ESMF_GridComp), pointer                :: TMPGCS(:)
  type(ESMF_State)   , pointer                :: TMPGIM(:)
  type(ESMF_State)   , pointer                :: TMPGEX(:)
  character(len=ESMF_MAXSTR), pointer         :: TMPNL(:)

  if (.not.associated(META%GCS)) then
     ! this is the first child to be added
     allocate(META%GCS(0), stat=status)
     VERIFY_(STATUS)
     allocate(META%GIM(0), stat=status)
     VERIFY_(STATUS)
     allocate(META%GEX(0), stat=status)
     VERIFY_(STATUS)
     allocate(META%GCNameList(0), stat=status)
     VERIFY_(STATUS)
  end if
  I = size(META%GCS) + 1
  MAPL_AddChildFromMeta = I
! realloc GCS, gcnamelist 
  allocate(TMPGCS(I), stat=status)
  VERIFY_(STATUS)
  allocate(TMPGIM(I), stat=status)
  VERIFY_(STATUS)
  allocate(TMPGEX(I), stat=status)
  VERIFY_(STATUS)
  allocate(TMPNL(I), stat=status)
  VERIFY_(STATUS)
  TMPGCS(1:I-1) = META%GCS
  TMPGIM(1:I-1) = META%GIM
  TMPGEX(1:I-1) = META%GEX
  TMPNL(1:I-1) = META%GCNameList
  deallocate(META%GCS)
  deallocate(META%GIM)
  deallocate(META%GEX)
  deallocate(META%GCNameList)
  META%GCS => TMPGCS
  META%GIM => TMPGIM
  META%GEX => TMPGEX
  META%GCNameList => TMPNL

  META%GCNameList(I) = trim(NAME)
  if (present(configfile)) then
     META%GCS(I) = ESMF_GridCompCreate   (     &
          NAME   = trim(NAME),                 &
          CONFIGFILE = configfile,             &
          grid = grid,                         &
          RC=STATUS )
     VERIFY_(STATUS)
  else
     META%GCS(I) = ESMF_GridCompCreate   (     &
          NAME   = trim(NAME),                 &
          CONFIG = META%CF,                    &
          grid = grid,                         &
          RC=STATUS )
     VERIFY_(STATUS)
  end if

! Create each child's import/export state
! ----------------------------------

  META%GIM(I) = ESMF_StateCreate (                         & 
       STATENAME = trim(META%GCNameList(I)) // '_Imports', &
       STATETYPE = ESMF_STATE_IMPORT,                      &
       RC=STATUS )
  VERIFY_(STATUS)


  META%GEX(I) = ESMF_StateCreate (                         &
       STATENAME = trim(META%GCNameList(I)) // '_Exports', &
       STATETYPE = ESMF_STATE_EXPORT,                      &
       RC=STATUS )
  VERIFY_(STATUS)

! create MAPL_Meta
  call MAPL_InternalStateCreate ( META%GCS(I), CHILD_META, RC=STATUS)
  VERIFY_(STATUS)

! put parentGC there
  if (present(parentGC)) then
     allocate(CHILD_META%parentGC, stat=status)
     VERIFY_(STATUS)
     CHILD_META%parentGC = parentGC
  end if

  call ESMF_GridCompSetServices ( META%GCS(I), SS, status )
  VERIFY_(STATUS)

  RETURN_(ESMF_SUCCESS)
end function MAPL_AddChildFromMeta

recursive integer function MAPL_AddChild(GC, NAME, SS, RC)
  type(ESMF_GridComp), intent(INOUT) :: GC
  character(len=*),    intent(IN   ) :: NAME
  external                           :: SS
  integer, optional  , intent(  OUT) :: rc
  
  
  character(len=ESMF_MAXSTR)                  :: IAm
  integer                                     :: STATUS

  type(MAPL_MetaComp), pointer                :: META

  Iam = "MAPL_AddChild"
  call MAPL_InternalStateRetrieve(GC, META, RC=status)
  VERIFY_(STATUS)

  MAPL_AddChild = MAPL_AddChildFromMeta(Meta, NAME, SS=SS, PARENTGC = GC, RC=status)
  VERIFY_(STATUS)

  RETURN_(ESMF_SUCCESS)
end function MAPL_AddChild

subroutine MAPL_DateStampGet (clock, DateStamp, rc)
  type (ESMF_Clock)                 :: clock
  character(len=*)        :: DateStamp
  integer, optional                 :: rc
  
  type(ESMF_Time)                   :: currentTime
  character(len=ESMF_MAXSTR)        :: TimeString
  character                         :: String(ESMF_MAXSTR)
  
  character(len=ESMF_MAXSTR)                  :: IAm
  character(len=ESMF_MAXSTR)                  :: COMP_NAME
  integer                                     :: STATUS

  character*4 year
  character*2 month
  character*2 day
  character*2 hour
  character*2 minute
  character*2 second

  equivalence ( string(01),TimeString )
  equivalence ( string(01),year       )
  equivalence ( string(06),month      )
  equivalence ( string(09),day        )
  equivalence ( string(12),hour       )
  equivalence ( string(15),minute     )
  equivalence ( string(18),second     )
  
  Iam = "MAPL_DateStampGet"
  
  call ESMF_ClockGet (clock, currTime=currentTime, rc=status)
  VERIFY_(STATUS)
  call ESMF_TimeGet  (currentTime, timeString=TimeString, rc=status)
  VERIFY_(STATUS)
  
  DateStamp = year // month // day // '_' // hour // minute // 'z.bin'

  RETURN_(ESMF_SUCCESS)
end subroutine MAPL_DateStampGet

!BOPI

! !IROUTINE: MAPL_GenericRunCouplers

! !INTERFACE:

subroutine MAPL_GenericRunCouplers( STATE, CHILD, CLOCK, RC )

! !ARGUMENTS:

  type (MAPL_MetaComp),    intent(INOUT) :: STATE
  integer,                 intent(IN   ) :: CHILD  ! Child Id
  type(ESMF_Clock),        intent(INOUT) :: CLOCK  ! The clock
  integer, optional,       intent(  OUT) :: RC     ! Error code:

!EOPI

! ErrLog Variables


  character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_GenericRunCouplers"
  integer                               :: STATUS
  integer                               :: J

  if(associated(STATE%GCS)) then
     do J=1,size(STATE%GCS)
        if(STATE%CCcreated(CHILD,J)) then
           call ESMF_CplCompRun (STATE%CCS(CHILD,J), STATE%GEX(CHILD), STATE%GIM(J), &
                CLOCK, RC=STATUS )
           VERIFY_(STATUS)
        endif
     enddo
  end if

  RETURN_(ESMF_SUCCESS)

end subroutine MAPL_GenericRunCouplers


!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

!BOPI

! !IROUTINE: MAPL_InternalStateGet

! !DESCRIPTION:
! This is the recommended way of getting the opaque MAPL Generic
! state object from the gridded component (GC). It can be called at any time
! {\em after} {\tt MAPL\_GenericSetServices} has been called on GC.
! Note that you get a pointer to the object. 
! \newline

! !INTERFACE:

  subroutine MAPL_InternalStateGet ( GC, INTERNAL_STATE, RC)

! !ARGUMENTS:
!  
   
    type(ESMF_GridComp),                  intent(IN   ) :: GC ! Gridded component
    type (MAPL_MetaComp),                       pointer :: INTERNAL_STATE
    integer,                    optional, intent(  OUT) :: RC ! Return code
!EOPI

! ErrLog Variables

    character(len=ESMF_MAXSTR)        :: IAm
    character(len=ESMF_MAXSTR)        :: COMP_NAME
    integer                           :: STATUS

! Local variables
! ---------------

    type (MAPL_GenericWrap )          :: WRAP
#if defined(ABSOFT) || defined(sysIRIX64)
    type(MAPL_MetaComp      ), target :: DUMMY
#endif

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "MAPL_InternalStateGet"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

#if defined(ABSOFT) || defined(sysIRIX64)
    WRAP%INTERNAL_STATE => DUMMY
#endif

    call ESMF_UserCompGetInternalState(GC, "M0APL_GenericInternalState", WRAP, STATUS)
    IF (STATUS /= ESMF_SUCCESS) then
       if (present(RC)) then
          RC = ESMF_FAILURE
       end if
       return
    END IF

    INTERNAL_STATE => WRAP%INTERNAL_STATE


    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_InternalStateGet

  subroutine MAPL_InternalStateCreate( GC, INTERNAL_STATE, RC)

    type(ESMF_GridComp),                  intent(INOUT) :: GC ! Gridded component
    type (MAPL_MetaComp),                       pointer :: INTERNAL_STATE
    integer,                    optional, intent(  OUT) :: RC ! Return code

! ErrLog Variables

    character(len=ESMF_MAXSTR)        :: IAm
    character(len=ESMF_MAXSTR)        :: COMP_NAME
    integer                           :: STATUS

! Local variables
! ---------------

    type (MAPL_GenericWrap )          :: WRAP
#if defined(ABSOFT) || defined(sysIRIX64)
    type(MAPL_MetaComp ), target      :: DUMMY
#endif

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "MAPL_InternalStateCreate"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

#if defined(ABSOFT) || defined(sysIRIX64)
    WRAP%INTERNAL_STATE => DUMMY
#endif
    call ESMF_UserCompGetInternalState(GC, "M0APL_GenericInternalState", WRAP, STATUS)
    ASSERT_(STATUS == ESMF_FAILURE)

! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate(INTERNAL_STATE, STAT=STATUS)
    VERIFY_(STATUS)
    WRAP%INTERNAL_STATE => INTERNAL_STATE

! Have ESMF save pointer to the wrapped internal state in the G.C.
! ----------------------------------------------------------------

    call ESMF_UserCompSetInternalState(GC, "M0APL_GenericInternalState", WRAP, STATUS)
    VERIFY_(STATUS)

! Initialize the config and grid in the generic state.
!-----------------------------------------------------

    call ESMF_GridCompGet( GC, CONFIG = INTERNAL_STATE%CF, &
         GRID = INTERNAL_STATE%GRID%ESMFGRID, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_InternalStateCreate


!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================


  subroutine MAPL_GenericStateDestroy (STATE,RC)
    type (MAPL_MetaComp), pointer              :: STATE
    integer, optional,           intent(  OUT) :: rc     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm = "MAPL_GenericStateDestroy"
    integer :: STATUS

    if(associated(STATE)) then
       call MAPL_SunOrbitDestroy    (STATE%ORBIT         ,RC=STATUS)
       VERIFY_(STATUS)

       call MAPL_GenericGridDestroy (STATE%GRID          ,RC=STATUS)
       VERIFY_(STATUS)

       call ESMF_StateValidate(STATE%INTERNAL, RC=STATUS)
       if (STATUS == ESMF_SUCCESS) then
          call ESMF_StateDestroy       (STATE%INTERNAL      ,RC=STATUS)
          VERIFY_(STATUS)
       end if

       call ESMF_StateValidate(STATE%FORCING, RC=STATUS)
       if (STATUS == ESMF_SUCCESS) then
          call ESMF_StateDestroy       (STATE%FORCING      ,RC=STATUS)
          VERIFY_(STATUS)
       end if

!       call MAPL_VarSpecDestroy     (STATE%IMPORT_SPEC   ,RC=STATUS)
!       VERIFY_(STATUS)

!       call MAPL_VarSpecDestroy     (STATE%EXPORT_SPEC   ,RC=STATUS)
!       VERIFY_(STATUS)

!       call MAPL_VarSpecDestroy     (STATE%INTERNAL_SPEC ,RC=STATUS)
!       VERIFY_(STATUS)

       call MAPL_VarSpecDestroy     (STATE%FRCSPEC   ,RC=STATUS)
       VERIFY_(STATUS)
       if(associated(STATE%IMPORT_SPEC)  ) deallocate(STATE%IMPORT_SPEC)
       if(associated(STATE%EXPORT_SPEC)  ) deallocate(STATE%EXPORT_SPEC)
       if(associated(STATE%INTERNAL_SPEC)) deallocate(STATE%INTERNAL_SPEC)
       if(associated(STATE%GCS          )) deallocate(STATE%GCS          )
       if(associated(STATE%GIM          )) deallocate(STATE%GIM          )
       if(associated(STATE%GEX          )) deallocate(STATE%GEX          )
       if(associated(STATE%GCNameList   )) deallocate(STATE%GCNameList   )
       if(associated(STATE%CCS          )) deallocate(STATE%CCS          )
       if(associated(STATE%CIM          )) deallocate(STATE%CIM          )
       if(associated(STATE%CEX          )) deallocate(STATE%CEX          )
       if(associated(STATE%CCCREATED    )) deallocate(STATE%CCCREATED    )
       if(associated(STATE%TIMES        )) deallocate(STATE%TIMES        )

!ALT: still to do: clean LINK, LOCSTREAM, EXCHANGEGRID, RECORD
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GenericStateDestroy



  subroutine MAPL_StateAddSpec(GC, INTERNAL, IMPORT, EXPORT, &
       SHORT_NAME, LONG_NAME,                                   &
       UNITS,  Dims, VLocation,DATATYPE, NUM_SUBTILES,              &
       REFRESH_INTERVAL, AVERAGING_INTERVAL, HALOWIDTH,         &
       PRECISION, DEFAULT, RC                 )

    type (ESMF_GridComp)            , intent(INOUT)   :: GC
    character (len=*)               , intent(IN)      :: SHORT_NAME
    character (len=*)  , optional   , intent(IN)      :: LONG_NAME
    character (len=*)  , optional   , intent(IN)      :: UNITS
    integer            , optional   , intent(IN)      :: DIMS
    integer            , optional   , intent(IN)      :: DATATYPE
    integer            , optional   , intent(IN)      :: VLOCATION
    integer            , optional   , intent(IN)      :: NUM_SUBTILES
    integer            , optional   , intent(IN)      :: REFRESH_INTERVAL
    integer            , optional   , intent(IN)      :: AVERAGING_INTERVAL
    integer            , optional   , intent(IN)      :: HALOWIDTH
    integer            , optional   , intent(IN)      :: PRECISION
    logical            , optional   , intent(IN)      :: INTERNAL, EXPORT, IMPORT
    real               , optional   , intent(IN)      :: DEFAULT
    integer            , optional   , intent(OUT)     :: RC

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_StateAddSpec"
    integer                               :: STATUS

    if(present(INTERNAL))then
       if(INTERNAL) then
          call MAPL_StateAddInternalSpec(GC, SHORT_NAME, LONG_NAME,             &
               UNITS,  Dims, VLocation,DATATYPE,NUM_SUBTILES,                          &
               REFRESH_INTERVAL, AVERAGING_INTERVAL, DEFAULT, PRECISION=PRECISION, RC=STATUS                     )
          VERIFY_(STATUS)
       endif
    endif

    if(present(IMPORT))then
       if(IMPORT) then
          call MAPL_StateAddImportSpec(GC, SHORT_NAME, LONG_NAME,               &
               UNITS,  Dims, VLocation,DATATYPE,NUM_SUBTILES,                          &
               REFRESH_INTERVAL, AVERAGING_INTERVAL, HALOWIDTH, DEFAULT, RC=STATUS         )
          VERIFY_(STATUS)
       endif
    endif

    if(present(EXPORT))then
       if(EXPORT) then
          call MAPL_StateAddExportSpec(GC, SHORT_NAME, LONG_NAME,               &
               UNITS,  Dims, VLocation,DATATYPE,NUM_SUBTILES,                          &
               REFRESH_INTERVAL, AVERAGING_INTERVAL, HALOWIDTH, DEFAULT, RC=STATUS        )
          VERIFY_(STATUS)
       endif
    endif
    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_StateAddSpec

  subroutine MAPL_StateSetSpecAttrib(STATE, NAME,    &
       INTERNAL, IMPORT, EXPORT, FORCING,            &
       REFRESH_INTERVAL, AVERAGING_INTERVAL, RC      )

    type (MAPL_MetaComp)            , intent(INOUT)   :: STATE
    character (len=*)               , intent(IN)      :: NAME
    integer            , optional   , intent(IN)      :: REFRESH_INTERVAL
    integer            , optional   , intent(IN)      :: AVERAGING_INTERVAL
    logical            , optional   , intent(IN)      :: INTERNAL, EXPORT, IMPORT, FORCING
    integer            , optional   , intent(OUT)     :: RC

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_StateSetSpecAttrib"
    integer                               :: STATUS
    type (MAPL_VarSpec ), pointer         :: SPEC=>null()


    if(present(INTERNAL))then
       if(INTERNAL) then
          SPEC=>STATE%INTERNAL_SPEC(MAPL_VarSpecGetIndex(STATE%INTERNAL_SPEC,NAME))
       endif
    endif

    if(present(IMPORT))then
       if(IMPORT) then
          SPEC=>STATE%IMPORT_SPEC(MAPL_VarSpecGetIndex(STATE%IMPORT_SPEC,NAME))
       endif
    endif

    if(present(EXPORT))then
       if(EXPORT) then
          SPEC=>STATE%EXPORT_SPEC(MAPL_VarSpecGetIndex(STATE%EXPORT_SPEC,NAME))
       endif
    endif

    if(present(FORCING))then
       if(FORCING) then
          SPEC=>STATE%FRCSPEC(MAPL_VarSpecGetIndex(STATE%FRCSPEC,NAME))
       endif
    endif

    ASSERT_(associated(SPEC))

    
    call MAPL_VarSpecSet(SPEC,            &
      ACCMLT_INTERVAL=AVERAGING_INTERVAL, &
      COUPLE_INTERVAL=REFRESH_INTERVAL,   &
                               RC=STATUS  )

    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_StateSetSpecAttrib


  subroutine MAPL_StateGetSpecAttrib(STATE, NAME,    &
       INTERNAL, IMPORT, EXPORT, FORCING,            &
       REFRESH_INTERVAL, AVERAGING_INTERVAL, RC      )

    type (MAPL_MetaComp)            , intent(IN)      :: STATE
    character (len=*)               , intent(IN )     :: NAME
    integer            , optional   , intent(OUT)     :: REFRESH_INTERVAL
    integer            , optional   , intent(OUT)     :: AVERAGING_INTERVAL
    logical            , optional   , intent(IN )     :: INTERNAL, EXPORT, IMPORT, FORCING
    integer            , optional   , intent(OUT)     :: RC

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_StateGetSpecAttrib"
    integer                               :: STATUS
    type (MAPL_VarSpec ), pointer         :: SPEC=>null()


    if(present(INTERNAL))then
       if(INTERNAL) then
          SPEC=>STATE%INTERNAL_SPEC(MAPL_VarSpecGetIndex(STATE%INTERNAL_SPEC,NAME))
       endif
    endif

    if(present(IMPORT))then
       if(IMPORT) then
          SPEC=>STATE%IMPORT_SPEC(MAPL_VarSpecGetIndex(STATE%IMPORT_SPEC,NAME))
       endif
    endif

    if(present(EXPORT))then
       if(EXPORT) then
          SPEC=>STATE%EXPORT_SPEC(MAPL_VarSpecGetIndex(STATE%EXPORT_SPEC,NAME))
       endif
    endif

    if(present(FORCING))then
       if(FORCING) then
          SPEC=>STATE%FRCSPEC(MAPL_VarSpecGetIndex(STATE%FRCSPEC,NAME))
       endif
    endif

    ASSERT_(associated(SPEC))

    
    call MAPL_VarSpecGet(SPEC,            &
      ACCMLT_INTERVAL=AVERAGING_INTERVAL, &
      COUPLE_INTERVAL=REFRESH_INTERVAL,   &
                               RC=STATUS  )

    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_StateGetSpecAttrib

  subroutine MAPL_InternalStateRetrieve(GC, INTERNAL_STATE, RC)

! !ARGUMENTS:
!  
    type(ESMF_GridComp),                  intent(INOUT) :: GC ! Gridded component
    type (MAPL_MetaComp),                       pointer :: INTERNAL_STATE
    integer,                    optional, intent(  OUT) :: RC ! Return code
!EOPI

! ErrLog Variables

    character(len=ESMF_MAXSTR)        :: IAm="MAPL_InternalStateRetrieve"
    integer                           :: STATUS

! Local variables
! ---------------

    call MAPL_InternalStateGet( GC, INTERNAL_STATE, RC=STATUS)
    if (STATUS /= ESMF_SUCCESS) then
       call MAPL_InternalStateCreate( GC, INTERNAL_STATE, RC=STATUS)
       VERIFY_(STATUS)
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_InternalStateRetrieve

!BOPI

! !IROUTINE: MAPL_StateAddImportSpec

! !DESCRIPTION:

!  Sets the specifications for an item in the {\tt IMPORT} state.
!  newline

! !INTERFACE:



  subroutine MAPL_StateAddImportSpec(GC, SHORT_NAME, LONG_NAME,             &
                                     UNITS,  Dims, VLocation,DATATYPE,NUM_SUBTILES,&
                                     REFRESH_INTERVAL, AVERAGING_INTERVAL, &
                                     HALOWIDTH, DEFAULT, RC)

! !ARGUMENTS:

    type (ESMF_GridComp)            , intent(INOUT)   :: GC
    character (len=*)               , intent(IN)      :: SHORT_NAME
    character (len=*)  , optional   , intent(IN)      :: LONG_NAME
    character (len=*)  , optional   , intent(IN)      :: UNITS
    integer            , optional   , intent(IN)      :: DIMS
    integer            , optional   , intent(IN)      :: DATATYPE
    integer            , optional   , intent(IN)      :: NUM_SUBTILES
    integer            , optional   , intent(IN)      :: VLOCATION
    integer            , optional   , intent(IN)      :: REFRESH_INTERVAL
    integer            , optional   , intent(IN)      :: AVERAGING_INTERVAL
    integer            , optional   , intent(IN)      :: HALOWIDTH
    real               , optional   , intent(IN)      :: DEFAULT
    integer            , optional   , intent(OUT)     :: RC

!EOPI

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_StateAddImportSpec"
    integer                               :: STATUS
    integer                               :: usable_AI
    integer                               :: usable_RI
    real                                  :: dt
    type (ESMF_Config)                    :: CF
    type (MAPL_MetaComp), pointer         :: STATE

    call MAPL_InternalStateRetrieve(GC, STATE, RC=status)
    VERIFY_(STATUS)

    CF = STATE%CF

!  Get the clock increment interval
    call ESMF_ConfigGetAttribute( CF, dt,  Label="RUN_DT:", RC=STATUS)
    VERIFY_(STATUS)

    if (present(REFRESH_INTERVAL)) then
       usable_RI = REFRESH_INTERVAL
    else
       usable_RI = nint(dt)
    endif

    if (present(AVERAGING_INTERVAL)) then
       usable_AI = AVERAGING_INTERVAL
    else
       usable_AI = nint(dt)
    endif

    call MAPL_VarSpecCreateInList(STATE%IMPORT_SPEC,                         &
       LONG_NAME  = LONG_NAME,                                               &
       UNITS      = UNITS,                                                   &
       SHORT_NAME = SHORT_NAME,                                              &
       DIMS       = DIMS,                                                    &
       STAT       = DATATYPE,                                                &
       NUM_SUBTILES=NUM_SUBTILES,                                            &
       ACCMLT_INTERVAL= usable_AI,                                           &
       COUPLE_INTERVAL= usable_RI,                                           &
       VLOCATION  = VLOCATION,                                               &
       HALOWIDTH  = HALOWIDTH,                                               &
       DEFAULT    = DEFAULT, RC=STATUS  )
    VERIFY_(STATUS)


    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_StateAddImportSpec

  subroutine MAPL_StateAddImportSpecFrmChld ( STATE, &
                                              SHORT_NAME, CHILD_ID, RC  )
    type (MAPL_MetaComp)            , intent(INOUT)   :: STATE
    character (len=*)               , intent(IN)      :: SHORT_NAME
    integer                         , intent(IN)      :: CHILD_ID
    integer            , optional   , intent(OUT)     :: RC


    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_StateAddImportSpecFrmChld"
    integer                               :: STATUS
    type (MAPL_VarSpec),      pointer     :: SPECS(:)
    integer                               :: I



    if (.not. associated(STATE%GCS)) then
       RETURN_(ESMF_FAILURE)
    end if

    call MAPL_GridCompGetVarSpecs(STATE%GCS(CHILD_ID), IMPORT=SPECS, RC=STATUS)
    VERIFY_(STATUS)

    I=MAPL_VarSpecGetIndex(SPECS, SHORT_NAME, RC=STATUS)

    if (I == -1) then
       RETURN_(ESMF_FAILURE)
    endif
    VERIFY_(STATUS)
    
    call MAPL_VarSpecAddRefToList(STATE%IMPORT_SPEC, SPECS(I), RC=STATUS)
    if (STATUS /= MAPL_DuplicateEntry) then
       VERIFY_(STATUS)
    else
       RETURN_(ESMF_FAILURE) ! ALT this needs to be revisited
    endif

!ALT: is reconnect needed
!    call MAPL_AddConnectivity ( STATE, SHORT_NAME=SHORT_NAME, &
!       FROM_IMPORT=CHILD_ID, TO_IMPORT=MAPL_Self, RC=STATUS  )
!    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_StateAddImportSpecFrmChld

!BOPI

! !IROUTINE: MAPL_StateAddExportSpec

! !DESCRIPTION:

!  Sets the specifications for an item in the {\tt EXPORT} state.
!  newline

! !INTERFACE:


  subroutine MAPL_StateAddExportSpec(GC, SHORT_NAME, LONG_NAME,             &
                                     UNITS,  Dims, VLocation,DATATYPE,NUM_SUBTILES,&
                                     REFRESH_INTERVAL, AVERAGING_INTERVAL, &
				     HALOWIDTH, DEFAULT, RC  )

! !ARGUMENTS:

    type (ESMF_GridComp)            , intent(INOUT)   :: GC
    character (len=*)               , intent(IN)      :: SHORT_NAME
    character (len=*)  , optional   , intent(IN)      :: LONG_NAME
    character (len=*)  , optional   , intent(IN)      :: UNITS
    integer            , optional   , intent(IN)      :: DIMS
    integer            , optional   , intent(IN)      :: DATATYPE
    integer            , optional   , intent(IN)      :: VLOCATION
    integer            , optional   , intent(IN)      :: NUM_SUBTILES
    integer            , optional   , intent(IN)      :: REFRESH_INTERVAL
    integer            , optional   , intent(IN)      :: AVERAGING_INTERVAL
    integer            , optional   , intent(IN)      :: HALOWIDTH
    real               , optional   , intent(IN)      :: DEFAULT
    integer            , optional   , intent(OUT)     :: RC


!EOPI


    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_StateAddExportSpec"
    integer                               :: STATUS
    integer                               :: usable_AI
    integer                               :: usable_RI
    real                                  :: usable_Default
    real                                  :: dt
    type (ESMF_Config)                    :: CF
    type (MAPL_MetaComp), pointer         :: STATE

    call MAPL_InternalStateRetrieve(GC, STATE, RC=status)
    VERIFY_(STATUS)

    CF = STATE%CF

!  Get the clock increment interval
    call ESMF_ConfigGetAttribute( CF, dt,  Label="RUN_DT:", RC=STATUS)
    VERIFY_(STATUS)

    if (present(REFRESH_INTERVAL)) then
       usable_RI = REFRESH_INTERVAL
    else
       usable_RI = nint(dt)
    endif

    if (present(AVERAGING_INTERVAL)) then
       usable_AI = AVERAGING_INTERVAL
    else
       usable_AI = nint(dt)
    endif

    call MAPL_VarSpecCreateInList(STATE%EXPORT_SPEC,                         &
       LONG_NAME  = LONG_NAME,                                               &
       UNITS      = UNITS,                                                   &
       SHORT_NAME = SHORT_NAME,                                              &
       DIMS       = DIMS,                                                    &
       STAT       = DATATYPE,                                                &
       NUM_SUBTILES=NUM_SUBTILES,                                            &
       ACCMLT_INTERVAL= usable_AI,                                           &
       COUPLE_INTERVAL= usable_RI,                                           &
       VLOCATION  = VLOCATION,                                               &
       HALOWIDTH  = HALOWIDTH,                                               &
       DEFAULT    = DEFAULT, RC=STATUS  )
    VERIFY_(STATUS)


    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_StateAddExportSpec



!BOPI

! !IROUTINE: MAPL_StateAddInternalSpec

! !DESCRIPTION:

!  Sets the specifications for an item in the {\tt INTERNAL} state.
!  newline

! !INTERFACE:

  subroutine MAPL_StateAddInternalSpec(GC,              &
                                       SHORT_NAME,         &
                                       LONG_NAME,          &
                                       UNITS,              &
                                       DIMS,               &
                                       VLOCATION,          &
                                       DATATYPE,           &
                                       NUM_SUBTILES,       &
                                       REFRESH_INTERVAL,   &
                                       AVERAGING_INTERVAL, &
                                       DEFAULT,            &
                                       NORESTART,          &
                                       HALOWIDTH,          &
                                       PRECISION,          &
                                       FRIENDLYTO,         &
                                       ADD2EXPORT,         &
                                       ATTR_RNAMES,        &
                                       ATTR_INAMES,        &
                                       ATTR_RVALUES,       &
                                       ATTR_IVALUES,       &
                                                       RC  )

! !ARGUMENTS:

    type (ESMF_GridComp)            , intent(INOUT)   :: GC
    character (len=*)               , intent(IN)      :: SHORT_NAME
    character (len=*)  , optional   , intent(IN)      :: LONG_NAME
    character (len=*)  , optional   , intent(IN)      :: UNITS
    integer            , optional   , intent(IN)      :: DIMS
    integer            , optional   , intent(IN)      :: DATATYPE
    integer            , optional   , intent(IN)      :: VLOCATION
    integer            , optional   , intent(IN)      :: NUM_SUBTILES
    integer            , optional   , intent(IN)      :: REFRESH_INTERVAL
    integer            , optional   , intent(IN)      :: AVERAGING_INTERVAL
    integer            , optional   , intent(IN)      :: PRECISION
    real               , optional   , intent(IN)      :: DEFAULT
    logical            , optional   , intent(IN)      :: NORESTART
    character (len=*)  , optional   , intent(IN)      :: HALOWIDTH
    character (len=*)  , optional   , intent(IN)      :: FRIENDLYTO
    logical            , optional   , intent(IN)      :: ADD2EXPORT
    character (len=*)  , optional   , intent(IN)      :: ATTR_INAMES(:)
    character (len=*)  , optional   , intent(IN)      :: ATTR_RNAMES(:)
    integer            , optional   , intent(IN)      :: ATTR_IVALUES(:)
    real               , optional   , intent(IN)      :: ATTR_RVALUES(:)
    integer            , optional   , intent(OUT)     :: RC

!EOPI

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_StateAddInternalSpec"
    integer                               :: STATUS
    logical                               :: usable_NR
    integer                               :: usable_HW
    integer                               :: I
    type (MAPL_MetaComp), pointer         :: STATE

    call MAPL_InternalStateRetrieve(GC, STATE, RC=status)
    VERIFY_(STATUS)

    if (present(NoRestart)) then
       usable_NR  = NoRestart
    else
       usable_NR = .false.
    endif

    if (present(HALOWIDTH)) then
       read(HALOWIDTH,'(I1)') usable_HW
    else
       usable_HW = 0
    endif

!ALT: NoRestart is still stubbed (more work is needed to support it)

    call MAPL_VarSpecCreateInList(STATE%INTERNAL_SPEC,                       &
       LONG_NAME  = LONG_NAME,                                               &
       UNITS      = UNITS,                                                   &
       SHORT_NAME = SHORT_NAME,                                              &
       DIMS       = DIMS,                                                    &
       STAT       = DATATYPE,                                                &
       NUM_SUBTILES=NUM_SUBTILES,                                            &
       ACCMLT_INTERVAL= AVERAGING_INTERVAL,                                  &
       COUPLE_INTERVAL= REFRESH_INTERVAL,                                    &
       VLOCATION  = VLOCATION,                                               &
       DEFAULT    = DEFAULT, FRIENDLYTO = FRIENDLYTO,                        &
       HALOWIDTH  = usable_HW, PRECISION=PRECISION,                          &
       ATTR_RNAMES=ATTR_RNAMES, ATTR_INAMES=ATTR_INAMES,                     &
       ATTR_RVALUES=ATTR_RVALUES, ATTR_IVALUES=ATTR_IVALUES,                 &
       RC=STATUS  )
    VERIFY_(STATUS)

!ALT: the next section is added here upon the request of Arlindo:
!     if FRIENDLYTO is set, we automatically 
!     add the spec/field to the export

    if (present(FRIENDLYTO) .or. present(ADD2EXPORT)) then

       I=MAPL_VarSpecGetIndex(STATE%INTERNAL_SPEC, SHORT_NAME, RC=STATUS)
       if (I == -1) then
          RETURN_(ESMF_FAILURE)
       endif
       VERIFY_(STATUS)
    
       call MAPL_VarSpecAddRefToList(STATE%EXPORT_SPEC, STATE%INTERNAL_SPEC(I), RC=STATUS)
       VERIFY_(STATUS)

    endif

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_StateAddInternalSpec




  subroutine MAPL_StateAddExportSpecFrmChld ( GC, &
                                              SHORT_NAME, CHILD_ID, RC  )
    type(ESMF_GridComp),            intent(INOUT) :: GC 
    character (len=*)               , intent(IN)      :: SHORT_NAME
    integer                         , intent(IN)      :: CHILD_ID
    integer            , optional   , intent(OUT)     :: RC


    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_StateAddExportSpecFrmChld"
    integer                               :: STATUS
    type (MAPL_VarSpec),      pointer     :: SPECS(:)
    integer                               :: I


    call MAPL_AddConnectivityE2E ( GC, SHORT_NAME, &
                                SRC_ID = CHILD_ID, &
                                TO_EXPORT = MAPL_Self, RC=STATUS  )
    VERIFY_(STATUS)


    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_StateAddExportSpecFrmChld

  subroutine MAPL_AddConnectivityOld ( GC, SHORT_NAME, TO_NAME, &
       FROM_IMPORT, FROM_EXPORT, TO_IMPORT, TO_EXPORT, RC  )

    type(ESMF_GridComp),            intent(INOUT) :: GC ! Gridded component
    character (len=*)             , intent(IN   ) :: SHORT_NAME
    character (len=*),    optional, intent(IN   ) :: TO_NAME
    integer,              optional, intent(IN   ) :: FROM_IMPORT
    integer,              optional, intent(IN   ) :: FROM_EXPORT
    integer,              optional, intent(IN   ) :: TO_IMPORT
    integer,              optional, intent(IN   ) :: TO_EXPORT
    integer,              optional, intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_AddConnectivity"
    integer                               :: STATUS
    type (MAPL_ConnectivityWrap)          :: connwrap
    type (MAPL_Connectivity), pointer     :: conn


    call ESMF_UserCompGetInternalState(gc, 'M0APL_Connectivity', &
                                       connwrap, status)
    if (STATUS == ESMF_FAILURE) then
       allocate(conn, STAT=STATUS)
       VERIFY_(STATUS)
       connwrap%ptr => conn
       call ESMF_UserCompSetInternalState(gc, 'M0APL_Connectivity', &
                                          connwrap, status)
       VERIFY_(STATUS)
    else
       conn => connwrap%ptr
    end if

    call MAPL_VarConnCreate(CONN%CONNECT, SHORT_NAME, TO_NAME=TO_NAME,        &
                            FROM_IMPORT=FROM_IMPORT, FROM_EXPORT=FROM_EXPORT, &
                            TO_IMPORT  =TO_IMPORT,   TO_EXPORT  =TO_EXPORT, RC=STATUS  )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_AddConnectivityOld

  subroutine MAPL_AddConnectivityE2E ( GC, SHORT_NAME, &
       SRC_ID, TO_EXPORT, RC  )

    type(ESMF_GridComp),            intent(INOUT) :: GC ! Gridded component
    character (len=*),              intent(IN   ) :: SHORT_NAME
    integer,                        intent(IN   ) :: SRC_ID !FROM_EXPORT
    integer,                        intent(IN   ) :: TO_EXPORT
    integer,              optional, intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_AddConnectivityE2E"
    integer                               :: STATUS
    type (MAPL_ConnectivityWrap)          :: connwrap
    type (MAPL_Connectivity), pointer     :: conn


    call ESMF_UserCompGetInternalState(gc, 'M0APL_Connectivity', &
                                       connwrap, status)
    if (STATUS == ESMF_FAILURE) then
       allocate(conn, STAT=STATUS)
       VERIFY_(STATUS)
       connwrap%ptr => conn
       call ESMF_UserCompSetInternalState(gc, 'M0APL_Connectivity', &
                                          connwrap, status)
       VERIFY_(STATUS)
    else
       conn => connwrap%ptr
    end if

    call MAPL_VarConnCreate(CONN%CONNECT, SHORT_NAME, &
                            FROM_EXPORT=SRC_ID, &
                            TO_EXPORT=TO_EXPORT, RC=STATUS  )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_AddConnectivityE2E

  subroutine MAPL_AddConnectivityRename ( GC, SRC_NAME, SRC_ID, &
       DST_NAME, DST_ID, RC  )

    type(ESMF_GridComp),            intent(INOUT) :: GC ! Gridded component
    character (len=*),              intent(IN   ) :: SRC_NAME !FROM_NAME==SHORT_NAME
    character (len=*),              intent(IN   ) :: DST_NAME !TO_NAME
    integer,                        intent(IN   ) :: SRC_ID !FROM_EXPORT
    integer,                        intent(IN   ) :: DST_ID !TO_IMPORT
    integer,              optional, intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_AddConnectivityRename"
    integer                               :: STATUS
    type (MAPL_ConnectivityWrap)          :: connwrap
    type (MAPL_Connectivity), pointer     :: conn


    call ESMF_UserCompGetInternalState(gc, 'M0APL_Connectivity', &
                                       connwrap, status)
    if (STATUS == ESMF_FAILURE) then
       allocate(conn, STAT=STATUS)
       VERIFY_(STATUS)
       connwrap%ptr => conn
       call ESMF_UserCompSetInternalState(gc, 'M0APL_Connectivity', &
                                          connwrap, status)
       VERIFY_(STATUS)
    else
       conn => connwrap%ptr
    end if

    call MAPL_VarConnCreate(CONN%CONNECT, SHORT_NAME=SRC_NAME, TO_NAME=DST_NAME,        &
                            FROM_EXPORT=SRC_ID, TO_IMPORT=DST_ID, RC=STATUS  )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_AddConnectivityRename

  subroutine MAPL_AddConnectivityRenameMany ( GC, SRC_NAME, SRC_ID, &
       DST_NAME, DST_ID, RC  )

    type(ESMF_GridComp),            intent(INOUT) :: GC ! Gridded component
    character (len=*),              intent(IN   ) :: SRC_NAME(:)
    character (len=*),              intent(IN   ) :: DST_NAME(:)
    integer,                        intent(IN   ) :: SRC_ID !FROM_EXPORT
    integer,                        intent(IN   ) :: DST_ID !TO_IMPORT
    integer,              optional, intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_AddConnectivityRename"
    integer                               :: STATUS
    integer                               :: I

    DO I = 1, size(SRC_NAME)
       call MAPL_AddConnectivity ( GC, SRC_NAME=SRC_NAME(I), DST_NAME=DST_NAME(I), &
                                   SRC_ID=SRC_ID, DST_ID=DST_ID, RC=status  )
       VERIFY_(STATUS)
    END DO

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_AddConnectivityRenameMany

  subroutine MAPL_AddConnectivityMany ( GC, SHORT_NAME, SRC_ID, DST_ID, RC  )

    type(ESMF_GridComp),            intent(INOUT) :: GC ! Gridded component
    character (len=*)             , intent(IN   ) :: SHORT_NAME(:)
    integer,                        intent(IN   ) :: SRC_ID
    integer,                        intent(IN   ) :: DST_ID
    integer,              optional, intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_AddConnectivityMany"
    integer                               :: STATUS
    integer                               :: I

    DO I = 1, size(SHORT_NAME)
!ALT attention: DST_NAME needs to be revisited once we remove the old style interface to
!               MAPL_AddConnectivity
       call MAPL_AddConnectivityRENAME ( GC, SRC_NAME=SHORT_NAME(I), DST_NAME=SHORT_NAME(I), &
                                   SRC_ID=SRC_ID, DST_ID=DST_ID, RC=status  )
       VERIFY_(STATUS)
    END DO


    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_AddConnectivityMany

  subroutine MAPL_DoNotConnect ( GC, SHORT_NAME, CHILD, RC )

    type(ESMF_GridComp),            intent(INOUT) :: GC ! Gridded component
    character (len=*)             , intent(IN   ) :: SHORT_NAME
    integer,                        intent(IN   ) :: CHILD
    integer,              optional, intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_DoNotConnect"
    integer                               :: STATUS
    type (MAPL_MetaComp), pointer         :: STATE
    type (MAPL_ConnectivityWrap)          :: connwrap
    type (MAPL_Connectivity), pointer     :: conn


    call ESMF_UserCompGetInternalState(gc, 'M0APL_Connectivity', &
                                       connwrap, status)
    if (STATUS == ESMF_FAILURE) then
       allocate(conn, STAT=STATUS)
       VERIFY_(STATUS)
       connwrap%ptr => conn
       call ESMF_UserCompSetInternalState(gc, 'M0APL_Connectivity', &
                                          connwrap, status)
       VERIFY_(STATUS)
    else
       conn => connwrap%ptr
    end if

    call MAPL_VarConnCreate(CONN%DONOTCONN, SHORT_NAME, &
       FROM_IMPORT=CHILD, RC=STATUS  )
    VERIFY_(STATUS)


    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_DoNotConnect

  subroutine MAPL_DoNotConnectMany ( GC, SHORT_NAME, CHILD, RC )

    type(ESMF_GridComp),            intent(INOUT) :: GC ! Gridded component
    character (len=*)             , intent(IN   ) :: SHORT_NAME(:)
    integer,                        intent(IN   ) :: CHILD
    integer,              optional, intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_DoNotConnectMany"
    integer                               :: STATUS
    integer                               :: I


    DO I = 1, size(SHORT_NAME)
       call MAPL_DoNotConnect(GC, SHORT_NAME(I), CHILD, RC=status)
       VERIFY_(STATUS)
    END DO

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_DoNotConnectMany

  subroutine MAPL_DoNotConnectAnyImport ( GC, CHILD, RC )

    type(ESMF_GridComp),            intent(INOUT) :: GC ! Gridded component
    integer,                        intent(IN   ) :: CHILD
    integer,              optional, intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_DoNotConnectAnyImport"
    integer                               :: STATUS

    call MAPL_DoNotConnect ( GC, SHORT_NAME="MAPL_AnyChildImport", &
                             CHILD=CHILD, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_DoNotConnectAnyImport


  subroutine MAPL_TerminateImportAll ( GC, ALL, RC )

    type(ESMF_GridComp),            intent(INOUT) :: GC ! Gridded component
    logical,                        intent(IN   ) :: ALL
    integer,              optional, intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_TerminateImportAll"
    integer                               :: STATUS
    type(MAPL_MetaComp), pointer          :: META
    integer                               :: I

    call MAPL_GetObjectFromGC(GC, META, RC=STATUS)
    VERIFY_(STATUS)

    if (associated(META%GCS)) then
       do I=1, size(META%GCS)
          call MAPL_TerminateImport ( GC, CHILD=I, RC=STATUS )
          VERIFY_(STATUS)
       end do
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_TerminateImportAll


  subroutine MAPL_StateAddExportSpecFrmAll ( STATE, RC  )
    type (MAPL_MetaComp)            , intent(INOUT)   :: STATE
    integer            , optional   , intent(OUT)     :: RC


    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_StateAddExportSpecFrmAll"
    integer                               :: STATUS
    type (MAPL_VarSpec),      pointer     :: SPECS(:)
    type (MAPL_VarSpec),      pointer     :: MYSPEC(:) => null()
    integer                               :: I
    integer                               :: N
    character(len=ESMF_MAXSTR)            :: NAME
    character(len=ESMF_MAXSTR)            :: GCNAME
    character (len=ESMF_MAXSTR)           :: LONG_NAME
    character (len=ESMF_MAXSTR)           :: UNITS
    integer                               :: DIMS
    integer                               :: VLOCATION
    integer                               :: NUM_SUBTILES
    integer                               :: ACCMLT_INTERVAL
    integer                               :: COUPLE_INTERVAL
    integer                               :: STAT
    type(ESMF_Field), pointer             :: FIELD



    if (.not. associated(STATE%GCS)) then
       RETURN_(ESMF_FAILURE)
    end if

    do N = 1, size(STATE%GCS)

       call MAPL_GridCompGetVarSpecs(STATE%GCS(N), EXPORT=SPECS, RC=STATUS)
       VERIFY_(STATUS)

       call ESMF_GridCompGet( STATE%GCS(N), name=GCNAME, RC=STATUS )
       VERIFY_(STATUS)

       
       do I = 1, size(SPECS)
    
          call MAPL_VarSpecGet(SPECS(I), SHORT_NAME=NAME,  &
                 LONG_NAME  = LONG_NAME,                                  &
                 UNITS      = UNITS,                                      &
                 DIMS       = DIMS,                                       &
                 VLOCATION  = VLOCATION,                                  &
                 NUM_SUBTILES=NUM_SUBTILES,                               &
                 STAT       = STAT,                                       &
                 ACCMLT_INTERVAL= ACCMLT_INTERVAL,                     &
                 COUPLE_INTERVAL= COUPLE_INTERVAL,                       &
                 RC=STATUS  )
          VERIFY_(STATUS)

          call MAPL_VarSpecGet(SPECS(I), FIELDPTR = FIELD, RC=STATUS  )
          VERIFY_(STATUS)


          call MAPL_VarSpecCreateInList(MYSPEC,                         &
                 SHORT_NAME=trim(NAME)//'_from_' // trim(GCNAME),         &
                 LONG_NAME  = LONG_NAME,                                  &
                 UNITS      = UNITS,                                      &
                 DIMS       = DIMS,                                       &
                 VLOCATION  = VLOCATION,                                  &
                 NUM_SUBTILES=NUM_SUBTILES,                               &
                 STAT       = STAT,                                       &
                 ACCMLT_INTERVAL= ACCMLT_INTERVAL,                     &
                 COUPLE_INTERVAL= COUPLE_INTERVAL,                       &
                 RC=STATUS  )
          VERIFY_(STATUS)

          call MAPL_VarSpecSet(MYSPEC(1), FIELDPTR = FIELD, RC=STATUS  )
          VERIFY_(STATUS)

          call MAPL_VarSpecAddRefToList(STATE%EXPORT_SPEC, MYSPEC(1), RC=STATUS)
          if (STATUS /= MAPL_DuplicateEntry) then
             VERIFY_(STATUS)
          else
             RETURN_(ESMF_FAILURE) ! ALT this needs to be revisited
          endif

          NULLIFY(MYSPEC)
!ALT specDestroy ?; otherwise we have small mem leak
!ALT    call MAPL_AddConnectivity ???


       end do

    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_StateAddExportSpecFrmAll



!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================



!BOPI

! !IROUTINE: MAPL_GenericStateGet

! !DESCRIPTION:
! This is the way of querying the opaque {\em MAPL\_Generic}
! state object. The arguments are:
!  \begin{description}
!   \item[STATE]
!    The \ggn\  object to be queried.     
!   \item[IM]
!      Size of the first horizontal dimension (X) of local arrays.
!   \item[JM]
!      Size of the second horizontal dimension (Y) of local arrays.
!   \item[LM]
!      Size of the vertical dimension.
!   \item[VERTDIM]
!      Position of the vertical dimension of 2 or higher dimensional arrays.
!   \item[NX]
!      Size of the DE array dimension aligned with the first horizontal dimension of arrays 
!   \item[NY]
!      Size of the DE array dimension aligned with the second horizontal dimension of arrays 
!   \item[NX0, NY0]
!      Coordinates of current DE.
!   \item[LONS]
!      X coordinates of array locations. Currently longitude in radians.
!   \item[LATS]
!      Y coordinates of array locations. Currently latitude in radians.
!   \item[INTERNAL\_ESMF\_STATE]
!      The gridded component's INTERNAL state.
!   \item[GCNames]
!      Names of the children.
!   \item[GCS]
!      The child gridded components.
!   \item[GIM]
!      The childrens' IMPORT states.
!   \item[GEX]
!      The childrens' EXPORT states.
!   \item[CCS]
!      Array of child-to-child couplers.
!   \end{description}

! !INTERFACE:


  subroutine MAPL_GenericStateGet (STATE, IM, JM, LM, VERTDIM,                &
                                   NX, NY, NX0, NY0, LAYOUT,                  &
                                   GCNames,                                   &
                                   LONS, LATS, ORBIT, RUNALARM,               &
                                   IMPORTspec, EXPORTspec, INTERNALspec,      &
                                   INTERNAL_ESMF_STATE,                       &
                                   TILETYPES, TILEKIND,                       &
                                   TILELATS,TILELONS,LOCSTREAM,               &
                                   EXCHANGEGRID,                               &
                                   CLOCK,                                      &
                                   GCS, CCS, GIM, GEX, CF, RC   )

! !ARGUMENTS:

    type (MAPL_MetaComp),           intent(INOUT) :: STATE
    type (ESMF_Alarm),    optional, intent(  OUT) :: RUNALARM
    type (MAPL_SunOrbit), optional, intent(  OUT) :: ORBIT
    integer,              optional, intent(  OUT) :: IM, JM, LM
    integer,              optional, intent(  OUT) :: VERTDIM
    integer,              optional, intent(  OUT) :: NX, NY, NX0, NY0
    type (ESMF_DELayout), optional, intent(  OUT) :: LAYOUT
    real, pointer,        optional                :: LONS(:,:)
    real, pointer,        optional                :: LATS(:,:)
    integer,              optional, intent(  OUT) :: RC     ! Error code:
    character(len=ESMF_MAXSTR),optional, pointer  :: GCNames(:)
    type (MAPL_VarSpec),  optional, pointer       :: IMPORTspec(:)
    type (MAPL_VarSpec),  optional, pointer       :: EXPORTspec(:)
    type (MAPL_VarSpec),  optional, pointer       :: INTERNALspec(:)
    type (ESMF_State),    optional, intent(  OUT) :: INTERNAL_ESMF_STATE
    integer,              optional, pointer       :: TILETYPES(:)
    integer,              optional, pointer       :: TILEKIND(:)
    real, pointer,        optional                :: TILELONS(:)
    real, pointer,        optional                :: TILELATS(:)
    type (MAPL_LocStream),optional, intent(  OUT) :: LOCSTREAM
    type (MAPL_LocStream),optional, intent(  OUT) :: EXCHANGEGRID
    type (ESMF_CLOCK)    ,optional, intent(  OUT) :: CLOCK
    type (ESMF_GridComp), optional, pointer       :: GCS(:)
    type (ESMF_CplComp),  optional, pointer       :: CCS(:,:)
    type (ESMF_State),    optional, pointer       :: GIM(:)
    type (ESMF_State),    optional, pointer       :: GEX(:)
    type (ESMF_Config),   optional, intent(  OUT) :: CF

!EOPI

    character(len=ESMF_MAXSTR), parameter :: IAm = "MAPL_GenericStateGet"
    integer                               :: STATUS

    real                                  :: ECC
    real                                  :: OB
    real                                  :: PER
    integer                               :: EQNX

     if(present(IM)) then
      IM=STATE%GRID%IM
     endif

     if(present(JM)) then
      JM=STATE%GRID%JM
     endif

     if(present(LM)) then
      LM=STATE%GRID%LM
     endif

     if(present(VERTDIM)) then
      VERTDIM=STATE%GRID%VERTDIM
     endif

     if(present(NX)) then
      NX=STATE%GRID%NX
     endif

     if(present(NY)) then
      NY=STATE%GRID%NY
     endif

     if(present(NX0)) then
      NX0=STATE%GRID%NX0
     endif

     if(present(NY0)) then
      NY0=STATE%GRID%NY0
     endif

     if(present(LAYOUT)) then
      LAYOUT=STATE%GRID%LAYOUT
     endif

     if(present(CF)) then
      CF=STATE%CF
     endif

     if(present(ORBIT)) then

        if(.not.MAPL_SunOrbitCreated(STATE%ORBIT)) then 

           call MAPL_GetResource(STATE, ECC, Label="ECCENTRICITY:", default=0.0167, &
                  RC=STATUS)
           VERIFY_(STATUS)

           call MAPL_GetResource(STATE, OB, Label="OBLIQUITY:"   , default=23.45 , &
                RC=STATUS)
           VERIFY_(STATUS)

           call MAPL_GetResource(STATE, PER, Label="PERIHELION:"  , default=102.0 , &
                RC=STATUS)
           VERIFY_(STATUS)

           call MAPL_GetResource(STATE, EQNX, Label="EQUINOX:"     , default=80    , &
                RC=STATUS)
           VERIFY_(STATUS)

           STATE%ORBIT = MAPL_SunOrbitCreate(STATE%CLOCK,ECC,OB,PER,EQNX,RC=STATUS)
           VERIFY_(STATUS)

        end if
        ORBIT=STATE%ORBIT
     end if

     if(present(RUNALARM)) then
      RUNALARM=STATE%ALARM(0)
     endif

     if(present(GCNames )) then
      GCNames=>STATE%GCNameList
     endif

     if(present(LONS    )) then
      LONS   =>STATE%GRID%LONS
     endif

     if(present(LATS    )) then
      LATS   =>STATE%GRID%LATS
     endif

     if(present(IMPORTspec)) then
      IMPORTspec =>STATE%IMPORT_SPEC
     endif

     if(present(EXPORTspec)) then
      EXPORTspec =>STATE%EXPORT_SPEC
     endif

     if(present(INTERNALspec)) then
      INTERNALspec =>STATE%INTERNAL_SPEC
     endif

     if(present(INTERNAL_ESMF_STATE)) then
      INTERNAL_ESMF_STATE = STATE%INTERNAL
     endif

     if(present(TILETYPES)) then
        call MAPL_LocStreamGet(STATE%LocStream, TILETYPE=TILETYPES, RC=STATUS)
        VERIFY_(STATUS)
     end if

     if(present(TILEKIND)) then
        call MAPL_LocStreamGet(STATE%LocStream, TILEKIND=TILEKIND, RC=STATUS)
        VERIFY_(STATUS)
     end if

     if(present(TILELONS)) then
        call MAPL_LocStreamGet(STATE%LocStream, TILELONS=TILELONS, RC=STATUS)
        VERIFY_(STATUS)
     end if

     if(present(TILELATS)) then
        call MAPL_LocStreamGet(STATE%LocStream, TILELATS=TILELATS, RC=STATUS)
        VERIFY_(STATUS)
     end if

     if(present(LOCSTREAM)) then
      LOCSTREAM = STATE%LOCSTREAM
     endif

     if(present(EXCHANGEGRID)) then
      EXCHANGEGRID = STATE%EXCHANGEGRID
     endif

     if(present(CLOCK)) then
      CLOCK = STATE%CLOCK
     endif

     if(present(GCS)) then
      GCS => STATE%GCS
     endif

     if(present(CCS)) then
      CCS => STATE%CCS
     endif

     if(present(GIM)) then
      GIM => STATE%GIM
     endif

     if(present(GEX)) then
      GEX => STATE%GEX
     endif

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GenericStateGet

  subroutine MAPL_GenericStateSet (STATE, ORBIT, LM, RUNALARM, CHILDINIT, &
                                   LOCSTREAM, EXCHANGEGRID, CLOCK, NAME, CF, RC)
    type (MAPL_MetaComp),            intent(INOUT) :: STATE
    type (ESMF_Alarm),     optional, intent(IN   ) :: RUNALARM
    type (MAPL_SunOrbit),  optional, intent(IN   ) :: ORBIT
    integer,               optional, intent(IN   ) :: LM
    logical,               optional, intent(IN   ) :: CHILDINIT
    type (MAPL_LocStream), optional, intent(IN   ) :: LOCSTREAM
    type (MAPL_LocStream), optional, intent(IN   ) :: EXCHANGEGRID
    type (ESMF_Clock)    , optional, intent(IN   ) :: CLOCK
    type (ESMF_Config)   , optional, intent(IN   ) :: CF
    character(len=*)     , optional, intent(IN   ) :: NAME
    integer,               optional, intent(  OUT) :: RC

    character(len=ESMF_MAXSTR), parameter :: IAm = "MAPL_GenericStateSet"
    integer :: STATUS


     if(present(LM)) then
      STATE%GRID%LM=LM
     endif

     if(present(ORBIT)) then
      STATE%ORBIT=ORBIT
     endif

     if(present(RUNALARM)) then
      STATE%ALARM(0)=RUNALARM
     endif

     if(present(CHILDINIT)) then
      STATE%CHILDINIT=CHILDINIT
     endif

     if(present(LOCSTREAM)) then
      STATE%LOCSTREAM=LOCSTREAM
     endif

     if(present(EXCHANGEGRID)) then
      STATE%EXCHANGEGRID=EXCHANGEGRID
     endif

     if(present(CLOCK)) then
      STATE%CLOCK=CLOCK
     endif

     if(present(NAME)) then
        STATE%COMPNAME=NAME
     endif

     if(present(Cf)) then
        STATE%CF=CF
     endif

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GenericStateSet

  subroutine MAPL_GenericStateSetFromGC (GC, ORBIT, LM, RUNALARM, CHILDINIT, &
                                   LOCSTREAM, EXCHANGEGRID, CLOCK, RC)
    type (ESMF_GridComp),            intent(IN   ) :: GC
    type (ESMF_Alarm),     optional, intent(IN   ) :: RUNALARM
    type (MAPL_SunOrbit),  optional, intent(IN   ) :: ORBIT
    integer,               optional, intent(IN   ) :: LM
    logical,               optional, intent(IN   ) :: CHILDINIT
    type (MAPL_LocStream), optional, intent(IN   ) :: LOCSTREAM
    type (MAPL_LocStream), optional, intent(IN   ) :: EXCHANGEGRID
    type (ESMF_Clock)    , optional, intent(IN   ) :: CLOCK
    integer,               optional, intent(  OUT) :: RC

    character(len=ESMF_MAXSTR), parameter :: IAm = "MAPL_GenericStateSetFromGC"
    integer :: STATUS

    type (MAPL_MetaComp), pointer         :: STATE

    call MAPL_InternalStateGet ( GC, STATE, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GenericStateSet (STATE, ORBIT, LM, RUNALARM, CHILDINIT, &
                               LOCSTREAM, EXCHANGEGRID, CLOCK, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GenericStateSetFromGC

  subroutine MAPL_GenericStateClockOn(STATE,NAME,RC)
    type (MAPL_MetaComp),        intent(INOUT) :: STATE
    character(len=*),            intent(IN   ) :: NAME
    integer, optional,           intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm = "MAPL_GenericStateClockOn"
    integer :: STATUS

    call MAPL_ProfClockOn(STATE%TIMES,NAME,RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GenericStateClockOn



  subroutine MAPL_StateAlarmAdd(STATE,ALARM,RC)
    type (MAPL_MetaComp),        intent(INOUT) :: STATE
    type (ESMF_Alarm),           intent(IN   ) :: ALARM
    integer, optional,           intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm = "MAPL_StateAlarmAdd"
    integer :: STATUS

    STATE%ALARMLAST = STATE%ALARMLAST + 1
    ASSERT_(STATE%ALARMLAST <= LAST_ALARM)

    STATE%ALARM(STATE%ALARMLAST) = ALARM

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_StateAlarmAdd

  subroutine MAPL_StateAlarmGet(STATE,ALARM,NAME,RC)
    type (MAPL_MetaComp),        intent(INOUT) :: STATE
    type (ESMF_Alarm),           intent(  OUT) :: ALARM
    character(len=*),            intent(IN   ) :: NAME
    integer, optional,           intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm = "MAPL_StateAlarmGet"
    integer :: STATUS, I
    character(len=ESMF_MAXSTR) :: ANAME

    do I=0,STATE%ALARMLAST
       call ESMF_AlarmGet(STATE%ALARM(I),ANAME,RC=STATUS)
       VERIFY_(STATUS)
       if(trim(NAME)/=trim(ANAME)) cycle
       ALARM=STATE%ALARM(I)
       RETURN_(ESMF_SUCCESS)
    end do

    RETURN_(ESMF_FAILURE)
  end subroutine MAPL_StateAlarmGet






!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================


  subroutine MAPL_GenericStateClockOff(STATE,NAME,RC)
    type (MAPL_MetaComp),        intent(INOUT) :: STATE
    character(len=*),            intent(IN   ) :: NAME
    integer, optional,           intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm = "MAPL_GenericStateClockOff"
    integer :: STATUS

    call MAPL_ProfClockOff(STATE%TIMES,NAME,RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GenericStateClockOff


!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

  subroutine MAPL_GenericStateClockAdd(GC, NAME, RC)

    type (ESMF_GridComp),        intent(INOUT) :: GC
    character(len=*),            intent(IN   ) :: NAME
    integer, optional,           intent(  OUT) :: RC     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm = "MAPL_GenericStateClockAdd"
    integer :: STATUS
    type (MAPL_MetaComp), pointer         :: STATE

    call MAPL_InternalStateRetrieve(GC, STATE, RC=status)
    VERIFY_(STATUS)

    call MAPL_ProfSet(STATE%TIMES,NAME=NAME,RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GenericStateClockAdd


!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================


  subroutine MAPL_GenericGridDestroy (GRID,RC)
    type (MAPL_GenericGrid),  intent(INOUT) :: GRID
    integer, optional,        intent(  OUT) :: rc     ! Error code:

    character(len=ESMF_MAXSTR), parameter :: IAm = "MAPL_GenericGridDestroy"
    integer :: STATUS

    if (associated(GRID%LATS)) deallocate(GRID%LATS)
    if (associated(GRID%LONS)) deallocate(GRID%LONS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GenericGridDestroy


!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================


  subroutine MAPL_ESMFStateWriteToFile(STATE,CLOCK,FILENAME,FILETYPE,RC)
    type(ESMF_State),                 intent(INOUT) :: STATE
    type(ESMF_Clock),                 intent(IN   ) :: CLOCK
    character(len=*),                 intent(IN   ) :: FILENAME
    character(LEN=*),                 intent(IN   ) :: FILETYPE
    integer, optional,                intent(  OUT) :: RC

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_ESMFStateWriteToFile"
    integer                               :: STATUS

    type (ESMF_StateItemType), pointer  :: ITEMTYPES(:)
    character(len=ESMF_MAXSTR ), pointer  :: ITEMNAMES(:)
    integer                               :: ITEMCOUNT
    integer                               :: UNIT
    integer                               :: I, J, N
!    integer                               :: IM_WORLD
!    integer                               :: JM_WORLD
!    integer                               :: MONTH, DAY, HOUR, MINUTE
!    integer                               :: YEAR, SECOND
!    type (ESMF_Time)                      :: CURRENTTIME
    type (ESMF_DELayout)                  :: LAYOUT

! Get information from state
!---------------------------

    call ESMF_StateGet(STATE,ITEMCOUNT=ITEMCOUNT,RC=STATUS)
    VERIFY_(STATUS)

    ASSERT_(ITEMCOUNT>0)

    allocate(ITEMNAMES(ITEMCOUNT),STAT=STATUS)
    VERIFY_(STATUS)
    allocate(ITEMTYPES(ITEMCOUNT),STAT=STATUS)
    VERIFY_(STATUS)

    call ESMF_StateGet(STATE,ITEMNAMELIST=ITEMNAMES,STATEITEMTYPELIST=ITEMTYPES,RC=STATUS)
    VERIFY_(STATUS)

! Open file
!----------
!I/O

    if (filetype == 'binary' .or. filetype == 'BINARY') then
       UNIT = GETFILE(FILENAME, form="unformatted", rc=status)
    elseif(filetype=="formatted".or.filetype=="FORMATTED") then
       UNIT = GETFILE(FILENAME, form="formatted", rc=status)
    else
       UNIT=0
    end if
    VERIFY_(STATUS)

! Write data
!-----------

    if(UNIT/=0) then
       do J = 1, ITEMCOUNT
          call MAPL_VarWrite(UNIT=UNIT, STATE=STATE, NAME=ITEMNAMES(J), rc=status)
          VERIFY_(STATUS)
       end do
       call FREE_FILE(UNIT)
    else
!!!    call MAPL_CFIOStateWrite(file=FILENAME, state=STATE, clock=CLOCK, rc=STATUS)
       STATUS = -1  ! not yet
       VERIFY_(STATUS)
    endif


    deallocate(ITEMNAMES) 
    deallocate(ITEMTYPES)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_ESMFStateWriteToFile

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================




!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================


  subroutine MAPL_ESMFStateReadFromFile(STATE,CLOCK,FILENAME,FILETYPE,RC)
    type(ESMF_State),                 intent(INOUT) :: STATE
    type(ESMF_Clock),                 intent(IN   ) :: CLOCK
    character(LEN=*),                 intent(IN   ) :: FILENAME
    character(LEN=*),                 intent(IN   ) :: FILETYPE
    integer, optional,                intent(  OUT) :: RC

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_ESMFStateReadFromFile"
    integer                               :: STATUS
    integer                               :: UNIT

! Open file
!----------

!I/O

    
    if (filetype == 'binary' .or. filetype == 'BINARY') then
       UNIT = GETFILE(FILENAME, form="unformatted", rc=status)
    elseif(filetype=="formatted".or.filetype=="FORMATTED") then
       UNIT = GETFILE(FILENAME, form="formatted", rc=status)
    else
       UNIT=0
    end if
    VERIFY_(STATUS)


! Read data
! ---------

    if(UNIT/=0) then
       call MAPL_VarRead(UNIT=UNIT, STATE=STATE, RC=STATUS)
       VERIFY_(STATUS)
    else
!!!       call MAPL_CFIOStateRead(FILE=FILENAME, CLOCK=CLOCK, STATE=STATE, RC=STATUS)
       STATUS = -1 ! not yet
       VERIFY_(STATUS)
    endif


    call ESMF_StateSetAttribute(STATE,'MAPL_Initialized',ESMF_TRUE,RC=STATUS)
    VERIFY_(STATUS)

!I/O
    call FREE_FILE(UNIT)


    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_ESMFStateReadFromFile

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================



  subroutine MAPL_StateCreateFromSpec(STATE,SPEC,GRID,TILEGRID,DEFER,RC)
    type(ESMF_State),                 intent(INOUT) :: STATE
    type(MAPL_VarSpec),               intent(INOUT) :: SPEC(:)
    type(ESMF_Grid),                  intent(IN   ) :: GRID
    logical, optional,                intent(IN   ) :: DEFER
    type(ESMF_Grid), optional,        intent(IN   ) :: TILEGRID
    integer, optional,                intent(  OUT) :: RC

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_StateCreateFromVarSpec"
    integer                               :: STATUS

    integer               :: COUNTS(3)
    integer               :: TILECOUNTS(2)
    integer               :: L
    type (ESMF_Array)     :: ARRAY
    type (ESMF_Field)     :: FIELD
    type (ESMF_Bundle)    :: BUNDLE
    type (ESMF_Field)     :: SPEC_FIELD
    type (ESMF_FieldDataMap) :: DATAMAP
    real(kind=ESMF_KIND_R4), pointer         :: VAR_1D(:), VAR_2D(:,:), VAR_3D(:,:,:)
    real(kind=ESMF_KIND_R8), pointer         :: VR8_1D(:), VR8_2D(:,:), VR8_3D(:,:,:)
    logical               :: usableDEFER
    integer               :: RANK
    integer               :: DIMS
    integer               :: STAT
    integer               :: KND
    integer               :: LOCATION
    character(ESMF_MAXSTR):: SHORT_NAME
    character(ESMF_MAXSTR):: LONG_NAME
    character(ESMF_MAXSTR):: UNITS
    character(ESMF_MAXSTR):: FRIENDLYTO
    integer               :: REFRESH
    integer               :: AVGINT
    real                  :: DEFAULT_VALUE
    type(ESMF_Grid)       :: GRIDPTR
    integer               :: I
    logical               ::  done
    integer               :: N, N1, N2, NE
    integer               :: HW
    integer               :: maplist(ESMF_MAXDIM)          ! mapping between array and grid
    character(len=ESMF_MAXSTR), pointer     :: ATTR_INAMES(:)
    character(len=ESMF_MAXSTR), pointer     :: ATTR_RNAMES(:)
    integer,                    pointer     :: ATTR_IVALUES(:)
    real,                       pointer     :: ATTR_RVALUES(:)

    call ESMF_GridGetDELocalInfo(GRID, &
         horzRelLoc=ESMF_CELL_CENTER, &
         vertRelLoc=ESMF_CELL_CELL, &
         localCellCountPerDim=COUNTS,RC=STATUS)
    VERIFY_(STATUS)

   if (present(TILEGRID)) then
      call ESMF_GridGetDELocalInfo(TILEGRID, horzRelLoc=ESMF_CELL_CENTER, localCellCountPerDim=TILECOUNTS,RC=STATUS)
      VERIFY_(STATUS)
   end if

   if (present(DEFER)) then
      usableDEFER = DEFER
   else
      usableDEFER = .false.
   end if

   do L=1,size(SPEC)

      call MAPL_VarSpecGet(SPEC(L),DIMS=DIMS,VLOCATION=LOCATION,   &
                           SHORT_NAME=SHORT_NAME, LONG_NAME=LONG_NAME, UNITS=UNITS,&
                           FIELD=SPEC_FIELD, &
                           STAT=STAT, DEFAULT = DEFAULT_VALUE, &
                           FRIENDLYTO=FRIENDLYTO, &
                           COUPLE_INTERVAL=REFRESH, &
                           ACCMLT_INTERVAL=AVGINT, &
                           HALOWIDTH=HW, &
                           PRECISION=KND, &
                           ATTR_RNAMES=ATTR_RNAMES, &
                           ATTR_INAMES=ATTR_INAMES, &
                           ATTR_RVALUES=ATTR_RVALUES, &
                           ATTR_IVALUES=ATTR_IVALUES, &
                           RC=STATUS )
      VERIFY_(STATUS)

      if(STAT /=  MAPL_CplNEEDED) then
!       cycle
      endif

      I=MAPL_VarSpecGetIndex(SPEC, SHORT_NAME, RC=STATUS)
      if (I /= L) then
         CALL WRITE_PARALLEL("===================>")
         CALL WRITE_PARALLEL(trim(Iam) //": var "// trim(SHORT_NAME) // " already exists. Skipping ...")
         cycle
      endif

      if (IAND(STAT, MAPL_BundleItem) /= 0) then

! Create an empty BUNDLE
! ----------------------
         bundle = ESMF_BundleCreate(NAME=SHORT_NAME, grid=grid, RC=STATUS)
         VERIFY_(STATUS)

! Put the BUNDLE in the state
! --------------------------
      
         call ESMF_StateAddBundle(STATE, BUNDLE, RC=STATUS )
         VERIFY_(STATUS)

         GOTO 10
!!!         cycle
      endif

      call ESMF_FieldGetArray(SPEC_FIELD, array, rc=status)
      if (status==ESMF_SUCCESS) then
         call MAPL_AllocateCoupling( SPEC_FIELD, RC=STATUS ) ! if 'DEFER' this allocates the data
         VERIFY_(STATUS)
         

!ALT we are creating new field so that we can optionally change the name of the field;
!    the important thing is that the data (ESMF_Array) is the SAME as the one in SPEC_Field

         FIELD = MAPL_FieldCreate(SPEC_FIELD, name=SHORT_NAME, RC=STATUS )
         VERIFY_(STATUS)

      else
!         call ESMF_FieldDestroy(SPEC_FIELD, RC=STATUS)
!         VERIFY_(STATUS)

! Create the appropriate ESMF FIELD
! ---------------------------------

         Dimensionality: select case(DIMS)

! Horizontal and vertical
! -----------------------

         case(MAPL_DimsHorzVert)
            GRIDPTR = GRID
            rank = 3
         
            select case(LOCATION)
            
            case(MAPL_VLocationCenter)
               if (knd == ESMF_KIND_R4) then
                  allocate(VAR_3D(1-HW:COUNTS(1)+HW, 1-HW:COUNTS(2)+HW, COUNTS(3)+2*HW), STAT=status)
                  VERIFY_(STATUS)
                  VAR_3D = DEFAULT_VALUE
                  ARRAY = ESMF_ArrayCreate(VAR_3D, ESMF_DATA_REF, HALOWIDTH=HW, RC=STATUS)
               else
                  allocate(VR8_3D(1-HW:COUNTS(1)+HW, 1-HW:COUNTS(2)+HW, COUNTS(3)+2*HW), STAT=status)
                  VERIFY_(STATUS)
                  VR8_3D = DEFAULT_VALUE
                  ARRAY = ESMF_ArrayCreate(VR8_3D, ESMF_DATA_REF, HALOWIDTH=HW, RC=STATUS)
               endif
               VERIFY_(STATUS)
               if (usableDEFER) then
                  call MAPL_ArrayF90Deallocate(array, rc=status)
                  VERIFY_(STATUS)
               endif
            case(MAPL_VLocationEdge  )
               if (knd == ESMF_KIND_R4) then
                  allocate(VAR_3D(1-HW:COUNTS(1)+HW, 1-HW:COUNTS(2)+HW, 0:COUNTS(3)), STAT=status)
                  VERIFY_(STATUS)
                  VAR_3D = DEFAULT_VALUE
                  ARRAY = ESMF_ArrayCreate(VAR_3D, ESMF_DATA_REF, RC=STATUS)
               else
                  allocate(VR8_3D(1-HW:COUNTS(1)+HW, 1-HW:COUNTS(2)+HW, 0:COUNTS(3)), STAT=status)
                  VERIFY_(STATUS)
                  VR8_3D = DEFAULT_VALUE
                  ARRAY = ESMF_ArrayCreate(VR8_3D, ESMF_DATA_REF, RC=STATUS)
               endif
               VERIFY_(STATUS)
               if (usableDEFER) then
                  call MAPL_ArrayF90Deallocate(array, rc=status)
                  VERIFY_(STATUS)
               endif

            case default
               RETURN_(ESMF_FAILURE)
            end select

! Horizontal only
! ---------------

         case(MAPL_DimsHorzOnly)
            GRIDPTR = GRID
            rank=2
            if (knd == ESMF_KIND_R4) then
               allocate(VAR_2D(1-HW:COUNTS(1)+HW, 1-HW:COUNTS(2)+HW), STAT=STATUS)
               VERIFY_(STATUS)
               VAR_2D = DEFAULT_VALUE
               ARRAY = ESMF_ArrayCreate(VAR_2D, ESMF_DATA_REF, HALOWIDTH=HW, RC=STATUS)
            else
               allocate(VR8_2D(1-HW:COUNTS(1)+HW, 1-HW:COUNTS(2)+HW), STAT=STATUS)
               VERIFY_(STATUS)
               VR8_2D = DEFAULT_VALUE
               ARRAY = ESMF_ArrayCreate(VR8_2D, ESMF_DATA_REF, HALOWIDTH=HW, RC=STATUS)
            end if
            VERIFY_(STATUS)
            if (usableDEFER) then
               call MAPL_ArrayF90Deallocate(array, rc=status)
               VERIFY_(STATUS)
            endif

! Vertical only
! -------------

         case(MAPL_DimsVertOnly)
            GRIDPTR = GRID
            rank=1
            
            select case(LOCATION)
               
            case(MAPL_VLocationCenter)
               if (knd == ESMF_KIND_R4) then
                  allocate(VAR_1D(COUNTS(3)), STAT=STATUS)
                  VERIFY_(STATUS)
                  VAR_1D = DEFAULT_VALUE
                  ARRAY = ESMF_ArrayCreate(VAR_1D, ESMF_DATA_REF, RC=STATUS)
               else
                  allocate(VR8_1D(COUNTS(3)), STAT=STATUS)
                  VERIFY_(STATUS)
                  VR8_1D = DEFAULT_VALUE
                  ARRAY = ESMF_ArrayCreate(VR8_1D, ESMF_DATA_REF, RC=STATUS)
               end if
               VERIFY_(STATUS)
               if (usableDEFER) then
                  call MAPL_ArrayF90Deallocate(array, rc=status)
                  VERIFY_(STATUS)
               endif

            case(MAPL_VLocationEdge  )
               if (knd == ESMF_KIND_R4) then
                  allocate(VAR_1D(0:COUNTS(3)), STAT=STATUS)
                  VERIFY_(STATUS)
                  VAR_1D = DEFAULT_VALUE
                  ARRAY = ESMF_ArrayCreate(VAR_1D, ESMF_DATA_REF, RC=STATUS)
               else
                  allocate(VR8_1D(0:COUNTS(3)), STAT=STATUS)
                  VERIFY_(STATUS)
                  VR8_1D = DEFAULT_VALUE
                  ARRAY = ESMF_ArrayCreate(VR8_1D, ESMF_DATA_REF, RC=STATUS)
               endif
               VERIFY_(STATUS)
               if (usableDEFER) then
                  call MAPL_ArrayF90Deallocate(array, rc=status)
                  VERIFY_(STATUS)
               endif
            end select

         case(MAPL_DimsTileOnly)
            GRIDPTR = TILEGRID
            rank=1

            if (knd == ESMF_KIND_R4) then
               allocate(VAR_1D(TILECOUNTS(1)), STAT=STATUS)
               VERIFY_(STATUS)
               VAR_1D = DEFAULT_VALUE
               ARRAY = ESMF_ArrayCreate(VAR_1D, ESMF_DATA_REF, RC=STATUS)
            else
               allocate(VR8_1D(TILECOUNTS(1)), STAT=STATUS)
               VERIFY_(STATUS)
               VR8_1D = DEFAULT_VALUE
               ARRAY = ESMF_ArrayCreate(VR8_1D, ESMF_DATA_REF, RC=STATUS)
            endif
            VERIFY_(STATUS)
            if (usableDEFER) then
               call MAPL_ArrayF90Deallocate(array, rc=status)
               VERIFY_(STATUS)
            endif

         case(MAPL_DimsTileTile)
            GRIDPTR = TILEGRID
            rank=2

            if (knd == ESMF_KIND_R4) then
               allocate(VAR_2D(TILECOUNTS(1), TILECOUNTS(2)), STAT=STATUS)
               VERIFY_(STATUS)
               VAR_2D = DEFAULT_VALUE
               ARRAY = ESMF_ArrayCreate(VAR_2D, ESMF_DATA_REF, RC=STATUS)
            else
               allocate(VR8_2D(TILECOUNTS(1), TILECOUNTS(2)), STAT=STATUS)
               VERIFY_(STATUS)
               VR8_2D = DEFAULT_VALUE
               ARRAY = ESMF_ArrayCreate(VR8_2D, ESMF_DATA_REF, RC=STATUS)
            endif
            VERIFY_(STATUS)
            if (usableDEFER) then
               call MAPL_ArrayF90Deallocate(array, rc=status)
               VERIFY_(STATUS)
            endif
            
! Invalid dimensionality
! ----------------------

         case default 
            RETURN_(ESMF_FAILURE)

         end select Dimensionality
         VERIFY_(STATUS)

! Create an ESMF FIELD with this array
! ------------------------------------

         call ESMF_FieldDataMapSetDefault(datamap, rank, rc=status)
         VERIFY_(STATUS)
         
         if (DIMS == MAPL_DimsVertOnly) then
            call ESMF_FieldDataMapGet(datamap, dataIndexList=maplist, rc=status)
            VERIFY_(STATUS)
            maplist(1) = 3 ! point to Z-axis
            call ESMF_FieldDataMapSet(datamap, dataIndexList=maplist, rc=status)
            VERIFY_(STATUS)
         end if

!         if ( rank == 3 .and. location == MAPL_VLocationEdge ) then
         if ( rank == 3 ) then

!ALT: the proper call should have been
#ifdef ESMF_2_2_WORKS
            FIELD = ESMF_FieldCreate(GRIDPTR, ARRAY,    &
                 horzRelloc = ESMF_CELL_CENTER,         &
                 vertRelloc = ESMF_CELL_VERTEX,         &
                 datamap=datamap,                       &
                 name   = SHORT_NAME,    RC=STATUS )
#else
            FIELD = ESMF_FieldCreate(GRIDPTR, ARRAY,    &
                 horzRelloc = ESMF_CELL_CENTER,         &
                 vertRelloc = ESMF_CELL_TOPFACE,        &
                 datamap=datamap,                       &
                 halowidth=hw,                          &
                 name   = SHORT_NAME,    RC=STATUS )
#endif
         else
            FIELD = ESMF_FieldCreate(GRIDPTR, ARRAY,    &
                 horzRelloc = ESMF_CELL_CENTER,         &
                 halowidth=hw,                          &
                 datamap=datamap,                       &
                 name   = SHORT_NAME,    RC=STATUS )
         end if
         VERIFY_(STATUS)

! Put the FIELD in the MAPL FIELD (VAR SPEC)
! --------------------------------

         call MAPL_VarSpecSet(SPEC(L),FIELD=FIELD,RC=STATUS)
         VERIFY_(STATUS)

      endif
! and in the FIELD in the state
! --------------------------
      
      call ESMF_StateAddField(STATE, FIELD, RC=STATUS )
      VERIFY_(STATUS)

      if (usableDEFER) then
         call ESMF_StateSetNeeded(STATE, SHORT_NAME,  &
              ESMF_NOTNEEDED, rc=status)
         VERIFY_(STATUS)
      end if

! Add SPECs to the FIELD

      call ESMF_FieldSetAttribute(FIELD, NAME='STAT', VALUE=STAT, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldSetAttribute(FIELD, NAME='DIMS', VALUE=DIMS, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldSetAttribute(FIELD, NAME='VLOCATION', VALUE=LOCATION, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldSetAttribute(FIELD, NAME='LONG_NAME', VALUE=LONG_NAME, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldSetAttribute(FIELD, NAME='UNITS', VALUE=UNITS, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldSetAttribute(FIELD, NAME='REFRESH_INTERVAL', VALUE=REFRESH, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldSetAttribute(FIELD, NAME='AVERAGING_INTERVAL', VALUE=AVGINT, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldSetAttribute(FIELD, NAME='HALOWIDTH', VALUE=HW, RC=STATUS)
      VERIFY_(STATUS)

      if (associated(ATTR_RNAMES)) then
         DO N = 1, size(ATTR_RNAMES) 
            call ESMF_FieldSetAttribute(FIELD, NAME=trim(ATTR_RNAMES(N)), &
                                        VALUE=ATTR_RVALUES(N), RC=STATUS)
            VERIFY_(STATUS)
         END DO
      end if

      if (associated(ATTR_INAMES)) then
         DO N = 1, size(ATTR_INAMES) 
            call ESMF_FieldSetAttribute(FIELD, NAME=trim(ATTR_INAMES(N)), &
                                        VALUE=ATTR_IVALUES(N), RC=STATUS)
            VERIFY_(STATUS)
         END DO
      end if

10    if (FRIENDLYTO /= "") then

! parse the string for ":" word delimiters
         done = .false.
         n1 = 1
         NE = len(FRIENDLYTO)

         DO WHILE(.not. DONE)
            N = INDEX(FRIENDLYTO(N1:NE), ':')
            IF (N == 0) then
               DONE = .TRUE.
               N2 = NE
            ELSE
               N2 = N1 + N - 2
            END IF
            if (N1 <= N2 .and. N2 > 0) then
               if (IAND(STAT, MAPL_BundleItem) /= 0) then
                  call ESMF_BundleSetAttribute(BUNDLE, &
                       NAME='FriendlyTo'//trim(FRIENDLYTO(N1:N2)), &
                       VALUE=ESMF_TRUE, RC=STATUS)
                  VERIFY_(STATUS)
               else
!print *,"DEBUG: setting FieldAttr:FriendlyTo"//trim(FRIENDLYTO(N1:N2))
                  call ESMF_FieldSetAttribute(FIELD, &
                       NAME='FriendlyTo'//trim(FRIENDLYTO(N1:N2)), &
                       VALUE=ESMF_TRUE, RC=STATUS)
                  VERIFY_(STATUS)
               end if
            end if

            N1 = N1 + N
         END DO

      end if
   enddo

   RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_StateCreateFromSpec


!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

!BOPI

! !IROUTINE: MAPL_GridCompGetVarSpec

! !INTERFACE:

subroutine MAPL_GridCompGetVarSpecs(GC,IMPORT,EXPORT,INTERNAL,RC)

! !ARGUMENTS:

    type(ESMF_GridComp),    intent(INOUT)  :: GC
    type(MAPL_VarSpec ), pointer, optional :: IMPORT(:)
    type(MAPL_VarSpec ), pointer, optional :: EXPORT(:)
    type(MAPL_VarSpec ), pointer, optional :: INTERNAL(:)
    integer,             optional, intent(OUT) :: RC

!EOPI

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)            :: IAm='MAPL_GridCompGetVarSpec'
    integer                               :: STATUS
    character(len=ESMF_MAXSTR)            :: COMP_NAME

    type (MAPL_MetaComp ), pointer        :: STATE

! Begin

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Get the private state
! ---------------------

    call MAPL_InternalStateRetrieve ( GC, STATE, RC=STATUS )
    VERIFY_(STATUS)

! Get the specs for the 3 ESMF states
! -----------------------------------

    call MAPL_StateGetVarSpecs(STATE,IMPORT,EXPORT,INTERNAL,RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

end subroutine MAPL_GridCompGetVarSpecs


!BOPI

! !IROUTINE: MAPL_StateGetVarSpec

! !INTERFACE:

subroutine MAPL_StateGetVarSpecs(STATE,IMPORT,EXPORT,INTERNAL,RC)

! !ARGUMENTS:

    type(MAPL_MetaComp),       intent(IN)  :: STATE
    type(MAPL_VarSpec ), pointer, optional :: IMPORT(:)
    type(MAPL_VarSpec ), pointer, optional :: EXPORT(:)
    type(MAPL_VarSpec ), pointer, optional :: INTERNAL(:)
    integer,             optional, intent(OUT) :: RC

!EOPI

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)            :: IAm='MAPL_StateGetVarSpec'
    integer                               :: STATUS

! Begin

! Get the specs for the 3 ESMF states
! -----------------------------------

    if(present(IMPORT)) then
     IMPORT => STATE%IMPORT_SPEC
    endif

    if(present(EXPORT)) then
     EXPORT => STATE%EXPORT_SPEC
    endif

    if(present(INTERNAL)) then
     INTERNAL => STATE%INTERNAL_SPEC
    endif

    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_StateGetVarSpecs



!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

!BOP

! !IROUTINE: MAPL_WireComponent

! !INTERFACE:

recursive subroutine MAPL_WireComponent(GC, RC)

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC
    integer,   optional                :: RC      ! return code

! !DESCRIPTION: This connects the child components, creates the couplers,
!               and adds child info to GC's import and export specs.

!EOP

!=============================================================================
!
! ErrLog Variables


    character(len=ESMF_MAXSTR)            :: IAm
    integer                               :: STATUS
    character(len=ESMF_MAXSTR)            :: COMP_NAME

! Local types
! -----------

    type SpecWrapper
     type (MAPL_VarSpec),              pointer :: SPEC(:)
    end type SpecWrapper

! Local variables
! ---------------

    type (MAPL_MetaComp),              pointer  :: STATE 
    type (ESMF_GridComp),              pointer  :: GCS(:)
    type (ESMF_CplComp),               pointer  :: CCS(:,:)
    type (MAPL_VarSpec),               pointer  :: IMPORT_SPECS(:)
    type (MAPL_VarSpec),               pointer  :: EXPORT_SPECS(:)
    type (MAPL_VarSpec),               pointer  :: IM_SPECS(:)
    type (MAPL_VarSpec),               pointer  :: EX_SPECS(:)
    type (MAPL_VarSpec),               pointer  :: SPECS(:)

  
    type (SpecWrapper), pointer, dimension(:,:) :: SRCS
    type (SpecWrapper), pointer, dimension(:,:) :: DSTS

    character(len=ESMF_MAXSTR)                  :: NAME
    character(len=ESMF_MAXSTR)                  :: SRCNAME
    character(len=ESMF_MAXSTR)                  :: DSTNAME
    character(len=ESMF_MAXSTR)                  :: CHILD_NAME
    character(len=ESMF_MAXSTR)                  :: SHORT_NAME
    character(len=ESMF_MAXSTR)                  :: ENAME
    integer                                     :: I, J, K, N
    integer                                     :: NC
    integer                                     :: NCPL
    integer                                     :: STAT
    logical                                     :: SATISFIED
    logical                                     :: PARENTIMPORT
    type (MAPL_ConnectivityWrap)                :: connwrap
    type (MAPL_Connectivity), pointer           :: conn
    type (MAPL_VarConn), pointer                :: CONNECT(:)
    type (MAPL_VarConn), pointer                :: DONOTCONN(:)

! Begin

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'MAPL_WireComponent'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the child components
! --------------------------------------------

    call MAPL_InternalStateRetrieve ( GC, STATE, RC=STATUS )
    VERIFY_(STATUS)

    GCS          => STATE%GCS
    IMPORT_SPECS => STATE%IMPORT_SPEC
    EXPORT_SPECS => STATE%EXPORT_SPEC

! First look for import couplings that are satisfied internally 
! -------------------------------------------------------------

    if (.not. associated(GCS)) then
       RETURN_(ESMF_SUCCESS)
    end if

    call ESMF_UserCompGetInternalState(gc, 'M0APL_Connectivity', &
                                       connwrap, status)
    if (STATUS == ESMF_FAILURE) then
       NULLIFY(CONNECT)
       NULLIFY(DONOTCONN)
    else
       conn => connwrap%ptr
       CONNECT => CONN%CONNECT
       DONOTCONN => CONN%DONOTCONN
    end if

    NC = size(GCS)

    allocate(SRCS(NC,NC), DSTS(NC,NC), STAT=STATUS)
    VERIFY_(STATUS)
    DO I=1,NC
       DO J=1,NC
          NULLIFY(SRCS(J,I)%SPEC)
          NULLIFY(DSTS(J,I)%SPEC)
       END DO
    END DO

! first check if we need to add Exports from children
    do I=1,NC       !  Cycle thru children's imports

       call MAPL_GridCompGetVarSpecs(GCS(I), EXPORT=EX_SPECS, RC=STATUS)
       VERIFY_(STATUS)

       if (.not. associated(EX_SPECS)) then
          cycle
       end if

       do K=1,size(EX_SPECS)
          call MAPL_VarSpecGet(EX_SPECS(K), SHORT_NAME=SHORT_NAME, RC=STATUS) 
          if (MAPL_VarIsConnected(CONNECT, SHORT_NAME=SHORT_NAME, &
                                  FROM_EXPORT=I, TO_EXPORT=MAPL_Self, &
                                  RC=STATUS)) then

             VERIFY_(STATUS)
             call MAPL_VarSpecAddRefToList(STATE%EXPORT_SPEC, EX_SPECS(K), &
                                           RC=STATUS)
             if (STATUS /= MAPL_DuplicateEntry) then
                VERIFY_(STATUS)
             else
                RETURN_(ESMF_FAILURE)
             endif
          end if
          VERIFY_(STATUS)
       end do

    end do

! try to satisfy imports internally
    do I=1,NC       !  Cycle thru children's imports

! check "do not connect" list for 
       PARENTIMPORT = .true.
       if (MAPL_VarIsListed(DONOTCONN, SHORT_NAME="MAPL_AnyChildImport", &
                            IMPORT=I, RC=STATUS)) then
          VERIFY_(STATUS)
          PARENTIMPORT = .false.
       end if
       
       call MAPL_GridCompGetVarSpecs(GCS(I), IMPORT=IM_SPECS, RC=STATUS)
       VERIFY_(STATUS)

       if (.not. associated(IM_SPECS)) then
          cycle
       end if

       do K=1,size(IM_SPECS)

          call MAPL_VarSpecGet(IM_SPECS(K), SHORT_NAME=SHORT_NAME, &
                               STAT=STAT, RC=STATUS) 
          VERIFY_(STATUS)

! do not connect Friendly bundles
          IF (IAND(STAT, MAPL_BundleItem) /= 0) then
             cycle
          end if

! check "do not connect" list 
          if (MAPL_VarIsListed(DONOTCONN, SHORT_NAME=SHORT_NAME, &
                               IMPORT=I, RC=STATUS)) then
             VERIFY_(STATUS)
             cycle
          end if
          VERIFY_(STATUS)

!  Cycle thru all exports
!------------------------
          call MAPL_VarSpecGet(IM_SPECS(K), STAT=STAT, RC=STATUS)
          VERIFY_(STATUS) 

          SATISFIED = .false.
          do J=1,NC      
             if(I==J) cycle
! then check if this is internally satisfied
             if (MAPL_VarIsConnected(CONNECT, IMPORT_NAME=SHORT_NAME, &
                                     IMPORT=I, EXPORT=J, RC=STATUS)) then
                
                VERIFY_(STATUS) 
                SATISFIED = .true.
                cycle
             end if
          end do

          if (SATISFIED) then
             STAT = ior(STAT,MAPL_CplSATISFIED)
             call MAPL_VarSpecSet(IM_SPECS(K), STAT=STAT, RC=STATUS)
             VERIFY_(STATUS)
          end if


          do J=1,NC      
             call MAPL_GridCompGetVarSpecs(GCS(J), EXPORT=EX_SPECS, RC=STATUS)
             VERIFY_(STATUS)

! Trying to satisfy I's imports from J's exports
! ----------------------------------------------


! then check if this is internally satisfied
             if (MAPL_VarIsConnected(CONNECT, &
                                     IMPORT_NAME=SHORT_NAME, EXPORT_NAME=ENAME, &
                                     IMPORT=I, EXPORT=J, RC=STATUS)) then
! If a match is found, add it to that coupler's src and dst specs
! ?? Mark the import satisfied and the export needed.
! -----------------------------------------------------------------
                VERIFY_(STATUS) 
                N =  MAPL_VarSpecGetIndex(EX_SPECS,ENAME,RC=STATUS) 
                if(N /= -1) then
                   VERIFY_(STATUS)
                else
                   RETURN_(ESMF_FAILURE)
                endif
!ALT: currently the function comparing the SPECS assumes SAME names;
!     so we temporarily change the SHORT_NAME, and restore the name after comparison
                call MAPL_VarSpecSet(EX_SPECS(N), SHORT_NAME=SHORT_NAME, RC=STATUS)
                VERIFY_(STATUS)
                if (EX_SPECS(N) == IM_SPECS(K)) then
                   call MAPL_VarSpecSet(EX_SPECS(N), SHORT_NAME=ENAME, RC=STATUS)
                   VERIFY_(STATUS)
! this a direct connection 
! SPECS are the same, no additional averaging is needed
                   call MAPL_Reconnect(STATE,       &
                        GCS(I), MAPL_Import, K,         &
                        GCS(J), MAPL_Export, N, RC=STATUS)
                   VERIFY_(STATUS)

                else
                   call MAPL_VarSpecSet(EX_SPECS(N), SHORT_NAME=ENAME, RC=STATUS)
                   VERIFY_(STATUS)
! coupler is needed
                   call MAPL_VarSpecAddRefToList(DSTS(J,I)%SPEC,IM_SPECS(K), RC=STATUS)
                   VERIFY_(STATUS)

                   call MAPL_VarSpecAddRefToList(SRCS(J,I)%SPEC,EX_SPECS(N), RC=STATUS)
                   VERIFY_(STATUS) 
                end if

             else
! Imports that are not internally satisfied have their specs put in the GC's
! import spec to be externally satisfied.  Their status is left unaltered. 
! --------------------------------------------------------------------------
                if (.not. SATISFIED .and. PARENTIMPORT) then
                   VERIFY_(STATUS) 
                   call MAPL_VarSpecGet(IM_SPECS(K), STAT=STAT, RC=STATUS)
                   VERIFY_(STATUS) 
                   if (iand(STAT,MAPL_CplSATISFIED) /= 0) then
                      cycle
                   end if
                   call MAPL_VarSpecAddRefToList(STATE%IMPORT_SPEC, IM_SPECS(K), RC=STATUS)
                   if (STATUS /= MAPL_DuplicateEntry) then
                      VERIFY_(STATUS)
                   else
                   
                      N =  MAPL_VarSpecGetIndex(STATE%IMPORT_SPEC, IM_SPECS(K),RC=STATUS) 
                      if(N /= -1) then
                         VERIFY_(STATUS)
                      else
                         RETURN_(ESMF_FAILURE)
                      endif

                      call MAPL_Reconnect(STATE,       &
                           GC, MAPL_Import, N,         &
                           GCS(I), MAPL_Import, K, RC=STATUS)
                      VERIFY_(STATUS)
                   end if
                endif
             end if
          enddo
       enddo
    enddo
 
    allocate(STATE%CCS(NC,NC),STAT=STATUS)
    VERIFY_(STATUS)
    allocate(STATE%CIM(NC,NC),STAT=STATUS)
    VERIFY_(STATUS)
    allocate(STATE%CEX(NC,NC),STAT=STATUS)
    VERIFY_(STATUS)
    allocate(STATE%CCcreated(NC,NC),STAT=STATUS)
    VERIFY_(STATUS)
    STATE%CCcreated = .false.

    CCS          => STATE%CCS

    do I=1,NC
       do J=1,NC

          if(associated(DSTS(J,I)%SPEC)) then
             if(I/=J) then
                call ESMF_GridCompGet( GCS(I), NAME=SRCNAME, RC=STATUS )
                VERIFY_(STATUS)

                call ESMF_GridCompGet( GCS(J), NAME=DSTNAME, RC=STATUS )
                VERIFY_(STATUS)

                CCS(J,I) = ESMF_CplCompCreate (                        &
                     NAME       = trim(SRCNAME)//'_2_'//trim(DSTNAME), & 
!                     LAYOUT     = STATE%GRID%LAYOUT,                   &
                     CONFIG     = STATE%CF,                  RC=STATUS )
                VERIFY_(STATUS)

!                STATE%CIM(J,I) = ESMF_StateCreate (     &
!                     STATENAME = trim(SRCNAME)//'_2_'//trim(DSTNAME) // '_Imports', &
!                     STATETYPE = ESMF_STATEEXPORT,               &
!                     RC=STATUS )
!                VERIFY_(STATUS)

!                STATE%CEX(J,I) = ESMF_StateCreate (     &
!                     STATENAME = trim(SRCNAME)//'_2_'//trim(DSTNAME) // '_Exports', &
!                     STATETYPE = ESMF_STATEEXPORT,               &
!                     RC=STATUS )
!                VERIFY_(STATUS)

                STATE%CCcreated(J,I) = .true.

                call WRITE_PARALLEL("Coupler needed for "//trim(SRCNAME)// ' and ' //&
                                    trim(DSTNAME))
                call ESMF_CplCompSetServices (CCS(J,I), GenericCplSetServices, STATUS )
                VERIFY_(STATUS)

                call MAPL_CplCompSetVarSpecs(CCS(J,I),SRCS(J,I)%SPEC,DSTS(J,I)%SPEC,RC=STATUS)
                VERIFY_(STATUS)

             endif
          endif
       enddo
    enddo


! Add my children's exports specs to mine
! ---------------------------------------

! currently disabled (they will be explicitly added to the EXPORT as nested)
!    call MAPL_StateAddExportSpecFrmAll ( STATE, RC=STATUS  )
!    VERIFY_(STATUS)
     

    deallocate(SRCS, DSTS, STAT=STATUS)
    VERIFY_(STATUS)

! Wire my children
! ---------------------------------------

    do I=1,NC
!!!ALT      call MAPL_WireComponent(GCS(I), RC=STATUS)
      VERIFY_(STATUS)
    end do

    RETURN_(ESMF_SUCCESS)
  
  end subroutine MAPL_WireComponent

  subroutine MAPL_BundleInit(GC,STATE,BUNDLENAME,RC)

! !ARGUMENTS:

    type(ESMF_GridComp),           intent(INOUT)  :: GC
    type(ESMF_State),              intent(IN   )  :: STATE
    character (len=*),             intent(IN   )  :: BUNDLENAME
    integer,             optional, intent(  OUT)  :: RC

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)            :: IAm='MAPL_BundleInit'
    integer                               :: STATUS
    character(len=ESMF_MAXSTR)            :: COMP_NAME

! Local variables
! ---------------

    type (MAPL_MetaComp),              pointer  :: INTERNAL_STATE 

    type (ESMF_Bundle)                          :: BUNDLE
    type (ESMF_Field)                           :: FIELD
    type (MAPL_VarSpec),               pointer  :: INTERNAL_SPEC(:)
    integer                                     :: I
    integer                                     :: STAT

!EOP

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the internal state
! --------------------------------------------

    call MAPL_InternalStateRetrieve(GC, INTERNAL_STATE, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_StateGetBundle(STATE, BUNDLENAME, BUNDLE, RC=STATUS)
    VERIFY_(STATUS)

    INTERNAL_SPEC => INTERNAL_STATE%INTERNAL_SPEC
    if (.not. associated(INTERNAL_SPEC)) then
       RETURN_(ESMF_FAILURE)
    end if
    do I = 1, size(INTERNAL_SPEC)
       call MAPL_VarSpecGet(INTERNAL_SPEC(I), FIELD=FIELD, STAT=STAT, RC=STATUS)
       VERIFY_(STATUS)

       if (ior(STAT, MAPL_FriendlyVariable) /= 0) then
          cycle
       end if

!ALT: alternatevly, we could get the field from the INTERNAL_ESMF_STATE

       call ESMF_BundleAddField(bundle, field, rc=status)
       VERIFY_(STATUS)
    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_BundleInit

!BOPI
! !IROUTINE: MAPL_StatePrintSpec

! !INTERFACE: 

  recursive subroutine MAPL_StatePrintSpec(GC, RC)

! !ARGUMENTS:

    type(ESMF_GridComp),           intent(INOUT)  :: GC
    integer,             optional, intent(  OUT)  :: RC

!EOPI

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)            :: IAm='MAPL_StatePrintSpec'
    integer                               :: STATUS
    character(len=ESMF_MAXSTR)            :: COMP_NAME

! Local variables
! ---------------

    type (MAPL_MetaComp),              pointer  :: INTERNAL_STATE 

    type (MAPL_VarSpec),               pointer  :: IMPORT_SPEC(:)
    type (MAPL_VarSpec),               pointer  :: EXPORT_SPEC(:)
    integer                                     :: I
    integer                                     :: STAT

!EOP

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the internal state
! --------------------------------------------

    call MAPL_InternalStateRetrieve ( GC, INTERNAL_STATE, RC=STATUS )
    VERIFY_(STATUS)

    IMPORT_SPEC => INTERNAL_STATE%IMPORT_SPEC
    EXPORT_SPEC => INTERNAL_STATE%EXPORT_SPEC

    if (associated(IMPORT_SPEC)) then
       call WRITE_PARALLEL("==========================")
       call WRITE_PARALLEL("IMPORT spec for " // trim(COMP_NAME))
       call WRITE_PARALLEL("==========================")
       if (associated(IMPORT_SPEC)) then
          do I = 1, size(IMPORT_SPEC)
             call MAPL_VarSpecPrint(IMPORT_SPEC(I), RC=STATUS)
             VERIFY_(STATUS)
          end do
       end if
    end if


    if (associated(EXPORT_SPEC)) then
       call WRITE_PARALLEL("==========================")
       call WRITE_PARALLEL("EXPORT spec for " // trim(COMP_NAME))
       call WRITE_PARALLEL("==========================")
       if (associated(EXPORT_SPEC)) then
          do I = 1, size(EXPORT_SPEC)
             call MAPL_VarSpecPrint(EXPORT_SPEC(I), RC=STATUS)
             VERIFY_(STATUS)
          end do
       end if
    end if

    if (associated(INTERNAL_STATE%GCS)) then
       do I = 1, size(INTERNAL_STATE%GCS)
          call MAPL_StatePrintSpec(INTERNAL_STATE%GCS(I), RC=STATUS)
          VERIFY_(STATUS)
       end do
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_StatePrintSpec


  recursive subroutine MAPL_GenericSpecEnum(GC, SPECS, RC)
! !ARGUMENTS:

    type(ESMF_GridComp),           intent(IN   )  :: GC
    type (MAPL_VarSpec),              pointer     :: SPECS(:)
    integer,             optional, intent(  OUT)  :: RC

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)            :: IAm
    integer                               :: STATUS
    character(len=ESMF_MAXSTR)            :: COMP_NAME

! Local variables
! ---------------

    type (MAPL_MetaComp),              pointer  :: STATE 
    type (MAPL_VarSpec),               pointer  :: IMPORT_SPEC(:)
    type (MAPL_VarSpec),               pointer  :: EXPORT_SPEC(:)
    integer                                     :: I, K
    integer                                     :: LBL

!EOP

! Get my name and set-up traceback handle
!----------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'MAPL_GenericSpecEnum'

! Retrieve the pointer to the internal state
!--------------------------------------------

    call MAPL_InternalStateGet ( GC, STATE, RC=STATUS )
    VERIFY_(STATUS)

    IMPORT_SPEC => STATE%IMPORT_SPEC
    EXPORT_SPEC => STATE%EXPORT_SPEC

! Add all import specs to the list
!  and label by place on list.
!---------------------------------

    if (associated(IMPORT_SPEC)) then
       do I = 1, size(IMPORT_SPEC)
          call MAPL_VarSpecAddRefToList(SPECS, IMPORT_SPEC(I), &
                                        ALLOW_DUPLICATES=.true., RC=STATUS)
          if (STATUS /= MAPL_DuplicateEntry) then
             VERIFY_(STATUS)
          end if
          LBL = size(SPECS)
          call MAPL_VarSpecSet(SPECS(LBL), LABEL=LBL, RC=STATUS )
          VERIFY_(STATUS)

       end do
    end if

! Add all export specs to the list
!  and label in a funny way.
!---------------------------------

    if (associated(EXPORT_SPEC)) then
       do I = 1, size(EXPORT_SPEC)
          call MAPL_VarSpecAddRefToList(SPECS, EXPORT_SPEC(I), &
                                        ALLOW_DUPLICATES=.true., RC=STATUS)
          if (STATUS /= MAPL_DuplicateEntry) then
             VERIFY_(STATUS)
          end if
          K = size(SPECS)
          call MAPL_VarSpecGet(SPECS(K), LABEL=LBL, RC=STATUS )
          VERIFY_(STATUS)
          if (LBL == 0) then
             call MAPL_VarSpecSet(SPECS(K), LABEL=K, RC=STATUS )
             VERIFY_(STATUS)
          end if
       end do
    end if

! Do the same for the children
!-----------------------------

    if (associated(STATE%GCS)) then
       do I = 1, size(STATE%GCS)
          call MAPL_GenericSpecEnum(STATE%GCS(I), SPECS, RC=STATUS)
          VERIFY_(STATUS)
       end do
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_GenericSpecEnum




  subroutine MAPL_LinkCreate(LINK,              &
       GC_FROM, STATETYPE_FROM, SPECINDEX_FROM, &
       GC_TO,   STATETYPE_TO,   SPECINDEX_TO,   &
       RC )

    type (MAPL_Link),              pointer     :: LINK(:)
    type(ESMF_GridComp),        intent(IN   )  :: GC_FROM
    integer,                    intent(IN   )  :: STATETYPE_FROM
    integer,                    intent(IN   )  :: SPECINDEX_FROM
    type(ESMF_GridComp),        intent(IN   )  :: GC_TO
    integer,                    intent(IN   )  :: STATETYPE_TO
    integer,                    intent(IN   )  :: SPECINDEX_TO
    integer,          optional, intent(  OUT)  :: RC     ! Error code:
    


    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_LinkCreate"
    integer                               :: STATUS

    integer                               :: I

    type (MAPL_Link ), pointer         :: TMP(:)
    type (MAPL_LinkType) :: FROM
    type (MAPL_LinkType) :: TO

      if(.not. associated(LINK)) then
       allocate(LINK(0),stat=STATUS)
       VERIFY_(STATUS)
      else
!ALT: check for duplicates ???
      endif
    
        
      I = size(LINK)

      allocate(TMP(I+1),stat=STATUS)
      VERIFY_(STATUS)
      
      TMP(1:I) = LINK
      deallocate(LINK, stat=STATUS)
      VERIFY_(STATUS)

      allocate(TMP(I+1)%Ptr,stat=STATUS)
      VERIFY_(STATUS)

      FROM = MAPL_LinkType(GC_FROM, STATETYPE_FROM, SPECINDEX_FROM)
      TO   = MAPL_LinkType(GC_TO,   STATETYPE_TO,   SPECINDEX_TO  )

      TMP(I+1)%Ptr = MAPL_LinkForm(FROM, TO)

      LINK => TMP

      RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_LinkCreate


  recursive subroutine MAPL_Reconnect(STATE,              &
       GC_FROM, STATETYPE_FROM, SPECINDEX_FROM, &
       GC_TO,   STATETYPE_TO,   SPECINDEX_TO,   &
       RC )

    type (MAPL_MetaComp)                       :: STATE
    type(ESMF_GridComp),        intent(IN   )  :: GC_FROM
    integer,                    intent(IN   )  :: STATETYPE_FROM
    integer,                    intent(IN   )  :: SPECINDEX_FROM
    type(ESMF_GridComp),        intent(IN   )  :: GC_TO
    integer,                    intent(IN   )  :: STATETYPE_TO
    integer,                    intent(IN   )  :: SPECINDEX_TO
    integer,          optional, intent(  OUT)  :: RC     ! Error code:
    


    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_Reconnect"
    integer                               :: STATUS


! Local variables
! ---------------

    type (MAPL_MetaComp),          pointer  :: PSTATE 


! Retrieve the pointer to the internal state of Root. 
! ----------------------------------------------------

    call MAPL_InternalStateRetrieve ( STATE%RootGC, PSTATE, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_LinkCreate(PSTATE%LINK,              &
       GC_FROM, STATETYPE_FROM, SPECINDEX_FROM, &
       GC_TO,   STATETYPE_TO,   SPECINDEX_TO,   &
       RC=STATUS )
    VERIFY_(STATUS)


    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_Reconnect

  integer function MAPL_LabelGet(LINK, RC)
    type (MAPL_LinkType)                       :: LINK
    integer,          optional, intent(  OUT)  :: RC     ! Error code:
    


    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_GetLabel"
    integer                               :: STATUS


! Local variables
! ---------------

    type (MAPL_MetaComp),              pointer  :: STATE 
    type (MAPL_VarSpec),               pointer  :: SPEC(:)


! Retrieve the pointer to the internal state of Root. 
! ----------------------------------------------------

    call MAPL_InternalStateRetrieve ( LINK%GC, STATE, RC=STATUS )
    VERIFY_(STATUS)

! Local aliases to the state
! ---------------------------------------------------

    if (LINK%StateType == MAPL_Import) then
       SPEC => STATE%IMPORT_SPEC
    else if (LINK%StateType == MAPL_Export) then
       SPEC => STATE%EXPORT_SPEC
    else
       RETURN_(ESMF_FAILURE)
    end if



    call MAPL_VarSpecGet(SPEC(LINK%SpecId), LABEL = MAPL_LabelGet, RC = STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end function MAPL_LabelGet

  subroutine MAPL_FriendlyGet ( GC, NAME, FIELD, REQUESTOR, RC )

! !ARGUMENTS:

    type(ESMF_GridComp),           intent(INOUT)  :: GC
    character(len=*),              intent(IN   )  :: NAME
    character(len=*),    optional, intent(IN   )  :: REQUESTOR !ALT (set to optional TEMPORARY)
    type(ESMF_Field),              intent(  OUT)  :: FIELD
    integer,             optional, intent(  OUT)  :: RC

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)            :: IAm='MAPL_FriendlyGet'
    integer                               :: STATUS
    character(len=ESMF_MAXSTR)            :: COMP_NAME

! Local variables
! ---------------

    type (MAPL_MetaComp),              pointer  :: STATE 
    type(ESMF_logical)                          :: FRIENDLY

    integer                                     :: N, STAT


! Retrieve the pointer to the internal state of Root. 
! ----------------------------------------------------

    call MAPL_InternalStateRetrieve ( GC, STATE, RC=STATUS )
    VERIFY_(STATUS)

    N =  MAPL_VarSpecGetIndex(STATE%INTERNAL_SPEC, NAME, RC=STATUS) 
    if(N /= -1) then
       VERIFY_(STATUS)
    else
       RETURN_(ESMF_FAILURE)
    endif

    call MAPL_VarSpecGet(STATE%INTERNAL_SPEC(N), STAT=STAT, RC=STATUS)
    VERIFY_(STATUS)

    ASSERT_(iand(STAT, MAPL_FriendlyVariable) /= 0)

    call ESMF_StateGetField(STATE%INTERNAL, NAME, FIELD, RC=STATUS)
    VERIFY_(STATUS)

    if (present(REQUESTOR)) then
       call ESMF_FieldGetAttribute  (FIELD, NAME="FriendlyTo"//trim(REQUESTOR),VALUE=FRIENDLY, RC=STATUS)
       VERIFY_(STATUS)
       ASSERT_(FRIENDLY==ESMF_TRUE)
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_FriendlyGet





 subroutine MAPL_CopyFriendlinessInState(STATEOUT, NAMEOUT,STATEIN,NAMEIN,RC)
    type(ESMF_STATE),      intent(INOUT) :: STATEOUT
    character(len=*),      intent(IN   ) :: nameOUT
    type(ESMF_STATE),      intent(IN   ) :: STATEIN
    character(len=*),      intent(IN   ) :: nameIN
    integer,    optional,  intent(  OUT) :: rc
    
    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_CopyFriendlinessInState"
    integer                               :: status 
    type(ESMF_FIELD)                      :: FIELDIN
    type(ESMF_FIELD)                      :: FIELDOUT

    call ESMF_StateGetField(STATEIN ,NAMEIN ,FIELDIN ,RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_StateGetField(STATEOUT,NAMEOUT,FIELDOUT,RC=STATUS)
    VERIFY_(STATUS)
    call  MAPL_CopyFriendlinessInField(FIELDOUT,FIELDIN,RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_CopyFriendlinessInState



 subroutine MAPL_CopyFriendlinessInField(FIELDOUT,FIELDIN,RC)
    type(ESMF_FIELD),      intent(INOUT) :: FIELDOUT
    type(ESMF_FIELD),      intent(IN   ) :: FIELDIN
    integer,    optional,  intent(  OUT) :: rc
    
    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_CopyFriendlinessInField"
    integer                               :: status 
    integer                               :: I, NF
    character(len=ESMF_MAXSTR)            :: NAME
    type(ESMF_LOGICAL)                    :: VALUE

    call ESMF_FieldGetAttributeCount(FIELDIN,NF      ,RC=STATUS)
    VERIFY_(STATUS)

    do I=1,NF
       call ESMF_FieldGetAttributeInfo(FIELDIN,I,NAME=NAME,RC=STATUS)
       VERIFY_(STATUS)
       NAME = trim(NAME)
       if(NAME(1:10)=='FriendlyTo') then
          call ESMF_FieldGetAttribute(FIELDIN , NAME=NAME, VALUE=VALUE, RC=STATUS)
          VERIFY_(STATUS)
          call ESMF_FieldSetAttribute(FIELDOUT, NAME=NAME, VALUE=VALUE, RC=STATUS)
          VERIFY_(STATUS)
       end if
    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_CopyFriendlinessInField


  subroutine MAPL_GridCompGetFriendlies0 ( GC, TO, BUNDLE, RC )

! !ARGUMENTS:

    type(ESMF_GridComp),           intent(INOUT)  :: GC
    character(len=*),              intent(IN   )  :: TO
    type(ESMF_Bundle  ),           intent(INOUT)  :: BUNDLE
    integer,             optional, intent(  OUT)  :: RC

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)            :: IAm='MAPL_GridCompGetFriendlies0'
    integer                               :: STATUS

! Local variables
! ---------------

    type (MAPL_MetaComp),        pointer  :: STATE 
    type (ESMF_State)                     :: INTERNAL
    type (ESMF_Field)                     :: FIELD
    type(ESMF_logical)                    :: FRIENDLY
    character (len=ESMF_MAXSTR), pointer  :: itemNameList(:)
    type(ESMF_StateItemType),    pointer  :: stateitemtypeList(:)
    type(ESMF_Bundle)                     :: B
 
    integer                               :: I, N
    integer                               :: J, NF

! Get my MAPL_Generic state
!--------------------------

    call MAPL_InternalStateGet ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    INTERNAL = STATE%INTERNAL

    call ESMF_StateGet(INTERNAL, ITEMCOUNT=N,  RC=STATUS)
    VERIFY_(STATUS)

    allocate(itemNameList(N)     ,STAT=STATUS)
    VERIFY_(STATUS)
    allocate(stateitemtypeList(N),STAT=STATUS)
    VERIFY_(STATUS)

    call ESMF_StateGet(INTERNAL,ITEMNAMELIST=itemNamelist,STATEITEMTYPELIST=stateitemtypeList,RC=STATUS)
    VERIFY_(STATUS)

    do I=1,N
       if(stateitemtypeList(I)==ESMF_STATEITEM_FIELD) then
          call ESMF_StateGetField(INTERNAL,itemNameList(I),FIELD,RC=STATUS)
          VERIFY_(STATUS)
          call ESMF_FieldGetAttribute  (FIELD, NAME="FriendlyTo"//trim(TO), VALUE=FRIENDLY, RC=STATUS)
          if(STATUS==ESMF_SUCCESS) then
             call ESMF_BundleAddField  (BUNDLE, FIELD, RC=STATUS )
             VERIFY_(STATUS)
          end if
       else if(stateitemtypeList(I)==ESMF_STATEITEM_BUNDLE) then
          call ESMF_StateGetBundle(INTERNAL,itemNameList(I), B, RC=STATUS)
          VERIFY_(STATUS)
          call ESMF_BundleGet(B,FieldCount=NF, RC=STATUS)
          VERIFY_(STATUS)
          call ESMF_BundleGetAttribute  (B, NAME="FriendlyTo"//trim(TO), VALUE=FRIENDLY, RC=STATUS)
          if(STATUS==ESMF_SUCCESS) then
! if the bundle is "friendly", copy every single field
             DO J=1,NF
                call ESMF_BundleGetField(B,   J,   FIELD,  RC=STATUS)
                VERIFY_(STATUS)
                call ESMF_BundleAddField  (BUNDLE, FIELD, RC=STATUS )
                VERIFY_(STATUS)
             END DO
          else
! check the fields for "friendliness"
             DO J=1,NF
                call ESMF_BundleGetField(B,   J,   FIELD,  RC=STATUS)
                VERIFY_(STATUS)
                call ESMF_FieldGetAttribute  (FIELD, NAME="FriendlyTo"//trim(TO), VALUE=FRIENDLY, RC=STATUS)
                if(STATUS==ESMF_SUCCESS) then
                   call ESMF_BundleAddField  (BUNDLE, FIELD, RC=STATUS )
                   VERIFY_(STATUS)
                END if
             END DO
          end if
       end if
    end do

    deallocate(itemNameList     ,STAT=STATUS)
    VERIFY_(STATUS)
    deallocate(stateitemtypeList,STAT=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GridCompGetFriendlies0


  subroutine MAPL_GridCompGetFriendlies1 ( GC, TO, BUNDLE, RC )

! !ARGUMENTS:

    type(ESMF_GridComp),           intent(INOUT)  :: GC(:)
    character(len=*),              intent(IN   )  :: TO
    type(ESMF_Bundle  ),           intent(INOUT)  :: BUNDLE
    integer,             optional, intent(  OUT)  :: RC

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)            :: IAm='MAPL_GridCompGetFriendlies1'
    integer                               :: STATUS, I

    do I=1,size(GC)
       call MAPL_GridCompGetFriendlies0(GC(I), TO, BUNDLE, RC=STATUS)
       VERIFY_(STATUS)
    end do

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GridCompGetFriendlies1








!================================
  recursive subroutine MAPL_GenericGridSet(gc,  grid, rc)
! !ARGUMENTS:

    type(ESMF_GridComp),           intent(INOUT)  :: GC
    type(ESMF_Grid),               intent(IN   )  :: GRID
    integer,             optional, intent(  OUT)  :: RC

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)            :: IAm='MAPL_GenericGridSet'
    integer                               :: STATUS
    character(len=ESMF_MAXSTR)            :: COMP_NAME

! Local variables
! ---------------

    type (MAPL_MetaComp),              pointer  :: STATE 
    type (ESMF_Grid)                            :: OWNGRID
    integer                                     :: I

! Retrieve the pointer to the internal state of Root. 
! ----------------------------------------------------

    call MAPL_InternalStateRetrieve ( GC, STATE, RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_GridCompGet(GC, GRID=OWNGRID, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_GridValidate(OWNGRID, RC=STATUS)
    if (STATUS /= ESMF_SUCCESS) then
       call ESMF_GridCompSet(GC, GRID=GRID, RC=STATUS)
       VERIFY_(STATUS)
    endif

    if (associated(STATE%GCS)) then
       do I = 1, size(STATE%GCS)
          call MAPL_GenericGridSet(STATE%GCS(I),  grid=grid, rc=status)
          VERIFY_(STATUS)
       end do
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GenericGridSet

   subroutine MAPL_SetVarSpecForCC(gcA, gcB, ccAxB, rc)
    type(ESMF_GridComp), intent(inout) :: GCA 
    type(ESMF_GridComp), intent(inout) :: GCB 
    type(ESMF_CplComp) , intent(inout) :: CCAxB 
    integer, optional,   intent(  out) :: RC     ! Error code:

! Local vars
    character(len=ESMF_MAXSTR)   :: Iam="MAPL_SetVarSpecForCC"
    character(len=ESMF_MAXSTR)   :: NAME
    integer                      :: STATUS
    integer                      :: I, N, STAT
    type (MAPL_VarSpec), pointer :: SRCS(:)
    type (MAPL_VarSpec), pointer :: DSTS(:)
    type (MAPL_VarSpec), pointer :: IM_SPECS(:), EX_SPECS(:)

! Begin

    NULLIFY(SRCS)
    NULLIFY(DSTS)

    call MAPL_GridCompGetVarSpecs(gcA, EXPORT=EX_SPECS, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompGetVarSpecs(gcB, IMPORT=IM_SPECS,  RC=STATUS)
    VERIFY_(STATUS)

    DO I = 1, size(IM_SPECS)
       call MAPL_VarSpecGet(IM_SPECS(I), STAT=STAT, RC=STATUS)
       VERIFY_(STATUS) 

       IF (IAND(STAT, MAPL_BundleItem) /= 0) then
          cycle
       END IF

       call MAPL_VarSpecAddRefToList(DSTS, IM_SPECS(I), RC=STATUS)
       VERIFY_(STATUS)
    END DO

    IF (.not. associated(DSTS)) then
       RETURN_(ESMF_FAILURE)
    END IF

    DO I = 1, size(DSTS)
       call MAPL_VarSpecGet(DSTS(I), STAT=STAT, SHORT_NAME=NAME, RC=STATUS)
       VERIFY_(STATUS) 

       N =  MAPL_VarSpecGetIndex(EX_SPECS, NAME, RC=STATUS) 
       if(N /= -1) then
          VERIFY_(STATUS)
       else
          call WRITE_PARALLEL("ERROR: cannot match spec:")
          call MAPL_VarSpecPrint(DSTS(I))
          RETURN_(ESMF_FAILURE)
       endif

       call MAPL_VarSpecAddRefToList(SRCS, DSTS(I), RC=STATUS)
       VERIFY_(STATUS)
    END DO

    call MAPL_CplCompSetVarSpecs(ccAxB, SRCS, DSTS, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_SetVarSpecForCC


  recursive subroutine MAPL_GenericConnCheck(GC, RC)
! !ARGUMENTS:

    type(ESMF_GridComp),           intent(INOUT)  :: GC
    integer,             optional, intent(  OUT)  :: RC

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)            :: IAm
    integer                               :: STATUS
    character(len=ESMF_MAXSTR)            :: COMP_NAME

! Local variables
! ---------------

    type (MAPL_MetaComp),              pointer  :: STATE 

    integer                                     :: I
    logical                                     :: err
    type (MAPL_ConnectivityWrap)          :: connwrap
    type (MAPL_Connectivity), pointer     :: conn



!EOP

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'MAPL_GenericConnCheck'

! Retrieve the pointer to the internal state
! --------------------------------------------

    call MAPL_InternalStateRetrieve ( GC, STATE, RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_UserCompGetInternalState(gc, 'M0APL_Connectivity', &
                                       connwrap, status)
    if (STATUS == ESMF_FAILURE) then
       allocate(conn, STAT=STATUS)
       VERIFY_(STATUS)
       connwrap%ptr => conn
       call ESMF_UserCompSetInternalState(gc, 'M0APL_Connectivity', &
                                          connwrap, status)
       VERIFY_(STATUS)
    else
       conn => connwrap%ptr
    end if

    err = .false.
    if (.not. MAPL_ConnCheckUnused(CONN%CONNECT)) then
       err = .true.
       CALL WRITE_PARALLEL("CONNECT ERRORS FOUND in " // trim(COMP_NAME))
    end if

    if (.not. MAPL_ConnCheckUnused(CONN%DONOTCONN)) then
       err = .true.
       CALL WRITE_PARALLEL("DO_NOT_CONNECT ERRORS FOUND in " // trim(COMP_NAME))
    end if

    if (associated(STATE%GCS)) then
       do I = 1, size(STATE%GCS)
          call MAPL_GenericConnCheck(STATE%GCS(I), RC=STATUS)
          if (status /= ESMF_SUCCESS) then
             err = .true.
          end if
       end do
    end if

    if (err) then
       RETURN_(ESMF_FAILURE)
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GenericConnCheck



  subroutine MAPL_GetResourceI41(STATE,VALUE,LABEL,DEFAULT,RC)
    type (MAPL_MetaComp),       intent(INOUT)    :: STATE
    character(len=*),           intent(IN   )    :: LABEL
    integer*4,                  intent(INOUT)    :: VALUE(:)
    integer*4, optional,        intent(IN   )    :: DEFAULT(:)
    integer  , optional,        intent(  OUT)    :: RC
    
    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_ResourceI4"
    integer                               :: STATUS
    integer                               :: I,J
    logical                               :: TEND
    character(len=ESMF_MAXSTR)            :: LBL(3)
    character(len=ESMF_MAXSTR)            :: TYPE 

    TYPE = STATE%COMPNAME(index(STATE%COMPNAME,":")+1:)

    LBL(1) = trim(STATE%COMPNAME)//'_'//trim(LABEL)
    LBL(2) = trim(TYPE)//'_'//trim(LABEL)
    LBL(3) = trim(LABEL)

    DO I=1,SIZE(LBL)
       call ESMF_ConfigFindLabel( STATE%CF, label=trim(LBL(I)), rc=status )
       IF (STATUS == ESMF_SUCCESS) then
          DO J=1,size(VALUE)
             call ESMF_ConfigGetAttribute(STATE%CF, VALUE(J), RC=STATUS )
             VERIFY_(STATUS)
          end DO
          call ESMF_ConfigNextLine  ( STATE%CF, tableEnd=tend, rc=STATUS )
          VERIFY_(STATUS)
          ASSERT_(TEND)
          RETURN_(ESMF_SUCCESS)
       END IF
    ENDDO

    if (present(DEFAULT)) then
       VALUE = DEFAULT
    else
       if (present(RC)) then
          RC = ESMF_FAILURE
          return
       end if
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GetResourceI41


  subroutine MAPL_GetResourceI4(STATE,VALUE,LABEL,DEFAULT,RC)
    type (MAPL_MetaComp),       intent(INOUT)    :: STATE
    character(len=*),           intent(IN   )    :: LABEL
    integer*4,                  intent(INOUT)    :: VALUE
    integer  , optional,        intent(IN   )    :: DEFAULT
    integer  , optional,        intent(  OUT)    :: RC
    
    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_ResourceI4"
    integer                               :: STATUS
    integer                               :: I
    character(len=ESMF_MAXSTR)            :: LBL(3)
    character(len=ESMF_MAXSTR)            :: TYPE 

    TYPE = STATE%COMPNAME(index(STATE%COMPNAME,":")+1:)

    LBL(1) = trim(STATE%COMPNAME)//'_'//trim(LABEL)
    LBL(2) = trim(TYPE)//'_'//trim(LABEL)
    LBL(3) = trim(LABEL)

    DO I=1,SIZE(LBL)
       call ESMF_ConfigFindLabel( STATE%CF, label=trim(LBL(I)), rc=status )
       IF (STATUS == ESMF_SUCCESS) then
          call ESMF_ConfigGetAttribute(STATE%CF, VALUE, &
               label = trim(LBL(I)), &
               default = DEFAULT, RC=STATUS )
          RETURN_(STATUS)
       END IF
    ENDDO


    if (present(DEFAULT)) then
       VALUE = DEFAULT
       RETURN_(ESMF_SUCCESS)
    end if
    RETURN_(ESMF_FAILURE)
  end subroutine MAPL_GetResourceI4


  subroutine MAPL_GetResourceI8(STATE,VALUE,LABEL,DEFAULT,RC)
    type (MAPL_MetaComp),            intent(INOUT)    :: STATE
    character(len=*),                intent(IN   )    :: LABEL
    integer(ESMF_KIND_I8),           intent(  OUT)    :: VALUE
    integer(ESMF_KIND_I8), optional, intent(IN   )    :: DEFAULT
    integer  , optional,             intent(  OUT)    :: RC

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_ResourceI8"
    integer                               :: STATUS
    integer                               :: I
    character(len=ESMF_MAXSTR)            :: LBL(3)
    character(len=ESMF_MAXSTR)            :: TYPE 

    TYPE = STATE%COMPNAME(index(STATE%COMPNAME,":")+1:)

    LBL(1) = trim(STATE%COMPNAME)//'_'//trim(LABEL)
    LBL(2) = trim(TYPE)//'_'//trim(LABEL)
    LBL(3) = trim(LABEL)

    DO I=1,SIZE(LBL)
       call ESMF_ConfigFindLabel( STATE%CF, label=trim(LBL(I)), rc=status )
       IF (STATUS == ESMF_SUCCESS) then
          call ESMF_ConfigGetAttribute(STATE%CF, VALUE, &
               label = trim(LBL(I)), &
               default = DEFAULT, RC=STATUS )
          RETURN_(STATUS)
       END IF
    ENDDO

    if (present(DEFAULT)) then
       VALUE = DEFAULT
       RETURN_(ESMF_SUCCESS)
    end if
    RETURN_(ESMF_FAILURE)
  end subroutine MAPL_GetResourceI8

  subroutine MAPL_GetResourceR4(STATE,VALUE,LABEL,DEFAULT,RC)
    type (MAPL_MetaComp),       intent(INOUT)    :: STATE
    character(len=*),           intent(IN   )    :: LABEL
    real*4,                     intent(INOUT)    :: VALUE
    real     , optional,        intent(IN   )    :: DEFAULT
    integer  , optional,        intent(  OUT)    :: RC
    
    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_ResourceR4"
    integer                               :: STATUS
    integer                               :: I
    character(len=ESMF_MAXSTR)            :: LBL(3)
    character(len=ESMF_MAXSTR)            :: TYPE 
    
    TYPE = STATE%COMPNAME(index(STATE%COMPNAME,":")+1:)
    
    LBL(1) = trim(STATE%COMPNAME)//'_'//trim(LABEL)
    LBL(2) = trim(TYPE)//'_'//trim(LABEL)
    LBL(3) = trim(LABEL)
    
    DO I=1,SIZE(LBL)
       call ESMF_ConfigFindLabel( STATE%CF, label=trim(LBL(I)), rc=status )
       IF (STATUS == ESMF_SUCCESS) then
          call ESMF_ConfigGetAttribute(STATE%CF, VALUE, &
               label = trim(LBL(I)), &
               default = DEFAULT, RC=STATUS )
          RETURN_(STATUS)
       END IF
    ENDDO
    
    if (present(DEFAULT)) then
       VALUE = DEFAULT
       RETURN_(ESMF_SUCCESS)
    end if
    RETURN_(ESMF_FAILURE)
  end subroutine MAPL_GetResourceR4


  subroutine MAPL_GetResourceR8(STATE,VALUE,LABEL,DEFAULT,RC)
    type (MAPL_MetaComp),         intent(INOUT)    :: STATE
    character(len=*),             intent(IN   )    :: LABEL
    real(ESMF_KIND_R8),           intent(  OUT)    :: VALUE
    real(ESMF_KIND_R8), optional, intent(IN   )    :: DEFAULT
    integer  , optional,          intent(  OUT)    :: RC

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_ResourceR8"
    integer                               :: STATUS
    integer                               :: I
    character(len=ESMF_MAXSTR)            :: LBL(3)
    character(len=ESMF_MAXSTR)            :: TYPE 

    TYPE = STATE%COMPNAME(index(STATE%COMPNAME,":")+1:)

    LBL(1) = trim(STATE%COMPNAME)//'_'//trim(LABEL)
    LBL(2) = trim(TYPE)//'_'//trim(LABEL)
    LBL(3) = trim(LABEL)

    DO I=1,SIZE(LBL)
       call ESMF_ConfigFindLabel( STATE%CF, label=trim(LBL(I)), rc=status )
       IF (STATUS == ESMF_SUCCESS) then
          call ESMF_ConfigGetAttribute(STATE%CF, VALUE, &
               label = trim(LBL(I)), &
               default = DEFAULT, RC=STATUS )
          RETURN_(STATUS)
       END IF
    ENDDO

    if (present(DEFAULT)) then
       VALUE = DEFAULT
       RETURN_(ESMF_SUCCESS)
    end if
    RETURN_(ESMF_FAILURE)
 end subroutine MAPL_GetResourceR8

 subroutine MAPL_GetResourceC(STATE,VALUE,LABEL,DEFAULT,RC)
   type (MAPL_MetaComp),       intent(INOUT)    :: STATE
   character(len=*),           intent(IN   )    :: LABEL
   character(len=*),           intent(INOUT)    :: VALUE
   character(len=*), optional, intent(IN   )    :: DEFAULT
   integer  , optional,        intent(  OUT)    :: RC

   character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_ResourceC"
   integer                               :: STATUS
   integer                               :: I
   character(len=ESMF_MAXSTR)            :: LBL(3)
   character(len=ESMF_MAXSTR)            :: TYPE 

   TYPE = STATE%COMPNAME(index(STATE%COMPNAME,":")+1:)
   
   LBL(1) = trim(STATE%COMPNAME)//'_'//trim(LABEL)
   LBL(2) = trim(TYPE)//'_'//trim(LABEL)
   LBL(3) = trim(LABEL)

   DO I=1,SIZE(LBL)
      call ESMF_ConfigFindLabel( STATE%CF, label=trim(LBL(I)), rc=status )
      IF (STATUS == ESMF_SUCCESS) then
         call ESMF_ConfigGetAttribute(STATE%CF, VALUE, &
              label = trim(LBL(I)), &
              default = DEFAULT, RC=STATUS )
         RETURN_(STATUS)
      END IF
   ENDDO

   if (present(DEFAULT)) then
      VALUE = DEFAULT
      RETURN_(ESMF_SUCCESS)
   end if

   if (present(RC)) then
      RC = ESMF_FAILURE
   end if

 end subroutine MAPL_GetResourceC

 integer function MAPL_GetNumSubtiles(STATE, RC)
    type (MAPL_MetaComp),       intent(INOUT)    :: STATE
    integer  , optional,        intent(  OUT)    :: RC

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_GetNumSubtiles"
    integer                               :: STATUS
    integer                               :: I
    integer                               :: DIMS
    integer                               :: NUM_SUBTILES

    MAPL_GetNumSubtiles = 1
    if (associated(STATE%INTERNAL_SPEC)) then
       DO I = 1, size(STATE%INTERNAL_SPEC)
          call MAPL_VarSpecGet(STATE%INTERNAL_SPEC(I), DIMS = DIMS,       &
                 NUM_SUBTILES=NUM_SUBTILES,                 &
                 RC=STATUS  )
          VERIFY_(STATUS)
          if (DIMS == MAPL_DimsTileTile) then
             MAPL_GetNumSubtiles = NUM_SUBTILES
             RETURN_(ESMF_SUCCESS)
          end if
       END DO
    end if
    RETURN_(ESMF_SUCCESS)
  end function MAPL_GetNumSubtiles


  recursive subroutine MAPL_AdjustIsNeeded ( GC, EXPORT, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component 
    type(ESMF_State),    intent(INOUT) :: EXPORT ! Export state
    integer, optional,   intent(  OUT) :: RC     ! Error code:

!EOPI

! ErrLog Variables

    character(len=ESMF_MAXSTR)    :: IAm
    character(len=ESMF_MAXSTR)    :: COMP_NAME
    integer                       :: STATUS

! Local derived type aliases

    type (MAPL_MetaComp),pointer     :: STATE 
    integer                          :: I
    integer                          :: ITEMCOUNT
    type (ESMF_StateItemType), pointer  :: ITEMTYPES(:)
    type (ESMF_Field)                   :: FIELD
    character(len=ESMF_MAXSTR ), pointer  :: ITEMNAMES(:)

    Iam = "MAPL_AdjustIsNeeded"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the internal state.
! -------------------------------------------

    call MAPL_InternalStateGet ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_StateGet(EXPORT,ITEMCOUNT=ITEMCOUNT,RC=STATUS)
    VERIFY_(STATUS)

    IF (ITEMCOUNT>0) then

       allocate(ITEMNAMES(ITEMCOUNT),STAT=STATUS)
       VERIFY_(STATUS)
       allocate(ITEMTYPES(ITEMCOUNT),STAT=STATUS)
       VERIFY_(STATUS)

       call ESMF_StateGet(EXPORT,ITEMNAMELIST=ITEMNAMES,STATEITEMTYPELIST=ITEMTYPES,RC=STATUS)
       VERIFY_(STATUS)

       do I = 1, ITEMCOUNT
          if(itemtypes(I)==ESMF_STATEITEM_FIELD) then
             call ESMF_StateGetField(EXPORT,itemNames(I),FIELD,RC=STATUS)
             VERIFY_(STATUS)
             if (MAPL_IsFieldAllocated(FIELD, RC=STATUS)) then
                VERIFY_(STATUS)
                call ESMF_StateSetNeeded(EXPORT, itemNames(I),  &
                     ESMF_NEEDED, rc=status)
                VERIFY_(STATUS)
             else
                VERIFY_(STATUS)
                call ESMF_StateSetNeeded(EXPORT, itemNames(I),  &
                     ESMF_NOTNEEDED, rc=status)
                VERIFY_(STATUS)
             end if

          end if
       end do

       deallocate(ITEMNAMES)
       deallocate(ITEMTYPES)
    end IF

    if (associated(STATE%GCS)) then
       do I = 1, size(STATE%GCS)
          call MAPL_AdjustIsNeeded(STATE%GCS(I), STATE%GEX(I), RC=STATUS)
          VERIFY_(STATUS)
       end do
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_AdjustIsNeeded

!BOPI

! !INTERFACE:

  logical function MAPL_IsFieldAllocated(FIELD, RC)

! !ARGUMENTS:

    type(ESMF_Field),    intent(IN   ) :: FIELD  ! Field
    integer, optional,   intent(  OUT) :: RC     ! Error code:

! !DESCRIPTION:  Shortcut for checking that field is allocated
 
!EOPI

! ErrLog Variables

    character(len=ESMF_MAXSTR), parameter :: IAm='MAPL_IsFieldAllocated'
    integer                               :: STATUS
    
    real(kind=4), dimension(:)        , pointer :: r4d1
    real(kind=4), dimension(:,:)      , pointer :: r4d2
    real(kind=4), dimension(:,:,:)    , pointer :: r4d3
    real(kind=4), dimension(:,:,:,:)  , pointer :: r4d4

    type(ESMF_Array)                        :: array
    type(ESMF_DataKind)                     :: dk
    type(ESMF_Pointer)                      :: base
    integer                                 :: rank

    call ESMF_FieldGetArray (FIELD, ARRAY, RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_ArrayGet(array, base=base, rc=status)
    VERIFY_(STATUS)

    MAPL_IsFieldAllocated = base /= ESMF_NULL_POINTER

#if 1

    MAPL_IsFieldAllocated = .true.

    call ESMF_ArrayGet(array, rank, kind=dk, rc=status)
    if (dk .eq. ESMF_R4) then
       if (rank .eq. 1) then
          call ESMF_ArrayGetData(array, r4d1, rc=status)
          VERIFY_(STATUS)
          if (.not. associated(r4d1)) then
             MAPL_IsFieldAllocated = .false.
             RETURN_(ESMF_SUCCESS)
          endif
       else if (rank .eq. 2) then
          call ESMF_ArrayGetData(array, r4d2, rc=status)
          VERIFY_(STATUS)
          if (.not. associated(r4d2)) then
             MAPL_IsFieldAllocated = .false.
             RETURN_(ESMF_SUCCESS)
          endif
       else if (rank .eq. 3) then
          call ESMF_ArrayGetData(array, r4d3, rc=status)
          VERIFY_(STATUS)
          if (.not. associated(r4d3)) then
             MAPL_IsFieldAllocated = .false.
             RETURN_(ESMF_SUCCESS)
          endif
       else if (rank .eq. 4) then
          call ESMF_ArrayGetData(array, r4d4, rc=status)
          VERIFY_(STATUS)
          if (.not. associated(r4d4)) then
             MAPL_IsFieldAllocated = .false.
             RETURN_(ESMF_SUCCESS)
          endif
       else
          RETURN_(ESMF_FAILURE)
       end if
    end if

#endif

    RETURN_(ESMF_SUCCESS)
  end function MAPL_IsFieldAllocated

  subroutine MAPL_ExchangeGridGet ( GC, EXCH, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(IN   ) :: GC     ! Gridded component 
    type(MAPL_LocStream),intent(INOUT) :: EXCH   ! Exchange grid
    integer, optional,   intent(  OUT) :: RC     ! Error code:

!EOPI

! ErrLog Variables

    character(len=ESMF_MAXSTR)    :: IAm
    character(len=ESMF_MAXSTR)    :: COMP_NAME
    integer                       :: STATUS

! Local derived type aliases

    type (MAPL_MetaComp),pointer     :: STATE 
    integer                          :: I

    Iam = "MAPL_ExchangeGridGet"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the internal state.
! -------------------------------------------

    call MAPL_InternalStateGet ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    EXCH = STATE%EXCHANGEGRID

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_ExchangeGridGet
    
  recursive subroutine MAPL_ExchangeGridSet ( GC, EXCH, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component 
    type(MAPL_LocStream),intent(IN   ) :: EXCH   ! Exchange grid
    integer, optional,   intent(  OUT) :: RC     ! Error code:

!EOPI

! ErrLog Variables

    character(len=ESMF_MAXSTR)    :: IAm
    character(len=ESMF_MAXSTR)    :: COMP_NAME
    integer                       :: STATUS

! Local derived type aliases

    type (MAPL_MetaComp),pointer     :: STATE 
    integer                          :: I

    Iam = "MAPL_ExchangeGridSet"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the internal state.
! -------------------------------------------

    call MAPL_InternalStateGet ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    STATE%EXCHANGEGRID=EXCH

    if (associated(STATE%GCS)) then
       do I = 1, size(STATE%GCS)
          call MAPL_ExchangeGridSet(STATE%GCS(I), exch, RC=STATUS)
          VERIFY_(STATUS)
       end do
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_ExchangeGridSet
    
  recursive subroutine MAPL_ExportStateGet(export, name, result, rc)
    type (ESMF_State), intent(IN   ) :: export(:)
    character(len=*),  intent(IN   ) :: name
    type (ESMF_State), intent(  OUT) :: result
    integer,           intent(  OUT) :: rc
    
    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_ExportStateGet"
    integer                               :: status 
    
    integer                               :: n, i, ni, k, j
    character(len=ESMF_MAXSTR)            :: sname
    type (ESMF_StateItemType), pointer    :: ITEMTYPES(:)
    character(len=ESMF_MAXSTR ), pointer  :: ITEMNAMES(:)
    type (ESMF_State), pointer            :: exptmp(:)
   
    n = size(export)
    
    do i = 1, n
       call ESMF_StateGet(export(i), name=sname, itemcount = ni, rc=status)
       VERIFY_(STATUS)
       if (sname == trim(name) // '_Exports') then
          result = export(i)
          RETURN_(ESMF_SUCCESS)
       end if
       
       allocate(itemtypes(ni), itemnames(ni), stat=status)
       VERIFY_(STATUS)
       
       call ESMF_StateGet(export(i), ITEMNAMELIST=ITEMNAMES,STATEITEMTYPELIST=ITEMTYPES,RC=STATUS)
       VERIFY_(STATUS)
      
       j = 0
       do k = 1, ni
          if (itemtypes(k) == ESMF_StateItem_State) then
             j = j+1
          end if
       end do
      
       allocate(exptmp(j), stat=status)
       VERIFY_(STATUS)
      
       j = 0
       do k = 1, ni
          if (itemtypes(k) == ESMF_StateItem_State) then
             j = j+1
             call ESMF_StateGetState(export(i), itemnames(k), exptmp(j) , rc=status)
             VERIFY_(STATUS)
          end if
       end do
      
       call MAPL_ExportStateGet(exptmp, name, result, rc=status)
       deallocate(exptmp)
       deallocate(itemtypes, itemnames)
       if (status == ESMF_SUCCESS) return
   end do

   rc = ESMF_FAILURE
   return
 end subroutine MAPL_ExportStateGet

 recursive subroutine MAPL_ImportStateGet ( GC, import, name, result, rc )
    type(ESMF_GridComp), intent(IN)  :: GC
    type(ESMF_State),    intent(IN)  :: import
    character(len=*),  intent(IN   ) :: name
    type (ESMF_State), intent(  OUT) :: result
    integer,           intent(  OUT) :: rc
    
    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_ImportStateGet"
    integer                               :: status 
    integer                               :: I
    character(len=ESMF_MAXSTR)            :: sname
    type (MAPL_MetaComp),pointer          :: STATE 
    

    call ESMF_StateGet(import, name=sname, rc=status)
    VERIFY_(STATUS)
    if (sname == trim(name) // '_Imports') then
       result = import
       RETURN_(ESMF_SUCCESS)
    end if

    call MAPL_InternalStateGet ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    if (associated(STATE%GCS)) then
       do I = 1, size(STATE%GCS)
          call MAPL_ImportStateGet(STATE%GCS(I), STATE%GIM(I), name, result, RC=STATUS)
          if (status == ESMF_SUCCESS) then
             RETURN_(ESMF_SUCCESS)
          end if
       end do
    end if
   rc = ESMF_FAILURE
   return
 end subroutine MAPL_ImportStateGet

 recursive subroutine MAPL_GetChildLocstream(GC, result, name, rc)
    type(ESMF_GridComp),   intent(IN   ) :: GC
    type (MAPL_LocStream), intent(  OUT) :: result
    character(len=*),      intent(IN   ) :: name
    integer,               intent(  OUT) :: rc
    
    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_GetChildLocstream"
    integer                               :: status 
    integer                               :: I
    character(len=ESMF_MAXSTR)            :: comp_name
    type (MAPL_MetaComp),pointer          :: STATE 
    
    call ESMF_GridCompGet(GC, name=comp_name, rc=status)
    VERIFY_(STATUS)

    call MAPL_InternalStateGet ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    if (name == comp_name) then
       result = state%locstream
       RETURN_(STATUS)
    end if


    if (associated(STATE%GCS)) then
       do I = 1, size(STATE%GCS)
          call MAPL_GetChildLocstream(STATE%GCS(I), result, name, rc=STATUS)
          if (status == ESMF_SUCCESS) then
             RETURN_(ESMF_SUCCESS)
          end if
       end do
    end if

    rc = ESMF_FAILURE
    return
  end subroutine MAPL_GetChildLocstream


  function MAPL_VerifyFriendlyInField(FIELD,FRIEND2COMP,RC) result(FRIENDLY)
    type (ESMF_field),          intent(IN)    :: FIELD
    character(len=*),           intent(IN)    :: FRIEND2COMP
    integer  , optional,        intent(OUT)   :: RC
    type(ESMF_logical)                        :: FRIENDLY

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_VerifyFriendlyField"
    integer                               :: STATUS

    call ESMF_FieldGetAttribute  (FIELD, NAME="FriendlyTo"//trim(FRIEND2COMP), &
                                        VALUE=FRIENDLY, RC=STATUS)
    if(STATUS/=ESMF_SUCCESS) FRIENDLY = ESMF_false

    RETURN_(ESMF_SUCCESS)

  end function MAPL_VerifyFriendlyInField

  function MAPL_VerifyFriendlyInState(STATE,NAME,FRIEND2COMP,RC) result(FRIENDLY)
    type (ESMF_State),          intent(IN)    :: STATE
    character(len=*),           intent(IN)    :: NAME
    character(len=*),           intent(IN)    :: FRIEND2COMP
    integer  , optional,        intent(OUT)   :: RC
    type(ESMF_logical)                        :: FRIENDLY

    character(len=ESMF_MAXSTR), parameter :: IAm="MAPL_VerifyFriendlyState"
    integer                               :: STATUS
    type (ESMF_field)                     :: FIELD

    call ESMF_StateGetField(STATE,NAME,FIELD,RC=STATUS)
    VERIFY_(STATUS)
    FRIENDLY=MAPL_VerifyFriendly(FIELD,FRIEND2COMP,RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end function MAPL_VerifyFriendlyInState

!==================================================================

subroutine MAPL_ReadForcing1(STATE,NAME,DATAFILE,CURRENTTIME,    &
                             FORCING,INIT_ONLY,ON_TILES,RC )

   type (MAPL_MetaComp),     intent(INOUT)   :: STATE
   character(len=*),         intent(IN   )   :: NAME
   character(len=*),         intent(IN   )   :: DATAFILE
   type (ESMF_Time),         intent(IN   )   :: CURRENTTIME
   real,                     intent(  OUT)   :: FORCING(:)   
   logical, optional,        intent(IN   )   :: INIT_ONLY
   logical, optional,        intent(IN   )   :: ON_TILES
   integer, optional,        intent(  OUT)   :: RC

   character(len=ESMF_MAXSTR)        :: IAm = "MAPL_ReadForcing1"
   integer                           :: STATUS

   call MAPL_ReadForcingX(STATE,NAME,DATAFILE,CURRENTTIME,      &
                          FORCING1=FORCING,INIT_ONLY=INIT_ONLY, &
                          ON_TILES=ON_TILES,                    &
                          RC=STATUS )

   RETURN_(STATUS)
end subroutine MAPL_ReadForcing1

!==================================================================

subroutine MAPL_ReadForcing2(STATE,NAME,DATAFILE,CURRENTTIME,    &
                             FORCING,INIT_ONLY,RC )

   type (MAPL_MetaComp),     intent(INOUT)   :: STATE
   character(len=*),         intent(IN   )   :: NAME
   character(len=*),         intent(IN   )   :: DATAFILE
   type (ESMF_Time),         intent(IN   )   :: CURRENTTIME
   real,                     intent(  OUT)   :: FORCING(:,:)   
   logical, optional,        intent(IN   )   :: INIT_ONLY
   integer, optional,        intent(  OUT)   :: RC

   character(len=ESMF_MAXSTR)        :: IAm = "MAPL_ReadForcing2"
   integer                           :: STATUS

   call MAPL_ReadForcingX(STATE,NAME,DATAFILE,CURRENTTIME,      &
                          FORCING2=FORCING,INIT_ONLY=INIT_ONLY, &
                          ON_TILES=.FALSE.,                     &
                          RC=STATUS )
   RETURN_(STATUS)

end subroutine MAPL_ReadForcing2

!==================================================================

subroutine MAPL_ReadForcingX(STATE,NAME,DATAFILE,CURRTIME,  &
                             FORCING1,FORCING2,INIT_ONLY,ON_TILES,RC )

   type (MAPL_MetaComp),     intent(INOUT)   :: STATE
   character(len=*),         intent(IN   )   :: NAME
   character(len=*),         intent(IN   )   :: DATAFILE
   type (ESMF_Time),         intent(IN   )   :: CURRTIME
   real, optional,           intent(  OUT)   :: FORCING1(:)   
   real, optional,           intent(  OUT)   :: FORCING2(:,:)   
   logical, optional,        intent(IN   )   :: INIT_ONLY
   logical, optional,        intent(IN   )   :: ON_TILES
   integer, optional,        intent(  OUT)   :: RC

! ErrLog Variables

   character(len=ESMF_MAXSTR)        :: IAm = "MAPL_ReadForcing"
   integer                           :: STATUS

! Locals

   type (MAPL_LocStream)             :: LOCSTREAM
   type (ESMF_Time)                  :: DATE1
   type (ESMF_Time)                  :: DATEN
   type (ESMF_Calendar)              :: CAL
   type (ESMF_DELayout)              :: LAYOUT
   type (ESMF_Time)                  :: MIDT1
   type (ESMF_Time)                  :: MIDT2
   type (ESMF_Time)                  :: CURRENTTIME
   type (ESMF_Grid)                  :: GRID
   type (ESMF_Field)                 :: FIELD
   real                              :: FAC
   integer                           :: IYR, IMM, IDD, IHR, IMN, ISC
   real, pointer                     :: PRE1(:  ), NEX1(:  )
   real, pointer                     :: PRE2(:,:), NEX2(:,:)
   real, pointer                     :: LONS(:,:), VAR2(:,:)
   type(ESMF_TimeInterval)           :: TIMEDIFF
   integer, parameter                :: FILE_HEADER_SIZE=14
   real                              :: REAL_HEADER(FILE_HEADER_SIZE)
   integer                           :: HEADER(FILE_HEADER_SIZE)
   integer                           :: UNIT
   integer                           :: DATE_PREV, YEAR
   logical                           :: USE_INIT_ONLY
   logical                           :: USE_ON_TILES
   logical                           :: TRANSFORM
   logical                           :: ONED
   logical                           :: CLIMATOLOGY
   integer                           :: YY, MM, DD
   integer                           :: H, M, S
   integer                           :: NSUBTILES
   type(ESMF_Grid)                   :: TILEGRID

! Process arguments
!------------------

   if(present(INIT_ONLY)) then
      USE_INIT_ONLY = INIT_ONLY
   else
      USE_INIT_ONLY = .false.
   end if

   if(present(ON_TILES)) then
      USE_ON_TILES = ON_TILES
   else
      USE_ON_TILES = .false.
   end if

   if    (present(FORCING1)) then
      ONED = .TRUE.
   elseif(present(FORCING2)) then
      ONED = .FALSE.
   else
      ASSERT_(.FALSE.)
   end if

! Get parameters from generic state.
!-----------------------------------

   call MAPL_GenericStateGet(STATE,             &
        LOCSTREAM=LOCSTREAM,                    &
        LONS=LONS,                              &
        LAYOUT=LAYOUT,                          &
                                      RC=STATUS )
   VERIFY_(STATUS)

! Get the calendar from the current time
!---------------------------------------

   call ESMF_TimeGet(CurrTime, calendar=cal, rc=status)
   VERIFY_(STATUS)

! Scratch space for reading data files
!-------------------------------------

    TRANSFORM = ONED .and. .not.USE_ON_TILES

    if(TRANSFORM) then
       allocate(VAR2(size(LONS,1),size(LONS,2)),STAT=STATUS)
       VERIFY_(STATUS)
    end if

! Get the grid form generic state 
!-------------------------------------------------------

   GRID=STATE%GRID%ESMFGRID

! Check if FRCSTATE already has previous and next vars. If not, create them
!---------------------------------------------------------------------------
   call ESMF_StateGetField(STATE%FORCING, trim(NAME)//'_PREV', FIELD, RC=STATUS)
   if (STATUS /= ESMF_SUCCESS) then
      if(ONED) then
         NSUBTILES = MAPL_GetNumSubtiles(STATE, RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_LocStreamGet(STATE%LocStream, TILEGRID=TILEGRID, RC=STATUS)
         VERIFY_(STATUS)
         
         ! create a spec
         call MAPL_VarSpecCreateInList(STATE%FRCSPEC,                   &
              SHORT_NAME = trim(NAME) // '_PREV',                       &
              LONG_NAME  = 'previos value of ' // trim(NAME),           &
              NUM_SUBTILES=NSUBTILES,                                   &
              DIMS       = MAPL_DimsTileOnly,                           &
              VLOCATION  = MAPL_VLocationNone,                          &
              COUPLE_INTERVAL = -1,                                     & ! time not set
              RC=STATUS  )
         VERIFY_(STATUS)
         
         call MAPL_VarSpecCreateInList(STATE%FRCSPEC,                   &
              SHORT_NAME = trim(NAME) // '_NEXT',                       &
              LONG_NAME  = 'next value of ' // trim(NAME),              &
              NUM_SUBTILES=NSUBTILES,                                   &
              DIMS       = MAPL_DimsTileOnly,                           &
              VLOCATION  = MAPL_VLocationNone,                          &
              COUPLE_INTERVAL = -1,                                     & ! time not set
              RC=STATUS  )
         VERIFY_(STATUS)

         ! create field and put it in FRCSTATE
         call MAPL_StateCreateFromSpec(STATE%FORCING,STATE%FRCSPEC,GRID,TILEGRID=TILEGRID,RC=STATUS)
      else
         call MAPL_VarSpecCreateInList(STATE%FRCSPEC,                   &
              SHORT_NAME = trim(NAME) // '_PREV',                       &
              LONG_NAME  = 'previos value of ' // trim(NAME),           &
              DIMS       = MAPL_DimsHorzOnly,                           &
              VLOCATION  = MAPL_VLocationNone,                          &
              COUPLE_INTERVAL = -1,                                     & ! time not set
              RC=STATUS  )
         VERIFY_(STATUS)
         
         call MAPL_VarSpecCreateInList(STATE%FRCSPEC,                   &
              SHORT_NAME = trim(NAME) // '_NEXT',                       &
              LONG_NAME  = 'next value of ' // trim(NAME),              &
              DIMS       = MAPL_DimsHorzOnly,                           &
              VLOCATION  = MAPL_VLocationNone,                          &
              COUPLE_INTERVAL = -1,                                     & ! time not set
              RC=STATUS  )
         VERIFY_(STATUS)

            ! create field and put it in FRCSTATE
         call MAPL_StateCreateFromSpec(STATE%FORCING,STATE%FRCSPEC,GRID,RC=STATUS)
      end if
   end if

! Get pointers to the endpoints in the forcing state
!----------------------------------------------------

   if(ONED) then
      call GET_POINTER(STATE%FORCING, PRE1, trim(NAME)//'_PREV', RC=STATUS)
      VERIFY_(STATUS)
      call GET_POINTER(STATE%FORCING, NEX1, trim(NAME)//'_NEXT', RC=STATUS)
      VERIFY_(STATUS)
   else
      call GET_POINTER(STATE%FORCING, PRE2, trim(NAME)//'_PREV', RC=STATUS)
      VERIFY_(STATUS)
      call GET_POINTER(STATE%FORCING, NEX2, trim(NAME)//'_NEXT', RC=STATUS)
      VERIFY_(STATUS)
   end if

! Set time endpoints. These are 0000-01-01:000000 if uninitialized
!----------------------------------------------------

   call MAPL_StateGetTimeStamp(STATE,trim(NAME)//'_NEXT',MIDT2, RC=STATUS )
   VERIFY_(STATUS)
   call MAPL_StateGetTimeStamp(STATE,trim(NAME)//'_PREV',MIDT1, RC=STATUS )
   VERIFY_(STATUS)

! Check if the input file is a climatology
!-----------------------------------------

   call MAPL_StateGetSpecAttrib(STATE,trim(NAME)//'_PREV',forcing=.true., &
        REFRESH_INTERVAL=DATE_PREV, rc=STATUS)
   VERIFY_(STATUS)

   if ( DATE_PREV < 0 ) then
      UNIT=GETFILE(DATAFILE,form='unformatted')
      call READ_PARALLEL(LAYOUT, REAL_HEADER, unit=UNIT)
      HEADER = nint(REAL_HEADER)
      call MAPL_BACKSPACE(UNIT,LAYOUT,RC=STATUS)
      VERIFY_(STATUS)
      CLIMATOLOGY = HEADER(1)==0
   else
      call ESMF_TimeGet(MIDT1, YY=YEAR, rc=status)
      VERIFY_(STATUS)
      CLIMATOLOGY = YEAR < 3
   end if

! Make a local copy of the current time
!--------------------------------------

!ALT clock are shallow objects!!!    CurrentTime = CurrTime
      call ESMF_TimeGet(CurrTime, &
                        YY=YY, MM=MM, DD=DD, &
                        H=H, M=M, S=S, rc=status)
      VERIFY_(STATUS)

      call ESMF_TimeSet(CurrentTime, &
                        YY=YY, MM=MM, DD=DD, &
                        H=H, M=M, S=S, rc=status)
      VERIFY_(STATUS)


   if(CLIMATOLOGY) then

!ALT: In Climatology mode we should not have leap year
!     so if the date is Feb.29, move it to Feb.28

      if (MM==2 .and. DD==29) DD=28

      call ESMF_TimeSet(CurrentTime, &
                        YY=1, MM=MM, DD=DD, &
                        H=H, M=M, S=S, rc=status)
      VERIFY_(STATUS)

   end if

! Make sure current time is between the current endpoints;
!  if not, refresh endpoints. This also works initially.
!-------------------------------------------------------

   if (CurrentTime < MIDT1 .or. CurrentTime > MIDT2) then
      call UPDATE_ENDPOINTS
   endif

! Interpolate ocean and ice temp and ice fraction
!   and update the forcing state on the tile grid.
!--------------------------------------------------

   if(USE_INIT_ONLY) then
      FAC=1.0
   else
      call MAPL_Interp_Fac (CurrentTime,MIDT1,MIDT2,FAC,RC=STATUS )
      VERIFY_(STATUS)
   end if

   ASSERT_(FAC >= 0.0)
   ASSERT_(FAC <= 1.0)

!  Update the friendly skin values
!---------------------------------

   if(ONED) then
      FORCING1 = FAC * PRE1  +  (1.0-FAC) * NEX1
   else
      FORCING2 = FAC * PRE2  +  (1.0-FAC) * NEX2
   end if

   if(TRANSFORM) deallocate(VAR2)

!  All done
!-----------

   RETURN_(ESMF_SUCCESS)

 contains

   subroutine UPDATE_ENDPOINTS

! Get forcing fortran unit
!-------------------------

    UNIT=GETFILE(DATAFILE,form='unformatted')

! Get the previous date
!----------------------

    call MAPL_StateGetSpecAttrib(STATE,trim(NAME)//'_PREV',forcing=.true., &
        REFRESH_INTERVAL=DATE_PREV, rc=STATUS)
    VERIFY_(STATUS)

! Check to see if forcing state buffers have been initialized
!-------------------------------------------------------------

    if ( DATE_PREV < 0 .or. CLIMATOLOGY) then

! If not, initialize them, by looping until correct times are found
!-----------------------------------------------------------------

       if(MAPL_AM_I_ROOT(LAYOUT)) then
          rewind(UNIT)
          do
             read(UNIT) REAL_HEADER
             HEADER = nint(REAL_HEADER)
             call ESMF_TimeSet(DATEN, &
                              YY=HEADER(7), MM=HEADER(8), DD=HEADER(9), &
                              H=HEADER(10), M=HEADER(11), S=HEADER(12), &
                              calendar=cal, rc=status)
             VERIFY_(STATUS)
    
! If the current time is beyond the item's final time, skip it.
!--------------------------------------------------------------

             if (DATEN < CurrentTime) then
                read(UNIT)
                cycle
             else
                backspace (UNIT)
                if(use_init_only) then
                   backspace (UNIT)
                   backspace (UNIT)
                end if
                exit
             end if
          end do
       end if

! We have found an interval that contains Current. Now get its initial time
!--------------------------------------------------------------------------

       call READ_PARALLEL(Layout, REAL_HEADER, unit=UNIT)
       HEADER = nint(REAL_HEADER)

       call ESMF_TimeSet(DATEN, &
                         YY=HEADER(7), MM=HEADER(8), DD=HEADER(9), &
                         H=HEADER(10), M=HEADER(11), S=HEADER(12), &
                         calendar=cal, rc=status)
       VERIFY_(STATUS)
    
       call ESMF_TimeSet(DATE1, &
                         YY=HEADER(1), MM=HEADER(2), DD=HEADER(3), &
                         H =HEADER(4), M =HEADER(5), S =HEADER(6), &
                         calendar=cal, rc=status)
       VERIFY_(STATUS)

! compute its central time
!-------------------------

       TimeDiff=DATEN-DATE1
       MIDT2=DATE1 + TimeDiff/2.0D0

       if(MIDT2<=CurrentTime) then ! The item we found is PREV
          MIDT1 = MIDT2
          
          call WRITE_PARALLEL("Previous time for"//trim(NAME)//"s is:", RC=STATUS)
          VERIFY_(STATUS)
          call WRITE_PARALLEL(HEADER(1:12), RC=STATUS)
          VERIFY_(STATUS)

! Read PREV
!----------

          call READIT('_PREV')

          call MAPL_StateSetTimeStamp(STATE,trim(NAME)//'_PREV',MIDT1, RC=STATUS )
          VERIFY_(STATUS)

! Read the header for NEXT
!-------------------------

          call READ_PARALLEL(Layout, REAL_HEADER, unit=UNIT)
          HEADER = nint(REAL_HEADER)

! Get NEXT's initial and final times
!-----------------------------------

          call ESMF_TimeSet(DATEN, &
                            YY=HEADER(7), MM=HEADER(8), DD=HEADER(9), &
                            H=HEADER(10), M=HEADER(11), S=HEADER(12), &
                            calendar=cal, rc=status)
          VERIFY_(STATUS)

          call ESMF_TimeSet(DATE1, &
                            YY=HEADER(1), MM=HEADER(2), DD=HEADER(3), &
                            H =HEADER(4), M =HEADER(5), S =HEADER(6), &
                            calendar=cal, rc=status)
          VERIFY_(STATUS)

! compute its central time
!-------------------------

          TimeDiff=DATEN-DATE1
          MIDT2=DATE1 + TimeDiff/2.0D0

! and read NEXT
!--------------

          call READIT('_NEXT')

          call MAPL_StateSetTimeStamp(STATE,trim(NAME)//'_NEXT',MIDT2, RC=STATUS )
          VERIFY_(STATUS)

       else

! Read NEXT
!----------

          call READIT('_NEXT')

          call MAPL_StateSetTimeStamp(STATE,trim(NAME)//'_NEXT',MIDT2, RC=STATUS )
          VERIFY_(STATUS)

! Back up to get PREV
!--------------------

          call MAPL_Backspace(UNIT, LAYOUT, COUNT=4, RC=STATUS); VERIFY_(STATUS)

! Read the header of PREV item
!-----------------------------

          call READ_PARALLEL(Layout, REAL_HEADER, unit=UNIT)
          HEADER = nint(REAL_HEADER)

! Get the item's initial and final times
!---------------------------------------

          call ESMF_TimeSet(DATEN, &
                            YY=HEADER(7), MM=HEADER(8), DD=HEADER(9), &
                            H=HEADER(10), M=HEADER(11), S=HEADER(12), &
                            calendar=cal, rc=status)
          VERIFY_(STATUS)

          call ESMF_TimeSet(DATE1, &
                            YY=HEADER(1), MM=HEADER(2), DD=HEADER(3), &
                            H =HEADER(4), M =HEADER(5), S =HEADER(6), &
                            calendar=cal, rc=status)
          VERIFY_(STATUS)

! compute its central time
!-------------------------

          TimeDiff=DATEN-DATE1
          MIDT1=DATE1 + TimeDiff/2.0D0

! Read PREV
!----------

          call READIT('_PREV')
          
          call MAPL_StateSetTimeStamp(STATE,trim(NAME)//'_PREV',MIDT1, RC=STATUS )
          VERIFY_(STATUS)

! Skip over NEXT to be positioned for subsequent reads
!-----------------------------------------------------

          call MAPL_Skip(UNIT, LAYOUT, COUNT=2, RC=STATUS ); VERIFY_(STATUS)

       end if
       
    elseif(.not.use_init_only) then ! Just need to update NEXT,  PREV is old NEXT

       if(ONED) then
          PRE1   = NEX1
       else
          PRE2   = NEX2
       end if

       MIDT1  = MIDT2

! Move time stamp from NEXT to PREV
!------------------------------------------

       call MAPL_StateGetTimeStamp(STATE,trim(NAME)//'_NEXT',MIDT2, RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_StateSetTimeStamp(STATE,trim(NAME)//'_PREV',MIDT2, RC=STATUS )
       VERIFY_(STATUS)

! Read the header of next item
!-----------------------------

       call READ_PARALLEL(Layout, REAL_HEADER, unit=UNIT, rc=status)
       VERIFY_(STATUS)
       HEADER = nint(REAL_HEADER)

! Get the item's initial and final time
!--------------------------------------
 
       call ESMF_TimeSet(DATE1, &
                        YY=HEADER(1), MM=HEADER(2), DD=HEADER(3), &
                        H =HEADER(4), M =HEADER(5), S =HEADER(6), &
                        calendar=cal, rc=status)
       VERIFY_(STATUS)
 
       call ESMF_TimeSet(DATEN, &
                        YY=HEADER(7), MM=HEADER(8), DD=HEADER(9), &
                        H =HEADER(10),M =HEADER(11),S =HEADER(12),&
                        calendar=cal, rc=status)
       VERIFY_(STATUS)

       TimeDiff=DATEN-DATE1
       MIDT2=DATE1 + TimeDiff/2.0D0

! Verify that it is the next item
!--------------------------------

       ASSERT_(MIDT2 >= CurrentTime)

! Read NEXT
!----------

       call READIT('_NEXT')

       call MAPL_StateSetTimeStamp(STATE,trim(NAME)//'_NEXT',MIDT2, RC=STATUS )
       VERIFY_(STATUS)

    endif

  end subroutine UPDATE_ENDPOINTS

  subroutine READIT(WHICH)
    character(len=5), intent(IN) :: WHICH
    real, pointer :: VAR1(:)

    if(TRANSFORM) then
       call MAPL_VarRead(UNIT, GRID, VAR2, RC=STATUS )
       VERIFY_(STATUS)
       call GET_POINTER(STATE%FORCING, VAR1, trim(NAME)//WHICH, RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_LocStreamTransform(LOCSTREAM, VAR1, VAR2, interp=.true.,RC=STATUS )
       VERIFY_(STATUS)
    else
       call MAPL_VarRead(UNIT, STATE%FORCING, trim(NAME)//WHICH, RC=STATUS)
       VERIFY_(STATUS)
    end if

  end subroutine READIT

end subroutine MAPL_READFORCINGX

!==================================================================

  subroutine MAPL_StateGetTimeStamp(STATE,NAME,TIME,RC)

    type (MAPL_MetaComp),     intent(INOUT)   :: STATE
    character(len=*),         intent(IN   )   :: NAME
    type (ESMF_Time),         intent(  OUT)   :: TIME
    integer, optional,        intent(  OUT)   :: RC

! ErrLog Variables

    character(len=ESMF_MAXSTR)        :: IAm = "MAPL_StateGetTimeStamp"
    integer                           :: STATUS

! Locals

    integer                           :: HOUR
    integer                           :: DATE
    integer                           :: IYR, IMM, IDD, IHR, IMN, ISC


    call MAPL_StateGetSpecAttrib(STATE,NAME,forcing=.true., &
        refresh_interval=DATE, averaging_interval=HOUR, RC=STATUS )
    VERIFY_(STATUS)

    if (DATE <= 0) then
       IYR=0; IMM = 1; IDD = 1
       IHR=0; IMN = 0; ISC = 0
    else
       call MAPL_UNPACKTIME(DATE,IYR,IMM,IDD)
       call MAPL_UNPACKTIME(HOUR,IHR,IMN,ISC)
    endif

    call ESMF_TimeSet(TIME, YY=IYR, MM=IMM, DD=IDD, H=IHR, M=IMN, S=ISC, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_StateGetTimeStamp

!==================================================================

  subroutine MAPL_StateSetTimeStamp(STATE,NAME,TIME,RC)

    type (MAPL_MetaComp),     intent(INOUT)   :: STATE
    character(len=*),         intent(IN   )   :: NAME
    type (ESMF_Time),         intent(IN   )   :: TIME
    integer, optional,        intent(  OUT)   :: RC

! ErrLog Variables

    character(len=ESMF_MAXSTR)        :: IAm = "MAPL_StateSetTimeStamp"
    integer                           :: STATUS

! Locals

    integer                           :: HOUR
    integer                           :: DATE
    integer                           :: IYR, IMM, IDD, IHR, IMN, ISC

    call ESMF_TimeGet(TIME ,YY=IYR, MM=IMM, DD=IDD, H=IHR, M=IMN, S=ISC, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_PackTime(DATE,IYR,IMM,IDD)
    call MAPL_PackTime(HOUR,IHR,IMN,ISC)

    call MAPL_StateSetSpecAttrib(STATE,NAME,forcing=.true., &
         refresh_interval=DATE, &
         averaging_interval=HOUR, rc=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine MAPL_StateSetTimeStamp

  subroutine MAPL_GenericMakeXchgNatural(STATE, RC)
    type (MAPL_MetaComp),     intent(INOUT)   :: STATE
    integer, optional,        intent(  OUT)   :: RC

! ErrLog Variables

    character(len=ESMF_MAXSTR)        :: IAm = "MAPL_GenericMakeXchgNatural"
    integer                           :: STATUS


    STATE%LOCSTREAM = STATE%ExchangeGrid

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GenericMakeXchgNatural

  
  subroutine MAPL_GridCreate(GC, rc)

    type(ESMF_GridComp),    intent(INOUT) :: GC
    integer, optional,      intent(OUT)   :: rc

    integer                               :: status
    character(len=ESMF_MAXSTR)            :: Comp_Name
    character(len=ESMF_MAXSTR)            :: IAm

    type (ESMF_VM)                        :: VM
    type (ESMF_DELayout)                  :: LAYOUT
    type (MAPL_MetaComp), pointer         :: STATE
    type (MAPL_LocStream)                 :: EXCH
    type (ESMF_Grid)                      :: GRID
    logical                               :: EXACT
    integer                               :: IM_WORLD
    integer                               :: JM_WORLD
    integer                               :: LM,L,NN
    integer                               :: NX, NY
    integer, allocatable                  :: IMS(:), JMS(:)
    character(len=ESMF_MAXSTR)            :: Gridname
    character(len=4)                      :: imsz
    character(len=4)                      :: jmsz
    character(len=2)                      :: date
    character(len=2)                      :: pole
    real(ESMF_KIND_R8)                    :: minCoord(2)
    real(ESMF_KIND_R8)                    :: deltaX, deltaY, deltaZ

    ! Query GC
    !---------

    Iam='AppGridCreate'
    call ESMF_GridCompGet( GC, vm=vm, name=Comp_Name,   rc = status )
    VERIFY_(STATUS)
    Iam = trim(Comp_Name)//Iam

    call ESMF_VMGetCurrent(vm, rc=status)
    VERIFY_(STATUS)

    ! Get MAPL object
    !----------------

    call MAPL_InternalStateGet(GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    ! Make layout  (Required in Config)
    !----------------------------------

    call MAPL_GetResource( STATE, NX, 'NX:', rc = status )
    VERIFY_(STATUS)
    call MAPL_GetResource( STATE, NY, 'NY:', rc = status )
    VERIFY_(STATUS)

    layout = ESMF_DELayoutCreate(vm, deCountList=(/NX, NY/), rc=status)
    VERIFY_(STATUS)

    ! Get the grid name (Required in Config)
    !---------------------------------------

    call MAPL_GetResource ( STATE, Gridname, 'GRIDNAME:', rc=status )
    VERIFY_(STATUS)

    ! Parse name for grid info 
    !-------------------------

    Gridname = AdjustL(Gridname)
    nn   = len_trim(Gridname)
    imsz = Gridname(3:index(Gridname,'x')-1)
    jmsz = Gridname(index(Gridname,'x')+1:nn-3)
    pole = Gridname(1:2)
    date = Gridname(nn-1:nn)

    ! Get grid sizes from grid name
    !------------------------------

    read(IMSZ,*) IM_WORLD
    read(JMSZ,*) JM_WORLD

    call MAPL_GetResource( STATE, LM, 'LM:', rc = status )
    VERIFY_(STATUS)

    ! Distribute the horizontal dimensions
    !-------------------------------------

    allocate( IMS(0:NX-1), stat=STATUS ); VERIFY_(STATUS)
    allocate( JMS(0:NY-1), stat=STATUS ); VERIFY_(STATUS)

    call MAPL_GetResource( STATE, ims, 'IMS:', rc = status )
    if(STATUS/=ESMF_SUCCESS) then
       call MAPL_DecomposeDim ( IM_WORLD,ims,nx )
    end if

    call MAPL_GetResource( STATE, jms, 'JMS:', rc = status )
    if(STATUS/=ESMF_SUCCESS) then
       call MAPL_DecomposeDim ( JM_WORLD,jms,ny )
    end if

    ! Get the type of grid from the grid name.
    !  Locations of uniform lat-lon grids are done exactly,
    !  all others are done approximately (within roundoff)
    !  from the exchange grid. 
    !-----------------------------------------------------


    EXACT = (pole=='PE' .or. pole=='PC' )                         .and. &
         (date=='DC' .or. date=='DE' .or. date=='GC' .or. date=='GC' )

    if(EXACT) then
       deltaX = (2.0*MAPL_PI)/IM_WORLD
       if (trim(POLE)=='PE') then
          deltaY = MAPL_PI/JM_WORLD
       else
          ASSERT_(JM_WORLD > 1)
          deltaY = MAPL_PI/(JM_WORLD-1)
       end if
       deltaZ = 1.0D0

       if(date(1:1)=='G') then
          MinCoord(1) = 0.0
       else
          MinCoord(1) = -MAPL_PI
       end if

       MinCoord(2) = -0.5*MAPL_PI

       if( date(2:2) /= 'E' )  MinCoord(1) = minCoord(1) - deltaX/2
       if( pole(2:2) /= 'E' )  MinCoord(2) = minCoord(2) - deltaY/2
    else
       deltaX   = 1.0
       deltaY   = 1.0
       deltaZ   = 1.0
       MinCoord = 0.0
    end if

    grid = ESMF_GridCreateHorzLatLonUni(           &
         counts      = (/IM_WORLD, JM_WORLD/),     &
         minGlobalCoordPerDim = minCoord,          &
         deltaPerDim = (/deltaX, deltaY /),        &
         horzStagger = ESMF_Grid_Horz_Stagger_A,   &
         periodic    = (/ESMF_TRUE, ESMF_FALSE/),  &
         name        = gridname,         rc=status )
    VERIFY_(STATUS)

    call ESMF_GridAddVertHeight(grid,              &
         delta       = (/(deltaZ, L=1,LM) /),      &
         vertStagger = ESMF_GRID_VERT_STAGGER_TOP, &
         rc=status )
    VERIFY_(STATUS)

    call ESMF_GridDistribute(grid,                 &
         deLayout        = layout,                 &
         countsPerDEDim1 = ims,                    &
         countsPerDEDim2 = jms,                    &
         rc=status )
    VERIFY_(STATUS)

    call ESMF_GridCompSet(GC, GRID=GRID, RC=STATUS)
    VERIFY_(STATUS)

    deallocate(ims)
    deallocate(jms)

    if(.not. EXACT) then

       ! get exchange
       EXCH = STATE%EXCHANGEGRID

       call MAPL_GridCoordAdjust(GRID, EXCH, RC=STATUS)
       VERIFY_(STATUS)
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GridCreate

  recursive subroutine MAPL_GetRootGC(GC, rootGC, RC)
    type(ESMF_GridComp),    intent(IN   ) :: GC
    type(ESMF_GridComp),    intent(  OUT) :: rootGC
    integer, optional,      intent(OUT)   :: rc

    integer                               :: status
    character(len=ESMF_MAXSTR)            :: IAm
    type (MAPL_MetaComp),     pointer     :: META 

    call MAPL_GetObjectFromGC(GC, META, RC=STATUS)
    VERIFY_(STATUS)

    if (.not. associated(META%parentGC)) then
       rootGC = GC
    else
       call MAPL_GetRootGC(META%parentGC, rootGC, RC=status)
       VERIFY_(STATUS)
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GetRootGC

  subroutine MAPL_GridCompSetEntryPoint(GC, subroutineType,  subroutineName, RC)
    type(ESMF_GridComp),                  intent(INOUT) :: GC         ! Gridded component
    character(len=*),                     intent(IN   ) :: subroutineType
    external                                            :: subroutineName
    integer,                    optional, intent(  OUT) :: RC         ! Return code

    integer                               :: status
    character(len=ESMF_MAXSTR)            :: IAm

    type (MAPL_MetaComp),     pointer     :: META 
    integer                               :: phase

    call MAPL_InternalStateRetrieve( GC, META, RC=STATUS)
    VERIFY_(STATUS)

    select case (subroutineType)
    case(ESMF_SETINIT)
       phase = MAPL_AddMethod(META%phase_init, RC=STATUS)
    case(ESMF_SETRUN)
       phase = MAPL_AddMethod(META%phase_run, RC=STATUS)
    case(ESMF_SETFINAL)
       phase = MAPL_AddMethod(META%phase_final, RC=STATUS)
    case("ESMF_WriteRestart") ! this should've been ESMF_SETWRITERESTART
                              ! it is used internally in ESMF but
                              ! it is not a public parameter
       phase = MAPL_AddMethod(META%phase_record, RC=STATUS)
! kludge
       call ESMF_GridCompSetEntryPoint(GC, ESMF_SETFINAL,  subroutineName, &
                                    MAPL_RecordPhase, status)
       RETURN_(STATUS)
    case("ESMF_ReadRestart")  ! this should've been ESMF_SETREADRESTART
                              ! it is used internally in ESMF but
                              ! it is not a public parameter
       phase = MAPL_AddMethod(META%phase_coldstart, RC=STATUS)
! kludge
       call ESMF_GridCompSetEntryPoint(GC, ESMF_SETINIT,  subroutineName, &
                                    MAPL_ColdstartPhase, status)
       RETURN_(STATUS)
    case default
       RETURN_(ESMF_FAILURE)
    end select
    VERIFY_(STATUS)

    call ESMF_GridCompSetEntryPoint(GC, subroutineType,  subroutineName, &
                                    phase, status)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_GridCompSetEntryPoint

  integer function MAPL_AddMethod(PHASE, RC)
    integer, pointer               :: PHASE(:)
    integer, optional, intent(out) :: rc

    integer :: I
    integer, pointer :: tmp(:)
    integer :: status
    character(len=ESMF_MAXSTR), parameter :: Iam="MAPL_AddMethod"

    MAPL_AddMethod = ESMF_SINGLEPHASE
    if (.not.associated(PHASE)) then
     ! this is the method to be added
       I = 1
       allocate(PHASE(I), stat=status)
       VERIFY_(STATUS)
       PHASE(I) = ESMF_SINGLEPHASE

    else
       I = size(PHASE) + 1
       allocate(TMP(I), stat=status)
       VERIFY_(STATUS)
       TMP(1:I-1) = PHASE
       TMP(I) = TMP(I-1)+1

       deallocate(PHASE)
       PHASE => TMP
    end if
    MAPL_AddMethod = PHASE(I)

    RETURN_(ESMF_SUCCESS)
  end function MAPL_AddMethod


end module MAPL_GenericMod
