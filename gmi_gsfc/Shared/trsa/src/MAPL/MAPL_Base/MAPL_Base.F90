! $Id: MAPL_Base.F90,v 1.2 2011-08-09 22:12:59 mrdamon Exp $
#include "MAPL_ErrLog.h"

module MAPL_BaseMod

  use ESMF_Mod

  implicit none
  private

!=============================================================================
!BOP

! !MODULE: -- A container module for global constants


! !PUBLIC PARAMETERS

integer, public, parameter :: MAPL_CplUNKNOWN        = 0
integer, public, parameter :: MAPL_CplSATISFIED      = 1
integer, public, parameter :: MAPL_CplNEEDED         = 2
integer, public, parameter :: MAPL_CplNOTNEEDED      = 4
integer, public, parameter :: MAPL_FriendlyVariable  = 8
integer, public, parameter :: MAPL_FieldItem         = 8
integer, public, parameter :: MAPL_BundleItem        = 16
integer, public, parameter :: MAPL_NoRestart         = 32

integer, public, parameter :: MAPL_VLocationNone   = 0
integer, public, parameter :: MAPL_VLocationEdge   = 1
integer, public, parameter :: MAPL_VLocationCenter = 2

integer, public, parameter :: MAPL_DimsUnknown     = 0
integer, public, parameter :: MAPL_DimsVertOnly    = 1
integer, public, parameter :: MAPL_DimsHorzOnly    = 2
integer, public, parameter :: MAPL_DimsHorzVert    = 3
integer, public, parameter :: MAPL_DimsTileOnly    = 4
integer, public, parameter :: MAPL_DimsTileTile    = 5

integer, public, parameter :: MAPL_DuplicateEntry  = -99
integer, public, parameter :: MAPL_Self = 0 
integer, public, parameter :: MAPL_Import = 1
integer, public, parameter :: MAPL_Export = 2
integer, public, parameter :: MAPL_ConnUnknown = -1
integer, public, parameter :: MAPL_RecordPhase  = 99
integer, public, parameter :: MAPL_ColdstartPhase  = 99
integer, public, parameter :: MAPL_FirstPhase   = 81
integer, public, parameter :: MAPL_SecondPhase  = MAPL_FirstPhase+1
integer, public, parameter :: MAPL_ThirdPhase   = MAPL_FirstPhase+2
integer, public, parameter :: MAPL_FourthPhase  = MAPL_FirstPhase+3
integer, public, parameter :: MAPL_FifthPhase   = MAPL_FirstPhase+4

real,    public, parameter :: MAPL_UNDEF              = 1.0e15  

integer, public, parameter :: MAPL_Ocean              = 0
integer, public, parameter :: MAPL_Lake               = 19
integer, public, parameter :: MAPL_LandIce            = 20

integer, public, parameter :: MAPL_BroadleafEvergreen = 1
integer, public, parameter :: MAPL_BroadleafDeciduous = 2
integer, public, parameter :: MAPL_Needleleaf         = 3
integer, public, parameter :: MAPL_GroundCover        = 4
integer, public, parameter :: MAPL_BroadleafShrubs    = 5
integer, public, parameter :: MAPL_Tundra             = 6
integer, public, parameter :: MAPL_BareSoil           = 7
integer, public, parameter :: MAPL_Desert             = 8
integer, public, parameter :: MAPL_NumVegTypes        = 8

integer, public, parameter :: MAPL_Land               = 100
integer, public, parameter :: MAPL_Vegetated          = 101


! !PUBLIC VARIABLES:

!EOP

public MAPL_ArrayF90Deallocate

type WRAP1R4
   real(kind=4), dimension(:)    , pointer :: ptr
end type WRAP1R4

type WRAP2R4
   real(kind=4), dimension(:,:)    , pointer :: ptr
end type WRAP2R4

type WRAP3R4
   real(kind=4), dimension(:,:,:)    , pointer :: ptr
end type WRAP3R4

type WRAP4R4
   real(kind=4), dimension(:,:,:,:)    , pointer :: ptr
end type WRAP4R4

type WRAP1R8
   real(kind=8), dimension(:)    , pointer :: ptr
end type WRAP1R8

type WRAP2R8
   real(kind=8), dimension(:,:)    , pointer :: ptr
end type WRAP2R8

type WRAP3R8
   real(kind=8), dimension(:,:,:)    , pointer :: ptr
end type WRAP3R8

type WRAP4R8
   real(kind=8), dimension(:,:,:,:)    , pointer :: ptr
end type WRAP4R8


public MAPL_RTRN
public MAPL_VRFY
public MAPL_ASRT

public MAPL_AllocateCoupling
public MAPL_ConnectCoupling
public MAPL_DecomposeDim
public MAPL_Interp_Fac
public MAPL_ClimInterpFac
public MAPL_PackTime
public MAPL_UnpackTime
public MAPL_TimeStringGet
public MAPL_FieldSetTime
public MAPL_FieldGetTime
public MAPL_tick
public MAPL_incymd
public MAPL_nhmsf
public MAPL_nsecf2
public MAPL_FieldCreate
public MAPL_RemapBounds

interface MAPL_FieldCreate
   module procedure MAPL_FieldCreateRename
   module procedure MAPL_FieldCreateRegrid
end interface

interface MAPL_FieldGetTime
   module procedure MAPL_GetFieldTimeFromField
   module procedure MAPL_GetFieldTimeFromState
end interface

interface MAPL_FieldSetTime
   module procedure MAPL_SetFieldTimeFromField
   module procedure MAPL_SetFieldTimeFromState
end interface

interface MAPL_AllocateCoupling
   module procedure MAPL_AllocateCouplingFromArray
   module procedure MAPL_AllocateCouplingFromField
end interface

interface MAPL_ConnectCoupling
   module procedure MAPL_ConnectCouplingFromArray
   module procedure MAPL_ConnectCouplingFromField
end interface

interface MAPL_RemapBounds
   module procedure MAPL_RemapBounds_3dr4
end interface

interface MAPL_VRFY
   module procedure MAPL_VRFY
   module procedure MAPL_VRFYt
end interface

interface MAPL_ASRT
   module procedure MAPL_ASRT
   module procedure MAPL_ASRTt
end interface

interface MAPL_RTRN
   module procedure MAPL_RTRN
   module procedure MAPL_RTRNt
end interface

contains
  subroutine MAPL_AllocateCouplingFromField(field, rc)
    type(ESMF_Field),  intent(INOUT) :: field
    integer, optional, intent(  OUT) :: rc             
    
    integer                                 :: status
    character(len=ESMF_MAXSTR), parameter   :: IAm='MAPL_AllocateCouplingFromField'
    type(ESMF_Array)                        :: array
    
#ifdef ESMF_1_0_4
    call ESMF_FieldGetData (FIELD, ARRAY, RC=STATUS)
#else
    call ESMF_FieldGetArray (FIELD, ARRAY, RC=STATUS)
#endif
    VERIFY_(STATUS)
    
    call MAPL_AllocateCouplingFromArray(array, rc=STATUS)
    VERIFY_(STATUS)
    
    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_AllocateCouplingFromField

  subroutine MAPL_AllocateCouplingFromArray(array, rc)
    type(ESMF_Array),  intent(INOUT) :: array
    integer, optional, intent(  OUT) :: rc             
    
    integer                                 :: rank
    type(ESMF_DataType)                     :: type
    type(ESMF_DataKind)                     :: dk
    integer, dimension(ESMF_MAXDIM)         :: counts, lbounds, ubounds
    type(ESMF_Pointer)                      :: base
    integer                                 :: status
    character(len=ESMF_MAXSTR), parameter   :: IAm='MAPL_AllocateCouplingFromArray'

    real(kind=4), dimension(:)        , pointer :: r4d1
    real(kind=4), dimension(:,:)      , pointer :: r4d2
    real(kind=4), dimension(:,:,:)    , pointer :: r4d3
    real(kind=4), dimension(:,:,:,:)  , pointer :: r4d4

    real(kind=8), dimension(:)        , pointer :: r8d1
    real(kind=8), dimension(:,:)      , pointer :: r8d2
    real(kind=8), dimension(:,:,:)    , pointer :: r8d3
    real(kind=8), dimension(:,:,:,:)  , pointer :: r8d4

    type (WRAP1R4) :: wrap1dr4
    type (WRAP2R4) :: wrap2dr4
    type (WRAP3R4) :: wrap3dr4
    type (WRAP4R4) :: wrap4dr4

    type (WRAP1R8) :: wrap1dr8
    type (WRAP2R8) :: wrap2dr8
    type (WRAP3R8) :: wrap3dr8
    type (WRAP4R8) :: wrap4dr8

    call ESMF_ArrayGet(array, rank, kind=dk, &
                       counts = counts, lbounds=lbounds, ubounds=ubounds, &
                       base=base, rc=status)
    VERIFY_(STATUS)
!ALT in case the counts=0 emsf keeps ubounds=lbounds
    where (counts==0) ubounds = lbounds + counts - 1
    if (dk .eq. ESMF_R4) then
       if (rank .eq. 1) then
          call ESMF_ArrayGetData(array, r4d1, rc=status)
          VERIFY_(STATUS)
          if (.not. associated(r4d1)) then
             allocate(r4d1(lbounds(1):ubounds(1)),  stat=status)
             VERIFY_(STATUS)

             r4d1 = 0.0
             call c_ESMC_ArraySetBaseAddr(array, r4d1, status) 
             VERIFY_(STATUS)

             wrap1dr4%ptr => r4d1
             call c_ESMC_ArraySetF90Ptr(array, wrap1dr4, status) 
             VERIFY_(STATUS)
          endif

       else if (rank .eq. 2) then
          call ESMF_ArrayGetData(array, r4d2, rc=status)
          VERIFY_(STATUS)
          if (.not. associated(r4d2)) then
             allocate(r4d2(lbounds(1):ubounds(1), &
                           lbounds(2):ubounds(2)),  stat=status)
             VERIFY_(STATUS)

             r4d2 = 0.0
             call c_ESMC_ArraySetBaseAddr(array, r4d2, status) 
             VERIFY_(STATUS)

             wrap2dr4%ptr => r4d2
             call c_ESMC_ArraySetF90Ptr(array, wrap2dr4, status) 
             VERIFY_(STATUS)
          endif

       else if (rank .eq. 3) then
          call ESMF_ArrayGetData(array, r4d3, rc=status)
          VERIFY_(STATUS)
          if (.not. associated(r4d3)) then
             allocate(r4d3(lbounds(1):ubounds(1), &
                           lbounds(2):ubounds(2), &
                           lbounds(3):ubounds(3)), stat=status)
             VERIFY_(STATUS)

             r4d3 = 0.0
             call c_ESMC_ArraySetBaseAddr(array, r4d3, status) 
             VERIFY_(STATUS)

             wrap3dr4%ptr => r4d3
             call c_ESMC_ArraySetF90Ptr(array, wrap3dr4, status) 
             VERIFY_(STATUS)
          endif

       else if (rank .eq. 4) then
          call ESMF_ArrayGetData(array, r4d4, rc=status)
          VERIFY_(STATUS)
          if (.not. associated(r4d4)) then
             allocate(r4d4(lbounds(1):ubounds(1), &
                           lbounds(2):ubounds(2), &
                           lbounds(3):ubounds(3), &
                           lbounds(4):ubounds(4)), stat=status)
             VERIFY_(STATUS)

             r4d4 = 0.0
             call c_ESMC_ArraySetBaseAddr(array, r4d4, status) 
             VERIFY_(STATUS)

             wrap4dr4%ptr => r4d4
             call c_ESMC_ArraySetF90Ptr(array, wrap4dr4, status) 
             VERIFY_(STATUS)
          end if

       else
          RETURN_(ESMF_FAILURE)
       endif
    else if (dk .eq. ESMF_R8) then
       if (rank .eq. 1) then
          call ESMF_ArrayGetData(array, r8d1, rc=status)
          VERIFY_(STATUS)
          if (.not. associated(r8d1)) then
             allocate(r8d1(lbounds(1):ubounds(1)),  stat=status)
             VERIFY_(STATUS)

             r8d1 = 0.0
             call c_ESMC_ArraySetBaseAddr(array, r8d1, status) 
             VERIFY_(STATUS)

             wrap1dr8%ptr => r8d1
             call c_ESMC_ArraySetF90Ptr(array, wrap1dr8, status) 
             VERIFY_(STATUS)
          endif

       else if (rank .eq. 2) then
          call ESMF_ArrayGetData(array, r8d2, rc=status)
          VERIFY_(STATUS)
          if (.not. associated(r8d2)) then
             allocate(r8d2(lbounds(1):ubounds(1), &
                           lbounds(2):ubounds(2)),  stat=status)
             VERIFY_(STATUS)

             r8d2 = 0.0
             call c_ESMC_ArraySetBaseAddr(array, r8d2, status) 
             VERIFY_(STATUS)

             wrap2dr8%ptr => r8d2
             call c_ESMC_ArraySetF90Ptr(array, wrap2dr8, status) 
             VERIFY_(STATUS)
          endif

       else if (rank .eq. 3) then
          call ESMF_ArrayGetData(array, r8d3, rc=status)
          VERIFY_(STATUS)
          if (.not. associated(r8d3)) then
             allocate(r8d3(lbounds(1):ubounds(1), &
                           lbounds(2):ubounds(2), &
                           lbounds(3):ubounds(3)), stat=status)
             VERIFY_(STATUS)

             r8d3 = 0.0
             call c_ESMC_ArraySetBaseAddr(array, r8d3, status) 
             VERIFY_(STATUS)

             wrap3dr8%ptr => r8d3
             call c_ESMC_ArraySetF90Ptr(array, wrap3dr8, status) 
             VERIFY_(STATUS)
          endif

       else if (rank .eq. 4) then
          call ESMF_ArrayGetData(array, r8d4, rc=status)
          VERIFY_(STATUS)
          if (.not. associated(r8d4)) then
             allocate(r8d4(lbounds(1):ubounds(1), &
                           lbounds(2):ubounds(2), &
                           lbounds(3):ubounds(3), &
                           lbounds(4):ubounds(4)), stat=status)
             VERIFY_(STATUS)

             r8d4 = 0.0
             call c_ESMC_ArraySetBaseAddr(array, r8d4, status) 
             VERIFY_(STATUS)

             wrap4dr8%ptr => r8d4
             call c_ESMC_ArraySetF90Ptr(array, wrap4dr8, status) 
             VERIFY_(STATUS)
          end if

       else
          RETURN_(ESMF_FAILURE)
       endif
    else
       RETURN_(ESMF_FAILURE)
    endif

    RETURN_(ESMF_SUCCESS)


  end subroutine MAPL_AllocateCouplingFromArray





  subroutine MAPL_ConnectCouplingFromField(field, from_field, rc)
    type(ESMF_Field),  intent(INOUT) :: field
    type(ESMF_Field),  intent(IN   ) :: from_field
    integer, optional, intent(  OUT) :: rc             
    
    integer                                 :: status
    character(len=ESMF_MAXSTR), parameter   :: IAm='MAPL_ConnectCouplingFromField'
    type(ESMF_Array)                        :: array
    type(ESMF_Array)                        :: from_array
    
#ifdef ESMF_1_0_4
    call ESMF_FieldGetData (FIELD, ARRAY, RC=STATUS)
#else
    call ESMF_FieldGetArray (FIELD, ARRAY, RC=STATUS)
#endif
    VERIFY_(STATUS)

#ifdef ESMF_1_0_4
    call ESMF_FieldGetData (FROM_FIELD, FROM_ARRAY, RC=STATUS)
#else
    call ESMF_FieldGetArray (FROM_FIELD, FROM_ARRAY, RC=STATUS)
#endif
    VERIFY_(STATUS)
    
    call MAPL_ConnectCouplingFromArray(array, from_array, rc=STATUS)
    VERIFY_(STATUS)
    
    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_ConnectCouplingFromField


  subroutine MAPL_ConnectCouplingFromArray(array, from_array, rc)
    type(ESMF_Array),  intent(INOUT) :: array
    type(ESMF_Array),  intent(IN   ) :: from_array
    integer, optional, intent(  OUT) :: rc             
    
    integer                                 :: rank
    type(ESMF_DataType)                     :: type
    type(ESMF_DataKind)                     :: dk
    integer, dimension(ESMF_MAXDIM)         :: counts, lbounds, ubounds
    type(ESMF_Pointer)                      :: base
    integer                                 :: status
    character(len=ESMF_MAXSTR), parameter   :: IAm='MAPL_ConnectCouplingFromArray'

    real(kind=4), dimension(:)        , pointer :: p4d1
    real(kind=4), dimension(:,:)      , pointer :: p4d2
    real(kind=4), dimension(:,:,:)    , pointer :: p4d3
    real(kind=4), dimension(:,:,:,:)  , pointer :: p4d4

    real(kind=4), dimension(:)        , pointer :: r4d1
    real(kind=4), dimension(:,:)      , pointer :: r4d2
    real(kind=4), dimension(:,:,:)    , pointer :: r4d3
    real(kind=4), dimension(:,:,:,:)  , pointer :: r4d4

    real(kind=8), dimension(:)        , pointer :: r8d1
    real(kind=8), dimension(:,:)      , pointer :: r8d2
    real(kind=8), dimension(:,:,:)    , pointer :: r8d3
    real(kind=8), dimension(:,:,:,:)  , pointer :: r8d4

    type (WRAP1R4) :: wrap1dr4
    type (WRAP2R4) :: wrap2dr4
    type (WRAP3R4) :: wrap3dr4
    type (WRAP4R4) :: wrap4dr4

    type (WRAP1R8) :: wrap1dr8
    type (WRAP2R8) :: wrap2dr8
    type (WRAP3R8) :: wrap3dr8
    type (WRAP4R8) :: wrap4dr8

    call ESMF_ArrayGet(array, rank, kind=dk, &
                       counts = counts, lbounds=lbounds, ubounds=ubounds, &
                       base=base, rc=status)
    VERIFY_(STATUS)
    where (counts==0) ubounds = lbounds + counts - 1
    if (dk .eq. ESMF_R4) then
       if (rank .eq. 1) then
          call ESMF_ArrayGetData(from_array, p4d1, rc=status)
          VERIFY_(STATUS)
          call ESMF_ArrayGetData(array, r4d1, rc=status)
          VERIFY_(STATUS)
          if (associated(r4d1)) then
             deallocate(r4d1,  stat=status)
             VERIFY_(STATUS)
          endif

          r4d1 => p4d1
          call c_ESMC_ArraySetBaseAddr(array, r4d1, status) 
          VERIFY_(STATUS)

          wrap1dr4%ptr => r4d1
          call c_ESMC_ArraySetF90Ptr(array, wrap1dr4, status) 
          VERIFY_(STATUS)

       else if (rank .eq. 2) then
          call ESMF_ArrayGetData(from_array, p4d2, rc=status)
          VERIFY_(STATUS)
          call ESMF_ArrayGetData(array, r4d2, rc=status)
          VERIFY_(STATUS)
          if (associated(r4d2)) then
             deallocate(r4d2,  stat=status)
             VERIFY_(STATUS)
          endif

          r4d2 => p4d2
          call c_ESMC_ArraySetBaseAddr(array, r4d2, status) 
          VERIFY_(STATUS)

          wrap2dr4%ptr => r4d2
          call c_ESMC_ArraySetF90Ptr(array, wrap2dr4, status) 
          VERIFY_(STATUS)

       else if (rank .eq. 3) then
          call ESMF_ArrayGetData(from_array, p4d3, rc=status)
          VERIFY_(STATUS)
          call ESMF_ArrayGetData(array, r4d3, rc=status)
          VERIFY_(STATUS)
          if (associated(r4d3)) then
             deallocate(r4d3, stat=status)
             VERIFY_(STATUS)
          end if
          r4d3 => p4d3
          call c_ESMC_ArraySetBaseAddr(array, r4d3, status) 
          VERIFY_(STATUS)

          wrap3dr4%ptr => r4d3
          call c_ESMC_ArraySetF90Ptr(array, wrap3dr4, status) 
          VERIFY_(STATUS)

       else if (rank .eq. 4) then
          call ESMF_ArrayGetData(from_array, p4d4, rc=status)
          VERIFY_(STATUS)
          call ESMF_ArrayGetData(array, r4d4, rc=status)
          VERIFY_(STATUS)
          if (associated(r4d4)) then
             deallocate(r4d4, stat=status)
             VERIFY_(STATUS)
          end if
          r4d4 => p4d4
          call c_ESMC_ArraySetBaseAddr(array, r4d4, status) 
          VERIFY_(STATUS)

          wrap4dr4%ptr => r4d4
          call c_ESMC_ArraySetF90Ptr(array, wrap4dr4, status) 
          VERIFY_(STATUS)

       else
          RETURN_(ESMF_FAILURE)
       endif
    else if (dk .eq. ESMF_R8) then
!ALT: temporaty set to FAIL; if compiles OK, copy and paste from above; replace r4=>r8
       RETURN_(ESMF_FAILURE)
    else
       RETURN_(ESMF_FAILURE)
    endif

    RETURN_(ESMF_SUCCESS)


  end subroutine MAPL_ConnectCouplingFromArray

  subroutine MAPL_ArrayF90Deallocate(array, rc)
    type(ESMF_Array),  intent(INOUT) :: array
    integer, optional, intent(  OUT) :: rc             
    
    integer                                 :: rank
    type(ESMF_DataType)                     :: type
    type(ESMF_DataKind)                     :: dk
    integer, dimension(ESMF_MAXDIM)         :: counts, lbounds, ubounds
    type(ESMF_Pointer)                      :: base
    integer                                 :: status
    character(len=ESMF_MAXSTR), parameter   :: IAm='MAPL_ArrayF90Deallocate'

    real(kind=4), dimension(:)        , pointer :: r4d1
    real(kind=4), dimension(:,:)      , pointer :: r4d2
    real(kind=4), dimension(:,:,:)    , pointer :: r4d3
    real(kind=4), dimension(:,:,:,:)  , pointer :: r4d4

    real(kind=8), dimension(:)        , pointer :: r8d1
    real(kind=8), dimension(:,:)      , pointer :: r8d2
    real(kind=8), dimension(:,:,:)    , pointer :: r8d3
    real(kind=8), dimension(:,:,:,:)  , pointer :: r8d4

    type (WRAP1R4) :: wrap1dr4
    type (WRAP2R4) :: wrap2dr4
    type (WRAP3R4) :: wrap3dr4
    type (WRAP4R4) :: wrap4dr4

    type (WRAP1R8) :: wrap1dr8
    type (WRAP2R8) :: wrap2dr8
    type (WRAP3R8) :: wrap3dr8
    type (WRAP4R8) :: wrap4dr8

    call ESMF_ArrayGet(array, rank, kind=dk, &
                       counts = counts, lbounds=lbounds, ubounds=ubounds, &
                       base=base, rc=status)
    VERIFY_(STATUS)
    if (dk .eq. ESMF_R4) then
       if (rank .eq. 1) then
          call ESMF_ArrayGetData(array, r4d1, rc=status)
          VERIFY_(STATUS)
          if (associated(r4d1)) then
             deallocate(r4d1,  stat=status)
             VERIFY_(STATUS)

             wrap1dr4%ptr => r4d1
             call c_ESMC_ArraySetF90Ptr(array, wrap1dr4, status) 
             VERIFY_(STATUS)
          endif

       else if (rank .eq. 2) then
          call ESMF_ArrayGetData(array, r4d2, rc=status)
          VERIFY_(STATUS)
          if (associated(r4d2)) then
             deallocate(r4d2,  stat=status)
             VERIFY_(STATUS)

             wrap2dr4%ptr => r4d2
             call c_ESMC_ArraySetF90Ptr(array, wrap2dr4, status) 
             VERIFY_(STATUS)
          endif

       else if (rank .eq. 3) then
          call ESMF_ArrayGetData(array, r4d3, rc=status)
          VERIFY_(STATUS)
          if (associated(r4d3)) then
             deallocate(r4d3, stat=status)
             VERIFY_(STATUS)

             wrap3dr4%ptr => r4d3
             call c_ESMC_ArraySetF90Ptr(array, wrap3dr4, status) 
             VERIFY_(STATUS)
          endif

       else if (rank .eq. 4) then
          call ESMF_ArrayGetData(array, r4d4, rc=status)
          VERIFY_(STATUS)
          if (associated(r4d4)) then
             deallocate(r4d4, stat=status)
             VERIFY_(STATUS)

             wrap4dr4%ptr => r4d4
             call c_ESMC_ArraySetF90Ptr(array, wrap4dr4, status) 
             VERIFY_(STATUS)
          end if

       else
          RETURN_(ESMF_FAILURE)
       endif
    else if (dk .eq. ESMF_R8) then
!ALT: temporaty set to FAIL; if compiles OK, copy and paste from above; replace r4=>r8
       RETURN_(ESMF_FAILURE)
    else
       RETURN_(ESMF_FAILURE)
    endif

    RETURN_(ESMF_SUCCESS)


  end subroutine MAPL_ArrayF90Deallocate

  subroutine MAPL_DecomposeDim ( dim_world,dim,NDEs )
      implicit   none
      integer    dim_world, NDEs
      integer    dim(0:NDEs-1)
      integer    n,im,rm,nbeg,nend
      im = dim_world/NDEs
      rm = dim_world-NDEs*im
      do n=0,NDEs-1
                      dim(n) = im
      if( n.le.rm-1 ) dim(n) = im+1
      enddo
  end subroutine MAPL_DecomposeDim

  subroutine MAPL_Interp_Fac (TIME0, TIME1, TIME2, FAC1, FAC2, RC)

!------------------------------------------------------------        

!  PURPOSE:
!  ========
!
!    Compute interpolation factors, fac, to be used 
!    in the calculation of the instantaneous boundary 
!    conditions, ie:
!
!     q(i,j) = fac1*q1(i,j) + (1.-fac1)*q2(i,j)
!
!    where:
!     q(i,j)  => Boundary Data valid    at time0
!     q1(i,j) => Boundary Data centered at time1
!     q2(i,j) => Boundary Data centered at time2

!  INPUT:
!  ======
!    time0    : Time of current timestep
!    time1    : Time of boundary data 1 
!    time2    : Time of boundary data 2 

!  OUTPUT:
!  =======
!     fac1    : Interpolation factor for Boundary Data 1
!
! ------------------------------------------------------------        
!               GODDARD LABORATORY FOR ATMOSPHERES            
! ------------------------------------------------------------        

    type(ESMF_Time),   intent(in ) :: TIME0, TIME1, TIME2
    real,              intent(out) :: FAC1
    real,    optional, intent(out) :: FAC2
    integer, optional, intent(out) :: RC
 
    type(ESMF_TimeInterval)        :: TimeDif1
    type(ESMF_TimeInterval)        :: TimeDif
 
    TimeDif1 = TIME2-TIME0
    TimeDif  = TIME2-TIME1
       
    FAC1 = TimeDif1/TimeDif

    if(present(FAC2)) FAC2 = 1.-FAC1
    if(present(RC  )) RC   = ESMF_SUCCESS
 
  end subroutine MAPL_Interp_Fac

  subroutine MAPL_ClimInterpFac (CLOCK,I1,I2,FAC, RC)

!------------------------------------------------------------        

    type(ESMF_CLOCK),  intent(in ) :: CLOCK
    integer,           intent(OUT) :: I1, I2
    real,              intent(out) :: FAC
    integer, optional, intent(out) :: RC
 
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR), parameter   :: IAm='MAPL_ClimInterpFac'

    type (ESMF_Time)                  :: CurrTime
    type (ESMF_Time)                  :: midMonth
    type (ESMF_Time)                  :: BEFORE, AFTER
    type (ESMF_TimeInterval)          :: oneMonth
    type (ESMF_Calendar)              :: cal

    call ESMF_ClockGet       ( CLOCK,    CurrTime=CurrTime, calendar=cal, rc=STATUS )
    VERIFY_(STATUS)
    call ESMF_TimeGet        ( CurrTime, midMonth=midMonth,               rc=STATUS )
    VERIFY_(STATUS)
    call ESMF_TimeIntervalSet( oneMonth, MM = 1, calendar=cal,            rc=status )
    VERIFY_(STATUS)

    if( CURRTIME < midMonth ) then
       AFTER    = midMonth
       midMonth = midMonth - oneMonth
       call ESMF_TimeGet (midMonth, midMonth=BEFORE, rc=STATUS )
       VERIFY_(STATUS)
    else
       BEFORE   = midMonth
       midMonth = midMonth + oneMonth
       call ESMF_TimeGet (midMonth, midMonth=AFTER , rc=STATUS )
       VERIFY_(STATUS)
    endif

    call MAPL_Interp_Fac( CURRTIME, BEFORE, AFTER, FAC, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_TimeGet (BEFORE, MM=I1, rc=STATUS )
    VERIFY_(STATUS)
    call ESMF_TimeGet (AFTER , MM=I2, rc=STATUS )
    VERIFY_(STATUS)
 

    RETURN_(ESMF_SUCCESS)
  end subroutine MAPL_ClimInterpFac


subroutine MAPL_TimeStringGet(TIMESTRING,YY,MM,DD,H,M,S)
  character(len=*),  intent (IN ) :: TIMESTRING
  integer, optional, intent (OUT) :: YY
  integer, optional, intent (OUT) :: MM
  integer, optional, intent (OUT) :: DD
  integer, optional, intent (OUT) :: H
  integer, optional, intent (OUT) :: M
  integer, optional, intent (OUT) :: S

  integer :: IYY, IMM, IDD, IHH, IMN, ISS

  read(TIMESTRING,'(I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2)') IYY,IMM,IDD,IHH,IMN,ISS
  
!ALT: SGI compiler does not like this format  read(TIMESTRING,'(I4,"-",I2,"-",I2,"T",I2,":",I2,":",I2)') IYY,IMM,IDD,IHH,IMN,ISS
  if(present(YY)) YY = IYY
  if(present(MM)) MM = IMM
  if(present(DD)) DD = IDD
  if(present(H )) H  = IHH
  if(present(M )) M  = IMN
  if(present(S )) S  = ISS

  return
end subroutine MAPL_TimeStringGet


subroutine MAPL_UnpackTime(TIME,IYY,IMM,IDD)
  integer, intent (IN ) :: TIME
  integer, intent (OUT) :: IYY
  integer, intent (OUT) :: IMM
  integer, intent (OUT) :: IDD
  IYY = TIME/10000
  IMM = mod(TIME/100,100)
  IDD = mod(TIME,100)
end subroutine MAPL_UnpackTime

subroutine MAPL_PackTime(TIME,IYY,IMM,IDD)
  integer, intent (OUT) :: TIME
  integer, intent (IN ) :: IYY
  integer, intent (IN ) :: IMM
  integer, intent (IN ) :: IDD
  TIME=IYY*10000+IMM*100+IDD
end subroutine MAPL_PackTime

subroutine MAPL_tick (nymd,nhms,ndt)
      integer nymd,nhms,ndt,nsec,nsecf
      nsecf(nhms) = nhms/10000*3600 + mod(nhms,10000)/100*60 + mod(nhms,100)
      IF(NDT.NE.0) THEN
      NSEC = NSECF(NHMS) + NDT
      IF (NSEC.GT.86400)  THEN
      DO WHILE (NSEC.GT.86400)
      NSEC = NSEC - 86400
      NYMD = MAPL_INCYMD (NYMD,1)
      ENDDO
      ENDIF   
      IF (NSEC.EQ.86400)  THEN
      NSEC = 0
      NYMD = MAPL_INCYMD (NYMD,1)
      ENDIF   
      IF (NSEC.LT.00000)  THEN
      DO WHILE (NSEC.LT.0)
      NSEC = 86400 + NSEC
      NYMD = MAPL_INCYMD (NYMD,-1)
      ENDDO
      ENDIF   
      NHMS = MAPL_NHMSF (NSEC)
      ENDIF   
      RETURN  
end subroutine MAPL_tick    

logical function MAPL_RTRN(A,iam,line,rc)
   integer,           intent(IN ) :: A
   character*(*),     intent(IN ) :: iam
   integer,           intent(IN ) :: line
   integer, optional, intent(OUT) :: RC

     MAPL_RTRN = .true.
     if(A/=ESMF_SUCCESS)print'(A40,I10)',Iam,line
     if(present(RC)) RC=A
end function MAPL_RTRN

logical function MAPL_VRFY(A,iam,line,rc)
   integer,           intent(IN ) :: A
   character*(*),     intent(IN ) :: iam
   integer,           intent(IN ) :: line
   integer, optional, intent(OUT) :: RC
     MAPL_VRFY = A/=ESMF_SUCCESS 
     if(MAPL_VRFY)then
       if(present(RC)) then
         print'(A40,I10)',Iam,line
         RC=A
       endif
     endif
end function MAPL_VRFY

logical function MAPL_ASRT(A,iam,line,rc)
   logical,           intent(IN ) :: A
   character*(*),     intent(IN ) :: iam
   integer,           intent(IN ) :: line
   integer, optional, intent(OUT) :: RC
     MAPL_ASRT = .not.A 
     if(MAPL_ASRT)then
       if(present(RC))then
         print'(A40,I10)',Iam,LINE
         RC=ESMF_FAILURE
       endif
     endif
end function MAPL_ASRT

logical function MAPL_RTRNt(A,text,iam,line,rc)
   integer,           intent(IN ) :: A
   character*(*),     intent(IN ) :: text,iam
   integer,           intent(IN ) :: line
   integer, optional, intent(OUT) :: RC

     MAPL_RTRNt = .true.
     if(A/=ESMF_SUCCESS)then
        print'(A40,I10)',Iam,line
        print *, text
     end if
     if(present(RC)) RC=A

end function MAPL_RTRNT

logical function MAPL_VRFYt(A,text,iam,line,rc)
   integer,           intent(IN ) :: A
   character*(*),     intent(IN ) :: iam,text
   integer,           intent(IN ) :: line
   integer, optional, intent(OUT) :: RC
     MAPL_VRFYt =  MAPL_VRFY(A,iam,line,rc)
     if(MAPL_VRFYt) print *, text
end function MAPL_VRFYT

logical function MAPL_ASRTt(A,text,iam,line,rc)
   logical,           intent(IN ) :: A
   character*(*),     intent(IN ) :: iam,text
   integer,           intent(IN ) :: line
   integer, optional, intent(OUT) :: RC
     MAPL_ASRTt =   MAPL_ASRT(A,iam,line,rc)
     if(MAPL_ASRTt) print *, text
end function MAPL_ASRTT

integer function MAPL_nsecf2 (nhhmmss,nmmdd,nymd)
      integer nhhmmss,nmmdd,nymd,nhms,nday,month
      integer nsday, ncycle,iday,iday2
      integer nsecf,i,nsegm,nsegd
      PARAMETER ( NSDAY  = 86400 )
      PARAMETER ( NCYCLE = 1461*24*3600 )
      INTEGER YEAR, DAY, SEC, YEAR0, DAY0, SEC0
      integer    MNDY(12,4), mnd48(48)
      DATA MND48/0,31,60,91,121,152,182,213,244,274,305,335,366,397,34*0 /
!     DATA MNDY /0,31,60,91,121,152,182,213,244,274,305,335,366,397,34*0 /
      equivalence ( mndy(1,1), mnd48(1) )
      nsecf(nhms) = nhms/10000*3600 + mod(nhms,10000)/100*60 + mod(nhms,100)
      MAPL_nsecf2 = nsecf( nhhmmss )
      if( nmmdd.eq.0 ) return
      DO 100 I=15,48
!     MNDY(I,1) = MNDY(I-12,1) + 365
      MND48(I) = MND48(I-12) + 365
100   CONTINUE
      nsegm =     nmmdd/100
      nsegd = mod(nmmdd,100)
      YEAR   = NYMD / 10000
      MONTH  = MOD(NYMD,10000) / 100
      DAY    = MOD(NYMD,100)
      SEC    = NSECF(nhhmmss)
      IDAY   = MNDY( MONTH ,MOD(YEAR ,4)+1 )
      month = month + nsegm
      If( month.gt.12 ) then
      month = month - 12
      year = year + 1
      endif
      IDAY2  = MNDY( MONTH ,MOD(YEAR ,4)+1 )
                    nday = iday2-iday
      if(nday.lt.0) nday = nday + 1461
                    nday = nday + nsegd
      MAPL_nsecf2 = MAPL_nsecf2 + nday*nsday
end function MAPL_nsecf2

integer function MAPL_nhmsf (nsec)
        implicit none
        integer  nsec
        MAPL_nhmsf =  nsec/3600*10000 + mod(nsec,3600)/60*100 + mod(nsec,60)
end function MAPL_nhmsf

integer function MAPL_incymd (NYMD,M)                                                  
      integer nymd,ny,nm,nd,m,ny00
      INTEGER NDPM(12)                                                          
      DATA    NDPM /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/             
      LOGICAL LEAP                                                              
      DATA    NY00     / 1900 /                                                 
      LEAP(NY) = MOD(NY,4).EQ.0 .AND. (NY.NE.0 .OR. MOD(NY00,400).EQ.0)         
      NY = NYMD / 10000                                                         
      NM = MOD(NYMD,10000) / 100                                                
      ND = MOD(NYMD,100) + M                                                    
      IF (ND.EQ.0) THEN                                                         
      NM = NM - 1                                                               
      IF (NM.EQ.0) THEN                                                         
          NM = 12                                                               
          NY = NY - 1                                                           
      ENDIF                                                                     
      ND = NDPM(NM)                                                             
      IF (NM.EQ.2 .AND. LEAP(NY))  ND = 29                                      
      ENDIF                                                                     
      IF (ND.EQ.29 .AND. NM.EQ.2 .AND. LEAP(NY))  GO TO 20                      
      IF (ND.GT.NDPM(NM)) THEN                                                  
      ND = 1                                                                    
      NM = NM + 1                                                               
      IF (NM.GT.12) THEN                                                        
          NM = 1                                                                
          NY = NY + 1                                                           
      ENDIF                                                                     
      ENDIF                                                                     
   20 CONTINUE                                                                  
      MAPL_INCYMD = NY*10000 + NM*100 + ND                                           
      RETURN                                                                    
end function MAPL_incymd


subroutine MAPL_PICKEM(II,JJ,IM,JM,COUNT)
integer, intent(IN ) :: IM, JM, COUNT
integer, intent(OUT) :: II(COUNT), JJ(COUNT)

integer, parameter :: NT=3

logical :: MASK(IM,JM)
integer :: L, NN, IX, JX
real    :: IIR(NT*COUNT), JJR(NT*COUNT)

   MASK=.true.

   NN=1

   call RANDOM_NUMBER(IIR)
   call RANDOM_NUMBER(JJR)


   do L=1, COUNT

      do
         IX=IIR(NN)*(IM-1)+2
         JX=JJR(NN)*(JM-2)+2

         NN = NN + 1

         if(MASK(IX,JX)) then
            II(L) = IX
            JJ(L) = JX
            MASK(IX-1:IX+1,JX-1:JX+1) = .false.
            exit
         endif

         if(NN>NT*COUNT) stop 222

      enddo
   enddo

!!$   DO L=1,JM
!!$      PRINT '(144L1)',MASK(:,L) 
!!$   ENDDO
!!$
!!$   PRINT *, COUNT, NN

   return
 end subroutine MAPL_PICKEM




    subroutine MAPL_GetFieldTimeFromField ( FIELD, TIME, RC )
      type(ESMF_Field),        intent(IN   ) :: FIELD
      type(ESMF_Time),         intent(  OUT) :: TIME
      integer, optional,       intent(  OUT) :: RC

      character(len=ESMF_MAXSTR),parameter   :: IAm=" MAPL_GetFieldTimeFromField"
      integer                                :: STATUS

      integer                                :: YEAR, MONTH, DAY
      integer                                :: HOUR, MINUTE, SCND
      character(len=ESMF_MAXSTR)             :: TIMESTAMP

      call ESMF_FieldGetAttribute(FIELD,     NAME="TimeStamp", VALUE=TIMESTAMP, RC=STATUS)
      if(STATUS/=0) then
         call ESMF_TimeSet          (TIME,      YY=0,                RC=STATUS)
      else
         call MAPL_TimeStringGet    (TIMESTAMP, YY=YEAR, MM=MONTH,  DD=DAY,   & 
                                                H =HOUR, M =MINUTE, S =SCND   )
         VERIFY_(STATUS)
         call ESMF_TimeSet          (TIME,      YY=YEAR, MM=MONTH,  DD=DAY,   &
                                                H =HOUR, M =MINUTE, S =SCND,  &
                                                                     RC=STATUS)
         VERIFY_(STATUS)
      end if

      RETURN_(ESMF_SUCCESS)
    end subroutine MAPL_GetFieldTimeFromField

! ------------------------------------------------------------------------------

    subroutine  MAPL_SetFieldTimeFromField (FIELD, TIME, RC )
      type(ESMF_FIELD),        intent(INOUT) :: FIELD
      type(ESMF_TIME),         intent(IN   ) :: TIME
      integer, optional,       intent(  OUT) :: RC

      character(len=ESMF_MAXSTR),parameter   :: IAm=" MAPL_SetFieldTimeFromField"
      integer                                :: STATUS

      integer                                :: YEAR, MONTH, DAY
      character(len=ESMF_MAXSTR)             :: TIMESTAMP

      call ESMF_TimeGet          (TIME,  timeString=TIMESTAMP,             RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldSetAttribute(FIELD, NAME="TimeStamp", VALUE=TIMESTAMP, RC=STATUS)
      VERIFY_(STATUS)

      RETURN_(ESMF_SUCCESS)
    end subroutine  MAPL_SetFieldTimeFromField


    subroutine  MAPL_GetFieldTimeFromState ( STATE, Fieldname, TIME, RC )
      type(ESMF_STATE),        intent(IN   ) :: STATE
      character(len=*),        intent(IN   ) :: Fieldname
      type(ESMF_Time),         intent(  OUT) :: TIME
      integer, optional,       intent(  OUT) :: RC

      character(len=ESMF_MAXSTR),parameter   :: IAm=" MAPL_GetFieldTimeFromState"
      integer                                :: STATUS

      type(ESMF_FIELD)                       :: FIELD
      integer                                :: YEAR, MONTH, DAY
      character(len=ESMF_MAXSTR)             :: TIMESTAMP

      call ESMF_StateGetField (STATE, FIELDNAME, FIELD, RC=STATUS )
      VERIFY_(STATUS)
      call MAPL_FieldGetTime  (FIELD, TIME,             RC=STATUS)
      VERIFY_(STATUS)

      RETURN_(ESMF_SUCCESS)
    end subroutine  MAPL_GetFieldTimeFromState

! ------------------------------------------------------------------------------

    subroutine  MAPL_SetFieldTimeFromState ( STATE, Fieldname, TIME, RC )
      type(ESMF_STATE),        intent(INOUT) :: STATE
      character(len=*),        intent(IN   ) :: Fieldname
      type(ESMF_Time),         intent(IN   ) :: TIME
      integer, optional,       intent(  OUT) :: RC

      character(len=ESMF_MAXSTR),parameter   :: IAm=" MAPL_SetFieldTimeFromState"
      integer                                :: STATUS

      type(ESMF_FIELD)                       :: FIELD
      integer                                :: YEAR, MONTH, DAY
      character(len=ESMF_MAXSTR)             :: TIMESTAMP

      call ESMF_StateGetField (STATE, FIELDNAME, FIELD, RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_FieldSetTime  (FIELD, TIME,             RC=STATUS)
      VERIFY_(STATUS)

      RETURN_(ESMF_SUCCESS)
    end subroutine  MAPL_SetFieldTimeFromState


    function MAPL_FieldCreateRename(FIELD, NAME, RC) RESULT(F)
      type (ESMF_Field), intent(IN   ) :: FIELD
      character(len=*),  intent(IN   ) :: NAME
      integer, optional, intent(  OUT) :: RC
      type (ESMF_Field)                :: F
      type (ESMF_DataType)             :: dt
      type (ESMF_DataKind)             :: dk

!   we are creating new field so that we can change the name of the field;
!   the important thing is that the data (ESMF_Array) and the grid (ESMF_Grid) 
!   are the SAME as the one in the original Field

      type(ESMF_FieldDataMap) :: datamap           
      type(ESMF_RelLoc)       :: relloc
      type(ESMF_Grid)         :: grid
      type(ESMF_Array)        :: array
      character(len=ESMF_MAXSTR)       :: attname
      integer                 :: status
      character(len=ESMF_MAXSTR), parameter :: Iam='MAPL_FieldCreateRename'

!ALT added kludge (next 6 lines)
      call ESMF_FieldGet(FIELD, name=attname, RC=STATUS)
      VERIFY_(STATUS)
      if (NAME == attname) then
         F = FIELD
         RETURN_(ESMF_SUCCESS)
      endif
      call ESMF_FieldGet(FIELD, grid=GRID, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldGet(FIELD, datamap=DATAMAP, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldGet(FIELD,  horzRelLoc=RELLOC, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldGetArray(FIELD, Array, RC=STATUS)
      VERIFY_(STATUS)
      
      F = ESMF_FieldCreate(GRID, ARRAY,          &
           datamap = datamap,                        &
           horzRelloc = relloc,                      &
           name  = NAME,            RC=STATUS )
      VERIFY_(STATUS)

      call MAPL_FieldCopyAttributes(FIELD_IN=field, FIELD_OUT=f, RC=status)
      VERIFY_(STATUS)

      RETURN_(ESMF_SUCCESS)
    end function MAPL_FieldCreateRename

    function MAPL_FieldCreateRegrid(FIELD, GRID, RC) RESULT(F)
      type (ESMF_Field), intent(IN   ) :: FIELD
      type (ESMF_Grid),  intent(IN   ) :: GRID
      integer, optional, intent(  OUT) :: RC
      type (ESMF_Field)                :: F
      type (ESMF_DataType)             :: dt
      type (ESMF_DataKind)             :: dk

!   we are creating new field so that we can change the grid of the field 
!   (and allocate array accordingly);

      type(ESMF_FieldDataMap) :: datamap           
      type(ESMF_RelLoc)       :: hrelloc
      type(ESMF_RelLoc)       :: vrelloc
      type(ESMF_Array)        :: array
      integer                 :: rank
      integer                 :: COUNTS(3)
      real, pointer           :: VAR_1D(:), VAR_2D(:,:), VAR_3D(:,:,:)
      character(len=ESMF_MAXSTR) :: NAME
      integer                 :: status
      integer                 :: DIMS
      character(len=ESMF_MAXSTR), parameter :: Iam='MAPL_FieldCreateRegrid'

      call ESMF_FieldGet(FIELD, name=name, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldGet(FIELD, datamap=DATAMAP, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldGet(FIELD,  horzRelLoc=HRELLOC, vertRelLoc=VRELLOC, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_FieldGetArray(FIELD, Array, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_GridGetDELocalInfo(GRID, &
           horzRelLoc=ESMF_CELL_CENTER, &
           vertRelLoc=ESMF_CELL_CELL, &
           localCellCountPerDim=COUNTS,RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet(array, rank=rank, rc=status)
      VERIFY_(STATUS)

      if (rank == 1) then
         rank = 2
         call ESMF_FieldDataMapSetDefault(datamap, rank, rc=status)
         VERIFY_(STATUS)
      end if

      if (rank == 2) then
!ALT halowidth assumed 0
         allocate(VAR_2D(COUNTS(1), COUNTS(2)), STAT=STATUS)
         VERIFY_(STATUS)
         ARRAY = ESMF_ArrayCreate(VAR_2D, ESMF_DATA_REF, RC=STATUS)
         VERIFY_(STATUS)
         DIMS = MAPL_DimsHorzOnly
      else
!ALT halowidth assumed 0
         allocate(VAR_3D(COUNTS(1), COUNTS(2), COUNTS(3)), STAT=STATUS)
         VERIFY_(STATUS)
         ARRAY = ESMF_ArrayCreate(VAR_3D, ESMF_DATA_REF, RC=STATUS)
         VERIFY_(STATUS)
         DIMS = MAPL_DimsHorzVert
      end if

      F = ESMF_FieldCreate(GRID, ARRAY,          &
           datamap = datamap,                        &
           horzRelloc = hrelloc,                     &
           vertRelloc = vrelloc,                     &
           name  = NAME,            RC=STATUS )
      VERIFY_(STATUS)

      call MAPL_FieldCopyAttributes(FIELD_IN=field, FIELD_OUT=f, RC=status)
      VERIFY_(STATUS)
! Overwrite DIMS attribute
      call ESMF_FieldSetAttribute(F, NAME='DIMS', VALUE=DIMS, RC=STATUS)
      VERIFY_(STATUS)

      RETURN_(ESMF_SUCCESS)
    end function MAPL_FieldCreateRegrid

    subroutine MAPL_FieldCopyAttributes(FIELD_IN, FIELD_OUT, RC)
      type (ESMF_Field), intent(IN   ) :: FIELD_IN
      type (ESMF_Field), intent(INOUT) :: FIELD_OUT
      integer, optional, intent(  OUT) :: RC

      type (ESMF_DataType)             :: dt
      type (ESMF_DataKind)             :: dk
      integer                          :: status
      character(len=ESMF_MAXSTR), parameter :: Iam='MAPL_FieldCopyAttributes'
      integer                          :: i, n, count
      character(len=ESMF_MAXSTR)       :: attname
      character(len=ESMF_MAXSTR)       :: att
      integer, pointer                 :: iptr(:)
      type(ESMF_Logical), pointer      :: lptr(:)
      real,    pointer                 :: rptr(:)

      call ESMF_FieldGetAttributeCount(field_in, count=n, rc=status)
      VERIFY_(STATUS)

      do i = 1, n
         call  ESMF_FieldGetAttributeInfo(field_in, attributeIndex=i, name=attname, &
                                          datatype=dt, datakind=dk, count=count, rc=status)
         VERIFY_(STATUS)

         if (dt == ESMF_Data_Integer) then
            allocate(iptr(count), stat=status)
            VERIFY_(STATUS)
            call ESMF_FieldGetAttribute(field_in,  NAME=attname, count=count, VALUELIST=iptr, RC=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldSetAttribute(field_out, NAME=attname, count=count, VALUELIST=iptr, RC=STATUS)
            VERIFY_(STATUS)
            deallocate(iptr)

         else if (dt == ESMF_Data_Logical) then
            allocate(lptr(count), stat=status)
            VERIFY_(STATUS)
            call ESMF_FieldGetAttribute(field_in,  NAME=attname, count=count, VALUELIST=lptr, RC=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldSetAttribute(field_out, NAME=attname, count=count, VALUELIST=lptr, RC=STATUS)
            VERIFY_(STATUS)
            deallocate(lptr)

         else if (dt == ESMF_Data_Real) then
            allocate(rptr(count), stat=status)
            VERIFY_(STATUS)
            call ESMF_FieldGetAttribute(field_in,  NAME=attname, count=count, VALUELIST=rptr, RC=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldSetAttribute(field_out, NAME=attname, count=count, VALUELIST=rptr, RC=STATUS)
            VERIFY_(STATUS)
            deallocate(rptr)

         else if (dt == ESMF_Data_Character) then
            call ESMF_FieldGetAttribute(field_in,  NAME=attname, VALUE=att, RC=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldSetAttribute(field_out, NAME=attname, VALUE=att, RC=STATUS)
            VERIFY_(STATUS)

         end if
      end do
      RETURN_(ESMF_SUCCESS)
    end subroutine MAPL_FieldCopyAttributes

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function MAPL_RemapBounds_3dr4(A,I1,IM,J1,JM,L1,LM)
        integer,      intent(IN) :: I1,IM,J1,JM,L1,LM
        real, target, intent(IN) :: A(I1:IM,J1:JM,L1:LM)
        real, pointer            :: MAPL_RemapBounds_3dr4(:,:,:)

        MAPL_RemapBounds_3dr4 => A
      end function MAPL_RemapBounds_3dr4

end module MAPL_BaseMod
