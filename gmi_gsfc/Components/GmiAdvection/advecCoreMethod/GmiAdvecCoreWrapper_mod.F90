#include "GmiESMF_ErrLog.h"
!--------------------------------------------------------------------------------
! NASA/GSFC - SIVO Code 610.3
!--------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiAdvecCoreWrapper_mod
!
      module GmiAdvecCoreWrapper_mod
!
! !USES:
      use ESMF_Mod
      use GmiESMF_ErrorChecking_mod
      use GmiArrayBundlePointer_mod    , only : t_GmiArrayBundle
      use GmiSpcConcentrationMethod_mod, only : t_SpeciesConcentration, &
     &               Get_concentration
      use GmiAdvectionMethod_mod       , only : t_Advection, &
     &       Get_numAdvectedSpecies, Get_advectedSpecies, Get_advectedSpeciesMap

!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: initializeAdvecCoreTracers
      public  :: fromGmiToAdvecCore
      public  :: fromAdvecCoreToGmi

#     include "GmiParameters.h"

! !DESCRIPTION:
!
! !AUTHOR:
! Jules.Kouatchou-1@nasa.gov
!EOP
!--------------------------------------------------------------------------------
      CONTAINS
!--------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initializeAdvecCoreTracers
!
! !INTERFACE:
!
      subroutine initializeAdvecCoreTracers (state, grid, Advection)
!
      implicit none
!
! !INPUT PARAMETERS:
      type(ESMF_Grid)  , intent(in) :: grid
      type(t_Advection), intent(in) :: Advection
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
! Created  in the bundle ESMF fields corresponding to advected tracers.
!
! !LOCAL VARIABLES:
!      integer :: STATUS, numVars, ib, rc, rank
!      integer :: counts(ESMF_MAXGRIDDIM)
!      real(ESMF_KIND_R4), pointer :: tempArray(:,:,:)
!      type(ESMF_Bundle)           :: tracerBun
!      type(ESMF_Field)            :: field
!      type (ESMF_FieldDataMap)    :: DATAMAP
!      integer :: numAdvectedSpecies
!      character(len=MAX_STRING_LENGTH), pointer :: advectedSpecies(:)
      character(len=ESMF_MAXSTR) :: IAm = "initializeAdvecCoreTracers"
!
!EOP
!--------------------------------------------------------------------------------
!BOC

!      ! Get species info from GMI Advection object
!
!      call Get_numAdvectedSpecies(Advection, numAdvectedSpecies)
!
!      allocate(advectedSpecies(numAdvectedSpecies))
!
!      call Get_advectedSpecies   (Advection, advectedSpecies)
!
!      ! GET ARRAY DIMENSIONS
!
!      call ESMF_GridGetDELocalInfo(grid, &
!     &          horzRelLoc=ESMF_CELL_CENTER, &
!     &          vertRelloc=ESMF_CELL_CELL, &
!     &          localCellCountPerDim=counts, RC=STATUS)
!      VERIFY_(STATUS)

!      ! Get the bundle from the state
!
!      call ESMF_StateGetBundle(state, "TRACERS", tracerBun, rc=STATUS)
!      VERIFY_(STATUS)
!
!      call ESMF_FieldDataMapSetDefault(datamap, 3, rc=status)
!      VERIFY_(STATUS)
!
!      SPECIES: do ib = 1, numAdvectedSpecies
!         allocate(tempArray(counts(1), counts(2), counts(3)))
!
!         field = ESMF_FieldCreate(grid, tempArray, &
!     &                        horzRelloc=ESMF_CELL_CENTER,     &
!     &                        datamap=datamap, &
!     &                        name=trim(advectedSpecies(ib)), &
!     &                        rc=STATUS)
!         VERIFY_(STATUS)
!
!         call ESMF_BundleAddField(tracerBun, field, rc=STATUS)
!         VERIFY_(STATUS)
!      end do SPECIES
!
!! Sanity check
!
!      call ESMF_BundleGet(tracerBun, fieldCount=numVars , rc=STATUS)
!      VERIFY_(STATUS)
!      ASSERT_(numAdvectedSpecies == numVars)
!
!! cleanup
!
!      deallocate(advectedSpecies)
!
      return

      end  subroutine initializeAdvecCoreTracers
!EOC
!--------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fromGmiToAdvecCore
!
! !INTERFACE:
!
      subroutine fromGmiToAdvecCore(state, grid, SpeciesConcentration, &
     &      Advection, UC, VC, MX_UR, MY_UR, MX, MY, MZ, DPEDT, DP, &
     &      i1, i2, ju1, j2, k1, k2)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer                      , intent(in) :: i1, i2, ju1, j2, k1, k2
      type(t_Advection)            , intent(in) :: Advection
      type(ESMF_Grid)              , intent(in) :: grid
      real(ESMF_KIND_R4), target   , intent(in) :: UC   (i1:i2,ju1:j2,k1:k2)
      real(ESMF_KIND_R4), target   , intent(in) :: VC   (i1:i2,ju1:j2,k1:k2)
      real(ESMF_KIND_R4), target   , intent(in) :: MX_UR(i1:i2,ju1:j2,k1:k2)
      real(ESMF_KIND_R4), target   , intent(in) :: MY_UR(i1:i2,ju1:j2,k1:k2)
      real(ESMF_KIND_R4), target   , intent(in) :: MX   (i1:i2,ju1:j2,k1:k2)
      real(ESMF_KIND_R4), target   , intent(in) :: MY   (i1:i2,ju1:j2,k1:k2)
      real(ESMF_KIND_R4), target   , intent(in) :: MZ   (i1:i2,ju1:j2,k1:k2)
      real(ESMF_KIND_R4), target   , intent(in) :: DPEDT(i1:i2,ju1:j2,k1:k2)
      real(ESMF_KIND_R4), target   , intent(in) :: DP   (i1:i2,ju1:j2,k1:k2)
      type (t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
!      integer :: STATUS, numVars, ib, rc
!      real(ESMF_KIND_R4), pointer, dimension(:,:,:) :: ptr3Edg
!      real(ESMF_KIND_R4), pointer, dimension(:,:,:) :: ptr3Dreal
!      type(ESMF_Field )                :: Field 
!      type(ESMF_Array )                :: Array
!      type(ESMF_Bundle)                :: tracerBun
!      type (t_GmiArrayBundle), pointer :: concentration(:)
!      integer                               :: numAdvectedSpecies
!      integer, pointer                      :: advectedSpeciesMap(:)
!      character(len=MAX_STRING_LENGTH), pointer :: advectedSpecies(:)
      character(len=ESMF_MAXSTR) :: IAm = "fromGmiToAdvecCore"
!
!EOP
!--------------------------------------------------------------------------------
!BOC
!      call Get_numAdvectedSpecies(Advection, numAdvectedSpecies)
!      
!      allocate(advectedSpecies(numAdvectedSpecies))
!      call Get_advectedSpecies(Advection, advectedSpecies)
!
!      allocate(advectedSpeciesMap(numAdvectedSpecies))
!      call Get_advectedSpeciesMap(Advection, advectedSpeciesMap)
!
!      ! Get the species concentration
!      call Get_concentration(SpeciesConcentration, concentration)
!
!      ! Get the bundle from the state
!
!      call ESMF_StateGetBundle(state, "TRACERS", tracerBun, rc=STATUS)
!      VERIFY_(STATUS)
!
!      call ESMF_BundleGet(tracerBun, fieldCount=numVars, rc=STATUS)
!      VERIFY_(STATUS)
!
!      ! Verify that the number of fields in the bundle is equal to the number
!      ! of advected species.
!
!      ASSERT_(numVars == numAdvectedSpecies)
!
!      do ib = 1, numVars
!         call ESMF_BundleGetDataPointer(tracerBun,advectedSpecies(ib), &
!                  ptr3Dreal, rc=status )
!
!         ptr3Dreal(:,:,ubound(ptr3Dreal,3):lbound(ptr3Dreal,3):-1) = &
!     &       concentration(advectedSpeciesMap(ib))%pArray3D(:,:,:)
!      end do
!
!      call copyData2StateField(state, 'UC', UC)
!      call copyData2StateField(state, 'VC', VC)
!      call copyData2StateField(state, 'MX_UR', MX_UR)
!      call copyData2StateField(state, 'MY_UR', MY_UR)
!      call copyData2StateField(state, 'MX', MX)
!      call copyData2StateField(state, 'MY', MY)
!      call copyData2StateField(state, 'MZ', MZ)
!      call copyData2StateField(state, 'DP', DP)
!      call copyData2StateField(state, 'DPEDT', DPEDT)
!
!      deallocate(advectedSpecies)
!      deallocate(advectedSpeciesMap)

      return

      end subroutine fromGmiToAdvecCore
!EOC
!--------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fromAdvecCoreToGmi
!
! !INTERFACE:
!
      subroutine fromAdvecCoreToGmi (state, grid, SpeciesConcentration, &
     &      Advection)
!
      implicit none
!
! !INPUT PARAMETERS:
      type(ESMF_Grid)  , intent(in) :: grid
      type(t_Advection), intent(in) :: Advection
      type(ESMF_State) , intent(in) :: state
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_SpeciesConcentration), intent(in) :: SpeciesConcentration
!
! !DESCRIPTION:
! Extract the species concentration from the advecCore state.
!
! !LOCAL VARIABLES:

!      integer :: ib, numVars, rc
!      type(ESMF_Bundle)                :: tracerBun
!      type (t_GmiArrayBundle), pointer :: concentration(:)
!      real(ESMF_KIND_R4),      pointer :: ptr3D(:,:,:)
!      integer                          :: numAdvectedSpecies
!      integer, pointer                 :: advectedSpeciesMap(:)
!      character(len=MAX_STRING_LENGTH), pointer :: advectedSpecies(:)
!
!      integer :: STATUS
      character(len=ESMF_MAXSTR) :: IAm = "fromAdvecCoreToGmi"
!
!EOP
!--------------------------------------------------------------------------------
!BOC  
!      call Get_numAdvectedSpecies(Advection, numAdvectedSpecies)
!      
!      allocate(advectedSpecies(numAdvectedSpecies))
!      call Get_advectedSpecies(Advection, advectedSpecies)
!
!      allocate(advectedSpeciesMap(numAdvectedSpecies))
!      call Get_advectedSpeciesMap(Advection, advectedSpeciesMap)
!
!      ! Get the species concentration
!      call Get_concentration(SpeciesConcentration, concentration)
!      
!      call ESMF_StateGetBundle (state, 'TRACERS', tracerBun, RC=STATUS )
!      VERIFY_(STATUS) 
!      
!      do ib=1,numAdvectedSpecies
!         call ESMF_BundleGetDataPointer(tracerBun,advectedSpecies(ib), &
!                  ptr3D, rc=status)
!
!         ASSERT_(size(concentration(advectedSpeciesMap(ib))%pArray3D,3)==size(ptr3D,3))
!         concentration(advectedSpeciesMap(ib))%pArray3D(:,:,:) = &
!     &           ptr3D(:,:,ubound(ptr3D,3):lbound(ptr3D,3):-1)
!      end do
!
!      deallocate(advectedSpecies)
!      deallocate(advectedSpeciesMap)
!
      return
      
      end subroutine fromAdvecCoreToGmi
!EOC  
!--------------------------------------------------------------------------------
!-------------------------------------------------------------------
!      subroutine copyData2StateField(state, name, ptr)
!!-------------------------------------------------------------------
!
!      implicit none
!
!      real(ESMF_KIND_R4),target :: ptr(:,:,:)
!      character(len=*) :: name
!      type(ESMF_State) :: state
!! 
!      type(ESMF_Field)  :: FIELD
!      real(ESMF_KIND_R4), pointer :: ptr3d(:,:,:)
!      integer :: STATUS, rc
!      character(len=ESMF_MAXSTR)       :: IAm
!      
!! start
!      IAm = "copyData2StateField"
!
!      ptr3d => ptr       
!      
!      call ESMF_StateGetField(state, name, field, rc=STATUS)
!      VERIFY_(STATUS)
!
!      call ESMF_FieldSetDataPointer(FIELD, ptr3d, rc=STATUS)
!      VERIFY_(STATUS)
!
!       return
!
!       end subroutine copyData2StateField 
!
!!-------------------------------------------------------------------
!   subroutine add2D2Bundle(ogrid,ptr,name,Bun)
!!-------------------------------------------------------------------
!
!   implicit none
!
!   type(ESMF_Grid)   :: ogrid
!   real,target :: ptr(:,:)
!   character(len=*) :: name
!   type(ESMF_Bundle) :: Bun
!! 
!   type(ESMF_Array)  :: EARRAY
!   type(ESMF_Field)  :: EFIELD
!   type(ESMF_FieldDataMap) :: DATAMAP           
!   real, pointer :: ptr2d(:,:)
!   integer :: STATUS, rc
!      character(len=ESMF_MAXSTR)       :: IAm
!      
!! start
!   IAm = "add2D2Bundle"
!
!   ptr2d => ptr       
!      
!   call ESMF_FieldDataMapSetDefault(datamap,2, rc=status)
!   VERIFY_(STATUS)
!   EARRAY = ESMF_ArrayCreate(ptr2d, ESMF_DATA_REF, RC=STATUS)
!   VERIFY_(STATUS) 
!   EFIELD = ESMF_FieldCreate(ogrid, EARRAY,            &
!        datamap    = datamap,                    &
!        horzRelloc = ESMF_CELL_CENTER,                     &
!        name       = name,             RC=STATUS )
!   VERIFY_(STATUS)
!   call ESMF_BundleAddField (Bun, EFIELD, RC=STATUS )
!   VERIFY_(STATUS)
!      
!   end subroutine add2D2Bundle 
!
!!-------------------------------------------------------------------
!   subroutine add3D2Bundle(ogrid,ptr,name,Bun)
!!-------------------------------------------------------------------
!   implicit none
!   type(ESMF_Grid)   :: ogrid
!   real(ESMF_KIND_R8), target :: ptr(:,:,:)
!   character(len=*) :: name
!   type(ESMF_Bundle) :: Bun
!!
!   type(ESMF_Array)  :: EARRAY
!   type(ESMF_Field)  :: EFIELD
!   type(ESMF_FieldDataMap) :: DATAMAP
!   real(ESMF_KIND_R8), pointer :: ptr3d(:,:,:)
!   integer :: STATUS, rc
!      character(len=ESMF_MAXSTR)       :: IAm
!
!! start
!   IAm = "add3D2Bundle"
!
!   ptr3d => ptr
!
!   call ESMF_FieldDataMapSetDefault(datamap,3, rc=status)
!   VERIFY_(STATUS)
!   EARRAY = ESMF_ArrayCreate(ptr3d, ESMF_DATA_REF, RC=STATUS)
!   VERIFY_(STATUS)
!   EFIELD = ESMF_FieldCreate(ogrid, EARRAY,            &
!        datamap    = datamap,                    &
!        horzRelloc = ESMF_CELL_CENTER,                     &
!        vertRelloc = ESMF_CELL_CELL,        &
!        name       = name,             RC=STATUS )
!   VERIFY_(STATUS)
!   call ESMF_BundleAddField (Bun, EFIELD, RC=STATUS )
!   VERIFY_(STATUS)
!
!   end subroutine add3D2Bundle
!!-------------------------------------------------------------------
!!-------------------------------------------------------------------
!   subroutine add3D2State(ogrid,ptr,name,state)
!!-------------------------------------------------------------------
!   implicit none
!   type(ESMF_Grid)   :: ogrid
!   real(ESMF_KIND_R8), target :: ptr(:,:,:)
!   character(len=*) :: name
!   type(ESMF_State) :: state
!!
!   type(ESMF_Array)  :: EARRAY
!   type(ESMF_Field)  :: EFIELD
!   type(ESMF_FieldDataMap) :: DATAMAP
!   real(ESMF_KIND_R8), pointer :: ptr3d(:,:,:)
!   integer :: STATUS, rc
!      character(len=ESMF_MAXSTR)       :: IAm
!
!! start
!   IAm = "add3D2State"
!
!   ptr3d => ptr
!
!   call ESMF_FieldDataMapSetDefault(datamap,3, rc=status)
!   VERIFY_(STATUS)
!   EARRAY = ESMF_ArrayCreate(ptr3d, ESMF_DATA_REF, RC=STATUS)
!   VERIFY_(STATUS)
!   EFIELD = ESMF_FieldCreate(ogrid, EARRAY,            &
!        datamap    = datamap,                    &
!        horzRelloc = ESMF_CELL_CENTER,                     &
!        vertRelloc = ESMF_CELL_CELL,        &
!        name       = name,             RC=STATUS )
!   VERIFY_(STATUS)
!   call ESMF_StateAddField (state, EFIELD, RC=STATUS )
!   VERIFY_(STATUS)
!
!   end subroutine add3D2State
!!-------------------------------------------------------------------
!
end Module GmiAdvecCoreWrapper_mod
