! $Id: FVAdvStateMod.F90,v 1.2 2011-08-09 22:12:59 mrdamon Exp $
!
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: FVAdvStateMod.F90 -- Manages the FVAdv component states
!
! !DESCRIPTION: 
!    This module contains functionality to allocate, define and deallocate
!    a test data set for the tracer advection component.  It takes as
!    arguments the import and export states of the tracer advection 
!    component.
!
!    \paragraph{Routines:}
!
!    \begin{tabular}{||l|l||}    \hline \hline
!      Routine                &  Description \\ \hline
!      FVAdvStateCreate   &  Allocate Tracer Advection states \\
!      FVAdvStateDestroy  &  Deallocate Tracer Advection states \\
!      FVAdvStateData     &  Fill Tracer Advection states with test data\\
!      \hline \hline
!    \end{tabular}
!      
!EOP
!------------------------------------------------------------------------------
!
!

#define ASSERT_(cond_) IF(.NOT. (cond_)) PRINT *, "assert failed in file: ", __FILE__, " at line: ",  __LINE__
  module FVAdvStateMod

    ! ESMF module, defines all ESMF data types and procedures
    use ESMF_Mod
    use FVadvcore_GridCompMod,    only :  T_ADV_STATE, T_ADV_STATE_WRAP

    implicit none
    private

    type(ESMF_Bundle), save :: tracer_bundle
    type(ESMF_Field), pointer, dimension(:), save :: tracer_field
    type(ESMF_Field), save  :: field_cx, field_cy, field_mfx, field_mfy, &
                               field_delp, field_newdelp
 
    ! States and DELayouts for the Subcomponents
    type(ESMF_DELayout), save :: layoutTop

    ! Public methods
    public FVAdvStateCreate, FVAdvStateDestroy, FVAdvStateData

!------------------------------------------------------------------------------

    contains
        
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOPI
! !IROUTINE: FVAdvStateCreate - allocate FVAdv component im/ex states

! !INTERFACE:
     subroutine FVAdvStateCreate(comp, importState, exportState)

! !USES:
    use FVAdvUtilsMod, only :  FVAdvUtilsInit
    use mod_comm, only : mp_sendirr, mp_recvirr
    implicit none

!
! !ARGUMENTS:
     type(ESMF_GridComp), intent(inout) :: comp
     type(ESMF_State), intent(inout) :: importState
     type(ESMF_State), intent(inout) :: exportState

!
! !DESCRIPTION:
!     This routine allocates an example import and export state
!     for the tracer advection component.
!     
!     The arguments are:
!     \begin{description}
!     \item[comp] 
!          Component.
!     \item[importState]
!          Import state.
!     \item[exportState]
!          Export state.
!     \end{description}
!
!EOPI
    integer, parameter            :: ntotq = 1    ! One tracer

    type (T_ADV_STATE), pointer :: state
    type (T_ADV_STATE_WRAP)     :: wrap

    integer :: rc

    integer :: npets, lpet
    integer :: ifirstxy, ilastxy, jfirstxy, jlastxy
    integer :: iq

    type(ESMF_Grid)            :: grid
    type(ESMF_Bundle)          :: bundle
    type(ESMF_Field)           :: field
    type(ESMF_VM)              :: vm
    type(ESMF_FieldDataMap)    :: datamap
    type(ESMF_ArraySpec)       :: arrayspec

    integer                    :: globalCounts(3)

    call ESMF_UserCompGetInternalState( comp, 'FVadvcoreState', wrap, rc )
    state => wrap%state
    ASSERT_(rc==0)

    call ESMF_GridCompGet(comp, vm=vm, grid=grid, rc=rc)
    ASSERT_(rc==0)
    call ESMF_VMGet(vm, LocalPet=lpet, rc=rc)
    ASSERT_(rc==0)

!
! Fetch the grid information
!
    call ESMF_GridGet( grid, horzRelloc=ESMF_CELL_CENTER,                   &
                       Vertrelloc=ESMF_CELL_CELL,                           &
                       globalCellCountPerDim=globalCounts )

!
! For 3D fields such as delta pressure
!
    call ESMF_ArraySpecSet(arrayspec, rank=3, type=ESMF_DATA_REAL, &
                           kind=ESMF_R8, rc=rc )
    ASSERT_( rc == 0 )

    tracer_bundle = ESMF_BundleCreate(NAME="Tracers", grid=grid, rc=rc)
    ASSERT_( rc == 0 )

    call ESMF_FieldDataMapSetDefault(datamap,3,(/1,2,0/),                  &
                                     counts=(/globalCounts(3)/),rc=rc)
    ASSERT_( rc == 0 )

!
! Add Tracer bundle to both import and export state
!
    allocate( tracer_field(ntotq) )
    do iq=1, ntotq
      tracer_field(iq) = ESMF_FieldCreate(grid, arrayspec, datamap=datamap,&
                                          horzRelloc=ESMF_CELL_CENTER,     &
                                          haloWidth=0, name="constituent", &
                                          rc=rc)
      ASSERT_( rc == 0 )
      call ESMF_BundleAddField(tracer_bundle, tracer_field(iq), rc=rc)
      ASSERT_( rc == 0 )
    enddo 
    call ESMF_StateAddBundle(importState, tracer_bundle, rc)
    ASSERT_( rc == 0 )
    call ESMF_StateAddBundle(exportState, tracer_bundle, rc)
    ASSERT_( rc == 0 )

!
! Add CX to import state
!
    field_cx = ESMF_FieldCreate(grid, arrayspec, datamap=datamap,         &
                                horzRelloc=ESMF_CELL_CENTER,              &
                                vertRelloc=ESMF_CELL_CELL,                &   ! Necessary?
                                name="CourantX", rc=rc)
    ASSERT_( rc == 0 )
    call ESMF_StateAddField(importState, field_cx, rc=rc)
    ASSERT_( rc == 0 )

!
! Add CY to import state
!
    field_cy = ESMF_FieldCreate(grid, arrayspec, datamap=datamap,         &
                                horzRelloc=ESMF_CELL_CENTER,              &
                                haloWidth=0, name="CourantY", rc=rc)
    ASSERT_( rc == 0 )
    call ESMF_StateAddField(importState, field_cy, rc)
    ASSERT_( rc == 0 )

!
! Add MFX to import state
!
    field_mfx= ESMF_FieldCreate(grid, arrayspec,  datamap=datamap,        &
                                horzRelloc=ESMF_CELL_CENTER,              &
                                haloWidth=0, name="MassFluxX", rc=rc)
    ASSERT_( rc == 0 )
    call ESMF_StateAddField(importState, field_mfx, rc)
    ASSERT_( rc == 0 )

!
! Add MFY to import state
!
    field_mfy= ESMF_FieldCreate(grid, arrayspec,  datamap=datamap,        &
                                horzRelloc=ESMF_CELL_CENTER,              &
                                haloWidth=0, name="MassFluxY", rc=rc)
    ASSERT_( rc == 0 )
    call ESMF_StateAddField(importState, field_mfy, rc)
    ASSERT_( rc == 0 )

!
! Add DELP to import state
!
    field_delp= ESMF_FieldCreate(grid, arrayspec,  datamap=datamap,       &
                                horzRelloc=ESMF_CELL_CENTER,              &
                                haloWidth=0, name="DeltaPressure", rc=rc)
    ASSERT_( rc == 0 )
    call ESMF_StateAddField(importState, field_delp, rc)
    ASSERT_( rc == 0 )

!
! Add DELP to import state
!
    field_newdelp= ESMF_FieldCreate(grid, arrayspec,  datamap=datamap,       &
                                   horzRelloc=ESMF_CELL_CENTER,              &
                                   haloWidth=0, name="NewDeltaPressure", rc=rc)
    ASSERT_( rc == 0 )
    call ESMF_StateAddField(importState, field_newdelp, rc)
    ASSERT_( rc == 0 )

    end subroutine FVAdvStateCreate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOPI
! !IROUTINE: FVAdvStateDestroy - destroy FVAdv component im/ex states
!
! !INTERFACE:
     subroutine FVAdvStateDestroy(comp, importState, exportState)

! !USES:
    use FVAdvUtilsMod, only :  FVAdvUtilsFinal
    implicit none

!
! !ARGUMENTS:
     type(ESMF_GridComp), intent(inout) :: comp
     type(ESMF_State), intent(inout) :: importState
     type(ESMF_State), intent(inout) :: exportState

!
! !DESCRIPTION:
!     This routine allocates an example import and export state
!     for the tracer advection component.
!     
!     The arguments are:
!     \begin{description}
!     \item[comp] 
!          Component.
!     \item[importState]
!          Import state.
!     \item[exportState]
!          Export state.
!     \end{description}
!
!EOPI
!------------------------------------------------------------------------------
    integer :: rc

! TODO:   call FVAdvUtilsFinal, and other deallocation routine

    end subroutine FVAdvStateDestroy


!------------------------------------------------------------------------------
!BOPI
! !IROUTINE: FVAdvStateData - defines FVAdv state with test data 

! !INTERFACE:
     subroutine FVAdvStateData(comp, importState, exportState)

! !USES:
    use FVAdvUtilsMod, only :  FVAdvUtilsInit, set_eta, &
                                   air_mass_flux, pudding, RemapBounds
    use TransposeMod, only  :  T_TRANSPOSE
    use parutilitiesmodule, only:  parpatterntype, parpatterncreate
    use mod_comm, only : mp_sendirr, mp_recvirr
    implicit none

!
! !ARGUMENTS:
     type(ESMF_GridComp), intent(inout) :: comp
     type(ESMF_State), intent(inout) :: importState
     type(ESMF_State), intent(inout) :: exportState  
!
! !DESCRIPTION:
!     This routine creates 'artificial' initial data for a test run.
!     In this case, the data correspond to Williamson Test 1 for 
!     simple advection.  
!     
!     The arguments are:
!     \begin{description}
!     \item[comp] 
!          Component.
!     \item[importState]
!          Import state.
!     \item[exportState]
!          Export state. (Not actually set; here only for symmetry
!     \item[rc] k
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors,
!          otherwise {\tt ESMF\_FAILURE}.
!     \end{description}
!
!EOPI
    type (T_ADV_STATE), pointer :: state
    type (T_ADV_STATE_WRAP)     :: wrap
    type (T_TRANSPOSE), pointer   :: transposes

    real(ESMF_KIND_R8), parameter :: HALF = 0.5
    real(ESMF_KIND_R8), parameter :: PI = 3.141592654
    real(ESMF_KIND_R8), parameter :: SECONDS_IN_DAY = 86400.0
    real(ESMF_KIND_R8), parameter :: EARTH_RADIUS = 6370000.0
    real(ESMF_KIND_R8), parameter :: ANGLE = HALF*PI - 0.05

    integer :: rc
    integer :: npets
    integer :: i, j, k, len
    integer :: dt, iord, jord, kord, ntotq, nq, lpet
    integer :: im, jm, km, ng
    integer :: ifirstxy, ilastxy, jfirstxy, jlastxy
    integer :: ks
    integer, parameter :: n_adj = 0   ! Number of adjustments to air_mass_flux
                                      ! 0 = none

    type(ESMF_Grid)   :: grid
    type(ESMF_Bundle) :: bundle
    type(ESMF_Field)  :: field
    type(ESMF_Array)  :: array
    type(ESMF_LocalArray) :: larray

    type (ParPatternType)           :: g2l3d

    real(ESMF_KIND_R8), allocatable :: clat(:)
    real(ESMF_KIND_R8), allocatable :: sine(:), cosp(:), sinp(:), cose(:)
    real(ESMF_KIND_R8), allocatable :: g_u(:,:,:), g_v(:,:,:)
    real(ESMF_KIND_R8), allocatable :: g_psg(:,:,:), g_ps(:,:)
    real(ESMF_KIND_R8), allocatable :: g_fx(:,:,:), g_fy(:,:,:)
    real(ESMF_KIND_R8), allocatable :: g_cx(:,:,:), g_cy(:,:,:)
    real(ESMF_KIND_R8), allocatable :: g_va(:,:,:)
    real(ESMF_KIND_R8), allocatable :: g_delp(:,:,:), g_delp1(:,:,:)
    real(ESMF_KIND_R8)              :: lat, lon, dp, dl, ph5, pint, ptop
    real(ESMF_KIND_R8)              :: velocityScale, deltaT
    real(ESMF_KIND_R8)              :: del_ak, del_bk

    real(ESMF_KIND_R8), allocatable :: ak(:), bk(:)
    real(ESMF_KIND_R8), allocatable :: ps(:,:)

    real(ESMF_KIND_R8), pointer     :: tmp3d(:,:,:)
    real(ESMF_KIND_R8), pointer     :: content(:,:,:)

    logical, allocatable            :: g_ffsl(:,:)

    type(ESMF_VM) :: vm

    velocityScale = 2 * PI * EARTH_RADIUS / (12 * SECONDS_IN_DAY)

    call ESMF_UserCompGetInternalState( comp, 'FVadvcoreState', wrap, rc )
    state => wrap%state
    ASSERT_(rc==0)
    transposes => state%transposes             ! pointer for convenience

    call ESMF_GridCompGet(comp, vm=vm, rc=rc)
    ASSERT_(rc==0)
    call ESMF_VMGet(vm, LocalPet=lpet, rc=rc)
    ASSERT_(rc==0)

    im = state%im
    jm = state%jm
    km = state%km

    dt = state%dt
    deltaT = REAL(dt, ESMF_KIND_R8)

    allocate( clat(jm) )

    dl  = (PI+PI)/im
    dp  = PI / (jm-1)

    lat = -HALF*PI
    do j = 1, jm
      clat(j) = lat
      lat     = lat + dp
    enddo

!
! Initialize the advection utilities
!
    call FVAdvUtilsInit( im, jm, km, deltaT, EARTH_RADIUS, clat )

!
! Retrieve the pointers to AK and BK and fill them with replicated data
!

    allocate( ak(km+1) )
    allocate( bk(km+1) )

    call set_eta(km, ks, ptop, pint, ak, bk)

!
! Update the pressure at top of atmosphere in state
!
    state%ptop = ptop

    if ( lpet == 0 ) then

!
! This code is being executed on the root PET because air_mass_flux
! is not available for an XY decomposition.  Since this is part of
! the initialization, it should not be an issue.
!
      allocate( g_u(im,jm,km), g_v(im,jm,km) )
      allocate( g_psg(im,jm,2), g_ps(im,jm) )
      allocate( g_fx(im,jm,km), g_fy(im,jm,km) )
      allocate( g_cx(im,jm,km), g_cy(im,jm,km) )
      allocate( g_va(im,jm,km), g_ffsl(jm,km) )
      allocate( g_delp(im,jm,km), g_delp1(im,jm,km) )

      do k = 1, km
        lat = -HALF*PI
        do j = 1, jm
          lon = -PI
          do i = 1, im
            lon = lon + dl
            g_u(i,j,k) =  velocityScale * ( COS(lat)*COS(ANGLE) +          &
                                          SIN(lat)*COS(lon)*SIN(ANGLE) )
            g_v(i,j,k) = -velocityScale * SIN(lon)*SIN(ANGLE)
          enddo
          lat = lat + dp
        enddo
      enddo

      iord = state%iord
      jord = state%jord

!  AK and BK read in from file, PSG is constant (no change in pressure)
!  U, V defined above.
    

      g_psg(:,:,:) = 100000.0    ! Standard surface pressure   (Pascals?)

      call air_mass_flux(im, jm, km, iord, jord,                         &
                         ak, bk, g_psg, g_ps, g_u, g_v,                  &
                         g_cx, g_cy, g_va, g_fx, g_fy,                   &
                         g_ffsl, g_delp1,  g_delp,                       &
                         deltaT, EARTH_RADIUS, n_adj)
    endif

    ifirstxy = state%ifirst
    ilastxy  = state%ilast 
    jfirstxy = state%jfirst
    jlastxy  = state%jlast 

    call ESMF_StateGetBundle(importState, "Tracers", bundle, rc=rc)
    ASSERT_( rc == 0 )
    call ESMF_BundleGetField( bundle, 1, field, rc=rc )
    ASSERT_( rc == 0 )
    call ESMF_FieldGetDataPointer(field, tmp3d, ESMF_DATA_REF, rc=rc)
    ASSERT_( rc == 0 )
!
! The bounds have to be remapped in order to use global indices subsequently
!
    content => RemapBounds(tmp3d, ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km)

!
! Define the constituent 'pudding' around LONG=1.5_r8*PI, LAT=0
!
    do k = 1, km
      lat = -HALF*PI + (jfirstxy-1)*dp
      do j = jfirstxy, jlastxy
        lon = -PI + (ifirstxy-1)*dl
        do i = ifirstxy, ilastxy
          content(i,j,k) = pudding(lon,lat)
          lon = lon + dl
        enddo
        lat = lat + dp
      enddo
    enddo

    allocate( ps( ifirstxy:ilastxy, jfirstxy:jlastxy ) )

    call mp_sendirr( g_ps, transposes%g2l_xy%SendDesc,&
                     transposes%g2l_xy%RecvDesc, ps )
    call mp_recvirr( ps, transposes%g2l_xy%RecvDesc )

!
! Create a 3D global-to-local scatter pattern from the 2D
!
    call parpatterncreate( transposes%g2l_xy, g2l3d, km )

! Get pointers to target arrays and scatter the global arrays

    call ESMF_StateGetField(importState, "CourantX", field, rc=rc)
    ASSERT_( rc == 0 )
    call ESMF_FieldGetDataPointer(field, tmp3d, ESMF_DATA_REF, rc=rc)
    ASSERT_( rc == 0 )
    call mp_sendirr( g_cx, g2l3d%SendDesc,&
                     g2l3d%RecvDesc, tmp3d )
    call mp_recvirr( tmp3d, g2l3d%RecvDesc )


    call ESMF_StateGetField(importState, "CourantY", field, rc=rc)
    ASSERT_( rc == 0 )
    call ESMF_FieldGetDataPointer(field, tmp3d, ESMF_DATA_REF, rc=rc)
    ASSERT_( rc == 0 )
    call mp_sendirr( g_cy, g2l3d%SendDesc,&
                     g2l3d%RecvDesc, tmp3d )
    call mp_recvirr( tmp3d, g2l3d%RecvDesc )

    call ESMF_StateGetField(importState, "MassFluxX", field, rc=rc)
    ASSERT_( rc == 0 )
    call ESMF_FieldGetDataPointer(field, tmp3d, ESMF_DATA_REF, rc=rc)
    ASSERT_( rc == 0 )
    call mp_sendirr( g_fx, g2l3d%SendDesc,&
                     g2l3d%RecvDesc, tmp3d )
    call mp_recvirr( tmp3d, g2l3d%RecvDesc )

    call ESMF_StateGetField(importState, "MassFluxY", field, rc=rc)
    ASSERT_( rc == 0 )
    call ESMF_FieldGetDataPointer(field, tmp3d, ESMF_DATA_REF, rc=rc)
    ASSERT_( rc == 0 )
    call mp_sendirr( g_fy, g2l3d%SendDesc,&
                     g2l3d%RecvDesc, tmp3d )
    call mp_recvirr( tmp3d, g2l3d%RecvDesc )

    call ESMF_StateGetField(importState, "DeltaPressure", field, rc=rc)
    ASSERT_( rc == 0 )
    call ESMF_FieldGetDataPointer(field, tmp3d, ESMF_DATA_REF, rc=rc)
    ASSERT_( rc == 0 )
    call mp_sendirr( g_delp, g2l3d%SendDesc,&
                     g2l3d%RecvDesc, tmp3d )
    call mp_recvirr( tmp3d, g2l3d%RecvDesc )


    call ESMF_StateGetField(importState, "NewDeltaPressure", field, rc=rc)
    ASSERT_( rc == 0 )
    call ESMF_FieldGetDataPointer(field, tmp3d, ESMF_DATA_REF, rc=rc)
    ASSERT_( rc == 0 )

    call setTargetDelP( ifirstxy, ilastxy, jfirstxy, jlastxy, km,  &
                        ak,       bk,      ps,       tmp3d )

    if ( lpet == 0 ) then
      deallocate( g_u, g_v )
      deallocate( g_psg, g_ps )
      deallocate( g_fx, g_fy )
      deallocate( g_cx, g_cy )
      deallocate( g_delp, g_delp1 )
      deallocate( g_ffsl, g_va )
    endif

    deallocate( ps   )
    deallocate( ak   )
    deallocate( bk   )
    deallocate( clat )

    print *, "Finished Tracer Advection Component Initial Data"

    contains

      subroutine setTargetDelP( ifirstxy, ilastxy, jfirstxy, jlastxy, km, &
                                ak,       bk,      ps,       delP  )
      
      integer, intent(in)  :: ifirstxy, ilastxy, jfirstxy, jlastxy, km
      real(ESMF_KIND_R8), intent(in) :: ak(km+1)
      real(ESMF_KIND_R8), intent(in) :: bk(km+1)
      real(ESMF_KIND_R8), intent(in) :: ps(ifirstxy:ilastxy, jfirstxy:jlastxy)
      real(ESMF_KIND_R8), intent(out):: delP(ifirstxy:ilastxy,jfirstxy:jlastxy,km)

      integer             :: i, j, k
      real(ESMF_KIND_R8)  :: del_ak, del_bk

      do k=1,km
        del_ak = ak(k+1) - ak(k)
        del_bk = bk(k+1) - bk(k)
        do j=jfirstxy,jlastxy
          do i=ifirstxy,ilastxy
            delp(i,j,k) = del_ak + del_bk*ps(i,j)
          enddo
        enddo
      enddo

      end subroutine setTargetDelP

    end subroutine FVAdvStateData
!------------------------------------------------------------------------------

  end module FVAdvStateMod
