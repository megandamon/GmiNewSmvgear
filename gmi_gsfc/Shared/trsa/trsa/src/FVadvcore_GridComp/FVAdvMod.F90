! $Id: FVAdvMod.F90,v 1.2 2011-08-09 22:13:00 mrdamon Exp $
!
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: FVAdvMod - Container for FVAdv_State type and its base methods
!
! !DESCRIPTION: 
!
!------------------------------------------------------------------------------
!
!
! !INTERFACE:


  module FVAdvMod

! !USES:

    use shr_kind_mod,       only : FVAdvReal => shr_kind_r8
    use mapz_module,        only : mapn_ppm_tracer
    use parutilitiesmodule, only : parinit
    use TracerMod,          only : T_Tracers, T_Tracer_Grid, Tracer_GridInit, &
                                   Tracer_GridGet, Tracer_GridFinal
    use TransposeMod


    implicit none
    private

! !PUBLIC TYPES:

    public T_FVadv_STATE

! !PUBLIC MEMBER FUNCTIONS:

    public FVadv_Final
    public FVadv_Run
    public FVadv_Init
    public FVadv_Get
    public FVadv_Remap
    public FVadv_PrepVars

! !PUBLIC DATA MEMBERS:
    public FVAdvReal

    real(FVAdvReal), parameter :: FVpi   = 3.141592653589793238462643383279_FVAdvReal
    real(FVAdvReal), parameter :: FVae   = 6376000.0_FVAdvReal
    real(FVAdvReal), parameter :: ZERO   = 0.0_FVAdvReal
    real(FVAdvReal), parameter :: ONE    = 1.0_FVAdvReal

!EOP

    type T_FVadv_STATE
      private
      type(T_TRACER_GRID) :: tracerGrid   !  tracer-specific grid information
      type(T_TRANSPOSE)   :: transposes   !  communication patterns
      integer             :: im, jm, km   !  problem dimensions
      integer             :: ifirst, ilast!  segment in X
      integer             :: jfirst, jlast!  segment in Y
      integer             :: kfirst, klast!  segment in Y
      integer, pointer    :: imxy(:)  
      integer, pointer    :: jmxy(:)
      integer             :: lpet
    end type T_FVadv_STATE

contains

!=================================================================

!=================================================================

    subroutine FVadv_Final(state,rc)
      type (T_FVadv_STATE),   intent(INOUT) :: state
      integer, optional,   intent(  out) :: rc

      rc = 0

      call Tracer_GridFinal( state%tracerGrid )
      call TransposeFinal( state%transposes )

      if(associated(state%imxy)) deallocate(state%imxy,stat=rc)
      if(rc /= 0) return
      if(associated(state%jmxy)) deallocate(state%jmxy)
      if(rc /= 0) return

    end subroutine FVadv_Final

!=================================================================

    subroutine FVadv_Init(state, im,jm,km,imxy,jmxy,lpet,comm,rc)

      type (T_FVadv_STATE),   intent(INOUT) :: state
      integer,    intent(in   ) :: im, jm, km
      integer,    intent(in   ) :: imxy(:), jmxy(:)
      integer,    intent(in   ) :: lpet
      integer,  optional,   intent(in   ) :: comm
      integer,  optional,   intent(  out) :: rc

      integer :: ng
      integer :: klow, ind
      integer :: nprxy_x, nprxy_y
      integer :: npryz_y, npryz_z
      integer :: myidxy_x, myidxy_y, myidyz_y, myidyz_z
      integer :: ifirstxy, ilastxy, jfirstxy, jlastxy
      integer :: jfirstyz, jlastyz, kfirstyz, klastyz

      integer, parameter :: JORD = 3

      integer, pointer :: kmyz(:), jmyz(:)

! Begin
!------

      rc = 0

! Initialize Pilgrim
!-------------------

      call parinit(comm)


      state%im   = im  
      state%jm   = jm  
      state%km   = km  
      state%lpet = lpet

      if(associated(state%imxy)) deallocate(state%imxy)
      allocate(state%imxy(size(imxy)),stat=rc)
      if(rc/=0) return
      state%imxy = imxy

      if(associated(state%jmxy)) deallocate(state%jmxy)
      allocate(state%jmxy(size(jmxy)),stat=rc)
      if(rc/=0) return
      state%jmxy = jmxy

! Calculate myid's (possibly use alternative method)
!---------------------------------------------------

      nprxy_x  = size(state%imxy)
      nprxy_y  = size(state%jmxy)
      npryz_y  = nprxy_y
      npryz_z  = nprxy_x

      myidyz_z = state%lpet / npryz_y
      myidxy_y = state%lpet / nprxy_x

      myidyz_y = state%lpet - myidyz_z*npryz_y
      myidxy_x = state%lpet - myidxy_y*nprxy_x

! For now jm is the same in both decompositions

      jmyz => state%jmxy

! Optimal decomposition for Z
!----------------------------

      allocate(kmyz(npryz_z),stat=rc)
      if(rc/=0) return

      klow         = km / npryz_z
      ind          = km - klow*npryz_z
      kmyz(:ind)   = klow+1
      kmyz(ind+1:) = klow

      if( km /= sum(kmyz) ) then
         rc = 1
         return
      endif

! Pick out the starting and finished indices on this DE
!------------------------------------------------------
      
      ifirstxy     = 1 + sum(state%imxy(1:myidxy_x))
      jfirstxy     = 1 + sum(state%jmxy(1:myidxy_y))

      ilastxy      = sum(state%imxy(1:myidxy_x+1))
      jlastxy      = sum(state%jmxy(1:myidxy_y+1))

      jfirstyz     = 1 + sum(jmyz(1:myidyz_y))
      kfirstyz     = 1 + sum(kmyz(1:myidyz_z))

      jlastyz      = sum(jmyz(1:myidyz_y+1))
      klastyz      = sum(kmyz(1:myidyz_z+1))

      state%ifirst = ifirstxy
      state%jfirst = jfirstxy
      state%kfirst = kfirstyz

      state%ilast  = ilastxy
      state%jlast  = jlastxy
      state%klast  = klastyz

      call Tracer_GridInit( FVae, FVpi, JORD,                    &
                            state%im, state%jm, state%km,        &
                            jfirstyz, jlastyz, kfirstyz, klastyz,&
                            myidyz_y, myidyz_z, npryz_y, npryz_z,&
                            state%tracerGrid )    

      call Tracer_GridGet( state%tracerGrid, ng_d = ng )

      call TransposeInit( state%im, state%jm, state%km, ng, 0,   &
                          ifirstxy, ilastxy, jfirstxy, jlastxy,  &
                          jfirstyz, jlastyz, kfirstyz, klastyz,  &
                          nprxy_x, nprxy_y, npryz_y, npryz_z,    &
                          state%imxy, state%jmxy, jmyz, kmyz,    &
                          state%transposes )
! All done
!---------

      deallocate( kmyz )

      rc = 0
  end subroutine FVadv_Init
!----------------------------------------------------------------


!----------------------------------------------------------------
!BOPI
! !IROUTINE: FVAdv_get
! !INTERFACE:
  subroutine FVAdv_get(state, tracerGrid, transposes, im, jm, km, &
                       ifirst, ilast, jfirst, jlast,              &
                       kfirst, klast, imxy, jmxy, lpet )

! !PARAMETERS:
    type(T_FVAdv_STATE), optional, intent(in   ) :: state
    type(T_TRACER_GRID), optional, intent(  out) :: tracerGrid   !  tracer-specific grid information
    type(T_TRANSPOSE),   optional, intent(  out) :: transposes   !  communication patterns
    integer,             optional, intent(  out) :: im, jm, km   !  problem dimensions
    integer,             optional, intent(  out) :: ifirst, ilast!  segment in X
    integer,             optional, intent(  out) :: jfirst, jlast!  segment in Y
    integer,             optional, intent(  out) :: kfirst, klast!  segment in Z
    integer,             optional, intent(  out) :: lpet
    integer,             optional, pointer       :: imxy(:)  
    integer,             optional, pointer       :: jmxy(:)

! !DESCRIPTION:
!    Perform the remapping on the tracers
!
! !REVISION HISTORY:
!    2007.06.01   Sawyer    Creation
!    2007.12.10   Sawyer    Bug fix:  CX, CY not staggered
!EOPI


    if ( present( tracerGrid ) ) tracerGrid = state%tracerGrid
    if ( present( transposes ) ) transposes = state%transposes
    if ( present( im ) )         im         = state%im
    if ( present( jm ) )         jm         = state%jm
    if ( present( km ) )         km         = state%km
    if ( present( ifirst ) )     ifirst     = state%ifirst
    if ( present( ilast  ) )     ilast      = state%ilast 
    if ( present( jfirst ) )     jfirst     = state%jfirst
    if ( present( jlast  ) )     jlast      = state%jlast 
    if ( present( kfirst ) )     kfirst     = state%kfirst
    if ( present( klast  ) )     klast      = state%klast 
    if ( present( lpet  ) )      lpet       = state%lpet
    if ( present( imxy ) )       imxy      => state%imxy
    if ( present( jmxy ) )       jmxy      => state%jmxy


  end subroutine FVAdv_get
!----------------------------------------------------------------


!----------------------------------------------------------------
!BOPI
! !IROUTINE: FVAdv_remap

! !INTERFACE:

  subroutine FVAdv_remap(state, tracer, pe1, pe2, kord)

    type(T_FVAdv_STATE), intent(in   ) :: state
    type (T_TRACERS),    intent(inout) :: tracer(:)
    real(FVAdvReal),     intent(in   ) :: pe1(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km+1)
    real(FVAdvReal),     intent(in   ) :: pe2(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km+1)
    integer,             intent(in   ) :: kord

! !DESCRIPTION:
!    Perform the remapping on the tracers
!
! !REVISION HISTORY:
!    06.08.03   Sawyer    Created from qmap (CTM) -- upgraded to 2D decomp
!    06.08.23   Sawyer    Revisions after 2006.08.22 walkthrough
!EOPI

   integer i,j,k,krd

   if(size(tracer)>0) then
      if(kord == 8) then
         krd = 8
      else
         krd = 7
      endif
      do j=state%jfirst,state%jlast
         call mapn_ppm_tracer ( state%km, pe1(:,j,:),                &
                                tracer, size(tracer),                &
                                state%km, pe2(:,j,:),                &
                                state%ifirst, state%ilast, j,        &
                                state%ifirst, state%ilast,           &
                                state%jfirst, state%jlast, 0, krd    )
      enddo
   endif

 end subroutine FVAdv_Remap
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !ROUTINE: FVAdv_PrepVars --- finalize the Lin tracer adv. scheme
!
! !INTERFACE:

 subroutine FVAdv_PrepVars(state,  dt,    ptop,   grav,                 &
                           uc,     vc,    mx_ur,  my_ur,                &
                           dp,     dpnew, mx,     my,     mz,           &
                           cx,     cy,    mfx,    mfy,    pe2      )

! !USES:
   use parutilitiesmodule, only: bcstop, sumop,  parcollective
   use mod_comm, only : mp_send3d, mp_recv3d

! !PARAMETERS:
    type(T_FVAdv_STATE), intent(in) :: state        ! Internal state
    real(FVAdvReal),     intent(in) :: dt           ! Long time step
    real, intent(in) :: ptop
    real, intent(in) :: grav
    real, intent(in) :: uc(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km)
    real, intent(in) :: vc(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km)
! mass flux (units: Pa m^2 s-1)
    real, intent(in) :: mx_ur(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km)
    real, intent(in) :: my_ur(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km)
    real, intent(in) :: dp(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km)
    real, intent(in) :: dpnew(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km)
    real, intent(in) :: mx(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km)
    real, intent(in) :: my(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km)
    real, intent(in) :: mz(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km)

! Courant numbers
    real(FVAdvReal), intent(out) :: cx(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km)
    real(FVAdvReal), intent(out) :: cy(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km)
! Normed mass fluxes (units: Pa)
    real(FVAdvReal), intent(out) :: mfx(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km)
    real(FVAdvReal), intent(out) :: mfy(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km)

! Final pressure (after remapping)
    real(FVAdvReal), intent(out) :: pe2(state%ifirst:state%ilast,state%jfirst:state%jlast,state%km+1)

! !DESCRIPTION:
!
!   Create the variables needed internally by the advection core,
!   in particular the Courant numbers CX and CY, the normed mass fluxes
!   MFX and MFY, and the final level pressure array PE2, which is
!   consistent with the given 3D mass fluxes, mx, my, and mz.
!   Note that dpnew is currently not used, but remains in the interface
!   as a possible alternative way to calculate PE2.
!
!   All arrays are unghosted, which means that ghosting is performed local
!   to this routine.
!
! !REVISION HISTORY:
!
!   2007.05.21   Sawyer     Creation
!   2007.12.07   Sawyer     Removed staggered grid calculation of CX, CY
!   2007.12.14   Sawyer     Added calculation of PE2
!EOP
!-----------------------------------------------------------------------
!BOC
!
  integer :: im, jm, km
  integer :: ifirst, ilast
  integer :: jfirst, jlast
  integer :: i, j, k
  integer :: js2g0, jn2g0, jn2m1

! Message passing variables
  integer :: dest, src
  integer :: npr_x, npr_y
  integer :: myid_x, myid_y
  integer :: iam
  integer :: commglobal

  real(FVAdvReal)              :: dlon, dphi, dtdy, dtgrav, rcap
  real(FVAdvReal)              :: sum1, sum2, maxdiff
  real(FVAdvReal), allocatable :: dtdx(:)        ! dt/dx

  real(FVAdvReal), pointer     :: cosp(:)        ! Volume cosine

  real(FVAdvReal), pointer     :: dlat(:)        ! D-grid latitudes

  real(FVAdvReal), allocatable :: mx_west(:,:)   ! Send buffer for MX
  real(FVAdvReal), allocatable :: mx_east(:,:)   ! Recv buffer for MX
  real(FVAdvReal), allocatable :: my_south(:,:)  ! Send buffer for MY
  real(FVAdvReal), allocatable :: my_north(:,:)  ! Recv buffer for MY
  real(FVAdvReal), allocatable :: my_sp(:,:)     ! South pole
  real(FVAdvReal), allocatable :: my_np(:,:)     ! South pole

  real(FVAdvReal) :: mfactor ! temporary value of normalization factor

  iam     = state%lpet

  im      = state%im
  jm      = state%jm
  km      = state%km
  ifirst  = state%ifirst
  ilast   = state%ilast
  jfirst  = state%jfirst
  jlast   = state%jlast

  dtgrav  = dt*grav

  call Tracer_GridGet( state%tracerGrid, commglobal = commglobal,           &
                       dlon = dlon, dphi = dphi, rcap = rcap, cosp = cosp )

  js2g0   = max(jfirst,2)
  jn2g0   = min(jlast,jm-1)
  jn2m1   = min(jlast-1,jm-1)

  allocate(mx_west(jfirst:jlast,km))
  allocate(mx_east(jfirst:jlast,km))
  allocate(my_south(ifirst:ilast,km))
  allocate(my_north(ifirst:ilast,km))
  allocate(my_sp(im,km), my_np(im,km) )

!
! Calculation of grid information (should this be in state or grid?)
! TODO: move all this stuff to tracerGrid

  allocate( dtdx(jm) )

  dtdy = dt / (FVae*dphi)
  do j=2,jm-1
    dtdx(j) = dt / (dlon*FVae*cosp(j))
  enddo

  npr_x  = size(state%imxy)
  npr_y  = size(state%jmxy)
  myid_y = iam / npr_x
  myid_x = iam - myid_y*npr_x


!
! MY must be ghosted on north
!
  if (npr_y > 1) then
!$omp  parallel do private(i,k)
     do k=1,km
        do i=ifirst,ilast
           my_south(i,k) = my(i,jfirst,k)
        enddo
     enddo
  
     call mp_send3d( commglobal,  iam-npr_x, iam+npr_x, im, jm, km,     &
                     ifirst, ilast, jfirst, jfirst, 1, km,              &
                     ifirst, ilast, jfirst, jfirst, 1, km, my_south )
  endif

! Convert horizontal mass fluxes from (Pa m^2 / s) 
! to (Pa per long time step) expected by FV

!$omp parallel do private(i,j,k,mfactor)
  do k = 1,km
    do j = js2g0,jn2g0
      mfactor = dt / ( dphi*FVae*dlon*FVae*cosp(j) )
      do i = ifirst,ilast
        mfx(i,j,k) = mx_ur(i,j,k) * mfactor
        cx(i,j,k) = dtdx(j)*uc(i,j,k)
      enddo
    enddo
    mfactor = dt / ( dphi*FVae*dlon*FVae )
    do j = jfirst,jlast
      do i = ifirst,ilast
        mfy(i,j,k) = my_ur(i,j,k) * mfactor             ! No cos factor!
        cy(i,j,k)  = dtdy*vc(i,j,k)
      enddo
    enddo
  enddo

  if (npr_y > 1) then
      call mp_recv3d( commglobal, iam+npr_x, im, jm, km,                 &
                      ifirst, ilast, jlast+1, jlast+1, 1, km,            &
                      ifirst, ilast, jlast+1, jlast+1, 1, km, my_north  )
  endif

!
! MX must be ghosted on east
!
  if (npr_x > 1) then
! Nontrivial x decomposition
     dest = myid_y*npr_x + MOD(iam+npr_x-1,npr_x)
     src  = myid_y*npr_x + MOD(iam+1,npr_x)
!$omp  parallel do private(j,k)
     do k = 1,km
        do j=jfirst,jlast
           mx_west(j,k) = mx(ifirst,j,k)
        enddo
     enddo
     call mp_send3d( commglobal, dest, src, im, jm, km,                &
                     ifirst, ifirst, jfirst, jlast, 1, km,             &
                     ifirst, ifirst, jfirst, jlast, 1, km, mx_west )
  endif

!
! Zero out undefined (j==1 and j==jm) parts of MX_UR and CX
!
  if (jfirst == 1) then
!$omp parallel do private(i,k)
    do k=1,km
      do i=ifirst,ilast
        mfx(i,1,k) = ZERO
        cx(i,1,k) = ZERO
      enddo
    enddo
  endif
  if (jlast == jm) then
!$omp parallel do private(i,k)
    do k=1,km
      do i=ifirst,ilast
        mfx(i,jm,k) = ZERO
        cx(i,jm,k) = ZERO
      enddo
    enddo
  endif

! Calculate consistent MY fluxes at poles

  my_sp = ZERO
  my_np = ZERO
!$omp  parallel do                  &
!$omp  default(shared)              &
!$omp  private(i,k)
  do k=1,km
     if ( jfirst == 1 ) then  ! SP
        do i=ifirst,ilast
           my_sp(i,k) = my(i,2,k)
        enddo
     endif
     if ( jlast == jm ) then  ! NP
        do i=ifirst,ilast
           my_np(i,k) = my(i,jm,k)
        enddo
     endif
  enddo

  if (npr_x > 1) then
     call mp_recv3d( commglobal, src, im, jm, km,                     &
                     ilast+1, ilast+1, jfirst, jlast, 1, km,          &
                     ilast+1, ilast+1, jfirst, jlast, 1, km, mx_east )
  else
!$omp  parallel do private(j,k)
     do k = 1,km
        do j=jfirst,jlast
           mx_east(j,k) = mx(ilast,j,k)
        enddo
     enddo
  endif

!
! Actually communicator commxy_x should used, but commglobal should work

  if (npr_x > 1) then
     call parcollective(commglobal, sumop, im, km, my_sp)
     call parcollective(commglobal, sumop, im, km, my_np)
  endif

!
! pe2 calculation
!
    pe2(:,:,1) = ptop
    maxdiff = ZERO
    do k=1,km-1
      do j=js2g0,jn2m1
        mfactor = dt / ( dphi*FVae*dlon*FVae*cosp(j) )
        do i=ifirst,ilast-1
           pe2(i,j,k+1) = dp(i,j,k) + mfactor*( -(mx(i+1,j,k)-mx(i,j,k) + &
              my(i,j+1,k)-my(i,j,k))) - dtgrav*(mz(i,j,k+1)-mz(i,j,k))
        enddo
        i=ilast
        pe2(i,j,k+1) = dp(i,j,k) + mfactor*( -(mx_east(j,k)-mx(i,j,k) +   &
            my(i,j+1,k)-my(i,j,k))) - dtgrav*(mz(i,j,k+1)-mz(i,j,k))
      enddo
      if (jlast /= jm) then
         j=jlast
         mfactor = dt / ( dphi*FVae*dlon*FVae*cosp(j) )
         do i=ifirst,ilast-1
            pe2(i,j,k+1) = dp(i,j,k) + mfactor*( -(mx(i+1,j,k)-mx(i,j,k) +&
               my_north(i,k)-my(i,j,k))) - dtgrav*(mz(i,j,k+1)-mz(i,j,k))
         enddo
         i=ilast
         pe2(i,j,k+1) = dp(i,j,k) + mfactor*( -(mx_east(j,k)-mx(i,j,k) +  &
               my_north(i,k)-my(i,j,k))) - dtgrav*(mz(i,j,k+1)-mz(i,j,k))
      endif
! Poles
      mfactor = dt * rcap / ( dphi*FVae*dlon*FVae )
      if ( jfirst == 1  ) then
         sum1 = -(SUM(my_sp(1:im,k)))*mfactor
         j=1
         do i=ifirst,ilast
            pe2(i,j,k+1) = dp(i,j,k) + sum1 - dtgrav*(mz(i,j,k+1)-mz(i,j,k))
         enddo
      endif
      if ( jlast  == jm ) then
         sum2 =  (SUM(my_np(1:im,k)))*mfactor
         j=jm
         do i=ifirst,ilast
            pe2(i,j,k+1) = dp(i,j,k) + sum2 - dtgrav*(mz(i,j,k+1)-mz(i,j,k))
         enddo
      endif

!
! Vertical integration
!
      do j=jfirst,jlast
        do i=ifirst,ilast
!!!           maxdiff = max(maxdiff,abs(pe2(i,j,k+1)-dpnew(i,j,k)))
           pe2(i,j,k+1) = pe2(i,j,k+1) + pe2(i,j,k)
!!!           pe2(i,j,k+1) = dpnew(i,j,k) + pe2(i,j,k)
        enddo
      enddo
    enddo
! Final level (earth surface)
    k = km
    do j=js2g0,jn2m1
      mfactor = dt / ( dphi*FVae*dlon*FVae*cosp(j) )
      do i=ifirst,ilast-1
         pe2(i,j,k+1) = dp(i,j,k) + mfactor*( -(mx(i+1,j,k)-mx(i,j,k) + &
            my(i,j+1,k)-my(i,j,k))) + dtgrav*mz(i,j,k)
      enddo
      i=ilast
      pe2(i,j,k+1) = dp(i,j,k) + mfactor*( -(mx_east(j,k)-mx(i,j,k) +   &
          my(i,j+1,k)-my(i,j,k))) + dtgrav*mz(i,j,k)
    enddo
    if (jlast /= jm) then
       j=jlast
       mfactor = dt / ( dphi*FVae*dlon*FVae*cosp(j) )
       do i=ifirst,ilast-1
          pe2(i,j,k+1) = dp(i,j,k) + mfactor*( -(mx(i+1,j,k)-mx(i,j,k) +&
             my_north(i,k)-my(i,j,k))) + dtgrav*mz(i,j,k)
       enddo
       i=ilast
       pe2(i,j,k+1) = dp(i,j,k) + mfactor*( -(mx_east(j,k)-mx(i,j,k) +  &
             my_north(i,k)-my(i,j,k))) + dtgrav*mz(i,j,k)
    endif
! Poles
    mfactor = dt * rcap / ( dphi*FVae*dlon*FVae )
    if ( jfirst == 1  ) then
       sum1 = -(SUM(my_sp(1:im,k)))*mfactor
       j=1
       do i=ifirst,ilast
          pe2(i,j,k+1) = dp(i,j,k) + sum1 + dtgrav*mz(i,j,k)
       enddo
    endif
    if ( jlast  == jm ) then
       sum2 =  (SUM(my_np(1:im,k)))*mfactor
       j=jm
       do i=ifirst,ilast
          pe2(i,j,k+1) = dp(i,j,k) + sum2 + dtgrav*mz(i,j,k)
       enddo
    endif

!
! Vertical integration
!
    do j=jfirst,jlast
      do i=ifirst,ilast
!           maxdiff = max(maxdiff,abs(pe2(i,j,k+1)-dpnew(i,j,k)))
         pe2(i,j,k+1) = pe2(i,j,k+1) + pe2(i,j,k)
!!!           pe2(i,j,k+1) = dpnew(i,j,k) + pe2(i,j,k)
      enddo
    enddo
!!!    print *, "MAXDIFF", maxdiff

  deallocate(mx_west)
  deallocate(mx_east)
  deallocate(my_south)
  deallocate(my_north)
  deallocate(my_sp, my_np)
  deallocate( dtdx )

  return

!EOC
end subroutine FVAdv_PrepVars
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !ROUTINE: FVAdv_Run --- Run tracer advection for one time step
!
! !INTERFACE:
 subroutine FVAdv_Run( state, tracerxy, delpxy, iord,    jord,  &
                       cxxy,  cyxy,      mfxxy,   mfyxy )
 
! !USES:

   use tp_core
   use fill_module

   use TransposeMod, only: T_TRANSPOSE
   use parutilitiesmodule, only: maxop, parcollective
   use mod_comm, only : mp_sendirr, mp_recvirr,            &
                        mp_sendirr_r4, mp_recvirr_r4,      &
                        mp_send4d_ns, mp_recv4d_ns,        &
                        mp_send4d_ns_r4, mp_recv4d_ns_r4,  &
                        mp_send3d_2, mp_recv3d_2, mp_barrier
   implicit none


   integer,  intent(in) :: iord,  jord
   type(T_FVAdv_State), target  :: state


! Delta pressure field initially, updated field on return
   real(FVAdvReal), intent(INOUT) :: delpxy(:,:,:)   !  Lagrangian change

! Courant Numbers
   real(FVAdvReal), intent(IN)  :: cxxy(:,:,:)
   real(FVAdvReal), intent(IN)  :: cyxy(:,:,:)

! Mass Fluxs
   real(FVAdvReal), intent(IN)  :: mfxxy(:,:,:)
   real(FVAdvReal), intent(IN)  :: mfyxy(:,:,:)

! Tracers, updated on return
   type(T_TRACERS), intent(inout):: tracerxy(:)

! !DESCRIPTION:
!
!  Perform large-time-step tracer transport using accumulated Courant
!  numbers (cx, cy) and the mass fluxes (mfx, mfy) within the Lagrangian
!  layers.  This routine is 100\% parallel in the vertical direction
!  (with SMP).  Merdional Courant number will be further split, if
!  necessary, to ensure stability.  Cy <= 1 away from poles; Cy $\le$
!  1/2 at the latitudes closest to the poles.
!
! !REVISION HISTORY:
!
!   SJL 99.04.13:  Delivery
!   WS  99.05.26:  Added jfirst:jlast concept; im, jm, km as parameters
!                  replaced IMR, JMR, JNP, NL with IM, JM-1, JM and KM
!   WS  99.09.27:  Documentation; indentation; jfirst:jlast 
!   WS  99.09.30:  Ghosting; loop limits; full parallelization; tested
!   SJL 99.10.15:  nsplt migrated to outermost loop to remove bug
!   SJL 99.12.19:  Local 2D arrays trimmed!
!   WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!   WS  00.07.13:  Changed PILGRIM API
!   AAM 00.08.29:  Added kfirst, klast
!   AAM 01.06.27:  Added y communicators
!   SJL 30.07.01:  MPI optimization/simplification
!   WS  02.04.24:  New mod_comm interfaces
!   WS  02.07.04:  Fixed 2D decomposition bug dest/src for mp_send3d
!   WS  03.11.19:  Merged in CAM changes by Mirin
!   WS  03.12.03:  Added GRID as argument, dynamics_vars removed
!   WS  04.08.25:  Simplification of interface with GRID
!   WS  04.10.07:  Removed dependency on spmd_dyn; info now in GRID
!   WS  05.04.04:  Transitioned to type T_TRACERS (supports r4 and FVAdvReal)
!   WS  05.04.09:  Each tracer now ghosted individually (save buffers)
!   WS  05.04.12:  Full support for either r4 or FVAdvReal tracers
!   WS  05.05.25:  Merged CAM and GEOS5, e.g. nsplt(k) opt. in CAM
!   PW  05.10.12:  Changes for Cray X1(E), alternative implementation
!                  of double buffering logic
!   WS  06.07.06:  Integrated into tracer advection module
!   WS  06.09.30:  Moved YZ decomposition inside this routine
!   WS  07.05.29:  Extensive changes based on suggestions by Suarez
!   WS  07.06.07:  Additional changes to support proper DP1/DP2
!
!EOP
!---------------------------------------------------------------------
!BOC


   type (T_TRACER_GRID), pointer :: grid
   type (T_TRANSPOSE),   pointer :: transposes

   real(FVAdvReal)  tiny
   parameter ( tiny = 1.e-10_FVAdvReal )

   type(T_TRACERS), pointer, dimension(:) :: tracer
   real(FVAdvReal), dimension(:,:,:), pointer    :: cx, cy, mfx, mfy, delp
   real(FVAdvReal), dimension(:,:,:), pointer    :: va, flx

! Local variables:
! 2d arrays

   real(FVAdvReal), allocatable  ::  a2(:,:)
   real(FVAdvReal), allocatable  :: fx(:,:)
   real(FVAdvReal), allocatable  :: fy(:,:)
   real(FVAdvReal), allocatable  :: cymax(:)

! 1d arrays
   real(FVAdvReal), dimension(:), pointer    :: cosp, sine, acosp


! Temporary FVAdvReal array for Q
   real(FVAdvReal), allocatable  :: q_r8(:,:,:,:)    ! Contains 2 3D ghosted arrays
   real(FVAdvReal), allocatable  :: dp1(:,:,:)

   logical :: ffsl(state%jm,state%kfirst:state%klast)
   integer :: nsplt(state%kfirst:state%klast)
   integer :: im, jm, km                    ! Dimensions
   integer :: ng                            ! Max number of ghost latitudes
   integer :: jfirst, jlast, kfirst, klast  ! YZ decomposition limits
   integer :: cur, nxt                      ! current and next q_r8 buffer indices

   integer :: i, j, k
   integer :: it, iq, kq, nq, max_nsplt
   integer :: k_courant, kend
   integer :: ktot
   integer :: js1gd, js2g0, js2gd, jn2g0,jn2gd,jn1g1,jn1gd
   integer :: dest, src
   integer :: iam, npryz_y
   integer :: comm_y, comm_z
   integer :: commglobal

   integer, parameter :: fill = .true.


   real(FVAdvReal) :: cy_global
   real(FVAdvReal) :: frac
   real(FVAdvReal) :: cmax
   real(FVAdvReal) :: sum1, sum2
   real(FVAdvReal) :: acap, rcap

   grid => state%TracerGrid
   transposes => state%transposes

   nq = size(tracerxy)

   if ( nq == 0 ) return    ! Nothing to advect

   cur     = 1
   nxt     = 2

   im      = state%im
   jm      = state%jm
   km      = state%km
   iam     = state%lpet

   call Tracer_GridGet( grid, ng_d = ng, acap = acap, rcap = rcap,            &
                        comm_y=comm_y, comm_z=comm_z, commglobal=commglobal,  &
                        jfirstyz = jfirst, jlastyz = jlast,                   &
                        kfirstyz = kfirst, klastyz = klast, npryz_y = npryz_y,&
                        cosp = cosp, sine = sine, acosp = acosp )

   ktot  = klast - kfirst + 1
   js2g0 = max(2,jfirst)
   jn2g0 = min(jm-1,jlast)
   jn1g1 = min(jm,jlast+1)
   js1gd = max(1,jfirst-ng)     ! NG latitudes on S (starting at 1)
   js2gd = max(2,jfirst-ng)     ! NG latitudes on S (starting at 2)
   jn2gd = min(jm-1,jlast+ng)   ! NG latitudes on S (ending at jm-1)
   jn1gd = min(jm,jlast+ng)     ! NG latitudes on N (ending at jm)

!
! Allocate the YZ-decomposed arrays needed by FVAdvRun
!
   allocate(q_r8(im,jfirst-ng:jlast+ng,kfirst:klast,1:2))
   allocate( tracer(nq) )
!
! If not twod_decomp, still allocate the ghosted arrays,		
! at least for the time being
!
   allocate( cx(im,jfirst-ng:jlast+ng,kfirst:klast) )
   allocate( cy(im,jfirst:jlast+1,kfirst:klast) )
   allocate( mfy(im,jfirst:jlast+1,kfirst:klast) )

   allocate( delp(im,jfirst:jlast,kfirst:klast) )
   allocate( mfx(im,jfirst:jlast,kfirst:klast) )
   allocate( dp1(im,jfirst:jlast,kfirst:klast) )

! First order upwinded value for courant number in Y on A-grid
   allocate( va(im,jfirst:jlast,kfirst:klast) )
! Scaled mass flux in X
   allocate( flx(im,jfirst:jlast,kfirst:klast) )

! 2D arrays
   allocate( a2(im,jfirst:jlast) )
   allocate( fx(im,jfirst:jlast) )
   allocate( fy(im,jfirst:jlast+1) )
   allocate( cymax(kfirst:klast) )

!
! Transpose to YZ decomposition
!
     call mp_sendirr( commglobal, cxxy, transposes%ng_xy2yz%SendDesc, &
                      transposes%ng_xy2yz%RecvDesc, cx )
     call mp_recvirr( commglobal, cx, transposes%ng_xy2yz%RecvDesc )

     call mp_sendirr( commglobal, cyxy, transposes%jp1_xy2yz%SendDesc,&
                      transposes%jp1_xy2yz%RecvDesc, cy )
     call mp_recvirr( commglobal, cy, transposes%jp1_xy2yz%RecvDesc )

     call mp_sendirr( commglobal, mfyxy, transposes%jp1_xy2yz%SendDesc,&
                      transposes%jp1_xy2yz%RecvDesc, mfy )
     call mp_recvirr( commglobal, mfy, transposes%jp1_xy2yz%RecvDesc )

!
! Transpose data to YZ decomposition
!
     call mp_sendirr( commglobal, mfxxy, transposes%xy2yz%SendDesc,    &
                      transposes%xy2yz%RecvDesc, mfx )
     call mp_recvirr( commglobal, mfx, transposes%xy2yz%RecvDesc )

     call mp_sendirr( commglobal, delpxy, transposes%xy2yz%SendDesc,   &
                      transposes%xy2yz%RecvDesc, delp )
     call mp_recvirr( commglobal, delp, transposes%xy2yz%RecvDesc )

     do i=1, nq
       if ( tracerxy(i)%is_r4 ) then
         tracer(i)%is_r4 = .TRUE.
         allocate( tracer(i)%content_r4(im,jfirst:jlast, kfirst:klast) )
         call mp_sendirr_r4( commglobal, tracerxy(i)%content_r4,       &
                             transposes%xy2yz%SendDesc,                &
                             transposes%xy2yz%RecvDesc,                &
                             tracer(i)%content_r4 )
         call mp_recvirr_r4( commglobal, tracer(i)%content_r4,         &
                             transposes%xy2yz%RecvDesc )
       else
         tracer(i)%is_r4 = .FALSE.
         allocate( tracer(i)%content(im,jfirst:jlast, kfirst:klast) )
         call mp_sendirr( commglobal, tracerxy(i)%content,             &
                          transposes%xy2yz%SendDesc,                   &
                          transposes%xy2yz%RecvDesc, tracer(i)%content )
         call mp_recvirr( commglobal, tracer(i)%content,               &
                          transposes%xy2yz%RecvDesc )
       endif
     enddo

!
! Perform the halo exchange
!
   call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast,        &
                      kfirst, klast, ng, ng, cx )
   call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast,        &
                      kfirst, klast, ng, ng, cx )

! Send one latitude of both cy and mfy to the south
   dest = iam-1
   src  = iam+1
   if ( mod(iam,npryz_y) == 0 ) dest = -1
   if ( mod(iam+1,npryz_y) == 0 ) src = -1
   call mp_send3d_2( commglobal, dest, src, im, jm, km,            &
                     1, im, jfirst, jlast+1, kfirst, klast,        &
                     1, im, jfirst, jfirst, kfirst, klast, cy, mfy)

!$omp parallel do default(shared) private(i,j,k,cmax)
   do k=kfirst,klast
        cymax(k) = ZERO
       do j=js2g0,jlast
            cmax = ZERO
          do i=1,im
            cmax = max( abs(cy(i,j,k)), cmax)
          enddo
            cymax(k) = max(cymax(k), cmax*(ONE + sine(j)**16) )
       enddo
   enddo

   call mp_recv3d_2( commglobal, src, im, jm, km,                     &
                     1, im, jfirst, jlast+1, kfirst, klast,           &
                     1, im, jlast+1, jlast+1, kfirst, klast, cy, mfy)

   call parcollective( comm_y, MAXOP, ktot, cymax )

!---------------------------------------------------------------------
! Determine the required value of nsplt for each level
!---------------------------------------------------------------------
    nsplt(:)  = int( ONE + cymax(:) )
    max_nsplt = maxval( nsplt(:) )		
    call parcollective( comm_z, MAXOP, max_nsplt )  ! Find global max
    nsplt(:)  = max_nsplt

    do k_courant = klast,kfirst,-1
       if( nsplt(k_courant) > 1 ) then
          exit
       end if
    end do
    k_courant = max( kfirst,k_courant )

!$omp  parallel do default(shared) private(i,j,k,frac) schedule(dynamic,1)

#if !defined(USE_OMP)
!CSD$ PARALLEL DO PRIVATE (I, J, K, FRAC)
#endif
 do k=kfirst,klast

    if( nsplt(k) .ne. 1 ) then
       frac = ONE / nsplt(k)
       do j=js2gd,jn2gd                  
          do i=1,im
            cx(i,j,k) =  cx(i,j,k) * frac      ! cx ghosted on N*ng S*ng
          enddo
       enddo

       do j=js2g0,jn2g0
          do i=1,im
            mfx(i,j,k) = mfx(i,j,k) * frac
          enddo
       enddo

       do j=js2g0,jn1g1                     
          do i=1,im
             cy(i,j,k) =  cy(i,j,k) * frac    ! cy ghosted on N
            mfy(i,j,k) = mfy(i,j,k) * frac    ! mfy ghosted on N
          enddo
       enddo
    endif

    do j=js2g0,jn2g0
       do i=1,im
          if (cy(i,j,k)*cy(i,j+1,k) > ZERO) then
             if ( cy(i,j,k) > ZERO) then
                va(i,j,k) = cy(i,j,k)
             else
                va(i,j,k) = cy(i,j+1,k)      ! cy ghosted on N
             endif
          else
             va(i,j,k) = ZERO
          endif
       enddo
    enddo

! Check if FFSL extension is needed.

    do j=js2gd,jn2gd             ! flux needed on N*ng S*ng
       ffsl(j,k) = .false.
       do i=1,im
          if ( abs(cx(i,j,k)) > ONE ) then  ! cx ghosted on N*ng S*ng
             ffsl(j,k) = .true.
             exit
          endif
       enddo
    enddo

! Scale E-W mass fluxes by CX if FFSL
    do j=js2g0,jn2g0
       if ( ffsl(j,k) ) then
          do i=1,im
             flx(i,j,k) = mfx(i,j,k) / sign( max(abs(cx(i,j,k)), tiny), &
                                        cx(i,j,k) )
          enddo
       else
          do i=1,im
             flx(i,j,k) = mfx(i,j,k)
          enddo
       endif
    enddo

 enddo
#if !defined(USE_OMP)
!CSD$ END PARALLEL DO
#endif

#ifdef TIMING_BARRIERS
 call mp_barrier(commglobal)
#endif

 do 6000 it=1, max_nsplt
    if ( it == 1 ) then
       kend = klast       !  The entire vertical slab needs to be considered
    else
       kend = k_courant   !  Only the subset including courant # > 1 considered
    endif
! WS 05.04.06:  send only the first tracer the rest at end of do iq loop
!               NOTE: there is per definition at least one tracer
    if ( tracer(1)%is_r4 ) then
      q_r8(1:im,jfirst:jlast,kfirst:kend,1) = &
                 tracer(1)%content_r4(1:im,jfirst:jlast,kfirst:kend)
    else
      q_r8(1:im,jfirst:jlast,kfirst:kend,1) = &
                 tracer(1)%content(1:im,jfirst:jlast,kfirst:kend)
    endif
    call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast,         &
                       kfirst, kend, ng, ng, q_r8(1,jfirst-ng,kfirst,1) )

!$omp parallel do default(shared) private(i,j,k,sum1,sum2)
    do k=kfirst,kend

       if (it <= nsplt(k)) then

          do j=js2g0,jn2g0
             do i=1,im-1
                dp1(i,j,k) =  delp(i,j,k) + mfx(i,j,k) - mfx(i+1,j,k) +  &
                              (mfy(i,j,k) - mfy(i,j+1,k)) * acosp(j)
             enddo
             dp1(im,j,k) = delp(im,j,k) + mfx(im,j,k) - mfx(1,j,k) +  &
                              (mfy(im,j,k) - mfy(im,j+1,k)) * acosp(j)
          enddo

          if ( jfirst == 1  ) then
             sum1 = ZERO
             do i=1,im
                sum1 = sum1 + mfy(i,2,k)
             end do

             sum1 = - sum1 * rcap
             do i=1,im
                dp1(i,1,k) = delp(i,1,k) +  sum1
             enddo
          endif

          if ( jlast == jm ) then
             sum2 = ZERO
             do i=1,im
                sum2 = sum2 + mfy(i,jm,k)
             end do

             sum2 = sum2 * rcap
             do i=1,im
                dp1(i,jm,k) = delp(i,jm,k) +  sum2
             enddo
          endif
       endif
    enddo

    do iq = 1, nq
!
! The buffer indices are exchanged, so that cur points to the current buffer,
! while nxt points to the one which will be used next.
!
       if ( mod(iq,2) == 0 ) then
          cur = 2
          nxt = 1
       else
          cur = 1
          nxt = 2
       endif
       call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast,        &
                          kfirst, kend, ng, ng, q_r8(1,jfirst-ng,kfirst,cur) )

!
! Pre-send the next tracer
!
       if ( iq < nq ) then
          if ( tracer(iq+1)%is_r4 ) then
             q_r8(1:im,jfirst:jlast,kfirst:kend,nxt) = &
                  tracer(iq+1)%content_r4(1:im,jfirst:jlast,kfirst:kend)
          else
             q_r8(1:im,jfirst:jlast,kfirst:kend,nxt) = &
                  tracer(iq+1)%content(1:im,jfirst:jlast,kfirst:kend)
          endif
          call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast,   &
                             kfirst, kend, ng, ng, q_r8(1,jfirst-ng,kfirst,nxt) )
       endif

!$omp parallel do default(shared)    &
!$omp private(i, j, k, kq, fx, fy, a2)
#if (!defined USE_OMP) && (!defined COUP_CSM)
!CSD$ PARALLEL DO PRIVATE (I, J, K, KQ, FX, FY, A2)
#endif
       do k=kfirst,kend
          if ( it <= nsplt(k) ) then
             call tp2c(a2, va(1:,jfirst:,k), q_r8(1:,jfirst-ng:,k,cur),&
                       cx(1:,jfirst-ng:,k) , cy(1:,jfirst:,k),       &
                       im, jm, iord, jord, ng,                       &
                       fx(1:,jfirst:), fy(1:,jfirst:), ffsl(1:,k),   &
                       rcap, acosp, flx(1:,jfirst:,k), mfy(1:,jfirst:,k),&
                       cosp, 1, jfirst, jlast )

             do j=jfirst,jlast
                do i=1,im
                   q_r8(i,j,k,cur) = q_r8(i,j,k,cur)*delp(i,j,k) + a2(i,j)
                enddo
             enddo

             if (fill) call fillxy (q_r8(1:,jfirst:,k,cur), im, jm, jfirst, &
                                    jlast, acap, cosp, acosp)
             if ( tracer(iq)%is_r4 ) then
                do j=jfirst,jlast
                   do i=1,im
                      tracer(iq)%content_r4(i,j,k) = q_r8(i,j,k,cur) / dp1(i,j,k)
                   enddo
                enddo
             else
                do j=jfirst,jlast
                   do i=1,im
                      tracer(iq)%content(i,j,k) = q_r8(i,j,k,cur) / dp1(i,j,k)
                   enddo
                enddo
             endif
          endif
       enddo
#if (!defined USE_OMP) && (!defined COUP_CSM)
!CSD$ END PARALLEL DO
#endif
    enddo  ! End of do iq=1, nq

!$omp parallel do private(i, j, k) schedule( dynamic,1 )
    do k=kfirst,kend
       if ( it <= nsplt(k) ) then   ! WS: apparent bug in trac2d, should be <= not <
          do j=jfirst,jlast
             do i=1,im
                delp(i,j,k) = dp1(i,j,k)
             enddo
          enddo
       endif
    enddo

6000  continue

!  Transpose back to XY decomposition
 call mp_sendirr( commglobal, delp, transposes%yz2xy%SendDesc,&
                  transposes%yz2xy%RecvDesc, delpxy )
 call mp_recvirr( commglobal, delpxy, transposes%yz2xy%RecvDesc )

 do i=1, nq
   if ( tracer(i)%is_r4 ) then
     call mp_sendirr_r4( commglobal, tracer(i)%content_r4,         &
                         transposes%yz2xy%SendDesc,                &
                         transposes%yz2xy%RecvDesc,                &
                         tracerxy(i)%content_r4 )
     call mp_recvirr_r4( commglobal, tracerxy(i)%content_r4,       &
                         transposes%yz2xy%RecvDesc )
     deallocate( tracer(i)%content_r4 )
   else
     call mp_sendirr( commglobal, tracer(i)%content,               &
                      transposes%yz2xy%SendDesc,                   &
                      transposes%yz2xy%RecvDesc,                   &
                      tracerxy(i)%content )
     call mp_recvirr( commglobal, tracerxy(i)%content,             &
                      transposes%yz2xy%RecvDesc )
     deallocate( tracer(i)%content )
   endif
 enddo

!
! YZ no longer needed:  deallocate space from heap
!
 deallocate( fx, fy, a2, cymax )
 deallocate( va, flx, dp1 )
 deallocate( cx, cy, mfy )
 deallocate( delp, mfx )
 deallocate( tracer )
 deallocate( q_r8 )

 return

 end subroutine FVAdv_Run
!-----------------------------------------------------------------------


end module FVAdvMod
