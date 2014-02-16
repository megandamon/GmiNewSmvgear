module FVAdvUtilsMod

 use shr_kind_mod, only: r4 => shr_kind_r4, r8 => shr_kind_r8

 private
 public   FVAdvUtilsTestCase
 real(r8), parameter :: zero=0.0_r8, half=0.5_r8, one=1.0_r8, two=2.0_r8

CONTAINS

!-------------------------------------------------------------------------
 subroutine FVAdvUtilsTestCase(nx, ny, im, jm, km, imxy, jmxy,            &
                               dt, pi, ae, velocityScale,                 &
                               uc, vc, mx_ur, my_ur, mx, my, mz,          &
                               delp, dpedt, content)
!-------------------------------------------------------------------------

 use mod_comm, only : gid, commglobal, mp_sendirr, mp_recvirr
 use decompmodule, only : DecompType, DecompCreate, DecompFree
 use parutilitiesmodule, only : CommGlobal, ParPatternType, ParInit,      &
                                ParPatternCreate, ParPatternFree

 implicit none

!-------
! Input
!-------

 integer, intent(in) :: nx      ! Comp grid X dim
 integer, intent(in) :: ny      ! Comp grid Y dim

 integer, intent(in) :: im      ! Global E-W dimension
 integer, intent(in) :: jm      ! Global N-S dimension
 integer, intent(in) :: km      ! Vertical dimension

 integer, intent(in) :: imxy(:) ! Distribution in X
 integer, intent(in) :: jmxy(:) ! Distribution in Y

 real(r8), intent(in):: dt      ! Time step in seconds
 real(r8), intent(in):: pi      ! Application's value of PI
 real(r8), intent(in):: ae      ! Earth's radius (m)

 real(r8), intent(in):: velocityScale          ! scale of velocity

 real,   intent(out) :: uc(:,:,:)              ! C-wind easterly
 real,   intent(out) :: vc(:,:,:)              ! C-wind northerly
 real,   intent(out) :: mx_ur(:,:,:)           ! Unremapped mass flux X (Pa m^2 s-1)
 real,   intent(out) :: my_ur(:,:,:)           ! Unremapped mass flux Y (Pa m^2 s-1)
 real,   intent(out) :: mx(:,:,:)              ! Mass flux X (Pa m^2 s-1)
 real,   intent(out) :: my(:,:,:)              ! Mass flux Y (Pa m^2 s-1)
 real,   intent(out) :: mz(:,:,:)              ! Mass flux Z (kg m^2 s-1)
 real,   intent(out) :: delp(:,:,:)            ! Delta pressure
 real,   intent(out) :: dpedt(:,:,:)           ! External increment to PE
 real,   intent(out) :: content(:,:,:)         ! Tracer content

!-----
! Local
!-----

 type (DecompType)               :: global2d, xy2d 
 type (ParPatternType)           :: g2l2d, g2l3d

 real(r8) :: clat(jm)    ! latitude in radian
 real(r8) :: dtdx5(jm)
 real(r8) :: dtdy5(jm)
 real(r8) :: cosp(jm)
 real(r8) :: cose(jm)
 real(r8) ::  gw(jm)
 real(r8) :: rgw(jm)

 real(r8) :: ps1(im,jm)

 real(r8) :: elat(jm+1)    ! cell edge latitude in radian
 real(r8) :: sine(jm+1)
 real(r8) :: dlat(jm)      ! delta-latitude in radian
 real(r8) :: dl, dp, lat, lon
 real(r8) :: angle

 integer  ::  i, j, k
 integer, allocatable :: xdist(:), ydist(:)

 real(r8), allocatable :: g_u(:,:,:), g_v(:,:,:), g_tr(:,:,:)
 real(r8), allocatable :: g_fx(:,:,:), g_fy(:,:,:), g_delp(:,:,:)

 real(r8), allocatable :: tmp(:,:,:)

!
! Due to a Bug (or Feature) in MAPL, this routine may be called
! before the initialization of FVadvcore (i.e., FVAdvInit).  Thus,
! PILGRIM needs to be initialized here too.

 call ParInit()   !  Initialize PILGRIM

! Temporary array has same size as output arrays but is R8
 allocate(tmp(size(delp,1),size(delp,2),size(delp,3)))

 allocate(xdist(nx), ydist(ny))
 xdist = 0;  ydist = 0;
 xdist(1) = im; ydist(1)=jm

 call decompcreate( nx, ny, xdist, ydist, global2d )
 deallocate(xdist, ydist)
 call decompcreate( nx, ny, imxy, jmxy, xy2d )
!
! BEWARE:  this next line causes the code to crash.  Not sure why...
!          Is PILGRIM initialized at this point?????
 call parpatterncreate(commglobal, global2d, xy2d, g2l2d, mod_method=0 )
 call parpatterncreate(commglobal, g2l2d, g2l3d, km )
 call parpatternfree( commglobal, g2l2d )

 call decompfree( global2d )
 call decompfree( xy2d )

 dl      = two*pi / im
 dp      = pi / (jm-1)
 angle   = half*pi - 0.05_r8   ! Slightly offset from pole

 lat     = -half*pi         ! S. Pole
 elat(1) =  lat
 sine(1) = -one
 cose(1) =  zero
 clat(1) =  lat

 do j=2,jm
    lat     = lat + dp
    clat(j) = lat
    elat(j) = half*(clat(j-1) + clat(j))
    sine(j) = sin(elat(j))
    cose(j) = cos(elat(j))
 enddo
 elat(jm+1) = half*pi       ! N. Pole
 sine(jm+1) = one

 dlat(1) = two*(elat(2) - elat(1))  ! Polar cap
 do j=2,jm-1
    dlat(j) = elat(j+1) - elat(j)
 enddo
 dlat(jm) = two*(elat(jm+1) - elat(jm))    ! Polar cap

 do j=1,jm
    gw(j)    = sine(j+1) - sine(j)
    cosp(j)  = gw(j) / dlat(j)
    rgw(j)   = one / gw(j)
    dtdx5(j) = half * dt / (dl*ae*cosp(j))
    dtdy5(j) = half * dt / (ae*dlat(j))
 enddo

 if ( gid == 0 ) then

!
! This code is being executed on the root PET because air_mass_flux
! is not available for an XY decomposition.  Since this is part of
! the initialization, it should not be an issue.
!
    allocate( g_u(im,jm,km), g_v(im,jm,km) )
    allocate( g_fx(im,jm,km), g_fy(im,jm,km), g_delp(im,jm,km) )
    allocate( g_tr(im,jm,km) )

    do k = 1, km
       lat = -HALF*PI
       do j = 1, jm
          lon = -PI
          do i = 1, im
             lon = lon + dl
             g_u(i,j,k) =  velocityScale * ( COS(lat)*COS(ANGLE) +          &
                                             SIN(lat)*COS(lon)*SIN(ANGLE) )
             g_v(i,j,k) = -velocityScale * SIN(lon)*SIN(ANGLE)
             g_tr(i,j,k) = pudding(PI,AE,lon,lat)
          enddo
          lat = lat + dp
       enddo
    enddo

!
! Define ps1  -- is there any reason not to use a constant?
!
    ps1 = 100000.0_r8
    call air_mass_flux(im, jm, km, ps1, g_u, g_v, g_fx, g_fy, g_delp, &
                       dtdx5, dtdy5, cosp, cose, rgw, gw,   &
                       dt, ae)
 endif

 call mp_sendirr( commglobal, g_tr, g2l3d%SendDesc, g2l3d%RecvDesc, tmp )
 call mp_recvirr( commglobal, tmp, g2l3d%RecvDesc )
 content = tmp

 call mp_sendirr( commglobal, g_u, g2l3d%SendDesc, g2l3d%RecvDesc, tmp )
 call mp_recvirr( commglobal, tmp, g2l3d%RecvDesc )
 uc = tmp

 call mp_sendirr( commglobal, g_v, g2l3d%SendDesc, g2l3d%RecvDesc, tmp )
 call mp_recvirr( commglobal, tmp, g2l3d%RecvDesc )
 vc = tmp

 call mp_sendirr( commglobal, g_fx, g2l3d%SendDesc, g2l3d%RecvDesc, tmp )
 call mp_recvirr( commglobal, tmp, g2l3d%RecvDesc )

!$omp parallel do private(i,j,k)
 do k = 1,km
    do j = 1,size(mx,2)
       do i = 1,size(mx,1)
!
! TEMPORARY:  initialize mx_ur and mx to same values
          mx_ur(i,j,k) = tmp(i,j,k) * (dl*ae*cosp(j)) * ae*dp / dt ! Pa m^2 / s
          mx(i,j,k) = tmp(i,j,k) * (dl*ae*cosp(j)) * ae*dp / dt ! Pa m^2 / s
       enddo
    enddo
 enddo


 call mp_sendirr( commglobal, g_fy, g2l3d%SendDesc, g2l3d%RecvDesc, tmp )
 call mp_recvirr( commglobal, tmp, g2l3d%RecvDesc )

!$omp parallel do private(i,j,k)
 do k = 1,km
    do j = 1,size(my,2)
       do i = 1,size(my,1)
!
! TEMPORARY:  initialize my_ur and my to same values
          my_ur(i,j,k) = tmp(i,j,k) * dl*ae * ae*dp / dt ! Pa m^2 / s  -- cos factor removed
          my(i,j,k) = tmp(i,j,k) * dl*ae * ae*dp / dt ! Pa m^2 / s  -- cos factor removed
       enddo
    enddo
 enddo

 call mp_sendirr( commglobal, g_delp, g2l3d%SendDesc, g2l3d%RecvDesc, tmp )
 call mp_recvirr( commglobal, tmp, g2l3d%RecvDesc )
 delp = tmp

!
! For the time being, leave these two zero.
!
 mz = 0.0
 dpedt = 0.0

 if ( gid == 0 ) then

    deallocate( g_u, g_v )
    deallocate( g_fx, g_fy, g_delp )
    deallocate( g_tr )

 endif

 end subroutine FVAdvUtilsTestCase
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  mfz_comp --- Calculate vertical mass flux
!
! !INTERFACE:
      subroutine mfz_comp( ae, grav, dt, commxy_x, im, jm, km,           &
                           ifirstxy, ilastxy, jfirstxy, jlastxy,         &
                           iam, myidxy_y, nprxy_x, nprxy_y,              &
                           dl, dp, acap, cosp, bk, mfx, mfy, mfz )
! !USES:
      use shr_kind_mod, only : r8 => shr_kind_r8
      use parutilitiesmodule, only: bcstop, sumop,  parcollective
      use mod_comm, only :  gid, commglobal, mp_send3d, mp_recv3d
      implicit none

! !INPUT PARAMETERS:
      real(r8), intent(in) :: ae                ! Radius of the Earth (m)
      real(r8), intent(in) :: grav              ! Gravity
      integer,  intent(in) :: dt                ! large time step in seconds
      integer,  intent(in) :: im, jm, km
      integer,  intent(in) :: iam, myidxy_y, nprxy_x, nprxy_y
      integer,  intent(in) :: ifirstxy, ilastxy, jfirstxy, jlastxy
      integer,  intent(in) :: commxy_x

      real(r8), intent(in) :: dl, dp, acap
      real(r8), intent(in) :: cosp(:)
      real(r8), intent(in) :: bk(:)

      real(r8), intent(inout) :: mfx(ifirstxy:ilastxy,jfirstxy:jlastxy,km) 
      real(r8), intent(inout) :: mfy(ifirstxy:ilastxy,jfirstxy:jlastxy,km) 

! !OUTPUT PARAMETERS:
      real(r8), intent(inout) :: mfz(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1)

! !DESCRIPTION:
!     Compute cell centered vertical mass flux
!        Vertical integration of the continuity equation using the 
!        divergence of the horizontal mass fluxes
!
! !REVISION HISTORY:
!   WMP  05.12.01   Created
!   WBS  07.07.30   Adapted for FVadvcore
!
! !BUGS:
!   Not yet tested...
!
!EOP
!---------------------------------------------------------------------
!BOC

      real(r8) mfy_north(ifirstxy:ilastxy,km)
      real(r8), allocatable :: mfx_east(:,:)     ! East halo
      real(r8), allocatable :: mfy_sp(:,:), mfy_np(:,:)
      real(r8) delp1xy(im,jm,km)
      real(r8) conv(ifirstxy:ilastxy,jfirstxy:jlastxy,km)
      real(r8) pit(ifirstxy:ilastxy,jfirstxy:jlastxy) ! pressure tendency
      real(r8) area, sum1, sum2
      integer  :: ierr, dest, src    ! SPMD related
      integer  :: i,j,k

      allocate(mfx_east(jfirstxy:jlastxy,km))     ! East halo
      allocate(mfy_sp(im,km), mfy_np(im,km) )

      call mp_send3d( commglobal, iam-nprxy_x, iam+nprxy_x, im, jm, km,  &
                      ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km,       &
                      ifirstxy, ilastxy, jfirstxy, jfirstxy, 1, km, mfy )
      call mp_recv3d( commglobal, iam+nprxy_x, im, jm, km,               &
                      ifirstxy, ilastxy, jlastxy+1, jlastxy+1, 1, km,    &
                      ifirstxy, ilastxy, jlastxy+1, jlastxy+1, 1, km,    &
                      mfy_north )
      if (nprxy_x > 1) then
! Nontrivial x decomposition
        dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
        src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
        call mp_send3d( commglobal, dest, src, im, jm, km,               &
                        ifirstxy, ilastxy, jfirstxy, jlastxy, 1,km,      &
                        ifirstxy, ifirstxy, jfirstxy, jlastxy, 1, km, mfx )
      endif

!         
! Prepare sum of mfy for divergence at poles
!
      mfy_sp = 0.0_r8
      mfy_np = 0.0_r8
!$omp  parallel do                  &
!$omp  default(shared)              &
!$omp  private(i,k)
      do k=1,km
        if ( jfirstxy == 1 ) then  ! SP
          do i=ifirstxy,ilastxy
            mfy_sp(i,k) = mfy(i,2,k)
          enddo
        endif
        if ( jlastxy == jm ) then  ! NP
          do i=ifirstxy,ilastxy
            mfy_np(i,k) = mfy(i,jm,k)
          enddo
        endif
! Periodic boundary  (for the case of no decomposition in X)
        do j=jfirstxy,jlastxy
           mfx_east(j,k) = mfx(ifirstxy,j,k)
        enddo
      enddo

      if ( nprxy_x > 1 ) then
! Non-trivial X decomposition
        call mp_recv3d( commglobal, src, im, jm, km,                             &
                        ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1, km,          &
                        ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1, km, mfx_east )
      endif
!                       
! Collect on all PEs the mfy at both poles
!
      if (nprxy_x > 1) then
        call parcollective(commxy_x, sumop, im, km, mfy_sp)
        call parcollective(commxy_x, sumop, im, km, mfy_np)
      endif

! Compute Convergence of the horizontal Mass flux
!$omp  parallel do                  &
!$omp  default(shared)              &
!$omp  private(i,j,k,sum1,sum2)
      do k=1,km
          do j=jfirstxy,jlastxy-1
             do i=ifirstxy,ilastxy-1
                conv(i,j,k) =  mfx(i,j,k) - mfx(i+1,j,k) +  &
                               mfy(i,j,k) - mfy(i,j+1,k)
             enddo
             conv(ilastxy,j,k) = mfx(ilastxy,j,k) - mfx_east(j,k) +     &
                                 mfy(ilastxy,j,k) - mfy(ilastxy,j+1,k)
          enddo
          j = jlastxy
          do i=ifirstxy,ilastxy-1
             conv(i,j,k) =  mfx(i,j,k) - mfx(i+1,j,k) +  &
                            mfy(i,j,k) - mfy_north(i,k)
          enddo
          conv(ilastxy,j,k) = mfx(ilastxy,j,k) - mfx_east(j,k) +     &
                              mfy(ilastxy,j,k) - mfy_north(ilastxy,k)

! Poles
          if ( jfirstxy == 1  ) then
            sum1 = -SUM(mfy_sp(1:im,k))
            do i=ifirstxy,ilastxy
              conv(i,1,k) = sum1
            enddo
          endif

          if ( jlastxy  == jm ) then
            sum2 =  SUM(mfy_np(1:im,k))
            do i=ifirstxy,ilastxy
              conv(i,jm,k) = sum2
            enddo
          endif
      enddo

! Surface pressure tendency
      pit(:,:) = 0.0
      do k=1,km
         do j=jfirstxy,jlastxy
            do i=ifirstxy,ilastxy
               pit(i,j) = pit(i,j) + conv(i,j,k)
            enddo
         enddo
      enddo

!  Sum over levels
      do k=2,km
         do j=jfirstxy,jlastxy
            do i=ifirstxy,ilastxy
               conv(i,j,k) = conv(i,j,k) + conv(i,j,k-1)
            enddo
         enddo
      enddo

      mfz(:,:,:) = 0.0
      do k=2,km
         do j=MAX(2,jfirstxy),MIN(jlastxy,jm-1)
            area = dl*cosp(j)*ae * dp*ae
            do i=ifirstxy,ilastxy
               mfz(i,j,k) = ( conv(i,j,k-1)  - bk(k)*pit(i,j) )/(grav*area)  ! Kg/m^2/s
            enddo
         enddo
! Poles
          if ( jfirstxy == 1  ) then
            j=1
            area = acap*(dl*ae * dp*ae) 
            do i=ifirstxy,ilastxy
               mfz(i,j,k) = ( conv(i,j,k-1)  - bk(k)*pit(i,j) ) / (grav*area)  ! Kg/m^2/s
            enddo
          endif
          if ( jlastxy  == jm ) then
            j=jm
            area = acap*(dl*ae * dp*ae)                      
            do i=ifirstxy,ilastxy
               mfz(i,j,k) = ( conv(i,j,k-1)  - bk(k)*pit(i,j) ) / (grav*area)  ! Kg/m^2/s
            enddo
          endif
      enddo

      deallocate(mfx_east)
      deallocate(mfy_sp,mfy_np)

      return
      end subroutine mfz_comp
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! Think about how to best realize this routine, needed only
! to get the mass fluxes in the example application
! It is outdated with respect to the FV dycore
!
 subroutine air_mass_flux(im, jm, km, ps1, u, v, fx, fy, delp, &
                          dtdx5, dtdy5, cosp, cose, rgw, gw,   &
                          dt, ae)

! !USES:
 use tp_core, only: tp2d

!------------------------------------------------------
! The hybrid ETA-coordinate:
! pressure at layer edges are defined as follows:
!
!        p(i,j,k) = ak(k) + bk(k)*ps(i,j)          (1)
!------------------------------------------------------

! Input from Data/Model:
! (u,v) is the time mean wind at Time=t+dt/2  
! delp1 is the layer thickness at Time=t

! Output:
! delp is the predicted thickness at Time=t+dt
! (fx,fy): background air mass flxues
! (cx,cy): CFL number

 implicit none

 integer, intent(in):: im
 integer, intent(in):: jm
 integer, intent(in):: km
 
 real(r8), intent(in):: dt
 real(r8), intent(in):: ae
 real(r8), intent(in):: ps1(im,jm)
 real(r8), intent(in):: u(im,jm,km)
 real(r8), intent(in):: v(im,jm,km)

 real(r8), intent(in):: dtdx5(:)
 real(r8), intent(in):: dtdy5(:)
 real(r8), intent(in):: cosp(:)
 real(r8), intent(in):: cose(:)
 real(r8), intent(in):: rgw(:)
 real(r8), intent(in):: gw(:)


! Output:
 real(r8), intent(out):: delp (im,jm,km)
 real(r8), intent(out):: fx(im,jm,km)
 real(r8), intent(out):: fy(im,jm,km)

! Local:
 integer, parameter   :: iord = 3
 integer, parameter   :: jord = 3

 real(r8)  :: ak(km+1)
 real(r8)  :: bk(km+1)

 real(r8)  :: delp1(im,jm,km)
 real(r8)  :: yms(im,jm,km)

 real(r8)  :: cx(im,jm,km)
 real(r8)  :: cy(im,jm,km)
 real(r8)  :: va(im,jm,km)

 logical   :: ffsl(jm,km)

 real(r8), parameter ::  tiny = 1.e-10
 real(r8)  :: dak, dbk
 real(r8)  :: dtoa, vt
 real(r8)  :: ptop, pint
 integer   :: i,j,k
 integer   :: ks

 dtoa = .5*dt/ae 

 call set_eta(km, ks, ptop, pint, ak, bk)

!$omp parallel do private(i, j, k, vt)

  do k=1,km
     do j=2, jm
         do i=1,im
            vt = v(i,j,k) + v(i,j-1,k)
            if ( vt > 0. ) then
               cy(i,j,k) = dtdy5(j-1)*vt
            else
               cy(i,j,k) = dtdy5(j)*vt
            endif
             yms(i,j,k) = dtoa*vt*cose(j)
         enddo
     enddo

     do j=2,jm-1
        do i=1,im
           if( cy(i,j,k)*cy(i,j+1,k) > 0. ) then
              if( cy(i,j,k) > 0. ) then
                  va(i,j,k) = cy(i,j,k)
              else
                  va(i,j,k) = cy(i,j+1,k)         
              endif
           else
              va(i,j,k) = 0.
          endif
        enddo
     enddo

     do j=2,jm-1
           cx(1,j,k) = dtdx5(j)*(u(1,j,k)+u(im,j,k))
        do i=2,im
           cx(i,j,k) = dtdx5(j)*(u(i,j,k)+u(i-1,j,k))
        enddo
     enddo

  enddo

!---------------------------------------------------
! Compute background mass-flux (fx, fy) and (cx, cy) 
!---------------------------------------------------

!$omp parallel do                             &
!$omp shared(im,jm,iord,jord) &
!$omp private(i, j, k, dak, dbk)

  do k=1,km

     do j=2,jm-1
        ffsl(j,k) = .false.
        do i=1,im
           if( abs(cx(i,j,k)) > 1. ) then
               ffsl(j,k) = .true.
               go to 2222
           endif
        enddo
2222  continue
     enddo

     dak = ak(k+1) - ak(k)
     dbk = bk(k+1) - bk(k)

     do j=1,jm
        do i=1,im
           delp1(i,j,k) = dak + dbk*ps1(i,j)
        enddo
     enddo 

     ! mg = 0 (9th argument):   No ghosting in this sequential version
     ! iv = 0 (15th argument):  
     ! jfirst = 1 (16th argument):   starting latitude (seq. version)
     ! jlast  = jm (17th argument):    ending latitude (seq. version)
     call tp2d(va(:,:,k), delp1(:,:,k), cx(:,:,k),                       &
               cy(:,:,k), im, jm, iord, jord, 0,                         &
               fx(:,:,k), fy(:,:,k),  ffsl(:,k),                         &
               cx(:,:,k), yms(:,:,k), cosp, 0, 1, jm)

     do j=2,jm-1
        do i=1,im-1
           delp(i,j,k) = delp1(i,j,k) + fx(i,j,k) - fx(i+1,j,k) +        &
                                        (fy(i,j,k)-fy(i,j+1,k))*rgw(j)
        enddo
        delp(im,j,k) = delp1(im,j,k) + fx(im,j,k) - fx(1,j,k) +          &
                                       (fy(im,j,k)-fy(im,j+1,k))*rgw(j)
      enddo

      do i=1,im
         delp(i,1,k) = delp1(i,1,k) - fy(i,2,k)*rgw(1)
      enddo
      call xpavg(delp(:,1,k), im)

      do i=1,im
         delp(i,jm,k) = delp1(i,jm,k) + fy(i,jm,k)*rgw(jm)
      enddo
      call xpavg(delp(:,jm,k), im)

     do j=2,jm-1
        if( ffsl(j,k) ) then
          do i=1,im
             fx(i,j,k) = fx(i,j,k)/sign(max(abs(cx(i,j,k)),tiny),cx(i,j,k))
          enddo
        endif
     enddo

  enddo

 end subroutine air_mass_flux
!--------------------------------------------------------------


!--------------------------------------------------------------
 subroutine xpavg(p, im)
! !USES:

      implicit none

! !INPUT PARAMETERS:
      integer im

! !INPUT/OUTPUT PARAMETERS:
      real(r8) p(im)

      integer i
      real(r8) sum1

      sum1 = 0.
      do i=1,im
         sum1 = sum1 + p(i)
      enddo
      sum1 = sum1 / im

      do i=1,im
         p(i) = sum1
      enddo
 end subroutine xpavg
!--------------------------------------------------------------

!--------------------------------------------------------------
 subroutine set_eta(plev, ks, ptop, pint, ak, bk)

 implicit none
      integer,  intent(in)  ::  plev
      integer,  intent(out) ::  ks
      real(r8), intent(out) ::  ak(plev+1), bk(plev+1)
      real(r8), intent(out) ::  ptop, pint

! NCAR specific
      real(r8) ::  a18(19),b18(19)              ! CCM3
      real(r8) ::  a26(27),b26(27)              ! CCM4
      real(r8) ::  a30(31),b30(31)              ! CCM4

! NASA only
      real(r8) ::  a30m(31),b30m(31)            ! smoothed CCM4 30-L
      real(r8) ::  a32(33),b32(33)
      real(r8) ::  a32old(33),b32old(33)
      real(r8) ::  a55(56),b55(56)
      real(r8) ::  a55old(56),b55old(56)
      real(r8) ::  a64(65),b64(65)
      real(r8) ::  a96(97),b96(97)

      integer k

! *** NCAR settings ***
      data a18 /291.70,  792.92,  2155.39,  4918.34,  8314.25, &
               7993.08, 7577.38,  7057.52,  6429.63,  5698.38, &
               4879.13, 3998.95,  3096.31,  2219.02,  1420.39, &
               754.13,  268.38,   0.0000,   0.0000 /

      data b18 /0.0000,    0.0000,    0.0000,   0.0000,   0.0000, &
                0.0380541, 0.0873088, 0.1489307, 0.2232996,       &
                0.3099406, 0.4070096, 0.5112977, 0.6182465,       &
                0.7221927, 0.8168173, 0.8957590, 0.9533137,       &
                0.9851122, 1.0  /
     
      data a26 /219.4067,  489.5209,   988.2418,   1805.201,      &
                2983.724,  4462.334,   6160.587,   7851.243,      &
                7731.271,  7590.131,   7424.086,   7228.744,      &
                6998.933,  6728.574,   6410.509,   6036.322,      &
                5596.111,  5078.225,   4468.96,    3752.191,      &
                2908.949,  2084.739,   1334.443,   708.499,       &
                252.136,   0.,         0. /

      data b26 /0.,         0.,         0.,         0.,           &
                0.,         0.,         0.,         0.,           &
                0.01505309, 0.03276228, 0.05359622, 0.07810627,   &
                0.1069411,  0.14086370, 0.180772,   0.227722,     &
                0.2829562,  0.3479364,  0.4243822,  0.5143168,    &
                0.6201202,  0.7235355,  0.8176768,  0.8962153,    &
                0.9534761,  0.9851122,  1.        /

      data a30 /225.523952394724, 503.169186413288, 1015.79474285245,  &
               1855.53170740604, 3066.91229343414,  4586.74766123295,  &
               6332.34828710556, 8070.14182209969,  9494.10423636436,  &
              11169.321089983,  13140.1270627975,  15458.6806893349,   &
              18186.3352656364, 17459.799349308,   16605.0657629967,   &
              15599.5160341263, 14416.541159153,   13024.8308181763,   &
              11387.5567913055,  9461.38575673103,  7534.44507718086,  &
               5765.89405536652, 4273.46378564835,  3164.26791250706,  &
               2522.12174236774, 1919.67375576496,  1361.80268600583,  &
                853.108894079924, 397.881818935275,    0.,             &
                  0.  /

      data b30 /0.,                 0.,                                 &
                0.,                 0.,                0.,              &
                0.,                 0.,                0.,              &
                0.,                 0.,                0.,              &
                0.,                 0.,                0.03935482725501,&
                0.085653759539127,  0.140122056007385, 0.20420117676258,&
                0.279586911201477,  0.368274360895157, 0.47261056303978,&
                0.576988518238068,  0.672786951065063, 0.75362843275070,&
                0.813710987567902,  0.848494648933411, 0.88112789392471,&
                0.911346435546875,  0.938901245594025, 0.96355980634689,&
                0.985112190246582,  1.   /

! *** NASA DAO settings ***

! Smoothed CCM4's 30-Level setup
      data a30m / 300.00000,     725.00000,    1500.00000,     &
             2600.00000,    3800.00000,    5050.00000,         &
             6350.00000,    7750.00000,    9300.00000,         &
            11100.00000,   13140.00000,   15458.00000,         &
            18186.33580,   20676.23761,   22275.23783,         &
            23025.65071,   22947.33569,   22038.21991,         &
            20274.24578,   17684.31619,   14540.98138,         &
            11389.69990,    8795.97971,    6962.67963,         &
             5554.86684,    4376.83633,    3305.84967,         &
             2322.63910,    1437.78398,     660.76994,         &
                0.00000 /

      data b30m / 0.00000,       0.00000,       0.00000,       &
                  0.00000,       0.00000,       0.00000,       &
                  0.00000,       0.00000,       0.00000,       &
                  0.00000,       0.00000,       0.00000,       &
                  0.00000,       0.00719,       0.02895,       &
                  0.06586,       0.11889,       0.18945,       &
                  0.27941,       0.38816,       0.50692,       &
                  0.61910,       0.70840,       0.77037,       &
                  0.81745,       0.85656,       0.89191,       &
                  0.92421,       0.95316,       0.97850,       &
                  1.00000 /

      data a32/40.00000,     100.00000,     200.00000,         &
            370.00000,     630.00000,    1000.00000,           &
           1510.00000,    2160.00000,    2900.00000,           &
           3680.00000,    4535.00000,    5505.00000,           &
           6607.26750,    7851.22980,    9236.56610,           &
          10866.34270,   12783.70000,   15039.30000,           &
          17693.00000,   20119.20876,   21686.49129,           &
          22436.28749,   22388.46844,   21541.75227,           &
          19873.78342,   17340.31831,   13874.44006,           &
          10167.16551,    6609.84274,    3546.59643,           &
           1270.49390,       0.00000,       0.00000   /

      data b32/0.00000,       0.00000,       0.00000,          &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00000,       0.00000,            &
             0.00000,       0.00696,       0.02801,            &
             0.06372,       0.11503,       0.18330,            &
             0.27033,       0.37844,       0.51046,            &
             0.64271,       0.76492,       0.86783,            &
             0.94329,       0.98511,       1.00000   /

      data a32old /300.0000,  454.1491,  652.5746,  891.9637, 1159.7102, &
             1492.8248, 1902.5026, 2400.4835, 2998.6740, 3708.6584,      &
             4541.1041, 5505.0739, 6607.2675, 7851.2298, 9236.5661,      &
            10866.3427, 12420.400, 13576.500, 14365.400, 14807.800,      &
             14915.500, 14691.400, 14129.400, 13214.800, 11923.200,      &
             10220.700,  8062.000,  5849.500,  3777.000,  2017.200,      &
               720.600,     0.000,     0.000 /

      data b32old /0.00, 0.0000000, 0.0000000, 0.0000000, 0.0000000,     &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,     &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,     &
              0.0000000, 0.003633 , 0.014628 , 0.033276 , 0.060071 ,     &
              0.095722 , 0.141171 , 0.197623 , 0.266571 , 0.349839 ,     &
              0.449632 , 0.568589 , 0.685887 , 0.793252 , 0.883128 ,     &
              0.948792 , 0.9851119, 1.0000000 /

      data a55/ 1.00000,       2.00000,       3.27000,              &
              4.75850,       6.60000,       8.93450,                &
             11.97030,      15.94950,      21.13490,                &
             27.85260,      36.50410,      47.58060,                &
             61.67790,      79.51340,     101.94420,                &
            130.05080,     165.07920,     208.49720,                &
            262.02120,     327.64330,     407.65670,                &
            504.68050,     621.68000,     761.98390,                &
            929.29430,    1127.68880,    1364.33920,                &
           1645.70720,    1979.15540,    2373.03610,                &
           2836.78160,    3380.99550,    4017.54170,                &
           4764.39320,    5638.79380,    6660.33770,                &
           7851.22980,    9236.56610,   10866.34270,                &
          12783.70000,   15039.30000,   17693.00000,                &
          20119.20876,   21686.49129,   22436.28749,                &
          22388.46844,   21541.75227,   19873.78342,                &
          17340.31831,   13874.44006,   10167.16551,                &
           6609.84274,    3546.59643,    1270.49390,                &
              0.00000,       0.00000   /

      data b55 / 0.00000,       0.00000,       0.00000,         &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00000,       0.00000,       0.00000,           &
               0.00696,       0.02801,       0.06372,           &
               0.11503,       0.18330,       0.27033,           &
               0.37844,       0.51046,       0.64271,           &
               0.76492,       0.86783,       0.94329,           &
               0.98511,       1.00000  /

      data a55old /1.0000,    2.0000,    3.2700,    4.7585,     6.6000, &
              8.9345,   11.9703,   15.9495,   21.1349,    27.8526,      &
             36.5041,   47.5806,   61.6779,   79.5134,   101.9442,      &
            130.0508,  165.0792,  208.4972,  262.0212,   327.6433,      &
            407.6567,  504.6805,  621.6800,  761.9839,   929.2943,      &
           1127.6888, 1364.3392, 1645.7072, 1979.1554,  2373.0361,      &
           2836.7816, 3380.9955, 4017.5417, 4764.3932,  5638.7938,      &
           6660.3377, 7851.2298, 9236.5661,10866.3427, 12420.400 ,      &
          13576.500 , 14365.400, 14807.800, 14915.500 , 14691.400,      &
          14129.400 , 13214.800, 11923.200, 10220.700 ,  8062.000,      &
           5849.500 ,  3777.000,  2017.200,   720.600,      0.000,      &
              0.000 /

      data b55old /   0.0000000, 0.0000000, 0.0000000,                  &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,    &
              0.0000000, 0.003633 , 0.014628 , 0.033276 , 0.060071 ,    &
              0.095722 , 0.141171 , 0.197623 , 0.266571 , 0.349839 ,    &
              0.449632 , 0.568589 , 0.685887 , 0.793252 , 0.883128 ,    &
              0.948792 , 0.9851119, 1.0000000 /

      data a64/1.00000,       1.80   ,       2.80086,     &
             3.93309,       5.20139,       6.77626,       &
             8.69654,      10.99483,      13.81736,       &
            17.26058,      21.43286,      26.45448,       &
            32.45730,      39.58402,      47.98678,       &
            57.82525,      69.26401,      82.46925,       &
            97.60468,     114.82686,     135.08787,       &
           158.92390,     186.96575,     219.95555,       &
           258.76633,     304.42522,     358.14053,       &
           421.33383,     495.67748,     583.13893,       &
           686.03282,     807.08215,     949.49044,       &
          1117.02644,    1314.12387,    1545.99882,       &
          1818.78771,    2139.70974,    2517.25793,       &
          2961.42386,    3483.96212,    4098.70138,       &
          4821.91034,    5672.72831,    6673.67169,       &
          7851.22983,    9236.56613,   10866.34270,       &
         12783.69059,   15039.35130,   17693.01955,       &
         20814.92310,   23609.16397,   25271.17281,       &
         25844.93368,   25345.63084,   23760.05052,       &
         21046.23129,   17132.35351,   12832.00555,       &
          8646.27815,    5012.23907,    2299.34286,       &
           781.15294,       0.00000  /

      data b64/0.00000,       0.00000,       0.00000,     &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00000,       0.00000,       &
             0.00000,       0.00879,       0.03537,       &
             0.08047,       0.14526,       0.23147,       &
             0.34138,       0.47789,       0.61606,       &
             0.74456,       0.85318,       0.93300,       &
             0.97730,       1.00000  /

      data a96/  1.00000,       2.32782,       3.34990,   &
               4.49484,       5.62336,       6.93048,     &
               8.41428,      10.06365,      11.97630,     &
              14.18138,      16.70870,      19.58824,     &
              22.84950,      26.52080,      30.62845,     &
              35.19588,      40.24273,      45.78375,     &
              51.82793,      58.43583,      65.62319,     &
              73.40038,      81.77154,      90.73373,     &
             100.27628,     110.82243,     122.47773,     &
             135.35883,     149.59464,     165.32764,     &
             182.71530,     201.93164,     223.16899,     &
             246.63988,     272.57922,     301.24661,     &
             332.92902,     367.94348,     406.64044,     &
             449.40720,     496.67181,     548.90723,     &
             606.63629,     670.43683,     740.94727,     &
             818.87329,     904.99493,    1000.17395,     &
            1105.36304,    1221.61499,    1350.09326,     &
            1492.08362,    1649.00745,    1822.43469,     &
            2014.10168,    2225.92627,    2460.02905,     &
            2718.75195,    3004.68530,    3320.69092,     &
            3669.93066,    4055.90015,    4482.46240,     &
            4953.88672,    5474.89111,    6050.68994,     &
            6687.04492,    7390.32715,    8167.57373,     &
            9026.56445,    9975.89648,   11025.06934,     &
           12184.58398,   13466.04785,   14882.28320,     &
           16447.46289,   18177.25781,   20088.97461,     &
           21886.89453,   23274.16602,   24264.66602,     &
           24868.31641,   25091.15430,   24935.41016,     &
           24399.52148,   23478.13281,   22162.01758,     &
           20438.00586,   18288.83984,   15693.01172,     &
           12624.54199,    9584.35352,    6736.55713,     &
            4231.34326,    2199.57910,     747.11890,     &
              0.00000 /

      data b96/0.00000,       0.00000,       0.00000,     &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00000,       0.00000,       0.00000,      &
              0.00315,       0.01263,       0.02853,      &
              0.05101,       0.08030,       0.11669,      &
              0.16055,       0.21231,       0.27249,      &
              0.34169,       0.42062,       0.51005,      &
              0.61088,       0.70748,       0.79593,      &
              0.87253,       0.93400,       0.97764,      &
              1.00000 /

      select case (plev)

! *** Original CCM3 18-Level setup ***
        case (18)
          ks = 4
          do k=1,plev+1
            ak(k) = a18(k)
            bk(k) = b18(k)
          enddo

        case (26)
! CCM4 26-Level setup ***
          ks = 7
          do k=1,plev+1
            ak(k) = a26(k)
            bk(k) = b26(k)
          enddo

        case (30)
! CCM4 30-Level setup ***
          ks = 12
          do k=1,plev+1
            ak(k) = a30(k)
            bk(k) = b30(k)
          enddo

! *** Revised 32-L setup with ptop at 0.4 mb ***
        case (32)
          ks = 18
          do k=1,plev+1
            ak(k) = a32(k)
            bk(k) = b32(k)
          enddo

! *** Revised 55-L setup with ptop at 0.01 mb ***
        case (55)
          ks = 41
          do k=1,plev+1
            ak(k) = a55(k)
            bk(k) = b55(k)
          enddo

! *** Others ***
        case (64)
          ks = 51
          do k=1,plev+1
            ak(k) = a64(k)
            bk(k) = b64(k)
          enddo

        case (96)
          ks = 77
          do k=1,plev+1
            ak(k) = a96(k)
            bk(k) = b96(k)
          enddo

        case DEFAULT

           stop "set_eta: This number of levels not supported"

      end select

          ptop = ak(1)
          pint = ak(ks+1) 

      end subroutine set_eta
!--------------------------------------------------------------

!--------------------------------------------------------------
  function pudding(pi,ae,lon,lat) result(q)
  implicit none
  real(r8), intent(in):: pi, ae, lon, lat
  real(r8)  :: q

  real(r8), parameter :: Q0 = 0.01
  real(r8), parameter :: LAT2  = 0._r8

  real(r8)  :: LONG2

  real(r8) ::  r1,r2

  Long2 = 1.5_r8*PI
  r2=ae/3._r8

  r1=ae*acos(sin(lat2)*sin(lat)+cos(lat2)*cos(lat)*cos(lon-long2))

  if( r1<r2) then
    q=0.5_r8*Q0*(1._r8 + cos(PI*r1/r2))
  else
    q=0._r8
  endif

  if ( abs(lon-long2) < 0.01 .and. abs(lat) < 0.01 ) then
    print *, q, "should be nearly", q0, "at (lat,lon)", lat, lon
  endif

  end function pudding
!--------------------------------------------------------------

!--------------------------------------------------------------
  function RemapBounds_3dr4(A,I1,IM,J1,JM,K1,KM)
    integer,          intent(IN) :: I1,IM,J1,JM,K1,KM
    real(r4), target, intent(IN) :: A(I1:IM,J1:JM,K1:KM)
    real(r4), pointer            :: RemapBounds_3dr4(:,:,:)

    RemapBounds_3dr4 => A
  end function RemapBounds_3dr4
!--------------------------------------------------------------
!--------------------------------------------------------------
  function RemapBounds_3dr8(A,I1,IM,J1,JM,K1,KM)
    integer,          intent(IN) :: I1,IM,J1,JM,K1,KM
    real(r8), target, intent(IN) :: A(I1:IM,J1:JM,K1:KM)
    real(r8), pointer            :: RemapBounds_3dr8(:,:,:)

    RemapBounds_3dr8 => A
  end function RemapBounds_3dr8
!--------------------------------------------------------------


end module FVAdvUtilsMod
!--------------------------------------------------------------

