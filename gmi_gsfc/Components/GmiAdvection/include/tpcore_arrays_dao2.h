
!=============================================================================
!
! $Id: tpcore_arrays_dao2.h,v 1.5 2011-08-09 22:12:57 mrdamon Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   tpcore_arrays_dao2.h
!
! DESCRIPTION
!   This include file tbd
!
!=============================================================================


!     ----------------
!     Integer commons.
!     ----------------

!     --------------------------------------------------------------------
!     jn : northward of latitude index = jn; Courant numbers could be > 1,
!          so use the flux-form semi-Lagrangian scheme
!     js : southward of latitude index = js; Courant numbers could be > 1,
!          so use the flux-form semi-Lagrangian scheme
!     --------------------------------------------------------------------

      integer, pointer :: jn(:)
      integer, pointer :: js(:)


!     =====================
      common  / gmitad_i1 /  &
!     =====================
     &  jn, js


!     -------------
!     Real commons.
!     -------------

!     -----------------------------------------------------------------
!     dps_ctm : CTM surface pressure tendency; sum over vertical of dpi
!               calculated from original mass fluxes (mb)
!     -----------------------------------------------------------------

      real*8, pointer :: dps_ctm(:,:)

!     ------------------------------------------------------------------------
!     dpi : divergence at a grid point; used to calculate vertical motion (mb)
!     wz  : large scale mass flux (per time step tdt) in the vertical
!           direction as diagnosed from the hydrostatic relationship (mb)
!     ------------------------------------------------------------------------

      real*8, pointer :: dpi(:,:,:)
      real*8, pointer :: wz (:,:,:)

!     ----------------------------------------------------------------------
!     adx    : cross term due to E-W advection (mixing ratio)
!     ady    : cross term due to N-S advection (mixing ratio)
!     crx    : Courant number in E-W direction
!     cry    : Courant number in N-S direction
!     delp1  : pressure thickness, the psudo-density in a hydrostatic system
!              at t1 (mb)
!     delp2  : vertical pressure change for each CTM box at t1+tdt (mb)
!     delpm  : pressure thickness, the psudo-density in a hydrostatic system
!              at t1+tdt/2 (approximate) (mb)
!     dq1    : species density (mb)
!     pu     : pressure at edges in "u"  (mb)
!     qqu    : concentration contribution from E-W advection (mixing ratio)
!     qqv    : concentration contribution from N-S advection (mixing ratio)
!     ua     : average of Courant numbers from il and il+1
!     va     : average of Courant numbers from ij and ij+1
!     ----------------------------------------------------------------------

      real*8, pointer :: adx  (:,:,:)
      real*8, pointer :: ady  (:,:,:)
      real*8, pointer :: crx  (:,:,:)
      real*8, pointer :: cry  (:,:,:)
      real*8, pointer :: delp1(:,:,:)
      real*8, pointer :: delp2(:,:,:)
      real*8, pointer :: delpm(:,:,:)
      real*8, pointer :: dq1  (:,:,:)
      real*8, pointer :: pu   (:,:,:)
      real*8, pointer :: qqu  (:,:,:)
      real*8, pointer :: qqv  (:,:,:)
      real*8, pointer :: ua   (:,:,:)
      real*8, pointer :: va   (:,:,:)

!     =====================
      common  / gmitad_r1 /  &
!     =====================
     &  dps_ctm,  &
!
     &  dpi,  &
     &  wz,  &
!
     &  adx, ady,  &
     &  crx, cry,  &
     &  delp1, delp2, delpm,  &
     &  dq1,  &
     &  pu,  &
     &  qqu, qqv,  &
     &  ua, va

