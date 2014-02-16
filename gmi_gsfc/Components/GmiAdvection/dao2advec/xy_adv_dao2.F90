!
!=============================================================================
!
! $Id: xy_adv_dao2.F90,v 1.9 2013-07-31 14:45:33 ssteenro Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Shian-Jiann Lin, DAO)
!   jrt@llnl.gov
!
! FILE
!   xy_adv_dao2.F
!
! ROUTINES
!   Xadv_Dao2
!   Yadv_Dao2
!
!=============================================================================
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Xadv_Dao2
!
! DESCRIPTION
!   This routine is the advective form E-W operator for computing the adx
!   (E-W) cross term.
!
! ARGUMENTS
!   iad : if iad = 1, use 1st order accurate scheme;
!         if iad = 2, use 2nd order accurate scheme
!   jn  : northward of latitude index = jn, Courant numbers could be > 1,
!         so use the flux-form semi-Lagrangian scheme
!   js  : southward of latitude index = js, Courant numbers could be > 1,
!         so use the flux-form semi-Lagrangian scheme
!   adx : cross term due to E-W advection (mixing ratio)
!   qqv : concentration contribution from N-S advection (mixing ratio)
!   ua  : average of Courant numbers from il and il+1
!
!-----------------------------------------------------------------------------
!
      subroutine Xadv_Dao2  &
     &  (iad, jn, js, adx, qqv, ua, &
     &   ilo, ihi, julo, jhi, ju1_gl, j2_gl, j1p, j2p, &
     &   i1, i2, ju1, j2, k1, k2)
!
      use GmiCalcDivergence_mod, only : Deter_Jrange
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: ju1_gl, j2_gl, j1p, j2p
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer :: iad
      integer :: jn (k1:k2)
      integer :: js (k1:k2)
      real*8, intent(out) :: adx(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: qqv(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: ua (ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik, iu
      integer :: jst, jend
!
      real*8  :: a1, b1, c1
      real*8  :: rdiff
      real*8  :: ril, riu
      real*8  :: ru
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      call Deter_Jrange (ju1, j2, j1p, j2p, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!
!     =============
      if (iad == 1) then
!     =============
!
!       ----------
!       1st order.
!       ----------
!
        do ik = k1, k2
          do ij = jst, jend
!
            if ((ij <= js(ik)) .or. (ij >= jn(ik))) then
!
!             --------------
!             In Polar area.
!             --------------
!
              do il = i1, i2
!
                iu  = ua(il,ij,ik)
                riu = iu
                ru  = ua(il,ij,ik) - riu
                iu  = il - iu
!
                if (ua(il,ij,ik) >= 0.0d0) then
                  rdiff = qqv(iu-1,ij,ik) - qqv(iu,ij,ik)
                else
                  rdiff = qqv(iu,ij,ik)   - qqv(iu+1,ij,ik)
                end if
!
                adx(il,ij,ik) = (qqv(iu,ij,ik) - qqv(il,ij,ik)) +  &
     &                          (ru * rdiff)
!
              end do
!
            else  ! js < ij < jn
!
!             ----------------
!             Eulerian upwind.
!             ----------------
!
              do il = i1, i2
!
                ril = il
                iu  = ril - ua(il,ij,ik)
!
                adx(il,ij,ik) = ua(il,ij,ik) *  &
     &                          (qqv(iu,ij,ik) - qqv(iu+1,ij,ik))
!
              end do
!
            end if
!
          end do
        end do
!
!
!     ==================
      else if (iad == 2) then
!     ==================
!
        do ik = k1, k2
          do ij = jst, jend
!
            if ((ij <= js(ik)) .or. (ij >= jn(ik))) then
!
!             --------------
!             In Polar area.
!             --------------
!
              do il = i1, i2
!
                iu  = Nint (ua(il,ij,ik))
                riu = iu
                ru  = riu - ua(il,ij,ik)
                iu  = il - iu
!
                a1 = 0.5d0 * (qqv(iu+1,ij,ik) + qqv(iu-1,ij,ik)) -  &
     &               qqv(iu,ij,ik)
!
                b1 = 0.5d0 * (qqv(iu+1,ij,ik) - qqv(iu-1,ij,ik))
!
                c1 = qqv(iu,ij,ik) - qqv(il,ij,ik)
!
                adx(il,ij,ik) = (ru * ((a1 * ru) + b1)) + c1
!
              end do
!
            else  ! js < ij < jn
!
!             ----------------
!             Eulerian upwind.
!             ----------------
!
              do il = i1, i2
!
                iu  = Nint (ua(il,ij,ik))
                riu = iu
                ru  = riu - ua(il,ij,ik)
                iu  = il - iu
!
                a1 = 0.5d0 * (qqv(iu+1,ij,ik) + qqv(iu-1,ij,ik)) -  &
     &               qqv(iu,ij,ik)
!
                b1 = 0.5d0 * (qqv(iu+1,ij,ik) - qqv(iu-1,ij,ik))
!
                c1 = qqv(iu,ij,ik) - qqv(il,ij,ik)
!
                adx(il,ij,ik) = (ru * ((a1 * ru) + b1)) + c1
!
              end do
!
            end if
!
          end do
        end do
!
!     ======
      end if
!     ======
!
!
      if (ju1 == ju1_gl) then
!
        adx(i1:i2,ju1,:) = 0.0d0
!
        if (j1p /= ju1_gl+1) then
!
          adx(i1:i2,ju1+1,:) = 0.0d0
!
        end if
!
      end if
!
!
      if (j2 == j2_gl) then
!
        adx(i1:i2,j2,:) = 0.0d0
!
        if (j1p /= ju1_gl+1) then
!
          adx(i1:i2,j2-1,:) = 0.0d0
!
        end if
!
      end if
!
!
      return
!
      end
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Yadv_Dao2
!
! DESCRIPTION
!   This routine is the advective form N-S operator for computing the ady
!   (N-S) cross term.
!
! ARGUMENTS
!   iad : if iad = 1, use 1st order accurate scheme;
!         if iad = 2, use 2nd order accurate scheme
!   ady : cross term due to N-S advection (mixing ratio)
!   qqu : concentration contribution from E-W advection (mixing ratio)
!   va  : average of Courant numbers from ij and ij+1
!
!-----------------------------------------------------------------------------
!
      subroutine Yadv_Dao2  &
     &  (iad, ady, qqu, va, &
     &   numLonDomains, gmi_nborder, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   j1p, j2p, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole, numDomains, mapi_all, ivert)
!
      use GmiCalcDivergence_mod, only : Deter_Jrange
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numLonDomains, gmi_nborder, numDomains
      integer, intent(in) :: commu_npole, commu_spole
      integer, intent(in) :: j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ivert
      integer, intent(in) :: mapi_all(2,numDomains)
      integer, intent(in) :: iad
      real*8, intent(out) :: ady(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: qqu(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: va (ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik
      integer :: jst, jend
      integer :: jv
!
      real*8  :: a1, b1, c1
      real*8  :: rij, rjv
      real*8  :: rv
!
      real*8  :: qquwk(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
!      do ij = julo, jhi
        qquwk(:,:,:) = qqu(:,:,:)
!      end do
!
!
!     ======================
      call Do_Yadv_Pole_I2d2 &
!     ======================
     &     (qqu, qquwk, &
     &   numLonDomains, gmi_nborder, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &    mapi_all, commu_npole, commu_spole, numDomains, ivert)
!
!
      call Deter_Jrange (ju1-1, j2+1, j1p-1, j2p+1, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!
!     =============
      if (iad == 1) then
!     =============
!
!       ----------
!       1st order.
!       ----------
!
        do ik = k1, k2
          do ij = jst, jend
            do il = i1, i2
!c?
              rij = ij
              jv  = rij - va(il,ij,ik)
!
              ady(il,ij,ik) = va(il,ij,ik) *  &
     &                        (qquwk(il,jv,ik) - qquwk(il,jv+1,ik))
!
            end do
          end do
        end do
!
!
!     ==================
      else if (iad == 2) then
!     ==================
!
        do ik = k1, k2
          do ij = jst, jend
            do il = i1, i2
!c?
              jv  = Nint (va(il,ij,ik))
              rjv = jv
              rv  = rjv - va(il,ij,ik)
              jv  = ij - jv
!
              a1 = 0.5d0 * (qquwk(il,jv+1,ik) + qquwk(il,jv-1,ik)) -  &
     &             qquwk(il,jv,ik)
!
              b1 = 0.5d0 * (qquwk(il,jv+1,ik) - qquwk(il,jv-1,ik))
!
              c1 = qquwk(il,jv,ik) - qquwk(il,ij,ik)
!
              ady(il,ij,ik) = (rv * ((a1 * rv) + b1)) + c1
!
            end do
          end do
        end do
!
      end if
!
!
!     =====================
      call Do_Yadv_Pole_Sum &
!     =====================
     &   ( ady, &
     &   numLonDomains, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole)
!
      return
!
      end
!
!
