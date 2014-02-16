!
!=============================================================================
!
! $Id: xyz_ppm_dao2.F90,v 1.12 2013-07-31 14:45:57 ssteenro Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Shian-Jiann Lin, DAO)
!   jrt@llnl.gov
!
! FILE
!   xyz_ppm_dao2.F
!
! ROUTINES
!   Xtp
!   Ytp
!   Xmist
!   Ymist
!   Fxppm
!   Fyppm
!   Fzppm
!   Lmtppm
!
!=============================================================================
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Xtp
!
! DESCRIPTION
!   This routine does horizontal advection in the E-W direction.
!
! ARGUMENTS
!   ilmt  : controls various options in E-W advection
!   jn    : northward of latitude index = jn, Courant numbers could be > 1,
!           so use the flux-form semi-Lagrangian scheme
!   js    : southward of latitude index = js, Courant numbers could be > 1,
!           so use the flux-form semi-Lagrangian scheme
!   pu    : pressure at edges in "u" (mb)
!   crx   : Courant number in E-W direction
!   dq1   : species density (mb)
!   qqv   : concentration contribution from N-S advection (mixing ratio)
!   xmass : horizontal mass flux in E-W direction (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Xtp  &
     &  (ilmt, jn, js, pu, crx, dq1, qqv, xmass, fx, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   j1p, j2p, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &   i1, i2, ju1, j2, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld)
!
      use GmiSubDomainsBC_mod, only : Gmi_Bc3du_For_Sub
      use GmiCalcDivergence_mod, only : Deter_Jrange
!
      implicit none
!
#     include "gem_msg_numbers.h"
#     include "tpcore_constants_dao2.h"
!
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numDomains, numLonDomains, numLatDomains, gmi_nborder
      integer, intent(in) :: i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: j1p, j2p
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ilmt
      integer, intent(in) :: communicatorWorld
      integer, intent(in) :: nb, sb, eb, wb
      integer, intent(in) :: ewflag(ilo:ihi)
      integer, intent(in) :: nspoleu(julo:jhi)
      integer, intent(in) :: north_domain, south_domain, east_domain, west_domain
      integer, intent(in) :: jn(k1:k2)
      integer, intent(in) :: js(k1:k2)
      real*8, intent(in)  :: pu   (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: crx  (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(inout) :: dq1  (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(inout) :: qqv  (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: xmass(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: fx (ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik, ic
      integer :: iu, ix
      integer :: jst, jend
      integer :: jvan
!
      integer :: isav(i1:i2)
!
      real*8  :: rc
      real*8  :: ric, ril
!
      real*8  :: dcx(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      dcx(:,:,:) = 0.0d0
      fx(:,:,:) = 0.0d0
!
!
      call Deter_Jrange (ju1, j2, j1p, j2p, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!
      if (IORD /= 1) then
!       ==========
        call Xmist &
!       ==========
     &   (dcx, qqv, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   j1p, j2p, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &   i1, i2, ju1, j2, k1, k2, ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld)
      end if
!
!
      jvan = Max (1, j2_gl / 18)
!
!
!     ==============
      do ik = k1, k2
        do ij = jst, jend
!       =================
!
!         ======================================
          if ((ij > js(ik)) .and. (ij < jn(ik))) then
!         ======================================
!
!           ------------------------------------------------------
!           Do horizontal Eulerian advection in the E-W direction.
!           ------------------------------------------------------
!
            if ((IORD == 1) .or.  &
     &          (ij == j1p) .or. (ij == j2p)) then
!
              do il = i1, i2
                ril = il
                iu  = ril - crx(il,ij,ik)
!
                fx(il,ij,ik) = qqv(iu,ij,ik)
              end do
!
            else
!
              if ((IORD == 2) .or.  &
     &            (ij <= (j1p+jvan)) .or. (ij >= (j2p-jvan))) then
!
                do il = i1, i2
                  ril = il
                  iu  = ril - crx(il,ij,ik)
!
                  fx(il,ij,ik) =  &
     &              qqv(iu,ij,ik) +  &
     &              (dcx(iu,ij,ik) *  &
     &               (Sign (1.0d0, crx(il,ij,ik)) - crx(il,ij,ik)))
                end do
!
              else
!
!               ==========
!DIR$           INLINE
                call Fxppm  &
     &            (ij, ik, ilmt, crx, dcx, fx, qqv, &
     &             ilo, ihi, julo, jhi, i1, i2, k1, k2)
!  qqv (inout) - can be updated
!DIR$           NOINLINE
!               ==========
!
              end if
!
            end if
!
            fx(i1:i2,ij,ik) = fx(i1:i2,ij,ik) * xmass(i1:i2,ij,ik)
!
!         ====
          else
!         ====
!
!           ------------------------------------------------------------
!           Do horizontal Conservative (flux-form) Semi-Lagrangian
!           advection in the E-W direction (van Leer at high latitudes).
!           ------------------------------------------------------------
!
            if ((IORD == 1) .or.  &
     &          (ij == j1p) .or. (ij == j2p)) then
!
              do il = i1, i2
                ic       = crx(il,ij,ik)
                isav(il) = il - ic
                ril      = il
                iu       = ril - crx(il,ij,ik)
                ric      = ic
                rc       = crx(il,ij,ik) - ric
!
                fx(il,ij,ik) = rc * qqv(iu,ij,ik)
              end do
!
            else
!
              do il = i1, i2
                ic       = crx(il,ij,ik)
                isav(il) = il - ic
                ril      = il
                iu       = ril - crx(il,ij,ik)
                ric      = ic
                rc       = crx(il,ij,ik) - ric
!
                fx(il,ij,ik) =  &
     &            rc *  &
     &            (qqv(iu,ij,ik) +  &
     &             (dcx(iu,ij,ik) * (Sign (1.0d0, rc) - rc)))
              end do
!
            end if
!
            do il = i1, i2
!
              if (crx(il,ij,ik) > 1.0d0) then
!
!#if (ARCH_OPTION == ARCH_CRAY)
!!               ========
!!DIR$           NOVECTOR
!!               ========
!#endif
!
                do ix = isav(il), il - 1
                  fx(il,ij,ik) = fx(il,ij,ik) + qqv(ix,ij,ik)
                end do
!
              else if (crx(il,ij,ik) < -1.0d0) then
!
                do ix = il, isav(il) - 1
                  fx(il,ij,ik) = fx(il,ij,ik) - qqv(ix,ij,ik)
                end do
!
!#if (ARCH_OPTION == ARCH_CRAY)
!!               ======
!!DIR$           VECTOR
!!               ======
!#endif
!
              end if
!
            end do
!
            fx(i1:i2,ij,ik) = pu(i1:i2,ij,ik) * fx(i1:i2,ij,ik)
!
!         ======
          end if
!         ======
!
!       ======
        end do
      end do
!     ======
!
!
!     ======================
      call Gmi_Bc3du_For_Sub &
!     ======================
     &  (fx, ACTM_BC_FX, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld, &
     &       numDomains)
!
!
      do ij = jst, jend
        do il = i1, i2
!
          dq1(il,ij,:) = dq1(il,ij,:) +  &
     &                   (fx(il,ij,:) - fx(il+1,ij,:))
!
        end do
      end do
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
!   Ytp
!
! DESCRIPTION
!   This routine does horizontal advection in the N-S direction.
!
! ARGUMENTS
!   jlmt      : controls various options in N-S advection
!   geofac_pc : special geometrical factor (geofac) for Polar cap
!   geofac    : geometrical factor for meridional advection; geofac uses
!               correct spherical geometry, and replaces acosp as the
!               meridional geometrical factor in tpcore
!   cry       : Courant number in N-S direction
!   dq1       : species density (mb)
!   qqu       : concentration contribution from E-W advection (mixing ratio)
!   qqv       : concentration contribution from N-S advection (mixing ratio)
!   ymass     : horizontal mass flux in N-S direction (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Ytp  &
     &  (jlmt, geofac_pc, geofac, cry, dq1, qqu, qqv, ymass, fy, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ivert, ilong, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld, &
     & commu_npole, commu_spole, mapi_all)
!
      use GmiCalcDivergence_mod, only : Deter_Jrange
!
      implicit none
!
#     include "tpcore_constants_dao2.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numDomains, numLonDomains, numLatDomains
      integer, intent(in) :: gmi_nborder
      integer, intent(in) :: j1p, j2p, ivert, ilong
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: jlmt
      integer, intent(in) :: mapi_all(2,numDomains)
      integer, intent(in) :: communicatorWorld
      integer, intent(in) :: nb, sb, eb, wb, commu_npole, commu_spole
      integer, intent(in) :: ewflag(ilo:ihi)
      integer, intent(in) :: nspoleu(julo:jhi)
      integer, intent(in) :: north_domain, south_domain, east_domain, west_domain
      real*8, intent(in)  :: geofac_pc
      real*8, intent(in)  :: geofac(ju1_gl:j2_gl)
      real*8, intent(in)  :: cry   (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(inout) :: dq1   (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: qqu   (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(inout)  :: qqv   (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: ymass (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: fy (ilo:ihi, julo:jhi, k1:k2)
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik
      integer :: jst, jend
      integer :: jv
!
      real*8  :: rj1p
!
      real*8  :: dcy(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      dcy(:,:,:) = 0.0d0
      fy(:,:,:) = 0.0d0
!
      rj1p = j1p
!
      call Deter_Jrange (ju1, j2+1, j1p, j2p+1, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!     ==============
      if (JORD == 1) then
!     ==============
!
        do ik = k1, k2
          do ij = jst, jend
            do il = i1, i2
!c?
              jv = rj1p - cry(il,ij,ik)
!
              qqv(il,ij,ik) = qqu(il,jv,ik)
!
            end do
          end do
        end do
!
!     ====
      else
!     ====
!
!       ==========
        call Ymist  &
!       ==========
     &    (4, dcy, qqu, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   i1_gl, i2_gl, ju1_gl, j2_gl, j1p, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld, &
     &       mapi_all, commu_npole, commu_spole, ivert)
!
!
        if ((JORD <= 0) .or. (JORD >= 3)) then
!
!         ==========
          call Fyppm  &
!         ==========
     &      (jlmt, cry, dcy, qqu, qqv, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ilong, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld, &
     &       mapi_all, commu_npole, commu_spole, ivert)
!
        else
!
          do ik = k1, k2
            do ij = jst, jend
              do il = i1, i2
!c?
                jv = rj1p - cry(il,ij,ik)
!
                qqv(il,ij,ik) =  &
     &            qqu(il,jv,ik) +  &
     &            ((Sign (1.0d0, cry(il,ij,ik)) - cry(il,ij,ik)) *  &
     &             dcy(il,jv,ik))
!
              end do
            end do
          end do
!
        end if
!
      end if
!
!
      do ij = jst, jend
        qqv(i1:i2,ij,:) = qqv(i1:i2,ij,:) * ymass(i1:i2,ij,:)
      end do
!
!
      call Deter_Jrange (ju1, j2, j1p, j2p, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!.sds.. save N-S species flux as diagnostic
      do ik = k1,k2
        do ij = i1,i2
          fy(ij,jst:jend+1,ik) = qqv(ij,jst:jend+1,ik) * geofac(jst:jend+1)
        enddo
      enddo
!
!... meridional flux update
      do ij = jst, jend
!
        dq1(i1:i2,ij,:) =  &
     &    dq1(i1:i2,ij,:) +  &
     &    (qqv(i1:i2,ij,:) - qqv(i1:i2,ij+1,:)) * geofac(ij)
!
      end do
!
!     ====================
      call Do_Ytp_Pole_Sum  &
!     ====================
     &  (geofac_pc, dq1, qqv, fy, &
     &   numLonDomains, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, j2p, ivert, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     & commu_npole, commu_spole)
!
      return
!
      end
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Xmist
!
! DESCRIPTION
!   This routine computes the linear tracer slope n the E-W direction.
!   It uses the Lin et. al. 1994 algorithm.
!
! ARGUMENTS
!   dcx : slope of concentration distribution in E-W direction (mixing ratio)
!   qqv : concentration contribution from N-S advection (mixing ratio)
!
!-----------------------------------------------------------------------------
!
      subroutine Xmist  &
     &  (dcx, qqv, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   j1p, j2p, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &   i1, i2, ju1, j2, k1, k2, ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld)
!
      use GmiSubDomainsBC_mod, only : Gmi_Bc3du_For_Sub
      use GmiCalcDivergence_mod, only : Deter_Jrange
!
      implicit none
!
#     include "gem_msg_numbers.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numDomains, numLonDomains, numLatDomains, gmi_nborder
      integer, intent(in) :: i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: j1p, j2p
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: communicatorWorld
      integer, intent(in) :: nb, sb, eb, wb
      integer, intent(in) :: ewflag(ilo:ihi)
      integer, intent(in) :: nspoleu(julo:jhi)
      integer, intent(in) :: north_domain, south_domain, east_domain, west_domain
      real*8, intent(out) :: dcx(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in) :: qqv(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik
      integer :: jst, jend
!
      real*8  :: pmax, pmin
      real*8  :: r24
      real*8  :: tmp
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      r24 = 1.0d0 / 24.0d0
!
!
      call Deter_Jrange (ju1, j2, j1p+1, j2p-1, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!
      do ik = k1, k2
        do ij = jst, jend
          do il = i1, i2
!
            tmp =  &
     &        ((8.0d0 * (qqv(il+1,ij,ik) - qqv(il-1,ij,ik))) +  &
     &         qqv(il-2,ij,ik) - qqv(il+2,ij,ik)) *  &
     &        r24
!
            pmax =  &
     &        Max (qqv(il-1,ij,ik), qqv(il,ij,ik), qqv(il+1,ij,ik)) -  &
     &        qqv(il,ij,ik)
!
            pmin =  &
     &        qqv(il,ij,ik) -  &
     &        Min (qqv(il-1,ij,ik), qqv(il,ij,ik), qqv(il+1,ij,ik))
!
            dcx(il,ij,ik) = Sign (Min (Abs (tmp), pmax, pmin), tmp)
!
          end do
        end do
      end do
!
!
!     ======================
      call Gmi_Bc3du_For_Sub &
!     ======================
     &    (dcx, ACTM_BC_DCX, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld, &
     &       numDomains)
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
!   Ymist
!
! DESCRIPTION
!   This routine computes the linear tracer slope n the N-S direction.
!   It uses the Lin et. al. 1994 algorithm.
!
! ARGUMENTS
!   id  : the "order" of the accuracy in the computed linear "slope"
!         (or mismatch, Lin et al. 1994); it is either 2 or 4.
!   dcy : slope of concentration distribution in N-S direction (mixing ratio)
!   qqu : concentration contribution from E-W advection (mixing ratio)
!
!-----------------------------------------------------------------------------
!
      subroutine Ymist  &
     &  (id, dcy, qqu, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   i1_gl, i2_gl, ju1_gl, j2_gl, j1p, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   ewflag, nspoleu, nb, sb, eb, wb,  &
     &   north_domain, south_domain, east_domain, west_domain,  &
     &   communicatorWorld, &
     &   mapi_all, commu_npole, commu_spole, ivert)
!
      use GmiSubDomainsBC_mod, only : Gmi_Bc3du_For_Sub
!
      implicit none
!
#     include "gem_msg_numbers.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numDomains, numLonDomains, numLatDomains
      integer, intent(in) :: gmi_nborder
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, j1p
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ivert
      integer, intent(in) :: id
      integer, intent(in) :: mapi_all(2,numDomains)
      integer, intent(in) :: communicatorWorld
      integer, intent(in) :: commu_npole, commu_spole
      integer, intent(in) :: nb, sb, eb, wb
      integer, intent(in) :: ewflag(ilo:ihi)
      integer, intent(in) :: nspoleu(julo:jhi)
      integer, intent(in) :: north_domain, south_domain, east_domain, west_domain
      real*8, intent(out) :: dcy(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in) :: qqu(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik
!
      real*8  :: pmax, pmin
      real*8  :: r24
      real*8  :: tmp
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      r24  = 1.0d0 / 24.0d0
!
!
!     ============
      if (id == 2) then
!     ============
!
        do ik = k1, k2
          do ij = ju1 - 1, j2 - 1
            do il = i1, i2
!
              tmp  = 0.25d0 * (qqu(il,ij+2,ik) - qqu(il,ij,ik))
!
              pmax =  &
     &          Max (qqu(il,ij,ik), qqu(il,ij+1,ik), qqu(il,ij+2,ik)) -  &
     &          qqu(il,ij+1,ik)
!
              pmin =  &
     &          qqu(il,ij+1,ik) -  &
     &          Min (qqu(il,ij,ik), qqu(il,ij+1,ik), qqu(il,ij+2,ik))
!
              dcy(il,ij+1,ik) = Sign (Min (Abs (tmp), pmin, pmax), tmp)
!
            end do
          end do
        end do
!
!     ====
      else
!     ====
!
!       ========================
        call Do_Ymist_Pole1_I2d2 &
!       ========================
     &     (dcy, qqu, &
     &   numLonDomains, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   mapi_all, commu_npole, commu_spole, numDomains, ivert)
!
        do ik = k1, k2
          do ij = ju1 - 2, j2 - 2
            do il = i1, i2
!
              tmp  = ((8.0d0 * (qqu(il,ij+3,ik) - qqu(il,ij+1,ik))) +  &
     &                qqu(il,ij,ik) - qqu(il,ij+4,ik)) *  &
     &               r24
!
              pmax =  &
     &          Max (qqu(il,ij+1,ik), qqu(il,ij+2,ik), qqu(il,ij+3,ik))  &
     &           - qqu(il,ij+2,ik)
!
              pmin =  &
     &          qqu(il,ij+2,ik) -  &
     &          Min (qqu(il,ij+1,ik), qqu(il,ij+2,ik), qqu(il,ij+3,ik))
!
              dcy(il,ij+2,ik) = Sign (Min (Abs (tmp), pmin, pmax), tmp)
!
            end do
          end do
        end do
!
      end if
!
!
!     ========================
      call Do_Ymist_Pole2_I2d2 &
!     ========================
     &    (dcy, qqu, &
     &   numLonDomains, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   mapi_all, commu_npole, commu_spole, numDomains, ivert)
!
!
!     ======================
      call Gmi_Bc3du_For_Sub &
!     ======================
     &  (dcy, ACTM_BC_DCY, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld, &
     &       numDomains)
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
!   Fxppm
!
! DESCRIPTION
!   This routine is the 1D "outer" flux form operator based on the Piecewise
!   Parabolic Method (PPM; see also Lin and Rood 1996) for computing the
!   fluxes in the E-W direction.
!
! ARGUMENTS
!   ij   : latitude
!   ik   : altitude
!   ilmt : controls various options in E-W advection
!   crx  : Courant number in E-W direction
!   dcx  : slope of concentration distribution in E-W direction (mixing ratio)
!   fx   : E-W flux (mixing ratio)
!   qqv  : concentration contribution from N-S advection (mixing ratio)
!
!-----------------------------------------------------------------------------
!
      subroutine Fxppm  &
     &  (ij, ik, ilmt, crx, dcx, fx, qqv, &
     &   ilo, ihi, julo, jhi, i1, i2, k1, k2)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: ilo, ihi, julo, jhi, i1, i2, k1, k2
      integer :: ij, ik
      integer, intent(in) :: ilmt
      real*8, intent(in)  :: crx(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: dcx(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: fx (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(inout) :: qqv(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      logical, save :: first = .true.
!
      integer :: il
      integer :: ilm1
      integer :: lenx
!
      real*8  :: r13, r23
      real*8  :: rval
!
      real*8, allocatable, save :: a6(:)
      real*8, allocatable, save :: al(:)
      real*8, allocatable, save :: ar(:)
!
      real*8, allocatable, save :: a61(:)
      real*8, allocatable, save :: al1(:)
      real*8, allocatable, save :: ar1(:)
!
      real*8, allocatable, save :: dcxi1(:)
      real*8, allocatable, save :: qqvi1(:)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
!     ==========
      if (first) then
!     ==========
!
        first = .false.
!
        Allocate (a6(ilo:ihi))
        Allocate (al(ilo:ihi))
        Allocate (ar(ilo:ihi))
        a6 = 0.0d0; al = 0.0d0; ar = 0.0d0
!
        Allocate (a61((ihi-1)-(ilo+1)+1))
        Allocate (al1((ihi-1)-(ilo+1)+1))
        Allocate (ar1((ihi-1)-(ilo+1)+1))
        a61 = 0.0d0; al1 = 0.0d0; ar1 = 0.0d0
!
        Allocate (dcxi1((ihi-1)-(ilo+1)+1))
        Allocate (qqvi1((ihi-1)-(ilo+1)+1))
        dcxi1 = 0.0d0; qqvi1 = 0.0d0
!
      end if
!
!
      r13 = 1.0d0 / 3.0d0
      r23 = 2.0d0 / 3.0d0
!
!
      do il = ilo + 1, ihi
!
        rval = 0.5d0 * (qqv(il-1,ij,ik) + qqv(il,ij,ik)) +  &
     &         (dcx(il-1,ij,ik) - dcx(il,ij,ik)) * r13
!
        al(il)   = rval
        ar(il-1) = rval
!
      end do
!
!
      do il = ilo + 1, ihi - 1
        a6(il) = 3.0d0 *  &
     &           (qqv(il,ij,ik) + qqv(il,ij,ik) - (al(il) + ar(il)))
      end do
!
!
!     ==============
      if (ilmt <= 2) then
!     ==============
!
        a61(:) = 0.0d0
        al1(:) = 0.0d0
        ar1(:) = 0.0d0
!
        dcxi1(:) = 0.0d0
        qqvi1(:) = 0.0d0
!
        lenx = 0
!
        do il = ilo + 1, ihi - 1
!
          lenx = lenx + 1
!
          a61(lenx)   = a6(il)
          al1(lenx)   = al(il)
          ar1(lenx)   = ar(il)
!
          dcxi1(lenx) = dcx(il,ij,ik)
          qqvi1(lenx) = qqv(il,ij,ik)
!
        end do
!
!       ===========
!DIR$   INLINE
        call Lmtppm  &
     &    (lenx, ilmt, a61, al1, ar1, dcxi1, qqvi1)
!DIR$   NOINLINE
!       ===========
!
        lenx = 0
!
        do il = ilo + 1, ihi - 1
!
          lenx = lenx + 1
!
          a6(il)   = a61(lenx)
          al(il)   = al1(lenx)
          ar(il)   = ar1(lenx)
!
          dcx(il,ij,ik) = dcxi1(lenx)
          qqv(il,ij,ik) = qqvi1(lenx)
!
        end do
!
      end if
!
!
      do il = i1, i2
!
        if (crx(il,ij,ik) > 0.0d0) then
!
          ilm1 = il - 1
!
          fx(il,ij,ik) =  &
     &      ar(ilm1) +  &
     &      0.5d0 * crx(il,ij,ik) *  &
     &      (al(ilm1) - ar(ilm1) +  &
     &       (a6(ilm1) * (1.0d0 - (r23 * crx(il,ij,ik)))))
!
        else
!
          fx(il,ij,ik) =  &
     &      al(il) -  &
     &      0.5d0 * crx(il,ij,ik) *  &
     &      (ar(il) - al(il) +  &
     &       (a6(il) * (1.0d0 + (r23 * crx(il,ij,ik)))))
!
        end if
!
      end do
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
!   Fyppm
!
! DESCRIPTION
!   This routine is the 1D "outer" flux form operator based on the Piecewise
!   Parabolic Method (PPM; see also Lin and Rood 1996) for computing the
!   fluxes in the N-S direction.
!
! ARGUMENTS
!   jlmt : controls various options in N-S advection
!   cry  : Courant number in N-S direction
!   dcy  : slope of concentration distribution in N-S direction (mixing ratio)
!   qqu  : concentration contribution from E-W advection (mixing ratio)
!   qqv  : concentration contribution from N-S advection (mixing ratio)
!
!-----------------------------------------------------------------------------
!
      subroutine Fyppm  &
     &  (jlmt, cry, dcy, qqu, qqv, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ilong, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld,&
     &       mapi_all, commu_npole, commu_spole, ivert)
!
      use GmiSubDomainsBC_mod, only : Gmi_Bc3du_For_Sub
      use GmiCalcDivergence_mod, only : Deter_Jrange
!
      implicit none
!
#     include "gem_msg_numbers.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numDomains, numLonDomains, numLatDomains
      integer, intent(in) :: gmi_nborder
      integer, intent(in) :: j1p, j2p, ilong
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ivert
      integer, intent(in) :: mapi_all(2,numDomains)
      integer, intent(in) :: communicatorWorld
      integer, intent(in) :: commu_npole, commu_spole
      integer, intent(in) :: nb, sb, eb, wb
      integer, intent(in) :: ewflag(ilo:ihi)
      integer, intent(in) :: nspoleu(julo:jhi)
      integer, intent(in) :: north_domain, south_domain, east_domain, west_domain
      integer, intent(in) :: jlmt
      real*8, intent(in)  :: cry(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: dcy(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: qqu(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: qqv(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: ijm1
      integer :: il, ij, ik
      integer :: jst, jend
      integer :: lenx
!
      real*8  :: r13, r23
!
      real*8  :: a61 (ilong*((jhi-1)-(julo+1)+1))
      real*8  :: al1 (ilong*((jhi-1)-(julo+1)+1))
      real*8  :: ar1 (ilong*((jhi-1)-(julo+1)+1))
      real*8  :: dcy1(ilong*((jhi-1)-(julo+1)+1))
      real*8  :: qqu1(ilong*((jhi-1)-(julo+1)+1))
!
      real*8  :: a6(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: al(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: ar(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      a6(:,:,:) = 0.0d0
      al(:,:,:) = 0.0d0
      ar(:,:,:) = 0.0d0
!
!
      r13 = 1.0d0 / 3.0d0
      r23 = 2.0d0 / 3.0d0
!
!
      do ij = julo + 1, jhi
        al(i1:i2,ij,:) =  &
     &    0.5d0 * (qqu(i1:i2,ij-1,:) + qqu(i1:i2,ij,:)) +  &
     &    (dcy(i1:i2,ij-1,:) - dcy(i1:i2,ij,:)) * r13
      end do
!
!
      do ij = julo, jhi - 1
        ar(i1:i2,ij,:) = al(i1:i2,ij+1,:)
      end do
!
!
!     =======================
      call Do_Fyppm_Pole_I2d2 &
!     =======================
     &   (al, ar, &
     &   numLonDomains, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   mapi_all, commu_npole, commu_spole, numDomains, ivert)
!
!
      do ij = julo + 1, jhi - 1
!
        a6(i1:i2,ij,:) =  &
     &    3.0d0 *  &
     &    (qqu(i1:i2,ij,:) + qqu(i1:i2,ij,:) -  &
     &     (al(i1:i2,ij,:) + ar(i1:i2,ij,:)))
!
      end do
!
!
!
!     ==============
      if (jlmt <= 2) then
!     ==============
!
        do ik = k1, k2
!
          lenx = 0
!
          do ij = julo + 1, jhi - 1
            do il = i1, i2
!
              lenx = lenx + 1
!
              a61 (lenx) = a6 (il,ij,ik)
              al1 (lenx) = al (il,ij,ik)
              ar1 (lenx) = ar (il,ij,ik)
              dcy1(lenx) = dcy(il,ij,ik)
              qqu1(lenx) = qqu(il,ij,ik)
!
            end do
          end do
!
!         ===========
!DIR$     INLINE
          call Lmtppm  &
     &      (lenx, jlmt, a61, al1, ar1, dcy1, qqu1)
!DIR$     NOINLINE
!         ===========
!
          lenx = 0
!
          do ij = julo + 1, jhi - 1
            do il = i1, i2
!
              lenx = lenx + 1
!
              a6(il,ij,ik) = a61(lenx)
              al(il,ij,ik) = al1(lenx)
              ar(il,ij,ik) = ar1(lenx)
!
            end do
          end do
!
        end do
!
      end if
!
!
      call Deter_Jrange (ju1, j2, j1p, j2p+1, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
      do ik = k1, k2
        do ij = jst, jend
!
          ijm1 = ij - 1
!
          do il = i1, i2
!
            if (cry(il,ij,ik) > 0.0d0) then
!
              qqv(il,ij,ik) =  &
     &          ar(il,ijm1,ik) +  &
     &          0.5d0 * cry(il,ij,ik) *  &
     &          (al(il,ijm1,ik) - ar(il,ijm1,ik) +  &
     &           (a6(il,ijm1,ik) * (1.0d0 - (r23 * cry(il,ij,ik)))))
!
            else
!
              qqv(il,ij,ik) =  &
     &          al(il,ij,ik) -  &
     &          0.5d0 * cry(il,ij,ik) *  &
     &          (ar(il,ij,ik) - al(il,ij,ik) +  &
     &           (a6(il,ij,ik) * (1.0d0 + (r23 * cry(il,ij,ik)))))
!
            end if
!
          end do
!
        end do
      end do
!
!
!     ======================
      call Gmi_Bc3du_For_Sub &
!     ======================
     &   (qqv, ACTM_BC_QQV, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld, &
     &       numDomains)
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
!   Fzppm
!
! DESCRIPTION
!   This routine is the 1D "outer" flux form operator based on the Piecewise
!   Parabolic Method (PPM; see also Lin and Rood 1996) for computing the
!   fluxes in the vertical direction.
!
!   Fzppm was modified by S.-J. Lin, 12/14/98, to allow the use of the KORD=7
!   (klmt=4) option.  KORD=7 enforces the 2nd monotonicity constraint of
!   Huynh (1996).  Note that in Huynh's original scheme, two constraints are
!   necessary for the preservation of monotonicity.  To use Huynh's
!   algorithm, it was modified as follows.  The original PPM is still used to
!   obtain the first guesses for the cell edges, and as such Huynh's 1st
!   constraint is no longer needed.  Huynh's median function is also replaced
!   by a simpler yet functionally equivalent in-line algorithm.
!
! ARGUMENTS
!   ispc  : species index
!   klmt  : controls various options in vertical advection
!   delp1 : pressure thickness, the pseudo-density in a hydrostatic system at
!           t1 (mb)
!c!?
!   dpi   : work array
!   wz    : large scale mass flux (per time step tdt) in the vertical
!           direction as diagnosed from the hydrostatic relationship (mb)
!   dq1   : species density (mb)
!   qq1   : species concentration (mixing ratio)
!
!-----------------------------------------------------------------------------
!
      subroutine Fzppm  &
     &  (ispc, klmt, delp1, dpi, wz, dq1, qq1, fz, &
     &   j1p, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &   ilong, ivert, i1, i2, ju1, j2, k1, k2)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in   ) :: ilong, ivert
      integer, intent(in   ) :: j1p, ju1_gl, j2_gl
      integer, intent(in   ) :: ilo, ihi, julo, jhi
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in)   :: ispc
      integer, intent(in)   :: klmt
      real*8, intent(in)    :: delp1(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)    :: wz  ( i1:i2,   ju1:j2,  k1:k2)
      real*8, intent(in)    :: qq1 (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(inout) :: dq1 (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out)   :: dpi ( i1:i2,   ju1:j2,  k1:k2)
      real*8, intent(out)   :: fz  (ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik
      integer :: k1p1, k1p2
      integer :: k2m1, k2m2
      integer :: lenx
!
      real*8  :: a1, a2
      real*8  :: aa, bb
      real*8  :: c0, c1, c2
      real*8  :: cm, cp
      real*8  :: fac1, fac2
      real*8  :: lac
      real*8  :: qmax, qmin
      real*8  :: qmp
      real*8  :: r13, r23
      real*8  :: tmp
!
      real*8  :: a61  (ilong*(ivert-4))
      real*8  :: al1  (ilong*(ivert-4))
      real*8  :: ar1  (ilong*(ivert-4))
      real*8  :: dca1 (ilong*(ivert-4))
      real*8  :: qq1a1(ilong*(ivert-4))
!
      real*8  :: a6   (i1:i2, k1:k2)
      real*8  :: al   (i1:i2, k1:k2)
      real*8  :: ar   (i1:i2, k1:k2)
      real*8  :: dca  (i1:i2, k1:k2)
      real*8  :: dlp1a(i1:i2, k1:k2)
      real*8  :: qq1a (i1:i2, k1:k2)
      real*8  :: wza  (i1:i2, k1:k2)
!
      real*8  :: dc   (i1:i2, ju1:j2, k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      a6(:,:) = 0.0d0
      al(:,:) = 0.0d0
      ar(:,:) = 0.0d0
      dc(:,:,:) = 0.0d0
!.sds... diagnostic vertical flux for species - set top to 0.0
      fz(:,:,:) = 0.0
!
!
      k1p1 = k1 + 1
      k1p2 = k1 + 2
!
      k2m1 = k2 - 1
      k2m2 = k2 - 2
!
      r13  = 1.0d0 / 3.0d0
      r23  = 2.0d0 / 3.0d0
!
!
!     -------------------
!     Compute dc for PPM.
!     -------------------
!
      do ik = k1, k2m1
        dpi(:,:,ik) = qq1(i1:i2,ju1:j2,ik+1) - qq1(i1:i2,ju1:j2,ik)
      end do
!
!
      do ik = k1p1, k2m1
        do ij = ju1, j2
          do il = i1, i2
!
            c0 = delp1(il,ij,ik) /  &
     &           (delp1(il,ij,ik-1) + delp1(il,ij,ik) +  &
     &            delp1(il,ij,ik+1))
!
            c1 = (delp1(il,ij,ik-1) + (0.5d0 * delp1(il,ij,ik))) /  &
     &           (delp1(il,ij,ik+1) + delp1(il,ij,ik))
!
            c2 = (delp1(il,ij,ik+1) + (0.5d0 * delp1(il,ij,ik))) /  &
     &           (delp1(il,ij,ik-1) + delp1(il,ij,ik))
!
            tmp = c0 *  &
     &            ((c1 * dpi(il,ij,ik)) +  &
     &             (c2 * dpi(il,ij,ik-1)))
!
            qmax =  &
     &        Max (qq1(il,ij,ik-1),  &
     &             qq1(il,ij,ik),  &
     &             qq1(il,ij,ik+1)) -  &
     &        qq1(il,ij,ik)
!
            qmin =  &
     &        qq1(il,ij,ik) -  &
     &        Min (qq1(il,ij,ik-1),  &
     &             qq1(il,ij,ik),  &
     &             qq1(il,ij,ik+1))
!
            dc(il,ij,ik) = Sign (Min (Abs (tmp), qmax, qmin), tmp)
!
          end do
        end do
      end do
!
!
!c?
!     -------------------------------------
!     Loop over latitudes (to save memory).
!     -------------------------------------
!
!     =======================
      ijloop: do ij = ju1, j2
!     =======================
!
        if (((ij == ju1_gl+1) .or. (ij == j2_gl-1)) .and.  &
     &      (j1p /= ju1_gl+1)) then
!         ============
          cycle ijloop
!         ============
        end if
!
!
        dca(:,:) = dc(:,ij,:)  ! the monotone slope
        wza(:,:) = wz(:,ij,:)
!
        dlp1a(:,:) = delp1(i1:i2,ij,:)
        qq1a (:,:) = qq1  (i1:i2,ij,:)
!
!
!       ----------------------------------------------------------------
!       Compute first guesses at cell interfaces.  First guesses are
!       required to be continuous.  Three-cell parabolic subgrid
!       distribution at model top; two-cell parabolic with zero gradient
!       subgrid distribution at the surface.
!       ----------------------------------------------------------------
!
!       ---------------------------
!       First guess top edge value.
!       ---------------------------
!
        do il = i1, i2
!
!         ------------------------------------------------------------
!         Three-cell PPM; compute a, b, & c of q = aP^2 + bP + c using
!         cell averages and dlp1a.
!         ------------------------------------------------------------
!
          fac1 = dpi(il,ij,k1p1) -  &
     &           dpi(il,ij,k1) * (dlp1a(il,k1p1) + dlp1a(il,k1p2)) /  &
     &           (dlp1a(il,k1) + dlp1a(il,k1p1))
!
          fac2 = (dlp1a(il,k1p1) + dlp1a(il,k1p2)) *  &
     &           (dlp1a(il,k1) + dlp1a(il,k1p1) + dlp1a(il,k1p2))
!
          aa = 3.0d0 * fac1 / fac2
!
          bb =  &
     &      2.0d0 * dpi(il,ij,k1) / (dlp1a(il,k1) + dlp1a(il,k1p1)) -  &
     &      r23 * aa * (2.0d0 * dlp1a(il,k1) + dlp1a(il,k1p1))
!
          al(il,k1) = qq1a(il,k1) -  &
     &                dlp1a(il,k1) *  &
     &                (r13 * aa * dlp1a(il,k1) +  &
     &                 0.5d0 * bb)
!
          al(il,k1p1) = dlp1a(il,k1) * (aa * dlp1a(il,k1) + bb) +  &
     &                  al(il,k1)
!
!         ---------------------
!         Check if change sign.
!         ---------------------
!
          if ((qq1a(il,k1) * al(il,k1)) <= 0.0d0) then
!
            al (il,k1) = 0.0d0
            dca(il,k1) = 0.0d0
!
          else
!
            dca(il,k1) = qq1a(il,k1) - al(il,k1)
!
          end if
!
        end do
!
!       -------
!       Bottom.
!       -------
!
        do il = i1, i2
!
!         ---------------------------------------------------
!         2-cell PPM with zero gradient right at the surface.
!         ---------------------------------------------------
!
!c        fac1 = dpi(il,ij,k2m1) * dlp1a(il,k2)**2 /
          fac1 = dpi(il,ij,k2m1) * (dlp1a(il,k2) * dlp1a(il,k2)) /  &
     &           ((dlp1a(il,k2) + dlp1a(il,k2m1)) *  &
     &            (2.0d0 * dlp1a(il,k2) + dlp1a(il,k2m1)))
!
          ar(il,k2) = qq1a(il,k2) + fac1
          al(il,k2) = qq1a(il,k2) - (fac1 + fac1)
!
          if ((qq1a(il,k2) * ar(il,k2)) <= 0.0d0) then
            ar(il,k2) = 0.0d0
          end if
!
          dca(il,k2) = ar(il,k2) - qq1a(il,k2)
!
        end do
!
!
!       ----------------------------------------
!       4th order interpolation in the interior.
!       ----------------------------------------
!
        do ik = k1p2, k2m1
          do il = i1, i2
!
            c1 = (dpi(il,ij,ik-1) * dlp1a(il,ik-1)) /  &
     &           (dlp1a(il,ik-1) + dlp1a(il,ik))
!
            c2 = 2.0d0 /  &
     &           (dlp1a(il,ik-2) + dlp1a(il,ik-1) +  &
     &            dlp1a(il,ik)   + dlp1a(il,ik+1))
!
            a1 = (dlp1a(il,ik-2) + dlp1a(il,ik-1)) /  &
     &           (2.0d0 * dlp1a(il,ik-1) + dlp1a(il,ik))
!
            a2 = (dlp1a(il,ik) + dlp1a(il,ik+1)) /  &
     &           (2.0d0 * dlp1a(il,ik) + dlp1a(il,ik-1))
!
            al(il,ik) =  &
     &        qq1a(il,ik-1) + c1 +  &
     &        c2 *  &
     &        (dlp1a(il,ik) * (c1 * (a1 - a2) + a2 * dca(il,ik-1)) -  &
     &         dlp1a(il,ik-1) * a1 * dca(il,ik))
!
          end do
        end do
!
!
        do ik = k1, k2m1
          ar(:,ik) = al(:,ik+1)
        end do
!
!
!       ---------------------------------------
!       Top & Bottom 2 layers always monotonic.
!       ---------------------------------------
!
        lenx = i2 - i1 + 1
!
        do ik = k1, k1p1
!
          do il = i1, i2
!
            a6(il,ik) =  &
     &        3.0d0 * (qq1a(il,ik) + qq1a(il,ik) -  &
     &                 (al(il,ik)  + ar(il,ik)))
          end do
!
!         ===========
!DIR$     INLINE
          call Lmtppm  &
     &      (lenx, 0, a6(i1,ik), al(i1,ik), ar(i1,ik),  &
     &       dca(i1,ik), qq1a(i1,ik))
!DIR$     NOINLINE
!         ===========
!
        end do
!
!
        do ik = k2m1, k2
!
          do il = i1, i2
!
            a6(il,ik) =  &
     &        3.0d0 * (qq1a(il,ik) + qq1a(il,ik) -  &
     &                 (al(il,ik)  + ar(il,ik)))
          end do
!
!         ===========
!DIR$     INLINE
          call Lmtppm  &
     &      (lenx, 0, a6(i1,ik), al(i1,ik), ar(i1,ik),  &
     &       dca(i1,ik), qq1a(i1,ik))
!DIR$     NOINLINE
!         ===========
!
        end do
!
!
!       ---------------------------
!       Interior depending on klmt.
!       ---------------------------
!
!       ==============
        if (klmt == 4) then
!       ==============
!
!         -------------------------------
!         KORD=7, Huynh's 2nd constraint.
!         -------------------------------
!
          do ik = k1p1, k2m1
            dca(:,ik) = dpi(:,ij,ik) - dpi(:,ij,ik-1)
          end do
!
!
          do ik = k1p2, k2m2
            do il = i1, i2
!
!             ------------
!             Right edges.
!             ------------
!
              qmp   = qq1a(il,ik) + (2.0d0 * dpi(il,ij,ik-1))
              lac   = qq1a(il,ik) +  &
     &                (1.5d0 * dca(il,ik-1)) + (0.5d0 * dpi(il,ij,ik-1))
              qmin  = Min (qq1a(il,ik), qmp, lac)
              qmax  = Max (qq1a(il,ik), qmp, lac)
!
              ar(il,ik) = Min (Max (ar(il,ik), qmin), qmax)
!
!             -----------
!             Left edges.
!             -----------
!
              qmp   = qq1a(il,ik) - (2.0d0 * dpi(il,ij,ik))
              lac   = qq1a(il,ik) +  &
     &                (1.5d0 * dca(il,ik+1)) - (0.5d0 * dpi(il,ij,ik))
              qmin  = Min (qq1a(il,ik), qmp, lac)
              qmax  = Max (qq1a(il,ik), qmp, lac)
!
              al(il,ik) = Min (Max (al(il,ik), qmin), qmax)
!
!             -------------
!             Recompute a6.
!             -------------
!
              a6(il,ik) =  &
     &          3.0d0 * (qq1a(il,ik) + qq1a(il,ik) -  &
     &                   (ar(il,ik)  + al(il,ik)))
            end do
          end do
!
!
!       ===================
        else if (klmt <= 2) then
!       ===================
!
          lenx = 0
!
          do ik = k1p2, k2m2
            do il = i1, i2
!
              lenx = lenx + 1
!
              al1  (lenx) = al  (il,ik)
              ar1  (lenx) = ar  (il,ik)
              dca1 (lenx) = dca (il,ik)
              qq1a1(lenx) = qq1a(il,ik)
!
              a61  (lenx) = 3.0d0 * (qq1a1(lenx) + qq1a1(lenx) -  &
     &                               (al1(lenx)  + ar1(lenx)))
            end do
          end do
!
!         ===========
!DIR$     INLINE
          call Lmtppm  &
     &      (lenx, klmt, a61, al1, ar1, dca1, qq1a1)
!DIR$     NOINLINE
!         ===========
!
          lenx = 0
!
          do ik = k1p2, k2m2
            do il = i1, i2
!
              lenx = lenx + 1
!
              a6  (il,ik) = a61  (lenx)
              al  (il,ik) = al1  (lenx)
              ar  (il,ik) = ar1  (lenx)
              dca (il,ik) = dca1 (lenx)
              qq1a(il,ik) = qq1a1(lenx)
!
            end do
          end do
!
!
        end if
!
!
        do ik = k1, k2m1
          do il = i1, i2
!
            if (wza(il,ik) > 0.0d0) then
!
              cm = wza(il,ik) / dlp1a(il,ik)
!
              dca(il,ik+1) =  &
     &          ar(il,ik) +  &
     &          0.5d0 * cm *  &
     &          (al(il,ik) - ar(il,ik) +  &
     &           a6(il,ik) * (1.0d0 - r23 * cm))
!
            else
!
              cp = wza(il,ik) / dlp1a(il,ik+1)
!
              dca(il,ik+1) =  &
     &          al(il,ik+1) +  &
     &          0.5d0 * cp *  &
     &          (al(il,ik+1) - ar(il,ik+1) -  &
     &           a6(il,ik+1) * (1.0d0 + r23 * cp))
!
            end if
!
          end do
        end do
!
!
        do ik = k1, k2m1
          dca(:,ik+1) = wza(:,ik) * dca(:,ik+1)
!.sds.. save vertical flux for species as diagnostic
          fz(i1:i2,ij,ik+1) = dca(:,ik+1)
        end do
!
!
        dq1(i1:i2,ij,k1) = dq1(i1:i2,ij,k1) - dca(:,k1p1)
        dq1(i1:i2,ij,k2) = dq1(i1:i2,ij,k2) + dca(:,k2)
!
!
        do ik = k1p1, k2m1
!
          dq1(i1:i2,ij,ik) =  &
     &      dq1(i1:i2,ij,ik) + dca(:,ik) - dca(:,ik+1)
!
        end do
!
!     =============
      end do ijloop
!     =============
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
!   Lmtppm
!
! DESCRIPTION
!   This routine enforces the full monotonic, semi-monotonic, or the
!   positive-definite constraint to the sub-grid parabolic distribution
!   of the Piecewise Parabolic Method (PPM).
!
! ARGUMENTS
!   lenx : vector length
!   lmt  : if 0 => full monotonicity;
!          if 1 => semi-monotonic constraint (no undershoots);
!          if 2 => positive-definite constraint
!   a6   : curvature of the test parabola
!   al   : left  edge value of the test parabola
!   ar   : right edge value of the test parabola
!   dc   : 0.5 * mismatch
!   qa   : cell-averaged value
!
!-----------------------------------------------------------------------------
!
      subroutine Lmtppm  &
     &  (lenx, lmt, a6, al, ar, dc, qa)
!
      implicit none
!
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer :: lenx
      integer, intent(in) :: lmt
      real*8  :: a6(lenx)
      real*8  :: al(lenx)
      real*8  :: ar(lenx)
      real*8  :: dc(lenx)
      real*8  :: qa(lenx)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il
!
      real*8  :: a6da
      real*8  :: da1, da2
      real*8  :: fmin, ftmp
      real*8  :: r12
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      r12 = 1.0d0 / 12.0d0
!
!
!     =============
      if (lmt == 0) then
!     =============
!
!       ----------------
!       Full constraint.
!       ----------------
!
        do il = 1, lenx
!
          if (dc(il) == 0.0d0) then
!
            a6(il) = 0.0d0
            al(il) = qa(il)
            ar(il) = qa(il)
!
          else
!
            da1  = ar(il) - al(il)
            da2  = da1    * da1
            a6da = a6(il) * da1
!
            if (a6da < -da2) then
!
              a6(il) = 3.0d0 * (al(il) - qa(il))
              ar(il) = al(il) - a6(il)
!
            else if (a6da > da2) then
!
              a6(il) = 3.0d0 * (ar(il) - qa(il))
              al(il) = ar(il) - a6(il)
!
            end if
!
          end if
!
        end do
!
!
!     ==================
      else if (lmt == 1) then
!     ==================
!
!       --------------------------
!       Semi-monotonic constraint.
!       --------------------------
!
        do il = 1, lenx
!
          if (Abs (ar(il) - al(il)) < -a6(il)) then
!
            if ((qa(il) < ar(il)) .and. (qa(il) < al(il))) then
!
              a6(il) = 0.0d0
              al(il) = qa(il)
              ar(il) = qa(il)
!
            else if (ar(il) > al(il)) then
!
              a6(il) = 3.0d0 * (al(il) - qa(il))
              ar(il) = al(il) - a6(il)
!
            else
!
              a6(il) = 3.0d0 * (ar(il) - qa(il))
              al(il) = ar(il) - a6(il)
!
            end if
!
          end if
!
        end do
!
!
!     ==================
      else if (lmt == 2) then
!     ==================
!
        do il = 1, lenx
!
          if (Abs (ar(il) - al(il)) < -a6(il)) then
!
            ftmp = ar(il) - al(il)
!
            fmin = qa(il) +  &
     &             0.25d0 * (ftmp * ftmp) / a6(il) +  &
     &             a6(il) * r12
!
            if (fmin < 0.0d0) then
!
              if ((qa(il) < ar(il)) .and. (qa(il) < al(il))) then
!
                a6(il) = 0.0d0
                al(il) = qa(il)
                ar(il) = qa(il)
!
              else if (ar(il) > al(il)) then
!
                a6(il) = 3.0d0 * (al(il) - qa(il))
                ar(il) = al(il) - a6(il)
!
              else
!
                a6(il) = 3.0d0 * (ar(il) - qa(il))
                al(il) = ar(il) - a6(il)
!
              end if
!
            end if
!
          end if
!
        end do
!
      end if
!
!
      return
!
      end
!
!
