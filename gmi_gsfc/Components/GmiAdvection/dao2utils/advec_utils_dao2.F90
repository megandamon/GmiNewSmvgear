!
!=============================================================================
!
! $Id: advec_utils_dao2.F90,v 1.10 2013-07-31 15:06:36 ssteenro Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Philip Cameron-Smith, LLNL)
!   jrt@llnl.gov
!
! FILE
!   advec_utils_dao2.F
!
! ROUTINES
!   Calc_Courant
!   Set_Press_Terms
!
!=============================================================================
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calc_Courant
!
! DESCRIPTION
!   This routine calculates courant numbers from the horizontal mass fluxes.
!
! ARGUMENTS
!   cose  : cosine of grid box edges
!   delpm : pressure thickness, the pseudo-density in a hydrostatic system
!           at t1+tdt/2 (approximate) (mb)
!   pu    : pressure at edges in "u"  (mb)
!   xmass : horizontal mass flux in E-W direction (mb)
!   ymass : horizontal mass flux in N-S direction (mb)
!   crx   : Courant number in E-W direction
!   cry   : Courant number in N-S direction
!
!-----------------------------------------------------------------------------
!
      subroutine Calc_Courant  &
     &  (cose, delpm, pu, xmass, ymass, crx, cry, &
     &   pr_diag, loc_proc, j1p, j2p, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)
!
      use GmiCalcDivergence_mod, only : Deter_Jrange
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: j1p, j2p
      integer, intent(in) :: ju1_gl, j2_gl
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      real*8, intent(in)  :: cose (ju1_gl:j2_gl+1)
      real*8, intent(in)  :: delpm(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: pu   (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: xmass(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: ymass(ilo:ihi, julo:jhi, k1:k2)
!
      real*8, intent(out) :: crx(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: cry(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: ij
      integer :: jst, jend
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Calc_Courant called by ', loc_proc
      end if
!
!
      crx(:,:,:) = 0.0d0
      cry(:,:,:) = 0.0d0
!
!
!     -----------------------------------
!     Calculate E-W horizontal mass flux.
!     -----------------------------------
!
      call Deter_Jrange (ju1, j2, j1p, j2p, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!
      do ij = jst, jend
!
        crx(i1:i2+1,ij,:) =  &
     &    xmass(i1:i2+1,ij,:) / pu(i1:i2+1,ij,:)
!
      end do
!
!
!     -----------------------------------
!     Calculate N-S horizontal mass flux.
!     -----------------------------------
!
      call Deter_Jrange (ju1, j2+1, j1p, j2p+1, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
      do ij = jst, jend
!
        cry(i1:i2,ij,:) =  &
     &    ymass(i1:i2,ij,:) /  &
     &    ((0.5d0 * cose(ij)) *  &
     &     (delpm(i1:i2,ij,:) + delpm(i1:i2,ij-1,:)))
!
      end do
!
!
      return
!
      end
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Set_Press_Terms
!
! DESCRIPTION
!   This routine sets the pressure terms:  delp1, delpm, pu.
!
! ARGUMENTS
!   dap   : pressure difference across layer from (ai * pt) term (mb)
!   dbk   : difference in bi across layer - the dSigma term
!   pres1 : surface pressure at t1     (mb)
!   pres2 : surface pressure at t1+tdt (mb)
!   delp1 : pressure thickness, the pseudo-density in a hydrostatic system
!           at t1 (mb)
!   delpm : pressure thickness, the pseudo-density in a hydrostatic system
!           at t1+tdt/2 (approximate)  (mb)
!   pu    : pressure at edges in "u"   (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Set_Press_Terms  &
     &  (dap, dbk, pres1, pres2, delp1, delpm, pu, &
     &   pr_diag, loc_proc, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &   j1p, j2p, i1, i2, ju1, j2, k1, k2)
!
      use GmiCalcDivergence_mod, only : Deter_Jrange
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: ju1_gl, j2_gl, j1p, j2p
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      real*8, intent(in)  :: dap(k1:k2)
      real*8, intent(in)  :: dbk(k1:k2)
      real*8, intent(in)  :: pres1(ilo:ihi, julo:jhi)
      real*8, intent(in)  :: pres2(ilo:ihi, julo:jhi)
!
      real*8, intent(out) :: delp1(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: delpm(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: pu   (ilo:ihi, julo:jhi, k1:k2)
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik
      integer :: jst, jend
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Set_Press_Terms called by ', loc_proc
      end if
!
      do ik = k1, k2
!
        delp1(:,:,ik) = dap(ik) + (dbk(ik) * pres1(:,:))
!
        delpm(:,:,ik) =  &
     &    dap(ik) +  &
     &    (dbk(ik) * 0.5d0 * (pres1(:,:) + pres2(:,:)))
!
      end do
!
      call Deter_Jrange (ju1, j2, j1p, j2p, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
      do ij = jst, jend
        do il = i1, i2
!
          pu(il,ij,:) = 0.5d0 * (delpm(il,ij,:) + delpm(il-1,ij,:))
!
        end do
      end do
!
      return
!
      end
!
!
