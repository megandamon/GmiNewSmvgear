!
!
!=============================================================================
!
! $Id: tpcore_dao2.F90,v 1.13 2013-07-31 15:37:39 ssteenro Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Shian-Jiann Lin, DAO)
!   jrt@llnl.gov
!
! FILE
!   tpcore_dao2.F
!
! ROUTINES
!   Do_Tpcore_Dao2
!   Set_Cross_Terms
!   Calc_Vert_Mass_Flux
!   Set_Jn_Js
!   Calc_Advec_Cross_Terms
!   Qckxyz
!   Set_Lmts
!   Alloc_Tpcore_Arrays
!   Dealloc_Tpcore_Arrays
!
!=============================================================================
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Tpcore_Dao2
!
! DESCRIPTION
!   This routine takes horizontal winds on sigma (or hybrid sigma-p)
!   surfaces and calculates mass fluxes, and then updates the 3D mixing
!   ratio fields one time step (tdt).  The basic scheme is a Multi-
!   Dimensional Flux Form Semi-Lagrangian (FFSL) based on the van Leer or
!   PPM (see Lin and Rood, 1995).
!
!   -------------------------------------------------------
!   See additional external documentation for more details.
!   -------------------------------------------------------
!
! ARGUMENTS
!   ----------------------------------------------------
!   More details on arguments in external documentation.
!   ----------------------------------------------------
!   first_this_tstp  : fist time this routine called this timestep
!   ispc             : species/group index
!   advec_consrv_opt : advection conserve option
!   geofac_pc : special geometrical factor (geofac) for Polar cap
!   geofac    : geometrical factor for meridional advection; geofac uses
!               correct spherical geometry, and replaces acosp as the
!               meridional geometrical factor in tpcore
!   cose      : cosine of grid box edges
!   dap       : pressure difference across layer from (ai * pt) term (mb)
!   dbk       : difference in bi across layer - the dSigma term
!   pres1     : surface pressure at t1     (mb)
!   pres2     : surface pressure at t1+tdt (mb)
!   xmass     : mass flux of air in E-W direction, at t1+tdt/2 (mb)
!   ymass     : mass flux of air in N-S direction, at t1+tdt/2 (mb)
!   zmass     : mass flux of air in vrt direction, at t1+tdt/2 (mb)
!   qq1       : single species/group concentration (mixing ratio)
!   fx        : constituent E-W mass flux
!   fy        : constituent N-S mass flux
!   fz        : constituent vertical mass flux
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Tpcore_Dao2  &
     &  (first_this_tstp, ispc, advec_consrv_opt,  &
     &   geofac_pc, geofac, cose, dap, dbk, pres1, pres2,  &
     &   xmass, ymass, zmass, qq1, fx, fy, fz, &
     &   pr_flux, pr_diag, &
     &   procID, num_species, numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &   i1, i2, ju1, j2, k1, k2, ilong, ivert, &
     &   ewflag, nspoleu, nb, sb, eb, wb,  &
     &   north_domain, south_domain, east_domain, west_domain,  &
     &   communicatorWorld, commu_npole, commu_spole, mapi_all)
!
      use GmiSubDomainsBC_mod, only : Gmi_Bc3du_For_Sub
      use GmiCalcDivergence_mod, only : calcDivergence
      use GmiMassFluxes_mod    , only : calcVerticalMassFlux
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
      logical, intent(in) :: pr_flux
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: gmi_nborder
      integer, intent(in) :: numDomains, numLonDomains, numLatDomains
      integer, intent(in) :: num_species
      integer, intent(in) :: ilong, ivert
      integer, intent(in) :: j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: communicatorWorld
      integer, intent(in) :: commu_npole, commu_spole
      integer, intent(in) :: nb, sb, eb, wb
      integer, intent(in) :: mapi_all(2,numDomains)
      integer, intent(in) :: ewflag(ilo:ihi)
      integer, intent(in) :: nspoleu(julo:jhi)
      integer, intent(in) :: north_domain, south_domain, east_domain, west_domain
!
!
#     include "tpcore_arrays_dao2.h"
!
      logical, intent(in) :: first_this_tstp
      integer, intent(in) :: ispc
      integer, intent(in) :: advec_consrv_opt
      real*8, intent(in)  :: geofac_pc
      real*8, intent(in)  :: geofac(ju1_gl:j2_gl)
      real*8  :: cose  (ju1_gl:j2_gl+1)
      real*8  :: dap   (k1:k2)
      real*8, intent(in) :: dbk   (k1:k2)
      real*8, intent(in) :: pres1 (ilo:ihi, julo:jhi)
      real*8, intent(in) :: pres2 (ilo:ihi, julo:jhi)
      real*8, intent(inout) :: xmass (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(inout) :: ymass (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent  (out) :: zmass (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(inout) :: qq1   (ilo:ihi, julo:jhi, k1:k2)
!.sds.. species fluxes
      real*8, intent(out) :: fx (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: fy (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: fz (ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      logical, save :: first = .true.
!
      integer :: ik
      integer :: num, k2m1
!
!     ----------------------------------------------------
!     ilmt : controls various options in E-W     advection
!     jlmt : controls various options in N-S     advection
!     klmt : controls various options in vertcal advection
!     ----------------------------------------------------
!
      integer, save :: ilmt, jlmt, klmt
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Do_Tpcore_Dao2 called by ', procID
      end if
!
!
      if (first) then
!
        first = .false.
!
!       =============
        call Set_Lmts  &
!       =============
     &    (ilmt, jlmt, klmt, pr_diag, procID, i2_gl, j2_gl)
!
      end if
!
!
      if (first_this_tstp) then
!
!       ====================
        call Set_Press_Terms  &
!       ====================
     &    (dap, dbk, pres1, pres2, delp1, delpm, pu, &
     &   pr_diag, procID, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &   j1p, j2p, i1, i2, ju1, j2, k1, k2)
!
!...intent(in)  dap - difference in ai across layer (mb)
!...intent(in)  dbk - difference in bi across layer (mb)
!...intent(in)  pres1 - surface pressure at t1 (mb)
!...intent(in)  pres2 - surface pressure at t1+tdt (mb)
!...intent(out) delp1 - pressure thickness at t1 (mb)
!...intent(out) delpm - pressure thickness at t1+tdt/2 (mb)
!...intent(out) pu - pressure at edges of box for "u" (mb)
!
!       ======================
        call Gmi_Bc3du_For_Sub &
!       ======================
     &   (pu, ACTM_BC_PU, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  & 
     &       communicatorWorld, numDomains)
!
!       =================
        call Calc_Courant  &
!       =================
     &    (cose, delpm, pu, xmass, ymass, crx, cry, &
     &   pr_diag, procID, j1p, j2p, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)
!
!       ======================
        call Gmi_Bc3du_For_Sub &
!       ======================
     &   (crx, ACTM_BC_CRX, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2, &
     &   ewflag, nspoleu, nb, sb, eb, wb,  &
     &   north_domain, south_domain, east_domain, west_domain,  & 
     &   communicatorWorld,  numDomains)
!
!       ======================
        call Gmi_Bc3du_For_Sub &
!       ======================
     &   (cry, ACTM_BC_CRY, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2, &
     &   ewflag, nspoleu, nb, sb, eb, wb,  &
     &   north_domain, south_domain, east_domain, west_domain,  & 
     &   communicatorWorld, numDomains)
!
!       ====================
        call calcDivergence  &
!       ====================
     &    (.true., geofac_pc, geofac, dpi, xmass, ymass, &
     &   pr_diag, procID, numLonDomains, j1p, j2p, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole)
!
        dps_ctm(:,:) = Sum (dpi(:,:,:), dim=3)
!
!       ====================
        call Set_Cross_Terms  &
!       ====================
     &    (crx, cry, ua, va, &
     &   pr_diag, procID, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld, mapi_all, commu_npole, commu_spole, ivert)
!
      call calcVerticalMassFlux (dbk, dps_ctm, dpi, wz, &
     &      pr_diag, procID, i1, i2, ju1, j2, k1, k2)
!
!.sds2.. have all mass flux here: east-west(xmass), north-south(ymass), vertical(wz)
!.sds2.. save omega (vertical flux) as diagnostic
        if (pr_flux) then
          k2m1 = k2 - 1
          zmass(i1:i2,ju1:j2,k1:k2m1) = -wz(i1:i2,ju1:j2,k2m1:k1:-1)
         endif
!
!
!       ==============
        call Set_Jn_Js  &
!       ==============
     &    (jn, js, crx, &
     &  pr_diag, procID, numDomains, &
     &  ilo, ihi, julo, jhi, ju1_gl, j2_gl, j1p, j2p, &
     &  i1, i2, ju1, j2, k1, k2, communicatorWorld)
!
      end if
!... end "first this time step only" block
!
!.sds.. convert to "mass"
      dq1(:,:,:) = qq1(:,:,:) * delp1(:,:,:)
!
!
!!DIR$ INLINE
!     ===========================
      call Calc_Advec_Cross_Terms  &
!     ===========================
     &  (jn, js, qq1, qqu, qqv, ua, va, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   j1p, j2p, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &   i1, i2, ju1, j2, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld)
!.sds.. notes on arrays
!  qq1 (in) - species mixing ratio
!  qqu (out)- concentration contribution from E-W advection cross terms(mixing ratio)
!  qqv (out)- concentration contribution from N-S advection cross terms(mixing ratio)
!  ua  (in) - average of Courant numbers from il and il+1
!  va  (in) - average of Courant numbers from ij and ij+1
!!DIR$ NOINLINE
!
!
!     ----------------------------------------------------
!     Add advective form E-W operator for E-W cross terms.
!     ----------------------------------------------------
!
!     ==============
      call Xadv_Dao2  &
!     ==============
     &  (2, jn, js, adx, qqv, ua, &
     &   ilo, ihi, julo, jhi, ju1_gl, j2_gl, j1p, j2p, &
     &   i1, i2, ju1, j2, k1, k2)
!.sds notes on output arrays
!  adx (out)- cross term due to E-W advection (mixing ratio)
!  qqv (in) - concentration contribution from N-S advection (mixing ratio)
!  ua  (in) - average of Courant numbers from il and il+1
!.sds
!
!     ----------------------------------------------------
!     Add advective form N-S operator for N-S cross terms.
!     ----------------------------------------------------
!
!     ==============
      call Yadv_Dao2  &
!     ==============
     &  (2, ady, qqu, va, &
     &   numLonDomains, gmi_nborder, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   j1p, j2p, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole, numDomains, mapi_all, ivert)
!.sds notes on output arrays
!  ady (out)- cross term due to N-S advection (mixing ratio)
!  qqu (in) - concentration contribution from N-S advection (mixing ratio)
!  va  (in) - average of Courant numbers from il and il+1
!.sds
!
!... update constituent array qq1 by adding in cross terms - use in fzppm
      qq1(i1:i2,ju1:j2,:) =  &
     &  qq1(i1:i2,ju1:j2,:) +  &
     &  ady(i1:i2,ju1:j2,:) + adx(i1:i2,ju1:j2,:)
!
!
!     ========
      call Xtp  &
!     ========
     &  (ilmt, jn, js, pu, crx, dq1, qqv, xmass, fx, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   j1p, j2p, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &   i1, i2, ju1, j2, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld)
!
!.sds notes on output arrays
!  pu  (in) - pressure at edges in "u" (mb)
!  crx (in) - Courant number in E-W direction
!  dq1 (inout) - species density (mb) - updated with the E-W flux (fx in Xtp)
!  qqv (inout) - concentration contribution from N-S advection (mixing ratio)
!  xmass(in) - horizontal mass flux in E-W direction (mb)
!  fx  (out) - species E-W mass flux
!.sds
!
!     ========
      call Ytp  &
!     ========
     &  (jlmt, geofac_pc, geofac, cry, dq1, qqu, qqv, ymass, fy, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ivert, ilong, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld, &
     & commu_npole, commu_spole, mapi_all)
!
!.sds notes on output arrays
!  cry (in) - Courant number in N-S direction
!  dq1 (inout) - species density (mb) - updated with the N-S flux (fy in Ytp)
!  qqu (in) - concentration contribution from E-W advection (mixing ratio)
!  qqv (inout) - concentration contribution from N-S advection (mixing ratio)  
!  ymass(in) - horizontal mass flux in E-W direction (mb)
!  fy  (out) - species N-S mass flux (need to mult by geofac)
!.sds
!
!     ==========
      call Fzppm  &
!     ==========
     &  (ispc, klmt, delp1, dpi, wz, dq1, qq1, fz, &
     &   j1p, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &   ilong, ivert, i1, i2, ju1, j2, k1, k2)
!
!.sds notes on output arrays
!   wz  (in) : vertical mass flux
!   dq1 (inout) : species density (mb)
!   qq1 (in) : species concentration (mixing ratio)
!.sds
!
!
      if (FILL) then
!       ===========
        call Qckxyz  &
!       ===========
     &    (ispc, dq1, &
     &   j1p, j2p, ju1_gl, j2_gl, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)
      end if
!
!
      if (advec_consrv_opt == 0) then
!
        do ik = k1, k2
!
          delp2(i1:i2,ju1:j2,ik) =  &
     &      dap(ik) +  &
     &      (dbk(ik) * (pres1(i1:i2,ju1:j2) +  &
     &      dps_ctm(i1:i2,ju1:j2)))
!
        end do
!
      else if ((advec_consrv_opt == 1) .or.  &
     &         (advec_consrv_opt == 2)) then
!
        do ik = k1, k2
!
          delp2(i1:i2,ju1:j2,ik) =  &
     &      dap(ik) +  &
     &      (dbk(ik) * pres2(i1:i2,ju1:j2))
!
        end do
!
      end if
!
!
      qq1(i1:i2,ju1:j2,:) =  &
     &  dq1(i1:i2,ju1:j2,:) / delp2(i1:i2,ju1:j2,:)
!
!
      if (j1p /= ju1_gl+1) then
!
        if (ju1 == ju1_gl) then
          qq1(i1:i2,ju1+1,:) = qq1(i1:i2,ju1,:)
        end if
!
        if (j2 == j2_gl) then
          qq1(i1:i2,j2-1,:)  = qq1(i1:i2,j2,:)
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
!   Set_Cross_Terms
!
! DESCRIPTION
!   This routine sets the cross terms for E-W horizontal advection.
!
! ARGUMENTS
!   crx : Courant number in E-W direction
!   cry : Courant number in N-S direction
!   ua  : average of Courant numbers from il and il+1
!   va  : average of Courant numbers from ij and ij+1
!
!-----------------------------------------------------------------------------
!
      subroutine Set_Cross_Terms  &
     &  (crx, cry, ua, va, &
     &   pr_diag, procID, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   ewflag, nspoleu, nb, sb, eb, wb,  &
     &   north_domain, south_domain, east_domain, west_domain,  &
     &   communicatorWorld, &
     &   mapi_all, commu_npole, commu_spole, ivert)
!
      use GmiSubDomainsBC_mod, only : Gmi_Bc3du_For_Sub
      use GmiCalcDivergence_mod, only : Deter_Jrange
!
      implicit none
!
#     include "gem_msg_numbers.h"
#     include "tpcore_constants_dao2.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: gmi_nborder
      integer, intent(in) :: numDomains, numLonDomains, numLatDomains
      integer, intent(in) :: j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ivert
      integer, intent(in) :: communicatorWorld, commu_npole, commu_spole
      integer, intent(in) :: nb, sb, eb, wb
      integer, intent(in) :: mapi_all(2,numDomains)
      integer, intent(in) :: ewflag(ilo:ihi)
      integer, intent(in) :: nspoleu(julo:jhi)
      integer, intent(in) :: north_domain, south_domain, east_domain, west_domain
      real*8, intent(in) :: crx(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in) :: cry(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: ua (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: va (ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij
      integer :: jst, jend
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Set_Cross_Terms called by ', procID
      end if
!
!
      if (.not. CROSS) then
!
        ua(:,:,:) = 0.0d0
        va(:,:,:) = 0.0d0
!
      else
!
        call Deter_Jrange (ju1, j2, j1p, j2p, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
        do ij = jst, jend
          do il = i1, i2
!
            ua(il,ij,:) = 0.5d0 * (crx(il,ij,:) + crx(il+1,ij,:))
!
          end do
        end do
!
!
        call Deter_Jrange (ju1, j2, ju1+1, j2-1, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
        do ij = jst, jend
          do il = i1, i2
!
            va(il,ij,:) = 0.5d0 * (cry(il,ij,:) + cry(il,ij+1,:))
!
          end do
        end do
!
!
!       =============================
        call Do_Cross_Terms_Pole_I2d2  &
!       =============================
     &    (cry, va, &
     &   numLonDomains, i1_gl, i2_gl, ju1_gl, j2_gl, j1p,  &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   mapi_all, commu_npole, commu_spole, numDomains, ivert)
!
!
!       ======================
!c?     call Gmi_Bc3du_For_Sub (ua, ACTM_BC_UA)
        call Gmi_Bc3du_For_Sub (va, ACTM_BC_VA, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2, &
     &   ewflag, nspoleu, nb, sb, eb, wb,  &
     &  north_domain, south_domain, east_domain, west_domain,  & 
     &  communicatorWorld, numDomains)
!       ======================
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
!   Set_Jn_Js
!
! DESCRIPTION
!   This routine determines jn and js, by looking where Courant # is > 1.
!
! ARGUMENTS
!   jn  : northward of latitude index = jn; Courant numbers could be > 1,
!         so use the flux-form semi-Lagrangian scheme
!   js  : southward of latitude index = js; Courant numbers could be > 1,
!         so use the flux-form semi-Lagrangian scheme
!   crx : Courant number in E-W direction
!
!-----------------------------------------------------------------------------
!
      subroutine Set_Jn_Js  &
     &  (jn, js, crx, &
     &  pr_diag, procID, numDomains, &
     &  ilo, ihi, julo, jhi, ju1_gl, j2_gl, j1p, j2p, &
     &  i1, i2, ju1, j2, k1, k2, communicatorWorld)
!
      use GmiReduce_mod, only : Gmi_Max_Reduce, Gmi_Min_Reduce
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: numDomains, communicatorWorld
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: ju1_gl, j2_gl, j1p, j2p
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer :: jn (k1:k2)
      integer :: js (k1:k2)
      real*8  :: crx(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik
      integer :: jn0, js0
      integer :: jst, jend
!
      real*8  :: rjn(k1:k2)
      real*8  :: rjs(k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Set_Jn_Js called by ', procID
      end if
!
!
      js0  = (j2_gl + 1 ) / 2
      jn0  = j2_gl - js0 + 1
!
      jst  = Max (ju1, j1p)
      jend = Min (j2,  js0)
!
!
      ikloop1: do ik = k1, k2
!
        js(ik) = j1p
!
        do ij = jend, jst, -1
          do il = i1, i2
!
            if (Abs (crx(il,ij,ik)) > 1.0d0) then
!
              js(ik) = ij
!
!             =============
              cycle ikloop1
!             =============
!
            end if
!
          end do
        end do
!
      end do ikloop1
!
!
      jst  = Max (ju1, jn0)
      jend = Min (j2,  j2p)
!
!
      ikloop2: do ik = k1, k2
!
        jn(ik) = j2p
!
        do ij = jst, jend
          do il = i1, i2
!
            if (Abs (crx(il,ij,ik)) > 1.0d0) then
!
              jn(ik) = ij
!
!             =============
              cycle ikloop2
!             =============
!
            end if
!
          end do
        end do
!
      end do ikloop2
!
!
      rjn(:) = jn(:)
      rjs(:) = js(:)
!
!     ======================
      call Gmi_Min_Reduce  &
!     ======================
     &  (k1, k2, rjn, numDomains, communicatorWorld)
!
!     ======================
      call Gmi_Max_Reduce  &
!     ======================
     &  (k1, k2, rjs, numDomains, communicatorWorld)
!
      jn(:) = Nint (rjn(:))
      js(:) = Nint (rjs(:))
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
!   Calc_Advec_Cross_Terms
!
! DESCRIPTION
!   This routine calculates the advective cross terms.
!
! ARGUMENTS
!   jn  : northward of latitude index = jn, Courant numbers could be > 1,
!         so use the flux-form semi-Lagrangian scheme
!   js  : southward of latitude index = js, Courant numbers could be > 1,
!         so use the flux-form semi-Lagrangian scheme
!   qq1 : species concentration (mixing ratio)
!   qqu : concentration contribution from E-W advection (mixing ratio)
!   qqv : concentration contribution from N-S advection (mixing ratio)
!   ua  : average of Courant numbers from il and il+1
!   va  : average of Courant numbers from ij and ij+1
!
!-----------------------------------------------------------------------------
!
      subroutine Calc_Advec_Cross_Terms  &
     &  (jn, js, qq1, qqu, qqv, ua, va, &
     &   numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   j1p, j2p, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &   i1, i2, ju1, j2, k1, k2,  ewflag, nspoleu, nb, sb, eb, wb,  &
     &   north_domain, south_domain, east_domain, west_domain,  & 
     &   communicatorWorld)
!
      use GmiSubDomainsBC_mod, only : Gmi_Bc3du_For_Sub
      use GmiCalcDivergence_mod, only : Deter_Jrange
!
      implicit none
!
#     include "gem_msg_numbers.h"
#     include "tpcore_constants_dao2.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numDomains, numLonDomains, numLatDomains, gmi_nborder
      integer, intent(in) :: j1p, j2p, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: communicatorWorld
      integer, intent(in) :: nb, sb, eb, wb
      integer, intent(in) :: ewflag(ilo:ihi)
      integer, intent(in) :: nspoleu(julo:jhi)
      integer, intent(in) :: north_domain, south_domain, east_domain, west_domain
      integer, intent(in) :: jn (k1:k2)
      integer, intent(in) :: js (k1:k2)
      real*8, intent(in)  :: qq1(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: qqu(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: qqv(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: ua (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: va (ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik, iu
      integer :: jst, jend
      integer :: jv
!
      real*8  :: ril, rij, riu
      real*8  :: ru
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
!     ================
      if (.not. CROSS) then
!     ================
!
        qqu(:,:,:) = qq1(:,:,:)
        qqv(:,:,:) = qq1(:,:,:)
!
!
!     ====
      else
!     ====
!
        qqu(:,:,:) = 0.0d0
        qqv(:,:,:) = 0.0d0
!
!
        call Deter_Jrange (ju1, j2, j1p, j2p, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!
        do ik = k1, k2
          do ij = jst, jend
!
            if ((ij <= js(ik)) .or. (ij >= jn(ik))) then
!
!             ----------------------------------------------------------
!             In Polar area, so need to deal with large courant numbers.
!             ----------------------------------------------------------
!
              do il = i1, i2
!
!c?
                iu  = ua(il,ij,ik)
                riu = iu
                ru  = ua(il,ij,ik) - riu
                iu  = il - iu
!
                if (ua(il,ij,ik) >= 0.0d0) then
!
                  qqu(il,ij,ik) =  &
     &              qq1(iu,ij,ik) +  &
     &              ru * (qq1(iu-1,ij,ik) - qq1(iu,ij,ik))
!
                else
!
                  qqu(il,ij,ik) =  &
     &              qq1(iu,ij,ik) +  &
     &              ru * (qq1(iu,ij,ik) - qq1(iu+1,ij,ik))
!
                end if
!
                qqu(il,ij,ik) = qqu(il,ij,ik) - qq1(il,ij,ik)
!
              end do
!
            else  ! js < ij < jn
!
!             ---------------------------
!             Do interior area (use PPM).
!             ---------------------------
!
              do il = i1, i2
!
                ril = il
                iu  = ril - ua(il,ij,ik)
!
                qqu(il,ij,ik) =  &
     &            ua(il,ij,ik) *  &
     &            (qq1(iu,ij,ik) - qq1(iu+1,ij,ik))
!
              end do
!
            end if
!
          end do
!
        end do
!
!
        call Deter_Jrange (ju1, j2, j1p, j2p, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!
        do ik = k1, k2
          do ij = jst, jend
            do il = i1, i2
!
!c?
              rij = ij
              jv  = rij - va(il,ij,ik)
!
              qqv(il,ij,ik) =  &
     &          va(il,ij,ik) *  &
     &          (qq1(il,jv,ik) - qq1(il,jv+1,ik))
!
            end do
          end do
        end do
!
!
        qqu(i1:i2,ju1:j2,:) =  &
     &    qq1(i1:i2,ju1:j2,:) + (0.5d0 * qqu(i1:i2,ju1:j2,:))
!
        qqv(i1:i2,ju1:j2,:) =  &
     &    qq1(i1:i2,ju1:j2,:) + (0.5d0 * qqv(i1:i2,ju1:j2,:))
!
!
!       ======================
        call Gmi_Bc3du_For_Sub (qqu, ACTM_BC_QQU, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  & 
     &       communicatorWorld, numDomains)
!
        call Gmi_Bc3du_For_Sub (qqv, ACTM_BC_QQV, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  & 
     &       communicatorWorld, numDomains)
!       ======================
!
!
!     ======
      end if
!     ======
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
!   Qckxyz
!
! DESCRIPTION
!   This routine checks for "filling".
!
! ARGUMENTS
!   ic  : species index
!   dq1 : species density (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Qckxyz  &
     &  (ic, dq1, &
     &   j1p, j2p, ju1_gl, j2_gl, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)
!
      use GmiCalcDivergence_mod, only : Deter_Jrange
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: j1p, j2p
      integer, intent(in) :: ju1_gl, j2_gl
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer :: ic
      real*8  :: dq1(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Parameter declarations.
!     ----------------------
!
      logical, parameter :: FILL_DIAG = .false.
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik
      integer :: ip
      integer :: jst, jend
      integer :: k1p1, k2m1
!
      real*8  :: dup, qup
      real*8  :: qly
      real*8  :: sum
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      ip = 0
!
!
      call Deter_Jrange (ju1, j2, j1p, j2p, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!
!     ----------
!     Top layer.
!     ----------
!
      k1p1 = k1 + 1
!
      do ij = jst, jend
        do il = i1, i2
!
          if (dq1(il,ij,k1) < 0.0d0) then
!
            ip = ip + 1
!
            dq1(il,ij,k1p1) = dq1(il,ij,k1p1) + dq1(il,ij,k1)
            dq1(il,ij,k1)   = 0.0d0
!
          end if
!
        end do
      end do
!
!
      do ik = k1 + 1, k2 - 1
        do ij = jst, jend
          do il = i1, i2
!
            if (dq1(il,ij,ik) < 0.0d0) then
!
              ip = ip + 1
!
!             -----------
!             From above.
!             -----------
!
              qup =  dq1(il,ij,ik-1)
              qly = -dq1(il,ij,ik)
              dup =  Min (qly, qup)
!
              dq1(il,ij,ik-1) = qup - dup
              dq1(il,ij,ik)   = dup - qly
!
!             -----------
!             From below.
!             -----------
!
              dq1(il,ij,ik+1) = dq1(il,ij,ik+1) + dq1(il,ij,ik)
              dq1(il,ij,ik)   = 0.0d0
!
            end if
!
          end do
        end do
      end do
!
!
!     -------------
!     Bottom layer.
!     -------------
!
      sum  = 0.0d0
!
      k2m1 = k2 - 1
!
      do ij = jst, jend
        do il = i1, i2
!
          if (dq1(il,ij,k2) < 0.0d0) then
!
            ip = ip + 1
!
!           -----------
!           From above.
!           -----------
!
            qup =  dq1(il,ij,k2m1)
            qly = -dq1(il,ij,k2)
            dup = Min (qly, qup)
!
            dq1(il,ij,k2m1) = qup - dup
!
!           -------------------------
!           From "below" the surface.
!           -------------------------
!
            sum = sum + qly - dup
!
            dq1(il,ij,k2) = 0.0d0
!
          end if
!
        end do
      end do
!
!
      if (FILL_DIAG) then
!
        if (ip > 0) then
          Write (6,*) 'WARNING: Negative tracer found and corrected.'
          Write (6,*)  &
     &      'spc = ', ic, ', vert fill pts      = ', ip
        end if
!
        if (sum > 1.0d-25) then
          Write (6,*)  &
     &      'spc = ', ic, ', mass src from grnd = ', sum
        end if
!
        do ik = k1, k2
          do ij = jst, jend
            do il = i1, i2
              if (dq1(il,ij,ik) < 1.0d-30) then
                Write (6,*) il, ij, ik, dq1(il,ij,ik)
              end if
            end do
          end do
        end do
!
      end if
!
!
!     =======================================
      where (dq1(i1:i2,jst:jend,:) < 1.0d-30)  &
     &  dq1(i1:i2,jst:jend,:) = 1.0d-30
!     =======================================
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
!   Set_Lmts
!
! DESCRIPTION
!   This routine sets ilmt, jlmt, klmt.
!
! ARGUMENTS
!   ilmt : controls various options in E-W     advection
!   jlmt : controls various options in N-S     advection
!   klmt : controls various options in vertcal advection
!
!-----------------------------------------------------------------------------
!
      subroutine Set_Lmts  &
     &  (ilmt, jlmt, klmt, &
     &   pr_diag, procID, i2_gl, j2_gl)
!
      implicit none
!
#     include "tpcore_constants_dao2.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: i2_gl, j2_gl
      integer :: ilmt
      integer :: jlmt
      integer :: klmt
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: j2_glm1
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Set_Lmts called by ', procID
      end if
!
      j2_glm1 = j2_gl - 1
!
!c?
      if (IORD <= 0) then
        if (i2_gl >= 144) then
          ilmt = 0
        else if (i2_gl >= 72) then
          ilmt = 1
        else
          ilmt = 2
        end if
      else
        ilmt = IORD - 3
      end if
!
!
!c?
      if (JORD <= 0) then
        if (j2_glm1 >= 90) then
          jlmt = 0
        else if (j2_glm1 >= 45) then
          jlmt = 1
        else
          jlmt = 2
        end if
      else
        jlmt = JORD - 3
      end if
!
      klmt = Max ((KORD-3), 0)
!
      if (pr_diag) then
        Write (6,*)  &
     &    'IORD = ', IORD, ', JORD = ', JORD, ', KORD = ', KORD,  &
     &    procID
        Write (6,*)  &
     &    'ilmt = ', ilmt, ', jlmt = ', jlmt, ', klmt = ', klmt,  &
     &    procID
      end if
!
      return
!
      end
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Alloc_Tpcore_Arrays
!
! DESCRIPTION
!   This routine allocates the tpcore arrays.
!
! ARGUMENTS
!   None
!
!-----------------------------------------------------------------------------
!
      subroutine Alloc_Tpcore_Arrays &
     & (pr_diag, procID, &
     &  ilo, ihi, julo, jhi, &
     &  i1, i2, ju1, j2, k1, k2)
!
      implicit none
!
#     include "tpcore_arrays_dao2.h"
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Alloc_Tpcore_Arrays called by ', procID
      end if
!
!
      Allocate (jn(k1:k2))
      Allocate (js(k1:k2))
      jn = 0; js = 0
!
      Allocate (dps_ctm(i1:i2, ju1:j2))
      dps_ctm = 0.0d0
!
      Allocate (dpi(i1:i2, ju1:j2, k1:k2))
      Allocate (wz (i1:i2, ju1:j2, k1:k2))
      dpi = 0.0d0; wz = 0.0d0
!
      Allocate (adx  (ilo:ihi, julo:jhi, k1:k2))
      Allocate (ady  (ilo:ihi, julo:jhi, k1:k2))
      adx = 0.0d0; ady = 0.0d0
!
      Allocate (crx  (ilo:ihi, julo:jhi, k1:k2))
      Allocate (cry  (ilo:ihi, julo:jhi, k1:k2))
      crx = 0.0d0; cry = 0.0d0
!
      Allocate (delp1(ilo:ihi, julo:jhi, k1:k2))
      Allocate (delp2(ilo:ihi, julo:jhi, k1:k2))
      Allocate (delpm(ilo:ihi, julo:jhi, k1:k2))
      delp1 = 0.0d0; delp2 = 0.0d0
      delpm = 0.0d0
!
      Allocate (dq1  (ilo:ihi, julo:jhi, k1:k2))
      Allocate (pu   (ilo:ihi, julo:jhi, k1:k2))
      dq1 = 0.0d0; pu = 0.0d0
!
      Allocate (qqu  (ilo:ihi, julo:jhi, k1:k2))
      Allocate (qqv  (ilo:ihi, julo:jhi, k1:k2))
      qqu = 0.0d0; qqv = 0.0d0
!
      Allocate (ua   (ilo:ihi, julo:jhi, k1:k2))
      Allocate (va   (ilo:ihi, julo:jhi, k1:k2))
      ua = 0.0d0; va = 0.0d0
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
!   Dealloc_Tpcore_Arrays
!
! DESCRIPTION
!   This routine deallocates the tpcore arrays.
!
! ARGUMENTS
!   None
!
!-----------------------------------------------------------------------------
!
      subroutine Dealloc_Tpcore_Arrays (pr_diag, procID)
!
      implicit none
!
#     include "tpcore_arrays_dao2.h"
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Dealloc_Tpcore_Arrays called by ', procID
      end if
!
!
      Deallocate (jn)
      Deallocate (js)
!
      Deallocate (dps_ctm)
!
      Deallocate (dpi)
      Deallocate (wz)
!
      Deallocate (adx)
      Deallocate (ady)
      Deallocate (crx)
      Deallocate (cry)
      Deallocate (delp1)
      Deallocate (delp2)
      Deallocate (delpm)
      Deallocate (dq1)
      Deallocate (pu)
      Deallocate (qqu)
      Deallocate (qqv)
      Deallocate (ua)
      Deallocate (va)
!
      return
!
      end
!
!
