  module GmiPressureFixer_mod
!
      use GmiCalcDivergence_mod, only : calcDivergence, Deter_Jrange
      use GmiMassFluxes_mod    , only : calcVerticalMassFlux, calcDelta_xyt
!
      implicit none
!
      private
      public  :: adjustPressFixer
      public  :: Init_Press_Fix, Calc_Delpm
      public  :: Do_Press_Fix_Llnl, Do_Press_Fix_Uci
!
#     include "GmiParameters.h"
!
  contains
!
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: adjustPressFixer
!
! !INTERFACE:
!
      subroutine adjustPressFixer &
     &  (metdata_name_org, met_grid_type, do_timinterp_winds,  &
     &   new_met_rec, advec_consrv_opt, pmet2_opt, press_fix_opt,  &
     &   tdt, geofac_pc, geofac, cose, cosp, rel_area, dap, dbk,  &
     &   pctm1, pctm2, pmet2, uux, vvx, xmass, ymass, zmass, &
     &   pr_diag, procID, numLonDomains, ilong, ilat, j1p, j2p, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jvlo, jhi, i1, i2, ju1, j2, k1, k2, &
     &     commu_npole, commu_spole, gmi_nborder)
!
      use GmiWrapMaster_mod, only : wrapMaster_2d, wrapMaster_3du
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID, numLonDomains, commu_npole, commu_spole
      integer, intent(in) :: ilong, ilat, gmi_nborder
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, j1p, j2p
      integer, intent(in) :: ilo, ihi, julo, jvlo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
                         ! first  part of metdata_name, e.g., "NCAR"
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
                         ! met grid type, 'A' or 'C'
      character (len=1)  :: met_grid_type
                         ! time interpolate wind fields?
      logical, intent(in) :: do_timinterp_winds
                         ! new met record?
      logical, intent(in) :: new_met_rec
                         ! advection_conserve option
      integer, intent(in) :: advec_consrv_opt
                         ! pmet2 option
      integer, intent(in) :: pmet2_opt
                         ! pressure fixer option
      integer, intent(in) :: press_fix_opt
                         ! model time step(s)
      real*8, intent(in) :: tdt
                         ! special geometrical factor (geofac) for Polar cap
      real*8  :: geofac_pc
                         ! geometrical factor for meridional advection; geofac uses
                         ! correct spherical geometry, and replaces acosp as the
                         ! meridional geometrical factor in tpcore
      real*8  :: geofac  (ju1_gl:j2_gl)
                         ! cosine of grid box edges
      real*8  :: cose    (ju1_gl:j2_gl+1)
                         ! cosine of grid box centers
      real*8  :: cosp    (ju1_gl:j2_gl)
                         ! relative surface area of grid box (fraction)
      real*8  :: rel_area(i1:i2, ju1:j2)
                         ! pressure difference across layer from (ai * pt) term (mb)
      real*8  :: dap     (k1:k2)
                         ! difference in bi across layer - the dSigma term
      real*8  :: dbk     (k1:k2)
                         !  CTM surface pressure at t1     (mb)
      real*8  :: pctm1   (ilo:ihi, julo:jhi)
                         ! CTM surface pressure at t1+tdt (mb)
      real*8  :: pctm2   (ilo:ihi, julo:jhi)
                         ! surface pressure     at t1+tdt (mb)
      real*8  :: pmet2   (ilo:ihi, julo:jhi)
                         ! current horizontal velocity in zonal (longitude) direction;
                         ! known at edges in longitude direction and centers in latitude
                         ! direction (m/s)
      real*8  :: uux     (ilo:ihi, julo:jhi, k1:k2)
                         ! current horizontal velocity in meridional(latitude) direction;
                         ! known at edges in latitude direction and centers in longitude
                         ! direction (m/s)
      real*8  :: vvx     (ilo:ihi, jvlo:jhi, k1:k2)
                         ! horizontal mass flux in E-W direction  (mb)
      real*8 , intent(inOut) :: xmass   (ilo:ihi, julo:jhi, k1:k2)
                         ! horizontal mass flux in N-S direction  (mb)
      real*8 , intent(inOut) :: ymass   (ilo:ihi, julo:jhi, k1:k2)
                         ! vertical mass flux  (mb)
      real*8 , intent(inOut) :: zmass   (i1:i2, ju1:j2, k1:k2)
!
! !DESCRIPTION:
! Initializes and calls the pressure fixer code.
!
! !LOCAL VARIABLES:
      logical, save :: DO_ADJUST_PRESS_DIAG = .false.
      logical, save :: first = .true.
!
!     --------------------------------------------------
!     dgpress   : global-pressure discrepancy
!     press_dev : RMS difference between pmet2 and pctm2
!                 (weighted by relative area)
!     --------------------------------------------------
!
      real*8  :: dgpress
      real*8  :: press_dev
!
!     -------------------------------------------------------------
!     dps : change of surface pressure from met field pressure (mb)
!     -------------------------------------------------------------
!
      real*8  :: dps(i1:i2, ju1:j2)
!
!     -----------------------------------------------------------------
!     dps_ctm : CTM surface pressure tendency; sum over vertical of dpi
!               calculated from original mass fluxes (mb)
!     -----------------------------------------------------------------
!
      real*8, allocatable, save :: dps_ctm(:,:)
!EOP
!-----------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6, *) 'adjustPressFixer called by ', procID
!
!     ==========
      if (first) then
!     ==========
!
        first = .false.
!
        Allocate (dps_ctm(i1:i2, ju1:j2))
        dps_ctm = 0.0d0
!
        call calcDelta_xyt (cosp, tdt, i2_gl, ju1_gl, j2_gl)
!
      end if
!
      dgpress =  Sum ((pmet2(i1_gl:i2_gl,ju1_gl:j2_gl) -  &
     &                 pctm1(i1_gl:i2_gl,ju1_gl:j2_gl)) *  &
     &                rel_area(i1_gl:i2_gl,ju1_gl:j2_gl))
!
      if (pmet2_opt == 1) then
        pmet2(:,:) = pmet2(:,:) - dgpress
      end if
!
!
      if (DO_ADJUST_PRESS_DIAG) then
        Write (6, *) 'Global mean surface pressure change (mb) = ',  &
     &                dgpress
      end if
!
!
      if (new_met_rec .or. do_timinterp_winds) then
!
!       ===================
        call Init_Press_Fix  &
!       ===================
     &    (metdata_name_org, met_grid_type, tdt, geofac_pc, geofac,  &
     &     cose, cosp, dap, dbk, dps, dps_ctm, rel_area, pctm1, pmet2,  &
     &     uux, vvx, xmass, ymass, zmass, pr_diag, procID, numLonDomains, &
     &     j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jvlo, jhi, &
     &     i1, i2, ju1, j2, k1, k2, commu_npole, commu_spole, gmi_nborder)
!
!
        if (DO_ADJUST_PRESS_DIAG) then
!
!         -------------------------------------------------------
!         Calculate the RMS pressure deviation (diagnostic only).
!         -------------------------------------------------------
!
          press_dev =  &
     &      Sqrt (Sum ((dps    (i1_gl:i2_gl,ju1_gl:j2_gl) -  &
     &                  dps_ctm(i1_gl:i2_gl,ju1_gl:j2_gl))**2 *  &
     &                  rel_area(i1_gl:i2_gl,ju1_gl:j2_gl)))
!
          Write (6,*) 'Before Pressure-Fixing:  Area weighted RMS deviation '
          Write (6,*) '  between pmet2 & pctm2 (mb) = ', press_dev
!
        end if
!
!
        if (press_fix_opt == 1) then
!
!         ======================
          call Do_Press_Fix_Llnl  &
!         ======================
     &      (geofac_pc, geofac, dbk, dps, dps_ctm, rel_area,  &
     &       xmass, ymass, &
     &       pr_diag, procID, numLonDomains, j1p, j2p, i1_gl, i2_gl, &
     &      ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &      commu_npole, commu_spole, gmi_nborder)
!
        else if (press_fix_opt == 2) then
!
!         =====================
          call Do_Press_Fix_Uci  &
!         =====================
     &      (geofac_pc, geofac, dbk, dps, dps_ctm, rel_area,  &
     &       xmass, ymass, &
     &       pr_diag, procID, numLonDomains, ilong, ilat, j1p, j2p, i1_gl, i2_gl, &
     &       ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &      commu_npole, commu_spole, gmi_nborder)
!
        end if
!
!
        if ((advec_consrv_opt == 0) .or.  &
     &      (advec_consrv_opt == 1)) then
!
          dps_ctm(i1:i2,ju1:j2) =  &
     &      pmet2(i1:i2,ju1:j2) - pctm1(i1:i2,ju1:j2)
!
!       -----------------------------------------------
!       else if (advec_consrv_opt == 2) then do nothing
!       -----------------------------------------------
!
        end if
!
      end if
!
      pctm2(i1:i2,ju1:j2) =  &
     &  pctm1(i1:i2,ju1:j2) + dps_ctm(i1:i2,ju1:j2)
!
!     ------------------------------------------------
!     Always done on Master, so only "wrap" necessary.
!     ------------------------------------------------
!
      call wrapMaster_2d  (pctm1, i1, i2, ju1, j2, ilo, ihi, julo, jhi, gmi_nborder)
      call wrapMaster_2d  (pctm2, i1, i2, ju1, j2, ilo, ihi, julo, jhi, gmi_nborder)
      call wrapMaster_3du (xmass, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, gmi_nborder)
      call wrapMaster_3du (ymass, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, gmi_nborder)
!
      if (DO_ADJUST_PRESS_DIAG) then
!
!       -------------------------------------------------------
!       Calculate the RMS pressure deviation (diagnostic only).
!       -------------------------------------------------------
!
        press_dev =  &
     &    Sqrt (Sum ((pmet2(i1_gl:i2_gl,ju1_gl:j2_gl) -  &
     &                 pctm2(i1_gl:i2_gl,ju1_gl:j2_gl))**2 *  &
     &                rel_area(i1_gl:i2_gl,ju1_gl:j2_gl)))
!
        Write (6,*) 'After Pressure-Fixing:  Area weighted RMS deviation '
        Write (6,*) '  between pmet2 & pctm2 (mb) = ', press_dev
!
      end if
!
      return
!
      end subroutine adjustPressFixer
!EOC
!-----------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calc_Delpm
!
! DESCRIPTION
!   This routine calculates the pressure term, delpm.
!
! ARGUMENTS
!   dap   : pressure difference across layer from (ai * pt) term (mb)
!   dbk   : difference in bi across layer - the dSigma term
!   pres1 : surface pressure at t1     (mb)
!   pres2 : surface pressure at t1+tdt (mb)
!   delpm : pressure thickness, the psudo-density in a hydrostatic system
!           at t1+tdt/2 (approximate)  (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Calc_Delpm  &
     &  (dap, dbk, pres1, pres2, delpm, &
     &   pr_diag, procID, k1, k2, ilo, ihi, julo, jhi)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: k1, k2, ilo, ihi, julo, jhi
      real*8, intent(in)  :: dap(k1:k2)
      real*8, intent(in)  :: dbk(k1:k2)
      real*8, intent(in)  :: pres1(ilo:ihi, julo:jhi)
      real*8, intent(in)  :: pres2(ilo:ihi, julo:jhi)
!
      real*8, intent(out) :: delpm(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: ik
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Calc_Delpm called by ', procID
      end if
!
      do ik = k1, k2
!
        delpm(:,:,ik) =  &
     &    dap(ik) +  &
     &    (dbk(ik) * 0.5d0 * (pres1(:,:) + pres2(:,:)))
!
      end do
!
!
      return
!
      end subroutine Calc_Delpm
!
!=============================================================================
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Philip Cameron-Smith, LLNL)
!   jrt@llnl.gov
!
! FILE
!   gmi_press_fixer.F
!
! ROUTINES
!   Init_Press_Fix
!   Do_Press_Fix_Llnl
!   Average_Press_Poles
!   Calc_Horiz_Mass_Flux
!   Calc_Horiz_Mass_Flux_Giss
!
!=============================================================================
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Init_Press_Fix
!
! DESCRIPTION
!   This routine initializes the pressure fixer.
!
! ARGUMENTS
!   metdata_name_org : first  part of metdata_name, e.g., "NCAR"
!   met_grid_type    : met grid type, 'A' or 'C'
!   tdt       : model time step (s)
!   geofac_pc : special geometrical factor (geofac) for Polar cap
!   geofac    : geometrical factor for meridional advection; geofac uses
!               correct spherical geometry, and replaces acosp as the
!               meridional geometrical factor in tpcore
!   cose      : cosine of grid box edges
!   cosp      : cosine of grid box centers
!   dap       : pressure difference across layer from (ai * pt) term (mb)
!   dbk       : difference in bi across layer - the dSigma term
!   dps       : change of surface pressure from met field pressure   (mb)
!   dps_ctm   : CTM surface pressure tendency; sum over vertical of dpi
!               calculated from original mass fluxes  (mb)
!   rel_area  : relative surface area of grid box (fraction)
!   pctm1     : CTM       surface pressure at t1      (mb)
!   pmet2     : met field surface pressure at t1+tdt  (mb)
!   uu        : wind velocity in E-W direction, at t1+tdt/2 (m/s)
!   vv        : wind velocity in N-S direction, at t1+tdt/2 (m/s)
!   vv        : wind velocity in N-S direction        (m/s)
!   xmass     : horizontal mass flux in E-W direction (mb)
!   ymass     : horizontal mass flux in N-S direction (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Init_Press_Fix  &
     &  (metdata_name_org, met_grid_type, tdt, geofac_pc, geofac,  &
     &   cose, cosp, dap, dbk, dps, dps_ctm, rel_area, pctm1, pmet2,  &
     &   uu, vv, xmass, ymass, zmass, pr_diag, procID, numLonDomains, &
     &   j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jvlo, jhi, i1, i2, &
     &   ju1, j2, k1, k2, commu_npole, commu_spole, gmi_nborder)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID, numLonDomains, commu_npole, commu_spole
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, j1p, j2p
      integer, intent(in) :: ilo, ihi, julo, jvlo, jhi, gmi_nborder
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      character (len=MAX_LENGTH_MET_NAME), intent(in) :: metdata_name_org
      character (len=1),  intent(in) :: met_grid_type
      real*8, intent(in)  :: geofac_pc
      real*8, intent(in)  :: tdt
      real*8, intent(in)  :: cose    (ju1_gl:j2_gl+1)
      real*8, intent(in)  :: cosp    (ju1_gl:j2_gl)
      real*8, intent(in)  :: geofac  (ju1_gl:j2_gl)
      real*8, intent(in)  :: dap     (k1:k2)
      real*8, intent(in)  :: dbk     (k1:k2)
      real*8, intent(in)  :: rel_area(i1:i2, ju1:j2)
      real*8, intent(in)  :: uu      (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: vv      (ilo:ihi, jvlo:jhi, k1:k2)
!
      real*8, intent(inout) :: pctm1(ilo:ihi, julo:jhi)
      real*8, intent(inout) :: pmet2(ilo:ihi, julo:jhi)
!
      real*8, intent(out) :: dps    (i1:i2, ju1:j2)
      real*8, intent(out) :: dps_ctm(i1:i2, ju1:j2)
      real*8, intent(out) :: xmass  (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: ymass  (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: zmass  (i1:i2, ju1:j2, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
!     --------------------------------------------------------------
!     dpi   : divergence at a grid point; used to calculate vertical
!             motion (mb)
!     --------------------------------------------------------------
!
      real*8  :: dpi(i1:i2, ju1:j2, k1:k2)
!
!     ---------------------------------------------------------------------
!     delpm : pressure thickness, the psudo-density in a hydrostatic system
!             at t1+tdt/2 (approximate) (mb)
!     ---------------------------------------------------------------------
!
      real*8  :: delpm(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) Write (6,*) 'Init_Press_Fix called by ', procID
!
!... average polar cap pressure
!
      if (j1p /= ju1_gl+1) then
!
!       ========================
        call Average_Press_Poles  &
!       ========================
     &    (rel_area, pctm1, &
     &     pr_diag, procID, i1, i2, ju1, j2, ilo, ihi, julo, jhi)
!
!       ========================
        call Average_Press_Poles  &
!       ========================
     &    (rel_area, pmet2, &
     &     pr_diag, procID, i1, i2, ju1, j2, ilo, ihi, julo, jhi)
!
      endif
!
!     -------------------------------------------------------------------
!     We need to calculate pressures at t1+tdt/2.  One ought to use pctm2
!     in the call to Calc_Delpm, but since we don't know it yet, we are
!     forced to use pmet2.  This should be good enough because it is only
!     used to convert the winds to the mass fluxes, which is done crudely
!     anyway and the mass fluxes will still get fixed OK.
!     -------------------------------------------------------------------
!
      dps(i1:i2,ju1:j2) = pmet2(i1:i2,ju1:j2) - pctm1(i1:i2,ju1:j2)
!
      if (metdata_name_org(1:4) == 'GISS') then
!
!       ==============================
        call Calc_Horiz_Mass_Flux_Giss  &
!       ==============================
     &    (tdt, cose, rel_area, uu, vv, xmass, ymass, &
     &     pr_diag, procID, ju1_gl, j2_gl, j1p, j2p, &
     &     ilo, ihi, julo, jvlo, jhi, i1, i2, ju1, j2, k1, k2, gmi_nborder)
!
      elseif (metdata_name_org(1:2) == 'EC') then
!
!       ==============================
        call Calc_Horiz_Mass_Flux_EC  &
!       ==============================
     &    (tdt, cose, rel_area, uu, vv, xmass, ymass, &
     &     pr_diag, procID, ju1_gl, j2_gl, j1p, j2p, &
     &     ilo, ihi, julo, jvlo, jhi, i1, i2, ju1, j2, k1, k2, gmi_nborder)
!
      else
!
!       ===============
        call Calc_Delpm  &
!       ===============
     &    (dap, dbk, pctm1, pmet2, delpm, &
     &     pr_diag, procID, k1, k2, ilo, ihi, julo, jhi)
!
!       =========================
        call Calc_Horiz_Mass_Flux  &
!       =========================
     &    (met_grid_type, tdt, cose, cosp, delpm,  &
     &     uu, vv, xmass, ymass, &
     &     pr_diag, procID, i2_gl, ju1_gl, j2_gl, &
     &     ilo, ihi, julo, jvlo, jhi, i1, i2, ju1, j2, k1, k2, gmi_nborder)
!
      end if
!
!
!     ====================
      call calcDivergence  &
!     ====================
     &  (.false., geofac_pc, geofac, dpi, xmass, ymass, &
     &   pr_diag, procID, numLonDomains, j1p, j2p, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole)
!
!
      dps_ctm(i1:i2,ju1:j2) = Sum (dpi(i1:i2,ju1:j2,:), dim=3)
!
      call calcVerticalMassFlux (dbk, dps_ctm, dpi, zmass, &
     &      pr_diag, procID, i1, i2, ju1, j2, k1, k2)
!
      return
!
      end subroutine Init_Press_Fix
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Press_Fix_Llnl
!
! DESCRIPTION
!   This routine fixes the mass fluxes to match the met field pressure
!   tendency.
!
! ARGUMENTS
!   geofac_pc   : special geometrical factor (geofac) for Polar cap
!   geofac      : geometrical factor for meridional advection; geofac uses
!                 correct spherical geometry, and replaces acosp as the
!                 meridional geometrical factor in tpcore
!   dbk         : difference in bi across layer - the dSigma term
!   dps         : change of surface pressure from met field pressure (mb)
!   dps_ctm     : CTM surface pressure tendency; sum over vertical of dpi
!                 calculated from original mass fluxes  (mb)
!   rel_area    : relative surface area of grid box (fraction)
!   xmass       : horizontal mass flux in E-W direction (mb)
!   ymass       : horizontal mass flux in N-S direction (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Press_Fix_Llnl  &
     &  (geofac_pc, geofac, dbk, dps, dps_ctm, rel_area, xmass, ymass, &
     &   pr_diag, procID, numLonDomains, j1p, j2p, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole, gmi_nborder)
!
      use GmiWrapMaster_mod, only : wrapMaster_3du
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID, numLonDomains, commu_npole, commu_spole
      integer, intent(in) :: j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilo, ihi, julo, jhi, gmi_nborder
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      real*8, intent(in)  :: geofac_pc
      real*8, intent(in)  :: geofac  (ju1_gl:j2_gl)
      real*8, intent(in)  :: dbk     (k1:k2)
      real*8, intent(in)  :: dps     (i1:i2, ju1:j2)
      real*8, intent(in)  :: rel_area(i1:i2, ju1:j2)
!
      real*8, intent(inout) :: dps_ctm(i1:i2, ju1:j2)
      real*8, intent(inout) :: xmass  (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(inout) :: ymass  (ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik
!
      real*8  :: dgpress
      real*8  :: fxmean
      real*8  :: ri2
!
      real*8  :: fxintegral(i1:i2+1)
!
      real*8  :: mmfd(ju1:j2)
      real*8  :: mmf (ju1:j2)
!
      real*8  :: ddps(i1:i2, ju1:j2)
!
!     ------------------------------------------------------------------------
!     dpi : divergence at a grid point; used to calculate vertical motion (mb)
!     ------------------------------------------------------------------------
!
      real*8  :: dpi(i1:i2, ju1:j2, k1:k2)
!
      real*8  :: xcolmass_fix(ilo:ihi, julo:jhi)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Do_Press_Fix_Llnl called by ', procID
      end if
!
!
      ri2 = i2_gl
!
      mmfd(:) = 0.0d0
!
      xcolmass_fix(:,:)   = 0.0d0
!
!
!     ------------------------------------------------------------
!     Calculate difference between GCM and LR predicted pressures.
!     ------------------------------------------------------------
!
      ddps(:,:) = dps(:,:) - dps_ctm(:,:)
!
!
!     --------------------------------------
!     Calculate global-pressure discrepancy.
!     --------------------------------------
!
      dgpress =  &
     &  Sum (ddps(i1:i2,ju1:j2) * rel_area(i1:i2,ju1:j2))
!
!
!     ----------------------------------------------------------
!     Calculate mean meridional flux divergence (df/dy).
!     Note that mmfd is actually the zonal mean pressure change,
!     which is related to df/dy by geometrical factors.
!     ----------------------------------------------------------
!
!     ------------------------
!     Handle non-Pole regions.
!     ------------------------
!
      do ij = j1p, j2p
        mmfd(ij) = -(sum(ddps(:,ij)) / ri2 - dgpress)
      end do
!
!     ---------------------------------------------
!     Handle poles.
!     Note that polar cap boxes have to all be averaged.
!     ---------------------------------------------
!
       do ij = ju1, j1p-1
         mmfd(ij) = -(ddps(1,ij) - dgpress)
       enddo
       do ij = j2p+1, j2
         mmfd(ij) = -(ddps(1,ij) - dgpress)
       enddo
!
!
!     ---------------------------------------------
!     Calculate mean meridional fluxes (cos(e)*fy).
!     ---------------------------------------------
!
      mmf(j1p) = mmfd(ju1) / geofac_pc
!
      do ij = j1p, j2p
        mmf(ij+1) = mmf(ij) + mmfd(ij) / geofac(ij)
      end do
!
!
!     ------------------------------------------------------------
!     Fix latitude bands.
!     Note that we don't need to worry about geometry here because
!     all boxes in a latitude band are identical.
!     Note also that fxintegral(i2+1) should equal fxintegral(i1),
!     i.e., zero.
!     ------------------------------------------------------------
!
      do ij = j1p, j2p
!
        fxintegral(:) = 0.0d0
!
        do il = i1, i2
!
          fxintegral(il+1) =  &
     &      fxintegral(il) -  &
     &      (ddps(il,ij) - dgpress) -  &
     &      mmfd(ij)
!
        end do
!
        fxmean = Sum (fxintegral(i1+1:i2+1)) / ri2
!
        do il = i1, i2+1
          xcolmass_fix(il,ij) = fxintegral(il) - fxmean
        end do
!
      end do
!
!
!     -------------------------------------
!     Distribute colmass_fix's in vertical.
!     -------------------------------------
!
      do ik = k1, k2
        do ij = j1p, j2p
          do il = i1, i2+1
!
            xmass(il,ij,ik) = xmass(il,ij,ik) +  &
     &                        xcolmass_fix(il,ij) * dbk(ik)
!
          end do
        end do
      end do
!
      do ik = k1, k2
        do ij = j1p, j2p+1
          do il = i1, i2
!
            ymass(il,ij,ik) = ymass(il,ij,ik) +  &
     &                        mmf(ij) * dbk(ik)
!
          end do
        end do
      end do
!
!
!     ------------------------------------------------
!     Always done on Master, so only "wrap" necessary.
!     ------------------------------------------------
!
      call wrapMaster_3du (xmass, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, gmi_nborder)
      call wrapMaster_3du (ymass, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, gmi_nborder)
!
!
!     ====================
      call calcDivergence  &
!     ====================
     &  (.false., geofac_pc, geofac, dpi, xmass, ymass, &
     &   pr_diag, procID, numLonDomains, j1p, j2p, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole)
!
!
      dps_ctm(i1:i2,ju1:j2) = Sum (dpi(i1:i2,ju1:j2,:), dim=3)
!
!
      return
!
      end subroutine Do_Press_Fix_Llnl
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Average_Press_Poles
!
! DESCRIPTION
!   This routine averages pressure at the Poles when the Polar cap is
!   enlarged.  It makes the last two latitudes equal.
!
! ARGUMENTS
!   rel_area : relative surface area of grid box (fraction)
!   press    : surface pressure (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Average_Press_Poles  &
     &  (rel_area, press, &
     &   pr_diag, procID, i1, i2, ju1, j2, ilo, ihi, julo, jhi)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: i1, i2, ju1, j2, ilo, ihi, julo, jhi
      real*8, intent(in)  :: rel_area(i1:i2, ju1:j2)
!
      real*8, intent(inout) :: press(ilo:ihi, julo:jhi)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      real*8  :: meanp
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Average_Press_Poles called by ', procID
      end if
!
      meanp =  &
     &  Sum (rel_area(i1:i2,ju1:ju1+1)  *  &
     &       press   (i1:i2,ju1:ju1+1)) /  &
     &  Sum (rel_area(i1:i2,ju1:ju1+1))
!
      press(i1:i2,ju1:ju1+1) = meanp
!
      meanp =  &
     &  Sum (rel_area(i1:i2,j2-1:j2)  *  &
     &       press   (i1:i2,j2-1:j2)) /  &
     &  Sum (rel_area(i1:i2,j2-1:j2))
!
      press(i1:i2,j2-1:j2) = meanp
!
      return
!
      end subroutine Average_Press_Poles
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calc_Horiz_Mass_Flux
!
! DESCRIPTION
!   This routine converts winds on A or C grid to mass flux on the C grid.
!
! ARGUMENTS
!   met_grid_type : met grid type, 'A' or 'C'
!   tdt   : model time step (s)
!   cose  : cosine of grid box edges
!   cosp  : cosine of grid box centers
!   delpm : pressure thickness, the psudo-density in a hydrostatic system
!           at t1+tdt/2 (approximate) (mb)
!   uu    : wind velocity  in E-W direction at t1+tdt/2 (m/s)
!   vv    : wind velocity  in N-S direction at t1+tdt/2 (m/s)
!   xmass : horizontal mass flux in E-W direction (mb)
!   ymass : horizontal mass flux in N-S direction (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Calc_Horiz_Mass_Flux  &
     &    (met_grid_type, tdt, cose, cosp, delpm,  &
     &     uu, vv, xmass, ymass, &
     &     pr_diag, procID, i2_gl, ju1_gl, j2_gl, &
     &     ilo, ihi, julo, jvlo, jhi, i1, i2, ju1, j2, k1, k2, gmi_nborder)
!
      use GmiWrapMaster_mod, only : wrapMaster_3du
!
      implicit none
!
#     include "gmi_phys_constants.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: i2_gl, ju1_gl, j2_gl, gmi_nborder
      integer, intent(in) :: ilo, ihi, julo, jvlo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      character (len=1), intent(in) :: met_grid_type
      real*8, intent(in)  :: tdt
      real*8, intent(in)  :: cose (ju1_gl:j2_gl+1)
      real*8, intent(in)  :: cosp (ju1_gl:j2_gl)
      real*8, intent(in)  :: delpm(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: uu   (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: vv   (ilo:ihi, jvlo:jhi, k1:k2)
!
      real*8, intent(out) :: xmass(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: ymass(ilo:ihi, julo:jhi, k1:k2)
!
!
!     -----------------------
!     Parameter declarations.
!     -----------------------
!
!     ----------------------------------------------------------------------
!     INTERP_MASS_FLUX_COMPONENTS :
!
!        If true, then interpolate winds and delta pressure separately from
!        A-grid to C-grid before multiplying to get mass flux on the C-grid.
!
!        If false, then multiply winds and delta pressure to get mass flux
!        on the A-grid, then interpolate mass flux on to the C-grid.
!     ----------------------------------------------------------------------
!
      logical, parameter :: INTERP_MASS_FLUX_COMPONENTS = .TRUE.
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      logical, save :: first = .true.
!
      integer :: il, ij, ik
!
!     -------------------------------
!     dl : spacing in longitude (rad)
!     dp : spacing in latitude  (rad)
!     -------------------------------
!
      real*8  :: dl
      real*8  :: dp
!
      real*8  :: ri2
      real*8  :: rj2m1
!
!     ------------------------
!     dtdy  : dt/dy      (s/m)
!     dtdy5 : 0.5 * dtdy (s/m)
!     ------------------------
!
      real*8, allocatable, save :: dtdy(:)
      real*8, allocatable, save :: dtdy5(:)
!
!     ------------------------
!     dtdx  : dt/dx      (s/m)
!     dtdx5 : 0.5 * dtdx (s/m)
!     ------------------------
!
      real*8, allocatable, save :: dtdx (:)
      real*8, allocatable, save :: dtdx5(:)
!
!     ---------------------------------------------------------------
!     tmp_pu1 : layer thickness on eastern edge of gridbox (mbar)
!     tmp_pu2 : layer thickness on southern edge of gridbox (mbar)
!     tmp_crx : zonal Courant number on eastern edge of gridbox
!     tmp_cry : meridional Courant number on southern edge of gridbox
!     ---------------------------------------------------------------
!
      real*8  :: tmp_pu1
      real*8  :: tmp_pu2
      real*8  :: tmp_crx
      real*8  :: tmp_cry
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Calc_Horiz_Mass_Flux called by ', procID
      end if
!
!
!     ==========
      if (first) then
!     ==========
!
        first = .false.
!
!
        Allocate (dtdx (ju1_gl:j2_gl))
        Allocate (dtdx5(ju1_gl:j2_gl))
        dtdx = 0.0d0
        dtdx5 = 0.0d0
!
        Allocate (dtdy (ju1_gl:j2_gl))
        Allocate (dtdy5(ju1_gl:j2_gl))
        dtdy = 0.0d0
        dtdy5 = 0.0d0
!
!
        ri2   = i2_gl
        rj2m1 = j2_gl - 1
!
        dl    = 2.0d0 * GMI_PI / ri2
        dp    = GMI_PI / rj2m1
!
!
!
!
!        dtdx (ju1_gl) = 0.0d0
!        dtdx5(ju1_gl) = 0.0d0
!
!        do ij = ju1_gl + 1, j2_gl - 1
        do ij = ju1_gl, j2_gl
          dtdx (ij) = tdt / (dl * RADEAR * cosp(ij))
          dtdx5(ij) = 0.5d0 * dtdx(ij)
          dtdy(ij)  = tdt / (RADEAR * dp)
          dtdy5(ij) = 0.5d0 * dtdy(ij)
        end do
!
!        dtdx (j2_gl)  = 0.0d0
!        dtdx5(j2_gl)  = 0.0d0
!... polar grid box is 1/2 size of rest of grid
        dtdy (ju1_gl) = dtdy (ju1_gl)*2
        dtdy5(ju1_gl) = dtdy5(ju1_gl)*2
        dtdy (j2_gl)  = dtdy (j2_gl)*2
        dtdy5(j2_gl)  = dtdy5(j2_gl)*2
!
!
      end if
!
!
      xmass(:,:,:) = 0.0d0
      ymass(:,:,:) = 0.0d0
!
!
      if (INTERP_MASS_FLUX_COMPONENTS .or. (met_grid_type == 'C')) then
!
!
        do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2
!
!             =========================
              if (met_grid_type == 'A') then  ! A grid
!             =========================
!
                tmp_crx = dtdx5(ij) * (uu(il,ij,ik) + uu(il-1,ij, ik))
!
                tmp_cry = dtdy5(ij) * (vv(il,ij,ik) + vv(il,  ij-1,ik))
!
!             ====
              else  ! C grid
!             ====
!
                tmp_crx = dtdx(ij) * uu(il-1,ij,  ik)
!.sds.. is this correct for C-grid met fields? Why shift by one?
                tmp_cry = dtdy(ij) * vv(il,  ij-1,ik)
!
              end if
!
!             -----------------------------------
!             Calculate E-W horizontal mass flux.
!             -----------------------------------
!
              tmp_pu1 = 0.5d0 * (delpm(il,ij,ik) + delpm(il-1,ij,ik))
!
              xmass(il,ij,ik) = tmp_crx * tmp_pu1
!
!
!             -----------------------------------
!             Calculate N-S horizontal mass flux.
!             -----------------------------------
!
              if(ij == ju1_gl) then
                tmp_pu2 = 0.0d0
              else
                tmp_pu2 = 0.5d0 * (delpm(il,ij,ik) + delpm(il,ij-1,ik))
              endif
!
              ymass(il,ij,ik) = tmp_cry * tmp_pu2 * cose(ij)
!
            end do
          end do
        end do
!
      else
!
        do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2
!
              xmass(il, ij, ik) = dtdx5(ij) *  &
     &                 ( uu(il-1,ij,ik) * delpm(il-1,ij,ik) +  &
     &                   uu(il  ,ij,ik) * delpm(il  ,ij,ik)  )
!
!.sds.. is this correct?? zero flux@pole by definition
              if(ij == ju1_gl) then
                ymass(il, ij, ik) = 0.0d0
              else
                ymass(il, ij, ik) = dtdy5(ij) * cose(ij) *  &
     &                 ( vv(il,ij-1,ik) * delpm(il,ij-1,ik) +  &
     &                   vv(il,ij  ,ik) * delpm(il,ij  ,ik)  )
              endif
!
            end do
          end do
        end do
!
      end if
!
!
!     ------------------------------------------------
!     Always done on Master, so only "wrap" necessary.
!     ------------------------------------------------
!
      call wrapMaster_3du (xmass, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, gmi_nborder)
      call wrapMaster_3du (ymass, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, gmi_nborder)
!
!
      return
!
      end subroutine Calc_Horiz_Mass_Flux
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calc_Horiz_Mass_Flux_Giss
!
! DESCRIPTION
!   This routine calculates the horizontal mass flux for GISS met data.
!
! ARGUMENTS
!   tdt      : model time step (s)
!   cose     : cosine of grid box edges
!   rel_area : relative surface area of grid box (fraction)
!   uu       : wind velocity in E-W direction (100N/s)
!   vv       : wind velocity in N-S direction (100N/s)
!   xmass    : horizontal mass flux in E-W direction (mb)
!   ymass    : horizontal mass flux in N-S direction (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Calc_Horiz_Mass_Flux_Giss  &
     &  (tdt, cose, rel_area, uu, vv, xmass, ymass, &
     &     pr_diag, procID, ju1_gl, j2_gl, j1p, j2p, &
     &     ilo, ihi, julo, jvlo, jhi, i1, i2, ju1, j2, k1, k2, gmi_nborder)
!
      use GmiWrapMaster_mod, only : wrapMaster_3du
!
      implicit none
!
#     include "gmi_phys_constants.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: ju1_gl, j2_gl, j1p, j2p, gmi_nborder
      integer, intent(in) :: ilo, ihi, julo, jvlo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
!
      real*8, intent(in)  :: tdt
      real*8, intent(in)  :: cose(ju1_gl:j2_gl+1)
      real*8, intent(in)  :: rel_area(i1:i2, ju1:j2)
      real*8, intent(in)  :: uu   (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: vv   (ilo:ihi, jvlo:jhi, k1:k2)
      real*8, intent(out) :: xmass(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: ymass(ilo:ihi, julo:jhi, k1:k2)
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
        Write (6,*) 'Calc_Horiz_Mass_Flux_Giss called by ', procID
      end if
!
!
!     -----------------------------------
!     Calculate E-W horizontal mass flux.
!     -----------------------------------
!
      call Deter_Jrange (ju1, j2, j1p, j2p, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
      do ij = jst, jend
        do il = i1, i2
!
          xmass(il,ij,:) =  &
     &      (uu(il,ij,:) * tdt) /  &
     &      (rel_area(il,ij) * SAREA_EARTH)
!
        end do
      end do
!
!
      call wrapMaster_3du (xmass, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, gmi_nborder)
!
!
!     -----------------------------------
!     Calculate N-S horizontal mass flux.
!     -----------------------------------
!
      call Deter_Jrange (ju1, j2+1, j1p, j2p+1, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!
      do ij = jst, jend
        do il = i1, i2
!
          ymass(il,ij,:) =  &
     &      (cose(ij) * vv(il,ij,:) * tdt) /  &
     &      (rel_area(il,ij) * SAREA_EARTH)
!
        end do
      end do
!
      call wrapMaster_3du (ymass, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, gmi_nborder)
!
      return
!
      end subroutine Calc_Horiz_Mass_Flux_Giss
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calc_Horiz_Mass_Flux_EC
!
! DESCRIPTION
!   This routine calculates the horizontal mass flux for EC-Oslo met data.
!
! ARGUMENTS
!   tdt      : model time step (s)
!   cose     : cosine of grid box edges
!   rel_area : relative surface area of grid box (fraction)
!   uu       : wind velocity in E-W direction (kg/s)
!   vv       : wind velocity in N-S direction (kg/s)
!   xmass    : horizontal mass flux in E-W direction (mb)
!   ymass    : horizontal mass flux in N-S direction (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Calc_Horiz_Mass_Flux_EC  &
     &  (tdt, cose, rel_area, uu, vv, xmass, ymass, &
     &     pr_diag, procID, ju1_gl, j2_gl, j1p, j2p, &
     &     ilo, ihi, julo, jvlo, jhi, i1, i2, ju1, j2, k1, k2, gmi_nborder)
!
      use GmiWrapMaster_mod, only : wrapMaster_3du
!
      implicit none
!
#     include "gmi_phys_constants.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: ju1_gl, j2_gl, j1p, j2p, gmi_nborder
      integer, intent(in) :: ilo, ihi, julo, jvlo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
!
      real*8, intent(in)  :: tdt
      real*8, intent(in)  :: cose(ju1_gl:j2_gl+1)
      real*8, intent(in)  :: rel_area(i1:i2, ju1:j2)
      real*8, intent(in)  :: uu   (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: vv   (ilo:ihi, jvlo:jhi, k1:k2)
      real*8, intent(out) :: xmass(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(out) :: ymass(ilo:ihi, julo:jhi, k1:k2)
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
        Write (6,*) 'Calc_Horiz_Mass_Flux_EC called by ', procID
      end if
!
!
!     -----------------------------------
!     Calculate E-W horizontal mass flux.
!     -----------------------------------
!
      call Deter_Jrange (ju1, j2, j1p, j2p, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
      do ij = jst, jend
        do il = i1, i2
!
          xmass(il,ij,:) =  &
     &      (uu(il,ij,:) * tdt) /  &
     &      (rel_area(il,ij) * SAREA_EARTH)
!
        end do
      end do
!
!
      call wrapMaster_3du (xmass, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, gmi_nborder)
!
!
!     -----------------------------------
!     Calculate N-S horizontal mass flux.
!     -----------------------------------
!
      call Deter_Jrange (ju1, j2+1, j1p, j2p+1, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!
      do ij = jst, jend
        do il = i1, i2
!
          ymass(il,ij,:) =  &
     &      (cose(ij) * vv(il,ij,:) * tdt) /  &
     &      (rel_area(il,ij) * SAREA_EARTH)
!
        end do
      end do
!
      call wrapMaster_3du (ymass, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, gmi_nborder)
!
      return
!
      end subroutine Calc_Horiz_Mass_Flux_EC
!
!
!=============================================================================
!
! CODE DEVELOPER
!   Peter Connell, LLNL (Original code from Michael Prather, UCI)
!     connell2@llnl.gov
!   Additional modifications by Philip Cameron-Smith, LLNL
!     cameronsmith1@llnl.gov
!   Additional modifications by John Tannahill, LLNL
!     jrt@llnl.gov
!
! FILE
!   gmi_press_fixer_uci.F
!
! ROUTINES
!   Do_Press_Fix_Uci
!   Do_Press_Filt_Uci
!   Do_Local_Filt_Uci
!   Do_Pole_Filt_Uci
!
!=============================================================================
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Press_Fix_Uci
!
! DESCRIPTION
!   This is an adaptation of M. Prather's UCI dyn0 filter for use in Gmimod.
!   The purpose of this code is to filter u and v mass fluxes and surface
!   pressures to generate a self-consistent set.
!
!   Adapted from generic CTM shell from UC Irvine (p-code 4.2, 9/01/00).
!
!   Horizontal pressure filter adjusts u & v to reduce error in (pctm - ps)
!     u + pressure filter ==> u#, v + filter ==> v#  (temporary)
!     The pressure filter does nearest neighbor flux (adjusting alfa/beta).
!
! ARGUMENTS
!   geofac_pc   : special geometrical factor (geofac) for Polar cap
!   geofac      : geometrical factor for meridional advection; geofac uses
!                 correct spherical geometry, and replaces acosp as the
!                 meridional geometrical factor in tpcore
!   dbk         : difference in bi across layer - the dSigma term
!   dps         : change of surface pressure from met field pressure (mb)
!   dps_ctm     : CTM surface pressure tendency; sum over vertical of dpi
!                 calculated from original mass fluxes  (mb)
!   rel_area    : relative surface area of grid box (fraction)
!   xmass       : horizontal mass flux in E-W direction (mb)
!   ymass       : horizontal mass flux in N-S direction (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Press_Fix_Uci  &
     &  (geofac_pc, geofac, dbk, dps, dps_ctm, rel_area, xmass, ymass, &
     &   pr_diag, procID, numLonDomains, ilong, ilat, j1p, j2p, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole, gmi_nborder)
!
      use GmiWrapMaster_mod, only : wrapMaster_3du
!
      implicit none
!
#     include "gmi_phys_constants.h"
!
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID, numLonDomains, commu_npole, commu_spole
      integer, intent(in) :: ilong, ilat
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, j1p, j2p
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, gmi_nborder
      real*8  :: geofac_pc
      real*8  :: geofac(ju1_gl:j2_gl)
      real*8  :: dbk   (k1:k2)
      real*8  :: dps   (i1:i2, ju1:j2)
      real*8  :: dps_ctm (i1:i2, ju1:j2)
      real*8  :: rel_area(i1:i2, ju1:j2)
      real*8  :: xmass   (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: ymass   (ilo:ihi, julo:jhi, k1:k2)
!
!     -----------------------
!     Parameter declarations.
!     -----------------------
!
      integer, parameter :: NTIMES_PARAM = 5
!
      real*8,  parameter :: G100  = 100.0d0 / GMI_G
!
      real*8,  parameter :: RG100 = 1.0d0 / G100
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik
      integer :: ntimes
!
!     --------------------------------
!     ax : perturbed alfa by merr (kg)
!     bx : perturbed beta by merr (kg)
!     --------------------------------
!
      real*8  :: ax(ilong+1, ilat)
      real*8  :: bx(ilong,   ilat+1)
!
!     -----------------------------------------------------------------
!     merr : mass error
!     perr : pressure error between CTM-GCM at new time (before filter)
!     -----------------------------------------------------------------
!
      real*8  :: merr(ilong, ilat)
      real*8  :: perr(ilong, ilat)
!
!     ------------------------------------------------------------
!     dpi : divergence at a grid point; used to calculate vertical
!           motion (mb)
!     ------------------------------------------------------------
!
      real*8  :: dpi(i1:i2, ju1:j2, k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Do_Press_Fix_Uci called.'
      end if
!
!
!     -----------------------------------------------------------
!     Begin filter of pressure errors.
!
!     Define the error in surface pressure, perr, expected at end
!     of time step.
!
!     Filter by error in adjacent boxes, weight by areas, adjust
!     ax & bx.
!     -----------------------------------------------------------
!
      perr(i1:i2,ju1:j2) =  &
     &  dps_ctm(i1:i2,ju1:j2) - dps(i1:i2,ju1:j2)
!
      merr(i1:i2,ju1:j2) =  &
     &  perr(i1:i2,ju1:j2) * rel_area(i1:i2,ju1:j2) * G100
!
      ntimes = NTIMES_PARAM
!
!     ======================
      call Do_Press_Filt_Uci  &
!     ======================
     &  (ntimes, merr, ax, bx, rel_area, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl, ilat, ilong)
!
!
!     ----------------------------------------
!     Distribute column mass flux in vertical.
!     ----------------------------------------
!
      do ij = ju1, j2
        do il = i1, i2
          ax(il,ij) = ax(il,ij) * RG100 / rel_area(il,ij)
        end do
        ax(i2+1,ij) = ax(i2+1,ij) * RG100 / rel_area(1,ij)
      end do
!
!
      do ik = k1, k2
        do ij = ju1, j2
          do il = i1, i2+1
!
            xmass(il,ij,ik) = xmass(il,ij,ik) +  &
     &                        (ax(il,ij) * dbk(ik))
!
          end do
        end do
      end do
!
!
      do ij = ju1, j2+1
        do il = i1, i2
!
          if (ij /= j2+1) then
!
            bx(il,ij) = bx(il,ij) * RG100 / rel_area(il,ij)
            bx(il,ij) = bx(il,ij) / geofac(ij)
!
          else
!
            bx(il,ij) = bx(il,ij) * RG100 / rel_area(il,ij-1)
            bx(il,ij) = bx(il,ij) / geofac(ij-1)
!
          end if
!
        end do
      end do
!
!
      do ik = k1, k2
        do ij = ju1, j2+1
          do il = i1, i2
!
            ymass(il,ij,ik) = ymass(il,ij,ik) +  &
     &                        (bx(il,ij) * dbk(ik))
!
          end do
        end do
      end do
!
!
      call wrapMaster_3du (xmass, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, gmi_nborder)
      call wrapMaster_3du (ymass, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, gmi_nborder)
!
!
!     ====================
      call calcDivergence  &
!     ====================
     &  (.false., geofac_pc, geofac, dpi, xmass, ymass, &
     &   pr_diag, procID, numLonDomains, j1p, j2p, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole)
!
!
      dps_ctm(i1:i2,ju1:j2) = Sum (dpi(i1:i2,ju1:j2,:), dim=3)
!
!
      return
!
      end subroutine Do_Press_Fix_Uci
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Press_Filt_Uci
!
! DESCRIPTION
!   This routine calculates coefficients ax and bx for corrections to alfa
!   and beta, the x and y mass flux.
!
!   Filter is smoothing presssure error, the error in pressure between
!   predicted ps(ctm) and ps(gcm).  The filtering can be done in multiple
!   steps:  local (del^2), longitude stripes, latitude belts.
!
!   Currently we only do local (del^2) filtering.
!
! ARGUMENTS
!   ntimes : number of iterations
!   merr   : mass error
!   ax     : perturbed alfa by merr (kg)
!   bx     : perturbed beta by merr (kg)
!   areaxy : area of grid box (m^2)
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Press_Filt_Uci  &
     &  (ntimes, merr, ax, bx, areaxy, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl, ilat, ilong)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilat, ilong
      integer :: ntimes
      real*8  :: merr(ilong,   ilat)
      real*8  :: ax  (ilong+1, ilat)
      real*8  :: bx  (ilong,   ilat+1)
      real*8  :: areaxy(i1_gl:i2_gl, ju1_gl:j2_gl)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      real*8  :: x0(ilong, ilat)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      ax(:,:) = 0.0d0
      bx(:,:) = 0.0d0
!
      x0(:,:) = merr(:,:)
!
!
!     ======================
      call Do_Local_Filt_Uci  &
!     ======================
     &  (ntimes, merr, ax, bx, areaxy, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl, ilat, ilong)
!
!
!     =====================
      call Do_Pole_Filt_Uci  &
!     =====================
     &  (1.0d0, merr, bx, areaxy, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl, ilat, ilong)
!
!
!     ------------------------------------------------------------------
!     It should be superfluous to recalculate merr here, but leave it in
!     anyway just in case.
!     ------------------------------------------------------------------
!
      merr(:,:) =  &
     &  x0(1:ilong,1:ilat) + ax(1:ilong,:) - ax(2:ilong+1,:) +  &
     &                       bx(:,1:ilat)  - bx(:,2:ilat+1)
!
!
      return
!
      end subroutine Do_Press_Filt_Uci
!
!
!-----------------------------------------------------------------------
!
! ROUTINE
!   Do_Local_Filt_Uci
!
! DESCRIPTION
!   This routine does the local filtering.
!
! ARGUMENTS
!   ntimes : number of iterations
!   merr   : mass error (kg)
!   ax     : perturbed alfa by merr (kg)
!   bx     : perturbed beta by merr (kg)
!   areaxy : area of grid box (m^2)
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Local_Filt_Uci  &
     &  (ntimes, merr, ax, bx, areaxy, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl, ilat, ilong)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilat, ilong
      integer :: ntimes
      real*8  :: merr(ilong,   ilat)
      real*8  :: ax  (ilong+1, ilat)
      real*8  :: bx  (ilong,   ilat+1)
      real*8  :: areaxy(i1_gl:i2_gl, ju1_gl:j2_gl)
!
!
!     -----------------------
!     Parameter declarations.
!     -----------------------
!
      real*8, parameter :: FNAZ8_INIT = 0.125d0
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: ij
      integer :: nfltr
!
      real*8  :: fnaz8
!
      real*8  :: x0(ilong, ilat)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
!     ------------------------------------------------------------------
!     Initialize corrective column mass flows (kg):  ax->alfa, bx->beta.
!     ------------------------------------------------------------------
!
      x0(:,:) = merr(:,:)
!
!
!     ------------------------------------------------------------------
!     Iterate over mass-error filter, accumulate corrections in ax & bx.
!     ------------------------------------------------------------------
!
      do nfltr = 1, ntimes
!
!       --------------------------
!       Calculate ax = E-W filter.
!       --------------------------
!
!       -----------------------------------------------------------
!       Calculate pressure filter E-W wind between boxes ij-1 & ij.
!       -----------------------------------------------------------
!
        fnaz8 = FNAZ8_INIT
!
        ax(2:ilong,2:ilat-1) =  &
     &    ax(2:ilong,2:ilat-1) +  &
     &    (fnaz8 * (merr(1:ilong-1,2:ilat-1) - merr(2:ilong,2:ilat-1)))
!
!
!       --------------------------------------------------------------
!       Calculate pressure filter E-W wind at edges ij=1 & ij=ilong+1.
!       --------------------------------------------------------------
!
        ax(ilong+1,2:ilat-1) =  &
     &    ax(ilong+1,2:ilat-1) +  &
     &    (fnaz8 * (merr(ilong,2:ilat-1) - merr(1,2:ilat-1)))
!
        ax(1,2:ilat-1) =  &
     &    ax(1,2:ilat-1) +  &
     &    (fnaz8 * (merr(ilong,2:ilat-1) - merr(1,2:ilat-1)))
!
!
!       ------------------------------------------------------------
!       Calculate bx = N-S filter, N-S wind between boxes ij-1 & ij.
!       ------------------------------------------------------------
!
        do ij = 3, ilat-1
!
          fnaz8 =  &
     &      (0.25d0 * areaxy(1,ij)) / (areaxy(1,ij-1) + areaxy(1,ij))
!
          bx(:,ij) =  &
     &      bx(:,ij) +  &
     &      (fnaz8 * (merr(:,ij-1) - merr(:,ij)))
!
        end do
!
!
!       -----------------------------------------------------------------
!       Enhance the filtering by factor of 2 ONLY into/out-of polar caps.
!       -----------------------------------------------------------------
!
        fnaz8 =  &
     &    (0.5d0 * areaxy(1,2)) / (areaxy(1,1) + areaxy(1,2))
!
        bx(:,2) =  &
     &    bx(:,2) +  &
     &    (fnaz8 * (merr(:,1) - merr(:,2)))
!
        fnaz8 =  &
     &    (0.5d0 * areaxy(1,ilat)) / (areaxy(1,ilat-1) + areaxy(1,ilat))
!
        bx(:,ilat) =  &
     &    bx(:,ilat) +  &
     &    (fnaz8 * (merr(:,ilat-1) - merr(:,ilat)))
!
!
!       ---------------------------------------------------------------
!       Need N-S flux across boundaries if window calc. (assume merr=0
!       outside).  ilat for optimal matrix/looping; it would be best to
!       define merr=0 for an oversized array merr(0:ilong+1,0:ilat+1).
!       ---------------------------------------------------------------
!
        merr(:,:) =  &
     &    x0(1:ilong,1:ilat) +  &
     &    ax(1:ilong,:) - ax(2:ilong+1,:) +  &
     &    bx(:,1:ilat)  - bx(:,2:ilat+1)
!
      end do
!
!
      return
!
      end subroutine Do_Local_Filt_Uci
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Pole_Filt_Uci
!
! DESCRIPTION
!   This routine does the Pole filtering.
!
! ARGUMENTS
!   coef   : correction scale factor
!   merr   : mass error (kg)
!   bx     : perturbed beta by merr (kg)
!   areaxy : area of grid box (m^2)
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Pole_Filt_Uci  &
     &  (coef, merr, bx, areaxy, &
     &  i1_gl, i2_gl, ju1_gl, j2_gl, ilat, ilong)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilat, ilong
!
      real*8  :: coef
      real*8  :: merr(ilong, ilat)
      real*8  :: bx  (ilong, ilat+1)
      real*8  :: areaxy(i1_gl:i2_gl, ju1_gl:j2_gl)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij
!
      real*8  :: erav
      real*8  :: sum
!
      real*8  :: bxj(ilat+1)
!
!
!     ----------------
!     Begin Execution.
!     ----------------
!
      bxj(:) = 0.0d0
!
!
!     ---------------------------------------------------------
!     Initialize corrective column mass flows (kg):  bxj->beta.
!     ---------------------------------------------------------
!
      do il = 1, ilong
!
        erav = 0.0d0
        sum  = 0.0d0
!
        do ij = 1, ilat
!
          erav = erav + merr (il,ij)
          sum  = sum  + areaxy(1,ij)
!
        end do
!
        erav = erav / sum
!
!
!       ------------------------------------------
!       Mass error filter, make corrections in bx.
!       ------------------------------------------
!
        do ij = 2, ilat
!
          bxj(ij) =  &
     &      bxj(ij-1) + merr(il,ij-1) -  &
     &      (areaxy(1,ij-1) * erav)
!
          bx(il,ij) =  &
     &      bx(il,ij) + (coef * bxj(ij))
!
        end do
!
        merr(il,:) = merr(il,:) +  &
     &               coef * (bxj(1:ilat) - bxj(2:ilat+1))
!
      end do
!
!
      return
!
      end subroutine Do_Pole_Filt_Uci
!
      end module GmiPressureFixer_mod
!
