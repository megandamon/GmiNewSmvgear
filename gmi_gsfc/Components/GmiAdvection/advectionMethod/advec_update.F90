!
!=============================================================================
!
! $Id: advec_update.F90,v 1.19 2013-07-31 14:59:27 ssteenro Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Shian-Jiann Lin, DAO)
!     jrt@llnl.gov
!   Additional modifications by Philip Cameron-Smith, LLNL
!     cameronsmith1@llnl.gov
!
! FILE
!   advec_update.F
!
! ROUTINES
!   Update_Advec
!   Average_Const_Poles
!   Write_Advec_Diag
!   Reset_Flux
!   Store_CFlux
!   Store_MZFlux
!   Store_MXYFlux
!
!=============================================================================
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update_Advec
!
! DESCRIPTION
!   This routine updates the advection.
!
! ARGUMENTS
!   do_chem_grp        : do_chemical groups?
!   num_adv_time_steps : number of advection time steps
!   advec_consrv_opt   : advection conserve option
!   advec_flag : array of flags that indicate whether or not to advect a
!                given species
!   pt         : pressure = (ai * pt) + (bi * psf) (mb)
!   geofac_pc  : special geometrical factor (geofac) for Polar cap
!   geofac     : geometrical factor for meridional advection; geofac uses
!                correct spherical geometry, and replaces acosp as the
!                meridional geometrical factor in tpcore
!   cose       : cosine of grid box edges
!   mcor       : area of grid box (m^2)
!   ai         : pressure = (ai * pt) + (bi * psf), ai at zone interface
!   bi         : pressure = (ai * pt) + (bi * psf), bi at zone interface
!   rel_area   : relative surface area of grid box (fraction)
!   pctm1      : CTM surface pressure at t1     (mb)
!   pctm2      : CTM surface pressure at t1+tdt (mb)
!   dap        : pressure difference across layer from (ai * pt) term (mb)
!   dbk        : difference in bi across layer - the dSigma term
!   xmass      : horizontal mass flux in E-W direction (mb)
!   ymass      : horizontal mass flux in N-S direction (mb)
!   const      : species concentration, known at zone centers (mixing ratio)
!
!-----------------------------------------------------------------------------
!
      subroutine Update_Advec  &
     &  (do_chem_grp, num_adv_time_steps, advec_consrv_opt, &
     &   advec_flag, pt, geofac_pc, geofac, cose, mcor, ai, bi, &
     &   rel_area, pctm1, pctm2, dap, dbk, xmass, ymass, const, &
     &   air_mass_flux, flux_x, flux_y, flux_z, psf_flux, &
     &   do_mean, pr_flux, pr_psf_flux, count_flux, &
     &   pr_diag, loc_proc, numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &   flux_species_num, flux_species,  mw, pr_const_flux, &
     &   j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &   i1, i2, ju1, j2, k1, k2, ilong, ivert, numSpecies, &
     &   commu_npole, commu_spole, ewflag, nspoleu, nb, sb, eb, wb,  &
     &   north_domain, south_domain, east_domain, west_domain,  &
     &   communicatorWorld, mapi_all)
!
      use GmiSpcConcentrationMethod_mod, only : isFixedConcentration
      use GmiSubDomainsBC_mod          , only : Gmi_Bc2d_For_Sub, Gmi_Bc3du_For_Sub
      use GmiUtilsMetFields_mod        , only : Check_Vert_Courant
!
      implicit none
!
!
#     include "gem_msg_numbers.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: do_mean
      logical, intent(in) :: pr_flux, pr_const_flux, pr_psf_flux
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer, intent(in) :: gmi_nborder
      integer, intent(in) :: numDomains, numLonDomains, numLatDomains
      integer, intent(in) :: numSpecies
      integer, intent(in) :: ilong, ivert
      integer, intent(in) :: j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: commu_npole, commu_spole
      integer, intent(in) :: communicatorWorld
      integer, intent(in) :: nb, sb, eb, wb
      integer, intent(in) :: ewflag(ilo:ihi)
      integer, intent(in) :: nspoleu(julo:jhi)
      integer, intent(in) :: mapi_all(2,numDomains)
      integer, intent(in) :: north_domain, south_domain, east_domain, west_domain
      integer, intent(inout) :: count_flux
#     include "setkin_group.h"
#     include "setkin_par.h"
!
      integer, intent(in) :: flux_species_num
      integer, intent(in) :: flux_species (numSpecies)
      real*8 , intent(in) :: mw (numSpecies)
!
      logical, intent(in) :: do_chem_grp
      integer, intent(in) :: num_adv_time_steps
      integer, intent(in) :: advec_consrv_opt
      integer, intent(in) :: advec_flag(numSpecies)
      real*8 , intent(in) :: pt
      real*8 , intent(in) :: geofac_pc
      real*8 , intent(in) :: geofac(ju1_gl:j2_gl)
      real*8 , intent(in) :: cose(ju1_gl:j2_gl+1)
      real*8 , intent(in) :: mcor(i1:i2, ju1:j2)
      real*8 , intent(in) :: ai(k1-1:k2)
      real*8 , intent(in) :: bi(k1-1:k2)
      real*8 , intent(in) :: rel_area(i1:i2, ju1:j2)
      real*8 , intent(in) :: pctm1(ilo:ihi, julo:jhi)
      real*8 , intent(in) :: pctm2(ilo:ihi, julo:jhi)
!
      real*8 , intent(inout) :: dap(k1:k2)
      real*8 , intent(inout) :: dbk(k1:k2)
      real*8 , intent(in) :: xmass(ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in) :: ymass(ilo:ihi, julo:jhi, k1:k2)
      real*8,  intent(inout) :: const(i1:i2, ju1:j2, k1:k2, numSpecies)
!
      real*8 , intent(inout) :: air_mass_flux(i1:i2,ju1:j2,k1:k2,3)
      real*8 , intent(inout) :: flux_x(i1:i2, ju1:j2, k1:k2, flux_species_num)
      real*8 , intent(inout) :: flux_y(i1:i2, ju1:j2, k1:k2, flux_species_num)
      real*8 , intent(inout) :: flux_z(i1:i2, ju1:j2, k1:k2, flux_species_num)
      real*8 , intent(inout) :: psf_flux(i1:i2,ju1:j2)
!
!
!     -----------------------
!     Parameter declarations.
!     -----------------------
!
!     ----------------------------------------------------------------------
!     Note that DO_ADVEC_DIAG only works when running on a single processor.
!     ----------------------------------------------------------------------
!
      logical, parameter :: DO_ADVEC_DIAG         = .false.
!
      logical, parameter :: DO_CHECK_VERT_COURANT = .false.
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      logical :: first_this_tstp
      logical :: is_top
!
      integer :: ic, idx_save
      integer :: icloopmax
      integer :: ig, im, imsgrp  ! used by chemical groups functionality
      integer :: istep
      integer :: k2m1
!
      real*8  :: alpha
      real*8  :: ristep
      real*8  :: rnum_adv_time_steps
!
      real*8  :: pvadv1(ilo:ihi, julo:jhi)
      real*8  :: pvadv2(ilo:ihi, julo:jhi)
!
      real*8  :: qq1       (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: xmass_vadv(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: ymass_vadv(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: zmass_vadv(ilo:ihi, julo:jhi, k1:k2)
!.sds.. species fluxes
      real*8 :: fx (ilo:ihi, julo:jhi, k1:k2)
      real*8 :: fy (ilo:ihi, julo:jhi, k1:k2)
      real*8 :: fz (ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Update_Advec called by ', loc_proc
      end if
!
!
      if (DO_ADVEC_DIAG) then
        is_top = .true.
!       =====================
        call Write_Advec_Diag  &
!       =====================
     &    (is_top, pt, ai, bi, mcor, pctm1, const, &
     &     pr_diag, loc_proc, numSpecies, numDomains, &
     &     i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
      end if
!
!
      if (j1p /= ju1_gl+1) then
!
        do ic = 1, numSpecies
!
!         ========================
          call Average_Const_Poles  &
!         ========================
     &      (dap, dbk, rel_area, pctm1, const(:,:,:,ic), &
     &   numLonDomains, ju1_gl, j2_gl, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, &
     &   commu_npole, commu_spole)
!
        end do
!
      end if
!
!
!     ----------------------------------------------------------------------
!     Temporarily reverse k's in several variables for use by the DAO2 code.
!     ----------------------------------------------------------------------
!
      dap(:) = dap(k2:k1:-1)
      dbk(:) = dbk(k2:k1:-1)
!
!.sds2      xmass(:,:,:) = xmass(:,:,k2:k1:-1)
!.sds2      ymass(:,:,:) = ymass(:,:,k2:k1:-1)
!
      const(:,:,:,:) = const(:,:,k2:k1:-1,:)
!
!
      rnum_adv_time_steps = num_adv_time_steps
!
      xmass_vadv(:,:,k1:k2) = xmass(:,:,k2:k1:-1) / rnum_adv_time_steps
      ymass_vadv(:,:,k1:k2) = ymass(:,:,k2:k1:-1) / rnum_adv_time_steps
!.sds2      xmass_vadv(:,:,:) = xmass(:,:,:) / rnum_adv_time_steps
!.sds2      ymass_vadv(:,:,:) = ymass(:,:,:) / rnum_adv_time_steps
!
!
      pvadv2(:,:) = pctm1(:,:)
!
!
!     ========================
      call Alloc_Tpcore_Arrays &
!     ========================
     & (pr_diag, loc_proc, &
     &  ilo, ihi, julo, jhi, &
     &  i1, i2, ju1, j2, k1, k2)
!
!
!     ==========
      isteploop: do istep = 1, num_adv_time_steps
!     ==========
!
        ristep = istep
        alpha  = ristep / rnum_adv_time_steps
!
        pvadv1(:,:) = pvadv2(:,:)
        pvadv2(:,:) =  &
     &    ((1.0d0 - alpha) * pctm1(:,:)) + (alpha * pctm2(:,:))
!
        if (DO_CHECK_VERT_COURANT) then
!
!         --------------------------------------------------------------------
!         Note:
!           * the answers may be slightly different here, rather than outside
!             the loop, because pvadv1 and pvadv2 have been interpolated
!             between pctm1 and pctm2, and hence we use slightly different
!             pressures; the Courant test here is less stringent, but more
!             correct;
!           * num_adv_time_steps is always 1 here, since xmass and ymass have
!             already been scaled by num_adv_time_steps;
!           * the levels and sign reported here are opposite to those reported
!             outside tpcore; this is because tpcore works from the top down.
!         --------------------------------------------------------------------
!
!         =======================
          call Check_Vert_Courant  &
!         =======================
     &      (.true., 1, geofac_pc, geofac, dap, dbk, pvadv1, pvadv2,  &
     &       xmass_vadv, ymass_vadv, &
     &   pr_diag, loc_proc, numLonDomains, j1p, j2p, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole)
!
        end if
!
!
        if (do_chem_grp) then
!
          icloopmax = numSpecies + NUMGRP
!
          qqgrp(:,:,:,:) = 0.0d0
!
          do ig = 1, NUMGRP
            do im = 1, MAXGRP_ELEM
!
              imsgrp = sgrp_elem_map(im,ig)
!
              if (imsgrp > 0) then
                qqgrp(:,:,:,ig) =  &
     &            qqgrp(:,:,:,ig) +  &
     &            const(:,:,:,imsgrp) * sgrp_fac(im,ig)
              end if
!
            end do
          end do
!
        else
!
          icloopmax = numSpecies
!
        end if
!
        first_this_tstp = .true.
!
!.sds.. reset flux species counter
        idx_save = 0
!
!       =======
        icloop: do ic = 1, icloopmax
!       =======
!
          if (ic <= numSpecies) then
!
            if ((advec_flag(ic) /= 1) .or. (isFixedConcentration (ic))) then
!             ============
              cycle icloop
!             ============
            end if
!
            qq1(:,:,:) = 0.0d0
!
            qq1(i1:i2,ju1:j2,:) = const(i1:i2,ju1:j2,:,ic)
!
          else
!
            ig = ic - numSpecies
!
            qq1(:,:,:) = 0.0d0
!
            qq1(i1:i2,ju1:j2,:) = qqgrp(:,:,:,ig)
!
          end if
!
!         ======================
          call Gmi_Bc3du_For_Sub &
!         ======================
     &      (qq1, ACTM_BC_QQ1, &
     &       numLonDomains, numLatDomains, gmi_nborder, &
     &       i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld, &
     &     numDomains)
!
!         ===================
          call Do_Tpcore_Dao2  &
!         ===================
     &      (first_this_tstp, ic, advec_consrv_opt,  &
     &       geofac_pc, geofac, cose, dap, dbk, pvadv1, pvadv2,  &
     &       xmass_vadv, ymass_vadv, zmass_vadv, qq1, fx, fy, fz, &
     &       pr_flux, pr_diag, &
     &       loc_proc, numSpecies, numDomains, numLonDomains, numLatDomains, gmi_nborder, &
     &       j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
     &       i1, i2, ju1, j2, k1, k2, ilong, ivert, &
     &       ewflag, nspoleu, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain,  &
     &       communicatorWorld, commu_npole, commu_spole, mapi_all)
!
          first_this_tstp = .false.
!
          if (ic <= numSpecies) then
            const(:,:,:,ic) = qq1(i1:i2,ju1:j2,:)
          else
!           ==================================
#           include "setkin_group_specifics.h"
!           ==================================
          end if
!
!.sds.. save the species fluxes as diagnostic
          if ((ic <= numSpecies) .and. pr_flux .and. pr_const_flux) then
!... are we to save this constiuents flux?
            if(flux_species(ic).eq.1) then
              idx_save = idx_save+1
!             =======================
              call Accum_Constit_Flux &
!             =======================
     &          (idx_save, ic, fx, fy, fz, mcor, mw, numSpecies, &
     &           flux_x, flux_y, flux_z, &
     &           i1, i2, ju1, j2, ilo, ihi, julo, jhi, &
     &           k1, k2, flux_species_num, pr_diag, loc_proc)
!
            endif
          endif
!
!       =============
        end do icloop
!       =============
!
!.sds
!.. accumulate, over time, the vertical air-mass flux for output to the flux.nc file
        if (pr_flux) then
!         ====================
          call Accum_MassZFlux &
!         ====================
     &      (zmass_vadv, air_mass_flux, &
     &       ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, pr_diag, loc_proc)
!
         endif
!.sds
!
!       =====================
        call Gmi_Bc2d_For_Sub  &
!       =====================
     &    (pvadv2, ACTM_BC_PVADV2, &
     &     numLonDomains, numLatDomains, gmi_nborder, &
     &     i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, &
     &     ewflag, nspoleu, nb, sb, eb, wb,  &
     &     north_domain, south_domain, east_domain, west_domain,  &
     &     communicatorWorld, &
     &     numDomains)
!
!     ================
      end do isteploop
!     ================
!
!.sds
!... save the horizontal fluxes for output to the flux.nc file
      if (pr_flux) then
!       =====================
        call Accum_MassXYFlux &
!       =====================
     &  (count_flux, air_mass_flux, psf_flux, &
     &   xmass, ymass, pctm2, pr_psf_flux, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, pr_diag, loc_proc)
!
       endif
!.sds
!
!     ==========================
      call Dealloc_Tpcore_Arrays (pr_diag, loc_proc)
!     ==========================
!
!     -------------------------------------------------------------
!     Reverse k's in several variables back the way they should be.
!     -------------------------------------------------------------
!
      dap(:) = dap(k2:k1:-1)
      dbk(:) = dbk(k2:k1:-1)
!
!.sds2      xmass(:,:,:) = xmass(:,:,k2:k1:-1)
!.sds2      ymass(:,:,:) = ymass(:,:,k2:k1:-1)
!
      const(:,:,:,:) = const(:,:,k2:k1:-1,:)
!
      if (DO_ADVEC_DIAG) then
        is_top = .false.
!       =====================
        call Write_Advec_Diag  &
!       =====================
     &    (is_top, pt, ai, bi, mcor, pctm2, const, &
     &     pr_diag, loc_proc, numSpecies, numDomains, &
     &     i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
        Write (6,'(a,es9.2)')  &
     &    'Max diff pctm2/pvadv2 after advection = ',  &
     &    Maxval (Abs (pctm2(:,:) - pvadv2(:,:)))
      end if
!
      return
!
      end
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Average_Const_Poles
!
! DESCRIPTION
!   This routine averages the species concentrations at the Poles when the
!   Polar cap is enlarged.  It makes the last two latitudes equal.
!
! ARGUMENTS
!   dap      : pressure difference across layer from (ai * pt) term (mb)
!   dbk      : difference in bi across layer - the dSigma term
!   rel_area : relative surface area of grid box (fraction)
!   pctm1    : CTM surface pressure at t1 (mb)
!   const1   : species concentration, known at zone center (mixing ratio)
!
!-----------------------------------------------------------------------------
!
      subroutine Average_Const_Poles  &
     &  (dap, dbk, rel_area, pctm1, const1, &
     &   numLonDomains, ju1_gl, j2_gl, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, &
     &   commu_npole, commu_spole)
!
      use GmiReduce_mod, only : Gmi_Sum1_Pole_Reduce
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numLonDomains
      integer, intent(in) :: ju1_gl, j2_gl
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: commu_npole, commu_spole
!
      real*8, intent(in)  :: dap(k1:k2)
      real*8, intent(in)  :: dbk(k1:k2)
      real*8, intent(in)  :: rel_area(i1:i2, ju1:j2)
      real*8, intent(in)  :: pctm1   (ilo:ihi, julo:jhi)
!
      real*8, intent(inout) :: const1(i1:i2, ju1:j2, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: ik
!
      real*8  :: meanq
      real*8  :: sum1, sum2
!
!     -----------------------------------------------------------------
!     delp1n : pressure thickness at North Pole, the psudo-density in a
!              hydrostatic system at t1 (mb)
!     delp1s : pressure thickness at South Pole, the psudo-density in a
!              hydrostatic system at t1 (mb)
!     -----------------------------------------------------------------
!
      real*8  :: delp1n(i1:i2, j2-1:j2,    k1:k2)
      real*8  :: delp1s(i1:i2,  ju1:ju1+1, k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
!     =================
      if (ju1 == ju1_gl) then
!     =================
!
        do ik = k1, k2
!
          delp1s(i1:i2,ju1:ju1+1,ik) =  &
     &      dap(ik) +  &
     &      (dbk(ik) * pctm1(i1:i2,ju1:ju1+1))
!
          sum1 = Sum (const1  (:,ju1:ju1+1,ik) *  &
     &                delp1s  (:,ju1:ju1+1,ik) *  &
     &                rel_area(:,ju1:ju1+1))
!
          sum2 = Sum (delp1s  (:,ju1:ju1+1,ik) *  &
     &                rel_area(:,ju1:ju1+1))
!
!         ============================
          call Gmi_Sum1_Pole_Reduce (sum1, numLonDomains, ju1, j2, ju1_gl, j2_gl, &
                                        commu_npole, commu_spole)
          call Gmi_Sum1_Pole_Reduce (sum2, numLonDomains, ju1, j2, ju1_gl, j2_gl, &
                                        commu_npole, commu_spole)
!         ============================
!
          meanq = sum1 / sum2
!
          const1(:,ju1:ju1+1,ik) = meanq
!
        end do
!
      end if
!
!
!     ================
      if (j2 == j2_gl) then
!     ================
!
        do ik = k1, k2
!
          delp1n(i1:i2,j2-1:j2,ik) =  &
     &      dap(ik) +  &
     &      (dbk(ik) * pctm1(i1:i2,j2-1:j2))
!
          sum1 = Sum (const1  (:,j2-1:j2,ik) *  &
     &                delp1n  (:,j2-1:j2,ik) *  &
     &                rel_area(:,j2-1:j2))
!
          sum2 = Sum (delp1n  (:,j2-1:j2,ik) *  &
     &                rel_area(:,j2-1:j2))
!
!         ============================
          call Gmi_Sum1_Pole_Reduce (sum1, numLonDomains, ju1, j2, ju1_gl, j2_gl, &
                                        commu_npole, commu_spole)
          call Gmi_Sum1_Pole_Reduce (sum2, numLonDomains, ju1, j2, ju1_gl, j2_gl, &
                                        commu_npole, commu_spole)
!         ============================
!
          meanq = sum1 / sum2
!
          const1(:,j2-1:j2,ik) = meanq
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
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Write_Advec_Diag
!
! DESCRIPTION
!   This routine writes out some advection diagnostics.
!
!   NOTE THAT CONST MUST BE INDEXED FROM THE GROUND UP HERE.
!
! ARGUMENTS
!   is_top : is call to this routine being made from the top or the bottom
!            of Update_Advec?
!   pt     : pressure = (ai * pt) + (bi * psf) (mb)
!   ai     : pressure = (ai * pt) + (bi * psf), ai at zone interface
!   bi     : pressure = (ai * pt) + (bi * psf), bi at zone interface
!   mcor   : area of grid box     (m^2)
!   pctm   : CTM surface pressure (mb)
!   const  : species concentration, known at zone centers (mixing ratio)
!
!-----------------------------------------------------------------------------
!
      subroutine Write_Advec_Diag  &
     &  (is_top, pt, ai, bi, mcor, pctm, const, &
     &   pr_diag, loc_proc, numSpecies, numDomains, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
      use GmiPrintError_mod, only : GmiPrintError
      use GmiTotalMass_mod , only : calcTotalMass
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
      integer, intent(in) :: loc_proc, numDomains
      integer, intent(in) :: numSpecies
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      logical,intent(in)  :: is_top
      real*8 ,intent(in)  :: pt
      real*8 ,intent(in)  :: ai(k1-1:k2)
      real*8 ,intent(in)  :: bi(k1-1:k2)
      real*8 ,intent(in)  :: mcor (i1:i2, ju1:j2)
      real*8 ,intent(in)  :: pctm (ilo:ihi, julo:jhi)
      real*8 ,intent(in)  :: const(i1:i2, ju1:j2, k1:k2, numSpecies)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      character (len=75) :: err_msg
!
      logical, save :: first1 = .true.
      logical, save :: first2 = .true.
!
      integer :: ic
!
      real*8  :: fac
      real*8  :: timestep_air_drift
      real*8  :: tot_air_drift, tot_air_moles
!
      real*8, save :: cum_adv_air_drift
      real*8, save :: initial_air_moles, start_air_moles
!
      real*8  :: rel_tracer_drift     (numSpecies)
      real*8  :: timestep_tracer_drift(numSpecies)
      real*8  :: total_tracer_drift   (numSpecies)
      real*8  :: tracer_moles         (numSpecies)
!
      real*8  :: mass (i1:i2, ju1:j2, k1:k2)
      real*8  :: moles(i1:i2, ju1:j2, k1:k2)
!
      real*8, allocatable, save :: cum_adv_tracer_drift(:)
      real*8, allocatable, save :: initial_tracer_moles(:)
      real*8, allocatable, save :: start_tracer_moles  (:)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Write_Advec_Diag called by ', loc_proc
      end if
!
!
      if (numDomains /= 1) then
        err_msg = 'Write_Advec_Diag not parallelized.'
        call GmiPrintError  &
     &    (err_msg, .true., 1, numDomains, 0, 0, 0.0d0, 0.0d0)
      end if
!
!
      if (first1) then
!
        first1 = .false.
!
        Allocate (cum_adv_tracer_drift(1:numSpecies))
        Allocate (initial_tracer_moles(1:numSpecies))
        Allocate (start_tracer_moles  (1:numSpecies))
        cum_adv_tracer_drift = 0.0d0
        initial_tracer_moles = 0.0d0
        start_tracer_moles   = 0.0d0
!
        cum_adv_air_drift = 0.0d0
!
      end if
!
!
      tracer_moles(:) = 0.0d0
!
!
!     ==================
      call calcTotalMass  &
!     ==================
     &  (pt, ai, bi, pctm, mcor, mass, &
     &   pr_diag, loc_proc, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
!
      fac = 1.0d0 / (MWTAIR * KGPG)
!
      moles(:,:,:) = mass(:,:,:) * fac
!
      tot_air_moles = Sum (moles(:,:,:))
!
!
      if (is_top) then
!
        Write (6,*) '-------- Top of Update_Advec  --------'
!
        Write (6,900) 'Moles of air     = ', tot_air_moles
!
        initial_air_moles = tot_air_moles
!
      else
!
        timestep_air_drift = tot_air_moles - initial_air_moles
!
        tot_air_drift      = tot_air_moles - start_air_moles
!
        cum_adv_air_drift  = cum_adv_air_drift + timestep_air_drift
!
        Write (6,*) '------- Bottom of Update_Advec -------'
!
        Write (6,900) 'Moles of Air     = ', tot_air_moles
!
        Write (6,905)  &
     &    '(Drift: Timestep = ', timestep_air_drift, ',  ',  &
     &    ' Tot. advect.    = ', cum_adv_air_drift,  ',  ',  &
     &    ' Total           = ', tot_air_drift, ')'
!
 900    format (a, es21.14)
 905    format (a, es9.2, a, a, es9.2, a, a, es9.2, a)
!
!
      end if
!
!
      do ic = 1, numSpecies
!
        tracer_moles(ic) =  &
     &    Sum (moles(:,:,:) * const(:,:,:,ic))
!
        if (is_top) then
!
          Write (6,910)  &
     &      'Moles of tracer (',ic,') = ', tracer_moles(ic)
!
          initial_tracer_moles(ic) = tracer_moles(ic)
!
        else
!
          total_tracer_drift(ic)    =  &
     &      tracer_moles(ic) - start_tracer_moles(ic)
!
          timestep_tracer_drift(ic) =  &
     &      tracer_moles(ic) - initial_tracer_moles(ic)
!
          cum_adv_tracer_drift(ic)  =  &
     &        cum_adv_tracer_drift(ic) + timestep_tracer_drift(ic)
!
          rel_tracer_drift(ic)      =  &
     &      timestep_tracer_drift(ic) /  &
     &      (initial_tracer_moles(ic) *  &
     &       Epsilon (initial_tracer_moles(ic)))
!
          if (Log10 (Abs (rel_tracer_drift(ic)) + 1d-9) < 3.0d0) then
!
            Write (6,920)  &
     &        'Moles of Tracer (', ic, ') = ', tracer_moles(ic),  &
     &        '  (Rel. drift = ', rel_tracer_drift(ic), ' * m.p.)'
!
          else
!
            Write (6,930)  &
     &        'Moles of Tracer (', ic, ') = ', tracer_moles(ic),  &
     &        '  (Rel. drift = ', rel_tracer_drift(ic), ' * m.p.)'
!
          end if
!
          Write (6,905)  &
     &      '(Drift: Timestep = ', timestep_tracer_drift(ic), ',  ',  &
     &      ' Tot. advect.    = ', cum_adv_tracer_drift (ic), ',  ',  &
     &      ' Total           = ', total_tracer_drift(ic), ')'
!
        end if
!
      end do
!
!
 910  format (a, i2, a, es21.14)
 920  format (a, i2, a, es21.14, a,  f7.1, a)
 930  format (a, i2, a, es21.14, a, es7.1, a)
!
!
      Write (6,*) '--------------------------------------'
!
!
      if (first2) then
!
        first2 = .false.
!
        start_air_moles       = tot_air_moles
        start_tracer_moles(:) = tracer_moles(:)
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
!   Reset_Flux
!
! DESCRIPTION
!   This routine resets the flux diagnostics.
!
! ARGUMENTS
!
!
!-----------------------------------------------------------------------------
!
      subroutine Reset_Flux &
     &  (pr_flux, pr_const_flux, pr_psf_flux, do_flux_reset, &
     &   air_mass_flux, count_flux, flux_x, flux_y, flux_z, psf_flux, &
     &   i1, i2, ju1, j2, k1, k2, flux_species_num, pr_diag, loc_proc)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in ) :: i1, i2, ju1, j2, k1, k2, loc_proc
      integer, intent(in ) :: flux_species_num
      logical, intent(in ) :: pr_flux, pr_const_flux, pr_psf_flux, pr_diag
      logical, intent(out) :: do_flux_reset
      integer, intent(out) :: count_flux
      real*8 , intent(out) :: air_mass_flux(i1:i2,ju1:j2,k1:k2,3)
      real*8 , intent(out) :: flux_x(i1:i2, ju1:j2, k1:k2, flux_species_num)
      real*8 , intent(out) :: flux_y(i1:i2, ju1:j2, k1:k2, flux_species_num)
      real*8 , intent(out) :: flux_z(i1:i2, ju1:j2, k1:k2, flux_species_num)
      real*8 , intent(out) :: psf_flux(i1:i2,ju1:j2)
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      logical, save :: first = .true.
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Reset_Flux called by ', loc_proc
      end if
!
      if(first) then
        first = .false.
        do_flux_reset = .true.
       endif
!... reset accumulated quantities
!      if(pr_flux .and. do_flux_reset) then
      if(pr_flux) then
!        print *,'reset flux2:',loc_proc
!... zero air mass flux
        air_mass_flux = 0.0d0
!... reset counter
        count_flux = 0
!
!... constituent component fluxes
        if(pr_const_flux) then
          flux_x = 0.0d0
          flux_y = 0.0d0
          flux_z = 0.0d0
        endif
!
!... constituent component fluxes
        if(pr_psf_flux) then
          psf_flux = 0.0d0
        endif
!... reset switch
        do_flux_reset = .false.
!
       endif
!
      return
!
      end
!
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Accum_Constit_Flux
!
! DESCRIPTION
!   This routine stores the species 3-D fluxes.
!
! ARGUMENTS
!   ispc : species flux index
!   ic   : species index in model
!   fx   : E-W mass flux for species ispc
!   fy   : N-S mass flux for species ispc
!   fz   : vertical mass flux for species ispc
!
!-----------------------------------------------------------------------------
!
      subroutine Accum_Constit_Flux &
     & (ispc, ic, fx, fy, fz, mcor, mw, numSpecies, &
     &  flux_x, flux_y, flux_z, &
     &  i1, i2, ju1, j2, ilo, ihi, julo, jhi, &
     &  k1, k2, flux_species_num, pr_diag, loc_proc)
!
      implicit none
!
!... need mwtair
#     include "gmi_phys_constants.h"
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: ispc, ic
!
      integer, intent(in) :: i1, i2, ju1, j2, ilo, ihi, julo, jhi
      integer, intent(in) :: k1, k2, flux_species_num, loc_proc
      logical, intent(in) :: pr_diag
      integer,intent(in) :: numSpecies
!... species fluxes
      real*8, intent(in) :: fx (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in) :: fy (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in) :: fz (ilo:ihi, julo:jhi, k1:k2)
!
      real*8, intent(in) :: mcor(i1:i2,ju1:j2)
      real*8, intent(in) :: mw(numSpecies)
!
      real*8, intent(inout) :: flux_x (i1:i2, ju1:j2, k1:k2, flux_species_num)
      real*8, intent(inout) :: flux_y (i1:i2, ju1:j2, k1:k2, flux_species_num)
      real*8, intent(inout) :: flux_z (i1:i2, ju1:j2, k1:k2, flux_species_num)
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: ik
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Accum_Constit_Flux called by ', loc_proc
      end if
!
!       ---------------------------------------------------
!       Extract mass fluxes for output to flux NetCDF file.
!       ---------------------------------------------------
!
!... Accumulate the mass fluxes (and weight wrt molecular weight to get kg)
      do ik = k1,k2
        flux_x(i1:i2,ju1:j2,ik,ispc) = flux_x(i1:i2,ju1:j2,ik,ispc)  &
     &   + (fx(i1:i2,ju1:j2,ik) * 100.d0/GMI_G * mw(ic)/MWTAIR) &
     &   * mcor(i1:i2,ju1:j2)
!
        flux_y(i1:i2,ju1:j2,ik,ispc) = flux_y(i1:i2,ju1:j2,ik,ispc)  &
     &   + (fy(i1:i2,ju1:j2,ik) * 100.d0/GMI_G * mw(ic)/MWTAIR) &
     &   * mcor(i1:i2,ju1:j2)
!
        flux_z(i1:i2,ju1:j2,ik,ispc) = flux_z(i1:i2,ju1:j2,ik,ispc)  &
     &   + (fz(i1:i2,ju1:j2,ik) * 100.d0/GMI_G * mw(ic)/MWTAIR) &
     &   * mcor(i1:i2,ju1:j2)
      enddo
!
!       -----------------------------------------------------
!       Special processing when using extended Polar boxes =>
!       duplicate the mass fluxes.
!       -----------------------------------------------------
!      if (j1p .ne. ju1_gl+1) then
!        flux_z(i1:i2,ju1_gl+1,k1:k2,ispc) = flux_z(i1:i2,ju1_gl,k1:k2,ispc)
!
!        flux_z(i1:i2,j2_gl-1,k1:k2,ispc) = flux_z(i1:i2,j2_gl,k1:k2,ispc)
!      endif
!
!
      return
!
      end
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Store_MZFlux
!
! DESCRIPTION
!   This routine stores the species 3-D fluxes.
!
! ARGUMENTS
!   zmass_vadv : vertical air mass flux for this sub step
!
!-----------------------------------------------------------------------------
!
      subroutine Accum_MassZFlux &
     &  (zmass_vadv, air_mass_flux, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, pr_diag, loc_proc)
!
      implicit none
!
!
!     ----------------------
!     Argument declarations.
!     ----------------------
      integer, intent(in) :: ilo, ihi, julo, jhi, loc_proc
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      logical, intent(in) :: pr_diag
!... species fluxes
      real*8, intent(in) :: zmass_vadv (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(inout) :: air_mass_flux (i1:i2,ju1:j2,k1:k2, 3)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
      integer :: k2m1
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Accum_MassZFlux called by ', loc_proc
      end if
!
!       ---------------------------------------------------
!       Store vertical air mass flux for output to flux NetCDF file.
!       ---------------------------------------------------
!
      k2m1 = k2 - 1
!... vertical mass flux, flip in vertical (zero at top)
      air_mass_flux(i1:i2,ju1:j2,k2,3) = 0.0d0
!
      air_mass_flux(i1:i2,ju1:j2,k1:k2m1,3) = air_mass_flux(i1:i2,ju1:j2,k1:k2m1,3) &
     &         + zmass_vadv(i1:i2,ju1:j2,k2m1:k1:-1)
!
      return
      end
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Store_MXYFlux
!
! DESCRIPTION
!   This routine stores the horizontal air mass fluxes.
!
! ARGUMENTS
!   xmass : E-W air mass flux for this sub step
!   ymass : N-S air mass flux for this sub step
!
!-----------------------------------------------------------------------------
!
      subroutine Accum_MassXYFlux &
     &  (count_flux, air_mass_flux, psf_flux, &
     &   xmass, ymass, pctm2, pr_psf_flux, &
     &   ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, pr_diag, loc_proc)
!
      implicit none
!
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2
      integer, intent(in) ::  loc_proc
      integer, intent(inout) :: count_flux
      logical, intent(in) :: pr_diag, pr_psf_flux
!
!... species fluxes
      real*8 , intent(in) :: xmass(ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in) :: ymass(ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in) :: pctm2(ilo:ihi, julo:jhi)
      real*8 , intent(inout) :: air_mass_flux(i1:i2, ju1:j2, k1:k2, 3)
      real*8 , intent(inout) :: psf_flux(i1:i2, ju1:j2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Accum_MassXYFlux called by ', loc_proc
      end if
!       --------------------------------------------------------------
!       Store horizontal air mass flux for output to flux NetCDF file.
!       --------------------------------------------------------------
      count_flux = count_flux + 1
!... E-W mass flux
      air_mass_flux(i1:i2,ju1:j2,k1:k2,1) = air_mass_flux(i1:i2,ju1:j2,k1:k2,1) &
     &       + xmass(i1:i2,ju1:j2,k2:k1:-1)
!... N-S mass flux
      air_mass_flux(i1:i2,ju1:j2,k1:k2,2) = air_mass_flux(i1:i2,ju1:j2,k1:k2,2) &
     &       + ymass(i1:i2,ju1:j2,k2:k1:-1)
!
      if(pr_psf_flux) then
        psf_flux(i1:i2,ju1:j2) = psf_flux(i1:i2,ju1:j2) &
     &       + pctm2(i1:i2,ju1:j2)
      endif
!
      return
      end
!
!
