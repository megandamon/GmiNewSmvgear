!=============================================================================
!
! $Id: sulfchem.F90,v 1.2 2011-08-09 22:12:58 mrdamon Exp $
!
! CODE DEVELOPER
!   originated from GRANTOUR programmed by John Walton
!   modified by Xiaohong Liu
!
! FILE
!   sulfchem.F
!
! ROUTINES
!   Do_Sulfchem
!
!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Sulfchem
!
! DESCRIPTION
!
!     *****************************************************************
!     * carries out the sulfur chemistry calculation
!     *****************************************************************
!
!     Performs a numerical integration of the full aqueous chemistry
!      equations for the species H2SO4, DMS, SO2 & H2O2 in cloudy air.
!     Performs an analytic calculation of the clear sky part of the
!      chemistry equations for the species H2SO4, DMS, SO2 & H2O2 in
!      clear air.
!     Combines the clear and cloudy concentrations to simulate mixing
!      between clear and cloudy regions on the basis of cloud fraction
!
!     ********** input:
!     *     dtchem  =  model time step in seconds
!     *     dtclc   =  time step for clear air chemistry
!     *     dtclmx  =  max. time step between mixing calculations in seconds
!     *     dtcl    =  time step between mixing calculations in seconds
!     *     dtaq    =  time step used in integration of the aqueous
!                      equations in seconds
!     *     dph     =  assumed pH of the water drops
!
!     *     tgcm       air temperature (K)
!     *     zmair      air density (# cm-3)
!     *     gmair      air density (gm cm-3)
!     *     oh         OH  concentration (cm-3)
!     *     ho2        HO2 concentration (cm-3)
!     *     no3        NO3 concentration (cm-3)
!     *     o3         O3  concentration (cm-3)
!     *     cfac       total cloud fraction (0-1)
!     *     cwac       cloud water mixing ratio for aq chem (g/g) (in-cloud)
!
!     *     qh1p       DMS + OH  reaction rate coefficient
!     *     qh2p       SO2 + OH     "      "        "
!     *     qh3p       HO2 + HO2    "      "        "
!     *     qj4p       H2O2 + hv photolysis rate coef.
!     *     qh5p       H2O2 + OH reaction rate coefficient
!     *     qh6p       SO2(aq) + H2O2(aq)  "        "
!     *     qh7p       SO2(aq) + O3(aq)    "        "
!     *     qh8p       HO2 + HO2 + H2O     "        "
!     *     qh9p       DMS + NO3 reaction rate coefficient
!
!     *     a1->a9     reaction coefficients
!
!
!     ********** output:
!
!     * Here is how I treat separate changes due to the aqueous
!       and clear sky chemical paths.
!
!       I assume that aqueous and clear sky chemistry act on the
!       total concentration, C, giving the concentration changes
!       dCa(C) & dCn(C).
!
!       The new values of the concentration at t+dt are
!
!           C(t+dt) = C(t) + cfp * dCa(C) + ( 1 - cfp ) * dCn(C)
!
!       In the case of H2SO4, these two terms contribute separately
!       to aqeuous and clear sky concentration.
!
!       In the case of all clear air or all cloud, these reduce to
!       what one would expect.
!
!       Definitions
!
!       C           species concentration
!
!       dCa, dCn    Concentration changes for aqueous & clear sky
!                   chemistry.
!
!       cfp         fractional cloud cover.
!
!     *----------------------------------------------------------------
!
! ARGUMENTS
!   itloop    : # of zones (ilong * ilat * ivert)
!   time      : time interval for each chemical calculation (seconds)
!   tgcm      : temperature (K)
!   pgcm      : atmospheric pressure at the center of each grid box (mb)
!   xr        : species concentration, known at zone centers
!               (molecules/cm^3 for gas-species, kg/kg for aerosol mass,
!                #particles/kg for sulfate aerosol number)
!   humidity  : specific humidity (g/kg)
!   semiss    : array of sulfur emissions (molecules/cm^3/s)
!   qjimp     : photolysis rate constants (s^-1)
!   qkimp     : thermal    rate constants (units vary)
!   cfac      : total cloud fraction [0 - 1]
!   cwac      : in-cloud liquid water in each grid box (g/g)
!   relh      : relative humidity [0 - 1]
!   aqua_infile_name : aquachem input file name
!
!-----------------------------------------------------------------------------

      subroutine Do_Sulfchem  &
     &  (itloop, time, tgcm, pgcm, xr, humidity,  &
     &   semiss, qjimp, qkimp, cfac, cwac, relh,  &
     &   aqua_infile_name, pr_diag, loc_proc, num_time_steps,  &
     &   massc, pr_sulf_src, dms_oh, dms_no3, so2_oh, so2_h2o2, so2_o3)

!
!---->load modules
!
      use modinit_umaer

      implicit none

#     include "gmi_phys_constants.h"
#     include "setkin_par.h"
#     include "sulfchem.h"
#     include "umaerosol.h"

!     ----------------------
!     IMPACT-related values.
!     ----------------------

      character*(*) :: aqua_infile_name

      integer :: itloop
      real*8  :: time
      real*8  :: tgcm    (itloop)
      real*8  :: pgcm    (itloop)
      real*8  :: xr      (itloop, NSP)
      real*8  :: humidity(itloop)
      real*8  :: semiss  (itloop, NSP)
      real*8  :: qjimp   (itloop, NUM_J)
      real*8  :: qkimp   (itloop, NUM_K)
      real*8  :: cfac    (itloop)
      real*8  :: cwac    (itloop)
      real*8  :: relh    (itloop)

      logical :: pr_diag
      integer :: loc_proc
      integer :: num_time_steps

      real*8  :: massc   (itloop)
      logical :: pr_sulf_src
      real*8  :: dms_oh  (itloop)
      real*8  :: dms_no3 (itloop)
      real*8  :: so2_oh  (itloop)
      real*8  :: so2_h2o2(itloop)
      real*8  :: so2_o3  (itloop)


! ----------------------------------------------------------
!  *** LOCAL  VARIABLES WERE DECLARED HERE.
! ----------------------------------------------------------

      integer :: numcl

      real*8   dtclmx, dtaqset,  &
     &         accset, dph

      real*8, parameter :: CMPM3 = CMPM * CMPM * CMPM

      parameter (dtclmx = 14400.d0)
      parameter (dtaqset = 600.d0)
      parameter (accset = 0.0001d0)
      parameter (dph = 4.5d0)

      real*8  &
     &  dtaq,    & ! time step used in integration of aqueous equations in seconds
     &  dtchem,  & ! model time step in seconds
     &  dtcl,    & ! time step between mixing calculations in seconds
     &  dtclc  ! time step for clear air chemistry

!     qh6       SO2(aq) + H2O2(aq) reaction rate coefficient
      real*8  qh6

!---->stuff for umaerosol
      logical, save :: first = .true.

      integer :: ng_i, nso4_i, nmomso4_i, nnon_i, nmomnon_i

      real*8 :: press (itloop)
      real*8 :: so4gas(itloop)
      real*8 :: srcgas(itloop)
      real*8 :: srcmas(itloop)
      real*8 :: srcpar(itloop)
      real*8 :: so4aer(nmomso4,nso4,itloop)
      real*8 :: so4non(nnon,itloop)
      real*8 :: aernon(nmomnon,nnon,itloop)
!---->end of stuff for umaerosol

!     air density (molecules/cm3).
      real*8 :: zmair(itloop)

!     air density (g/cm3).
      real*8 :: gmair(itloop)

!     air pressure (bar)
      real*8 :: pgcmb(itloop)

!     H2O concentration in molecules/cm3.
      real*8 :: h2o(itloop)

!     concentrations of O3, OH, HO2 and NO3 in molecules/cm3.
      real*8 :: o3(itloop)
      real*8 :: oh(itloop)
      real*8 :: ho2(itloop)
      real*8 :: no3(itloop)

!     rtcd (s-1) rate coef. for SO2 -> dust transfer
!     rtcs (s-1) rate coef. for SO2 -> sslt transfer
      real*8 :: rtcd(itloop,4)
      real*8 :: rtcs(itloop,4)

!     cp and dcp have units [molecules/cm3 for gas-species and
!     gm spec./gm air for aerosols]
      real*8 :: cp (itloop,NSP)
      real*8 :: dcp(itloop,NSP)

!     H2O2 photolysi rate coefficient in 1/s
      real*8 :: qj4p(itloop)

!     qh1p       DMS + OH  reaction rate coefficient
!     qh2p       SO2 + OH     "      "        "
!     qh3p       HO2 + HO2    "      "        "
!     qh5p       H2O2 + OH reaction rate coefficient
!     qh7p       SO2(aq) + O3(aq)    "        "
!     qh8p       HO2 + HO2 + H2O     "        "
!     qh9p       DMS + NO3    "      "        "

      real*8 :: qh1p(itloop)
      real*8 :: qh2p(itloop)
      real*8 :: qh3p(itloop)
      real*8 :: qh5p(itloop)
      real*8 :: qh7p(itloop)
      real*8 :: qh8p(itloop)
      real*8 :: qh9p(itloop)

!     henry's law coefficients (hwp, hyp, h3p) for SO2, H2O2, O3
      real*8 :: hwp(itloop)
      real*8 :: hyp(itloop)
      real*8 :: h3p(itloop)

!     array of grid numbers for grids in clouds
      integer :: idxcl(itloop)

!     a1-9 coefficients in aqueous equations (only for cloudy grids)
      real*8 :: acoef(itloop,9)

      zmair = 0.0d0
      gmair = 0.0d0
      h2o   = 0.0d0
      o3    = 0.0d0
      oh    = 0.0d0
      ho2   = 0.0d0
      no3   = 0.0d0
      rtcd  = 0.0d0
      rtcs  = 0.0d0
      cp    = 0.0d0
      dcp   = 0.0d0

      qj4p  = 0.0d0
      qh1p  = 0.0d0
      qh2p  = 0.0d0
      qh3p  = 0.0d0
      qh5p  = 0.0d0
      qh7p  = 0.0d0
      qh8p  = 0.0d0
      qh9p  = 0.0d0
      hwp   = 0.0d0
      hyp   = 0.0d0
      h3p   = 0.0d0
      idxcl = 0.0d0
      acoef = 0.0d0

      so4gas = 0.0d0
      srcgas = 0.0d0
      srcmas = 0.0d0
      srcpar = 0.0d0
      so4aer = 0.0d0
      so4non = 0.0d0
      aernon = 0.0d0

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Do_Sulfchem called by ', loc_proc
      end if

      if (first) then

        first = .false.

        nso4_i   =nso4
        nmomso4_i=nmomso4
        nnon_i   =nnon
        nmomnon_i=nmomnon
        ng_i     =ng

        call umaerinit(nso4_i,nmomso4_i,nnon_i,nmomnon_i,ng_i,  &
     &                 sigmod,radmerg,  &
     &                 radgsrc,siggsrc,  &
     &                 radgso4,siggso4,fracso4,  &
     &                 radgnon,siggnon,fracnon,  &
     &                 rhonon,xmolnon,nrh_flg,  &
     &                 crh1,crh2,crh3,crh4)

      end if


!     *****************************************************************
!     * Set chemistry parameters
!     *****************************************************************

      dtchem = time
      dtcl   = dtclmx
      dtaq   = dtaqset

      hpp  = 10**(-dph)
      rhpp = 1./hpp

!     *****************************************************************
!     * Compute air density, "zmair" [#/cm3] and "gmair" [gm/cm3] at
!       model grid points.
!     *****************************************************************

      do ijk = 1, itloop
        zmair(ijk) = xr(ijk, IMGAS)
        gmair(ijk) = ( MWTAIR / AVOGAD ) * zmair(ijk)
      end do

!     *****************************************************************
!     * Calculate water vapor concentration in molecules/cm^3
!     *****************************************************************

      do ijk = 1, itloop
        h2o(ijk) = humidity(ijk) * MWTAIR / MWTH2O / GPKG  &
     &             * xr(ijk, IMGAS)
      end do

      do ijk = 1, itloop
        o3(ijk)  = xr(ijk, IO3)
        oh(ijk)  = xr(ijk, IOH)
        ho2(ijk) = xr(ijk, IHO2)
        no3(ijk) = xr(ijk, INO3)
      end do

!     change unit of pressure at the grid points from millibar to bar
      do ijk = 1, itloop
        pgcmb(ijk) = pgcm(ijk) / MBPB      ! bar
      end do

      do ijk = 1, itloop
        press(ijk) = pgcm(ijk) * PASPMB    ! unit: Pa
      end do

!     *****************************************************************
!     * Define the H2O2 photolysis rate coefficient
!     *****************************************************************

      do ijk = 1, itloop
        qj4p(ijk) = qjimp(ijk, 1)
      end do

!     *****************************************************************
!     * Given by IMPACT concentrations for species changed by
!       sulfur chemistry
!     *****************************************************************

      do ic = 1, NSP
        cp(:, ic) = xr(:, ic)
      end do


!     *****************************************************************
!     * Calculate the SO2 mass transfer rate coefficients for
!       dust and sea salt
!     *****************************************************************

#ifndef NONON
!          =================
!     call Do_So2toPar_Rates
!          =================
!    &  (itloop, pgcmb, tgcm,
!    &   relh, aqua_infile_name,
!    &   gmair, cp, rtcd, rtcs)
#endif

!     *****************************************************************
!     * Calculate gas-phase reaction rates
!       (qh1p, qh2p, qh3p, qh5p, qh8p, qh9p)
!     *****************************************************************

!          ===============
      call Do_QK_Gas_Rates  &
!          ===============
     &  (itloop, tgcm, zmair,  &
     &   qh1p, qh2p, qh3p, qh5p, qh8p, qh9p)


!     ================================================================
      if (Mod (num_time_steps, (Nint(dtcl) / Nint(dtchem))) == 0) then
!     ================================================================

!     *****************************************************************
!     * Set up for aqueous chemistry
!       (generate numcl, idxcl)
!     *****************************************************************

!          ================
      call Do_Setup_Aquchem  &
!          ================
     &  (itloop, numcl, idxcl, cfac, cwac)


      if (numcl .gt. 0) then

!     *******************************************************************
!     * There are some grid points in clouds. Set up for aqueous chem.
!       calculate the aqueous-phase reaction rates, effective Henry's law
!       coefficients, and aqueous equation coeff. (a1-a9)
!     *******************************************************************

!            ===============
        call Do_QK_Aqu_Rates  &
!            ===============
     &    (itloop, tgcm, hpp, rhpp, numcl, idxcl,  &
     &     oh, ho2, h2o, no3, o3,  &
     &     cwac, zmair,  &
     &     qh1p, qh2p, qh3p, qj4p, qh5p, qh8p, qh9p,  &
     &     qh6, qh7p, hwp, hyp, h3p,  &
     &     acoef)

        dtclc = dtcl

      else

!     ******************************************************************
!     * There are no cloudy grid points, advance the entire dtchem step.
!     * Advance the clear air chemistry only.
!     ******************************************************************

        dtclc = dtcl   ! dtchem

      endif

!     =====
      endif
!     =====

!.... time step over cloud mixing time (dtclc)

!     tcl = 0.0d0

!     ========
! 20  continue
!     ========


!     ***********************************************************
!     * Calculate the species changes due to clear sky chemistry.
!     ***********************************************************

!          ===============
      call Do_Gas_Sulfchem  &
!          ===============
!    &  (dtclc, itloop, cp, dcp, cfac, zmair, semiss,
     &  (dtchem, itloop, cp, dcp, cfac, zmair, semiss, srcgas,  &
     &   oh, ho2, h2o, no3,  &
     &   massc, pr_sulf_src, dms_oh, dms_no3, so2_oh,  &
     &   qh1p, qh2p, qh3p, qj4p, qh5p, qh8p, qh9p)


!     =================================================================
      if (Mod (num_time_steps, (Nint(dtcl) / Nint(dtchem))) == 0) then
!     =================================================================

!     *****************************************************************
!     * Integrate aqueous equations over dtclc and add the aqueous
!       change to the clear sky changes already computed.
!     *****************************************************************

      if (numcl .gt. 0) then

!            ===============
        call Do_Aqu_Sulfchem  &
!            ===============
     &    (accset, dtaq, dtchem, dtcl, dtclc,  &
     &     itloop, numcl, idxcl, semiss, srcmas,  &
     &     massc, pr_sulf_src,  &
     &     dms_oh, dms_no3, so2_oh, so2_h2o2, so2_o3,  &
     &     cp, dcp, cfac, zmair, acoef, loc_proc)

      end if

!     ======
      end if
!     ======

!
!---->call umaerosol
!
      where(cp(:, ISO4G)    < 1.0d-30) cp(:, ISO4G)  = 1.0d-30

      where(cp(:, ISO4M1)   < 1.0d-30) cp(:, ISO4M1) = 1.0d-30
      where(cp(:, ISO4N1)   < 1.0d-30) cp(:, ISO4N1) = 1.0d-30
      where(cp(:, ISO4M2)   < 1.0d-30) cp(:, ISO4M2) = 1.0d-30
      where(cp(:, ISO4N2)   < 1.0d-30) cp(:, ISO4N2) = 1.0d-30

      where(cp(:, ISO4NOC:ISO4S4) < 1.0d-30)  &
     &      cp(:, ISO4NOC:ISO4S4) = 1.0d-30

!     assume 2% of anthropogenic sulfur emission as H2SO4, 98% as SO2
      srcpar(:)     = semiss(:,IFSO2) * 0.02d0 * CMPM3      ! #molecules/m3/s
!     srcpar(:)     = 1.0d-30                               ! no direct emission

      where(srcgas(:)       < 1.0d-30) srcgas(:)     = 1.0d-30
      where(srcmas(:)       < 1.0d-30) srcmas(:)     = 1.0d-30
      where(srcpar(:)       < 1.0d-30) srcpar(:)     = 1.0d-30

      so4gas(:)     = cp(:,ISO4G) * CMPM3                   ! #molecules/m3

      so4aer(1,1,:) = cp(:,ISO4M1)                            & ! #molecules/m3
     &              * wtmair/wtmso4*zmair(:)*CMPM3
      so4aer(2,1,:) = cp(:,ISO4N1)                            & ! #particles/m3
     &              * gmair(:)*CMPM3/GPKG
      so4aer(1,2,:) = cp(:,ISO4M2)                            & ! #molecules/m3
     &              * wtmair/wtmso4*zmair(:)*CMPM3
      so4aer(2,2,:) = cp(:,ISO4N2)                            & ! #particles/m3
     &              * gmair(:)*CMPM3/GPKG

      do ijk = 1, itloop

        so4non(1:nnon,ijk) = cp(ijk,ISO4NOC:ISO4S4)           & ! #molecules/m3
     &                     * wtmair/wtmso4*zmair(ijk)*CMPM3

        aernon(1,1:nnon,ijk) = cp(ijk,INOC:ISSLT4)            & ! kg/m3
     &                       * gmair(ijk)*CMPM3/GPKG

      end do

!          =========
      call umaerosol  &
!          =========
     &     (dtchem, itloop, tgcm, press, relh, cfac,  &
     &      so4gas, srcgas, srcmas, srcpar,  &
     &      so4aer, so4non, aernon)


!     give back concentrations

      cp(:,ISO4G)  = so4gas(:) / CMPM3                      ! molecules/cm3

      cp(:,ISO4M1) = so4aer(1,1,:)                            & ! kg/kg
     &             / (wtmair/wtmso4*zmair(:)*CMPM3)
      cp(:,ISO4N1) = so4aer(2,1,:)                            & ! #particles/kg
     &             / (gmair(:)*CMPM3/GPKG)
      cp(:,ISO4M2) = so4aer(1,2,:)                            & ! kg/kg
     &             / (wtmair/wtmso4*zmair(:)*CMPM3)
      cp(:,ISO4N2) = so4aer(2,2,:)                            & ! #particles/kg
     &             / (gmair(:)*CMPM3/GPKG)

      do ijk = 1, itloop

        cp(ijk,ISO4NOC:ISO4S4) = so4non(1:nnon,ijk)           & ! kg/kg
     &                 / (wtmair/wtmso4*zmair(ijk)*CMPM3)

      end do


!     *****************************************************************
!     * Impose the clear plus cloudy sky concentration changes on
!       the grid cells.
!     *****************************************************************

!          ====================
      call Do_Const_Sulf_Update  &
!          ====================
     &  (itloop, cp, dcp)


!     *****************************************************************
!     * Calculate the mass transfer of SO2 vapor to water on the
!       surface of particles.
!     * The SO2 will appear on the particles as SO4.
!     *****************************************************************
#ifndef NONON
!          =================
!     call Do_Const_So2toPar
!          =================
!    &  (dtchem, itloop, zmair,
!    &   massc, pr_sulf_src,
!    &   rtcd, rtcs, cp)
#endif

!.... advance the time

!     tcl = tcl + dtclc
!                                 ========
!     if (tcl .lt. dtchem - 0.01d0) go to 20
!                                 ========

!     *****************************************************************
!     * Give back to IMPACT data
!     *****************************************************************

      do ic = 1, NSP
        xr(:, ic) = cp(:, ic)
      end do


      return

      end

