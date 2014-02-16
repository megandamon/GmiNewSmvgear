
!=============================================================================
!
! $Id: emiss_soil.F90,v 1.12 2013-08-28 20:54:52 jkouatch Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
!   Original code from:
!     Harvard tropospheric emissions module for 3D applications;
!       by Yuhang Wang, Gerry Gardner, and Prof. Daniel Jacob
!       of Harvard University (Release V2.1)
!
! FILE
!   emiss_soil.F
!
! ROUTINES
!   Calc_Canopy_Nox
!   Precip_Frac
!   Soil_Nox
!   Pulseit
!   Soil_Base
!   Soil_Can_Red_Fac
!   Soil_Fert_Fac
!   Soil_Temp_Fac
!   Soil_Type
!
! HISTORY
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calc_Canopy_Nox
!
! DESCRIPTION
!   This routine calculates canopynox using a resistance-in-series model
!   (Harvard Model / Depvel Version 3.2:  5/27/97).
!
!   Literature cited:  see drydep.F.
!
! ARGUMENTS
!   ino_num   : index of NO in const
!   lsnow     : array of flags that indicate land, water, or ice
!               (1 => water, 2 => land, 3 => ice)
!   ireg      : number of land types in a grid square
!   iland     : land type id in grid square for ireg land types
!   iuse      : fraction of grid box area occupied by land type (mil^-1?)
!   mw        : array of species' molecular weights (g/mol)
!   cfrac     : fractional cloud cover
!   radiat    : solar radiation (W/m^2)
!   suncos    : cosines of the solar zenith angle
!   tempk     : surface air temperature (degK)
!   xlai      : leaf area index of land type for current month
!   canopynox : deposition rate constant for NOx
!
!-----------------------------------------------------------------------------

      subroutine Calc_Canopy_Nox  &
     &  (ino_num, lsnow, ireg, iland, iuse, mw, cfrac, radiat,  &
     &   suncos, tempk, xlai, canopynox, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, num_species)

      use GmiResistance_mod, only : CanopyResistance, SurfaceResistance

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_emiss_constants.h"
#     include "gmi_drydep_data.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, num_species
      integer :: ino_num
      integer :: lsnow    (i1:i2, ju1:j2)
      integer :: ireg     (i1:i2, ju1:j2)
      integer :: iland    (i1:i2, ju1:j2, NTYPE)
      integer :: iuse     (i1:i2, ju1:j2, NTYPE)
      real*8  :: mw       (num_species)
      real*8  :: cfrac    (i1:i2, ju1:j2)
      real*8  :: radiat   (i1:i2, ju1:j2)
      real*8  :: suncos   (i1:i2, ju1:j2)
      real*8  :: tempk    (i1:i2, ju1:j2)
      real*8  :: xlai     (i1:i2, ju1:j2, NTYPE)
      real*8  :: canopynox(i1:i2, ju1:j2, NTYPE)

!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter ::  &
     &  F0_NO    = 0.0d0,     & ! reactivity factor for oxidation of biological
                            ! substances
     &  HSTAR_NO = 1.9d-3,    & ! Henry's Law constant for NO
     &  PRESS    = 1.5d5


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: idep1
      integer :: il, ij
      integer :: iolson
      integer :: ldt

      real*8  :: rdc
      real*8  :: rix
      real*8  :: rsurfc
      real*8  :: rt
      real*8  :: tempc1, tempk1

      real*8  :: rac (NTYPE)
      real*8  :: rclo(NTYPE)
      real*8  :: rcls(NTYPE)
      real*8  :: rgso(NTYPE)
      real*8  :: rgss(NTYPE)
      real*8  :: ri  (NTYPE)
      real*8  :: rlu (NTYPE)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Calc_Canopy_Nox called by ', loc_proc
      end if


      do ij = ju1, j2
        do il = i1, i2

          tempk1 = tempk(il,ij)
          tempc1 = tempk1 + ABS_ZERO

!         --------------------------------------------------------------------
!         Compute bulk surface resistance for gases.  Adjust external surface
!         resistances for temperature; from Wesely [1989], expression given in
!         text on p. 1296.  There is no evidence that the resistance continues
!         to increase at temperatures below -18 C, so at colder temperatures
!         hold the resistance fixed.
!         --------------------------------------------------------------------

          if (tempc1 < -18.0d0) then
            rt = 1.2d9
          else
            rt = 1000.0d0 * Exp (-tempc1 - 4.0d0)
          end if

!         --------------------------------------------------------------------
!         Get surface resistances; loop over land types, ldt.  The land types
!         within each grid square are defined using the Olson land-type
!         database.  Each of the Olson land types is assigned a corresponding
!         "deposition land type" with characteristic values of surface
!         resistance components.  There are 74 Olson land-types but only 11
!         deposition land-types (i.e., many of the Olson land types share the
!         same deposition characteristics).  Surface resistance components for
!         the "deposition land types" are from Wesely [1989] except for
!         tropical forests [Jacob and Wofsy, 1990] and for tundra [Jacob,
!         et al., 1992].  All surface resistance components are normalized to
!         a leaf area index of unity.  Olson land types, deposition land
!         types, and surface resistance components are read from an input data
!         file; check that file for further details.
!         --------------------------------------------------------------------

!         ================================
          LDTLOOP: do ldt = 1, ireg(il,ij)
!         ================================

!                                     =============
            if (iuse(il,ij,ldt) == 0) cycle LDTLOOP
!                                     =============

            if (lsnow(il,ij) == 3) then  ! snow or ice
              idep1  = 1
            else
              iolson = iland(il,ij,ldt) + 1
              idep1  = idep(iolson)
            end if

!           =======================
            call SurfaceResistance  &
!           =======================
     &        (idep1, rac(ldt), rclo(ldt), rcls(ldt), rgso(ldt),  &
     &         rgss(ldt), ri(ldt), rlu(ldt), rt, tempc1,  &
     &         cfrac(il,ij), radiat(il,ij), suncos(il,ij),  &
     &         xlai(il,ij,ldt), rix)

!           ----------------------------------------------------------------
!           Compute aerodynamic resistance to lower elements in lower part
!           of the canopy or structure, assuming level terrain; equation (5)
!           of Wesely [1989].
!           ----------------------------------------------------------------

            rdc =  &
     &        100.0d0 *  &
     &        (1.0d0 + (1000.0d0 / (radiat(il,ij) + 10.0d0)))

!           ---------------------------------------------------------------
!           Species-dependent corrections to resistances are from equations
!           (6)-(9) of Wesely [1989].
!           ---------------------------------------------------------------

!           ======================
            call CanopyResistance  &
!           ======================
     &        (rdc, rix, PRESS, tempk1, F0_NO, HSTAR_NO, mw(ino_num),  &
     &         rac(ldt), rclo(ldt), rcls(ldt), rgso(ldt),  &
     &         rgss(ldt), rlu(ldt), rsurfc)

            canopynox(il,ij,ldt) =  1.0d0 / rsurfc

!         ==============
          end do LDTLOOP
!         ==============

        end do
      end do


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Precip_Frac
!
! DESCRIPTION
!   This routine computes the fraction of a grid box that is actually
!   precipitating, along with the precipitation rate.
!
!   This version of Precip_Frac replaces Yuhang Wang's original version, as
!   used in the Harvard code prior to 10/18/99.
!
!  Reference:
!    Liu, H. Y., D. J. Jacob, I. Bey, R. M. Yantosca, and D. M. Koch,
!    Three-dimensional simulation of $210Pb$ and $7Be$ in the Harvard-DAO
!    tropospheric chemistry model, Eos Trans. AGU, 80 (17), S32, 1999a.
!
! ARGUMENTS
!   INPUT:
!     preacc : DAO total      precipitation at ground (mm/day)
!     precon : DAO convective precipitation at ground (mm/day)
!   OUTPUT:
!     frac   : fraction of grid box undergoing precipitation (unitless)
!     rate   : rate of precipitation for grid box (i,j)      (mm/day)
!
!-----------------------------------------------------------------------------

      subroutine Precip_Frac  &
     &  (preacc, precon, frac, rate, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, ju1_gl, j2_gl)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, ju1_gl, j2_gl
      real*8  :: preacc(i1:i2, ju1:j2)
      real*8  :: precon(i1:i2, ju1:j2)
      real*8  :: frac  (i1:i2, ju1:j2)
      real*8  :: rate  (i1:i2, ju1:j2)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij
      integer :: jstart, jend

      real*8  :: frac_convec
      real*8  :: frac_lsprecip


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Precip_Frac called by ', loc_proc
      end if


!     ----------------------------------------------
!     For the polar boxes there is no precipitation;
!     set rate = 0, frac = 0 for those boxes.
!     ----------------------------------------------

      if (ju1 == ju1_gl) then  ! South Pole
        frac(:,ju1) = 0.0d0
        rate(:,ju1) = 0.0d0
        jstart = ju1 + 1
      else
        jstart = ju1
      end if

      if (j2 == j2_gl) then    ! North Pole
        frac(:,j2) = 0.0d0
        rate(:,j2) = 0.0d0
        jend = j2 - 1
      else
        jend = j2
      end if


      do ij = jstart, jend
        do il = i1, i2

!         ---------------------------------------------------------------
!         Large-scale precipitation at (il,ij) =
!           preacc(il,ij) - precon(il,ij).
!         If there is large-scale precipitation at grid box (il,ij), then
!         assume that it covers 7% of the area of grid box(il,ij).
!         ---------------------------------------------------------------

          if ((preacc(il,ij) - precon(il,ij)) > 0.0d0) then
            frac_lsprecip = 7.0d-2
          else
            frac_lsprecip = 0.0d0
          end if

!         ----------------------------------------------------------------
!         Convective precipitation at (il,ij) = precon(il,ij); if there is
!         convective precipitation at (il,ij), then assume that it covers
!         0.3% of the area of grid box (il,ij).
!         ----------------------------------------------------------------

          if (precon(il,ij) > 0.0d0) then
            frac_convec = 3.0d-3
          else
            frac_convec = 0.0d0
          end if

!         -------------------------------------------------------------------
!         frac = total fraction of grid box (il,ij) covered by precipitation;
!         the possible values of frac are:  0.0%, 0.3%, 7.0%, or 7.3%.
!         -------------------------------------------------------------------

          frac(il,ij) = frac_lsprecip + frac_convec

!         -----------------------------------------------------------
!         rate = total precipitation rate in mm/day, adjusted for the
!         fraction of the grid box that is precipitating.
!         -----------------------------------------------------------

          if (frac(il,ij) > 0.0d0) then
            rate(il,ij) = preacc(il,ij) / frac(il,ij)
          else
            rate(il,ij) = 0.0d0
          end if

        end do
      end do


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Soil_Nox
!
! DESCRIPTION
!   This routine calculates the total NOx emission from the soil.
!
!   Based on Yienger and Levy [1995];
!   see Wang et al [1998]: Global Simulation of Tropospheric
!   O3-NOx-hydrocarbon; JGR Vol 103, pages 10713-10725.
!
! ARGUMENTS
!   nymd      : year/month/day  (YYYYMMDD)
!   ireg      : number of land types in a grid square
!   iland     : land type id in grid square for ireg land types
!   iuse      : fraction of grid box area occupied by land type (mil^-1?)
!   indexsoil : i,j of the grid
!   nconsoil  : Olson -> soil type
!   tdt       : model time step (s)
!   ylmid     : latitude        (deg)
!   frac      : precipitating fraction of grid square
!   rate      : local precipitation rate, i.e., grid-scale precipitation
!               rate divided by frac (mm/day)
!   radiat    : solar radiation (W/m^2)
!   tempk     : temperature     (degK)
!   windsqr   : surface wind speed squared ((m/s)^2)
!   canopynox : deposition rate constant for NOx
!   xlai      : leaf area index of land type for current month
!   soilfert  : fertilizers (ng N/m^2/s)
!   soilprep  : two months of observed precipitation (mm/day/box)
!   soilpuls  : tracking of wet/dry & three types of pulsing (Y&L, 94);
!     soilpuls(1,mm)    => flag for wet (-1) or dry (+1) soil; only dry soil
!                          is subject to pulsing
!     soilpuls(nn+1,mm) => fraction of grid square affected by equivalent
!                          fresh pulsing (e.g., sprinkle, shower, heavy rain)
!   xsoilnox  : soil NOx (molec/cm^2/s)
!
!-----------------------------------------------------------------------------

      subroutine Soil_Nox  &
     &  (nymd, ireg, iland, iuse, indexsoil, nconsoil, tdt, ylmid,  &
     &   frac, rate, radiat, tempk, windsqr, canopynox, xlai,  &
     &   soilfert, soilprep, soilpuls, xsoilnox, soil_day_save, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, ju1_gl, j2_gl)

      use GmiTimeControl_mod, only : GmiSplitDateTime
      use m_set_NLANDHAR    , only : NLANDHAR

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_emiss_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, ju1_gl, j2_gl
      integer :: nymd
      integer :: ireg     (i1:i2, ju1:j2)
      integer :: iland    (i1:i2, ju1:j2, NTYPE)
      integer :: iuse     (i1:i2, ju1:j2, NTYPE)
      integer :: indexsoil(2, NLANDHAR)
      integer :: nconsoil (NVEGTYPE)
      real*8  :: tdt
      real*8  :: ylmid    (ju1_gl:j2_gl)
      real*8  :: frac     (i1:i2, ju1:j2)
      real*8  :: rate     (i1:i2, ju1:j2)
      real*8  :: radiat   (i1:i2, ju1:j2)
      real*8  :: tempk    (i1:i2, ju1:j2)
      real*8  :: windsqr  (i1:i2, ju1:j2)
      real*8  :: canopynox(i1:i2, ju1:j2, NTYPE)
      real*8  :: xlai     (i1:i2, ju1:j2, NTYPE)
      real*8  :: soilfert (NLANDHAR)
      real*8  :: soilprep (2, NLANDHAR)
      real*8  :: soilpuls (NPULSE+1, NLANDHAR)
      real*8  :: xsoilnox (i1:i2, ju1:j2)
      integer, intent(inOut) :: soil_day_save


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8, external  :: Pulseit
      real*8, external  :: Soil_Base
      real*8, external  :: Soil_Can_Red_Fac
      real*8, external  :: Soil_Fert_Fac
      real*8, external  :: Soil_Temp_Fac


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: iday, imonth
      integer :: il, ij, it
      integer :: mm, nn
      integer :: ydummy

      real*8  :: riuse
      real*8  :: rpulse
      real*8  :: tempc1


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Soil_Nox called by ', loc_proc
      end if


!     ====================
      call GmiSplitDateTime  &
!     ====================
     &  (nymd, ydummy, imonth, iday)


!     ==============
      call Soil_Type  &
!     ==============
     &  (iday, soilprep, soilpuls, soil_day_save, pr_diag, loc_proc)


      do mm = 1, NLANDHAR

        il = indexsoil(1,mm)
        ij = indexsoil(2,mm)

        if ((il >=  i1) .and. (il <= i2) .and.  &
     &      (ij >= ju1) .and. (ij <= j2)) then

!         ---------------------------------------------------
!         Pulsing factor "function Pulseit ( )";
!         ECO system dependent;
!         Temperature factor "function Soil_Temp_Fac ( )";
!         Base emission with fertilization;
!         Canopy reduction;
!         Soil NOx emissions (watch out for trop. evergreen).
!         ---------------------------------------------------

          rpulse =  &
     &      Pulseit (mm, tdt, frac(il,ij), rate(il,ij), soilpuls)

          tempc1 = tempk(il,ij) + ABS_ZERO

          do it = 1, ireg(il,ij)

            nn    = nconsoil(iland(il,ij,it)+1)
            riuse = iuse(il,ij,it)

            xsoilnox(il,ij) =  &
     &        xsoilnox(il,ij) +  &
     &        ((Soil_Temp_Fac (mm, nn, tempc1, soilpuls) *  &
     &          Soil_Base     (mm, nn, soilprep, soilpuls, rpulse)) +  &
     &         Soil_Fert_Fac  (ij, mm, nn, imonth, ylmid, soilfert, ju1_gl, j2_gl)) *  &
     &        (1.0d0 -  &
     &         Soil_Can_Red_Fac  &
     &           (il, ij, it, nn, radiat, windsqr, canopynox, xlai, i1, i2, ju1, j2)) *  &
     &        riuse /  &
     &        1000.0d0

          end do

        end if

      end do


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Pulseit
!
! DESCRIPTION
!   This routine calculates the increase (or "pulse") of soil NOx
!   emission due to precipitation falling over a dry grid square and
!   JGR 100, 11,447-11,464, 1995.
!
! ARGUMENTS
!   mm       : land grid-square index
!   tdt      : model time step (s)
!   frac1    : precipitating fraction of grid square
!   rate1    : local precipitation rate, i.e., grid-scale precipitation
!              rate divided by frac1 (mm/day)
!   soilpuls : tracking of wet/dry & three types of pulsing (Y&L, 94);
!     soilpuls(1,mm)    => flag for wet (-1) or dry (+1) soil; only dry soil
!                          is subject to pulsing
!     soilpuls(nn+1,mm) => fraction of grid square affected by equivalent
!                          fresh pulsing (e.g., sprinkle, shower, heavy rain)
!
!-----------------------------------------------------------------------------

      function Pulseit (mm, tdt, frac1, rate1, soilpuls)

      use m_set_NLANDHAR, only : NLANDHAR

      implicit none

#     include "gmi_emiss_constants.h"
#     include "gmi_time_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: mm
      real*8  :: tdt
      real*8  :: frac1
      real*8  :: rate1
      real*8  :: soilpuls(NPULSE+1, NLANDHAR)

!     ----------------------
!     Function declarations.
!     ----------------------

      real*8  :: Pulseit

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: nn

      real*8  :: area      ! fraction of grid square affected by pulsing
      real*8  :: rfac
      real*8  :: rsecpdy
      real*8  :: tdt_days

      real*8  ::  exp_fac(NPULSE)    ! pulsing decay per time step (day^-1)
      !logical, save ::  first = .true.
      !real*8, save  ::  exp_fac(NPULSE) ! pulsing decay per time step (day^-1)

!     ----------------
!     Begin execution.
!     ----------------

      Pulseit = 0.0d0

      rsecpdy  = SECPDY

      tdt_days = tdt / rsecpdy

!!     ==========
!      if (first) then
!!     ==========
!!
!        first = .false.

        exp_fac(1) = Exp (-PULSE_DECAY(1) * tdt_days)
        exp_fac(2) = Exp (-PULSE_DECAY(2) * tdt_days)
        exp_fac(3) = Exp (-PULSE_DECAY(3) * tdt_days)

!      end if

      if (soilpuls(1,mm) <= 0.0d0) then  ! wet, no pulsing

        Pulseit = 1.0d0

      else                               ! dry, subject to pulsing

        do nn = 1, NPULSE

          if (soilpuls(nn+1,mm) < 1.0d-3) then  ! no pulsing, assume evap.

            soilpuls(nn+1,mm) = 0.0d0

          else        ! pulse from previous time step decays exponentially

            soilpuls(nn+1,mm) = soilpuls(nn+1,mm) * exp_fac(nn)

          end if

        end do

!       -------------------------------------------------------------------
!       Determine if a new pulse is to be applied to the grid square due
!       to precipitation over the current time step.  The pulse is applied
!       to the grid square fraction, frac1, experiencing precipitation.
!       Assume a characteristic 1-day duration for precipitation in a given
!       subgrid area of the grid square, so that the full extent of pulsing
!       (PULSE_FAC) is realized over 24 hours; for a model time step of
!       tdt seconds, reduce the pulsing by a factor tdt/SECPDY.
!       -------------------------------------------------------------------

        rfac = frac1 * tdt_days

        if      ((rate1 >=  1.0d0) .and.  &
     &           (rate1 <   5.0d0)) then  ! sprinkle

          soilpuls(2,mm) = soilpuls(2,mm) + rfac

        else if ((rate1 >=  5.0d0) .and.  &
     &           (rate1 <  15.0d0)) then  ! shower

          soilpuls(3,mm) = soilpuls(3,mm) + rfac

        else if (rate1  >= 15.0d0) then   ! heavy rain

          soilpuls(4,mm) = soilpuls(4,mm) + rfac

        end if

!       ------
!       Scale.
!       ------

        area = 0.0d0

!       ---------------------------------------------------------------------
!       Add up the contributions of the different types of pulses to obtain
!       the total pulsing multiplicative factor Pulseit; PULSE_FAC is the
!       multiplicative factor for fresh pulsing of each type.  Also determine
!       the fractional grid square area, area, affected by pulsing.  Assume
!       that the area occupied by the different pulses is additive, i.e.,
!       that successive pulses apply to different areas of the grid square
!       and that the area co-occupied by a pulse decreases as the pulsing
!       decays.  If the resulting area is in excess of unity then the pulsing
!       must be scaled back to the grid square area.  If the area is less
!       than unity then we have to account for non-pulsing emissions from the
!       (1-area) non-pulsing fraction of the grid square.
!       ---------------------------------------------------------------------

        do nn = 1, NPULSE

          Pulseit = Pulseit + (PULSE_FAC(nn) * soilpuls(nn+1,mm))

          area    = area + soilpuls(nn+1,mm)

        end do

        if (area < 1.0d0) then

          Pulseit = Pulseit + 1.0d0 - area

        else

          Pulseit = Pulseit / area

          do nn = 1, NPULSE

            soilpuls(nn+1,mm) = soilpuls(nn+1,mm) / area

          end do

        end if

      end if


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Soil_Base
!
! DESCRIPTION
!   This routine updates the land fractions of three types of pulsings;
!   updates only if within the calculated window.
!
! ARGUMENTS
!   mm       : land grid-square index
!   nn       : soil type
!   soilprep : two months of observed precipitation (mm/day/box)
!   soilpuls : tracking of wet/dry & three types of pulsing (Y&L, 94);
!     soilpuls(1,mm)    => flag for wet (-1) or dry (+1) soil; only dry soil
!                          is subject to pulsing
!     soilpuls(nn+1,mm) => fraction of grid square affected by equivalent
!                          fresh pulsing (e.g., sprinkle, shower, heavy rain)
!   rpulse   : pulsing rate (units?)
!
!-----------------------------------------------------------------------------

      function Soil_Base  &
     &  (mm, nn, soilprep, soilpuls, rpulse)

      use m_set_NLANDHAR    , only : NLANDHAR

      implicit none

#     include "gmi_emiss_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: mm, nn
      real*8  :: soilprep(2, NLANDHAR)
      real*8  :: soilpuls(NPULSE+1, NLANDHAR)
      real*8  :: rpulse


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8  :: Soil_Base


!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter ::  &
     &  UNITCONV = 4.3d9  ! ng N/m^2/s -> molec/cm^2/s


!     ----------------
!     Begin execution.
!     ----------------

      Soil_Base = 0.0d0


      if (nn == 1) then       ! desert

        Soil_Base = 0.0d0

      else if (nn == 2) then  ! tropical rain forest

        if (soilprep(2,mm) > 1.0d0) then  ! wet

          Soil_Base = SOIL_AW(2)

        else                              ! dry

          Soil_Base = SOIL_AD(2)

        end if

      else if ((nn ==8) .or. (nn == 9)) then

        Soil_Base = SOIL_AW(nn)

        if (nn == 9) then

          Soil_Base = Soil_Base / 30.0d0

        end if

      else  ! other

        if (soilpuls(1,mm) > 0.0d0) then  ! dry

          Soil_Base = SOIL_AD(nn) * rpulse

        else                              ! wet

          Soil_Base = SOIL_AW(nn)

        end if

      end if


!     --------------
!     Convert units.
!     --------------

      Soil_Base = Soil_Base * UNITCONV


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Soil_Can_Red_Fac
!
! DESCRIPTION
!   This routine calculates the canopy reduction factor.
!
!   Wang et al.: [1998] JGR vol. 103 p10713-10725
!
! ARGUMENTS
!   il : longitude index
!   ij : latitude  index
!   it : number in vegetation type of the grid
!   nn : soil type
!   radiat    : solar radiation (W/m^2)
!   windsqr   : surface wind speed squared ((m/s)^2)
!   canopynox : deposition rate constant for NOx
!   xlai      : leaf area index of land type for current month
!
!-----------------------------------------------------------------------------

      function Soil_Can_Red_Fac  &
     &  (il, ij, it, nn, radiat, windsqr, canopynox, xlai, &
     &  i1, i2, ju1, j2)

      implicit none

#     include "gmi_emiss_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2
      integer :: il, ij
      integer :: it, nn
      real*8  :: radiat   (i1:i2, ju1:j2)
      real*8  :: windsqr  (i1:i2, ju1:j2)
      real*8  :: canopynox(i1:i2, ju1:j2, NTYPE)
      real*8  :: xlai     (i1:i2, ju1:j2, NTYPE)


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8  :: Soil_Can_Red_Fac


!     -----------------------
!     Parameter declarations.
!     -----------------------

!c?
!     -----------------------------------------------------------------
!     Coefficient ALPHA (2.8d-2, 5.6d-3) day, night canopy ventilation;
!     time of 1 hour day, 5 hour night;
!     VFDAY, VFNIGHT -> ALPHA scaled.
!     -----------------------------------------------------------------

      real*8, parameter ::  &
     &  VFDAY   = 1.0d-2,    & ! ventilation velocity in day   (m/s)
     &  VFNIGHT = 0.2d-2   ! ventilation velocity in night (m/s)


!     ----------------------
!     Variable declarations.
!     ----------------------

      real*8  ::  vfnew    ! ventilation rate constant for NOx (m/s)


!     ----------------
!     Begin execution.
!     ----------------

      Soil_Can_Red_Fac = 0.0d0


      if (radiat(il,ij) > 0.0d0) then
        vfnew = VFDAY
      else
        vfnew = VFNIGHT
      end if


      if ((xlai     (il,ij,it) > 0.0d0) .and.  &
     &    (canopynox(il,ij,it) > 0.0d0)) then

        vfnew =  &
     &    vfnew *  &
     &    Sqrt (windsqr(il,ij) / 9.0d0 * 7.0d0 / xlai(il,ij,it)) *  &
     &    (SOIL_EXT(2) / SOIL_EXT(nn))

        Soil_Can_Red_Fac =  &
     &    canopynox(il,ij,it) / (canopynox(il,ij,it) + vfnew)

      end if


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Soil_Fert_Fac
!
! DESCRIPTION
!   This routine calculates the NOx emission from fertilizer.
!
! ARGUMENTS
!   ij            : latitude index
!   mm            : land grid-square index
!   nn            : soil type
!   current_month : current month
!   ylmid         : latitude    (deg)
!   soilfert      : fertilizers (ng N/m^2/s)
!
!-----------------------------------------------------------------------------

      function Soil_Fert_Fac  &
     &  (ij, mm, nn, current_month, ylmid, soilfert, ju1_gl, j2_gl)

      use m_set_NLANDHAR    , only : NLANDHAR

      implicit none

#     include "gmi_emiss_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: ju1_gl, j2_gl
      integer :: ij
      integer :: mm, nn
      integer :: current_month
      real*8  :: ylmid   (ju1_gl:j2_gl)
      real*8  :: soilfert(NLANDHAR)


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8  :: Soil_Fert_Fac


!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter ::  &
     &  UNITCONV = 4.3d9  ! ng N/m^2/s -> molec/cm^2/s


!     ----------------
!     Begin execution.
!     ----------------

      Soil_Fert_Fac = 0.0d0


!                                    ======
      if ((nn /= 8) .and. (nn /= 9)) return
!                                    ======


!c?
!     -------------------------------------------------
!     Agriculture/fertilizer has to be in "ng N/m^2/s";
!     including the boxes on 30 degree lines.
!     -------------------------------------------------

      if (ylmid(ij) >  28.0d0) then

        if ((current_month >= 5) .and.  &
     &      (current_month <= 8)) then    ! summer in N. Hemis.

          Soil_Fert_Fac = soilfert(mm)

        else                              ! winter in N. Hemis.

          Soil_Fert_Fac = 0.0d0

        end if

      else if (ylmid(ij) > -28.0d0) then  ! tropics

        Soil_Fert_Fac = soilfert(mm)

      else

        if ((current_month <=  2) .or.  &
     &      (current_month >= 11)) then   ! summer in S. Hemis.

          Soil_Fert_Fac = soilfert(mm)

        else                              ! winter in S. Hemis.

          Soil_Fert_Fac = 0.0d0

        end if

      end if


      if (nn == 9) then
        Soil_Fert_Fac = Soil_Fert_Fac / 30.0d0
      end if


!     --------------
!     Convert units.
!     --------------

      Soil_Fert_Fac = Soil_Fert_Fac * UNITCONV
!     This is only for preindustrial runs - tag GMI_PREIND_CORRECT     Soil_Fert_Fac = 0.


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Soil_Temp_Fac
!
! DESCRIPTION
!   This routine calculates the soil temperature factor.
!
!   Yienger and Levy [1995] JGR 100, 11447-11464
!
! ARGUMENTS
!   mm       : land grid-square index
!   nn       : soil type
!   tempc1   : surface air temperature (degC)
!   soilpuls : tracking of wet/dry & three types of pulsing (Y&L, 94);
!     soilpuls(1,mm)    => flag for wet (-1) or dry (+1) soil; only dry soil
!                          is subject to pulsing
!     soilpuls(nn+1,mm) => fraction of grid square affected by equivalent
!                          fresh pulsing (e.g., sprinkle, shower, heavy rain)
!
!-----------------------------------------------------------------------------

      function Soil_Temp_Fac  &
     &  (mm, nn, tempc1, soilpuls)

      use m_set_NLANDHAR    , only : NLANDHAR

      implicit none

#     include "gmi_emiss_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: mm, nn
      real*8  :: tempc1
      real*8  :: soilpuls(NPULSE+1, NLANDHAR)


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8  :: Soil_Temp_Fac


!     ----------------------
!     Variable declarations.
!     ----------------------

      real*8  :: xtempc1


!     ----------------
!     Begin execution.
!     ----------------

      Soil_Temp_Fac = 0.0d0


      xtempc1 = tempc1


      if (nn <= 2) then  ! desert and rain forest

        Soil_Temp_Fac = 1.0d0

      else if ((soilpuls(1,mm) > 0.0d0) .and.  &
     &         (nn /= 8) .and. (nn /= 9)) then

!       -------------------------------------------------------
!       Dry:
!         surface temperature -> soil temperature;
!         convert the lowest model level air temp to soil temp;
!         based on observations of Johansson et. al. [1988];
!         add 5 degrees C to model temperature.
!       -------------------------------------------------------

        xtempc1 = xtempc1 + 5.0d0

        if (xtempc1 > 30.0d0) then       ! optimal

          Soil_Temp_Fac = 1.0d0

        else if (xtempc1 >  0.0d0) then  ! cold-linear

          Soil_Temp_Fac = xtempc1 / 30.0d0

        else

          Soil_Temp_Fac = 0.0d0

        end if

      else

!       ---------------------------------------------------------------------
!       Wet:
!         surface temperature -> soil temperature;
!         convert the lowest model level air temp to soil temp;
!         use the empirical relationships derived by Williams et al. [1992b];
!         ECO system dependent.
!       ---------------------------------------------------------------------

        xtempc1 = SOIL_T2(nn) + (SOIL_T1(nn) * xtempc1)

        if (xtempc1 >= 30.0d0) then       ! optimal

          Soil_Temp_Fac = 21.97d0

        else if (xtempc1 >= 10.0d0) then  ! exponential

          Soil_Temp_Fac = Exp (0.103d0 * xtempc1)

        else if (xtempc1 >   0.0d0) then  ! cold-linear

          Soil_Temp_Fac = 0.28d0 * xtempc1

        else

          Soil_Temp_Fac = 0.0d0

        end if

      end if


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Soil_Type
!
! DESCRIPTION
!   This routine determines whether soil is dry or wet for all land grid
!   boxes; updated daily.
!
! ARGUMENTS
!   iday     : day of the current month
!   soilprep : two months of observed precipitation (mm/day/box)
!   soilpuls : tracking of wet/dry & three types of pulsing (Y&L, 94);
!     soilpuls(1,mm)    => flag for wet (-1) or dry (+1) soil; only dry soil
!                          is subject to pulsing
!     soilpuls(nn+1,mm) => fraction of grid square affected by equivalent
!                          fresh pulsing (e.g., sprinkle, shower, heavy rain)
!
!-----------------------------------------------------------------------------

      subroutine Soil_Type (iday, soilprep, soilpuls, soil_day_save, &
                            pr_diag, loc_proc)

      use m_set_NLANDHAR    , only : NLANDHAR

      implicit none

#     include "gmi_emiss_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc
      integer :: iday
      real*8  :: soilprep(2, NLANDHAR)
      real*8  :: soilpuls(NPULSE+1, NLANDHAR)
      integer, intent(inOut) ::  soil_day_save


!     -----------------------
!     Parameter declarations.
!     -----------------------

      integer, parameter ::  IDAYS_TO_TEST = 14      ! number of days for pulse

      real*8,  parameter ::  WETSOIL       = 10.0d0  ! criteria for wet soil (mm);
                                ! above 10 mm for two weeks

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: mm, nn
      integer :: ncurr          ! number of days in current  month
      integer :: nprev          ! number of days in previous month

      real*8  :: rain           ! total rain (units?)
      real*8  :: rncurr
      real*8  :: rnprev


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Soil_Type called by ', loc_proc
      end if


      if (iday /= soil_day_save) then

        soil_day_save  = iday


        if (iday >= IDAYS_TO_TEST) then
          ncurr = IDAYS_TO_TEST
          nprev = 0
        else
          ncurr = iday
          nprev = IDAYS_TO_TEST - iday
        end if

        rncurr = ncurr
        rnprev = nprev


        do mm = 1, NLANDHAR

          rain = (soilprep(1,mm) * rnprev) + (soilprep(2,mm) * rncurr)

          if (rain > WETSOIL) then  ! wet

            soilpuls(1,mm) = -1.0d0

            do nn = 1, NPULSE

              soilpuls(nn+1,mm) = 0.0d0

            end do

          else                      ! dry

            soilpuls(1,mm) = 1.0d0

          end if

        end do


      end if


      return

      end

