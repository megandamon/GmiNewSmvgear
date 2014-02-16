
!=============================================================================
!
! $Id: emiss_harvard.F90,v 1.22 2013-08-28 20:54:52 jkouatch Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   emiss_harvard.F
!
! ROUTINES
!   Update_Emiss_Harvard
!   Do_Biogenic_Emiss
!   Do_Soil_Emiss
!   Add_Emiss_Harvard
!   Add_Semiss_Harvard_Inchem
!
! HISTORY
!   - November 23, 2004 * Jules Kouatchou
!     Added the variable "surf_emiss_out2" as argument of the routine
!     Add_Emiss_Harvard. It is used to produce surface emission diagnostics
!     for soil NO, monoterpenes CO, methanol CO and biogenic propene.
!   - August 19, 2005 * Jules Kouatchou
!     Set here the value of NLANDHAR using the function
!     NLANDHAR_expected.
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update_Emiss_Harvard
!
! DESCRIPTION
!   This routine updates the Harvard biogenic (isoprene, CO, & propene) and
!   soil emissions (NOx).
!
! HISTORY
!   * August 12, 2005 * Jules Kouatchou
!     Allocate the variables index_soil, soil_fert, soil_precip, and
!     soil_pulse the first time this routine is called.
!
! ARGUMENTS
!   iisoprene_num : index of isoprene in const
!   ino_num       : index of NO       in const
!   nhms          : hour/minute/sec  (HHMMSS)
!   nymd          : year/month/day   (YYYYMMDD)
!   lwi_flags     : array of flags that indicate land, water, or ice
!   tdt           : model time step  (s)
!   mw            : array of species' molecular weights (g/mol)
!   latdeg        : latitude         (deg)
!   londeg        : longitude        (deg)
!   mcor          : area of grid box (m^2)
!   tempk         : surface air temperature  (degK)
!   pardif        : Diffuse photosynthetically active radiation (0.35-0.70 um)
!   pardir        : Direct  photosynthetically active radiation (0.35-0.70 um)
!   radswg        : net downward shortwave radiation at ground (W/m^2)
!   surf_rough    : surface roughness        (m)
!   con_precip    : convective precipitation (mm/day)
!   tot_precip    : total precipitation      (mm/day)
!   ustar         : ustar                    (m/s)
!   max_cloud     : maximum overlap cloud fraction for LW
!   ran_cloud     : random  overlap cloud fraction for LW
!   emiss_isop    : isoprene    emissions    (kg/s)
!   emiss_monot   : monoterpene emissions    (kg/s)
!   emiss_nox     : NOx         emissions    (kg/s)
!
!-----------------------------------------------------------------------------

      subroutine Update_Emiss_Harvard  &
     &  (iisoprene_num, ino_num, lwi_flags, tdt, mw, &
     &   cosSolarZenithAngle, latdeg, nymd, mcor, &
     &   aefIsop, aefMbo, aefMonot, isoLai, isoLaiCurr, isoLaiPrev, &
     &   T_15_AVG, tempk, pardif, pardir, radswg, surf_rough, con_precip,  &
     &   tot_precip, ustar, fracCloudCover, &
     &   emiss_isop, emiss_monot, emiss_nox, &
     &   index_soil, ncon_soil, soil_fert, soil_precip, soil_pulse, &
     &   ireg, iland, iuse, convert_isop, convert_monot, &
     &   coeff_isop, base_isop, base_monot, xlai, &
     &   doMEGANemission, days_btw_m, soil_day_save, pr_diag, loc_proc, &
     &   i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl, ilong, num_species)

      use m_set_NLANDHAR, only : NLANDHAR
      use GmiEmissionMEGAN_mod, only : calcBiogenicMEGANemission

      implicit none

#     include "gmi_emiss_constants.h"
#     include "gmi_time_constants.h"
#     include "gmi_phys_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

     logical, intent(in) :: pr_diag, doMEGANemission
     integer, intent(in) :: loc_proc, nymd
     integer, intent(in) :: i1, i2, ju1, j2, k1, k2
     integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, ilong, num_species
     integer :: iisoprene_num, ino_num
     real*8  :: latdeg     (ju1_gl:j2_gl), cosSolarZenithAngle (i1:i2, ju1:j2)
     integer :: lwi_flags  (i1:i2, ju1:j2)
     real*8  :: tdt
      integer, intent(in) :: days_btw_m ! ! days between midmonths in the LAI data
      real*8 , intent(in) :: T_15_AVG            (i1:i2, ju1:j2)
      real*8 , intent(in) :: aefIsop             (i1:i2, ju1:j2)
      real*8 , intent(in) :: aefMbo              (i1:i2, ju1:j2)
      real*8 , intent(in) :: aefMonot            (i1:i2, ju1:j2)
      real*8 , intent(in) :: isoLai              (i1:i2, ju1:j2)
      real*8 , intent(in) :: isoLaiCurr          (i1:i2, ju1:j2)
      real*8 , intent(in) :: isoLaiPrev          (i1:i2, ju1:j2)
     real*8  :: mw         (num_species)  , mcor       (i1:i2, ju1:j2)
     real*8  :: tempk      (i1:i2, ju1:j2), radswg     (i1:i2, ju1:j2)
     real*8  :: pardif     (i1:i2, ju1:j2), pardir     (i1:i2, ju1:j2)
     real*8  :: surf_rough (i1:i2, ju1:j2), con_precip (i1:i2, ju1:j2)
     real*8  :: tot_precip (i1:i2, ju1:j2), ustar      (i1:i2, ju1:j2)
     real*8  :: fracCloudCover  (i1:i2, ju1:j2)
     real*8  :: emiss_isop (i1:i2, ju1:j2), emiss_monot(i1:i2, ju1:j2)
     real*8  :: emiss_nox  (i1:i2, ju1:j2)
     integer :: iland(i1:i2, ju1:j2, NTYPE), ireg(i1:i2, ju1:j2)
     integer :: iuse (i1:i2, ju1:j2, NTYPE)
     real*8  :: convert_isop (NVEGTYPE), convert_monot(NVEGTYPE)
     real*8  :: base_isop (i1:i2, ju1:j2, NTYPE), base_monot(i1:i2, ju1:j2, NTYPE)
     real*8  :: xlai      (i1:i2, ju1:j2, NTYPE)
     real*8  :: soil_fert(NLANDHAR), soil_precip(2,NLANDHAR), soil_pulse(NPULSE+1,NLANDHAR)
     integer :: ncon_soil (NVEGTYPE), index_soil(2, NLANDHAR)
     real*8  :: coeff_isop   (NPOLY)
     integer, intent(inOut) :: soil_day_save

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ij
      integer :: nsec_jan1

      real*8  :: decl, tdtinv
      real*8  :: rdistsq
      real*8  :: days, time

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_Emiss_Harvard called by ', loc_proc
      end if


       if (doMEGANemission) then
!         ======================
          call calcBiogenicMEGANemission  &
!         ======================
     &      (emiss_isop, emiss_monot, days_btw_m, &
     &       aefIsop, aefMbo, aefMonot, isoLai, isoLaiCurr, isoLaiPrev, &
     &       tempk, T_15_AVG, pardif, pardir, cosSolarZenithAngle, &
     &       pr_diag, loc_proc, i1, i2, ju1, j2)

          ! Perform unit conversions
          tdtinv = 1.0d0 / tdt

          emiss_isop = emiss_isop * tdtinv / ATOMSC_PER_MOLECISOP *  &
     &                        (mw(iisoprene_num) / AVOGAD) * KGPG

          emiss_monot = emiss_monot * tdtinv / ATOMSC_PER_MOLECMONOT *  &
     &                          (MWTMONOT / AVOGAD) * KGPG
       else
!         ======================
          call Do_Biogenic_Emiss  &
!         ======================
     &      (iisoprene_num, ireg, iland, iuse, tdt, mw, cosSolarZenithAngle, mcor,  &
     &       tempk, fracCloudCover, coeff_isop, convert_isop, convert_monot,  &
     &       xlai, base_isop, base_monot, emiss_isop, emiss_monot, &
     &       pr_diag, loc_proc, i1, i2, ju1, j2, num_species)
      end if


!     ==================
      call Do_Soil_Emiss  &
!     ==================
     &  (ino_num, nymd, lwi_flags, ireg, iland, iuse, index_soil,  &
     &   ncon_soil, tdt, mw, latdeg, cosSolarZenithAngle, mcor, tempk, radswg,  &
     &   surf_rough, con_precip, tot_precip, ustar, fracCloudCover,  &
     &   soil_fert, soil_precip, soil_pulse, xlai, emiss_nox, soil_day_save, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, ju1_gl, j2_gl, num_species)


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Biogenic_Emiss
!
! DESCRIPTION
!   This routine updates the Harvard biogenic emissions (isoprene, CO, &
!   propene).
!
! ARGUMENTS
!   iisoprene_num : index of isoprene in const
!   ireg   : number of land types in a grid square
!   iland  : land type id in grid square for ireg land types
!   iuse   : fraction of grid box area occupied by land type (mils?)
!   tdt    : model time step  (s)
!   mw     : array of species' molecular weights (g/mol)
!   cossza : cosines of the solar zenith angle
!   mcor   : area of grid box (m^2)
!   tempk  : surface air temperature (degK)
!   fracCloudCover    : fractional cloud cover
!   coeff_isop    : coefficients used for polynomial fit
!   convert_isop  : isoprene    emissions by landtype   (atomsC/cm^2 leaf/s)
!   convert_monot : monoterpene emissions by landtype   (atomsC/cm^2 leaf/s)
!   xlai          : leaf area index of land type for current month
!   base_isop     : baseline emissions for isoprene     (kgC/box/step?)
!   base_monot    : baseline emissions for monoterpenes (kgC/box/step?)
!   emiss_isop    : isoprene    emissions (kg/s)
!   emiss_monot   : monoterpene emissions (kg/s)
!
!-----------------------------------------------------------------------------

      subroutine Do_Biogenic_Emiss  &
     &  (iisoprene_num, ireg, iland, iuse, tdt, mw, cossza, mcor,  &
     &   tempk, fracCloudCover, coeff_isop, convert_isop, convert_monot,  &
     &   xlai, base_isop, base_monot, emiss_isop, emiss_monot, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, num_species)

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_emiss_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, num_species
      integer :: iisoprene_num
      integer :: ireg (i1:i2, ju1:j2)
      integer :: iland(i1:i2, ju1:j2, NTYPE)
      integer :: iuse (i1:i2, ju1:j2, NTYPE)
      real*8  :: tdt
      real*8  :: mw   (num_species)
      real*8  :: cossza       (i1:i2, ju1:j2)
      real*8  :: mcor         (i1:i2, ju1:j2)
      real*8  :: tempk        (i1:i2, ju1:j2)
      real*8  :: fracCloudCover   (i1:i2, ju1:j2)
      real*8  :: coeff_isop   (NPOLY)
      real*8  :: convert_isop (NVEGTYPE)
      real*8  :: convert_monot(NVEGTYPE)
      real*8  :: xlai         (i1:i2, ju1:j2, NTYPE)
      real*8  :: base_isop    (i1:i2, ju1:j2, NTYPE)
      real*8  :: base_monot   (i1:i2, ju1:j2, NTYPE)
      real*8  :: emiss_isop   (i1:i2, ju1:j2)
      real*8  :: emiss_monot  (i1:i2, ju1:j2)


!     ----------------------
!     Function declarations.
!     ----------------------

      real*8, external  :: Biogenic_Isop
      real*8, external  :: Biogenic_Monot


!     ----------------------
!     Variable declarations.
!     ----------------------

!      logical, save :: first = .true.

      integer :: il, ij

      integer :: iuse1(NTYPE)

      real*8  :: biemiss, bmemiss
      real*8  :: tdtinv

      real*8  :: base_isop1 (NTYPE)
      real*8  :: base_monot1(NTYPE)
      real*8  :: xlai1      (NTYPE)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Do_Biogenic_Emiss called by ', loc_proc
      end if


!!     ==========
!      if (first) then
!!     ==========
!
!        first = .false.
!
!!       ==================
!        call Biogenic_Base  &
!!       ==================
!     &    (ireg, iland, tdt, convert_isop, convert_monot, mcor,  &
!     &     base_isop, base_monot, pr_diag, loc_proc, i1, i2, ju1, j2)
!
!      end if


      tdtinv = 1.0d0 / tdt


      do ij = ju1, j2
        do il = i1, i2

          iuse1(:) = iuse(il,ij,:)

          base_isop1 (:) = base_isop (il,ij,:)
          base_monot1(:) = base_monot(il,ij,:)
          xlai1      (:) = xlai      (il,ij,:)


          biemiss =  &
!           =============
     &      Biogenic_Isop  &
!           =============
     &        (ireg(il,ij), iuse1(:), fracCloudCover(il,ij), cossza(il,ij),  &
     &         tempk(il,ij), coeff_isop(:), base_isop1(:), xlai1(:))

          emiss_isop(il,ij) =  &
     &      biemiss * tdtinv / ATOMSC_PER_MOLECISOP *  &
     &      (mw(iisoprene_num) / AVOGAD) * KGPG


          bmemiss =  &
!           ==============
     &      Biogenic_Monot  &
!           ==============
     &        (ireg(il,ij), iuse1(:), tempk(il,ij), base_monot1(:),  &
     &         xlai1(:))

          emiss_monot(il,ij) =  &
     &      bmemiss * tdtinv / ATOMSC_PER_MOLECMONOT *  &
     &      (MWTMONOT / AVOGAD) * KGPG

        end do
      end do


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Soil_Emiss
!
! DESCRIPTION
!   This routine updates the Harvard soil emissions (NOx)
!
! ARGUMENTS
!   ino_num       : index of NO in const
!   nymd          : year/month/day  (YYYYMMDD)
!   lwi_flags     : array of flags that indicate land, water, or ice
!   ireg          : number of land types in a grid square
!   iland         : land type id in grid square for ireg land types
!   iuse          : fraction of grid box area occupied by land type (mil^-1?)
!   index_soil    : i,j of the grid
!   ncon_soil     : Olson -> soil type
!   tdt           : model time step (s)
!   mw            : array of species' molecular weights (g/mol)
!   latdeg        : latitude (deg)
!   cossza        : cosines of the solar zenith angle
!   mcor          : area of grid box         (m^2)
!   tempk         : surface air temperature  (degK)
!   radswg        : net downward shortwave radiation at ground (W/m^2)
!   surf_rough    : surface roughness        (m)
!   con_precip    : convective precipitation (mm/day)
!   tot_precip    : total precipitation      (mm/day)
!   ustar         : ustar                    (m/s)
!   fracCloudCover    : fractional cloud cover
!   soil_fert     : fertilizers   (ng N/m^2/s)
!   soil_precip   : two months of observed precipitation (mm/day/box)
!   soil_pulse    : tracking of wet/dry & three types of pulsing (Y&L, 94)
!   xlai          : leaf area index of land type for month #1
!   emiss_nox     : NOx emissions (kg/s)
!
!-----------------------------------------------------------------------------

      subroutine Do_Soil_Emiss  &
     &  (ino_num, nymd, lwi_flags, ireg, iland, iuse, index_soil,  &
     &   ncon_soil, tdt, mw, latdeg, cossza, mcor, tempk, radswg,  &
     &   surf_rough, con_precip, tot_precip, ustar, fracCloudCover,  &
     &   soil_fert, soil_precip, soil_pulse, xlai, emiss_nox, soil_day_save,   &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, ju1_gl, j2_gl, num_species)

      use m_set_NLANDHAR, only : NLANDHAR

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_emiss_constants.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, ju1_gl, j2_gl, num_species
      integer :: ino_num
      integer :: nymd
      integer :: lwi_flags  (i1:i2, ju1:j2)
      integer :: ireg       (i1:i2, ju1:j2)
      integer :: iland      (i1:i2, ju1:j2, NTYPE)
      integer :: iuse       (i1:i2, ju1:j2, NTYPE)
      integer :: index_soil (2, NLANDHAR)
      integer :: ncon_soil  (NVEGTYPE)
      real*8  :: tdt
      real*8  :: mw         (num_species)
      real*8  :: latdeg     (ju1_gl:j2_gl)
      real*8  :: cossza     (i1:i2, ju1:j2)
      real*8  :: mcor       (i1:i2, ju1:j2)
      real*8  :: tempk      (i1:i2, ju1:j2)
      real*8  :: radswg     (i1:i2, ju1:j2)
      real*8  :: surf_rough (i1:i2, ju1:j2)
      real*8  :: con_precip (i1:i2, ju1:j2)
      real*8  :: tot_precip (i1:i2, ju1:j2)
      real*8  :: ustar      (i1:i2, ju1:j2)
      real*8  :: fracCloudCover (i1:i2, ju1:j2)
      real*8  :: soil_fert  (NLANDHAR)
      real*8  :: soil_precip(2, NLANDHAR)
      real*8  :: soil_pulse (NPULSE+1, NLANDHAR)
      real*8  :: xlai       (i1:i2, ju1:j2, NTYPE)
      real*8  :: emiss_nox  (i1:i2, ju1:j2)
      integer, intent(inOut) :: soil_day_save


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij

      real*8  :: cmpm2
      real*8  :: rfac
      real*8  :: windsp

!     ----------------------------------------------------------------------
!     frac_precip : fraction of grid box undergoing precipitation (unitless)
!     rate_precip : rate of precipitation for grid box (i,j) (mm/day)
!     windsp2     : surface wind speed squared ((m/s)^2)
!     xsoil_nox   : soil NOx (molec/cm^2/s)
!     ----------------------------------------------------------------------

      real*8  :: frac_precip(i1:i2, ju1:j2)
      real*8  :: rate_precip(i1:i2, ju1:j2)
      real*8  :: windsp2    (i1:i2, ju1:j2)
      real*8  :: xsoil_nox  (i1:i2, ju1:j2)

!     ---------------------------------------------
!     canopy_nox : deposition rate constant for NOx
!     ---------------------------------------------

      real*8  :: canopy_nox (i1:i2, ju1:j2, NTYPE)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Do_Soil_Emiss called by ', loc_proc
      end if


      frac_precip (:,:) = 0.0d0
      rate_precip (:,:) = 0.0d0
      windsp2     (:,:) = 0.0d0
      xsoil_nox   (:,:) = 0.0d0

      canopy_nox(:,:,:) = 0.0d0


      do ij = ju1, j2
        do il = i1, i2

!         --------------------------------------------------------------
!         Calculate the wind speed (m/s) (see Seinfeld, page 494, 1986).
!         --------------------------------------------------------------

          windsp = (ustar(il,ij) / 0.4d0) *  &
     &             Log (10.0d0 / surf_rough(il,ij))

          windsp2(il,ij) = windsp * windsp

        end do
      end do


!     ====================
      call Calc_Canopy_Nox  &
!     ====================
     &  (ino_num, lwi_flags, ireg, iland, iuse, mw, fracCloudCover,  &
     &   radswg, cossza, tempk, xlai, canopy_nox, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, num_species)


!soi
!soi  ----------------------------------------------------------------
!soi  If you want to run the original soil nox test case that was sent
!soi  to us, search for all the "csoi" comments in this file and make
!soi  the "necessary" changes.
!soi
!soi  Also change the following in dao2no_dao46_02.in:
!soi    tfinal_days = 1.0d0,
!soi    start_ymd   = 970601,
!soi    tdt         = 14400.0d0,
!soi    mw(4)       = 46.01d0,
!soi
!soi  Also comment out the error checks between the "cemi" markers in
!soi  gmi_namelist.F.
!soi  ----------------------------------------------------------------

!soi  con_precip(:,:) = 0.0d0
!soi  radswg    (:,:) = 0.0d0
!soi  tot_precip(:,:) = 0.0d0

!soi  tempk     (62,22)   = 300.25839d0
!soi  windsp2   (62,22)   =  68.315167385258405d0

!soi  canopy_nox(62,22,:) =   0.0d0
!soi  canopy_nox(62,22,1) =   7.16978099301216181d-4
!soi  canopy_nox(62,22,2) =   3.787423305420255d-4
!soi  canopy_nox(62,22,3) =   3.787423305420255d-4
!soi


!     ================
      call Precip_Frac  &
!     ================
     &  (tot_precip, con_precip, frac_precip, rate_precip, &
     &  pr_diag, loc_proc, i1, i2, ju1, j2, ju1_gl, j2_gl)


!     =============
      call Soil_Nox  &
!     =============
     &  (nymd, ireg, iland, iuse, index_soil, ncon_soil, tdt, latdeg,  &
     &   frac_precip, rate_precip, radswg, tempk, windsp2, canopy_nox,  &
     &   xlai, soil_fert, soil_precip, soil_pulse, xsoil_nox, soil_day_save, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, ju1_gl, j2_gl)


!soi
!soi  Write (6,*)
!soi  Write (6,*) 'On J90s =>'
!soi  Write (6,*)
!soi  Write (6,*)
!soi &  'Result should be => xsoil_nox(62,22) =  12721337142.50'
!soi  Write (6,*)
!soi  Write (6,*)
!soi &  'Result is        => xsoil_nox(62,22) = ', xsoil_nox(62,22)
!soi  Write (6,*)
!soi  call Stopcode (" ")
!soi


      cmpm2 = CMPM * CMPM

      do ij = ju1, j2
        do il = i1, i2

          rfac = (mw(ino_num) * KGPG) * (mcor(il,ij) * cmpm2) / AVOGAD

          emiss_nox(il,ij) = xsoil_nox(il,ij) * rfac

        end do
      end do


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Add_Emiss_Harvard
!
! DESCRIPTION
!   This routine adds the Harvard emissions to const:
!
!     1) Take emissions of kg/s and multiply by the time step to get total
!        kg of emission over time step.
!
!     2) Divide by mass of the zone to obtain emissions in terms of mixing
!        ratio; also multiply by the ratio of molecular weight of air to
!        molecular weight of the chemical emission to get volume mixing
!        ratio from mass mixing ratio.
!
!     3) Add emitted mixing ratio amount to existing mixing ratio of const.
!
! ARGUMENTS
!   pr_surf_emiss : should surface emissions be accumulated for output?
!   iisoprene_num : index of isoprene in const
!   ico_num       : index of CO       in const
!   ipropene_num  : index of propene  in const
!   ino_num       : index of NO       in const
!   tdt           : model time step    (s)
!   mw            : array of species' molecular weights (g/mol)
!   mcor          : surface area of grid box (m^2)
!   surf_emiss_out: accumulated surface emissions (kg/m^2/delta t output)
!   mass          : total mass of the atmosphere within each grid box (kg)
!   const         : species concentration, known at zone centers
!                   (mixing ratio)
!   emiss_isop    : isoprene   emissions (kg/s)
!   emiss_monot   : monterpene emissions (kg/s)
!   emiss_nox     : NOx        emissions (kg/s)
!
!-----------------------------------------------------------------------------

      subroutine Add_Emiss_Harvard  &
     &  (pr_surf_emiss, pr_emiss_3d, iisoprene_num, ico_num, ipropene_num,  &
     &   ino_num, tdt, mw, mcor, surf_emiss_out, surf_emiss_out2,  &
     &   emiss_3d_out, mass, concentration, emiss_isop, emiss_monot, emiss_nox, &
     &   do_ShipEmission, emiss_hno3, emiss_o3, ihno3_num, io3_num, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, num_species)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      implicit none

#     include "gmi_phys_constants.h"
#     include "gmi_emiss_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, i1, i2, ju1, j2, k1, k2, num_species
      logical :: pr_surf_emiss, do_ShipEmission
      logical :: pr_emiss_3d
      integer :: iisoprene_num, ico_num, ipropene_num, ino_num, io3_num, ihno3_num
      real*8  :: tdt
      real*8  :: mw   (num_species)
      real*8  :: mcor           (i1:i2, ju1:j2)
      real*8  :: surf_emiss_out (i1:i2, ju1:j2, num_species)
      real*8  :: surf_emiss_out2(i1:i2, ju1:j2, 6)
      real*8  :: emiss_3d_out   (i1:i2, ju1:j2, k1:k2, num_species)
      real*8  :: mass           (i1:i2, ju1:j2, k1:k2)
      real*8  :: emiss_isop     (i1:i2, ju1:j2)
      real*8  :: emiss_monot    (i1:i2, ju1:j2)
      real*8  :: emiss_nox      (i1:i2, ju1:j2)
      real*8  :: emiss_o3       (i1:i2, ju1:j2)
      real*8  :: emiss_hno3     (i1:i2, ju1:j2)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)

      integer :: ik

!     ----------------------
!     Variable declarations.
!     ----------------------

      real*8  :: emass(i1:i2, ju1:j2)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Add_Emiss_Harvard called by ', loc_proc
      end if


      if (iisoprene_num > 0) then

!       ----------------------------
!       Biogenic source of isoprene.
!       ----------------------------

        emass(:,:) =  &
     &    emiss_isop(:,:) * tdt

        concentration(iisoprene_num)%pArray3D(:,:,1) =  &
     &    concentration(iisoprene_num)%pArray3D(:,:,1) +  &
     &    ((emass(:,:) / mass(:,:,1)) *  &
     &     (MWTAIR / mw(iisoprene_num)))

       if (pr_surf_emiss)  &
     &   surf_emiss_out(:,:,iisoprene_num) =  &
     &     surf_emiss_out(:,:,iisoprene_num) +  &
     &     (emass(:,:) / mcor(:,:))


       if (pr_emiss_3d) then
          emiss_3d_out(:,:,1,iisoprene_num) =  &
     &      emiss_3d_out(:,:,1,iisoprene_num) +  &
     &       (emass(:,:) / mcor(:,:))

       end if
      end if

      if (ico_num > 0) then

!       --------------------------------------------------------------------
!       Biogenic source of CO from methanol oxidation, scaled from isoprene.
!       --------------------------------------------------------------------

        emass(:,:) =  &
     &    emiss_isop(:,:) * tdt *  &
     &    ICO_FAC_ISOP

        concentration(ico_num)%pArray3D(:,:,1) =  &
     &    concentration(ico_num)%pArray3D(:,:,1) +  &
     &    ((emass(:,:) / mass(:,:,1)) *  &
     &     (MWTAIR / mw(ico_num)))

       if (pr_surf_emiss)  &
     &   surf_emiss_out(:,:,ico_num) =  &
     &     surf_emiss_out(:,:,ico_num) +  &
     &     (emass(:,:) / mcor(:,:))
!!!!!!!!!!!!!!!!CO_methanol
       if (pr_surf_emiss)  &
     &   surf_emiss_out2(:,:,1) =  &
     &     surf_emiss_out2(:,:,1) +  &
     &     (emass(:,:) / mcor(:,:))
!!!!!!!!!!!!!!!!


       if (pr_emiss_3d) then
          emiss_3d_out(:,:,1,ico_num) =  &
     &      emiss_3d_out(:,:,1,ico_num) +  &
     &       (emass(:,:) / mcor(:,:))

       end if

!       ---------------------------------------------------------------------
!       Biogenic source of CO from monoterpene oxidation.
!       ---------------------------------------------------------------------

        emass(:,:) =  &
     &    emiss_monot(:,:) * tdt *  &
     &    ICO_FAC_MONOT

        concentration(ico_num)%pArray3D(:,:,1) =  &
     &    concentration(ico_num)%pArray3D(:,:,1) +  &
     &    ((emass(:,:) / mass(:,:,1)) *  &
     &     (MWTAIR / mw(ico_num)))

       if (pr_surf_emiss)  &
     &   surf_emiss_out(:,:,ico_num) =  &
     &     surf_emiss_out(:,:,ico_num) +  &
     &     (emass(:,:) / mcor(:,:))

!!!!!!!!!!!!!!!!CO_monoterpene
       if (pr_surf_emiss)  &
     &   surf_emiss_out2(:,:,2) =  &
     &     surf_emiss_out2(:,:,2) +  &
     &     (emass(:,:) / mcor(:,:))
!!!!!!!!!!!!!!!!

       if (pr_emiss_3d)  &
     &    emiss_3d_out(:,:,1,ico_num) =  &
     &      emiss_3d_out(:,:,1,ico_num) +  &
     &       (emass(:,:) / mcor(:,:))


      end if

!     ---------------------
      if (ipropene_num > 0) then
!     ---------------------

!       -------------------------------------------------
!       Biogenic source of propene, scaled from isoprene.
!       -------------------------------------------------

        emass(:,:) =  &
     &    emiss_isop(:,:) * tdt *  &
     &    BIOSCAL *  &
     &    ((ATOMSC_PER_MOLECISOP / ATOMSC_PER_MOLECPRPE) *  &
     &     (mw(ipropene_num)     / mw(iisoprene_num)))

        concentration(ipropene_num)%pArray3D(:,:,1) =  &
     &    concentration(ipropene_num)%pArray3D(:,:,1) +  &
     &    ((emass(:,:) / mass(:,:,1)) *  &
     &     (MWTAIR / mw(ipropene_num)))

       if (pr_surf_emiss)  &
     &   surf_emiss_out(:,:,ipropene_num) =  &
     &     surf_emiss_out(:,:,ipropene_num) +  &
     &     (emass(:,:) / mcor(:,:))

!!!!!!!!!!!!!!!!Biogenic_propene
       if (pr_surf_emiss)  &
     &   surf_emiss_out2(:,:,3) =  &
     &     surf_emiss_out2(:,:,3) +  &
     &     (emass(:,:) / mcor(:,:))
!!!!!!!!!!!!!!!!

       if (pr_emiss_3d) then
          emiss_3d_out(:,:,1,ipropene_num) =  &
     &      emiss_3d_out(:,:,1,ipropene_num) +  &
     &       (emass(:,:) / mcor(:,:))
       end if
!!!!!!!!!!!!!!!!!!

      end if


!     ----------------
      if (ino_num > 0) then
!     ----------------

!       -------------------
!       Soil source of NOx.
!       -------------------

        emass(:,:) =  &
     &    emiss_nox(:,:) * tdt

        concentration(ino_num)%pArray3D(:,:,1) =  &
     &    concentration(ino_num)%pArray3D(:,:,1) +  &
     &    ((emass(:,:) / mass(:,:,1)) *  &
     &     (MWTAIR / mw(ino_num)))

       if (pr_surf_emiss)  &
     &   surf_emiss_out(:,:,ino_num) =  &
     &     surf_emiss_out(:,:,ino_num) +  &
     &     (emass(:,:) / mcor(:,:))

!!!!!!!!!!!!!!!!Soil_NOx
       if (pr_surf_emiss)  &
     &   surf_emiss_out2(:,:,4) =  &
     &     surf_emiss_out2(:,:,4) +  &
     &     (emass(:,:) / mcor(:,:))
!!!!!!!!!!!!!!!!
       if (pr_emiss_3d) then
          emiss_3d_out(:,:,1,ino_num) =  &
     &      emiss_3d_out(:,:,1,ino_num) +  &
     &       (emass(:,:) / mcor(:,:))

       end if
!!!!!!!!!!!!!!!!!!

      end if

      if (pr_surf_emiss) then
        if (do_ShipEmission) then
            ! --------------
            ! Source of HNO3.
            ! --------------
          if (ihno3_num > 0) then
            emass(:,:) =  emiss_hno3(:,:) * tdt
            surf_emiss_out2(:,:,5) = surf_emiss_out2(:,:,5) +  &
     &                                  (emass(:,:) / mcor(:,:))
          end if
            ! -------------
            ! Source of O3.
            ! -------------
          if (io3_num > 0) then
            emass(:,:) =  emiss_o3(:,:) * tdt
            surf_emiss_out2(:,:,6) = surf_emiss_out2(:,:,6) +  &
     &                                  (emass(:,:) / mcor(:,:))
          end if
        else
          surf_emiss_out2(:,:,5) = 0.0d0
          surf_emiss_out2(:,:,6) = 0.0d0
        end if
      end if

      return

      end
