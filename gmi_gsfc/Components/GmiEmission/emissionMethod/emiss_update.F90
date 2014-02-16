
!=============================================================================
!
! $Id: emiss_update.F90,v 1.28 2013-08-28 20:54:52 jkouatch Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   emiss_update.F
!
! ROUTINES
!   Update_Emiss
!   Update_Semiss_Inchem
!
! HISTORY
!   - November 23, 2004 * Jules Kouatchou
!     Modified the routine Update_Emiss by adding the variable
!     "surf_emiss_out2" as argument of Update_Emiss and Add_Emiss_Harvard.
!   - December 8, 2005 * Bigyani Das
!     Added changes to update DMS, dust, sea salt emissions daily introducing
!     medt = 86400.0d0 and changing mdt to medt for the aerocom runs
!     when "do_aerocom" is true.
!=============================================================================
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update_Emiss
!
! DESCRIPTION
!   This routine updates const based on emissions.
!
! ARGUMENTS
!   lwi_flags     : array of flags that indicate land, water, or ice
!   latdeg        : latitude  (deg)
!   londeg        : longitude (deg)
!   mcor          : area of grid box (m^2)
!   emiss_isop    : isoprene    emissions (kg/s)
!   emiss_monot   : monoterpene emissions (kg/s)
!   emiss_nox     : NOx         emissions (kg/s)
!   radswg        : net downward shortwave radiation at ground (W/m^2)
!   surf_air_temp : surface air temperature  (degK)
!   surf_rough    : surface roughness   (m)
!   con_precip    : convective precipitation (mm/day)
!   tot_precip    : total precipitation (mm/day)
!   ustar         : ustar (m/s)
!   mass          : total mass of the atmosphere within each grid box (kg)
!   max_cloud     : maximum overlap cloud fraction for LW
!   ran_cloud     : random  overlap cloud fraction for LW
!   kel           : temperature (degK)
!   surf_emiss_out: surface emissions to be written to output (kg/m^2/time)
!   emiss_3d_out: surface emissions to be written to output (kg/m^2/time)
!   surf_emiss_out2: surface emissions to be written to output (kg/m^2/time)
!   const         : species concentration, known at zone centers (mixing ratio)
!   emissionArray : array of emissions      (kg/s)
!   emiss_dust_t  : array of dust emissions (kg/s)
!   emiss_dust    : array of dust emissions for 6 hours (kg/s)
!   emiss_aero    : array of aerosol emissions (kg/s)
!   pbl           : boundary layer height (m)
!   humidity      : specific humidity (g/kg)
!   pctm1         : surface pressure at t1 (mb)
!-----------------------------------------------------------------------------

      subroutine Update_Emiss  &
     &  (lwi_flags, cosSolarZenithAngle, latdeg, mcor, emiss_isop, emiss_monot,  &
     &   emiss_nox, do_ShipEmission, emiss_hno3, emiss_o3, ihno3_num, io3_num, &
     &   doMEGANemission, days_btw_m, &
     &   T_15_AVG, aefIsop, aefMbo, aefMonot, isoLai, isoLaiCurr, isoLaiPrev, &
     &   radswg, surf_air_temp, pardif, pardir, surf_rough, con_precip,  &
     &   tot_precip, ustar, mass, fracCloudCover, kel,  &
     &   surf_emiss_out, surf_emiss_out2, emiss_3d_out, &
     &   aerosolEmiss3D, aerosolSurfEmiss, aerosolSurfEmissMap, &
     &   concentration, emissionArray, emiss_dust_t, emiss_dust,  &
     &   emiss_aero_t, emiss_aero,pbl, gridBoxHeight, &
     &   index_soil, ncon_soil, soil_fert, soil_precip, soil_pulse, &
     &   ireg, iland, iuse, convert_isop, convert_monot, &
     &   coeff_isop, base_isop, base_monot, xlai, &
     &   IBOC, IBBC, INOC, IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, ISSLT4, &
     &   IFSO2, INSO2, INDMS, IDUST1, IDUST2, IDUST3, IDUST4, IDUST5, IAN, IMGAS, INO, &
     &   iisoprene_num, ino_num, ico_num, ipropene_num,   &
     &   pr_surf_emiss, pr_emiss_3d, pr_diag, loc_proc, met_opt, &
     &   emiss_opt, chem_opt, emiss_aero_opt, emiss_dust_opt, do_aerocom, &
     &   do_semiss_inchem, do_drydep, emiss_map, emiss_map_dust, emiss_map_aero, &
     &   ndust, nst_dust, nt_dust, naero, nymd, num_time_steps, mw, tdt, ndt, &
     &   emiss_timpyr, num_emiss,  isop_scale, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, &
     &   i1_gl, i2_gl, ju1_gl, j2_gl, ilong, num_species, soil_day_save)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiTimeControl_mod       , only : GmiSplitDateTime
      use GmiHeightSigmaLevel_mod  , only : CalcHeightSigmaLevel
      use CalcAerosolEmissDiagn_mod, only : calcAerosolEmissDiagn
      use m_set_NLANDHAR           , only : NLANDHAR

      implicit none

#     include "gmi_emiss_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in   ) :: IBOC, IBBC, INOC, IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, ISSLT4
      integer, intent(in   ) :: IFSO2, INSO2, INDMS, IDUST1, IDUST2, IDUST3, IDUST4, IDUST5, IAN, IMGAS, INO
      logical, intent(in   ) :: pr_surf_emiss, pr_emiss_3d, pr_diag
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in   ) :: i1_gl, i2_gl, ju1_gl, j2_gl, ilong
      integer, intent(in   ) :: num_species, emiss_timpyr, num_emiss
      integer, intent(in   ) :: ndust, nst_dust, nt_dust, naero
      integer, intent(in   ) :: nymd, num_time_steps, ndt
      integer, intent(in   ) :: iisoprene_num, ino_num, ico_num, ipropene_num, io3_num, ihno3_num
      real*8 , intent(in   ) :: mw(num_species)
      real*8 , intent(in   ) :: isop_scale(12)
      real*8 , intent(in   ) :: tdt
      integer, intent(in   ) :: met_opt, emiss_opt, chem_opt, emiss_aero_opt, emiss_dust_opt
      integer, intent(in   ) :: emiss_map(num_emiss), emiss_map_dust(num_species)
      integer, intent(in   ) :: emiss_map_aero(num_species)
      logical, intent(in   ) :: do_semiss_inchem, do_drydep, do_aerocom, do_ShipEmission
      integer, intent(in   ) :: lwi_flags(i1:i2, ju1:j2)
      real*8 , intent(in   ) :: cosSolarZenithAngle(i1:i2, ju1:j2)
      real*8 , intent(in   ) :: latdeg   (ju1_gl:j2_gl)
      real*8 , intent(in   ) :: mcor     (i1:i2, ju1:j2)
      real*8 , intent(  out) :: emiss_isop   (i1:i2, ju1:j2)
      real*8 , intent(  out) :: emiss_monot  (i1:i2, ju1:j2)
      real*8 , intent(  out) :: emiss_nox    (i1:i2, ju1:j2)
      logical, intent(in) :: doMEGANemission
      integer, intent(in) :: days_btw_m
      real*8 , intent(in) :: T_15_AVG            (i1:i2, ju1:j2)
      real*8 , intent(in) :: aefIsop             (i1:i2, ju1:j2)
      real*8 , intent(in) :: aefMbo              (i1:i2, ju1:j2)
      real*8 , intent(in) :: aefMonot            (i1:i2, ju1:j2)
      real*8 , intent(in) :: isoLai              (i1:i2, ju1:j2)
      real*8 , intent(in) :: isoLaiCurr          (i1:i2, ju1:j2)
      real*8 , intent(in) :: isoLaiPrev          (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: emiss_o3     (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: emiss_hno3     (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: radswg       (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: surf_air_temp(i1:i2, ju1:j2)
      real*8 , intent(in   ) :: pardif(i1:i2, ju1:j2)
      real*8 , intent(in   ) :: pardir(i1:i2, ju1:j2)
      real*8 , intent(in   ) :: surf_rough   (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: con_precip   (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: tot_precip   (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: ustar        (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: mass         (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in   ) :: kel          (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(inout) :: surf_emiss_out (i1:i2, ju1:j2, num_species)
      real*8 , intent(inout) :: surf_emiss_out2 (i1:i2, ju1:j2, 6)
      real*8 , intent(inout) :: emiss_3d_out (i1:i2, ju1:j2, k1:k2, num_species)
      real*8 , intent(inout) :: aerosolSurfEmiss   (i1:i2, ju1:j2, ndust+naero+5)
      real*8 , intent(inout) :: aerosolEmiss3D   (i1:i2, ju1:j2, k1:k2, 5)
      integer, intent(in   ) :: aerosolSurfEmissMap(ndust+naero+5)
      real*8 , intent(in   ) :: emiss_dust_t (i1:i2, ju1:j2, ndust, nst_dust:nst_dust+nt_dust-1)
      real*8 , intent(  out) :: emiss_dust   (i1:i2, ju1:j2, ndust)
      real*8 , intent(  out) :: emiss_aero   (i1:i2, ju1:j2, naero)
      real*8 , intent(in   ) :: emiss_aero_t (i1:i2, ju1:j2, naero, emiss_timpyr)
      real*8 , intent(in   ) :: pbl          (i1:i2, ju1:j2)
      real*8 , intent(in   ) :: gridBoxHeight (i1:i2, ju1:j2, k1:k2) ! for Llnl emission
      real*8 , intent(in   ) :: fracCloudCover (i1:i2, ju1:j2) ! cloud fraction
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)
      type (t_GmiArrayBundle), intent(in   ) :: emissionArray(num_emiss)
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

      integer :: idumyr, month, idumday
      integer :: num_new_emiss_dust
      real*8  :: medt

      real*8  :: tempk(i1:i2, ju1:j2)

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Update_Emiss called by ', loc_proc
      end if

!     -----------------------------------
!     Update dust emission every 6 hours.
!     -----------------------------------

      if (emiss_dust_opt == 1) then    ! GMI dust emissions

        if (do_aerocom) then
           medt = 86400.0d0
        else
           medt = 21600.0d0
        end if

        num_new_emiss_dust =  &
     &    (num_time_steps / (Nint (medt) / ndt)) + 1

        if (Mod (num_time_steps, (Nint (medt) / ndt)) == 0) then

          emiss_dust(i1:i2,ju1:j2,1:ndust) =  &
     &      emiss_dust_t(i1:i2,ju1:j2,1:ndust,num_new_emiss_dust)

        end if

      end if


      if (btest(emiss_opt,0) .or. btest(emiss_opt,1)) then

!       ===================
        call Add_Emiss_Llnl  &
!       ===================
     &    (do_aerocom,  &
     &     pr_surf_emiss, pr_emiss_3d, mcor, surf_emiss_out, emiss_3d_out,  &
     &     mass, concentration, emissionArray, emiss_dust,  emiss_aero_t, emiss_aero, pbl, gridBoxHeight, &
     &   IBOC, IBBC, INOC, IFOC, IFBC, ISSLT1, ISSLT2, ISSLT3, ISSLT4, &
     &   IFSO2, INSO2, INDMS, IDUST1, IDUST2, IDUST3, IDUST4, IDUST5, &
     &   pr_diag, loc_proc, &
     &   chem_opt, emiss_aero_opt, emiss_dust_opt, &
     &   do_semiss_inchem, emiss_map, emiss_map_dust, emiss_map_aero, &
     &   ndust, naero, nymd, mw, tdt, &
     &   emiss_timpyr, num_emiss, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species)
 
        ! Compute surface emission diagnostics for aerosols.

        if ((pr_surf_emiss) .or. (pr_emiss_3d)) then
           if ((emiss_aero_opt > 0) .or. (emiss_dust_opt > 0)) then
              call calcAerosolEmissDiagn(aerosolEmiss3D, aerosolSurfEmiss, &
                       aerosolSurfEmissMap, &
                       emissionArray, emiss_dust, emiss_aero, tdt, mcor, &
                       pr_surf_emiss, pr_emiss_3d, &
                       emiss_aero_opt, emiss_dust_opt, emiss_map, num_species, &
                       ndust, naero, i1, i2, ju1, j2, k1, k2, num_emiss)
           end if
        end if

        if (btest(emiss_opt,1)) then

          emiss_isop (:,:) = 0.0d0
          emiss_monot(:,:) = 0.0d0
          emiss_nox  (:,:) = 0.0d0

          if (met_opt == 3) then
            tempk(i1:i2,ju1:j2) = surf_air_temp(i1:i2,ju1:j2)
          else
            tempk(i1:i2,ju1:j2) = kel(i1:i2,ju1:j2,1)
          end if

!         =========================
          call Update_Emiss_Harvard  &
!         =========================
     &      (iisoprene_num, ino_num, lwi_flags, tdt, mw, cosSolarZenithAngle, &
     &       latdeg, nymd, mcor, &
     &       aefIsop, aefMbo, aefMonot, isoLai, isoLaiCurr, isoLaiPrev, &
     &       T_15_AVG, tempk, pardif, pardir, radswg, surf_rough,  &
     &       con_precip, tot_precip, ustar, fracCloudCover,  &
     &       emiss_isop, emiss_monot, emiss_nox, &
     &       index_soil, ncon_soil, soil_fert, soil_precip, soil_pulse, &
     &       ireg, iland, iuse, convert_isop, convert_monot, &
     &       coeff_isop, base_isop, base_monot, xlai, &
     &       doMEGANemission, days_btw_m, soil_day_save, pr_diag, loc_proc, &
     &       i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl, ilong, num_species)

!         ------------------------------------------------
!         Scale the isoprene emissions on a monthly basis.
!         Default is no scaling (isop_scale = 1.0d0).
!         ------------------------------------------------

          call GmiSplitDateTime (nymd, idumyr, month, idumday)

          emiss_isop(:,:) = emiss_isop(:,:) * isop_scale(month)

          if (.not. do_semiss_inchem) then
!           ======================
            call Add_Emiss_Harvard  &
!           ======================
     &        (pr_surf_emiss, pr_emiss_3d, iisoprene_num, ico_num,  &
     &         ipropene_num,  &
     &         ino_num, tdt, mw, mcor, surf_emiss_out, surf_emiss_out2,  &
     &         emiss_3d_out, mass, concentration, emiss_isop, emiss_monot, emiss_nox, &
     &   do_ShipEmission, emiss_hno3, emiss_o3, ihno3_num, io3_num, &   
     &         pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, num_species)

          end if


        end if

      end if


      return

      end
