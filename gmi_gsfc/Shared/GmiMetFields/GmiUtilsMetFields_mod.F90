!-----------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!-----------------------------------------------------------------------------
!BOP
!
! !IMODULE: GmiUtilsMetFields_mod
!
! !INTERFACE:
!
      module GmiUtilsMetFields_mod
!
! !USES:
      use ESMF_Mod, only : ESMF_MAXSTR
      use GmiPrintError_mod, only : GmiPrintError
      use GmiWrapMaster_mod, only : wrapMaster_2d, wrapMaster_3du, wrapMaster_3dv
      use GmiPressureFixer_mod, only : Calc_Delpm
      use m_netcdf_io_read , only : Ncrd_1d, Ncrd_3d, Ncrd_4d
      use m_netcdf_io_get_dimlen, only : Ncget_Unlim_Dimlen, Ncget_Dimlen
      use m_netcdf_io_open      , only : Ncop_Rd
      use m_netcdf_io_close     , only : Nccl
      use m_netcdf_io_checks    , only : Ncdoes_Var_Exist
      use GmiCalcDivergence_mod, only : calcDivergence, Deter_Jrange
      use GmiMassFluxes_mod    , only : calcVerticalMassFlux
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: Read_Met1, Set_Met1_Simple
      public  :: Read_Met2, Calc_Var_Adv_Tstp, Calc_Courant_Met
      public  :: Get_Next_Met_File, Check_Met1_Range, Check_Vert_Courant
      public  :: Calc_Humidity, Convert_Pbl
!
#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
!
!EOP
!-----------------------------------------------------------------------------
      contains
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Met1_Simple
!
! !INTERFACE:
!
! !DESCRIPTION:
!   This routine sets up a simple set of met1 data (ps, u, v, & kel).
!
!   Note that velocity fields are defined on p. 37 of Williamson & Rasch's
!   1988 paper "Two-dimensional Semi-lagrangian Transport with Space
!   Preserving Interpolation".
!
! ARGUMENTS
!   psx     : surface pressure at beginning of met time period (mb)
!   del_psx : surface pressure at end       of met time period (mb)
!   uux     : current horizontal velocity in zonal (longitude) direction;
!             known at edges in longitude direction and centers in latitude
!             direction (m/s)
!   del_uux : horizontal velocity in zonal (longitude) direction at end of
!             met time period; known at edges in longitude direction and
!             centers in latitude direction  (m/s)
!   vvx     : current horizontal velocity in meridional (latitude) direction;
!             known at edges in latitude direction and centers in longitude
!             direction (m/s)
!   del_vvx : horizontal velocity in meridional (latitude) direction at end
!             of met time period; known at edges in latitude direction and
!             centers in longitude direction (m/s)
!   kel     : current temperature (degK)
!   del_kel : temperature at end of met time period (degK)
!   coscen  : cosine of latitude of zone centers = cos(dlatr)
!   dlatr   : latitude of zone center in latitude direction  (rad)
!   do_wind_pole : make winds blow over Pole?
!
!-----------------------------------------------------------------------------
!
      subroutine Set_Met1_Simple (psx, del_psx, uux, del_uux, vvx, del_vvx, &
     &              kel, del_kel, coscen, dlatr, do_wind_pole, pr_diag, &
     &              procID, i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, ilo_gl, &
     &              ihi_gl, julo_gl, jvlo_gl, jhi_gl, k1, k2)
!
      implicit none
!
! !INTPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, k1, k2
      integer, intent(in) :: ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl
      real*8 , intent(in) :: coscen (ju1_gl:j2_gl)
      real*8 , intent(in) :: dlatr  (ju1_gl:j2_gl)
      logical, intent(in) :: do_wind_pole
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: psx    (ilo_gl:ihi_gl, julo_gl:jhi_gl)
      real*8 , intent(out) :: del_psx(ilo_gl:ihi_gl, julo_gl:jhi_gl)
      real*8 , intent(out) :: uux    (ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)
      real*8 , intent(out) :: del_uux(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)
      real*8 , intent(out) :: vvx    (ilo_gl:ihi_gl, jvlo_gl:jhi_gl, k1:k2)
      real*8 , intent(out) :: del_vvx(ilo_gl:ihi_gl, jvlo_gl:jhi_gl, k1:k2)
      real*8 , intent(out) :: kel    (ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)
      real*8 , intent(out) :: del_kel(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)
!
! !LOCAL VARIABLES:
      integer :: il, ij, ik
      real*8  :: angle
      real*8  :: d2
      real*8  :: f0
      real*8  :: ri2
      real*8  :: uw1
      real*8  :: dl (i1_gl:i2_gl)
      real*8  :: dl2(i1_gl:i2_gl)
      real*8  :: uw2(i1_gl:i2_gl, ju1_gl:j2_gl)
      real*8  :: vw2(i1_gl:i2_gl, jv1_gl:j2_gl)
      character(len=ESMF_MAXSTR) :: IAm
!
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = " Set_Met1_Simple"
!
      if (pr_diag) Write (6,*) Iam, 'Set_Met1_Simple called by ', procID
!
      ri2 = i2_gl
      d2  = 360.0d0 / ri2
      f0  = -d2 / 2.0d0
!
      do il = i1_gl, i2_gl
        dl(il)  =       (il - 1) * d2  * RADPDEG
        dl2(il) = (f0 + (il - 1) * d2) * RADPDEG
      end do
!
      if (do_wind_pole) then
        angle = 0.5d0 * GMI_PI
      else
        angle = 0.0d0
      end if
!
      uw1 = 2.0d0 * GMI_PI * RADEAR  / (12.0d0 * SECPDY)
!
      do ij = ju1_gl, j2_gl
!
        uw2(i1_gl:i2_gl,ij) =  &
     &    uw1 *  &
     &    (Cos (angle) * coscen(ij) +  &
     &     Sin (angle) * Sin (dlatr(ij)) * Cos (dl2(i1_gl:i2_gl)))
!
      end do
!
      do ij = jv1_gl, j2_gl
        vw2(i1_gl:i2_gl,ij) = -uw1 * Sin (angle) * Sin (dl(i1_gl:i2_gl))
      end do
!
      psx    (i1_gl:i2_gl,ju1_gl:j2_gl) = 1000.0d0
      del_psx(i1_gl:i2_gl,ju1_gl:j2_gl) = 1000.0d0
!
      do ik = k1, k2
        uux(i1_gl:i2_gl,ju1_gl:j2_gl,ik) = uw2(i1_gl:i2_gl,ju1_gl:j2_gl)
        vvx(i1_gl:i2_gl,jv1_gl:j2_gl,ik) = vw2(i1_gl:i2_gl,jv1_gl:j2_gl)
      end do
!
      where (Abs (uux(i1_gl:i2_gl,ju1_gl:j2_gl,:)) < 0.000001d0)
        uux(i1_gl:i2_gl,ju1_gl:j2_gl,:) = 0.0d0
      end where
!
      where (Abs (vvx(i1_gl:i2_gl,jv1_gl:j2_gl,:)) < 0.000001d0)
        vvx(i1_gl:i2_gl,jv1_gl:j2_gl,:) = 0.0d0
      end where
!
      del_uux(i1_gl:i2_gl,ju1_gl:j2_gl,:) = uux(i1_gl:i2_gl,ju1_gl:j2_gl,:)
      del_vvx(i1_gl:i2_gl,jv1_gl:j2_gl,:) = vvx(i1_gl:i2_gl,jv1_gl:j2_gl,:)
!
      kel    (i1_gl:i2_gl,ju1_gl:j2_gl,:) = 273.0d0
      del_kel(i1_gl:i2_gl,ju1_gl:j2_gl,:) = 273.0d0
!
      return
!
      end subroutine Set_Met1_Simple
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read_Met1
!
! !INTERFACE:
!
! !DESCRIPTION:
!   This routine reads in a set of met1 data (ps, u, v, & kel).
!
! ARGUMENTS
!   ncid_met1 : NetCDF file id to read met1 data from
!   m1rnum_in : next NetCDF met1 record to read
!   ps  : array of pressure    data (mb)
!   u   : array of u wind      data (m/s)
!   v   : array of v wind      data (m/s)
!   kel : array of temperature data (degK)
!
!-----------------------------------------------------------------------------
!
      subroutine Read_Met1 (metdata_name_model, &
     &               ncid_met1, m1rnum_in, ps, uux, vvx, kel, &
     &               pr_diag, procID, i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, &
     &               k1, k2, ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      character (len=*), intent(in) :: metdata_name_model
      integer, intent(in) :: ncid_met1
      integer, intent(in) :: m1rnum_in
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, k1, k2
      integer, intent(in) :: ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl
      real*8  :: ps (ilo_gl:ihi_gl, julo_gl:jhi_gl)
      real*8  :: uux(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)
      real*8  :: vvx(ilo_gl:ihi_gl, jvlo_gl:jhi_gl, k1:k2)
      real*8  :: kel(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)
!
!
!     -----------------------
!     Parameter declarations.
!     -----------------------
!
      character (len=MAX_LENGTH_VAR_NAME) :: KEL_VNAM
      character (len=MAX_LENGTH_VAR_NAME) :: PS_VNAM
      character (len=MAX_LENGTH_VAR_NAME) :: U_VNAM
      character (len=MAX_LENGTH_VAR_NAME) :: V_VNAM
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
      real*8 , allocatable :: uux_temp(:,:,:)
      real*8 , allocatable :: vvx_temp(:,:,:)
      real*8 , allocatable :: kel_temp(:,:,:)
!
      integer :: isize, jusize, jvsize, ksize
!
      integer :: cnt3d (3), cnt4d (4)
      integer :: strt3d(3), strt4d(4)
      character(len=ESMF_MAXSTR) :: IAm
!
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = "Read_Met1"
!
      if (pr_diag) Write (6,*) Iam, ' called by ', procID
!
      if (metdata_name_model(1:5) .eq. "GEOS5") then
         KEL_VNAM = 'T'
         PS_VNAM  = 'PS'
         U_VNAM   = 'U'
         V_VNAM   = 'V'
      else
         KEL_VNAM = 'kel'
         PS_VNAM  = 'psf'
         U_VNAM   = 'u'
         V_VNAM   = 'v'
      end if
!
      isize  = i2_gl - i1_gl + 1
      jusize = j2_gl - ju1_gl + 1
!      jvsize = j2_gl - jv1_gl + 1
      ksize  = k2 - k1 + 1
!
      strt3d(1) = 1
      strt3d(2) = 1
      strt3d(3) = m1rnum_in
!
      cnt3d (:) = (/ isize, jusize, 1 /)
!
      call Ncrd_3d (ps(i1_gl:i2_gl,ju1_gl:j2_gl), ncid_met1, PS_VNAM, strt3d, cnt3d)
!
      strt4d(1) = 1
      strt4d(2) = 1
      strt4d(3) = k1
      strt4d(4) = m1rnum_in
!
      cnt4d (:) = (/ isize, jusize, ksize, 1 /)
!
      call Ncrd_4d (uux(i1_gl:i2_gl,ju1_gl:j2_gl,:), ncid_met1, U_VNAM,  &
     &              strt4d, cnt4d)
      call Ncrd_4d (kel(i1_gl:i2_gl,ju1_gl:j2_gl,:), ncid_met1, KEL_VNAM,  &
     &              strt4d, cnt4d)
!
      strt4d(2) = 1
      cnt4d (:) = (/ isize, jusize, ksize, 1 /)
!
      call Ncrd_4d (vvx(i1_gl:i2_gl,ju1_gl:j2_gl,:), ncid_met1, V_VNAM, strt4d, cnt4d)
!
      if (metdata_name_model(1:5) .eq. "GEOS5") then
!
         ps = ps / 100.0d0
!
         ! GEOS5 level 72 is surface - invert the k level
!
         allocate(uux_temp(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2))
         allocate(vvx_temp(ilo_gl:ihi_gl, jvlo_gl:jhi_gl, k1:k2))
         allocate(kel_temp(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2))
!
         uux_temp(:,:,k1:k2) = uux(:,:,k2:k1:-1)
         vvx_temp(:,:,k1:k2) = vvx(:,:,k2:k1:-1)
         kel_temp(:,:,k1:k2) = kel(:,:,k2:k1:-1)
!
         uux(:,:,:) = uux_temp(:,:,:)
         vvx(:,:,:) = vvx_temp(:,:,:)
         kel(:,:,:) = kel_temp(:,:,:)
!
         deallocate(uux_temp)
         deallocate(vvx_temp)
         deallocate(kel_temp)
!
      end if
!
      return
!
      end subroutine Read_Met1
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read_Met2
!
! !INTERFACE:
!
! !DESCRIPTION:
!   This routine reads in a set of met2 data (everything except ps, u, v,
!   & kel).
!
! ARGUMENTS
!   metdata_name_org   : first  part of metdata_name, e.g., "NCAR"
!   metdata_name_model : second part of metdata_name, e.g., "MATCH"
!   do_wetdep     : do wet deposition?
!   ncid_met2     : NetCDF file id to read met2 data from
!   m2rnum_in     : next NetCDF met2 record to read
!   chem_opt      : chemistry      option
!   convec_opt    : convection     option
!   emiss_in_opt  : emissions      option
!   sfalbedo_opt  : surface albedo option
!   uvalbedo_opt  : UV albedo      option
!   lwi_flags     : Land, water, ice flags (1: water, 2: land, 3: ice)
!   con_precip    : convective precipitation (mm/day)
!   grnd_temp     : ground temperature (degK)
!   pbl           : planetary boundary layer depth (m)
!   radswg        : net downward shortwave radiation at ground (W/m^2)
!   saldif        : surface albedo for diffuse light (nr IR)  (frac. 0-1)
!   saldir        : surface albedo for direct  light (nr IR)  (frac. 0-1)
!   sasdif        : surface albedo for diffuse light (uv/vis) (frac. 0-1)
!   sasdir        : surface albedo for direct  light (uv/vis) (frac. 0-1)
!   surf_air_temp : surface air temperature  (degK)
!   surf_alb_uv   : bulk surface albedo (fraction 0-1)
!   surf_rough    : surface roughness   (m)
!   tot_precip    : total precipitation (mm/day)
!   ustar         : ustar (m/s)
!   md            : convective mass flux in downdraft (kg/m^2*s)
!   cmf           : convective mass flux              (kg/m^2*s)
!   dtrn          : detrainment rate  (DAO:kg/m^2*s, NCAR:s^-1)
!   ed            : entrainment into convective downdraft (s^-1)
!   eu            : entrainment into convective updraft   (s^-1)
!   humidity      : specific humidity (g/kg)
!   kzz           : array of vertical diffusion coefficients (m^2/s)
!   max_cloud     : maximum overlap cloud fraction for LW
!   moistq        : moisture changes due to wet processes (g/kg/day)
!   rain          : rainfall across cell edges (mm/day)
!   ran_cloud     : random  overlap cloud fraction for LW
!   tau_cloud     : cloud optical thickness in visible
!   aconv         : convective cloud fraction (from CCM3)
!   astrat        : stratiform cloud fraction (from CCM3)
!   clwc          : cloud liq water, grid box avg (g/m^3) (from CCM3)
!   relhum        : relative humidity (from CCM3)
!   rhclear       : relative humidity of clear region (from CCM3)
!   cldmas        : convective mass flux in updraft (Pa/s) (from CCM3)
!   omega         : vertical pressure velocity      (Pa/s) (from CCM3)
!   rain_zm       : Z-M Rain Production (kg/kg/s) (from GEOS4)
!   rain_hk       : Hack Rain Production (kg/kg/s) (from GEOS4)
!   rain_ls       : Laarge-Scale Rain Production (kg/kg/s) (from GEOS4)
!   zmdu          : Z-M convective detrainment flux in updraft   (Pa/s)
!   zmeu          : Z-M convective entrainment flux in updraft   (Pa/s)
!   zmed          : Z-M convective entrainment flux in downdraft (Pa/s)
!   zmmd          : Z-M convective mass flux in downdraft        (Pa/s)
!   zmmu          : Z-M convective mass flux in updraft          (Pa/s)
!   hkdu          : Hack convective detrainment flux in updraft   (Pa/s)
!   hkeu          : Hack convective entrainment flux in updraft   (Pa/s)
!   hkmu          : Hack convective mass flux in updraft          (Pa/s)
!
!   u10m          : 10 meter Zonal Wind      (m/s)
!   v10m          : 10 meter Meridional Wind (m/s)
!   gwet          : Root zone soil wetness   (fraction)
!   pardif        : Diffuse photosynthetically active radiation (0.35-0.70 um)
!   pardir        : Direct  photosynthetically active radiation (0.35-0.70 um)
!
!-----------------------------------------------------------------------------
!
      subroutine Read_Met2  &
     &  (mcor, mass, metdata_name_org, metdata_name_model, do_wetdep, ncid_met2,  &
     &   m2rnum_in, chem_opt, convec_opt, emiss_in_opt,  &
     &   gwet_opt, sfalbedo_opt, uvalbedo_opt, &
     &   lwi_flags, con_precip, grnd_temp, pbl, radswg, saldif, saldir, &
     &   sasdif, sasdir, surf_air_temp, surf_alb_uv, surf_rough, tot_precip, &
     &   ustar, md, cmf, dtrn, ed, eu, humidity, kzz, max_cloud, moistq, rain, &
     &   ran_cloud, tau_cloud, aconv, astrat, clwc, relhum, rhclear, cldmas, &
     &   omega, rain_zm, rain_hk, rain_ls, taucli, tauclw, &
     &   zmdu, zmeu, zmed, zmmd, zmmu, hkdu, &
     &   hkeu, hkmu, u10m, v10m, gwet, pardif, pardir, pr_diag, procID, &
     &   i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, j2_gl)
!
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, i1_gl, ju1_gl, j2_gl
      logical, intent(in) :: do_wetdep
      integer :: kSave
      integer :: ncid_met2
      integer :: m2rnum_in
      integer :: chem_opt
      integer :: convec_opt
      integer :: gwet_opt
      integer :: emiss_in_opt
      integer :: sfalbedo_opt
      integer :: uvalbedo_opt
      real*8 , intent(in) :: mcor (i1:i2,ju1:j2)       ! area of grid box (m^2)
      real*8 , intent(in) :: mass (i1:i2,ju1:j2,k1:k2)
      integer :: lwi_flags    (i1:i2, ju1:j2)
      real*8  :: con_precip   (i1:i2, ju1:j2)
      real*8  :: grnd_temp    (i1:i2, ju1:j2)
      real*8  :: pbl          (i1:i2, ju1:j2)
      real*8  :: radswg       (i1:i2, ju1:j2)
      real*8  :: saldif       (i1:i2, ju1:j2)
      real*8  :: saldir       (i1:i2, ju1:j2)
      real*8  :: sasdif       (i1:i2, ju1:j2)
      real*8  :: sasdir       (i1:i2, ju1:j2)
      real*8  :: surf_air_temp(i1:i2, ju1:j2)
      real*8  :: surf_alb_uv  (i1:i2, ju1:j2)
      real*8  :: surf_rough   (i1:i2, ju1:j2)
      real*8  :: tot_precip   (i1:i2, ju1:j2)
      real*8  :: ustar        (i1:i2, ju1:j2)
      real*8  :: md           (i1:i2, ju1:j2, k1:k2)
      real*8  :: cmf          (i1:i2, ju1:j2, k1:k2)
      real*8  :: dtrn         (i1:i2, ju1:j2, k1:k2)
      real*8  :: ed           (i1:i2, ju1:j2, k1:k2)
      real*8  :: eu           (i1:i2, ju1:j2, k1:k2)
!... GEOS4 cloud flux fields
      real*8  :: zmdu         (i1:i2, ju1:j2, k1:k2)
      real*8  :: zmeu         (i1:i2, ju1:j2, k1:k2)
      real*8  :: zmed         (i1:i2, ju1:j2, k1:k2)
      real*8  :: zmmd         (i1:i2, ju1:j2, k1:k2)
      real*8  :: zmmu         (i1:i2, ju1:j2, k1:k2)
      real*8  :: hkdu         (i1:i2, ju1:j2, k1:k2)
      real*8  :: hkeu         (i1:i2, ju1:j2, k1:k2)
      real*8  :: hkmu         (i1:i2, ju1:j2, k1:k2)
!
      real*8  :: humidity     (i1:i2, ju1:j2, k1:k2)
      real*8  :: kzz          (i1:i2, ju1:j2, k1:k2)
      real*8  :: max_cloud    (i1:i2, ju1:j2, k1:k2)
      real*8  :: moistq       (i1:i2, ju1:j2, k1:k2)
      real*8  :: rain         (i1:i2, ju1:j2, k1:k2)
!... GEOS4 rain fields
      real*8  :: rain_zm      (i1:i2, ju1:j2, k1:k2)
      real*8  :: rain_hk      (i1:i2, ju1:j2, k1:k2)
      real*8  :: rain_ls      (i1:i2, ju1:j2, k1:k2)
      real*8  :: taucli       (i1:i2, ju1:j2, k1:k2)
      real*8  :: tauclw       (i1:i2, ju1:j2, k1:k2)
      real*8  :: pficu        (i1:i2, ju1:j2, k1:k2+1)
      real*8  :: pflcu        (i1:i2, ju1:j2, k1:k2+1)
      real*8  :: pfilsan      (i1:i2, ju1:j2, k1:k2+1)
      real*8  :: pfllsan      (i1:i2, ju1:j2, k1:k2+1)
!
      real*8  :: ran_cloud    (i1:i2, ju1:j2, k1:k2)
      real*8  :: tau_cloud    (i1:i2, ju1:j2, k1:k2)
      real*8  :: aconv        (i1:i2, ju1:j2, k1:k2)
      real*8  :: astrat       (i1:i2, ju1:j2, k1:k2)
      real*8  :: clwc         (i1:i2, ju1:j2, k1:k2)
      real*8  :: relhum       (i1:i2, ju1:j2, k1:k2)
      real*8  :: rhclear      (i1:i2, ju1:j2, k1:k2)
      real*8  :: cldmas       (i1:i2, ju1:j2, k1:k2)
      real*8  :: omega        (i1:i2, ju1:j2, k1:k2)
!... GEOS4
      real*8  :: u10m         (i1:i2, ju1:j2)
      real*8  :: v10m         (i1:i2, ju1:j2)
      real*8  :: gwet         (i1:i2, ju1:j2)
      real*8  :: pardif       (i1:i2, ju1:j2)
      real*8  :: pardir       (i1:i2, ju1:j2)
!... GEOS5
      real*8  :: delp         (i1:i2, ju1:j2, k1:k2)
!... EC_Oslo
      real*8  :: ciwc         (i1:i2, ju1:j2, k1:k2)
      real*8  :: psfc         (i1:i2, ju1:j2)
      real*8  :: ec_ai        (k1:k2+1)
      real*8  :: ec_bi        (k1:k2+1)
!
!.sds... used for EC only
!.sds... value for each layer, but 37 and 40 level EC fields have same bottom levels,
!.sds...   split top 3 levels in half, so I can have one array
!.sds...   Calculate tau for liquid water clouds. Use marine stratus for levels 1-13,
!.sds...       marine cumulus for levels 20-40, interpolate between.
      real*8, dimension(40), save :: REFF = ( (/  &
     &   9.6010000d0, 9.6010000d0, 9.6010000d0, 9.6010000d0, 9.6010000d0 &
     &  ,9.6010000d0, 9.6010000d0, 9.6010000d0, 9.6010000d0, 9.6010000d0 &
     & ,10.0408571d0,10.4807143d0,10.9205714d0,11.3604286d0,11.8002857d0 &
     & ,12.6800000d0,12.6800000d0,12.6800000d0,12.6800000d0,12.6800000d0 &
     & ,12.6800000d0,12.6800000d0,12.6800000d0,12.6800000d0,12.6800000d0 &
     & ,12.6800000d0,12.6800000d0,12.6800000d0,12.6800000d0,12.6800000d0 &
     & ,12.6800000d0,12.6800000d0,12.6800000d0,12.6800000d0,12.6800000d0 &
     & ,12.6800000d0,12.6800000d0,12.6800000d0,12.6800000d0,12.6800000d0 /) )
      real*8, dimension(40), save :: QEXT = (  (/ &
     &   2.0994d0, 2.0994d0, 2.0994d0, 2.0994d0, 2.0994d0, 2.0994d0 &
     &  ,2.0994d0, 2.0994d0, 2.0994d0, 2.0994d0, 2.0994d0, 2.0994d0 &
     &  ,2.0994d0, 2.0810d0, 2.0810d0, 2.0810d0, 2.0810d0, 2.0810d0 &
     &  ,2.0810d0, 2.0810d0, 2.0810d0, 2.0810d0, 2.0810d0, 2.0810d0 &
     &  ,2.0810d0, 2.0810d0, 2.0810d0, 2.0810d0, 2.0810d0, 2.0810d0 &
     &  ,2.0810d0, 2.0810d0, 2.0810d0, 2.0810d0, 2.0810d0, 2.0810d0 &
     &  ,2.0810d0, 2.0810d0, 2.0810d0, 2.0810d0 /) )
!
!     -----------------------
!     Parameter declarations.
!     -----------------------
!
!...ok
      character (len=MAX_LENGTH_VAR_NAME) :: CON_PRECIP_VNAM    = 'precon'
      character (len=MAX_LENGTH_VAR_NAME) :: GRND_TEMP_VNAM     = 'tg'
      character (len=MAX_LENGTH_VAR_NAME) :: LWI_VNAM           = 'lwi'
      character (len=MAX_LENGTH_VAR_NAME) :: PBL_VNAM           = 'pbl'
      character (len=MAX_LENGTH_VAR_NAME) :: RADSWG_VNAM        = 'radswg'
!... in GEOS4, not in DAO
      character (len=MAX_LENGTH_VAR_NAME) :: SF_ALB_LDIF_VNAM   = 'saldif'
      character (len=MAX_LENGTH_VAR_NAME) :: SF_ALB_LDIR_VNAM   = 'saldir'
      character (len=MAX_LENGTH_VAR_NAME) :: SF_ALB_SDIF_VNAM   = 'sasdif'
      character (len=MAX_LENGTH_VAR_NAME) :: SF_ALB_SDIR_VNAM   = 'sasdir'
!... in GEOS4, not in DAO
      character (len=MAX_LENGTH_VAR_NAME), parameter :: VERT_PRSE_VEL_VNAM = 'omega'
      character (len=MAX_LENGTH_VAR_NAME) :: SNOW_VNAM          = 'snowh'
!...ok
      character (len=MAX_LENGTH_VAR_NAME), parameter :: SURF_AIR_GEOS_VNAM = 't10m'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: SURF_AIR_GEO4_VNAM = 't2m'
      character (len=MAX_LENGTH_VAR_NAME) :: SURF_AIR_TEMP_VNAM = 'ts'
      character (len=MAX_LENGTH_VAR_NAME) :: SURF_ROUGH_VNAM    = 'z0'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: TOT_PRECIP_VNAM    = 'preacc'
      character (len=MAX_LENGTH_VAR_NAME) :: USTAR_VNAM         = 'ustar'
      character (len=MAX_LENGTH_VAR_NAME) :: HUMIDITY_VNAM      = 'sphu'
      character (len=MAX_LENGTH_VAR_NAME) :: KZZ_VNAM           = 'kzz'
      character (len=MAX_LENGTH_VAR_NAME) :: CMF_VNAM           = 'cldmas'
!... new units
      character (len=MAX_LENGTH_VAR_NAME) :: DTRN_DAO_VNAM      = 'dtrain'
!... GEOS4 names for convection properties, including rain
      character (len=MAX_LENGTH_VAR_NAME), parameter :: DTRNUP_ZM_GEO4     = 'zmdu'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: ETRNUP_ZM_GEO4     = 'zmeu'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: ETRNDN_ZM_GEO4     = 'zmed'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: MFLXDN_ZM_GEO4     = 'zmmd'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: MFLXUP_ZM_GEO4     = 'zmmu'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: DTRNUP_HK_GEO4     = 'cmfdtr'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: ETRNUP_HK_GEO4     = 'cmfetr'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: MFLXUP_HK_GEO4     = 'cmfmc2'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: RAIN_ZM_GEO4       = 'zmdqr'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: RAIN_HK_GEO4       = 'cmfdqr2'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: RAIN_LS_GEO4       = 'dqrl'
!... non-GEOS4 names for convection properties, including rain
      character (len=MAX_LENGTH_VAR_NAME) :: DTRN_NCAR_VNAM     = 'du'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: ED_VNAM            = 'ed'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: EU_VNAM            = 'eu'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: MD_VNAM            = 'md'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: RAIN_VNAM          = 'rain'
!... cloud optical depths
      character (len=MAX_LENGTH_VAR_NAME), parameter :: TAU_CLOUD_NCAR_VNAM= 'CLDVISOT'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: TAU_CLOUD_GISS_VNAM= 'cldod'
      character (len=MAX_LENGTH_VAR_NAME) :: TAU_CLDIC_GEO_VNAM= 'taucli'
      character (len=MAX_LENGTH_VAR_NAME) :: TAU_CLDWA_GEO_VNAM= 'tauclw'
!... moistq became dqcond in GEOS4 - not saved in 5 year gcm run
      character (len=MAX_LENGTH_VAR_NAME), parameter :: MOISTQ_GEOS4_VNAM  = 'dqcond'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: MOISTQ_VNAM        = 'moistq'
!...ok
      character (len=MAX_LENGTH_VAR_NAME), parameter :: CONV_MASS_UPD_VNAM = 'cldmas'
!... MAX and RAN seem to be used to calculate cloud fraction
      character (len=MAX_LENGTH_VAR_NAME) :: CLOUD_FRAC_VNAM    = 'cloud'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: MAX_CLOUD_VNAM     = 'clmolw'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: RAN_CLOUD_VNAM     = 'clrolw'
!.. no equiv right now... look
      character (len=MAX_LENGTH_VAR_NAME), parameter :: SURF_ALB_UV_VNAM   = 'albd'
!.. no equiv - don't seem to be used
      character (len=MAX_LENGTH_VAR_NAME), parameter :: CLD_LIQ_WATER_VNAM = 'clwc'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: CONV_CLD_FRAC_VNAM = 'aconv'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: REL_HUM_VNAM       = 'relhum'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: REL_HUM_CLEAR_VNAM = 'rhclear'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: STRT_CLD_FRAC_VNAM = 'astrat'
!... GEOS4 names
      character (len=MAX_LENGTH_VAR_NAME) :: U10m_VNAM   = 'u10m'
      character (len=MAX_LENGTH_VAR_NAME) :: V10m_VNAM   = 'v10m'
      character (len=MAX_LENGTH_VAR_NAME) :: gwet_vnam1  = 'gwet1'
      character (len=MAX_LENGTH_VAR_NAME) :: gwet_vnam2  = 'gwettop'
      character (len=MAX_LENGTH_VAR_NAME) :: PARDIF_VNAM = 'pardif'
      character (len=MAX_LENGTH_VAR_NAME) :: PARDIR_VNAM = 'pardir'
!... GEOS5 names
      character (len=MAX_LENGTH_VAR_NAME), parameter ::  SURF_AIR_GEO5_VNAM = 'T2M'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: PRECANV_VNAM = 'PRECANV'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: PRECSNO_VNAM = 'PRECSNO'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: PRECLSC_VNAM = 'PRECLSC'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: DQLDTMST_VNAM = 'DQLDTMST'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: DQIDTMST_VNAM = 'DQIDTMST'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: DQVDTMST_VNAM = 'DQVDTMST'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: RAIN_ZM_GEO5 = 'DQRCON'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: RAIN_LS_GEO5 = 'DQRLSC'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: OPTDEPTH_VNAM = 'OPTDEPTH'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: TAUCLI_VNAM = 'TAUCLI'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: TAUCLW_VNAM = 'TAUCLW'
!... GEOS5 MERRA 300 names
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: PFICU_VNAM = 'PFICU'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: PFLCU_VNAM = 'PFLCU'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: PFILSAN_VNAM = 'PFILSAN'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: PFLLSAN_VNAM = 'PFLLSAN'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: DELP_VNAM = 'DELP'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: FROCEAN_VNAM       = 'FROCEAN'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: FRLAKE_VNAM       = 'FRLAKE'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: FRSEAICE_VNAM       = 'FRSEAICE'
      character (len=MAX_LENGTH_VAR_NAME), parameter  :: FRLANDICE_VNAM       = 'FRLANDICE'
!
!... EC_Oslo names still unused
!.EC        float SLH(rec, jju, ii) ;
!.EC                SLH:long_name = "Surface Latent Heat Flux" ;
!.EC                SLH:units = "W/m2" ;
!.EC        float SHF(rec, jju, ii) ;
!.EC                SHF:long_name = "Surface Sensible Heat Flux" ;
!.EC                SHF:units = "W/m2" ;
!.EC        float SMF(rec, jju, ii) ;
!.EC                SMF:long_name = "Surface Momentum Flux(from N-S & E-W stress)" ;
!.EC                SMF:units = "dens*m^2/s^2" ;
!.EC        float SFQ(rec, jju, ii) ;
!.EC                SFQ:long_name = "Surface(2m) Water Vapor" ;
!.EC                SFQ:units = "kg-H2O/kg-air" ;
!.EC        float MSLP(rec, jju, ii) ;
!.EC                MSLP:long_name = "Mean Sea-Level Pressure  (estimate only)" ;
!.EC                MSLP:units = "hPa" ;
!.EC        float SWVL1(rec, jju, ii) ;
!.EC                SWVL1:long_name = "Volumetric soil water layer 1" ;
!.EC                SWVL1:units = "fraction?" ;
!.EC        float SWVL2(rec, jju, ii) ;
!.EC                SWVL2:long_name = "Volumetric soil water layer 2" ;
!.EC                SWVL2:units = "fraction?" ;
!.EC        float SWVL3(rec, jju, ii) ;
!.EC                SWVL3:long_name = "Volumetric soil water layer 3" ;
!.EC                SWVL3:units = "fraction?" ;
!.EC        float SWVL4(rec, jju, ii) ;
!.EC                SWVL4:long_name = "Volumetric soil water layer 4" ;
!.EC                SWVL4:units = "fraction?" ;
!.EC        float LSPRECFR(rec, jju, ii) ;
!.EC                LSPRECFR:long_name = "Large-scale precipitation fraction (accumulated)" ;
!.EC                LSPRECFR:units = "fraction" ;
!.EC        float UVRAD(rec, jju, ii) ;
!.EC                UVRAD:long_name = "Downward UV radiation at the surface" ;
!.EC                UVRAD:units = "W/m^2?" ;
!.EC        float PHIS(rec, jju, ii) ;
!.EC                PHIS:long_name = "Surface Geopotential (Z)" ;
!.EC                PHIS:units = "m?" ;
!.EC        float LSPREC(rec, jju, ii) ;
!.EC                LSPREC:long_name = "Large Scale Precipitation (stratiform) (accumulated)" ;
!.EC                LSPREC:units = "m?" ;
!.EC        float CONVPREC(rec, jju, ii) ;
!.EC                CONVPREC:long_name = "Convective Precipitation (accumulated)" ;
!.EC                CONVPREC:units = "m?" ;
!.EC        float LSM(rec, jju, ii) ;
!.EC                LSM:long_name = "Land/Sea mask (LSM), NOTE: {0-water,1-land}" ;
!.EC                LSM:units = "flag" ;
!.EC        float SUND(rec, jju, ii) ;
!.EC                SUND:long_name = "Sun Shine Duration (SUND) (accumulated)" ;
!.EC                SUND:units = "s" ;
!.EC        float SNSR(rec, jju, ii) ;
!.EC                SNSR:long_name = "Surface net solar radiation, clear sky (accumulated)" ;
!.EC                SNSR:units = "W/m^2?" ;
!.EC        float FSR(rec, jju, ii) ;
!.EC                FSR:long_name = "Forecast Surface Roughness (FSR) " ;
!.EC                FSR:units = "m?" ;
!
      character (len=MAX_LENGTH_VAR_NAME), parameter :: CONV_MASS_UPD_EC_VNAM = 'CWETE'  ! kg/s
      character (len=MAX_LENGTH_VAR_NAME), parameter :: CONV_MASS_DWN_EC_VNAM = 'CWETD'  ! kg/s
      character (len=MAX_LENGTH_VAR_NAME), parameter :: ENTRAIN_UPD_EC_VNAM = 'CENTU'  ! kg/s
      character (len=MAX_LENGTH_VAR_NAME), parameter :: ENTRAIN_DWN_EC_VNAM = 'CENTD'  ! kg/s
!
      character (len=MAX_LENGTH_VAR_NAME), parameter :: RAIN_LS_EC_VNAM     = 'PRECLS'  !kg/s
      character (len=MAX_LENGTH_VAR_NAME), parameter :: RAIN_CON_EC_VNAM     = 'PRECCON'  !kg/s
      character (len=MAX_LENGTH_VAR_NAME), parameter :: HUMIDITY_EC_VNAM = 'Q'  ! kg/kg
!... must convert to optical depth
      character (len=MAX_LENGTH_VAR_NAME), parameter :: TAU_CLDIC_EC_VNAM= 'CLDLWC'  !
!... must convert to optical depth
      character (len=MAX_LENGTH_VAR_NAME), parameter :: TAU_CLDWA_EC_VNAM= 'CLDIWC'
      character (len=MAX_LENGTH_VAR_NAME), parameter :: CLOUD_FRAC_EC_VNAM = 'CLDFR'  ! 0-1
!
      character (len=MAX_LENGTH_VAR_NAME), parameter :: U10m_EC_VNAM     = 'SFU'      ! m/s
      character (len=MAX_LENGTH_VAR_NAME), parameter :: V10m_EC_VNAM     = 'SFV'      ! m/s
      character (len=MAX_LENGTH_VAR_NAME), parameter :: SURF_AIR_EC_VNAM = 'SFT'      ! K
      character (len=MAX_LENGTH_VAR_NAME), parameter :: SEAICE_EC_VNAM = 'SeaIce'      !
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      character (len=MAX_LENGTH_VAR_NAME) :: dtrn_vnam
      character (len=MAX_LENGTH_VAR_NAME) :: gwet_VNAM
      character (len=75) :: err_msg
!
      integer :: ij, ik, ilong, ilat, ivert
!
      integer :: cnt3d (3), cnt4d (4)
      integer :: strt3d(3), strt4d(4)
      integer :: cnt1d (1), strt1d(1)
!
      real*8  :: cos_correction, rij, rtotj, geos4moistq_conv
!
      real*8  :: lwi_flags_real (i1:i2, ju1:j2)
      real*8  :: ci_flags_real (i1:i2, ju1:j2)
      real*8  :: snow_cover (i1:i2, ju1:j2)
      real*8  :: frlandice (i1:i2, ju1:j2)
      real*8  :: frOcean (i1:i2, ju1:j2)
      real*8  :: frLake (i1:i2, ju1:j2)
      real*8  :: frSeaice (i1:i2, ju1:j2)
!
!... GEOS4 and EC have cloud optical depth split between ice and water - add for total
      real*8  :: tauice_cloud    (i1:i2, ju1:j2, k1:k2)
      real*8  :: tauwater_cloud    (i1:i2, ju1:j2, k1:k2)
!
!... work arrays for making sure precip is (+)ve
      logical :: mask (i1:i2, ju1:j2)
      real*8  :: zero (i1:i2, ju1:j2)
!
      !... GEOS5 precipitation fields
      real*8  :: precanv      (i1:i2, ju1:j2)
      real*8  :: precsno      (i1:i2, ju1:j2)
      real*8  :: preclsc      (i1:i2, ju1:j2)
!
      !... GEOS5 arrays for inverting the vertical levels
      real*8 , allocatable :: tempVar(:,:,:)
      character(len=ESMF_MAXSTR) :: IAm
!
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = "Read_Met2"
!
      if (pr_diag) Write (6,*) IAm, ' called by ', procID
!
      if (metdata_name_model(1:5) .eq. "GEOS5") then
         kzz_vnam           = 'KH'
         tau_cldic_geo_vnam = 'TAUCLI'
         tau_cldwa_geo_vnam = 'TAUCLW'
         humidity_vnam      = 'QV'
         grnd_temp_vnam     = 'TSKIN'
         pbl_vnam           = 'PBLH'
         radswg_vnam        = 'SWGNET'
         surf_rough_vnam    = 'Z0M'
         ustar_vnam         = 'USTAR'
         con_precip_vnam    = 'PRECCON'
         u10m_vnam          = 'U10M'
         v10m_vnam          = 'V10M'
         gwet_vnam1         = 'GWETROOT'
         gwet_vnam2         = 'GWETTOP'
         pardif_vnam        = 'PARDF'
         pardir_vnam        = 'PARDR'
         lwi_vnam           = 'LWI'
         snow_vnam          = 'SNODP'
         sf_alb_ldif_vnam   = 'ALBNIRDF'
         sf_alb_ldir_vnam   = 'ALBNIRDR'
         sf_alb_sdif_vnam   = 'ALBVISDF'
         sf_alb_sdir_vnam   = 'ALBVISDR'
         cmf_vnam           = 'CMFMC'
         dtrn_DAO_vnam      = 'DTRAIN'
         cloud_frac_vnam    = 'CLOUD'
!
         allocate (tempVar(i1:i2,ju1:j2,k1:k2))
      else
         KZZ_VNAM           = 'kzz'
         TAU_CLDIC_GEO_VNAM = 'taucli'
         TAU_CLDWA_GEO_VNAM = 'tauclw'
         HUMIDITY_VNAM      = 'sphu'
         GRND_TEMP_VNAM     = 'tg'
         PBL_VNAM           = 'pbl'
         RADSWG_VNAM        = 'radswg'
         SURF_ROUGH_VNAM    = 'z0'
         USTAR_VNAM         = 'ustar'
         CON_PRECIP_VNAM    = 'precon'
         U10m_VNAM          = 'u10m'
         V10m_VNAM          = 'v10m'
         gwet_vnam1         = 'gwet1'
         gwet_vnam2         = 'gwettop'
         PARDIF_VNAM        = 'pardif'
         PARDIR_VNAM        = 'pardir'
         LWI_VNAM           = 'lwi'
         SNOW_VNAM          = 'snowh'
         SF_ALB_LDIF_VNAM   = 'saldif'
         SF_ALB_LDIR_VNAM   = 'saldir'
         SF_ALB_SDIF_VNAM   = 'sasdif'
         SF_ALB_SDIR_VNAM   = 'sasdir'
         CMF_VNAM           = 'cldmas'
         DTRN_DAO_VNAM      = 'dtrain'
         CLOUD_FRAC_VNAM    = 'cloud'
      end if
!
!.sds... EC met fields
      if (metdata_name_org(1:2) == 'EC') then
         PBL_VNAM           = 'BLH'
         SF_ALB_LDIF_VNAM   = 'SA'
         SF_ALB_LDIR_VNAM   = 'SA'
         SF_ALB_SDIF_VNAM   = 'SA'
         SF_ALB_SDIR_VNAM   = 'SA'
         CMF_VNAM           = 'CWETE'
         GRND_TEMP_VNAM     = 'SFT'
         SURF_AIR_TEMP_VNAM = 'SFT'
         LWI_VNAM           = 'LSM'
         DTRN_NCAR_VNAM     = 'CENTU'
         snow_vnam          = 'SnowDepth'
      endif
!
!
      ilong = i2 - i1  + 1
      ilat  = j2 - ju1 + 1
      ivert = k2 - k1  + 1
!
      if (gwet_opt == 0) then
         gwet_VNAM = gwet_vnam1
      else
         gwet_VNAM = gwet_vnam2
      end if
!
!... conversion factor for moistq field to put into units of g/kg/day
      geos4moistq_conv = (1000.d0 * 86400.d0)
!
      strt3d(1) = i1  -  i1_gl + 1
      strt3d(2) = ju1 - ju1_gl + 1
      strt3d(3) = m2rnum_in
!
      cnt3d (:) = (/ ilong, ilat, 1 /)
!
      call Ncrd_3d (grnd_temp, ncid_met2,GRND_TEMP_VNAM, strt3d,cnt3d)
      call Ncrd_3d (pbl,       ncid_met2,PBL_VNAM,       strt3d,cnt3d)
!.sds...
!.bad fix for thse fields, but what to use?
      if (metdata_name_org(1:2) == 'EC') then
        radswg = 1
        surf_rough = 1
        ustar = 1
      else
        call Ncrd_3d (radswg,    ncid_met2,RADSWG_VNAM,    strt3d,cnt3d)
        call Ncrd_3d (surf_rough, ncid_met2, SURF_ROUGH_VNAM,strt3d,cnt3d)
        call Ncrd_3d (ustar,      ncid_met2, USTAR_VNAM,     strt3d,cnt3d)
      endif
!
      !... use rain production rate from large scale and convection
      if (metdata_name_model(1:5) == 'GEOS5') then
         call ncrd_3d (con_precip, ncid_met2, CON_PRECIP_VNAM,strt3d,cnt3d)
         call ncrd_3d (precanv,ncid_met2, PRECANV_VNAM,strt3d,cnt3d)
         call ncrd_3d (precsno,ncid_met2, PRECSNO_VNAM,strt3d,cnt3d)
         call ncrd_3d (preclsc,ncid_met2, PRECLSC_VNAM,strt3d,cnt3d)
         tot_precip(:,:) = precanv(:,:) + precsno(:,:) + preclsc(:,:) + con_precip(:,:)
!
!.sds.. convert from kg/m^2/s to mm/day divide by 1000kg/m^3 and mult by 1000 mm/m and mult by 86400 s/day
         con_precip(:,:) = con_precip(:,:) * 86400.0d0
         tot_precip(:,:) = tot_precip(:,:) * 86400.0d0
      elseif (metdata_name_org(1:2) /= 'EC') then
         call Ncrd_3d (tot_precip, ncid_met2, TOT_PRECIP_VNAM,strt3d,cnt3d)
         call Ncrd_3d (con_precip, ncid_met2, CON_PRECIP_VNAM,strt3d,cnt3d)
      end if
!
!... use 10m temp for GEOS3
      if (metdata_name_model(1:5) == 'GEOS3') then
        call Ncrd_3d (surf_air_temp, ncid_met2, SURF_AIR_GEOS_VNAM,strt3d,cnt3d)
!... use 2m temp for GEOS4
      elseif (metdata_name_model(1:5) == 'GEOS4') then
        call Ncrd_3d (surf_air_temp, ncid_met2, SURF_AIR_GEO4_VNAM,strt3d,cnt3d)
      elseif (metdata_name_model(1:5) == 'GEOS5') then
        call ncrd_3d (surf_air_temp, ncid_met2, SURF_AIR_GEO5_VNAM,strt3d,cnt3d)
!... use 'ts' for others
      else
        call Ncrd_3d (surf_air_temp, ncid_met2, SURF_AIR_TEMP_VNAM,strt3d,cnt3d)
      end if
!
!
      if ((metdata_name_model(1:5) == 'GEOS4') .or. &
     &    (metdata_name_model(1:5) == 'GEOS5')) then
         call Ncrd_3d(gwet  , ncid_met2, gwet_VNAM,strt3d,cnt3d)
         call Ncrd_3d(u10m  , ncid_met2, U10m_VNAM,strt3d,cnt3d)
         call Ncrd_3d(v10m  , ncid_met2, V10m_VNAM,strt3d,cnt3d)
         call Ncrd_3d(pardif, ncid_met2, PARDIF_VNAM,strt3d,cnt3d)
         call Ncrd_3d(pardir, ncid_met2, PARDIR_VNAM,strt3d,cnt3d)
      end if
!
!.sds... EC 10m winds
      if (metdata_name_org(1:2) == 'EC') then
         call Ncrd_3d(u10m  , ncid_met2, U10m_EC_VNAM,strt3d,cnt3d)
         call Ncrd_3d(v10m  , ncid_met2, V10m_EC_VNAM,strt3d,cnt3d)
      endif
!
!     -----------------------------------------------
!     Test to see if Land water ice flags are needed.
!     -----------------------------------------------
!
      if ((chem_opt     /= 0) .or.  &
     &    (convec_opt   /= 0) .or.  &
     &    (emiss_in_opt /= 0)) then
!
!       --------------------------------------------------------------------
!       Include here any data sets that you would like lwi (land water ice)
!       information to be read from the met data rather than simply set with
!       hardwired data staements. For 1x1 degree resolution it must be read
!       from the met data.
!       --------------------------------------------------------------------
!
        if ((metdata_name_org  (1:3) == 'DAO') .and.  &
     &      (metdata_name_model(1:5) == 'GEOS3')) then
!
          call Ncrd_3d (lwi_flags_real, ncid_met2, LWI_VNAM, strt3d, cnt3d)
!
          lwi_flags(:,:) = lwi_flags_real(:,:)
!
!         -------------------------------------------------------------
!         Convert DAO GEOS3 LWI to GMIMOD's LWI (1 water 2 land 3 ice).
!         -------------------------------------------------------------
!
          where (lwi_flags <= 8)   lwi_flags = 2
          where (lwi_flags == 10)  lwi_flags = 2
          where (lwi_flags == 9)   lwi_flags = 3
          where (lwi_flags == 101) lwi_flags = 3
          where (lwi_flags == 100) lwi_flags = 1
!
        end if
!
!
!... get land-water-ice flags for GEOS4/GEOS5 - already 1-3
        if (metdata_name_org  (1:4) == 'GMAO') then
!
           if (metdata_name_model(1:5) == 'GEOS4') then
              call Ncrd_3d (lwi_flags_real, ncid_met2, LWI_VNAM, strt3d,cnt3d)
           endif
!
           if (metdata_name_model(1:5) == 'GEOS5') then
              call Ncrd_3d (frlandice, ncid_met2, 'FRLANDICE', strt3d, cnt3d)
!
              ! Megan Damon Sept 2010
              ! LWI not in MERRA output stream
              ! FROCEAN, FRLAKE, and FRLANDICE are constant
              if (metdata_name_model(6:10) == "MERRA") then
                 call Ncrd_3d (frOcean, ncid_met2, FROCEAN_VNAM, strt3d, cnt3d)
                 call Ncrd_3d (frLake, ncid_met2, FRLAKE_VNAM, strt3d, cnt3d)
                 call Ncrd_3d (frSeaice, ncid_met2, FRSEAICE_VNAM, strt3d, cnt3d)
                 lwi_flags_real = 1.0 ! land
                 where ( frOcean+frLake >= 0.6)  lwi_flags_real = 0.0 ! water
                 where ( lwi_flags_real == 0 .and. frSeaice > 0.5)  lwi_flags_real = 2.0 ! ice
                 where ( lwi_flags_real == 0 .and. grnd_temp < 271.40) lwi_flags_real = 2.0 ! ice
                 where ( frLandIce > 0.5) lwi_flags_real = 2.0 ! land ice
              endif
           endif
!
           !... Convert DAO GEOS4/GEOS5 LWI to GMIMOD's LWI (1 water 2 land 3 sea-ice)
           lwi_flags(:,:) = lwi_flags_real(:,:) + 1
!
           !... Read in snow cover to update ice flag to include land-ice
           call Ncrd_3d (snow_cover, ncid_met2, SNOW_VNAM, strt3d, cnt3d)
          if (metdata_name_model(1:5) == 'GEOS4') then
             !... set ice flag where 1 inch of snow or more (assuming 1 to 10 ratio)
             !where (snow_cover >= 0.00254) lwi_flags = 3
             !... set ice flag where 10 inch of snow or more (assuming 1 to 10 ratio)
             !where (snow_cover >= 0.0254) lwi_flags = 3
             !... set ice flag where 10 cm of water equivalent ice
             !... GEOS4 is water equivalent
             where (snow_cover(:,:) >= 0.1) lwi_flags(:,:) = 3
          else if (metdata_name_model(1:5) == 'GEOS5') then
             !... GEOS5 is snow depth
             where (snow_cover(:,:) >= 0.35 .and. snow_cover(:,:) < 1000.) lwi_flags(:,:) = 3
             !... GEOS5 has a land ice parameter for glaciated areas (don't seem to have any snow)
             where (frlandice(:,:) >= 0.1) lwi_flags(:,:) = 3
          endif
!
       endif
!
!... get land-water-ice flags for EC use LSM and CI
        if (metdata_name_org (1:2) == 'EC') then
          call Ncrd_3d (lwi_flags_real, ncid_met2, LWI_VNAM, strt3d, cnt3d)
          call Ncrd_3d (ci_flags_real, ncid_met2, SEAICE_EC_VNAM, strt3d, cnt3d)
          !... Convert EC LSM (0 water 1 land) and CI to LWI (1 water 2 land 3 sea-ice)
          lwi_flags(:,:) = lwi_flags_real(:,:) + 1
          where(ci_flags_real(:,:) .gt. 0.5) lwi_flags(:,:) = 3
          !... get EC snow depth
          call Ncrd_3d (snow_cover, ncid_met2, snow_vnam, strt3d, cnt3d)
          !... EC is snow depth? water equivalent? this empirically determined
          where (snow_cover(:,:) >= 0.11) lwi_flags(:,:) = 3
        endif
!
!... get land-water-ice flags for EC use LSM and CI
        if (metdata_name_org (1:2) == 'EC') then
          call Ncrd_3d (lwi_flags_real, ncid_met2, LWI_VNAM, strt3d, cnt3d)
          call Ncrd_3d (ci_flags_real, ncid_met2, SEAICE_EC_VNAM, strt3d, cnt3d)
          !... Convert EC LSM (0 water 1 land) and CI to LWI (1 water 2 land 3 sea-ice)
          lwi_flags(:,:) = lwi_flags_real(:,:) + 1
          where(ci_flags_real(:,:) .gt. 0.5) lwi_flags(:,:) = 3
          !... get EC snow depth
          call Ncrd_3d (snow_cover, ncid_met2, snow_vnam, strt3d, cnt3d)
          !... EC is snow depth? water equivalent? this empirically determined
          where (snow_cover(:,:) >= 0.11) lwi_flags(:,:) = 3
        endif
      end if
!
      if (sfalbedo_opt == 3) then
!
        if ((metdata_name_org  (1:4) == 'NCAR') .or.  &
     &      (metdata_name_model(1:4) == 'CCM3') .or.  &
     &      (metdata_name_model(1:5) == 'GEOS4') .or. &
     &      (metdata_name_org  (1:2) == 'EC') .or. &
     &      (metdata_name_model(1:5) == 'GEOS5')) then
!
          call Ncrd_3d (saldif, ncid_met2, SF_ALB_LDIF_VNAM, strt3d, cnt3d)
          call Ncrd_3d (saldir, ncid_met2, SF_ALB_LDIR_VNAM, strt3d, cnt3d)
          call Ncrd_3d (sasdif, ncid_met2, SF_ALB_SDIF_VNAM, strt3d, cnt3d)
          call Ncrd_3d (sasdir, ncid_met2, SF_ALB_SDIR_VNAM, strt3d, cnt3d)
!
        else
!
          err_msg = 'No sfalbedo met data problem in Read_Met2.'
          call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
!
        end if
!
      end if
!
!
      if (uvalbedo_opt == 3) then
!
        if ((metdata_name_org  (1:4) == 'GISS') .or.  &
     &      (metdata_name_model(1:2) == 'GS')) then
!
          call Ncrd_3d (surf_alb_uv, ncid_met2, SURF_ALB_UV_VNAM, strt3d, cnt3d)
!
        else
!
          err_msg = 'No surf_alb_uv met data problem in Read_Met2.'
          call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
!
        end if
!
      end if
!
      strt4d(1) = i1  -  i1_gl + 1
      strt4d(2) = ju1 - ju1_gl + 1
      strt4d(3) = k1
      strt4d(4) = m2rnum_in
!
      cnt4d (:) = (/ ilong, ilat, ivert, 1 /)
!
      if (convec_opt /= 0 .or. do_wetdep) then
!
        call Ncrd_4d (cmf, ncid_met2, CMF_VNAM, strt4d, cnt4d)
!
!... Inverting vertical levels for GEOS5
        if (metdata_name_model(1:5) == 'GEOS5') then
           tempVar(:,:,k1:k2) = cmf    (:,:,k2:k1:-1)
           cmf    (:,:,k1:k2) = tempVar(:,:,k1:k2)
        end if
!
!... convert to kg/m^2/s from Pa/s for GEOS4
        if (metdata_name_org(1:4) == 'GMAO' .and. metdata_name_model(1:5) == 'GEOS4') then
           cmf = cmf/GMI_G
        end if
!
!... convert to kg/m^2/s from kg/s for EC-Oslo
        if (metdata_name_org(1:2) == 'EC') then
           cmf(:,:,k1:k2) = cmf(:,:,k1:k2) / spread(mcor(:,:),3,k2-k1+1)
        end if
!
        if (convec_opt /= 0) then
!
          if (metdata_name_org(1:3) == 'DAO' .or.  &
     &        metdata_name_org(1:4) == 'GMAO') then
            dtrn_vnam = DTRN_DAO_VNAM
          else
            dtrn_vnam = DTRN_NCAR_VNAM
          endif
!
!... convert to kg/m^2/s from kg/s for EC-Oslo has no detrainment - calc as residual
          if (metdata_name_org(1:2) == 'EC') then
             do ik=k1+1,k2
                dtrn(i1:i2,ju1:j2,ik) = cmf(i1:i2,ju1:j2,ik-1) - cmf(i1:i2,ju1:j2,ik)
             enddo
          else
             call Ncrd_4d (dtrn, ncid_met2, dtrn_vnam, strt4d, cnt4d)
          end if
!
          ! Inverting vertical levels
          if (metdata_name_model(1:5) == 'GEOS5') then
             tempVar(:,:,k1:k2) = dtrn   (:,:,k2:k1:-1)
             dtrn   (:,:,k1:k2) = tempVar(:,:,k1:k2)
          end if
!
!... convert to kg/m^2/s from Pa/s for GEOS4
          if (metdata_name_org(1:4) == 'GMAO' .and. metdata_name_model(1:5) == 'GEOS4') dtrn = dtrn/GMI_G
!
        endif
!
      endif
!
      if (convec_opt == 2) then
!
         if (metdata_name_org(1:3) == 'DAO' .or. metdata_name_model(1:5) == 'GEOS5') then
            do ik = k1, k2-1
               eu(:,:,ik) = cmf(:,:,ik+1) - cmf(:,:,ik) + dtrn(:,:,ik)
            enddo
            ed = 0.0
            md = 0.0
!
         else if (metdata_name_org(1:2) == 'EC') then
            !... read and convert downdraft flux to /s from kg/s for EC-Oslo
               call Ncrd_4d (eu, ncid_met2, ENTRAIN_UPD_EC_VNAM, strt4d, cnt4d)
               eu(:,:,k1:k2) = eu(:,:,k1:k2) / spread(mcor(:,:),3,k2-k1+1)
            !... read and convert downdraft flux to /s from kg/s for EC-Oslo
               call Ncrd_4d (ed, ncid_met2, ENTRAIN_DWN_EC_VNAM, strt4d, cnt4d)
               ed(:,:,k1:k2) = ed(:,:,k1:k2) / spread(mcor(:,:),3,k2-k1+1)
            !... read and convert downdraft flux to /s from kg/s for EC-Oslo
               call Ncrd_4d (md, ncid_met2, CONV_MASS_DWN_EC_VNAM, strt4d, cnt4d)
               md(:,:,k1:k2) = md(:,:,k1:k2) / spread(mcor(:,:),3,k2-k1+1)
            !... for convec=2 add in updraft entrain to total detrain
               dtrn(i1:i2,ju1:j2,k1:k2) = dtrn(i1:i2,ju1:j2,k1:k2) + eu(i1:i2,ju1:j2,k1:k2)
!
         else if (metdata_name_org(1:4) == 'GMAO' .and. metdata_name_model(1:5) == 'GEOS4') then
!
!... Zhang-McFarland and hack fluxes and convert to kg/m2/s from Pa/s
               call Ncrd_4d (zmed, ncid_met2, ETRNDN_ZM_GEO4, strt4d, cnt4d)
               zmed = zmed/GMI_G
               call Ncrd_4d (zmeu, ncid_met2, ETRNUP_ZM_GEO4, strt4d, cnt4d)
               zmeu = zmeu/GMI_G
               call Ncrd_4d (hkeu, ncid_met2, ETRNUP_HK_GEO4, strt4d, cnt4d)
               hkeu = hkeu/GMI_G
               call Ncrd_4d (zmmd, ncid_met2, MFLXDN_ZM_GEO4, strt4d, cnt4d)
               zmmd = zmmd/GMI_G
!
!... for convec_opt=2 collapse entrain for separate convection
!...  to appropriate arrays
               ed = zmed
               eu = zmeu + hkeu
               md = zmmd
!... all others
          else
               call Ncrd_4d (eu, ncid_met2, EU_VNAM, strt4d, cnt4d)
               call Ncrd_4d (ed, ncid_met2, ED_VNAM, strt4d, cnt4d)
               call Ncrd_4d (md, ncid_met2, MD_VNAM, strt4d, cnt4d)
          endif
!
      endif
!
!
!!!!!!! For the full GEOS4 convection
      if (convec_opt == 3 .and. metdata_name_org(1:4) == 'GMAO' .and.  &
     &    metdata_name_model(1:5) == 'GEOS4') then
!
!... Zhang-McFarland fluxes and convert to kg/m2/s from Pa/s
        call Ncrd_4d (zmeu, ncid_met2, ETRNUP_ZM_GEO4, strt4d, cnt4d)
        zmeu = zmeu/GMI_G
        call Ncrd_4d (zmed, ncid_met2, ETRNDN_ZM_GEO4, strt4d, cnt4d)
        zmed = zmed/GMI_G
        call Ncrd_4d (zmmd, ncid_met2, MFLXDN_ZM_GEO4, strt4d, cnt4d)
        zmmd = zmmd/GMI_G
        call Ncrd_4d (zmmu, ncid_met2, MFLXUP_ZM_GEO4, strt4d, cnt4d)
        zmmu = zmmu/GMI_G
!... don't read the detrainment, but calculate as residual
!        call Ncrd_4d (zmdu, ncid_met2, DTRNUP_ZM_GEO4, strt4d, cnt4d)
!        zmdu = zmdu/GMI_G
        do ik=k1+1,k2
          zmdu(i1:i2,ju1:j2,ik) = zmeu(i1:i2,ju1:j2,ik) +  &
     &       zmmu(i1:i2,ju1:j2,ik-1) - zmmu(i1:i2,ju1:j2,ik)
        enddo
!
!... hack fluxes and convert to kg/m2/s from Pa/s
        call Ncrd_4d (hkeu, ncid_met2, ETRNUP_HK_GEO4, strt4d, cnt4d)
        hkeu = hkeu/GMI_G
        call Ncrd_4d (hkmu, ncid_met2, MFLXUP_HK_GEO4, strt4d, cnt4d)
        hkmu = hkmu/GMI_G
!... don't read the detrainment, but calculate as residual
!        call Ncrd_4d (hkdu, ncid_met2, DTRNUP_HK_GEO4, strt4d, cnt4d)
!        hkdu = hkdu/GMI_G
        do ik=k1+1,k2
          hkdu(i1:i2,ju1:j2,ik) = hkeu(i1:i2,ju1:j2,ik) +  &
     &       hkmu(i1:i2,ju1:j2,ik-1) - hkmu(i1:i2,ju1:j2,ik)
        enddo
      endif
!
!     --------------------------------------------------------------------
!     CORRECT ERRORS IN THE GISS 2PRIME CONVECTION UNITS.
!
!     The LLNL GISS data was multiplied by g*0.01; to remove this error we
!     multiply by 100./g.  Also the data was assummed to be at a fixed 45
!     degree latitude grid box;  to correct that error multiply by
!     cos(45)/cos(lat).
!     At this point the units of the downdraft entrainment is still fraction of box/s
!     --------------------------------------------------------------------
!
      if ((convec_opt /= 0) .and.  &
     &    (metdata_name_org  (1:4) == 'GISS' ) .and.  &
     &    (metdata_name_model(1:6) == '2prime')) then
!
        cmf(:,:,:) = cmf(:,:,:) * 100.0d0 / GMI_G
        eu (:,:,:) = eu (:,:,:) * 100.0d0 / GMI_G
        ed (:,:,:) = ed (:,:,:) * 100.0d0 / GMI_G
        md (:,:,:) = md (:,:,:) * 100.0d0 / GMI_G
!
        do ij = ju1, j2
!
          rij = ij
          rtotj = (j2_gl - ju1_gl)
          cos_correction = cos(45.0d0 * RADPDEG) /  &
     &       cos(((((rij-1) / rtotj) * 180.0d0) - 90.0d0) * RADPDEG)
!
          cmf(:,ij,:) = cmf(:,ij,:) * cos_correction
!
        enddo
!
!
! jhk cmf layer lifting
!
        do ik = k2, k1+1, -1
           cmf(:,:,ik) = cmf(:,:,ik-1)
        enddo
!
        cmf(:,:,k1) = cmf(:,:,k1) * 0.2d0
        cmf(:,:,13:k2) = 0.0d0
!
!
! end of cmf layer lifting
!
      end if
!
!... get precipitation and moistq
      if (((metdata_name_org  (1:4) == 'NCAR' ) .and.  &
     &     (metdata_name_model(1:5) == 'MATCH')) .or.  &
     &    (metdata_name_org(1:4) == 'GISS')) then
!
         call Ncrd_4d (rain, ncid_met2, RAIN_VNAM, strt4d, cnt4d)
!
!       ---------------------------------------------------------------
!       Convective and total precipitation were incorrect in the GISS
!       data set.  The ratio was OK.  We use the bottom layer of the 3D
!       rain variable for the total precipitation hitting the ground.
!       ---------------------------------------------------------------
!
         if ((metdata_name_org  (1:4) == 'GISS') .and.  &
     &       (metdata_name_model(1:6) == '2prime')) then
!
            where(tot_precip > 0.0d0)
               con_precip(:,:) = rain(:,:,1) * (con_precip(:,:) / tot_precip(:,:))
            elsewhere
               con_precip(:,:) = 0.0d0
            endwhere
!
            tot_precip(:,:) = rain(:,:,1)
!
         endif
!
      elseif (metdata_name_org  (1:4) == 'GMAO') then
!
         ! Read in convective and large scale rain fields
!
         if (metdata_name_model(1:5) == 'GEOS4') then
!... get GEOS4 3D rain fields
            call Ncrd_4d (rain_zm, ncid_met2, RAIN_ZM_GEO4, strt4d, cnt4d)
            call Ncrd_4d (rain_hk, ncid_met2, RAIN_HK_GEO4, strt4d, cnt4d)
            call Ncrd_4d (rain_ls, ncid_met2, RAIN_LS_GEO4, strt4d, cnt4d)
!
!... convert units from kg/kg/s to g/kg/day
            rain_zm(:,:,:) = rain_zm(:,:,:) * geos4moistq_conv
            rain_hk(:,:,:) = rain_hk(:,:,:) * geos4moistq_conv
            rain_ls(:,:,:) = rain_ls(:,:,:) * geos4moistq_conv
!
!... put total rain production into moistq for GEOS4 fields
            moistq(:,:,:) = rain_ls(:,:,:)+rain_zm(:,:,:)+rain_hk(:,:,:)
!
            moistq(:,:,:) = -moistq(:,:,:)
         endif
!
         if (metdata_name_model(1:5) == 'GEOS5') then
!
         ! Megan Damon September 2010
         ! DQRCON and DQRLSC not saved on native grid in MERRA
         ! Moved calculation here to save processing time
            if (metdata_name_model(1:10) == "GEOS5MERRA") then
               kSave = cnt4d(3)
               cnt4d(3) = kSave + 1
               call ncrd_4d (pficu, ncid_met2, PFICU_VNAM, strt4d, cnt4d)
               call ncrd_4d (pflcu, ncid_met2, PFLCU_VNAM, strt4d, cnt4d)
               call ncrd_4d (pfilsan, ncid_met2, PFILSAN_VNAM, strt4d, cnt4d)
               call ncrd_4d (pfllsan, ncid_met2, PFLLSAN_VNAM, strt4d, cnt4d)
               cnt4d(3) = kSave
               do ik = k1, k2
                  rain_zm(:,:,ik) = pficu(:,:,ik+1)-pficu(:,:,ik)+pflcu(:,:,ik+1)-pflcu(:,:,ik)
                  rain_ls(:,:,ik) = pfilsan(:,:,ik+1)-pfilsan(:,:,ik)+pfllsan(:,:,ik+1)-pfllsan(:,:,ik)
               enddo
               call Ncrd_4d (delp, ncid_met2, DELP_VNAM, strt4d, cnt4d)
               rain_zm(:,:,:) = rain_zm(:,:,:)/delp(:,:,:) * 9.81
               rain_ls(:,:,:) = rain_ls(:,:,:)/delp(:,:,:) * 9.81
               ! no hack rain provided by GEOS5
               rain_hk(:,:,:) =  0.0d0
!... rain production must be positive
               where(rain_zm(:,:,:) <= 0.0d0) rain_zm(:,:,:) = 0.0
               where(rain_ls(:,:,:) <= 0.0d0) rain_ls(:,:,:) = 0.0
!... convert units from kg/kg/s to g/kg/day
               rain_zm(:,:,:) = rain_zm(:,:,:) * geos4moistq_conv
               rain_ls(:,:,:) = rain_ls(:,:,:) * geos4moistq_conv
               ! Inverting the vertical levels
               tempVar(:,:,k1:k2) = rain_zm(:,:,k2:k1:-1)
               rain_zm(:,:,k1:k2) = tempVar(:,:,k1:k2)
               tempVar(:,:,k1:k2) = rain_ls(:,:,k2:k1:-1)
               rain_ls(:,:,k1:k2) = tempVar(:,:,k1:k2)
            else
               call ncrd_4d (rain_zm, ncid_met2, RAIN_ZM_GEO5, strt4d, cnt4d)
               call ncrd_4d (rain_ls, ncid_met2, RAIN_LS_GEO5, strt4d, cnt4d)
!
               ! Inverting the vertical levels
               tempVar(:,:,k1:k2) = rain_zm(:,:,k2:k1:-1)
               rain_zm(:,:,k1:k2) = tempVar(:,:,k1:k2)
               tempVar(:,:,k1:k2) = rain_ls(:,:,k2:k1:-1)
               rain_ls(:,:,k1:k2) = tempVar(:,:,k1:k2)
               ! no hack rain provided by GEOS5
               rain_hk(:,:,:) =  0.0d0
!... convert units from kg/kg/s to g/kg/day
               rain_zm(:,:,:) = rain_zm(:,:,:) * geos4moistq_conv
               rain_ls(:,:,:) = rain_ls(:,:,:) * geos4moistq_conv
            endif
!
            if (metdata_name_model(1:8) == "GEOS5DAS" .or. &
     &          metdata_name_model(1:8) == "GEOS5CCM") then
               ! get GEOS5 tendencies from moist physics fields
               call ncrd_4d (moistq, ncid_met2, DQVDTMST_VNAM, strt4d, cnt4d)
               ! Inverting the vertical levels
               tempVar(:,:,k1:k2) = moistq (:,:,k2:k1:-1)
               moistq (:,:,k1:k2) = tempVar(:,:,k1:k2)
               ! and convert units from kg/kg/s to g/kg/day
               moistq(:,:,:) = moistq(:,:,:) * GEOS4moistq_conv
!
            elseif (metdata_name_model(1:10) == "GEOS5MERRA") then
               !... put total rain production into moistq for GEOS5merra fields
               moistq(:,:,:) = rain_ls(:,:,:)+rain_zm(:,:,:)
               !... convert units and reverse sign
               moistq(:,:,:) = -moistq(:,:,:)
            else
               err_msg = "The metdata name is not supported"
               call GmiPrintError (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
            end if
         endif
!
      !... EC_Oslo fields
      elseif (metdata_name_org(1:2) == 'EC') then
         call ncrd_4d (rain_ls, ncid_met2, RAIN_LS_EC_VNAM, strt4d, cnt4d)
         call ncrd_4d (rain_zm, ncid_met2, RAIN_CON_EC_VNAM, strt4d, cnt4d)
         rain = rain_ls + rain_zm
         !... difference rain to get moistq
         moistq(:,:,k1:k2-1) = (rain(:,:,k1:k2-1) - rain(:,:,k1+1:k2))
         !... set top level
         moistq(:,:,k1:k2  ) = rain(:,:,k1:k2)
         !... convert from kg/s to g/kg/day
         moistq(:,:,k1:k2-1) = moistq(:,:,k1:k2) * GPKG * SECPDY / mass(:,:,k1:k2)
         !... fill in precip at ground (units?)
         con_precip(:,:) = moistq(:,:,k1) * 86400.0d0 * 0.0
         tot_precip(:,:) = moistq(:,:,k1) * 86400.0d0
!
      else
        call Ncrd_4d (moistq, ncid_met2, MOISTQ_VNAM, strt4d, cnt4d)
!
      endif
!
      if (((metdata_name_org  (1:4) == 'NCAR') .and.  &
     &     (metdata_name_model(1:4) == 'CCM3'))) then
!
        call Ncrd_4d (astrat, ncid_met2,STRT_CLD_FRAC_VNAM,strt4d,cnt4d)
        call Ncrd_4d (clwc,   ncid_met2,CLD_LIQ_WATER_VNAM,strt4d,cnt4d)
        call Ncrd_4d (aconv,  ncid_met2,CONV_CLD_FRAC_VNAM,strt4d,cnt4d)
        call Ncrd_4d (relhum, ncid_met2,REL_HUM_VNAM      ,strt4d,cnt4d)
        call Ncrd_4d (rhclear,ncid_met2,REL_HUM_CLEAR_VNAM,strt4d,cnt4d)
!
        call Ncrd_4d (cldmas, ncid_met2,CONV_MASS_UPD_VNAM,strt4d,cnt4d)
        call Ncrd_4d (omega,  ncid_met2,VERT_PRSE_VEL_VNAM,strt4d,cnt4d)
!
        max_cloud(:,:,:) = 0.0d0
        ran_cloud(:,:,:) = aconv(:,:,:) + astrat(:,:,:)
!
!       --------------------------------------------------------------------
!       Adjust the cloud liquid water to make it an average of the grid box.
!       --------------------------------------------------------------------
!
        clwc(:,:,:) = clwc(:,:,:) * (aconv(:,:,:) + astrat(:,:,:))
!
      else if ((metdata_name_org  (1:4) == 'GMAO') .and.  &
     &        ((metdata_name_model(1:5) == 'GEOS4') .or.  &
     &         (metdata_name_model(1:5) == 'GEOS5'))) then
!.sds         if (metdata_name_model(1:13) == "GEOS5MERRA300") then
!.sds            ! read in three clouds
!.sds            ! use three clouds to get cloud area fraction
!.sds            ! sum tau clouds and then multiple by cloud area fraction
!.sds            call ncrd_4d (cfan, ncid_met2,CFAN_VNAM, strt4d, cnt4d)
!.sds            call ncrd_4d (cfcu, ncid_met2,CFCU_VNAM, strt4d, cnt4d)
!.sds            call ncrd_4d (cfls, ncid_met2,CFLS_VNAM, strt4d, cnt4d)
!.sds            max_cloud(:,:,:) = cfan(:,:,:)+cfcu(:,:,:)+ cfls(:,:,:)
!.sds!... why this instead of above?            max_cloud(:,:,:) = real(0)
!.sds         else
!.sds            !... looks like max_cloud and ran_cloud are used to calc cloud fraction
        call Ncrd_4d (max_cloud,ncid_met2,CLOUD_FRAC_VNAM,strt4d,cnt4d)
!.sds         endif
        ran_cloud(:,:,:) = 0.0d0
!... get specific humidity and convert to g/kg
        call Ncrd_4d (humidity,  ncid_met2, HUMIDITY_VNAM, strt4d,cnt4d)
!
        ! Inverting the vertical levels
        if (metdata_name_model(1:5) == 'GEOS5') then
           tempVar  (:,:,k1:k2) = max_cloud(:,:,k2:k1:-1)
           max_cloud(:,:,k1:k2) = tempVar  (:,:,k1:k2)
!
           tempVar (:,:,k1:k2) = humidity(:,:,k2:k1:-1)
           humidity(:,:,k1:k2) = tempVar (:,:,k1:k2)
        end if
!
        humidity(:,:,:) = humidity(:,:,:) * GPKG
!
      elseif (metdata_name_org(1:2) == 'EC') then
        !... EC specific humidity and convert from kg/kg to g/kg
        call Ncrd_4d (humidity,  ncid_met2, HUMIDITY_EC_VNAM, strt4d,cnt4d)
        humidity(:,:,:) = humidity(:,:,:)*GPKG
        !... assuming EC cloud fraction is maximum overlap
        call Ncrd_4d (max_cloud, ncid_met2, CLOUD_FRAC_EC_VNAM,strt4d,cnt4d)
!
        ran_cloud(:,:,:) = 0.0d0
!
      else
!
        call Ncrd_4d (humidity,  ncid_met2, HUMIDITY_VNAM, strt4d,cnt4d)
        call Ncrd_4d (max_cloud, ncid_met2, MAX_CLOUD_VNAM,strt4d,cnt4d)
        call Ncrd_4d (ran_cloud, ncid_met2, RAN_CLOUD_VNAM,strt4d,cnt4d)
!
        if ((metdata_name_org  (1:4) == 'GISS') .and.  &
     &      (metdata_name_model(1:6) == '2prime')) then
!
!         --------------------------------------------------------
!         Correct the humidity units.  The max_cloud and ran_cloud
!         variables are flawed in the GISS 2prime data.  It was
!         determined that it is best to just set them to zero.
!         --------------------------------------------------------
!
          humidity (:,:,:) = humidity(:,:,:) * GPKG
          max_cloud(:,:,:) = 0.0d0
          ran_cloud(:,:,:) = 0.0d0
        endif
!
      endif
!
!
!     -------------------------------------------------------------------
!     Read in optical thickness from any meteorlogical data that includes
!     it.  At this point we have optical thickness from GISS and from
!     MATCH3.  The MATCH data needs to be multiplied by cloud fraction
!     to generate a grid box average value.
!     -------------------------------------------------------------------
!
      if (metdata_name_org(1:4) == 'GISS') then
!
        call Ncrd_4d  &
     &    (tau_cloud, ncid_met2, TAU_CLOUD_GISS_VNAM, strt4d, cnt4d)
!
! Bryan Duncan 6-2005
! Set the 1st layer of tau_cloud to 0. for GISS
            tau_cloud(:,:,1) = 0.
!
! Bryan Duncan 6-2005
! Sensitivity test : temporary coding
! This is used for sensitivity test related to cvs tag GMI_GISS_PHASEII
!        tau_cloud(:,:,:) = tau_cloud(:,:,:)/2.
!
!
      else if (metdata_name_org  (1:4) == 'GMAO') then
         call Ncrd_4d  &
              &          (tauice_cloud, ncid_met2, TAU_CLDIC_GEO_VNAM, strt4d,cnt4d)
!
         call Ncrd_4d  &
              &          (tauwater_cloud, ncid_met2, TAU_CLDWA_GEO_VNAM, strt4d,cnt4d)
         ! Needed for fastJX? Megan Damon
         taucli(:,:,:) = tauice_cloud(:,:,:)
         tauclw(:,:,:) = tauwater_cloud(:,:,:)
!
!
         ! Megan Damon Sept 2010
         ! CCM support added for ARRA Aviations work
         if (metdata_name_model(1:5) == 'GEOS4' .or. &
              &     metdata_name_model(1:8) == 'GEOS5CCM') then
!
            tau_cloud(:,:,:) = tauice_cloud(:,:,:) + tauwater_cloud(:,:,:)
!
         else if (metdata_name_model(1:5) == 'GEOS5') then
            call Ncrd_4d (tau_cloud, ncid_met2, OPTDEPTH_VNAM, strt4d,cnt4d)
         endif
!
         if (metdata_name_model(1:5) == 'GEOS5') then
            ! Inverting vertical levels
            tempVar  (:,:,k1:k2) = tau_cloud(:,:,k2:k1:-1)
            tau_cloud(:,:,:) = tempVar  (:,:,:)
            !.sds.. since these are in-cloud opt depths, we need to mult by
            !  cloud fraction for linear overlap
            !    tau_cloud(:,:,:) = tau_cloud(:,:,:)*max_cloud(:,:,:)
            ! Since these are in-cloud opt depths, we need to mult by
            !        cloud fraction**1.5 for approx random overlap (Hongyu Liu)(hyl, bmy, 10/24/08)
            tau_cloud(:,:,:) = tau_cloud(:,:,:)*max_cloud(:,:,:)**1.5d0
!
            ! Jules (01/09/2012): This was added for FastJX 6.5
            IF (metdata_name_model(1:10) == "GEOS5MERRA") THEN
               call ncrd_4d (taucli, ncid_met2, TAUCLI_VNAM, strt4d, cnt4d)
               tempVar(:,:,k1:k2) = taucli(:,:,k2:k1:-1)
               taucli(:,:,k1:k2)  = tempVar(:,:,k1:k2)*max_cloud(:,:,:)**1.5d0
!
               call ncrd_4d (tauclw, ncid_met2, TAUCLW_VNAM, strt4d, cnt4d)
               tempVar(:,:,k1:k2) = tauclw(:,:,k2:k1:-1)
               tauclw(:,:,k1:k2)  = tempVar(:,:,k1:k2)*max_cloud(:,:,:)**1.5d0
!
               tau_cloud(:,:,:) = taucli(:,:,:) + tauclw(:,:,:)
            END IF
!
         end if
!
      else if (metdata_name_model(1:5) == 'MATCH') then
!
        call Ncrd_4d  &
     &    (tau_cloud, ncid_met2, TAU_CLOUD_NCAR_VNAM, strt4d, cnt4d)
!
        tau_cloud(:,:,:) =  &
     &    tau_cloud(:,:,:) *  &
     &    (max_cloud(:,:,:) + ran_cloud(:,:,:))
!
!
      else if (metdata_name_org(1:2) == 'EC') then
!
        !... need cloud optical thickness - have cloud liquid/ice water content in grid box
        call Ncrd_4d (clwc, ncid_met2,TAU_CLDWA_EC_VNAM,strt4d,cnt4d)
        call Ncrd_4d (ciwc, ncid_met2,TAU_CLDIC_EC_VNAM,strt4d,cnt4d)
!
        call Ncrd_3d (psfc,  ncid_met2,'psf',strt3d,cnt3d)
        strt1d(:) = (/ 1 /)
        cnt1d(:) = (/ k2-k1+2 /)
        call Ncrd_1d (ec_ai, ncid_met2,'ai',strt1d,cnt1d)
        call Ncrd_1d (ec_bi, ncid_met2,'bi',strt1d,cnt1d)
!
!-----------------cloud visible opacity used in fast-JX-----------------
!
! tau=(q/(4/3)*reff(microns)*1e-6*1000 kg/m3)*
!     (lwc(kg/kg)/cfr)*(delp(hPa)*100/9.8(m/s2))
        do ik = k1,k2
!.sds... don't think I need this, but what is clwc/ciwc really - in cloud or in box - I assume in box
!          if (CLDFR(L) .gt. 0.d0)  then
!            LWCFR(L) = clwc(L)/max_cloud(L)
!            IWCFR(L) = ciwc(L)/max_cloud(L)
!          endif
!
           tauwater_cloud(:,:,ik) = (QEXT(ik) / ( (4.d0/3.d0) * REFF(ik)*1.d-3 )) *  &
     &        clwc(:,:,ik) * ( ((psfc(:,:)*ec_bi(ik)+ec_ai(ik)) -  &
     &        (psfc(:,:)*ec_bi(ik+1)+ec_ai(ik+1)) ) * (100.d0 / GMI_G) )
!
! Calculate beta(ext) and tau for ice water clouds using
! Sun and Shine expression for beta
! currently no correction for small (<20 micron) crystals
           tauice_cloud(:,:,ik) = ( 1.d3 * ciwc(:,:,ik) /   &
     &        (3.06d-2 + 2.5481d-1*1.d3*ciwc(:,:,ik) )) *  &
     &        (7.d0 * dlog((psfc(:,:)*ec_bi(ik)+ec_ai(ik))/(psfc(:,:)*ec_bi(ik+1)+ec_ai(ik+1)) ) )
        enddo
!
!  Add together taus
        tau_cloud(:,:,:) = (tauice_cloud(:,:,:) + tauwater_cloud(:,:,:)) *  &
     &    (max_cloud(:,:,:) + ran_cloud(:,:,:))
!
!
!!.sds... convert cloud liquid/ice water content from kg/kg per box to g/m^3
!        clwc(:,:,:) = clwc(:,:,:) * GPKG *  mass(:,:,:) / (mcor(:,:,:)*box height)
!        ciwc(:,:,:) = ciwc(:,:,:) * GPKG *  mass(:,:,:) / (mcor(:,:,:)*box height)
!
      endif
!
!
      if (Ncdoes_Var_Exist (ncid_met2, KZZ_VNAM)) then
!
!       ------------------------------------------------------------------
!       Read kzz if it is in the data file.  Also shift the kzz's down one
!       level, and make sure that kzz is positive (some negative values
!       were found in some DAO data).
!       ------------------------------------------------------------------
!
        call Ncrd_4d (kzz, ncid_met2, KZZ_VNAM, strt4d, cnt4d)
!
!.sds for GEOS5, which has all edge values of kzz and is top-down in file, need
!.sds  to make sure we preserve the bottom 0.0 value, not the top one
!
        if (metdata_name_model(1:5) == 'GEOS5') then
           tempVar(:,:,k1+1:k2) = kzz(:,:,k2:k1+1:-1)
           tempVar(:,:,k1) = kzz(:,:,k1)
           kzz(:,:,k1:k2) = tempVar(:,:,k1:k2)
        end if
!
        if (metdata_name_org(1:3) == 'DAO' .or.  &
     &      metdata_name_org(1:4) == 'GMAO') then
!
          do ik = k1, k2-1
            kzz(:,:,ik) = kzz(:,:,ik+1)
          end do
!
          kzz(:,:,k2) = 0.0d0
!
        endif
!
        where (kzz(:,:,:) < 0.0d0)
          kzz(:,:,:) = 0.0d0
        end where
!
      else
!
        kzz(:,:,:) = -1.0d50
!
      endif
!
      return
!
      end subroutine Read_Met2
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_Met1_Range
!
! !INTERFACE:
!
! !DESCRIPTION:
!   This routine does a gross check on the met1 data (ps, u, v, & kel)
!   values.
!
! ARGUMENTS
!   metdata_name_org : first part of metdata_name, e.g., "NCAR"
!   ps  : array of pressure    data (mb)
!   u   : array of u wind      data (m/s)
!   v   : array of v wind      data (m/s)
!   kel : array of temperature data (degK)
!
!-----------------------------------------------------------------------------
!
      subroutine Check_Met1_Range (metdata_name_org, ps, u, v, kel, &
     &   procID, i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, ilo_gl, ihi_gl, &
     &   julo_gl, jvlo_gl, jhi_gl, k1, k2)
!
      use GmiCheckRange_mod, only : checkRange2d, checkRange3d
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: procID
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl
      integer, intent(in) :: ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl, k1, k2
      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
      real*8  :: ps (ilo_gl:ihi_gl, julo_gl:jhi_gl)
      real*8  :: u  (ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)
      real*8  :: v  (ilo_gl:ihi_gl, jvlo_gl:jhi_gl, k1:k2)
      real*8  :: kel(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)
!EOP
!-----------------------------------------------------------------------------
!BOC
      call checkRange2d  &
     &  ('ps', procID, i1_gl, i2_gl, ju1_gl, j2_gl, ps(i1_gl:i2_gl,ju1_gl:j2_gl),  &
     &   0.0d0, 5000.0d0)
!
!
      if (metdata_name_org(1:4) /= 'GISS' .and. metdata_name_org(1:2) /= 'EC') then
!
        call checkRange3d  &
     &    ('u', procID, i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2,  &
     &     u(i1_gl:i2_gl,ju1_gl:j2_gl,:), -400.0d0, 400.0d0)
!
        call checkRange3d  &
     &    ('v', procID, i1_gl, i2_gl, jv1_gl, j2_gl, k1, k2,  &
     &     v(i1_gl:i2_gl,jv1_gl:j2_gl,:), -400.0d0, 400.0d0)
!
      else
!
!       --------------------------------------------------------------
!       Note:  GISS met fields contain mass fluxes in units of 100N/s
!              and EC met fields contain mass fluxes in units of kg/s.
!       --------------------------------------------------------------
!
        call checkRange3d  &
     &    ('u', procID, i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2,  &
     &     u(i1_gl:i2_gl,ju1_gl:j2_gl,:), -4d10, 4d10)
!
        call checkRange3d  &
     &    ('v', procID, i1_gl, i2_gl, jv1_gl, j2_gl, k1, k2,  &
     &     v(i1_gl:i2_gl,jv1_gl:j2_gl,:), -4d10, 4d10)
!
      end if
!
      call checkRange3d  &
     &  ('kel', procID, i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2,  &
     &   kel(i1_gl:i2_gl,ju1_gl:j2_gl,:), 100.0d0, 400.0d0)
!
      return
!
      end subroutine Check_Met1_Range
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Next_Met_File
!
! !INTERFACE:
!
! !DESCRIPTION:
!   This routine gets the name of the next met data file and the number
!   of records in it.
!
! ARGUMENTS
!   do_cycle_met     : cycle met data?
!   out_of_data      : out of met data?
!   met_infile_num   : index of current met input file
!   met_num_infiles  : total number of  met input files
!   mnum_recs        : number of NetCDF met records in file
!   mrnum_in         : next NetCDF met record to read
!   ncid_met         : NetCDF met input file id
!   met_infile_names : met input file names
!
!-----------------------------------------------------------------------------
!
      subroutine Get_Next_Met_File  &
     &  (do_cycle_met, out_of_data, met_infile_num, met_num_infiles,  &
     &   mnum_recs, mrnum_in, ncid_met, met_infile_names, &
     &   pr_diag, procID)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      logical :: do_cycle_met
      logical :: out_of_data
      integer :: met_infile_num
      integer :: met_num_infiles
      integer :: mnum_recs
      integer :: mrnum_in
      integer :: ncid_met
      character (len=MAX_LENGTH_FILE_NAME) :: met_infile_names(*)
      character(len=ESMF_MAXSTR) :: IAm
!
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = "Get_Next_Met_File"
!
      if (pr_diag) Write (6,*) IAm, ' called by ', procID
!
      out_of_data = .false.
!
!
      call Nccl (ncid_met)
!
      met_infile_num = met_infile_num + 1
!
!
      if (met_infile_num > met_num_infiles) then
!
        if (do_cycle_met) then
!
          met_infile_num = 1
!
          Write (6,*)
          Write (6,*) 'Cycling on met data.'
          Write (6,*)
!
        else
!
          out_of_data = .true.
!
        end if
!
      end if
!
!
      if (.not. out_of_data) then
!
        call Ncop_Rd (ncid_met, met_infile_names(met_infile_num))
!
        call Ncget_Unlim_Dimlen (ncid_met, mnum_recs)
!
        mrnum_in = 1
!
      else
!
        Write (6,*)
        Write (6,*) 'Out of met data in Get_Next_Met_File.'
        Write (6,*)
!
      end if
!
      return
!
      end subroutine Get_Next_Met_File
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Var_Adv_Tstp
!
! !INTERFACE:
!
! !DESCRIPTION:
!   This routine calculates the number of advection time steps.
!
!   Note that delpm inside tpcore is calculated using pctm2, not pmet2, so
!   we need to use pctm2 to calculate the delpm for calculating the Courant
!   numbers.  The difference should be small, but if we aren't this careful,
!   we could violate the Courant limit inside tpcore, and hence cause an
!   array out of bounds error.
!
! ARGUMENTS
!   num_adv_time_stps : number of advection time steps
!   cose      : cosine of grid box edges
!   dap       : pressure difference across layer from (ai * pt) term (mb)
!   dbk       : difference in bi across layer - the dSigma term
!   geofac_pc : special geometrical factor (geofac) for Polar cap
!   geofac    : geometrical factor for meridional advection; geofac uses
!               correct spherical geometry, and replaces acosp as the
!               meridional geometrical factor in tpcore
!   pctm1     : CTM surface pressure at t1     (mb)
!   pctm2     : CTM surface pressure at t1+tdt (mb)
!   xmass     : horizontal mass flux in E-W direction (mb)
!   ymass     : horizontal mass flux in N-S direction (mb)
!   zmass     : vertical   mass flux                  (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Calc_Var_Adv_Tstp  &
     &  (num_adv_time_steps, cose, dap, dbk, geofac_pc, geofac,  &
     &   pctm1, pctm2, xmass, ymass, zmass, &
     &   commu_npole, commu_spole, pr_diag, procID, numLonDomains, &
     &   gmi_nborder, j1p, j2p, i1, i2, ju1, j2, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, k1, k2, ilo, ihi, julo, jhi)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in)  :: pr_diag
      integer, intent(in)  :: procID, numLonDomains
      integer, intent(in)  :: gmi_nborder, j1p, j2p, i1, i2, ju1, j2, i1_gl, i2_gl
      integer, intent(in)  :: ju1_gl, j2_gl, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in)  :: commu_npole, commu_spole
      real*8,  intent(in)  :: cose(ju1_gl:j2_gl+1)
      real*8,  intent(in)  :: dap (k1:k2)
      real*8,  intent(in)  :: dbk (k1:k2)
      real*8,  intent(in)  :: geofac_pc
      real*8,  intent(in)  :: geofac(ju1_gl:j2_gl)
      real*8,  intent(in)  :: pctm1(ilo:ihi, julo:jhi)
      real*8,  intent(in)  :: pctm2(ilo:ihi, julo:jhi)
      real*8,  intent(in)  :: xmass(ilo:ihi, julo:jhi, k1:k2)
      real*8,  intent(in)  :: ymass(ilo:ihi, julo:jhi, k1:k2)
      real*8,  intent(in)  :: zmass(i1_gl:i2_gl, ju1_gl:j2_gl, k1:k2)
!
      integer, intent(out) :: num_adv_time_steps
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik
!
      real*8  :: courantx, couranty
      real*8  :: max_crx1, max_cry1
      real*8  :: max_crx2, max_cry2
      real*8  :: eff_courant
      real*8  :: max_vert_courant_down, max_vert_courant_up
      real*8  :: rgmi_nborder
!
!     -----------------------------------------------------------------------
!     delpm   : pressure thickness, the psudo-density in a hydrostatic system
!               at t1+tdt/2 (approximate) (mb)
!     dps_ctm : CTM surface pressure tendency; sum over vertical of dpi
!               calculated from original mass fluxes (mb)
!     dpi     : divergence at a grid point; used to calculate vertical motion
!               (mb)
!     wz      : large scale mass flux (per time step tdt) in the vertical
!               direction as diagnosed from the hydrostatic relationship (mb)
!     -----------------------------------------------------------------------
!
      real*8  :: delpm(ilo:ihi, julo:jhi, k1:k2)
!
      real*8  :: dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl)
      real*8  :: minpctm(i1_gl:i2_gl, ju1_gl:j2_gl)
!
      real*8  :: dpi    (i1_gl:i2_gl, ju1_gl:j2_gl, k1:k2)
      character(len=ESMF_MAXSTR) :: IAm
!
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = "CalC_Var_Adv_Tstp"
!
      if (pr_diag) Write (6,*) IAm, ' called by ', procID
!
!     --------------------------------------------------------------
!     Calculate Horizontal Courant numbers assuming first pctm1, and
!     then pctm2, to find the largest values.
!     --------------------------------------------------------------
!
!     ===============
      call Calc_Delpm  &
!     ===============
     &  (dap, dbk, pctm1, pctm1, delpm, &
     &   pr_diag, procID, k1, k2, ilo, ihi, julo, jhi)
!
!     =====================
      call Calc_Courant_Met  &
!     =====================
     &  (cose, delpm, xmass, ymass, max_crx1, max_cry1, &
     &   i1, i2, ju1, j2, k1, k2, j1p, j2p, ju1_gl, j2_gl, &
     &   ilo, ihi, julo, jhi, pr_diag, procID)
!
!
!     ===============
      call Calc_Delpm  &
!     ===============
     &  (dap, dbk, pctm2, pctm2, delpm, &
     &   pr_diag, procID, k1, k2, ilo, ihi, julo, jhi)
!
!     =====================
      call Calc_Courant_Met  &
!     =====================
     &  (cose, delpm, xmass, ymass, max_crx2, max_cry2, &
     &   i1, i2, ju1, j2, k1, k2, j1p, j2p, ju1_gl, j2_gl, &
     &   ilo, ihi, julo, jhi, pr_diag, procID)
!
!
      courantx  = Max (max_crx1, max_crx2)
      couranty  = Max (max_cry1, max_cry2)
!
      if (pr_diag) then
        Write (6,*) "courantx = ", courantx
        Write (6,*) "couranty = ", couranty
      end if
!
!
      rgmi_nborder = gmi_nborder
!
      eff_courant  = Max (courantx/rgmi_nborder, couranty)
!
!
!     ----------------------------------
!     Calculate vertical Courant values.
!     ----------------------------------
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
      dps_ctm(:,:) = Sum (dpi(:,:,:), dim=3)
!
!!     ========================
!      call calcVerticalMassFlux  &
!!     ========================
!     &  (dbk, dps_ctm, dpi, wz, &
!     &   pr_diag, procID, i1_gl, i2_gl, ju1_gl, j2_gl, k1, k2)
!
!
      minpctm(:,:) =  &
     &  Min (pctm1(i1_gl:i2_gl,ju1_gl:j2_gl),  &
     &       pctm2(i1_gl:i2_gl,ju1_gl:j2_gl))
!
!
      max_vert_courant_up   = 0.0d0
!
      do ik = k1, k2
        do ij = ju1_gl, j2_gl
          do il = i1_gl, i2_gl
!
            max_vert_courant_up =  &
     &        Max (max_vert_courant_up,  &
     &             Abs (zmass(il,ij,ik)   /  &
     &                  (dap(ik) + (dbk(ik) * minpctm(il,ij)))))
!
          end do
        end do
      end do
!
!
      max_vert_courant_down = 0.0d0
!
      do ik = k1, k2-1
        do ij = ju1_gl, j2_gl
          do il = i1_gl, i2_gl
!
            max_vert_courant_down =  &
     &        Max (max_vert_courant_down,  &
     &             Abs (zmass(il,ij,ik+1) /  &
     &                  (dap(ik) + (dbk(ik) * minpctm(il,ij)))))
!
          end do
        end do
      end do
!
!
      if (pr_diag) then
        Write (6,*) "max_vert_courant_up   = ", max_vert_courant_up
        Write (6,*) "max_vert_courant_down = ", max_vert_courant_down
      end if
!
!
      eff_courant =  &
     &   Max (eff_courant, max_vert_courant_up, max_vert_courant_down)
!
!
!     ---------------------------------------------------------------------
!     Note that 1.01 has been used below, instead of 1.00, in order to
!     protect against any numerical precision events whereby an eff_courant
!     of 0.99999 gets calculated as 1.00001 inside tpcore.  It is unlikely
!     that any problem would ever occur, but it is better to be safe than
!     sorry.
!     ---------------------------------------------------------------------
!
      num_adv_time_steps = eff_courant + 1.01d0
!
      return
!
      end subroutine Calc_Var_Adv_Tstp
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Courant_Met
!
! !INTERFACE:
!
! !DESCRIPTION:
!   This routine calculates courant numbers from the horizontal mass fluxes.
!
! ARGUMENTS
!   cose  : cosine of grid box edges
!   delpm : pressure thickness, the psudo-density in a hydrostatic system
!           at t1+tdt/2 (approximate) (mb)
!   xmass : horizontal mass flux in E-W direction (mb)
!   ymass : horizontal mass flux in N-S direction (mb)
!   max_crx : Maximum Courant number in E-W direction
!   max_cry : Maximum Courant number in N-S direction
!
!-----------------------------------------------------------------------------
!
      subroutine Calc_Courant_Met  &
     &  (cose, delpm, xmass, ymass, max_crx, max_cry, &
     &   i1, i2, ju1, j2, k1, k2, j1p, j2p, ju1_gl, j2_gl, &
     &   ilo, ihi, julo, jhi, pr_diag, procID)
!
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: ilo, ihi, julo, jhi, ju1_gl, j2_gl
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, j1p, j2p
      real*8, intent(in)  :: cose (ju1_gl:j2_gl+1)
      real*8, intent(in)  :: delpm(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: xmass(ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: ymass(ilo:ihi, julo:jhi, k1:k2)
!
      real*8, intent(out) :: max_crx
      real*8, intent(out) :: max_cry
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik
      integer :: jst, jend
!
      real*8  :: tmp_pu
!
      real*8  :: tmp_crx
      real*8  :: tmp_cry
      character(len=ESMF_MAXSTR) :: IAm
!
!EOP
!-----------------------------------------------------------------------------
!BOC
      Iam = "Calc_Courant_Met"
!
      if (pr_diag) Write (6,*) IAm, ' called by ', procID
!
      max_crx = 0.0d0
      max_cry = 0.0d0
!
      tmp_pu  = 0.0d0
!
      tmp_crx = 0.0d0
      tmp_cry = 0.0d0
!
!     -----------------------------------
!     Calculate E-W horizontal mass flux.
!     -----------------------------------
!
      call Deter_Jrange (ju1, j2, j1p, j2p, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!
      do ik = k1, k2
        do ij = jst, jend
          do il = i1, i2+1
!
            tmp_pu = 0.5d0 * (delpm(il,ij,ik) + delpm(il-1,ij,ik))
!
            tmp_crx = xmass(il,ij,ik) / tmp_pu
!
            max_crx = MAX( ABS(tmp_crx), max_crx )
!
          end do
        end do
      end do
!
!
!     -----------------------------------
!     Calculate N-S horizontal mass flux.
!     -----------------------------------
!
      call Deter_Jrange (ju1, j2+1, j1p, j2p+1, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!
      do ik = k1, k2
        do ij = jst, jend
          do il = i1, i2
!
            tmp_pu = 0.5d0 * (delpm(il,ij,ik) + delpm(il,ij-1,ik))
!
            tmp_cry = ymass(il,ij,ik) / (tmp_pu * cose(ij))
!
            max_cry = MAX( ABS(tmp_cry), max_cry )
!
          end do
        end do
      end do
!
!
      return
!
      end subroutine Calc_Courant_Met
!EOC
!------------------------------------------------------------------------
!
! ROUTINE
!   Check_Vert_Courant
!
! DESCRIPTION
!   This routine checks that the vertical courant condition isn't broken.
!
!   For reasonable met fields the vertical velocity is generally so small
!   that the vertical Courant limit is not a concern.  However, if the winds
!   are really noisy it could be a problem, so it is worth checking.
!
!   NOTE: This routine has been written to go outside Tpcore.  This means
!         that the vertical indexing is from the bottom, and the vertical
!         flux is calculated starting from the bottom.  This also means that
!         positive flux is up, rather than down.  However, this routine
!         should work without change inside Tpcore (where everything is
!         reversed), but the sign and indexing of any problem found will be
!         reported differently.
!
! ARGUMENTS
!   do_reduction : set to False if called on Master;
!                  set to True  if called by Slaves
!   num_adv_time_stps : number of advection time steps
!   geofac_pc    : special geometrical factor (geofac) for Polar cap
!   geofac       : geometrical factor for meridional advection; geofac uses
!                  correct spherical geometry, and replaces acosp as the
!                  meridional geometrical factor in tpcore
!   dap   : pressure difference across layer from (ai * pt) term (mb)
!   dbk   : difference in bi across layer - the dSigma term
!   pctm1 : CTM surface pressure at t1     (mb)
!   pctm2 : CTM surface pressure at t1+tdt (mb)
!   xmass : horizontal mass flux in E-W direction (mb)
!   ymass : horizontal mass flux in N-S direction (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Check_Vert_Courant  &
     &  (do_reduction, num_adv_time_steps, geofac_pc, geofac,  &
     &   dap, dbk, pctm1, pctm2, xmass, ymass, &
     &   pr_diag, procID, numLonDomains, j1p, j2p, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole)
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
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      logical, intent(in) :: do_reduction
      integer, intent(in) :: num_adv_time_steps
      real*8,  intent(in) :: geofac_pc
      real*8,  intent(in) :: geofac(ju1_gl:j2_gl)
      real*8,  intent(in) :: dap(k1:k2)
      real*8,  intent(in) :: dbk(k1:k2)
      real*8,  intent(in) :: pctm1(ilo:ihi, julo:jhi)
      real*8,  intent(in) :: pctm2(ilo:ihi, julo:jhi)
      real*8,  intent(in) :: xmass(ilo:ihi, julo:jhi, k1:k2)
      real*8,  intent(in) :: ymass(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij, ik
!
      real*8  :: r_num_adv_time_steps
!
      real*8  :: dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl)
!
      real*8  :: dp
      real*8  :: dpi(i1_gl:i2_gl, ju1_gl:j2_gl, k1:k2)
      real*8  :: wz (i1_gl:i2_gl, ju1_gl:j2_gl, k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      r_num_adv_time_steps = num_adv_time_steps
!
!
!     ====================
      call calcDivergence  &
!     ====================
     &  (do_reduction, geofac_pc, geofac, dpi, xmass, ymass, &
     &   pr_diag, procID, numLonDomains, j1p, j2p, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole)
!
      dps_ctm(:,:) = Sum (dpi(:,:,:), dim=3)
!
!
!     ========================
      call calcVerticalMassFlux  &
!     ========================
     &  (dbk, dps_ctm, dpi, wz, &
     &   pr_diag, procID, i1, i2, ju1, j2, k1, k2)
!
!
      wz(:,:,:) = wz(:,:,:) / r_num_adv_time_steps
!
!
      do ik = k1, k2
        do ij = ju1_gl, j2_gl
          do il = i1_gl, i2_gl
!
            dp = dap(ik) +  &
     &           dbk(ik) * Min (pctm1(il,ij), pctm2(il,ij))
!
            if (Abs (wz(il,ij,ik)) >= dp) then
!
              Write (6,*) "Vertical courant number exceeded:"
              Write (6,*) "  At (", il, ij, ik, ")"
              Write (6,*) "  wz(k) = ", wz(il,ij,ik)
              Write (6,*) "  dp(k) = ", dp
!
            end if
!
            if (ik /= k2) then  ! wz not defined for k2+1
              if (Abs (wz(il,ij,ik+1)) >= dp) then
!
                Write (6,*) "Vertical courant number exceeded:"
                Write (6,*) "  At (", il, ij, ik, ")"
                Write (6,*) "  wz(k+1) = ", wz(il,ij,ik+1)
                Write (6,*) "  dp(k)   = ", dp
!
              end if
            end if
!
          end do
        end do
      end do
!
      return
!
      end subroutine Check_Vert_Courant
!
!!-----------------------------------------------------------------------------
!!
!! ROUTINE
!!   Finish_Met
!!
!! DESCRIPTION
!!   This routine does some final calculations/conversions of humidity and pbl,
!!   if necessary.
!!
!! ARGUMENTS
!!   metdata_name_org   : first  part of metdata_name, e.g., "NCAR"
!!   metdata_name_model : second part of metdata_name, e.g., "MATCH"
!!   new_met_rec   : new met record?
!!   t_cloud_ice   : temperature for ice formation in clouds (degK)
!!   kel           : current temperature (degK)
!!   del_kel       : temperature at end of met time period   (degK)
!!   psx           : surface pressure at beginning of met time period  (mb)
!!   del_psx       : surface pressure at end       of met time period  (mb)
!!   press3e       : atmospheric pressure at the edge of each grid box (mb)
!!   pbl           : planetary boundary layer depth (in:mb, out:m)
!!   surf_air_temp : surface air temperature (degK)
!!   humidity      : specific humidity (g/kg)
!!   relhum        : relative humidity (fraction)
!!
!!-----------------------------------------------------------------------------
!
!      subroutine Finish_Met  &
!     &  (metdata_name_org, metdata_name_model, new_met_rec,  &
!     &   t_cloud_ice, kel, del_kel, psx, del_psx, press3e,  &
!     &   pbl, surf_air_temp, humidity, relhum,  &
!     &   del_pbl, del_surf_air_temp, del_humidity, del_relhum, &
!     &   pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
!      implicit none
!
!
!!     ----------------------
!!     Argument declarations.
!!     ----------------------
!
!      logical, intent(in) :: pr_diag
!      integer, intent(in) :: procID
!      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
!      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
!      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_model
!      logical :: new_met_rec
!      real*8  :: t_cloud_ice
!      real*8  :: kel    (ilo:ihi, julo:jhi, k1:k2)
!      real*8  :: del_kel(ilo:ihi, julo:jhi, k1:k2)
!      real*8  :: psx    (ilo:ihi, julo:jhi)
!      real*8  :: del_psx(ilo:ihi, julo:jhi)
!      real*8  :: press3e(ilo:ihi, julo:jhi, k1-1:k2)
!      real*8  :: pbl          (i1:i2, ju1:j2)
!      real*8  :: surf_air_temp(i1:i2, ju1:j2)
!      real*8  :: humidity     (i1:i2, ju1:j2, k1:k2)
!      real*8  :: relhum       (i1:i2, ju1:j2, k1:k2)
!      real*8  :: del_pbl          (i1:i2, ju1:j2)
!      real*8  :: del_surf_air_temp(i1:i2, ju1:j2)
!      real*8  :: del_humidity     (i1:i2, ju1:j2, k1:k2)
!      real*8  :: del_relhum       (i1:i2, ju1:j2, k1:k2)
!
!!     ----------------------
!!     Variable declarations.
!!     ----------------------
!
!      logical, save :: first_pbl = .true.
!
!
!!     ----------------
!!     Begin execution.
!!     ----------------
!
!      if (pr_diag) then
!        Write (6,*) 'Finish_Met called by ', procID
!      end if
!
!
!      if (new_met_rec) then
!
!        if ((metdata_name_org  (1:4) == 'NCAR') .and.  &
!     &      (metdata_name_model(1:4) == 'CCM3')) then
!
!!         ==================
!          call Calc_Humidity  &
!!         ==================
!     &      (t_cloud_ice, kel, press3e, relhum, humidity, &
!     &       pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
!!         ==================
!          call Calc_Humidity  &
!!         ==================
!     &      (t_cloud_ice, del_kel, press3e, del_relhum, del_humidity, &
!     &       pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
!        end if
!
!        if (metdata_name_org(1:4) /= 'GISS' .and.  &
!     &      (metdata_name_org  (1:4) /= 'GMAO' .and.  &
!     &       metdata_name_model(1:5) /= 'GEOS4')) then
!
!          if (first_pbl) then
!
!            first_pbl = .false.
!
!!           ================
!            call Convert_Pbl  &
!!           ================
!     &        (psx, surf_air_temp, humidity, pbl, &
!     &         pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
!          end if
!
!!         ================
!          call Convert_Pbl  &
!!         ================
!     &      (del_psx, del_surf_air_temp, del_humidity, del_pbl, &
!     &       pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
!        end if
!
!      end if
!
!      return
!
!      end subroutine Finish_Met
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calc_Humidity
!
! DESCRIPTION
!   This routine calculates specific humidity (mb) using a formula from
!   Seinfeld+Pandis (p. 781).
!
! ARGUMENTS
!   t_cloud_ice : temperature of ice formation in clouds (degK)
!   kel         : temperature (degK)
!   press3e     : atmospheric pressure at the edge of each grid box (mb)
!   relhum      : relative humidity (from CCM3)
!   humidity    : specific humidity (in:mb, out:g/kg)
!
!-----------------------------------------------------------------------------
!
      subroutine Calc_Humidity  &
     &  (t_cloud_ice, kel, press3e, relhum, humidity, &
     &   pr_diag, procID, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      real*8  :: t_cloud_ice
      real*8  :: kel     (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: press3e (ilo:ihi, julo:jhi, k1-1:k2)
      real*8  :: relhum  (i1:i2, ju1:j2, k1:k2)
      real*8  :: humidity(i1:i2, ju1:j2, k1:k2)
!
!
!     ----------------------
!     Parameter declarations.
!     ----------------------
!
      real*8, parameter :: A0W = 6.107799961d+00
      real*8, parameter :: A1W = 4.436518521d-01
      real*8, parameter :: A2W = 1.428945805d-02
      real*8, parameter :: A3W = 2.650648471d-04
      real*8, parameter :: A4W = 3.031240396d-06
      real*8, parameter :: A5W = 2.034080948d-08
      real*8, parameter :: A6W = 6.136820929d-11
!
      real*8, parameter :: A0I = 6.109177956d+00
      real*8, parameter :: A1I = 5.034698970d-01
      real*8, parameter :: A2I = 1.886013408d-02
      real*8, parameter :: A3I = 4.176223716d-04
      real*8, parameter :: A4I = 5.824720280d-06
      real*8, parameter :: A5I = 4.838803174d-08
      real*8, parameter :: A6I = 1.838826904d-10
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Calc_Humidity called by ', procID
      end if
!
!
      where (kel(i1:i2,ju1:j2,:) > t_cloud_ice)
!
        humidity(:,:,:) =  &
     &    relhum(:,:,:) *  &
     &    (A0W +  &
     &     A1W * (kel(i1:i2,ju1:j2,:) + ABS_ZERO)    +  &
     &     A2W * (kel(i1:i2,ju1:j2,:) + ABS_ZERO)**2 +  &
     &     A3W * (kel(i1:i2,ju1:j2,:) + ABS_ZERO)**3 +  &
     &     A4W * (kel(i1:i2,ju1:j2,:) + ABS_ZERO)**4 +  &
     &     A5W * (kel(i1:i2,ju1:j2,:) + ABS_ZERO)**5 +  &
     &     A6W * (kel(i1:i2,ju1:j2,:) + ABS_ZERO)**6)
!
      elsewhere
!
        humidity(:,:,:) =  &
     &    relhum(:,:,:) *  &
     &    (A0I +  &
     &     A1I * (kel(i1:i2,ju1:j2,:) + ABS_ZERO)    +  &
     &     A2I * (kel(i1:i2,ju1:j2,:) + ABS_ZERO)**2 +  &
     &     A3I * (kel(i1:i2,ju1:j2,:) + ABS_ZERO)**3 +  &
     &     A4I * (kel(i1:i2,ju1:j2,:) + ABS_ZERO)**4 +  &
     &     A5I * (kel(i1:i2,ju1:j2,:) + ABS_ZERO)**5 +  &
     &     A6I * (kel(i1:i2,ju1:j2,:) + ABS_ZERO)**6)
!
      end where
!
!
!     ---------------------------------
!     Convert humidity from mb to g/kg.
!     ---------------------------------
!
      humidity(:,:,:) =  &
     &  humidity(:,:,:) * GPKG /  &
     &  (press3e(i1:i2,ju1:j2,k1-1:k2-1) -  &
     &   press3e(i1:i2,ju1:j2,k1:k2))
!
!
      return
!
      end subroutine Calc_Humidity
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Convert_Pbl
!
! DESCRIPTION
!   This routine converts planetary boundary layer depth (pbl) from units
!   of millibars to meters.
!
! ARGUMENTS
!   psf           : surface pressure  (hPa)
!   surf_air_temp : surface air temperature (degK)
!   humidity      : specific humidity (g/kg)
!   pbl           : planetary boundary layer depth (in:mb, out:m)
!
!-----------------------------------------------------------------------------
!
      subroutine Convert_Pbl  &
     &  (psf, surf_air_temp, humidity, pbl, &
     &   pr_diag, procID, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      real*8  :: psf          (ilo:ihi, julo:jhi)
      real*8  :: surf_air_temp(i1:i2, ju1:j2)
      real*8  :: humidity     (i1:i2, ju1:j2, k1:k2)
      real*8  :: pbl          (i1:i2, ju1:j2)
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Convert_Pbl called by ', procID
      end if
!
      pbl(:,:) =  &
     &  29.3d0 * surf_air_temp(:,:) *  &
     &  (1.0d0 +  &
     &   (0.61d-3 * humidity(:,:,1))) *  &
     &  Log (psf(i1:i2,ju1:j2) / (psf(i1:i2,ju1:j2) - pbl(:,:)))
!
      return
!
      end subroutine Convert_Pbl
!EOC
!------------------------------------------------------------------------------
!
      end module GmiUtilsMetFields_mod
!
