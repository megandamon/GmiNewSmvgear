!-------------------------------------------------------------------------
 module GmiTracerMethod_mod
!
!  use GmiGrid_mod, only : t_gmiGrid, Get_i1, Get_i2, Get_ju1, Get_j2,      &
!         Get_k1, Get_k2, Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl,      &
!         Get_ilo, Get_ihi, Get_julo, Get_jhi, Get_ilo_gl, Get_ihi_gl,      &
!         Get_julo_gl, Get_jhi_gl, Get_ilong, Get_ilat, Get_ivert,          &
!         Get_numSpecies
!
  public  :: t_SF6
  type t_SF6
   character (len=132) :: SF6_infile_name
   real*8 , pointer :: j_SF6 (:,:,:) => null()        ! used in the SF6 module
   real*8 , pointer :: k_SF6 (:,:,:) => null()        ! used in the SF6 module
   real*8 , pointer :: AD_SF6 (:,:,:) => null()       ! used in the SF6 module
   real*8 , pointer :: O1D_SF6 (:,:,:) => null()      ! used in the SF6 module
  end type t_SF6
!
  public  :: t_Linoz
  type t_Linoz
   character (len=132) :: linoz_infile_name
   real*8 , pointer :: O3_linoz (:,:) => null()          ! used in the linoz module
   real*8 , pointer :: T_linoz (:,:) => null()           ! used in the linoz module
   real*8 , pointer :: colO3_linoz (:,:) => null()       ! used in the linoz module
   real*8 , pointer :: PmL_linoz (:,:) => null()         ! used in the linoz module
   real*8 , pointer :: dPmLdO3_linoz (:,:) => null()     ! used in the linoz module
   real*8 , pointer :: dPmLdT_linoz (:,:) => null()      ! used in the linoz module
   real*8 , pointer :: dPmLdcolO3_linoz (:,:) => null()  ! used in the linoz module
   real*8 , pointer :: Cariolle_Loss_linoz (:) => null() ! used in the linoz module
   real*8 :: O3SSFC_linoz
   real*8 :: O3STAU_linoz
   real*8 :: O3PAUZ_linoz
   real*8 :: CLLOAD_linoz
   real*8 :: TPSC_linoz
  end type t_Linoz
!
  public  :: t_stratO3
  type t_stratO3
   character (len=132) :: SO3daily_infile_name      ! used in the Strat_O3 module
   character (len=132) :: SO3monthly_infile_name    ! used in the Strat_O3 module
   real*8 , pointer :: daily_o3    (:,:,:) => null()  ! used in the Strat_O3 module
   real*8 , pointer :: monthly_o1d (:,:,:) => null()  ! used in the Strat_O3 module
   real*8 , pointer :: monthly_oh  (:,:,:) => null()  ! used in the Strat_O3 module
   real*8 , pointer :: monthly_ho2 (:,:,:) => null()  ! used in the Strat_O3 module
  end type t_stratO3
!
!------------------------------------------------------------------
!  Arrays and dimensions needed for calculation of the linearized
!   ozone, "linoz".
!------------------------------------------------------------------
!
  CONTAINS
!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update_Radon_Lead
!
! !INTERFACE:
!
      subroutine Update_Radon_Lead  &
     &  (lwi_flags, mcor, mass, surf_air_temp, concentration, IRn, IPb, &
     &   pr_diag, loc_proc, &
     &   tdt, i1, i2, ju1, j2, k1, k2, num_species)
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!
      implicit none
!
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
!
! !INPUT PARAMETERS:
      logical, intent(in   ) :: pr_diag
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: num_species
      integer, intent(in   ) :: IRn, IPb
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2
      real*8 , intent(in   ) :: tdt
      integer, intent(in   ) :: lwi_flags     (i1:i2, ju1:j2)
                              ! array of flags that indicate land, water, or ice
      real*8 , intent(in   ) :: mcor          (i1:i2, ju1:j2)
                              ! area of grid box (m^2)
      real*8 , intent(in   ) :: mass          (i1:i2, ju1:j2, k1:k2)
                              ! total mass of the atmosphere within each grid box (kg)
      real*8 , intent(in   ) :: surf_air_temp (i1:i2, ju1:j2)
                              ! surface area temperature (degK)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
                              ! species concentration, known at zone centers (mixing ratio)
!
! !DESCRIPTION:
!   This routine updates the Radon_Lead chemistry.
!
! !REVISION HISTORY:
!   Initial code.
!EOP
!-------------------------------------------------------------------------
!
!     -----------------------
!     Parameter declarations.
!     -----------------------
!
      real*8, parameter :: CMPM2 = CMPM * CMPM
!
!     -----------------------------------------------
!     Radon land scaling factor for cold temperatures
!     (Rind+Lerner, JGR 1996).
!     -----------------------------------------------
!
      real*8, parameter :: COLD_TEMP_LAND     = 273.0d0    ! degK
      real*8, parameter :: COLD_TEMP_LAND_FAC =   0.6d0
!
      real*8, parameter :: EFOLD_TIME_RN222   =   5.5d0    ! days
!
      real*8, parameter :: RN_SOURCE_LAND     =   1.000d0  ! atoms/cm^2/s
      real*8, parameter :: RN_SOURCE_OCEAN    =   0.005d0  ! atoms/cm^2/s
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij
!
      real*8  :: exp_fac
      real*8  :: rn_source
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Update_Radon_Lead called by ', loc_proc
      end if
!
!     ----------------------------------
!     Lead increases due to Radon decay.
!     ----------------------------------
!
      exp_fac = Exp (-tdt / (EFOLD_TIME_RN222 * SECPDY))
!
      concentration(IPb)%pArray3D(:,:,:) =  &
     &  concentration(IPb)%pArray3D(:,:,:) +  &
     &  (1.0d0 - exp_fac) * concentration(IRn)%pArray3D(:,:,:)
!
!     -----------------------------
!     Radon decreases due to decay.
!     -----------------------------
!
      concentration(IRn)%pArray3D(:,:,:) = exp_fac * concentration(IRn)%pArray3D(:,:,:)
!
!
!     ---------------------------------------------------------------------
!     The radon comes from the land and the ocean.  This was originally
!     scaled in the code to match the lwi flags in the various met fields,
!     but this is no longer done.  Instead, if the radon needs to be scaled
!     to get exactly 16 kg/yr, the scaling must be done to the output,
!     after the run is completed.
!     ---------------------------------------------------------------------
!
      do ij = ju1, j2
        do il = i1, i2
!
          if (lwi_flags(il,ij) == 1) then  ! ocean
!
            rn_source = RN_SOURCE_OCEAN
          else if (lwi_flags(il,ij) == 3) then  ! ice
             rn_source = 0.0d0
!
          else if (lwi_flags(il,ij) == 2) then  ! land
!
            rn_source = RN_SOURCE_LAND
!
            if (surf_air_temp(il,ij) < COLD_TEMP_LAND) then
              rn_source = rn_source * COLD_TEMP_LAND_FAC
            end if
!
          end if
!
          if ((lwi_flags(il,ij) == 1) .or.  &
     &        (lwi_flags(il,ij) == 2)) then
!
            concentration(IRn)%pArray3D(il,ij,1) =  &
     &        concentration(IRn)%pArray3D(il,ij,1) +  &
     &        (rn_source * mcor(il,ij) * CMPM2 * tdt) /  &
     &        (mass(il,ij,1) * (AVOGAD * GPKG / MWTAIR))
!
          end if
!
        end do
      end do
!
!
      return
!
      end subroutine Update_Radon_Lead
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update_CH3I
!
! !INTERFACE:
!
      subroutine Update_CH3I  &
     &  (lwi_flags, mcor, mass, concentration, ICH3I, &
     &   pr_diag, loc_proc, tdt, i1, i2, ju1, j2, ilo, ihi, julo, jhi, k1, k2, num_species)
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!
      implicit none
!
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
!
! !INPUT PARAMETERS:
      logical, intent(in   ) :: pr_diag
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: num_species
      integer, intent(in   ) :: ICH3I
      integer, intent(in   ) :: i1, i2, ju1, j2, ilo, ihi, julo, jhi, k1, k2
      real*8 , intent(in   ) :: tdt
      integer, intent(in   ) :: lwi_flags (i1:i2, ju1:j2) ! array of flags that indicate land, water, or ice
      real*8 , intent(in   ) :: mcor (i1:i2, ju1:j2)      ! area of grid box (m^2)
      real*8 , intent(in   ) :: mass (i1:i2, ju1:j2, k1:k2)  ! total mass of the atmosphere within each grid box (kg)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
                              ! species concentration, known at zone centers (mixing ratio)
!
! !DESCRIPTION:
!   This routine updates the CH3I "chemistry".
!
! !REVISION HISTORY:
!   Initial code.
!EOP
!-------------------------------------------------------------------------
!
!     -----------------------
!     Parameter declarations.
!     -----------------------
!
      real*8, parameter :: EFOLD_TIME_CH3I    =   5.0d0    ! days
!
      real*8, parameter :: CH3I_SOURCE_LAND   =   0.000d0  ! molecules/cm^2/s
      real*8, parameter :: CH3I_SOURCE_ICE    =   0.000d0  ! molecules/cm^2/s
      real*8, parameter :: CH3I_SOURCE_OCEAN  =   1.000d0  ! molecules/cm^2/s
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ij
!
      real*8  :: exp_fac
      real*8  :: CH3I_source
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Update_CH3I called by ', loc_proc
      end if
!
!     -----------------------------
!     CH3I decreases due to decay.
!     -----------------------------
!
      exp_fac = Exp (-tdt / (EFOLD_TIME_CH3I * SECPDY))
      concentration(ICH3I)%pArray3D(:,:,:) = exp_fac * concentration(ICH3I)%pArray3D(:,:,:)
!
!
!     ---------------------------------------------------------------------
!     The methyl iodide comes from the ocean.
!     ---------------------------------------------------------------------
!
      do ij = ju1, j2
         do il = i1, i2
!
            if (lwi_flags(il,ij) == 1) then  ! ocean
               CH3I_source = CH3I_SOURCE_OCEAN
            else if (lwi_flags(il,ij) == 3) then  ! ice
               CH3I_source = CH3I_SOURCE_ICE
            else if (lwi_flags(il,ij) == 2) then  ! land
               CH3I_source = CH3I_SOURCE_LAND
            end if
!
            concentration(ICH3I)%pArray3D(il,ij,1) =  &
     &        concentration(ICH3I)%pArray3D(il,ij,1) +  &
     &        (CH3I_source * mcor(il,ij) * CMPM * CMPM * tdt) /  &
     &        (mass(il,ij,1) * (AVOGAD * GPKG / MWTAIR))
!
         end do
      end do
!
!
      return
!
      end subroutine Update_CH3I
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update_Beryl
!
! !INTERFACE:
!
      subroutine Update_Beryl (latdeg, press3c, concentration, IBe7, IBe10, &
     &                         pr_diag, loc_proc, &
     &                         tdt, be_opt, t_half_be7, t_half_be10, &
     &                         yield_be7, yield_be10, &
     &                         i1, i2, ju1, j2, k1, k2, &
     &                         ilo, ihi, julo, jhi, ju1_gl, j2_gl, num_species)
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiInterpolation_mod     , only : Interp_Bilinear
!
      implicit none
!
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
#     include "gmi_cosmogenic.h"
!
! !INPUT PARAMETERS:
      logical, intent(in   ) :: pr_diag
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: num_species
      integer, intent(in   ) :: IBe7, IBe10
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in   ) :: ilo, ihi, julo, jhi, ju1_gl, j2_gl
      integer, intent(in   ) :: be_opt
      real*8 , intent(in   ) :: tdt
      real*8 , intent(in   ) :: yield_be7, yield_be10
      real*8 , intent(in   ) :: t_half_be7
!                             ! half life of Beryllium-7,  or other cosmogenic radionuclide (days)
      real*8 , intent(in   ) :: t_half_be10
!                             ! half life of Beryllium-10, or other cosmogenic radionuclide (days)
      real*8 , intent(in   ) :: latdeg (ju1_gl:j2_gl)
!                             ! latitude of grid box center (deg)
      real*8 , intent(in   ) :: press3c(ilo:ihi, julo:jhi, k1:k2)
!                             ! atmospheric pressure at the center of each grid box (mb)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!                             ! species concentration (mixing ratio)
!
! !DESCRIPTION:
!   This routine updates the Beryllium chemistry.
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      logical, save :: first = .true.
!
      integer :: il, ij, ik
!
      real*8, allocatable, save :: be7_source_atmo (:,:,:)
      real*8, allocatable, save :: be10_source_atmo(:,:,:)
!
      real*8, allocatable, save :: star_be7 (:,:,:)
      real*8, allocatable, save :: star_be10(:,:,:)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Update_Beryl called by ', loc_proc
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
        Allocate (be7_source_atmo (i1:i2, ju1:j2, k1:k2))
        Allocate (be10_source_atmo(i1:i2, ju1:j2, k1:k2))
        be7_source_atmo = 0.0d0; be10_source_atmo = 0.0d0
!
        Allocate (star_be7 (i1:i2, ju1:j2, k1:k2))
        Allocate (star_be10(i1:i2, ju1:j2, k1:k2))
        star_be7 = 0.0d0; star_be10 = 0.0d0
!
!
!       --------------------------------------------------
!       Interpolate the star values onto the current grid.
!       --------------------------------------------------
!
        if (be_opt == 1) then
!
          do ik = k1, k2
            do ij = ju1, j2
              do il = i1, i2
!
                call Interp_Bilinear  &
     &            (latdeg(ij), press3c(il,ij,ik),  &
     &             star_be7(il,ij,ik), CR1_LAT, CR1_PRES, NUM_CR1_LAT,  &
     &             NUM_CR1_PRES, star_cr1_be)
!
                call Interp_Bilinear  &
     &            (latdeg(ij), press3c(il,ij,ik),  &
     &             star_be10(il,ij,ik), CR1_LAT, CR1_PRES, NUM_CR1_LAT,  &
     &             NUM_CR1_PRES, star_cr1_be)
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
                call Interp_Bilinear  &
     &            (latdeg(ij), press3c(il,ij,ik),  &
     &             star_be7(il,ij,ik), CR2_LAT, CR2_PRES, NUM_CR2_LAT,  &
     &             NUM_CR2_PRES, star_cr2_be7)
!
                call Interp_Bilinear  &
     &            (latdeg(ij), press3c(il,ij,ik),  &
     &             star_be10(il,ij,ik), CR2_LAT, CR2_PRES, NUM_CR2_LAT,  &
     &             NUM_CR2_PRES, star_cr2_be10)
!
              end do
            end do
          end do
!
        end if
!
!
        be7_source_atmo (:,:,:) = star_be7 (:,:,:) * yield_be7
        be10_source_atmo(:,:,:) = star_be10(:,:,:) * yield_be10
!
!     ======
      end if
!     ======
!
!
!     --------------------------------------------------------
!     Beryllium-7, with a half-life of t_half_be7, decreases,
!     due to decay to Lithium-7.  We do not model the Lithium.
!     t_half_be7 can be specified in the namelist file,
!     otherwise it will default to the half-life of Be-7.
!     --------------------------------------------------------
!
      concentration(IBe7)%pArray3D(:,:,:) =  &
     &  Exp (-GMI_LN2 * tdt / (t_half_be7 * SECPDY)) *  &
     &  concentration(IBe7)%pArray3D(:,:,:)
!
!
!     ---------------------------------------------------------
!     Beryllium-10, with a half-life of t_half_be10, decreases,
!     due to decay to Boron-10.  We do not model the Boron.
!     t_half_be10 can be specified in the namelist file,
!     otherwise it will default to the half-life of Be-10.
!     ---------------------------------------------------------
!
      concentration(IBe10)%pArray3D(:,:,:) =  &
     &  Exp (-GMI_LN2 * tdt / (t_half_be10 * SECPDY)) *  &
     &  concentration(IBe10)%pArray3D(:,:,:)
!
!
!     ----------------------------------------
!     Beryllium is produced in the atmosphere.
!     ----------------------------------------
!
      concentration(IBe7)%pArray3D(:,:,:) =  &
     &  concentration(IBe7)%pArray3D(:,:,:) +  &
     &  be7_source_atmo(:,:,:) * MWTAIR * tdt / AVOGAD
!
      concentration(IBe10)%pArray3D(:,:,:) =  &
     &  concentration(IBe10)%pArray3D(:,:,:) +  &
     &  be10_source_atmo(:,:,:) * MWTAIR * tdt / AVOGAD
!
!
      return
!
      end subroutine Update_Beryl
!
!!-------------------------------------------------------------------------
!!BOP
!!
!! !ROUTINE: Update_Beryl_2
!!
!! !INTERFACE:
!!
!      subroutine Update_Beryl_2 (latdeg, press3c, concentration, IBe7S, IBe10S, IBe7T, IBe10T, &
!     &                         pr_diag, loc_proc, tdt, be_opt, t_half_be7, t_half_be10, &
!     &                         yield_be7, yield_be10, Ie90, &
!     &                         i1, i2, ju1, j2, k1, k2, &
!     &                         ilo, ihi, julo, jhi, ju1_gl, j2_gl, num_species)
!!
!! !USES:
!      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!      use GmiInterpolation_mod     , only : Interp_Bilinear
!!
!      implicit none
!
!#     include "gmi_phys_constants.h"
!#     include "gmi_time_constants.h"
!#     include "gmi_cosmogenic.h"
!!
!! !INPUT PARAMETERS:
!      logical, intent(in   ) :: pr_diag
!      integer, intent(in   ) :: loc_proc
!      integer, intent(in   ) :: num_species
!      integer, intent(in   ) :: IBe7S, IBe10S, IBe7T, IBe10T, Ie90
!      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2
!      integer, intent(in   ) :: ilo, ihi, julo, jhi, ju1_gl, j2_gl
!      integer, intent(in   ) :: be_opt
!      real*8 , intent(in   ) :: tdt
!      real*8 , intent(in   ) :: yield_be7, yield_be10
!      real*8 , intent(in   ) :: t_half_be7
!!                             ! half life of Beryllium-7,  or other cosmogenic radionuclide (days)
!      real*8 , intent(in   ) :: t_half_be10
!!                             ! half life of Beryllium-10, or other cosmogenic radionuclide (days)
!      real*8 , intent(in   ) :: latdeg (ju1_gl:j2_gl)
!!                             ! latitude of grid box center (deg)
!      real*8 , intent(in   ) :: press3c(ilo:ihi, julo:jhi, k1:k2)
!!                             ! atmospheric pressure at the center of each grid box (mb)
!!
!! !INPUT/OUTPUT PARAMETERS:
!      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!!                             ! species concentration (mixing ratio)
!!
!! !DESCRIPTION:
!!   This routine updates the Beryllium chemistry.
!!
!! !REVISION HISTORY:
!!  Initial code.
!!
!!EOP
!!-------------------------------------------------------------------------
!!     ----------------------
!!     Variable declarations.
!!     ----------------------
!
!      logical, save :: first = .true.
!
!      integer :: il, ij, ik
!
!      real*8, allocatable, save :: be7_source_atmo (:,:,:)
!      real*8, allocatable, save :: be10_source_atmo(:,:,:)
!
!      real*8, allocatable, save :: star_be7 (:,:,:)
!      real*8, allocatable, save :: star_be10(:,:,:)
!
!      real*8, save :: e90_trpps_val
!
!!     ----------------
!!     Begin execution.
!!     ----------------
!
!      if (pr_diag) then
!        Write (6,*) 'Update_Beryl_2 called by ', loc_proc
!      end if
!
!
!!     ==========
!      if (first) then
!!     ==========
!
!        first = .false.
!
!        e90_trpps_val = 90.0e-9
!
!        Allocate (be7_source_atmo (i1:i2, ju1:j2, k1:k2))
!        Allocate (be10_source_atmo(i1:i2, ju1:j2, k1:k2))
!        be7_source_atmo = 0.0d0
!        be10_source_atmo = 0.0d0
!
!        Allocate (star_be7 (i1:i2, ju1:j2, k1:k2))
!        Allocate (star_be10(i1:i2, ju1:j2, k1:k2))
!        star_be7 = 0.0d0
!        star_be10 = 0.0d0
!
!
!!       --------------------------------------------------
!!       Interpolate the star values onto the current grid.
!!       --------------------------------------------------
!
!        if (be_opt == 1) then
!
!          do ik = k1, k2
!            do ij = ju1, j2
!              do il = i1, i2
!
!                call Interp_Bilinear  &
!     &            (latdeg(ij), press3c(il,ij,ik),  &
!     &             star_be7(il,ij,ik), CR1_LAT, CR1_PRES, NUM_CR1_LAT,  &
!     &             NUM_CR1_PRES, star_cr1_be)
!
!                call Interp_Bilinear  &
!     &            (latdeg(ij), press3c(il,ij,ik),  &
!     &             star_be10(il,ij,ik), CR1_LAT, CR1_PRES, NUM_CR1_LAT,  &
!     &             NUM_CR1_PRES, star_cr1_be)
!
!              end do
!            end do
!          end do
!
!        else
!
!          do ik = k1, k2
!            do ij = ju1, j2
!              do il = i1, i2
!
!                call Interp_Bilinear  &
!     &            (latdeg(ij), press3c(il,ij,ik),  &
!     &             star_be7(il,ij,ik), CR2_LAT, CR2_PRES, NUM_CR2_LAT,  &
!     &             NUM_CR2_PRES, star_cr2_be7)
!
!                call Interp_Bilinear  &
!     &            (latdeg(ij), press3c(il,ij,ik),  &
!     &             star_be10(il,ij,ik), CR2_LAT, CR2_PRES, NUM_CR2_LAT,  &
!     &             NUM_CR2_PRES, star_cr2_be10)
!
!              end do
!            end do
!          end do
!
!        end if
!
!
!        be7_source_atmo (:,:,:) = star_be7 (:,:,:) * yield_be7
!        be10_source_atmo(:,:,:) = star_be10(:,:,:) * yield_be10
!
!!     ======
!      end if
!!     ======
!
!
!!     --------------------------------------------------------
!!     Beryllium-7, with a half-life of t_half_be7, decreases,
!!     due to decay to Lithium-7.  We do not model the Lithium.
!!     t_half_be7 can be specified in the namelist file,
!!     otherwise it will default to the half-life of Be-7.
!!     --------------------------------------------------------
!
!      concentration(IBe7S)%pArray3D(:,:,:) =  &
!     &  Exp (-GMI_LN2 * tdt / (t_half_be7 * SECPDY)) *  &
!     &  concentration(IBe7S)%pArray3D(:,:,:)
!      concentration(IBe7T)%pArray3D(:,:,:) =  &
!     &  Exp (-GMI_LN2 * tdt / (t_half_be7 * SECPDY)) *  &
!     &  concentration(IBe7T)%pArray3D(:,:,:)
!
!
!!     ---------------------------------------------------------
!!     Beryllium-10, with a half-life of t_half_be10, decreases,
!!     due to decay to Boron-10.  We do not model the Boron.
!!     t_half_be10 can be specified in the namelist file,
!!     otherwise it will default to the half-life of Be-10.
!!     ---------------------------------------------------------
!
!      concentration(IBe10S)%pArray3D(:,:,:) =  &
!     &  Exp (-GMI_LN2 * tdt / (t_half_be10 * SECPDY)) *  &
!     &  concentration(IBe10S)%pArray3D(:,:,:)
!      concentration(IBe10T)%pArray3D(:,:,:) =  &
!     &  Exp (-GMI_LN2 * tdt / (t_half_be10 * SECPDY)) *  &
!     &  concentration(IBe10T)%pArray3D(:,:,:)
!
!
!!     ----------------------------------------
!!     Beryllium is produced in the atmosphere.
!!     ----------------------------------------
!
!      where(concentration(Ie90)%pArray3D(:,:,:) .le. e90_trpps_val)
!        concentration(IBe7S)%pArray3D(:,:,:) =  &
!          concentration(IBe7S)%pArray3D(:,:,:) +  &
!          be7_source_atmo(:,:,:) * MWTAIR * tdt / AVOGAD
!      elsewhere
!        concentration(IBe7T)%pArray3D(:,:,:) =  &
!          concentration(IBe7T)%pArray3D(:,:,:) +  &
!          be7_source_atmo(:,:,:) * MWTAIR * tdt / AVOGAD
!      endwhere
!
!      where(concentration(Ie90)%pArray3D(:,:,:) .le. e90_trpps_val)
!        concentration(IBe10S)%pArray3D(:,:,:) =  &
!          concentration(IBe10S)%pArray3D(:,:,:) +  &
!          be10_source_atmo(:,:,:) * MWTAIR * tdt / AVOGAD
!      elsewhere
!        concentration(IBe10T)%pArray3D(:,:,:) =  &
!          concentration(IBe10T)%pArray3D(:,:,:) +  &
!          be10_source_atmo(:,:,:) * MWTAIR * tdt / AVOGAD
!      endwhere
!
!
!      return
!
!      end subroutine Update_Beryl_2
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update_Radon_Lead_Beryl  (hyl, 08/30/06, 11/23/10)
!
! !INTERFACE:
!
!
      subroutine Update_Radon_Lead_Beryl  &
     &  (latdeg, tropp, press3c, lwi_flags, mcor, mass, mw, &
     &   surf_air_temp, concentration, IRn, IPb, IPbS, IBe7, IBe7S, IBe10, IBe10S,  &
     &   decay_3d_out, &                                                 !jules,hyl
     &   pr_surf_emiss, pr_emiss_3d, surf_emiss_out, emiss_3d_out, &     !hyl
     &   pr_diag, loc_proc, &
     &   tdt, be_opt, t_half_be7, t_half_be10, &
     &   yield_be7, yield_be10, &
     &   i1, i2, ju1, j2, k1, k2, &
     &   ilo, ihi, julo, jhi, ju1_gl, j2_gl, num_species)
!
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiInterpolation_mod     , only : Interp_Bilinear
!
      implicit none
!
#     include "gmi_phys_constants.h"
#     include "gmi_time_constants.h"
#     include "gmi_cosmogenic.h"
!
! !INPUT PARAMETERS:
      logical, intent(in   ) :: pr_surf_emiss, pr_emiss_3d    !hyl
      logical, intent(in   ) :: pr_diag
      integer, intent(in   ) :: loc_proc
      integer, intent(in   ) :: num_species, IRn, IPb, IPbS, IBe7, IBe7S, IBe10, IBe10S
      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in   ) :: ilo, ihi, julo, jhi, ju1_gl, j2_gl
      real*8 , intent(in   ) :: mw(num_species)                                   !hyl
      real*8 , intent(inout) :: surf_emiss_out (i1:i2, ju1:j2, num_species)       !hyl
                              ! accumulation of surface emission for output (kg/m^2/constOutputFrequency)
!.. emiss_3d_out, unit: kg/m2/constOutputFrequency (hyl,1/5/2011)
      real*8 , intent(inout) :: emiss_3d_out (i1:i2, ju1:j2, k1:k2, num_species)  !hyl
                              ! accumulation of 3d emissions for output (kg/m^2/constOutputFrequency)
      real*8 , intent(inout) :: decay_3d_out (i1:i2, ju1:j2, k1:k2, num_species)  !jules,hyl
      integer, intent(in   ) :: be_opt
      real*8 , intent(in   ) :: tdt
      integer, intent(in   ) :: lwi_flags     (i1:i2, ju1:j2)
                              ! array of flags that indicate land, water, or ice
      real*8 , intent(in   ) :: mcor          (i1:i2, ju1:j2)
                              ! area of grid box (m^2)
      real*8 , intent(in   ) :: mass          (i1:i2, ju1:j2, k1:k2)
                              ! total mass of the atmosphere within each grid box (kg)
      real*8 , intent(in   ) :: surf_air_temp (i1:i2, ju1:j2)
                              ! surface area temperature (degK)
      real*8 , intent(in   ) :: yield_be7, yield_be10
      real*8 , intent(in   ) :: t_half_be7
!                             ! half life of Beryllium-7,  or other cosmogenic radionuclide (days)
      real*8 , intent(in   ) :: t_half_be10
!                             ! half life of Beryllium-10, or other cosmogenic radionuclide (days)
      real*8 , intent(in   ) :: latdeg (ju1_gl:j2_gl)
!                             ! latitude of grid box center (deg)
      real*8 , intent(in   ) :: tropp (i1:i2, ju1:j2)
!                             ! tropopause pressure (mb)
      real*8 , intent(in   ) :: press3c(ilo:ihi, julo:jhi, k1:k2)
!                             ! atmospheric pressure at the center of each grid box (mb)
! !INPUT/OUTPUT PARAMETERS:
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!                             ! species concentration (mixing ratio)
!
! !DESCRIPTION:
!   This routine updates the Radon_Lead_Beryl chemistry.
!
! !REVISION HISTORY:
!   Initial code.
!   Merged Update_Radon_Lead and Update_Beryl to Update_Radon_Lead_Beryl (hyl, 08/30/06)
!EOP
!-------------------------------------------------------------------------
!
!     -----------------------
!     Parameter declarations.
!     -----------------------
!
      real*8, parameter :: CMPM2 = CMPM * CMPM
!
!     -----------------------------------------------
!     Radon land scaling factor for cold temperatures
!     (Rind+Lerner, JGR 1996).
!     -----------------------------------------------
!
      real*8, parameter :: COLD_TEMP_LAND     = 273.0d0    ! degK
!     real*8, parameter :: COLD_TEMP_LAND_FAC =   0.6d0    !orig
      real*8, parameter :: COLD_TEMP_LAND_FAC =   0.333d0  !hyl,10/07/06
!
      real*8, parameter :: EFOLD_TIME_RN222   =   5.5d0
      real*8, parameter :: EFOLD_TIME_PB210   =11742.8d0   !half-life 22.3years, e-fold 22.3/ln2=22.3/alog(2)=32.17 years
!
      real*8, parameter :: RN_SOURCE_LAND     =   1.000d0  ! atoms/cm^2/s
      real*8, parameter :: RN_SOURCE_OCEAN    =   0.005d0  ! atoms/cm^2/s
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      logical, save :: first = .true.   !hyl
!
      integer :: il, ij, ik             !hyl
      real*8  :: ste_factor
      integer :: ltpause(i1:i2, ju1:j2)
!
      real*8  :: rn_exp_fac, pb_exp_fac  !hyl
      real*8  :: rn_source
!
      real*8, allocatable, save :: be7_source_atmo (:,:,:)
      real*8, allocatable, save :: be10_source_atmo(:,:,:)
!
      real*8, allocatable, save :: star_be7 (:,:,:)
      real*8, allocatable, save :: star_be10(:,:,:)
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (pr_diag) then
        Write (6,*) 'Update_Radon_Lead_Beryl called by ', loc_proc
      end if
!
!... locate the model layer of tropopause (hyl, 12/02/07)
       do ik=k2,k1,-1
         do ij=ju1,j2
           do il=i1,i2
             if (press3c(il,ij,ik) < tropp(il,ij)) then
               !tropp is grid center pressure in calcTropopausePress_WMO (hyl,3/9/11)
                ltpause(il,ij) = ik-1  !standard (hyl,3/9/11)
               !ltpause(il,ij) = ik    !orig
             endif
           enddo
         enddo
       enddo
!
!     ----------------------------------
!     Lead increases due to Radon decay.
!     ----------------------------------
!
      rn_exp_fac = Exp (-tdt / (EFOLD_TIME_RN222 * SECPDY))
      pb_exp_fac = Exp (-tdt / (EFOLD_TIME_Pb210 * SECPDY))
!
      concentration(IPb)%pArray3D(i1:i2,ju1:j2,k1:k2) =  &
        concentration(IPb)%pArray3D(i1:i2,ju1:j2,k1:k2) +  &
        (1.0d0 - rn_exp_fac) * concentration(IRn)%pArray3D(i1:i2,ju1:j2,k1:k2)
!
!... diagnostic
      if (pr_emiss_3d) then
        do ik = k1, k2
          emiss_3d_out(i1:i2,ju1:j2,ik,IPb) =  emiss_3d_out(i1:i2,ju1:j2,ik,IPb) +         &
             (1.0d0 - rn_exp_fac) * concentration(IRn)%pArray3D(i1:i2,ju1:j2,ik)   &
             * mw(IPb)/MWTAIR * mass(i1:i2,ju1:j2,ik) / mcor (i1:i2,ju1:j2)
        enddo
      endif
!
!... stratospheric Lead-210 due to the Source Gas Injection (SGI) pathway (i.e.
!... Lead-210 resulting from the decay of Radon-222 that is in the
!... stratosphere at the time of its decay. See Considine et al. ACP2005
!... for details (hyl,09/01/06)
!
!... this section of code causes "IOT trap signal" !hyl, 2/7/11:
       do ij=ju1,j2
         do il=i1,i2
           do ik=k1,k2
             if (ik > ltpause(il,ij) ) then
               concentration(IPbS)%pArray3D(il,ij,ik) =  &
                 concentration(IPbS)%pArray3D(il,ij,ik) +  &
                 (1.0d0 - rn_exp_fac) * concentration(IRn)%pArray3D(il,ij,ik)
!
!... diagnostic
               if (pr_emiss_3d) then
                 emiss_3d_out(il,ij,ik,IPbS) =  emiss_3d_out(il,ij,ik,IPbS)             &
                   + (1.0d0 - rn_exp_fac) * concentration(IRn)%pArray3D(il,ij,ik)   &
                   * mw(IPbS)/MWTAIR * mass(il,ij,ik) / mcor(il,ij)
               endif
             endif
           enddo
         enddo
       enddo
!----------------------------------------------------------------------
!
!     ----------------------------------
!     Radon and Lead decrease due to radioactive decay.
!     ----------------------------------
      concentration(IRn)%pArray3D(i1:i2,ju1:j2,k1:k2) = rn_exp_fac * concentration(IRn)%pArray3D(i1:i2,ju1:j2,k1:k2)
!
      concentration(IPb)%pArray3D(i1:i2,ju1:j2,k1:k2) = pb_exp_fac * concentration(IPb)%pArray3D(i1:i2,ju1:j2,k1:k2)
!
!... diagnostic
      do ik = k1, k2
!... Rn
        decay_3d_out(i1:i2,ju1:j2,ik,IRn) = decay_3d_out(i1:i2,ju1:j2,ik,IRn) +     &
           (1 - rn_exp_fac) * concentration(IRn)%pArray3D(i1:i2,ju1:j2,ik)  &
           * mw(IRn)/MWTAIR * mass(i1:i2,ju1:j2,ik) / mcor(i1:i2,ju1:j2)
!... 210Pb
        decay_3d_out(i1:i2,ju1:j2,ik,IPb) = decay_3d_out(i1:i2,ju1:j2,ik,IPb) +     &
           (1 - pb_exp_fac) * concentration(IPb)%pArray3D(i1:i2,ju1:j2,ik)  &
           * mw(IPb)/MWTAIR * mass(i1:i2,ju1:j2,ik) / mcor(i1:i2,ju1:j2)
!... strat 210Pb
        decay_3d_out(i1:i2,ju1:j2,ik,IPbS) = decay_3d_out(i1:i2,ju1:j2,ik,IPbS) +     &
           (1 - pb_exp_fac) * concentration(IPbS)%pArray3D(i1:i2,ju1:j2,ik)  &
           * mw(IPbS)/MWTAIR * mass(i1:i2,ju1:j2,ik) / mcor(i1:i2,ju1:j2)
      enddo
!
!     ---------------------------------------------------------------------
!     The radon comes from the land and the ocean.  This was originally
!     scaled in the code to match the lwi flags in the various met fields,
!     but this is no longer done.  Instead, if the radon needs to be scaled
!     to get exactly 16 kg/yr, the scaling must be done to the output,
!     after the run is completed.
!     ---------------------------------------------------------------------
!
      do ij = ju1, j2
        do il = i1, i2
!
          if (lwi_flags(il,ij) == 1) then  ! ocean
!
            rn_source = RN_SOURCE_OCEAN
            ! radon emissions at high latitudes following Jacob et al. [1997] (hyl,10/9/2006)
            if (abs(latdeg(ij)).ge.70.0d0) then
              rn_source = 0.0d0
            endif
            if (abs(latdeg(ij)).ge.60.0d0 .and. &
              abs(latdeg(ij)).lt.70.0d0) then
              rn_source = 0.005
            endif
          else if (lwi_flags(il,ij) == 3) then  ! ice
            rn_source = 0.0d0
!
          else if (lwi_flags(il,ij) == 2) then  ! land
!
           if (abs(latdeg(ij)).ge.70.0d0) then
              rn_source = 0.0d0
            endif
            if (abs(latdeg(ij)).ge.60.0d0 .and. &
              abs(latdeg(ij)).lt.70.0d0) then
              rn_source = 0.005
            endif
            if (abs(latdeg(ij)).lt.60.0d0) then
              rn_source = RN_SOURCE_LAND
!
              if (surf_air_temp(il,ij) < COLD_TEMP_LAND) then
                rn_source = rn_source * COLD_TEMP_LAND_FAC
              end if
            endif
!
          end if
          if ((lwi_flags(il,ij) == 1) .or.  &
              (lwi_flags(il,ij) == 2)) then
!
            concentration(IRn)%pArray3D(il,ij,1) =  &
              concentration(IRn)%pArray3D(il,ij,1) +  &
              (rn_source * mcor(il,ij) * CMPM2 * tdt) /  &
              (mass(il,ij,1) * (AVOGAD * GPKG / MWTAIR))
!
!... surf_emiss_out: unit = "kg/m2/constOutputFrequency" (hyl,1/5/11)
            if (pr_surf_emiss) then
              surf_emiss_out(il,ij,IRn) =  surf_emiss_out(il,ij,IRn)   &
                + (rn_source * CMPM2 * tdt) / AVOGAD * mw(IRn) / GPKG
            endif
!
!... Use emiss_3d_out to diagnose 222Rn emissions. This should produce the same
!... results as surf_emiss_out (hyl 1/20/11 - tested OK)
        !emiss_3d_out, unit: kg/m2/constOutputFrequency (hyl,1/5/2011)
            if (pr_emiss_3d) then
              emiss_3d_out(il,ij,1,IRn) =  emiss_3d_out(il,ij,1,IRn)  &
                + (rn_source * CMPM2 * tdt) / AVOGAD * mw(IRn) / GPKG
            endif
!
          endif
        enddo
      enddo
!
!============= Beryllium ==============
!
!     ==========
      if (first) then
!     ==========
!
        first = .false.
!
        Allocate (be7_source_atmo (i1:i2, ju1:j2, k1:k2))
        Allocate (be10_source_atmo(i1:i2, ju1:j2, k1:k2))
        be7_source_atmo = 0.0d0
        be10_source_atmo = 0.0d0
!
        Allocate (star_be7 (i1:i2, ju1:j2, k1:k2))
        Allocate (star_be10(i1:i2, ju1:j2, k1:k2))
        star_be7 = 0.0d0
        star_be10 = 0.0d0
!
!       --------------------------------------------------
!       Interpolate the star values onto the current grid.
!       --------------------------------------------------
!
        if (be_opt == 1) then
!
          do ik = k1, k2
            do ij = ju1, j2
              do il = i1, i2
!
                call Interp_Bilinear  &
                  (latdeg(ij), press3c(il,ij,ik),  &
                   star_be7(il,ij,ik), CR1_LAT, CR1_PRES, NUM_CR1_LAT,  &
                   NUM_CR1_PRES, star_cr1_be)
!
                call Interp_Bilinear  &
                  (latdeg(ij), press3c(il,ij,ik),  &
                   star_be10(il,ij,ik), CR1_LAT, CR1_PRES, NUM_CR1_LAT,  &
                   NUM_CR1_PRES, star_cr1_be)
!
              end do
            end do
          end do
!
        else
          do ik = k1, k2
            do ij = ju1, j2
              do il = i1, i2
!
                call Interp_Bilinear  &
                  (latdeg(ij), press3c(il,ij,ik),  &
                   star_be7(il,ij,ik), CR2_LAT, CR2_PRES, NUM_CR2_LAT,  &
                   NUM_CR2_PRES, star_cr2_be7)
!
                call Interp_Bilinear  &
                  (latdeg(ij), press3c(il,ij,ik),  &
                   star_be10(il,ij,ik), CR2_LAT, CR2_PRES, NUM_CR2_LAT,  &
                   NUM_CR2_PRES, star_cr2_be10)
!
!
              enddo
            enddo
          enddo
!
        endif
!
!... be7_source_atmo and be10_source_atmo unit: atoms/g/s (hyl,1/12/2011)
        be7_source_atmo (i1:i2,ju1:j2,k1:k2) = star_be7 (i1:i2,ju1:j2,k1:k2) * yield_be7
        be10_source_atmo(i1:i2,ju1:j2,k1:k2) = star_be10(i1:i2,ju1:j2,k1:k2) * yield_be10
!
      endif   ! "first"
!
!--------------------------------------------------------
!... Beryllium-7, with a half-life of t_half_be7, decreases,
!...  due to decay to Lithium-7.  We do not model the Lithium.
!...  t_half_be7 can be specified in the namelist file,
!...  otherwise it will default to the half-life of Be-7.
!--------------------------------------------------------
!
!... Beryllium-7 decay (hyl,09/01/06)
      concentration(IBe7)%pArray3D(i1:i2,ju1:j2,k1:k2) =  &
        Exp (-GMI_LN2 * tdt / (t_half_be7 * SECPDY)) *  &
        concentration(IBe7)%pArray3D(i1:i2,ju1:j2,k1:k2)
!
!... stratospheric Beryllium-7 decay (hyl,09/01/06)
      concentration(IBe7S)%pArray3D(i1:i2,ju1:j2,k1:k2) =  &
        Exp (-GMI_LN2 * tdt / (t_half_be7 * SECPDY)) *  &
        concentration(IBe7S)%pArray3D(i1:i2,ju1:j2,k1:k2)
!
!---------------------------------------------------------
!... Beryllium-10, with a half-life of t_half_be10, decreases,
!... due to decay to Boron-10.  We do not model the Boron.
!... t_half_be10 can be specified in the namelist file,
!... otherwise it will default to the half-life of Be-10.
!---------------------------------------------------------
!
!... Beryllium-10 decay (hyl,09/01/06)
      concentration(IBe10)%pArray3D(i1:i2,ju1:j2,k1:k2) =  &
        Exp (-GMI_LN2 * tdt / (t_half_be10 * SECPDY)) *  &
        concentration(IBe10)%pArray3D(i1:i2,ju1:j2,k1:k2)
!
!... stratospheric Beryllium-10 decay (hyl,09/01/06)
      concentration(IBe10S)%pArray3D(i1:i2,ju1:j2,k1:k2) =  &
        Exp (-GMI_LN2 * tdt / (t_half_be10 * SECPDY)) *  &
        concentration(IBe10S)%pArray3D(i1:i2,ju1:j2,k1:k2)
!
!... diagnostic
      do ik = k1, k2
!... Beryllium-7 (hyl,09/01/06)
        decay_3d_out (i1:i2,ju1:j2,ik,IBe7) = decay_3d_out (i1:i2,ju1:j2,ik,IBe7)    &
         + ( 1 - Exp (-GMI_LN2 * tdt / (t_half_be7 * SECPDY)) )  &
         * concentration(IBe7)%pArray3D(i1:i2,ju1:j2,ik)                   &
         * mw(IBe7)/MWTAIR * mass(i1:i2,ju1:j2,ik) / mcor(i1:i2,ju1:j2)
!
!... stratospheric Beryllium-7 (hyl,09/01/06)
        decay_3d_out (i1:i2,ju1:j2,ik,IBe7S) = decay_3d_out (i1:i2,ju1:j2,ik,IBe7S)    &
         + ( 1 - Exp (-GMI_LN2 * tdt / (t_half_be7 * SECPDY)) )  &
         * concentration(IBe7S)%pArray3D(i1:i2,ju1:j2,ik)                   &
         * mw(IBe7S)/MWTAIR * mass(i1:i2,ju1:j2,ik) / mcor(i1:i2,ju1:j2)
        decay_3d_out (i1:i2,ju1:j2,ik,IBe10) = decay_3d_out (i1:i2,ju1:j2,ik,IBe10)    &
         + ( 1 - Exp (-GMI_LN2 * tdt / (t_half_be10 * SECPDY)) ) &
         * concentration(IBe10)%pArray3D(i1:i2,ju1:j2,ik)                   &
         * mw(IBe10)/MWTAIR * mass(i1:i2,ju1:j2,ik) / mcor(i1:i2,ju1:j2)
!
        !stratospheric Beryllium-10 (hyl,09/01/06)
        decay_3d_out (i1:i2,ju1:j2,ik,IBe10S) = decay_3d_out (i1:i2,ju1:j2,ik,IBe10S)    &
         + ( 1 - Exp (-GMI_LN2 * tdt / (t_half_be10 * SECPDY)) ) &
         * concentration(IBe10S)%pArray3D(i1:i2,ju1:j2,ik)                   &
         * mw(IBe10S)/MWTAIR * mass(i1:i2,ju1:j2,ik) / mcor(i1:i2,ju1:j2)
!
      enddo
!
!     ----------------------------------------
!     Beryllium is produced in the atmosphere.
!     ----------------------------------------
!
!... adjustment of cross-tropopause transport of Beryl (hyl, 11/26/07)
      ste_factor = 1.0d0
!
!... Reduce Beryl-7 and Beryl-10 emissions in the stratosphere
!... to adjust cross-tropopause transport. This is only for the
!... simulation of tropospheric Beryl-7 and Beryl-10. (hyl, 11/27/07)
!
       do ik=k1,k2
          do ij=ju1,j2
             do il=i1,i2
                if ( ik > ltpause(il,ij) ) then
!
                   concentration(IBe7)%pArray3D(il,ij,ik) =  &
                     concentration(IBe7)%pArray3D(il,ij,ik) +  &
                     be7_source_atmo(il,ij,ik) * MWTAIR * tdt / AVOGAD  &
                     / ste_factor   !hyl
!
                   concentration(IBe10)%pArray3D(il,ij,ik) =  &
                     concentration(IBe10)%pArray3D(il,ij,ik) +  &
                     be10_source_atmo(il,ij,ik) * MWTAIR * tdt / AVOGAD  &
                     / ste_factor   !hyl
!
                   !emiss_3d_out unit: kg/m2/constOutputFrequency;
                   !be7_source_atmo unit: atoms/g/s (hyl,1/12/2011)
                   if (pr_emiss_3d) then
                     emiss_3d_out(il,ij,ik,IBe7) =  emiss_3d_out(il,ij,ik,IBe7) &
                      + be7_source_atmo(il,ij,ik) * tdt / AVOGAD  &
                      * mw(IBe7) * mass(il,ij,ik) / mcor(il,ij)   &
                      / ste_factor   !hyl
                     emiss_3d_out(il,ij,ik,IBe10) =  emiss_3d_out(il,ij,ik,IBe10) &
                      + be10_source_atmo(il,ij,ik) * tdt / AVOGAD &
                      * mw(IBe10) * mass(il,ij,ik) / mcor(il,ij)  &
                      / ste_factor   !hyl
                   endif
!
                else
!
                   concentration(IBe7)%pArray3D(il,ij,ik) =  &
                     concentration(IBe7)%pArray3D(il,ij,ik) +  &
                     be7_source_atmo(il,ij,ik) * MWTAIR * tdt / AVOGAD
!
                   concentration(IBe10)%pArray3D(il,ij,ik) =  &
                     concentration(IBe10)%pArray3D(il,ij,ik) +  &
                     be10_source_atmo(il,ij,ik) * MWTAIR * tdt / AVOGAD
!
                   !be7_source_atmo unit: atoms/g/s (hyl,1/12/2011)
                   if (pr_emiss_3d) then
                     emiss_3d_out(il,ij,ik,IBe7) =  emiss_3d_out(il,ij,ik,IBe7) &
                      + be7_source_atmo(il,ij,ik) * tdt / AVOGAD   &
                      * mw(IBe7) * mass(il,ij,ik) / mcor(il,ij)
                     emiss_3d_out(il,ij,ik,IBe10) =  emiss_3d_out(il,ij,ik,IBe10) &
                      + be10_source_atmo(il,ij,ik) * tdt / AVOGAD  &
                      * mw(IBe10) * mass(il,ij,ik) / mcor(il,ij)
                   endif
!
                endif
!
             enddo
          enddo
       enddo
!
      do ij=ju1,j2
         do il=i1,i2
            do ik=k1,k2
               if ( ik > ltpause(il,ij) ) then
                  concentration(IBe7S)%pArray3D(il,ij,ik) =    &
                    concentration(IBe7S)%pArray3D(il,ij,ik) +  &
                    be7_source_atmo(il,ij,ik) * MWTAIR * tdt / AVOGAD  &
                    / ste_factor   !hyl
!
                  concentration(IBe10S)%pArray3D(il,ij,ik) =    &
                    concentration(IBe10S)%pArray3D(il,ij,ik) +  &
                    be10_source_atmo(il,ij,ik) * MWTAIR * tdt / AVOGAD &
                    / ste_factor   !hyl
!
                  !be7_source_atmo unit: atoms/g/s (hyl,1/12/2011)
                  if (pr_emiss_3d) then
                    emiss_3d_out(il,ij,ik,IBe7S) =  emiss_3d_out(il,ij,ik,IBe7S) &
                      +  be7_source_atmo(il,ij,ik) * tdt / AVOGAD  &
                      * mw(IBe7S) * mass(il,ij,ik) / mcor(il,ij)   &
                      / ste_factor   !hyl
                    emiss_3d_out(il,ij,ik,IBe10S) =  emiss_3d_out(il,ij,ik,IBe10S) &
                      + be10_source_atmo(il,ij,ik) * tdt / AVOGAD  &
                      * mw(IBe10S) * mass(il,ij,ik) / mcor(il,ij)  &
                      / ste_factor   !hyl
                  endif
!
               endif
            enddo
         enddo
      enddo
!
      return
!
      end subroutine Update_Radon_Lead_Beryl
!
!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialize_Linoz
!
! !INTERFACE:
!
!
  subroutine Initialize_Linoz (self, month, latdeg, averagePressEdge,  &
   i1, i2, ju1, j2, k1, k2, ju1_gl, j2_gl)
!
! USES:
  use m_netcdf_io_open,       only : Ncop_Rd
  use m_netcdf_io_close,      only : Nccl
  use m_netcdf_io_get_dimlen, only : Ncget_Dimlen
  use GmiInterpolation_mod,   only : Interp, Interp_Bilinear
!   read in parameters necessary for linoz calculation
!
  implicit none
!
# include "netcdf.inc"
!
! !INPUT PARAMETERS:
  integer, intent(in) :: month
  real*8 , intent(in) :: latdeg(ju1_gl:j2_gl)
  real*8 , intent(in) :: averagePressEdge(k1-1:k2)
  integer, intent(in) :: i1, i2, ju1, j2, k1, k2
  integer, intent(in) :: ju1_gl, j2_gl
!
! !OUTPUT PARAMETERS:
  type (t_Linoz), intent(inout)   :: self
!
! !DESCRIPTION:
!  Open a NetCDF file for reading of linoz parameters and interpolate to current grid
!
! !LOCAL VARIABLES:
  integer, parameter :: MAX_LENGTH_VAR_NAME = 132
  character (len=MAX_LENGTH_VAR_NAME) :: err_msg
  integer :: ierr
  integer :: ncid
  integer :: ii_num, jj_num, kk_num
  integer :: cnt1d (1)
  integer :: strt1d(1)
  integer :: cnt3d (3)
  integer :: strt3d(3)
  integer :: varid
  integer :: ij, ik, inum
  real*8  :: temp_out
!
!
! !AUTHOR:
!  Stephen Steenrod (GSFC/USRA)
!
!
!... dimension parameter names for Linoz input file
  character (len=MAX_LENGTH_VAR_NAME), parameter :: LAT_DNAM  = 'latitude_dim'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: PRES_DNAM = 'height_dim'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: TIME_DNAM = 'time_dim'
!... variable parameter names for Linoz input file
  character (len=MAX_LENGTH_VAR_NAME), parameter :: LAT_VNAM = 'latitude_dim'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: HGT_VNAM = 'height_dim'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: PRS_VNAM = 'pressure'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: TIME_VNAM = 'time_dim'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: CAR_VNAM = 'Cariolle_Loss'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: O3_VNAM = 'O3'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: T_VNAM = 'Temperature'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: COLO3_VNAM = 'colO3'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: PML_VNAM = 'PmL'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: DPMLDO3_VNAM = 'dPmLdO3'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: DPMLDT_VNAM = 'dPmLdT'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: DPMLDCOLO3_VNAM = 'dPmLdcolO3'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: O3SSFC_VNAM = 'O3SSFC'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: O3STAU_VNAM = 'O3STAU'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: O3PAUZ_VNAM = 'O3PAUZ'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: CLLOAD_VNAM = 'CLLOAD'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: TPSC_VNAM = 'TPSC'
!
  real, allocatable :: lat_linoz_in (:)
  real, allocatable :: height_linoz_in (:)
  real, allocatable :: pressure_linoz_in (:)
  real, allocatable :: Cariolle_Loss_linoz_in (:)
  real, allocatable :: O3_linoz_in (:,:,:)
  real, allocatable :: T_linoz_in (:,:,:)
  real, allocatable :: colO3_linoz_in (:,:,:)
  real, allocatable :: PmL_linoz_in (:,:,:)
  real, allocatable :: dPmLdO3_linoz_in (:,:,:)
  real, allocatable :: dPmLdT_linoz_in (:,:,:)
  real, allocatable :: dPmLdcolO3_linoz_in (:,:,:)
  real :: O3SSFC_linoz
  real :: O3STAU_linoz
  real :: O3PAUZ_linoz
  real :: CLLOAD_linoz
  real :: TPSC_linoz
!... work array
  real*8, allocatable :: temp (:,:)
  real*8, allocatable :: temp1d (:)
  real*8, allocatable :: alogplev (:)
  real*8, allocatable :: lat_linoz (:)
  real*8, allocatable :: logpressure_linoz (:)
!
! !REVISION HISTORY:
!  Initial code.
!
  call Ncop_Rd (ncid, self%linoz_infile_name)
!
  call Ncget_Dimlen (ncid, LAT_DNAM, ii_num)
  call Ncget_Dimlen (ncid, PRES_DNAM, jj_num)
  call Ncget_Dimlen (ncid, TIME_DNAM, kk_num)
!
!... get latitudes of input data
  strt1d(1) = 1
  cnt1d(1) = ii_num
!
  allocate (lat_linoz_in(ii_num))
!  call Ncrd_1d (lat_linoz_in, ncid, LAT_VNAM, strt1d, cnt1d)
  ierr = Nf_Inq_Varid (ncid, LAT_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, lat_linoz_in)
!
!... get vertical coordinates of input data
  strt1d(1) = 1
  cnt1d(1) = jj_num
!
  allocate (height_linoz_in (jj_num))
!  call Ncrd_1d (height_linoz_in, ncid, HGT_VNAM, strt1d, cnt1d)
  ierr = Nf_Inq_Varid (ncid, HGT_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, height_linoz_in)
!
  allocate (pressure_linoz_in (jj_num))
!  call Ncrd_1d (pressure_linoz_in, ncid, PRS_VNAM, strt1d, cnt1d)
  ierr = Nf_Inq_Varid (ncid, PRS_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, pressure_linoz_in)
!
  allocate (Cariolle_Loss_linoz_in (jj_num))
!  call Ncrd_1d (Cariolle_Loss_linoz_in, ncid, CAR_VNAM, strt3d, cnt3d)
  ierr = Nf_Inq_Varid (ncid, CAR_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, Cariolle_Loss_linoz_in)
!
!... get 3D coordinates of input data
  strt3d = (/1,1,1/)
  cnt3d = (/ii_num,jj_num,kk_num/)
!... Ozone climotology
  allocate (O3_linoz_in (ii_num,jj_num,kk_num))
!  call Ncrd_3d (O3_linoz_in, ncid, O3_VNAM, strt3d, cnt3d)
  ierr = Nf_Inq_Varid (ncid, O3_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt3d, cnt3d, O3_linoz_in)
!
!... Temperature climotology
  allocate (T_linoz_in (ii_num,jj_num,kk_num))
!  call Ncrd_3d (T_linoz_in, ncid, T_VNAM, strt3d, cnt3d)
  ierr = Nf_Inq_Varid (ncid, T_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt3d, cnt3d, T_linoz_in)
!
!... Column Ozone climotology
  allocate (colO3_linoz_in (ii_num,jj_num,kk_num))
!  call Ncrd_3d (colO3_linoz_in, ncid, COLO3_VNAM, strt3d, cnt3d)
  ierr = Nf_Inq_Varid (ncid, COLO3_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt3d, cnt3d, colO3_linoz_in)
!
!... Ozone Production - Loss climotology
  allocate (PmL_linoz_in (ii_num,jj_num,kk_num))
!  call Ncrd_3d (PmL_linoz_in, ncid, PML_VNAM, strt3d, cnt3d)
  ierr = Nf_Inq_Varid (ncid, PML_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt3d, cnt3d, PmL_linoz_in)
!
!... d(Ozone Production - Loss)/d(Ozone) climotology
  allocate (dPmLdO3_linoz_in (ii_num,jj_num,kk_num))
!  call Ncrd_3d (dPmLdO3_linoz_in, ncid, DPMLDO3_VNAM, strt3d, cnt3d)
  ierr = Nf_Inq_Varid (ncid, DPMLDO3_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt3d, cnt3d, dPmLdO3_linoz_in)
!
!... d(Ozone Production - Loss)/d(Temperature) climotology
  allocate (dPmLdT_linoz_in (ii_num,jj_num,kk_num))
!  call Ncrd_3d (dPmLdT_linoz_in, ncid, DPMLDT_VNAM, strt3d, cnt3d)
  ierr = Nf_Inq_Varid (ncid, DPMLDT_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt3d, cnt3d, dPmLdT_linoz_in)
!
!... d(Ozone Production - Loss)/d(Column Ozone) climotology
  allocate (dPmLdcolO3_linoz_in (ii_num,jj_num,kk_num))
!  call Ncrd_3d (dPmLdcolO3_linoz_in, ncid, DPMLDCOLO3_VNAM, strt3d, cnt3d)
  ierr = Nf_Inq_Varid (ncid, DPMLDCOLO3_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt3d, cnt3d, dPmLdcolO3_linoz_in)
!
!... Linoz parameters
!  O3SSFC
  ierr = Nf_Inq_Varid (ncid, O3SSFC_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, 1, 1, O3SSFC_linoz)
!  O3STAU
  ierr = Nf_Inq_Varid (ncid, O3STAU_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, 1, 1, O3STAU_linoz)
!  O3SPAUZ
  ierr = Nf_Inq_Varid (ncid, O3PAUZ_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, 1, 1, O3PAUZ_linoz)
!  CLLOAD
  ierr = Nf_Inq_Varid (ncid, CLLOAD_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, 1, 1, CLLOAD_linoz)
!  TPSC
  ierr = Nf_Inq_Varid (ncid, TPSC_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, 1, 1, TPSC_linoz)
!
  call Nccl (ncid)
!
!
!... Need to interpolate 3D fields to current model grid
  allocate (lat_linoz (ii_num))
  lat_linoz (:) = lat_linoz_in (:)
!... need to set end points to make sure pole is covered
  lat_linoz(1) = -91
  lat_linoz(ii_num) = 91
!
  allocate (logpressure_linoz (jj_num))
  logpressure_linoz (:) = log(pressure_linoz_in (:))
!
!... get average log of pressures in vertical of model grid
  allocate (alogplev (k1:k2))
  do ik = k1, k2
    alogplev(ik) = (log(averagePressEdge(ik-1))+log(averagePressEdge(ik)))*0.5d0
  enddo
!
  allocate (temp (ii_num,jj_num))
!... loop over all 3d input fields and put on model grid
  do inum = 1,7
!
    select case (inum)
      case (1)
        temp(:,:) = O3_linoz_in(:,:,month)
      case (2)
        temp(:,:) = T_linoz_in(:,:,month)
      case (3)
        temp(:,:) = colO3_linoz_in(:,:,month)
      case (4)
        temp(:,:) = PmL_linoz_in(:,:,month)
      case (5)
        temp(:,:) = dPmLdO3_linoz_in(:,:,month)
      case (6)
        temp(:,:) = dPmLdT_linoz_in(:,:,month)
      case (7)
        temp(:,:) = dPmLdcolO3_linoz_in(:,:,month)
      case default
    end select
!
    do ik = k1, k2
!
      do ij = ju1, j2
!
!... interpolate linoz input to current grid
        call Interp_Bilinear(latdeg(ij), alogplev(ik), temp_out,       &
                   lat_linoz, logpressure_linoz, ii_num, jj_num, temp)
!
!.. put result in appropriate array
        select case (inum)
          case (1)
            self%O3_linoz(ij,ik) = temp_out
          case (2)
            self%T_linoz(ij,ik) = temp_out
          case (3)
            self%colO3_linoz(ij,ik) = temp_out
          case (4)
            self%PmL_linoz(ij,ik) = temp_out
          case (5)
            self%dPmLdO3_linoz(ij,ik) = temp_out
          case (6)
            self%dPmLdT_linoz(ij,ik) = temp_out
          case (7)
            self%dPmLdcolO3_linoz(ij,ik) = temp_out
          case default
        end select
!
      enddo
    enddo
!
  enddo
!
!... Interpolate Cariolle Loss to current vertical grid
!... approximate the grid pt pressure by using sfc pressure of 1000 hPa
  allocate (temp1d (jj_num))
  temp1d(:) = Cariolle_Loss_linoz_in(:)
!
!... interpolate to current grid
  call Interp (logpressure_linoz, temp1d, jj_num, &
               alogplev, self%Cariolle_Loss_linoz, k2-k1+1 )
!
  self%O3SSFC_linoz = O3SSFC_linoz
  self%O3STAU_linoz = O3STAU_linoz
  self%O3PAUZ_linoz = O3PAUZ_linoz
  self%CLLOAD_linoz = CLLOAD_linoz
  self%TPSC_linoz   = TPSC_linoz
!
  deallocate (temp)
  deallocate (temp1d)
  deallocate (alogplev)
  deallocate (lat_linoz)
  deallocate (logpressure_linoz)
  deallocate (lat_linoz_in)
  deallocate (height_linoz_in)
  deallocate (pressure_linoz_in)
  deallocate (Cariolle_Loss_linoz_in)
  deallocate (O3_linoz_in)
  deallocate (T_linoz_in)
  deallocate (colO3_linoz_in)
  deallocate (PmL_linoz_in)
  deallocate (dPmLdO3_linoz_in)
  deallocate (dPmLdT_linoz_in)
  deallocate (dPmLdcolO3_linoz_in)
!
  return
!
  end subroutine Initialize_Linoz
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update_LINOZ
!
! !INTERFACE:
!
!-------------------------------------------------------------------------
!
  subroutine Update_LINOZ  &
   (linoz_infile_name, ILINOZ, concentration, nymd, nhms, tdt, londeg, latdeg,  &
    mcor, pr_diag, procID, kel, mass, gridBoxHeight, press3c,  &
    averagePressEdge, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi,  &
    i1_gl, i2_gl, ju1_gl, j2_gl, num_species)
!
!
! !USES:
  use GmiArrayBundlePointer_mod    , only : t_GmiArrayBundle
  use GmiTimeControl_mod, only : GmiSplitDateTime, GetSecondsFromJanuary1
  use GmiSolar_mod, only : CalcCosSolarZenithAngle
!
  implicit none
!
# include "gmi_phys_constants.h"
# include "gmi_time_constants.h"
# include "GmiParameters.h"
!# include "gmi_time_constants.h"
!
! !DESCRIPTION:
!   This routine updates the LINOZ "chemistry".
!
! !REVISION HISTORY:
!   Initial code.
!EOP
!-------------------------------------------------------------------------
!
!
!INPUT PARAMETERS:
!
  character (len=MAX_LENGTH_FILE_NAME),intent(in) :: linoz_infile_name
  integer, intent(in) :: ILINOZ
  integer, intent(in) :: nymd, nhms
  real*8,  intent(in) :: tdt
  real*8,  intent(in) :: londeg(i1_gl:i2_gl)
  real*8,  intent(in) :: latdeg(ju1_gl:j2_gl)
  real*8,  intent(in) :: mcor(i1:i2,ju1:j2)
  integer, intent(in) :: procID
  logical, intent(in) :: pr_diag
!
!  real*8,  intent(in) :: am(k1:k2)
!  real*8,  intent(in) :: bm(k1:k2)
!  real*8,  intent(in) :: pt
  real*8,  intent(in) :: kel(ilo:ihi,julo:jhi,k1:k2)
  real*8,  intent(in) :: mass(i1:i2,ju1:j2,k1:k2)
  real*8,  intent(in) :: gridBoxHeight(i1:i2,ju1:j2,k1:k2)
  real*8,  intent(in) :: press3c (ilo:ihi, julo:jhi, k1:k2) ! atmospheric pressure at the center of each grid box (mb)
  real*8,  intent(in) :: averagePressEdge (k1-1:k2)
!
  integer, intent(in) :: i1, i2, ju1, j2, k1, k2
  integer, intent(in) :: ilo, ihi, julo, jhi
  integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
  integer, intent(in) :: num_species
!
!
! !INPUT/OUTPUT PARAMETERS:
  type(t_GmiArrayBundle),       intent(inOut) :: concentration(num_species)
!
!
! !LOCAL VARIABLES:
  !... work arrays
  type(t_Linoz), save :: self
!
  real*8, allocatable :: cosSolarZenithAngle(:,:)
!
  real*8, allocatable :: linoz(:,:,:)
  real*8, allocatable :: col_o3 (:,:)
  real*8  :: new_linoz, sso3, pmlnet, dco3, dtmp, derco3, f0, zo3ssfc
  real*8  :: dertmp, dero3, climpml, climo3, clpsc, carloss, maxSZAPSCLoss, fxyzs
  integer :: i, j, k, k_10km, k_58km
!
  real*8, parameter :: DUCONST = 0.5 * AVOGAD / (MWTAIR * 2.687d+16 * 10.0d0)
!
!... pressure limits for validity of linoz "chemistry"
  real*8, parameter :: pmin = 0.2371374
  real*8, parameter :: pmax = 273.842  ! half way between 316.2278 and 237.1374
!... relax to O3SSFC_linoz @ 2day efold below this
  real*8, parameter :: psink = 486.967 ! half way between 562.3414 and 421.6965
!  integer, parameter :: levo3bl = 4
!
  real*8 :: time, days
  integer :: day, month, idumyear, nsec_jan1
  integer, save :: month_save = -999
  logical, save :: first = .true.
!
!---DOBF = 0.5*6.022d23*1000/(48.*1.d4*2.687d16) = weighting for half-layer
!  real*8,parameter ::  DOBF = 23345.45d0
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
  if (pr_diag) Write (6,*) 'Update_LINOZ called by ', procID
!
!... get date info
  call GmiSplitDateTime (nymd, idumyear, month, day)
!
!... calc sun angle for polar loss
  call GetSecondsFromJanuary1(nsec_jan1, nymd, nhms)
  time = nsec_jan1
  days = time / SECPDY
  allocate(cosSolarZenithAngle(i1:i2,ju1:j2))
  call CalcCosSolarZenithAngle(days, latDeg, lonDeg, cosSolarZenithAngle, &
         i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl)
!
!... need to read in appropriate month's climo values
  if(month_save .ne. month) then
    if(first) then
      call Allocate_O3_linoz (self, ju1, j2, k1, k2)
      call Allocate_T_linoz (self, ju1, j2, k1, k2)
      call Allocate_colO3_linoz (self, ju1, j2, k1, k2)
      call Allocate_PmL_linoz (self, ju1, j2, k1, k2)
      call Allocate_dPmLdO3_linoz (self, ju1, j2, k1, k2)
      call Allocate_dPmLdT_linoz (self, ju1, j2, k1, k2)
      call Allocate_dPmLdcolO3_linoz (self, ju1, j2, k1, k2)
      call Allocate_Cariolle_Loss_linoz (self, k1, k2)
      self%linoz_infile_name = linoz_infile_name
      first = .false.
    endif
!
    call Initialize_Linoz(self, month, latdeg, averagePressEdge,  &
      i1, i2, ju1, j2, k1, k2, ju1_gl, j2_gl)
!
    month_save = month
!
  endif
!
!... set up work space
  allocate (linoz(i1:i2,ju1:j2,k1:k2))
  allocate (col_o3(ju1:j2,k1:k2))
!
!... convert from vmr to mass
  linoz(:,:,:) = concentration(ILINOZ)%pArray3D(i1:i2,ju1:j2,k1:k2) * mass(i1:i2,ju1:j2,k1:k2)
!
!... get surface linoz mmr
  zo3ssfc = self%O3SSFC_linoz / (MWTAIR/(16.d0*3.d0))
!
!... relax to climo values below 237.1374 hPa
  if (self%O3STAU_linoz .gt. 300.d0) then
    fxyzs = 1.d0 - exp(-tdt/self%O3STAU_linoz)
  else
    fxyzs = 1.d0
  endif
!
!... loop over longitudes
  do i=i1,i2
!... calculate the column ozone of our calculated linoz, DUCONST has built in 1/2
    do j=ju1,j2
      col_o3(j,k2) = linoz(i, j, k2) * DUCONST / mcor(i,j)
    enddo
    do k=k2-1,k1,-1
      do j=ju1,j2
        col_o3(j,k) = col_o3(j,k+1) &
                     + (linoz(i,j,k) + linoz(i,j,k+1)) * DUCONST / mcor(i,j)
      enddo
    enddo
!
!... LINOZ chemistry over altitudes LZMIN:LM, boundary over 1:LZLBO3
    do k=k1,k2
      do j=ju1,j2
!
!... relax to climo in troposphere
        if(press3c(i,j,k).gt.psink) then
          linoz(i,j,k) = linoz(i,j,k) - (linoz(i,j,k)-zo3ssfc*mass(i,j,k))*fxyzs
!
!... calculate parameterization over UT and stratosphere
        elseif(press3c(i,j,k).le.pmax .and. press3c(i,j,k).ge.pmin) then
!
!---climatological ozone (v/v = mole fraction mixing ratio)
          climo3  = self%O3_linoz(j,k) * mass(i,j,k)
!---climatological p-l
          climpml = self%PmL_linoz(j,k)
!---partial derivative: d(P-L)/dO3 < 0
          dero3   = self%dPmLdO3_linoz(j,k)
!---partial derivative: d(P-L)/dT
          dertmp  = self%dPmLdT_linoz(j,k)
!---partial derivative: d(P-L)/dcol-O3
          derco3  = self%dPmLdcolO3_linoz(j,k)
!---differences from climatology
          dtmp = kel(i,j,k) - self%T_linoz(j,k)
          dco3 = col_o3(j,k) - self%colO3_linoz (j,k)
!
!---steady-state ozone:  can be negative in lower strat, but timescale is long
          sso3 = climo3 - mass(i,j,k) * (climpml+dco3*derco3+dtmp*dertmp)/dero3
!---use dero3 (<0) timescale to decay to Steady-State.
          pmlnet = (sso3 - linoz(i,j,k))*(1.d0 - exp(dero3*tdt))
          new_linoz = linoz(i,j,k) + pmlnet
!
!---PSC-activ loss depends on  T (incl. <195K & lat>40), sun above horizon(L)
!---enhanced polar O3 loss (Cariolle 1990) with loss freq TPSCLZ(LR) < 0
          if (self%CLLOAD_linoz .gt. 0.0d0) then
            if (abs(latdeg(j)) .gt. 40.0d0) then
              if (kel(i,j,k) .lt. self%TPSC_linoz) then
!... calculate maximum SZA for PSC loss (= tangent ht as sunset)
                f0 = max(16.d0*log10(1000.d0/(gridBoxHeight(i,j,k))),0.d0)
                f0 = (90.d0 + sqrt(f0)) * RADPDEG
                maxSZAPSCLoss = cos(f0)
                if (cosSolarZenithAngle(i,j) .gt. maxSZAPSCLoss) then
                  clpsc = (self%CLLOAD_linoz/2.75d0)**2
                  carloss = 1.d0 - exp(self%Cariolle_Loss_linoz(k)*clpsc*tdt)
                  new_linoz = linoz(i,j,k) + linoz(i,j,k)*carloss
                endif
              endif
            endif
          endif
!
          linoz(i,j,k) = new_linoz
!
!... linoz chemistry tables only valid to 58km
        elseif(press3c(i,j,k) .lt. pmin) then
          linoz(i,j,k) = 0.0d0
        endif
!
      enddo
    enddo
!
  enddo
!
!... convert to back to vmr
!
  concentration(ILINOZ)%pArray3D(i1:i2,ju1:j2,k1:k2) = linoz(i1:i2,ju1:j2,k1:k2) / mass(i1:i2,ju1:j2,k1:k2)
!
  deallocate (linoz)
  deallocate (col_o3)
  deallocate(cosSolarZenithAngle)
!
  return
!
  end subroutine Update_LINOZ
!
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialize_SF6
!
! !INTERFACE:
!
!
  subroutine Initialize_SF6 (self, month, londeg, latdeg, averagePressEdge,  &
   i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl)
!
! USES:
  use m_netcdf_io_open,       only : Ncop_Rd
  use m_netcdf_io_close,      only : Nccl
  use m_netcdf_io_get_dimlen, only : Ncget_Dimlen
  use GmiInterpolation_mod,   only : Interp, Interp_Trilinear
!   read in parameters necessary for SF6 calculation
!
  implicit none
!
# include "netcdf.inc"
# include "GmiParameters.h"
!
! !INPUT PARAMETERS:
  integer, intent(in) :: month
  real*8 , intent(in) :: londeg(i1_gl:i2_gl)
  real*8 , intent(in) :: latdeg(ju1_gl:j2_gl)
  real*8 , intent(in) :: averagePressEdge(k1-1:k2)
  integer, intent(in) :: i1, i2, ju1, j2, k1, k2
  integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
!
! !OUTPUT PARAMETERS:
!... There are 4 variables output:
!...  AD (air density, molec/cm3),
!...  j (photolysis rate, 1/s),
!...  k (rate constant for k*[O1D][SF6]; cm6/molec2/s,
!...  O1D (O1D; mixing ratio).
!
  type (t_SF6), intent(inout)   :: self
!
! !DESCRIPTION:
!  Open a NetCDF file for reading of SF6 parameters and interpolate to current grid
!
! !LOCAL VARIABLES:
  character (len=MAX_LENGTH_VAR_NAME) :: err_msg
  integer :: ierr
  integer :: ncid
  integer :: ii_num, jj_num, kk_num
  integer :: cnt1d (1)
  integer :: strt1d(1)
  integer :: cnt3d (4)
  integer :: strt3d(4)
  integer :: varid
  integer :: ii, ij, ik, inum
  real*8  :: temp_out
!
!
! !AUTHOR:
!  Stephen Steenrod (GSFC/USRA)
!
!
!... dimension parameter names for SF6 input file
  character (len=MAX_LENGTH_VAR_NAME), parameter :: LON_DNAM  = 'lon'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: LAT_DNAM  = 'lat'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: PRES_DNAM = 'level'
!... variable parameter names for SF6 input file
  character (len=MAX_LENGTH_VAR_NAME), parameter :: J_VNAM = 'j'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: K_VNAM = 'k'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: O1D_VNAM = 'O1D'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: AD_VNAM = 'AD'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: LON_VNAM = 'lon'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: LAT_VNAM = 'lat'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: PRS_VNAM = 'lev'
!
  real, allocatable :: lon_sf6_in (:)
  real, allocatable :: lat_sf6_in (:)
  real, allocatable :: pressure_sf6_in (:)
  real, allocatable :: ad_sf6_in (:,:,:)
  real, allocatable :: j_sf6_in (:,:,:)
  real, allocatable :: k_sf6_in (:,:,:)
  real, allocatable :: o1d_sf6_in (:,:,:)
!... work array
  real*8, allocatable :: temp (:,:,:)
  real*8, allocatable :: alogplev (:)
  real*8, allocatable :: lon_sf6 (:)
  real*8, allocatable :: lat_sf6 (:)
  real*8, allocatable :: logpressure_sf6 (:)
!
! !REVISION HISTORY:
!  Initial code.
!
  call Ncop_Rd (ncid, self%sf6_infile_name)
!
  call Ncget_Dimlen (ncid, LON_DNAM, ii_num)
  call Ncget_Dimlen (ncid, LAT_DNAM, jj_num)
  call Ncget_Dimlen (ncid, PRES_DNAM, kk_num)
!
!... get latitudes of input data
  strt1d(1) = 1
!
  cnt1d(1) = ii_num
  allocate (lon_sf6_in(ii_num))
!  call Ncrd_1d (lat_sf6_in, ncid, LAT_VNAM, strt1d, cnt1d)
  ierr = Nf_Inq_Varid (ncid, LON_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, lon_sf6_in)
!
  cnt1d(1) = jj_num
  allocate (lat_sf6_in(jj_num))
!  call Ncrd_1d (lat_sf6_in, ncid, LAT_VNAM, strt1d, cnt1d)
  ierr = Nf_Inq_Varid (ncid, LAT_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, lat_sf6_in)
!
!... get vertical coordinates of input data
  cnt1d(1) = kk_num
  allocate (pressure_sf6_in (jj_num))
!  call Ncrd_1d (pressure_sf6_in, ncid, PRS_VNAM, strt1d, cnt1d)
  ierr = Nf_Inq_Varid (ncid, PRS_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, pressure_sf6_in)
!
!... get 3D coordinates of input data
  strt3d = (/1,1,1,month/)
!
  cnt3d = (/ii_num,jj_num,kk_num,1/)
!... j(SF6) climotology
  allocate (j_sf6_in (ii_num,jj_num,kk_num))
!  call Ncrd_3d (j_sf6_in, ncid, O3_VNAM, strt3d, cnt3d)
  ierr = Nf_Inq_Varid (ncid, j_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt3d, cnt3d, j_sf6_in)
!
!... k(SF6) climotology
  allocate (k_sf6_in (ii_num,jj_num,kk_num))
!  call Ncrd_3d (k_sf6_in, ncid, T_VNAM, strt3d, cnt3d)
  ierr = Nf_Inq_Varid (ncid, k_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt3d, cnt3d, k_sf6_in)
!
!... AD climotology
  allocate (AD_sf6_in (ii_num,jj_num,kk_num))
!  call Ncrd_3d (colO3_sf6_in, ncid, COLO3_VNAM, strt3d, cnt3d)
  ierr = Nf_Inq_Varid (ncid, AD_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt3d, cnt3d, ad_sf6_in)
!
!... O1D climotology
  allocate (o1d_sf6_in (ii_num,jj_num,kk_num))
!  call Ncrd_3d (PmL_sf6_in, ncid, PML_VNAM, strt3d, cnt3d)
  ierr = Nf_Inq_Varid (ncid, O1D_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt3d, cnt3d, o1d_sf6_in)
!
  call Nccl (ncid)
!
!
!... Need to interpolate 3D fields to current model grid
  allocate (lon_sf6 (ii_num+2))
  lon_sf6 (2:ii_num+1) = lon_sf6_in (:)
!... need to set end points to wrap
  lon_sf6(1) = lon_sf6_in(1) - (lon_sf6_in(2) - lon_sf6_in(1))
  lon_sf6(ii_num+2) = lon_sf6_in(1) + (lon_sf6_in(ii_num) - lon_sf6_in(ii_num-1)) + 360
!
!... Need to interpolate 3D fields to current model grid
  allocate (lat_sf6 (jj_num))
  lat_sf6 (:) = lat_sf6_in (:)
!... need to set end points to make sure pole is covered
  lat_sf6(1) = -91
  lat_sf6(jj_num) = 91
!
  allocate (logpressure_sf6 (kk_num))
  logpressure_sf6 (:) = log(pressure_sf6_in (:))
!
!... get average log of pressures in vertical of model grid
  allocate (alogplev (k1:k2))
  do ik = k1, k2
    alogplev(ik) = (log(averagePressEdge(ik-1))+log(averagePressEdge(ik)))*0.5d0
  enddo
  allocate (temp (ii_num+2,jj_num,kk_num))
!... loop over all 3d input fields and put on model grid
  do inum = 1,4
!
    select case (inum)
      case (1)
        temp(2:ii_num+1,:,:) = j_sf6_in(:,:,:)
        temp(1,:,:)          = j_sf6_in(ii_num,:,:)
        temp(ii_num+2,:,:)   = j_sf6_in(1,:,:)
      case (2)
        temp(2:ii_num+1,:,:) = k_sf6_in(:,:,:)
        temp(1,:,:)          = k_sf6_in(ii_num,:,:)
        temp(ii_num+2,:,:)   = k_sf6_in(1,:,:)
      case (3)
        temp(2:ii_num+1,:,:) = AD_sf6_in(:,:,:)
        temp(1,:,:)          = AD_sf6_in(ii_num,:,:)
        temp(ii_num+2,:,:)   = AD_sf6_in(1,:,:)
      case (4)
        temp(2:ii_num+1,:,:) = O1D_sf6_in(:,:,:)
        temp(1,:,:)          = O1D_sf6_in(ii_num,:,:)
        temp(ii_num+2,:,:)   = O1D_sf6_in(1,:,:)
      case default
    end select
!
    do ik = k1, k2
!
      do ij = ju1, j2
!
        do ii = i1, i2
!
!... interpolate 3D sf6 input to current 3D grid
           call Interp_Trilinear(londeg(ii), latdeg(ij), alogplev(ik), temp_out,       &
                   lon_sf6, lat_sf6, logpressure_sf6, ii_num+2, jj_num, kk_num, temp)
!
!.. put result in appropriate array
           select case (inum)
              case (1)
                self%j_sf6(ii,ij,ik) = temp_out
              case (2)
                self%k_sf6(ii,ij,ik) = temp_out
              case (3)
                self%AD_sf6(ii,ij,ik) = temp_out
              case (4)
                self%O1D_sf6(ii,ij,ik) = temp_out
              case default
            end select
!
        enddo
      enddo
    enddo
!
  enddo
!
!
  deallocate (temp)
  deallocate (lon_sf6)
  deallocate (lon_sf6_in)
  deallocate (lat_sf6)
  deallocate (lat_sf6_in)
  deallocate (alogplev)
  deallocate (logpressure_sf6)
  deallocate (pressure_sf6_in)
  deallocate (j_sf6_in)
  deallocate (k_sf6_in)
  deallocate (AD_sf6_in)
  deallocate (O1D_sf6_in)
!
  return
!
  end subroutine Initialize_sf6
!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update_SF6
!
! !INTERFACE:
!
!-------------------------------------------------------------------------
!
  subroutine Update_SF6  &
   (SF6_infile_name, ISF6, concentration, nymd, tdt, londeg, latdeg,  &
    averagePressEdge, pr_diag, procID, i1, i2, ju1, j2, k1, k2,  &
    ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl, num_species)
!
!
! !USES:
  use GmiArrayBundlePointer_mod    , only : t_GmiArrayBundle
  use GmiTimeControl_mod, only : GmiSplitDateTime, GetSecondsFromJanuary1
!
  implicit none
!
# include "gmi_phys_constants.h"
# include "gmi_time_constants.h"
# include "GmiParameters.h"
!# include "gmi_time_constants.h"
!
! !DESCRIPTION:
!   This routine updates the SF6 "chemistry".
!
! !REVISION HISTORY:
!   Initial code.
!EOP
!-------------------------------------------------------------------------
!
!
!INPUT PARAMETERS:
!
  character (len=MAX_LENGTH_FILE_NAME),intent(in) :: SF6_infile_name
  integer, intent(in) :: ISF6
  integer, intent(in) :: nymd
  real*8,  intent(in) :: tdt
  real*8,  intent(in) :: londeg(i1_gl:i2_gl)
  real*8,  intent(in) :: latdeg(ju1_gl:j2_gl)
  real*8,  intent(in) :: averagePressEdge (k1-1:k2)
  logical, intent(in) :: pr_diag
  integer, intent(in) :: procID
!
  integer, intent(in) :: i1, i2, ju1, j2, k1, k2
  integer, intent(in) :: ilo, ihi, julo, jhi
  integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
  integer, intent(in) :: num_species
!
!
! !INPUT/OUTPUT PARAMETERS:
  type(t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!
!
! !LOCAL VARIABLES:
  type(t_SF6), save :: self
!
  integer :: i, j, k
!
  integer :: day, month, idumyear
  integer, save :: month_save = -999
  logical, save :: first = .true.
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
  if (pr_diag) Write (6,*) 'Update_SF6 called by ', procID
!
!... get date info
  call GmiSplitDateTime (nymd, idumyear, month, day)
!
!... need to read in appropriate month's climo values
  if(month_save .ne. month) then
    if(first) then
      call Allocate_J_SF6 (self, i1, i2, ju1, j2, k1, k2)
      call Allocate_K_SF6 (self, i1, i2, ju1, j2, k1, k2)
      call Allocate_AD_SF6 (self, i1, i2, ju1, j2, k1, k2)
      call Allocate_O1D_SF6 (self, i1, i2, ju1, j2, k1, k2)
      self%SF6_infile_name = SF6_infile_name
      first = .false.
    endif
!
    call Initialize_SF6(self, month, londeg, latdeg, averagePressEdge,  &
      i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl)
!
    month_save = month
!
  endif
!
!... SF6 chemistry involves losses by photolysis and O1D
  concentration(ISF6)%pArray3D(i1:i2,ju1:j2,k1:k2) = concentration(ISF6)%pArray3D(i1:i2,ju1:j2,k1:k2)  &
    - ( self%j_SF6(i1:i2,ju1:j2,k1:k2) * concentration(ISF6)%pArray3D(i1:i2,ju1:j2,k1:k2)  &
      + ( self%k_SF6(i1:i2,ju1:j2,k1:k2) * concentration(ISF6)%pArray3D(i1:i2,ju1:j2,k1:k2)  &
         * self%O1D_SF6(i1:i2,ju1:j2,k1:k2) * self%AD_SF6(i1:i2,ju1:j2,k1:k2)) ) * tdt
!
  return
!
  end subroutine Update_SF6
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update_Uniform
!
! !INTERFACE:
!
!-------------------------------------------------------------------------
!
 subroutine Update_Uniform  &
    (mcor, mass, concentration, IUniform, tdt,  &
     pr_diag, loc_proc, i1, i2, ju1, j2, ilo, ihi, julo, jhi, k1, k2, num_species)
!
! !USES:
  use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!
  implicit none
!
# include "gmi_phys_constants.h"
# include "gmi_time_constants.h"
!
! !INPUT PARAMETERS:
  real*8 , intent(in   ) :: mcor (i1:i2, ju1:j2)            ! area of grid box  (m^2)
  real*8 , intent(in   ) :: mass (i1:i2, ju1:j2, k1:k2)   ! total mass of the atmosphere within each grid box (kg)
!
  integer, intent(in   ) :: IUniform
  real*8 , intent(in   ) :: tdt
  logical, intent(in   ) :: pr_diag
  integer, intent(in   ) :: loc_proc
  integer, intent(in   ) :: i1, i2, ju1, j2, ilo, ihi, julo, jhi, k1, k2
  integer, intent(in   ) :: num_species
!
! !INPUT/OUTPUT PARAMETERS:
  type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
                              ! species concentration, known at zone centers (mixing ratio)
!
! !DESCRIPTION:
!   This routine updates the Uniform tracer - uniform surface emission of .
!
! !REVISION HISTORY:
!   Initial code.
!EOP
!-------------------------------------------------------------------------
!
! -----------------------
! Parameter declarations.
! -----------------------
!
  real*8, parameter :: EFOLD_TIME_Uniform = 25.0d0
!
  real*8, parameter :: Uniform_SOURCE  =  2400d12  ! g/y
  real*8, parameter :: Uniform_MW      =  28.0d0   ! CO like
!
!
! ----------------------
! Variable declarations.
! ----------------------
!
  integer :: il, ij
!
  real*8  :: exp_fac, source
!
!
! ----------------
! Begin execution.
! ----------------
!
  if (pr_diag) then
    Write (6,*) 'Update_Uniform called by ', loc_proc
  end if
!
! -----------------------------
! "Uniform" decreases due to decay.
! -----------------------------
!
  exp_fac = Exp (-tdt / (EFOLD_TIME_Uniform * SECPDY))
  concentration(IUniform)%pArray3D(:,:,:) = exp_fac * concentration(IUniform)%pArray3D(:,:,:)
!
!
! ---------------------------------------------------------------------
! The "Uniform" tracer is emitted uniformily from the surface - 2400 Tg/y
! ---------------------------------------------------------------------
! convert source to molecules/cm^2/s
  source = Uniform_source/Uniform_MW/(60.d0*60.d0*24.d0*365.25d0)*6.022d23/510.d16
!
  concentration(IUniform)%pArray3D(i1:i2,ju1:j2,1) =           &
         concentration(IUniform)%pArray3D(i1:i2,ju1:j2,1) +    &
         (source * mcor(i1:i2,ju1:j2) * CMPM * CMPM * tdt) /  &
         (mass(i1:i2,ju1:j2,1) * (AVOGAD * GPKG / MWTAIR))
!
  return
!
  end subroutine Update_Uniform
!EOC
!.sds!-------------------------------------------------------------------------
!.sds!BOP
!.sds!
!.sds! !ROUTINE: Initialize_stratO3
!.sds!
!.sds! !INTERFACE:
!.sds!
!.sds!
!.sds  subroutine Initialize_stratO3 (self, month, day, month_save, londeg, latdeg, press3c,  &
!.sds   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl)
!.sds!
!.sds! USES:
!.sds  use m_netcdf_io_open,       only : Ncop_Rd
!.sds  use m_netcdf_io_close,      only : Nccl
!.sds  use m_netcdf_io_get_dimlen, only : Ncget_Dimlen
!.sds  use GmiInterpolation_mod,   only : Interp, Interp_Trilinear
!.sds!   read in parameters necessary for Strat_O3 calculation
!.sds!
!.sds  implicit none
!.sds!
!.sds# include "netcdf.inc"
!.sds# include "GmiParameters.h"
!.sds!
!.sds! !INPUT PARAMETERS:
!.sds  integer, intent(in) :: month, day, month_save
!.sds  real*8 , intent(in) :: londeg(i1_gl:i2_gl)
!.sds  real*8 , intent(in) :: latdeg(ju1_gl:j2_gl)
!.sds  real*8,  intent(in) :: press3c (ilo:ihi, julo:jhi, k1:k2) ! atmospheric pressure at the center of each grid box (mb)
!.sds  integer, intent(in) :: i1, i2, ju1, j2, k1, k2
!.sds  integer, intent(in) :: ilo, ihi, julo, jhi
!.sds  integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
!.sds!
!.sds! !OUTPUT PARAMETERS:
!.sds!... There are 4 variables output:
!.sds!...  AD (air density, molec/cm3),
!.sds!...  j (photolysis rate, 1/s),
!.sds!...  k (rate constant for k*[O1D][SF6]; cm6/molec2/s,
!.sds!...  O1D (O1D; mixing ratio).
!.sds!
!.sds  type (t_stratO3), intent(inout)   :: self
!.sds!
!.sds! !DESCRIPTION:
!.sds!  Open a NetCDF file for reading of Strat_O3 parameters and interpolate to current grid
!.sds!
!.sds! !LOCAL VARIABLES:
!.sds  character (len=MAX_LENGTH_VAR_NAME) :: err_msg
!.sds  integer :: ierr
!.sds  integer :: ncid
!.sds  integer :: ii_num, jj_num, kk_num
!.sds  integer :: cnt1d (1)
!.sds  integer :: strt1d(1)
!.sds  integer :: cnt4d (4)
!.sds  integer :: strt4d(4)
!.sds  integer :: cnt5d (5)
!.sds  integer :: strt5d(5)
!.sds  integer :: varid
!.sds  integer :: ii, ij, ik, inum
!.sds  real*8  :: temp_out
!.sds  real*8  :: alogplev
!.sds!
!.sds!
!.sds! !AUTHOR:
!.sds!  Stephen Steenrod (GSFC/USRA)
!.sds!
!.sds!
!.sds!... dimension parameter names for Strat_O3 input file
!.sds  character (len=MAX_LENGTH_VAR_NAME), parameter :: LON_DNAM  = 'longitude'
!.sds  character (len=MAX_LENGTH_VAR_NAME), parameter :: LAT_DNAM  = 'latitude'
!.sds  character (len=MAX_LENGTH_VAR_NAME), parameter :: PRES_DNAM = 'pres'
!.sds!... variable parameter names for Strat_O3 input file
!.sds  character (len=MAX_LENGTH_VAR_NAME), parameter :: O3_VNAM  = 'FULLCHEM_O3'
!.sds  character (len=MAX_LENGTH_VAR_NAME), parameter :: O1D_VNAM = 'FULLCHEM_O1D'
!.sds  character (len=MAX_LENGTH_VAR_NAME), parameter :: OH_VNAM  = 'FULLCHEM_OH'
!.sds  character (len=MAX_LENGTH_VAR_NAME), parameter :: HO2_VNAM = 'FULLCHEM_HO2'
!.sds  character (len=MAX_LENGTH_VAR_NAME), parameter :: LON_VNAM = 'longitude'
!.sds  character (len=MAX_LENGTH_VAR_NAME), parameter :: LAT_VNAM = 'latitude'
!.sds  character (len=MAX_LENGTH_VAR_NAME), parameter :: PRS_VNAM = 'pres'
!.sds!
!.sds  real, allocatable :: lon_so3_in (:)
!.sds  real, allocatable :: lat_so3_in (:)
!.sds  real, allocatable :: pressure_so3_in (:)
!.sds  real, allocatable :: o3_in (:,:,:)
!.sds  real, allocatable :: o1d_in (:,:,:)
!.sds  real, allocatable :: oh_in (:,:,:)
!.sds  real, allocatable :: ho2_in (:,:,:)
!.sds!... work array
!.sds  real*8, allocatable :: lon_so3 (:)
!.sds  real*8, allocatable :: lat_so3 (:)
!.sds  real*8, allocatable :: logpressure_so3 (:)
!.sds  real*8, allocatable :: temp (:,:,:)
!.sds!
!.sds! !REVISION HISTORY:
!.sds!  Initial code.
!.sds!
!.sds!
!.sds!... monthly file
!.sds  if(month .ne. month_save) then
!.sds    call Ncop_Rd (ncid, self%SO3monthly_infile_name)
!.sds!
!.sds    call Ncget_Dimlen (ncid, LON_DNAM, ii_num)
!.sds    call Ncget_Dimlen (ncid, LAT_DNAM, jj_num)
!.sds    call Ncget_Dimlen (ncid, PRES_DNAM, kk_num)
!.sds!
!.sds!... get latitudes of input data
!.sds    strt1d(1) = 1
!.sds!
!.sds    cnt1d(1) = ii_num
!.sds    allocate (lon_so3_in(ii_num))
!.sds    ierr = Nf_Inq_Varid (ncid, LON_VNAM, varid)
!.sds    ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, lon_so3_in)
!.sds!
!.sds    cnt1d(1) = jj_num
!.sds    allocate (lat_so3_in(jj_num))
!.sds    ierr = Nf_Inq_Varid (ncid, LAT_VNAM, varid)
!.sds    ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, lat_so3_in)
!.sds!
!.sds!... get vertical coordinates of input data
!.sds    cnt1d(1) = kk_num
!.sds    allocate (pressure_so3_in (jj_num))
!.sds    ierr = Nf_Inq_Varid (ncid, PRS_VNAM, varid)
!.sds    ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, pressure_so3_in)
!.sds!
!.sds!... get 3D coordinates of input data
!.sds    strt4d = (/1,1,1,month/)
!.sds!
!.sds    cnt4d = (/ii_num,jj_num,kk_num,1/)
!.sds!
!.sds!... O1D monthly average field
!.sds    allocate (o1d_in (ii_num,jj_num,kk_num))
!.sds    ierr = Nf_Inq_Varid (ncid, O1D_VNAM, varid)
!.sds    ierr = Nf_Get_Vara_Real (ncid, varid, strt4d, cnt4d, o1d_in)
!.sds!
!.sds!... OH monthly average field
!.sds    allocate (oh_in (ii_num,jj_num,kk_num))
!.sds    ierr = Nf_Inq_Varid (ncid, OH_VNAM, varid)
!.sds    ierr = Nf_Get_Vara_Real (ncid, varid, strt4d, cnt4d, oh_in)
!.sds!
!.sds!... HO2 monthly average field
!.sds    allocate (ho2_in (ii_num,jj_num,kk_num))
!.sds    ierr = Nf_Inq_Varid (ncid, HO2_VNAM, varid)
!.sds    ierr = Nf_Get_Vara_Real (ncid, varid, strt4d, cnt4d, ho2_in)
!.sds!
!.sds    call Nccl (ncid)
!.sds!
!.sds!... Need to interpolate 3D fields to current model grid
!.sds    allocate (lon_so3 (ii_num+2))
!.sds    lon_so3(2:ii_num+1) = lon_so3_in (:)
!.sds!... need to set end points to wrap
!.sds    lon_so3(1) = lon_so3_in(1) - (lon_so3_in(2) - lon_so3_in(1))
!.sds    lon_so3(ii_num+2) = lon_so3_in(1) + (lon_so3_in(ii_num) - lon_so3_in(ii_num-1)) + 360
!.sds!
!.sds!... Need to interpolate 3D fields to current model grid
!.sds    allocate (lat_so3 (jj_num))
!.sds    lat_so3(:) = lat_so3_in (:)
!.sds!... need to set end points to make sure pole is covered
!.sds    lat_so3(1) = -91
!.sds    lat_so3(jj_num) = 91
!.sds!
!.sds    allocate (logpressure_so3 (kk_num))
!.sds    logpressure_so3(:) = log(pressure_so3_in (:))
!.sds!
!.sds    allocate (temp (ii_num+2,jj_num,kk_num))
!.sds!
!.sds!... loop over all 3d input fields and put on model grid
!.sds    do inum = 1,3
!.sds!
!.sds      select case (inum)
!.sds        case (1)
!.sds          temp(2:ii_num+1,:,:) = O1D_in(:,:,:)
!.sds          temp(1,:,:)          = O1D_in(ii_num,:,:)
!.sds          temp(ii_num+2,:,:)   = O1D_in(1,:,:)
!.sds        case (2)
!.sds          temp(2:ii_num+1,:,:) = OH_in(:,:,:)
!.sds          temp(1,:,:)          = OH_in(ii_num,:,:)
!.sds          temp(ii_num+2,:,:)   = OH_in(1,:,:)
!.sds        case (3)
!.sds          temp(2:ii_num+1,:,:) = HO2_in(:,:,:)
!.sds          temp(1,:,:)          = HO2_in(ii_num,:,:)
!.sds          temp(ii_num+2,:,:)   = HO2_in(1,:,:)
!.sds        case default
!.sds      end select
!.sds!
!.sds      do ik = k1, k2
!.sds        do ij = ju1, j2
!.sds          do ii = i1, i2
!.sds!.. log pressure of box
!.sds            alogplev = log(press3c(ii,ij,ik))
!.sds!
!.sds!... interpolate 3D so3 input to current 3D grid
!.sds             call Interp_Trilinear(londeg(ii), latdeg(ij), alogplev, temp_out,       &
!.sds                     lon_so3, lat_so3, logpressure_so3, ii_num+2, jj_num, kk_num, temp)
!.sds!
!.sds!... put result in appropriate array
!.sds             select case (inum)
!.sds                case (1)
!.sds                  self%monthly_o1d(ii,ij,ik) = temp_out
!.sds                case (2)
!.sds                  self%monthly_oh(ii,ij,ik) = temp_out
!.sds                case (3)
!.sds                  self%monthly_ho2(ii,ij,ik) = temp_out
!.sds                case default
!.sds              end select
!.sds!
!.sds          enddo
!.sds        enddo
!.sds      enddo
!.sds!
!.sds    enddo
!.sds!... deallocate work variables
!.sds    deallocate (lon_so3_in)
!.sds    deallocate (lat_so3_in)
!.sds    deallocate (pressure_so3_in)
!.sds    deallocate (o1d_in)
!.sds    deallocate (oh_in)
!.sds    deallocate (ho2_in)
!.sds    deallocate (lon_so3)
!.sds    deallocate (lat_so3)
!.sds    deallocate (logpressure_so3)
!.sds    deallocate (temp)
!.sds!... end monthly file
!.sds  endif
!.sds!
!.sds!... daily file
!.sds  call Ncop_Rd (ncid, self%SO3daily_infile_name)
!.sds!
!.sds  call Ncget_Dimlen (ncid, LON_DNAM, ii_num)
!.sds  call Ncget_Dimlen (ncid, LAT_DNAM, jj_num)
!.sds  call Ncget_Dimlen (ncid, PRES_DNAM, kk_num)
!.sds!
!.sds!. get latitudes of input data
!.sds  strt1d(1) = 1
!.sds!
!.sds  cnt1d(1) = ii_num
!.sds  allocate (lon_so3_in(ii_num))
!.sds  ierr = Nf_Inq_Varid (ncid, LON_VNAM, varid)
!.sds  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, lon_so3_in)
!.sds!
!.sds  cnt1d(1) = jj_num
!.sds  allocate (lat_so3_in(jj_num))
!.sds  ierr = Nf_Inq_Varid (ncid, LAT_VNAM, varid)
!.sds  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, lat_so3_in)
!.sds!
!.sds!. get vertical coordinates of input data
!.sds  cnt1d(1) = kk_num
!.sds  allocate (pressure_so3_in (jj_num))
!.sds  ierr = Nf_Inq_Varid (ncid, PRS_VNAM, varid)
!.sds  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, pressure_so3_in)
!.sds!
!.sds!. get 3D coordinates of input data
!.sds  strt5d = (/1,1,1,day,month/)
!.sds!
!.sds  cnt5d = (/ii_num,jj_num,kk_num,1,1/)
!.sds!
!.sds!. O3 daily average field
!.sds  allocate (o3_in (ii_num,jj_num,kk_num))
!.sds  ierr = Nf_Inq_Varid (ncid, O3_VNAM, varid)
!.sds  ierr = Nf_Get_Vara_Real (ncid, varid, strt5d, cnt5d, o3_in)
!.sds!
!.sds  call Nccl (ncid)
!.sds!
!.sds!
!.sds!... Need to interpolate 3D fields to current model grid
!.sds  allocate (lon_so3 (ii_num+2))
!.sds  lon_so3(2:ii_num+1) = lon_so3_in (:)
!.sds!... need to set end points to wrap
!.sds  lon_so3(1) = lon_so3_in(1) - (lon_so3_in(2) - lon_so3_in(1))
!.sds  lon_so3(ii_num+2) = lon_so3_in(1) + (lon_so3_in(ii_num) - lon_so3_in(ii_num-1)) + 360
!.sds!
!.sds!... Need to interpolate 3D fields to current model grid
!.sds  allocate (lat_so3 (jj_num))
!.sds  lat_so3(:) = lat_so3_in (:)
!.sds!... need to set end points to make sure pole is covered
!.sds  lat_so3(1) = -91
!.sds  lat_so3(jj_num) = 91
!.sds!
!.sds  allocate (logpressure_so3 (kk_num))
!.sds  logpressure_so3(:) = log(pressure_so3_in (:))
!.sds!
!.sds  allocate (temp (ii_num+2,jj_num,kk_num))
!.sds!
!.sds!... put 3d input field on model grid
!.sds  temp(2:ii_num+1,:,:) = O3_in(:,:,:)
!.sds  temp(1,:,:)          = O3_in(ii_num,:,:)
!.sds  temp(ii_num+2,:,:)   = O3_in(1,:,:)
!.sds!
!.sds  do ik = k1, k2
!.sds    do ij = ju1, j2
!.sds      do ii = i1, i2
!.sds!... log pressure of box
!.sds        alogplev = log(press3c(ii,ij,ik))
!.sds!
!.sds!. interpolate 3D so3 input to current 3D grid
!.sds        call Interp_Trilinear(londeg(ii), latdeg(ij), alogplev, temp_out,       &
!.sds                 lon_so3, lat_so3, logpressure_so3, ii_num+2, jj_num, kk_num, temp)
!.sds!
!.sds!.put result in appropriate array
!.sds        self%daily_o3(ii,ij,ik) = temp_out
!.sds!
!.sds      enddo
!.sds    enddo
!.sds  enddo
!.sds!
!.sds!
!.sds  deallocate (lon_so3_in)
!.sds  deallocate (lat_so3_in)
!.sds  deallocate (pressure_so3_in)
!.sds  deallocate (o3_in)
!.sds  deallocate (lon_so3)
!.sds  deallocate (lat_so3)
!.sds  deallocate (logpressure_so3)
!.sds  deallocate (temp)
!.sds!
!.sds  return
!.sds!
!.sds  end subroutine Initialize_stratO3
!.sds!
!.sds!-------------------------------------------------------------------------
!.sds!BOP
!.sds!
!.sds! !ROUTINE: Update_stratO3
!.sds!
!.sds! !INTERFACE:
!.sds!
!.sds!-------------------------------------------------------------------------
!.sds!
!.sds  subroutine Update_stratO3  &
!.sds   (SO3daily_infile_name, SO3monthly_infile_name, Istrat_O3, Ie90,  &
!.sds    concentration, humidity, kel, press3c, nymd, tdt, londeg, latdeg,   &
!.sds    pr_diag, procID, i1, i2, ju1, j2, k1, k2,  &
!.sds    ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl, num_species)
!.sds!
!.sds!
!.sds! !USES:
!.sds  use GmiArrayBundlePointer_mod    , only : t_GmiArrayBundle
!.sds  use GmiTimeControl_mod, only : GmiSplitDateTime, GetSecondsFromJanuary1
!.sds!
!.sds  implicit none
!.sds!
!.sds# include "gmi_phys_constants.h"
!.sds# include "gmi_time_constants.h"
!.sds# include "GmiParameters.h"
!.sds!# include "gmi_time_constants.h"
!.sds!
!.sds! !DESCRIPTION:
!.sds!   This routine updates the "Stratospheric Ozone" "tropospheric chemistry".
!.sds!
!.sds! !REVISION HISTORY:
!.sds!   Initial code: 29 August 2012
!.sds!EOP
!.sds!-------------------------------------------------------------------------
!.sds!
!.sds!
!.sds!INPUT PARAMETERS:
!.sds!
!.sds  character (len=MAX_LENGTH_FILE_NAME),intent(in) :: SO3daily_infile_name, SO3monthly_infile_name
!.sds  integer, intent(in) :: Istrat_O3
!.sds  integer, intent(in) :: Ie90
!.sds  real*8 , intent(in) :: humidity(i1:i2,   ju1:j2,   k1:k2)
!.sds  real*8,  intent(in) :: kel(ilo:ihi,julo:jhi,k1:k2)
!.sds  real*8,  intent(in) :: press3c (ilo:ihi, julo:jhi, k1:k2) ! atmospheric pressure at the center of each grid box (mb)
!.sds  integer, intent(in) :: nymd
!.sds  real*8,  intent(in) :: tdt
!.sds  real*8,  intent(in) :: londeg(i1_gl:i2_gl)
!.sds  real*8,  intent(in) :: latdeg(ju1_gl:j2_gl)
!.sds  logical, intent(in) :: pr_diag
!.sds  integer, intent(in) :: procID
!.sds!
!.sds  integer, intent(in) :: i1, i2, ju1, j2, k1, k2
!.sds  integer, intent(in) :: ilo, ihi, julo, jhi
!.sds  integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
!.sds  integer, intent(in) :: num_species
!.sds!
!.sds!
!.sds! !INPUT/OUTPUT PARAMETERS:
!.sds  type(t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!.sds!
!.sds!
!.sds! !LOCAL VARIABLES:
!.sds  type(t_StratO3), save :: self
!.sds!
!.sds  integer :: i, j, k
!.sds!
!.sds  integer :: day, month, idumyear
!.sds  integer, save :: month_save = -999
!.sds  integer, save :: day_save = -999
!.sds  logical, save :: first = .true.
!.sds  real*8,  save :: fac
!.sds!
!.sds  real*8, parameter :: k_o1d = 1.200D-10
!.sds  real*8, allocatable :: k_h2o (:,:,:)
!.sds  real*8, allocatable :: k_oh  (:,:,:)
!.sds  real*8, allocatable :: k_ho2 (:,:,:)
!.sds  real*8, allocatable :: o1d_term (:,:,:)
!.sds  real*8, allocatable :: h2o_term (:,:,:)
!.sds  real*8, allocatable :: oh_term  (:,:,:)
!.sds  real*8, allocatable :: ho2_term (:,:,:)
!.sds  real*8, allocatable :: mgas (:,:,:)
!.sds!
!.sds!
!.sds!     ----------------
!.sds!     Begin execution.
!.sds!     ----------------
!.sds!
!.sds  if (pr_diag) Write (6,*) 'Update_stratO3 called by ', procID
!.sds!
!.sds!
!.sds!... get date info
!.sds  call GmiSplitDateTime (nymd, idumyear, month, day)
!.sds!
!.sds!... need to read in appropriate month's climo values
!.sds  if(day_save .ne. day) then
!.sds    if(first) then
!.sds!... convert specific humidity to vol mixing ratio
!.sds      fac = MWTAIR / (MWTH2O * GPKG)
!.sds!...
!.sds      self%SO3daily_infile_name = SO3daily_infile_name
!.sds      self%SO3monthly_infile_name = SO3monthly_infile_name
!.sds      call Allocate_daily_o3 (self, i1, i2, ju1, j2, k1, k2)
!.sds      call Allocate_monthly_o1d (self, i1, i2, ju1, j2, k1, k2)
!.sds      call Allocate_monthly_oh  (self, i1, i2, ju1, j2, k1, k2)
!.sds      call Allocate_monthly_ho2 (self, i1, i2, ju1, j2, k1, k2)
!.sds      first = .false.
!.sds    endif
!.sds!
!.sds    call Initialize_StratO3(self, month, day, month_save, londeg, latdeg, press3c,  &
!.sds      i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl)
!.sds!
!.sds    day_save = day
!.sds    month_save = month
!.sds!
!.sds  endif
!.sds!
!.sds!... local work arrays
!.sds  Allocate(k_h2o(i1:i2, ju1:j2, k1:k2))
!.sds  k_h2o(i1:i2, ju1:j2, k1:k2) = 0.0d0
!.sds  Allocate(k_oh(i1:i2, ju1:j2, k1:k2))
!.sds  k_oh(i1:i2, ju1:j2, k1:k2) = 0.0d0
!.sds  Allocate(k_ho2(i1:i2, ju1:j2, k1:k2))
!.sds  k_ho2(i1:i2, ju1:j2, k1:k2) = 0.0d0
!.sds!
!.sds  Allocate(o1d_term(i1:i2, ju1:j2, k1:k2))
!.sds  o1d_term(i1:i2, ju1:j2, k1:k2) = 0.0d0
!.sds  Allocate(h2o_term(i1:i2, ju1:j2, k1:k2))
!.sds  h2o_term(i1:i2, ju1:j2, k1:k2) = 0.0d0
!.sds  Allocate(oh_term (i1:i2, ju1:j2, k1:k2))
!.sds  oh_term(i1:i2, ju1:j2, k1:k2) = 0.0d0
!.sds  Allocate(ho2_term(i1:i2, ju1:j2, k1:k2))
!.sds  ho2_term(i1:i2, ju1:j2, k1:k2) = 0.0d0
!.sds  Allocate(mgas(i1:i2, ju1:j2, k1:k2))
!.sds  mgas(i1:i2, ju1:j2, k1:k2) = 0.0d0
!.sds!
!.sds!... fix stratosphere to input field (daily output from fullchem model
!.sds!...  - tropopause is e90=90ppb (80.6+-0.5% of total mass) (Prather)
!.sds  where(concentration(Ie90)%pArray3D(i1:i2,ju1:j2,k1:k2) <= 90.0d-9 )
!.sds!
!.sds    concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) = self%daily_o3(i1:i2,ju1:j2,k1:k2)
!.sds!
!.sds!... tropospheric loss
!.sds  elsewhere
!.sds    mgas(i1:i2,ju1:j2,k1:k2) =   &
!.sds      ((press3c(i1:i2,ju1:j2,k1:k2) * MB2CGS) / (kel (i1:i2,ju1:j2,k1:k2) * BOLTZMN_E))
!.sds!... convert mixing ratio to molecules/cm^3
!.sds    concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) =  &
!.sds      concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) * mgas(i1:i2,ju1:j2,k1:k2)
!.sds!
!.sds!... loss terms
!.sds!.... O1D + O3 = 2 O2 and O1D + O3 = 2O + O2 (same rates)
!.sds    o1d_term(i1:i2,ju1:j2,k1:k2) = k_o1d  &
!.sds      * concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2)  &
!.sds      * self%monthly_o1d(i1:i2,ju1:j2,k1:k2)  &
!.sds      * mgas(i1:i2,ju1:j2,k1:k2)
!.sds!
!.sds!.... H2O + O1D = 2 OH
!.sds    k_h2o(i1:i2,ju1:j2,k1:k2) = 1.630D-10 * exp(-60.0D+00 / kel(i1:i2,ju1:j2,k1:k2))
!.sds    h2o_term(i1:i2,ju1:j2,k1:k2) = k_h2o(i1:i2,ju1:j2,k1:k2) &
!.sds      * (humidity(i1:i2,ju1:j2,k1:k2) * fac)  &
!.sds      * self%monthly_o1d(i1:i2,ju1:j2,k1:k2)  &
!.sds      * mgas(i1:i2,ju1:j2,k1:k2)
!.sds!
!.sds!.... O3 + OH = HO2 + O2
!.sds    k_oh(i1:i2,ju1:j2,k1:k2) = 1.700D-12 * exp(-940.0D+00 / kel(i1:i2,ju1:j2,k1:k2))
!.sds    oh_term(i1:i2,ju1:j2,k1:k2) = k_oh(i1:i2,ju1:j2,k1:k2)  &
!.sds      * concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2)  &
!.sds      * self%monthly_oh(i1:i2,ju1:j2,k1:k2)  &
!.sds      * mgas(i1:i2,ju1:j2,k1:k2)
!.sds!
!.sds!.... HO2 + O3 = 2 O2 + OH
!.sds    k_ho2(i1:i2,ju1:j2,k1:k2) = 1.000D-14 * exp(-490.0D+00 / kel(i1:i2,ju1:j2,k1:k2))
!.sds    ho2_term(i1:i2,ju1:j2,k1:k2) = k_ho2(i1:i2,ju1:j2,k1:k2)  &
!.sds      * concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2)  &
!.sds      * self%monthly_ho2(i1:i2,ju1:j2,k1:k2)  &
!.sds      * mgas(i1:i2,ju1:j2,k1:k2)
!.sds!
!.sds!... apply loss to field
!.sds    concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) = &
!.sds      concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) &
!.sds      - ( o1d_term(i1:i2,ju1:j2,k1:k2) + h2o_term(i1:i2,ju1:j2,k1:k2) &
!.sds      + oh_term(i1:i2,ju1:j2,k1:k2)  + ho2_term(i1:i2,ju1:j2,k1:k2) ) * tdt
!.sds!...
!.sds!... convert from molecules to mixing ratio
!.sds    concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) =  &
!.sds      concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) / mgas(i1:i2,ju1:j2,k1:k2)
!.sds!
!.sds  end where
!.sds!
!.sds  deallocate(o1d_term)
!.sds  deallocate(h2o_term)
!.sds  deallocate(oh_term )
!.sds  deallocate(ho2_term)
!.sds  deallocate(mgas)
!.sds!
!.sds  return
!.sds!
!.sds  end subroutine Update_stratO3
!.sds!EOC
!.sds!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialize_stratO3
!
! !INTERFACE:
!
!
  subroutine Initialize_stratO3 (self, month, day, month_save, londeg, latdeg, press3c,  &
   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl)
!
! USES:
  use m_netcdf_io_open,       only : Ncop_Rd
  use m_netcdf_io_close,      only : Nccl
  use m_netcdf_io_get_dimlen, only : Ncget_Dimlen
  use GmiInterpolation_mod,   only : Interp, Interp_Trilinear
!   read in parameters necessary for Strat_O3 calculation
!
  implicit none
!
# include "netcdf.inc"
# include "GmiParameters.h"
!
! !INPUT PARAMETERS:
  integer, intent(in) :: month, day, month_save
  real*8 , intent(in) :: londeg(i1_gl:i2_gl)
  real*8 , intent(in) :: latdeg(ju1_gl:j2_gl)
  real*8,  intent(in) :: press3c (ilo:ihi, julo:jhi, k1:k2) ! atmospheric pressure at the center of each grid box (mb)
  integer, intent(in) :: i1, i2, ju1, j2, k1, k2
  integer, intent(in) :: ilo, ihi, julo, jhi
  integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
!
! !OUTPUT PARAMETERS:
!
  type (t_stratO3), intent(inout)   :: self
!
! !DESCRIPTION:
!  Open a NetCDF file for reading of Strat_O3 parameters and interpolate to current grid
!
! !LOCAL VARIABLES:
  character (len=MAX_LENGTH_VAR_NAME) :: err_msg
  integer :: ierr
  integer :: ncid
  integer :: ii_num, jj_num, kk_num
  integer :: cnt1d (1)
  integer :: strt1d(1)
  integer :: cnt2d (2)
  integer :: strt2d(2)
  integer :: cnt3d (3)
  integer :: strt3d(3)
  integer :: cnt4d (4)
  integer :: strt4d(4)
  integer :: cnt5d (5)
  integer :: strt5d(5)
  integer :: varid
  integer :: ii, ij, ik, inum
  real*8  :: temp_out
  real*8  :: alogplev
!
!
! !AUTHOR:
!  Stephen Steenrod (GSFC/USRA)
!
!
!... dimension parameter names for Strat_O3 input files
  character (len=MAX_LENGTH_VAR_NAME), parameter :: LON_DNAM  = 'longitude'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: LAT_DNAM  = 'latitude'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: PRES_DNAM = 'pres'
!... variable parameter names for Strat_O3 input files
  character (len=MAX_LENGTH_VAR_NAME), parameter :: O3_VNAM  = 'FULLCHEM_O3'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: O1D_VNAM = 'FULLCHEM_k_O1D_H2O'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: OH_VNAM  = 'FULLCHEM_k_OH_O3'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: HO2_VNAM = 'FULLCHEM_k_HO2_O3'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: LON_VNAM = 'longitude'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: LAT_VNAM = 'latitude'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: PT_VNAM  = 'pt'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: AM_VNAM  = 'am'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: BM_VNAM  = 'bm'
  character (len=MAX_LENGTH_VAR_NAME), parameter :: PRS_VNAM = 'psfc'
!
  real, allocatable :: lon_so3_in (:)
  real, allocatable :: lat_so3_in (:)
  real :: pt_so3_in
  real, allocatable :: am_so3_in (:)
  real, allocatable :: bm_so3_in (:)
  real, allocatable :: psfc_so3_in (:,:)
  real, allocatable :: o3_in (:,:,:)
  real, allocatable :: o1d_in (:,:,:)
  real, allocatable :: oh_in (:,:,:)
  real, allocatable :: ho2_in (:,:,:)
!... work array
  real*8, allocatable :: lon_so3 (:)
  real*8, allocatable :: lat_so3 (:)
  real*8, allocatable :: logpressure_so3 (:)
  real*8, allocatable :: temp (:,:,:)
!
! !REVISION HISTORY:
!  Initial code.
!
!
!... monthly file
  if(month .ne. month_save) then
    call Ncop_Rd (ncid, self%SO3monthly_infile_name)
!
    call Ncget_Dimlen (ncid, LON_DNAM, ii_num)
    call Ncget_Dimlen (ncid, LAT_DNAM, jj_num)
    call Ncget_Dimlen (ncid, PRES_DNAM, kk_num)
!
!... get latitudes of input data
    strt1d(1) = 1
!
    cnt1d(1) = ii_num
    allocate (lon_so3_in(ii_num))
    ierr = Nf_Inq_Varid (ncid, LON_VNAM, varid)
    ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, lon_so3_in)
!
    cnt1d(1) = jj_num
    allocate (lat_so3_in(jj_num))
    ierr = Nf_Inq_Varid (ncid, LAT_VNAM, varid)
    ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, lat_so3_in)
!
!... get a's, b's, pt and sfc press to calc vertical coordinates of input data
!
    cnt1d(1) = 1
    ierr = Nf_Inq_Varid (ncid, PT_VNAM, varid)
    ierr = Nf_Get_Var1_Real (ncid, varid, strt1d, pt_so3_in)
!
    cnt1d(1) = kk_num
    allocate (am_so3_in(kk_num))
    ierr = Nf_Inq_Varid (ncid, AM_VNAM, varid)
    ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, am_so3_in)
!
    cnt1d(1) = kk_num
    allocate (bm_so3_in(kk_num))
    ierr = Nf_Inq_Varid (ncid, BM_VNAM, varid)
    ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, bm_so3_in)
!
    strt3d = (/1,1,month/)
    cnt3d = (/ii_num,jj_num,1/)
    allocate (psfc_so3_in (ii_num, jj_num))
    ierr = Nf_Inq_Varid (ncid, PRS_VNAM, varid)
    ierr = Nf_Get_Vara_Real (ncid, varid, strt3d, cnt3d, psfc_so3_in)
!
!... get 3D coordinates of input data
    strt4d = (/1,1,1,month/)
!
    cnt4d = (/ii_num,jj_num,kk_num,1/)
!
!... O1D monthly average field
    allocate (o1d_in (ii_num,jj_num,kk_num))
    ierr = Nf_Inq_Varid (ncid, O1D_VNAM, varid)
    ierr = Nf_Get_Vara_Real (ncid, varid, strt4d, cnt4d, o1d_in)
!
!... OH monthly average field
    allocate (oh_in (ii_num,jj_num,kk_num))
    ierr = Nf_Inq_Varid (ncid, OH_VNAM, varid)
    ierr = Nf_Get_Vara_Real (ncid, varid, strt4d, cnt4d, oh_in)
!
!... HO2 monthly average field
    allocate (ho2_in (ii_num,jj_num,kk_num))
    ierr = Nf_Inq_Varid (ncid, HO2_VNAM, varid)
    ierr = Nf_Get_Vara_Real (ncid, varid, strt4d, cnt4d, ho2_in)
!
    call Nccl (ncid)
!
!... Need to interpolate 3D fields to current model grid
    allocate (lon_so3 (ii_num+2))
    lon_so3(2:ii_num+1) = lon_so3_in (:)
!... need to set end points to wrap
    lon_so3(1) = lon_so3_in(1) - (lon_so3_in(2) - lon_so3_in(1))
    lon_so3(ii_num+2) = lon_so3_in(1) + (lon_so3_in(ii_num) - lon_so3_in(ii_num-1)) + 360
!
!... Need to interpolate 3D fields to current model grid
    allocate (lat_so3 (jj_num))
    lat_so3(:) = lat_so3_in (:)
!... need to set end points to make sure pole is covered
    lat_so3(1) = -91
    lat_so3(jj_num) = 91
!
    allocate (logpressure_so3 (kk_num))
!
    allocate (temp (ii_num+2,jj_num,kk_num))
!
!... loop over all 3d input fields and put on model grid
    do inum = 1,3
!... no interpolation needed
      if(ii_num.eq.(i2_gl-i1_gl+1).and.jj_num.eq.(j2_gl-ju1_gl+1).and.kk_num.eq.(k2-k1+1)) then
        select case (inum)
          case (1)
            self%monthly_o1d(:,:,:) = O1D_in(i1:i2,ju1:j2,k1:k2)
          case (2)
            self%monthly_oh(:,:,:) = OH_in(i1:i2,ju1:j2,k1:k2)
          case (3)
            self%monthly_ho2(:,:,:) = HO2_in(i1:i2,ju1:j2,k1:k2)
          case default
        end select
!
      else
!... interpolate 
        select case (inum)
          case (1)
            temp(2:ii_num+1,:,:) = O1D_in(:,:,:)
            temp(1,:,:)          = O1D_in(ii_num,:,:)
            temp(ii_num+2,:,:)   = O1D_in(1,:,:)
          case (2)
            temp(2:ii_num+1,:,:) = OH_in(:,:,:)
            temp(1,:,:)          = OH_in(ii_num,:,:)
            temp(ii_num+2,:,:)   = OH_in(1,:,:)
          case (3)
            temp(2:ii_num+1,:,:) = HO2_in(:,:,:)
            temp(1,:,:)          = HO2_in(ii_num,:,:)
            temp(ii_num+2,:,:)   = HO2_in(1,:,:)
          case default
        end select
!
        do ik = k1, k2
          do ij = ju1, j2
            do ii = i1, i2
!.. log pressure of box
                alogplev = log(press3c(ii,ij,ik))
!.. log pressure of input data
                logpressure_so3 = log(pt_so3_in*am_so3_in(ik) + bm_so3_in(ik)*psfc_so3_in(ii,ij))
!
!... interpolate 3D so3 input to current 3D grid
                 call Interp_Trilinear(londeg(ii), latdeg(ij), alogplev, temp_out,       &
                       lon_so3, lat_so3, logpressure_so3, ii_num+2, jj_num, kk_num, temp)
!
!... put result in appropriate array
               select case (inum)
                  case (1)
                    self%monthly_o1d(ii,ij,ik) = temp_out
                  case (2)
                    self%monthly_oh(ii,ij,ik) = temp_out
                  case (3)
                    self%monthly_ho2(ii,ij,ik) = temp_out
                  case default
                end select
!
            enddo
          enddo
        enddo
      endif
!
    enddo
!... deallocate work variables
    deallocate (lon_so3_in)
    deallocate (lat_so3_in)
    deallocate (am_so3_in)
    deallocate (bm_so3_in)
    deallocate (psfc_so3_in)
    deallocate (o1d_in)
    deallocate (oh_in)
    deallocate (ho2_in)
    deallocate (lon_so3)
    deallocate (lat_so3)
    deallocate (logpressure_so3)
    deallocate (temp)
!... end monthly file
  endif
!
!... daily file
  call Ncop_Rd (ncid, self%SO3daily_infile_name)
!
  call Ncget_Dimlen (ncid, LON_DNAM, ii_num)
  call Ncget_Dimlen (ncid, LAT_DNAM, jj_num)
  call Ncget_Dimlen (ncid, PRES_DNAM, kk_num)
!
!. get latitudes of input data
  strt1d(1) = 1
!
  cnt1d(1) = ii_num
  allocate (lon_so3_in(ii_num))
  ierr = Nf_Inq_Varid (ncid, LON_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, lon_so3_in)
!
  cnt1d(1) = jj_num
  allocate (lat_so3_in(jj_num))
  ierr = Nf_Inq_Varid (ncid, LAT_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, lat_so3_in)
!
!... get a's, b's, pt and sfc press to calc vertical coordinates of input data
!
  cnt1d(1) = 1
  ierr = Nf_Inq_Varid (ncid, PT_VNAM, varid)
  ierr = Nf_Get_Var1_Real (ncid, varid, strt1d, pt_so3_in)
!
  cnt1d(1) = kk_num
  allocate (am_so3_in(kk_num))
  ierr = Nf_Inq_Varid (ncid, AM_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, am_so3_in)
!
  cnt1d(1) = kk_num
  allocate (bm_so3_in(kk_num))
  ierr = Nf_Inq_Varid (ncid, BM_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt1d, cnt1d, bm_so3_in)
!
  strt4d = (/1,1,day,month/)
  cnt4d = (/ii_num,jj_num,1,1/)
  allocate (psfc_so3_in (ii_num, jj_num))
  ierr = Nf_Inq_Varid (ncid, PRS_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt4d, cnt4d, psfc_so3_in)
!
!. get 3D coordinates of input data
  strt5d = (/1,1,1,day,month/)
!
  cnt5d = (/ii_num,jj_num,kk_num,1,1/)
!
!. O3 daily average field
  allocate (o3_in (ii_num,jj_num,kk_num))
  ierr = Nf_Inq_Varid (ncid, O3_VNAM, varid)
  ierr = Nf_Get_Vara_Real (ncid, varid, strt5d, cnt5d, o3_in)
!
  call Nccl (ncid)
!
!
!... no interpolation needed
  if(ii_num.eq.(i2_gl-i1_gl+1).and.jj_num.eq.(j2_gl-ju1_gl+1).and.kk_num.eq.(k2-k1+1)) then
    self%daily_o3(:,:,:) = O3_in(i1:i2,ju1:j2,k1:k2)
!
  else
!... interpolate 
!... Need to interpolate 3D fields to current model grid
    allocate (lon_so3 (ii_num+2))
    lon_so3(2:ii_num+1) = lon_so3_in (:)
!... need to set end points to wrap
    lon_so3(1) = lon_so3_in(1) - (lon_so3_in(2) - lon_so3_in(1))
    lon_so3(ii_num+2) = lon_so3_in(1) + (lon_so3_in(ii_num) - lon_so3_in(ii_num-1)) + 360
!
!... Need to interpolate 3D fields to current model grid
    allocate (lat_so3 (jj_num))
    lat_so3(:) = lat_so3_in (:)
!... need to set end points to make sure pole is covered
    lat_so3(1) = -91
    lat_so3(jj_num) = 91
!
    allocate (logpressure_so3 (kk_num))
!
    allocate (temp (ii_num+2,jj_num,kk_num))
!
!... put 3d input field on model grid
    temp(2:ii_num+1,:,:) = O3_in(:,:,:)
    temp(1,:,:)          = O3_in(ii_num,:,:)
    temp(ii_num+2,:,:)   = O3_in(1,:,:)
!
    do ik = k1, k2
      do ij = ju1, j2
        do ii = i1, i2
!... log pressure of box
          alogplev = log(press3c(ii,ij,ik))
!... log pressure of input data
          logpressure_so3 = log(pt_so3_in*am_so3_in(ik) + bm_so3_in(ik)*psfc_so3_in(ii,ij))
!
!... interpolate 3D so3 input to current 3D grid
          call Interp_Trilinear(londeg(ii), latdeg(ij), alogplev, temp_out,       &
                   lon_so3, lat_so3, logpressure_so3, ii_num+2, jj_num, kk_num, temp)
!
!...  put result in appropriate array
          self%daily_o3(ii,ij,ik) = temp_out
!
        enddo
      enddo
    enddo
    deallocate (lon_so3)
    deallocate (lat_so3)
    deallocate (logpressure_so3)
    deallocate (temp)
  endif
!
!
  deallocate (lon_so3_in)
  deallocate (lat_so3_in)
  deallocate (am_so3_in)
  deallocate (bm_so3_in)
  deallocate (psfc_so3_in)
  deallocate (o3_in)
!
  return
!
  end subroutine Initialize_stratO3
!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update_stratO3
!
! !INTERFACE:
!
!-------------------------------------------------------------------------
!
  subroutine Update_stratO3  &
   (SO3daily_infile_name, SO3monthly_infile_name, Istrat_O3, Ie90,  &
    concentration, humidity, kel, press3c, nymd, tdt, londeg, latdeg,   &
    pr_diag, procID, i1, i2, ju1, j2, k1, k2,  &
    ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl, num_species)
!
!
! !USES:
  use GmiArrayBundlePointer_mod    , only : t_GmiArrayBundle
  use GmiTimeControl_mod, only : GmiSplitDateTime, GetSecondsFromJanuary1
!
  implicit none
!
# include "gmi_phys_constants.h"
# include "gmi_time_constants.h"
# include "GmiParameters.h"
!# include "gmi_time_constants.h"
!
! !DESCRIPTION:
!   This routine updates the "Stratospheric Ozone" "tropospheric chemistry".
!
! !REVISION HISTORY:
!   Initial code: 29 August 2012
!EOP
!-------------------------------------------------------------------------
!
!
!INPUT PARAMETERS:
!
  character (len=MAX_LENGTH_FILE_NAME),intent(in) :: SO3daily_infile_name, SO3monthly_infile_name
  integer, intent(in) :: Istrat_O3
  integer, intent(in) :: Ie90
  real*8 , intent(in) :: humidity(i1:i2,   ju1:j2,   k1:k2)
  real*8,  intent(in) :: kel(ilo:ihi,julo:jhi,k1:k2)
  real*8,  intent(in) :: press3c (ilo:ihi, julo:jhi, k1:k2) ! atmospheric pressure at the center of each grid box (mb)
  integer, intent(in) :: nymd
  real*8,  intent(in) :: tdt
  real*8,  intent(in) :: londeg(i1_gl:i2_gl)
  real*8,  intent(in) :: latdeg(ju1_gl:j2_gl)
  logical, intent(in) :: pr_diag
  integer, intent(in) :: procID
!
  integer, intent(in) :: i1, i2, ju1, j2, k1, k2
  integer, intent(in) :: ilo, ihi, julo, jhi
  integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
  integer, intent(in) :: num_species
!
!
! !INPUT/OUTPUT PARAMETERS:
  type(t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!
!
! !LOCAL VARIABLES:
  type(t_StratO3), save :: self
!
  integer :: i, j, k
!
  integer :: day, month, idumyear
  integer, save :: month_save = -999
  integer, save :: day_save = -999
  logical, save :: first = .true.
  real*8,  save :: fac
  real*8,  save :: e90_thresh = 90.0d-9 
!  real*8,  save :: e90_thresh = 85.0d-9 
!  real*8,  save :: e90_thresh = 80.0d-9 
!  real*8,  save :: e90_thresh = 75.0d-9 
!
  real*8, allocatable :: loss_term (:,:,:)
!  real*8, allocatable :: mgas (:,:,:)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
  if (pr_diag) Write (6,*) 'Update_stratO3 called by ', procID
!
!
!  print *,'sdssds 00: ',procID,maxval(concentration(Istrat_O3)%pArray3D(:,:,:)),minval(concentration(Istrat_O3)%pArray3D(:,:,:))
!... get date info
  call GmiSplitDateTime (nymd, idumyear, month, day)
!
!  print *,'sdssds 01: ',procID,maxval(concentration(Istrat_O3)%pArray3D(:,:,:)),minval(concentration(Istrat_O3)%pArray3D(:,:,:))
!... need to read in appropriate month's climo values
  if(day_save .ne. day) then
    if(first) then
!... convert specific humidity to vol mixing ratio
      fac = MWTAIR / (MWTH2O * GPKG)
!...
      self%SO3daily_infile_name = SO3daily_infile_name
      self%SO3monthly_infile_name = SO3monthly_infile_name
      call Allocate_daily_o3 (self, i1, i2, ju1, j2, k1, k2)
      call Allocate_monthly_o1d (self, i1, i2, ju1, j2, k1, k2)
      call Allocate_monthly_oh  (self, i1, i2, ju1, j2, k1, k2)
      call Allocate_monthly_ho2 (self, i1, i2, ju1, j2, k1, k2)
      first = .false.
    endif
!
    call Initialize_StratO3(self, month, day, month_save, londeg, latdeg, press3c,  &
      i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl)
!
    day_save = day
    month_save = month
!
  endif
!
!  print *,'sdssds 02: ',procID,maxval(concentration(Istrat_O3)%pArray3D(:,:,:)),minval(concentration(Istrat_O3)%pArray3D(:,:,:))
!... local work arrays
  Allocate(loss_term(i1:i2, ju1:j2, k1:k2))
  loss_term(i1:i2, ju1:j2, k1:k2) = 0.0d0
!
!  Allocate(mgas(i1:i2, ju1:j2, k1:k2))
!  mgas(i1:i2, ju1:j2, k1:k2) = 0.0d0
!
!... fix stratosphere to input field (daily output from fullchem model
!...  - tropopause is e90=90ppb (80.6+-0.5% of total mass) (Prather)
!  print *,'sdssds 03: ',procID,maxval(concentration(Istrat_O3)%pArray3D(:,:,:)),minval(concentration(Istrat_O3)%pArray3D(:,:,:))
  where(concentration(Ie90)%pArray3D(i1:i2,ju1:j2,k1:k2) <= e90_thresh )
!
    concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) = self%daily_o3(i1:i2,ju1:j2,k1:k2)
!
!... tropospheric loss
  elsewhere
!
!... loss terms
!... L(O3, t0) = k1[O1D][H2O] + k2[HO2][O3] + k3[OH][O3} + other terms?
!... The above loss has units of cm-3s-1, and all the densities are those calculated at t0.
!...  We can define an "instantaneous" first-order loss time constant as:
!...   ? (t0) = [O3](t0)/L(O3,t0)
!... We then calculate O3,strat(t1) by assuming a linear loss for  O3,strat(t0) with
!...  constant ? (t0) over the interval [t0, t1].
!... We could also get ? as the average of the time constants at t0 and  t1.
!... A simple solution would be:
!...  [O3_strat](t1) = [O3_strat](t0) exp (-(t1-t0) / (O3](t0)/L(O3,t0)) )
!
!.... O1D + H2O = 2 OH
!.... HO2 + O3 = 2 O2 + OH
!.... OH + O3 = HO2 + O2
    loss_term(i1:i2,ju1:j2,k1:k2) = concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) &
      / ( self%monthly_o1d(i1:i2,ju1:j2,k1:k2) + self%monthly_ho2(i1:i2,ju1:j2,k1:k2) &
        + self%monthly_oh(i1:i2,ju1:j2,k1:k2) )
!!... calc m
!    mgas(i1:i2,ju1:j2,k1:k2) =   &
!      ((press3c(i1:i2,ju1:j2,k1:k2) * MB2CGS) / (kel (i1:i2,ju1:j2,k1:k2) * BOLTZMN_E))
!!... convert mixing ratio to molecules/cm^3
!    concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) =  &
!      concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) * mgas(i1:i2,ju1:j2,k1:k2)
!
!... apply loss to field
    concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) = &
      concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) &
       * exp(-tdt/loss_term(i1:i2,ju1:j2,k1:k2))
!!...
!!... convert from molecules cm-3 back to volume mixing ratio
!    concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) =  &
!      concentration(Istrat_O3)%pArray3D(i1:i2,ju1:j2,k1:k2) / mgas(i1:i2,ju1:j2,k1:k2)
!
  end where
!  print *,'sdssds 04: ',procID,maxval(concentration(Istrat_O3)%pArray3D(:,:,:)),minval(concentration(Istrat_O3)%pArray3D(:,:,:))
!
  deallocate(loss_term)
!  deallocate(mgas)
!
  return
!
  end subroutine Update_stratO3
!EOC
!
!-------------------------------------------------------------------------
   subroutine Allocate_O3_linoz (self, ju1, j2, k1, k2)
    implicit none
    integer       , intent(in   ) :: ju1, j2, k1, k2
    type (t_Linoz), intent(inOut) :: self
    allocate(self%O3_linoz(ju1:j2,k1:k2))
    self%O3_linoz(:,:) = 0.0d0
    return
  end subroutine Allocate_O3_linoz
!-------------------------------------------------------------------------
   subroutine Allocate_T_linoz (self, ju1, j2, k1, k2)
    implicit none
    integer       , intent(in   ) :: ju1, j2, k1, k2
    type (t_Linoz), intent(inOut) :: self
    allocate(self%T_linoz(ju1:j2,k1:k2))
    self%T_linoz(:,:) = 0.0d0
    return
  end subroutine Allocate_T_linoz
!-------------------------------------------------------------------------
   subroutine Allocate_colO3_linoz (self, ju1, j2, k1, k2)
    implicit none
    integer       , intent(in   ) :: ju1, j2, k1, k2
    type (t_Linoz), intent(inOut) :: self
    allocate(self%colO3_linoz(ju1:j2,k1:k2))
    self%colO3_linoz(:,:) = 0.0d0
    return
  end subroutine Allocate_colO3_linoz
!-------------------------------------------------------------------------
   subroutine Allocate_PmL_linoz (self, ju1, j2, k1, k2)
    implicit none
    integer       , intent(in   ) :: ju1, j2, k1, k2
    type (t_Linoz), intent(inOut) :: self
    allocate(self%PmL_linoz(ju1:j2,k1:k2))
    self%PmL_linoz(:,:) = 0.0d0
    return
  end subroutine Allocate_PmL_linoz
!-------------------------------------------------------------------------
   subroutine Allocate_dPmLdO3_linoz (self, ju1, j2, k1, k2)
    implicit none
    integer       , intent(in   ) :: ju1, j2, k1, k2
    type (t_Linoz), intent(inOut) :: self
    allocate(self%dPmLdO3_linoz(ju1:j2,k1:k2))
    self%dPmLdO3_linoz(:,:) = 0.0d0
    return
  end subroutine Allocate_dPmLdO3_linoz
!-------------------------------------------------------------------------
   subroutine Allocate_dPmLdT_linoz (self, ju1, j2, k1, k2)
    implicit none
    integer       , intent(in   ) :: ju1, j2, k1, k2
    type (t_Linoz), intent(inOut) :: self
    allocate(self%dPmLdT_linoz(ju1:j2,k1:k2))
    self%dPmLdT_linoz(:,:) = 0.0d0
    return
  end subroutine Allocate_dPmLdT_linoz
!-------------------------------------------------------------------------
   subroutine Allocate_dPmLdcolO3_linoz (self, ju1, j2, k1, k2)
    implicit none
    integer       , intent(in   ) :: ju1, j2, k1, k2
    type (t_Linoz), intent(inOut) :: self
    allocate(self%dPmLdcolO3_linoz(ju1:j2,k1:k2))
    self%dPmLdcolO3_linoz(:,:) = 0.0d0
    return
  end subroutine Allocate_dPmLdcolO3_linoz
!-------------------------------------------------------------------------
   subroutine Allocate_Cariolle_Loss_linoz (self, k1, k2)
    implicit none
    integer       , intent(in   ) :: k1, k2
    type (t_Linoz), intent(inOut) :: self
    allocate(self%Cariolle_Loss_linoz(k1:k2))
    self%Cariolle_Loss_linoz(:) = 0.0d0
    return
  end subroutine Allocate_Cariolle_Loss_linoz
!-------------------------------------------------------------------------
   subroutine Allocate_J_SF6 (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer     , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_SF6), intent(inOut) :: self
    allocate(self%j_SF6(i1:i2,ju1:j2,k1:k2))
    self%j_SF6(:,:,:) = 0.0d0
    return
  end subroutine Allocate_J_SF6
!-------------------------------------------------------------------------
   subroutine Allocate_K_SF6 (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer     , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_SF6), intent(inOut) :: self
    allocate(self%k_SF6(i1:i2,ju1:j2,k1:k2))
    self%k_SF6(:,:,:) = 0.0d0
    return
  end subroutine Allocate_K_SF6
!-------------------------------------------------------------------------
   subroutine Allocate_AD_SF6 (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer     , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_SF6), intent(inOut) :: self
    allocate(self%AD_SF6(i1:i2,ju1:j2,k1:k2))
    self%AD_SF6(:,:,:) = 0.0d0
    return
  end subroutine Allocate_AD_SF6
!-------------------------------------------------------------------------
   subroutine Allocate_O1D_SF6 (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer     , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_SF6), intent(inOut) :: self
    allocate(self%O1D_SF6(i1:i2,ju1:j2,k1:k2))
    self%O1D_SF6(:,:,:) = 0.0d0
    return
  end subroutine Allocate_O1D_SF6
!-------------------------------------------------------------------------
   subroutine Allocate_daily_o3 (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer         , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_stratO3), intent(inOut) :: self
    allocate(self%daily_o3(i1:i2,ju1:j2,k1:k2))
    self%daily_o3(:,:,:) = 0.0d0
    return
  end subroutine Allocate_daily_o3
!-------------------------------------------------------------------------
   subroutine Allocate_monthly_o1d (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer         , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_stratO3), intent(inOut) :: self
    allocate(self%monthly_o1d(i1:i2,ju1:j2,k1:k2))
    self%monthly_o1d(:,:,:) = 0.0d0
    return
  end subroutine Allocate_monthly_o1d
!-------------------------------------------------------------------------
   subroutine Allocate_monthly_oh (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer         , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_stratO3), intent(inOut) :: self
    allocate(self%monthly_oh(i1:i2,ju1:j2,k1:k2))
    self%monthly_oh(:,:,:) = 0.0d0
    return
  end subroutine Allocate_monthly_oh
!-------------------------------------------------------------------------
   subroutine Allocate_monthly_ho2 (self, i1, i2, ju1, j2, k1, k2)
    implicit none
    integer         , intent(in   ) :: i1, i2, ju1, j2, k1, k2
    type (t_stratO3), intent(inOut) :: self
    allocate(self%monthly_ho2(i1:i2,ju1:j2,k1:k2))
    self%monthly_ho2(:,:,:) = 0.0d0
    return
  end subroutine Allocate_monthly_ho2
!
 end module GmiTracerMethod_mod
!
