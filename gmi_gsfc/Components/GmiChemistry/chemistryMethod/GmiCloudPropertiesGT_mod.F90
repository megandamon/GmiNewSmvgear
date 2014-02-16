!------------------------------------------------------------------------------
! NASA GSFC - SIVO code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiCloudPropertiesGT_mod
!
! !INTERFACE:
!
      module GmiCloudPropertiesGT_mod
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use  GmiSolver_SavedVariables_mod
!      
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: calcCloudPropertiesGT
!
#     include "gmi_time_constants.h"
#     include "gmi_phys_constants.h"
!#     include "gmiCloudParametersGT.h"
!
! !DESCRIPTION:
! Provides routines to calculate cloud properties.
!
! !AUTHOR:
! Nicholas Meskhidze,       GA TECH,   nmeskhidze@eas.gatech.edu
! Rafaella P. Sotiropoulou, GA TECH,   rafaella,sotiropoulou@eas.gatech.edu
! Jules Kouatchou,          NASA GSFC, Jules.Kouatchou-1@nasa.gov
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcCloudPropertiesGT
!
! !INTERFACE:
!
      subroutine calcCloudPropertiesGT (cloudGT, londeg, latdeg, kel, concentration, &
     &               humidity, moistq, cmf1, press3c, gridBoxHeight, mass,  &
     &               lwi_flags, totalCloudFraction, tau_cloud, radswg,      &
     &               cloud_param, cloud_tau, cloudDroplet_opt, pr_diag,     &
     &               loc_proc, num_species, imgas_num, i1, i2, ju1, j2, k1, &
     &               k2, ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl)
!

      implicit none

#     include "setkin_par.h"
#     include "setkin_depos.h"
#     include "phot_lookup.h"

! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: num_species, imgas_num, loc_proc
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ilo, ihi, julo, jhi, i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: cloudDroplet_opt   
                             ! 1:BL Correlation;
                             ! 2:NS Parameterization;
                             ! 3:AG Parameterization;
                             ! 4:SK Correlation;
      real*8 , intent(in) :: londeg       (i1_gl:i2_gl)
                             ! longitude (deg)
      real*8 , intent(in) :: latdeg       (ju1_gl:j2_gl)
                             ! latitude (deg)
      real*8 , intent(in) :: press3c      (ilo:ihi, julo:jhi, k1:k2)
                             ! atmospheric pressure at the center of each grid box (mb)
      real*8 , intent(in) :: mass         (i1:i2, ju1:j2, k1:k2)
                             ! total mass of the atmosphere within each grid box (kg)
      integer, intent(in) :: lwi_flags    (i1:i2, ju1:j2)
                             ! array of flags that indicate land-2, water-1, or ice-3
      real*8 , intent(in) :: kel          (ilo:ihi, julo:jhi, k1:k2)
                             ! temperature (degK)
      real*8 , intent(in) :: humidity     (i1:i2, ju1:j2, k1:k2)
                             ! specific humidity  (g/kg)
      real*8 , intent(in) :: moistq       (i1:i2, ju1:j2, k1:k2)
                             ! moisture changes due to wet processes (g/kg/day)
      real*8 , intent(in) :: cmf1         (i1:i2,  ju1:j2, k1:k2)
                             ! convective mass flux used for cloud fraction calculation (kg/m^2/s)
      real*8 , intent(in) :: radswg       (i1:i2, ju1:j2)
                             ! net downward shortwave radiation at ground (W/m^2)
      real*8 , intent(in) :: tau_cloud    (i1:i2, ju1:j2, k1:k2)
                             ! optical depth (dimensionless)
      real*8 , intent(in) :: gridBoxHeight(i1:i2, ju1:j2, k1:k2)
                             ! height of each grid box (m)
      real*8 , intent(in) :: totalCloudFraction(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_CloudParametersGT), intent(inOut) :: cloudGT
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
                              ! species concentration, known at zone centers (mixing ratio)
      real*8                 , intent(inOut) :: cloud_param(i1:i2, ju1:j2, k1:k2,num_species)
      real*8                 , intent(inOut) :: cloud_tau  (i1:i2, ju1:j2, k1:k2)
!
! !DESCRIPTION:
! This routine calculates some cloud properties.
!
! !DEFINED PARAMETERS:
      real*8, parameter :: EXT_EFF = 2.0d0        ! extinction efficiency of a water droplet
      real*8, parameter :: LWD     = 1.0d+06      ! liquid water densitry (g/m^3)
      real*8, parameter :: PI      = 3.14159260d0
      real*8, parameter :: ONETH   = 0.33333333d0
      real*8, parameter :: THRFTH  = 0.750d0
      real*8, parameter :: ONEOM   = 1.0d-6
      real*8, parameter :: MAX_D   = 2.5d+3
      real*8, parameter :: MIN_D   = 40.0d0
      real*8, parameter :: MLN     = 1.0d+6
      real*8, parameter :: BETAcon = 0.67d0
      real*8, parameter :: BETAmar = 0.80d0
      real*8, parameter :: GRAV    = 9.81d0         ! g constant 
!
! !LOCAL VARIABLES:
      real*8  :: hl      (ju1_gl:j2_gl)            !rsot
      real*8  :: height_m(i1:i2, ju1:j2, k1:k2)    !rsot
      real*8  :: QLWC  (i1:i2, ju1:j2, k1:k2)
      real*8  :: QCRIT (i1:i2, ju1:j2, k1:k2)
      real*8  :: Rv    (i1:i2, ju1:j2, k1:k2)
      real*8  :: RHOAIR(i1:i2, ju1:j2, k1:k2)
      real*8  :: FCLD  (i1:i2, ju1:j2, k1:k2)
      real*8  :: res   (i1:i2, ju1:j2)
      real*8  :: resLWP(i1:i2, ju1:j2)
      real*8  :: afrac (k1:k2)
      real*8  :: LWPfrac(k1:k2)

      !logical, save :: first = .true.
      !logical, save :: cloud_top = .true.

      integer :: il, ij, ik
      integer :: ic, icx
      integer :: idumday, idumyear
      integer :: imon
      integer :: month
      integer :: aero_typ

      real*8  :: fac

      real*8, allocatable :: c0 (:)
      real*8, allocatable :: c1 (:)
      real*8, allocatable :: c2 (:)
      real*8, allocatable :: cfc(:)
      real*8, allocatable :: cfs(:)
      real*8, allocatable :: rhc(:)

!     -------------------------------------------
!     cmfcol     : convective mass flux   (mb/hr)
!     height_asl : height above sea level (m)
!     -------------------------------------------

      real*8, allocatable :: cmfcol    (:)
      real*8, allocatable :: height_asl(:)

!     -----------------------------------------------------
!     clfrac     : total cloud fraction [0 - 1]
!     cloud_frac : Cloud fraction calculated from maximum - random
!                  overalap scheme  
!     lwccol     : in-cloud liquid water in each grid box (g/m3)
!     lwpath     : liquid water path (g/m2)
!     lwppath    : alternative way to calculate liquid water path (g/m2)
!     relhum     : relative humidity [0 - 1]
!     -----------------------------------------------------

      real*8, allocatable :: clfrac (:,:,:)
      real*8, allocatable :: cloud_frac (:,:,:)
      real*8, allocatable :: lwccol (:,:,:)
      real*8, allocatable :: lwpath (:,:,:)
      real*8, allocatable :: lwppath(:,:,:)
      real*8, allocatable :: relhum (:,:,:)
!
!     ----------------------------------------------------
!     eff_rad     : cloud drop effective radius (m)
!     eff_rad_n   : Cloud Effective radius calculated from
!                   total SO4 (m)
!     eff_rad_par : Cloud effective radius calculated from
!                   Fountoukis & Nenes, 2005
!     ----------------------------------------------------
      real*8, allocatable :: eff_rad      (:,:,:)
      real*8, allocatable :: eff_rad_n    (:,:,:)
      real*8, allocatable :: eff_rad_par  (:,:,:)
!
!     ----------------------------------------------
!     cloud_opt_dep     : cloud optical depth (unitless)
!     cloud_depth       : depth of cloud (m)
!     cloud_opt_dep_n   : new cloud optical depth calculated
!                         using efffective radius as a function
!                         of total SO4 (unitless)
!     cloud_opt_dep_par : Cloud optical depth  using efffective radius
!                         calculated from Fountoukis & Nenes, 2005
!     ----------------------------------------------
!
      real*8, allocatable :: cloud_depth      (:,:,:)
      real*8, allocatable :: cloud_opt_dep    (:,:,:)
      real*8, allocatable :: cloud_opt_dep_n  (:,:,:)
      real*8, allocatable :: cloud_opt_dep_par(:,:,:)
!
!     ----------------------------------------------------
!
!      cdnc : Cloud Droplet Number Concentration  (cm-3)
!             for Continental Stratiform & Maritime clouds
!      cdnc :  (cm-3) calculated from Fountoukis & Nenes, 2005
!    ------------------------------------------------------
!
       real*8, allocatable :: cdnc(:,:,:)
       real*8, allocatable :: cdnc_par(:,:,:)
!
!     -----------------------------------------------------
!     alb     : Cloud albedo
!     alb_n   : Claud albedo using calculated effective radius
!     alb_par : Claud albedo using effective radius calculated
!             from Fountoukis & Nenes, 2005
!     -------------------------------------------------------
!
      real*8, allocatable :: alb(:,:,:)
      real*8, allocatable :: alb_n(:,:,:)
      real*8, allocatable :: alb_par(:,:,:)
!
!    __________________________________________________________
!
!     w_cloud      : Cloud updraft velocity (m/s)
!     super_cloud   : Maximum supersaturation
!     tot_so4      : Total SO4 used in NS parameterization
!     __________________________________________________________
!
      real*8, allocatable :: w_cloud(:,:,:)
      real*8, allocatable :: super_cloud(:,:,:)
      real*8,     pointer :: tot_so4 (:,:,:)
!     ----------------------------------------------------------
      real*8, allocatable :: tau_cloud_t(:,:,:)
      real*8, allocatable :: tauro(:,:,:)

      real*8, allocatable :: AIE_1(:,:)
      real*8, allocatable :: AIE_2(:,:)
      real*8, allocatable :: AIE_3(:,:)

!     -----------------------------------------------------
!     QautNS_KK1 : Autoconversion by Khairoutdinov and Kogan, [2000] Eq. 29
!     QautNS_KK2 : Autoconversion by Khairoutdinov and Kogan, [2000] Eq. 30
!     QautNS_R   : Autoconversion by Rostayn, [2000] Eq. 15
!     -------------------------------------------------------
       real*8, allocatable :: QautNS_KK1(:,:,:)
       real*8, allocatable :: QautNS_KK2(:,:,:)
       real*8, allocatable :: QautNS_R(:,:,:)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) then
        Write (6,*) 'calcCloudPropertiesGT called by ', loc_proc
      end if

      Allocate (c0 (k1:k2))
      Allocate (c1 (k1:k2))
      Allocate (c2 (k1:k2))
      c0 = 0.0d0; c1 = 0.0d0; c2 = 0.0d0

      Allocate (cfc(k1:k2))
      Allocate (cfs(k1:k2))
      cfc = 0.0d0; cfs = 0.0d0

      Allocate (rhc(k1:k2))
      rhc = 0.0d0

      Allocate (cmfcol    (k1:k2))
      Allocate (height_asl(k1:k2))
      cmfcol = 0.0d0; height_asl = 0.0d0

      Allocate (clfrac (i1:i2, ju1:j2, k1:k2))
      Allocate (cloud_frac (i1:i2, ju1:j2,k1:k2))
      clfrac = 0.0d0; cloud_frac=0.0d0

      Allocate (lwccol (i1:i2, ju1:j2, k1:k2))
      Allocate (lwpath (i1:i2, ju1:j2, k1:k2))
      Allocate (lwppath(i1:i2, ju1:j2, k1:k2))
      Allocate (relhum (i1:i2, ju1:j2, k1:k2))
      lwccol = 0.0d0; lwpath = 0.0d0;lwppath = 0.0d0;relhum = 0.0d0

      Allocate (eff_rad    (i1:i2, ju1:j2, k1:k2))
      eff_rad = 0.0d0

      Allocate (cloud_opt_dep    (i1:i2, ju1:j2, k1:k2))
      Allocate (cloud_opt_dep_n  (i1:i2, ju1:j2, k1:k2))
      Allocate (cloud_opt_dep_par(i1:i2, ju1:j2, k1:k2))
      cloud_opt_dep = 0.0d0; cloud_opt_dep_n = 0.0d0
      cloud_opt_dep_par = 0.0d0

      Allocate (cloud_depth(i1:i2, ju1:j2, k1:k2))
      cloud_depth = 0.0d0

      Allocate (cdnc(i1:i2, ju1:j2, k1:k2))
      Allocate (cdnc_par(i1:i2, ju1:j2, k1:k2))
      Allocate (eff_rad_n(i1:i2, ju1:j2, k1:k2))
      Allocate (eff_rad_par(i1:i2, ju1:j2, k1:k2))
      cdnc = 0.0d0; cdnc_par = 0.0d0; eff_rad_n=0.0d0
      eff_rad_par=0.0d0

      Allocate (alb(i1:i2, ju1:j2, k1:k2))
      Allocate (alb_n(i1:i2, ju1:j2,k1:k2))
      Allocate (alb_par(i1:i2, ju1:j2, k1:k2))
      alb = 0.0d0; alb_n = 0.0d0; alb_par=0.0d0

      Allocate (super_cloud(i1:i2, ju1:j2, k1:k2))
      Allocate (w_cloud   (i1:i2, ju1:j2, k1:k2))
      Allocate (tot_so4   (i1:i2, ju1:j2, k1:k2))
      super_cloud = 0.0d0; w_cloud = 0.0d0; tot_so4 = 0.0d0
 
      Allocate(tau_cloud_t(i1:i2, ju1:j2, k1:k2))
      Allocate(tauro(i1:i2, ju1:j2, k1:k2))
      tau_cloud_t =0.0d0
      tauro=0.0d0
   
      Allocate(AIE_1(i1:i2, ju1:j2))
      Allocate(AIE_2(i1:i2, ju1:j2))
      Allocate(AIE_3(i1:i2, ju1:j2))
      AIE_1 = 0.0d0; AIE_2 = 0.0d0;
      AIE_3 = 0.0d0

!!!! Begin rsot
      Allocate (QautNS_KK1(i1:i2, ju1:j2, k1:k2))
      Allocate (QautNS_KK2(i1:i2, ju1:j2, k1:k2))
      Allocate (QautNS_R(i1:i2, ju1:j2, k1:k2))
      QautNS_KK1=0.0d0; QautNS_KK2=0.0d0; QautNS_R=0.0d0 
!!!! End rsot


!     if (first) then
 
!       first = .false.
     
!     -----------------------------------------------------------
!     Calculate relative humidity from Seinfeld (1986) p. 181.
!     The first rh is the temperature dependent parameter a.
!     The second rh is the saturation vapor pressure of water.
!     The third rh is the actual relative humidity as a fraction.
!     Then make sure rh is between 0 and 0.95.
!     -----------------------------------------------------------


      do ij = ju1, j2
        do il = i1, i2

          relhum(il,ij,:) =  &
     &      1.0d0 - (373.15d0 / kel(il,ij,:))

          relhum(il,ij,:) =  &
     &      1013.25d0 * Exp (13.3185d0 * relhum(il,ij,:)    -  &
     &                        1.9760d0 * relhum(il,ij,:)**2 -  &
     &                        0.6445d0 * relhum(il,ij,:)**3 -  &
     &                        0.1299d0 * relhum(il,ij,:)**4)

          relhum(il,ij,:) =  &
     &      humidity(il,ij,:) * MWTAIR / 18.0d0 / GPKG *  &
     &      press3c(il,ij,:) / relhum(il,ij,:)

          relhum(il,ij,:) =  &
     &      Max (Min (relhum(il,ij,:), 0.95d0), 0.0d0)

        end do
      end do

!
       cloud_param(:,:,:,1)=relhum(:,:,:)  !uncommented by rsot
!      cloud_param(:,:,:,1) = flux_gt(:,:,:)  ! flux_gt_AIE(:,:,:)
!       -------------------------------------------------------------------
!       Calculate cloud water from a parameterization of Kiehl found in
!       "Sensitivity of the Simulated Climate to a Diagnostic Formulation
!       for Cloud Liquid Water" by James Hack, Journal of Climate, Vol. 11,
!       July 1998, p 1499.
!       -------------------------------------------------------------------

        do ij = ju1, j2
          hl(ij) = 1080.d0 + 2000.d0 * (COS(latdeg(ij) * RADPDEG)) ** 2.d0
!          hl(ij) = 1080.d0 + 2000.d0 * (DCOS(latdeg(ij) * RADPDEG)) ** 2.d0
          
          do il = i1, i2

            height_asl(:) =  &
     &        -29.3d0 * kel(il,ij,:) *  &
     &        Log (press3c(il,ij,:) / 1013.25d0)

            height_asl(:) = Max (height_asl(:), 0.0d0)
            height_m(il,ij,:) = height_asl(:)


            lwccol(il,ij,:) =  0.21d0 * dexp(-height_asl(:)/hl(ij))  ! 0.18d0 *  & this # was revised in a later paper

!!!!!!!     &             (1080.0d0 +  &
!!!!!!!     &              (2000.0d0 * Cos (latdeg(ij) * RADPDEG) ** 2.0d0))) !!!!!commented by rsot on 03/26/07

                
!           ---------------------------------------------------------------
!           Calculate stratiform cloud fractions from a parameterization of
!           Sundqvist et al. (1988).
!           ---------------------------------------------------------------

!           -----------------------------------------------------
!           Now set up the constants to be used in the stratiform
!           parameterization; from Xu and Krueger.
!           -----------------------------------------------------

            where ((height_asl(:) < 15000.0d0) .and.  &
     &             (height_asl(:)> 5900.0d0))
!            c1 (:) =  0.012d0
!            c2 (:) =  0.321d0
              rhc(:) =  0.445d0
!!! rsot 03/26/07              rhc(:) =  0.460d0
            end where

            where ((height_asl(:) <= 5900.0d0) .and.  &
     &             (height_asl(:)> 2350.0d0))
!              c1 (:) = -0.102d0
!              c2 (:) =  7.473d0
              rhc(:) =  0.715d0
            end where

            where  (height_asl(:) <= 2350.0d0)
!              c1 (:) =  0.144d0
!              c2 (:) =  8.050d0
              rhc(:) =  0.780d0
            end where

            where (relhum(il,ij,:)>rhc(:) .and.  &
     &             height_asl(:)   < 15000.0d0)

              cfs(:) =  &
     &          1.0d0 -  &
     &          Sqrt (1.0d0 - (relhum(il,ij,:) - rhc(:)) /  &
     &                         (1.0d0 - rhc(:)))
              cfs(:) = Max (cfs(:), 0.0d0)

            elsewhere

              cfs(:) = 0.0d0

            end where
!
!           ------------------------------------------------------------
!           Calculate convective cloud fractions from a parameterization
!           of Xu and Krueger, (1991) Monthly Weather Review, Vol. 119 p 342.
!           ------------------------------------------------------------

!           --------------------------------------------------
!           Now set up the constants to be used in the cumulus
!           parameterization.
!           --------------------------------------------------

            where ((height_asl(:) < 15000.0d0) .and.  &
     &             (height_asl(:)> 5900.0d0))
	     c0(:) = 0.0131d0
	     c1(:) = 0.0123d0
	     c2(:) = 0.0030d0
!!! rsot 03/26/07              c0(:) = 0.0708d0  ! rsot 03/26/07
!!! rsot 03/26/07              c1(:) = 0.0941d0  ! rsot 03/26/07
!!! rsot 03/26/07              c2(:) = 0.0233d0  ! rsot 03/26/07
            end where

            where ((height_asl(:) <= 5900.0d0) .and.  &
     &             (height_asl(:)> 2350.0d0))
              c0(:) = 0.0722d0
              c1(:) = 0.0760d0
              c2(:) = 0.0254d0
            end where

            where  (height_asl(:) <= 2350.0d0)
              c0(:) = 0.0337d0
              c1(:) = 0.0453d0
              c2(:) = 0.0237d0
            end where


!           -------------------------------------
!           Calculate the cumulus cloud fraction.
!           -------------------------------------

            cmfcol(:) = cmf1(il,ij,:) * 3600.0d0 * 9.8d0 / 100.0d0

            where (cmfcol(:)>0.01d0)

              cfc(:) = c0(:) +  &
     &                 c1(:) *  Log10 (cmfcol(:)) +  &
     &                 c2(:) * (Log10 (cmfcol(:)))**2
              cfc(:) = Max (cfc(:), 0.0d0)

            elsewhere

              cfc(:) = 0.0d0

            end where


!           ---------------------------------------------------------
!           Compute total cloud fraction for combined large-scale and
!           convective clouds (GEOS-1 GCM descriptive documentation).
!           ---------------------------------------------------------

            clfrac(il,ij,:) =  &
     &        1.0d0 - (1.0d0 - cfc(:)) * (1.0d0 - cfs(:))

          end do
        end do


!      -------------------------
!      Added by Xhliu, 02/05/03.
!      -------------------------

        where ((clfrac(:,:,:) <= 1.0d-2) .or.  &
     &         (lwccol(:,:,:) <= 1.0d-9))
          clfrac(:,:,:) = 0.0d0
          lwccol(:,:,:) = 0.0d0
        end where

!     ======
!      end if
!     ======
!
!
       cloud_param(:,:,:,2) = clfrac(:,:,:)

!!!!added by rsot to calculate total cloud fraction

        do il = i1, i2
          do ij = ju1, j2
            do ik = k1, k2
		afrac(ik)= clfrac(il, ij, ik)
            end do
	        res(il, ij)=maxval(afrac)
          end do
        end do
       cloud_param(:,:,1,2) = res(:,:)

       afrac(:) = 0.d0
       res(:,:) = 0.d0

!!!!end of calculations
!             
!
!       --------------------------------------------------------------------
!       Compute total cloud fraction using max_cloud, ran_cloud
!       --------------------------------------------------------------------

        cloud_frac(i1:i2, ju1:j2, k1:k2) = totalCloudFraction(i1:i2, ju1:j2, k1:k2)

!        do ik = k1, k2
!          do ij = ju1, j2
!            do il = i1, i2
!                   cloud_frac(il, ij, ik) =  &
!     &                  1.0d0 -  &
!     &                 (1.0d0 - max_cloud(il, ij, ik)) *  &
!     &                 (1.0d0 - ran_cloud(il, ij, ik))
!            end do
!          end do
!        end do

!!!!added by rsot to calculate total cloud fraction in an alternate way

!	 do il = i1, i2
!	   do ij = ju1, j2
!	     do ik = k1, k2
!		 afrac(ik)= cloud_frac(il, ij, ik)
!	     end do
!		 res(il, ij)=maxval(afrac)
!	   end do
!        end do
!       cloud_param(:,:,2,2) = res(:,:)

!!!!end of calculations


!         cloud_param(:,:,:,1)=cloud_frac(:,:,:)
!
!                 ===============================
!!!!!!!!!         clfrac(:,:,:)=cloud_frac(:,:,:)      !!!!!!!!!!!!>>>>>>>>>>>>>>>>>>>>>>>>
!                 ===============================
!          cloud_param(:,:,:,15)=clfrac(:,:,:) ! commented by rsot 071006
!rsot
          FCLD(:,:,:) = clfrac(:,:,:)
!rsot
!
!       -----------------------------------------
!       In-cloud liquid water mixing ratio (g/g).
!       -----------------------------------------

        do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2
   
              if (moistq(il,ij,ik) < -1.0d-5) then

                if (clfrac(il,ij,ik)>1.0d-2) then

                  lwccol(il,ij,ik) =  &
     &              lwccol(il,ij,ik) * 1.0d-6 /  &
     &              (concentration(imgas_num)%pArray3D(il,ij,ik) * MWTAIR / AVOGAD)

!c                -----------------------------------------------------------
!c                Do not allow cloud in the lowest level according to NCAR
!c                CCM3 because LWC is prescribed as the maximum in the lowest
!c                level, Xhliu, 01/17/03.
!c                -----------------------------------------------------------

                  if (ik == k1) lwccol(il,ij,ik) = 0.0d0

!c                ------------------------------------------------------------
!c                Do not allow cloud liquid water when temperature below -25C.
!c                ------------------------------------------------------------

                  if (kel(il,ij,ik) < 248.0d0) then
                    lwccol(il,ij,ik) = 0.0d0
                  end if

                else

                  lwccol(il,ij,ik) = 0.0d0

                end if

              else

                lwccol(il,ij,ik) = 0.0d0

              end if

            end do
          end do
        end do
!
		QLWC(:,:,:)=lwccol(:,:,:)   !rsot

!
!      ================================================================
!      Alternate way of calculating cloud fraction
!      ----------------------------------------------------------------
!        The cloud overlap scheme used is a maximum-random overlap
!         scheme, i.e., continuous cloud layers are assumed to be maximally
!        overlapped, while discontinuous cloud layers are randomly
!        overlapped (Geleyn and Hollingsworth, 1979; Feng
!        et al., 2004).
!        ----------------------------------------------------------------

!      cloud_frac(:,:,:) =  &
!    &  1.0d0 -  &
!    &  ((1.0d0 - Maxval (max_cloud(:,:,:))) *  &
!    &   (Product ((1.0d0 - ran_cloud(:,:,:)))))

!       ----------------------------------------------------------
!       Compute cloud drop effective radius from Kiehl, JGR, 1998.
!       ----------------------------------------------------------

!       -----------------------------------------
!       Ocean, sea ice, and land (kel < 243.15d0).
!
!       -----------------------------------------

        eff_rad(:,:,:) = 1.0d-5

        do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2

              if (lwi_flags(il,ij) == 2) then  ! land

!rsot                if ((kel(il,ij,ik) >= 243.15d0) .and.  &
!rsot     &              (kel(il,ij,ik) <= 263.15d0)) then

!rsot                  eff_rad(il,ij,ik) =  &
!rsot     &              5.0d-6 -  &
!rsot     &              5.0d-6 *  &
!rsot     &              ((kel(il,ij,ik) + ABS_ZERO + 10.d0) / 20.d0)

!rsot                else if (kel(il,ij,ik) >= 263.15d0) then
                if (kel(il,ij,ik) >= 263.15d0) then

                  eff_rad(il,ij,ik) = 5.0d-6

                end if

              end if

            end do
          end do
        end do

!
!!!!!!!!!!!!!!!!!!!!!!!!        cloud_param(:,:,:,4) = eff_rad(:,:,:)

!       ----------------------------------------------------
!       Estimate depth of a cloud.
!       ----------------------------------------------------

       do il = i1, i2
         do ij = ju1, j2
           do ik = k2, k1, -1
!
             if (moistq(il,ij,ik) < -1.0d-5) then
               cloud_depth(il,ij,ik) = gridBoxHeight(il,ij,ik)
             end if
!
           end do
         end do
       end do

!  
!      ---------------------------------------------------------
!       In-cloud liquid water mixing ratio- lwccol is in (g/g)
!       Convert to g/m^3 to be used in cloud optical depth
!       calculations
!      ---------------------------------------------------------
!
        do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2
                lwccol(il,ij,ik)=  &
     &              lwccol(il,ij,ik) / 1.0d-6 *  &
     &              (concentration(imgas_num)%pArray3D(il,ij,ik) * MWTAIR / AVOGAD)
             end do
          end do
        end do
!
!       -----------------------------------------
!       Grid-averaged cloud liquid water (g/m^3).
!       -----------------------------------------

!       do ik = k1, k2
!         do ij = ju1, j2
!           do il = i1, i2

!             lwccol(il,ij,ik) =  &
!    &          lwccol(il,ij,ik) *clfrac(il,ij,ik)

!           end do
!         end do
!       end do


!         cloud_param(:,:,:,10) =lwccol(:,:,:)   !! commented by rsot
          cloud_param(:,:,:,3)  =lwccol(:,:,:)   !! added by rsot

!!!!Begin: additions by rsot

!        ---------------------------------------------
!        Grid-averaged cloud liquid water path (g/m^2)     !rsot
!        ---------------------------------------------

       do ik = k1, k2
         do ij = ju1, j2
           do il = i1, i2

             if (ik == 1) then
		lwpath(il,ij,ik)=0.d0
	     else
		lwpath(il,ij,ik)=-0.21d0*hl(ij)*clfrac(il,ij,ik) &
     &               *(dexp(-height_m(il,ij,ik)/hl(ij))-dexp(-height_m(il,ij,ik-1)/hl(ij)))
	     endif

           end do
         end do
       end do

!        ------------------------------------------------------------------------
!        Alternate way to calculate grid-averaged cloud liquid water path (g/m^2)     !rsot
!        ------------------------------------------------------------------------
       do ik = k1, k2
         do ij = ju1, j2
           do il = i1, i2
	     if (ik == 1) then
		lwppath(il,ij,ik)=0.d0
	     else
                lwppath(il,ij,ik)=lwccol(il,ij,ik)*cloud_depth(il,ij,ik)*clfrac(il,ij,ik)
	     endif	
           end do
         end do
       end do

        do il = i1, i2
          do ij = ju1, j2
            do ik = k1, k2
		LWPfrac(ik)= lwppath(il, ij, ik)
            end do
	        resLWP(il, ij)=sum(LWPfrac)
          end do
        end do
       
	lwpath(:,:,1)= resLWP(:,:)


!!!!End: additions by rsot
!
!  ========================================================================
!  For GISS and FVGCM use cloud optical depth provided    
!  ========================================================================

          tau_cloud_t(:,:,:)= tau_cloud(:,:,:)

!Beg Jules
!!  For DAO
!!    ------------------------------------------------------------------------
!           if (metdata_name_org_c(1:3) == 'DAO') then
!!         ---------------------------------------------------------------
!!         The met data did not have optical depth, so calculate it from
!!         the cloud fractions, temperature and delta pressure.  This
!!         algorithm comes from the GEOS-CHEM model and was given to us by
!!         Bob Yantosca.
!!         ---------------------------------------------------------------
!              where (kel(i1:i2,ju1:j2,:)  <= 190.66d0)  &
!     &           tauro(:,:,:) = 0.0d0
!
!              where ((kel(i1:i2,ju1:j2,:)> 190.66d0) .and.  &
!     &             (kel(i1:i2,ju1:j2,:) <= 263.16d0))  &
!     &           tauro(:,:,:) =  &
!     &             (kel(i1:i2,ju1:j2,:) - 190.66d0)**2 * 2.0d-6
!
!              where ((kel(i1:i2,ju1:j2,:)> 263.16d0) .and.  &
!     &              (kel(i1:i2,ju1:j2,:) <= 273.38d0))  &
!     &           tauro(:,:,:) =  &
!     &              (kel(i1:i2,ju1:j2,:) * 6.95d-3) - 1.82d0
!
!              where (kel(i1:i2,ju1:j2,:) > 273.38d0)  &
!     &           tauro(:,:,:) = 0.08d0
!
!             tau_cloud_t(:,:,:) =  &
!     &          (0.16d0       * max_cloud(:,:,:)) +  &
!     &          (tauro(:,:,:) * ran_cloud(:,:,:))
!
!              do ij = ju1, j2
!                do il = i1, i2
! 
!                   tau_cloud_t(il,ij,k1:k2) =  &
!     &              tau_cloud_t(il,ij,k1:k2) *  &
!     &              (  press3e(il,ij,k1-1:k2-1) -  press3e(il,ij,k1:k2) )
!
!!     &              ((ai(k1-1:k2-1)*pt + bi(k1-1:k2-1) * pctm2(il,ij))-  &
!!     &               (ai(k1:k2)    *pt + bi(k1:k2)     * pctm2(il,ij)))
!
!                end do
!              end do
!
!                cloud_param(:,:,:,4) = tau_cloud_t(:,:,:)
!           end if
!End Jules

                cloud_param(:,:,:,4) = tau_cloud_t(:,:,:)
!!
!       ---------------------------------------------------------------
!       Estimate cloud optical depth (Seinfeld & Pandis, p.1173, 1998).
!       ---------------------------------------------------------------

        do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2

               if ((lwi_flags(il,ij) == 2) .and. (kel(il,ij,ik) >= 263.15d0)) then     !Use only for liquid clouds

                   cloud_opt_dep(il,ij,ik) =  (clfrac(il,ij,ik)**1.5) *  &
     &                     ((3.0d0 * cloud_depth(il,ij,ik) *  &
     &                       lwccol(il,ij,ik) * EXT_EFF) /  &
     &                      (4.0d0 * LWD * eff_rad(il,ij,ik)))
               
               else if ((lwi_flags(il,ij) == 1) .and. (kel(il,ij,ik) >= 269.15d0)) then     !Use only for liquid clouds

                   cloud_opt_dep(il,ij,ik) =  (clfrac(il,ij,ik)**1.5) *  &
     &                     ((3.0d0 * cloud_depth(il,ij,ik) *  &
     &                       lwccol(il,ij,ik) * EXT_EFF) /  &
     &                      (4.0d0 * LWD * eff_rad(il,ij,ik)))
               else
                   cloud_opt_dep(il,ij,ik)=tau_cloud_t(il,ij,ik)  ! Use what is provided from meteorology
               end if

            end do
          end do
        end do


              cloud_param(:,:,:,5) = cloud_opt_dep(:,:,:)
!
!     ----------------------------------------------------------------------------------
!     Use cloud optical depth to estimate cloud albedo (Seinfeld & Pandis, p.1174, 1998).
!
!      -----------------------------------------------------------------------------------
 
        do ik = k1, k2      
           do ij = ju1, j2
             do il = i1, i2
                 alb(il,ij,ik)=cloud_opt_dep(il,ij,ik)/(cloud_opt_dep(il,ij,ik)+7.7)
            end do
          end do
       end do 
        alb(:,:,:) =  &
     &      Max (Min (alb(:,:,:), 0.95d0), 0.0d0)

!         cloud_param(:,:,:,10) = alb(:,:,:)    ! droebit mieci lwccol
          cloud_param(:,:,:,10) = alb(:,:,:)    ! uncommented by rsot 071006
!
!!!!Start rsot modifications; Most part originally by nika
!
!     =======================================================================
!
!                    Calculate Cloud Microphysical Properties
!
!     =======================================================================
!
!  ----------------------------------------------------------------------------------
!                  Boucher and Lohmann, 1995 empirical formulation
!     sulafate is in ug/m3
!  ----------------------------------------------------------------------------------
       do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2
!
              if (cloud_param(il,ij,ik,7).gt.0.d0) then
! For land:
              if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>= 263.15d0).and.(lwi_flags(il,ij) == 2)) then
                      aero_typ = 2                            ! Polluted continental
     	  	      call NdBL (cloud_param(il,ij,ik,7),aero_typ, cdnc(il,ij,ik))
              cdnc(il,ij,ik) = MIN(MAX(cdnc(il,ij,ik),MIN_D),MAX_D)
! For ocean:
              else if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>= 269.15d0).and.(lwi_flags(il,ij) == 1)) then
                      aero_typ = 1                          ! Polluted Marine
     	  	      call NdBL (cloud_param(il,ij,ik,7),aero_typ, cdnc(il,ij,ik))
              cdnc(il,ij,ik) = MIN(MAX(cdnc(il,ij,ik),MIN_D),MAX_D)
	      else
	      cdnc(il,ij,ik) = 0.d0
              end if
 
	      else
	      cdnc(il,ij,ik) = 0.d0
              end if

            end do
          end do
        end do

        IF (cloudDroplet_opt == 1) THEN
             cdnc_par(:,:,:) = cdnc(:,:,:)
	ENDIF
!
!  ----------------------------------------------------------------------------------
!         Nenes and Seinfeld 2003 & Fountoukis and Nenes, 2005 parameterizaton
!  ----------------------------------------------------------------------------------
      IF (cloudDroplet_opt == 2) THEN
          do ik = k1, k2
            do ij = ju1, j2
              do il = i1, i2
!
              if (cloud_param(il,ij,ik,7).gt.0.d0) then
! For land:
              if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>= 263.15d0).and.(lwi_flags(il,ij) == 2)) then
                  aero_typ = 2                            ! Polluted continental
                  call paramdriver  &
     &                    (cloudGT, kel(il,ij,ik),press3c(il,ij,ik),w_cloud(il,ij,ik),  &
     &                     cloud_param(il,ij,ik,7), cdnc_par(il,ij,ik),  &
     &                     super_cloud(il,ij,ik), aero_typ)
              cdnc_par(il,ij,ik) = MIN(MAX(cdnc_par(il,ij,ik)/MLN,MIN_D),MAX_D)   !convert from  #/m3 to #/cm3

! For ocean:
              else if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>= 269.15d0).and.(lwi_flags(il,ij) == 1)) then
                  aero_typ = 1                          ! Polluted Marine
                  call paramdriver  &
     &                    (cloudGT, kel(il,ij,ik),press3c(il,ij,ik),w_cloud(il,ij,ik),  &
     &                     cloud_param(il,ij,ik,7), cdnc_par(il,ij,ik),  &
     &                     super_cloud(il,ij,ik),aero_typ)
              cdnc_par(il,ij,ik) = MIN(MAX(cdnc_par(il,ij,ik)/MLN,MIN_D),MAX_D)   !convert from  #/m3 to #/cm3
	      else
	      cdnc_par(il,ij,ik) = 0.d0
              endif

	      else
	      cdnc_par(il,ij,ik) = 0.d0
              endif

            end do
         end do
        end do
!                                       
      ENDIF
!
!! Start rsot
!  ----------------------------------------------------------------------------------
!                     Abdul-Razzak and Ghan, 1999 parameterization
!  ----------------------------------------------------------------------------------
      IF (cloudDroplet_opt == 3) THEN
          do ik = k1, k2
            do ij = ju1, j2
              do il = i1, i2
		
              if (cloud_param(il,ij,ik,7).gt.0.d0) then
! For land:
              if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>= 263.15d0).and.(lwi_flags(il,ij) == 2)) then
                  aero_typ = 2                            ! Polluted continental
                  call paramAG  &
     &                    (cloudGT, kel(il,ij,ik),press3c(il,ij,ik),w_cloud(il,ij,ik),  &
     &                     cloud_param(il,ij,ik,7), cdnc_par(il,ij,ik),  &
     &                     super_cloud(il,ij,ik), aero_typ)
              cdnc_par(il,ij,ik) = MIN(MAX(cdnc_par(il,ij,ik)/MLN,MIN_D),MAX_D)   !convert from  #/m3 to #/cm3

! For ocean:
              else if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>= 269.15d0).and.(lwi_flags(il,ij) == 1)) then
                  aero_typ = 1                          ! Polluted Marine
                  call paramAG  &
     &                    (cloudGT, kel(il,ij,ik),press3c(il,ij,ik),w_cloud(il,ij,ik),  &
     &                     cloud_param(il,ij,ik,7), cdnc_par(il,ij,ik),  &
     &                     super_cloud(il,ij,ik),aero_typ)
              cdnc_par(il,ij,ik) = MIN(MAX(cdnc_par(il,ij,ik)/MLN,MIN_D),MAX_D)   !convert from  #/m3 to #/cm3
! For the rest:
	      else
	      cdnc_par(il,ij,ik) = 0.d0
              end if

	      else
	      cdnc_par(il,ij,ik) = 0.d0
              end if

            end do
         end do
        end do
!
      ENDIF
!
!  ----------------------------------------------------------------------------------
!                     Segal and Khain, 2006 parameterization
!  ----------------------------------------------------------------------------------
      IF (cloudDroplet_opt == 4) THEN
          do ik = k1, k2
            do ij = ju1, j2
              do il = i1, i2
!
		
              if (cloud_param(il,ij,ik,7).gt.0.d0) then
! For land:
              if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>= 263.15d0).and.(lwi_flags(il,ij) == 2)) then
                  aero_typ = 2                            ! Polluted continental
                  call paramSK  &
     &                    (w_cloud(il,ij,ik), cloud_param(il,ij,ik,7), &
     &                     cdnc_par(il,ij,ik), aero_typ)
              cdnc_par(il,ij,ik) = MIN(MAX(cdnc_par(il,ij,ik),MIN_D),MAX_D)  

! For ocean:
              else if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>= 269.15d0).and.(lwi_flags(il,ij) == 1)) then
                  aero_typ = 1                          ! Polluted Marine
                  call paramSK  &
     &                    (w_cloud(il,ij,ik), cloud_param(il,ij,ik,7), &
     &                     cdnc_par(il,ij,ik), aero_typ)
              cdnc_par(il,ij,ik) = MIN(MAX(cdnc_par(il,ij,ik),MIN_D),MAX_D)   
! For the rest:
	      else
	      cdnc_par(il,ij,ik) = 0.d0
	
              end if

	      else
	      cdnc_par(il,ij,ik) = 0.d0
              end if

            end do
         end do
        end do
!
      ENDIF
!
!! End rsot
!     ----------------------------------------------------------------------
!       Calculate Droplet Effective Redius  eff_rad_n  (m)
!                                            using Boucher and Lohmann, 1995
!
!     ----------------------------------------------------------------------
         do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2

              if  (cdnc(il,ij,ik).gt.0.d0) then

! For land:
              if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>= 263.15d0).and.(lwi_flags(il,ij) == 2)) then
                                 eff_rad_n(il,ij,ik)=(((lwccol(il,ij,ik)*THRFTH*ONEOM)/  &
     &                                     (PI*LWD *BETAcon*cdnc(il,ij,ik)))**ONETH)
              eff_rad_n(il,ij,ik) = MIN(MAX(eff_rad_n(il,ij,ik),5.0d-6),20.d-6)
! For ocean:
              else if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>= 269.15d0).and.(lwi_flags(il,ij) == 1)) then
                                 eff_rad_n(il,ij,ik)=(((lwccol(il,ij,ik)*THRFTH*ONEOM)/  &
     &                                     (PI*LWD *BETAmar*cdnc(il,ij,ik)))**ONETH)
              eff_rad_n(il,ij,ik) = MIN(MAX(eff_rad_n(il,ij,ik),5.0d-6),20.d-6)
! For the rest:
              else
              eff_rad_n(il,ij,ik) = 0.0d0
              end if

              else
              eff_rad_n(il,ij,ik) = 0.0d0
!
!rsot                       end if
               end if
            end do
          end do
        end do
!

!
!     ----------------------------------------------------------------------
!         Calculate Droplet Effective Redius using other parameteizations
!                                 eff_rad_par  (m)
!     ----------------------------------------------------------------------
!        eff_rad_par(:,:,:) = 1.0d-5  ! set for  the ice clouds

         do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2

              if  (cdnc_par(il,ij,ik).gt.0.d0) then
!
! For land:
              if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>= 263.15d0).and.(lwi_flags(il,ij) == 2)) then
                                    eff_rad_par(il,ij,ik)=(((lwccol(il,ij,ik)*THRFTH*ONEOM)/    &
     &                                           (PI*LWD *BETAcon*cdnc_par(il,ij,ik)))**ONETH)
              eff_rad_par(il,ij,ik) = MIN(MAX(eff_rad_par(il,ij,ik),5.0d-6),20.d-6)
! For ocean:
              else if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>= 269.15d0).and.(lwi_flags(il,ij) == 1)) then
                                     eff_rad_par(il,ij,ik)=(((lwccol(il,ij,ik)*THRFTH*ONEOM)/  &
     &                                            (PI*LWD *BETAmar*cdnc_par(il,ij,ik)))**ONETH)
              eff_rad_par(il,ij,ik) = MIN(MAX(eff_rad_par(il,ij,ik),5.0d-6),20.d-6)
! For the rest:
              else
              eff_rad_par(il,ij,ik) = 0.0d0
              end if

              else
              eff_rad_par(il,ij,ik) = 0.0d0
!
!rsot                     end if
               end if
            end do
          end do
        end do
!
!    ------------------------------------------------------------------------
!    Calculate cloud optical depth using new Droplet Effective Redius from BL
!                               cloud_opt_dep_n        
!    ------------------------------------------------------------------------

      do ik = k1, k2
         do ij = ju1, j2
            do il = i1, i2

              if  (cdnc(il,ij,ik).gt.0.d0)  then
! For land:
                 if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>=263.15d0).and.(lwi_flags(il,ij)==2)) then    ! Continental cloud
                   cloud_opt_dep_n(il,ij,ik) = (clfrac(il,ij,ik)**1.5)* &
     &                       ((3.0d0 * cloud_depth(il,ij,ik) *  &
     &                         lwccol(il,ij,ik) * EXT_EFF) /  &
     &                         (4.0d0 * LWD * eff_rad_n(il,ij,ik)))
! For ocean:
                 else if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>=269.15d0).and.(lwi_flags(il,ij)==1)) then    ! Marine cloud
                   cloud_opt_dep_n(il,ij,ik) = (clfrac(il,ij,ik)**1.5)* &
     &                       ((3.0d0 * cloud_depth(il,ij,ik) *  &
     &                         lwccol(il,ij,ik) * EXT_EFF) /  &
     &                         (4.0d0 * LWD * eff_rad_n(il,ij,ik)))
                 else
!!!                   cloud_opt_dep_n(il,ij,ik)= tau_cloud_t(il,ij,ik)
                    cloud_opt_dep_n(il,ij,ik)= 0.d0
                 end if

              else
                 cloud_opt_dep_n(il,ij,ik)= 0.d0
              end if

           end do
        end do
      end do
!
!    ----------------------------------------------------------------------
!    Calculate cloud optical depth using new Droplet Effective Redius from
!    other parameterizations
!                             cloud_opt_dep_par
!    ----------------------------------------------------------------------

      do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2

              if  (cdnc_par(il,ij,ik).gt.0.d0) then
! For land:
              if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>=263.15d0).and.(lwi_flags(il,ij)==2)) then  ! continental cloud
                    cloud_opt_dep_par(il,ij,ik) = (clfrac(il,ij,ik)**1.5)* &
     &                   ((3.0d0 * cloud_depth(il,ij,ik) *  &
     &                   lwccol(il,ij,ik) * EXT_EFF) /  &
     &                   (4.0d0 * LWD * eff_rad_par(il,ij,ik)))
! For ocean:
                 else if ((lwccol(il,ij,ik)>0.0d0).and.(kel(il,ij,ik)>=269.15d0).and.(lwi_flags(il,ij)==1)) then    ! Marine cloud
                    cloud_opt_dep_par(il,ij,ik) = (clfrac(il,ij,ik)**1.5)* &
     &                   ((3.0d0 * cloud_depth(il,ij,ik) *  &
     &                   lwccol(il,ij,ik) * EXT_EFF) /  &
     &                   (4.0d0 * LWD * eff_rad_par(il,ij,ik)))
                 else
                    cloud_opt_dep_par(il,ij,ik)= 0.d0
!!!                   cloud_opt_dep_par(il,ij,ik)= tau_cloud_t(il,ij,ik)
                 end if

              else
                 cloud_opt_dep_par(il,ij,ik)= 0.d0
              end if

            end do
          end do
        end do
!
!    -------------------------------------------------------------------------
!     Use cloud optical depth determined from BL to estimate new cloud albedo
!                       (Seinfeld & Pandis, p.1174, 1998)
!    -------------------------------------------------------------------------
       do ik = k1, k2
          do ij = ju1, j2
             do il = i1, i2

             alb_n(il,ij,ik)=cloud_opt_dep_n(il,ij,ik)/(cloud_opt_dep_n(il,ij,ik)+7.7)
             alb_n(il,ij,ik) = Max (Min (alb_n(il,ij,ik), 0.95d0), 0.0d0)

             end do
          end do
       end do

!
!    -----------------------------------------------------------------------------------
!     Use cloud optical depth from other parameterizations to estimate new cloud albedo
!                           (Seinfeld & Pandis, p.1174, 1998).
!    -----------------------------------------------------------------------------------
       do ik = k1, k2
          do ij = ju1, j2
             do il = i1, i2

             alb_par(il,ij,ik)=cloud_opt_dep_par(il,ij,ik)/(cloud_opt_dep_par(il,ij,ik)+7.7)
             alb_par(il,ij,ik) = Max (Min (alb_par(il,ij,ik), 0.95d0), 0.0d0)

             end do
          end do
       end do
!
!     -----------------------------------------------------------------------------------
!     calculate AIE by calculating  AIE= 1/3*F(TOA)*T^2*Ac*Rc(1-Rc) . Monthly averages of 
!     this variable will be multiplied on a (N(CD)-N(PI))/N(CD) to get Indirect Forcing
!     -----------------------------------------------------------------------------------
 
!     >>>>>>>>>>>>>   Using GMI outputs <<<<<<<<<<<<<<<<<<<     
 
!      do ij = ju1, j2
!        do il = i1, i2
!          AIE_1(il,ij)=ONETH *(flux_gt_AIE(il,ij,2)*clfrac(il,ij,3)*   &
!    &                alb(il,ij)*(1.0d0- alb(il,ij)))

!        end do
!      end do        
!     cloud_param(:,:,1,14) = AIE_1(:,:)

!     >>>>>>>>>>>>>  Using BL outputs   <<<<<<<<<<<<<<<<
!       do ij = ju1, j2
!         do il = i1, i2
!           AIE_2(il,ij)=ONETH*(flux_gt_AIE(il,ij,2)*clfrac(il,ij,3)*   &
!    &                alb_n(il,ij)*(1.0d0- alb_n(il,ij)))
!        end do
!      end do

!       cloud_param(:,:,1,9)= AIE_2(:,:)

!    >>>>>>>>>>>>  Using NS outputs  <<<<<<<<<<<<<<<<<<<
    
!      do ij = ju1, j2
!        do il = i1, i2
!           AIE_3(il,ij)=ONETH*(flux_gt_AIE(il,ij,2)*clfrac(il,ij,3)*   &
!    &                alb_par(il,ij)*(1.0d0- alb_par(il,ij)))
!        end do
!      end do

!       cloud_param(:,:,1,4)= AIE_3(:,:)


!----------------------------------------------------------------------------------
!   Assign cloud optical depth that will be later used in the photolysis section
!----------------------------------------------------------------------------------


!                             use   cloud_opt_dep(:,:,:)      for GMI default
!                             use   cloud_opt_dep_n(:,:,:)    for BL
!                             use   cloud_opt_dep_par(:,:,:)  for all others including BL when cloudDroplet_opt=1 

!             ===========================================
              cloud_tau(:,:,:)= cloud_opt_dep_par(:,:,:)    			!!! rsot
!             ===========================================

!
!!!!Start rsot
!  --------------------------------------------------------------------------------
!                        Autoconversion Rate Calculations
!  --------------------------------------------------------------------------------

       do ik = k1, k2
          do ij = ju1, j2
             do il = i1, i2

!  Estimate Air density (kg m^-3)
    	RHOAIR(il,ij,ik)=1.d3*concentration(imgas_num)%pArray3D(il,ij,ik)*MWTAIR/AVOGAD

        if  (cdnc_par(il,ij,ik).gt.0.d0) then

!  --------------------------------------------------------------------------------
!               Estimation of autoconversion rate from Droplet Number
!                    [Kharoutdinov and Kogan, 2000 (Eq. 29)]
!  --------------------------------------------------------------------------------
	QLWC(il,ij,ik)=MIN(QLWC(il,ij,ik), 3.d-03)  		!(upper limit for the QLWC)
	QautNS_KK1(il,ij,ik) = FCLD(il,ij,ik)*1350.d0*(QLWC(il,ij,ik)**2.47d0)* &
     &                (MAX(CDNC_par(il,ij,ik), 40d0))**(-1.79d0)

!  --------------------------------------------------------------------------------
!            Estimation of autoconversion rate from Mean Volume Radius, Rv 
!                    [Kharoutdinov and Kogan, 2000 (Eq. 30)]
!  --------------------------------------------------------------------------------
	Rv(il,ij,ik)=1.d6*eff_rad_par(il,ij,ik)/1.1d0  !Rv in um
        Rv(il,ij,ik)=MIN(MAX(Rv(il,ij,ik),4.5d0),18d0)
	IF (Rv(il,ij,ik) .GT. 7.d0) THEN		! Set critical radius for onset of Autoconversion
	  QautNS_KK2(il,ij,ik)  = FCLD(il,ij,ik)*4.1d-15*(Rv(il,ij,ik))**5.67d0
	ELSE
	  QautNS_KK2(il,ij,ik)  = 0.d0
	ENDIF

!  --------------------------------------------------------------------------------
!             Estimation of autoconversion rate from Number of Droplets 
!                           [Rotstayn, 1999 (Eq. 15)]
!  --------------------------------------------------------------------------------
	QautNS_R(il,ij,ik) = FCLD(il,ij,ik)*0.104d0*GRAV*0.55d0*RHOAIR(il,ij,ik)**(4d0/3d0) &
     &         /(1.8d-5*(MAX(CDNC_par(il,ij,ik)*1d6, 40d6)*1d3)**(1d0/3d0))*QLWC(il,ij,ik)**(7d0/3d0)

!	Heaviside function
	QCRIT(il,ij,ik)= 4d0/3d0*PI*1d3*(7.5d-6)**3*MAX(CDNC_par(il,ij,ik)*1d6,40d6)/RHOAIR(il,ij,ik)
	IF (QCRIT(il,ij,ik) .LE. 1d-25) QCRIT(il,ij,ik)= 1d-25
	IF ((QLWC(il,ij,ik) - QCRIT(il,ij,ik)) .LT. 0.d0) THEN
	   QautNS_R(il,ij,ik)  = 0d0
	ELSE IF ((QLWC(il,ij,ik) - QCRIT(il,ij,ik)) .EQ. 0.d0) THEN
	   QautNS_R(il,ij,ik) = 0.5d0 * QautNS_R(il,ij,ik)
	ELSE IF ((QLWC(il,ij,ik) - QCRIT(il,ij,ik)) .GT. 0.d0) THEN
	   QautNS_R(il,ij,ik)  = QautNS_R(il,ij,ik)
	ENDIF

!  --------------------------------------------------------------------------------
!
! For the rest of the cases:
         else if  (cdnc_par(il,ij,ik) .le. 0.d0) then

	 QautNS_KK1(il,ij,ik) = 0.d0
	 QautNS_KK2(il,ij,ik) = 0.d0
	 QautNS_R(il,ij,ik)   = 0.d0

!  --------------------------------------------------------------------------------
       endif
             end do
          end do
       end do

!  ----------------------------------------------------------------------------------
!                                  Save Diagnostics
!  ----------------------------------------------------------------------------------
       cloud_param(:,:,:,8)  = cdnc(:,:,:)
       cloud_param(:,:,:,9)  = eff_rad_n(:,:,:)
       cloud_param(:,:,:,11) = cloud_opt_dep_n(:,:,:)
       cloud_param(:,:,:,12) = lwpath(:,:,:)
       cloud_param(:,:,:,13) = cdnc_par  (:,:,:)
       cloud_param(:,:,:,14) = w_cloud   (:,:,:)
       cloud_param(:,:,1,14) = radswg(:,:)         !!!!!!!!w_cloud   (:,:,:)
       cloud_param(:,:,:,15) = super_cloud(:,:,:)   !uncomment rsot 071006
       cloud_param(:,:,:,16) = eff_rad_par(:,:,:)
       cloud_param(:,:,:,17) = cloud_opt_dep_par(:,:,:)
       cloud_param(:,:,:,18) = alb_par(:,:,:)
       cloud_param(:,:,:,19) = QautNS_KK1  (:,:,:)
       cloud_param(:,:,:,20) = QautNS_KK2  (:,:,:)
       cloud_param(:,:,:,21) = QautNS_R    (:,:,:)

      Deallocate (AIE_1)
      Deallocate (AIE_2)
      Deallocate (AIE_3)
      Deallocate (c0)
      Deallocate (c1)
      Deallocate (c2)
      Deallocate (cfc)
      Deallocate (cfs)
      Deallocate (clfrac)
      Deallocate (cmfcol)
      Deallocate (height_asl)
      Deallocate (lwccol)
      Deallocate (lwpath)
      Deallocate (relhum)
      Deallocate (rhc)
      Deallocate (eff_rad)
      Deallocate (eff_rad_n)
      Deallocate (cloud_depth)
      Deallocate (cloud_opt_dep)
      Deallocate (cloud_opt_dep_n)
      Deallocate (cdnc)
      Deallocate (alb)
      Deallocate (alb_n)
      Deallocate (cdnc_par)
      Deallocate (w_cloud)
      Deallocate (super_cloud)
      Deallocate (eff_rad_par)
      Deallocate (cloud_opt_dep_par)
      Deallocate (alb_par)
      Deallocate (tot_so4)
      Deallocate (tau_cloud_t)
      Deallocate (tauro)
      Deallocate(cloud_frac)
      Deallocate (QautNS_KK1)
      Deallocate (QautNS_KK2)
      Deallocate (QautNS_R)

      return

      end subroutine calcCloudPropertiesGT
!EOC
!------------------------------------------------------------------------------
!
!*************************************************************************
!
! DROPLET ACTIVATION PARAMETERIZATION FILES START HERE
!
! ************************************************************************
!
!
!=======================================================================
!======================================================================
!
!       Code Developer
!       Nicholas Meskhidze, GA TECH
!       nmeskhidze@eas.gatech.edu
!    -------------------------------
!     DESCRIPTION
!
! *** PROGRAM PARAMExample
! *** THIS SUBROUTINE PROVIDES AN EXAMPLE FOR USING THE PARAMETERIZATION
!
!     The calling sequence is as follows:
!     1. CALL GAULEG   (once, only if integration over a PDF of updrafts)
!     2. CALL CCNSPEC  (to convert aerosol/chemistry info to CCN)
!     3. CALL PDFACTIV (to calculate Ndroplet and Smax for known W)
!            or
!        CALL WACTIV   (to calculate W for known Ndroplet)
!
!     Notes:
!     GAULEG calculates the GAUSS integration points and weighting factors
!     This routine needs to be called only if PDF up updraft velocity is
!     employed; if so, GAULEG needs to be called only once, in the INIT
!     routine of the GCM, as the xgs_par, wgs_par values are saved in the COMMON
!     block in 'gmi_par.h'. Npgauss is a PARAMETER defined in the
!     INCLUDE file 'gmi_par.h' and is the # of integration points used
!     to integrate over the PDF of updrafts. When used, the PDF is assumed
!     to be Gaussian.
!
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      subroutine paramdriver (cloudGT, tparc,p_parc, wparc,sulf_p,nact,smax,ityp)
!
      implicit none
!
      type (t_CloudParametersGT), intent(inOut) :: cloudGT
      integer, parameter ::nmdm=3     ! max # of lognormal modes.
      integer  :: nmodes,ityp, nmd
      real*8 :: tpi(nmdm), dpgi(nmdm),  sigi(nmdm),  vhfi(nmdm),  &
     &          amsi(nmdm),densi(nmdm), denii(nmdm), amfsi(nmdm),  &
     &          tpart, npt,nact, nacti, wparc,tparc,pparc,p_parc,  &
     &          accom,sigw,smax,sulfi,sulf_p
!cc
! ** initialize: calculate gauss quadrature points *******************
!
!!      call gauleg (xgs_par, wgs_par, npgauss)
!
! ** specify input for calculating ccn spectrum **********************
!
      if (ityp == 1) wparc = 0.35d0     ! Specify updraft velocity of
      if (ityp == 2) wparc = 1.0d0      ! air parcel
!
!     tparc  = 298.0d0    ! temperature (k)
      pparc  = p_parc*100  !convert mb to (Pa)
      accom  = 0.042d0    ! accommodation coefficient (common for all ccn)
!
      nmodes = 3       ! # of aerosol modes
!     ityp   = 1       ! 1=marine (na), 2=continental, 3=remote marine
      sulfi  = sulf_p  !-1 ! <0, take default conc's; >0, use as scaling factor
!
      call aertyp (tpi, dpgi,  sigi,   vhfi, amsi, densi, denii,  &
     &                  amfsi, nmodes, nmd,  sulfi,ityp)

!
! ** calculate ccn spectrum ****************************************
!
      call ccnspec (cloudGT, tpi,dpgi,sigi,amfsi,vhfi,amsi,densi,denii,  &
     &              tparc,pparc,nmodes)
!
      sigw  = 0.0   ! 0=calculate for single updraft (no pdf)


      call pdfactiv (cloudGT, wparc,sigw,tparc,pparc,nact,smax) ! calculate ndrop

!
! *** calculations complete ***********************************************
!
      return
      end subroutine paramdriver
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!=======================================================================
!
!       Code Developer
!       Nicholas Meskhidze, GA TECH
!       nmeskhidze@eas.gatech.edu
!    -------------------------------
!     DESCRIPTION
!
! *** block data blkpar
! *** this subroutine provides initial (default) values to program
!     parameters via data statements
!
! *** written by athanasios nenes
!
!=======================================================================
!
!      subroutine setData()
!!!
!        implicit none
!!!
!!
!      amw_par     = 18d-3                   ! water molecular weight
!      ama_par     = 29d-3                   ! air molecular weight
!      grav_par    = 9.81d0                  ! g constant
!      rgas_par    = 8.31d0                  ! universal gas constant
!!
!      accom_par = 0.06                     ! default accommodation coef
!!
!      maxit_par   = 100                     ! max iterations for solution
!      eps_par     = 1d-6                    ! convergence criterion
!!
!      pi_par      = 3.1415927d0             ! some constants
!      zero_par    = 0d0
!      great_par   = 1d30
!      sq2pi_par   = 2.5066282746d0
!!
!!     ccnspst    = .false.                  ! internal consistency check
!!
!! *** end of block data subprogram *************************************
!!
!      end subroutine setData
!
!=======================================================================
!
!       Code Developer
!       Nicholas Meskhidze, GA TECH
!       nmeskhidze@eas.gatech.edu
!    -------------------------------
!     DESCRIPTION
!
! *** subroutine ccnspec
! *** this subroutine calculates the ccn spectrum of the aerosol using
!     the appropriate form of kohler theory
!
! *** written by athanasios nenes
!
!=======================================================================
!
      subroutine ccnspec (cloudGT, tpi,dpgi,sigi,amfsi,vhfi,amsi,densi,denii,  &
     &                     tparc,pparc,nmodes)
!
        implicit none
!
      type (t_CloudParametersGT), intent(inOut) :: cloudGT
      integer  :: nmodes,i,j,k
      real*8   :: dpgi(nmodes), vhfi(nmodes), amsi(nmodes),  &
     &            densi(nmodes),sigi(nmodes), tpi(nmodes),  &
     &            amfsi(nmodes),denii(nmodes),  &
     &            tparc, pparc, amfi,denp,vlfs,par1,  &
     &            par2, wparc,smax
!
      cloudGT%nmd_par  = nmodes                ! save aerosol params in common
      do i=1,cloudGT%nmd_par
         cloudGT%dpg_par(i) = dpgi(i)
         cloudGT%vhf_par(i) = vhfi(i)
         cloudGT%ams_par(i) = amsi(i)
         cloudGT%dens_par(i)= densi(i)
         cloudGT%sig_par(i) = sigi(i)
         cloudGT%tp_par(i)  = tpi(i)
         cloudGT%amfs_par(i)= amfsi(i)
         cloudGT%deni_par(i)= denii(i)
      enddo
!
      cloudGT%temp_par = tparc             ! save parcel props in common
      cloudGT%pres_par = pparc
!
      call props(cloudGT)                   ! calculate thermophysical properties
!
! *** calculate critical properties
!
      cloudGT%akoh_par = 4d0*amw_par*cloudGT%surt_par/rgas_par/cloudGT%temp_par/denw_par
!
           ! curvature param
      do k=1,cloudGT%nmd_par
         amfi = max(1.0-cloudGT%amfs_par(k),0.0)                 ! insoluble mass.frac.
         denp = cloudGT%amfs_par(k)*cloudGT%dens_par(k) + amfi*cloudGT%deni_par(k)  ! particle density
         vlfs = cloudGT%amfs_par(k)/cloudGT%dens_par(k)/(cloudGT%amfs_par(k)/cloudGT%dens_par(k)+  &
     &                              amfi/cloudGT%deni_par(k)) ! vol.fr.salt
         par1 = 4d0*denw_par*cloudGT%ams_par(k)/27d0/cloudGT%vhf_par(k)/denp/amw_par/  &
     &                                                 cloudGT%dpg_par(k)**3
         par1 = par1/vlfs                    ! adjust for vol.frac. of salt
         par2 = sqrt(par1*cloudGT%akoh_par**3)
         cloudGT%sg_par(k)= exp(par2) - 1d0              ! sc of dpg
      enddo
!
! *** end of subroutine ccnspec ****************************************
!
      return
      end subroutine ccnspec


!=======================================================================
!
!
!       Code Developer
!       Nicholas Meskhidze, GA TECH
!       nmeskhidze@eas.gatech.edu
!    -------------------------------
!     DESCRIPTION
!
! *** subroutine pdfactiv
! *** this subroutine calculates the ccn activation fraction according
!     to the nenes and seinfeld (2003) parameterization, with
!     modification for non-contunuum effects as proposed by fountoukis
!     and nenes (in preparation). this routine calculates for a pdf of
!     updraft velocities.
!
! *** written by athanasios nenes
!
!=======================================================================
!
      subroutine pdfactiv (cloudGT, wparc,sigw,tparc,pparc,nact,smax)
!
        implicit none
!
!
      type (t_CloudParametersGT), intent(inOut) :: cloudGT
      integer  :: i, isec
      real*8   :: tpart, nact, nacti,wparc,smax
      real*8   :: pdf, dpnmx,sigw,plimt,probi,whi,wlo,  &
     &            scal, wpi,smaxi,tparc,pparc
!
! *** case where updraft is very small
!
      if (wparc.le.1d-6) then
         smax  = 0d0
         nact  = 0d0
         isec  = 1
         dpnmx = great_par
         return
      endif
!
! *** single updraft case
!
      if (sigw.lt.1e-10) then
         call activate (cloudGT, wparc,nact,smax)
         cloudGT%wpdbg(1) = wparc                      ! save debug info
         cloudGT%pddbg(1) = 1.0
         cloudGT%nadbg(1) = nact
         cloudGT%smdbg(1) = smax
!
! *** pdf of updrafts
!
      else
         nact  = zero_par
         smax  = zero_par
         plimt = 1e-3     ! probability of high updraft limit
         probi = sqrt(-2.0*log(plimt*sigw*sq2pi_par))
         whi   = wparc + sigw*probi             ! upper updrft limit
         wlo   = 0.05  ! wparc - sigw*probi     ! low updrft limit
         scal  = 0.5*(whi-wlo)                  ! scaling for updrafts
         do i=1,npgauss
            wpi  = wlo + scal*(1.0-cloudGT%xgs_par(i))      ! updraft
            call activate (cloudGT, wpi,nacti,smaxi)     ! # of drops
            pdf  = (1.0/sq2pi_par/sigw)*exp(-0.5*((wpi-wparc)/  &
     &                                  sigw)**2) !
            nact = nact + cloudGT%wgs_par(i)*(pdf*nacti)    ! integral for drops
            smax = smax + cloudGT%wgs_par(i)*(pdf*smaxi)    ! integral for smax
            cloudGT%wpdbg(i) = wpi                      ! save debug info
            cloudGT%pddbg(i) = pdf
            cloudGT%nadbg(i) = nacti
            cloudGT%smdbg(i) = smaxi
            if (pdf.lt.plimt) goto 100
         enddo
 100     nact = nact*scal                       ! scale integrals
         smax = smax*scal
      endif
!
      return
!
! *** end of subroutine pdfactiv ****************************************
!
      end subroutine pdfactiv


!=======================================================================
!
!
!       Code Developer
!       Nicholas Meskhidze, GA TECH
!       nmeskhidze@eas.gatech.edu
!    -------------------------------
!     DESCRIPTION
!
! *** subroutine wactiv
! *** this subroutine calculates the updraft necessary to achieve a drop
!     concentration.
!
! *** written by athanasios nenes
!
!=======================================================================
!
      subroutine wactiv (cloudGT, nact, wparc, smax)
!
      implicit none
!
      type (t_CloudParametersGT), intent(inOut) :: cloudGT
        integer  :: i,niter
      real*8   :: tpart, nact, nact1, nact2, nact3, wparc, y1, x1,smax,  &
     &            x2,y2,x3,y3,sign,descr
!
! *** initial values for bisection **************************************
!
      x1   = 1e-3          ! low value of updraft
      call activate (cloudGT, x1, nact1, smax)
      y1   = nact1/nact - 1d0
!
      x2   = 20           ! high value of updraft
      call activate (cloudGT, x2, nact2, smax)
      y2   = nact2/nact - 1d0
!
! *** perform bisection *************************************************
!
20    do 30 i=1,maxit_par
         x3   = 0.5*(x1+x2)
         call activate (cloudGT, x3, nact3, smax)
         y3   = nact3/nact - 1d0
!
         if (sign(1.d0,y1)*sign(1.d0,y3) .le. zero_par) then
!                                          ! (y1*y3 .le. zero)
            y2    = y3
            x2    = x3
         else
            y1    = y3
            x1    = x3
         endif
!
         if (abs(x2-x1) .le. eps_par*x1) goto 40
         niter = i
!
30    continue
!
! *** converged ; return ************************************************
!
40    x3    = 0.5*(x1+x2)
      call activate (cloudGT, x3, nact3, smax)
      y3    = nact3/nact - 1d0
      wparc = x3
!
      return
!
! *** end of subroutine wactiv ****************************************
!
      end subroutine wactiv


!=======================================================================
!
!
!       Code Developer
!       Nicholas Meskhidze, GA TECH
!       nmeskhidze@eas.gatech.edu
!    -------------------------------
!     DESCRIPTION
!
! *** subroutine activate
! *** this subroutine calculates the ccn activation fraction according
!     to the nenes and seinfeld (2003) parameterization, with
!     modification for non-contunuum effects as proposed by fountoukis
!     and nenes (in preparation).
!
! *** written by athanasios nenes
!
!=======================================================================
!
      subroutine activate (cloudGT, wparc,ndrpl,smax)

        implicit none
!
!
      type (t_CloudParametersGT), intent(inOut) :: cloudGT
      integer :: i,niter
      real*8  :: ndrpl, wparc, beta,cf1, cf2,x1,sinteg1,sinteg2,y1,  &
     &           x2,y2,x3,y3, sign,smax
!
! *** setup common block variables
!
      cloudGT%wparcel = wparc
!
! *** setup constants
!
      cloudGT%alfa_par = grav_par*amw_par*dhv_par/cpair_par/rgas_par/cloudGT%temp_par  &
     &       /cloudGT%temp_par -grav_par*ama_par/rgas_par/cloudGT%temp_par
      cloudGT%bet1_par = cloudGT%pres_par*ama_par/cloudGT%psat_par/amw_par + amw_par*dhv_par*  &
     &       dhv_par/cpair_par/rgas_par/cloudGT%temp_par/cloudGT%temp_par
      cloudGT%bet2_par = rgas_par*cloudGT%temp_par*denw_par/cloudGT%psat_par/cloudGT%dv_par/amw_par/  &
     &       4d0 +dhv_par*denw_par/4d0/cloudGT%aka_par/cloudGT%temp_par*(dhv_par*  &
     &       amw_par/rgas_par/cloudGT%temp_par - 1d0)
      beta = 0.5d0*pi_par*cloudGT%bet1_par*denw_par/cloudGT%bet2_par/cloudGT%alfa_par/  &
     &       wparc/cloudGT%dair_par
      cf1  = 0.5*(((1/cloudGT%bet2_par)/(cloudGT%alfa_par*wparc))**0.5)
      cf2  = cloudGT%akoh_par/3d0
!
! *** initial values for bisection **************************************
!
      x1   = 1.0d-5   ! min cloud supersaturation -> 0

      call sintegral (cloudGT, x1,ndrpl,sinteg1,sinteg2)
      y1   = (sinteg1*cf1+sinteg2*cf2)*beta*x1 - 1d0
!
      x2   = 1d0      ! max cloud supersaturation = 100%
      call sintegral (cloudGT, x2,ndrpl,sinteg1,sinteg2)
      y2   = (sinteg1*cf1+sinteg2*cf2)*beta*x2 - 1d0
!
! *** perform bisection *************************************************
!
20    do 30 i=1,maxit_par
         x3   = 0.5*(x1+x2)
!
         call sintegral (cloudGT, x3,ndrpl,sinteg1,sinteg2)
         y3 = (sinteg1*cf1+sinteg2*cf2)*beta*x3 - 1d0
!
         if (sign(1.d0,y1)*sign(1.d0,y3) .le. zero_par) then
!                                          ! (y1*y3 .le. zero)
            y2    = y3
            x2    = x3
         else
            y1    = y3
            x1    = x3
         endif
!
         if (abs(x2-x1) .le. eps_par*x1) goto 40
         niter = i
!
30    continue
!
! *** converged ; return ************************************************
!
40    x3   = 0.5*(x1+x2)
!
      call sintegral (cloudGT, x3,ndrpl,sinteg1,sinteg2)
      y3   = (sinteg1*cf1+sinteg2*cf2)*beta*x3 - 1d0
      smax = x3
!
      return
!
! *** end of subroutine activate ****************************************
!
      end subroutine activate


!=======================================================================
!
!
!       Code Developer
!       Nicholas Meskhidze, GA TECH
!       nmeskhidze@eas.gatech.edu
!    -------------------------------
!     DESCRIPTION
!
! *** subroutine sintegral
! *** this subroutine calculates the condensation integrals, according
!     to the population splitting algorithm of nenes and seinfeld (2003)
!     modal formulation according to fountoukis and nenes (2004)
!
! *** written by athanasios nenes
!
!=======================================================================
!
      subroutine sintegral (cloudGT, spar, summa, sum, summat)

        implicit none
!
      type (t_CloudParametersGT), intent(inOut) :: cloudGT
      integer ::j,i,k
      real*8 :: sum, summat, summa, nd(nsmx_par),  &
     &          integ1(nsmx_par),integ2(nsmx_par)
      real*8 :: erf, descr,spar,ratio, ssplt2, ssplt1,  &
     &          sqrt, ssplt, sqtwo, dlgsg,dlgsp, ekth
      real*8 :: orism1, orism2, orism3, orism4, orism5
!
! *** here is where the criterion with the descriminant is put. when it
!     is < 0, then set crit2 = .true. otherwise, set the two values of
!     ssplt and continue.
!
      descr  = 1d0 - (16d0/9d0)*cloudGT%alfa_par*cloudGT%wparcel*cloudGT%bet2_par*  &
     &         (cloudGT%akoh_par/spar**2)**2
!
      if (descr.le.0d0) then
         cloudGT%crit2  = .true.             ! ssplt1,ssplt2 do not exist
         ratio  = (2.0d7/3.0)*cloudGT%akoh_par*spar**(-0.3824)
         if (ratio.gt.1.0) ratio = 1.0
         ssplt2 = spar*ratio
      else
         cloudGT%crit2  = .false.
         ssplt1 = 0.5d0*(1d0-sqrt(descr)) ! min root of both
         ssplt2 = 0.5d0*(1d0+sqrt(descr)) ! max root of both
         ssplt1 = sqrt(ssplt1)*spar       ! multiply ratios with smax
         ssplt2 = sqrt(ssplt2)*spar
      endif
!
      ssplt = ssplt2  ! store ssplit in common
!
! *** calculate integrals
!
      sum       = 0
      summat    = 0
      summa     = 0
!
      sqtwo     = sqrt(2d0)
!
      do 999 j = 1, cloudGT%nmd_par
!
      dlgsg     = dlog(cloudGT%sig_par(j))
      dlgsp     = dlog(cloudGT%sg_par(j)/spar)
      orism1    = 2.d0*dlog(cloudGT%sg_par(j)/ssplt2)/(3.d0*sqtwo*dlgsg )
      orism2    = orism1 - 3.d0*dlgsg/(2.d0*sqtwo)
      orism3    = 2.d0*dlgsp/(3.d0*sqtwo*dlgsg)-3.d0*dlgsg/(2.d0*sqtwo)
      orism4    = orism1 + 3.d0*dlgsg/sqtwo
      orism5    = 2.d0*dlgsp/(3*sqtwo*dlgsg)
      ekth      = exp(9d0/2d0*dlgsg*dlgsg)
      integ1(j) = cloudGT%tp_par(j)*spar*((1-erf(orism1)) -  &
     &            0.5d0*((cloudGT%sg_par(j)/spar)**2)*ekth*(1-erf(orism4)))
      integ2(j) = (exp(9d0/8d0*dlgsg*dlgsg)*cloudGT%tp_par(j)/cloudGT%sg_par(j))*  &
     &            (erf(orism2) - erf(orism3))
!
! *** calculate number of drops
!
      nd(j)     = (cloudGT%tp_par(j)/2.0)*(1.0-erf(orism5))
!
      sum       = sum    + integ1(j)
      summat    = summat + integ2(j)
      summa     = summa  + nd(j)
 999  continue
!
      return
      end subroutine sintegral

!=======================================================================
!
!
!       Code Developer
!       Nicholas Meskhidze, GA TECH
!       nmeskhidze@eas.gatech.edu
!    -------------------------------
!     DESCRIPTION
!
! *** subroutine props
! *** this subroutine calculates the thermophysical properties
!
! *** written by athanasios nenes
!
!=======================================================================
!
      subroutine props(cloudGT)
!
      implicit none
!  
      type (t_CloudParametersGT), intent(inOut) :: cloudGT
!      real*8  :: vpres, sft,presa,dbig,dlow,coef
      real*8  :: presa,dbig,dlow,coef
!
!      denw_par  = 1d3                         ! water density
!      dhv_par   = 2.25d6                      ! water enthalpy of vaporization
!      cpair_par = 1.0061d3                    ! air cp
      presa = cloudGT%pres_par/1.013d5                ! pressure (pa)
      cloudGT%dair_par  = cloudGT%pres_par*ama_par/rgas_par/cloudGT%temp_par            ! air density
!
      cloudGT%aka_par   = (4.39+0.071*cloudGT%temp_par)*1d-3  ! air thermal conductivity
!
      cloudGT%dv_par    = (0.211d0/presa)*(cloudGT%temp_par/273d0)**1.94
      cloudGT%dv_par    = cloudGT%dv_par*1d-4                 ! water vapor diffusivity in air
      dbig  = 5.0d-6
      dlow  = 0.207683*((accom_par)**(-0.33048))
      dlow  = dlow*1d-6
!
! dv average
!
      coef  = ((2*pi_par*amw_par/(rgas_par*cloudGT%temp_par))**0.5)
!
      cloudGT%dv_par    = (cloudGT%dv_par/(dbig-dlow))*((dbig-dlow)-(2*cloudGT%dv_par/accom_par)  &
     &        *coef*(dlog((dbig+(2*cloudGT%dv_par/accom_par)*coef)/(dlow+  &
     &        (2*cloudGT%dv_par/accom_par)*coef))))             ! non-continuum effects
!
      cloudGT%psat_par  = vpres(cloudGT%temp_par)*(1e5/1.0d3) ! saturation vapor pressure
!
      cloudGT%surt_par  = sft(cloudGT%temp_par)             ! surface tension for water (j m-2)
!
      return
!
! *** end of subroutine props *******************************************
!
      end subroutine props



!=======================================================================
!
! *** function vpres
! *** this function calculates saturated water vapour pressure as a
!     function of temperature. valid for temperatures between -50 and
!     50 c.
!
! ======================== arguments / usage ===========================
!
!  input:
!     [t]
!     real variable.
!     ambient temperature expressed in kelvin.
!
!  output:
!     [vpres]
!     real variable.
!     saturated vapor pressure expressed in mbar.
!
!=======================================================================
!
      real*8 function vpres (t)
!
      implicit none
!
      integer  ::i
      real*8 :: a(0:6), t,ttemp
      data a/6.107799610e+0, 4.436518521e-1, 1.428945805e-2,  &
     &       2.650648471e-4, 3.031240396e-6, 2.034080948e-8,  &
     &       6.136820929e-11/
!
! calculate polynomial (without exponentiation).
!
      ttemp = t-273.0d0
      vpres = a(6)*ttemp
      do i=5,1,-1
         vpres = (vpres + a(i))*ttemp
      enddo
      vpres = vpres + a(0)
!
! end of function vpres
!
      return
      end function vpres



!=======================================================================
!
! *** function sft
! *** this function calculates water surface tension as a
!     function of temperature. valid for temperatures between -40 and
!     40 c.
!
! ======================== arguments / usage ===========================
!
!  input:
!     [t]
!     real variable.
!     ambient temperature expressed in kelvin.
!
!  output:
!     [sft]
!     real variable.
!     surface tension expressed in j m-2.
!
!=======================================================================
!
      real*8 function sft (t)
!
      implicit none
!
      real*8 :: t,tpars
!
      tpars = t-273.15d0
      sft   = 0.0761-1.55e-4*tpars
!
      return
      end function sft


! ***********************************************************************
!
      subroutine gauleg (x,w,n)
!
! calculation of points and weights for n point gauss integration
! ***********************************************************************
      implicit none
!
      integer           :: n,m,i,j
      real*8            :: x(n), w(n),xm,xl,z,p1,p2,p3,pp,z1
      real*8, parameter :: eps_par1=1.e-6
      real*8, parameter :: x1=-1.0, x2=1.0
!
! calculation
!
      m=(n+1)/2d0
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.eps_par1)go to 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      end subroutine gauleg

!=======================================================================
!
! *** real function erf
! *** this subroutine calculates the error function
!
! *** obtained from numerical recipies
!
!=======================================================================
!
      real*8 function erf(x)
!
      implicit none
!
        real*8  :: x
!
      if(x.lt.0.)then
        erf=-gammp(.5d0,x**2)
      else
        erf=gammp(.5d0,x**2)
      endif
      return
      end function erf

!
!=======================================================================
!
      real*8 function gammln(xx)
!
!=======================================================================
!
      implicit none
!
      integer  :: j
      real*8 :: cof(6),stp,half,one,fpf,x,tmp,ser,xx
!
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,  &
     &    -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do 11 j=1,6
        x=x+one
        ser=ser+cof(j)/x
11    continue
      gammln=tmp+log(stp*ser)
      return
      end function gammln


!
!=======================================================================
!
      real*8 function gammp(a,x)
!
!=======================================================================
!
      implicit none
!
      real*8 a,x,gln,gamser,gammcf

      if(x.lt.0.d0.or.a.le.0.d0)pause
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.d0-gammcf
      endif
      return
      end function gammp

!
!=======================================================================
!
      subroutine gcf(gammcf,a,x,gln)
!
!=======================================================================
!
      implicit none
!
      integer            :: n
      real*8, parameter  :: itmax=100,eps_par1=3.e-7
!        real*8             :: gln, gammln,gold,a0,a1,x,b0,b1,fac,an,  &
        real*8             :: gln, gold,a0,a1,x,b0,b1,fac,an,  &
     &                      float,ana,a,anf,g,gammcf
      gln=gammln(a)
      gold=0.
      a0=1.
      a1=x
      b0=0.
      b1=1.
      fac=1.
      do 11 n=1,itmax
        an=float(n)
        ana=an-a
        a0=(a1+a0*ana)*fac
        b0=(b1+b0*ana)*fac
        anf=an*fac
        a1=x*a0+anf*a1
        b1=x*b0+anf*b1
        if(a1.ne.0.)then
          fac=1./a1
          g=b1*fac
          if(abs((g-gold)/g).lt.eps_par1)go to 1
          gold=g
        endif
11    continue
      pause 'a too large, itmax too small'
1     gammcf=exp(-x+a*log(x)-gln)*g
      return
      end subroutine gcf


!
!=======================================================================
!
      subroutine gser(gamser,a,x,gln)
!
!=======================================================================
!
      implicit none
!
      integer  :: n
      real*8, parameter ::  itmax=100,eps_par1=3.e-7
!        real*8  :: gln, gammln,x,gamser,a,ap,sum,del,abs
        real*8  :: gln, x,gamser,a,ap,sum,del,abs

      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,itmax
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*eps_par1)go to 1
11    continue
      pause 'a too large, itmax too small'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      end subroutine gser

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!=======================================================================
!
!
!       Code Developer
!       Nicholas Meskhidze, GA TECH
!       nmeskhidze@eas.gatech.edu
!    -------------------------------
!     DESCRIPTION
!
! *** subroutine aertyp
! *** lognormal distribution info, scaled to sulfate mass
!
! *** written by athanasios nenes
!
!=======================================================================
!
      subroutine aertyp (tpi,   dpgi, sigi, vhfi, amsi, densi, denii,  &
     &                  amfsi, nmdm, nmd,  sulfi,ityp)
!
      implicit none
!
      integer  :: i, ityp, nmdm,nmd
      real*8   :: tpi(nmdm), dpgi(nmdm),  sigi(nmdm),  vhfi(nmdm),  &
     &                 amsi(nmdm),densi(nmdm), denii(nmdm), amfsi(nmdm),  &
     &                 sulfi, amts
!
! *********************************************************************
! *** "polluted" (north atlantic) marine aerosol
! *********************************************************************
!
      if (ityp.eq.1) then
!
! nucleation mode
!
      tpi  (1) = 230d6      ! total concentration (# m-3)
      dpgi (1) = 0.02d-6    ! modal diameter (m)
      sigi (1) = 1.47d0     ! geometric dispersion (sigma_g)
      densi(1) = 1760.0     ! density of soluble fraction (kg m-3)
      denii(1) = 2100.0     ! density of insoluble fraction (kg m-3)
      amfsi(1) = 0.33       ! soluble mass fraction
      amsi (1) = 1.3200e-01 ! molar mass of soluble fraction (kg mol-1)
      vhfi (1) = 3.0        ! van't hoff factor for soluble frac.(ions molec-1)
!
! accumulation mode
!
      tpi  (2) = 176.7d6    ! total concentration (# m-3)
      dpgi (2) = 0.092d-6   ! modal diameter (m)
      sigi (2) = 1.6d0      ! geometric dispersion (sigma_g)
      densi(2) = 1760.0     ! density of soluble fraction (kg m-3)
      denii(2) = 2100.0     ! density of insoluble fraction (kg m-3)
      amfsi(2) = 0.33       ! soluble mass fraction
      amsi (2) = 1.3200e-01 ! molar mass of soluble fraction (kg mol-1)
      vhfi (2) = 3.0        ! van't hoff factor for soluble frac.(ions molec-1)
!
! coarse mode
!
      tpi  (3) = 3.1d6      ! total concentration (# m-3)
      dpgi (3) = 0.58d-6    ! modal diameter (m)
      sigi (3) = 2.49d0     ! geometric dispersion (sigma_g)
      densi(3) = 1760.0     ! density of soluble fraction (kg m-3)
      denii(3) = 2100.0     ! density of insoluble fraction (kg m-3)
      amfsi(3) = 0.95       ! soluble mass fraction
      amsi (3) = 1.3200e-01 ! molar mass of soluble fraction (kg mol-1)
      vhfi (3) = 3.0        ! van't hoff factor for soluble frac.(ions molec-1)
!
! cumulative properties
!
      amts     = 1.75
      nmd      = 3
!
! *********************************************************************
! *** continental aerosol
! *********************************************************************
!
      else if (ityp.eq.2) then
!
! nucleation mode
!
      tpi  (1) = 1000d6     ! total concentration (# m-3)
      dpgi (1) = 0.016d-6   ! modal diameter (m)
      sigi (1) = 1.6d0      ! geometric dispersion (sigma_g)
      densi(1) = 1760.0     ! density of soluble fraction (kg m-3)
      denii(1) = 2100.0     ! density of insoluble fraction (kg m-3)
      amfsi(1) = 0.5        ! soluble mass fraction
      amsi (1) = 1.3200e-01 ! molar mass of soluble fraction (kg mol-1)
      vhfi (1) = 3.0        ! van't hoff factor for soluble frac.(ions molec-1)
!
! accumulation mode
!
      tpi  (2) = 800d6      ! total concentration (# m-3)
      dpgi (2) = 0.067d-6   ! modal diameter (m)
      sigi (2) = 2.1d0      ! geometric dispersion (sigma_g)
      densi(2) = 1760.0     ! density of soluble fraction (kg m-3)
      denii(2) = 2100.0     ! density of insoluble fraction (kg m-3)
      amfsi(2) = 0.5        ! soluble mass fraction
      amsi (2) = 1.3200e-01 ! molar mass of soluble fraction (kg mol-1)
      vhfi (2) = 3.0        ! van't hoff factor for soluble frac.(ions molec-1)
!
! coarse mode
!
      tpi  (3) = 0.72d6     ! total concentration (# m-3)
      dpgi (3) = 0.93d-6    ! modal diameter (m)
      sigi (3) = 2.2d0      ! geometric dispersion (sigma_g)
      densi(3) = 1760.0     ! density of soluble fraction (kg m-3)
      denii(3) = 2100.0     ! density of insoluble fraction (kg m-3)
      amfsi(3) = 0.5        ! soluble mass fraction
      amsi (3) = 1.3200e-01 ! molar mass of soluble fraction (kg mol-1)
      vhfi (3) = 3.0        ! van't hoff factor for soluble frac.(ions molec-1)
!
! cumulative properties
!
      amts     = 4.46
      nmd      = 3
!
! *********************************************************************
! *** "clean" (southern oceans) marine aerosol
! *********************************************************************
!
      else if (ityp.eq.3) then
!
! nucleation mode
!
      tpi  (1) = 310d6      ! total concentration (# m-3)
      dpgi (1) = 0.018d-6   ! modal diameter (m)
      sigi (1) = 1.4d0      ! geometric dispersion (sigma_g)
      densi(1) = 1760.0     ! density of soluble fraction (kg m-3)
      denii(1) = 2100.0     ! density of insoluble fraction (kg m-3)
      amfsi(1) = 0.33       ! soluble mass fraction
      amsi (1) = 1.3200e-01 ! molar mass of soluble fraction (kg mol-1)
      vhfi (1) = 3.0        ! van't hoff factor for soluble frac.(ions molec-1)
!
! accumulation mode
!
      tpi  (2) = 70d6       ! total concentration (# m-3)
      dpgi (2) = 0.075d-6   ! modal diameter (m)
      sigi (2) = 1.6d0      ! geometric dispersion (sigma_g)
      densi(2) = 1760.0     ! density of soluble fraction (kg m-3)
      denii(2) = 2100.0     ! density of insoluble fraction (kg m-3)
      amfsi(2) = 0.33       ! soluble mass fraction
      amsi (2) = 1.3200e-01 ! molar mass of soluble fraction (kg mol-1)
      vhfi (2) = 3.0        ! van't hoff factor for soluble frac.(ions molec-1)
!
! coarse mode
!
      tpi  (3) = 3.1d6      ! total concentration (# m-3)
      dpgi (3) = 0.62d-6    ! modal diameter (m)
      sigi (3) = 2.7d0      ! geometric dispersion (sigma_g)
      densi(3) = 1760.0     ! density of soluble fraction (kg m-3)
      denii(3) = 2100.0     ! density of insoluble fraction (kg m-3)
      amfsi(3) = 0.95       ! soluble mass fraction
      amsi (3) = 1.3200e-01 ! molar mass of soluble fraction (kg mol-1)
      vhfi (3) = 3.0        ! van't hoff factor for soluble frac.(ions molec-1)
!
! cumulative properties
!
      amts     = 1.75
      nmd      = 3
      endif
!
! *********************************************************************
! *** process distributions and return
! *********************************************************************
!
      if (sulfi.lt.0d0) then
         sulfi = amts
      else
         do i=1,nmd
            tpi(i) = tpi(i)*sulfi/amts    ! scale dist. by sulfate mass
         enddo
      endif
!
      return
!
      end subroutine aertyp

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!=======================================================================

!
!       Code Developer
!       Nicholas Meskhidze, GA TECH
!       nmeskhidze@eas.gatech.edu
!    -------------------------------
!     DESCRIPTION
! *** subroutine wcalcbl
! *** calculate updraft velocities using lance et al. (2004) algorithm
!
! *** written by athanasios nenes
!
!=======================================================================
!
      subroutine wcalcbl (sulfi,ityp,wparc)
        implicit none
!
!     ------------------------
!     Variable declaration
!     ------------------------
      integer  :: i
      real*8   :: sulfi, wparc
      integer  :: ityp
      real*8   :: so4(20), wmarplt(20), wcont(20), wmarcln(20),  &
     &            delso4, s4coef, delw

      data so4/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0,4.0,  &
     &         5.0,6.0,7.0,8.0,9.0,10.0,20.0/

      data wmarplt/19.99999046,19.99999046,7.962497371,4.928205238,  &
     &             3.918289601,3.334171832,2.914480751,2.572751031,  &
     &             2.270527083,1.987856653,0.602901416,0.44409889,  &
     &             0.387454362,0.357389792,0.338374179,0.325068108,  &
     &             0.315116418,0.307310984,0.300964896,0.268002568/

      data wcont  /19.99999046,19.99999046,9.073062339,5.772469968,  &
     &             4.37407545,3.533317682,2.944610585,2.498520734,  &
     &             2.146460258,1.863475131,0.770942349,0.53355026,  &
     &             0.433430788,0.37678775,0.339495286,0.312558016,  &
     &             0.291826546,0.27509919,0.261080878,0.147503164/

      data wmarcln/19.99999046,19.99999046,16.76764915,8.263419023,  &
     &             6.85000131,6.189548113,5.795380848,5.529094916,  &
     &             5.334755344,5.185088412,4.499296164,4.159116087,  &
     &             3.857820137,3.536240548,3.169981167,2.755074907,  &
     &             2.332597448,1.981314774,1.730071109,1.025660195/
!
      if (sulfi.le.so4(1)) then   ! so4 < minimum limit
         if (ityp.eq.1) then
            wparc = wmarplt(1)
         else if (ityp.eq.2) then
            wparc = wcont(1)
         else
            wparc = wmarcln(1)
         endif
!
      elseif (sulfi.ge.so4(20)) then  ! so4>maximum limit
         if (ityp.eq.1) then
            wparc = wmarplt(20)
         else if (ityp.eq.2) then
            wparc = wcont(20)
         else
            wparc = wmarcln(20)
         endif
!
      else                           ! retrieve the w with interpolation
         do i=2,20
            if (so4(i-1).lt.sulfi .and. sulfi.le.so4(i)) then
               delso4 = so4(i)-so4(i-1)
               s4coef = sulfi -so4(i-1)
               if (ityp.eq.1) then
                  delw  = wmarplt(i)-wmarplt(i-1)
                  wparc = wmarplt(i-1) + s4coef*delw/delso4
               else if (ityp.eq.2) then
                  delw  = wcont(i)-wcont(i-1)
                  wparc = wcont(i-1) + s4coef*delw/delso4
               else
                  delw  = wmarcln(i)-wmarcln(i-1)
                  wparc = wmarcln(i-1) + s4coef*delw/delso4
               endif
               goto 10
            endif
         enddo
      endif
!
  10  continue
      return
      end subroutine wcalcbl


!  **********************============================**************************
!
!  END NIKA
!c*************************************************************************
!
!
! Begin rsot
!=======================================================================
!
!       Code Developer
!       Rafaella Sotiropoulou, GATECH
!       rsot@eas.gatech.edu
!    -------------------------------
! *** SUBROUTINE NdBL
! *** CALCULATE DROPLET CONCENTRATIONS USING BOUCHER AND LOHMANN CORRELATION
!
! *** WRITTEN BY ATHANASIOS NENES
!
!=======================================================================
!
      subroutine NdBL (sulfi,ityp, ndrop)
!
      implicit none
!
      integer  :: ityp
      real*8   :: sulfi, ndrop
!
      if (ityp.eq.1 .or. ityp.eq.3) then       ! Marine
         ndrop = 0.48d0*dlog10(sulfi) + 2.06d0
      else                                     ! Continental (average)
         ndrop = 0.41d0*dlog10(sulfi) + 2.21d0
      endif
!
      ndrop = 10d0**ndrop         ! cm-3
      ndrop = MAX(ndrop,40d0)     ! Restrict Minimum droplet number
!
      return
      end subroutine NdBL

!*************************************************************************
!
! DROPLET ACTIVATION PARAMETERIZATION: Abdul-Razzak and Ghan, 1998
!
! ************************************************************************
!
!
!       Code Developer
!       Rafaella Sotiropoulou, GATECH
!       rsot@eas.gatech.edu
!    -------------------------------
! *** CALCULATE DROPLET CONCENTRATIONS USING Abdul-Razzak & Ghan parameterization
!
! *** PROGRAM PARAMAG
!
!=======================================================================
!
      subroutine paramAG (cloudGT, tparc,p_parc, wparc,sulf_p,ndrops,SmaxAG,ityp)
!
      implicit none
!
      type (t_CloudParametersGT), intent(inOut) :: cloudGT
      integer, parameter ::nmdm=3     ! max # of lognormal modes.
      integer  :: nmodes,ityp, nmd
      real*8 :: tpi(nmdm), dpgi(nmdm),  sigi(nmdm),  vhfi(nmdm),  &
     &          amsi(nmdm),densi(nmdm), denii(nmdm), amfsi(nmdm),  &
     &          tpart, npt, ndrops, nacti, wparc,tparc,pparc,p_parc,  &
     &          accom,sigw,SmaxAG,sulfi,sulf_p
!
! ** specify input for calculating ccn spectrum **********************
!
      if (ityp == 1) wparc = 0.35d0     ! Specify updraft velocity of
      if (ityp == 2) wparc = 1.0d0      ! air parcel
!
      pparc  = p_parc*100.d0  !convert mb to (Pa)
!
      nmodes = 3       ! # of aerosol modes
      sulfi  = sulf_p  !-1 ! <0, take default conc's; >0, use as scaling factor
!
!
      call aertyp (tpi, dpgi,  sigi,   vhfi, amsi, densi, denii,  &
     &                  amfsi, nmodes, nmd,  sulfi,ityp)
!
! ** calculate ccn spectrum ****************************************
!
      call ccnspcAG (cloudGT, tpi,dpgi,sigi,amfsi,vhfi,amsi,densi,denii,  &
     &              tparc,pparc,nmodes)
!
! ** calculate Ndroplet and Smax for known W ****************************
!
      call dropform(cloudGT, WPARC,TPARC,PPARC,ndrops,SmaxAG)

!
!!rsot	ndrops=ndrops*1d-6 !Convert to cm-3
!
! *** calculations complete ***********************************************
!
      return
      end subroutine paramAG
!=======================================================================
!
!       Code Developer
!       Rafaella Sotiropoulou, GATECH
!       rsot@eas.gatech.edu
!    -------------------------------
!
!     DESCRIPTION
!
! *** subroutine ccnspcAG
! *** this subroutine calculates the ccn spectrum of the aerosol using
!     the appropriate form of kohler theory
!
!
!=======================================================================
!
      subroutine ccnspcAG (cloudGT, tpi,dpgi,sigi,amfsi,vhfi,amsi,densi,denii,  &
     &                     tparc,pparc,nmodes)
!
        implicit none
!
      type (t_CloudParametersGT), intent(inOut) :: cloudGT
      integer  :: nmodes,i,j,k
      real*8   :: dpgi(nmodes), vhfi(nmodes), amsi(nmodes),  &
     &            densi(nmodes),sigi(nmodes), tpi(nmodes),  &
     &            amfsi(nmodes),denii(nmodes), beta(nmodes),  &
     &            tparc, pparc, amfi,denp,vlfs,par1,  &
     &            par2, wparc

!
      cloudGT%nmd_par  = nmodes                ! save aerosol params in common
      do i=1,cloudGT%nmd_par
         cloudGT%dpg_AG(i) = dpgi(i)
         cloudGT%vhf_AG(i) = vhfi(i)
         cloudGT%ams_AG(i) = amsi(i)
         cloudGT%dens_AG(i)= densi(i)
         cloudGT%sig_AG(i) = sigi(i)
         cloudGT%tp_AG(i)  = tpi(i)
         cloudGT%amfs_AG(i)= amfsi(i)
         cloudGT%deni_AG(i)= denii(i)
      enddo
!
      cloudGT%temp_par = tparc               ! save parcel propsAG in common
      cloudGT%pres_par = pparc
!
      call propsAG(cloudGT)                   ! calculate thermophysical properties
!
! *** calculate critical properties
!
      cloudGT%akoh_AG = 2d0*amw_par*cloudGT%surt_par/rgas_par/cloudGT%temp_par/denw_par   ! A
!
           ! curvature param
      do k=1,cloudGT%nmd_par
         amfi = max(1.d0-cloudGT%amfs_AG(k),0.d0)                 ! insoluble mass.frac.
         denp = cloudGT%amfs_AG(k)*cloudGT%dens_AG(k) + amfi*cloudGT%deni_AG(k)   ! particle density
   	 beta(k) = vhfi(k)*amw_par*cloudGT%amfs_AG(k)*denp/amsi(k)/denw_par
	 cloudGT%sg_AG(k) = 2.d0/sqrt(beta(k))*(2.d0*cloudGT%akoh_AG/3.d0/dpgi(k))**1.5d0
      enddo

!
! *** end of subroutine ccnspcAG ****************************************
!
      return
      end subroutine ccnspcAG
!=======================================================================
!
!       Code Developer
!       Rafaella Sotiropoulou, GATECH
!       rsot@eas.gatech.edu
!    -------------------------------
!
!     DESCRIPTION
!
! *** subroutine propsAG
! *** this subroutine calculates the thermophysical properties
!
! *** written by athanasios nenes
!
!=======================================================================
!
      subroutine propsAG(cloudGT)
!
      implicit none
!
      type (t_CloudParametersGT), intent(inOut) :: cloudGT
!      real*8  :: vpres, sft,presa,coef,arga,arlhv
      real*8  :: presa,coef,arga,arlhv
!
!      denw_par  = 1d3						 ! water density
!
      ARGA  = 0.167d0 + 3.67d-4*cloudGT%temp_par
      ARLHV = 597.3d0*(273.15d0/cloudGT%temp_par)**ARGA
      cloudGT%DHV_AG = ARLHV/2.38844d-8*1d-4				 ! Enthalpy of evaporization (J/kg)
!!!!!!	  dhv_par   = 2.25d6					 ! water enthalpy of vaporization
!      cpair_par = 1.0061d3					 ! air cp
      presa = cloudGT%pres_par/1.013d5					 ! pressure (pa)
      cloudGT%dair_par  = cloudGT%pres_par*ama_par/rgas_par/cloudGT%temp_par             ! air density
!
      cloudGT%aka_par   = (4.39d0+0.071d0*cloudGT%temp_par)*1d-3		 ! air thermal conductivity
!
      cloudGT%dv_AG     = (0.211d0/presa)*(cloudGT%temp_par/273d0)**1.94d0
      cloudGT%dv_AG     = cloudGT%dv_AG*1d-4					 ! water vapor diffusivity in air
      cloudGT%psat_par  = vpres(cloudGT%temp_par)*(1d5/1.0d3)		         ! saturation vapor pressure
!
      cloudGT%surt_par  = sft(cloudGT%temp_par)					 ! surface tension for water (j m-2)
!
      return
!
! *** end of subroutine propsAG *******************************************
!
      end subroutine propsAG

!=======================================================================
!
!       Code Developer
!       Rafaella Sotiropoulou, GATECH
!       rsot@eas.gatech.edu
!    -------------------------------
!     DESCRIPTION
!
! *** subroutine dropform
! *** this subroutine calculates the ccn activation fraction according
!     to the Abdul-Razzak and Ghan (1998) parameterization
!
!=======================================================================
!
      subroutine dropform (cloudGT, WPARC,TPARC,PPARC,ndrops,SmaxAG)

        implicit none
!
!
      type (t_CloudParametersGT), intent(inOut) :: cloudGT
      integer :: i
      real*8  :: ndrops, wparc, SmaxAG, erf, NdAG(cloudGT%nmd_par)
      real*8  :: zetaAG, etaAG(cloudGT%nmd_par), fAG(cloudGT%nmd_par), gAG(cloudGT%nmd_par), u(cloudGT%nmd_par)
      real*8  :: tparc, pparc, sf1(cloudGT%nmd_par), sf2(cloudGT%nmd_par), sf3(cloudGT%nmd_par), sf4(cloudGT%nmd_par),sf5, smax
!
! *** setup common block variables
!
      cloudGT%wparcel = wparc
!
! *** setup constants
!
      cloudGT%alfa_AG = grav_par*amw_par*cloudGT%dhv_AG/cpair_par/rgas_par/(cloudGT%temp_par**2d0)  &
     &         - grav_par*ama_par/rgas_par/cloudGT%temp_par						! alpha Eq. 11
      cloudGT%bet1_AG = rgas_par*cloudGT%temp_par/cloudGT%psat_par/amw_par + amw_par*(cloudGT%dhv_AG**2.d0) &
     &         / cpair_par/cloudGT%pres_par/ama_par/cloudGT%temp_par				        	! gamma Eq. 12
      cloudGT%bet2_AG = denw_par*rgas_par*cloudGT%temp_par/cloudGT%psat_par/cloudGT%dv_AG/amw_par/  &
     &           1d0 + cloudGT%dhv_AG*denw_par/1d0/cloudGT%aka_par/cloudGT%temp_par*(cloudGT%dhv_AG*  &
     &           amw_par/rgas_par/cloudGT%temp_par - 1d0)		                                ! 1/G Eq. 16

 ! *** calculate Smax for the Abdul-Razzak and Ghan (1998) parameterization

          zetaAG= 2.d0*cloudGT%akoh_AG/3.d0*(cloudGT%alfa_AG*wparc*cloudGT%bet2_AG)**0.5d0

	  SmaxAG = 0.d0
	  fAG(i) = 0.d0
	  gAG(i) = 0.d0
	  u(i)   = 0.d0
    	  sf1(i) = 0.d0
    	  sf2(i) = 0.d0
    	  sf3(i) = 0.d0
    	  sf4(i) = 0.d0
    	  sf5    = 0.d0
	  smax   = 0.d0

	do i=1,cloudGT%nmd_par
	    etaAG(i)  = ((cloudGT%alfa_AG*wparc*cloudGT%bet2_AG)**1.5d0)/2.d0/pi_par/denw_par/cloudGT%bet1_AG/cloudGT%tp_AG(i)
	    fAG(i) = 0.5d0*dexp(2.5d0*dlog(cloudGT%sig_AG(i))*dlog(cloudGT%sig_AG(i)))
	    gAG(i) = 1.d0 + 0.25d0*dlog(cloudGT%sig_AG(i))


!		  sf1(i)= fAG(i)*sqrt(zetaAG/etaAG(i))**3.d0

!		  sf2(i)= etaAG(i)+3d0*zetaAG
!		  sf3(i)= sg_par(i)*sg_par(i)

!		  sf4(i)=1.d0/sf3(i)*(sf1(i)+gAG(i)*(sf3(i)/sf2(i))**(3.d0/4.d0))
!		  sf5=sf5+sf4(i)


	    SmaxAG = SmaxAG + 1.d0/cloudGT%sg_AG(i)/cloudGT%sg_AG(i)* (fAG(i)*(zetaAG/etaAG(i))**1.5d0 + &
     &	    gAG(i)*(cloudGT%sg_AG(i)*cloudGT%sg_AG(i)/(etaAG(i)+3.d0*zetaAG))**(3.d0/4.d0))

        enddo
!	     smax = smax + 1.d0/sqrt(sf5)
	     SmaxAG = 1d0/sqrt(SmaxAG)
      
      do i=1,cloudGT%nmd_par
	     u(i)= 2.d0*dlog(cloudGT%sg_AG(i)/SmaxAG)/3.d0/sqrt(2.d0)/dlog(cloudGT%sig_AG(i))
      enddo

	Ndrops = 0.d0

! *** calculate the ccn activation fraction according
!     to the Abdul-Razzak and Ghan (1998) parameterization
!
      do i=1,cloudGT%nmd_par
	  NdAG(i) = cloudGT%tp_AG(i)/2.d0*(1.d0-erf(u(i)))
	    if (NdAG(i) .gt. cloudGT%tp_AG(i)) NdAG(i) = cloudGT%tp_AG(i)
	  Ndrops = Ndrops + NdAG(i)
      enddo

      return
!
! *** end of subroutine dropform ****************************************
!
      end subroutine dropform
!=======================================================================
!*************************************************************************
!
! DROPLET ACTIVATION PARAMETERIZATION: Segal and Khain, 2006
!
! ************************************************************************
!       Code Developer
!       Rafaella Sotiropoulou, GATECH
!       rsot@eas.gatech.edu
!    -------------------------------
! *** CALCULATE DROPLET CONCENTRATIONS USING Segal & Khain parameterization
!
! *** PROGRAM PARAMSK
!
!=======================================================================
!
      subroutine paramSK (wparc,sulf_p,dropSK,ityp)
!
      implicit none
!
      integer, parameter ::nmdm=3     ! max # of lognormal modes.
      integer  :: nmodes,ityp, nmd, i
      real*8 :: tpi(nmdm), dpgi(nmdm),  sigi(nmdm),  vhfi(nmdm),  &
     &          amsi(nmdm),densi(nmdm), denii(nmdm), amfsi(nmdm),  &
     &          tpart, npt, dropSK, nacti, wparc,tparc,pparc,p_parc,  &
     &          sigw,sulfi,sulf_p,wparcav,sigwpar,avgNcn,sigNcn,Wi, &
     &          avgWi,sigWi,b0,b1,b2,b3,b11,b22,b33,b12,b13,b23,gNcn, &
     &          AWcb, ANcn, AWi
!
      pparc  = p_parc*100.d0  !convert mb to (Pa)
!
      nmodes = 3       ! # of aerosol modes
      sulfi  = sulf_p  !-1 ! <0, take default conc's; >0, use as scaling factor
!      sulfi  = -1 ! <0, take default conc's; >0, use as scaling factor
!
!
      call aertyp (tpi, dpgi,  sigi,   vhfi, amsi, densi, denii,  &
     &                  amfsi, nmodes, nmd,  sulfi,ityp)
!
! ******************** calculate Ndroplet for known W *********************
!
! **  Calculate the aerosol number
      gNcn =0.d0
!
      do i=1, nmodes
         gNcn=gNcn+amfsi(i)*tpi(i)*1.0d-6            !convert to cm-3
      enddo

	IF (ITYP == 1) THEN
        IF (gNcn .le. 40.d0) goto 1000
	ELSE IF (ITYP == 2) THEN
        IF (gNcn .le. 75.d0) goto 1000
	ENDIF

      IF (ityp == 1) THEN   ! Marine

! ** specify input for calculating droplet number **********
	 IF (gNcn .gt. 100.0d0) THEN        ! Polluted Marine (Region 2)
             wparc    =   0.35d0
             wparcav  =    1.0d0
             sigwpar  =    0.5d0
             avgNcn   =   400.d0
             sigNcn   =0.30103d0
             Wi       =dlog10(1.6d0)
             avgWi    =    0.4d0
             sigWi    =  0.075d0
             b0       =  252.1d0
             b1       =   46.8d0
             b2       =  145.2d0
             b3       =  -29.5d0
             b11      =  -13.0d0
             b22      =   38.1d0
             b33      =   -6.1d0
             b12      =   30.7d0
             b13      =   -0.9d0
             b23      =  -18.7d0
	 ELSE IF (gNcn .gt. 40.0d0)   THEN
              IF (gNcn .le. 100.0d0)  THEN   ! Clean Marine (Region 1)
             wparc    =   0.35d0
             wparcav  =    0.5d0
             sigwpar  =   0.25d0
             avgNcn   =   100.d0
             sigNcn   =0.30103d0
             Wi       =dlog10(1.6d0)
             avgWi    =  0.225d0
             sigWi    =  0.075d0
             b0       =   67.2d0
             b1       =   10.7d0
             b2       =   41.7d0
             b3       =   -6.1d0
             b11      =   -2.8d0
             b22      =   12.0d0
             b33      =   -1.4d0
             b12      =    7.2d0
             b13      =    0.1d0
             b23      =   -3.9d0
	      ENDIF
	 ENDIF
!
         AWcb	= (wparc-wparcav)/sigwpar
         ANcn	= (dlog10(gNcn)-dlog10(avgNcn))/sigNcn
         AWi	= (Wi-avgWi)/sigWi
!
         dropSK =  b0 + b1*AWcb + b2*ANcn + b3*AWi+    &
     &             b11*AWcb*AWcb + b22*ANcn*ANcn + b33*AWi*AWi + &
     &             b12*AWcb*ANcn + b13*AWcb*AWi  + b23*ANcn*AWi
!
!
      ELSE IF (ityp == 2) THEN    !Continental 

!
! ** specify input for calculating droplet number **********
	 IF (gNcn .ge. 450.0d0) THEN        ! Polluted Continental (Region 5)
             wparc    =    1.0d0
             wparcav  =    0.5d0
             sigwpar  =   0.25d0
             avgNcn   =  1700.d0
             sigNcn   =0.30103d0
             Wi       =dlog10(2.1d0)
             avgWi    =  0.225d0
             sigWi    =  0.075d0
             b0       =  647.8d0
             b1       =   55.0d0
             b2       =  371.0d0
             b3       =  -55.1d0
             b11      =  -23.7d0
             b22      =   90.0d0
             b33      =  -23.2d0
             b12      =  -25.5d0
             b13      =   -0.5d0
             b23      =  -24.8d0
	 ELSE IF (gNcn .lt. 450.0d0) THEN
              IF (gNcn .ge. 75.0d0)  THEN   ! Clean Continental (Region 6)
             wparc    =    1.0d0
             wparcav  =    0.5d0
             sigwpar  =   0.25d0
             avgNcn   =   400.d0
             sigNcn   =0.30103d0
             Wi       =dlog10(2.1d0)
             avgWi    =    0.3d0
             sigWi    =  0.075d0
             b0       =  298.8d0
             b1       =   26.5d0
             b2       =  189.6d0
             b3       =  -35.3d0
             b11      =   -5.1d0
             b22      =   49.8d0
             b33      =    6.9d0
             b12      =   19.2d0
             b13      =   -3.2d0
             b23      =  -20.4d0
	      ENDIF
	 ENDIF
!
         AWcb	= (wparc-wparcav)/sigwpar
         ANcn	= (dlog10(gNcn)-dlog10(avgNcn))/sigNcn
         AWi	= (Wi-avgWi)/sigWi
!
         dropSK =  b0 + b1*AWcb + b2*ANcn + b3*AWi+    &
     &             b11*AWcb*AWcb + b22*ANcn*ANcn + b33*AWi*AWi + &
     &             b12*AWcb*ANcn + b13*AWcb*AWi  + b23*ANcn*AWi

      ENDIF
!
! *** calculations complete ***********************************************
!
1000 continue
      return
      end subroutine paramSK 
!=======================================================================

      end module GmiCloudPropertiesGT_mod
