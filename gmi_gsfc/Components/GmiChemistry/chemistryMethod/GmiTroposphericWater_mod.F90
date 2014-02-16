!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiTroposphericWater_mod
!
! !INTERFACE:
!
      module GmiTroposphericWater_mod
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: addTroposphericWater, removeTroposphericWater
!
#     include "gmi_phys_constants.h"
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: addTroposphericWater
!
! !INTERFACE:
!
      subroutine addTroposphericWater  &
     &  (chem_mecha, ih2o_num, press3c, concentration, tropp, humidity,  &
     &   pchem_water, strat_water, pr_diag, loc_proc, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species)
!
      implicit none
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: chem_mecha
      logical         , intent(in) :: pr_diag
      integer         , intent(in) :: loc_proc
      integer         , intent(in) :: i1, i2, ju1, j2, k1, k2
      integer         , intent(in) :: ilo, ihi, julo, jhi
      integer         , intent(in) :: num_species
      ! species concentration array index for water
      integer         , intent(in) :: ih2o_num
      ! atmospheric pressure at the center of each grid box (mb)
      real*8          , intent(in) :: press3c(ilo:ihi, julo:jhi, k1:k2)
      ! tropopause pressure (mb)
      real*8          , intent(in) :: tropp  (i1:i2, ju1:j2)
      ! specific humidity   (g/kg)
      real*8          , intent(in) :: humidity   (i1:i2, ju1:j2, k1:k2)
!
! !OUTPUT PARAMETERS:
      ! stratospheric and tropospheric water vapor before chemistry
      ! operator (mixing ratio)
      real*8 , intent(out) :: pchem_water(i1:i2, ju1:j2, k1:k2)
      ! stratospheric water vapor before chemistry operator (mixing ratio)
      real*8 , intent(out) :: strat_water(i1:i2, ju1:j2, k1:k2)
!
! !INPUT/OUTPUT PARAMETERS:
      ! species concentration, known at zone centers (mixing ratio)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!
! !DESCRIPTION:
!   This routine adds the tropospheric water from the humidity variable into
!   the water vapor before doing chemistry; only do this below the calculated
!   tropopause pressure (tropp).  Make sure the water is never less than
!   2 ppm.
!
! !DEFINED PARAMETER:
      real*8, parameter :: MIN_STRAT_WATER = 3.0d-6
!
! !LOCAL VARIABLES:
      real*8  :: fac
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'addTroposphericWater called by ', loc_proc

      IF (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) THEN
         strat_water(:,:,:)    =  MIN_STRAT_WATER
         concentration(ih2o_num)%pArray3D(:,:,:) = strat_water(:,:,:)
      ELSE
         strat_water(:,:,:) = concentration(ih2o_num)%pArray3D(:,:,:)
      ENDIF

      fac = MWTAIR / (MWTH2O * GPKG)

      where (press3c(i1:i2,ju1:j2,:) > Spread (tropp(:,:), 3, k2))

        concentration(ih2o_num)%pArray3D(:,:,:) = &
     &      concentration(ih2o_num)%pArray3D(:,:,:) + (humidity(:,:,:) * fac)

      end where

      pchem_water(:,:,:) = concentration(ih2o_num)%pArray3D(:,:,:)

      where (concentration(ih2o_num)%pArray3D(:,:,:) < MIN_STRAT_WATER)
        concentration(ih2o_num)%pArray3D(:,:,:) = MIN_STRAT_WATER
      end where

      return

      end subroutine addTroposphericWater
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: removeTroposphericWater
!
! !INTERFACE:
!
      subroutine removeTroposphericWater  &
     &  (chem_mecha, ih2o_num, concentration, pchem_water, strat_water, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, num_species)
!
      implicit none
!
! !INPUT PARAMETER:
      character(len=*), intent(in) :: chem_mecha
      logical         , intent(in) :: pr_diag
      integer         , intent(in) :: loc_proc
      integer         , intent(in) :: i1, i2, ju1, j2, k1, k2
      integer         , intent(in) :: num_species
      ! species concentration array index for water
      integer         , intent(in) :: ih2o_num
      ! stratospheric and tropospheric water vapor before chemistry
      ! operator (mixing ratio)
      real*8          , intent(in) :: pchem_water(i1:i2, ju1:j2, k1:k2)
      ! stratospheric water vapor before chemistry operator (mixing ratio)
      real*8          , intent(in) :: strat_water(i1:i2, ju1:j2, k1:k2)
!
! !INPUT/OUTPUT PARAMETER:
      ! species concentration, known at zone centers (mixing ratio)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(num_species)
!
! !DESCRIPTION:
!  Removes the tropospheric water from the water vapor in the
!  species concentration array, now that chemistry is finished.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'removeTroposphericWater called by ', loc_proc

      IF (chem_mecha == 'strat_trop' .or. chem_mecha == 'strat_trop_aerosol' ) THEN
         concentration(ih2o_num)%pArray3D(:,:,:) = strat_water(:,:,:)
      ELSE
         concentration(ih2o_num)%pArray3D(:,:,:) =  &
     &        strat_water(:,:,:) *  &
     &        (1.0d0 + ((concentration(ih2o_num)%pArray3D(:,:,:) - &
     &                   pchem_water(:,:,:)) / pchem_water(:,:,:)))
      ENDIF

      return

      end subroutine removeTroposphericWater
!EOC
!------------------------------------------------------------------------------
      end module GmiTroposphericWater_mod
