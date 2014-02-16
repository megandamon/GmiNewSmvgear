  module GmiRelativeHumidity_mod

      implicit none
      
      private 
      public   :: CalcRelativeHumidity

 contains

!-----------------------------------------------------------------------------
!
! ROUTINE
!   CalcRelativeHumidity
!
! DESCRIPTION
!   This routine calculates relative humidity  using a formula from
!   Seinfeld (1986) (p. 181).
!
! ARGUMENTS
!   kel         : temperature (degK)
!   press3c     : atmospheric pressure at the center of each grid box (mb)
!   humidity    : specific humidity (g/kg)
!   relhum      : relative humidity (fraction)
!
!-----------------------------------------------------------------------------

      subroutine CalcRelativeHumidity  &
     &  (kel, press3c, humidity, relhum, &
     &   pr_diag, loc_proc, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)


      implicit none

#     include "gmi_phys_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in ) :: pr_diag
      integer, intent(in ) :: loc_proc
      integer, intent(in ) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      real*8 , intent(in ) :: kel     (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in ) :: press3c (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in ) :: humidity(i1:i2, ju1:j2, k1:k2)
      real*8 , intent(out) :: relhum  (i1:i2, ju1:j2, k1:k2)

!     ----------------------
!     Parameter declarations.
!     ----------------------


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Calc_Humidity_Rel called by ', loc_proc
      end if

!     -----------------------------------------------------------------
!     First calculate relative humidity from Seinfeld (1986) p. 181.
!     The first  relhume is the temperature dependent parameter a.
!     The second relhume is the saturation vapor pressure of water.
!     The third  relhume is the actual relative humidity as a fraction.
!     Then make sure relhume is between 0 and 0.95.
!     -----------------------------------------------------------------


      relhum(:,:,:) = 1.0d0 - (373.15d0 / kel(i1:i2,ju1:j2,:))

      relhum(:,:,:) =  &
     &  1013.25d0 * Exp (13.3185d0 * relhum(:,:,:)    -  &
     &                    1.9760d0 * relhum(:,:,:)**2 -  &
     &                    0.6445d0 * relhum(:,:,:)**3 -  &
     &                    0.1299d0 * relhum(:,:,:)**4)

      relhum(:,:,:) =  &
     &  humidity(:,:,:) * MWTAIR / 18.0d0 /  &
     &  GPKG * press3c(i1:i2,ju1:j2,:) / relhum(:,:,:)

!      relhum(:,:,:) = Max (Min (relhum(:,:,:), 1.0d0), 1.0d-30)
      relhum(:,:,:) = Max (Min (relhum(:,:,:), 0.95d0), 0.0d0)

      return

      end subroutine CalcRelativeHumidity

      end module GmiRelativeHumidity_mod
